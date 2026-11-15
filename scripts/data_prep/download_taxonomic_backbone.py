"""Download CoL/POWO/BHL/NatureServe taxonomic and subnational caches into flat `data/taxonomy/` CSVs and SQLite Tier A."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union
from urllib import error as _uerr, parse as _uparse

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import load_dotenv_defaults, resolve_target_genus
from scripts.core.db import (
    get_db_connection, is_canonical_or_explicit_db_context, resolve_db_path,
    upsert_natureserve_subnational_ranges, upsert_taxon_backbone_profiles,
)
from scripts.core.epa_ecoregions import _CANADIAN_PROV_ABBR_TO_FULL, _US_STATE_ABBR_TO_FULL, _normalize_us_state_name
from scripts.core.io_utils import atomic_write_csv_rows, atomic_write_text, default_http_get, to_text

PROFILE_CSV_NAME = "public_taxon_backbone_profiles.csv"
RANGE_CSV_NAME = "public_natureserve_subnational_ranges.csv"
MANIFEST_JSON_NAME = "public_taxonomic_backbone_manifest.json"
COL_WCVP_DATASET_ID = 2232
_USER_AGENT = "Plantmatics-TaxonomicBackbone/1.0 (UNC-Chapel Hill Herbarium)"

PROFILE_FIELDS: Tuple[str, ...] = (
    "species_standardized", "target_genus", "accepted_scientific_name", "accepted_authorship", "basionym_name",
    "basionym_authorship", "homotypic_synonyms_json", "heterotypic_synonyms_json", "is_hybrid_taxon", "hybrid_formula",
    "publication_citation", "publication_year", "protologue_url", "locus_classicus_text", "tdwg_native_codes_json",
    "tdwg_introduced_codes_json", "natureserve_element_id", "natureserve_global_rank", "nwpl_indicator_code",
    "nwpl_wetland_score", "backbone_source", "updated_at",
)
RANGE_FIELDS: Tuple[str, ...] = ("species_standardized", "country_code", "subnation_code", "s_rank", "is_native", "is_exotic", "is_historical_or_extirpated")
_US_ABBR_SET: Set[str] = {k.upper() for k in _US_STATE_ABBR_TO_FULL}
_US_FULL_TO_ABBR: Dict[str, str] = {v: k.upper() for k, v in _US_STATE_ABBR_TO_FULL.items()}
_CA_FULL_TO_ABBR: Dict[str, str] = {
    "ontario": "ON", "quebec": "QC", "new brunswick": "NB", "nova scotia": "NS", "prince edward island": "PE",
    "newfoundland and labrador": "NL", "manitoba": "MB", "saskatchewan": "SK", "alberta": "AB",
    "british columbia": "BC", "yukon": "YT", "northwest territories": "NT", "nunavut": "NU",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _validate_flat_taxonomy_dir(tax_dir: Path) -> None:
    if tax_dir.exists() and not tax_dir.is_dir():
        raise ValueError(f"Taxonomy path is not a directory: {tax_dir}")
    subdirs = [p.name for p in tax_dir.iterdir() if p.is_dir()] if tax_dir.exists() else []
    if subdirs:
        raise ValueError(f"Flat staging invariant violated under {tax_dir}; found subdirs: {subdirs}")
    tax_dir.mkdir(parents=True, exist_ok=True)


def _http_request_json(url: str, *, http_get_fn: Optional[Callable[[str], Any]] = None, post_body: Optional[Mapping[str, Any]] = None, timeout: float = 20.0) -> Optional[Any]:
    try:
        if http_get_fn is not None:
            raw = http_get_fn(url)
            if isinstance(raw, (dict, list)):
                return raw
            txt = to_text(raw) if isinstance(raw, (bytes, bytearray, str)) else ""
            return json.loads(txt) if txt.strip() else None
        headers = {"Accept": "application/json"}
        data = json.dumps(post_body).encode("utf-8") if post_body is not None else None
        if data is not None:
            headers["Content-Type"] = "application/json"
        txt = to_text(default_http_get(url, timeout=timeout, user_agent=_USER_AGENT, headers=headers, data=data)).strip()
        return json.loads(txt) if txt else None
    except (_uerr.URLError, TimeoutError, json.JSONDecodeError, OSError, ValueError) as exc:
        print(f"[WARN] HTTP request failed for {url}: {exc}")
        return None


def _parse_species_name(raw_name: Any, target_genus: str) -> Tuple[Optional[str], str, bool]:
    """Normalize a botanical binomial into `(species_standardized, display_name, is_hybrid)`."""
    if not raw_name:
        return None, "", False
    txt = re.sub(r"\s+", " ", str(raw_name).strip())
    if re.search(r"\b(?:var|subsp|ssp|f|cv)\.?\s+", txt, flags=re.IGNORECASE):
        return None, "", False
    m = re.match(r"^([A-Za-z-]+)(?:\s*(×)\s*|\s+[xX]\s+|\s+)([a-z][a-z-]+)\b", txt)
    if not m:
        return None, "", False
    genus, epithet = m.group(1).capitalize(), m.group(3).lower()
    if genus.lower() != target_genus.strip().lower() or epithet in {"sp", "spp", "cf", "aff", "x"}:
        return None, "", False
    is_hybrid = bool(m.group(2) or re.match(r"^[A-Za-z-]+\s+[xX×]\s*[a-z]", txt))
    sp_std = f"{genus} {epithet}"
    return sp_std, (f"{genus} × {epithet}" if is_hybrid else sp_std), is_hybrid


def _year_from_text(val: Any) -> Optional[int]:
    if isinstance(val, int) and 1700 <= val <= 2100:
        return val
    m = re.search(r"\b(1[789]\d{2}|20\d{2})\b", str(val or ""))
    return int(m.group(1)) if m else None


def _syn_label(syn: Any) -> str:
    if not isinstance(syn, Mapping):
        return str(syn or "").strip()
    if syn.get("label"):
        return str(syn["label"]).strip()
    nm = syn.get("name")
    n_str = nm.get("scientificName") if isinstance(nm, Mapping) else nm
    auth = syn.get("author") or syn.get("authorship") or (nm.get("authorship") if isinstance(nm, Mapping) else "")
    return f"{n_str or ''} {auth or ''}".strip()


def _new_profile(sp_std: str, target_genus: str, source: str) -> Dict[str, Any]:
    return {
        **dict.fromkeys(PROFILE_FIELDS, None), "species_standardized": sp_std, "target_genus": target_genus,
        "accepted_scientific_name": sp_std, "is_hybrid_taxon": 0, "backbone_source": source, "updated_at": _now_iso(),
        "_homotypic": set(), "_heterotypic": set(), "_tdwg_native": set(), "_tdwg_intro": set(),
    }


def _ingest_distributions(prof: Dict[str, Any], dist_obj: Any) -> None:
    if isinstance(dist_obj, dict):
        for bucket, key in (("natives", "_tdwg_native"), ("introduced", "_tdwg_intro")):
            for item in dist_obj.get(bucket) or []:
                if code := str((item or {}).get("tdwgCode") or (item or {}).get("area") or "").strip().upper():
                    prof[key].add(code)
    elif isinstance(dist_obj, list):
        for item in dist_obj:
            if isinstance(item, Mapping) and (code := str(item.get("tdwgCode") or item.get("areaId") or item.get("area") or "").strip().upper()):
                status = str(item.get("status") or item.get("establishmentMeans") or "").lower()
                prof["_tdwg_intro" if any(k in status for k in ("introduc", "exotic", "alien")) else "_tdwg_native"].add(code)


def _apply_basionym(prof: Dict[str, Any], bas: Any, fallback_auth: Any = None) -> None:
    if isinstance(bas, Mapping):
        prof["basionym_name"] = str(bas.get("scientificName") or bas.get("name") or "").strip() or prof["basionym_name"]
        prof["basionym_authorship"] = str(bas.get("authorship") or "").strip() or prof["basionym_authorship"]
    elif isinstance(bas, str) and bas.strip():
        prof["basionym_name"] = bas.strip()
    if fallback_auth and not prof["basionym_authorship"]:
        prof["basionym_authorship"] = str(fallback_auth).strip() or None


def _parse_col_entries(entries: Sequence[Mapping[str, Any]], target_genus: str, http_get_fn: Optional[Callable[[str], Any]] = None) -> Dict[str, Dict[str, Any]]:
    profiles: Dict[str, Dict[str, Any]] = {}
    syn_queue: List[Tuple[str, str, bool]] = []
    probe_info = True
    for item in entries:
        usage = item.get("usage") if isinstance(item.get("usage"), Mapping) else item
        name_obj = usage.get("name") if isinstance(usage.get("name"), Mapping) else usage
        sci = name_obj.get("scientificName") or name_obj.get("canonicalName") or item.get("scientificName") or ""
        rank = str(name_obj.get("rank") or item.get("rank") or "").lower()
        status = str(usage.get("status") or item.get("status") or "accepted").lower()
        auth = str(name_obj.get("authorship") or item.get("authorship") or "").strip() or None
        acc_raw = usage.get("accepted") if usage.get("accepted") is not None else item.get("accepted")
        if "synonym" in status or isinstance(acc_raw, (Mapping, str)):
            acc_obj = acc_raw if isinstance(acc_raw, Mapping) else {}
            acc_nm = acc_obj.get("name") if isinstance(acc_obj.get("name"), Mapping) else acc_obj
            acc_sci = acc_raw if isinstance(acc_raw, str) else (acc_nm.get("scientificName") or acc_obj.get("label") or "")
            acc_std, _, _ = _parse_species_name(acc_sci, target_genus)
            if acc_std and (full_lbl := f"{sci} {auth}".strip() if auth else str(sci).strip()):
                syn_queue.append((acc_std, full_lbl, ("homotypic" in status) or bool(usage.get("homotypic") or item.get("homotypic"))))
            continue
        if rank and rank not in {"species", "nothospecies"}:
            continue
        sp_std, disp_name, is_hyb = _parse_species_name(sci, target_genus)
        if not sp_std:
            continue
        prof = profiles.setdefault(sp_std, _new_profile(sp_std, target_genus, "CoL_WCVP_2232"))
        prof["accepted_scientific_name"] = disp_name
        prof["accepted_authorship"] = auth or prof["accepted_authorship"]
        prof["is_hybrid_taxon"] = 1 if (is_hyb or rank == "nothospecies" or name_obj.get("notho")) else prof["is_hybrid_taxon"]
        if prof["is_hybrid_taxon"] and not prof["hybrid_formula"]:
            prof["hybrid_formula"] = str(item.get("hybridFormula") or disp_name)
        _apply_basionym(prof, item.get("basionym") or usage.get("basionym") or name_obj.get("basionym"), name_obj.get("basionymAuthorship"))
        if pub_cit := str(item.get("publicationCitation") or name_obj.get("publishedIn") or "").strip():
            prof["publication_citation"] = pub_cit
        prof["publication_year"] = prof["publication_year"] or _year_from_text(pub_cit or name_obj.get("publishedInYear") or item.get("publicationYear"))
        prof["protologue_url"] = str(item.get("protologueUrl") or prof["protologue_url"] or "").strip() or None
        prof["locus_classicus_text"] = str(item.get("locusClassicus") or prof["locus_classicus_text"] or "").strip() or None
        for bucket, key in (("_homotypic", "homotypicSynonyms"), ("_heterotypic", "heterotypicSynonyms")):
            prof[bucket].update(lbl for lbl in (_syn_label(s) for s in item.get(key) or []) if lbl)
        _ingest_distributions(prof, item.get("distributions"))
        if probe_info and (uid := usage.get("id")) and not prof["_homotypic"] and not prof["_heterotypic"]:
            info = _http_request_json(f"https://api.checklistbank.org/dataset/{COL_WCVP_DATASET_ID}/taxon/{_uparse.quote(str(uid))}/info", http_get_fn=http_get_fn)
            if not isinstance(info, Mapping) or "result" in info:
                probe_info = False
            else:
                syns = info.get("synonyms") or {}
                if isinstance(syns, Mapping):
                    for bucket, key in (("_homotypic", "homotypic"), ("_heterotypic", "heterotypic")):
                        prof[bucket].update(lbl for lbl in (_syn_label(s) for s in syns.get(key) or []) if lbl)
                _ingest_distributions(prof, info.get("distributions"))
    for acc_std, syn_lbl, is_hom in syn_queue:
        if acc_std in profiles:
            profiles[acc_std]["_homotypic" if is_hom else "_heterotypic"].add(syn_lbl)
    return profiles


def _fetch_powo_fallback(target_genus: str, http_get_fn: Optional[Callable[[str], Any]] = None) -> Dict[str, Dict[str, Any]]:
    payload = _http_request_json(f"https://powo.science.kew.org/api/2/search?{_uparse.urlencode({'q': f'genus:{target_genus}'})}", http_get_fn=http_get_fn)
    if not isinstance(payload, Mapping):
        return {}
    profiles: Dict[str, Dict[str, Any]] = {}
    probe_detail = True
    for item in payload.get("results") or []:
        if not isinstance(item, Mapping) or item.get("accepted") is False:
            continue
        sp_std, disp_name, is_hyb = _parse_species_name(item.get("name"), target_genus)
        if not sp_std:
            continue
        prof = profiles.setdefault(sp_std, _new_profile(sp_std, target_genus, "Kew_POWO_v2"))
        prof["accepted_scientific_name"] = disp_name
        prof["accepted_authorship"] = str(item.get("author") or "").strip() or None
        prof["is_hybrid_taxon"] = 1 if (is_hyb or item.get("hybrid")) else 0
        detail = item
        if probe_detail and item.get("fqId") and "synonyms" not in item and "distribution" not in item:
            resp = _http_request_json(f"https://powo.science.kew.org/api/2/taxon/{_uparse.quote(str(item['fqId']))}", http_get_fn=http_get_fn)
            if isinstance(resp, Mapping) and "results" not in resp:
                detail = resp
            else:
                probe_detail = False
        _apply_basionym(prof, detail.get("basionym"))
        prof["publication_citation"] = str(detail.get("reference") or "").strip() or prof["publication_citation"]
        prof["publication_year"] = _year_from_text(detail.get("publicationYear") or prof["publication_citation"])
        prof["hybrid_formula"] = str(detail.get("hybridFormula") or "").strip() or (disp_name if prof["is_hybrid_taxon"] else None)
        for syn in detail.get("synonyms") or []:
            if isinstance(syn, Mapping) and (lbl := _syn_label(syn)):
                prof["_homotypic" if (bool(syn.get("homotypic")) or str(syn.get("taxonomicStatus") or "").lower().startswith("hom")) else "_heterotypic"].add(lbl)
        _ingest_distributions(prof, detail.get("distribution"))
    return profiles


def fetch_col_powo_backbone(target_genus: str, http_get_fn: Optional[Callable[[str], Any]] = None) -> Dict[str, Dict[str, Any]]:
    """Query CoL WCVP dataset 2232 with Kew POWO v2 fallback for accepted species, basionyms, synonyms, and TDWG codes."""
    all_entries: List[Mapping[str, Any]] = []
    for page in range(5):
        url = f"https://api.checklistbank.org/dataset/{COL_WCVP_DATASET_ID}/nameusage/search?" + _uparse.urlencode({"q": target_genus, "content": "SCIENTIFIC_NAME", "limit": 1000, "offset": page * 1000})
        payload = _http_request_json(url, http_get_fn=http_get_fn)
        if not isinstance(payload, Mapping) or not payload.get("result"):
            break
        batch = [r for r in payload.get("result") or [] if isinstance(r, Mapping)]
        all_entries.extend(batch)
        if payload.get("last", True) or len(batch) < 1000:
            break
    raw_profiles = (_parse_col_entries(all_entries, target_genus, http_get_fn=http_get_fn) if all_entries else {}) or _fetch_powo_fallback(target_genus, http_get_fn=http_get_fn)
    return {
        sp: {
            **{k: v for k, v in rec.items() if not k.startswith("_")},
            "homotypic_synonyms_json": json.dumps(sorted(rec["_homotypic"]), ensure_ascii=False),
            "heterotypic_synonyms_json": json.dumps(sorted(rec["_heterotypic"]), ensure_ascii=False),
            "tdwg_native_codes_json": json.dumps(sorted(rec["_tdwg_native"]), ensure_ascii=False),
            "tdwg_introduced_codes_json": json.dumps(sorted(rec["_tdwg_intro"]), ensure_ascii=False),
        } for sp, rec in sorted(raw_profiles.items())
    }


def enrich_bhl_protologues(profiles: Dict[str, Dict[str, Any]], bhl_api_key: Optional[str] = None, http_get_fn: Optional[Callable[[str], Any]] = None) -> None:
    """Enrich taxon profiles with Biodiversity Heritage Library protologue URLs and locus/year metadata."""
    load_dotenv_defaults()
    eff_key = (bhl_api_key if bhl_api_key is not None else os.environ.get("BHL_API_KEY", "")).strip()
    if not eff_key:
        print("[INFO] BHL_API_KEY not configured; skipping BHL protologue URL enrichment.")
        return
    for sp, prof in profiles.items():
        if prof.get("protologue_url") or not (q_name := str(prof.get("basionym_name") or prof.get("accepted_scientific_name") or sp).strip()):
            continue
        url = "https://www.biodiversitylibrary.org/api3?" + _uparse.urlencode({"op": "NameGetDetail", "name": q_name, "apikey": eff_key, "format": "json"})
        payload = _http_request_json(url, http_get_fn=http_get_fn)
        if not isinstance(payload, Mapping) or str(payload.get("Status") or "").lower() == "error":
            continue
        res = payload.get("Result")
        titles = (res.get("Titles") if isinstance(res, Mapping) else res) or []
        if not isinstance(titles, list) or not titles or not isinstance(titles[0], Mapping):
            continue
        t0, items = titles[0], titles[0].get("Items") or []
        p_url = t0.get("PageUrl") or t0.get("TitleUrl") or (res.get("ProtologueUrl") if isinstance(res, Mapping) else None)
        if not p_url and isinstance(items, list) and items and isinstance(items[0], Mapping):
            pages = items[0].get("Pages") or []
            p_url = (pages[0].get("PageUrl") if pages and isinstance(pages[0], Mapping) else None) or items[0].get("ItemUrl")
        if p_url:
            prof["protologue_url"] = str(p_url).strip()
        prof["publication_year"] = prof.get("publication_year") or _year_from_text(t0.get("PublicationDate"))
        if not prof.get("locus_classicus_text") and t0.get("PublicationPlace"):
            prof["locus_classicus_text"] = str(t0["PublicationPlace"]).strip()


def _normalize_subnation(country_raw: Any, sub_raw: Any) -> Tuple[Optional[str], Optional[str]]:
    c_tok = str(country_raw or "").strip().upper()
    c_code = "US" if c_tok in {"US", "USA", "UNITED STATES"} else ("CA" if c_tok in {"CA", "CAN", "CANADA"} else None)
    s_raw = str(sub_raw or "").strip()
    if re.match(r"^(?:US|CA)-([A-Za-z]{2})$", s_raw, flags=re.IGNORECASE):
        cp, sc = s_raw.split("-", 1)
        return (c_code or cp.upper()), sc.upper()
    s_up = s_raw.upper()
    if c_code in {None, "US"} and s_up in _US_ABBR_SET:
        return "US", s_up
    norm = _normalize_us_state_name(s_raw) or ""
    full = _CANADIAN_PROV_ABBR_TO_FULL.get(norm, norm)
    if full in _CA_FULL_TO_ABBR and c_code in {None, "CA"}:
        return "CA", _CA_FULL_TO_ABBR[full]
    if full in _US_FULL_TO_ABBR and c_code in {None, "US"}:
        return "US", _US_FULL_TO_ABBR[full]
    return ("CA", s_up) if (c_code == "CA" and s_up in _CA_FULL_TO_ABBR.values()) else (None, None)


def fetch_natureserve_ranges(species_list: Sequence[str], target_genus: str, http_get_fn: Optional[Callable[[str], Any]] = None) -> Tuple[Dict[str, Dict[str, str]], List[Dict[str, Any]]]:
    """Query NatureServe Explorer API for global G-ranks and US/CA subnational distributions."""
    ranks_by_taxon: Dict[str, Dict[str, str]] = {}
    range_map: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
    target_set = {sp for sp in species_list if sp}
    ns_body = {"criteriaType": "species", "textCriteria": [{"paramType": "quickSearch", "searchToken": target_genus, "matchAgainst": "allNames", "operator": "contains"}], "pagingOptions": {"page": 0, "recordsPerPage": 250}}
    search_payload = (
        _http_request_json("https://explorer.natureserve.org/api/data/speciesSearch?" + _uparse.urlencode({"text": target_genus}), http_get_fn=http_get_fn)
        if http_get_fn is not None else _http_request_json("https://explorer.natureserve.org/api/data/speciesSearch", post_body=ns_body)
    )
    probe_taxon = True
    for item in (list(search_payload.get("results") or []) if isinstance(search_payload, Mapping) else []):
        if not isinstance(item, Mapping):
            continue
        sp_std, _, _ = _parse_species_name(item.get("scientificName") or item.get("primaryCommonName") or "", target_genus)
        if not sp_std or (target_set and sp_std not in target_set):
            continue
        uid, g_rank = str(item.get("uniqueId") or item.get("elementGlobalId") or "").strip(), str(item.get("roundedGRank") or item.get("gRank") or "").strip()
        detail = item
        if probe_taxon and uid and "nations" not in item:
            resp = _http_request_json(f"https://explorer.natureserve.org/api/data/taxon/{_uparse.quote(uid)}", http_get_fn=http_get_fn)
            if isinstance(resp, Mapping) and ("nations" in resp or "roundedGRank" in resp or "gRank" in resp):
                detail, g_rank = resp, str(resp.get("roundedGRank") or resp.get("gRank") or g_rank).strip()
            else:
                probe_taxon = False
        if uid or g_rank:
            ranks_by_taxon[sp_std] = {"natureserve_element_id": uid, "natureserve_global_rank": g_rank}
        for nat in (n for n in (detail.get("nations") or []) if isinstance(n, Mapping)):
            n_code = nat.get("nationCode") or nat.get("isoCode") or nat.get("countryCode")
            for sub in (s for s in (nat.get("subnations") or []) if isinstance(s, Mapping)):
                c_code, sub_code = _normalize_subnation(n_code, sub.get("subnationCode") or sub.get("name"))
                if not c_code or not sub_code:
                    continue
                s_rank = str(sub.get("sRank") or sub.get("roundedSRank") or "").strip()
                exotic = bool(sub.get("exotic")) or ("exotic" in str(sub.get("origin") or "").lower()) or s_rank.upper().startswith("SE")
                range_map[(sp_std, c_code, sub_code)] = {
                    "species_standardized": sp_std, "country_code": c_code, "subnation_code": sub_code, "s_rank": s_rank,
                    "is_native": 1 if ((not exotic) if sub.get("native") is None else bool(sub.get("native"))) else 0,
                    "is_exotic": int(exotic), "is_historical_or_extirpated": int(bool(sub.get("extirpated")) or s_rank.upper().startswith(("SH", "SX"))),
                }
    return ranks_by_taxon, [range_map[k] for k in sorted(range_map)]


def _read_csv_safe(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def sync_backbone_csvs_to_db(taxonomy_dir: Path, db_path: Union[str, Path], demo: bool = False) -> Dict[str, int]:
    """Read flat `data/taxonomy/` CSVs and atomically upsert them into SQLite Tier A tables."""
    _validate_flat_taxonomy_dir(taxonomy_dir)
    if demo or os.environ.get("PLANTMATICS_DEMO_MODE") == "1" or os.environ.get("PACKERA_DEMO_MODE") == "1":
        return {"taxon_backbone_profiles": 0, "natureserve_subnational_ranges": 0}
    prof_rows, range_rows = _read_csv_safe(taxonomy_dir / PROFILE_CSV_NAME), _read_csv_safe(taxonomy_dir / RANGE_CSV_NAME)
    conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=PROJECT_ROOT)
    try:
        return {
            "taxon_backbone_profiles": upsert_taxon_backbone_profiles(conn, prof_rows, demo=False, preserve_existing_on_null=True),
            "natureserve_subnational_ranges": upsert_natureserve_subnational_ranges(conn, range_rows, demo=False, preserve_existing_on_null=False),
        }
    finally:
        conn.close()


def download_taxonomic_backbone(
    taxonomy_dir: Union[str, Path] = PROJECT_ROOT / "data" / "taxonomy", target_genus: Optional[str] = None, db_path: Optional[Union[str, Path]] = None,
    *, sync_only: bool = False, force: bool = False, bhl_api_key: Optional[str] = None,
    col_http_fn: Optional[Callable[[str], Any]] = None, bhl_http_fn: Optional[Callable[[str], Any]] = None, natureserve_http_fn: Optional[Callable[[str], Any]] = None,
) -> Dict[str, Any]:
    """Download or reload CoL/POWO/BHL/NatureServe taxonomy caches and sync them to SQLite."""
    eff_dir, eff_genus = Path(taxonomy_dir), resolve_target_genus(explicit_genus=target_genus)
    _validate_flat_taxonomy_dir(eff_dir)
    eff_db = resolve_db_path(db_path=db_path, workspace_root=PROJECT_ROOT)
    prof_csv, range_csv, manifest_json = eff_dir / PROFILE_CSV_NAME, eff_dir / RANGE_CSV_NAME, eff_dir / MANIFEST_JSON_NAME
    if not sync_only and (force or not prof_csv.is_file()):
        profiles = fetch_col_powo_backbone(eff_genus, http_get_fn=col_http_fn)
        enrich_bhl_protologues(profiles, bhl_api_key=bhl_api_key, http_get_fn=bhl_http_fn)
        ranks, ns_rows = fetch_natureserve_ranges(sorted(profiles), eff_genus, http_get_fn=natureserve_http_fn)
        prof_list = [{**profiles[sp], **{k: v for k, v in ranks.get(sp, {}).items() if v}} for sp in sorted(profiles)]
        atomic_write_csv_rows(prof_csv, PROFILE_FIELDS, prof_list)
        atomic_write_csv_rows(range_csv, RANGE_FIELDS, ns_rows)
        atomic_write_text(manifest_json, json.dumps({"target_genus": eff_genus, "generated_at": _now_iso(), "n_profiles": len(prof_list), "n_natureserve_ranges": len(ns_rows), "files": [PROFILE_CSV_NAME, RANGE_CSV_NAME]}, indent=2) + "\n")
    sync_ok = db_path is not None or bool(os.environ.get("PLANTMATICS_DB_PATH") or os.environ.get("PACKERA_DB_PATH")) or is_canonical_or_explicit_db_context(eff_dir, workspace_root=PROJECT_ROOT)
    if not sync_ok:
        print("[WARN] Standalone taxonomy dir without explicit --db; skipping canonical DB sync.")
    counts = sync_backbone_csvs_to_db(eff_dir, eff_db) if sync_ok else {"taxon_backbone_profiles": 0, "natureserve_subnational_ranges": 0}
    return {"taxonomy_dir": str(eff_dir), "target_genus": eff_genus, "db_path": str(eff_db), **counts}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Download CoL/POWO/BHL/NatureServe taxonomic backbone into SQLite.")
    parser.add_argument("--taxonomy-dir", type=Path, default=PROJECT_ROOT / "data" / "taxonomy")
    parser.add_argument("--genus", type=str, default=None)
    parser.add_argument("--db", type=Path, default=None)
    parser.add_argument("--sync-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(download_taxonomic_backbone(taxonomy_dir=args.taxonomy_dir, target_genus=args.genus, db_path=args.db, sync_only=args.sync_only, force=args.force), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
