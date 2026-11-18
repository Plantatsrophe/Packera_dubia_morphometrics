"""Non-circular 4-layer taxonomic specialist discovery (OpenAlex, WCVP/Wikidata P428, Bionomia, overrides)."""

from __future__ import annotations

import csv, json, os, re, unicodedata
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple
from urllib import error as _uerr, parse as _uparse

from scripts.core.config import SPECIALIST_PATTERNS, resolve_target_family, resolve_target_genus
from scripts.core.db import get_db_connection, is_canonical_or_explicit_db_context, load_taxonomic_specialists, resolve_db_path, upsert_taxonomic_specialists
from scripts.core.io_utils import atomic_write_csv_rows, default_http_get, to_text

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SPECIALISTS_CSV_NAME = "public_taxonomic_specialists.csv"
SPECIALIST_FIELDS: Tuple[str, ...] = (
    "specialist_id", "scope_type", "scope_taxon", "canonical_name", "qualified_regex", "ipni_abbr", "orcid",
    "wikidata_qid", "openalex_works_count", "bionomia_det_count", "tier_cap", "source", "user_verified", "updated_at",
)
_AMBIGUOUS_COMMON_SURNAMES: Set[str] = set("anderson atwood baker bell britton brooks brown bush campbell carter chapman clark collins cook cooper cox davis diamond dorn duncan estes evans fisher freeman fuller golden gray green hall heller hill horn howard jackson james jeffrey johnson jones kelly king lee lewis long love martin miller mitchell moore morgan morris murphy nelson packer parker perry phillips powell price reed richardson roberts robinson rogers ross russell scott smith stewart taylor thomas thompson turner walker ward watson weber welsh white williams wilson wood wright young".split())
_TAXONOMIC_CONCEPTS: Tuple[str, ...] = (
    "taxonomy", "systematics", "phylogen", "revision", "monograph", "nomenclat", "nov.", "sp. nov", "comb. nov",
    "phytokeys", "taxon", "systematic botany", "brittonia", "novon", "phytoneuron", "sida", "rhodora",
    "harvard pap", "blumea", "kew bull", "ann. missouri", "nordic j. bot", "telopea", "phytotaxa", "madrono", "madroño", "castanea",
)
_INVALID_INSTS: Set[str] = {"", "UNKNOWN", "UNKNOWN_INST", "NONE", "NAN", "NULL", "NA"}
_PATTERN_CACHE: Dict[Tuple[Any, ...], Tuple[List[str], List[str]]] = {}


def _ascii_fold(text: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFKD", str(text or "")) if not unicodedata.combining(c)).lower()

def _to_int_safe(val: Any) -> int:
    try: return int(float(str(val or 0).strip() or 0))
    except (TypeError, ValueError): return 0

def _norm_orcid(val: Any) -> Optional[str]:
    s = re.sub(r"^https?://orcid\.org/", "", str(val or "").strip(), flags=re.IGNORECASE).strip()
    return s if re.match(r"^\d{4}-\d{4}-\d{4}-[\dX]{4}$", s, flags=re.IGNORECASE) else None

def _year_from_val(val: Any) -> Optional[int]:
    if isinstance(val, int) and 1600 <= val <= 2100: return val
    m = re.search(r"\b(1[6789]\d{2}|20\d{2})\b", str(val or ""))
    return int(m.group(1)) if m else None

def _http_json(url: str, http_fn: Optional[Callable[[str], Any]] = None, timeout: float = 12.0) -> Optional[Any]:
    try:
        raw = http_fn(url) if http_fn is not None else default_http_get(url, timeout=timeout, headers={"Accept": "application/json"})
        if isinstance(raw, (dict, list)): return raw
        txt = to_text(raw).strip() if isinstance(raw, (bytes, bytearray, str)) else ""
        return json.loads(txt) if txt else None
    except (_uerr.URLError, TimeoutError, json.JSONDecodeError, OSError, ValueError):
        return None

def _read_csv_rows(csv_path: Path) -> List[Dict[str, Any]]:
    if not csv_path.is_file() or csv_path.stat().st_size == 0: return []
    with open(csv_path, "r", encoding="utf-8", newline="") as fh: return [dict(r) for r in csv.DictReader(fh)]

def _parse_name_parts(raw_name: str) -> Tuple[List[str], str]:
    s = re.sub(r"\bet\s+al\.?|\b(?:dr|prof|mr|mrs|ms|jr|sr|ii|iii|iv)\.?\b", "", str(raw_name or ""), flags=re.IGNORECASE).strip(" ,;.")
    if not s: return [], ""
    if "," in s:
        sp, gp = [p.strip() for p in s.split(",", 1)]
        st = [re.sub(r"[^A-Za-zÀ-ÿ'’-]", "", t) for t in sp.split()]
        gt = [re.sub(r"[^A-Za-zÀ-ÿ]", "", t) for t in re.sub(r"\.", ". ", gp).replace(".", " ").split()]
        return [t for t in gt if t], (st[-1] if st else "")
    toks = [t for t in (re.sub(r"[^A-Za-zÀ-ÿ'’-]", "", x.replace(".", "")) for x in re.sub(r"\.(?=[A-Za-zÀ-ÿ])", ". ", s).split()) if t]
    return (toks[:-1], toks[-1]) if len(toks) >= 2 else ([], toks[0] if toks else "")

def _build_qualified_regex(canonical_name: str, ipni_abbr: Optional[str] = None, allow_unambiguous_bare: bool = False) -> Optional[str]:
    giv, surname = _parse_name_parts(canonical_name)
    sur_low = _ascii_fold(surname).strip("-'’")
    if len(sur_low) < 2 or not re.search(r"[a-z]{2,}", sur_low): return None
    sur_vars = {re.escape(surname), re.escape(sur_low)}
    sur_rx = f"(?:{'|'.join(sorted(sur_vars))})" if len(sur_vars) > 1 else re.escape(surname)
    branches: List[str] = []
    if giv:
        first, f0, f0a = giv[0], giv[0][0], _ascii_fold(giv[0][0]).upper()
        fi = f"[{re.escape(f0)}{re.escape(f0a)}]" if f0.lower() != f0a.lower() else re.escape(f0)
        mid_init = "".join(rf"(?:{re.escape(m[0])}\.?\s*)?" for m in giv[1:])
        branches.append(rf"\b{fi}\.?\s*{mid_init}{sur_rx}\b")
        if len(first) > 1:
            mid_full = "".join(rf"(?:{re.escape(m)}\.?\s*|{re.escape(m[0])}\.?\s*)?" for m in giv[1:])
            branches.extend([rf"\b{re.escape(first)}\s+{mid_full}{sur_rx}\b", rf"\b{sur_rx}\s*,\s*(?:{re.escape(first)}|{fi})\b"])
        else: branches.append(rf"\b{sur_rx}\s*,\s*{fi}\b")
    elif allow_unambiguous_bare and sur_low not in _AMBIGUOUS_COMMON_SURNAMES and len(sur_low) >= 4:
        branches.append(rf"\b{sur_rx}\b")
    else: return None
    if ipni_abbr and (ig := _parse_name_parts(ipni_abbr))[0] and _ascii_fold(ig[1]) == sur_low:
        branches.append(rf"\b{re.escape(ig[0][0][0])}\.?\s*{sur_rx}\b")
    pat = "|".join(dict.fromkeys(branches))
    try: rx = re.compile(pat, re.IGNORECASE)
    except re.error: return None
    return None if ((giv or sur_low in _AMBIGUOUS_COMMON_SURNAMES) and rx.search(surname)) else pat

def _author_key(canonical_name: str, orcid: Optional[str] = None) -> str:
    if orcid:
        return f"orcid:{orcid.lower()}"
    giv, sur = _parse_name_parts(canonical_name)
    return f"name:{_ascii_fold(sur)}:{_ascii_fold(giv[0])[0] if giv else ''}"

def _new_rec(scope_taxon: str, key: str, cname: str, qreg: str, source: str, *, tier_cap: str = "Tier_1_Gold",
             ipni_abbr: Optional[str] = None, orcid: Optional[str] = None, wikidata_qid: Optional[str] = None,
             works: int = 0, dets: int = 0) -> Dict[str, Any]:
    sid = re.sub(r"[^a-z0-9:_-]+", "_", f"genus:{scope_taxon}:{key}".lower()).strip("_")
    return {
        "specialist_id": sid, "scope_type": "genus", "scope_taxon": scope_taxon, "canonical_name": cname,
        "qualified_regex": qreg, "ipni_abbr": ipni_abbr, "orcid": orcid, "wikidata_qid": wikidata_qid,
        "openalex_works_count": works, "bionomia_det_count": dets, "tier_cap": tier_cap, "source": source, "user_verified": 0,
    }

def _query_openalex_specialists(target_genus: str, target_family: Optional[str], http_fn: Optional[Callable[[str], Any]]) -> List[Dict[str, Any]]:
    payload = _http_json("https://api.openalex.org/works?" + _uparse.urlencode({"search": f"{target_genus} taxonomy systematics revision", "per-page": 50}), http_fn=http_fn)
    if not isinstance(payload, Mapping):
        return []
    g_low, f_low, agg = target_genus.strip().lower(), (target_family or "").strip().lower(), {}
    for idx, work in enumerate(payload.get("results") or []):
        if not isinstance(work, Mapping):
            continue
        loc = work.get("primary_location") or work.get("host_venue") or {}
        src = str(((loc.get("source") or {}) if isinstance(loc, Mapping) else {}).get("display_name") or loc.get("display_name") or "")
        concepts = " ".join(str(c.get("display_name") or "") for c in (work.get("concepts") or []) + (work.get("topics") or []) if isinstance(c, Mapping))
        blob = f"{work.get('title') or work.get('display_name') or ''} {src} {concepts}".lower()
        if not ((g_low in blob or (f_low and f_low in blob)) and any(k in blob for k in _TAXONOMIC_CONCEPTS)):
            continue
        wid = str(work.get("id") or f"w_{idx}")
        for rel in work.get("authorships") or []:
            auth = (rel.get("author") if isinstance(rel.get("author"), Mapping) else rel) if isinstance(rel, Mapping) else {}
            cname, orcid = str(auth.get("display_name") or "").strip(), _norm_orcid(auth.get("orcid"))
            giv, sur = _parse_name_parts(cname)
            if not (cname and sur and giv):
                continue
            rec = agg.setdefault(_author_key(cname, orcid), {"cname": cname, "orcid": orcid, "works": set(), "exp": 0})
            rec["works"].add(wid); rec["exp"] = max(rec["exp"], _to_int_safe(auth.get("works_count")))
    return [_new_rec(target_genus, k, r["cname"], qreg, "openalex", orcid=r["orcid"], works=w_cnt)
            for k, r in agg.items() if (w_cnt := max(len(r["works"]), int(r["exp"]))) >= 2 and (qreg := _build_qualified_regex(r["cname"]))]

def _split_wcvp_authors(accepted_authorship: str) -> List[str]:
    stripped = re.sub(r"\bet\s+al\.?", "", re.sub(r"\([^)]*\)", " ", str(accepted_authorship or "")), flags=re.IGNORECASE).strip(" ,;.")
    tokens: List[str] = []
    for p in re.split(r"\s*(?:&|;|\band\b|\bex\b)\s*", stripped, flags=re.IGNORECASE):
        if not (part := p.strip(" ,;.")):
            continue
        if "," in part and not re.match(r"^[A-Za-zÀ-ÿ'’-]+\s*,\s*[A-Za-z]\.?(?:\s*[A-Za-z]\.?)*$", part):
            tokens.extend(sub.strip(" ,.") for sub in part.split(",") if sub.strip(" ,."))
        else:
            tokens.append(part)
    return [t for t in tokens if len(_ascii_fold(re.sub(r"[^A-Za-z]", "", t))) >= 2]

def _query_wikidata_p428(tokens: Sequence[str], http_fn: Optional[Callable[[str], Any]]) -> Dict[str, Dict[str, Any]]:
    if not (toks := sorted({t.strip() for t in tokens if t.strip()})[:40]):
        return {}
    vals = " ".join(f'"{t.replace(chr(34), "")}"' for t in toks)
    sparql = (f"SELECT ?item ?itemLabel ?p428 ?orcid ?birth ?death ?floruit WHERE {{ VALUES ?p428 {{ {vals} }} ?item wdt:P428 ?p428 . "
              f"OPTIONAL {{ ?item wdt:P496 ?orcid . }} OPTIONAL {{ ?item wdt:P569 ?birth . }} OPTIONAL {{ ?item wdt:P570 ?death . }} "
              f'OPTIONAL {{ ?item wdt:P1317 ?floruit . }} SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }} }}')
    payload = _http_json("https://query.wikidata.org/sparql?" + _uparse.urlencode({"format": "json", "query": sparql}), http_fn=http_fn)
    rows = (((payload.get("results") or {}).get("bindings")) if isinstance(payload.get("results"), Mapping) else (payload.get("bindings") or payload.get("authors") or [])) if isinstance(payload, Mapping) else []
    out: Dict[str, Dict[str, Any]] = {}
    for b in rows or []:
        if not isinstance(b, Mapping):
            continue
        gv = lambda k: str(b[k].get("value") if isinstance(b.get(k), Mapping) else (b.get(k) or "")).strip()
        if not (p428 := gv("p428") or gv("ipni_abbr")):
            continue
        out[p428.lower()] = {
            "ipni_abbr": p428, "label": gv("itemLabel") or gv("canonical_name") or p428, "orcid": _norm_orcid(gv("orcid")),
            "wikidata_qid": (gv("item") or gv("wikidata_qid")).rsplit("/", 1)[-1] or None,
            "birth_year": _year_from_val(gv("birth") or b.get("birth_year")), "death_year": _year_from_val(gv("death") or b.get("death_year")),
            "floruit_year": _year_from_val(gv("floruit") or b.get("floruit_year")),
        }
    return out

def _extract_wcvp_wikidata(target_genus: str, tax_dir: Optional[Path], wcvp_profiles: Optional[Mapping[str, Any]],
                           sync_only: bool, wikidata_http_fn: Optional[Callable[[str], Any]]) -> List[Dict[str, Any]]:
    profs: List[Mapping[str, Any]] = ([v for v in wcvp_profiles.values() if isinstance(v, Mapping)] if isinstance(wcvp_profiles, Mapping)
                                      else ([v for v in wcvp_profiles if isinstance(v, Mapping)] if isinstance(wcvp_profiles, Sequence) else []))
    if not profs and tax_dir is not None:
        profs = next((r for fn in ("public_taxon_backbone_profiles.csv", "public_taxon_backbone.csv") if (r := _read_csv_rows(tax_dir / fn))), [])
    g_low, cands = target_genus.strip().lower(), []
    for prof in profs:
        tg, sp = str(prof.get("target_genus") or "").strip().lower(), str(prof.get("species_standardized") or prof.get("accepted_scientific_name") or "").strip().lower()
        if (tg and tg != g_low) or (not tg and sp and not sp.startswith(f"{g_low} ")):
            continue
        if (pub_y := _year_from_val(prof.get("publication_year") or prof.get("publication_citation"))) is not None and pub_y < 1950:
            continue
        bas_a = str(prof.get("basionym_authorship") or "").strip().lower()
        cands.extend((tok, pub_y) for tok in _split_wcvp_authors(str(prof.get("accepted_authorship") or "")) if not (bas_a and tok.lower() == bas_a))
    wd_map = _query_wikidata_p428([c[0] for c in cands], wikidata_http_fn) if (cands and (wikidata_http_fn is not None or not sync_only)) else {}
    out: Dict[str, Dict[str, Any]] = {}
    for tok, pub_y in cands:
        wd = wd_map.get(tok.lower(), {})
        by, dy, fy = wd.get("birth_year"), wd.get("death_year"), wd.get("floruit_year")
        if (dy is not None and dy < 1950) or (by is not None and by < 1900) or not ((pub_y is not None and pub_y >= 1950) or (by is not None and by >= 1900) or (fy is not None and fy >= 1950)):
            continue
        cname = str(wd.get("label") or tok).strip()
        if qreg := _build_qualified_regex(cname, ipni_abbr=tok, allow_unambiguous_bare=True):
            key = _author_key(cname, orcid=wd.get("orcid"))
            out[key] = _new_rec(target_genus, key, cname, qreg, "wcvp_wikidata", ipni_abbr=wd.get("ipni_abbr") or tok, orcid=wd.get("orcid"), wikidata_qid=wd.get("wikidata_qid"))
    return list(out.values())

def _enrich_bionomia(target_genus: str, merged: Dict[str, Dict[str, Any]], http_fn: Optional[Callable[[str], Any]]) -> None:
    if http_fn is None or not isinstance((payload := _http_json(f"https://api.bionomia.net/users/search.json?{_uparse.urlencode({'q': target_genus})}", http_fn=http_fn)), (dict, list)):
        return
    items = payload if isinstance(payload, list) else (payload.get("records") or payload.get("results") or payload.get("determiners") or payload.get("specimens") or [])
    occ_agg, summaries = {}, []
    for it in items:
        if not isinstance(it, Mapping):
            continue
        if "institutionCode" in it or "identifiedBy" in it:
            if not (det := str(it.get("identifiedBy") or "").strip()):
                continue
            r = occ_agg.setdefault(_author_key(det, _norm_orcid(it.get("orcid"))), {"cname": det, "orcid": _norm_orcid(it.get("orcid")), "insts": set(), "vars": set(), "cnt": 0, "dated": 0})
            r["cnt"] += 1; r["vars"].add(det); r["dated"] += int(bool(_year_from_val(it.get("dateIdentified") or it.get("eventDate"))))
            if (inst := str(it.get("institutionCode") or "").strip().upper()) not in _INVALID_INSTS:
                r["insts"].add(inst)
        elif cname := str(it.get("canonical_name") or it.get("fullname") or it.get("name") or "").strip():
            insts = {str(x).strip().upper() for x in (it.get("institutionCodes") or it.get("institutions") or []) if str(x).strip().upper() not in _INVALID_INSTS}
            n_i, cnt = int(it.get("institution_count") or len(insts)), int(it.get("bionomia_det_count") or it.get("determinations_count") or it.get("identified_count") or it.get("count") or 0)
            summaries.append({"cname": cname, "orcid": _norm_orcid(it.get("orcid")), "n_inst": n_i, "cnt": cnt, "vars": [str(v).strip() for v in (it.get("variants") or it.get("name_variants") or [cname]) if str(v).strip()], "mass": bool(it.get("is_mass_transfer") or it.get("undated_batch")) or (n_i <= 1 and cnt >= 20 and it.get("dated") is False)})
    summaries.extend({"cname": r["cname"], "orcid": r["orcid"], "n_inst": len(r["insts"]), "cnt": r["cnt"], "vars": sorted(r["vars"]), "mass": len(r["insts"]) <= 1 and r["dated"] == 0} for r in occ_agg.values())
    for s in summaries:
        if not (pats := [p for v in s["vars"] if (p := _build_qualified_regex(v))]):
            continue
        key, gold = _author_key(s["cname"], s["orcid"]), (s["n_inst"] >= 3) and not s["mass"]
        if key in merged:
            ex = merged[key]
            ex["bionomia_det_count"] = max(_to_int_safe(ex.get("bionomia_det_count")), int(s["cnt"]))
            ex["qualified_regex"] = "|".join(dict.fromkeys([ex["qualified_regex"], *pats]))
            ex["source"] = "+".join(sorted((set(str(ex.get("source") or "").split("+")) | {"bionomia"}) - {""}))
            if not gold and _to_int_safe(ex.get("user_verified")) != 1 and not str(ex.get("source", "")).startswith("seed_"):
                ex["tier_cap"] = "Tier_2_Silver"
        else:
            merged[key] = _new_rec(target_genus, key, s["cname"], "|".join(dict.fromkeys(pats)), "bionomia", tier_cap=("Tier_1_Gold" if gold else "Tier_2_Silver"), orcid=s["orcid"], dets=int(s["cnt"]))

def _load_db_specialists(db_path: Optional[Path], eff_dir: Path, *, read_only: bool) -> List[Dict[str, Any]]:
    if not (db_path is not None or is_canonical_or_explicit_db_context(eff_dir, workspace_root=PROJECT_ROOT)): return []
    res_db = resolve_db_path(db_path=db_path, workspace_root=PROJECT_ROOT)
    if str(res_db) == ":memory:" or not Path(res_db).is_file(): return []
    conn = get_db_connection(db_path=res_db, read_only=read_only, init_schema=(not read_only), workspace_root=PROJECT_ROOT)
    try: return load_taxonomic_specialists(conn)
    finally: conn.close()

def _merge_record(merged: Dict[str, Dict[str, Any]], row: Mapping[str, Any]) -> None:
    cname, qreg = str(row.get("canonical_name") or "").strip(), str(row.get("qualified_regex") or "").strip()
    if not (cname and qreg): return
    sid = str(row.get("specialist_id") or "").strip()
    key = f"seed:{sid.rsplit(':seed_', 1)[-1]}" if ":seed_" in sid else _author_key(cname, _norm_orcid(row.get("orcid")))
    if key not in merged:
        merged[key] = {**dict(row), "openalex_works_count": _to_int_safe(row.get("openalex_works_count")), "bionomia_det_count": _to_int_safe(row.get("bionomia_det_count")), "user_verified": _to_int_safe(row.get("user_verified"))}
        return
    ex = merged[key]
    ex["openalex_works_count"] = max(_to_int_safe(ex.get("openalex_works_count")), _to_int_safe(row.get("openalex_works_count")))
    ex["bionomia_det_count"] = max(_to_int_safe(ex.get("bionomia_det_count")), _to_int_safe(row.get("bionomia_det_count")))
    for k in ("ipni_abbr", "orcid", "wikidata_qid"): ex[k] = ex.get(k) or row.get(k)
    ex["qualified_regex"] = "|".join(dict.fromkeys([ex["qualified_regex"], qreg]))
    ex["source"] = "+".join(sorted((set(str(ex.get("source") or "").split("+")) | set(str(row.get("source") or "").split("+"))) - {""}))
    if (uv := _to_int_safe(row.get("user_verified"))) != 0:
        ex["user_verified"] = uv
        if uv == 1: ex["tier_cap"] = "Tier_1_Gold"

def discover_and_sync_specialists(
    target_genus: str, target_family: Optional[str] = None, taxonomy_dir: Optional[Path] = None,
    db_path: Optional[Path] = None, wcvp_profiles: Optional[Mapping[str, Any]] = None, force: bool = False,
    sync_only: bool = False, openalex_http_fn: Optional[Callable[[str], Any]] = None,
    wikidata_http_fn: Optional[Callable[[str], Any]] = None, bionomia_http_fn: Optional[Callable[[str], Any]] = None,
) -> List[Dict[str, Any]]:
    """Run 4-layer non-circular specialist discovery and sync flat CSV + Tier A `taxonomic_specialists`."""
    eff_g, eff_f = (target_genus or resolve_target_genus() or "Packera").strip(), (target_family or resolve_target_family() or "Asteraceae").strip()
    eff_dir = Path(taxonomy_dir) if taxonomy_dir is not None else (PROJECT_ROOT / "data" / "taxonomy")
    eff_dir.mkdir(parents=True, exist_ok=True)
    csv_path = eff_dir / SPECIALISTS_CSV_NAME
    sync_db = db_path is not None or is_canonical_or_explicit_db_context(eff_dir, workspace_root=PROJECT_ROOT)
    is_demo = os.environ.get("PLANTMATICS_DEMO_MODE") == "1" or os.environ.get("PACKERA_DEMO_MODE") == "1"
    existing = _read_csv_rows(csv_path) + ([] if is_demo else _load_db_specialists(db_path, eff_dir, read_only=False))
    merged, overrides = {}, []
    for row in existing:
        if str(row.get("scope_taxon") or "").strip().lower() not in {eff_g.lower(), eff_f.lower()}: continue
        if _to_int_safe(row.get("user_verified")) in (1, -1): overrides.append(dict(row))
        if not force: _merge_record(merged, row)
    if not sync_only or openalex_http_fn is not None:
        for row in _query_openalex_specialists(eff_g, eff_f, openalex_http_fn): _merge_record(merged, row)
    for row in _extract_wcvp_wikidata(eff_g, eff_dir, wcvp_profiles, sync_only, wikidata_http_fn): _merge_record(merged, row)
    _enrich_bionomia(eff_g, merged, bionomia_http_fn)
    if eff_g.lower() == "packera":
        for idx, pat in enumerate(SPECIALIST_PATTERNS):
            merged.setdefault(f"seed:{idx:03d}", _new_rec("Packera", f"seed_{idx:03d}", f"Packera Seed {idx + 1}", pat, "seed_packera"))
    for ov in overrides: _merge_record(merged, ov)
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    out_rows = sorted([{**{k: r.get(k) for k in SPECIALIST_FIELDS}, "updated_at": str(r.get("updated_at") or now)} for r in merged.values()],
                      key=lambda r: (-_to_int_safe(r.get("user_verified")), str(r.get("tier_cap") or ""), str(r.get("canonical_name") or "")))
    other_rows = [r for r in _read_csv_rows(csv_path) if str(r.get("scope_taxon") or "").strip().lower() not in {eff_g.lower(), eff_f.lower()}]
    atomic_write_csv_rows(csv_path, SPECIALIST_FIELDS, other_rows + out_rows)
    if sync_db and not is_demo and out_rows:
        conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=PROJECT_ROOT)
        try: upsert_taxonomic_specialists(conn, out_rows, demo=False, preserve_existing_on_null=False)
        finally: conn.close()
    _PATTERN_CACHE.clear()
    return out_rows

def _load_patterns_uncached(eff_g: str, eff_f: str, db_path: Optional[Path], eff_dir: Path) -> Tuple[List[str], List[str]]:
    rows = _read_csv_rows(eff_dir / SPECIALISTS_CSV_NAME) + _load_db_specialists(db_path, eff_dir, read_only=True)
    by_id = {str(r["specialist_id"]): r for r in rows if str(r.get("scope_taxon") or "").strip().lower() in {eff_g.lower(), eff_f.lower()} and str(r.get("specialist_id") or "").strip()}
    active, excluded, excl_names = [], [], []
    for r in by_id.values():
        try: re.compile(pat := str(r.get("qualified_regex") or "").strip())
        except re.error: continue
        if not pat: continue
        uv = _to_int_safe(r.get("user_verified"))
        if uv == -1:
            excluded.append(pat)
            if cn := str(r.get("canonical_name") or "").strip(): excl_names.append(cn)
        elif uv == 1 or str(r.get("tier_cap") or "") == "Tier_1_Gold": active.append(pat)
    if eff_g.lower() == "packera":
        for pat in SPECIALIST_PATTERNS:
            if pat not in excluded and not any(re.search(pat, nm, re.IGNORECASE) for nm in excl_names): active.append(pat)
    return list(dict.fromkeys(active)), list(dict.fromkeys(excluded))

def load_specialist_pattern_sets(
    target_genus: Optional[str] = None, target_family: Optional[str] = None,
    db_path: Optional[Path] = None, taxonomy_dir: Optional[Path] = None,
) -> Tuple[List[str], List[str]]:
    """Return `(active_gold_patterns, excluded_patterns)` cached by DB/CSV mtimes."""
    eff_g, eff_f = (target_genus or resolve_target_genus() or "Packera").strip(), (target_family or resolve_target_family() or "Asteraceae").strip()
    eff_dir = Path(taxonomy_dir) if taxonomy_dir is not None else Path(os.environ.get("PLANTMATICS_TAXONOMY_DIR") or os.environ.get("PACKERA_TAXONOMY_DIR") or (PROJECT_ROOT / "data" / "taxonomy"))
    res_db = resolve_db_path(db_path=db_path, workspace_root=PROJECT_ROOT)
    if str(res_db) == ":memory:":
        return _load_patterns_uncached(eff_g, eff_f, db_path, eff_dir)
    db_p, csv_p, wal_p = Path(res_db), eff_dir / SPECIALISTS_CSV_NAME, Path(f"{res_db}-wal")
    key = (eff_g.lower(), eff_f.lower(), str(db_p), str(eff_dir), db_p.stat().st_mtime_ns if db_p.is_file() else 0, wal_p.stat().st_mtime_ns if wal_p.is_file() else 0, csv_p.stat().st_mtime_ns if csv_p.is_file() else 0)
    if key not in _PATTERN_CACHE:
        _PATTERN_CACHE[key] = _load_patterns_uncached(eff_g, eff_f, db_path, eff_dir)
    return _PATTERN_CACHE[key]

def load_active_specialist_patterns(
    target_genus: Optional[str] = None, target_family: Optional[str] = None,
    db_path: Optional[Path] = None, taxonomy_dir: Optional[Path] = None,
) -> List[str]:
    """Load compiled-safe active Tier_1_Gold specialist regex patterns for `target_genus`."""
    return list(load_specialist_pattern_sets(target_genus, target_family, db_path, taxonomy_dir)[0])
