#!/usr/bin/env python3
"""
Module: scripts/data_prep/step06_voucher_linker_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  SQLite voucher & cytology prior loaders, exsiccatae duplicate indexing,
  and 4-tier Step 06 sample-to-voucher linking (`G1` through `G4`).
"""

from __future__ import annotations

import re, sqlite3
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from scripts.data_prep.step06_sample_metadata_utils import (
    _NON_TYPE_STATUS_TOKENS, _canonical_taxon_equiv_key, _normalize_event_date_part,
    _normalize_exsiccatae_key, _safe_text, sanitize_step06_sample_coordinates, strip_smart_quotes,
)


def _load_db_vouchers(conn: Optional[sqlite3.Connection]) -> List[Dict[str, Any]]:
    if conn is None:
        return []
    try:
        cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()}
    except sqlite3.Error:
        return []
    if not cols:
        return []
    select_cols = [c for c in ("catalogNumber", "recordNumber", "recordedBy", "scientificName", "species_standardized", "decimalLatitude", "decimalLongitude", "stateProvince", "county", "eventDate", "determiner_tier", "type_status", "image_path", "exsiccatae_key", "is_primary_duplicate") if c in cols]
    return [{c: r[idx] for idx, c in enumerate(select_cols)} for r in conn.execute(f"SELECT {', '.join(select_cols)} FROM vouchers;").fetchall()]


def _normalize_exsiccatae_lookup_key(raw: Any) -> str:
    return " ".join(str(raw or "").strip().split())


def _parse_int_ploidy(val: Any) -> Optional[int]:
    if val is None or isinstance(val, bool):
        return None
    txt = str(val).strip()
    if not txt or txt.upper() in ("NA", "NAN", "NONE", "NULL"):
        return None
    txt = txt[:-1].strip() if txt.lower().endswith("x") else txt
    try:
        iv = int(round(float(txt)))
    except ValueError:
        return None
    return iv if 2 <= iv <= 16 else None


def _load_db_cytology_priors(conn_or_db_path: Optional[Any]) -> Dict[str, Dict[str, Any]]:
    """Safely load Step 06 empirical Bayes cytology priors from SQLite tables."""
    empty: Dict[str, Dict[str, Any]] = {"by_catalog": {}, "by_exsiccatae": {}, "by_genbank": {}, "by_taxon": {}}
    if conn_or_db_path is None:
        return empty
    conn, owns_conn = None, False
    if isinstance(conn_or_db_path, (str, Path)):
        db_p = Path(conn_or_db_path)
        if not db_p.exists():
            return empty
        try:
            conn, owns_conn = sqlite3.connect(f"file:{db_p}?mode=ro", uri=True), True
        except sqlite3.Error:
            return empty
    else:
        conn = conn_or_db_path
    try:
        by_catalog, by_exsiccatae, by_genbank, by_taxon = {}, {}, {}, {}
        try:
            cyto_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(cytological_records);").fetchall()}
        except sqlite3.Error:
            cyto_cols = set()
        ploidy_col = "ploidy_level" if "ploidy_level" in cyto_cols else ("voucher_ploidy_level" if "voucher_ploidy_level" in cyto_cols else None)
        if ploidy_col:
            tier_col = "cytology_linkage_tier" if "cytology_linkage_tier" in cyto_cols else ("linkage_tier" if "linkage_tier" in cyto_cols else None)
            wanted = [c for c in ("catalogNumber", "exsiccatae_key", "recordedBy", "recordNumber", "eventDate", "genbank_accession", ploidy_col, tier_col) if c and c in cyto_cols]
            col_pos = {c: idx for idx, c in enumerate(wanted)}
            cat_buckets: Dict[str, Tuple[int, set]] = {}
            ex_buckets: Dict[str, Tuple[str, set]] = {}
            gb_buckets: Dict[str, set] = {}
            for row in conn.execute(f"SELECT {', '.join(wanted)} FROM cytological_records;").fetchall():
                if (pl := _parse_int_ploidy(row[col_pos[ploidy_col]])) is None:
                    continue
                tier_up = (str(row[col_pos[tier_col]] or "").strip().upper() if tier_col and tier_col in col_pos else "")
                if (not tier_up or tier_up.startswith("C1_DIRECT") or tier_up.startswith("C2_EXSICCATAE")) and "catalogNumber" in col_pos:
                    if cat := str(row[col_pos["catalogNumber"]] or "").strip():
                        rank = 0 if tier_up.startswith("C1_") else (1 if tier_up.startswith("C2_") else 2)
                        cur_rank, cur_set = cat_buckets.get(cat, (99, set()))
                        if rank < cur_rank: cat_buckets[cat] = (rank, {pl})
                        elif rank == cur_rank: cur_set.add(pl); cat_buckets[cat] = (cur_rank, cur_set)
                ex_norm = _normalize_exsiccatae_lookup_key(row[col_pos["exsiccatae_key"]] if "exsiccatae_key" in col_pos else "")
                if not ex_norm and "recordedBy" in col_pos and "recordNumber" in col_pos:
                    try:
                        ex_norm = _normalize_exsiccatae_lookup_key(_normalize_exsiccatae_key(row[col_pos["recordedBy"]], row[col_pos["recordNumber"]], row[col_pos["eventDate"]] if "eventDate" in col_pos else None))
                    except Exception:
                        ex_norm = ""
                if ex_norm:
                    lk = ex_norm.lower(); canon, pl_set = ex_buckets.get(lk, (ex_norm, set())); pl_set.add(pl); ex_buckets[lk] = (canon, pl_set)
                if "genbank_accession" in col_pos and (gb := str(row[col_pos["genbank_accession"]] or "").strip()):
                    gb_buckets.setdefault(gb, set()).add(pl)
            for cat, (_, pl_set) in cat_buckets.items():
                if len(pl_set) == 1: by_catalog[cat] = next(iter(pl_set))
            for canon, pl_set in ex_buckets.values():
                if len(pl_set) == 1: by_exsiccatae[canon] = next(iter(pl_set))
            for gb, pl_set in gb_buckets.items():
                if len(pl_set) == 1: by_genbank[gb] = next(iter(pl_set))
        try:
            tax_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(taxon_cytology_summary);").fetchall()}
        except sqlite3.Error:
            tax_cols = set()
        if "species_standardized" in tax_cols:
            t_wanted = [c for c in ("species_standardized", "modal_ploidy_level", "is_multi_cytotype_taxon") if c in tax_cols]
            t_pos = {c: idx for idx, c in enumerate(t_wanted)}
            for row in conn.execute(f"SELECT {', '.join(t_wanted)} FROM taxon_cytology_summary;").fetchall():
                if not (sp := str(row[t_pos["species_standardized"]] or "").strip()):
                    continue
                modal = _parse_int_ploidy(row[t_pos["modal_ploidy_level"]]) if "modal_ploidy_level" in t_pos else None
                mraw = row[t_pos["is_multi_cytotype_taxon"]] if "is_multi_cytotype_taxon" in t_pos else 0
                is_multi = (1 if mraw.strip().lower() in ("1", "true", "yes") else 0) if isinstance(mraw, str) else (1 if bool(mraw) else 0)
                by_taxon[sp] = {"modal_ploidy_level": modal, "is_multi_cytotype_taxon": is_multi}
        return {"by_catalog": by_catalog, "by_exsiccatae": by_exsiccatae, "by_genbank": by_genbank, "by_taxon": by_taxon}
    except sqlite3.Error:
        return empty
    finally:
        if owns_conn and conn is not None:
            conn.close()


def _parse_pop_aliases(raw_val: Any) -> List[str]:
    return [t.strip().upper() for t in re.split(r"[,;|]+", str(raw_val or "")) if t.strip()]


def _match_field_record(sample_id: str, parsed: Dict[str, Optional[str]], field_rows: Sequence[Mapping[str, Any]]) -> Optional[Mapping[str, Any]]:
    pop, code, sp_std = (parsed.get("population_id") or "").upper(), (parsed.get("taxon_code") or "").upper(), parsed.get("species_standardized")
    if not pop:
        return None
    if exact_pid := [r for r in field_rows if r.get("project_id") and strip_smart_quotes(r["project_id"]).upper() == sample_id.upper()]:
        return exact_pid[0]
    candidates = [r for r in field_rows if (r.get("population_id") or "").upper() == pop or pop in _parse_pop_aliases(r.get("population_alias"))]
    if not candidates:
        return None
    if code == "NOV":
        if novel_cands := [r for r in candidates if r.get("is_putative_novel_or_hybrid_query")]:
            return novel_cands[0]
    if same_code := [r for r in candidates if (r.get("taxon_code") or "").upper() == code and code]:
        return same_code[0]
    sp_eq = _canonical_taxon_equiv_key(sp_std)
    if same_sp := [r for r in candidates if sp_eq and _canonical_taxon_equiv_key(r.get("species_standardized")) == sp_eq]:
        return same_sp[0]
    return None


def _coerce_flag_int(val: Any) -> int:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return 0
    if isinstance(val, bool):
        return 1 if val else 0
    s = str(val).strip().lower()
    if not s or s in {"false", "f", "no", "n", "nan", "none", "null", "na"}:
        return 0
    if s in {"true", "t", "yes", "y"}:
        return 1
    try:
        return int(float(s))
    except ValueError:
        return 0


def _authority_tier_rank(voucher: Mapping[str, Any]) -> int:
    det = str(voucher.get("determiner_tier") or "").strip().lower()
    type_stat = str(voucher.get("type_status") or "").strip().upper()
    is_type = bool(type_stat and type_stat not in _NON_TYPE_STATUS_TOKENS and not re.search(r"\b(?:NOT\s+A?\s*TYPE|NON[-\s]?TYPE)\b", type_stat))
    if is_type or det.startswith("tier_0"): return 3
    if det.startswith("tier_1") or det == "gold": return 2
    if det.startswith("tier_2") or det == "silver": return 1
    return 0


def _duplicate_rank_key(voucher: Mapping[str, Any]) -> Tuple[int, int, int, str]:
    return (-int(bool(_safe_text(voucher.get("image_path")))), -_coerce_flag_int(voucher.get("is_primary_duplicate")), -_authority_tier_rank(voucher), str(voucher.get("catalogNumber") or ""))


def _voucher_exsiccatae_keys(voucher: Mapping[str, Any]) -> Tuple[str, str, str]:
    stored = _safe_text(voucher.get("exsiccatae_key"))
    col, num, dt = str(voucher.get("recordedBy") or ""), str(voucher.get("recordNumber") or ""), str(voucher.get("eventDate") or "")
    return stored, _normalize_exsiccatae_key(col, num, dt), _normalize_exsiccatae_key(col, num, "")


def _build_exsiccatae_index(vouchers: Sequence[Mapping[str, Any]]) -> Dict[str, List[Mapping[str, Any]]]:
    index: Dict[str, List[Mapping[str, Any]]] = {}
    seen_per_key: Dict[str, set] = {}
    def _add(key: str, v: Mapping[str, Any]) -> None:
        if not key: return
        cat = str(v.get("catalogNumber") or "")
        bucket = seen_per_key.setdefault(key, set())
        if cat and cat in bucket: return
        if cat: bucket.add(cat)
        index.setdefault(key, []).append(v)
    for v in vouchers:
        stored, norm, dateless = _voucher_exsiccatae_keys(v)
        if stored: _add(stored, v); _add(stored.lower(), v)
        if norm: _add(norm, v)
        if dateless: _add(f"__dateless__:{dateless}", v)
    return index


def _extract_meta_exsiccatae_keys(meta: Optional[Mapping[str, Any]]) -> Tuple[str, str, str]:
    if not meta:
        return "", "", ""
    explicit = _safe_text(meta.get("exsiccatae_key") or meta.get("Exsiccatae_Key"))
    col = str(meta.get("recordedBy") or meta.get("recorded_by") or meta.get("collector") or "")
    num = str(meta.get("recordNumber") or meta.get("record_number") or meta.get("collector_number") or meta.get("collection_number") or meta.get("collector_id") or meta.get("Collector_ID") or "")
    dt = str(meta.get("eventDate") or meta.get("event_date") or meta.get("date_collected") or "")
    return explicit, _normalize_exsiccatae_key(col, num, dt), _normalize_exsiccatae_key(col, num, "")


def _lookup_exsiccatae_candidates(explicit_key: str, norm_key: str, dateless_key: str, exsiccatae_index: Mapping[str, Sequence[Mapping[str, Any]]]) -> List[Mapping[str, Any]]:
    hits: List[Mapping[str, Any]] = []
    seen_cats: set = set()
    def _extend(rows: Sequence[Mapping[str, Any]]) -> None:
        for v in rows:
            cat = str(v.get("catalogNumber") or "")
            if cat and cat in seen_cats: continue
            if cat: seen_cats.add(cat)
            hits.append(v)
    if explicit_key:
        _extend(exsiccatae_index.get(explicit_key) or exsiccatae_index.get(explicit_key.lower()) or ())
    if norm_key:
        _extend(exsiccatae_index.get(norm_key) or ())
        if dateless_key and norm_key != dateless_key:
            for cand in exsiccatae_index.get(f"__dateless__:{dateless_key}", ()):
                if not _normalize_event_date_part(cand.get("eventDate")): _extend((cand,))
        elif dateless_key:
            dbuck = list(exsiccatae_index.get(f"__dateless__:{dateless_key}", ()))
            if len({_normalize_event_date_part(c.get("eventDate")) for c in dbuck if _normalize_event_date_part(c.get("eventDate"))}) <= 1:
                _extend(dbuck)
    return hits


def _resolve_output_exsiccatae_key(voucher: Optional[Mapping[str, Any]], user_meta: Optional[Mapping[str, Any]], field_row: Optional[Mapping[str, Any]]) -> Optional[str]:
    if voucher:
        v_stored, v_norm, _ = _voucher_exsiccatae_keys(voucher)
        if v_stored: return v_stored
        if v_norm: return v_norm
    for meta in (user_meta, field_row):
        exp_k, norm_k, _ = _extract_meta_exsiccatae_keys(meta)
        if exp_k: return exp_k
        if norm_k: return norm_k
    return None


def _score_voucher_match(voucher: Mapping[str, Any], code: str, sp_std: Optional[str], pop_tok: str, is_alias_fallback: bool) -> Optional[Tuple[Any, ...]]:
    rec, cat = strip_smart_quotes(voucher.get("recordNumber") or ""), strip_smart_quotes(voucher.get("catalogNumber") or "")
    v_sp = str(voucher.get("species_standardized") or "").strip()
    pat = re.compile(rf"(?:^|[-_]){re.escape(pop_tok)}(?:[-_]|$)", re.IGNORECASE)
    if not ((rec and pat.search(rec)) or (cat and pat.search(cat))):
        return None
    rec_up, code_up = rec.upper(), code.upper()
    prefix_compatible = prefix_incompatible = False
    if rec_up.startswith("PAC") and len(rec_up) >= 6:
        rec_code = "NAUR" if rec_up.startswith("PACNAUR") else rec_up[3:6]
        if rec_code == code_up or {rec_code, code_up}.issubset({"TOM", "DUB"}) or (rec_code == "TX" and code_up in ("TOM", "DUB")):
            prefix_compatible = True
        else:
            prefix_incompatible = True
    if prefix_incompatible:
        return None
    sp_eq, v_sp_eq = _canonical_taxon_equiv_key(sp_std), _canonical_taxon_equiv_key(v_sp)
    sp_match = bool(sp_eq and v_sp_eq == sp_eq)
    if not prefix_compatible and sp_eq and v_sp_eq and v_sp_eq != sp_eq:
        return None
    has_coords = 0 if (voucher.get("decimalLatitude") is not None and voucher.get("decimalLongitude") is not None) else 1
    return (1 if is_alias_fallback else 0, 0 if prefix_compatible else 1, 0 if sp_match else 1, 3 - _authority_tier_rank(voucher), has_coords, 0 if _safe_text(voucher.get("image_path")) else 1, -_coerce_flag_int(voucher.get("is_primary_duplicate")), str(voucher.get("catalogNumber") or ""))


def _link_sample_to_voucher(
    sample_id: str, parsed: Dict[str, Optional[str]], user_meta: Optional[Mapping[str, Any]],
    field_row: Optional[Mapping[str, Any]], vouchers: Sequence[Mapping[str, Any]],
    exsiccatae_index: Optional[Mapping[str, Sequence[Mapping[str, Any]]]] = None,
) -> Dict[str, Any]:
    """Resolve 4-tier voucher linkage (`G1_Direct_Voucher` .. `G4_Ungeoreferenced_Tip`) with coordinate hygiene."""
    voucher_by_cat = {str(v.get("catalogNumber") or ""): v for v in vouchers if v.get("catalogNumber")}
    voucher_by_cat_ci = {k.upper(): v for k, v in voucher_by_cat.items()}
    eff_index = exsiccatae_index if exsiccatae_index is not None else _build_exsiccatae_index(vouchers)

    def _siblings(v: Mapping[str, Any]) -> List[Mapping[str, Any]]:
        return _lookup_exsiccatae_candidates(*_voucher_exsiccatae_keys(v), eff_index)

    def _upgrade(pv: Mapping[str, Any]) -> Mapping[str, Any]:
        if _safe_text(pv.get("image_path")):
            return pv
        sibs = list(_siblings(pv))
        u_exp, u_norm, u_dl = _extract_meta_exsiccatae_keys(user_meta)
        if u_exp or u_norm:
            sibs.extend(_lookup_exsiccatae_candidates(u_exp, u_norm, u_dl, eff_index))
        imaged = [s for s in sibs if _safe_text(s.get("image_path"))]
        return min(imaged, key=_duplicate_rank_key) if imaged else pv

    def _ret_g1(v: Mapping[str, Any], tier: str = "G1_Direct_Voucher") -> Dict[str, Any]:
        return {"catalogNumber": str(v.get("catalogNumber")), "voucher": v, "linkage_tier": tier, "exsiccatae_key": _resolve_output_exsiccatae_key(v, user_meta, field_row)}

    valid_tiers = {"G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G3_Georeferenced_Unvouchered", "G4_Ungeoreferenced_Tip"}
    u_tier = str(user_meta.get("linkage_tier") or "").strip() if user_meta else ""
    u_tier = u_tier if u_tier in valid_tiers else ""
    if user_meta and user_meta.get("catalogNumber"):
        u_cat = str(user_meta["catalogNumber"]).strip()
        if matched := voucher_by_cat.get(u_cat) or voucher_by_cat_ci.get(u_cat.upper()):
            return _ret_g1(_upgrade(matched), u_tier or "G1_Direct_Voucher")

    u_exp, u_norm, u_dl = _extract_meta_exsiccatae_keys(user_meta)
    if (u_exp or u_norm) and u_tier not in ("G3_Georeferenced_Unvouchered", "G4_Ungeoreferenced_Tip"):
        if u_cands := _lookup_exsiccatae_candidates(u_exp, u_norm, u_dl, eff_index):
            return _ret_g1(min(u_cands, key=_duplicate_rank_key))

    if u_tier not in ("G3_Georeferenced_Unvouchered", "G4_Ungeoreferenced_Tip"):
        if field_row and field_row.get("project_id") and strip_smart_quotes(field_row["project_id"]).upper() == sample_id.upper():
            f_exp, f_norm, f_dl = _extract_meta_exsiccatae_keys(field_row)
            if (f_exp or f_norm) and (f_cands := _lookup_exsiccatae_candidates(f_exp, f_norm, f_dl, eff_index)):
                return _ret_g1(min(f_cands, key=_duplicate_rank_key))
        sp_eq_main = _canonical_taxon_equiv_key(parsed.get("species_standardized"))
        probe_ids = [sample_id] + ([str(user_meta["project_id"])] if user_meta and user_meta.get("project_id") else [])
        if field_row and field_row.get("project_id"):
            fpid = strip_smart_quotes(field_row["project_id"])
            fsp_eq = _canonical_taxon_equiv_key(field_row.get("species_standardized"))
            if not (sp_eq_main and fsp_eq and sp_eq_main != fsp_eq):
                probe_ids.extend([fpid, re.sub(r"-\d{8}$", "", fpid)])
        for pid in probe_ids:
            if not (pid_c := strip_smart_quotes(pid).upper()):
                continue
            cat_hits, rec_hits = [], []
            for v in vouchers:
                cat_c, rec_c = strip_smart_quotes(v.get("catalogNumber") or "").upper(), strip_smart_quotes(v.get("recordNumber") or "").upper()
                vsp_eq = _canonical_taxon_equiv_key(v.get("species_standardized"))
                if pid_c == cat_c: cat_hits.append(v)
                elif rec_c and (pid_c == rec_c or pid_c.startswith(rec_c + "-")) and not (sp_eq_main and vsp_eq and sp_eq_main != vsp_eq): rec_hits.append(v)
            if cat_hits:
                return _ret_g1(_upgrade(cat_hits[0]))
            if rec_hits:
                exp_v: List[Mapping[str, Any]] = []; seen_c: set = set()
                for rv in rec_hits:
                    for cand in [rv, *_siblings(rv)]:
                        cid_v = str(cand.get("catalogNumber") or "")
                        if cid_v and cid_v in seen_c: continue
                        if cid_v: seen_c.add(cid_v)
                        exp_v.append(cand)
                return _ret_g1(min(exp_v, key=_duplicate_rank_key))

        pop, code, sp_std = (parsed.get("population_id") or "").upper(), (parsed.get("taxon_code") or "").upper(), parsed.get("species_standardized")
        if pop and vouchers:
            aliases = list(dict.fromkeys([pop, *_parse_pop_aliases((user_meta or {}).get("population_alias")), *_parse_pop_aliases((field_row or {}).get("population_alias"))]))
            scored: List[Tuple[Tuple[Any, ...], Mapping[str, Any]]] = []
            for idx, tok in enumerate(aliases):
                if tok_hits := [(sc, v) for v in vouchers if (sc := _score_voucher_match(v, code=code, sp_std=sp_std, pop_tok=tok, is_alias_fallback=(idx > 0))) is not None]:
                    scored.extend(tok_hits); break
            if scored:
                scored.sort(key=lambda item: item[0])
                return _ret_g1(scored[0][1], "G2_Population_Voucher_Proxy")

    out_ex = _resolve_output_exsiccatae_key(None, user_meta, field_row)
    c_lat = (user_meta or {}).get("decimalLatitude") if (user_meta and user_meta.get("decimalLatitude") is not None) else (field_row or {}).get("decimalLatitude")
    c_lon = (user_meta or {}).get("decimalLongitude") if (user_meta and user_meta.get("decimalLongitude") is not None) else (field_row or {}).get("decimalLongitude")
    c_state = (user_meta or {}).get("stateProvince") or (field_row or {}).get("stateProvince") or ""
    c_county = (user_meta or {}).get("county") or (field_row or {}).get("county") or ""
    c_loc = (user_meta or {}).get("locality") or (field_row or {}).get("notes") or ""
    s_lat, s_lon, was_nullified = sanitize_step06_sample_coordinates(c_lat, c_lon, c_state, county=c_county, locality=c_loc, is_vouchered=False)
    if u_tier == "G4_Ungeoreferenced_Tip" or was_nullified:
        tier_out = "G4_Ungeoreferenced_Tip"
    elif u_tier == "G3_Georeferenced_Unvouchered":
        tier_out = "G3_Georeferenced_Unvouchered"
    else:
        tier_out = "G3_Georeferenced_Unvouchered" if (s_lat is not None and s_lon is not None) else "G4_Ungeoreferenced_Tip"
    return {"catalogNumber": None, "voucher": None, "linkage_tier": tier_out, "exsiccatae_key": out_ex}

