"""step01_concept_authority_utils.py - Step 01 taxonomic concept tracking, hybrid/open-nomenclature
gating, typified basionym extraction, determination year parsing, specialist constants, and segregate caps."""

from __future__ import annotations

import datetime
import json
import math
from pathlib import Path
import re
import sqlite3
from typing import Any, Dict, FrozenSet, List, Mapping, Optional, Set, Tuple, Union

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]

# Post-1960 Packera revisionists and segregate authors qualifying for Tier_1_Gold.
CORE_PACKERA_GOLD_SURNAMES: FrozenSet[str] = frozenset({
    "barkley", "kowal", "trock", "bain", "mahoney", "soderholm", "mandel",
    "moore-pollard", "weakley", "poindexter", "gramling", "boufford", "nesom", "uttal",
})

# Initial/given-name-qualified patterns for core specialists with common surnames
# (bare "estes", "freeman", "golden", "moore", "weber", "packer" remain Bronze).
CORE_PACKERA_GOLD_REGEXES: Tuple[re.Pattern, ...] = tuple(re.compile(p, re.IGNORECASE) for p in (
    r"\b(?:c\.?\s*c\.?|craig(?:\s+c\.?)?)\s+freeman\b|\bfreeman\s*,\s*(?:c\.?\s*c\.?|craig)\b",
    r"\b(?:j\.?\s*l\.?|judith\s+l\.?)\s+golden\b|\bgolden\s*,\s*(?:j\.?\s*l\.?|judith)\b",
    r"\b(?:e\.?\s*r\.?|erika\s+r\.?)\s+moore\b|\bmoore\s*,\s*(?:e\.?\s*r\.?|erika)\b",
    r"\b(?:l\.?\s*d\.?|dwayne)\s+estes\b|\bestes\s*,\s*(?:l\.?\s*d\.?|dwayne)\b",
    r"\b(?:w\.?\s*a\.?|william\s+a\.?)\s+weber\b|\bweber\s*,\s*(?:w\.?\s*a\.?|william)\b",
    r"\b(?:j\.?\s*g\.?|john\s+g\.?)\s+packer\b|\bpacker\s*,\s*(?:j\.?\s*g\.?|john)\b",
))

# Regional flora authorities and historical Senecioneae authors routed to Tier_2_Silver.
REGIONAL_AND_HISTORICAL_SILVER_SURNAMES: FrozenSet[str] = frozenset({
    "cronquist", "greenman", "pippen", "semple", "schilling", "strother", "brouillet",
    "kral", "godfrey", "sorrie", "floden", "wieboldt", "mcmillan", "judziewicz",
    "okennon", "okane", "ruyle-douglas", "yeatts", "spellenberg", "beauchamp",
    "mabberley", "fernald", "rydberg", "boivin",
})

_SILVER_PATTERN_STRS: Tuple[str, ...] = (
    r"\b(?:c\.?|charles)\s+jeffrey\b|\bjeffrey\s*,\s*(?:c\.?|charles)\b",
    r"\b(?:b\.?\s*l\.?|billie\s*(?:l(?:ee)?\.?)?)\s+turner\b|\bturner\s*,\s*(?:b\.?\s*l\.?|billie)\b",
    r"(?:^|[\s|,;(])[áa]\.?\s*l[öo]ve\b|\bd\.?\s*l[öo]ve\b|\b[áa]skell\s+l[öo]ve\b|\bdoris\s+l[öo]ve\b",
    r"\b(?:g\.?\s*w\.?|george\s+w\.?)\s+douglas\b|\bdouglas\s*,\s*(?:g\.?\s*w\.?|george)\b",
    r"\b(?:e\.?\s*l\.?|edward\s+lee)\s+greene\b|\bgreene\s*,\s*(?:e\.?\s*l\.?|edward)\b",
    r"\b(?:a\.?|asa)\s+gray\b|\bgray\s*,\s*(?:a\.?|asa)\b",
    r"\b(?:j\.?\s*k\.?|john\s+kunkel)\s+small\b|\bsmall\s*,\s*(?:j\.?\s*k\.?|john)\b",
    r"\b(?:alph\.?|alphonso)\s+wood\b|\bwood\s*,\s*(?:alph\.?|alphonso)\b",
    r"\b(?:s\.?\s*f\.?|sidney\s+fay)\s+blake\b|\bblake\s*,\s*(?:s\.?\s*f\.?|sidney)\b",
    r"\b(?:s\.?|sereno)\s+watson\b|\bwatson\s*,\s*(?:s\.?|sereno)\b",
    r"\b(?:s\.?\s*l\.?|stanley\s+l\.?)\s+welsh\b|\bwelsh\s*,\s*(?:s\.?\s*l\.?|stanley)\b",
    r"\b(?:n\.?\s*d\.?|duane)\s+atwood\b|\batwood\s*,\s*(?:n\.?\s*d\.?|duane)\b",
    r"\b(?:r\.?\s*d\.?|robert\s+d\.?)\s+dorn\b|\bdorn\s*,\s*(?:r\.?\s*d\.?|robert)\b",
    r"\b(?:n\.?\s*l\.?|nathaniel\s+l(?:ord)?\.?)\s+britton\b|\bbritton\s*,\s*(?:n\.?\s*l\.?|nathaniel)\b",
    r"\b(?:a\.?\s*a\.?|amos\s+a(?:rthur)?\.?)\s+heller\b|\bheller\s*,\s*(?:a\.?\s*a\.?|amos)\b",
    r"\b(?:j\.?\s*j\.?\s*n\.?|julian(?:\s+j\.?\s*n\.?)?)\s+campbell\b|\bcampbell\s*,\s*(?:j\.?\s*j\.?\s*n\.?|julian)\b",
    r"\bo['’]?kennon\b|\br\.?\s*j\.?\s*o['’]?kennon\b|\bo['’]?kane\b|\bs\.?\s*l\.?\s*o['’]?kane\b",
)
REGIONAL_AND_HISTORICAL_SILVER_REGEXES: Tuple[re.Pattern, ...] = tuple(re.compile(p, re.IGNORECASE) for p in _SILVER_PATTERN_STRS)

_DEFAULT_GENERA: FrozenSet[str] = frozenset({"Packera", "Senecio", "Cineraria"})
_NON_EPITHET: FrozenSet[str] = frozenset({"sp", "spp", "cf", "aff", "vel", "indet", "hybrid", "unknown", "na", "x", "var", "subsp", "ssp", "f", "forma", "nothosubsp", "nothovar", "sensu", "auct", "ex", "in", "sl"})
_TYPE_OF_RE = re.compile(r"\b(?:holo|iso|lecto|neo|epi|syn|isolecto|isoneo|isosyn|para|isopara|co)?type\s+of\s+([A-Z][a-z]+(?:\s+[a-z][a-z-]+).*)$", re.IGNORECASE)
_OPEN_QUAL_RE = re.compile(r"(?:\b(?:cf|aff|vel\s+aff)\b\.?|\?|\bs\.?\s*l\.(?!\s+[A-Z][a-z])|\bsensu\s+lato\b)", re.IGNORECASE)
_HYBRID_REMARK_RE = re.compile(r"\b(?:putative\s+hybrid|hybrid\s+with|hybrid\s+between|intermediate\s+(?:between|toward)|intergrade\s+with|likely\s+hybrid)\b", re.IGNORECASE)
_KNOWN_HOMOTYPIC_PAIRS: FrozenSet[Tuple[str, str]] = frozenset({("senecio tomentosus", "packera dubia"), ("packera tomentosa", "packera dubia"), ("cineraria dubia", "packera dubia"), ("senecio dubius", "packera dubia")})
_FALLBACK_POST1980_SEGREGATES: Dict[str, Tuple[int, Set[str], Set[str]]] = {"Packera serpenticola": (2014, {"NC"}, {"burke"}), "Packera insulae-regalis": (2011, {"MI"}, {"keweenaw"})}
_US_STATE_TO_ABBR: Dict[str, str] = {
    "ALABAMA": "AL", "ARKANSAS": "AR", "CONNECTICUT": "CT", "DELAWARE": "DE", "DISTRICT OF COLUMBIA": "DC",
    "FLORIDA": "FL", "GEORGIA": "GA", "ILLINOIS": "IL", "INDIANA": "IN", "IOWA": "IA", "KENTUCKY": "KY",
    "LOUISIANA": "LA", "MAINE": "ME", "MARYLAND": "MD", "MASSACHUSETTS": "MA", "MICHIGAN": "MI",
    "MINNESOTA": "MN", "MISSISSIPPI": "MS", "MISSOURI": "MO", "NEW HAMPSHIRE": "NH", "NEW JERSEY": "NJ",
    "NEW YORK": "NY", "NORTH CAROLINA": "NC", "OHIO": "OH", "PENNSYLVANIA": "PA", "RHODE ISLAND": "RI",
    "SOUTH CAROLINA": "SC", "TENNESSEE": "TN", "TEXAS": "TX", "VERMONT": "VT", "VIRGINIA": "VA",
    "WEST VIRGINIA": "WV", "WISCONSIN": "WI",
}


def _clean_str(val: Any) -> str:
    """Normalize Darwin Core scalars to single-spaced text; map null markers to ''."""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return ""
    if pd is not None:
        try:
            if pd.isna(val):
                return ""
        except (TypeError, ValueError):
            pass
    s = re.sub(r"\s+", " ", str(val).strip())
    return "" if s.lower() in {"none", "nan", "null", "na", "n/a"} else s


def is_interspecific_hybrid_formula(raw_text: Any, allowed_genera: Optional[Set[str]] = None) -> Tuple[bool, str, Optional[str]]:
    """Detect two-epithet hybrid formulas while never flagging single-epithet nothospecies."""
    txt = _clean_str(raw_text)
    if not txt:
        return False, "", None
    allowed = {g.capitalize() for g in (allowed_genera or _DEFAULT_GENERA)}
    # Require full left parent + right epithet around × or spaced x/X so ×memmingeri stays valid.
    for m_sep in re.finditer(r"(?:\s*×\s*|\s+[xX]\s+)", txt):
        left, right = txt[: m_sep.start()].strip(), txt[m_sep.end() :].strip()
        m_left = re.match(r"^([A-Za-z]+)\s+([a-z][a-z-]+)\b", left) if (left and right) else None
        if not m_left:
            continue
        g1, ep1 = m_left.group(1).capitalize(), m_left.group(2).lower().strip("-")
        if g1 not in allowed or len(ep1) < 2 or ep1 in _NON_EPITHET:
            continue
        g2, ep2 = None, None
        m_full, m_abbr, m_bare = re.match(r"^([A-Za-z]+)\s+([a-z][a-z-]+)\b", right), re.match(r"^([A-Za-z])\.\s*([a-z][a-z-]+)\b", right), re.match(r"^([a-z][a-z-]+)(?:\s|$|\()", right)
        if m_full and m_full.group(1).capitalize() in allowed:
            g2, ep2 = m_full.group(1).capitalize(), m_full.group(2).lower().strip("-")
        elif m_abbr:
            hits = [g for g in allowed if g.upper().startswith(m_abbr.group(1).upper())]
            g2, ep2 = (hits[0] if len(hits) == 1 else None), m_abbr.group(2).lower().strip("-")
        elif m_bare:
            g2, ep2 = g1, m_bare.group(1).lower().strip("-")
        if g2 and ep2 and len(ep2) >= 2 and ep2 not in _NON_EPITHET:
            return True, f"{g1} {ep1} × {g2} {ep2}", f"{g1} {ep1}"
    return False, "", None


def has_open_nomenclature_qualifier(raw_name: Any = None, identification_qualifier: Any = None, verbatim_identification: Any = None, identification_remarks: Any = None) -> bool:
    """Return True when fields carry cf./aff./vel aff./?/s.l. or explicit hybrid remarks."""
    iq = _clean_str(identification_qualifier)
    if iq and (_OPEN_QUAL_RE.search(iq) or _HYBRID_REMARK_RE.search(iq) or iq.lower().strip(".") in {"cf", "aff", "vel aff", "s.l", "sensu lato", "?", "hybrid"}):
        return True
    if any((v := _clean_str(f)) and (_OPEN_QUAL_RE.search(v) or _HYBRID_REMARK_RE.search(v)) for f in (raw_name, verbatim_identification)):
        return True
    rem = _clean_str(identification_remarks)
    return bool(rem and (_HYBRID_REMARK_RE.search(rem) or re.search(r"(?:\bs\.l\.|\bsensu\s+lato\b)", rem)))


def _extract_taxon_concept_str(raw_txt: str, allowed: Optional[Set[str]] = None) -> Tuple[Optional[str], str]:
    """Extract canonical 'Genus epithet' and non-autonym var./subsp. rank (ignoring f./forma)."""
    cleaned = re.sub(r"\s+", " ", _OPEN_QUAL_RE.sub(" ", raw_txt)).strip()
    m_bin = re.search(r"\b([A-Z][a-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-z][a-z-]+)\b", cleaned)
    if not m_bin:
        return None, ""
    g, ep = m_bin.group(1).capitalize(), m_bin.group(2).lower().strip("-")
    if len(ep) < 2 or ep in _NON_EPITHET or (allowed is not None and g not in allowed):
        return None, ""
    infra = ""
    if (m_inf := re.search(r"\b(subsp|ssp|var)\.?\s+([a-z][a-z-]+)\b", cleaned[m_bin.end() :], flags=re.IGNORECASE)):
        rk, iep = ("subsp." if m_inf.group(1).lower() in {"subsp", "ssp"} else "var."), m_inf.group(2).lower().strip("-")
        if len(iep) >= 2 and iep not in _NON_EPITHET and iep != ep:  # Strip botanical autonyms.
            infra = f"{rk} {iep}"
    return f"{g} {ep}", infra


def extract_typified_basionym(rec: Mapping[str, Any], is_type_specimen: bool, allowed_genera: Optional[Set[str]] = None) -> str:
    """Resolve verbatim basionym/trinomial on verified types via typifiedName -> typeStatus -> verbatim fields."""
    if not is_type_specimen or not isinstance(rec, Mapping):
        return ""
    allowed = {g.capitalize() for g in (allowed_genera or _DEFAULT_GENERA)}
    cands: List[str] = [s for s in (_clean_str(rec.get("typifiedName")),) if s]
    for tfld in ("typeStatus", "type_status"):
        if (m_of := _TYPE_OF_RE.search(_clean_str(rec.get(tfld)))):
            cands.append(m_of.group(1))
    for entry in (rec.get("identifications") if isinstance(rec.get("identifications"), list) else ()):
        if isinstance(entry, Mapping) and (m_of := _TYPE_OF_RE.search(_clean_str(entry.get("typeStatus")))):
            cands.append(m_of.group(1))
    cands.extend(v for k in ("verbatimScientificName", "verbatimIdentification", "scientificName", "species_raw") if (v := _clean_str(rec.get(k))))
    for cand in cands:
        binom, infra = _extract_taxon_concept_str(cand, allowed=allowed)
        if binom:
            return f"{binom} {infra}".strip()
    return ""


def _parse_year_token(raw_val: Any, max_year: int) -> Optional[int]:
    """Return the first standalone 4-digit year within [1753, max_year], else None."""
    for m_y in re.finditer(r"(?<!\d)(1[7-9]\d{2}|20\d{2}|2100)(?!\d)", _clean_str(raw_val)):
        if 1753 <= (yr := int(m_y.group(1))) <= max_year:
            return yr
    return None


def extract_identification_year(date_identified_raw: Any = None, identified_by_raw: Any = None, identifications: Any = None, max_year: int = 2100) -> Optional[int]:
    """Parse 4-digit determination year (1753..current_year) from dateIdentified, identifications, or identifiedBy."""
    hi = min(int(max_year), datetime.date.today().year)
    if (yr := _parse_year_token(date_identified_raw, hi)) is not None:
        return yr
    if isinstance(identifications, list):
        ordered = sorted([e for e in identifications if isinstance(e, Mapping)], key=lambda e: (1 if e.get("isCurrent") is True else 0, _clean_str(e.get("dateIdentified"))), reverse=True)
        for entry in ordered:
            if (yr := _parse_year_token(entry.get("dateIdentified"), hi) or _parse_year_token(entry.get("identifiedBy"), hi)) is not None:
                return yr
    return _parse_year_token(identified_by_raw, hi)


def _is_homotypic_pair(v_bin: str, r_bin: str, backbone_profiles: Optional[Mapping[str, Mapping[str, Any]]]) -> bool:
    """Return True when verbatim_bin and resolved_bin share a nomenclatural basionym."""
    v_low, r_low = v_bin.lower().strip(), r_bin.lower().strip()
    if v_low == r_low or (v_low, r_low) in _KNOWN_HOMOTYPIC_PAIRS:
        return True
    stem = lambda b: re.sub(r"(?:us|a|um|is|e|ii|i)$", "", re.sub(r"[^a-z-]", "", (b.lower().split() or [""])[-1]))
    if len(v_low.split()) >= 2 and len(r_low.split()) >= 2 and stem(v_low) and stem(v_low) == stem(r_low):
        return True
    if backbone_profiles and r_bin in backbone_profiles:
        prof = backbone_profiles[r_bin]
        if _clean_str(prof.get("basionym_name")).lower() == v_low:
            return True
        if (raw_homo := _clean_str(prof.get("homotypic_synonyms_json"))):
            try:
                return any(_clean_str(x).lower().startswith(v_low) for x in json.loads(raw_homo))
            except Exception:
                pass
    return False


def parse_voucher_taxonomic_concept(
    rec: Mapping[str, Any], resolved_species: Optional[str], is_type_specimen: bool,
    synonym_genera: Set[str], backbone_profiles: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> Dict[str, Any]:
    """Resolve pre-lumping verbatim concept, typified basionym, synonymy status, and open/hybrid flags."""
    allowed = {g.capitalize() for g in (synonym_genera or _DEFAULT_GENERA)}
    id_year = extract_identification_year(rec.get("dateIdentified"), rec.get("identifiedBy"), rec.get("identifications"))
    typ_bas = extract_typified_basionym(rec, is_type_specimen=is_type_specimen, allowed_genera=allowed)
    raw_fields = [v for k in ("verbatimScientificName", "verbatimIdentification", "scientificName", "species_raw") if (v := _clean_str(rec.get(k)))]
    raw_primary = raw_fields[0] if raw_fields else ""
    for cand in raw_fields + [_clean_str(rec.get("identificationQualifier"))]:
        is_hyb, norm_f, p1 = is_interspecific_hybrid_formula(cand, allowed_genera=allowed)
        if is_hyb:
            return {"species_standardized_hint": resolved_species or (p1 if (p1 and p1.startswith("Packera ")) else None), "infraspecific_taxon": "", "typified_basionym": typ_bas, "verbatim_taxon_concept": norm_f, "concept_synonymy_status": "Hybrid_Formula", "identification_year": id_year, "is_open_or_hybrid_concept": True}
    if has_open_nomenclature_qualifier(raw_primary, rec.get("identificationQualifier"), rec.get("verbatimIdentification"), rec.get("identificationRemarks")):
        hint_bin, verb_bin = None, ""
        for cand in raw_fields:
            b, inf = _extract_taxon_concept_str(cand, allowed=allowed)
            if b:
                verb_bin, hint_bin = f"{b} {inf}".strip(), (b if b.startswith("Packera ") else None)
                break
        return {"species_standardized_hint": resolved_species or hint_bin, "infraspecific_taxon": "", "typified_basionym": typ_bas, "verbatim_taxon_concept": verb_bin or raw_primary, "concept_synonymy_status": "Open_Nomenclature", "identification_year": id_year, "is_open_or_hybrid_concept": True}
    src = typ_bas if (is_type_specimen and typ_bas) else raw_primary
    v_bin, infra = _extract_taxon_concept_str(src, allowed=allowed)
    if not v_bin and resolved_species:
        v_bin, infra = _extract_taxon_concept_str(resolved_species, allowed=allowed)
    eff_res = resolved_species or (v_bin if (v_bin and v_bin.startswith("Packera ")) else None)
    status = ("Infraspecific_Taxon" if infra else ("Accepted_Sensu_Stricto" if not (v_bin and eff_res) or v_bin.lower() == eff_res.lower() else ("Homotypic_Synonym" if _is_homotypic_pair(v_bin, eff_res, backbone_profiles) else "Heterotypic_Synonym")))
    return {"species_standardized_hint": eff_res, "infraspecific_taxon": infra, "typified_basionym": typ_bas, "verbatim_taxon_concept": f"{v_bin} {infra}".strip() if v_bin else (resolved_species or ""), "concept_synonymy_status": status, "identification_year": id_year, "is_open_or_hybrid_concept": False}


def _norm_county(val: Any) -> str:
    return re.sub(r"\s+(?:county|parish|borough|co\.?)$", "", _clean_str(val).lower()).strip()


def _norm_state(val: Any) -> str:
    s = _clean_str(val).upper()
    return s if len(s) == 2 else _US_STATE_TO_ABBR.get(s, s)


def _load_segregate_profiles(db_path: Optional[Union[str, Path]], workspace_dir: Optional[Path]) -> Dict[str, Tuple[int, Set[str], Set[str], List[Tuple[float, float]]]]:
    """Load post-1980 segregate publication years alongside county, subnational, and point footprints."""
    segs = {sp: (yr, set(sts), set(cnts), []) for sp, (yr, sts, cnts) in _FALLBACK_POST1980_SEGREGATES.items()}
    cand_db = Path(db_path) if (db_path and str(db_path) != ":memory:") else (Path(workspace_dir) / "data" / "plantmatics.sqlite" if workspace_dir else None)
    if cand_db is None or not cand_db.is_file():
        return segs
    try:
        with sqlite3.connect(str(cand_db)) as conn:
            tbls = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()}
            if "taxon_backbone_profiles" in tbls:
                for sp, pyr in conn.execute("SELECT species_standardized, publication_year FROM taxon_backbone_profiles WHERE publication_year >= 1980;").fetchall():
                    if sp and pyr:
                        ent = segs.setdefault(str(sp), (int(pyr), set(), set(), []))
                        segs[str(sp)] = (min(ent[0], int(pyr)), ent[1], ent[2], ent[3])
            if "usda_county_ranges" in tbls:
                for sp, st, cnt in conn.execute("SELECT canonical_packera, state_norm, county_norm FROM usda_county_ranges;").fetchall():
                    if str(sp) in segs:
                        segs[str(sp)][1].add(_norm_state(st))
                        if _norm_county(cnt):
                            segs[str(sp)][2].add(_norm_county(cnt))
            if "natureserve_subnational_ranges" in tbls:
                for sp, subn in conn.execute("SELECT species_standardized, subnation_code FROM natureserve_subnational_ranges WHERE is_native = 1;").fetchall():
                    if str(sp) in segs and subn:
                        segs[str(sp)][1].add(_norm_state(subn))
            if "vouchers" in tbls:
                for sp, lat, lon in conn.execute("SELECT species_standardized, decimalLatitude, decimalLongitude FROM vouchers WHERE decimalLatitude IS NOT NULL AND decimalLongitude IS NOT NULL;").fetchall():
                    if str(sp) in segs and lat is not None and lon is not None:
                        segs[str(sp)][3].append((float(lat), float(lon)))
    except Exception:
        pass
    return segs


def apply_concept_and_segregate_recency_caps(
    df: "pd.DataFrame", db_path: Optional[Union[str, Path]] = None, workspace_dir: Optional[Path] = None, county_buffer_km: float = 50.0,
) -> "pd.DataFrame":
    """Cap non-type determiner tiers for open/hybrid concepts, heterotypic synonyms, and pre-segregate sheets."""
    if pd is None or df is None or df.empty or "determiner_tier" not in df.columns:
        return df
    out = df.copy()
    segs = _load_segregate_profiles(db_path, workspace_dir)
    for idx in out.index:
        ts_val = _clean_str(out.at[idx, "type_status"] if "type_status" in out.columns else "")
        if ts_val and ts_val.upper() not in {"NONE", "NOT A TYPE", "NON-TYPE"}:
            continue  # Never demote verified nomenclatural types.
        v_oh = out.at[idx, "is_open_or_hybrid_concept"] if "is_open_or_hybrid_concept" in out.columns else False
        c_stat = _clean_str(out.at[idx, "concept_synonymy_status"] if "concept_synonymy_status" in out.columns else "")
        if bool(v_oh is True or str(v_oh).strip().lower() in {"1", "true", "yes"}) or c_stat in {"Open_Nomenclature", "Hybrid_Formula"}:
            out.at[idx, "determiner_tier"] = "Tier_3_Bronze"
            continue
        if _clean_str(out.at[idx, "determiner_tier"]) != "Tier_1_Gold" or not (sp_row := _clean_str(out.at[idx, "species_standardized"] if "species_standardized" in out.columns else out.at[idx, "scientificName"])):
            continue
        if c_stat == "Heterotypic_Synonym":
            out.at[idx, "determiner_tier"] = "Tier_2_Silver"
            continue
        eff_yr = next((int(float(out.at[idx, c])) for c in ("identification_year", "year") if c in out.columns and out.at[idx, c] is not None and not pd.isna(out.at[idx, c]) and str(out.at[idx, c]).strip()), None)
        if eff_yr is None:
            continue
        st_row = _norm_state(out.at[idx, "stateProvince"] if "stateProvince" in out.columns else "")
        cnt_row = _norm_county(out.at[idx, "county"] if "county" in out.columns else "")
        lat_r = out.at[idx, "decimalLatitude"] if "decimalLatitude" in out.columns else None
        lon_r = out.at[idx, "decimalLongitude"] if "decimalLongitude" in out.columns else None
        for seg_sp, (pub_yr, seg_sts, seg_cnts, seg_pts) in segs.items():
            if sp_row == seg_sp or eff_yr >= pub_yr:
                continue
            in_cnt = bool(cnt_row and cnt_row in seg_cnts and (not seg_sts or st_row in seg_sts))
            in_sub = bool(not seg_cnts and st_row and st_row in seg_sts)
            in_buf = bool(seg_pts and lat_r is not None and lon_r is not None and not pd.isna(lat_r) and not pd.isna(lon_r) and any(
                12742.0 * math.asin(min(1.0, math.sqrt(math.sin((math.radians(sl) - math.radians(float(lat_r))) / 2) ** 2 + math.cos(math.radians(float(lat_r))) * math.cos(math.radians(sl)) * math.sin((math.radians(sn) - math.radians(float(lon_r))) / 2) ** 2))) <= county_buffer_km
                for sl, sn in seg_pts
            ))
            if in_cnt or in_sub or in_buf:
                out.at[idx, "determiner_tier"] = "Tier_2_Silver"
                break
    return out
