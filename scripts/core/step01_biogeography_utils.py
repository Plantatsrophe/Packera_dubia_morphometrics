"""
===============================================================================
Module: step01_biogeography_utils.py
Project: Plantmatics Integrative Morphometrics and Phylogenomics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Config-driven Eastern North American taxon residency resolution,
    orthographic canonicalization, and ecoregion-bounded transcontinental
    biogeographic screening helpers for Step 01 ingestion and cohort pruning.
===============================================================================
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import sqlite3
from typing import Any, FrozenSet, Optional, Set, Tuple, Union

from scripts.core.config import (
    CIS_MISSISSIPPI_STATES_PROVINCES,
    GREAT_PLAINS_STRADDLING_STATES,
    PROJECT_ROOT,
    WESTERN_LONGITUDE_THRESHOLD,
    _DEFAULT_CONFIG,
    resolve_target_genus,
)
from scripts.core.db import resolve_db_path
from scripts.core.epa_ecoregions import (
    DEFAULT_CEC_NA_L3_GEOJSON,
    DEFAULT_EPA_L3_GEOJSON,
    _is_canadian_province_label,
    _is_epa_download_mocked,
    _point_in_epa_l3,
    infer_regional_group,
    load_cec_na_level3_spatial_index,
    standardize_epa_l3_label,
)
from scripts.core.taxonomy_authority import (
    _clean_upper_state_province,
    is_excluded_western_region,
)

EASTERN_NA_TDWG_L3_CODES: FrozenSet[str] = frozenset({
    "ALA", "CNT", "DEL", "DWC", "FLO", "GEO", "ILL", "INI", "KTY", "MAI",
    "MAS", "MRY", "MIC", "MSI", "NCA", "NWH", "NWJ", "NWY", "OHI", "PEN",
    "RHO", "SCA", "TEN", "VER", "VRG", "WDC", "WIS", "WVA", "LAB", "NBR",
    "NFL", "NSC", "ONT", "PEI", "QUE",
})

PLAINS_CORDILLERAN_ECOREGIONS: FrozenSet[str] = frozenset({
    "Northwestern_Great_Plains", "Northwestern_Glaciated_Plains",
    "Northern_Glaciated_Plains", "Western_High_Plains", "Central_Great_Plains",
    "South_Central_Semi_Arid_Prairies", "Texas_Blackland_Prairies",
    "Edwards_Plateau", "East_Central_Texas_Plains", "Southern_Texas_Plains",
    "Cross_Timbers", "High_Plains", "Aspen_Parkland_Northern_Glaciated_Plains",
    "Boreal_Plain", "Canadian_Rockies", "Middle_Rockies", "Southern_Rockies",
    "Wyoming_Basin", "Colorado_Plateaus", "Chihuahuan_Deserts",
})

_DEFAULT_SYNONYM_GENERA: FrozenSet[str] = frozenset({"Packera", "Senecio", "Cineraria"})
_DEFAULT_TRANSCONTINENTAL_TAXA: FrozenSet[str] = frozenset(
    {"Packera paupercula", "Packera indecora", "Packera pauciflora"}
)
_CIS_MS_UPPER_SET: Set[str] = {str(s).strip().upper() for s in CIS_MISSISSIPPI_STATES_PROVINCES}


def normalize_packera_orthography(name: Any) -> str:
    """Normalize orthographic and hybrid-marker variants to canonical binomial form."""
    if name is None:
        return ""
    txt = re.sub(r"\s+", " ", str(name).strip())
    if not txt or txt.lower() in {"nan", "none", "null", "na"}:
        return ""
    low = txt.lower()
    if low == "packera millefolia":
        return "Packera millefolium"
    if re.fullmatch(r"packera\s*(?:×|x)?\s*memmingeri", low):
        return "Packera memmingeri"
    return txt


def _canonical_taxon_equiv_key(name: Any) -> str:
    """Return a lowercase canonical matching key normalizing millefolia and hybrid memmingeri."""
    norm = normalize_packera_orthography(name).lower()
    if not norm:
        return ""
    m_hyb = re.match(r"^([a-z]+)(?:\s*×\s*|\s+[x×]\s+|\s+)([a-z][a-z-]*)$", norm)
    return f"{m_hyb.group(1)} {m_hyb.group(2)}" if m_hyb else norm


def resolve_synonym_genera(cfg: Optional[Any] = None) -> Set[str]:
    """Read configured `taxa.synonym_genera`, always including `resolve_target_genus(cfg)`."""
    active_cfg = cfg if cfg is not None else _DEFAULT_CONFIG
    raw_list = getattr(getattr(active_cfg, "taxa", None), "synonym_genera", None) or []
    genera = {str(g).strip().capitalize() for g in raw_list if str(g).strip()} or set(_DEFAULT_SYNONYM_GENERA)
    target_genus = (resolve_target_genus(active_cfg) or "").strip().capitalize()
    if target_genus:
        genera.add(target_genus)
    return genera


def resolve_transcontinental_taxa(cfg: Optional[Any] = None) -> Set[str]:
    """Read configured `taxa.transcontinental_boreal_taxa` with canonical fallback."""
    active_cfg = cfg if cfg is not None else _DEFAULT_CONFIG
    raw_list = getattr(getattr(active_cfg, "taxa", None), "transcontinental_boreal_taxa", None) or []
    taxa = {normalize_packera_orthography(t) for t in raw_list if normalize_packera_orthography(t)}
    return taxa if taxa else set(_DEFAULT_TRANSCONTINENTAL_TAXA)


def _read_sqlite_eastern_residents(
    db_path: Optional[Union[str, Path]] = None,
    workspace_root: Optional[Path] = None,
) -> Set[str]:
    """Query Tier A SQLite taxonomic backbone and NatureServe tables read-only when present."""
    found: Set[str] = set()
    eff_ws = Path(workspace_root).resolve() if workspace_root is not None else PROJECT_ROOT
    try:
        resolved = resolve_db_path(db_path=db_path, workspace_root=eff_ws)
    except Exception:
        return found
    if not resolved or str(resolved) == ":memory:":
        return found
    db_file = Path(resolved)
    if not db_file.is_file() or db_file.stat().st_size == 0:
        return found

    conn: Optional[sqlite3.Connection] = None
    try:
        conn = sqlite3.connect(f"file:{db_file}?mode=ro", uri=True)
    except sqlite3.Error:
        try:
            conn = sqlite3.connect(str(db_file))
        except sqlite3.Error:
            return found

    try:
        tables = {str(r[0]) for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()}
        if "taxon_backbone_profiles" in tables:
            cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(taxon_backbone_profiles);").fetchall()}
            if {"species_standardized", "tdwg_native_codes_json"}.issubset(cols):
                for sp_raw, tdwg_json in conn.execute(
                    "SELECT species_standardized, tdwg_native_codes_json FROM taxon_backbone_profiles;"
                ).fetchall():
                    if not sp_raw or not tdwg_json:
                        continue
                    try:
                        codes = json.loads(str(tdwg_json))
                    except Exception:
                        codes = [c.strip() for c in str(tdwg_json).split(",") if c.strip()]
                    if {str(c).strip().upper() for c in (codes or []) if str(c).strip()} & EASTERN_NA_TDWG_L3_CODES:
                        sp_norm = normalize_packera_orthography(sp_raw)
                        if sp_norm:
                            found.add(sp_norm)

        if "natureserve_subnational_ranges" in tables:
            cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(natureserve_subnational_ranges);").fetchall()}
            if {"species_standardized", "subnation_code", "is_native"}.issubset(cols):
                for sp_raw, sub_code in conn.execute(
                    "SELECT species_standardized, subnation_code FROM natureserve_subnational_ranges WHERE is_native = 1;"
                ).fetchall():
                    if sp_raw and sub_code and str(sub_code).strip().upper() in _CIS_MS_UPPER_SET:
                        sp_norm = normalize_packera_orthography(sp_raw)
                        if sp_norm:
                            found.add(sp_norm)
    except sqlite3.Error:
        pass
    finally:
        if conn is not None:
            conn.close()
    return found


def resolve_eastern_resident_whitelist(
    cfg: Optional[Any] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_root: Optional[Path] = None,
) -> Set[str]:
    """
    Union `taxa.eastern_resident_whitelist` from config with native Eastern NA taxa
    documented in SQLite `taxon_backbone_profiles` (TDWG Level 3) and
    `natureserve_subnational_ranges` (`is_native = 1`), normalizing orthography.
    """
    active_cfg = cfg if cfg is not None else _DEFAULT_CONFIG
    taxa_cfg = getattr(active_cfg, "taxa", None)
    raw_wl = list(getattr(taxa_cfg, "eastern_resident_whitelist", None) or [])
    if not raw_wl and _DEFAULT_CONFIG is not active_cfg:
        raw_wl = list(getattr(_DEFAULT_CONFIG.taxa, "eastern_resident_whitelist", None) or [])

    whitelist: Set[str] = {norm for sp in raw_wl for norm in [normalize_packera_orthography(sp)] if norm}
    whitelist.update(_read_sqlite_eastern_residents(db_path=db_path, workspace_root=workspace_root))

    cfg_syns = getattr(taxa_cfg, "synonyms", None) or {}
    if isinstance(cfg_syns, dict):
        wl_keys = {_canonical_taxon_equiv_key(s) for s in whitelist}
        for canon_k, syn_vals in cfg_syns.items():
            c_norm = normalize_packera_orthography(canon_k)
            if c_norm and _canonical_taxon_equiv_key(c_norm) in wl_keys:
                whitelist.add(c_norm)
                for sv in syn_vals or []:
                    s_clean = re.sub(r"\s+", " ", str(sv).strip())
                    if s_clean:
                        whitelist.add(s_clean)
    return whitelist


def _resolve_ecoregion_offline_safe(
    lat: Optional[float],
    lon: Optional[float],
    state_province: Optional[str],
    regional_group: Optional[str],
) -> Optional[str]:
    """Resolve EPA/CEC Level III ecoregion without triggering multi-hundred-MB downloads."""
    eff_region = standardize_epa_l3_label(regional_group)
    if eff_region:
        return eff_region
    if DEFAULT_CEC_NA_L3_GEOJSON.is_file() or _is_epa_download_mocked():
        return infer_regional_group(lat, lon, state_province=state_province)
    if DEFAULT_EPA_L3_GEOJSON.is_file():
        hit = infer_regional_group(lat, lon, state_province=state_province, ecoregion_path=DEFAULT_EPA_L3_GEOJSON)
        if hit is not None:
            return hit
    if lat is not None and lon is not None:
        try:
            lat_f, lon_f = float(lat), float(lon)
            if math.isfinite(lat_f) and math.isfinite(lon_f):
                if lat_f >= 41.6 or _is_canadian_province_label(state_province):
                    cec_bundle = load_cec_na_level3_spatial_index(allow_download=False)
                    return _point_in_epa_l3(lon_f, lat_f, cec_bundle)
        except (ValueError, TypeError):
            pass
    return None


def evaluate_transcontinental_ecoregion_gate(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    regional_group: Optional[str] = None,
    is_type_specimen: bool = False,
    western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
) -> Tuple[bool, str]:
    """
    Screen transcontinental boreal taxa using EPA/CEC Level III ecoregion membership
    instead of an arbitrary -95.5W meridian cliff.
    """
    if is_type_specimen:
        return False, "type_exempt"
    if is_excluded_western_region(
        state_province,
        lat=lat,
        lon=lon,
        western_longitude_threshold=western_longitude_threshold,
        exclude_western=True,
        is_type_specimen=False,
    ):
        return True, "excluded_plains_cordilleran"

    eff_region = _resolve_ecoregion_offline_safe(
        lat=lat, lon=lon, state_province=state_province, regional_group=regional_group
    )
    up_state = _clean_upper_state_province(state_province)
    if eff_region in PLAINS_CORDILLERAN_ECOREGIONS or (
        eff_region is None and up_state in GREAT_PLAINS_STRADDLING_STATES
    ):
        return True, "excluded_plains_cordilleran"

    if lon is not None:
        try:
            lon_f = float(lon)
            if math.isfinite(lon_f) and lon_f < -95.5:
                return False, "boreal_ecoregion_extension"
        except (ValueError, TypeError):
            pass
    return False, "core_eastern"
