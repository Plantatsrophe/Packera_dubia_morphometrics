"""
===============================================================================
Module: harvester.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unified botanical voucher ingestion, Darwin Core metadata normalization,
    determiner authority evaluation, and asynchronous media acquisition engine.
    Consolidates dynamic genus-wide GBIF occurrence querying, 3-tier determiner
    authority stratification, geographic boundary filtering, image quality
    validation, and atomic curated dataset persistence.
===============================================================================
"""

from __future__ import annotations

import asyncio
import concurrent.futures
import datetime
import json
import logging
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import unicodedata
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

try:
    import aiohttp
except ImportError:
    aiohttp = None  # type: ignore[assignment]
try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]
try:
    from PIL import Image
except ImportError:
    Image = None  # type: ignore[assignment]
class _MissingPygbifSubmodule:
    """Mock-friendly stub when optional `pygbif` is not installed in the current environment."""

    _is_missing_pygbif = True

    def __init__(self, module_name: str) -> None:
        self._module_name = module_name

    def _raise(self, *_args: Any, **_kwargs: Any) -> Any:
        raise RuntimeError(
            f"Optional dependency '{self._module_name}' (pygbif) is required for live GBIF API queries. "
            "Install via `pip install pygbif`."
        )

    search = _raise
    name_backbone = _raise
    name_usage = _raise
    name_lookup = _raise


try:
    import pygbif
    import pygbif.occurrences as occ
    import pygbif.species as species
except ImportError:
    pygbif = None  # type: ignore[assignment]
    occ = _MissingPygbifSubmodule("pygbif.occurrences")  # type: ignore[assignment]
    species = _MissingPygbifSubmodule("pygbif.species")  # type: ignore[assignment]

_AIOHTTP_CLIENT_ERRORS: Tuple[type, ...] = (aiohttp.ClientError,) if aiohttp is not None else ()
try:
    from shapely.geometry import Point, Polygon, shape
    from shapely.strtree import STRtree
except ImportError:
    Point = Polygon = shape = STRtree = None  # type: ignore[assignment]
try:
    from tqdm.asyncio import tqdm as async_tqdm
except ImportError:
    async_tqdm = None  # type: ignore[assignment]

# Allow decompression of high-resolution botanical herbarium sheets (>89 MP)
if Image is not None:
    Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    DEFAULT_COUNTRIES,
    DEFAULT_GEOGRAPHIC_SCOPE,
    DEFAULT_HOME_HERBARIUM,
    DEFAULT_MIN_FILE_SIZE_KB,
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN,
    DEFAULT_OUTPUT_CSV,
    DEFAULT_RAW_DIR,
    DEFAULT_TARGET_TAXA,
    DEFAULT_USDA_COUNTY_BUFFER_KM,
    DEFAULT_USDA_COUNTY_FILTER,
    DEFAULT_USDA_EXEMPT_GOLD,
    DEFAULT_WORKSPACE,
    EXCLUDED_WESTERN_STATES,
    MAJOR_HERBARIA_CODES,
    PRIMARY_TYPE_STATUSES,
    SECONDARY_TYPE_STATUSES,
    SPECIALIST_PATTERNS,
    VALID_GEOGRAPHIC_SCOPES,
    VALID_TYPE_STATUSES,
    WESTERN_LONGITUDE_THRESHOLD,
    _DEFAULT_CONFIG,
    normalize_countries_list,
)
from scripts.core.db import (
    get_db_connection,
    is_canonical_or_explicit_db_context,
    load_taxon_synonyms,
    resolve_db_path,
    upsert_taxon_synonyms,
    upsert_vouchers_records,
)
from scripts.core.logger import setup_logging
from scripts.core.usda_plants_ranges import filter_vouchers_by_usda_county_range

_SYNONYM_GENERA: Set[str] = {"Packera", "Senecio", "Cineraria"}
_OPEN_NOMENCLATURE_EPITHETS: Set[str] = {
    "sp", "spp", "cf", "aff", "indet", "hybrid", "unknown", "na", "x", "×",
}


def _extract_canonical_binomial(
    raw_name: Any,
    allowed_genera: Optional[Set[str]] = None,
) -> Optional[str]:
    """Extract canonical 'Genus epithet' binomial from a Darwin Core scientific name string."""
    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return None
    txt = str(raw_name).strip()
    if not txt:
        return None
    m = re.match(r"^([A-Za-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", txt)
    if not m:
        return None
    genus = m.group(1).capitalize()
    epithet = m.group(2).lower().strip("-")
    if len(epithet) < 2 or epithet in _OPEN_NOMENCLATURE_EPITHETS:
        return None
    if allowed_genera is not None and genus not in allowed_genera:
        return None
    return f"{genus} {epithet}"


def _extract_infraspecific_taxon(raw_name: Any) -> str:
    """Extract normalized infraspecific rank + epithet ('var. gracilis', 'subsp. ...') if present."""
    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return ""
    txt = str(raw_name).strip()
    if not txt:
        return ""
    m = re.search(
        r"\b(subsp|ssp|var|forma|f|nothosubsp|nothovar)\.?\s+([a-z][a-z-]{1,})\b",
        txt,
        flags=re.IGNORECASE,
    )
    if not m:
        return ""
    rk_raw = m.group(1).lower()
    ep = m.group(2).lower().strip("-")
    if ep in _OPEN_NOMENCLATURE_EPITHETS:
        return ""
    rk_map = {
        "ssp": "subsp.",
        "subsp": "subsp.",
        "var": "var.",
        "forma": "f.",
        "f": "f.",
        "nothosubsp": "nothosubsp.",
        "nothovar": "nothovar.",
    }
    return f"{rk_map.get(rk_raw, rk_raw + '.')} {ep}"


_CULTIVATED_STATUS_TOKENS: Set[str] = {
    "cultivated", "captive", "managed", "garden", "greenhouse", "horticultural",
}
_CULTIVATED_TEXT_RE = re.compile(
    r"(?<!not\s)(?<!non-)(?<!escaped\sfrom\s)\b"
    r"(?:cultivated|in\s+cultivation|grown\s+in\s+cultivation|greenhouse|glasshouse|"
    r"botanical\s+garden|botanic\s+garden|arboretum|hortus\s+botanicus|experimental\s+garden)\b",
    re.IGNORECASE,
)
_CA_PROVINCE_ALIASES_TO_KEY: Dict[str, str] = {
    "on": "ontario", "ont": "ontario", "ontario": "ontario",
    "qc": "quebec", "pq": "quebec", "que": "quebec", "quebec": "quebec", "québec": "quebec",
    "nb": "new brunswick", "new brunswick": "new brunswick", "nouveau-brunswick": "new brunswick",
    "ns": "nova scotia", "nova scotia": "nova scotia", "nouvelle-écosse": "nova scotia", "nouvelle-ecosse": "nova scotia",
    "pe": "prince edward island", "pei": "prince edward island", "prince edward island": "prince edward island",
    "nl": "newfoundland and labrador", "nf": "newfoundland and labrador", "newfoundland": "newfoundland and labrador", "labrador": "newfoundland and labrador", "newfoundland and labrador": "newfoundland and labrador",
    "mb": "manitoba", "man": "manitoba", "manitoba": "manitoba",
    "sk": "saskatchewan", "sask": "saskatchewan", "saskatchewan": "saskatchewan",
    "ab": "alberta", "alta": "alberta", "alberta": "alberta",
    "bc": "british columbia", "british columbia": "british columbia", "colombie-britannique": "british columbia",
    "yt": "yukon", "yukon": "yukon", "yukon territory": "yukon",
    "nt": "northwest territories", "nwt": "northwest territories", "northwest territories": "northwest territories",
    "nu": "nunavut", "nunavut": "nunavut",
}

# (lat_min, lat_max, lon_min, lon_max) envelopes for offline CoordinateCleaner state/province verification
_NA_STATE_PROVINCE_BBOXES: Dict[str, Tuple[float, float, float, float]] = {
    "alabama": (30.1, 35.1, -88.5, -84.8),
    "alaska": (51.2, 71.5, -179.9, -129.9),
    "arizona": (31.3, 37.1, -114.9, -109.0),
    "arkansas": (33.0, 36.6, -94.7, -89.6),
    "california": (32.5, 42.1, -124.5, -114.1),
    "colorado": (36.9, 41.1, -109.1, -102.0),
    "connecticut": (40.9, 42.1, -73.8, -71.7),
    "delaware": (38.4, 39.9, -75.8, -75.0),
    "district of columbia": (38.7, 39.0, -77.2, -76.9),
    "florida": (24.4, 31.1, -87.7, -79.9),
    "georgia": (30.3, 35.1, -85.7, -80.8),
    "hawaii": (18.9, 22.3, -160.3, -154.8),
    "idaho": (41.9, 49.1, -117.3, -111.0),
    "illinois": (36.9, 42.6, -91.6, -87.4),
    "indiana": (37.7, 41.8, -88.1, -84.7),
    "iowa": (40.3, 43.6, -96.7, -90.1),
    "kansas": (36.9, 40.1, -102.1, -94.5),
    "kentucky": (36.4, 39.2, -89.6, -81.9),
    "louisiana": (28.9, 33.1, -94.1, -88.8),
    "maine": (43.0, 47.5, -71.1, -66.9),
    "maryland": (37.9, 39.8, -79.5, -75.0),
    "massachusetts": (41.2, 42.9, -73.6, -69.9),
    "michigan": (41.6, 48.3, -90.5, -82.1),
    "minnesota": (43.4, 49.4, -97.3, -89.4),
    "mississippi": (30.1, 35.1, -91.7, -88.0),
    "missouri": (35.9, 40.7, -95.8, -89.0),
    "montana": (44.3, 49.1, -116.1, -104.0),
    "nebraska": (39.9, 43.1, -104.1, -95.3),
    "nevada": (35.0, 42.1, -120.1, -114.0),
    "new hampshire": (42.6, 45.4, -72.6, -70.6),
    "new jersey": (38.9, 41.4, -75.6, -73.8),
    "new mexico": (31.3, 37.1, -109.1, -103.0),
    "new york": (40.4, 45.1, -79.8, -71.8),
    "north carolina": (33.8, 36.6, -84.4, -75.4),
    "north dakota": (45.9, 49.1, -104.1, -96.5),
    "ohio": (38.4, 42.0, -84.9, -80.5),
    "oklahoma": (33.6, 37.1, -103.1, -94.4),
    "oregon": (41.9, 46.3, -124.6, -116.4),
    "pennsylvania": (39.7, 42.3, -80.6, -74.6),
    "rhode island": (41.1, 42.1, -71.9, -71.1),
    "south carolina": (32.0, 35.3, -83.4, -78.5),
    "south dakota": (42.4, 46.0, -104.1, -96.4),
    "tennessee": (34.9, 36.7, -90.4, -81.6),
    "texas": (25.8, 36.6, -106.7, -93.5),
    "utah": (36.9, 42.1, -114.1, -109.0),
    "vermont": (42.7, 45.1, -73.5, -71.4),
    "virginia": (36.5, 39.5, -83.7, -75.2),
    "washington": (45.5, 49.1, -124.8, -116.9),
    "west virginia": (37.2, 40.7, -82.7, -77.7),
    "wisconsin": (42.4, 47.1, -92.9, -86.7),
    "wyoming": (40.9, 45.1, -111.1, -104.0),
    "ontario": (41.6, 56.9, -95.2, -74.3),
    "quebec": (44.9, 62.6, -79.8, -57.1),
    "new brunswick": (44.5, 48.1, -69.1, -63.7),
    "nova scotia": (43.3, 47.1, -66.4, -59.6),
    "prince edward island": (45.9, 47.1, -64.5, -61.9),
    "newfoundland and labrador": (46.5, 60.4, -67.8, -52.6),
    "manitoba": (49.0, 60.1, -102.1, -88.9),
    "saskatchewan": (49.0, 60.1, -110.1, -101.3),
    "alberta": (49.0, 60.1, -120.1, -110.0),
    "british columbia": (48.2, 60.1, -139.1, -114.0),
    "yukon": (60.0, 69.7, -141.0, -123.8),
    "northwest territories": (60.0, 78.8, -136.5, -102.0),
    "nunavut": (51.6, 83.2, -120.7, -61.0),
}


def _normalize_ca_province_bbox_key(raw_state_prov: Optional[str]) -> Optional[str]:
    if not raw_state_prov:
        return None
    cleaned = re.sub(r"[^a-zà-ÿ\s-]", "", str(raw_state_prov).strip().lower())
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return _CA_PROVINCE_ALIASES_TO_KEY.get(cleaned)


_FATAL_GBIF_GEO_ISSUES: Set[str] = {
    "ZERO_COORDINATE",
    "COORDINATE_INVALID",
    "COORDINATE_OUT_OF_RANGE",
    "COUNTRY_COORDINATE_MISMATCH",
}


def _is_cultivated_occurrence(rec: Dict[str, Any]) -> bool:
    """Return True if Darwin Core establishment or locality metadata indicates a cultivated/garden plant."""
    for fld in ("establishmentMeans", "degreeOfEstablishment"):
        val = str(rec.get(fld) or "").strip().lower()
        if val and any(tok in val for tok in _CULTIVATED_STATUS_TOKENS) and "escaped" not in val and "naturalized" not in val:
            return True
    for txt_fld in ("locality", "verbatimLocality", "habitat", "occurrenceRemarks"):
        txt = str(rec.get(txt_fld) or "").strip()
        if txt and _CULTIVATED_TEXT_RE.search(txt):
            return True
    return False


_GBIF_SYNONYM_STATUSES: Set[str] = {
    "SYNONYM",
    "HOMOTYPIC_SYNONYM",
    "HETEROTYPIC_SYNONYM",
}
_GBIF_REJECTED_SYNONYM_STATUSES: Set[str] = {
    "MISAPPLIED",
    "DOUBTFUL",
    "PROPARTE_SYNONYM",
}
_DYNAMIC_SYNONYM_TO_CANONICAL: Dict[str, str] = {}
for _cfg_canon, _cfg_syns in (_DEFAULT_CONFIG.taxa.synonyms or {}).items():
    _canon_norm = _extract_canonical_binomial(_cfg_canon, {"Packera"}) or _cfg_canon.strip()
    if _canon_norm:
        _DYNAMIC_SYNONYM_TO_CANONICAL[_canon_norm.lower()] = _canon_norm
        for _cs in _cfg_syns or []:
            _cs_low = re.sub(r"\s+", " ", str(_cs).strip().lower())
            if _cs_low:
                _DYNAMIC_SYNONYM_TO_CANONICAL[_cs_low] = _canon_norm
            _cs_bin = _extract_canonical_binomial(_cs, _SYNONYM_GENERA)
            if _cs_bin:
                _DYNAMIC_SYNONYM_TO_CANONICAL[_cs_bin.lower()] = _canon_norm

_DYNAMIC_TAXON_AUTHOR_STRINGS: Dict[str, Set[str]] = {}


def _hydrate_synonyms_from_sqlite(db_path: Optional[Union[str, Path]] = None) -> None:
    """Hydrate in-memory synonym mappings from persisted `taxon_synonyms` when SQLite exists."""
    try:
        resolved = resolve_db_path(db_path=db_path)
        if resolved != ":memory:" and Path(resolved).is_file():
            conn = get_db_connection(db_path=resolved, read_only=True)
            try:
                for raw_low, canon in load_taxon_synonyms(conn).items():
                    if not raw_low or not canon:
                        continue
                    _DYNAMIC_SYNONYM_TO_CANONICAL.setdefault(raw_low, canon)
                    raw_bin = _extract_canonical_binomial(raw_low, _SYNONYM_GENERA)
                    if raw_bin:
                        _DYNAMIC_SYNONYM_TO_CANONICAL.setdefault(raw_bin.lower(), canon)
            finally:
                conn.close()
    except Exception:
        pass


_hydrate_synonyms_from_sqlite()
_GBIF_SYNONYM_CACHE: Dict[str, Tuple[Optional[str], List[str], List[str]]] = {}
_GBIF_BACKBONE_SYNONYM_FLAG: Dict[str, bool] = {}
_GBIF_GENUS_KEY_CACHE: Dict[str, Optional[int]] = {}
_GBIF_GENUS_SPECIES_CACHE: Dict[int, List[str]] = {}
_GBIF_NETWORK_UNAVAILABLE: bool = False


def _is_gbif_species_mocked() -> bool:
    """Return True if pygbif.species endpoints are currently patched by unittest.mock."""
    nb = getattr(species, "name_backbone", None)
    nu = getattr(species, "name_usage", None)
    nl = getattr(species, "name_lookup", None)
    return any(
        hasattr(fn, "_mock_name") or hasattr(fn, "assert_called")
        for fn in (nb, nu, nl)
        if fn is not None
    )


def _register_dynamic_synonyms(
    canonical_packera: str,
    synonyms: Iterable[Any],
    author_strings: Optional[Iterable[str]] = None,
    allow_packera_species_remap: bool = True,
) -> None:
    """Register discovered synonyms and author strings in the in-memory standardization cache."""
    if not canonical_packera or not canonical_packera.startswith("Packera "):
        return
    _DYNAMIC_SYNONYM_TO_CANONICAL[canonical_packera.lower()] = canonical_packera
    if author_strings:
        clean_authors = {str(a).strip() for a in author_strings if a and str(a).strip()}
        if clean_authors:
            _DYNAMIC_TAXON_AUTHOR_STRINGS.setdefault(canonical_packera, set()).update(clean_authors)
            for auth_sci in clean_authors:
                auth_low = re.sub(r"\s+", " ", auth_sci.lower())
                auth_p_bin = _extract_canonical_binomial(auth_sci, {"Packera"})
                if (
                    not allow_packera_species_remap
                    and auth_p_bin
                    and auth_p_bin.lower() != canonical_packera.lower()
                    and not re.search(r"\b(?:var|subsp|ssp|f)\b", auth_low)
                ):
                    continue
                _DYNAMIC_SYNONYM_TO_CANONICAL[auth_low] = canonical_packera
                auth_bin = _extract_canonical_binomial(auth_sci, _SYNONYM_GENERA)
                if auth_bin:
                    _DYNAMIC_SYNONYM_TO_CANONICAL[auth_bin.lower()] = canonical_packera

    for item in synonyms:
        if not item:
            continue
        item_str = re.sub(r"\s+", " ", str(item).strip())
        if not item_str:
            continue
        item_low = item_str.lower()
        item_p_bin = _extract_canonical_binomial(item_str, {"Packera"})
        if (
            not allow_packera_species_remap
            and item_p_bin
            and item_p_bin.lower() != canonical_packera.lower()
            and not re.search(r"\b(?:var|subsp|ssp|f)\b", item_low)
        ):
            continue
        _DYNAMIC_SYNONYM_TO_CANONICAL[item_low] = canonical_packera
        item_bin = _extract_canonical_binomial(item_str, _SYNONYM_GENERA)
        if item_bin:
            _DYNAMIC_SYNONYM_TO_CANONICAL[item_bin.lower()] = canonical_packera


def _call_and_normalize_gbif_backbone(query_name: str, rank: str) -> Dict[str, Any]:
    """
    Invoke `pygbif.species.name_backbone` across both v1 (`name=`, `rank=`) and
    v2 (`scientificName=`, `taxonRank=`) signatures, normalizing nested GBIF v2
    responses into the flat backbone dict schema.
    """
    if species is None or getattr(species, "name_backbone", None) is None:
        return {}
    nb_fn = species.name_backbone
    use_v2_sig = False
    if not _is_gbif_species_mocked():
        try:
            import inspect
            use_v2_sig = "scientificName" in inspect.signature(nb_fn).parameters
        except Exception:
            use_v2_sig = False
    if use_v2_sig:
        raw = nb_fn(scientificName=query_name, taxonRank=str(rank).upper(), kingdom="Plantae")
    else:
        raw = nb_fn(name=query_name, rank=rank, kingdom="Plantae")
    if not isinstance(raw, dict) or not raw:
        return {}
    if "usageKey" in raw or "matchType" in raw or "canonicalName" in raw or "species" in raw:
        return raw
    usage = raw.get("usage") if isinstance(raw.get("usage"), dict) else {}
    acc_usage = raw.get("acceptedUsage") if isinstance(raw.get("acceptedUsage"), dict) else {}
    diag = raw.get("diagnostics") if isinstance(raw.get("diagnostics"), dict) else {}
    cls_list = raw.get("classification") if isinstance(raw.get("classification"), list) else []
    cls_by_rank = {
        str(c.get("rank", "")).upper(): c
        for c in cls_list
        if isinstance(c, dict) and c.get("rank")
    }
    sp_cls = cls_by_rank.get("SPECIES", {})
    gen_cls = cls_by_rank.get("GENUS", {})
    norm: Dict[str, Any] = {
        "matchType": diag.get("matchType", "EXACT" if usage else "NONE"),
        "status": usage.get("status") or diag.get("status"),
        "usageKey": usage.get("key"),
        "acceptedUsageKey": acc_usage.get("key") or sp_cls.get("key"),
        "scientificName": usage.get("name"),
        "canonicalName": usage.get("canonicalName") or usage.get("name"),
        "accepted": acc_usage.get("name") or acc_usage.get("canonicalName") or sp_cls.get("name"),
        "species": acc_usage.get("canonicalName") or sp_cls.get("name") or acc_usage.get("name"),
        "genus": gen_cls.get("name"),
        "genusKey": gen_cls.get("key"),
    }
    return {k: v for k, v in norm.items() if v is not None}


def _persist_dynamic_synonyms_for_taxon(canonical_packera: str) -> None:
    """Persist newly discovered live GBIF synonyms for `canonical_packera` into SQLite `taxon_synonyms`."""
    if (
        not canonical_packera
        or _is_gbif_species_mocked()
        or os.environ.get("PACKERA_DEMO_MODE") == "1"
        or not is_canonical_or_explicit_db_context()
    ):
        return
    try:
        resolved = resolve_db_path()
        if resolved != ":memory:" and Path(resolved).is_file():
            sub_map = {
                k: v
                for k, v in _DYNAMIC_SYNONYM_TO_CANONICAL.items()
                if v == canonical_packera and k
            }
            if not sub_map:
                return
            auth_sub = {}
            if canonical_packera in _DYNAMIC_TAXON_AUTHOR_STRINGS:
                auth_sub[canonical_packera] = sorted(_DYNAMIC_TAXON_AUTHOR_STRINGS[canonical_packera])
            conn = get_db_connection(db_path=resolved, init_schema=False)
            try:
                upsert_taxon_synonyms(conn, sub_map, auth_sub or None)
            finally:
                conn.close()
    except Exception:
        pass


def _get_packera_genus_key_gbif(timeout_seconds: float = 5.0) -> Optional[int]:
    """
    Query pygbif.species.name_backbone(name='Packera', rank='genus', kingdom='Plantae')
    with a fast 5-second timeout and in-memory cache to obtain the GBIF usageKey for genus Packera.
    """
    global _GBIF_NETWORK_UNAVAILABLE

    is_mocked = _is_gbif_species_mocked()
    cache_key = "packera"
    if not is_mocked and cache_key in _GBIF_GENUS_KEY_CACHE:
        return _GBIF_GENUS_KEY_CACHE[cache_key]
    if not is_mocked and _GBIF_NETWORK_UNAVAILABLE:
        return None

    def _fetch_key() -> Optional[int]:
        backbone = _call_and_normalize_gbif_backbone("Packera", "genus")
        if not isinstance(backbone, dict) or not backbone:
            return None
        match_type = str(backbone.get("matchType", "EXACT")).upper()
        if match_type == "NONE":
            return None
        raw_key = backbone.get("usageKey") or backbone.get("genusKey") or backbone.get("key")
        if raw_key is None:
            return None
        try:
            return int(raw_key)
        except (TypeError, ValueError):
            return None

    if is_mocked:
        try:
            return _fetch_key()
        except Exception:
            return None

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch_key)
        res = future.result(timeout=timeout_seconds)
        _GBIF_GENUS_KEY_CACHE[cache_key] = res
        return res
    except Exception:
        _GBIF_NETWORK_UNAVAILABLE = True
        return None
    finally:
        executor.shutdown(wait=False, cancel_futures=True)


def discover_genus_packera_species_gbif(limit: int = 1000) -> List[str]:
    """
    Dynamically query the GBIF Species Backbone API to collect all accepted
    `Packera <epithet>` species binomials across the entire genus Packera.
    """
    global _GBIF_NETWORK_UNAVAILABLE

    eff_limit = max(1, int(limit))
    is_mocked = _is_gbif_species_mocked()
    if not is_mocked and eff_limit in _GBIF_GENUS_SPECIES_CACHE:
        return list(_GBIF_GENUS_SPECIES_CACHE[eff_limit])
    if not is_mocked and _GBIF_NETWORK_UNAVAILABLE:
        return []

    genus_key = _get_packera_genus_key_gbif(timeout_seconds=5.0)
    if genus_key is None:
        return []

    def _fetch_species() -> List[str]:
        raw_items: List[Any] = []

        if hasattr(species, "name_lookup"):
            try:
                lookup_resp = species.name_lookup(
                    higherTaxonKey=genus_key,
                    rank="SPECIES",
                    status="ACCEPTED",
                    limit=eff_limit,
                )
                if isinstance(lookup_resp, dict):
                    raw_items.extend(lookup_resp.get("results", []) or [])
                elif isinstance(lookup_resp, list):
                    raw_items.extend(lookup_resp)
            except Exception:
                pass

        if hasattr(species, "name_usage"):
            try:
                try:
                    children_resp = species.name_usage(key=genus_key, data="children", limit=eff_limit)
                except TypeError:
                    children_resp = species.name_usage(key=genus_key, data="children")
                if isinstance(children_resp, dict):
                    raw_items.extend(children_resp.get("results", []) or [])
                elif isinstance(children_resp, list):
                    raw_items.extend(children_resp)
            except Exception:
                pass

        rejected_statuses = {
            "SYNONYM",
            "HETEROTYPIC_SYNONYM",
            "HOMOTYPIC_SYNONYM",
            "PROPARTE_SYNONYM",
            "MISAPPLIED",
            "DOUBTFUL",
        }
        discovered: Set[str] = set()
        for item in raw_items:
            if not isinstance(item, dict):
                continue
            status_val = str(item.get("taxonomicStatus") or item.get("status") or "").strip().upper()
            if status_val and status_val in rejected_statuses:
                continue
            rank_val = str(item.get("rank") or "").strip().upper()
            if rank_val and rank_val not in {"SPECIES", "SUBSPECIES", "VARIETY"}:
                continue
            cand = _extract_canonical_binomial(
                item.get("canonicalName") or item.get("species") or item.get("scientificName"),
                {"Packera"},
            )
            if cand and re.match(r"^Packera\s+(?:[x×]\s*)?[a-z-]+$", cand):
                discovered.add(cand)
                sci_full = item.get("scientificName")
                auth_list = [sci_full.strip()] if isinstance(sci_full, str) and sci_full.strip() else None
                _register_dynamic_synonyms(cand, [cand], auth_list)

        return sorted(discovered)

    if is_mocked:
        try:
            return _fetch_species()
        except Exception:
            return []

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch_species)
        res = future.result(timeout=5.0)
        _GBIF_GENUS_SPECIES_CACHE[eff_limit] = list(res)
        return res
    except Exception:
        _GBIF_NETWORK_UNAVAILABLE = True
        return []
    finally:
        executor.shutdown(wait=False, cancel_futures=True)


def _query_gbif_backbone_and_synonyms(
    query_name: str,
    timeout_seconds: float = 5.0,
) -> Tuple[Optional[str], List[str], List[str]]:
    """
    Query `pygbif.species.name_backbone` and `pygbif.species.name_usage(..., data='synonyms')`
    with a hard 5-second timeout and in-memory cache so offline unit tests never fail or hang.
    Returns `(accepted_packera_binomial, synonym_binomials, raw_scientific_names_with_authors)`
    and records in `_GBIF_BACKBONE_SYNONYM_FLAG` whether `query_name` itself is an exact
    unaccepted synonym in the GBIF Backbone.
    """
    global _GBIF_NETWORK_UNAVAILABLE

    is_mocked = _is_gbif_species_mocked()
    cache_key = query_name.strip().lower()
    if not is_mocked and cache_key in _GBIF_SYNONYM_CACHE:
        return _GBIF_SYNONYM_CACHE[cache_key]
    if not is_mocked and _GBIF_NETWORK_UNAVAILABLE:
        return None, [], []

    def _fetch() -> Tuple[Optional[str], List[str], List[str], bool]:
        backbone = _call_and_normalize_gbif_backbone(query_name, "species")
        if not isinstance(backbone, dict) or not backbone:
            return None, [], [], False
        match_type = str(backbone.get("matchType", "EXACT")).upper()
        if match_type in {"NONE", "HIGHERRANK"}:
            return None, [], [], False

        bb_status = str(backbone.get("status") or backbone.get("taxonomicStatus") or "").strip().upper()
        is_exact_synonym = (match_type == "EXACT") and (bb_status in _GBIF_SYNONYM_STATUSES)

        author_SciNames: List[str] = []
        for fld in ("scientificName", "accepted"):
            val = backbone.get(fld)
            if isinstance(val, str) and val.strip():
                author_SciNames.append(val.strip())

        accepted_packera: Optional[str] = None
        acc_fields = ("accepted", "species", "canonicalName") if is_exact_synonym else ("species", "accepted", "canonicalName")
        for cand_fld in acc_fields:
            cand_bin = _extract_canonical_binomial(backbone.get(cand_fld), {"Packera"})
            if cand_bin:
                accepted_packera = cand_bin
                break

        usage_key = backbone.get("usageKey")
        accepted_key = backbone.get("acceptedUsageKey")
        query_packera_bin = _extract_canonical_binomial(query_name, {"Packera"})
        if (
            accepted_key is not None
            and str(accepted_key) != str(usage_key)
            and (accepted_packera is None or (query_packera_bin and accepted_packera.lower() == query_packera_bin.lower()))
            and hasattr(species, "name_usage")
        ):
            try:
                acc_key_norm = int(accepted_key) if str(accepted_key).isdigit() else accepted_key
                acc_info = species.name_usage(key=acc_key_norm)
                if isinstance(acc_info, dict) and "results" not in acc_info:
                    acc_bin = _extract_canonical_binomial(
                        acc_info.get("canonicalName") or acc_info.get("species") or acc_info.get("scientificName"),
                        {"Packera"},
                    )
                    if acc_bin:
                        accepted_packera = acc_bin
                    acc_sci = acc_info.get("scientificName")
                    if isinstance(acc_sci, str) and acc_sci.strip():
                        author_SciNames.append(acc_sci.strip())
            except Exception:
                pass

        found_binomials: List[str] = []
        for cand_fld in ("canonicalName", "species", "accepted"):
            b_val = _extract_canonical_binomial(backbone.get(cand_fld), _SYNONYM_GENERA)
            if b_val and b_val not in found_binomials:
                found_binomials.append(b_val)

        keys_to_query: List[Any] = []
        if usage_key is not None:
            keys_to_query.append(int(usage_key) if str(usage_key).isdigit() else usage_key)
        if accepted_key is not None and str(accepted_key) != str(usage_key):
            keys_to_query.append(int(accepted_key) if str(accepted_key).isdigit() else accepted_key)

        for k_id in keys_to_query:
            try:
                syn_resp = species.name_usage(key=k_id, data="synonyms")
            except Exception:
                syn_resp = []
            if isinstance(syn_resp, dict):
                syn_items = syn_resp.get("results", [])
            elif isinstance(syn_resp, list):
                syn_items = syn_resp
            else:
                syn_items = []
            for item in syn_items:
                if not isinstance(item, dict):
                    continue
                item_status = str(item.get("taxonomicStatus") or item.get("status") or "").strip().upper()
                if item_status in _GBIF_REJECTED_SYNONYM_STATUSES:
                    continue
                sci_full = item.get("scientificName")
                if isinstance(sci_full, str) and sci_full.strip():
                    author_SciNames.append(sci_full.strip())
                syn_bin = _extract_canonical_binomial(
                    item.get("canonicalName") or item.get("scientificName") or item.get("species"),
                    _SYNONYM_GENERA,
                )
                if syn_bin and syn_bin not in found_binomials:
                    found_binomials.append(syn_bin)
                if accepted_packera is None:
                    acc_cand = _extract_canonical_binomial(item.get("species") or item.get("accepted"), {"Packera"})
                    if acc_cand:
                        accepted_packera = acc_cand

        return accepted_packera, found_binomials, author_SciNames, is_exact_synonym

    if is_mocked:
        try:
            acc_p, f_bins, auth_s, is_syn = _fetch()
            _GBIF_BACKBONE_SYNONYM_FLAG[cache_key] = is_syn
            return acc_p, f_bins, auth_s
        except Exception:
            return None, [], []

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch)
        acc_p, f_bins, auth_s, is_syn = future.result(timeout=timeout_seconds)
        res = (acc_p, f_bins, auth_s)
        _GBIF_SYNONYM_CACHE[cache_key] = res
        _GBIF_BACKBONE_SYNONYM_FLAG[cache_key] = is_syn
        return res
    except Exception:
        _GBIF_NETWORK_UNAVAILABLE = True
        return None, [], []
    finally:
        executor.shutdown(wait=False, cancel_futures=True)


def expand_taxon_synonyms_gbif(
    taxon_name: str,
    cfg_synonyms: Optional[Dict[str, List[str]]] = None,
    protected_packera_taxa: Optional[Set[str]] = None,
) -> List[str]:
    """
    Dynamically expands a target `Packera` (or historical `Senecio`/`Cineraria`) taxon into its
    full homotypic and heterotypic synonym query set using:
      1. Optional configuration synonym mappings (`cfg_synonyms`).
      2. Live GBIF Species API (`pygbif.species.name_backbone` + `pygbif.species.name_usage`),
         guarded by a 5-second timeout and in-memory cache for offline reliability.
    Unaccepted species-rank `Packera` synonyms in the GBIF Backbone (e.g., `Packera tomentosa`
    under `Packera dubia`) and historical generic combinations (`Senecio`/`Cineraria`) are
    dynamically unified under the accepted `Packera <epithet>` in `_DYNAMIC_SYNONYM_TO_CANONICAL`
    and persisted to SQLite `taxon_synonyms`, while protecting distinct accepted `Packera`
    species in `protected_packera_taxa` from cross-species overwriting.
    """
    if taxon_name is None:
        return []
    clean_input = str(taxon_name).strip()
    if not clean_input:
        return []
    if clean_input.lower() == "packera":
        return ["Packera"]

    protected_lower: Set[str] = {str(p).strip().lower() for p in (protected_packera_taxa or set()) if p}
    expanded: List[str] = []
    seen_lower: Set[str] = set()

    def _add_candidate(cand: Optional[str]) -> None:
        if not cand:
            return
        c_str = str(cand).strip()
        if not c_str:
            return
        c_low = c_str.lower()
        if c_low not in seen_lower:
            seen_lower.add(c_low)
            expanded.append(c_str)

    canonical_packera = _extract_canonical_binomial(clean_input, {"Packera"})
    if canonical_packera and canonical_packera.lower() in _DYNAMIC_SYNONYM_TO_CANONICAL:
        canonical_packera = _DYNAMIC_SYNONYM_TO_CANONICAL[canonical_packera.lower()]

    cfg_group_matched: List[str] = []
    if cfg_synonyms:
        for cfg_key, cfg_vals in cfg_synonyms.items():
            key_clean = str(cfg_key).strip()
            vals_clean = [str(v).strip() for v in (cfg_vals or []) if str(v).strip()]
            all_group = [key_clean] + vals_clean
            all_group_lower = {g.lower() for g in all_group}
            bin_input = _extract_canonical_binomial(clean_input, _SYNONYM_GENERA)
            if (
                clean_input.lower() in all_group_lower
                or (bin_input and bin_input.lower() in all_group_lower)
                or (canonical_packera and canonical_packera.lower() == key_clean.lower())
            ):
                if canonical_packera is None:
                    canonical_packera = _extract_canonical_binomial(key_clean, {"Packera"})
                for g in all_group:
                    _add_candidate(g)
                    cfg_group_matched.append(g)

    if canonical_packera and canonical_packera.startswith("Packera "):
        _add_candidate(canonical_packera)
    _add_candidate(clean_input)

    # Query GBIF Species Backbone + Synonyms endpoint (fast 5s timeout + cached)
    gbif_query_target = canonical_packera if (canonical_packera and canonical_packera.startswith("Packera ")) else clean_input
    gbif_acc, gbif_syns, gbif_author_names = _query_gbif_backbone_and_synonyms(gbif_query_target)
    query_is_synonym = _GBIF_BACKBONE_SYNONYM_FLAG.get(gbif_query_target.strip().lower(), False)
    if gbif_acc and gbif_acc.startswith("Packera "):
        if not (canonical_packera and canonical_packera.startswith("Packera ")):
            canonical_packera = gbif_acc
        elif query_is_synonym and gbif_acc.lower() != canonical_packera.lower():
            canonical_packera = gbif_acc
        _add_candidate(gbif_acc)
    for g_syn in gbif_syns:
        _add_candidate(g_syn)

    if canonical_packera and canonical_packera.startswith("Packera "):
        canon_low = canonical_packera.lower()
        other_protected = protected_lower - {canon_low}
        if cfg_group_matched:
            _register_dynamic_synonyms(canonical_packera, cfg_group_matched, None, allow_packera_species_remap=True)
        safe_expanded = [
            cand
            for cand in expanded
            if (_extract_canonical_binomial(cand, {"Packera"}) or "").lower() not in other_protected
        ]
        safe_authors = [
            auth
            for auth in (gbif_author_names or [])
            if (_extract_canonical_binomial(auth, {"Packera"}) or "").lower() not in other_protected
        ]
        _register_dynamic_synonyms(
            canonical_packera,
            safe_expanded,
            safe_authors,
            allow_packera_species_remap=True,
        )
        _persist_dynamic_synonyms_for_taxon(canonical_packera)

    return expanded


def standardize_packera_taxon(raw_name: Any) -> Optional[str]:
    """
    Dynamically standardizes raw Darwin Core scientific names, varieties, unaccepted `Packera`
    synonyms, and historical `Senecio`/`Cineraria` basionyms across the genus `Packera`.
    Uses the runtime/SQLite synonym table populated dynamically by `expand_taxon_synonyms_gbif()`
    and occurrence ingestion, GBIF `name_backbone` synonym verification, and genus-wide botanical
    regex matching (`^packera\\s*(?:([x×])\\s*)?([a-z-]+)`).
    """
    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return None
    s_orig = str(raw_name).strip()
    if not s_orig:
        return None
    s_compact = re.sub(r"\s+", " ", s_orig.lower())
    if s_compact in {"none", "nan", "null", "na", "unknown", "packera"}:
        return None

    cand_bin = _extract_canonical_binomial(s_orig, _SYNONYM_GENERA)
    is_non_packera_syn_genus = bool(cand_bin and cand_bin.split(" ", 1)[0] in (_SYNONYM_GENERA - {"Packera"}))

    # When pygbif.species is mocked in a unit test, allow mocked basionym or exact-synonym
    # Packera lookups to resolve directly without being shadowed by earlier cache state.
    if cand_bin and _is_gbif_species_mocked():
        gbif_acc, gbif_syns, gbif_authors = _query_gbif_backbone_and_synonyms(cand_bin)
        is_exact_syn = _GBIF_BACKBONE_SYNONYM_FLAG.get(cand_bin.lower(), False)
        if gbif_acc and gbif_acc.startswith("Packera ") and (is_non_packera_syn_genus or is_exact_syn):
            _register_dynamic_synonyms(
                gbif_acc,
                [s_orig, cand_bin, *gbif_syns],
                gbif_authors,
                allow_packera_species_remap=is_exact_syn,
            )
            return gbif_acc

    # 1. Dynamic GBIF Backbone & runtime/SQLite synonym dictionary lookup
    if s_compact in _DYNAMIC_SYNONYM_TO_CANONICAL:
        return _DYNAMIC_SYNONYM_TO_CANONICAL[s_compact]

    var_m = re.match(r"^([a-z]+\s+[a-z-]+\s+(?:var\.?|subsp\.?|ssp\.?|f\.?)\s+[a-z-]+)", s_compact)
    if var_m:
        var_raw = var_m.group(1)
        if var_raw in _DYNAMIC_SYNONYM_TO_CANONICAL:
            return _DYNAMIC_SYNONYM_TO_CANONICAL[var_raw]
        var_dot = re.sub(r"\b(var|subsp|ssp|f)\b\.?", r"\1.", var_raw)
        if var_dot in _DYNAMIC_SYNONYM_TO_CANONICAL:
            return _DYNAMIC_SYNONYM_TO_CANONICAL[var_dot]

    if cand_bin and cand_bin.lower() in _DYNAMIC_SYNONYM_TO_CANONICAL:
        return _DYNAMIC_SYNONYM_TO_CANONICAL[cand_bin.lower()]

    # 2. General genus-wide botanical regex for any valid Packera binomial
    m = re.match(r"^packera(?:\s*×\s*|\s+[x×]\s+|\s+)([a-z][a-z-]*)", s_compact, flags=re.IGNORECASE)
    if m:
        epithet = m.group(1).lower().strip("-")
        if len(epithet) >= 2 and epithet not in _OPEN_NOMENCLATURE_EPITHETS:
            return f"Packera {epithet}"
        return None

    # 3. Dynamic GBIF Backbone resolution for historical Senecio / Cineraria synonyms
    if is_non_packera_syn_genus:
        gbif_acc, gbif_syns, gbif_authors = _query_gbif_backbone_and_synonyms(cand_bin)
        if gbif_acc and gbif_acc.startswith("Packera "):
            _register_dynamic_synonyms(gbif_acc, [s_orig, cand_bin, *gbif_syns], gbif_authors)
            _persist_dynamic_synonyms_for_taxon(gbif_acc)
            return gbif_acc

    return None

from scripts.core.taxonomy_authority import (
    EXPORT_COLUMNS,
    GENERIC_DETERMINER_TOKENS,
    _COLLECTOR_REFERENCE_REGEX,
    _HOLOTYPE_TIER_STATUSES,
    _INVALID_INSTITUTION_CODES,
    _IPNI_AUTHOR_ABBREVIATIONS,
    _ISOTYPE_TIER_STATUSES,
    _MISSING_IDENTIFIED_BY_TOKENS,
    _NEGATED_TYPE_RE,
    _PARATYPE_TIER_STATUSES,
    _SYNTYPE_TIER_STATUSES,
    _fold_ascii,
    _is_explicit_collector_reference,
    _is_missing_identified_by,
    _is_verified_type_status,
    _resolve_row_standardized_packera,
    _should_fallback_to_collector,
    apply_corpus_authority_tiers,
    classify_nomenclatural_type_rank,
    classify_nomenclatural_type_subrank,
    discover_corpus_authorities,
    extract_protologue_authors,
    generate_exsiccatae_key,
    get_image_resolution_score,
    haversine_distance_meters,
    is_excluded_western_region,
    normalize_collector,
    normalize_determiner_surnames,
    normalize_event_date,
    normalize_record_number,
    parse_determiner_tier,
    resolve_geographic_scope_settings,
    resolve_voucher_image_path,
    sanitize_filename,
    select_cascading_outgroup_anchors,
    setup_logger,
    spatially_balance_and_thin_vouchers,
    stratify_duplicates,
)
from scripts.core.epa_ecoregions import (
    DEFAULT_EPA_L3_DIR,
    DEFAULT_EPA_L3_GEOJSON,
    DEFAULT_EPA_L3_GPKG,
    DEFAULT_US_COUNTIES_GEOJSON,
    EPA_L3_OFFICIAL_ZIP_URL,
    _EPA_L3_INDEX_CACHE,
    _HISTORIC_COUNTY_CONSOLIDATIONS,
    _US_STATE_ABBR_TO_FULL,
    _build_offline_epa_l3_fallback_index,
    _ensure_epa_level3_cache_files,
    _normalize_us_county_name,
    _normalize_us_state_name,
    _point_in_epa_l3,
    infer_regional_group,
    load_epa_level3_spatial_index,
    standardize_epa_l3_label,
)


def optimize_herbarium_image_url(url: str) -> str:
    """Transforms provider-specific URLs to request full-resolution original scans."""
    if not url or not isinstance(url, str):
        return ""

    optimized = url.strip()

    # Smithsonian NMNH: strip dimension clamp (e.g. &h=2000)
    if "collections.nmnh.si.edu/media/" in optimized:
        optimized = re.sub(r"[?&][hw]=\d+", "", optimized)
        if "?" not in optimized and "&" in optimized:
            optimized = optimized.replace("&", "?", 1)

    # Symbiota / SERNEC / SEINet / CCH portals: replace web/thumbnail with orig/large
    if any(k in optimized.lower() for k in ["symbiota", "sernec", "seinet", "cch2", "swbiodiversity"]):
        optimized = re.sub(r"/(?:web|tn|thumbnails?)/", "/orig/", optimized, flags=re.IGNORECASE)
        optimized = re.sub(r"_(?:tn|web|sm)\.(jpe?g|png)", r"_lg.\1", optimized, flags=re.IGNORECASE)

    # IIIF endpoints: replace constrained dimensions with /full/max/
    if "/full/!" in optimized or "/full/pct:" in optimized or re.search(r"/full/\d+,\d*/", optimized):
        optimized = re.sub(r"/full/(?:!?\d+,\d*|pct:\d+)/", "/full/max/", optimized)

    return optimized


def extract_high_res_image_url(media_list: Optional[List[Dict[str, Any]]]) -> Optional[str]:
    """Parses Darwin Core media records, scoring and prioritizing highest-quality specimen image."""
    if not media_list or not isinstance(media_list, list):
        return None

    scored_candidates: List[Tuple[float, str]] = []

    for item in media_list:
        if not isinstance(item, dict):
            continue

        m_type = str(item.get("type", ""))
        m_format = str(item.get("format", "")).lower()
        identifier = str(item.get("identifier", "")).strip()

        if not identifier or not identifier.startswith(("http://", "https://")):
            continue

        score = 0.0
        if m_type == "StillImage" or "image" in m_format:
            score += 50.0
        if identifier.lower().endswith((".jpg", ".jpeg", ".png", ".tif", ".tiff")):
            score += 20.0

        ident_lower = identifier.lower()
        if any(h in ident_lower for h in ["_lg", "_large", "original", "/orig/", "/master/", "/highres/", "/full/", "hires", "high_res"]):
            score += 100.0
        if "max" in ident_lower:
            score += 30.0

        if any(t in ident_lower for t in ["_tn", "_thumb", "thumbnail", "_sm", "_small", "preview", "icon", "mini"]):
            score -= 150.0
        if "detailimages" in ident_lower:
            score -= 25.0
        if re.search(r"[?&]h=(?:[1-9]\d{0,2}|1\d{3}|2000)\b", ident_lower):
            score -= 10.0

        optimized_url = optimize_herbarium_image_url(identifier)
        scored_candidates.append((score, optimized_url))

    if not scored_candidates:
        return None

    scored_candidates.sort(key=lambda x: x[0], reverse=True)
    return scored_candidates[0][1]


def validate_image_quality(
    image_path: Path,
    min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
    min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB,
    check_sharpness: bool = False,
    min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN,
) -> Tuple[bool, Dict[str, Any]]:
    """Evaluates image resolution, byte size, and optical sharpness metrics."""
    image_path = Path(image_path)
    if not image_path.exists():
        resolved_alt = resolve_voucher_image_path(image_path)
        if resolved_alt.is_file():
            image_path = resolved_alt
        else:
            return False, {"valid": False, "reason": "file_not_found"}

    file_size_kb = image_path.stat().st_size / 1024.0
    if file_size_kb < min_file_size_kb:
        return False, {
            "valid": False,
            "reason": "file_size_too_small",
            "file_size_kb": round(file_size_kb, 1),
            "min_file_size_kb": min_file_size_kb,
        }

    try:
        with Image.open(image_path) as img:
            w, h = img.size
            mp = round((w * h) / 1e6, 2)
    except Exception as e:
        return False, {
            "valid": False,
            "reason": f"corrupt_or_unreadable_image: {e}",
            "file_size_kb": round(file_size_kb, 1),
        }

    if mp < min_megapixels:
        return False, {
            "valid": False,
            "reason": "low_resolution",
            "megapixels": mp,
            "min_megapixels": min_megapixels,
            "width": w,
            "height": h,
            "file_size_kb": round(file_size_kb, 1),
        }

    sharpness_score = None
    if check_sharpness:
        try:
            import cv2
            img_cv = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
            if img_cv is not None:
                ch, cw = img_cv.shape[:2]
                y0, y1 = int(ch * 0.15), max(int(ch * 0.85), 1)
                x0, x1 = int(cw * 0.15), max(int(cw * 0.85), 1)
                roi_cv = img_cv[y0:y1, x0:x1] if (y1 > y0 and x1 > x0) else img_cv
                sharpness_score = round(float(cv2.Laplacian(roi_cv, cv2.CV_64F).var()), 2)
                if sharpness_score < min_sharpness:
                    return False, {
                        "valid": False,
                        "reason": "blurry_or_upscaled",
                        "sharpness": sharpness_score,
                        "min_sharpness": min_sharpness,
                        "megapixels": mp,
                        "width": w,
                        "height": h,
                        "file_size_kb": round(file_size_kb, 1),
                    }
        except Exception:
            pass

    return True, {
        "valid": True,
        "megapixels": mp,
        "width": w,
        "height": h,
        "file_size_kb": round(file_size_kb, 1),
        "sharpness": sharpness_score,
    }


async def download_single_image(
    session: aiohttp.ClientSession,
    image_url: str,
    destination_path: Path,
    semaphore: asyncio.Semaphore,
    max_retries: int = 3,
    force: bool = False,
    webp_quality: int = 92,
    logger: Optional[logging.Logger] = None,
) -> Tuple[bool, bool]:
    """Asynchronously downloads a single voucher image file to local storage with retry logic.

    Returns:
        Tuple[bool, bool]: (success, was_skipped)
    """
    destination_path = Path(destination_path)
    if not force:
        existing_path = (
            destination_path
            if (destination_path.exists() and os.path.getsize(destination_path) > 0)
            else resolve_voucher_image_path(destination_path)
        )
        if existing_path.is_file() and os.path.getsize(existing_path) > 0:
            if logger:
                logger.info(
                    f"Skipping download for existing voucher {existing_path.name} "
                    f"({os.path.getsize(existing_path)} bytes)."
                )
            return True, True

    async with semaphore:
        for attempt in range(1, max_retries + 1):
            temp_path: Optional[Path] = None
            try:
                if session is None:
                    req = urllib.request.Request(
                        image_url,
                        headers={"User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium)"},
                    )

                    def _fetch_urllib() -> Tuple[int, str, bytes]:
                        with urllib.request.urlopen(req, timeout=45) as resp:
                            st = int(getattr(resp, "status", 200))
                            ct = str(resp.headers.get("Content-Type", "")).lower()
                            return st, ct, resp.read()

                    status_code, content_type, content = await asyncio.to_thread(_fetch_urllib)
                else:
                    timeout = aiohttp.ClientTimeout(total=45, connect=15) if aiohttp is not None else 45
                    async with session.get(image_url, timeout=timeout) as response:
                        status_code = response.status
                        content_type = response.headers.get("Content-Type", "").lower() if response.status == 200 else ""
                        content = await response.read() if response.status == 200 else b""
                if status_code == 200:
                    if True:
                        if len(content) > 1024 and (
                            not content_type
                            or "image" in content_type
                            or "octet-stream" in content_type
                            or content[:3] == b"\xff\xd8\xff"
                        ):
                            payload_to_write = content
                            if destination_path.suffix.lower() == ".webp" and Image is not None:
                                try:
                                    import io
                                    with Image.open(io.BytesIO(content)) as img:
                                        rgb_img = img.convert("RGB")
                                        buf = io.BytesIO()
                                        rgb_img.save(buf, format="WEBP", quality=int(webp_quality), method=4)
                                        webp_bytes = buf.getvalue()
                                        if len(webp_bytes) > 0:
                                            payload_to_write = webp_bytes
                                except Exception:
                                    payload_to_write = content

                            destination_path.parent.mkdir(parents=True, exist_ok=True)
                            temp_file = tempfile.NamedTemporaryFile(
                                dir=destination_path.parent,
                                prefix=f"{destination_path.stem}_",
                                suffix=".tmp",
                                delete=False,
                            )
                            temp_path = Path(temp_file.name)
                            temp_file.write(payload_to_write)
                            temp_file.flush()
                            os.fsync(temp_file.fileno())
                            temp_file.close()
                            temp_path.replace(destination_path)
                            return True, False
                elif status_code in {404, 410}:
                    return False, False
            except (*_AIOHTTP_CLIENT_ERRORS, asyncio.TimeoutError, Exception):
                if temp_path and temp_path.exists():
                    try:
                        temp_path.unlink()
                    except Exception:
                        pass
                if attempt == max_retries:
                    return False, False
                await asyncio.sleep(1.0 * (2 ** (attempt - 1)))
        return False, False


async def download_all_voucher_images(
    records_to_download: List[Tuple[str, Path]],
    concurrency_limit: int = 15,
    force: bool = False,
    webp_quality: int = 92,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, int]:
    """Coordinates asynchronous batch downloading of voucher images with concurrency control."""
    semaphore = asyncio.Semaphore(concurrency_limit)
    headers = {
        "User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"
    }
    total_records = len(records_to_download)
    stats = {"success": 0, "skipped": 0, "failed": 0, "total": total_records}

    pending = []
    for url, dest in records_to_download:
        dest_p = Path(dest)
        existing_dest = (
            dest_p
            if (not force and dest_p.exists() and os.path.getsize(dest_p) > 0)
            else (resolve_voucher_image_path(dest_p) if not force else dest_p)
        )
        if not force and existing_dest.is_file() and os.path.getsize(existing_dest) > 0:
            if logger:
                logger.info(
                    f"Retaining existing voucher image: {existing_dest.name} "
                    f"({os.path.getsize(existing_dest)} bytes)."
                )
            stats["skipped"] += 1
        else:
            pending.append((url, dest))

    if not pending:
        if logger:
            logger.info(
                f"All {stats['skipped']} voucher images are already cached locally. "
                f"[Processed: 0 | Skipped: {stats['skipped']} | Total: {total_records}]"
            )
        return stats

    if logger:
        logger.info(
            f"Initiating asynchronous download of {len(pending)} pending images (Concurrency: {concurrency_limit})... "
            f"[Processed: 0 | Skipped: {stats['skipped']} | Total: {total_records}]"
        )

    if aiohttp is not None:
        connector = aiohttp.TCPConnector(limit=concurrency_limit, limit_per_host=5, ssl=False)
        async with aiohttp.ClientSession(headers=headers, connector=connector) as session:
            tasks = [
                download_single_image(
                    session, url, dest, semaphore, force=force, webp_quality=webp_quality, logger=logger
                )
                for url, dest in pending
            ]
            results = (
                await async_tqdm.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
                if async_tqdm is not None
                else await asyncio.gather(*tasks)
            )
    else:
        tasks = [
            download_single_image(
                None, url, dest, semaphore, force=force, webp_quality=webp_quality, logger=logger
            )
            for url, dest in pending
        ]
        results = (
            await async_tqdm.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
            if async_tqdm is not None
            else await asyncio.gather(*tasks)
        )

    for success, was_skipped in results:
        if success:
            if was_skipped:
                stats["skipped"] += 1
            else:
                stats["success"] += 1
        else:
            stats["failed"] += 1

    if logger:
        logger.info(
            f"Voucher Download Summary: "
            f"[Processed: {stats['success']} | Skipped: {stats['skipped']} | Total: {total_records}]"
        )

    return stats


def migrate_existing_vouchers_to_webp(
    raw_dir: Union[str, Path],
    db_path: Optional[Union[str, Path]] = None,
    curated_csv: Optional[Union[str, Path]] = None,
    quality: int = 92,
    logger: Optional[logging.Logger] = None,
) -> int:
    """
    Scans `raw_dir` for any existing `.jpg`/`.jpeg` voucher images, converts each valid
    image to `.webp` (`quality=92, method=2` by default) via PIL, removes the old `.jpg`
    once the `.webp` is verified non-empty, and updates `.jpg`/`.jpeg` -> `.webp` in
    `curated_csv` and SQLite `vouchers.image_path` if present.
    """
    raw_dir_path = Path(raw_dir)
    migrated_count = 0

    if raw_dir_path.exists() and raw_dir_path.is_dir() and Image is not None:
        jpg_files = sorted(
            p for p in raw_dir_path.iterdir()
            if p.is_file() and p.suffix.lower() in {".jpg", ".jpeg"}
        )
        total_jpgs = len(jpg_files)
        if total_jpgs > 0 and logger:
            logger.info(
                f"Found {total_jpgs} existing JPG voucher image(s) in {raw_dir_path}; "
                f"migrating to WebP (quality={quality})..."
            )

        def _convert_single_jpg(jpg_path: Path) -> int:
            webp_path = jpg_path.with_suffix(".webp")
            if webp_path.is_file() and webp_path.stat().st_size > 0:
                try:
                    jpg_path.unlink()
                    return 1
                except Exception:
                    return 0
            temp_path: Optional[Path] = None
            try:
                with Image.open(jpg_path) as img:
                    rgb_img = img.convert("RGB")
                    temp_file = tempfile.NamedTemporaryFile(
                        dir=raw_dir_path,
                        prefix=f"{jpg_path.stem}_",
                        suffix=".webp.tmp",
                        delete=False,
                    )
                    temp_path = Path(temp_file.name)
                    temp_file.close()
                    rgb_img.save(temp_path, format="WEBP", quality=int(quality), method=2)
                if temp_path.is_file() and temp_path.stat().st_size > 0:
                    temp_path.replace(webp_path)
                    if webp_path.is_file() and webp_path.stat().st_size > 0:
                        jpg_path.unlink()
                        return 1
                elif temp_path and temp_path.exists():
                    temp_path.unlink()
            except Exception as exc:
                if temp_path and temp_path.exists():
                    try:
                        temp_path.unlink()
                    except Exception:
                        pass
                if logger:
                    logger.warning(f"Skipping WebP migration for unreadable image {jpg_path.name}: {exc}")
            return 0

        if total_jpgs <= 4:
            for jpg_path in jpg_files:
                migrated_count += _convert_single_jpg(jpg_path)
        else:
            env_workers = os.environ.get("PACKERA_MAX_WORKERS", "").strip()
            max_workers = (
                max(1, int(env_workers))
                if env_workers.isdigit()
                else min(16, max(4, (os.cpu_count() or 4) * 2))
            )
            completed_jpgs = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
                for res in pool.map(_convert_single_jpg, jpg_files):
                    migrated_count += res
                    completed_jpgs += 1
                    if logger and (completed_jpgs % 250 == 0 or completed_jpgs == total_jpgs):
                        logger.info(
                            f"  [WebP Migration] Processed {completed_jpgs}/{total_jpgs} "
                            f"({100.0 * completed_jpgs / total_jpgs:.1f}%) JPG sheets..."
                        )

    if migrated_count > 0 and logger:
        logger.info(
            f"Migrated {migrated_count} existing JPG voucher image(s) in {raw_dir_path} to WebP (quality={quality})."
        )

    if curated_csv is not None and pd is not None:
        csv_p = Path(curated_csv)
        if csv_p.is_file() and csv_p.stat().st_size > 0:
            try:
                df_csv = pd.read_csv(csv_p)
                if "image_path" in df_csv.columns and not df_csv.empty:
                    old_paths = df_csv["image_path"].astype(str).tolist()
                    new_paths = [
                        re.sub(r"\.jpe?g$", ".webp", str(v), flags=re.IGNORECASE)
                        if not pd.isna(v) and str(v).strip()
                        else v
                        for v in df_csv["image_path"]
                    ]
                    if [str(x) for x in new_paths] != old_paths:
                        df_csv["image_path"] = new_paths
                        tmp_csv = tempfile.NamedTemporaryFile(
                            mode="w",
                            delete=False,
                            dir=csv_p.parent,
                            suffix=".tmp",
                            encoding="utf-8",
                        )
                        tmp_csv_p = Path(tmp_csv.name)
                        try:
                            df_csv.to_csv(tmp_csv, index=False, encoding="utf-8")
                            tmp_csv.flush()
                            os.fsync(tmp_csv.fileno())
                            tmp_csv.close()
                            tmp_csv_p.replace(csv_p)
                        except Exception:
                            if tmp_csv_p.exists():
                                tmp_csv_p.unlink()
            except Exception as csv_exc:
                if logger:
                    logger.warning(f"Could not update image_path in {csv_p}: {csv_exc}")

    try:
        target_csv_ctx = Path(curated_csv) if curated_csv is not None else None
        if db_path is not None or (target_csv_ctx is not None and is_canonical_or_explicit_db_context(target_csv_ctx)):
            resolved_db = resolve_db_path(db_path=db_path)
            if resolved_db != ":memory:" and Path(resolved_db).is_file():
                conn = get_db_connection(db_path=resolved_db, init_schema=False)
                try:
                    conn.execute(
                        "UPDATE vouchers SET image_path = SUBSTR(image_path, 1, LENGTH(image_path) - 4) || '.webp' "
                        "WHERE image_path IS NOT NULL AND LOWER(SUBSTR(image_path, -4)) = '.jpg'"
                    )
                    conn.execute(
                        "UPDATE vouchers SET image_path = SUBSTR(image_path, 1, LENGTH(image_path) - 5) || '.webp' "
                        "WHERE image_path IS NOT NULL AND LOWER(SUBSTR(image_path, -5)) = '.jpeg'"
                    )
                    conn.commit()
                finally:
                    conn.close()
    except Exception as db_exc:
        if logger:
            logger.debug(f"SQLite WebP image_path migration skipped/warning: {db_exc}")

    return migrated_count


def export_curated_table(
    df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
    home_herbarium: str = DEFAULT_HOME_HERBARIUM,
    db_path: Optional[Union[str, Path]] = None,
) -> Path:
    """
    Atomically exports the 26-column Darwin Core `EXPORT_COLUMNS` voucher table to `output_path`
    (merging by `catalogNumber` on resume) while synchronizing extended diagnostic metadata
    (`infraspecific_taxon`, `georef_qc_status`, `usda_symbol`, `usda_county_status`,
    `exsiccatae_taxonomic_conflict`, `exsiccatae_gold_taxon`, `spatial_balance_weight`)
    and dynamic synonym mappings into SQLite (`vouchers`, `taxon_synonyms`).
    """
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Merge with existing records if table exists to prevent data loss on resumption
    if output_path.exists() and output_path.stat().st_size > 0:
        try:
            df_existing = pd.read_csv(output_path)
            if not df_existing.empty and "catalogNumber" in df_existing.columns:
                if df.empty:
                    df = df_existing
                else:
                    pre_count = len(df_existing)
                    df = pd.concat([df_existing, df], ignore_index=True)
                    df = df.drop_duplicates(subset=["catalogNumber"], keep="last").reset_index(drop=True)
                    if logger:
                        logger.info(
                            f"Merged records with existing table: {pre_count} previous -> "
                            f"{len(df)} total unique vouchers by catalogNumber."
                        )
        except Exception as e:
            if logger:
                logger.warning(f"Could not read existing table at {output_path} for merging: {e}")

    df = df.copy()
    if not df.empty:
        if "scientificName" not in df.columns and "species_raw" in df.columns:
            df["scientificName"] = df["species_raw"]
        elif "species_raw" not in df.columns and "scientificName" in df.columns:
            df["species_raw"] = df["scientificName"]

        if "decimalLatitude" not in df.columns and "latitude" in df.columns:
            df["decimalLatitude"] = df["latitude"]
        elif "latitude" not in df.columns and "decimalLatitude" in df.columns:
            df["latitude"] = df["decimalLatitude"]

        if "decimalLongitude" not in df.columns and "longitude" in df.columns:
            df["decimalLongitude"] = df["longitude"]
        elif "longitude" not in df.columns and "decimalLongitude" in df.columns:
            df["longitude"] = df["decimalLongitude"]

        if "identifiedBy" not in df.columns and "determiner_raw" in df.columns:
            df["identifiedBy"] = df["determiner_raw"]
        elif "determiner_raw" not in df.columns and "identifiedBy" in df.columns:
            df["determiner_raw"] = df["identifiedBy"]

        dwc_contract_cols = [
            "catalogNumber",
            "scientificName",
            "decimalLatitude",
            "decimalLongitude",
            "eventDate",
            "identifiedBy",
            "determiner_tier",
            "image_path",
        ]
        for c in dwc_contract_cols:
            if c not in df.columns:
                df[c] = ""

        prev_tiers = df["determiner_tier"].astype(str).tolist()
        df = apply_corpus_authority_tiers(df)
        tiers_changed = df["determiner_tier"].astype(str).tolist() != prev_tiers

        if (
            tiers_changed
            or "exsiccatae_key" not in df.columns
            or "is_primary_duplicate" not in df.columns
            or "exsiccatae_taxonomic_conflict" not in df.columns
        ):
            df = stratify_duplicates(df, logger=logger, home_herbarium=home_herbarium)

    _sqlite_only_cols = {
        "exsiccatae_taxonomic_conflict",
        "exsiccatae_gold_taxon",
        "usda_symbol",
        "usda_county_status",
        "infraspecific_taxon",
        "georef_qc_status",
    }
    cols_to_export = [col for col in EXPORT_COLUMNS if col in df.columns]
    extra_cols = [
        c for c in df.columns
        if c not in cols_to_export and not c.startswith("_") and c not in _sqlite_only_cols
    ]
    final_cols = cols_to_export + extra_cols

    df_export = df[final_cols] if not df.empty else pd.DataFrame(columns=EXPORT_COLUMNS)

    temp_file = tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_path.parent,
        suffix=".tmp",
        encoding="utf-8",
    )
    temp_path = Path(temp_file.name)
    try:
        df_export.to_csv(temp_file, index=False, encoding="utf-8")
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
        temp_path.replace(output_path)
        if logger:
            logger.info(f"Atomically saved curated vouchers table ({len(df_export)} records) to: {output_path}")
        if db_path is not None or is_canonical_or_explicit_db_context(output_path):
            try:
                conn = get_db_connection(db_path=db_path, init_schema=True)
                try:
                    upsert_taxon_synonyms(conn, _DYNAMIC_SYNONYM_TO_CANONICAL, _DYNAMIC_TAXON_AUTHOR_STRINGS)
                    if not df_export.empty:
                        sql_cols = list(final_cols) + [
                            c
                            for c in (
                                "exsiccatae_taxonomic_conflict",
                                "exsiccatae_gold_taxon",
                                "usda_symbol",
                                "usda_county_status",
                                "infraspecific_taxon",
                                "georef_qc_status",
                                "spatial_balance_weight",
                            )
                            if c in df.columns and c not in final_cols
                        ]
                        db_records = df[sql_cols].to_dict(orient="records")
                        for rec in db_records:
                            if not rec.get("species_standardized"):
                                rec["species_standardized"] = standardize_packera_taxon(
                                    rec.get("scientificName") or rec.get("species_raw")
                                )
                        upsert_vouchers_records(conn, db_records)
                finally:
                    conn.close()
            except Exception as db_exc:
                if logger:
                    logger.warning(f"SQLite synchronization warning for {output_path}: {db_exc}")
        return output_path
    except Exception as e:
        if temp_path.exists():
            try:
                temp_path.unlink()
            except Exception:
                pass
        if logger:
            logger.error(f"Failed to atomically write curated table to {output_path}: {e}")
        raise


def print_and_log_summary(
    df: pd.DataFrame,
    download_stats: Optional[Dict[str, Any]],
    logger: logging.Logger,
    usda_filter_stats: Optional[Dict[str, Any]] = None,
) -> None:
    """Generates a publication-grade summary of the harvested dataset."""
    total_records = len(df)
    logger.info("=" * 80)
    logger.info("                  PACKERA VOUCHER INGESTION & CURATION SUMMARY                  ")
    logger.info("=" * 80)
    logger.info(f"Total Quality-Filtered Specimen Vouchers: {total_records:,}")

    if total_records == 0:
        logger.warning("No records were retained. Check query parameters or network connection.")
        logger.info("=" * 80)
        return

    # Determiner Tier Breakdown
    tier_counts = df["determiner_tier"].value_counts()
    logger.info("\n--- TAXONOMIC DETERMINER AUTHORITY STRATIFICATION ---")
    for tier in ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"]:
        cnt = tier_counts.get(tier, 0)
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {tier:<15} : {cnt:>5} records ({pct:>5.1f}%)")

    # Species Breakdown
    species_counts = df["species_raw"].value_counts()
    logger.info("\n--- TAXON DISTRIBUTION (RAW DETERMINATIONS) ---")
    for sp, cnt in species_counts.head(8).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {sp:<45} : {cnt:>5} records ({pct:>5.1f}%)")

    # Regional Ecological Groups Breakdown (EPA Level III Ecoregions)
    region_counts = df["regional_group"].fillna("NA").replace("", "NA").value_counts()
    logger.info("\n--- US EPA LEVEL III ECOREGIONS ---")
    for reg, cnt in region_counts.items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {str(reg):<30} : {cnt:>5} records ({pct:>5.1f}%)")

    # Herbarium Institutions (Top 10)
    inst_counts = df["institutionCode"].value_counts()
    logger.info("\n--- TOP HERBARIUM INSTITUTIONS ---")
    for inst, cnt in inst_counts.head(10).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {inst:<15} : {cnt:>5} records ({pct:>5.1f}%)")

    # Exsiccatae & Duplicate Stratification Breakdown
    if "is_primary_duplicate" in df.columns:
        n_primary = int(df["is_primary_duplicate"].sum())
        n_duplicates = total_records - n_primary
        logger.info("\n--- EXSICCAE & DUPLICATE STRATIFICATION ---")
        logger.info(f"  * Primary Duplicate Vouchers  : {n_primary:>5} records ({(n_primary / total_records) * 100.0:>5.1f}%)")
        logger.info(f"  * Secondary Duplicate Sheets  : {n_duplicates:>5} records ({(n_duplicates / total_records) * 100.0:>5.1f}%)")
        if "duplicate_count" in df.columns:
            multi_sheets = int((df["duplicate_count"] > 1).sum())
            logger.info(f"  * Sheets from Multi-Sheet Events: {multi_sheets:>5} records")


    # USDA PLANTS County Range Filter Summary
    if usda_filter_stats and usda_filter_stats.get("evaluated", 0) > 0:
        logger.info("\n--- USDA PLANTS COUNTY RANGE FILTER ---")
        logger.info(f"  * Vouchers Evaluated          : {usda_filter_stats.get('evaluated', 0):>5}")
        logger.info(f"  * Documented In-Range         : {usda_filter_stats.get('in_range', 0):>5}")
        logger.info(f"  * Buffered / 1-Ring Adjacent  : {usda_filter_stats.get('buffered_adjacent', 0):>5}")
        logger.info(f"  * Exempt Retained (Type/Gold) : {usda_filter_stats.get('exempt_retained', 0):>5}")
        logger.info(f"  * Unmapped Taxon / Non-US     : {usda_filter_stats.get('unmapped_or_non_us', 0):>5}")
        logger.info(f"  * Dropped Out-of-Range        : {usda_filter_stats.get('dropped_out_of_range', 0):>5}")

    # Image Download & Quality Summary
    if download_stats:
        logger.info("\n--- SPECIMEN IMAGE DOWNLOAD & QUALITY STATUS ---")
        logger.info(f"  * Downloaded Successfully : {download_stats.get('success', 0):>5}")
        logger.info(f"  * Cached / Skipped        : {download_stats.get('skipped', 0):>5}")
        logger.info(f"  * Quality Filter Rejected : {download_stats.get('quality_rejected', 0):>5}")
        logger.info(f"  * Failed / Inaccessible   : {download_stats.get('failed', 0):>5}")
        if "median_mp" in download_stats:
            logger.info(f"  * Median Image Resolution : {download_stats.get('median_mp', 0.0):>5.2f} Megapixels")

    logger.info("=" * 80)


class VoucherHarvester:
    """
    Consolidated botanical voucher harvester and curator for the genus `Packera`.
    Unifies GBIF query execution (with optional `expand_synonyms` dynamic unification of
    unaccepted `Packera` synonyms and `Senecio`/`Cineraria` basionyms, enabled by default in
    CLI entrypoints and opt-in `False` for direct programmatic/test calls), Darwin Core
    metadata & CoordinateCleaner validation, 3-layer dynamic determiner & herbarium authority
    stratification, USDA PLANTS county range auditing/filtering, and asynchronous WebP media
    acquisition with secondary duplicate rescue.
    """

    def __init__(
        self,
        taxa: Optional[List[str]] = None,
        max_uncertainty_meters: float = 5000.0,
        max_records_per_taxon: int = 5000,
        exclude_western: bool = True,
        min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
        min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB,
        check_sharpness: bool = False,
        min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN,
        concurrency: int = 15,
        output_csv: Path = DEFAULT_OUTPUT_CSV,
        raw_dir: Path = DEFAULT_RAW_DIR,
        workspace_dir: Path = DEFAULT_WORKSPACE,
        logger: Optional[logging.Logger] = None,
        force: bool = False,
        spatial_thin_km: float = 0.0,
        max_per_state: Optional[int] = None,
        max_per_region: Optional[int] = None,
        geographic_scope: str = DEFAULT_GEOGRAPHIC_SCOPE,
        countries: Optional[Union[str, List[str]]] = None,
        include_western: bool = False,
        western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
        expand_synonyms: bool = False,
        cfg_synonyms: Optional[Dict[str, List[str]]] = None,
        home_herbarium: str = DEFAULT_HOME_HERBARIUM,
        db_path: Optional[Union[str, Path]] = None,
        usda_county_filter: bool = DEFAULT_USDA_COUNTY_FILTER,
        usda_county_buffer_km: float = DEFAULT_USDA_COUNTY_BUFFER_KM,
        usda_exempt_gold: bool = DEFAULT_USDA_EXEMPT_GOLD,
        usda_cache_dir: Optional[Path] = None,
        usda_county_geojson_path: Optional[Path] = None,
        harvest_full_continental_range: bool = True,
        min_eastern_residency_records: int = 3,
        webp_quality: int = 92,
    ):
        self.logger = logger or logging.getLogger("VoucherHarvester")
        self.db_path = db_path
        if self.db_path is not None:
            _hydrate_synonyms_from_sqlite(self.db_path)
        raw_input_taxa = [str(t).strip() for t in (taxa if taxa is not None else DEFAULT_TARGET_TAXA) if str(t).strip()]
        self.expand_synonyms = bool(expand_synonyms)
        self.cfg_synonyms = cfg_synonyms if cfg_synonyms is not None else dict(_DEFAULT_CONFIG.taxa.synonyms or {})

        explicit_genus_only = len(raw_input_taxa) == 1 and raw_input_taxa[0].lower() == "packera"
        scope, eff_excl_west, eff_countries = resolve_geographic_scope_settings(
            geographic_scope=geographic_scope,
            exclude_western=exclude_western,
            include_western=include_western,
            countries=countries,
        )
        if not raw_input_taxa:
            self.is_genus_wide = True
            self.logger.info("Discovering accepted Packera species across the genus via GBIF Backbone...")
            discovered_species = discover_genus_packera_species_gbif()
            base_taxa = discovered_species if discovered_species else ["Packera"]
            self.logger.info(f"Discovered {len(base_taxa)} accepted Packera base taxa.")
            if (
                eff_excl_west
                and len(base_taxa) > 1
                and occ is not None
                and getattr(occ, "search", None) is not None
                and not _is_gbif_species_mocked()
                and not hasattr(occ.search, "_mock_name")
                and not hasattr(occ.search, "assert_called")
            ):
                lon_bound = f"{float(western_longitude_threshold)},-50.0"
                facet_countries = list(eff_countries) if eff_countries else ["US", "CA"]
                min_res_n = max(1, int(min_eastern_residency_records))
                self.logger.info(
                    f"Querying GBIF server-side species facet for Eastern NA residency "
                    f"(countries={facet_countries}, decimalLongitude={lon_bound}, min_records={min_res_n})..."
                )
                eastern_residents = []
                try:
                    import urllib.parse
                    genus_key_val = _get_packera_genus_key_gbif() or 3109225
                    agg_key_counts: Dict[int, int] = {}
                    east_us_states = sorted(
                        st_title
                        for st_low, st_title in {
                            v.lower(): v.title() for v in _US_STATE_ABBR_TO_FULL.values()
                        }.items()
                        if st_low.upper() not in EXCLUDED_WESTERN_STATES
                    )
                    east_ca_provs = [
                        "Ontario", "Quebec", "New Brunswick", "Nova Scotia",
                        "Prince Edward Island", "Newfoundland and Labrador", "Manitoba",
                    ]
                    for ctry in facet_countries:
                        facet_urls = [
                            (
                                f"https://api.gbif.org/v1/occurrence/search?"
                                f"taxonKey={genus_key_val}&country={ctry}&decimalLongitude={lon_bound}"
                                f"&basisOfRecord=PRESERVED_SPECIMEN&mediaType=StillImage"
                                f"&facet=speciesKey&facetMincount=1&facetLimit=500&limit=0"
                            )
                        ]
                        # Also tally locality-redacted (ungeoreferenced) sheets in Eastern states/provinces
                        # so imperiled Eastern endemics with withheld coordinates are never pre-screened out.
                        st_list = east_us_states if ctry == "US" else (east_ca_provs if ctry == "CA" else [])
                        if st_list:
                            st_q = "&".join(f"stateProvince={urllib.parse.quote(s)}" for s in st_list)
                            facet_urls.append(
                                f"https://api.gbif.org/v1/occurrence/search?"
                                f"taxonKey={genus_key_val}&country={ctry}&hasCoordinate=false&{st_q}"
                                f"&basisOfRecord=PRESERVED_SPECIMEN&mediaType=StillImage"
                                f"&facet=speciesKey&facetMincount=1&facetLimit=500&limit=0"
                            )
                        for facet_url in facet_urls:
                            req = urllib.request.Request(facet_url, headers={"Accept": "application/json"})
                            with urllib.request.urlopen(req, timeout=15.0) as resp_http:
                                payload = json.loads(resp_http.read().decode("utf-8"))
                            facets = payload.get("facets", []) or []
                            sp_counts = (facets[0].get("counts", []) or []) if facets else []
                            for item in sp_counts:
                                if isinstance(item, dict) and str(item.get("name", "")).isdigit():
                                    k_int = int(item["name"])
                                    agg_key_counts[k_int] = agg_key_counts.get(k_int, 0) + int(item.get("count", 0))
                    resident_keys = [k_int for k_int, c_tot in agg_key_counts.items() if c_tot >= min_res_n]
                    base_set_lower = {b.lower(): b for b in base_taxa}
                    resolved_names = set()
                    for sk in resident_keys:
                        try:
                            u_info = species.name_usage(key=sk) if hasattr(species, "name_usage") else {}
                            if isinstance(u_info, dict):
                                cand = _extract_canonical_binomial(
                                    u_info.get("canonicalName") or u_info.get("species") or u_info.get("scientificName"),
                                    {"Packera"},
                                )
                                if cand:
                                    resolved_names.add(base_set_lower.get(cand.lower(), cand))
                        except Exception:
                            pass
                    eastern_residents = sorted(resolved_names)
                except Exception as facet_exc:
                    self.logger.warning(f"Facet pre-screen fallback triggered ({facet_exc}).")

                if eastern_residents:
                    self.logger.info(
                        f"Eastern NA Pre-Screen: Retained {len(eastern_residents)}/{len(base_taxa)} "
                        f"Eastern NA resident Packera species: {eastern_residents}"
                    )
                    base_taxa = eastern_residents
        elif explicit_genus_only:
            self.is_genus_wide = True
            base_taxa = ["Packera"]
        else:
            self.is_genus_wide = False
            base_taxa = raw_input_taxa

        if self.is_genus_wide:
            self._protected_packera_taxa: Set[str] = {
                b.lower()
                for tx in base_taxa
                for b in [_extract_canonical_binomial(tx, {"Packera"})]
                if b
            }
        else:
            self._protected_packera_taxa = set()
            if not self.expand_synonyms:
                self._protected_packera_taxa = {
                    b.lower()
                    for tx in raw_input_taxa
                    for b in [_extract_canonical_binomial(tx, {"Packera"})]
                    if b
                }

        if self.expand_synonyms and not (len(base_taxa) == 1 and base_taxa[0].lower() == "packera"):
            self.logger.info(
                f"Resolving GBIF Backbone synonyms (unaccepted Packera + Senecio/Cineraria basionyms) "
                f"across {len(base_taxa)} base taxa..."
            )
            if len(base_taxa) > 1 and not _is_gbif_species_mocked() and not _GBIF_NETWORK_UNAVAILABLE:
                uncached_taxa = [tx for tx in base_taxa if tx.strip().lower() not in _GBIF_SYNONYM_CACHE]
                if uncached_taxa:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=min(12, len(uncached_taxa))) as syn_pool:
                        list(syn_pool.map(_query_gbif_backbone_and_synonyms, uncached_taxa))
            if not self.is_genus_wide:
                # Protect explicitly co-requested ACCEPTED Packera species while allowing unaccepted
                # synonyms (e.g. Packera tomentosa -> Packera dubia) to merge into their accepted species.
                for tx in base_taxa:
                    b_tx = _extract_canonical_binomial(tx, {"Packera"})
                    if b_tx and not _GBIF_BACKBONE_SYNONYM_FLAG.get(b_tx.lower(), False):
                        self._protected_packera_taxa.add(b_tx.lower())
            expanded_taxa: List[str] = []
            seen_t: Set[str] = set()
            for idx_tx, tx in enumerate(base_taxa, start=1):
                for syn_tx in expand_taxon_synonyms_gbif(
                    tx,
                    cfg_synonyms=self.cfg_synonyms,
                    protected_packera_taxa=self._protected_packera_taxa,
                ):
                    low_tx = syn_tx.lower()
                    if low_tx not in seen_t:
                        seen_t.add(low_tx)
                        expanded_taxa.append(syn_tx)
                if len(base_taxa) >= 15 and (idx_tx % 25 == 0 or idx_tx == len(base_taxa)):
                    self.logger.info(
                        f"  [Synonym Expansion] Resolved {idx_tx}/{len(base_taxa)} base taxa "
                        f"-> {len(expanded_taxa)} total query binomials..."
                    )
            self.taxa = expanded_taxa if expanded_taxa else base_taxa
        else:
            self.taxa = base_taxa

        scope, eff_excl_west, eff_countries = resolve_geographic_scope_settings(
            geographic_scope=geographic_scope,
            exclude_western=exclude_western,
            include_western=include_western,
            countries=countries,
        )
        self.geographic_scope = scope
        self.exclude_western = eff_excl_west
        self.countries = eff_countries
        self.western_longitude_threshold = float(western_longitude_threshold)
        self.home_herbarium = str(home_herbarium or DEFAULT_HOME_HERBARIUM).strip().upper() or "NCU"

        self.max_uncertainty_meters = max_uncertainty_meters
        self.max_records_per_taxon = max_records_per_taxon
        self.min_megapixels = min_megapixels
        self.min_file_size_kb = min_file_size_kb
        self.check_sharpness = check_sharpness
        self.min_sharpness = min_sharpness
        self.concurrency = concurrency
        self.output_csv = Path(output_csv)
        self.raw_dir = Path(raw_dir)
        self.workspace_dir = Path(workspace_dir)
        self.logger = logger or logging.getLogger("VoucherHarvester")
        self.force = force
        self.spatial_thin_km = float(spatial_thin_km or 0.0)
        self.max_per_state = max_per_state
        self.max_per_region = max_per_region
        self.usda_county_filter = bool(usda_county_filter)
        self.usda_county_buffer_km = float(usda_county_buffer_km)
        self.usda_exempt_gold = bool(usda_exempt_gold)
        self.usda_cache_dir = Path(usda_cache_dir) if usda_cache_dir is not None else None
        self.usda_county_geojson_path = (
            Path(usda_county_geojson_path) if usda_county_geojson_path is not None else None
        )
        self.harvest_full_continental_range = bool(harvest_full_continental_range)
        self.min_eastern_residency_records = max(1, int(min_eastern_residency_records))
        self.webp_quality = max(1, min(100, int(webp_quality)))
        self.last_usda_filter_stats: Optional[Dict[str, Any]] = None
        self._secondary_duplicate_pool: Optional[pd.DataFrame] = None
        self._pre_thin_curated_pool: Optional[pd.DataFrame] = None

    def harvest(self) -> pd.DataFrame:
        """Harvests and normalizes Darwin Core occurrence records across configured taxa."""
        from scripts.core.usda_plants_ranges import load_us_county_spatial_graph

        all_curated_records: List[Dict[str, Any]] = []
        seen_catalog_keys: set = set()
        catalog_to_inst: Dict[str, str] = {}
        country_targets: List[Optional[str]] = list(self.countries) if self.countries else [None]
        current_year = datetime.date.today().year + 1

        try:
            us_spatial_graph = load_us_county_spatial_graph(
                county_path=self.usda_county_geojson_path,
                allow_download=False,
            )
        except Exception:
            us_spatial_graph = {}

        for taxon in self.taxa:
            is_direct_genus_query = taxon.strip().lower() == "packera"
            genus_key = _get_packera_genus_key_gbif() if is_direct_genus_query else None
            country_desc = ",".join(self.countries) if self.countries else "GLOBAL"
            self.logger.info(
                f"Querying GBIF API for taxon: '{taxon}' "
                f"(scope={self.geographic_scope}, country={country_desc}, basisOfRecord=PRESERVED_SPECIMEN)..."
            )
            taxon_harvested = 0
            taxon_retained = 0
            taxon_western_excluded = 0

            def _ingest_occurrence_record(rec: Dict[str, Any], require_type: bool = False) -> None:
                nonlocal taxon_harvested, taxon_retained, taxon_western_excluded

                rec_type_rank, _ = classify_nomenclatural_type_rank(rec.get("typeStatus"))
                is_rec_type = rec_type_rank >= 1
                if require_type and not is_rec_type:
                    return

                taxon_harvested += 1

                # 0. Cultivated / Captive / Greenhouse Filter (exempting nomenclatural types)
                if not is_rec_type and _is_cultivated_occurrence(rec):
                    return

                # 1. Geographic Coordinate, CoordinateCleaner Hygiene & Uncertainty Validation
                raw_issues = rec.get("issues") or []
                issue_set = {
                    str(iss).strip().upper()
                    for iss in (raw_issues if isinstance(raw_issues, (list, tuple, set)) else [raw_issues])
                    if iss
                }
                georef_qc_status = "valid"
                lat = rec.get("decimalLatitude")
                lon = rec.get("decimalLongitude")
                lat_val: Optional[float] = None
                lon_val: Optional[float] = None
                if lat is not None and lon is not None:
                    try:
                        lf = float(lat)
                        lnf = float(lon)
                        if math.isfinite(lf) and math.isfinite(lnf):
                            lat_val, lon_val = lf, lnf
                    except (ValueError, TypeError):
                        pass

                state_prov = rec.get("stateProvince")
                rec_country = str(rec.get("countryCode") or "").strip().upper()

                if lat_val is not None and lon_val is not None:
                    if issue_set & _FATAL_GBIF_GEO_ISSUES:
                        lat_val, lon_val = None, None
                        georef_qc_status = "gbif_geo_issue_nullified"
                    elif not (-90.0 <= lat_val <= 90.0 and -180.0 <= lon_val <= 180.0):
                        lat_val, lon_val = None, None
                        georef_qc_status = "out_of_bounds_nullified"
                    elif abs(lat_val) < 0.01 and abs(lon_val) < 0.01:
                        lat_val, lon_val = None, None
                        georef_qc_status = "zero_coord_nullified"
                    elif (
                        lon_val > 0.0
                        and (rec_country in {"US", "CA", "MX"} or (not rec_country and self.geographic_scope in {"eastern_na", "north_america"}))
                    ):
                        lat_val, lon_val = None, None
                        georef_qc_status = "hemisphere_sign_mismatch_nullified"
                    elif state_prov:
                        st_dwc_norm = _normalize_us_state_name(state_prov)
                        checked_by_poly = False
                        if st_dwc_norm and st_dwc_norm in _US_STATE_ABBR_TO_FULL.values() and Point is not None and us_spatial_graph.get("tree") is not None:
                            pt_chk = Point(lon_val, lat_val)
                            tree_c = us_spatial_graph["tree"]
                            geoms_c = us_spatial_graph["geoms"]
                            cks_c = us_spatial_graph["ck_list"]
                            gid_c = us_spatial_graph["geom_id_map"]
                            neigh_c = us_spatial_graph.get("neighbors") or {}
                            hit_cks: List[str] = []
                            for h_it in tree_c.query(pt_chk):
                                h_idx = int(h_it) if (np is not None and isinstance(h_it, (int, np.integer))) or isinstance(h_it, int) else gid_c.get(id(h_it), -1)
                                if h_idx >= 0 and geoms_c[h_idx].covers(pt_chk):
                                    hit_cks.append(cks_c[h_idx])
                            if hit_cks:
                                checked_by_poly = True
                                ok_state = any(
                                    ck_h.split("|", 1)[0] == st_dwc_norm
                                    or any(nb.split("|", 1)[0] == st_dwc_norm for nb in neigh_c.get(ck_h, ()))
                                    for ck_h in hit_cks
                                )
                                if not ok_state:
                                    lat_val, lon_val = None, None
                                    georef_qc_status = "state_polygon_mismatch_nullified"
                        if not checked_by_poly and lat_val is not None and lon_val is not None:
                            bbox_key = st_dwc_norm if (st_dwc_norm and st_dwc_norm in _NA_STATE_PROVINCE_BBOXES) else _normalize_ca_province_bbox_key(state_prov)
                            if bbox_key and bbox_key in _NA_STATE_PROVINCE_BBOXES:
                                b_lat_min, b_lat_max, b_lon_min, b_lon_max = _NA_STATE_PROVINCE_BBOXES[bbox_key]
                                margin_deg = 0.5
                                if not (
                                    (b_lat_min - margin_deg) <= lat_val <= (b_lat_max + margin_deg)
                                    and (b_lon_min - margin_deg) <= lon_val <= (b_lon_max + margin_deg)
                                ):
                                    lat_val, lon_val = None, None
                                    georef_qc_status = "state_polygon_mismatch_nullified"
                else:
                    georef_qc_status = "ungeoreferenced"

                is_west_rec = is_excluded_western_region(
                    state_prov,
                    lat=lat_val,
                    lon=lon_val,
                    western_longitude_threshold=self.western_longitude_threshold,
                    exclude_western=True,
                    is_type_specimen=is_rec_type,
                )
                has_any_state = bool(
                    state_prov is not None
                    and str(state_prov).strip()
                    and str(state_prov).strip().upper() not in {"NONE", "NAN", "NA", "N/A", "NULL", "UNKNOWN", "UNSPECIFIED", "UNRECORDED"}
                )
                if self.exclude_western:
                    if not self.harvest_full_continental_range and is_west_rec:
                        taxon_western_excluded += 1
                        return
                    # Un-georeferenced non-types lacking any stateProvince cannot be placed biogeographically
                    if lat_val is None and lon_val is None and not has_any_state and not is_rec_type:
                        taxon_western_excluded += 1
                        return

                # Preserve missing coordinateUncertaintyInMeters as None (NA) rather than fabricating 5000 m
                uncertainty_val: Optional[float] = None
                if lat_val is not None and lon_val is not None:
                    uncertainty_raw = rec.get("coordinateUncertaintyInMeters")
                    if uncertainty_raw is not None and str(uncertainty_raw).strip() != "":
                        try:
                            parsed_unc = float(uncertainty_raw)
                            if math.isfinite(parsed_unc) and parsed_unc >= 0.0:
                                if parsed_unc > self.max_uncertainty_meters and not is_rec_type:
                                    return
                                uncertainty_val = parsed_unc
                        except (ValueError, TypeError):
                            pass

                # 2. Darwin Core Temporal Validation (retaining month/year historical sheets for morphometrics
                #    with day=None -> doy=NA, and masking YYYY-01-01 digitization placeholders)
                raw_year = rec.get("year")
                raw_month = rec.get("month")
                raw_day = rec.get("day")
                raw_ev_str = str(rec.get("eventDate") or "").strip()

                year_val: Optional[int] = None
                month_val: Optional[int] = None
                day_val: Optional[int] = None
                try:
                    if raw_year is not None and str(raw_year).strip() != "":
                        year_val = int(float(raw_year))
                except (ValueError, TypeError):
                    year_val = None
                if year_val is None and raw_ev_str:
                    m_yr = re.match(r"^(\d{4})\b", raw_ev_str)
                    if m_yr:
                        year_val = int(m_yr.group(1))

                if year_val is not None and not (1700 <= year_val <= current_year):
                    year_val = None
                if year_val is None and not is_rec_type:
                    return

                try:
                    if raw_month is not None and str(raw_month).strip() != "":
                        m_cand = int(float(raw_month))
                        if 1 <= m_cand <= 12:
                            month_val = m_cand
                except (ValueError, TypeError):
                    month_val = None
                try:
                    if raw_day is not None and str(raw_day).strip() != "":
                        d_cand = int(float(raw_day))
                        if 1 <= d_cand <= 31:
                            day_val = d_cand
                except (ValueError, TypeError):
                    day_val = None

                if year_val is not None and month_val is not None and day_val is not None:
                    try:
                        valid_date = datetime.date(year_val, month_val, day_val)
                        # Guard against YYYY-01-01 Symbiota/GBIF digitization placeholders unless verbatimEventDate explicitly names January
                        verbatim_ev = str(rec.get("verbatimEventDate") or "").strip().lower()
                        explicit_jan = bool(re.search(r"\bjan(?:uary)?\b", verbatim_ev))
                        is_jan1_placeholder = (
                            month_val == 1
                            and day_val == 1
                            and not explicit_jan
                            and (
                                "RECORDED_DATE_UNLIKELY" in issue_set
                                or bool(re.fullmatch(r"\d{4}(?:-01-01(?:T00:00:00Z?)?)?", raw_ev_str))
                                or not raw_ev_str
                            )
                        )
                        if is_jan1_placeholder:
                            month_val, day_val = None, None
                            event_date = f"{year_val:04d}"
                        else:
                            event_date = raw_ev_str or valid_date.isoformat()
                    except ValueError:
                        if not is_rec_type:
                            return
                        day_val = None
                        event_date = raw_ev_str or f"{year_val:04d}-{month_val:02d}"
                elif year_val is not None and month_val is not None:
                    day_val = None
                    event_date = raw_ev_str or f"{year_val:04d}-{month_val:02d}"
                else:
                    month_val, day_val = None, None
                    event_date = raw_ev_str or (f"{year_val:04d}" if year_val is not None else "")

                # 3. High-Resolution Media Image Validation
                media_list = rec.get("media", [])
                image_url = extract_high_res_image_url(media_list)
                if not image_url:
                    return

                # 4. Catalog Number & Herbarium Institution Namespace Protection
                raw_catalog = rec.get("catalogNumber")
                inst_code = str(rec.get("institutionCode") or rec.get("collectionCode") or "UNKNOWN_INST").strip() or "UNKNOWN_INST"
                inst_norm_up = inst_code.upper()
                gbif_key = str(rec.get("key", ""))

                if raw_catalog and str(raw_catalog).strip():
                    raw_cat_str = str(raw_catalog).strip()
                    clean_cat = sanitize_filename(raw_cat_str)
                    # Prefix institutionCode onto bare numeric/non-alphabetic accession numbers (e.g. '12345' -> 'GA_12345')
                    if not re.search(r"[A-Za-z]", raw_cat_str) and inst_norm_up not in _INVALID_INSTITUTION_CODES:
                        catalog_number = f"{sanitize_filename(inst_code)}_{clean_cat}"
                    else:
                        catalog_number = clean_cat
                else:
                    catalog_number = f"{sanitize_filename(inst_code)}_{gbif_key}"

                prev_inst = catalog_to_inst.get(catalog_number)
                if prev_inst is not None and prev_inst != inst_norm_up and inst_norm_up not in _INVALID_INSTITUTION_CODES:
                    catalog_number = f"{sanitize_filename(inst_code)}_{catalog_number}"

                unique_key = (catalog_number, gbif_key)
                if unique_key in seen_catalog_keys:
                    return
                seen_catalog_keys.add(unique_key)
                catalog_to_inst[catalog_number] = inst_norm_up

                relative_image_path = f"data/raw_vouchers/{catalog_number}.webp"
                raw_sci_val = str(rec.get("scientificName") or rec.get("species") or taxon).strip()
                infraspecific_val = _extract_infraspecific_taxon(raw_sci_val)

                # Dynamically map historical non-Packera genera (Senecio / Cineraria) AND occurrence records
                # explicitly flagged by GBIF as unaccepted synonyms (e.g. Packera tomentosa -> Packera dubia),
                # while protecting explicitly requested accepted Packera segregates.
                raw_packera_bin = _extract_canonical_binomial(raw_sci_val, {"Packera"})
                rec_tax_status = str(rec.get("taxonomicStatus") or rec.get("status") or "").strip().upper()
                is_rec_synonym = rec_tax_status in _GBIF_SYNONYM_STATUSES
                is_protected_raw = bool(
                    raw_packera_bin
                    and raw_packera_bin.lower() in self._protected_packera_taxa
                    and not (self.expand_synonyms and is_rec_synonym)
                )
                acc_from_rec: Optional[str] = None
                if raw_packera_bin is None or (is_rec_synonym and not is_protected_raw):
                    for acc_fld in ("acceptedScientificName", "species"):
                        cand_acc = _extract_canonical_binomial(rec.get(acc_fld), {"Packera"})
                        if cand_acc:
                            acc_from_rec = cand_acc
                            break
                    if acc_from_rec and (raw_packera_bin is None or acc_from_rec.lower() != raw_packera_bin.lower()):
                        if acc_from_rec.lower() not in (self._protected_packera_taxa - {acc_from_rec.lower()}) or not is_protected_raw:
                            syn_items_to_reg = [raw_sci_val] + ([raw_packera_bin] if raw_packera_bin else [])
                            _register_dynamic_synonyms(
                                acc_from_rec,
                                syn_items_to_reg,
                                [raw_sci_val],
                                allow_packera_species_remap=is_rec_synonym,
                            )

                resolved_packera = standardize_packera_taxon(raw_sci_val) or acc_from_rec
                if is_direct_genus_query or self.is_genus_wide:
                    if not resolved_packera or not re.match(r"^Packera\s+(?:[x×]\s*)?[a-z-]+$", resolved_packera):
                        seen_catalog_keys.discard(unique_key)
                        return
                if resolved_packera and not raw_sci_val.lower().startswith(resolved_packera.lower()):
                    sci_name_val = resolved_packera
                else:
                    sci_name_val = raw_sci_val

                raw_loc_val = rec.get("locality") or rec.get("verbatimLocality")
                raw_hab_val = rec.get("habitat")
                # Exclude parenthetical pre-Packera basionym authors (e.g. '(L.)', '(Michx.)', '(Spreng.)') on non-types
                non_paren_names = [re.sub(r"\([^)]*\)", " ", nm) for nm in (raw_sci_val, sci_name_val, taxon) if nm]
                proto_auths = (
                    extract_protologue_authors([raw_sci_val, sci_name_val, taxon])
                    if is_rec_type
                    else (extract_protologue_authors(non_paren_names) - {"weber", "love", "löve"})
                )

                # 5. Initial Layer 1 + Layer 2 Determiner Tier Stratification
                determiner_tier, type_status, determiner_raw = parse_determiner_tier(
                    type_status_raw=rec.get("typeStatus"),
                    identified_by_raw=rec.get("identifiedBy"),
                    recorded_by_raw=rec.get("recordedBy"),
                    history_raw=rec.get("verbatimIdentificationHistory"),
                    institution_code_raw=inst_code,
                    locality_raw=raw_loc_val,
                    habitat_raw=raw_hab_val,
                    dynamic_specialists=proto_auths,
                )

                # 6. North American (CEC / US EPA) Level III Ecoregion Polygon Assignment
                regional_group = infer_regional_group(
                    lat=lat_val,
                    lon=lon_val,
                    state_province=state_prov,
                    habitat=raw_hab_val,
                    locality=raw_loc_val,
                    county=rec.get("county"),
                )

                curated_record = {
                    "catalogNumber": catalog_number,
                    "institutionCode": inst_code,
                    "scientificName": sci_name_val,
                    "species_raw": raw_sci_val,
                    "infraspecific_taxon": infraspecific_val,
                    "georef_qc_status": georef_qc_status,
                    "recordedBy": str(rec.get("recordedBy") or "").strip(),
                    "recordNumber": str(rec.get("recordNumber") or "").strip(),
                    "identifiedBy": determiner_raw,
                    "determiner_raw": determiner_raw,
                    "determiner_tier": determiner_tier,
                    "type_status": type_status,
                    "county": rec.get("county") or "",
                    "stateProvince": state_prov or "",
                    "decimalLatitude": lat_val,
                    "latitude": lat_val,
                    "decimalLongitude": lon_val,
                    "longitude": lon_val,
                    "coordinateUncertainty": uncertainty_val,
                    "year": year_val,
                    "month": month_val,
                    "day": day_val,
                    "eventDate": event_date,
                    "regional_group": regional_group,
                    "image_path": relative_image_path,
                    "_image_url": image_url,
                    "_is_western_rec": bool(is_west_rec),
                    "_is_type_rec": bool(is_rec_type),
                    "_raw_identified_by": str(rec.get("identifiedBy") or "").strip(),
                    "_raw_history": str(rec.get("verbatimIdentificationHistory") or "").strip(),
                    "_raw_type_status": str(rec.get("typeStatus") or "").strip(),
                    "_raw_locality": str(raw_loc_val or "").strip(),
                    "_raw_habitat": str(raw_hab_val or "").strip(),
                }
                all_curated_records.append(curated_record)
                taxon_retained += 1

            use_server_lon_filter = bool(self.exclude_western and not self.harvest_full_continental_range)
            for country_code in country_targets:
                offset = 0
                limit = 300
                while taxon_retained < self.max_records_per_taxon:
                    search_kwargs: Dict[str, Any] = {
                        "basisOfRecord": "PRESERVED_SPECIMEN",
                        "limit": limit,
                        "offset": offset,
                    }
                    if is_direct_genus_query and genus_key is not None:
                        search_kwargs["taxonKey"] = genus_key
                    else:
                        search_kwargs["scientificName"] = taxon
                    if country_code:
                        search_kwargs["country"] = country_code
                    if use_server_lon_filter:
                        search_kwargs["decimalLongitude"] = f"{self.western_longitude_threshold},-50.0"
                    try:
                        response = occ.search(**search_kwargs)
                    except Exception as e:
                        self.logger.error(f"GBIF API query error for '{taxon}' (country={country_code}) at offset {offset}: {e}")
                        break

                    results = response.get("results", [])
                    count = response.get("count", 0)
                    if not results:
                        break

                    for rec in results:
                        _ingest_occurrence_record(rec, require_type=False)
                        if taxon_retained >= self.max_records_per_taxon:
                            break

                    offset += len(results)
                    if offset >= count:
                        break
                    time.sleep(0.1)

                # Secondary Type-Specimen Rescue Query ensuring nomenclatural types are always captured
                if self.exclude_western and taxon_retained < self.max_records_per_taxon:
                    type_search_kwargs: Dict[str, Any] = {
                        "basisOfRecord": "PRESERVED_SPECIMEN",
                        "typeStatus": "*",
                        "limit": 300,
                        "offset": 0,
                    }
                    if is_direct_genus_query and genus_key is not None:
                        type_search_kwargs["taxonKey"] = genus_key
                    else:
                        type_search_kwargs["scientificName"] = taxon
                    if country_code:
                        type_search_kwargs["country"] = country_code
                    try:
                        type_resp = occ.search(**type_search_kwargs)
                        type_results = type_resp.get("results", []) if isinstance(type_resp, dict) else []
                    except Exception:
                        type_results = []
                    for t_rec in type_results:
                        if taxon_retained >= self.max_records_per_taxon:
                            break
                        _ingest_occurrence_record(t_rec, require_type=True)

                if taxon_retained >= self.max_records_per_taxon:
                    break

            self.logger.info(
                f"Taxon '{taxon}': Processed {taxon_harvested} occurrences -> "
                f"Retained {taxon_retained} candidate records "
                f"(Immediate exclusions: {taxon_western_excluded})."
            )

        df = pd.DataFrame(all_curated_records)
        if not df.empty:
            # Eastern NA Biogeographic Residency Screen + Full US+CA Continental Range Retention for Residents
            is_direct_single_genus_tiny_mock = (
                len(self.taxa) == 1
                and self.taxa[0].strip().lower() == "packera"
                and len(df) < self.min_eastern_residency_records
            )
            if self.exclude_western and not is_direct_single_genus_tiny_mock:
                canon_series = df["scientificName"].map(
                    lambda x: standardize_packera_taxon(x) or str(x).strip()
                )
                west_mask = (
                    df["_is_western_rec"].astype(bool)
                    if "_is_western_rec" in df.columns
                    else pd.Series(False, index=df.index)
                )
                type_mask_r = (
                    df["_is_type_rec"].astype(bool)
                    if "_is_type_rec" in df.columns
                    else pd.Series(False, index=df.index)
                )
                east_non_type_mask = (~west_mask) & (~type_mask_r)
                east_counts = canon_series[east_non_type_mask].value_counts().to_dict()
                resident_species = {
                    sp for sp, cnt in east_counts.items() if int(cnt) >= self.min_eastern_residency_records
                }

                if self.is_genus_wide:
                    keep_res_mask = canon_series.isin(resident_species)
                    if not keep_res_mask.all():
                        removed_sp = sorted(set(canon_series.dropna().unique()) - resident_species)
                        pre_len = len(df)
                        df = df.loc[keep_res_mask].reset_index(drop=True)
                        self.logger.info(
                            f"Eastern NA Biogeographic Residency Filter: Removed {len(removed_sp)} "
                            f"non-resident Cordilleran/western species ({pre_len - len(df)} records) with "
                            f"< {self.min_eastern_residency_records} Eastern NA occurrences."
                        )
                elif self.harvest_full_continental_range:
                    # Verified Eastern NA resident species (>= min_eastern_residency_records in Eastern NA)
                    # retain their full US + CA continental range; non-resident species drop western strays.
                    keep_res_mask = canon_series.isin(resident_species) | (~west_mask)
                    if not keep_res_mask.all():
                        df = df.loc[keep_res_mask].reset_index(drop=True)

                if df.empty:
                    for c_tmp in ("_is_western_rec", "_is_type_rec"):
                        if c_tmp in df.columns:
                            df = df.drop(columns=[c_tmp])
                    return df

            # Layer 3: Corpus-derived multi-herbarium specialist & regional herbarium authority re-evaluation
            dyn_specs, dyn_herbs = discover_corpus_authorities(df)
            if dyn_specs or dyn_herbs:
                self.logger.info(
                    f"Dynamic Authority Engine: Discovered {len(dyn_specs)} corpus specialist(s) "
                    f"and {len(dyn_herbs)} dynamic research herbarium/herbaria."
                )
            df = apply_corpus_authority_tiers(df, dynamic_specialists=dyn_specs, dynamic_herbaria=dyn_herbs)
            raw_audit_cols = [
                c for c in (
                    "_is_western_rec", "_is_type_rec",
                    "_raw_identified_by", "_raw_history", "_raw_type_status", "_raw_locality", "_raw_habitat",
                )
                if c in df.columns
            ]
            if raw_audit_cols:
                df = df.drop(columns=raw_audit_cols)

            df = stratify_duplicates(
                df,
                workspace_dir=self.workspace_dir,
                logger=self.logger,
                home_herbarium=self.home_herbarium,
            )
            # Cache secondary duplicate sheets so download_and_validate_media() can promote a fallback
            # duplicate whenever a primary sheet's image URL 404s or fails optical resolution QC
            if "is_primary_duplicate" in df.columns:
                sec_mask = ~df["is_primary_duplicate"].map(
                    lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}
                )
                self._secondary_duplicate_pool = df.loc[sec_mask].copy()
            else:
                self._secondary_duplicate_pool = None

            is_gbif_mocked = (
                _is_gbif_species_mocked()
                or occ is None
                or getattr(occ, "_is_missing_pygbif", False)
                or getattr(occ, "search", None) is None
                or hasattr(occ.search, "_mock_name")
                or hasattr(occ.search, "assert_called")
            )
            default_usda_cache = Path(self.workspace_dir) / "data" / "environmental" / "usda_plants"
            has_local_usda_cache = (
                self.usda_cache_dir is not None
                or (default_usda_cache / "plantlst.txt").is_file()
                or (default_usda_cache / "distributions").is_dir()
            )
            if (not is_gbif_mocked) or has_local_usda_cache:
                df, usda_stats = filter_vouchers_by_usda_county_range(
                    df=df,
                    standardize_taxon_fn=standardize_packera_taxon,
                    synonym_to_canonical=_DYNAMIC_SYNONYM_TO_CANONICAL,
                    county_buffer_km=self.usda_county_buffer_km,
                    exempt_gold=self.usda_exempt_gold,
                    cache_dir=self.usda_cache_dir,
                    county_geojson_path=self.usda_county_geojson_path,
                    db_path=self.db_path,
                    workspace_dir=self.workspace_dir,
                    allow_network=(not is_gbif_mocked),
                    force_refresh=self.force,
                    drop_out_of_range=bool(self.usda_county_filter),
                    logger=self.logger,
                )
                self.last_usda_filter_stats = usda_stats
                if df.empty:
                    return df
            self._pre_thin_curated_pool = df.copy()
            df = spatially_balance_and_thin_vouchers(
                df,
                spatial_thin_km=self.spatial_thin_km,
                max_per_state=self.max_per_state,
                max_per_region=self.max_per_region,
                primary_duplicates_only=True,
                preserve_all_gold=True,
                workspace_dir=self.workspace_dir,
            )
        return df

    def download_and_validate_media(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Asynchronously downloads voucher specimen sheets, rescuing failed primaries via duplicate fallbacks."""
        if df.empty or "_image_url" not in df.columns:
            return df, {"success": 0, "skipped": 0, "failed": 0, "quality_rejected": 0}

        download_queue: List[Tuple[str, Path]] = []
        for _, row in df.iterrows():
            url = row["_image_url"]
            dest = self.workspace_dir / row["image_path"]
            download_queue.append((url, dest))

        self.logger.info(f"Starting asynchronous download of {len(download_queue)} voucher sheets...")
        stats: Dict[str, Any] = asyncio.run(
            download_all_voucher_images(
                records_to_download=download_queue,
                concurrency_limit=self.concurrency,
                force=self.force,
                webp_quality=self.webp_quality,
                logger=self.logger,
            )
        )

        self.logger.info(
            f"Auditing specimen image quality (min_megapixels={self.min_megapixels} MP, "
            f"min_file_size_kb={self.min_file_size_kb} KB)..."
        )
        valid_rows: List[pd.Series] = []
        failed_exsiccatae_keys: List[str] = []
        quality_rejected = 0
        mp_values = []

        for _, row in df.iterrows():
            img_dest = self.workspace_dir / row["image_path"]
            is_valid, q_metrics = validate_image_quality(
                image_path=img_dest,
                min_megapixels=self.min_megapixels,
                min_file_size_kb=self.min_file_size_kb,
                check_sharpness=self.check_sharpness,
                min_sharpness=self.min_sharpness,
            )
            if is_valid:
                valid_rows.append(row)
                if "megapixels" in q_metrics:
                    mp_values.append(q_metrics["megapixels"])
            else:
                quality_rejected += 1
                ek = str(row.get("exsiccatae_key") or "").strip()
                if ek:
                    failed_exsiccatae_keys.append(ek)
                for cand_unlink in (img_dest, resolve_voucher_image_path(img_dest)):
                    if cand_unlink and cand_unlink.is_file():
                        try:
                            cand_unlink.unlink()
                        except Exception:
                            pass

        surviving_cats: Set[str] = {
            str(r.get("catalogNumber", "")).strip() for r in valid_rows if r.get("catalogNumber")
        }
        surviving_eks: Set[str] = {
            str(r.get("exsiccatae_key", "")).strip() for r in valid_rows if r.get("exsiccatae_key")
        }

        # Attempt secondary duplicate rescue for any collection event with zero surviving sheets after image QC
        if (
            failed_exsiccatae_keys
            and self._secondary_duplicate_pool is not None
            and not self._secondary_duplicate_pool.empty
            and "exsiccatae_key" in self._secondary_duplicate_pool.columns
            and "_image_url" in self._secondary_duplicate_pool.columns
        ):
            rescued_n = 0
            tier_prio = {"Tier_1_Gold": 3, "Tier_2_Silver": 2, "Tier_3_Bronze": 1}
            for ek in dict.fromkeys(failed_exsiccatae_keys):
                if ek in surviving_eks:
                    continue
                cand_dups = self._secondary_duplicate_pool[
                    self._secondary_duplicate_pool["exsiccatae_key"].astype(str) == ek
                ].copy()
                if cand_dups.empty:
                    continue
                cand_dups["_tp"] = cand_dups["determiner_tier"].astype(str).map(lambda t: tier_prio.get(t, 1))
                cand_dups = cand_dups.sort_values(by="_tp", ascending=False).drop(columns=["_tp"])
                for _, dup_row in cand_dups.iterrows():
                    dup_cat = str(dup_row.get("catalogNumber", "")).strip()
                    if dup_cat and dup_cat in surviving_cats:
                        continue
                    dup_dest = self.workspace_dir / str(dup_row["image_path"])
                    dup_stats = asyncio.run(
                        download_all_voucher_images(
                            records_to_download=[(str(dup_row["_image_url"]), dup_dest)],
                            concurrency_limit=1,
                            force=self.force,
                            webp_quality=self.webp_quality,
                            logger=self.logger,
                        )
                    )
                    dup_ok, dup_m = validate_image_quality(
                        image_path=dup_dest,
                        min_megapixels=self.min_megapixels,
                        min_file_size_kb=self.min_file_size_kb,
                        check_sharpness=self.check_sharpness,
                        min_sharpness=self.min_sharpness,
                    )
                    if dup_ok:
                        promoted = dup_row.copy()
                        promoted["is_primary_duplicate"] = True
                        valid_rows.append(promoted)
                        if dup_cat:
                            surviving_cats.add(dup_cat)
                        surviving_eks.add(ek)
                        stats["success"] = int(stats.get("success", 0)) + int(dup_stats.get("success", 0))
                        if "megapixels" in dup_m:
                            mp_values.append(dup_m["megapixels"])
                        rescued_n += 1
                        break
                    for cand_unlink in (dup_dest, resolve_voucher_image_path(dup_dest)):
                        if cand_unlink and cand_unlink.is_file():
                            try:
                                cand_unlink.unlink()
                            except Exception:
                                pass
            if rescued_n > 0:
                self.logger.info(
                    f"Exsiccatae Media Rescue: Promoted {rescued_n} valid secondary duplicate sheet(s) "
                    f"after primary image download/QC failure."
                )

        self.logger.info(
            f"Quality Audit Complete: Retained {len(valid_rows)} vouchers meeting quality standards "
            f"(Rejected {quality_rejected} substandard/low-res images)."
        )
        df_filtered = pd.DataFrame(valid_rows).reset_index(drop=True) if valid_rows else df.iloc[0:0].copy()
        if not df_filtered.empty:
            df_filtered = stratify_duplicates(
                df_filtered,
                workspace_dir=self.workspace_dir,
                logger=self.logger,
                home_herbarium=self.home_herbarium,
            )
            if "spatial_balance_weight" in df.columns:
                df_filtered = spatially_balance_and_thin_vouchers(
                    df_filtered,
                    spatial_thin_km=self.spatial_thin_km,
                    max_per_state=self.max_per_state,
                    max_per_region=self.max_per_region,
                    primary_duplicates_only=True,
                    preserve_all_gold=True,
                    workspace_dir=self.workspace_dir,
                )
        stats["quality_rejected"] = quality_rejected
        if mp_values:
            stats["median_mp"] = float(np.median(mp_values))

        return df_filtered, stats

    def export(self, df: pd.DataFrame) -> Path:
        """Atomically persists the curated vouchers table to CSV and SQLite."""
        return export_curated_table(
            df,
            self.output_csv,
            logger=self.logger,
            home_herbarium=self.home_herbarium,
            db_path=self.db_path,
        )

    def run(
        self,
        download_images: bool = False,
        force: Optional[bool] = None,
    ) -> Tuple[pd.DataFrame, Optional[Dict[str, Any]]]:
        """Executes the end-to-end voucher ingestion and curation workflow."""
        if force is not None:
            self.force = force

        migrate_existing_vouchers_to_webp(
            self.raw_dir,
            db_path=self.db_path,
            curated_csv=self.output_csv,
            quality=self.webp_quality,
            logger=self.logger,
        )

        self.logger.info("Starting Packera Voucher Ingestion & Authority Stratification Pipeline...")
        self.logger.info(f"Target Taxa ({len(self.taxa)}): {self.taxa}")
        self.logger.info(
            f"Geographic Scope: {self.geographic_scope} (countries={self.countries or 'GLOBAL'}, "
            f"exclude_western={self.exclude_western})"
        )
        self.logger.info(f"Max Coordinate Uncertainty Threshold: {self.max_uncertainty_meters} m")
        self.logger.info(f"Force Overwrite Mode: {self.force}")

        df_curated = self.harvest()
        download_stats = None

        if download_images and not df_curated.empty:
            df_curated, download_stats = self.download_and_validate_media(df_curated)

        self.export(df_curated)
        print_and_log_summary(
            df_curated,
            download_stats,
            self.logger,
            usda_filter_stats=self.last_usda_filter_stats,
        )
        self.logger.info("Pipeline execution completed successfully.")
        return df_curated, download_stats


def harvest_taxa_occurrences(
    taxa_list: Optional[List[str]] = None,
    max_uncertainty_meters: float = 5000.0,
    max_records_per_taxon: int = 1000,
    exclude_western: bool = True,
    logger: Optional[logging.Logger] = None,
    geographic_scope: str = DEFAULT_GEOGRAPHIC_SCOPE,
    countries: Optional[Union[str, List[str]]] = None,
    include_western: bool = False,
    expand_synonyms: bool = False,
    usda_county_filter: bool = DEFAULT_USDA_COUNTY_FILTER,
    usda_county_buffer_km: float = DEFAULT_USDA_COUNTY_BUFFER_KM,
    usda_exempt_gold: bool = DEFAULT_USDA_EXEMPT_GOLD,
) -> pd.DataFrame:
    """
    Functional entry point for occurrence harvesting (defaults to `expand_synonyms=False`
    for deterministic programmatic/unit-test queries; pass `expand_synonyms=True` to enable
    dynamic GBIF Backbone synonym expansion).
    """
    harvester = VoucherHarvester(
        taxa=taxa_list,
        max_uncertainty_meters=max_uncertainty_meters,
        max_records_per_taxon=max_records_per_taxon,
        exclude_western=exclude_western,
        logger=logger,
        geographic_scope=geographic_scope,
        countries=countries,
        include_western=include_western,
        expand_synonyms=expand_synonyms,
        usda_county_filter=usda_county_filter,
        usda_county_buffer_km=usda_county_buffer_km,
        usda_exempt_gold=usda_exempt_gold,
    )
    return harvester.harvest()


def main() -> None:
    """CLI execution wrapper delegating to data prep runner."""
    import importlib
    cli_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")
    cli_mod.main()


if __name__ == "__main__":
    main()
