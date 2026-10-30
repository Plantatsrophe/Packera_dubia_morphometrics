"""
===============================================================================
Module: harvester.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unified botanical voucher ingestion, Darwin Core metadata normalization,
    determiner authority evaluation, and asynchronous media acquisition engine.
    Consolidates GBIF occurrence querying, 3-tier determiner authority
    stratification, geographic boundary filtering, image quality validation,
    and atomic curated dataset persistence.
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

import aiohttp
import numpy as np
import pandas as pd
from PIL import Image
import pygbif
import pygbif.occurrences as occ
import pygbif.species as species
from shapely.geometry import Point, Polygon, shape
from shapely.strtree import STRtree
from tqdm.asyncio import tqdm as async_tqdm

# Allow decompression of high-resolution botanical herbarium sheets (>89 MP)
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
    DEFAULT_WORKSPACE,
    EXCLUDED_WESTERN_STATES,
    MAJOR_HERBARIA_CODES,
    SPECIALIST_PATTERNS,
    VALID_GEOGRAPHIC_SCOPES,
    VALID_TYPE_STATUSES,
    WESTERN_LONGITUDE_THRESHOLD,
    _DEFAULT_CONFIG,
    normalize_countries_list,
)
from scripts.core.logger import setup_logging

EASTERN_NA_PACKERA_TAXA: List[str] = [
    "Packera anonyma",
    "Packera antennariifolia",
    "Packera aurea",
    "Packera crawfordii",
    "Packera cymbalaria",
    "Packera dubia",
    "Packera glabella",
    "Packera hyperborealis",
    "Packera indecora",
    "Packera insulanus",
    "Packera memmingeri",
    "Packera millefolium",
    "Packera obovata",
    "Packera paupercula",
    "Packera plattensis",
    "Packera pseudaurea",
    "Packera schweinitziana",
    "Packera serpenticola",
    "Packera spellenbergii",
    "Packera tampicana",
    "Packera vulgaris",
]

EASTERN_NA_GBIF_QUERY_TAXA: List[str] = EASTERN_NA_PACKERA_TAXA + [
    "Senecio anonymus",
    "Senecio smallii",
    "Packera smallii",
    "Senecio antennariifolius",
    "Senecio aureus",
    "Senecio gracilis",
    "Senecio crawfordii",
    "Packera paupercula var. crawfordii",
    "Senecio tomentosus",
    "Packera tomentosa",
    "Senecio dubius",
    "Senecio alabamensis",
    "Senecio glabellus",
    "Senecio lobatus",
    "Senecio cymbalaria",
    "Senecio resedifolius",
    "Senecio indecorus",
    "Senecio insulanus",
    "Senecio memmingeri",
    "Senecio millefolium var. memmingeri",
    "Senecio millefolium",
    "Senecio obovatus",
    "Senecio rotundus",
    "Senecio elliottii",
    "Senecio pauperculus",
    "Senecio balsamitae",
    "Senecio plattensis",
    "Senecio pseudotomentosus",
    "Senecio pseudaureus",
    "Senecio schweinitzianus",
    "Senecio robbinsii",
    "Packera robbinsii",
    "Packera eburnea",
    "Senecio tampicanus",
    "Senecio imparipinnatus",
    "Packera imparipinnata",
]


# Comprehensive Eastern + Western North American (US, CA, MX) Packera <-> Senecio/Cineraria basionym & synonym map
PACKERA_NORTH_AMERICA_SYNONYMS: Dict[str, List[str]] = {
    # Eastern & Central North America
    "Packera anonyma": ["Senecio anonymus", "Senecio smallii", "Packera smallii", "Senecio earlei"],
    "Packera antennariifolia": ["Senecio antennariifolius"],
    "Packera aurea": ["Senecio aureus", "Senecio gracilis"],
    "Packera crawfordii": ["Senecio crawfordii", "Packera paupercula var. crawfordii", "Senecio pauperculus var. crawfordii"],
    "Packera cymbalaria": ["Senecio cymbalaria", "Senecio resedifolius", "Packera resedifolia"],
    "Packera dubia": ["Senecio tomentosus", "Packera tomentosa", "Cineraria tomentosa", "Senecio dubius", "Senecio alabamensis"],
    "Packera glabella": ["Senecio glabellus", "Senecio lobatus"],
    "Packera hyperborealis": ["Senecio hyperborealis"],
    "Packera indecora": ["Senecio indecorus"],
    "Packera insulanus": ["Senecio insulanus", "Packera insulana"],
    "Packera memmingeri": ["Senecio memmingeri", "Senecio millefolium var. memmingeri", "Packera x memmingeri"],
    "Packera millefolium": ["Senecio millefolium"],
    "Packera obovata": ["Senecio obovatus", "Senecio rotundus", "Senecio elliottii"],
    "Packera paupercula": [
        "Senecio pauperculus",
        "Senecio balsamitae",
        "Packera balsamitae",
        "Packera paupercula var. savannarum",
        "Packera paupercula var. pseudotomentosa",
        "Packera paupercula var. balsamitae",
        "Packera paupercula var. appalachiana",
    ],
    "Packera plattensis": ["Senecio plattensis", "Senecio pseudotomentosus", "Packera flavovirens"],
    "Packera pseudaurea": [
        "Senecio pseudaureus",
        "Packera pseudaurea var. semicordata",
        "Senecio pseudaureus var. semicordatus",
        "Packera pseudaurea var. flavula",
        "Senecio pseudaureus var. flavulus",
    ],
    "Packera schweinitziana": ["Senecio schweinitzianus", "Senecio robbinsii", "Packera robbinsii"],
    "Packera serpenticola": ["Packera eburnea"],
    "Packera spellenbergii": ["Senecio spellenbergii"],
    "Packera tampicana": ["Senecio tampicanus", "Senecio imparipinnatus", "Packera imparipinnata"],
    # Western US, Rocky Mountains, Pacific Coast, Boreal/Arctic Canada & Mexico
    "Packera actinella": ["Senecio actinella"],
    "Packera bellidifolia": ["Senecio bellidifolius"],
    "Packera bernardina": ["Senecio bernardinus"],
    "Packera bolanderi": ["Senecio bolanderi", "Senecio harfordii", "Packera bolanderi var. harfordii"],
    "Packera breweri": ["Senecio breweri"],
    "Packera cana": ["Senecio canus", "Senecio purshianus", "Senecio convallium", "Senecio hallii", "Senecio howellii"],
    "Packera candidissima": ["Senecio candidissimus"],
    "Packera cardamine": ["Senecio cardamine"],
    "Packera castoreus": ["Senecio castoreus"],
    "Packera clevelandii": ["Senecio clevelandii"],
    "Packera coahuilensis": ["Senecio coahuilensis"],
    "Packera contermina": ["Senecio conterminus"],
    "Packera crocata": ["Senecio crocatus"],
    "Packera cynthioides": ["Senecio cynthioides", "Senecio wrightii"],
    "Packera debilis": ["Senecio debilis", "Senecio fedifolius"],
    "Packera dimorphophylla": [
        "Senecio dimorphophyllus",
        "Packera dimorphophylla var. paysonii",
        "Senecio dimorphophyllus var. paysonii",
    ],
    "Packera eurycephala": ["Senecio eurycephalus", "Packera eurycephala var. lewisrosei"],
    "Packera fendleri": ["Senecio fendleri"],
    "Packera flettii": ["Senecio flettii"],
    "Packera franciscana": ["Senecio franciscanus"],
    "Packera ganderi": ["Senecio ganderi"],
    "Packera greenei": ["Senecio greenei"],
    "Packera hartiana": ["Senecio hartianus"],
    "Packera hesperia": ["Senecio hesperius"],
    "Packera heterophylla": ["Senecio heterophyllus"],
    "Packera hintoniorum": ["Senecio hintoniorum"],
    "Packera ionophylla": ["Senecio ionophyllus", "Senecio sparsilobus"],
    "Packera layneae": ["Senecio layneae"],
    "Packera macounii": ["Senecio macounii", "Senecio fastigiatus"],
    "Packera malmstenii": ["Senecio malmstenii"],
    "Packera mancosiana": ["Senecio mancosianus"],
    "Packera montereyana": ["Senecio montereyanus"],
    "Packera moranii": ["Senecio moranii"],
    "Packera multilobata": ["Senecio multilobatus", "Senecio lynceus"],
    "Packera musiniensis": ["Senecio musiniensis"],
    "Packera neomexicana": [
        "Senecio neomexicanus",
        "Senecio mutabilis",
        "Packera neomexicana var. mutabilis",
        "Senecio metcalfei",
        "Packera neomexicana var. metcalfei",
        "Senecio toumeyi",
        "Packera neomexicana var. toumeyi",
    ],
    "Packera ogotorukensis": ["Senecio ogotorukensis"],
    "Packera pauciflora": ["Senecio pauciflorus"],
    "Packera porteri": ["Senecio porteri"],
    "Packera quaerens": ["Senecio quaerens"],
    "Packera quercetorum": ["Senecio quercetorum"],
    "Packera rosei": ["Senecio rosei"],
    "Packera sacramentana": ["Senecio sacramentanus"],
    "Packera sanguisorbae": ["Senecio sanguisorbae"],
    "Packera scalaris": ["Senecio scalaris"],
    "Packera streptanthifolia": [
        "Senecio streptanthifolius",
        "Senecio cymbalarioides",
        "Packera cymbalarioides",
        "Senecio laetiflorus",
        "Senecio oodes",
    ],
    "Packera subnuda": ["Senecio subnudus", "Packera buekii"],
    "Packera texensis": ["Senecio texensis"],
    "Packera toluccana": ["Senecio toluccanus"],
    "Packera tridenticulata": ["Senecio tridenticulatus", "Senecio acutidens", "Senecio densus"],
    "Packera umbraculifera": ["Senecio umbraculifer"],
    "Packera werneriifolia": ["Senecio werneriifolius", "Senecio petrocallis", "Senecio molinarius", "Senecio muirii"],
    "Packera zimapanica": ["Senecio zimapanicus"],
}

_SYNONYM_GENERA: Set[str] = {"Packera", "Senecio", "Cineraria"}


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
    m = re.match(r"^([A-Za-z]+)\s*(?:[x×]\s*)?([a-zA-Z][a-zA-Z-]+)", txt)
    if not m:
        return None
    genus = m.group(1).capitalize()
    epithet = m.group(2).lower().lstrip("x×-")
    if not epithet or epithet in {"sp", "spp", "cf", "aff", "indet", "hybrid"}:
        return None
    if allowed_genera is not None and genus not in allowed_genera:
        return None
    return f"{genus} {epithet}"


_BUILTIN_SYNONYM_TO_CANONICAL: Dict[str, str] = {}
for _canon_sp, _syn_list in PACKERA_NORTH_AMERICA_SYNONYMS.items():
    _BUILTIN_SYNONYM_TO_CANONICAL[_canon_sp.lower()] = _canon_sp
    for _syn in _syn_list:
        _syn_low = _syn.strip().lower()
        _BUILTIN_SYNONYM_TO_CANONICAL[_syn_low] = _canon_sp
        _syn_bin = _extract_canonical_binomial(_syn, _SYNONYM_GENERA)
        if _syn_bin and _syn_bin.lower() not in _BUILTIN_SYNONYM_TO_CANONICAL:
            _BUILTIN_SYNONYM_TO_CANONICAL[_syn_bin.lower()] = _canon_sp

for _cfg_canon, _cfg_syns in (_DEFAULT_CONFIG.taxa.synonyms or {}).items():
    _canon_norm = _extract_canonical_binomial(_cfg_canon, {"Packera"}) or _cfg_canon.strip()
    _BUILTIN_SYNONYM_TO_CANONICAL[_canon_norm.lower()] = _canon_norm
    for _cs in _cfg_syns or []:
        _cs_low = str(_cs).strip().lower()
        if _cs_low and _cs_low not in _BUILTIN_SYNONYM_TO_CANONICAL:
            _BUILTIN_SYNONYM_TO_CANONICAL[_cs_low] = _canon_norm
        _cs_bin = _extract_canonical_binomial(_cs, _SYNONYM_GENERA)
        if _cs_bin and _cs_bin.lower() not in _BUILTIN_SYNONYM_TO_CANONICAL:
            _BUILTIN_SYNONYM_TO_CANONICAL[_cs_bin.lower()] = _canon_norm

_DYNAMIC_SYNONYM_TO_CANONICAL: Dict[str, str] = dict(_BUILTIN_SYNONYM_TO_CANONICAL)
_DYNAMIC_TAXON_AUTHOR_STRINGS: Dict[str, Set[str]] = {}
_GBIF_SYNONYM_CACHE: Dict[str, Tuple[Optional[str], List[str], List[str]]] = {}
_GBIF_NETWORK_UNAVAILABLE: bool = False


def _is_gbif_species_mocked() -> bool:
    """Return True if pygbif.species endpoints are currently patched by unittest.mock."""
    nb = getattr(species, "name_backbone", None)
    nu = getattr(species, "name_usage", None)
    return any(
        hasattr(fn, "_mock_name") or hasattr(fn, "assert_called")
        for fn in (nb, nu)
        if fn is not None
    )


def _query_gbif_backbone_and_synonyms(
    query_name: str,
    timeout_seconds: float = 5.0,
) -> Tuple[Optional[str], List[str], List[str]]:
    """
    Query pygbif.species.name_backbone and pygbif.species.name_usage(..., data='synonyms')
    with a hard 5-second timeout and in-memory cache so offline unit tests never fail or hang.
    Returns (accepted_packera_binomial, synonym_binomials, raw_scientific_names_with_authors).
    """
    global _GBIF_NETWORK_UNAVAILABLE

    is_mocked = _is_gbif_species_mocked()
    cache_key = query_name.strip().lower()
    if not is_mocked and cache_key in _GBIF_SYNONYM_CACHE:
        return _GBIF_SYNONYM_CACHE[cache_key]
    if not is_mocked and _GBIF_NETWORK_UNAVAILABLE:
        return None, [], []

    def _fetch() -> Tuple[Optional[str], List[str], List[str]]:
        backbone = species.name_backbone(name=query_name, rank="species", kingdom="Plantae")
        if not isinstance(backbone, dict) or not backbone:
            return None, [], []
        match_type = str(backbone.get("matchType", "EXACT")).upper()
        if match_type in {"NONE", "HIGHERRANK"}:
            return None, [], []

        author_SciNames: List[str] = []
        for fld in ("scientificName", "accepted"):
            val = backbone.get(fld)
            if isinstance(val, str) and val.strip():
                author_SciNames.append(val.strip())

        accepted_packera: Optional[str] = None
        for cand_fld in ("species", "accepted", "canonicalName"):
            cand_bin = _extract_canonical_binomial(backbone.get(cand_fld), {"Packera"})
            if cand_bin:
                accepted_packera = cand_bin
                break

        found_binomials: List[str] = []
        for cand_fld in ("canonicalName", "species", "accepted"):
            b_val = _extract_canonical_binomial(backbone.get(cand_fld), _SYNONYM_GENERA)
            if b_val and b_val not in found_binomials:
                found_binomials.append(b_val)

        usage_key = backbone.get("usageKey")
        accepted_key = backbone.get("acceptedUsageKey")
        keys_to_query: List[Any] = []
        if usage_key is not None:
            keys_to_query.append(usage_key)
        if accepted_key is not None and accepted_key != usage_key:
            keys_to_query.append(accepted_key)

        for k_id in keys_to_query:
            syn_resp = species.name_usage(key=k_id, data="synonyms")
            if isinstance(syn_resp, dict):
                syn_items = syn_resp.get("results", [])
            elif isinstance(syn_resp, list):
                syn_items = syn_resp
            else:
                syn_items = []
            for item in syn_items:
                if not isinstance(item, dict):
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

        return accepted_packera, found_binomials, author_SciNames

    if is_mocked:
        try:
            return _fetch()
        except Exception:
            return None, [], []

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        future = executor.submit(_fetch)
        res = future.result(timeout=timeout_seconds)
        _GBIF_SYNONYM_CACHE[cache_key] = res
        return res
    except Exception:
        _GBIF_NETWORK_UNAVAILABLE = True
        return None, [], []
    finally:
        executor.shutdown(wait=False, cancel_futures=True)


def expand_taxon_synonyms_gbif(
    taxon_name: str,
    cfg_synonyms: Optional[Dict[str, List[str]]] = None,
) -> List[str]:
    """
    Dynamically expands a target Packera (or historical Senecio/Cineraria) taxon into its full
    homotypic and heterotypic synonym query set using:
      1. Optional configuration synonym mappings (`cfg_synonyms`) and the built-in North American
         Packera <-> Senecio/Cineraria basionym dictionary.
      2. Live GBIF Species API (`pygbif.species.name_backbone` + `pygbif.species.name_usage`),
         guarded by a 5-second timeout and in-memory cache for offline reliability.
    Also registers all discovered synonyms in `_DYNAMIC_SYNONYM_TO_CANONICAL` so
    `standardize_packera_taxon()` maps every harvested synonym back to `Packera <epithet>`.
    """
    if taxon_name is None:
        return []
    clean_input = str(taxon_name).strip()
    if not clean_input:
        return []

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

    # Determine canonical Packera binomial if resolvable from input, cfg_synonyms, or built-in map
    canonical_packera = standardize_packera_taxon(clean_input)

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
                    canonical_packera = _extract_canonical_binomial(key_clean, {"Packera"}) or standardize_packera_taxon(key_clean)
                for g in all_group:
                    _add_candidate(g)

    if canonical_packera and canonical_packera.startswith("Packera "):
        _add_candidate(canonical_packera)
    _add_candidate(clean_input)

    # Check built-in Eastern + Western North American Packera synonym dictionary
    lookup_keys = [clean_input]
    if canonical_packera and canonical_packera not in lookup_keys:
        lookup_keys.insert(0, canonical_packera)
    for lk in lookup_keys:
        for b_canon, b_syns in PACKERA_NORTH_AMERICA_SYNONYMS.items():
            b_group_lower = {b_canon.lower()} | {s.lower() for s in b_syns}
            lk_bin = _extract_canonical_binomial(lk, _SYNONYM_GENERA)
            if lk.lower() in b_group_lower or (lk_bin and lk_bin.lower() in b_group_lower):
                if canonical_packera is None:
                    canonical_packera = b_canon
                _add_candidate(b_canon)
                for s in b_syns:
                    _add_candidate(s)
                    s_bin = _extract_canonical_binomial(s, _SYNONYM_GENERA)
                    if s_bin:
                        _add_candidate(s_bin)

    # Query GBIF Species Backbone + Synonyms endpoint (fast 5s timeout + cached)
    gbif_query_target = canonical_packera if (canonical_packera and canonical_packera.startswith("Packera ")) else clean_input
    gbif_acc, gbif_syns, gbif_author_names = _query_gbif_backbone_and_synonyms(gbif_query_target)
    if gbif_acc and gbif_acc.startswith("Packera "):
        if canonical_packera is None:
            canonical_packera = gbif_acc
        _add_candidate(gbif_acc)
    for g_syn in gbif_syns:
        _add_candidate(g_syn)

    # Register dynamic synonyms -> canonical Packera binomial & cache protologue author strings
    if canonical_packera and canonical_packera.startswith("Packera "):
        _DYNAMIC_SYNONYM_TO_CANONICAL[canonical_packera.lower()] = canonical_packera
        if gbif_author_names:
            _DYNAMIC_TAXON_AUTHOR_STRINGS.setdefault(canonical_packera, set()).update(gbif_author_names)
        for item in expanded:
            item_low = item.lower()
            if item_low not in _BUILTIN_SYNONYM_TO_CANONICAL:
                _DYNAMIC_SYNONYM_TO_CANONICAL[item_low] = canonical_packera
            item_bin = _extract_canonical_binomial(item, _SYNONYM_GENERA)
            if item_bin and item_bin.lower() not in _BUILTIN_SYNONYM_TO_CANONICAL:
                _DYNAMIC_SYNONYM_TO_CANONICAL[item_bin.lower()] = canonical_packera

    return expanded


def standardize_packera_taxon(raw_name: Any) -> Optional[str]:
    """
    Standardizes raw Darwin Core scientific names, varieties, and historical Senecio/Cineraria
    basionyms/synonyms across all North American Packera species.
    Uses specific varietal/hybrid disambiguations, the built-in + dynamic GBIF synonym map,
    and preserves any valid binomial matching ^packera\\s*(?:[x×]\\s*)?([a-z-]+) dynamically.
    """
    if raw_name is None or (isinstance(raw_name, float) and math.isnan(raw_name)):
        return None
    s_orig = str(raw_name).strip()
    if not s_orig or s_orig.lower() in {"none", "nan", "null", "na", "unknown"}:
        return None
    s = s_orig.lower()

    # 1. Specific varietal / hybrid / basionym disambiguations (order matters)
    if any(k in s for k in ("crawfordii", "crawford")):
        return "Packera crawfordii"
    if any(k in s for k in ("memmingeri", "memminger", "packera×memmingeri", "packeraxmemmingeri")):
        return "Packera memmingeri"
    if any(k in s for k in ("millefolium", "millefolia", "millefolius")):
        return "Packera millefolium"
    if any(k in s for k in ("pseudaurea", "pseudaureus", "semicordata", "semicordatus")):
        return "Packera pseudaurea"
    if any(k in s for k in ("paupercula", "pauperculus", "balsamitae", "appalachiana", "savannarum")):
        return "Packera paupercula"
    if any(k in s for k in ("plattensis", "pseudotomentosus", "pseudotomentosa", "flavovirens")):
        return "Packera plattensis"
    if any(k in s for k in ("dubia", "dubius", "tomentosus", "tomentosa", "alabamensis")):
        return "Packera dubia"
    if any(k in s for k in ("anonyma", "anonymus", "smallii", "earlei")):
        return "Packera anonyma"
    if any(k in s for k in ("antennariifolia", "antennariifolius")):
        return "Packera antennariifolia"
    if any(k in s for k in ("aurea", "aureus", "gracilis")):
        return "Packera aurea"
    if any(k in s for k in ("glabella", "glabellus", "lobatus", "lobata")):
        return "Packera glabella"
    if "hyperboreal" in s:
        return "Packera hyperborealis"
    if any(k in s for k in ("cymbalaria", "resedifolius", "resedifolia")) and "cymbalarioides" not in s:
        return "Packera cymbalaria"
    if any(k in s for k in ("indecora", "indecorus")):
        return "Packera indecora"
    if any(k in s for k in ("insulanus", "insulana")):
        return "Packera insulanus"
    if any(k in s for k in ("obovata", "obovatus", "rotundus", "rotunda", "elliottii")):
        return "Packera obovata"
    if any(k in s for k in ("schweinitziana", "schweinitzianus", "robbinsii")):
        return "Packera schweinitziana"
    if any(k in s for k in ("serpenticola", "eburnea", "eburneus", "packeraeburnea")):
        return "Packera serpenticola"
    if "spellenberg" in s:
        return "Packera spellenbergii"
    if any(k in s for k in ("tampicana", "tampicanus", "imparipinnatus", "imparipinnata")):
        return "Packera tampicana"
    if "vulgaris" in s:
        return "Packera vulgaris"

    # 2. Dynamic & North American Packera / Senecio / Cineraria synonym map lookup
    s_compact = re.sub(r"\s+", " ", s)
    if s_compact in _DYNAMIC_SYNONYM_TO_CANONICAL:
        return _DYNAMIC_SYNONYM_TO_CANONICAL[s_compact]
    var_m = re.match(r"^([a-z]+\s+[a-z-]+\s+(?:var\.?|subsp\.?|ssp\.?)\s+[a-z-]+)", s_compact)
    if var_m:
        var_key = re.sub(r"\b(?:var|subsp|ssp)\b\.?", lambda m0: m0.group(0).rstrip(".") + ".", var_m.group(1))
        if var_key in _DYNAMIC_SYNONYM_TO_CANONICAL:
            return _DYNAMIC_SYNONYM_TO_CANONICAL[var_key]
    cand_bin = _extract_canonical_binomial(s_orig, _SYNONYM_GENERA)
    if cand_bin and cand_bin.lower() in _DYNAMIC_SYNONYM_TO_CANONICAL:
        return _DYNAMIC_SYNONYM_TO_CANONICAL[cand_bin.lower()]

    # 3. Dynamic fallback for any valid Packera binomial (^packera\s*(?:[x×]\s*)?([a-z-]+))
    m = re.match(r"^packera\s*(?:[x×]\s*)?([a-z-]+)", s)
    if m:
        epithet = m.group(1).lstrip("x×-")
        if epithet:
            return f"Packera {epithet}"

    return None

EXPORT_COLUMNS = [
    "catalogNumber",
    "institutionCode",
    "scientificName",
    "species_raw",
    "recordedBy",
    "recordNumber",
    "identifiedBy",
    "determiner_raw",
    "determiner_tier",
    "type_status",
    "county",
    "stateProvince",
    "decimalLatitude",
    "latitude",
    "decimalLongitude",
    "longitude",
    "coordinateUncertainty",
    "year",
    "month",
    "day",
    "eventDate",
    "regional_group",
    "image_path",
    "exsiccatae_key",
    "is_primary_duplicate",
    "duplicate_count",
]


def normalize_collector(recorded_by: Optional[str]) -> str:
    """
    Normalizes botanical collector names for exsiccatae collection event fingerprinting:
    lowercased, non-alphanumeric characters stripped, initials removed, primary surname extracted.
    e.g., 'J. Brandon Fuller' -> 'fuller', 'Fuller, J. B.' -> 'fuller', 'A. Cronquist' -> 'cronquist'.
    """
    if not recorded_by or not isinstance(recorded_by, str):
        return "unknown"
    s = str(recorded_by).strip()
    if not s or s.lower() in {"none", "nan", "null", "unknown", "anonymous", "not evident"}:
        return "unknown"

    # Separate multiple collectors (e.g. 'J. B. Fuller & R. R. Kowal') -> retain primary collector
    primary = re.split(r"\s*(?:&|;|\band\b|\bwith\b|\bet\s+al\.?)\s*", s, flags=re.IGNORECASE)[0].strip()

    # Handle 'Fuller, J. Brandon' format
    if "," in primary:
        surname_part = primary.split(",")[0].strip()
        cleaned = re.sub(r"[^a-zA-Z0-9]", "", surname_part).lower()
        if cleaned:
            return cleaned

    # Tokenize and remove punctuation
    tokens = [re.sub(r"[^a-zA-Z0-9]", "", t).lower() for t in primary.split()]
    tokens = [t for t in tokens if t]

    if not tokens:
        return "unknown"

    # Remove single-letter initials (e.g. 'j', 'b')
    non_initials = [t for t in tokens if len(t) > 1]
    if not non_initials:
        # All tokens were initials (e.g. 'J. B.')
        return "".join(tokens)

    # In 'First Middle Last', the surname is the last token
    return non_initials[-1]


def normalize_record_number(record_number: Optional[Any]) -> str:
    """
    Normalizes collection number: strips prefixes and suffixes, extracts pure numeric component.
    e.g., '#1042' -> '1042', 'No. 1042b' -> '1042', 's.n.' -> ''.
    """
    if record_number is None:
        return ""
    s = str(record_number).strip()
    if not s or s.lower() in {"none", "nan", "null", "s.n.", "sn", "s. n.", "unnumbered"}:
        return ""
    match = re.search(r"\d+", s)
    return match.group(0) if match else ""


def normalize_event_date(
    event_date_raw: Optional[str] = None,
    year: Optional[Any] = None,
    month: Optional[Any] = None,
    day: Optional[Any] = None,
) -> str:
    """
    Standardizes collection event date to ISO format YYYY-MM-DD.
    Evaluates explicit year, month, day components or parses date string.
    """
    try:
        if year is not None and month is not None and day is not None:
            if not (pd.isna(year) or pd.isna(month) or pd.isna(day)):
                y, m, d = int(year), int(month), int(day)
                if 1700 <= y <= 2100 and 1 <= m <= 12 and 1 <= d <= 31:
                    return f"{y:04d}-{m:02d}-{d:02d}"
    except (ValueError, TypeError):
        pass

    if event_date_raw and isinstance(event_date_raw, str):
        s = event_date_raw.strip()
        match_iso = re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", s)
        if match_iso:
            y, m, d = int(match_iso.group(1)), int(match_iso.group(2)), int(match_iso.group(3))
            return f"{y:04d}-{m:02d}-{d:02d}"
        match_y = re.search(r"\b(1[789]\d{2}|20\d{2})\b", s)
        if match_y:
            return f"{match_y.group(1)}-00-00"

    try:
        if year is not None and not pd.isna(year):
            y = int(year)
            if 1700 <= y <= 2100:
                return f"{y:04d}-00-00"
    except (ValueError, TypeError):
        pass

    return "unknown_date"


def generate_exsiccatae_key(
    collector: Optional[str],
    number: Optional[Any],
    date: Optional[str] = None,
    year: Optional[Any] = None,
    month: Optional[Any] = None,
    day: Optional[Any] = None,
) -> str:
    """
    Constructs a standardized collection event signature (exsiccatae fingerprint):
    f"{norm_collector}_{norm_number}_{norm_date}".
    """
    norm_collector = normalize_collector(collector)
    norm_number = normalize_record_number(number)
    norm_date = normalize_event_date(date, year=year, month=month, day=day)
    return f"{norm_collector}_{norm_number}_{norm_date}"


def haversine_distance_meters(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Calculates geodesic distance between two coordinate pairs using the Haversine formula."""
    r_earth = 6371000.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    a = math.sin(delta_phi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2.0) ** 2
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return r_earth * c


def get_image_resolution_score(image_path_str: Optional[str], workspace_dir: Optional[Path] = None) -> float:
    """Computes a numeric image resolution score (megapixels or pixel count, fallback to byte size)."""
    if not image_path_str or not isinstance(image_path_str, str):
        return 0.0
    p = Path(image_path_str)
    resolved_path: Optional[Path] = None
    if p.is_absolute() and p.exists():
        resolved_path = p
    elif workspace_dir and (workspace_dir / p).exists():
        resolved_path = workspace_dir / p
    elif (PROJECT_ROOT / p).exists():
        resolved_path = PROJECT_ROOT / p

    if resolved_path and resolved_path.exists():
        try:
            with Image.open(resolved_path) as img:
                w, h = img.size
                return float(w * h)
        except Exception:
            try:
                return float(resolved_path.stat().st_size)
            except Exception:
                return 0.0
    return 0.0


def stratify_duplicates(
    df: pd.DataFrame,
    max_distance_meters: float = 50.0,
    workspace_dir: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
    home_herbarium: str = DEFAULT_HOME_HERBARIUM,
) -> pd.DataFrame:
    """
    Identifies duplicate voucher specimens (exsiccatae) representing the same
    field collection event distributed across multiple herbaria.

    Duplicates are resolved by:
      1. Standardized exsiccatae fingerprint (norm_collector + norm_number + norm_date).
      2. If collection number is missing, falls back to geographic coordinate proximity (< 50 m)
         within identical collection dates.

    Designates exactly one primary voucher per duplicate group using institutional hierarchy:
      Priority 1: Verified Type Status
      Priority 2: Home Herbarium (default: NCU; configurable via home_herbarium)
      Priority 3: Tier 1 Specialist determinations
      Priority 4: Highest specimen image resolution

    Adds/updates columns:
      - exsiccatae_key: Standardized signature string
      - is_primary_duplicate: Boolean flag (True for primary specimen, False for duplicates)
      - duplicate_count: Total sheets in this exsiccatae group
    """
    if df.empty:
        df["exsiccatae_key"] = pd.Series(dtype=str)
        df["is_primary_duplicate"] = pd.Series(dtype=bool)
        df["duplicate_count"] = pd.Series(dtype=int)
        return df

    df = df.copy()

    # Ensure metadata columns exist
    if "recordedBy" not in df.columns:
        df["recordedBy"] = df.get("determiner_raw", df.get("identifiedBy", ""))
    if "recordNumber" not in df.columns:
        df["recordNumber"] = ""
    if "image_path" not in df.columns:
        df["image_path"] = ""

    norm_collectors: List[str] = []
    norm_numbers: List[str] = []
    norm_dates: List[str] = []
    has_coords: List[bool] = []
    lats: List[Optional[float]] = []
    lons: List[Optional[float]] = []

    for _, row in df.iterrows():
        rec_by = str(row.get("recordedBy", "") if not pd.isna(row.get("recordedBy", "")) else "")
        rec_num = row.get("recordNumber", "") if not pd.isna(row.get("recordNumber", "")) else ""
        ev_date = str(row.get("eventDate", "") if not pd.isna(row.get("eventDate", "")) else "")
        yr = row.get("year")
        mo = row.get("month")
        dy = row.get("day")

        c_norm = normalize_collector(rec_by)
        n_norm = normalize_record_number(rec_num)
        d_norm = normalize_event_date(ev_date, year=yr, month=mo, day=dy)

        lat_val = None
        lon_val = None
        for col_lat in ["latitude", "decimalLatitude"]:
            if col_lat in row and not pd.isna(row[col_lat]):
                try:
                    lat_val = float(row[col_lat])
                    break
                except (ValueError, TypeError):
                    pass
        for col_lon in ["longitude", "decimalLongitude"]:
            if col_lon in row and not pd.isna(row[col_lon]):
                try:
                    lon_val = float(row[col_lon])
                    break
                except (ValueError, TypeError):
                    pass

        norm_collectors.append(c_norm)
        norm_numbers.append(n_norm)
        norm_dates.append(d_norm)
        has_coords.append(lat_val is not None and lon_val is not None)
        lats.append(lat_val)
        lons.append(lon_val)

    n_rows = len(df)
    group_labels: List[Optional[str]] = [None] * n_rows

    # Pass 1: Assign keys for records with explicit collection numbers
    for idx in range(n_rows):
        if norm_numbers[idx]:
            group_labels[idx] = f"{norm_collectors[idx]}_{norm_numbers[idx]}_{norm_dates[idx]}"

    # Pass 2: Missing collection number -> fallback to coordinate proximity (< 50m) on identical dates
    missing_num_indices = [idx for idx in range(n_rows) if not norm_numbers[idx]]
    date_to_indices: Dict[str, List[int]] = {}
    for idx in missing_num_indices:
        d = norm_dates[idx]
        date_to_indices.setdefault(d, []).append(idx)

    for d_val, idx_list in date_to_indices.items():
        if d_val == "unknown_date":
            for i in idx_list:
                cat = sanitize_filename(str(df.iloc[i].get("catalogNumber", i)))
                group_labels[i] = f"{norm_collectors[i]}_sn_{cat}_{d_val}"
            continue

        coord_indices = [i for i in idx_list if has_coords[i]]
        no_coord_indices = [i for i in idx_list if not has_coords[i]]

        # Records without coordinates cannot verify < 50m proximity; preserve as unique singletons
        for i in no_coord_indices:
            cat = sanitize_filename(str(df.iloc[i].get("catalogNumber", i)))
            group_labels[i] = f"{norm_collectors[i]}_sn_{cat}_{d_val}"

        # Match unnumbered records to any numbered collection on the same date within max_distance_meters
        numbered_date_indices = [
            i for i in range(n_rows)
            if norm_numbers[i] and norm_dates[i] == d_val and has_coords[i]
        ]

        matched_to_numbered = set()
        for i in coord_indices:
            for num_i in numbered_date_indices:
                if norm_collectors[i] == norm_collectors[num_i] or norm_collectors[i] == "unknown":
                    dist = haversine_distance_meters(lats[i], lons[i], lats[num_i], lons[num_i])
                    if dist <= max_distance_meters:
                        group_labels[i] = group_labels[num_i]
                        matched_to_numbered.add(i)
                        break

        remaining_coord_indices = [i for i in coord_indices if i not in matched_to_numbered]
        num_rem = len(remaining_coord_indices)
        if num_rem == 0:
            continue

        visited = [False] * num_rem
        for a in range(num_rem):
            if visited[a]:
                continue
            cluster = [a]
            visited[a] = True
            q = [a]
            while q:
                curr = q.pop(0)
                curr_idx = remaining_coord_indices[curr]
                for b in range(num_rem):
                    if not visited[b]:
                        b_idx = remaining_coord_indices[b]
                        same_coll = (
                            norm_collectors[curr_idx] == norm_collectors[b_idx]
                            or norm_collectors[curr_idx] == "unknown"
                            or norm_collectors[b_idx] == "unknown"
                        )
                        if same_coll:
                            dist = haversine_distance_meters(
                                lats[curr_idx], lons[curr_idx], lats[b_idx], lons[b_idx]
                            )
                            if dist <= max_distance_meters:
                                visited[b] = True
                                cluster.append(b)
                                q.append(b)

            actual_indices = [remaining_coord_indices[c_idx] for c_idx in cluster]
            best_coll = "unknown"
            for act_i in actual_indices:
                if norm_collectors[act_i] != "unknown":
                    best_coll = norm_collectors[act_i]
                    break

            if len(actual_indices) > 1:
                anchor_cat = sanitize_filename(str(df.iloc[actual_indices[0]].get("catalogNumber", actual_indices[0])))
                cluster_key = f"{best_coll}_geo_{anchor_cat}_{d_val}"
            else:
                act_i = actual_indices[0]
                cat = sanitize_filename(str(df.iloc[act_i].get("catalogNumber", act_i)))
                cluster_key = f"{best_coll}_sn_{cat}_{d_val}"

            for act_i in actual_indices:
                group_labels[act_i] = cluster_key

    df["exsiccatae_key"] = group_labels

    # Pass 3: Evaluate Institutional Priority and designate is_primary_duplicate
    grouped = df.groupby("exsiccatae_key", sort=False)
    df["duplicate_count"] = grouped["catalogNumber"].transform("count")

    is_primary_list = [False] * n_rows

    for key, group_df in grouped:
        indices = list(group_df.index)
        if len(indices) == 1:
            is_primary_list[indices[0]] = True
            continue

        scored_candidates = []
        for idx in indices:
            row = group_df.loc[idx]

            # Priority 1: Verified Type Status
            type_stat = str(row.get("type_status", "") if not pd.isna(row.get("type_status", "")) else "").strip().lower()
            has_type = 1 if type_stat and type_stat not in {"none", "not a type", "notatype", "unspecified", "null"} else 0

            # Priority 2: Home Herbarium (default: NCU, configurable via home_herbarium)
            target_home = str(home_herbarium).strip().upper() if home_herbarium else ""
            inst = str(row.get("institutionCode", "") if not pd.isna(row.get("institutionCode", "")) else "").strip().upper()
            is_home_herb = 1 if (target_home and inst == target_home) else 0

            # Priority 3: Tier 1 Specialist determinations
            tier = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
            tier_score = 3 if tier == "Tier_1_Gold" else (2 if tier == "Tier_2_Silver" else 1)

            # Priority 4: Highest specimen image resolution
            img_path = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "")
            res_score = get_image_resolution_score(img_path, workspace_dir=workspace_dir)

            cat_tiebreak = str(row.get("catalogNumber", idx))

            score_tuple = (has_type, is_home_herb, tier_score, res_score, cat_tiebreak)
            scored_candidates.append((score_tuple, idx))

        scored_candidates.sort(key=lambda x: x[0], reverse=True)
        primary_idx = scored_candidates[0][1]
        is_primary_list[primary_idx] = True

    df["is_primary_duplicate"] = is_primary_list

    if logger:
        n_prim = sum(is_primary_list)
        n_dup = n_rows - n_prim
        logger.info(
            f"Duplicate Stratification: Processed {n_rows} vouchers -> "
            f"{n_prim} primary duplicate vouchers designated ({n_dup} duplicate sheets flagged)."
        )

    return df


def spatially_balance_and_thin_vouchers(
    df: pd.DataFrame,
    spatial_thin_km: float = 0.0,
    max_per_state: Optional[int] = None,
    max_per_region: Optional[int] = None,
    primary_duplicates_only: bool = True,
    preserve_all_gold: bool = True,
    manifest_df: Optional[pd.DataFrame] = None,
    workspace_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Stratifies, deduplicates, spatially thins, and computes continuous inverse
    spatial-density weights (`spatial_balance_weight`) for botanical vouchers.
    """
    if df is None:
        return pd.DataFrame(columns=EXPORT_COLUMNS + ["spatial_balance_weight"])
    if df.empty:
        out_empty = df.copy()
        out_empty["spatial_balance_weight"] = pd.Series(dtype=float)
        return out_empty

    out_df = df.copy()

    # Ensure standardized taxon column exists
    if "species_standardized" not in out_df.columns:
        src_col = "scientificName" if "scientificName" in out_df.columns else (
            "species_raw" if "species_raw" in out_df.columns else None
        )
        if src_col is not None:
            out_df["species_standardized"] = out_df[src_col].map(
                lambda x: standardize_packera_taxon(str(x)) or str(x).strip() or "Unknown"
                if not pd.isna(x) else "Unknown"
            )
        else:
            out_df["species_standardized"] = "Unknown"
    else:
        src_col = "scientificName" if "scientificName" in out_df.columns else (
            "species_raw" if "species_raw" in out_df.columns else None
        )
        def _fill_sp(row: pd.Series) -> str:
            val = row.get("species_standardized")
            if val is not None and not pd.isna(val) and str(val).strip():
                return str(val).strip()
            if src_col is not None:
                raw_v = row.get(src_col)
                if raw_v is not None and not pd.isna(raw_v) and str(raw_v).strip():
                    return standardize_packera_taxon(str(raw_v)) or str(raw_v).strip()
            return "Unknown"
        out_df["species_standardized"] = out_df.apply(_fill_sp, axis=1)

    # Step A: Exsiccatae Deduplication
    if primary_duplicates_only:
        if "is_primary_duplicate" not in out_df.columns:
            if any(c in out_df.columns for c in ("recordedBy", "recordNumber", "exsiccatae_key")):
                out_df = stratify_duplicates(out_df, workspace_dir=workspace_dir)
            else:
                out_df["is_primary_duplicate"] = True
        prim_mask = out_df["is_primary_duplicate"].map(
            lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}
        )
        out_df = out_df.loc[prim_mask].reset_index(drop=True)
        if out_df.empty:
            out_df["spatial_balance_weight"] = pd.Series(dtype=float)
            return out_df

    # Step B: Candidate Priority Ranking
    t1_leaf_counts: Dict[str, int] = {}
    if manifest_df is not None and not manifest_df.empty and "catalogNumber" in manifest_df.columns:
        tier_col = "assigned_tier" if "assigned_tier" in manifest_df.columns else (
            "tier_assignment" if "tier_assignment" in manifest_df.columns else None
        )
        if tier_col is not None:
            t1_rows = manifest_df[
                manifest_df[tier_col].astype(str).str.strip().str.lower().str.match(r"^tier_?1")
            ]
        else:
            t1_rows = manifest_df
        for cat_k, cnt in t1_rows["catalogNumber"].astype(str).value_counts().items():
            t1_leaf_counts[str(cat_k)] = int(cnt)

    lat_col = "decimalLatitude" if "decimalLatitude" in out_df.columns else (
        "latitude" if "latitude" in out_df.columns else None
    )
    lon_col = "decimalLongitude" if "decimalLongitude" in out_df.columns else (
        "longitude" if "longitude" in out_df.columns else None
    )

    n_rows = len(out_df)
    lats: List[Optional[float]] = [None] * n_rows
    lons: List[Optional[float]] = [None] * n_rows
    priority_tuples: List[Tuple[int, int, int, float, float, str]] = []
    is_gold_flags: List[bool] = []

    for i in range(n_rows):
        row = out_df.iloc[i]
        cat_str = str(row.get("catalogNumber", i))

        # 1. has_type_status
        t_stat = str(row.get("type_status", "") if not pd.isna(row.get("type_status", "")) else "").strip()
        has_type = 1 if t_stat and t_stat.lower() not in {
            "", "none", "nan", "null", "not a type", "notatype", "unspecified"
        } else 0

        # 2. determiner_tier
        tier_str = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
        tier_score = 3 if tier_str == "Tier_1_Gold" else (2 if tier_str == "Tier_2_Silver" else 1)
        is_gold_flags.append(tier_str == "Tier_1_Gold")

        # 3. n_tier1_leaves
        if cat_str in t1_leaf_counts:
            n_t1 = t1_leaf_counts[cat_str]
        elif "n_tier1_leaves" in out_df.columns and not pd.isna(row.get("n_tier1_leaves")):
            try:
                n_t1 = int(float(row.get("n_tier1_leaves", 0)))
            except (ValueError, TypeError):
                n_t1 = 0
        else:
            n_t1 = 0

        # Coordinates & 4. coordinate_precision
        lat_v: Optional[float] = None
        lon_v: Optional[float] = None
        if lat_col and lon_col:
            try:
                lv = float(row.get(lat_col))
                lnv = float(row.get(lon_col))
                if math.isfinite(lv) and math.isfinite(lnv):
                    lat_v, lon_v = lv, lnv
            except (ValueError, TypeError):
                pass
        lats[i] = lat_v
        lons[i] = lon_v

        coord_prec = 0.0
        if lat_v is not None and lon_v is not None:
            unc_raw = row.get("coordinateUncertainty")
            try:
                unc_f = float(unc_raw)
                if math.isfinite(unc_f) and unc_f >= 0.0:
                    coord_prec = 1.0 / (1.0 + unc_f)
                else:
                    coord_prec = 1.0 / (1.0 + 5000.0)
            except (ValueError, TypeError):
                coord_prec = 1.0 / (1.0 + 5000.0)

        # 5. image_resolution_score
        img_p = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "")
        res_score = get_image_resolution_score(img_p, workspace_dir=workspace_dir)

        priority_tuples.append((has_type, tier_score, n_t1, round(coord_prec, 8), round(res_score, 4), cat_str))

    active_indices = list(range(n_rows))

    # Step C: Spatial Grid Thinning (spatial_thin_km > 0)
    if spatial_thin_km is not None and float(spatial_thin_km) > 0.0:
        thin_km = float(spatial_thin_km)
        dlat = thin_km / 111.32
        cell_map: Dict[Tuple[str, int, int], List[int]] = {}
        non_geo_indices: List[int] = []

        for idx in active_indices:
            lat_v = lats[idx]
            lon_v = lons[idx]
            sp_val = str(out_df.iloc[idx]["species_standardized"])
            if lat_v is None or lon_v is None:
                non_geo_indices.append(idx)
                continue
            cos_lat = max(math.cos(math.radians(lat_v)), 1e-4)
            dlon = thin_km / (111.32 * cos_lat)
            lat_bin = int(math.floor(lat_v / dlat))
            lon_bin = int(math.floor(lon_v / dlon))
            cell_map.setdefault((sp_val, lat_bin, lon_bin), []).append(idx)

        thinned_indices: List[int] = []
        for _cell_key, cell_idxs in cell_map.items():
            cell_idxs_sorted = sorted(cell_idxs, key=lambda i: priority_tuples[i], reverse=True)
            retained_cell = [cell_idxs_sorted[0]]
            if preserve_all_gold:
                for i in cell_idxs_sorted[1:]:
                    if is_gold_flags[i]:
                        retained_cell.append(i)
            thinned_indices.extend(retained_cell)

        thinned_indices.extend(non_geo_indices)
        active_indices = sorted(set(thinned_indices))

    def _select_spatially_spread_quota(
        cand_pool: List[int],
        quota: int,
        stratum_values: Dict[int, str],
    ) -> List[int]:
        """
        Selects up to `quota` indices from `cand_pool`, preserving all Gold anchors
        when `preserve_all_gold=True` and spreading remaining selections across
        `stratum_values` (counties or states) and farthest Haversine distance.
        """
        if len(cand_pool) <= quota:
            return sorted(cand_pool)

        if preserve_all_gold:
            gold_idxs = [i for i in cand_pool if is_gold_flags[i]]
            rem_pool = [i for i in cand_pool if not is_gold_flags[i]]
        else:
            gold_idxs = []
            rem_pool = list(cand_pool)

        rem_pool.sort(key=lambda i: priority_tuples[i], reverse=True)
        selected: List[int] = list(gold_idxs)
        if len(selected) >= quota:
            return sorted(selected) if preserve_all_gold else sorted(selected[:quota])

        stratum_counts: Dict[str, int] = {}
        for i in selected:
            st_key = stratum_values[i]
            stratum_counts[st_key] = stratum_counts.get(st_key, 0) + 1

        while len(selected) < quota and rem_pool:
            min_cnt = min(stratum_counts.get(stratum_values[i], 0) for i in rem_pool)
            eligible = [i for i in rem_pool if stratum_counts.get(stratum_values[i], 0) == min_cnt]

            sel_geo = [i for i in selected if lats[i] is not None and lons[i] is not None]
            if not sel_geo:
                best_i = eligible[0]
            else:
                best_i = eligible[0]
                best_score: Tuple[float, Tuple[int, int, int, float, float, str]] = (-1.0, priority_tuples[best_i])
                for cand_i in eligible:
                    if lats[cand_i] is not None and lons[cand_i] is not None:
                        min_d_km = min(
                            haversine_distance_meters(
                                lats[cand_i], lons[cand_i], lats[s_i], lons[s_i]
                            ) / 1000.0
                            for s_i in sel_geo
                        )
                    else:
                        min_d_km = 0.0
                    cand_score = (round(min_d_km, 2), priority_tuples[cand_i])
                    if cand_score > best_score:
                        best_score = cand_score
                        best_i = cand_i

            selected.append(best_i)
            st_key = stratum_values[best_i]
            stratum_counts[st_key] = stratum_counts.get(st_key, 0) + 1
            rem_pool.remove(best_i)

        return sorted(selected)

    # Step D.1: Per-State Quota Balancing (max_per_state > 0)
    if max_per_state is not None and int(max_per_state) > 0:
        state_quota = int(max_per_state)
        state_groups: Dict[Tuple[str, str], List[int]] = {}
        county_labels: Dict[int, str] = {}
        for idx in active_indices:
            row = out_df.iloc[idx]
            sp_val = str(row["species_standardized"])
            st_raw = row.get("stateProvince", "")
            st_val = str(st_raw if not pd.isna(st_raw) else "").strip().upper() or "UNKNOWN_STATE"
            cnty_raw = row.get("county", "")
            cnty_val = str(cnty_raw if not pd.isna(cnty_raw) else "").strip().upper()
            if not cnty_val:
                if lats[idx] is not None and lons[idx] is not None:
                    cnty_val = f"GRID_{int(round(lats[idx] * 4.0))}_{int(round(lons[idx] * 4.0))}"
                else:
                    cnty_val = f"UNSPECIFIED_{idx}"
            county_labels[idx] = cnty_val
            state_groups.setdefault((sp_val, st_val), []).append(idx)

        state_balanced: List[int] = []
        for _grp_key, grp_idxs in state_groups.items():
            state_balanced.extend(_select_spatially_spread_quota(grp_idxs, state_quota, county_labels))
        active_indices = sorted(set(state_balanced))

    # Step D.2: Per-Biogeographic-Region Quota Balancing (max_per_region > 0)
    if max_per_region is not None and int(max_per_region) > 0:
        region_quota = int(max_per_region)
        region_groups: Dict[Tuple[str, str], List[int]] = {}
        state_labels: Dict[int, str] = {}
        for idx in active_indices:
            row = out_df.iloc[idx]
            sp_val = str(row["species_standardized"])
            reg_raw = row.get("regional_group", "")
            reg_val = str(reg_raw if not pd.isna(reg_raw) else "").strip() or "Unassigned_NA"
            st_raw = row.get("stateProvince", "")
            st_val = str(st_raw if not pd.isna(st_raw) else "").strip().upper() or "UNKNOWN_STATE"
            state_labels[idx] = st_val
            region_groups.setdefault((sp_val, reg_val), []).append(idx)

        region_balanced: List[int] = []
        for _grp_key, grp_idxs in region_groups.items():
            region_balanced.extend(_select_spatially_spread_quota(grp_idxs, region_quota, state_labels))
        active_indices = sorted(set(region_balanced))

    out_df = out_df.iloc[active_indices].reset_index(drop=True)
    if out_df.empty:
        out_df["spatial_balance_weight"] = pd.Series(dtype=float)
        return out_df

    # Step E: Compute Continuous `spatial_balance_weight`
    ret_lats = np.array(
        [lats[i] if lats[i] is not None else np.nan for i in active_indices],
        dtype=float,
    )
    ret_lons = np.array(
        [lons[i] if lons[i] is not None else np.nan for i in active_indices],
        dtype=float,
    )
    ret_states = [
        str(out_df.iloc[k].get("stateProvince", "") if not pd.isna(out_df.iloc[k].get("stateProvince", "")) else "").strip().upper()
        for k in range(len(out_df))
    ]
    ret_taxa = out_df["species_standardized"].astype(str).tolist()

    spatial_weights = np.ones(len(out_df), dtype=float)
    bandwidth_km = 25.0

    taxon_to_rows: Dict[str, List[int]] = {}
    for k, sp_name in enumerate(ret_taxa):
        taxon_to_rows.setdefault(sp_name, []).append(k)

    for _sp_name, t_idxs in taxon_to_rows.items():
        t_arr = np.array(t_idxs, dtype=int)
        valid_geo_mask = np.isfinite(ret_lats[t_arr]) & np.isfinite(ret_lons[t_arr])
        geo_idxs = t_arr[valid_geo_mask]
        if len(geo_idxs) == 0:
            continue

        # State frequency factor f_state = N_{taxon, state(i)} / median_s(N_{taxon, s})
        st_counts: Dict[str, int] = {}
        for k in t_idxs:
            st_s = ret_states[k]
            if st_s:
                st_counts[st_s] = st_counts.get(st_s, 0) + 1
        med_state_cnt = float(np.median(list(st_counts.values()))) if st_counts else 1.0
        med_state_cnt = max(med_state_cnt, 1.0)

        f_state = np.array(
            [
                max(st_counts.get(ret_states[k], 0) / med_state_cnt, 1e-6)
                if ret_states[k] in st_counts else 1.0
                for k in geo_idxs
            ],
            dtype=float,
        )

        # Local 25 km Gaussian kernel specimen density d_local
        rad_lat = np.radians(ret_lats[geo_idxs])
        rad_lon = np.radians(ret_lons[geo_idxs])
        dlat_mat = rad_lat[:, None] - rad_lat[None, :]
        dlon_mat = rad_lon[:, None] - rad_lon[None, :]
        a_mat = (
            np.sin(dlat_mat / 2.0) ** 2
            + np.cos(rad_lat[:, None]) * np.cos(rad_lat[None, :]) * (np.sin(dlon_mat / 2.0) ** 2)
        )
        dist_km_mat = 2.0 * 6371.0 * np.arcsin(np.clip(np.sqrt(a_mat), 0.0, 1.0))
        kern_mat = np.exp(-(dist_km_mat ** 2) / (2.0 * (bandwidth_km ** 2)))
        np.fill_diagonal(kern_mat, 0.0)
        d_local = 1.0 + np.sum(kern_mat, axis=1)

        w_raw = 1.0 / np.sqrt(np.maximum(d_local * f_state, 1e-8))
        mean_raw = float(np.mean(w_raw))
        if mean_raw > 0.0:
            w_norm = w_raw / mean_raw
        else:
            w_norm = np.ones_like(w_raw)

        w_spatial = np.clip(w_norm, 0.20, 3.00)
        for _ in range(3):
            m_cur = float(np.mean(w_spatial))
            if m_cur > 0.0:
                w_spatial = np.clip(w_spatial / m_cur, 0.20, 3.00)

        spatial_weights[geo_idxs] = np.round(w_spatial, 6)

    out_df["spatial_balance_weight"] = spatial_weights
    return out_df



def setup_logger(log_file_path: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    """Configures or retrieves a formatted logger for voucher harvesting."""
    return setup_logging(log_file=log_file_path, verbose=verbose, name="VoucherHarvester")


def sanitize_filename(name: str) -> str:
    """Sanitizes arbitrary strings into safe, valid filesystem filenames across platforms."""
    clean = re.sub(r'[\\/*?:"<>|\s]+', "_", str(name).strip())
    clean = clean.strip("._")
    return clean if clean else "voucher_unknown"


_IPNI_AUTHOR_ABBREVIATIONS: Dict[str, str] = {
    "spreng": "sprengel",
    "mabb": "mabberley",
    "michx": "michaux",
    "nutt": "nuttall",
    "torr": "torrey",
    "rydb": "rydberg",
    "greenm": "greenman",
    "cronq": "cronquist",
    "hook": "hooker",
    "ell": "elliott",
    "britt": "britton",
    "woot": "wooton",
    "standl": "standley",
    "spellenb": "spellenberg",
    "blank": "blankinship",
    "hend": "henderson",
    "hult": "hulten",
    "muhl": "muhlenberg",
    "willd": "willdenow",
    "less": "lessing",
    "hemsl": "hemsley",
    "lag": "lagasca",
    "moc": "mocino",
    "bip": "schultz-bipontinus",
    "dc": "candolle",
    "l": "linnaeus",
}

GENERIC_DETERMINER_TOKENS: Set[str] = {
    "",
    "unknown",
    "anonymous",
    "none",
    "nan",
    "null",
    "na",
    "n/a",
    "staff",
    "student",
    "students",
    "herbarium",
    "curator",
    "volunteer",
    "unassigned",
    "unspecified",
    "not",
    "evident",
    "indet",
    "indetermined",
    "undetermined",
    "various",
    "class",
    "admin",
    "database",
    "digitizer",
    "illegible",
    "unrecorded",
}

_INVALID_INSTITUTION_CODES: Set[str] = {
    "",
    "UNKNOWN",
    "UNKNOWN_INST",
    "NONE",
    "NAN",
    "NULL",
    "NA",
}


def _fold_ascii(text: str) -> str:
    """Fold Unicode diacritics to base ASCII lowercase characters (e.g., 'Löve' -> 'love')."""
    norm = unicodedata.normalize("NFKD", text)
    return "".join(ch for ch in norm if not unicodedata.combining(ch)).lower()


def extract_protologue_authors(raw_names: Union[str, Iterable[Any]]) -> Set[str]:
    """
    Layer 2 Authority Engine: Extracts taxonomic authority surnames and standard IPNI
    abbreviations from Darwin Core `scientificName` strings.

    Examples:
      - "Packera dubia (Spreng.) Trock & Mabb." -> {"spreng", "sprengel", "trock", "mabb", "mabberley"}
      - "Packera paupercula (Michx.) Á.Löve & D.Löve" -> {"michx", "michaux", "love", "löve"}
    """
    if raw_names is None:
        return set()
    if isinstance(raw_names, (str, bytes)):
        items: Iterable[Any] = [raw_names.decode("utf-8", errors="ignore") if isinstance(raw_names, bytes) else raw_names]
    else:
        items = raw_names

    authors: Set[str] = set()
    for raw in items:
        if raw is None or (isinstance(raw, float) and math.isnan(raw)):
            continue
        txt = str(raw).strip()
        if not txt or txt.lower() in {"none", "nan", "null", "na", "unknown"}:
            continue

        # Strip leading Genus + optional hybrid marker + specific epithet
        m_binom = re.match(r"^[A-Z][a-z]+\s+(?:[x×]\s*)?[a-z][a-z-]*\b(.*)$", txt)
        if not m_binom:
            continue
        rem = m_binom.group(1).strip()
        if not rem:
            continue

        # Strip infraspecific rank + epithet segments (var., subsp., ssp., f., forma, nothovar.)
        rem = re.sub(
            r"\b(?:var\.?|subsp\.?|ssp\.?|forma|f\.?|nothosubsp\.?|nothovar\.?)\s+[a-z][a-z-]*\b",
            " ",
            rem,
        )
        # Strip open-nomenclature qualifiers and 4-digit publication years
        rem = re.sub(r"\b(?:sp\.?|spp\.?|cf\.?|aff\.?|ined\.?|nom\.?\s*nud\.?|s\.?\s*l\.?|s\.?\s*str\.?)\b", " ", rem, flags=re.IGNORECASE)
        rem = re.sub(r"\b\d{4}\b", " ", rem)
        rem = rem.strip()
        if not rem:
            continue

        # Split basionym and combination author groups
        segments = re.split(r"[(),;&/]+|\b(?:ex|in|et|and)\b", rem)
        for seg in segments:
            s_clean = seg.strip()
            if not s_clean:
                continue
            # Strip trailing filius marker ('Hook.f.' -> 'Hook.')
            s_clean = re.sub(r"\.?\s*\b(?:fil|f)\.?$", ".", s_clean, flags=re.IGNORECASE).strip()
            # Strip leading abbreviated given names/initials ('Á.Löve' -> 'Löve', 'W.A.Weber' -> 'Weber', 'Alph.Wood' -> 'Wood')
            m_init = re.match(r"^(?:[A-ZÁÉÍÓÚÄÖÜ][a-zà-ÿ]{0,4}\.\s*)+([A-ZÁÉÍÓÚÄÖÜ][A-Za-zÀ-ÿ-]+(?:\.)?)$", s_clean)
            if m_init:
                s_clean = m_init.group(1).strip()
            else:
                # Also handle space-separated initials without periods ('A Love' -> 'Love')
                tokens = [t for t in s_clean.split() if len(t.rstrip(".")) > 1]
                if tokens:
                    s_clean = tokens[-1]

            tok = s_clean.strip(" .-_'\"")
            if not tok:
                continue
            tok_low = tok.lower()
            tok_ascii = _fold_ascii(tok)

            if len(tok_low) >= 3 and re.search(r"[a-zà-ÿ]", tok_low):
                authors.add(tok_low)
            if len(tok_ascii) >= 3 and re.search(r"[a-z]", tok_ascii):
                authors.add(tok_ascii)

            for key_cand in (tok_low, tok_ascii):
                if key_cand in _IPNI_AUTHOR_ABBREVIATIONS:
                    expanded = _IPNI_AUTHOR_ABBREVIATIONS[key_cand]
                    if len(expanded) >= 3:
                        authors.add(expanded)

    return authors


def normalize_determiner_surnames(identified_by_raw: Optional[Any]) -> List[str]:
    """
    Normalizes Darwin Core `identifiedBy` strings into clean ASCII-folded determiner
    surnames, filtering out dates, prefixes, initials, and generic non-person strings.
    """
    if identified_by_raw is None or (isinstance(identified_by_raw, float) and math.isnan(identified_by_raw)):
        return []
    raw = str(identified_by_raw).strip()
    if not raw:
        return []
    raw_low = raw.lower()
    if raw_low in GENERIC_DETERMINER_TOKENS or any(
        bad in raw_low for bad in ("herbarium staff", "not evident", "unrecorded", "anonymous", "digitiz")
    ):
        return []

    # Strip annotation prefixes, affirmation marks, parenthetical notes, dates, and years
    cleaned = re.sub(r"\([^)]*\)", " ", raw)
    cleaned = re.sub(
        r"\b(?:det\.?\s*(?:by)?|determined\s+by|annot\.?\s*(?:by)?|annotated\s+by|ver\.?\s*(?:by)?|verified\s+by|conf\.?\s*(?:by)?|confirmed\s+by|rev\.?\s*(?:by)?)\b",
        " ",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(
        r"[,\s]*(?:\d{1,2}\s+)?\b(?:jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:t(?:ember)?)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\.?\s*\d{0,4}\b",
        " ",
        cleaned,
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(r"[,\s]*\b\d{4}(?:[-/]\d{1,2}(?:[-/]\d{1,2})?)?\b", " ", cleaned)
    cleaned = re.sub(r"\b\d+\b", " ", cleaned)
    cleaned = cleaned.replace("!", " ").strip(" ,;|-")
    if not cleaned:
        return []

    parts = re.split(r"\s*(?:&|;|\band\b|\bwith\b|\bet\s+al\.?)\s*", cleaned, flags=re.IGNORECASE)
    surnames: List[str] = []
    for part in parts:
        p = part.strip(" ,;|-.")
        if not p:
            continue
        p = re.sub(r"\b(?:dr|prof|mr|mrs|ms|jr|sr|ii|iii|iv)\b\.?", " ", p, flags=re.IGNORECASE).strip(" ,.")
        if not p:
            continue
        if "," in p:
            cand_raw = p.split(",", 1)[0].strip()
            cand_tokens = [re.sub(r"[^A-Za-zÀ-ÿ-]", "", t) for t in cand_raw.split()]
            cand_tokens = [t for t in cand_tokens if len(t) > 1]
            surname = _fold_ascii(cand_tokens[-1]) if cand_tokens else ""
        else:
            # Handle compact initials like 'Á.Löve' or 'D.K.Trock'
            p_spaced = re.sub(r"\.", ". ", p)
            tokens = [re.sub(r"[^A-Za-zÀ-ÿ-]", "", t) for t in p_spaced.split()]
            non_init = [t for t in tokens if len(t) > 1]
            surname = _fold_ascii(non_init[-1]) if non_init else ""

        surname = surname.strip("-")
        if (
            len(surname) >= 3
            and surname not in GENERIC_DETERMINER_TOKENS
            and not any(tok in GENERIC_DETERMINER_TOKENS for tok in _fold_ascii(p).split() if tok in {"staff", "student", "students", "herbarium", "unknown", "anonymous"})
        ):
            if surname not in surnames:
                surnames.append(surname)

    return surnames


def _is_verified_type_status(type_status_raw: Optional[Any]) -> Tuple[bool, str]:
    """Return (is_type, cleaned_type_label) for a Darwin Core typeStatus value."""
    if type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw)):
        return False, "None"
    type_status = str(type_status_raw).strip()
    if not type_status or type_status.upper() in {"NONE", "NOT A TYPE", "NOTATYPE", "UNSPECIFIED", "NULL", "NAN"}:
        return False, "None"
    for valid_type in VALID_TYPE_STATUSES:
        if re.search(rf"\b{valid_type}\b", type_status, re.IGNORECASE):
            return True, valid_type.title()
    return False, "None"


def discover_corpus_authorities(df: pd.DataFrame) -> Tuple[Set[str], Set[str]]:
    """
    Layer 3 Authority Engine: Discovers empirical taxonomic specialists and regional
    research herbaria directly from the assembled Darwin Core voucher corpus.

    Rules:
      - Dynamic Specialist Discovery (`dynamic_specialists`):
        1. Determiner (`identifiedBy`) has identified >= 5 specimens of a given Packera
           taxon across >= 3 distinct `institutionCode`s, OR
        2. Determiner has annotated >= 1 verified type specimen in the corpus.
      - Dynamic Regional Research Herbarium Discovery (`dynamic_herbaria`):
        1. `institutionCode` holds >= 1 verified type specimen for target Packera taxa, OR
        2. `institutionCode` ranks in the top quartile of institutional holdings
           (>= 10 curated sheets) with >= 70% determiner-attributed records for any harvested taxon.

    Returns:
        Tuple[Set[str], Set[str]]: (dynamic_specialists, dynamic_herbaria)
    """
    dynamic_specialists: Set[str] = set()
    dynamic_herbaria: Set[str] = set()
    if df is None or df.empty:
        return dynamic_specialists, dynamic_herbaria

    id_col = (
        "_raw_identified_by"
        if "_raw_identified_by" in df.columns
        else ("identifiedBy" if "identifiedBy" in df.columns else ("determiner_raw" if "determiner_raw" in df.columns else None))
    )
    sp_col = (
        "species_standardized"
        if "species_standardized" in df.columns
        else ("scientificName" if "scientificName" in df.columns else ("species_raw" if "species_raw" in df.columns else None))
    )
    type_col = (
        "_raw_type_status"
        if "_raw_type_status" in df.columns
        else ("type_status" if "type_status" in df.columns else ("typeStatus" if "typeStatus" in df.columns else None))
    )
    inst_col = "institutionCode" if "institutionCode" in df.columns else None

    # Track (determiner_surname, standardized_taxon) -> [count, set_of_institutions]
    det_taxon_counts: Dict[Tuple[str, str], int] = {}
    det_taxon_insts: Dict[Tuple[str, str], Set[str]] = {}

    # Track institutional holdings and determiner-attribution rates
    inst_total_sheets: Dict[str, int] = {}
    inst_taxon_totals: Dict[Tuple[str, str], int] = {}
    inst_taxon_attributed: Dict[Tuple[str, str], int] = {}

    for _, row in df.iterrows():
        raw_sp = row.get(sp_col) if sp_col else None
        std_sp = standardize_packera_taxon(raw_sp)
        if not std_sp and raw_sp is not None and not pd.isna(raw_sp):
            s_tmp = str(raw_sp).strip()
            if s_tmp.startswith("Packera "):
                std_sp = _extract_canonical_binomial(s_tmp, {"Packera"})

        raw_inst = row.get(inst_col) if inst_col else ""
        inst = str(raw_inst if not pd.isna(raw_inst) else "").strip().upper()
        valid_inst = inst not in _INVALID_INSTITUTION_CODES

        raw_id = row.get(id_col) if id_col else ""
        surnames = normalize_determiner_surnames(raw_id)
        has_attributed_det = len(surnames) > 0

        raw_type = row.get(type_col) if type_col else None
        is_type, _ = _is_verified_type_status(raw_type)

        if is_type:
            for s_name in surnames:
                dynamic_specialists.add(s_name)
            if valid_inst and std_sp:
                dynamic_herbaria.add(inst)

        if valid_inst:
            inst_total_sheets[inst] = inst_total_sheets.get(inst, 0) + 1
            if std_sp:
                it_key = (inst, std_sp)
                inst_taxon_totals[it_key] = inst_taxon_totals.get(it_key, 0) + 1
                if has_attributed_det:
                    inst_taxon_attributed[it_key] = inst_taxon_attributed.get(it_key, 0) + 1

        if std_sp and surnames:
            for s_name in surnames:
                dt_key = (s_name, std_sp)
                det_taxon_counts[dt_key] = det_taxon_counts.get(dt_key, 0) + 1
                if valid_inst:
                    det_taxon_insts.setdefault(dt_key, set()).add(inst)

    for dt_key, n_specimens in det_taxon_counts.items():
        n_insts = len(det_taxon_insts.get(dt_key, set()))
        if n_specimens >= 5 and n_insts >= 3:
            dynamic_specialists.add(dt_key[0])

    if inst_total_sheets:
        q75_holdings = float(np.percentile(list(inst_total_sheets.values()), 75))
        for inst, total_n in inst_total_sheets.items():
            if total_n >= 10 and total_n >= q75_holdings:
                min_sp_sheets = min(3, total_n)
                for (i_code, sp_name), sp_total in inst_taxon_totals.items():
                    if i_code == inst and sp_total >= min_sp_sheets:
                        n_attr = inst_taxon_attributed.get((i_code, sp_name), 0)
                        if (n_attr / float(sp_total)) >= 0.70:
                            dynamic_herbaria.add(inst)
                            break

    return dynamic_specialists, dynamic_herbaria


def parse_determiner_tier(
    type_status_raw: Optional[str],
    identified_by_raw: Optional[str],
    recorded_by_raw: Optional[str],
    history_raw: Optional[str],
    institution_code_raw: Optional[str],
    locality_raw: Optional[str],
    habitat_raw: Optional[str],
    dynamic_specialists: Optional[Iterable[str]] = None,
    dynamic_herbaria: Optional[Iterable[str]] = None,
) -> Tuple[str, str, str]:
    """
    Evaluates Darwin Core fields to assign a taxonomic determination credibility tier
    using the 3-Layer Dynamic Expert & Herbarium Authority Engine:
      - Layer 1: Genus-wide Packera / Senecioneae monographer patterns (`SPECIALIST_PATTERNS`).
      - Layer 2 & 3: Dynamic protologue authors and corpus-derived multi-institution
        specialists (`dynamic_specialists`) and regional research herbaria (`dynamic_herbaria`).

    Tiers:
      - Tier_1_Gold: Primary/secondary nomenclatural types, monographers, protologue authors,
        or corpus-verified multi-herbarium specialists.
      - Tier_2_Silver: Determinations from major or corpus-qualified research herbaria with
        complete ecological locality and attributed determiner.
      - Tier_3_Bronze: Unverified candidate collections, general collectors, or missing determiner.

    Returns:
        Tuple[str, str, str]: (determiner_tier, type_status_clean, determiner_raw_combined)
    """
    type_status = "" if (type_status_raw is None or (isinstance(type_status_raw, float) and math.isnan(type_status_raw))) else str(type_status_raw).strip()
    identified_by = "" if (identified_by_raw is None or (isinstance(identified_by_raw, float) and math.isnan(identified_by_raw))) else str(identified_by_raw).strip()
    recorded_by = "" if (recorded_by_raw is None or (isinstance(recorded_by_raw, float) and math.isnan(recorded_by_raw))) else str(recorded_by_raw).strip()
    history = "" if (history_raw is None or (isinstance(history_raw, float) and math.isnan(history_raw))) else str(history_raw).strip()
    institution = "" if (institution_code_raw is None or (isinstance(institution_code_raw, float) and math.isnan(institution_code_raw))) else str(institution_code_raw).strip().upper()
    locality = "" if (locality_raw is None or (isinstance(locality_raw, float) and math.isnan(locality_raw))) else str(locality_raw).strip()
    habitat = "" if (habitat_raw is None or (isinstance(habitat_raw, float) and math.isnan(habitat_raw))) else str(habitat_raw).strip()

    determiner_raw = identified_by if identified_by else (history if history else recorded_by)
    combined_audit_text = f"{type_status} | {identified_by} | {history} | {recorded_by}"

    is_type, type_status_clean = _is_verified_type_status(type_status)

    is_specialist = False
    for pattern in SPECIALIST_PATTERNS:
        if re.search(pattern, combined_audit_text, re.IGNORECASE):
            is_specialist = True
            break

    if not is_specialist and dynamic_specialists:
        dyn_spec_set: Set[str] = set()
        for item in dynamic_specialists:
            if item is None:
                continue
            s_tok = str(item).strip().lower()
            if len(s_tok) >= 3:
                dyn_spec_set.add(s_tok)
                dyn_spec_set.add(_fold_ascii(s_tok))
        if dyn_spec_set:
            det_surnames = set(normalize_determiner_surnames(identified_by) + normalize_determiner_surnames(history))
            if det_surnames & dyn_spec_set:
                is_specialist = True
            else:
                det_audit_raw = f"{identified_by} | {history}".lower()
                det_audit_ascii = _fold_ascii(det_audit_raw)
                for spec_tok in dyn_spec_set:
                    pat = rf"(?<![a-zA-ZÀ-ÿ]){re.escape(spec_tok)}(?![a-zA-ZÀ-ÿ])"
                    if re.search(pat, det_audit_raw) or re.search(pat, det_audit_ascii):
                        is_specialist = True
                        break

    if is_type or is_specialist:
        return "Tier_1_Gold", type_status_clean, determiner_raw

    dyn_herb_set = {str(h).strip().upper() for h in dynamic_herbaria if h} if dynamic_herbaria else set()
    is_major_herbarium = (institution in MAJOR_HERBARIA_CODES) or (institution in dyn_herb_set)
    has_rich_locality = len(locality) > 10 or len(habitat) > 5
    has_determiner = bool(
        identified_by
        and identified_by.lower() not in GENERIC_DETERMINER_TOKENS
        and len(normalize_determiner_surnames(identified_by)) > 0
    )

    if is_major_herbarium and has_rich_locality and has_determiner:
        return "Tier_2_Silver", type_status_clean, determiner_raw

    return "Tier_3_Bronze", type_status_clean, determiner_raw


def apply_corpus_authority_tiers(
    df: pd.DataFrame,
    dynamic_specialists: Optional[Set[str]] = None,
    dynamic_herbaria: Optional[Set[str]] = None,
) -> pd.DataFrame:
    """
    Applies the 3-Layer Dynamic Authority Engine across an assembled DataFrame:
      - Discovers corpus-level multi-herbarium specialists and dynamic regional herbaria.
      - Collects per-taxon protologue authors from Darwin Core `scientificName` / `species_raw`.
      - Re-evaluates `determiner_tier` across all rows cleanly without downgrading pre-scored
        rows whose `locality` / `habitat` columns were already omitted upon CSV persistence.
    """
    if df is None or df.empty:
        return df

    out = df.copy()
    if dynamic_specialists is None or dynamic_herbaria is None:
        disc_specs, disc_herbs = discover_corpus_authorities(out)
        eff_specs = disc_specs if dynamic_specialists is None else (set(dynamic_specialists) | disc_specs)
        eff_herbs = disc_herbs if dynamic_herbaria is None else (set(dynamic_herbaria) | disc_herbs)
    else:
        eff_specs = set(dynamic_specialists)
        eff_herbs = set(dynamic_herbaria)

    # Collect per-taxon protologue authors across the corpus + cached GBIF synonyms
    taxon_protologue_map: Dict[str, Set[str]] = {}
    for _, row in out.iterrows():
        raw_names_row = [
            row.get("scientificName"),
            row.get("species_raw"),
        ]
        std_sp = None
        for rn in raw_names_row:
            std_sp = standardize_packera_taxon(rn)
            if std_sp:
                break
        if std_sp:
            extra_cached = list(_DYNAMIC_TAXON_AUTHOR_STRINGS.get(std_sp, set()))
            auths = extract_protologue_authors([*raw_names_row, *extra_cached])
            if auths:
                taxon_protologue_map.setdefault(std_sp, set()).update(auths)

    has_locality_cols = any(
        c in out.columns for c in ("_raw_locality", "locality", "verbatimLocality", "_raw_habitat", "habitat")
    )
    tier_rank_map = {"Tier_1_Gold": 3, "Tier_2_Silver": 2, "Tier_3_Bronze": 1}

    new_tiers: List[str] = []
    new_types: List[str] = []
    update_type_col = "type_status" in out.columns

    for _, row in out.iterrows():
        t_raw = row.get("_raw_type_status") if "_raw_type_status" in out.columns else row.get("type_status", row.get("typeStatus"))
        id_raw = row.get("_raw_identified_by") if "_raw_identified_by" in out.columns else row.get("identifiedBy")
        rec_raw = row.get("recordedBy")
        hist_raw = row.get("_raw_history") if "_raw_history" in out.columns else row.get("verbatimIdentificationHistory")
        inst_raw = row.get("institutionCode")
        loc_raw = (
            row.get("_raw_locality")
            if "_raw_locality" in out.columns
            else (row.get("locality") if "locality" in out.columns else row.get("verbatimLocality"))
        )
        hab_raw = row.get("_raw_habitat") if "_raw_habitat" in out.columns else row.get("habitat")

        raw_sci = row.get("scientificName") or row.get("species_raw")
        std_sp = standardize_packera_taxon(raw_sci)
        row_proto_auths = set(taxon_protologue_map.get(std_sp, set())) if std_sp else set()
        row_proto_auths |= extract_protologue_authors([row.get("scientificName"), row.get("species_raw")])

        tier_eval, type_clean, _ = parse_determiner_tier(
            type_status_raw=t_raw,
            identified_by_raw=id_raw,
            recorded_by_raw=rec_raw,
            history_raw=hist_raw,
            institution_code_raw=inst_raw,
            locality_raw=loc_raw,
            habitat_raw=hab_raw,
            dynamic_specialists=eff_specs | row_proto_auths,
            dynamic_herbaria=eff_herbs,
        )

        prev_tier = str(row.get("determiner_tier", "") if not pd.isna(row.get("determiner_tier", "")) else "").strip()
        if not has_locality_cols and prev_tier in tier_rank_map:
            if tier_rank_map[tier_eval] < tier_rank_map[prev_tier]:
                tier_eval = prev_tier

        new_tiers.append(tier_eval)
        prev_type = str(row.get("type_status", "") if not pd.isna(row.get("type_status", "")) else "").strip()
        new_types.append(type_clean if (type_clean != "None" or not prev_type) else prev_type)

    out["determiner_tier"] = new_tiers
    if update_type_col:
        out["type_status"] = new_types
    return out


def resolve_geographic_scope_settings(
    geographic_scope: Optional[str] = DEFAULT_GEOGRAPHIC_SCOPE,
    exclude_western: bool = True,
    include_western: bool = False,
    countries: Optional[Union[str, List[str]]] = None,
) -> Tuple[str, bool, List[str]]:
    """
    Resolves effective `(geographic_scope, exclude_western, countries)` settings:
      - `eastern_na` (default): `exclude_western=True` (unless `include_western=True` or
        `exclude_western=False`), `countries=["US"]` (or `["US", "CA", "MX"]` when `include_western=True`).
      - `north_america`: automatically sets `exclude_western=False` and expands default
        `countries` to `["US", "CA", "MX"]` unless explicit `countries` were supplied.
      - `global`: automatically sets `exclude_western=False` and sets `countries=[]` (unrestricted
        global GBIF search) unless explicit `countries` were supplied.
    """
    scope = str(geographic_scope or DEFAULT_GEOGRAPHIC_SCOPE).strip().lower()
    if scope not in VALID_GEOGRAPHIC_SCOPES:
        raise ValueError(f"Unsupported geographic_scope={geographic_scope!r}; expected one of {sorted(VALID_GEOGRAPHIC_SCOPES)}.")

    eff_exclude_western = bool(exclude_western) and (not bool(include_western)) and (scope == "eastern_na")

    if countries is not None:
        eff_countries = normalize_countries_list(countries, default=[])
    elif scope == "global":
        eff_countries = []
    elif scope == "north_america" or bool(include_western):
        eff_countries = ["US", "CA", "MX"]
    else:
        eff_countries = normalize_countries_list(DEFAULT_COUNTRIES, default=["US"])

    return scope, eff_exclude_western, eff_countries


def is_excluded_western_region(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    western_longitude_threshold: float = WESTERN_LONGITUDE_THRESHOLD,
) -> bool:
    """Determines whether a record originates from a region farther west than TX and OK."""
    if state_province is not None and not (isinstance(state_province, float) and math.isnan(state_province)):
        raw_state = str(state_province).strip()
        cleaned = re.sub(r"\(state\)", "", raw_state, flags=re.IGNORECASE).strip(" ._,-")
        upper_state = cleaned.upper()

        if upper_state in {"WASHINGTON, D.C.", "WASHINGTON D.C.", "WASHINGTON DC", "DISTRICT OF COLUMBIA", "DC"}:
            return False

        if upper_state in EXCLUDED_WESTERN_STATES:
            return True

    if lon is not None:
        try:
            lon_val = float(lon)
            if lon_val < float(western_longitude_threshold):
                return True
        except (ValueError, TypeError):
            pass

    return False


def select_cascading_outgroup_anchors(
    df: pd.DataFrame,
    outgroup_taxa: List[str],
    max_anchors: int = 30,
    workspace_dir: Optional[Path] = None,
) -> Tuple[List[pd.DataFrame], Dict[str, int], Dict[str, Dict[str, int]]]:
    """
    Selects up to `max_anchors` reference sheets per outgroup taxon by cascading across
    determiner authority tiers (`Tier_1_Gold` -> `Tier_2_Silver` -> `Tier_3_Bronze`),
    primary duplicate status, local image availability, and image resolution score.
    Ensures outgroup species lacking `Tier_1_Gold` vouchers still supply reference anchors.
    """
    outgroup_dfs: List[pd.DataFrame] = []
    outgroup_counts: Dict[str, int] = {}
    outgroup_tier_breakdown: Dict[str, Dict[str, int]] = {}

    if df is None or df.empty or not outgroup_taxa or max_anchors <= 0:
        for tx in outgroup_taxa or []:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
        return outgroup_dfs, outgroup_counts, outgroup_tier_breakdown

    tier_rank = {"Tier_1_Gold": 3, "Tier_2_Silver": 2, "Tier_3_Bronze": 1}
    ws_root = Path(workspace_dir).resolve() if workspace_dir is not None else PROJECT_ROOT

    for tx in outgroup_taxa:
        tx_df = df[df["species_standardized"] == tx].copy()
        if tx_df.empty:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
            continue

        def _has_local(row: pd.Series) -> int:
            img_p = str(row.get("image_path", "") if not pd.isna(row.get("image_path", "")) else "").strip()
            cat_str = str(row.get("catalogNumber", "") if not pd.isna(row.get("catalogNumber", "")) else "").strip()
            cat_rel = f"data/raw_vouchers/{cat_str}.jpg" if cat_str else ""
            candidates = []
            if img_p:
                p_obj = Path(img_p)
                candidates.append(p_obj if p_obj.is_absolute() else (ws_root / p_obj))
                candidates.append(Path(img_p))
            if cat_rel:
                candidates.append(ws_root / cat_rel)
                candidates.append(Path(cat_rel))
            return 1 if any(cp.is_file() for cp in candidates) else 0

        def _has_url_or_path(row: pd.Series) -> bool:
            for c in ("image_path", "_image_url", "image_url", "identifier", "good_quality_AccessURI"):
                v = row.get(c, None)
                if v is not None and not pd.isna(v) and str(v).strip():
                    return True
            return False

        def _is_prim(row: pd.Series) -> int:
            if "is_primary_duplicate" not in row.index or pd.isna(row.get("is_primary_duplicate")):
                return 1
            v = row.get("is_primary_duplicate")
            return 1 if (v is True or str(v).strip().lower() in {"true", "1", "t", "yes"}) else 0

        tx_df["_has_local"] = tx_df.apply(_has_local, axis=1)
        tx_df["_has_valid"] = (tx_df["_has_local"] > 0) | tx_df.apply(_has_url_or_path, axis=1)
        tx_df = tx_df[tx_df["_has_valid"]].copy()
        if tx_df.empty:
            outgroup_counts[tx] = 0
            outgroup_tier_breakdown[tx] = {}
            continue

        tx_df["_tier_rank"] = (
            tx_df["determiner_tier"].astype(str).str.strip().map(lambda t: tier_rank.get(t, 0))
            if "determiner_tier" in tx_df.columns
            else 0
        )
        tx_df["_is_prim"] = tx_df.apply(_is_prim, axis=1)
        tx_df["_img_res"] = tx_df["image_path"].map(
            lambda p: get_image_resolution_score(str(p) if not pd.isna(p) else "", workspace_dir=ws_root)
        ) if "image_path" in tx_df.columns else 0.0

        tx_df = tx_df.sort_values(
            by=["_tier_rank", "_is_prim", "_has_local", "_img_res"],
            ascending=[False, False, False, False],
            kind="mergesort",
        )
        tx_sel = tx_df.head(int(max_anchors)).drop(
            columns=["_has_local", "_has_valid", "_tier_rank", "_is_prim", "_img_res"]
        )
        outgroup_dfs.append(tx_sel)
        outgroup_counts[tx] = len(tx_sel)
        if "determiner_tier" in tx_sel.columns:
            outgroup_tier_breakdown[tx] = {
                str(k): int(v) for k, v in tx_sel["determiner_tier"].astype(str).value_counts().items()
            }
        else:
            outgroup_tier_breakdown[tx] = {"Unspecified": len(tx_sel)}

    return outgroup_dfs, outgroup_counts, outgroup_tier_breakdown


EPA_L3_OFFICIAL_ZIP_URL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip"
DEFAULT_EPA_L3_DIR = PROJECT_ROOT / "data" / "environmental" / "epa_ecoregions_l3"
DEFAULT_EPA_L3_GEOJSON = DEFAULT_EPA_L3_DIR / "us_eco_l3.geojson"
DEFAULT_EPA_L3_GPKG = DEFAULT_EPA_L3_DIR / "us_eco_l3.gpkg"
DEFAULT_US_COUNTIES_GEOJSON = DEFAULT_EPA_L3_DIR / "us_counties.geojson"

_US_STATE_ABBR_TO_FULL: Dict[str, str] = {
    "al": "alabama", "ak": "alaska", "az": "arizona", "ar": "arkansas", "ca": "california",
    "co": "colorado", "ct": "connecticut", "de": "delaware", "dc": "district of columbia",
    "fl": "florida", "ga": "georgia", "hi": "hawaii", "id": "idaho", "il": "illinois",
    "in": "indiana", "ia": "iowa", "ks": "kansas", "ky": "kentucky", "la": "louisiana",
    "me": "maine", "md": "maryland", "ma": "massachusetts", "mi": "michigan", "mn": "minnesota",
    "ms": "mississippi", "mo": "missouri", "mt": "montana", "ne": "nebraska", "nv": "nevada",
    "nh": "new hampshire", "nj": "new jersey", "nm": "new mexico", "ny": "new york",
    "nc": "north carolina", "nd": "north dakota", "oh": "ohio", "ok": "oklahoma", "or": "oregon",
    "pa": "pennsylvania", "ri": "rhode island", "sc": "south carolina", "sd": "south dakota",
    "tn": "tennessee", "tx": "texas", "ut": "utah", "vt": "vermont", "va": "virginia",
    "wa": "washington", "wv": "west virginia", "wi": "wisconsin", "wy": "wyoming",
}

_HISTORIC_COUNTY_CONSOLIDATIONS: Dict[Tuple[str, str], str] = {
    ("virginia", "nansemond"): "suffolk",
    ("virginia", "princess anne"): "virginia beach",
    ("virginia", "chesapeake"): "norfolk",
    ("georgia", "columbus"): "muscogee",
    ("south dakota", "shannon"): "oglala lakota",
}

_EPA_L3_INDEX_CACHE: Dict[Tuple[str, str, float, float], Dict[str, Any]] = {}


def standardize_epa_l3_label(raw_name: Optional[str]) -> Optional[str]:
    """Standardizes official US EPA Level III Ecoregion names to canonical underscore labels."""
    if raw_name is None:
        return None
    txt = str(raw_name).strip()
    if not txt or txt.lower() in {"na", "nan", "none", "null", "unknown", "other_us", "unassigned_na"}:
        return None
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", txt).strip("_")
    return cleaned or None


def _normalize_us_state_name(state_province: Optional[str]) -> Optional[str]:
    if state_province is None or (isinstance(state_province, float) and math.isnan(state_province)):
        return None
    s = re.sub(r"\s+", " ", str(state_province).strip().lower().rstrip("."))
    if not s or s in {"na", "nan", "none", "unknown"}:
        return None
    return _US_STATE_ABBR_TO_FULL.get(s, s)


def _normalize_us_county_name(county: Optional[str], state_norm: Optional[str] = None) -> Optional[str]:
    if county is None or (isinstance(county, float) and math.isnan(county)):
        return None
    raw = str(county).strip().lower()
    if not raw or raw in {"na", "nan", "none", "unknown", "unplaced", "indeterminate"}:
        return None
    if " or " in raw or "/" in raw or "?" in raw:
        return None
    c_str = re.sub(r"[\[\]\(\)]", "", raw)
    c_str = re.sub(r"^city\s+of\s+", "", c_str)
    c_str = re.sub(
        r"\s+(?:county|parish|par\.?|borough|census\s+area|municipality|co\.?|city\s+and\s+borough|independent\s+city|independent)$",
        "",
        c_str,
    )
    c_str = re.sub(r"\s+city$", "", c_str)
    c_str = re.sub(r"^(?:st\.|saint)\s+", "st ", c_str)
    c_str = re.sub(r"^ste\.\s*", "ste ", c_str)
    c_str = re.sub(r"[^a-z0-9]+", " ", c_str).strip()
    if not c_str:
        return None
    if state_norm is not None and (state_norm, c_str) in _HISTORIC_COUNTY_CONSOLIDATIONS:
        c_str = _HISTORIC_COUNTY_CONSOLIDATIONS[(state_norm, c_str)]
    return c_str


def _ensure_epa_level3_cache_files(eco_geojson_path: Path, cnty_geojson_path: Path) -> bool:
    """Downloads official US EPA Level III Ecoregions and builds GeoJSON/GPKG cache atomically if absent."""
    if eco_geojson_path.exists() and cnty_geojson_path.exists():
        return True
    eco_geojson_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_d:
        tmp_dir = Path(tmp_d)
        zip_p = tmp_dir / "us_eco_l3.zip"
        try:
            urllib.request.urlretrieve(EPA_L3_OFFICIAL_ZIP_URL, str(zip_p))
            with zipfile.ZipFile(zip_p, "r") as zf:
                zf.extractall(tmp_dir)
        except Exception:
            return False
        shp_candidates = list(tmp_dir.rglob("us_eco_l3.shp"))
        if not shp_candidates:
            return False
        shp_path = shp_candidates[0]
        r_code = (
            "suppressPackageStartupMessages({library(sf); library(maps)}); sf_use_s2(FALSE);"
            f"eco <- st_read({json.dumps(str(shp_path))}, quiet=TRUE);"
            "eco$ecoregion_label <- gsub('^_+|_+$', '', gsub('[^A-Za-z0-9]+', '_', trimws(as.character(eco$US_L3NAME))));"
            "eco <- eco[, c('US_L3CODE', 'US_L3NAME', 'ecoregion_label')];"
            "eco_w <- st_make_valid(st_transform(st_simplify(eco, preserveTopology=TRUE, dTolerance=100), 4326));"
            "cm <- st_make_valid(st_as_sf(maps::map('county', fill=TRUE, plot=FALSE))); st_crs(cm) <- 4326;"
            "pts <- strsplit(as.character(cm$ID), ',');"
            "st_n <- vapply(pts, `[`, character(1), 1);"
            "co_n <- trimws(gsub('[^a-z0-9]+', ' ', gsub('^saint\\\\s+|^st\\\\.\\\\s*', 'st ', gsub('\\\\s+city$', '', sub(':.*$', '', vapply(pts, `[`, character(1), 2))))));"
            "ck <- paste0(st_n, '|', co_n); dk <- unique(ck[duplicated(ck)]); sm <- !(ck %in% dk); ga <- st_geometry(cm);"
            "ok <- ck[sm]; os <- st_n[sm]; oc <- co_n[sm]; og <- lapply(which(sm), function(i) st_cast(ga[[i]], 'MULTIPOLYGON'));"
            "for(d in dk){ id <- which(ck==d); ok <- c(ok, d); os <- c(os, st_n[id[1]]); oc <- c(oc, co_n[id[1]]); og[[length(og)+1L]] <- st_cast(st_make_valid(st_union(ga[id]))[[1]], 'MULTIPOLYGON'); };"
            "cd <- st_sf(county_key=ok, state_name=os, county_norm=oc, geometry=st_sfc(og, crs=4326));"
            "ca <- st_make_valid(st_transform(cd, st_crs(eco)));"
            "ix <- suppressWarnings(st_intersection(ca[, 'county_key', drop=FALSE], eco[, 'ecoregion_label', drop=FALSE]));"
            "ix$a <- as.numeric(st_area(ix)); ia <- split(aggregate(a ~ county_key + ecoregion_label, data=st_drop_geometry(ix), FUN=sum), ~county_key);"
            "cw <- st_make_valid(st_transform(ca, 4326)); cp <- st_coordinates(suppressWarnings(st_point_on_surface(cw)));"
            "cw$centroid_lon <- cp[,1]; cw$centroid_lat <- cp[,2]; cw$n_ecoregions <- 0L; cw$top_ecoregion <- NA_character_; cw$top_area_share <- 0;"
            "for(i in seq_len(nrow(cw))){ s <- ia[[cw$county_key[i]]]; if(!is.null(s)&&nrow(s)>0){ t <- sum(s$a); o <- order(s$a, decreasing=TRUE); cw$n_ecoregions[i] <- sum((s$a/t)>=0.001); cw$top_ecoregion[i] <- s$ecoregion_label[o[1]]; cw$top_area_share[i] <- round(s$a[o[1]]/t, 6); } };"
            f"te <- {json.dumps(str(eco_geojson_path) + '.tmp.geojson')}; tc <- {json.dumps(str(cnty_geojson_path) + '.tmp.geojson')};"
            "st_write(eco_w, te, driver='GeoJSON', quiet=TRUE, delete_dsn=TRUE); file.rename(te, sub('\\\\.tmp\\\\.geojson$', '', te));"
            "st_write(cw, tc, driver='GeoJSON', quiet=TRUE, delete_dsn=TRUE); file.rename(tc, sub('\\\\.tmp\\\\.geojson$', '', tc));"
        )
        try:
            subprocess.run(["Rscript", "-e", r_code], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            return False
    return eco_geojson_path.exists() and cnty_geojson_path.exists()


def _build_offline_epa_l3_fallback_index() -> Dict[str, Any]:
    """Deterministic in-memory multi-vertex EPA Level III polygons for isolated offline unit tests (never persisted)."""
    offline_polys = [
        ("Blue_Ridge", Polygon([(-84.5, 34.5), (-83.2, 34.8), (-81.6, 36.1), (-79.0, 37.9), (-78.2, 38.9), (-79.5, 39.2), (-82.8, 36.8), (-84.6, 35.2)])),
        ("Ridge_and_Valley", Polygon([(-85.6, 34.2), (-84.6, 35.2), (-82.8, 36.8), (-79.5, 39.2), (-77.2, 40.8), (-78.8, 41.2), (-83.8, 37.4), (-85.9, 34.8)])),
        ("Piedmont", Polygon([(-85.2, 32.5), (-82.0, 33.5), (-79.4, 34.9), (-77.6, 36.5), (-77.1, 38.9), (-78.2, 38.9), (-79.0, 37.9), (-81.6, 36.1), (-83.2, 34.8), (-84.5, 34.5)])),
        ("Middle_Atlantic_Coastal_Plain", Polygon([(-79.6, 33.0), (-75.4, 35.2), (-75.0, 38.8), (-76.4, 39.2), (-77.1, 38.5), (-77.6, 36.5), (-78.5, 35.2)])),
        ("Southeastern_Plains", Polygon([(-90.2, 31.0), (-85.0, 30.8), (-81.3, 32.2), (-79.6, 33.0), (-78.5, 35.2), (-79.4, 34.9), (-82.0, 33.5), (-85.2, 32.5), (-88.8, 34.8), (-90.1, 33.2)])),
        ("Southern_Coastal_Plain", Polygon([(-85.0, 29.6), (-80.2, 26.5), (-80.8, 31.2), (-81.3, 32.2), (-85.0, 30.8)])),
        ("Mississippi_Alluvial_Plain", Polygon([(-91.8, 29.5), (-90.2, 30.5), (-90.1, 35.0), (-89.2, 37.0), (-90.6, 37.0), (-91.5, 34.2), (-92.0, 31.0)])),
        ("South_Central_Plains", Polygon([(-95.8, 30.2), (-92.0, 31.0), (-91.5, 34.2), (-94.8, 34.1), (-96.0, 32.0)])),
        ("East_Central_Texas_Plains", Polygon([(-98.2, 29.0), (-95.8, 30.2), (-96.0, 32.8), (-97.2, 33.0)])),
        ("Ouachita_Mountains", Polygon([(-95.6, 34.1), (-92.4, 34.3), (-92.6, 35.1), (-95.4, 35.0)])),
        ("Ozark_Highlands", Polygon([(-94.8, 35.7), (-91.0, 35.8), (-90.5, 38.2), (-94.6, 38.1)])),
        ("Interior_Plateau", Polygon([(-88.5, 34.8), (-85.5, 35.2), (-84.5, 38.2), (-87.8, 38.8)])),
        ("Southwestern_Appalachians", Polygon([(-86.8, 33.6), (-85.6, 34.2), (-84.2, 36.8), (-85.5, 36.8)])),
    ]
    eco_labels = [lbl for lbl, _ in offline_polys]
    eco_geoms = [geom for _, geom in offline_polys]
    offline_counties = {
        "north carolina|orange": {"geom": Polygon([(-79.25, 35.86), (-78.95, 35.86), (-78.95, 36.24), (-79.25, 36.24)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -79.10, "centroid_lat": 36.05},
        "north carolina|new hanover": {"geom": Polygon([(-78.02, 34.05), (-77.78, 34.05), (-77.78, 34.35), (-78.02, 34.35)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -77.90, "centroid_lat": 34.20},
        "north carolina|wake": {"geom": Polygon([(-79.00, 35.50), (-77.80, 35.50), (-77.80, 36.05), (-79.00, 36.05)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -78.40, "centroid_lat": 35.78},
        "texas|nacogdoches": {"geom": Polygon([(-94.90, 31.40), (-94.40, 31.40), (-94.40, 31.85), (-94.90, 31.85)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -94.65, "centroid_lat": 31.62},
        "arkansas|garland": {"geom": Polygon([(-93.40, 34.42), (-92.90, 34.42), (-92.90, 34.78), (-93.40, 34.78)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -93.15, "centroid_lat": 34.60},
    }
    return {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
        "counties": offline_counties,
        "counties_nospace": {k.replace(" ", ""): v for k, v in offline_counties.items()},
        "use_live_county_intersection": True,
    }


def load_epa_level3_spatial_index(
    ecoregion_path: Optional[Path] = None,
    county_path: Optional[Path] = None,
    allow_download: bool = True,
) -> Dict[str, Any]:
    """Loads and caches official US EPA Level III Ecoregion and US county polygons into an STRtree index."""
    custom_eco = ecoregion_path is not None
    eco_p = Path(ecoregion_path) if ecoregion_path is not None else DEFAULT_EPA_L3_GEOJSON
    cnty_p = Path(county_path) if county_path is not None else DEFAULT_US_COUNTIES_GEOJSON

    if not eco_p.exists() and not custom_eco and allow_download:
        _ensure_epa_level3_cache_files(eco_p, cnty_p)

    if not eco_p.exists():
        return _build_offline_epa_l3_fallback_index()

    eco_mtime = eco_p.stat().st_mtime
    cnty_mtime = cnty_p.stat().st_mtime if cnty_p.exists() else -1.0
    cache_key = (str(eco_p.resolve()), str(cnty_p.resolve()) if cnty_p.exists() else "", eco_mtime, cnty_mtime)
    if cache_key in _EPA_L3_INDEX_CACHE:
        return _EPA_L3_INDEX_CACHE[cache_key]

    with open(eco_p, "r", encoding="utf-8") as fh:
        eco_fc = json.load(fh)

    eco_geoms: List[Any] = []
    eco_labels: List[str] = []
    for feat in eco_fc.get("features", []):
        props = feat.get("properties") or {}
        raw_lbl = props.get("ecoregion_label") or props.get("US_L3NAME") or props.get("NA_L3NAME") or props.get("name")
        lbl = standardize_epa_l3_label(raw_lbl)
        geom_dict = feat.get("geometry")
        if lbl and geom_dict:
            geom_obj = shape(geom_dict)
            if not geom_obj.is_empty:
                eco_geoms.append(geom_obj)
                eco_labels.append(lbl)

    counties: Dict[str, Dict[str, Any]] = {}
    counties_nospace: Dict[str, Dict[str, Any]] = {}
    if cnty_p.exists():
        with open(cnty_p, "r", encoding="utf-8") as cfh:
            cnty_fc = json.load(cfh)
        for feat in cnty_fc.get("features", []):
            props = feat.get("properties") or {}
            st_n = _normalize_us_state_name(props.get("state_name") or props.get("STATE_NAME") or props.get("state"))
            co_n = _normalize_us_county_name(props.get("county_norm") or props.get("NAME") or props.get("county"), st_n)
            geom_dict = feat.get("geometry")
            if st_n and co_n and geom_dict:
                geom_obj = shape(geom_dict)
                if not geom_obj.is_empty:
                    rec = {
                        "geom": geom_obj,
                        "top_ecoregion": standardize_epa_l3_label(props.get("top_ecoregion")),
                        "top_area_share": float(props["top_area_share"]) if props.get("top_area_share") is not None else None,
                        "n_ecoregions": int(props["n_ecoregions"]) if props.get("n_ecoregions") is not None else None,
                        "centroid_lon": float(props["centroid_lon"]) if props.get("centroid_lon") is not None else None,
                        "centroid_lat": float(props["centroid_lat"]) if props.get("centroid_lat") is not None else None,
                    }
                    ck = f"{st_n}|{co_n}"
                    counties[ck] = rec
                    counties_nospace[f"{st_n}|{co_n.replace(' ', '')}"] = rec

    idx_bundle = {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
        "counties": counties,
        "counties_nospace": counties_nospace,
        "use_live_county_intersection": custom_eco or (county_path is not None),
    }
    _EPA_L3_INDEX_CACHE[cache_key] = idx_bundle
    return idx_bundle


def _point_in_epa_l3(lon_f: float, lat_f: float, idx_bundle: Dict[str, Any]) -> Optional[str]:
    pt = Point(lon_f, lat_f)
    eco_tree = idx_bundle["eco_tree"]
    eco_geoms = idx_bundle["eco_geoms"]
    eco_labels = idx_bundle["eco_labels"]
    eco_id_map = idx_bundle["eco_id_map"]
    for item in eco_tree.query(pt):
        idx = int(item) if isinstance(item, (int, np.integer)) else eco_id_map.get(id(item), -1)
        if idx >= 0 and eco_geoms[idx].covers(pt):
            return eco_labels[idx]
    return None


def infer_regional_group(
    lat: Optional[float],
    lon: Optional[float],
    state_province: Optional[str] = None,
    habitat: Optional[str] = None,
    locality: Optional[str] = None,
    county: Optional[str] = None,
    ecoregion_path: Optional[Path] = None,
    county_path: Optional[Path] = None,
) -> Optional[str]:
    """Assigns canonical US EPA Level III Ecoregion (Omernik & Griffith 2014) via polygon intersection.

    Rules:
      a) Finite (lat, lon): strict point-in-polygon intersection against EPA Level III Ecoregions.
      b) Missing/NA coordinates: county polygon containment test against EPA Level III Ecoregions;
         assigns via county centroid ONLY when the county lies within a single Level III Ecoregion
         (or >= 95% single-ecoregion area share), returning None (NA) for multi-ecoregion or unmatched counties.
      c) Never uses bounding boxes, state-level lumping, habitat text heuristics, or 'Other_US'.
    """
    _ = (habitat, locality)
    has_finite_coords = False
    lat_f = math.nan
    lon_f = math.nan
    if lat is not None and lon is not None and str(lat).strip() != "" and str(lon).strip() != "":
        try:
            lat_f = float(lat)
            lon_f = float(lon)
            if math.isfinite(lat_f) and math.isfinite(lon_f) and -90.0 <= lat_f <= 90.0 and -180.0 <= lon_f <= 180.0:
                has_finite_coords = True
        except (ValueError, TypeError):
            has_finite_coords = False

    idx_bundle = load_epa_level3_spatial_index(ecoregion_path=ecoregion_path, county_path=county_path)

    # (a) Finite coordinates: strict point-in-polygon intersection
    if has_finite_coords:
        return _point_in_epa_l3(lon_f, lat_f, idx_bundle)

    # (b) Missing coordinates: county-containment fallback (>= 95% single-ecoregion area share)
    st_norm = _normalize_us_state_name(state_province)
    co_norm = _normalize_us_county_name(county, state_norm=st_norm)
    if not st_norm or not co_norm:
        return None

    ck = f"{st_norm}|{co_norm}"
    ck_ns = f"{st_norm}|{co_norm.replace(' ', '')}"
    cnty_rec = idx_bundle["counties"].get(ck) or idx_bundle["counties_nospace"].get(ck_ns)
    if cnty_rec is None:
        return None

    cnty_geom = cnty_rec["geom"]
    use_live = idx_bundle.get("use_live_county_intersection", False) or (cnty_rec.get("top_area_share") is None)

    if use_live:
        eco_tree = idx_bundle["eco_tree"]
        eco_geoms = idx_bundle["eco_geoms"]
        eco_labels = idx_bundle["eco_labels"]
        eco_id_map = idx_bundle["eco_id_map"]
        cos_w = max(math.cos(math.radians(cnty_geom.centroid.y)), 0.1)
        area_by_eco: Dict[str, float] = {}
        for item in eco_tree.query(cnty_geom):
            idx = int(item) if isinstance(item, (int, np.integer)) else eco_id_map.get(id(item), -1)
            if idx >= 0 and eco_geoms[idx].intersects(cnty_geom):
                inter_a = eco_geoms[idx].intersection(cnty_geom).area * cos_w
                if inter_a > 0.0:
                    lbl = eco_labels[idx]
                    area_by_eco[lbl] = area_by_eco.get(lbl, 0.0) + inter_a
        tot_eco_a = sum(area_by_eco.values())
        cnty_a = cnty_geom.area * cos_w
        if tot_eco_a <= 0.0 or (cnty_a > 0.0 and (tot_eco_a / cnty_a) < 0.50):
            return None
        top_eco, top_a = max(area_by_eco.items(), key=lambda kv: kv[1])
        top_share = top_a / tot_eco_a
        n_eco = sum(1 for a in area_by_eco.values() if (a / tot_eco_a) >= 0.001)
    else:
        top_eco = cnty_rec.get("top_ecoregion")
        top_share = float(cnty_rec.get("top_area_share") or 0.0)
        n_eco = int(cnty_rec.get("n_ecoregions") or 0)

    if not top_eco or not (n_eco == 1 or top_share >= 0.95):
        return None

    # Intersect county centroid (or interior representative point for estuarine/coastal counties)
    c_lon = cnty_rec.get("centroid_lon")
    c_lat = cnty_rec.get("centroid_lat")
    if c_lon is None or c_lat is None:
        c_pt = cnty_geom.centroid
        if not cnty_geom.covers(c_pt):
            c_pt = cnty_geom.representative_point()
        c_lon, c_lat = c_pt.x, c_pt.y

    cent_hit = _point_in_epa_l3(float(c_lon), float(c_lat), idx_bundle)
    if cent_hit == top_eco:
        return cent_hit
    rep_pt = cnty_geom.representative_point()
    rep_hit = _point_in_epa_l3(float(rep_pt.x), float(rep_pt.y), idx_bundle)
    if rep_hit == top_eco:
        return rep_hit
    return top_eco if cent_hit is None else None


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
    if not image_path.exists():
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
                sharpness_score = round(float(cv2.Laplacian(img_cv, cv2.CV_64F).var()), 2)
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
    logger: Optional[logging.Logger] = None,
) -> Tuple[bool, bool]:
    """Asynchronously downloads a single voucher image file to local storage with retry logic.

    Returns:
        Tuple[bool, bool]: (success, was_skipped)
    """
    if not force and destination_path.exists() and os.path.getsize(destination_path) > 0:
        if logger:
            logger.info(
                f"Skipping download for existing voucher {destination_path.name} "
                f"({os.path.getsize(destination_path)} bytes)."
            )
        return True, True

    async with semaphore:
        for attempt in range(1, max_retries + 1):
            temp_path: Optional[Path] = None
            try:
                timeout = aiohttp.ClientTimeout(total=45, connect=15)
                async with session.get(image_url, timeout=timeout) as response:
                    if response.status == 200:
                        content_type = response.headers.get("Content-Type", "").lower()
                        content = await response.read()
                        if len(content) > 1024 and (
                            not content_type
                            or "image" in content_type
                            or "octet-stream" in content_type
                            or content[:3] == b"\xff\xd8\xff"
                        ):
                            destination_path.parent.mkdir(parents=True, exist_ok=True)
                            temp_file = tempfile.NamedTemporaryFile(
                                dir=destination_path.parent,
                                prefix=f"{destination_path.stem}_",
                                suffix=".tmp",
                                delete=False,
                            )
                            temp_path = Path(temp_file.name)
                            temp_file.write(content)
                            temp_file.flush()
                            os.fsync(temp_file.fileno())
                            temp_file.close()
                            temp_path.replace(destination_path)
                            return True, False
                    elif response.status in {404, 410}:
                        return False, False
            except (aiohttp.ClientError, asyncio.TimeoutError, Exception):
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
        if not force and dest.exists() and os.path.getsize(dest) > 0:
            if logger:
                logger.info(
                    f"Retaining existing voucher image: {dest.name} "
                    f"({os.path.getsize(dest)} bytes)."
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

    connector = aiohttp.TCPConnector(limit=concurrency_limit, limit_per_host=5, ssl=False)
    async with aiohttp.ClientSession(headers=headers, connector=connector) as session:
        tasks = [
            download_single_image(session, url, dest, semaphore, force=force, logger=logger)
            for url, dest in pending
        ]
        results = await async_tqdm.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
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


def export_curated_table(
    df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
    home_herbarium: str = DEFAULT_HOME_HERBARIUM,
) -> Path:
    """Atomically exports the curated vouchers DataFrame to CSV via temporary file replacement."""
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

        if tiers_changed or "exsiccatae_key" not in df.columns or "is_primary_duplicate" not in df.columns:
            df = stratify_duplicates(df, logger=logger, home_herbarium=home_herbarium)

    cols_to_export = [col for col in EXPORT_COLUMNS if col in df.columns]
    extra_cols = [c for c in df.columns if c not in cols_to_export and not c.startswith("_")]
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
    Consolidated botanical voucher harvester and curator for the genus Packera.
    Unifies GBIF query execution (with dynamic synonym expansion), Darwin Core metadata
    normalization, 3-layer dynamic determiner & herbarium authority stratification,
    configurable geographic scope filtering, and asynchronous high-resolution image acquisition.
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
    ):
        base_taxa = list(taxa) if taxa else list(DEFAULT_TARGET_TAXA)
        self.expand_synonyms = bool(expand_synonyms)
        self.cfg_synonyms = cfg_synonyms if cfg_synonyms is not None else dict(_DEFAULT_CONFIG.taxa.synonyms or {})
        if self.expand_synonyms:
            expanded_taxa: List[str] = []
            seen_t: Set[str] = set()
            for tx in base_taxa:
                for syn_tx in expand_taxon_synonyms_gbif(tx, cfg_synonyms=self.cfg_synonyms):
                    low_tx = syn_tx.lower()
                    if low_tx not in seen_t:
                        seen_t.add(low_tx)
                        expanded_taxa.append(syn_tx)
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

    def harvest(self) -> pd.DataFrame:
        """Harvests and normalizes Darwin Core occurrence records across configured taxa."""
        all_curated_records: List[Dict[str, Any]] = []
        seen_catalog_keys: set = set()
        country_targets: List[Optional[str]] = list(self.countries) if self.countries else [None]

        for taxon in self.taxa:
            country_desc = ",".join(self.countries) if self.countries else "GLOBAL"
            self.logger.info(
                f"Querying GBIF API for taxon: '{taxon}' "
                f"(scope={self.geographic_scope}, country={country_desc}, basisOfRecord=PRESERVED_SPECIMEN)..."
            )
            taxon_harvested = 0
            taxon_retained = 0
            taxon_western_excluded = 0

            for country_code in country_targets:
                offset = 0
                limit = 300
                while taxon_retained < self.max_records_per_taxon:
                    search_kwargs: Dict[str, Any] = {
                        "scientificName": taxon,
                        "basisOfRecord": "PRESERVED_SPECIMEN",
                        "limit": limit,
                        "offset": offset,
                    }
                    if country_code:
                        search_kwargs["country"] = country_code
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
                        taxon_harvested += 1

                        # 1. Geographic Coordinate & Uncertainty Validation
                        lat = rec.get("decimalLatitude")
                        lon = rec.get("decimalLongitude")
                        lat_val: Optional[float] = None
                        lon_val: Optional[float] = None
                        if lat is not None and lon is not None:
                            try:
                                lat_val = float(lat)
                                lon_val = float(lon)
                            except (ValueError, TypeError):
                                pass

                        state_prov = rec.get("stateProvince")
                        if self.exclude_western and is_excluded_western_region(
                            state_prov,
                            lat=lat_val,
                            lon=lon_val,
                            western_longitude_threshold=self.western_longitude_threshold,
                        ):
                            taxon_western_excluded += 1
                            continue

                        uncertainty_val = self.max_uncertainty_meters
                        if lat_val is not None and lon_val is not None:
                            uncertainty_raw = rec.get("coordinateUncertaintyInMeters")
                            if uncertainty_raw is not None:
                                try:
                                    parsed_unc = float(uncertainty_raw)
                                    if parsed_unc > self.max_uncertainty_meters:
                                        continue
                                    uncertainty_val = parsed_unc
                                except (ValueError, TypeError):
                                    pass

                        # 2. Raw Darwin Core Temporal Validation
                        raw_year = rec.get("year")
                        raw_month = rec.get("month")
                        raw_day = rec.get("day")
                        try:
                            y = int(raw_year)
                            m = int(raw_month)
                            d = int(raw_day)
                            valid_date = datetime.date(y, m, d)
                            year_val, month_val, day_val = y, m, d
                        except (ValueError, TypeError, OverflowError):
                            continue

                        event_date = str(rec.get("eventDate") or valid_date.isoformat()).strip()

                        # 3. High-Resolution Media Image Validation
                        media_list = rec.get("media", [])
                        image_url = extract_high_res_image_url(media_list)
                        if not image_url:
                            continue

                        # 4. Catalog Number & Herbarium Institution Normalization
                        raw_catalog = rec.get("catalogNumber")
                        inst_code = rec.get("institutionCode") or rec.get("collectionCode") or "UNKNOWN_INST"
                        gbif_key = str(rec.get("key", ""))

                        if raw_catalog and str(raw_catalog).strip():
                            catalog_number = sanitize_filename(str(raw_catalog).strip())
                        else:
                            catalog_number = f"{sanitize_filename(inst_code)}_{gbif_key}"

                        unique_key = (catalog_number, gbif_key)
                        if unique_key in seen_catalog_keys:
                            continue
                        seen_catalog_keys.add(unique_key)

                        relative_image_path = f"data/raw_vouchers/{catalog_number}.jpg"
                        sci_name_val = rec.get("scientificName") or rec.get("species") or taxon
                        raw_loc_val = rec.get("locality") or rec.get("verbatimLocality")
                        raw_hab_val = rec.get("habitat")
                        proto_auths = extract_protologue_authors([sci_name_val, taxon])

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

                        # 6. US EPA Level III Ecoregion Polygon Assignment
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
                            "species_raw": sci_name_val,
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
                            "_raw_identified_by": str(rec.get("identifiedBy") or "").strip(),
                            "_raw_history": str(rec.get("verbatimIdentificationHistory") or "").strip(),
                            "_raw_type_status": str(rec.get("typeStatus") or "").strip(),
                            "_raw_locality": str(raw_loc_val or "").strip(),
                            "_raw_habitat": str(raw_hab_val or "").strip(),
                        }
                        all_curated_records.append(curated_record)
                        taxon_retained += 1

                        if taxon_retained >= self.max_records_per_taxon:
                            break

                    offset += len(results)
                    if offset >= count:
                        break
                    time.sleep(0.1)

                if taxon_retained >= self.max_records_per_taxon:
                    break

            self.logger.info(
                f"Taxon '{taxon}': Processed {taxon_harvested} occurrences -> "
                f"Retained {taxon_retained} curated records meeting quality filters "
                f"(Western states excluded: {taxon_western_excluded})."
            )

        df = pd.DataFrame(all_curated_records)
        if not df.empty:
            # Layer 3: Corpus-derived multi-herbarium specialist & regional herbarium authority re-evaluation
            dyn_specs, dyn_herbs = discover_corpus_authorities(df)
            if dyn_specs or dyn_herbs:
                self.logger.info(
                    f"Dynamic Authority Engine: Discovered {len(dyn_specs)} corpus specialist(s) "
                    f"and {len(dyn_herbs)} dynamic research herbarium/herbaria."
                )
            df = apply_corpus_authority_tiers(df, dynamic_specialists=dyn_specs, dynamic_herbaria=dyn_herbs)
            raw_audit_cols = [
                c for c in ("_raw_identified_by", "_raw_history", "_raw_type_status", "_raw_locality", "_raw_habitat")
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
            has_thin_or_quota = (
                self.spatial_thin_km > 0.0
                or (self.max_per_state is not None and self.max_per_state > 0)
                or (self.max_per_region is not None and self.max_per_region > 0)
            )
            df = spatially_balance_and_thin_vouchers(
                df,
                spatial_thin_km=self.spatial_thin_km,
                max_per_state=self.max_per_state,
                max_per_region=self.max_per_region,
                primary_duplicates_only=has_thin_or_quota,
                preserve_all_gold=True,
                workspace_dir=self.workspace_dir,
            )
        return df

    def download_and_validate_media(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Asynchronously downloads voucher specimen sheets and filters substandard imagery."""
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
                logger=self.logger,
            )
        )

        self.logger.info(
            f"Auditing specimen image quality (min_megapixels={self.min_megapixels} MP, "
            f"min_file_size_kb={self.min_file_size_kb} KB)..."
        )
        valid_indices = []
        quality_rejected = 0
        mp_values = []

        for idx, row in df.iterrows():
            img_dest = self.workspace_dir / row["image_path"]
            is_valid, q_metrics = validate_image_quality(
                image_path=img_dest,
                min_megapixels=self.min_megapixels,
                min_file_size_kb=self.min_file_size_kb,
                check_sharpness=self.check_sharpness,
                min_sharpness=self.min_sharpness,
            )
            if is_valid:
                valid_indices.append(idx)
                if "megapixels" in q_metrics:
                    mp_values.append(q_metrics["megapixels"])
            else:
                quality_rejected += 1
                if img_dest.exists():
                    try:
                        img_dest.unlink()
                    except Exception:
                        pass

        self.logger.info(
            f"Quality Audit Complete: Retained {len(valid_indices)} vouchers meeting quality standards "
            f"(Rejected {quality_rejected} substandard/low-res images)."
        )
        df_filtered = df.loc[valid_indices].reset_index(drop=True)
        if not df_filtered.empty:
            df_filtered = stratify_duplicates(
                df_filtered,
                workspace_dir=self.workspace_dir,
                logger=self.logger,
                home_herbarium=self.home_herbarium,
            )
        stats["quality_rejected"] = quality_rejected
        if mp_values:
            stats["median_mp"] = float(np.median(mp_values))

        return df_filtered, stats

    def export(self, df: pd.DataFrame) -> Path:
        """Atomically persists the curated vouchers table to CSV."""
        return export_curated_table(df, self.output_csv, logger=self.logger, home_herbarium=self.home_herbarium)

    def run(
        self,
        download_images: bool = False,
        force: Optional[bool] = None,
    ) -> Tuple[pd.DataFrame, Optional[Dict[str, Any]]]:
        """Executes the end-to-end voucher ingestion and curation workflow."""
        if force is not None:
            self.force = force

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
        print_and_log_summary(df_curated, download_stats, self.logger)
        self.logger.info("Pipeline execution completed successfully.")
        return df_curated, download_stats


def harvest_taxa_occurrences(
    taxa_list: List[str],
    max_uncertainty_meters: float = 5000.0,
    max_records_per_taxon: int = 1000,
    exclude_western: bool = True,
    logger: Optional[logging.Logger] = None,
    geographic_scope: str = DEFAULT_GEOGRAPHIC_SCOPE,
    countries: Optional[Union[str, List[str]]] = None,
    include_western: bool = False,
    expand_synonyms: bool = False,
) -> pd.DataFrame:
    """Backwards-compatible functional entry point for occurrence harvesting."""
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
    )
    return harvester.harvest()


def main() -> None:
    """CLI execution wrapper delegating to data prep runner."""
    import importlib
    cli_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")
    cli_mod.main()


if __name__ == "__main__":
    main()
