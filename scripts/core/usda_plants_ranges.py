"""
===============================================================================
Module: usda_plants_ranges.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Programmatic USDA NRCS PLANTS Database county-level distribution fetcher,
    synonym-aware symbol resolver, local CSV + SQLite (`usda_county_ranges`)
    cache manager, and buffered US county polygon range filter for Step 01
    voucher harvesting.

    Official Public-Domain Endpoints:
      1. Master Checklist (`plantlst.txt`):
         https://plants.sc.egov.usda.gov/DocumentLibrary/Txt/plantlst.txt
      2. Symbol Search (`POST`):
         https://plantsservices.sc.egov.usda.gov/api/PlantSearch
      3. Per-Symbol Distribution CSV (`GET`):
         https://plantsservices.sc.egov.usda.gov/api/PlantProfileDataDownload?symbol=<SYMBOL>
===============================================================================
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import io
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tempfile
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union
import urllib.error
import urllib.parse
import urllib.request

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]
try:
    from shapely.geometry import Point, shape
    from shapely.ops import nearest_points, unary_union
    from shapely.strtree import STRtree
except ImportError:
    Point = shape = nearest_points = unary_union = STRtree = None  # type: ignore[assignment]

from scripts.core.db import (
    get_db_connection,
    is_canonical_or_explicit_db_context,
    resolve_db_path,
)
from scripts.core.epa_ecoregions import (
    DEFAULT_EPA_L3_GEOJSON,
    DEFAULT_US_COUNTIES_GEOJSON,
    _US_STATE_ABBR_TO_FULL,
    _ensure_epa_level3_cache_files,
    _normalize_us_county_name,
    _normalize_us_state_name,
    load_epa_level3_spatial_index,
)
from scripts.core.taxonomy_authority import (
    classify_nomenclatural_type_rank,
    haversine_distance_meters,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_USDA_CACHE_DIR = PROJECT_ROOT / "data" / "environmental" / "usda_plants"
USDA_PLANTLST_URL = "https://plants.sc.egov.usda.gov/DocumentLibrary/Txt/plantlst.txt"
USDA_PLANT_SEARCH_API_URL = "https://plantsservices.sc.egov.usda.gov/api/PlantSearch"
USDA_PLANT_PROFILE_API_URL = "https://plantsservices.sc.egov.usda.gov/api/PlantProfile"
USDA_DISTRIBUTION_DOCS_API_URL = (
    "https://plantsservices.sc.egov.usda.gov/api/PlantProfile/getDownloadDistributionDocumentation"
)
USDA_DISTRIBUTION_CSV_API_URL = "https://plantsservices.sc.egov.usda.gov/api/PlantProfileDataDownload"
USDA_USER_AGENT = "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"

_US_FULL_STATE_NAMES: Set[str] = set(_US_STATE_ABBR_TO_FULL.values())
_EXCLUDED_STATUS_TOKENS = ("excluded", "(x)", "false report", "erroneous")

# Process-level caches
_PLANTLST_INDEX_CACHE: Dict[Tuple[str, float], Dict[str, Any]] = {}
_COUNTY_SPATIAL_GRAPH_CACHE: Dict[Tuple[str, float], Dict[str, Any]] = {}
_USDA_NETWORK_UNAVAILABLE: bool = False


def reset_usda_network_state() -> None:
    """Reset cached network-unavailability state (used by unit tests)."""
    global _USDA_NETWORK_UNAVAILABLE
    _USDA_NETWORK_UNAVAILABLE = False


def _extract_clean_binomial(raw_name: Any) -> Tuple[Optional[str], bool]:
    """
    Extract lowercase `'genus epithet'` binomial and whether the string specifies
    an infraspecific rank (`var.`, `subsp.`, `ssp.`, `f.`).
    """
    if raw_name is None:
        return None, False
    txt = str(raw_name).strip()
    if not txt:
        return None, False
    m = re.match(r"^([A-Za-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", txt)
    if not m:
        return None, False
    genus = m.group(1).lower()
    epithet = m.group(2).lower().strip("-")
    if len(epithet) < 2 or epithet in {"sp", "spp", "cf", "aff", "indet", "hybrid", "unknown", "na", "x"}:
        return None, False
    is_infraspecific = bool(re.search(r"\b(?:var|subsp|ssp|f)\.?\s+[a-z-]{2,}", txt.lower()))
    return f"{genus} {epithet}", is_infraspecific


def _atomic_write_text(target_path: Path, text: str) -> None:
    """Atomically write UTF-8 text via a sibling temporary file + fsync + replace."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(target_path.parent),
        prefix=f".{target_path.stem}_",
        suffix=f"{target_path.suffix}.tmp",
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        tmp_path.replace(target_path)
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise


def _http_get_text(url: str, timeout: float = 12.0) -> Optional[str]:
    """Perform an unauthenticated HTTP GET returning decoded UTF-8 text or None on error."""
    global _USDA_NETWORK_UNAVAILABLE
    if _USDA_NETWORK_UNAVAILABLE:
        return None
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": USDA_USER_AGENT,
            "Accept": "application/json, text/csv, text/plain, */*",
        },
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if getattr(resp, "status", 200) != 200:
                return None
            raw_bytes = resp.read()
            return raw_bytes.decode("utf-8-sig", errors="replace")
    except urllib.error.HTTPError:
        return None
    except Exception:
        _USDA_NETWORK_UNAVAILABLE = True
        return None


def _http_post_text(url: str, payload: Mapping[str, Any], timeout: float = 15.0) -> Optional[str]:
    """Perform an unauthenticated HTTP POST with JSON payload returning raw UTF-8 text (CSV or JSON)."""
    global _USDA_NETWORK_UNAVAILABLE
    if _USDA_NETWORK_UNAVAILABLE:
        return None
    body_bytes = json.dumps(dict(payload)).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=body_bytes,
        headers={
            "User-Agent": USDA_USER_AGENT,
            "Content-Type": "application/json;charset=UTF-8",
            "Accept": "text/csv, application/json, text/plain, */*",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if getattr(resp, "status", 200) != 200:
                return None
            return resp.read().decode("utf-8-sig", errors="replace")
    except urllib.error.HTTPError:
        return None
    except Exception:
        _USDA_NETWORK_UNAVAILABLE = True
        return None


def _http_post_json(url: str, payload: Mapping[str, Any], timeout: float = 10.0) -> Optional[Dict[str, Any]]:
    """Perform an unauthenticated HTTP POST with JSON payload returning parsed dict or None."""
    txt = _http_post_text(url, payload=payload, timeout=timeout)
    if not txt:
        return None
    try:
        parsed = json.loads(txt)
        return parsed if isinstance(parsed, dict) else None
    except Exception:
        return None


def _parse_plantlst_text(text: str) -> Dict[str, Any]:
    """
    Parse USDA PLANTS `plantlst.txt` (`"Symbol","Synonym Symbol","Scientific Name with Author",...`)
    into species-rank and infraspecific binomial -> accepted USDA `Symbol` lookup tables.
    """
    species_binom_to_symbols: Dict[str, List[str]] = {}
    infra_binom_to_symbols: Dict[str, List[str]] = {}
    symbol_to_accepted_name: Dict[str, str] = {}

    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        if not isinstance(row, dict):
            continue
        sym = str(row.get("Symbol") or "").strip().upper()
        syn_sym = str(row.get("Synonym Symbol") or "").strip().upper()
        sci_auth = str(
            row.get("Scientific Name with Author")
            or row.get("Scientific Name with Authors")
            or row.get("Scientific Name")
            or ""
        ).strip()
        if not sym or not sci_auth:
            continue
        binom, is_infra = _extract_clean_binomial(sci_auth)
        if not binom:
            continue
        if not syn_sym and not is_infra and sym not in symbol_to_accepted_name:
            symbol_to_accepted_name[sym] = sci_auth
        target_map = infra_binom_to_symbols if is_infra else species_binom_to_symbols
        bucket = target_map.setdefault(binom, [])
        if sym not in bucket:
            # Prioritize accepted rows (empty Synonym Symbol) at front of list
            if not syn_sym:
                bucket.insert(0, sym)
            else:
                bucket.append(sym)

    return {
        "species_map": species_binom_to_symbols,
        "infra_map": infra_binom_to_symbols,
        "symbol_names": symbol_to_accepted_name,
    }


def load_or_fetch_plantlst_index(
    cache_dir: Optional[Path] = None,
    allow_network: bool = True,
    timeout: float = 20.0,
    logger: Optional[logging.Logger] = None,
) -> Optional[Dict[str, Any]]:
    """Load cached `plantlst.txt` index or download once from USDA PLANTS if allowed."""
    eff_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_USDA_CACHE_DIR
    plantlst_path = eff_dir / "plantlst.txt"

    if not plantlst_path.is_file() and allow_network and not _USDA_NETWORK_UNAVAILABLE:
        if logger:
            logger.info(f"Downloading USDA PLANTS master checklist to cache: {plantlst_path} ...")
        fetched = _http_get_text(USDA_PLANTLST_URL, timeout=timeout)
        if fetched and "Symbol" in fetched[:256]:
            try:
                _atomic_write_text(plantlst_path, fetched)
            except Exception as exc:
                if logger:
                    logger.debug(f"Could not cache plantlst.txt at {plantlst_path}: {exc}")
                return _parse_plantlst_text(fetched)

    if not plantlst_path.is_file() or plantlst_path.stat().st_size == 0:
        return None

    mtime = plantlst_path.stat().st_mtime
    cache_key = (str(plantlst_path.resolve()), mtime)
    if cache_key in _PLANTLST_INDEX_CACHE:
        return _PLANTLST_INDEX_CACHE[cache_key]

    try:
        text = plantlst_path.read_text(encoding="utf-8-sig", errors="replace")
        parsed = _parse_plantlst_text(text)
        _PLANTLST_INDEX_CACHE[cache_key] = parsed
        return parsed
    except Exception:
        return None


def resolve_usda_symbols_for_taxon(
    canonical_packera: str,
    synonym_names: Optional[Iterable[str]] = None,
    cache_dir: Optional[Path] = None,
    allow_network: bool = True,
    logger: Optional[logging.Logger] = None,
) -> List[str]:
    """
    Resolve one or more USDA PLANTS accepted Symbols for `canonical_packera`
    across its homotypic and heterotypic synonyms (e.g., mapping `Packera dubia`
    through `Packera tomentosa` / `Senecio tomentosus` -> `PATO4`).
    """
    candidates: List[str] = []
    seen_cand: Set[str] = set()

    def _add_cand(raw: Any) -> None:
        b, _ = _extract_clean_binomial(raw)
        if b and b not in seen_cand:
            seen_cand.add(b)
            candidates.append(b)

    _add_cand(canonical_packera)
    for syn in synonym_names or ():
        _add_cand(syn)
    fallback_map = {
        "packera dubia": ("Packera tomentosa", "Senecio tomentosus", "Senecio dubius", "Cineraria tomentosa"),
        "packera crawfordii": ("Senecio crawfordii",),
        "packera anonyma": ("Packera anonyma", "Senecio anonymus", "Senecio smallii", "Packera smallii"),
        "packera smallii": ("Packera anonyma", "Senecio smallii", "Senecio anonymus"),
    }
    for fallback_syn in fallback_map.get(canonical_packera.strip().lower(), ()):
        _add_cand(fallback_syn)

    if not candidates:
        return []

    eff_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_USDA_CACHE_DIR
    symbols: List[str] = []

    # 1. Check local / downloaded plantlst.txt index first (O(1), zero per-species API latency)
    idx = load_or_fetch_plantlst_index(cache_dir=eff_dir, allow_network=allow_network, logger=logger)
    if idx is not None:
        sp_map: Dict[str, List[str]] = idx["species_map"]
        infra_map: Dict[str, List[str]] = idx["infra_map"]
        for cand in candidates:
            for sym in sp_map.get(cand, ()):
                if sym not in symbols:
                    symbols.append(sym)
        if not symbols:
            for cand in candidates:
                for sym in infra_map.get(cand, ()):
                    if sym not in symbols:
                        symbols.append(sym)
        if symbols:
            return symbols

    # 2. Fallback: POST /api/PlantSearch per candidate binomial
    if allow_network and not _USDA_NETWORK_UNAVAILABLE:
        for cand in candidates:
            title_cand = f"{cand.split(' ', 1)[0].capitalize()} {cand.split(' ', 1)[1]}"
            resp = _http_post_json(
                USDA_PLANT_SEARCH_API_URL,
                {"SearchCriterion": title_cand, "SearchType": 1},
            )
            if not resp:
                continue
            results = resp.get("PlantSearchResults") or []
            if not isinstance(results, list):
                continue
            for item in results:
                if not isinstance(item, dict):
                    continue
                sym = str(item.get("AcceptedSymbol") or item.get("Symbol") or "").strip().upper()
                sci = str(item.get("ScientificName") or "").strip()
                rank = str(item.get("Rank") or "").strip().lower()
                item_binom, is_infra = _extract_clean_binomial(sci)
                if sym and item_binom in seen_cand and not is_infra and rank in {"", "species"}:
                    if sym not in symbols:
                        symbols.append(sym)
            if symbols:
                break

    return symbols


def _slice_csv_from_header(raw_text: str) -> str:
    """Strip any leading metadata/title lines prior to the actual CSV header row."""
    cleaned = raw_text.lstrip("\ufeff").strip()
    if not cleaned:
        return ""
    lines = cleaned.splitlines()
    for idx, ln in enumerate(lines[:12]):
        low = ln.lower()
        if "," in low and any(h in low for h in ("state", "county", "fips", "symbol")):
            return "\n".join(lines[idx:])
    return cleaned


def _convert_usda_json_records_to_csv(records: Sequence[Mapping[str, Any]], default_symbol: str = "") -> str:
    """Serialize USDA JSON distribution dicts into canonical 7-column CSV text."""
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["Symbol", "Scientific Name", "Common Name", "State", "County", "FIPS Code", "Distribution Status"])
    wrote = 0
    for rec in records:
        if not isinstance(rec, Mapping):
            continue
        norm = {re.sub(r"[\s_]+", " ", str(k or "").strip().lower()): v for k, v in rec.items()}
        sym = str(norm.get("symbol") or norm.get("acceptedsymbol") or default_symbol).strip().upper()
        sci = str(norm.get("scientific name") or norm.get("scientificname") or "").strip()
        com = str(norm.get("common name") or norm.get("commonname") or "").strip()
        st = str(
            norm.get("state")
            or norm.get("state and province")
            or norm.get("statename")
            or norm.get("state name")
            or norm.get("stateprovince")
            or norm.get("location")
            or ""
        ).strip()
        co = str(
            norm.get("county")
            or norm.get("county name")
            or norm.get("countyname")
            or ""
        ).strip()
        st_fip = str(norm.get("state fip") or norm.get("statefip") or "").strip()
        co_fip = str(norm.get("county fip") or norm.get("countyfip") or "").strip()
        fips = str(norm.get("fips code") or norm.get("fips") or norm.get("fipscode") or "").strip()
        if not fips and (st_fip or co_fip):
            fips = f"{st_fip.zfill(2)}{co_fip.zfill(3)}" if (st_fip and co_fip) else (co_fip or st_fip)
        status = str(
            norm.get("distribution status")
            or norm.get("distributionstatus")
            or norm.get("native status")
            or norm.get("nativestatus")
            or norm.get("status")
            or ""
        ).strip()
        if st:
            writer.writerow([sym, sci, com, st, co, fips, status])
            wrote += 1
    return buf.getvalue() if wrote > 0 else ""


def _normalize_usda_distribution_payload(raw_text: Optional[str], default_symbol: str = "") -> Optional[str]:
    """Accept raw CSV text, JSON-wrapped CSV text, or JSON distribution objects/arrays from USDA PLANTS."""
    if not raw_text:
        return None
    s = raw_text.lstrip("\ufeff").strip()
    if not s:
        return None

    if s[0] in ('{', '[', '"') and (s[0] != '"' or s.endswith('"')):
        try:
            obj = json.loads(s)
        except Exception:
            obj = None
        if isinstance(obj, str):
            return _normalize_usda_distribution_payload(obj, default_symbol=default_symbol)
        if isinstance(obj, list):
            csv_conv = _convert_usda_json_records_to_csv(obj, default_symbol=default_symbol)
            return csv_conv or None
        if isinstance(obj, dict):
            for str_key in ("Csv", "CsvData", "csv", "csvData", "FileContents", "Content", "Data", "data"):
                val = obj.get(str_key)
                if isinstance(val, str) and "," in val:
                    res = _normalize_usda_distribution_payload(val, default_symbol=default_symbol)
                    if res:
                        return res
            for list_key in (
                "PlantsDistributionResults",
                "plantsDistributionResults",
                "DistributionDocumentation",
                "Distributions",
                "distributions",
                "Counties",
                "counties",
                "Data",
                "data",
                "Results",
                "results",
            ):
                val = obj.get(list_key)
                if isinstance(val, list) and val:
                    csv_conv = _convert_usda_json_records_to_csv(val, default_symbol=default_symbol)
                    if csv_conv:
                        return csv_conv
            return None

    sliced = _slice_csv_from_header(s)
    parsed_probe = parse_usda_distribution_csv(sliced, default_symbol=default_symbol)
    if parsed_probe["counties"] or parsed_probe["state_only_states"]:
        return sliced
    return None


def parse_usda_distribution_csv(
    csv_text: str,
    default_symbol: str = "",
) -> Dict[str, Any]:
    """
    Parse USDA PLANTS distribution CSV across both USDA PLANTS 3.0 export schemas:
      - `Symbol,Scientific Name,Common Name,State and Province,County Name,FIPS Code,Distribution Status`
      - `Symbol,Synonym Symbol,Scientific Name,State FIP,State,County FIP,County,Native Status`

    Distinguishes county-resolved states from state-only (`No County Data`) states:
    a state is marked `state_only` ONLY when USDA records non-excluded presence for
    that state with blank `County` and zero county-level rows in that state.
    """
    counties: Set[str] = set()
    counties_by_state: Dict[str, Set[str]] = {}
    present_us_states: Set[str] = set()
    db_rows: List[Dict[str, Any]] = []

    if not csv_text or not csv_text.strip():
        return {
            "counties": counties,
            "state_only_states": set(),
            "documented_states": set(),
            "rows": db_rows,
        }

    sliced_csv = _slice_csv_from_header(csv_text)
    reader = csv.DictReader(io.StringIO(sliced_csv))
    raw_parsed_rows: List[Tuple[str, str, Optional[str], str, str]] = []

    for row in reader:
        if not isinstance(row, dict):
            continue
        norm_row = {
            re.sub(r"[\s_]+", " ", str(k or "").strip().lower()): str(v or "").strip()
            for k, v in row.items()
            if k is not None
        }
        sym = (norm_row.get("symbol") or default_symbol).strip().upper()
        raw_state = (
            norm_row.get("state")
            or norm_row.get("state and province")
            or norm_row.get("stateprovince")
            or norm_row.get("state name")
            or norm_row.get("statename")
            or ""
        )
        raw_county = (
            norm_row.get("county")
            or norm_row.get("county name")
            or norm_row.get("countyname")
            or ""
        )
        st_fip = norm_row.get("state fip") or norm_row.get("statefip") or ""
        co_fip = norm_row.get("county fip") or norm_row.get("countyfip") or ""
        fips = norm_row.get("fips code") or norm_row.get("fips") or norm_row.get("fipscode") or ""
        if not fips and (st_fip or co_fip):
            fips = f"{st_fip.zfill(2)}{co_fip.zfill(3)}" if (st_fip and co_fip) else (co_fip or st_fip)
        native_status = (
            norm_row.get("native status")
            or norm_row.get("distribution status")
            or norm_row.get("nativestatus")
            or norm_row.get("distributionstatus")
            or norm_row.get("status")
            or ""
        )

        low_status = native_status.lower()
        if any(tok in low_status for tok in _EXCLUDED_STATUS_TOKENS):
            continue

        st_norm = _normalize_us_state_name(raw_state)
        if not st_norm or st_norm not in _US_FULL_STATE_NAMES:
            continue

        present_us_states.add(st_norm)
        co_norm = _normalize_us_county_name(raw_county, state_norm=st_norm) if raw_county else None
        raw_parsed_rows.append((sym, st_norm, co_norm, fips, native_status))
        if co_norm:
            ck = f"{st_norm}|{co_norm}"
            counties.add(ck)
            counties_by_state.setdefault(st_norm, set()).add(co_norm)

    state_only_states: Set[str] = {
        st for st in present_us_states if not counties_by_state.get(st)
    }
    documented_states: Set[str] = set(counties_by_state.keys()) | state_only_states

    seen_pk: Set[Tuple[str, str]] = set()
    for sym, st_norm, co_norm, fips, native_status in raw_parsed_rows:
        if co_norm:
            pk = (st_norm, co_norm)
            if pk in seen_pk:
                continue
            seen_pk.add(pk)
            db_rows.append({
                "usda_symbol": sym,
                "state_norm": st_norm,
                "county_norm": co_norm,
                "fips_code": fips,
                "native_status": native_status,
                "is_state_only": 0,
            })
        elif st_norm in state_only_states:
            pk = (st_norm, "")
            if pk in seen_pk:
                continue
            seen_pk.add(pk)
            db_rows.append({
                "usda_symbol": sym,
                "state_norm": st_norm,
                "county_norm": "",
                "fips_code": fips,
                "native_status": native_status,
                "is_state_only": 1,
            })

    return {
        "counties": counties,
        "state_only_states": state_only_states,
        "documented_states": documented_states,
        "rows": db_rows,
    }


def _load_taxon_range_from_sqlite(
    conn: sqlite3.Connection,
    canonical_packera: str,
) -> Optional[Dict[str, Any]]:
    """Load persisted USDA county distribution for `canonical_packera` from SQLite if present."""
    try:
        cur = conn.execute(
            """
            SELECT usda_symbol, state_norm, county_norm, fips_code, native_status, is_state_only
            FROM usda_county_ranges
            WHERE canonical_packera = ?;
            """,
            (canonical_packera,),
        )
        fetched = cur.fetchall()
    except sqlite3.Error:
        return None

    if not fetched:
        return None

    counties: Set[str] = set()
    state_only_states: Set[str] = set()
    documented_states: Set[str] = set()
    symbols: List[str] = []

    for r in fetched:
        sym = str(r[0] or "").strip().upper()
        st_norm = str(r[1] or "").strip().lower()
        co_norm = str(r[2] or "").strip().lower()
        is_st_only = bool(int(r[5] or 0))
        if sym and sym not in symbols:
            symbols.append(sym)
        if not st_norm:
            continue
        documented_states.add(st_norm)
        if is_st_only or not co_norm:
            state_only_states.add(st_norm)
        else:
            counties.add(f"{st_norm}|{co_norm}")

    state_only_states = {
        st for st in state_only_states
        if not any(ck.startswith(f"{st}|") for ck in counties)
    }
    return {
        "canonical_packera": canonical_packera,
        "symbols": symbols,
        "counties": counties,
        "state_only_states": state_only_states,
        "documented_states": documented_states,
    }


def _persist_taxon_range_to_sqlite(
    conn: sqlite3.Connection,
    canonical_packera: str,
    rows: Sequence[Mapping[str, Any]],
) -> int:
    """Persist normalized USDA county presence rows for `canonical_packera` into `usda_county_ranges`."""
    if not rows:
        return 0
    now_iso = datetime.now(timezone.utc).isoformat()
    payload: List[Tuple[str, str, str, str, Optional[str], Optional[str], int, str]] = []
    for r in rows:
        st_norm = str(r.get("state_norm") or "").strip().lower()
        co_norm = str(r.get("county_norm") or "").strip().lower()
        sym = str(r.get("usda_symbol") or "").strip().upper()
        if not st_norm or not sym:
            continue
        payload.append((
            canonical_packera,
            sym,
            st_norm,
            co_norm,
            str(r.get("fips_code") or "").strip() or None,
            str(r.get("native_status") or "").strip() or None,
            1 if r.get("is_state_only") else 0,
            now_iso,
        ))
    if not payload:
        return 0

    sql = """
    INSERT INTO usda_county_ranges (
        canonical_packera, usda_symbol, state_norm, county_norm,
        fips_code, native_status, is_state_only, updated_at
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(canonical_packera, state_norm, county_norm) DO UPDATE SET
        usda_symbol = excluded.usda_symbol,
        fips_code = COALESCE(excluded.fips_code, usda_county_ranges.fips_code),
        native_status = COALESCE(excluded.native_status, usda_county_ranges.native_status),
        is_state_only = excluded.is_state_only,
        updated_at = excluded.updated_at;
    """
    with conn:
        conn.executemany(sql, payload)
    return len(payload)


def _fetch_usda_symbol_distribution_csv(sym: str, timeout: float = 15.0) -> Optional[str]:
    """
    Download county distribution CSV for USDA symbol `sym` using:
      1. `GET /api/PlantProfile?symbol=<SYM>` (extracts numeric `Id`/`MasterId` or inline records,
         and short-circuits if mocked `_http_get_text` returns CSV directly).
      2. `POST /api/PlantProfile/getDownloadDistributionDocumentation` with `MasterId` + `Symbol`.
      3. `GET /api/PlantProfileDataDownload?symbol=<SYM>` legacy fallback.
    """
    master_id: Optional[int] = None
    prof_url = f"{USDA_PLANT_PROFILE_API_URL}?{urllib.parse.urlencode({'symbol': sym})}"
    prof_raw = _http_get_text(prof_url, timeout=timeout)
    if prof_raw:
        direct_norm = _normalize_usda_distribution_payload(prof_raw, default_symbol=sym)
        if direct_norm:
            return direct_norm
        try:
            prof_obj = json.loads(prof_raw)
            if isinstance(prof_obj, dict):
                raw_id = (
                    prof_obj.get("Id")
                    or prof_obj.get("MasterId")
                    or prof_obj.get("PlantId")
                    or prof_obj.get("id")
                    or prof_obj.get("masterId")
                )
                if raw_id is not None and str(raw_id).strip().isdigit():
                    master_id = int(str(raw_id).strip())
        except Exception:
            pass

    if master_id is None and not _USDA_NETWORK_UNAVAILABLE:
        srch = _http_post_json(USDA_PLANT_SEARCH_API_URL, {"SearchCriterion": sym, "SearchType": 0}, timeout=10.0)
        if srch and isinstance(srch.get("PlantSearchResults"), list):
            for item in srch["PlantSearchResults"]:
                if isinstance(item, dict):
                    item_sym = str(item.get("Symbol") or item.get("AcceptedSymbol") or "").strip().upper()
                    raw_id = item.get("Id") or item.get("MasterId") or item.get("PlantId")
                    if item_sym == sym and raw_id is not None and str(raw_id).strip().isdigit():
                        master_id = int(str(raw_id).strip())
                        break

    if not _USDA_NETWORK_UNAVAILABLE:
        post_payload: Dict[str, Any] = {
            "Field": "Symbol",
            "SortBy": "sortSciName",
            "Offset": None,
            "MasterId": master_id,
            "Symbol": sym,
            "Value": sym,
        }
        docs_raw = _http_post_text(USDA_DISTRIBUTION_DOCS_API_URL, payload=post_payload, timeout=timeout)
        docs_norm = _normalize_usda_distribution_payload(docs_raw, default_symbol=sym)
        if docs_norm:
            return docs_norm

    if not _USDA_NETWORK_UNAVAILABLE:
        dl_url = f"{USDA_DISTRIBUTION_CSV_API_URL}?{urllib.parse.urlencode({'symbol': sym})}"
        legacy_raw = _http_get_text(dl_url, timeout=timeout)
        legacy_norm = _normalize_usda_distribution_payload(legacy_raw, default_symbol=sym)
        if legacy_norm:
            return legacy_norm

    return None


def fetch_or_load_usda_taxon_range(
    canonical_packera: str,
    synonym_names: Optional[Iterable[str]] = None,
    cache_dir: Optional[Path] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Path = PROJECT_ROOT,
    allow_network: bool = True,
    force_refresh: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """
    Load or download the unioned USDA PLANTS county range for `canonical_packera`
    across all mapped USDA symbols and synonyms, persisting to CSV cache and SQLite.
    """
    eff_cache_dir = Path(cache_dir) if cache_dir is not None else (Path(workspace_dir) / "data" / "environmental" / "usda_plants")
    dist_dir = eff_cache_dir / "distributions"

    # 1. Check SQLite cache first unless force_refresh=True
    use_sqlite = db_path is not None or (
        os.environ.get("PACKERA_DEMO_MODE") != "1"
        and is_canonical_or_explicit_db_context(eff_cache_dir, workspace_root=workspace_dir)
    )
    resolved_db: Optional[Union[str, Path]] = None
    if use_sqlite:
        try:
            resolved_db = resolve_db_path(db_path=db_path, workspace_root=workspace_dir)
        except Exception:
            resolved_db = None

    if not force_refresh and resolved_db is not None and resolved_db != ":memory:" and Path(resolved_db).is_file():
        try:
            conn = get_db_connection(db_path=resolved_db, read_only=True, workspace_root=workspace_dir)
            try:
                cached_sql = _load_taxon_range_from_sqlite(conn, canonical_packera)
                if cached_sql is not None and (cached_sql["counties"] or cached_sql["state_only_states"]):
                    return cached_sql
            finally:
                conn.close()
        except Exception:
            pass

    # 2. Check taxon-named local CSV override/cache `<Packera_epithet>.csv` if present
    slug_csv = dist_dir / f"{re.sub(r'[^A-Za-z0-9_-]+', '_', canonical_packera.strip())}.csv"
    if not force_refresh and slug_csv.is_file() and slug_csv.stat().st_size > 0:
        parsed_slug = parse_usda_distribution_csv(
            slug_csv.read_text(encoding="utf-8-sig", errors="replace"),
            default_symbol=canonical_packera.upper(),
        )
        syms = sorted({str(r["usda_symbol"]) for r in parsed_slug["rows"] if r.get("usda_symbol")})
        result_slug = {
            "canonical_packera": canonical_packera,
            "symbols": syms,
            "counties": parsed_slug["counties"],
            "state_only_states": parsed_slug["state_only_states"],
            "documented_states": parsed_slug["documented_states"],
        }
        if resolved_db is not None and parsed_slug["rows"]:
            try:
                conn_w = get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=workspace_dir)
                try:
                    _persist_taxon_range_to_sqlite(conn_w, canonical_packera, parsed_slug["rows"])
                finally:
                    conn_w.close()
            except Exception:
                pass
        return result_slug

    # 3. Resolve USDA symbol(s) via checklist / PlantSearch API
    symbols = resolve_usda_symbols_for_taxon(
        canonical_packera=canonical_packera,
        synonym_names=synonym_names,
        cache_dir=eff_cache_dir,
        allow_network=allow_network,
        logger=logger,
    )
    if not symbols:
        return {
            "canonical_packera": canonical_packera,
            "symbols": [],
            "counties": set(),
            "state_only_states": set(),
            "documented_states": set(),
        }

    all_counties: Set[str] = set()
    all_state_only: Set[str] = set()
    all_doc_states: Set[str] = set()
    all_rows: List[Dict[str, Any]] = []

    for sym in symbols:
        sym_csv_path = dist_dir / f"{sym}.csv"
        csv_text: Optional[str] = None
        if not force_refresh and sym_csv_path.is_file() and sym_csv_path.stat().st_size > 0:
            try:
                csv_text = sym_csv_path.read_text(encoding="utf-8-sig", errors="replace")
            except Exception:
                csv_text = None

        if csv_text is not None:
            probe_cached = parse_usda_distribution_csv(csv_text, default_symbol=sym)
            if not probe_cached["counties"] and not probe_cached["state_only_states"] and allow_network:
                csv_text = None

        if csv_text is None and allow_network and not _USDA_NETWORK_UNAVAILABLE:
            csv_text = _fetch_usda_symbol_distribution_csv(sym, timeout=15.0)
            if csv_text:
                try:
                    _atomic_write_text(sym_csv_path, csv_text)
                except Exception as exc:
                    if logger:
                        logger.debug(f"Could not write cached USDA CSV {sym_csv_path}: {exc}")

        if not csv_text:
            continue

        parsed = parse_usda_distribution_csv(csv_text, default_symbol=sym)
        all_counties.update(parsed["counties"])
        all_state_only.update(parsed["state_only_states"])
        all_doc_states.update(parsed["documented_states"])
        all_rows.extend(parsed["rows"])

    # Any state with >=1 county record across synonyms is county-resolved rather than state-only
    all_state_only = {
        st for st in all_state_only
        if not any(ck.startswith(f"{st}|") for ck in all_counties)
    }
    filtered_rows = [
        r for r in all_rows
        if not (r.get("is_state_only") and r.get("state_norm") not in all_state_only)
    ]

    if resolved_db is not None and filtered_rows:
        try:
            conn_w = get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=workspace_dir)
            try:
                _persist_taxon_range_to_sqlite(conn_w, canonical_packera, filtered_rows)
            finally:
                conn_w.close()
        except Exception as db_exc:
            if logger:
                logger.debug(f"SQLite usda_county_ranges persistence skipped: {db_exc}")

    return {
        "canonical_packera": canonical_packera,
        "symbols": symbols,
        "counties": all_counties,
        "state_only_states": all_state_only,
        "documented_states": all_doc_states,
    }


def _ensure_standalone_us_counties_geojson(cnty_geojson_path: Path) -> bool:
    """Build `us_counties.geojson` locally via R `maps` + `sf` if absent and Rscript is installed."""
    if cnty_geojson_path.is_file() and cnty_geojson_path.stat().st_size > 0:
        return True
    rscript_bin = shutil.which("Rscript")
    if not rscript_bin:
        return False
    cnty_geojson_path.parent.mkdir(parents=True, exist_ok=True)
    r_code = (
        "suppressPackageStartupMessages({library(sf); library(maps)}); sf_use_s2(FALSE);"
        "cm <- st_make_valid(st_as_sf(maps::map('county', fill=TRUE, plot=FALSE))); st_crs(cm) <- 4326;"
        "pts <- strsplit(as.character(cm$ID), ',');"
        "st_n <- vapply(pts, `[`, character(1), 1);"
        "co_n <- trimws(gsub('[^a-z0-9]+', ' ', gsub('^saint\\\\s+|^st\\\\.\\\\s*', 'st ', gsub('\\\\s+city$', '', sub(':.*$', '', vapply(pts, `[`, character(1), 2))))));"
        "ck <- paste0(st_n, '|', co_n); dk <- unique(ck[duplicated(ck)]); sm <- !(ck %in% dk); ga <- st_geometry(cm);"
        "ok <- ck[sm]; os <- st_n[sm]; oc <- co_n[sm]; og <- lapply(which(sm), function(i) st_cast(ga[[i]], 'MULTIPOLYGON'));"
        "for(d in dk){ id <- which(ck==d); ok <- c(ok, d); os <- c(os, st_n[id[1]]); oc <- c(oc, co_n[id[1]]); og[[length(og)+1L]] <- st_cast(st_make_valid(st_union(ga[id]))[[1]], 'MULTIPOLYGON'); };"
        "cd <- st_sf(county_key=ok, state_name=os, county_norm=oc, geometry=st_sfc(og, crs=4326));"
        "cp <- st_coordinates(suppressWarnings(st_point_on_surface(cd)));"
        "cd$centroid_lon <- cp[,1]; cd$centroid_lat <- cp[,2];"
        f"tc <- {json.dumps(str(cnty_geojson_path) + '.tmp.geojson')};"
        "st_write(cd, tc, driver='GeoJSON', quiet=TRUE, delete_dsn=TRUE);"
        "file.rename(tc, sub('\\\\.tmp\\\\.geojson$', '', tc));"
    )
    try:
        subprocess.run([rscript_bin, "-e", r_code], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
    except Exception:
        return False
    return cnty_geojson_path.is_file() and cnty_geojson_path.stat().st_size > 0


def load_us_county_spatial_graph(
    county_path: Optional[Path] = None,
    ecoregion_path: Optional[Path] = None,
    allow_download: bool = True,
) -> Dict[str, Any]:
    """
    Build and cache a spatial index (`STRtree`) and 1-ring topological county
    adjacency map over `us_counties.geojson`.
    """
    cnty_p = Path(county_path) if county_path is not None else DEFAULT_US_COUNTIES_GEOJSON
    eco_p = Path(ecoregion_path) if ecoregion_path is not None else DEFAULT_EPA_L3_GEOJSON

    if not cnty_p.exists() and county_path is None:
        if allow_download:
            _ensure_epa_level3_cache_files(eco_p, cnty_p)
        if not cnty_p.exists():
            _ensure_standalone_us_counties_geojson(cnty_p)

    mtime = cnty_p.stat().st_mtime if cnty_p.exists() else -1.0
    cache_key = (str(cnty_p.resolve()) if cnty_p.exists() else "OFFLINE_FALLBACK", mtime)
    if cache_key in _COUNTY_SPATIAL_GRAPH_CACHE:
        return _COUNTY_SPATIAL_GRAPH_CACHE[cache_key]

    idx_bundle = load_epa_level3_spatial_index(
        ecoregion_path=eco_p if eco_p.exists() else ecoregion_path,
        county_path=cnty_p if cnty_p.exists() else county_path,
        allow_download=False,
    )
    counties_dict: Dict[str, Dict[str, Any]] = dict(idx_bundle.get("counties") or {})
    # If eco_p did not exist yet, load cnty_p directly if it exists
    if cnty_p.exists() and len(counties_dict) <= 10:
        try:
            cnty_fc = json.loads(cnty_p.read_text(encoding="utf-8"))
            counties_dict = {}
            for feat in cnty_fc.get("features", []):
                props = feat.get("properties") or {}
                st_n = _normalize_us_state_name(props.get("state_name") or props.get("STATE_NAME") or props.get("state"))
                co_n = _normalize_us_county_name(props.get("county_norm") or props.get("NAME") or props.get("county"), st_n)
                geom_dict = feat.get("geometry")
                if st_n and co_n and geom_dict and shape is not None:
                    geom_obj = shape(geom_dict)
                    if not geom_obj.is_empty:
                        c_pt = geom_obj.representative_point()
                        counties_dict[f"{st_n}|{co_n}"] = {
                            "geom": geom_obj,
                            "centroid_lon": float(props["centroid_lon"]) if props.get("centroid_lon") is not None else float(c_pt.x),
                            "centroid_lat": float(props["centroid_lat"]) if props.get("centroid_lat") is not None else float(c_pt.y),
                        }
        except Exception:
            pass

    counties_nospace: Dict[str, str] = {
        f"{ck.split('|', 1)[0]}|{ck.split('|', 1)[1].replace(' ', '')}": ck
        for ck in counties_dict
        if "|" in ck
    }

    ck_list: List[str] = []
    geoms: List[Any] = []
    for ck, rec in counties_dict.items():
        g = rec.get("geom")
        if g is not None and not g.is_empty:
            ck_list.append(ck)
            geoms.append(g)

    tree = STRtree(geoms) if (STRtree is not None and geoms) else None
    geom_id_map = {id(g): i for i, g in enumerate(geoms)}
    neighbors: Dict[str, Set[str]] = {ck: set() for ck in ck_list}

    if tree is not None:
        for i, g_i in enumerate(geoms):
            ck_i = ck_list[i]
            # 1e-3 deg (~100 m) tolerance bridges simplified polygon vertices across county lines
            buffered_i = g_i.buffer(1e-3)
            for hit in tree.query(buffered_i):
                j = int(hit) if (np is not None and isinstance(hit, (int, np.integer))) or isinstance(hit, int) else geom_id_map.get(id(hit), -1)
                if j >= 0 and j != i and buffered_i.intersects(geoms[j]):
                    ck_j = ck_list[j]
                    neighbors[ck_i].add(ck_j)
                    neighbors[ck_j].add(ck_i)

    graph = {
        "counties": counties_dict,
        "counties_nospace": counties_nospace,
        "ck_list": ck_list,
        "geoms": geoms,
        "tree": tree,
        "geom_id_map": geom_id_map,
        "neighbors": neighbors,
        "is_complete_national_layer": len(counties_dict) >= 500,
    }
    _COUNTY_SPATIAL_GRAPH_CACHE[cache_key] = graph
    return graph


def resolve_voucher_us_county(
    lat: Optional[float],
    lon: Optional[float],
    state_province: Optional[str],
    county: Optional[str],
    spatial_graph: Mapping[str, Any],
) -> Tuple[Optional[str], Optional[str], bool]:
    """
    Resolve `(state_norm, county_norm, is_non_us)` for a voucher using:
      1. Point-in-polygon against `us_counties.geojson` when finite coordinates exist
         (with a 0.02 deg coastal barrier-island snap when state matches or is blank).
      2. Normalized Darwin Core `(stateProvince, county)` fallback.
    """
    st_dwc = _normalize_us_state_name(state_province)
    is_explicit_non_us = bool(st_dwc and st_dwc not in _US_FULL_STATE_NAMES)

    has_coords = False
    lat_f = math.nan
    lon_f = math.nan
    if lat is not None and lon is not None:
        try:
            lat_f = float(lat)
            lon_f = float(lon)
            if math.isfinite(lat_f) and math.isfinite(lon_f) and -90.0 <= lat_f <= 90.0 and -180.0 <= lon_f <= 180.0:
                has_coords = True
        except (TypeError, ValueError):
            has_coords = False

    tree = spatial_graph.get("tree")
    geoms = spatial_graph.get("geoms") or []
    ck_list = spatial_graph.get("ck_list") or []
    geom_id_map = spatial_graph.get("geom_id_map") or {}

    if has_coords and tree is not None and Point is not None:
        pt = Point(lon_f, lat_f)
        for hit in tree.query(pt):
            idx = int(hit) if (np is not None and isinstance(hit, (int, np.integer))) or isinstance(hit, int) else geom_id_map.get(id(hit), -1)
            if idx >= 0 and geoms[idx].covers(pt):
                st_hit, co_hit = ck_list[idx].split("|", 1)
                return st_hit, co_hit, False

        # Coastal / estuarine tolerance (~2 km)
        pt_buf = pt.buffer(0.02)
        best_idx = -1
        best_dist = float("inf")
        for hit in tree.query(pt_buf):
            idx = int(hit) if (np is not None and isinstance(hit, (int, np.integer))) or isinstance(hit, int) else geom_id_map.get(id(hit), -1)
            if idx >= 0 and geoms[idx].intersects(pt_buf):
                st_cand = ck_list[idx].split("|", 1)[0]
                if st_dwc and st_dwc in _US_FULL_STATE_NAMES and st_cand != st_dwc:
                    continue
                d_val = float(geoms[idx].distance(pt))
                if d_val < best_dist:
                    best_dist = d_val
                    best_idx = idx
        if best_idx >= 0:
            st_hit, co_hit = ck_list[best_idx].split("|", 1)
            return st_hit, co_hit, False

    if is_explicit_non_us:
        return st_dwc, None, True

    if st_dwc and st_dwc in _US_FULL_STATE_NAMES:
        co_dwc = _normalize_us_county_name(county, state_norm=st_dwc)
        if co_dwc:
            ns_key = f"{st_dwc}|{co_dwc.replace(' ', '')}"
            canon_ck = (spatial_graph.get("counties_nospace") or {}).get(ns_key)
            if canon_ck:
                _, co_dwc = canon_ck.split("|", 1)
        return st_dwc, co_dwc, False

    # Coordinates outside US without US state label (e.g., Canada/Mexico)
    if has_coords and spatial_graph.get("is_complete_national_layer"):
        return None, None, True

    return None, None, False


def _is_within_geodesic_buffer(
    lat_f: Optional[float],
    lon_f: Optional[float],
    voucher_ck: Optional[str],
    doc_counties: Set[str],
    doc_union_geom: Any,
    spatial_graph: Mapping[str, Any],
    buffer_km: float,
) -> bool:
    """Return True if the voucher lies within `buffer_km` of any documented county."""
    if buffer_km <= 0.0 or not doc_counties:
        return False
    buffer_m = float(buffer_km) * 1000.0
    counties_dict: Mapping[str, Mapping[str, Any]] = spatial_graph.get("counties") or {}

    q_lat = lat_f
    q_lon = lon_f
    if (q_lat is None or q_lon is None) and voucher_ck and voucher_ck in counties_dict:
        c_rec = counties_dict[voucher_ck]
        if c_rec.get("centroid_lat") is not None and c_rec.get("centroid_lon") is not None:
            q_lat = float(c_rec["centroid_lat"])
            q_lon = float(c_rec["centroid_lon"])

    if q_lat is None or q_lon is None:
        return False

    # Exact polygon boundary geodesic check when documented county geometries are available
    if doc_union_geom is not None and Point is not None and nearest_points is not None and not doc_union_geom.is_empty:
        pt = Point(q_lon, q_lat)
        near_geom_pt, _ = nearest_points(doc_union_geom, pt)
        dist_m = haversine_distance_meters(q_lat, q_lon, float(near_geom_pt.y), float(near_geom_pt.x))
        if dist_m <= buffer_m:
            return True

    # Centroid-to-centroid fallback for counties not in a partial/offline polygon layer
    for d_ck in doc_counties:
        c_rec = counties_dict.get(d_ck)
        if not c_rec:
            continue
        clat = c_rec.get("centroid_lat")
        clon = c_rec.get("centroid_lon")
        if clat is not None and clon is not None:
            if haversine_distance_meters(q_lat, q_lon, float(clat), float(clon)) <= buffer_m:
                return True

    return False


def filter_vouchers_by_usda_county_range(
    df: "pd.DataFrame",
    standardize_taxon_fn: Callable[[Any], Optional[str]],
    synonym_to_canonical: Optional[Mapping[str, str]] = None,
    county_buffer_km: float = 50.0,
    exempt_gold: bool = True,
    cache_dir: Optional[Path] = None,
    county_geojson_path: Optional[Path] = None,
    db_path: Optional[Union[str, Path]] = None,
    workspace_dir: Path = PROJECT_ROOT,
    allow_network: bool = True,
    force_refresh: bool = False,
    drop_out_of_range: bool = True,
    logger: Optional[logging.Logger] = None,
) -> Tuple["pd.DataFrame", Dict[str, int]]:
    """
    Evaluate harvested vouchers against USDA PLANTS county distributions (with 1-ring
    county adjacency and `county_buffer_km` geodesic buffer) and filter out non-exempt
    `out_of_range` misidentifications prior to Step 01 image downloading.

    Exemptions:
      - Verified nomenclatural type specimens (`type_rank >= 1`) are unconditionally retained.
      - `Tier_1_Gold` specialist determinations are retained when `exempt_gold=True`.
      - Taxa with zero mapped US counties in USDA PLANTS (`unmapped_taxon`) and non-US
        occurrences (`non_us_skip`) fail open so valid records are never dropped due to
        missing checklist coverage.
    """
    stats: Dict[str, int] = {
        "evaluated": 0,
        "in_range": 0,
        "buffered_adjacent": 0,
        "exempt_retained": 0,
        "unmapped_or_non_us": 0,
        "dropped_out_of_range": 0,
    }
    if df is None or df.empty:
        return df, stats

    out_df = df.copy()
    stats["evaluated"] = len(out_df)

    # Build reverse synonym map: canonical_packera -> sorted synonym strings
    canon_to_syns: Dict[str, List[str]] = {}
    for raw_low, canon_val in (synonym_to_canonical or {}).items():
        if canon_val:
            canon_to_syns.setdefault(canon_val, []).append(raw_low)

    # Resolve canonical Packera taxon per row
    resolved_taxa: List[str] = []
    for _, row in out_df.iterrows():
        std_tx: Optional[str] = None
        for col in ("species_standardized", "scientificName", "species_raw"):
            if col in row and row[col] is not None and not pd.isna(row[col]) and str(row[col]).strip():
                std_tx = standardize_taxon_fn(row[col])
                if std_tx:
                    break
        resolved_taxa.append(std_tx or "")

    unique_taxa = sorted({tx for tx in resolved_taxa if tx.startswith("Packera ")})
    if not unique_taxa:
        out_df["usda_symbol"] = ""
        out_df["usda_county_status"] = "unmapped_taxon"
        stats["unmapped_or_non_us"] = len(out_df)
        return out_df, stats

    taxon_ranges: Dict[str, Dict[str, Any]] = {}
    for tx in unique_taxa:
        taxon_ranges[tx] = fetch_or_load_usda_taxon_range(
            canonical_packera=tx,
            synonym_names=canon_to_syns.get(tx, ()),
            cache_dir=cache_dir,
            db_path=db_path,
            workspace_dir=workspace_dir,
            allow_network=allow_network,
            force_refresh=force_refresh,
            logger=logger,
        )

    has_any_mapped = any(
        bool(r["counties"] or r["state_only_states"])
        for r in taxon_ranges.values()
    )
    if not has_any_mapped:
        out_df["usda_symbol"] = ""
        out_df["usda_county_status"] = "unmapped_taxon"
        stats["unmapped_or_non_us"] = len(out_df)
        return out_df, stats

    spatial_graph = load_us_county_spatial_graph(
        county_path=county_geojson_path,
        allow_download=allow_network,
    )
    counties_dict = spatial_graph["counties"]
    neighbors_map: Mapping[str, Set[str]] = spatial_graph["neighbors"]

    # Precompute 1-ring neighbor sets and unioned Shapely geometries per mapped taxon
    taxon_expanded_neighbors: Dict[str, Set[str]] = {}
    taxon_doc_union_geom: Dict[str, Any] = {}
    for tx, r_info in taxon_ranges.items():
        doc_cks: Set[str] = r_info["counties"]
        adj_set: Set[str] = set()
        doc_geoms_tx: List[Any] = []
        for ck in doc_cks:
            adj_set.update(neighbors_map.get(ck, ()))
            c_rec = counties_dict.get(ck)
            if c_rec and c_rec.get("geom") is not None:
                doc_geoms_tx.append(c_rec["geom"])
        adj_set -= doc_cks
        taxon_expanded_neighbors[tx] = adj_set
        taxon_doc_union_geom[tx] = (
            unary_union(doc_geoms_tx) if (unary_union is not None and doc_geoms_tx) else None
        )

    keep_mask: List[bool] = []
    usda_symbols_col: List[str] = []
    usda_status_col: List[str] = []
    dropped_by_taxon: Dict[str, int] = {}

    for idx_pos, (_, row) in enumerate(out_df.iterrows()):
        tx = resolved_taxa[idx_pos]
        r_info = taxon_ranges.get(tx)
        sym_str = ";".join(r_info["symbols"]) if (r_info and r_info["symbols"]) else ""
        usda_symbols_col.append(sym_str)

        if not r_info or (not r_info["counties"] and not r_info["state_only_states"]):
            usda_status_col.append("unmapped_taxon")
            stats["unmapped_or_non_us"] += 1
            keep_mask.append(True)
            continue

        raw_ts = row.get("type_status")
        if raw_ts is None or pd.isna(raw_ts) or not str(raw_ts).strip():
            raw_ts = row.get("typeStatus")
        type_rank, _ = classify_nomenclatural_type_rank(raw_ts)
        is_type = type_rank >= 1
        det_tier = str(row.get("determiner_tier") or "").strip()
        is_exempt = is_type or (exempt_gold and det_tier == "Tier_1_Gold")

        lat_v: Optional[float] = None
        lon_v: Optional[float] = None
        for c_lat in ("decimalLatitude", "latitude"):
            if c_lat in row and row[c_lat] is not None and not pd.isna(row[c_lat]):
                try:
                    fv = float(row[c_lat])
                    if math.isfinite(fv):
                        lat_v = fv
                        break
                except (TypeError, ValueError):
                    pass
        for c_lon in ("decimalLongitude", "longitude"):
            if c_lon in row and row[c_lon] is not None and not pd.isna(row[c_lon]):
                try:
                    fv = float(row[c_lon])
                    if math.isfinite(fv):
                        lon_v = fv
                        break
                except (TypeError, ValueError):
                    pass

        st_raw = None if pd.isna(row.get("stateProvince")) else str(row.get("stateProvince") or "")
        co_raw = None if pd.isna(row.get("county")) else str(row.get("county") or "")
        st_norm, co_norm, is_non_us = resolve_voucher_us_county(
            lat=lat_v,
            lon=lon_v,
            state_province=st_raw,
            county=co_raw,
            spatial_graph=spatial_graph,
        )

        if is_non_us or (st_norm is None and lat_v is None and lon_v is None):
            usda_status_col.append("non_us_or_unresolved")
            stats["unmapped_or_non_us"] += 1
            keep_mask.append(True)
            continue

        doc_cks = r_info["counties"]
        st_only = r_info["state_only_states"]
        doc_states = r_info["documented_states"]
        voucher_ck = f"{st_norm}|{co_norm}" if (st_norm and co_norm) else None

        if (st_norm and st_norm in st_only) or (voucher_ck and voucher_ck in doc_cks):
            usda_status_col.append("in_range")
            stats["in_range"] += 1
            keep_mask.append(True)
            continue

        # Ungeoreferenced state-only label without county name: retain if state has documented presence
        if st_norm and co_norm is None and lat_v is None and st_norm in doc_states:
            usda_status_col.append("in_range")
            stats["in_range"] += 1
            keep_mask.append(True)
            continue

        is_adj = bool(voucher_ck and voucher_ck in taxon_expanded_neighbors.get(tx, set()))
        if not is_adj and county_buffer_km > 0.0:
            is_adj = _is_within_geodesic_buffer(
                lat_f=lat_v,
                lon_f=lon_v,
                voucher_ck=voucher_ck,
                doc_counties=doc_cks,
                doc_union_geom=taxon_doc_union_geom.get(tx),
                spatial_graph=spatial_graph,
                buffer_km=county_buffer_km,
            )

        if is_adj:
            usda_status_col.append("buffered_adjacent")
            stats["buffered_adjacent"] += 1
            keep_mask.append(True)
            continue

        usda_status_col.append("out_of_range")
        if is_exempt:
            stats["exempt_retained"] += 1
            keep_mask.append(True)
        elif not drop_out_of_range:
            keep_mask.append(True)
        else:
            stats["dropped_out_of_range"] += 1
            dropped_by_taxon[tx] = dropped_by_taxon.get(tx, 0) + 1
            keep_mask.append(False)

    out_df["usda_symbol"] = usda_symbols_col
    out_df["usda_county_status"] = usda_status_col
    filtered_df = out_df.loc[keep_mask].reset_index(drop=True)

    if logger:
        logger.info(
            f"USDA PLANTS County Range Filter (buffer={county_buffer_km:.1f} km, exempt_gold={exempt_gold}): "
            f"Evaluated {stats['evaluated']} -> Retained {len(filtered_df)} "
            f"[in_range={stats['in_range']}, buffered_adjacent={stats['buffered_adjacent']}, "
            f"exempt={stats['exempt_retained']}, unmapped/non_us={stats['unmapped_or_non_us']}] | "
            f"Dropped out-of-range={stats['dropped_out_of_range']}."
        )
        if dropped_by_taxon:
            top_drops = ", ".join(f"{k}: {v}" for k, v in sorted(dropped_by_taxon.items()))
            logger.info(f"  Out-of-range vouchers dropped by taxon: {top_drops}")

    return filtered_df, stats
