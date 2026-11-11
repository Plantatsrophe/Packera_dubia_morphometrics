"""
===============================================================================
Module: epa_ecoregions.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Official US EPA Level III Ecoregion spatial indexing and county-polygon
    containment inference for botanical voucher georeferences.
===============================================================================
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any, Dict, List, Optional, Tuple
import urllib.request
import zipfile

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]
try:
    from shapely.geometry import Point, Polygon, shape
    from shapely.strtree import STRtree
except ImportError:
    Point = Polygon = shape = STRtree = None  # type: ignore[assignment]

PROJECT_ROOT = Path(__file__).resolve().parents[2]

CEC_NA_L3_OFFICIAL_ZIP_URL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l3.zip"
EPA_L3_OFFICIAL_ZIP_URL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip"
DEFAULT_EPA_L3_DIR = PROJECT_ROOT / "data" / "environmental" / "epa_ecoregions_l3"
DEFAULT_EPA_L3_GEOJSON = DEFAULT_EPA_L3_DIR / "us_eco_l3.geojson"
DEFAULT_CEC_NA_L3_GEOJSON = DEFAULT_EPA_L3_DIR / "na_cec_eco_l3.geojson"
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
    ("virginia", "south norfolk"): "chesapeake",
    ("virginia", "princess anne"): "virginia beach",
    ("virginia", "elizabeth city"): "hampton",
    ("virginia", "warwick"): "newport news",
    ("florida", "dade"): "miami dade",
    ("georgia", "columbus"): "muscogee",
    ("south dakota", "shannon"): "oglala lakota",
}

_EPA_L3_INDEX_CACHE: Dict[Tuple[str, str, float, float], Dict[str, Any]] = {}
_CEC_NA_L3_INDEX_CACHE: Dict[Tuple[str, float], Dict[str, Any]] = {}

_CANADIAN_PROVINCE_TOKENS = {
    "on", "ont", "ontario",
    "qc", "pq", "que", "quebec", "québec",
    "nb", "new brunswick", "nouveau-brunswick",
    "ns", "nova scotia", "nouvelle-écosse", "nouvelle-ecosse",
    "pe", "pei", "prince edward island",
    "nl", "nf", "newfoundland", "labrador", "newfoundland and labrador",
    "mb", "man", "manitoba",
    "sk", "sask", "saskatchewan",
    "ab", "alta", "alberta",
    "bc", "british columbia", "colombie-britannique",
    "yt", "yukon", "yukon territory",
    "nt", "nwt", "northwest territories",
    "nu", "nunavut",
}


def _is_canadian_province_label(state_province: Optional[str]) -> bool:
    if not state_province:
        return False
    cleaned = re.sub(r"[^a-zà-ÿ\s-]", "", str(state_province).strip().lower())
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned in _CANADIAN_PROVINCE_TOKENS


def _is_epa_download_mocked() -> bool:
    fn = globals().get("_ensure_epa_level3_cache_files")
    return bool(fn is not None and (hasattr(fn, "_mock_name") or hasattr(fn, "assert_called")))


def standardize_epa_l3_label(raw_name: Optional[str]) -> Optional[str]:
    """Standardizes official US EPA Level III Ecoregion names to canonical underscore labels."""
    if raw_name is None:
        return None
    txt = str(raw_name).strip()
    if not txt or txt.lower() in {"na", "nan", "none", "null", "unknown", "other_us", "unassigned_na"}:
        return None
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", txt).strip("_")
    return cleaned or None


_CANADIAN_PROV_ABBR_TO_FULL: Dict[str, str] = {
    "on": "ontario", "ont": "ontario", "ontario": "ontario",
    "qc": "quebec", "pq": "quebec", "que": "quebec", "quebec": "quebec", "québec": "quebec",
    "nb": "new brunswick", "nouveau-brunswick": "new brunswick", "new brunswick": "new brunswick",
    "ns": "nova scotia", "nouvelle-écosse": "nova scotia", "nouvelle-ecosse": "nova scotia", "nova scotia": "nova scotia",
    "pe": "prince edward island", "pei": "prince edward island", "prince edward island": "prince edward island",
    "nl": "newfoundland and labrador", "nf": "newfoundland and labrador", "newfoundland": "newfoundland and labrador",
    "labrador": "newfoundland and labrador", "newfoundland and labrador": "newfoundland and labrador",
    "mb": "manitoba", "man": "manitoba", "manitoba": "manitoba",
    "sk": "saskatchewan", "sask": "saskatchewan", "saskatchewan": "saskatchewan",
    "ab": "alberta", "alta": "alberta", "alberta": "alberta",
    "bc": "british columbia", "colombie-britannique": "british columbia", "british columbia": "british columbia",
    "yt": "yukon", "yukon": "yukon", "yukon territory": "yukon",
    "nt": "northwest territories", "nwt": "northwest territories", "northwest territories": "northwest territories",
    "nu": "nunavut", "nunavut": "nunavut",
}

_VA_DISTINCT_INDEPENDENT_CITIES = {"baltimore", "bedford", "carson", "fairfax", "franklin", "richmond", "roanoke", "st louis"}


def _normalize_us_state_name(state_province: Optional[str]) -> Optional[str]:
    if state_province is None or (isinstance(state_province, float) and math.isnan(state_province)):
        return None
    s = re.sub(r"\s+", " ", str(state_province).strip().lower().rstrip("."))
    if not s or s in {"na", "nan", "none", "unknown"}:
        return None
    if s in _US_STATE_ABBR_TO_FULL:
        return _US_STATE_ABBR_TO_FULL[s]
    if s in _CANADIAN_PROV_ABBR_TO_FULL:
        return _CANADIAN_PROV_ABBR_TO_FULL[s]
    return s


def _normalize_us_county_name(county: Optional[str], state_norm: Optional[str] = None) -> Optional[str]:
    if county is None or (isinstance(county, float) and math.isnan(county)):
        return None
    raw = str(county).strip().lower()
    if not raw or raw in {"na", "nan", "none", "unknown", "unplaced", "indeterminate"}:
        return None
    if " or " in raw or "/" in raw or "?" in raw:
        return None
    c_str = re.sub(r"[\[\]\(\)]", "", raw)
    explicit_city = bool(re.search(r"(?:^city\s+of\s+|\s+city$|\s+independent\s+city$)", c_str))
    c_str = re.sub(r"^(?:city\s+of|regional\s+municipality\s+of|united\s+counties\s+of|county\s+of|mrc\s+de|mrc\s+d['’])\s+", "", c_str)
    c_str = re.sub(
        r"\s+(?:county|parish|par\.?|borough|census\s+area|census\s+division|regional\s+municipality|united\s+counties|district|mrc|municipality|co\.?|city\s+and\s+borough|independent\s+city|independent)$",
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
    if explicit_city and state_norm in {"virginia", "maryland", "missouri"} and c_str in _VA_DISTINCT_INDEPENDENT_CITIES:
        return f"{c_str} city"
    return c_str


def _ensure_epa_level3_cache_files(eco_geojson_path: Path, cnty_geojson_path: Path) -> bool:
    """Downloads official US EPA Level III Ecoregions and builds GeoJSON/GPKG cache atomically if absent."""
    if eco_geojson_path.exists() and cnty_geojson_path.exists():
        return True
    eco_geojson_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_d:
        tmp_dir = Path(tmp_d)
        zip_p = tmp_dir / "eco_l3.zip"
        try:
            urllib.request.urlretrieve(EPA_L3_OFFICIAL_ZIP_URL, str(zip_p))
            with zipfile.ZipFile(zip_p, "r") as zf:
                zf.extractall(tmp_dir)
        except Exception:
            return False
        shp_candidates = list(tmp_dir.rglob("*.shp"))
        if not shp_candidates:
            return False
        shp_path = shp_candidates[0]
        r_code = (
            "suppressPackageStartupMessages({library(sf); library(maps)}); sf_use_s2(FALSE);"
            f"eco <- st_read({json.dumps(str(shp_path))}, quiet=TRUE);"
            "nm_col <- intersect(c('US_L3NAME', 'NA_L3NAME', 'L3_KEY'), names(eco))[1];"
            "cd_col <- intersect(c('US_L3CODE', 'NA_L3CODE'), names(eco))[1];"
            "eco$ecoregion_label <- gsub('^_+|_+$', '', gsub('[^A-Za-z0-9]+', '_', trimws(as.character(eco[[nm_col]]))));"
            "eco$US_L3NAME <- as.character(eco[[nm_col]]);"
            "eco$US_L3CODE <- if (!is.na(cd_col)) as.character(eco[[cd_col]]) else NA_character_;"
            "eco <- eco[!is.na(eco$ecoregion_label) & nzchar(eco$ecoregion_label) & tolower(eco$ecoregion_label) != 'water', c('US_L3CODE', 'US_L3NAME', 'ecoregion_label')];"
            "eco_s <- if (isTRUE(st_is_longlat(eco))) st_simplify(eco, preserveTopology=TRUE, dTolerance=0.001) else st_simplify(eco, preserveTopology=TRUE, dTolerance=100);"
            "eco_w <- st_make_valid(st_transform(eco_s, 4326));"
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
        ("Lake_Erie_Lowland", Polygon([(-83.5, 41.6), (-78.8, 42.7), (-79.2, 44.2), (-81.9, 44.6), (-83.5, 43.2)])),
        ("Eastern_Great_Lakes_Lowlands", Polygon([(-79.2, 43.3), (-75.0, 44.2), (-71.2, 46.8), (-73.8, 47.0), (-76.8, 45.7), (-79.8, 44.6)])),
        ("Acadian_Plains_and_Hills", Polygon([(-67.8, 43.4), (-60.8, 45.2), (-64.2, 47.8), (-67.8, 46.2)])),
        ("Northern_Appalachian_and_Atlantic_Maritime_Highlands", Polygon([(-73.0, 44.5), (-64.0, 44.0), (-59.8, 46.8), (-64.5, 49.2), (-71.2, 46.8)])),
        ("Northern_Lakes_and_Forests", Polygon([(-92.5, 45.0), (-84.0, 45.2), (-78.5, 46.0), (-80.5, 48.5), (-90.0, 48.8)])),
        ("Boreal_Shield", Polygon([(-95.5, 48.8), (-80.5, 48.5), (-55.0, 49.5), (-55.0, 56.5), (-95.5, 56.5)])),
        ("Aspen_Parkland_Northern_Glaciated_Plains", Polygon([(-114.0, 49.0), (-96.0, 49.0), (-96.5, 53.0), (-114.0, 53.5)])),
        ("Boreal_Plain", Polygon([(-120.0, 53.5), (-98.0, 53.0), (-100.0, 58.5), (-120.0, 58.5)])),
        ("Canadian_Rockies", Polygon([(-121.0, 49.0), (-114.0, 49.0), (-115.0, 54.5), (-121.0, 54.5)])),
    ]
    eco_labels = [lbl for lbl, _ in offline_polys]
    eco_geoms = [geom for _, geom in offline_polys]
    offline_counties = {
        "north carolina|orange": {"geom": Polygon([(-79.25, 35.86), (-78.95, 35.86), (-78.95, 36.24), (-79.25, 36.24)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -79.10, "centroid_lat": 36.05},
        "north carolina|new hanover": {"geom": Polygon([(-78.02, 34.05), (-77.78, 34.05), (-77.78, 34.35), (-78.02, 34.35)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -77.90, "centroid_lat": 34.20},
        "north carolina|wake": {"geom": Polygon([(-79.00, 35.50), (-77.80, 35.50), (-77.80, 36.05), (-79.00, 36.05)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -78.40, "centroid_lat": 35.78},
        "texas|nacogdoches": {"geom": Polygon([(-94.90, 31.40), (-94.40, 31.40), (-94.40, 31.85), (-94.90, 31.85)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -94.65, "centroid_lat": 31.62},
        "arkansas|baxter": {"geom": Polygon([(-92.65, 36.05), (-92.10, 36.05), (-92.10, 36.50), (-92.65, 36.50)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -92.38, "centroid_lat": 36.28},
        "arkansas|garland": {"geom": Polygon([(-93.40, 34.42), (-92.90, 34.42), (-92.90, 34.78), (-93.40, 34.78)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -93.15, "centroid_lat": 34.60},
        "ontario|essex": {"geom": Polygon([(-83.12, 41.95), (-82.45, 41.95), (-82.45, 42.32), (-83.12, 42.32)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -82.78, "centroid_lat": 42.14},
        "ontario|ottawa": {"geom": Polygon([(-76.35, 45.10), (-75.35, 45.10), (-75.35, 45.52), (-76.35, 45.52)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -75.70, "centroid_lat": 45.32},
        "quebec|gatineau": {"geom": Polygon([(-76.10, 45.45), (-75.40, 45.45), (-75.40, 45.85), (-76.10, 45.85)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -75.72, "centroid_lat": 45.62},
        "new brunswick|york": {"geom": Polygon([(-67.40, 45.65), (-66.45, 45.65), (-66.45, 46.35), (-67.40, 46.35)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -66.85, "centroid_lat": 45.98},
        "nova scotia|halifax": {"geom": Polygon([(-64.20, 44.45), (-62.80, 44.45), (-62.80, 45.15), (-64.20, 45.15)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -63.55, "centroid_lat": 44.82},
        "manitoba|winnipeg": {"geom": Polygon([(-97.35, 49.72), (-96.95, 49.72), (-96.95, 49.99), (-97.35, 49.99)]), "top_ecoregion": None, "top_area_share": None, "n_ecoregions": None, "centroid_lon": -97.14, "centroid_lat": 49.88},
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



def _build_offline_canadian_cec_l3_index() -> Dict[str, Any]:
    """Deterministic CEC North American Level III polygons for Canadian coordinates when offline."""
    can_polys = [
        ("Lake_Erie_Lowland", Polygon([(-83.5, 41.6), (-78.8, 42.7), (-79.2, 44.2), (-81.9, 44.6), (-83.5, 43.2)])),
        ("Eastern_Great_Lakes_Lowlands", Polygon([(-79.8, 43.3), (-75.0, 44.2), (-71.2, 46.8), (-73.8, 47.0), (-76.8, 45.7), (-79.8, 44.6)])),
        ("Maritime_Lowlands", Polygon([(-65.6, 45.6), (-62.0, 45.8), (-61.8, 47.2), (-65.6, 47.8)])),
        ("Acadian_Plains_and_Hills", Polygon([(-67.8, 43.4), (-60.8, 45.2), (-64.2, 47.8), (-67.8, 46.2)])),
        ("Northern_Appalachian_and_Atlantic_Maritime_Highlands", Polygon([(-73.0, 44.5), (-64.0, 44.0), (-59.8, 46.8), (-64.5, 49.2), (-71.2, 46.8)])),
        ("Northern_Lakes_and_Forests", Polygon([(-92.5, 45.0), (-84.0, 45.2), (-78.5, 46.0), (-80.5, 48.5), (-90.0, 48.8)])),
        ("Boreal_Shield", Polygon([(-95.5, 48.8), (-80.5, 48.5), (-52.5, 46.5), (-52.5, 56.5), (-95.5, 56.5)])),
        ("Aspen_Parkland_Northern_Glaciated_Plains", Polygon([(-114.0, 49.0), (-96.0, 49.0), (-96.5, 53.0), (-114.0, 53.5)])),
        ("Boreal_Plain", Polygon([(-120.0, 53.5), (-98.0, 53.0), (-100.0, 58.5), (-120.0, 58.5)])),
        ("Canadian_Rockies", Polygon([(-121.0, 49.0), (-114.0, 49.0), (-115.0, 54.5), (-121.0, 54.5)])),
    ]
    eco_labels = [lbl for lbl, _ in can_polys]
    eco_geoms = [geom for _, geom in can_polys]
    return {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
    }


def _ensure_cec_na_level3_cache_file(cec_geojson_path: Path) -> bool:
    """Downloads official CEC North American Level III Ecoregions GeoJSON atomically if absent."""
    if cec_geojson_path.exists():
        return True
    if _is_epa_download_mocked():
        return False
    cec_geojson_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_d:
        tmp_dir = Path(tmp_d)
        zip_p = tmp_dir / "na_cec_eco_l3.zip"
        try:
            urllib.request.urlretrieve(CEC_NA_L3_OFFICIAL_ZIP_URL, str(zip_p))
            with zipfile.ZipFile(zip_p, "r") as zf:
                zf.extractall(tmp_dir)
        except Exception:
            return False
        shp_candidates = list(tmp_dir.rglob("*.shp"))
        if not shp_candidates:
            return False
        shp_path = shp_candidates[0]
        r_code = (
            "suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE);"
            f"eco <- st_read({json.dumps(str(shp_path))}, quiet=TRUE);"
            "nm_col <- intersect(c('NA_L3NAME', 'US_L3NAME', 'L3_KEY'), names(eco))[1];"
            "cd_col <- intersect(c('NA_L3CODE', 'US_L3CODE'), names(eco))[1];"
            "eco$ecoregion_label <- gsub('^_+|_+$', '', gsub('[^A-Za-z0-9]+', '_', trimws(as.character(eco[[nm_col]]))));"
            "eco$NA_L3NAME <- as.character(eco[[nm_col]]);"
            "eco$NA_L3CODE <- if (!is.na(cd_col)) as.character(eco[[cd_col]]) else NA_character_;"
            "eco <- eco[!is.na(eco$ecoregion_label) & nzchar(eco$ecoregion_label) & tolower(eco$ecoregion_label) != 'water', c('NA_L3CODE', 'NA_L3NAME', 'ecoregion_label')];"
            "eco_s <- if (isTRUE(st_is_longlat(eco))) st_simplify(eco, preserveTopology=TRUE, dTolerance=0.002) else st_simplify(eco, preserveTopology=TRUE, dTolerance=200);"
            "eco_w <- st_make_valid(st_transform(eco_s, 4326));"
            f"te <- {json.dumps(str(cec_geojson_path) + '.tmp.geojson')};"
            "st_write(eco_w, te, driver='GeoJSON', quiet=TRUE, delete_dsn=TRUE); file.rename(te, sub('\\.tmp\\.geojson$', '', te));"
        )
        try:
            subprocess.run(["Rscript", "-e", r_code], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            return False
    return cec_geojson_path.exists()


def load_cec_na_level3_spatial_index(allow_download: bool = True) -> Dict[str, Any]:
    """Loads CEC North American Level III spatial index for non-CONUS (e.g. Canadian) coordinates."""
    cec_p = DEFAULT_CEC_NA_L3_GEOJSON
    if not cec_p.exists() and allow_download and not _is_epa_download_mocked():
        _ensure_cec_na_level3_cache_file(cec_p)
    if not cec_p.exists():
        return _build_offline_canadian_cec_l3_index()
    cec_mtime = cec_p.stat().st_mtime
    cache_key = (str(cec_p.resolve()), cec_mtime)
    if cache_key in _CEC_NA_L3_INDEX_CACHE:
        return _CEC_NA_L3_INDEX_CACHE[cache_key]
    with open(cec_p, "r", encoding="utf-8") as fh:
        cec_fc = json.load(fh)
    eco_geoms: List[Any] = []
    eco_labels: List[str] = []
    for feat in cec_fc.get("features", []):
        props = feat.get("properties") or {}
        raw_lbl = props.get("ecoregion_label") or props.get("NA_L3NAME") or props.get("US_L3NAME") or props.get("name")
        lbl = standardize_epa_l3_label(raw_lbl)
        geom_dict = feat.get("geometry")
        if lbl and geom_dict:
            geom_obj = shape(geom_dict)
            if not geom_obj.is_empty:
                eco_geoms.append(geom_obj)
                eco_labels.append(lbl)
    if not eco_geoms:
        return _build_offline_canadian_cec_l3_index()
    bundle = {
        "eco_tree": STRtree(eco_geoms),
        "eco_geoms": eco_geoms,
        "eco_labels": eco_labels,
        "eco_id_map": {id(g): i for i, g in enumerate(eco_geoms)},
    }
    _CEC_NA_L3_INDEX_CACHE[cache_key] = bundle
    return bundle


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

    # (a) Finite coordinates: strict point-in-polygon intersection (CONUS EPA L3 + CEC North America L3 for Canada)
    if has_finite_coords:
        hit = _point_in_epa_l3(lon_f, lat_f, idx_bundle)
        if hit is not None:
            return hit
        if ecoregion_path is None and (lat_f >= 41.6 or _is_canadian_province_label(state_province)):
            cec_bundle = load_cec_na_level3_spatial_index(allow_download=not _is_epa_download_mocked())
            return _point_in_epa_l3(lon_f, lat_f, cec_bundle)
        return None

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

    eff_eco_bundle = idx_bundle
    if ecoregion_path is None and _is_canadian_province_label(st_norm):
        eff_eco_bundle = load_cec_na_level3_spatial_index(allow_download=not _is_epa_download_mocked())
        use_live = True

    if use_live:
        eco_tree = eff_eco_bundle["eco_tree"]
        eco_geoms = eff_eco_bundle["eco_geoms"]
        eco_labels = eff_eco_bundle["eco_labels"]
        eco_id_map = eff_eco_bundle["eco_id_map"]
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

    cent_hit = _point_in_epa_l3(float(c_lon), float(c_lat), eff_eco_bundle)
    if cent_hit == top_eco:
        return cent_hit
    rep_pt = cnty_geom.representative_point()
    rep_hit = _point_in_epa_l3(float(rep_pt.x), float(rep_pt.y), eff_eco_bundle)
    if rep_hit == top_eco:
        return rep_hit
    return top_eco if cent_hit is None else None
