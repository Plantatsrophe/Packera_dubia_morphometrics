"""
===============================================================================
Module: coord_hygiene_utils.py
Project: Plantmatics Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Stdlib-only coordinate hygiene quarantine classifier flagging herbarium
    institution rooftops, US county centroids, and coarse 0.01-degree gridded
    coordinates without mutating raw Darwin Core latitude/longitude fields.
===============================================================================
"""

from __future__ import annotations

import math
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple
import unicodedata

_EARTH_RADIUS_M: float = 6371008.8

HERBARIUM_AND_BOTANICAL_INSTITUTION_COORDS: Dict[str, Tuple[float, float]] = {
    "NCU": (35.9108, -79.0505), "DUKE": (36.0014, -78.9382), "NCSC": (35.7877, -78.6703),
    "NY": (40.8623, -73.8803), "GH": (42.3785, -71.1160), "US": (38.8889, -77.0260),
    "MO": (38.6128, -90.2594), "F": (41.8663, -87.6170), "PH": (39.9571, -75.1714),
    "BRIT": (32.7397, -97.3631), "TEX": (30.2850, -97.7340), "WIS": (43.0742, -89.4012),
    "MIN": (44.9845, -93.1832), "TENN": (35.9544, -83.9295), "GA": (33.9480, -83.3773),
    "WILLI": (37.2707, -76.7153), "LSU": (30.4133, -91.1800), "CM": (40.4433, -79.9496),
    "MU": (39.5097, -84.7345), "CSCN": (42.8186, -103.0005), "DAO": (45.3920, -75.7160),
    "CAN": (45.3920, -75.7160), "MT": (45.5569, -73.5575), "TRT": (43.6677, -79.3948),
}

_WILD_HABITAT_ESCAPE_RE = re.compile(
    r"\b(?:\d+(?:\.\d+)?\s*(?:mi|miles?|km|m)\.?\s+(?:[NSEW]{1,3}|north|south|east|west)\b"
    r"|outcrop|flatrock|serpentine|ultramafic|diabase|gabbro|limestone|dolomite|bluff"
    r"|cliff|gorge|seep|fen|bog|swamp|creek|river|branch|ridge|summit|gap"
    r"|national\s+forest|state\s+park|wildlife\s+management|nature\s+preserve"
    r"|biological\s+reserve|mason\s+farm)\b",
    re.IGNORECASE,
)

_VERBATIM_COORD_DIGITS_RE = re.compile(
    r"(?:\b\d{2,3}\.\d{3,}\b|\b\d{1,3}\s*[°º]\s*\d{1,2}|\b(?:lat|lon|gps)\b[:\s]*[-+]?\d{2})",
    re.IGNORECASE,
)

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

_CA_PROV_ABBR_TO_FULL: Dict[str, str] = {
    "on": "ontario", "ont": "ontario", "qc": "quebec", "pq": "quebec", "que": "quebec",
    "nb": "new brunswick", "ns": "nova scotia", "pe": "prince edward island",
    "pei": "prince edward island", "nl": "newfoundland and labrador",
    "nf": "newfoundland and labrador", "mb": "manitoba", "man": "manitoba",
    "sk": "saskatchewan", "sask": "saskatchewan", "ab": "alberta", "alta": "alberta",
    "bc": "british columbia", "yt": "yukon", "nt": "northwest territories",
    "nwt": "northwest territories", "nu": "nunavut",
}

_ADMIN_LOCALITY_TOKENS: Set[str] = {
    "county", "parish", "borough", "co", "usa", "us", "united", "states", "canada",
}

_STATE_PROV_TOKENS: Set[str] = set()
for _abbr, _full in (*_US_STATE_ABBR_TO_FULL.items(), *_CA_PROV_ABBR_TO_FULL.items()):
    _STATE_PROV_TOKENS.add(_abbr)
    _STATE_PROV_TOKENS.update(re.findall(r"[a-z]+", _full))


def haversine_distance_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Compute great-circle distance in meters using mean Earth radius 6,371,008.8 m."""
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi, dlam = math.radians(lat2 - lat1), math.radians(lon2 - lon1)
    a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2.0) ** 2
    return _EARTH_RADIUS_M * 2.0 * math.atan2(math.sqrt(a), math.sqrt(max(0.0, 1.0 - a)))


def _clean_text(val: Any) -> str:
    if val is None or isinstance(val, bool) or (isinstance(val, float) and math.isnan(val)):
        return ""
    s = str(val).strip()
    return "" if s.lower() in {"", "na", "nan", "none", "null"} else s


def _parse_optional_finite_float(val: Any) -> Optional[float]:
    if val is None or isinstance(val, bool):
        return None
    try:
        f = float(val)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _fold_ascii_lower(text: str) -> str:
    folded = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()
    return re.sub(r"\b(?:[a-z]\.){2,}", lambda m: m.group(0).replace(".", ""), folded)


def _normalize_state(state_prov: Any) -> str:
    s = re.sub(r"\s+", " ", _fold_ascii_lower(_clean_text(state_prov)).rstrip("."))
    if not s or s in {"unknown", "unspecified", "unrecorded"}:
        return ""
    return _US_STATE_ABBR_TO_FULL.get(s) or _CA_PROV_ABBR_TO_FULL.get(s) or s


def _normalize_county(county: Any) -> str:
    raw = _fold_ascii_lower(_clean_text(county))
    if not raw or raw in {"unknown", "unplaced", "indeterminate"}:
        return ""
    s = re.sub(r"[\[\]\(\)]", "", raw).strip()
    s = re.sub(r"^(?:county\s+of|parish\s+of|city\s+of)\s+", "", s)
    s = re.sub(r"\s+(?:county|parish|par\.?|borough|census\s+area|municipality|co\.?)$", "", s)
    s = re.sub(r"\s+city$", "", s)
    s = re.sub(r"^(?:st\.|saint)\s+", "st ", s)
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def _extract_county_centroid(rec: Mapping[str, Any]) -> Optional[Tuple[float, float]]:
    clat = _parse_optional_finite_float(rec.get("centroid_lat"))
    clon = _parse_optional_finite_float(rec.get("centroid_lon"))
    if clat is not None and clon is not None:
        return clat, clon
    geom = rec.get("geom")
    if geom is not None and hasattr(geom, "centroid"):
        try:
            c = geom.centroid
            cy = _parse_optional_finite_float(getattr(c, "y", None))
            cx = _parse_optional_finite_float(getattr(c, "x", None))
            if cy is not None and cx is not None:
                return cy, cx
        except Exception:
            return None
    return None


def locality_only_county_or_state(locality: Any, county: Any = "", state_prov: Any = "") -> bool:
    """Return True if all alphabetic tokens in `locality` only name county/state/country terms."""
    loc_str = _clean_text(locality)
    if not loc_str:
        return False
    loc_tokens = re.findall(r"[a-z]+", _fold_ascii_lower(loc_str))
    if not loc_tokens:
        return False
    county_words = set(re.findall(r"[a-z]+", _fold_ascii_lower(_clean_text(county))))
    state_words = set(re.findall(r"[a-z]+", _fold_ascii_lower(_clean_text(state_prov))))
    allowed = _ADMIN_LOCALITY_TOKENS | _STATE_PROV_TOKENS | county_words | state_words
    return all(tok in allowed for tok in loc_tokens)


def evaluate_coordinate_centroid_quarantine(
    lat: Any,
    lon: Any,
    locality: str = "",
    habitat: str = "",
    county: str = "",
    state_prov: str = "",
    uncertainty_m: Any = None,
    issues: Optional[Iterable[str]] = None,
    is_type_specimen: bool = False,
    us_spatial_graph: Optional[Mapping[str, Any]] = None,
    hit_county_keys: Optional[Iterable[str]] = None,
    inst_radius_m: float = 250.0,
    county_centroid_radius_m: float = 750.0,
) -> Tuple[bool, str]:
    """Classify non-type coordinates matching herbarium rooftops, county centroids, or coarse grids."""
    if is_type_specimen:
        return False, ""
    lat_f = _parse_optional_finite_float(lat)
    lon_f = _parse_optional_finite_float(lon)
    if lat_f is None or lon_f is None:
        return False, ""

    loc_str, hab_str = _clean_text(locality), _clean_text(habitat)
    loc_hab = f"{loc_str} {hab_str}".strip()
    has_wild_escape = bool(_WILD_HABITAT_ESCAPE_RE.search(loc_hab))

    if not has_wild_escape:
        for ilat, ilon in HERBARIUM_AND_BOTANICAL_INSTITUTION_COORDS.values():
            if haversine_distance_m(lat_f, lon_f, ilat, ilon) <= float(inst_radius_m):
                return True, "herbarium_institution_rooftop"

    unc_val = _parse_optional_finite_float(uncertainty_m)
    unc_ok = (unc_val is None) or (unc_val >= 2000.0)
    counties = (us_spatial_graph or {}).get("counties") or {}

    if counties and unc_ok and not has_wild_escape and not _VERBATIM_COORD_DIGITS_RE.search(loc_hab):
        norm_state = _normalize_state(state_prov)
        norm_county = _normalize_county(county)
        hit_list = [str(k) for k in (hit_county_keys or ()) if k and str(k) in counties]
        ns_map = (us_spatial_graph or {}).get("counties_nospace") or {}

        if not norm_county and not hit_list:
            cand_iter = (
                rec for ck, rec in counties.items()
                if not (norm_state and "|" in ck and ck.split("|", 1)[0] != norm_state)
            )
            bbox_guard = True
        else:
            cand_keys = list(dict.fromkeys(hit_list))
            if norm_state and norm_county:
                exact_ck = f"{norm_state}|{norm_county}"
                mapped_ck = exact_ck if exact_ck in counties else ns_map.get(f"{norm_state}|{norm_county.replace(' ', '')}")
                if mapped_ck and mapped_ck in counties and mapped_ck not in cand_keys:
                    cand_keys.append(mapped_ck)
            elif norm_county and not norm_state and not hit_list:
                cand_keys.extend(ck for ck in counties if ck.endswith(f"|{norm_county}") and ck not in cand_keys)
            cand_iter = (counties[ck] for ck in cand_keys if ck in counties)
            bbox_guard = False

        for rec in cand_iter:
            centroid = _extract_county_centroid(rec) if isinstance(rec, Mapping) else None
            if centroid is None:
                continue
            clat, clon = centroid
            if bbox_guard and (abs(lat_f - clat) > 0.02 or abs(lon_f - clon) > 0.02):
                continue
            if haversine_distance_m(lat_f, lon_f, clat, clon) <= float(county_centroid_radius_m):
                return True, "county_centroid_suspect"

    if (
        not has_wild_escape
        and abs(lat_f * 100.0 - round(lat_f * 100.0)) < 1e-6
        and abs(lon_f * 100.0 - round(lon_f * 100.0)) < 1e-6
    ):
        issue_set = {str(x).strip().upper() for x in (issues or ()) if x}
        if "COORDINATE_ROUNDED" in issue_set:
            return True, "coarse_grid_rounded"
        if unc_val is None and bool(loc_str) and locality_only_county_or_state(loc_str, county, state_prov):
            return True, "coarse_grid_rounded"

    return False, ""
