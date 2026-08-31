"""
===============================================================================
Module: harvester_metadata.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Taxonomic, geographic, and phenological metadata evaluation utilities for
    herbarium voucher processing. Provides Determiner Authority stratification
    (Tiers 1–3), circular harmonic phenology transformations, regional eco-
    geographic classifications, and western boundary exclusion filters.
===============================================================================
"""

from __future__ import annotations

import datetime
import math
import re
from typing import Any, Optional, Tuple

from scripts.core.config import (
    EXCLUDED_WESTERN_STATES,
    MAJOR_HERBARIA_CODES,
    SPECIALIST_PATTERNS,
    VALID_TYPE_STATUSES,
    WESTERN_LONGITUDE_THRESHOLD,
)


def sanitize_filename(name: str) -> str:
    """
    Sanitizes arbitrary strings into safe, valid filesystem filenames across OS platforms.

    Args:
        name: Raw identifier or catalog number string.

    Returns:
        str: Cleaned alphanumeric filename string safe for Windows/Linux filesystems.
    """
    clean = re.sub(r'[\\/*?:"<>|\s]+', '_', str(name).strip())
    clean = clean.strip("._")
    return clean if clean else "voucher_unknown"


def parse_determiner_tier(
    type_status_raw: Optional[str],
    identified_by_raw: Optional[str],
    recorded_by_raw: Optional[str],
    history_raw: Optional[str],
    institution_code_raw: Optional[str],
    locality_raw: Optional[str],
    habitat_raw: Optional[str]
) -> Tuple[str, str, str]:
    """
    Evaluates Darwin Core fields to assign a taxonomic determination credibility tier.

    Tiers:
      - Tier_1_Gold: Primary/secondary nomenclatural types or verified specialist annotations
      - Tier_2_Silver: Determinations from major research herbaria with complete ecological locality
      - Tier_3_Bronze: Unverified candidate collections, general collectors, or missing determiner

    Args:
        type_status_raw: Darwin Core typeStatus.
        identified_by_raw: Darwin Core identifiedBy.
        recorded_by_raw: Darwin Core recordedBy.
        history_raw: Darwin Core verbatimIdentificationHistory.
        institution_code_raw: Darwin Core institutionCode.
        locality_raw: Darwin Core locality / verbatimLocality.
        habitat_raw: Darwin Core habitat.

    Returns:
        Tuple[str, str, str]: (determiner_tier, type_status_clean, determiner_raw_combined)
    """
    type_status = str(type_status_raw).strip() if type_status_raw is not None else ""
    identified_by = str(identified_by_raw).strip() if identified_by_raw is not None else ""
    recorded_by = str(recorded_by_raw).strip() if recorded_by_raw is not None else ""
    history = str(history_raw).strip() if history_raw is not None else ""
    institution = str(institution_code_raw).strip().upper() if institution_code_raw is not None else ""
    locality = str(locality_raw).strip() if locality_raw is not None else ""
    habitat = str(habitat_raw).strip() if habitat_raw is not None else ""

    determiner_raw = identified_by if identified_by else (history if history else recorded_by)
    combined_audit_text = f"{type_status} | {identified_by} | {history} | {recorded_by}"

    type_status_clean = "None"
    is_type = False
    if type_status and type_status.upper() not in {"NONE", "NOT A TYPE", "NOTATYPE", "UNSPECIFIED", "NULL"}:
        for valid_type in VALID_TYPE_STATUSES:
            if re.search(rf"\b{valid_type}\b", type_status, re.IGNORECASE):
                type_status_clean = valid_type.title()
                is_type = True
                break

    # 1. Tier 1 (Gold Standard): Type specimen OR specialist determination
    is_specialist = False
    for pattern in SPECIALIST_PATTERNS:
        if re.search(pattern, combined_audit_text, re.IGNORECASE):
            is_specialist = True
            break

    if is_type or is_specialist:
        return "Tier_1_Gold", type_status_clean, determiner_raw

    # 2. Tier 2 (Silver Standard): Major herbarium with rich locality/habitat metadata and non-empty determiner
    is_major_herbarium = institution in MAJOR_HERBARIA_CODES
    has_rich_locality = len(locality) > 10 or len(habitat) > 5
    has_determiner = bool(identified_by and identified_by.lower() not in {"unknown", "anonymous", "none", "null"})

    if is_major_herbarium and has_rich_locality and has_determiner:
        return "Tier_2_Silver", type_status_clean, determiner_raw

    # 3. Tier 3 (Bronze Standard): General collections, unverified determinations, or blank determiner
    return "Tier_3_Bronze", type_status_clean, determiner_raw


def calculate_circular_phenology(year: Any, month: Any, day: Any) -> Optional[Tuple[int, float, float]]:
    """
    Computes circular harmonic phenological coordinates (doy, sin, cos) from collection date.

    Transforms day-of-year into continuous periodic trigonometric functions on the unit circle
    to prevent boundary discontinuities between December 31 (DOY 365) and January 1 (DOY 1).

    Args:
        year: Collection year integer.
        month: Collection month integer (1-12).
        day: Collection day integer (1-31).

    Returns:
        Optional[Tuple[int, float, float]]: (doy, pheno_sin, pheno_cos) or None if invalid date.
    """
    try:
        y = int(year)
        m = int(month)
        d = int(day)
        collection_date = datetime.date(y, m, d)
        doy = collection_date.timetuple().tm_yday

        theta = 2.0 * math.pi * (doy / 365.25)
        pheno_sin = round(math.sin(theta), 6)
        pheno_cos = round(math.cos(theta), 6)

        return doy, pheno_sin, pheno_cos
    except (ValueError, TypeError, OverflowError):
        return None


def is_excluded_western_region(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None
) -> bool:
    """
    Determines whether a record originates from a US state or territory farther west than Texas and Oklahoma.

    Args:
        state_province: State or province name / abbreviation string (or None/NaN).
        lat: Optional decimal latitude float.
        lon: Optional decimal longitude float.

    Returns:
        bool: True if the record should be excluded as a western locality; False otherwise.
    """
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
            if lon_val < WESTERN_LONGITUDE_THRESHOLD:
                return True
        except (ValueError, TypeError):
            pass

    return False


def infer_regional_group(
    lat: Optional[float],
    lon: Optional[float],
    state_province: Optional[str] = None,
    habitat: Optional[str] = None,
    locality: Optional[str] = None
) -> str:
    """
    Assigns an ecological / physiographic regional group based on coordinates, state, and habitat.

    Regional Divisions:
      - Coastal_Plain_Sandhills: Atlantic/Gulf Coastal Plain, Sandhills, Maritime strand
      - Piedmont_Granite_Flatrocks: Piedmont province, granitic outcrops, flatrocks, diabase
      - Appalachian_Highlands: Blue Ridge, Ridge and Valley, Appalachian Plateau
      - Interior_Prairie_Midwest: Interior Low Plateaus, Ozarks, Tallgrass & Mixedgrass Prairies
      - Other_US: Geographic fallback for peripheral occurrences

    Args:
        lat: Decimal latitude (float or None if unrecorded).
        lon: Decimal longitude (float or None if unrecorded).
        state_province: State or province name.
        habitat: Habitat text description.
        locality: Locality text description.

    Returns:
        str: Standardized regional group identifier.
    """
    text_context = f"{state_province or ''} {habitat or ''} {locality or ''}".lower()

    if any(k in text_context for k in ["sandhill", "sand hill", "longleaf", "coastal plain", "dune", "maritime", "pocosin"]):
        return "Coastal_Plain_Sandhills"
    if any(k in text_context for k in ["flatrock", "granite outcrop", "granite", "diabase", "piedmont", "monadnock"]):
        return "Piedmont_Granite_Flatrocks"
    if any(k in text_context for k in ["blue ridge", "appalachian", "balds", "cove", "ridge and valley", "smoky", "high elevation"]):
        return "Appalachian_Highlands"
    if any(k in text_context for k in ["prairie", "glade", "limestone glade", "cedar glade", "ozark", "interior low plateau", "barren"]):
        return "Interior_Prairie_Midwest"

    state = (state_province or "").upper().strip()

    coastal_states = {"FL", "FLORIDA", "LA", "LOUISIANA", "MS", "MISSISSIPPI"}
    if state in coastal_states:
        return "Coastal_Plain_Sandhills"

    if state in {"NC", "NORTH CAROLINA", "SC", "SOUTH CAROLINA", "GA", "GEORGIA", "VA", "VIRGINIA"}:
        if lon is not None:
            if lon > -78.0:
                return "Coastal_Plain_Sandhills"
            elif -81.0 <= lon <= -78.0:
                return "Piedmont_Granite_Flatrocks"
            else:
                return "Appalachian_Highlands"
        else:
            return "Piedmont_Granite_Flatrocks"

    if state in {"TN", "TENNESSEE", "KY", "KENTUCKY", "WV", "WEST VIRGINIA", "PA", "PENNSYLVANIA"}:
        if lon is not None:
            if lon > -84.0:
                return "Appalachian_Highlands"
            else:
                return "Interior_Prairie_Midwest"
        else:
            return "Appalachian_Highlands"

    if state in {"MO", "MISSOURI", "AR", "ARKANSAS", "IL", "ILLINOIS", "IN", "INDIANA", "OH", "OHIO", "IA", "IOWA", "KS", "KANSAS", "NE", "NEBRASKA", "OK", "OKLAHOMA", "TX", "TEXAS"}:
        return "Interior_Prairie_Midwest"

    if lat is not None and lon is not None:
        if 24.0 <= lat <= 38.0 and -85.0 <= lon <= -75.0:
            return "Piedmont_Granite_Flatrocks"
        elif 34.0 <= lat <= 45.0 and -84.0 <= lon <= -70.0:
            return "Appalachian_Highlands"
        elif 28.0 <= lat <= 49.0 and -102.0 <= lon <= -84.0:
            return "Interior_Prairie_Midwest"
        elif 25.0 <= lat <= 35.0 and -98.0 <= lon <= -80.0:
            return "Coastal_Plain_Sandhills"

    return "Other_US"
