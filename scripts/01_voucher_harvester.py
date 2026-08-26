#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: Senior Bioinformatician & Biodiversity Data Engineer
Date: August 2026

Description:
    Automated Darwin Core voucher harvester and data-curation pipeline querying
    the GBIF API (via pygbif) for Packera dubia and related taxa in the US.
    Implements multi-tiered taxonomic authority scoring to mitigate 20-40%
    herbarium misidentifications, calculates harmonic circular phenology,
    assigns regional ecological groups, asynchronously downloads high-res
    specimen sheets, and exports standardized curated metadata.

Target Taxa:
    - Packera dubia (Spreng.) Trock & Mabb.
    - Packera tomentosa (Michx.) C. Jeffrey / Senecio tomentosus Michx.
    - Packera anonyma (Alph. Wood) W.A. Weber & Á. Löve
    - Packera plattensis (Nutt.) W.A. Weber & Á. Löve
    - Packera paupercula (Michx.) Á. Löve & D. Löve

Determiner Credibility Tiers:
    - Tier_1_Gold: Nomenclatural types (Holotypes, Isotypes, Lectotypes, etc.)
                   OR determined/annotated by recognized monographic specialists:
                   Barkley, Trock, Kowal, Weakley, Bain, Mahoney, Fuller.
    - Tier_2_Silver: Determined by botanists affiliated with major regional
                     research herbaria (NCU, GA, US, NY, BRIT, MO, WIS, VDB, FLAS)
                     with rich locality and ecological habitat metadata.
    - Tier_3_Bronze: General floristic collections, unverified collector IDs,
                     or blank determiner fields.

Usage:
    py scripts/01_voucher_harvester.py --download-images --max-records 1000
===============================================================================
"""

import os
import sys
import re
import math
import time
import datetime
import argparse
import asyncio
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

import numpy as np
import pandas as pd
import aiohttp
from tqdm.asyncio import tqdm as async_tqdm
from tqdm import tqdm
import pygbif.occurrences as occ
import pygbif.species as species

# -----------------------------------------------------------------------------
# Configuration and Constants
# -----------------------------------------------------------------------------

# Target focal taxa in the Packera dubia complex and allied North American clades
DEFAULT_TARGET_TAXA = [
    "Packera dubia",
    "Packera tomentosa",
    "Senecio tomentosus",
    "Packera anonyma",
    "Packera plattensis",
    "Packera paupercula",
]

# Recognized monographic taxonomic specialists for Packera / Senecioneae (Tier 1 Gold)
SPECIALIST_PATTERNS = [
    r"\bBarkley\b",
    r"\bT\.?\s*M\.?\s*Barkley\b",
    r"\bTheodore\s+M\.?\s+Barkley\b",
    r"\bTrock\b",
    r"\bD\.?\s*K\.?\s*Trock\b",
    r"\bDebra\s+K\.?\s+Trock\b",
    r"\bKowal\b",
    r"\bR\.?\s*R\.?\s*Kowal\b",
    r"\bRobert\s+R\.?\s+Kowal\b",
    r"\bWeakley\b",
    r"\bA\.?\s*S\.?\s*Weakley\b",
    r"\bAlan\s+S\.?\s+Weakley\b",
    r"\bBain\b",
    r"\bJ\.?\s*F\.?\s*Bain\b",
    r"\bJohn\s+F\.?\s+Bain\b",
    r"\bMahoney\b",
    r"\bA\.?\s*M\.?\s*Mahoney\b",
    r"\bAlison\s+M\.?\s+Mahoney\b",
    r"\bFuller\b",
    r"\bJ\.?\s*B\.?\s*Fuller\b",
    r"\bBrandon\s+Fuller\b",
]

# Major regional research herbaria with high curatorial standards (Tier 2 Silver)
MAJOR_HERBARIA_CODES = {
    "NCU", "GA", "US", "NY", "BRIT", "MO", "WIS", "VDB", "FLAS", "TEX", "LL", "TENN", "F"
}

# Recognized nomenclatural type status designations (Tier 1 Gold)
VALID_TYPE_STATUSES = {
    "HOLOTYPE", "ISOTYPE", "LECTOTYPE", "ISOLECTOTYPE", "SYNTYPE", "ISOSYNTYPE",
    "NEOTYPE", "ISONEOTYPE", "PARATYPE", "ISOPARATYPE", "EPITYPE", "TYPE", "COTYPE"
}

# Base directories
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_VOUCHERS_DIR = DATA_DIR / "raw_vouchers"
TABLES_DIR = DATA_DIR / "tables"
OUTPUTS_DIR = BASE_DIR / "outputs"
REPORTS_DIR = OUTPUTS_DIR / "reports"

# Default Output paths
DEFAULT_OUTPUT_CSV = TABLES_DIR / "curated_vouchers.csv"
DEFAULT_SUMMARY_LOG = REPORTS_DIR / "voucher_ingestion_summary.log"

# -----------------------------------------------------------------------------
# Logging Configuration
# -----------------------------------------------------------------------------

def setup_logger(log_file_path: Optional[Path] = None) -> logging.Logger:
    """
    Configures a thread-safe and formatted logger outputting to both console and disk.
    
    Args:
        log_file_path: Optional destination Path for writing log messages.
        
    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger("VoucherHarvester")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # Formatter for consistent timestamps and log levels
    formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console Handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File Handler (if requested and directory created)
    if log_file_path:
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file_path, mode="w", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger

# -----------------------------------------------------------------------------
# Darwin Core Metadata Parsing & Quality Scoring Functions
# -----------------------------------------------------------------------------

def sanitize_filename(name: str) -> str:
    """
    Sanitizes arbitrary strings into safe, valid filesystem filenames across OS platforms.
    
    Args:
        name: Raw identifier or catalog number string.
        
    Returns:
        str: Cleaned alphanumeric filename string safe for Windows/Linux filesystems.
    """
    # Replace slashes, colons, spaces, and invalid filesystem characters with underscores
    clean = re.sub(r'[\\/*?:"<>|\s]+', '_', str(name).strip())
    # Remove leading/trailing underscores and dots
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
    # Clean strings
    type_status = str(type_status_raw).strip() if type_status_raw is not None else ""
    identified_by = str(identified_by_raw).strip() if identified_by_raw is not None else ""
    recorded_by = str(recorded_by_raw).strip() if recorded_by_raw is not None else ""
    history = str(history_raw).strip() if history_raw is not None else ""
    institution = str(institution_code_raw).strip().upper() if institution_code_raw is not None else ""
    locality = str(locality_raw).strip() if locality_raw is not None else ""
    habitat = str(habitat_raw).strip() if habitat_raw is not None else ""

    # Combined determiner text for regex audit
    determiner_raw = identified_by if identified_by else (history if history else recorded_by)
    combined_audit_text = f"{type_status} | {identified_by} | {history} | {recorded_by}"

    # Normalize type status representation
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
        # Validate calendar date
        collection_date = datetime.date(y, m, d)
        doy = collection_date.timetuple().tm_yday  # 1 to 366
        
        # Calculate circular coordinates based on astronomical solar year (365.25 days)
        theta = 2.0 * math.pi * (doy / 365.25)
        pheno_sin = round(math.sin(theta), 6)
        pheno_cos = round(math.cos(theta), 6)
        
        return doy, pheno_sin, pheno_cos
    except (ValueError, TypeError, OverflowError):
        return None


def infer_regional_group(
    lat: float,
    lon: float,
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
        lat: Decimal latitude.
        lon: Decimal longitude.
        state_province: State or province name.
        habitat: Habitat text description.
        locality: Locality text description.
        
    Returns:
        str: Standardized regional group identifier.
    """
    text_context = f"{state_province or ''} {habitat or ''} {locality or ''}".lower()
    
    # Keyword checks from ecological descriptions
    if any(k in text_context for k in ["sandhill", "sand hill", "longleaf", "coastal plain", "dune", "maritime", "pocosin"]):
        return "Coastal_Plain_Sandhills"
    if any(k in text_context for k in ["flatrock", "granite outcrop", "granite", "diabase", "piedmont", "monadnock"]):
        return "Piedmont_Granite_Flatrocks"
    if any(k in text_context for k in ["blue ridge", "appalachian", "balds", "cove", "ridge and valley", "smoky", "high elevation"]):
        return "Appalachian_Highlands"
    if any(k in text_context for k in ["prairie", "glade", "limestone glade", "cedar glade", "ozark", "interior low plateau", "barren"]):
        return "Interior_Prairie_Midwest"

    # Coordinate and state bounding heuristic for North America
    state = (state_province or "").upper().strip()
    
    # Southeastern Atlantic & Gulf Coastal Plain
    coastal_states = {"FL", "FLORIDA", "LA", "LOUISIANA", "MS", "MISSISSIPPI"}
    if state in coastal_states:
        return "Coastal_Plain_Sandhills"
        
    # State-based and coordinate classification
    if state in {"NC", "NORTH CAROLINA", "SC", "SOUTH CAROLINA", "GA", "GEORGIA", "VA", "VIRGINIA"}:
        if lon > -78.0:
            return "Coastal_Plain_Sandhills"
        elif -81.0 <= lon <= -78.0:
            return "Piedmont_Granite_Flatrocks"
        else:
            return "Appalachian_Highlands"
            
    if state in {"TN", "TENNESSEE", "KY", "KENTUCKY", "WV", "WEST VIRGINIA", "PA", "PENNSYLVANIA"}:
        if lon > -84.0:
            return "Appalachian_Highlands"
        else:
            return "Interior_Prairie_Midwest"

    if state in {"MO", "MISSOURI", "AR", "ARKANSAS", "IL", "ILLINOIS", "IN", "INDIANA", "OH", "OHIO", "IA", "IOWA", "KS", "KANSAS", "NE", "NEBRASKA", "OK", "OKLAHOMA", "TX", "TEXAS"}:
        return "Interior_Prairie_Midwest"

    # Geographic bounding box fallback
    if 24.0 <= lat <= 38.0 and -85.0 <= lon <= -75.0:
        return "Piedmont_Granite_Flatrocks"
    elif 34.0 <= lat <= 45.0 and -84.0 <= lon <= -70.0:
        return "Appalachian_Highlands"
    elif 28.0 <= lat <= 49.0 and -102.0 <= lon <= -84.0:
        return "Interior_Prairie_Midwest"
    elif 25.0 <= lat <= 35.0 and -98.0 <= lon <= -80.0:
        return "Coastal_Plain_Sandhills"

    return "Other_US"


def extract_high_res_image_url(media_list: Optional[List[Dict[str, Any]]]) -> Optional[str]:
    """
    Parses Darwin Core media records to extract the highest-quality specimen sheet image URL.
    
    Args:
        media_list: List of media dictionaries from the GBIF occurrence payload.
        
    Returns:
        Optional[str]: Verified HTTP/HTTPS image URL, or None if unavailable.
    """
    if not media_list or not isinstance(media_list, list):
        return None

    candidate_urls = []
    for item in media_list:
        if not isinstance(item, dict):
            continue
        
        # Check media type and format
        m_type = item.get("type", "")
        m_format = item.get("format", "").lower()
        identifier = item.get("identifier", "").strip()

        if not identifier or not identifier.startswith(("http://", "https://")):
            continue

        # Prioritize JPEG/PNG or StillImage
        if m_type == "StillImage" or "image" in m_format or identifier.lower().endswith((".jpg", ".jpeg", ".png", ".tif", ".tiff")):
            candidate_urls.append(identifier)
        else:
            candidate_urls.append(identifier)

    return candidate_urls[0] if candidate_urls else None

# -----------------------------------------------------------------------------
# Asynchronous Image Downloader
# -----------------------------------------------------------------------------

async def download_single_image(
    session: aiohttp.ClientSession,
    image_url: str,
    destination_path: Path,
    semaphore: asyncio.Semaphore,
    max_retries: int = 3
) -> bool:
    """
    Asynchronously downloads a single voucher image file to local storage with retry logic.
    
    Args:
        session: Active aiohttp ClientSession.
        image_url: Remote image URL.
        destination_path: Local filesystem destination Path.
        semaphore: Asyncio Semaphore for concurrency throttling.
        max_retries: Number of exponential backoff retry attempts.
        
    Returns:
        bool: True if image downloaded successfully or already existed, False otherwise.
    """
    # Skip if file already exists with non-zero size
    if destination_path.exists() and destination_path.stat().st_size > 1024:
        return True

    async with semaphore:
        for attempt in range(1, max_retries + 1):
            try:
                # Set reasonable connection and read timeouts
                timeout = aiohttp.ClientTimeout(total=45, connect=15)
                async with session.get(image_url, timeout=timeout) as response:
                    if response.status == 200:
                        content_type = response.headers.get("Content-Type", "").lower()
                        # Verify we received image binary data
                        content = await response.read()
                        if len(content) > 1024 and (not content_type or "image" in content_type or "octet-stream" in content_type or content[:3] == b'\xff\xd8\xff'):
                            # Ensure parent directory exists
                            destination_path.parent.mkdir(parents=True, exist_ok=True)
                            # Write atomically using temporary file
                            temp_path = destination_path.with_suffix(".tmp")
                            with open(temp_path, "wb") as f:
                                f.write(content)
                            temp_path.replace(destination_path)
                            return True
                    elif response.status in {404, 410}:
                        # Permanent client error; do not retry
                        return False
            except (aiohttp.ClientError, asyncio.TimeoutError, Exception):
                if attempt == max_retries:
                    return False
                await asyncio.sleep(1.0 * (2 ** (attempt - 1)))
        return False


async def download_all_voucher_images(
    records_to_download: List[Tuple[str, Path]],
    concurrency_limit: int = 15,
    logger: Optional[logging.Logger] = None
) -> Dict[str, int]:
    """
    Coordinates asynchronous batch downloading of voucher images with a visual progress bar.
    
    Args:
        records_to_download: List of tuples (image_url, destination_path).
        concurrency_limit: Maximum concurrent HTTP requests.
        logger: Logger instance.
        
    Returns:
        Dict[str, int]: Download statistics (success, skipped, failed).
    """
    semaphore = asyncio.Semaphore(concurrency_limit)
    headers = {
        "User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"
    }
    
    stats = {"success": 0, "skipped": 0, "failed": 0}
    
    # Filter items that need actual downloading vs already cached
    pending = []
    for url, dest in records_to_download:
        if dest.exists() and dest.stat().st_size > 1024:
            stats["skipped"] += 1
        else:
            pending.append((url, dest))

    if not pending:
        if logger:
            logger.info(f"All {stats['skipped']} voucher images are already cached locally.")
        return stats

    if logger:
        logger.info(f"Initiating asynchronous download of {len(pending)} pending images (Concurrency: {concurrency_limit})...")

    connector = aiohttp.TCPConnector(limit=concurrency_limit, limit_per_host=5, ssl=False)
    async with aiohttp.ClientSession(headers=headers, connector=connector) as session:
        tasks = [
            download_single_image(session, url, dest, semaphore)
            for url, dest in pending
        ]
        
        # Track progress with async tqdm
        results = await async_tqdm.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
        
        for success in results:
            if success:
                stats["success"] += 1
            else:
                stats["failed"] += 1

    return stats

# -----------------------------------------------------------------------------
# GBIF Query & Ingestion Pipeline
# -----------------------------------------------------------------------------

def harvest_taxa_occurrences(
    taxa_list: List[str],
    max_uncertainty_meters: float = 5000.0,
    max_records_per_taxon: int = 1000,
    logger: Optional[logging.Logger] = None
) -> pd.DataFrame:
    """
    Executes paginated queries to the GBIF Occurrence API and filters records according to DwC standards.
    
    Filters applied:
      - Basis of record: PRESERVED_SPECIMEN
      - Geographic scope: United States (country='US')
      - Spatial validation: Non-null coordinates, coordinateUncertaintyInMeters <= max_uncertainty_meters
      - Temporal validation: Valid collection dates (year, month, day)
      - Media validation: High-resolution sheet image URL present
      
    Args:
        taxa_list: List of botanical species names to query.
        max_uncertainty_meters: Maximum allowed georeferencing uncertainty in meters.
        max_records_per_taxon: Maximum records to harvest per taxon query.
        logger: Logger instance.
        
    Returns:
        pd.DataFrame: Curated dataframe containing parsed metadata and authority scores.
    """
    if logger is None:
        logger = logging.getLogger("VoucherHarvester")

    all_curated_records = []
    seen_catalog_keys = set()

    for taxon in taxa_list:
        logger.info(f"Querying GBIF API for taxon: '{taxon}' (country=US, basisOfRecord=PRESERVED_SPECIMEN)...")
        
        offset = 0
        limit = 300  # GBIF page size limit
        taxon_harvested = 0
        taxon_retained = 0

        while taxon_harvested < max_records_per_taxon:
            fetch_limit = min(limit, max_records_per_taxon - taxon_harvested)
            try:
                # Query GBIF occurrence search endpoint
                response = occ.search(
                    scientificName=taxon,
                    country="US",
                    basisOfRecord="PRESERVED_SPECIMEN",
                    hasCoordinate=True,
                    limit=fetch_limit,
                    offset=offset
                )
            except Exception as e:
                logger.error(f"GBIF API query error for '{taxon}' at offset {offset}: {e}")
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
                if lat is None or lon is None:
                    continue

                try:
                    lat_val = float(lat)
                    lon_val = float(lon)
                except (ValueError, TypeError):
                    continue

                # Coordinate uncertainty check (must be <= 5000 meters)
                uncertainty_raw = rec.get("coordinateUncertaintyInMeters")
                if uncertainty_raw is None:
                    # Filter out null uncertainty as per strict pipeline requirements
                    continue
                
                try:
                    uncertainty_val = float(uncertainty_raw)
                    if uncertainty_val > max_uncertainty_meters:
                        continue
                except (ValueError, TypeError):
                    continue

                # 2. Temporal & Phenology Validation (year, month, day)
                year = rec.get("year")
                month = rec.get("month")
                day = rec.get("day")
                pheno_res = calculate_circular_phenology(year, month, day)
                if pheno_res is None:
                    continue
                doy, pheno_sin, pheno_cos = pheno_res

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

                # Ensure unique catalog identifier across duplicate uploads
                unique_key = (catalog_number, gbif_key)
                if unique_key in seen_catalog_keys:
                    continue
                seen_catalog_keys.add(unique_key)

                # Local image path destination
                relative_image_path = f"data/raw_vouchers/{catalog_number}.jpg"

                # 5. Taxonomic Authority & Determiner Tier Stratification
                type_status_raw = rec.get("typeStatus")
                identified_by_raw = rec.get("identifiedBy")
                recorded_by_raw = rec.get("recordedBy")
                history_raw = rec.get("verbatimIdentificationHistory")
                locality_raw = rec.get("locality") or rec.get("verbatimLocality")
                habitat_raw = rec.get("habitat")
                species_raw = rec.get("scientificName") or rec.get("species") or taxon

                determiner_tier, type_status, determiner_raw = parse_determiner_tier(
                    type_status_raw=type_status_raw,
                    identified_by_raw=identified_by_raw,
                    recorded_by_raw=recorded_by_raw,
                    history_raw=history_raw,
                    institution_code_raw=inst_code,
                    locality_raw=locality_raw,
                    habitat_raw=habitat_raw
                )

                # 6. Regional Ecological Group Assignment
                state_prov = rec.get("stateProvince")
                regional_group = infer_regional_group(
                    lat=lat_val,
                    lon=lon_val,
                    state_province=state_prov,
                    habitat=habitat_raw,
                    locality=locality_raw
                )

                # Append standardized record
                curated_record = {
                    "catalogNumber": catalog_number,
                    "institutionCode": inst_code,
                    "species_raw": species_raw,
                    "determiner_raw": determiner_raw,
                    "determiner_tier": determiner_tier,
                    "type_status": type_status,
                    "latitude": lat_val,
                    "longitude": lon_val,
                    "coordinateUncertainty": uncertainty_val,
                    "doy": doy,
                    "pheno_sin": pheno_sin,
                    "pheno_cos": pheno_cos,
                    "regional_group": regional_group,
                    "image_path": relative_image_path,
                    "_image_url": image_url,  # Temporary internal column for downloading
                }
                all_curated_records.append(curated_record)
                taxon_retained += 1

            offset += len(results)
            if offset >= count:
                break
            time.sleep(0.1)  # Respectful GBIF rate limiting

        logger.info(f"Taxon '{taxon}': Processed {taxon_harvested} occurrences -> Retained {taxon_retained} curated records meeting quality filters.")

    # Convert to DataFrame
    df = pd.DataFrame(all_curated_records)
    return df

# -----------------------------------------------------------------------------
# Summary Reporting & Execution Flow
# -----------------------------------------------------------------------------

def print_and_log_summary(df: pd.DataFrame, download_stats: Optional[Dict[str, int]], logger: logging.Logger) -> None:
    """
    Generates a publication-grade terminal and file log summary of the harvested dataset.
    
    Args:
        df: Curated vouchers DataFrame.
        download_stats: Optional download counters dictionary.
        logger: Logger instance.
    """
    total_records = len(df)
    logger.info("=" * 80)
    logger.info("                  PACKERA VOUCHER INGESTION & CURATION SUMMARY                  ")
    logger.info("=" * 80)
    logger.info(f"Total Quality-Filtered Specimen Vouchers: {total_records:,}")

    if total_records == 0:
        logger.warning("No records were retained. Check query parameters or network connection.")
        logger.info("=" * 80)
        return

    # 1. Determiner Tier Breakdown
    tier_counts = df["determiner_tier"].value_counts()
    logger.info("\n--- TAXONOMIC DETERMINER AUTHORITY STRATIFICATION ---")
    for tier in ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"]:
        cnt = tier_counts.get(tier, 0)
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {tier:<15} : {cnt:>5} records ({pct:>5.1f}%)")

    # 2. Species Breakdown
    species_counts = df["species_raw"].value_counts()
    logger.info("\n--- TAXON DISTRIBUTION (RAW DETERMINATIONS) ---")
    for sp, cnt in species_counts.head(8).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {sp:<45} : {cnt:>5} records ({pct:>5.1f}%)")

    # 3. Regional Ecological Groups Breakdown
    region_counts = df["regional_group"].value_counts()
    logger.info("\n--- REGIONAL ECO-GEOGRAPHIC GROUPS ---")
    for reg, cnt in region_counts.items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {reg:<30} : {cnt:>5} records ({pct:>5.1f}%)")

    # 4. Herbarium Institutions (Top 10)
    inst_counts = df["institutionCode"].value_counts()
    logger.info("\n--- TOP HERBARIUM INSTITUTIONS ---")
    for inst, cnt in inst_counts.head(10).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {inst:<15} : {cnt:>5} records ({pct:>5.1f}%)")

    # 5. Image Download Summary
    if download_stats:
        logger.info("\n--- HIGH-RESOLUTION IMAGE DOWNLOAD STATUS ---")
        logger.info(f"  * Downloaded Successfully : {download_stats.get('success', 0):>5}")
        logger.info(f"  * Cached / Skipped        : {download_stats.get('skipped', 0):>5}")
        logger.info(f"  * Failed / Inaccessible   : {download_stats.get('failed', 0):>5}")

    logger.info("=" * 80)


def main() -> None:
    """
    Main entry point for command-line execution of the voucher harvester.
    """
    parser = argparse.ArgumentParser(
        description="Automated GBIF Voucher Harvester & Determiner Authority Scorer for Packera dubia."
    )
    parser.add_argument(
        "--taxa",
        nargs="+",
        default=DEFAULT_TARGET_TAXA,
        help="List of scientific binomials or names to harvest from GBIF."
    )
    parser.add_argument(
        "--max-uncertainty",
        type=float,
        default=5000.0,
        help="Maximum coordinate uncertainty in meters (default: 5000.0)."
    )
    parser.add_argument(
        "--max-records-per-taxon",
        type=int,
        default=1000,
        help="Maximum records to harvest per taxon (default: 1000)."
    )
    parser.add_argument(
        "--download-images",
        action="store_true",
        default=False,
        help="Flag to enable asynchronous high-resolution specimen image downloading."
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=15,
        help="Max concurrent asynchronous image downloads (default: 15)."
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=str(DEFAULT_OUTPUT_CSV),
        help=f"Target path for curated metadata CSV (default: {DEFAULT_OUTPUT_CSV})."
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default=str(DEFAULT_SUMMARY_LOG),
        help=f"Path for summary log file (default: {DEFAULT_SUMMARY_LOG})."
    )

    args = parser.parse_args()

    # Create destination directories
    output_csv_path = Path(args.output_csv)
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    RAW_VOUCHERS_DIR.mkdir(parents=True, exist_ok=True)
    log_file_path = Path(args.log_file)
    log_file_path.parent.mkdir(parents=True, exist_ok=True)

    # Initialize logger
    logger = setup_logger(log_file_path)
    logger.info("Starting Packera Voucher Ingestion & Authority Stratification Pipeline...")
    logger.info(f"Target Taxa: {args.taxa}")
    logger.info(f"Max Coordinate Uncertainty Threshold: {args.max_uncertainty} m")

    # Step 1: Harvest and curate metadata records from GBIF
    df_curated = harvest_taxa_occurrences(
        taxa_list=args.taxa,
        max_uncertainty_meters=args.max_uncertainty,
        max_records_per_taxon=args.max_records_per_taxon,
        logger=logger
    )

    download_stats = None

    # Step 2: Asynchronously download specimen images if enabled
    if args.download_images and not df_curated.empty and "_image_url" in df_curated.columns:
        download_queue = []
        for _, row in df_curated.iterrows():
            url = row["_image_url"]
            dest = BASE_DIR / row["image_path"]
            download_queue.append((url, dest))

        logger.info(f"Starting asynchronous download of {len(download_queue)} voucher sheets...")
        download_stats = asyncio.run(
            download_all_voucher_images(
                records_to_download=download_queue,
                concurrency_limit=args.concurrency,
                logger=logger
            )
        )

    # Step 3: Export standardized CSV (excluding internal temporary columns)
    export_columns = [
        "catalogNumber",
        "institutionCode",
        "species_raw",
        "determiner_raw",
        "determiner_tier",
        "type_status",
        "latitude",
        "longitude",
        "coordinateUncertainty",
        "doy",
        "pheno_sin",
        "pheno_cos",
        "regional_group",
        "image_path"
    ]
    
    if not df_curated.empty:
        df_export = df_curated[[col for col in export_columns if col in df_curated.columns]]
        df_export.to_csv(output_csv_path, index=False, encoding="utf-8")
        logger.info(f"Successfully exported {len(df_export)} curated records to: {output_csv_path}")
    else:
        # Create empty table with standardized headers
        pd.DataFrame(columns=export_columns).to_csv(output_csv_path, index=False, encoding="utf-8")
        logger.warning(f"No records met all filtering criteria. Created empty table at: {output_csv_path}")

    # Step 4: Output comprehensive summary log
    print_and_log_summary(df_curated, download_stats, logger)
    logger.info("Pipeline execution completed successfully.")


if __name__ == "__main__":
    main()
