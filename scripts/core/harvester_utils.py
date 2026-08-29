import os
import sys
import re
import datetime
import argparse
import asyncio
import aiohttp
import logging
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
import numpy as np
import pandas as pd
import pygbif.occurrences as occ
from tqdm.asyncio import tqdm as async_tqdm

# Botanical & Morphometric Configuration Constants
from scripts.core.config import (
    DEFAULT_WORKSPACE,
    DEFAULT_RAW_DIR,
    DEFAULT_TARGET_TAXA,
    DEFAULT_OUTPUT_CSV,
    DEFAULT_SUMMARY_LOG,
    VALID_TYPE_STATUSES,
    SPECIALIST_PATTERNS,
    MAJOR_HERBARIA_CODES,
)


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
        # Guard against unrecorded/None coordinates before numeric comparison
        if lon is not None:
            if lon > -78.0:
                return "Coastal_Plain_Sandhills"
            elif -81.0 <= lon <= -78.0:
                return "Piedmont_Granite_Flatrocks"
            else:
                return "Appalachian_Highlands"
        else:
            # Fallback default assignment for Southeast/Mid-Atlantic states when coordinates are missing
            return "Piedmont_Granite_Flatrocks"
            
    if state in {"TN", "TENNESSEE", "KY", "KENTUCKY", "WV", "WEST VIRGINIA", "PA", "PENNSYLVANIA"}:
        # Guard against unrecorded/None coordinates before numeric comparison
        if lon is not None:
            if lon > -84.0:
                return "Appalachian_Highlands"
            else:
                return "Interior_Prairie_Midwest"
        else:
            # Fallback default assignment for Appalachian states when coordinates are missing
            return "Appalachian_Highlands"

    if state in {"MO", "MISSOURI", "AR", "ARKANSAS", "IL", "ILLINOIS", "IN", "INDIANA", "OH", "OHIO", "IA", "IOWA", "KS", "KANSAS", "NE", "NEBRASKA", "OK", "OKLAHOMA", "TX", "TEXAS"}:
        return "Interior_Prairie_Midwest"

    # Geographic bounding box fallback (safely execute only when both coordinates are non-None)
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


# =============================================================================
# Asynchronous Image Downloader Utilities
# =============================================================================

async def download_single_image(
    session: aiohttp.ClientSession,
    image_url: str,
    destination_path: Path,
    semaphore: asyncio.Semaphore,
    max_retries: int = 3
) -> bool:
    """
    Asynchronously downloads a single voucher image file to local storage with retry logic.
    
    Implements exponential backoff on transient network faults, checks file size integrity,
    and performs atomic writes via temporary files to avoid corrupting image assets.
    
    Args:
        session: Active aiohttp ClientSession.
        image_url: Remote image URL string.
        destination_path: Local filesystem destination Path.
        semaphore: Asyncio Semaphore for concurrency throttling.
        max_retries: Number of exponential backoff retry attempts (default: 3).
        
    Returns:
        bool: True if image downloaded successfully or already existed on disk, False otherwise.
    """
    # Skip download if file already exists locally with non-zero byte size (> 1KB)
    if destination_path.exists() and destination_path.stat().st_size > 1024:
        return True

    async with semaphore:
        for attempt in range(1, max_retries + 1):
            try:
                # Set reasonable connection and read timeouts to prevent hanging sockets
                timeout = aiohttp.ClientTimeout(total=45, connect=15)
                async with session.get(image_url, timeout=timeout) as response:
                    if response.status == 200:
                        content_type = response.headers.get("Content-Type", "").lower()
                        content = await response.read()
                        # Validate binary payload is non-trivial and matches typical image characteristics
                        if len(content) > 1024 and (not content_type or "image" in content_type or "octet-stream" in content_type or content[:3] == b'\xff\xd8\xff'):
                            # Ensure target directory structure exists
                            destination_path.parent.mkdir(parents=True, exist_ok=True)
                            # Perform atomic write via temporary file to prevent partial writes
                            temp_path = destination_path.with_suffix(".tmp")
                            with open(temp_path, "wb") as f:
                                f.write(content)
                            temp_path.replace(destination_path)
                            return True
                    elif response.status in {404, 410}:
                        # Permanent HTTP 404/410 errors indicate unrecoverable missing assets; skip retrying
                        return False
            except (aiohttp.ClientError, asyncio.TimeoutError, Exception):
                if attempt == max_retries:
                    return False
                # Exponential backoff delay: 1s, 2s, 4s...
                await asyncio.sleep(1.0 * (2 ** (attempt - 1)))
        return False


async def download_all_voucher_images(
    records_to_download: List[Tuple[str, Path]],
    concurrency_limit: int = 15,
    logger: Optional[logging.Logger] = None
) -> Dict[str, int]:
    """
    Coordinates asynchronous batch downloading of voucher images with concurrency control and progress tracking.
    
    Args:
        records_to_download: List of tuples where each tuple contains (image_url, destination_path).
        concurrency_limit: Maximum concurrent HTTP connections (default: 15).
        logger: Optional Logger instance for diagnostics.
        
    Returns:
        Dict[str, int]: Download metrics dictionary with 'success', 'skipped', and 'failed' counts.
    """
    semaphore = asyncio.Semaphore(concurrency_limit)
    headers = {
        "User-Agent": "PackeraResearchBot/1.0 (UNC Chapel Hill Herbarium; Evolutionary Morphometrics Lab)"
    }
    
    stats = {"success": 0, "skipped": 0, "failed": 0}
    
    # Filter items: distinguish cached vouchers vs pending downloads
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

    # Establish TCP connector with host and total connection limits
    connector = aiohttp.TCPConnector(limit=concurrency_limit, limit_per_host=5, ssl=False)
    async with aiohttp.ClientSession(headers=headers, connector=connector) as session:
        tasks = [
            download_single_image(session, url, dest, semaphore)
            for url, dest in pending
        ]
        
        # Track asynchronous task completion with async_tqdm progress bar
        results = await async_tqdm.gather(*tasks, desc="Downloading Voucher Sheets", unit="img")
        
        for success in results:
            if success:
                stats["success"] += 1
            else:
                stats["failed"] += 1

    return stats


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


