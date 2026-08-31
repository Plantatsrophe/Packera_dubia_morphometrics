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
from PIL import Image
Image.MAX_IMAGE_PIXELS = None  # Allow decompression of high-resolution botanical sheets (>89 MP)

# Botanical & Morphometric Configuration Constants
from scripts.core.config import (
    DEFAULT_WORKSPACE,
    DEFAULT_RAW_DIR,
    DEFAULT_TARGET_TAXA,
    DEFAULT_OUTPUT_CSV,
    DEFAULT_SUMMARY_LOG,
    DEFAULT_QUARANTINE_DIR,
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_FILE_SIZE_KB,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN,
    VALID_TYPE_STATUSES,
    SPECIALIST_PATTERNS,
    MAJOR_HERBARIA_CODES,
    EXCLUDED_WESTERN_STATES,
    WESTERN_LONGITUDE_THRESHOLD,
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


def is_excluded_western_region(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None
) -> bool:
    """
    Determines whether a record originates from a US state or territory farther west than Texas and Oklahoma.
    
    Excluded regions include the Rocky Mountains, Intermountain West, Southwest, Pacific Northwest,
    West Coast, and Pacific states (e.g. CO, NM, WY, MT, UT, AZ, NV, ID, WA, OR, CA, AK, HI).
    
    Disambiguates 'Washington' (State, excluded) from 'Washington, D.C.' / 'District of Columbia' (East, retained).
    
    Args:
        state_province: State or province name / abbreviation string (or None/NaN).
        lat: Optional decimal latitude float.
        lon: Optional decimal longitude float.
        
    Returns:
        bool: True if the record should be excluded as a western locality; False otherwise.
    """
    if state_province is not None and not (isinstance(state_province, float) and math.isnan(state_province)):
        raw_state = str(state_province).strip()
        # Clean state string (remove trailing punctuation, qualifiers like '(State)', etc.)
        cleaned = re.sub(r"\(state\)", "", raw_state, flags=re.IGNORECASE).strip(" ._,-")
        upper_state = cleaned.upper()
        
        # Disambiguation: Preserve District of Columbia / Washington, D.C.
        if upper_state in {"WASHINGTON, D.C.", "WASHINGTON D.C.", "WASHINGTON DC", "DISTRICT OF COLUMBIA", "DC"}:
            return False
            
        if upper_state in EXCLUDED_WESTERN_STATES:
            return True
            
    # Coordinate fallback: If state is unrecorded or not in excluded list, check western longitude bound
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


def optimize_herbarium_image_url(url: str) -> str:
    """
    Transforms provider-specific URLs to request full-resolution original scans rather than
    dimension-clamped web previews or downscaled thumbnails.
    
    Args:
        url: Raw image URL string.
        
    Returns:
        str: Optimized image URL requesting maximum resolution.
    """
    if not url or not isinstance(url, str):
        return ""

    optimized = url.strip()

    # 1. Smithsonian NMNH URLs: strip dimension clamp (e.g. &h=2000) to request native resolution
    if "collections.nmnh.si.edu/media/" in optimized:
        optimized = re.sub(r"[?&][hw]=\d+", "", optimized)
        if "?" not in optimized and "&" in optimized:
            optimized = optimized.replace("&", "?", 1)

    # 2. Symbiota / SERNEC / SEINet / CCH portals: replace web/thumbnail paths with original / large
    if any(k in optimized.lower() for k in ["symbiota", "sernec", "seinet", "cch2", "swbiodiversity"]):
        optimized = re.sub(r"/(?:web|tn|thumbnails?)/", "/orig/", optimized, flags=re.IGNORECASE)
        optimized = re.sub(r"_(?:tn|web|sm)\.(jpe?g|png)", r"_lg.\1", optimized, flags=re.IGNORECASE)

    # 3. IIIF Image API Endpoints: replace constrained dimensions with /full/max/
    if "/full/!" in optimized or "/full/pct:" in optimized or re.search(r"/full/\d+,\d*/", optimized):
        optimized = re.sub(r"/full/(?:!?\d+,\d*|pct:\d+)/", "/full/max/", optimized)

    return optimized


def extract_high_res_image_url(media_list: Optional[List[Dict[str, Any]]]) -> Optional[str]:
    """
    Parses Darwin Core media records, scoring and prioritizing the highest-quality specimen sheet image URL.
    
    Scores each candidate:
      + High resolution indicators (_lg, original, /orig/, master, highres, full)
      - Low resolution / thumbnail indicators (_tn, thumb, preview, _sm, small, icon)
    Applies institution-specific URL optimizations to ensure full-resolution asset retrieval.
    
    Args:
        media_list: List of media dictionaries from the GBIF occurrence payload.
        
    Returns:
        Optional[str]: Verified and optimized HTTP/HTTPS image URL, or None if unavailable.
    """
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

        # Media type and format baseline
        if m_type == "StillImage" or "image" in m_format:
            score += 50.0
        if identifier.lower().endswith((".jpg", ".jpeg", ".png", ".tif", ".tiff")):
            score += 20.0

        ident_lower = identifier.lower()

        # High-resolution indicators (+ points)
        if any(h in ident_lower for h in ["_lg", "_large", "original", "/orig/", "/master/", "/highres/", "/full/", "hires", "high_res"]):
            score += 100.0
        if "max" in ident_lower:
            score += 30.0

        # Low-resolution / thumbnail penalties (- points)
        if any(t in ident_lower for t in ["_tn", "_thumb", "thumbnail", "_sm", "_small", "preview", "icon", "mini"]):
            score -= 150.0
        if "detailimages" in ident_lower:
            score -= 25.0
        if re.search(r"[?&]h=(?:[1-9]\d{0,2}|1\d{3}|2000)\b", ident_lower):
            score -= 10.0

        # Apply URL optimization rewrite
        optimized_url = optimize_herbarium_image_url(identifier)
        scored_candidates.append((score, optimized_url))

    if not scored_candidates:
        return None

    # Sort descending by quality score and select highest-scoring candidate URL
    scored_candidates.sort(key=lambda x: x[0], reverse=True)
    return scored_candidates[0][1]


def validate_image_quality(
    image_path: Path,
    min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
    min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB,
    check_sharpness: bool = False,
    min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN,
) -> Tuple[bool, Dict[str, Any]]:
    """
    Evaluates image resolution, byte size, and optical sharpness metrics.
    
    Args:
        image_path: Path to the local image file.
        min_megapixels: Minimum resolution threshold in Megapixels (W x H / 1e6).
        min_file_size_kb: Minimum compressed file size in Kilobytes.
        check_sharpness: Whether to compute Laplacian variance edge sharpness.
        min_sharpness: Minimum acceptable Laplacian variance score.
        
    Returns:
        Tuple[bool, Dict[str, Any]]: (is_valid, metrics_dict)
    """
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

    # 5. Image Download & Quality Summary
    if download_stats:
        logger.info("\n--- SPECIMEN IMAGE DOWNLOAD & QUALITY STATUS ---")
        logger.info(f"  * Downloaded Successfully : {download_stats.get('success', 0):>5}")
        logger.info(f"  * Cached / Skipped        : {download_stats.get('skipped', 0):>5}")
        logger.info(f"  * Quality Filter Rejected : {download_stats.get('quality_rejected', 0):>5}")
        logger.info(f"  * Failed / Inaccessible   : {download_stats.get('failed', 0):>5}")
        if "median_mp" in download_stats:
            logger.info(f"  * Median Image Resolution : {download_stats.get('median_mp', 0.0):>5.2f} Megapixels")

    logger.info("=" * 80)


