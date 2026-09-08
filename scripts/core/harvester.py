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
import datetime
import logging
import math
import os
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import aiohttp
import numpy as np
import pandas as pd
from PIL import Image
import pygbif.occurrences as occ
from tqdm.asyncio import tqdm as async_tqdm

# Allow decompression of high-resolution botanical herbarium sheets (>89 MP)
Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
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
    VALID_TYPE_STATUSES,
    WESTERN_LONGITUDE_THRESHOLD,
)
from scripts.core.logger import setup_logging

EXPORT_COLUMNS = [
    "catalogNumber",
    "institutionCode",
    "scientificName",
    "species_raw",
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
]


def setup_logger(log_file_path: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    """Configures or retrieves a formatted logger for voucher harvesting."""
    return setup_logging(log_file=log_file_path, verbose=verbose, name="VoucherHarvester")


def sanitize_filename(name: str) -> str:
    """Sanitizes arbitrary strings into safe, valid filesystem filenames across platforms."""
    clean = re.sub(r'[\\/*?:"<>|\s]+', "_", str(name).strip())
    clean = clean.strip("._")
    return clean if clean else "voucher_unknown"


def parse_determiner_tier(
    type_status_raw: Optional[str],
    identified_by_raw: Optional[str],
    recorded_by_raw: Optional[str],
    history_raw: Optional[str],
    institution_code_raw: Optional[str],
    locality_raw: Optional[str],
    habitat_raw: Optional[str],
) -> Tuple[str, str, str]:
    """
    Evaluates Darwin Core fields to assign a taxonomic determination credibility tier.

    Tiers:
      - Tier_1_Gold: Primary/secondary nomenclatural types or verified specialist annotations.
      - Tier_2_Silver: Determinations from major research herbaria with complete ecological locality.
      - Tier_3_Bronze: Unverified candidate collections, general collectors, or missing determiner.

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

    is_specialist = False
    for pattern in SPECIALIST_PATTERNS:
        if re.search(pattern, combined_audit_text, re.IGNORECASE):
            is_specialist = True
            break

    if is_type or is_specialist:
        return "Tier_1_Gold", type_status_clean, determiner_raw

    is_major_herbarium = institution in MAJOR_HERBARIA_CODES
    has_rich_locality = len(locality) > 10 or len(habitat) > 5
    has_determiner = bool(identified_by and identified_by.lower() not in {"unknown", "anonymous", "none", "null"})

    if is_major_herbarium and has_rich_locality and has_determiner:
        return "Tier_2_Silver", type_status_clean, determiner_raw

    return "Tier_3_Bronze", type_status_clean, determiner_raw


def is_excluded_western_region(
    state_province: Optional[str],
    lat: Optional[float] = None,
    lon: Optional[float] = None,
) -> bool:
    """Determines whether a record originates from a US state farther west than TX and OK."""
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
    locality: Optional[str] = None,
) -> str:
    """Assigns an ecological / physiographic regional group based on coordinates, state, and habitat."""
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
        return "Piedmont_Granite_Flatrocks"

    if state in {"TN", "TENNESSEE", "KY", "KENTUCKY", "WV", "WEST VIRGINIA", "PA", "PENNSYLVANIA"}:
        if lon is not None:
            return "Appalachian_Highlands" if lon > -84.0 else "Interior_Prairie_Midwest"
        return "Appalachian_Highlands"

    if state in {"MO", "MISSOURI", "AR", "ARKANSAS", "IL", "ILLINOIS", "IN", "INDIANA", "OH", "OHIO", "IA", "IOWA", "KS", "KANSAS", "NE", "NEBRASKA", "OK", "OKLAHOMA", "TX", "TEXAS"}:
        return "Interior_Prairie_Midwest"

    if lat is not None and lon is not None:
        if 24.0 <= lat <= 38.0 and -85.0 <= lon <= -75.0:
            return "Piedmont_Granite_Flatrocks"
        if 34.0 <= lat <= 45.0 and -84.0 <= lon <= -70.0:
            return "Appalachian_Highlands"
        if 28.0 <= lat <= 49.0 and -102.0 <= lon <= -84.0:
            return "Interior_Prairie_Midwest"
        if 25.0 <= lat <= 35.0 and -98.0 <= lon <= -80.0:
            return "Coastal_Plain_Sandhills"

    return "Other_US"


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

    # Regional Ecological Groups Breakdown
    region_counts = df["regional_group"].value_counts()
    logger.info("\n--- REGIONAL ECO-GEOGRAPHIC GROUPS ---")
    for reg, cnt in region_counts.items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {reg:<30} : {cnt:>5} records ({pct:>5.1f}%)")

    # Herbarium Institutions (Top 10)
    inst_counts = df["institutionCode"].value_counts()
    logger.info("\n--- TOP HERBARIUM INSTITUTIONS ---")
    for inst, cnt in inst_counts.head(10).items():
        pct = (cnt / total_records) * 100.0
        logger.info(f"  * {inst:<15} : {cnt:>5} records ({pct:>5.1f}%)")

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
    Consolidated botanical voucher harvester and curator for the Packera dubia species complex.
    Unifies GBIF query execution, Darwin Core metadata normalization, 3-tier determiner authority
    stratification, geographic filtering, and asynchronous high-resolution image acquisition.
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
    ):
        self.taxa = taxa or DEFAULT_TARGET_TAXA
        self.max_uncertainty_meters = max_uncertainty_meters
        self.max_records_per_taxon = max_records_per_taxon
        self.exclude_western = exclude_western
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

    def harvest(self) -> pd.DataFrame:
        """Harvests and normalizes Darwin Core occurrence records across configured taxa."""
        all_curated_records: List[Dict[str, Any]] = []
        seen_catalog_keys: set = set()

        for taxon in self.taxa:
            self.logger.info(f"Querying GBIF API for taxon: '{taxon}' (country=US, basisOfRecord=PRESERVED_SPECIMEN)...")
            offset = 0
            limit = 300
            taxon_harvested = 0
            taxon_retained = 0
            taxon_western_excluded = 0

            while taxon_retained < self.max_records_per_taxon:
                try:
                    response = occ.search(
                        scientificName=taxon,
                        country="US",
                        basisOfRecord="PRESERVED_SPECIMEN",
                        limit=limit,
                        offset=offset,
                    )
                except Exception as e:
                    self.logger.error(f"GBIF API query error for '{taxon}' at offset {offset}: {e}")
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
                    if self.exclude_western and is_excluded_western_region(state_prov, lat=lat_val, lon=lon_val):
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

                    # 5. Taxonomic Authority & Determiner Tier Stratification
                    determiner_tier, type_status, determiner_raw = parse_determiner_tier(
                        type_status_raw=rec.get("typeStatus"),
                        identified_by_raw=rec.get("identifiedBy"),
                        recorded_by_raw=rec.get("recordedBy"),
                        history_raw=rec.get("verbatimIdentificationHistory"),
                        institution_code_raw=inst_code,
                        locality_raw=rec.get("locality") or rec.get("verbatimLocality"),
                        habitat_raw=rec.get("habitat"),
                    )

                    # 6. Regional Ecological Group Assignment
                    regional_group = infer_regional_group(
                        lat=lat_val,
                        lon=lon_val,
                        state_province=state_prov,
                        habitat=rec.get("habitat"),
                        locality=rec.get("locality") or rec.get("verbatimLocality"),
                    )

                    curated_record = {
                        "catalogNumber": catalog_number,
                        "institutionCode": inst_code,
                        "scientificName": rec.get("scientificName") or rec.get("species") or taxon,
                        "species_raw": rec.get("scientificName") or rec.get("species") or taxon,
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
                    }
                    all_curated_records.append(curated_record)
                    taxon_retained += 1

                    if taxon_retained >= self.max_records_per_taxon:
                        break

                offset += len(results)
                if offset >= count:
                    break
                time.sleep(0.1)

            self.logger.info(
                f"Taxon '{taxon}': Processed {taxon_harvested} occurrences -> "
                f"Retained {taxon_retained} curated records meeting quality filters "
                f"(Western states excluded: {taxon_western_excluded})."
            )

        return pd.DataFrame(all_curated_records)

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
        stats["quality_rejected"] = quality_rejected
        if mp_values:
            stats["median_mp"] = float(np.median(mp_values))

        return df_filtered, stats

    def export(self, df: pd.DataFrame) -> Path:
        """Atomically persists the curated vouchers table to CSV."""
        return export_curated_table(df, self.output_csv, logger=self.logger)

    def run(
        self,
        download_images: bool = False,
        force: Optional[bool] = None,
    ) -> Tuple[pd.DataFrame, Optional[Dict[str, Any]]]:
        """Executes the end-to-end voucher ingestion and curation workflow."""
        if force is not None:
            self.force = force

        self.logger.info("Starting Packera Voucher Ingestion & Authority Stratification Pipeline...")
        self.logger.info(f"Target Taxa: {self.taxa}")
        self.logger.info(f"Max Coordinate Uncertainty Threshold: {self.max_uncertainty_meters} m")
        self.logger.info(f"Exclude Western States (> TX & OK): {self.exclude_western}")
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
) -> pd.DataFrame:
    """Backwards-compatible functional entry point for occurrence harvesting."""
    harvester = VoucherHarvester(
        taxa=taxa_list,
        max_uncertainty_meters=max_uncertainty_meters,
        max_records_per_taxon=max_records_per_taxon,
        exclude_western=exclude_western,
        logger=logger,
    )
    return harvester.harvest()


def main() -> None:
    """CLI execution wrapper delegating to data prep runner."""
    import importlib
    cli_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")
    cli_mod.main()


if __name__ == "__main__":
    main()
