"""
===============================================================================
Module: harvester_utils.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unified utility facade for botanical voucher harvesting, re-exporting metadata
    evaluation routines, media processing, and diagnostic logging summaries.
===============================================================================
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

# Re-export metadata processing functions for 100% backwards compatibility
from scripts.core.harvester_metadata import (
    calculate_circular_phenology,
    infer_regional_group,
    is_excluded_western_region,
    parse_determiner_tier,
    sanitize_filename,
)

# Re-export media handling functions for 100% backwards compatibility
from scripts.core.harvester_media import (
    download_all_voucher_images,
    download_single_image,
    extract_high_res_image_url,
    optimize_herbarium_image_url,
    validate_image_quality,
)

__all__ = [
    "sanitize_filename",
    "parse_determiner_tier",
    "calculate_circular_phenology",
    "is_excluded_western_region",
    "infer_regional_group",
    "optimize_herbarium_image_url",
    "extract_high_res_image_url",
    "validate_image_quality",
    "download_single_image",
    "download_all_voucher_images",
    "setup_logger",
    "print_and_log_summary",
]


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

    formatter = logging.Formatter(
        fmt="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_file_path:
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file_path, mode="w", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


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
