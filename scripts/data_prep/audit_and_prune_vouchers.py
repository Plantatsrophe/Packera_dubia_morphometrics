#!/usr/bin/env python3
"""
===============================================================================
Script: audit_and_prune_vouchers.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Audits the optical resolution, file size, and image quality metrics across
    harvested botanical voucher sheets in data/raw_vouchers/. Flags or prunes
    substandard (< 8 MP or blurry/upscaled) sheets and synchronizes the curated
    vouchers metadata table (data/tables/curated_vouchers.csv).
===============================================================================
"""

import sys
import shutil
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    DEFAULT_WORKSPACE,
    DEFAULT_RAW_DIR,
    DEFAULT_CURATED_CSV,
    DEFAULT_QUARANTINE_DIR,
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_FILE_SIZE_KB,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN,
)
from scripts.core.harvester_utils import (
    setup_logger,
    validate_image_quality,
)

Image.MAX_IMAGE_PIXELS = None  # Allow decompression of high-resolution botanical sheets (>89 MP)


def audit_and_prune_dataset(
    raw_dir: Path,
    csv_path: Path,
    quarantine_dir: Path,
    min_megapixels: float = DEFAULT_MIN_MEGAPIXELS,
    min_file_size_kb: float = DEFAULT_MIN_FILE_SIZE_KB,
    check_sharpness: bool = False,
    min_sharpness: float = DEFAULT_MIN_SHARPNESS_LAPLACIAN,
    mode: str = "dry-run",
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    """
    Scans local voucher sheets, computes optical metrics, and optionally prunes
    substandard files while synchronizing curated_vouchers.csv.
    
    Args:
        raw_dir: Directory containing raw voucher images.
        csv_path: Path to curated metadata CSV table.
        quarantine_dir: Destination directory for quarantined substandard files.
        min_megapixels: Minimum resolution threshold in Megapixels.
        min_file_size_kb: Minimum compressed file size in KB.
        check_sharpness: Whether to calculate Laplacian variance sharpness.
        min_sharpness: Minimum acceptable Laplacian variance score.
        mode: Operation mode: 'dry-run', 'quarantine', or 'delete'.
        logger: Logger instance.
        
    Returns:
        Dict[str, Any]: Detailed audit summary dictionary.
    """
    if logger is None:
        logger = logging.getLogger("VoucherAuditor")

    if not raw_dir.exists():
        logger.error(f"Raw vouchers directory not found at: {raw_dir}")
        return {"error": "raw_dir_not_found"}

    df = pd.read_csv(csv_path) if csv_path.exists() else pd.DataFrame()
    inst_by_cat = dict(zip(df["catalogNumber"], df["institutionCode"])) if not df.empty and "catalogNumber" in df.columns else {}

    image_files = sorted(list(raw_dir.glob("*.jpg")) + list(raw_dir.glob("*.jpeg")) + list(raw_dir.glob("*.png")))
    total_images = len(image_files)

    logger.info(f"Auditing {total_images:,} voucher images in {raw_dir}...")
    logger.info(f"Quality Criteria: min_megapixels={min_megapixels} MP, min_file_size_kb={min_file_size_kb} KB, check_sharpness={check_sharpness}")

    accepted_catalogs = set()
    rejected_records = []
    herbarium_stats = defaultdict(lambda: {"count": 0, "accepted": 0, "rejected": 0, "mp_list": []})

    for img_path in image_files:
        cat_num = img_path.stem
        inst = inst_by_cat.get(cat_num, "UNKNOWN")
        herbarium_stats[inst]["count"] += 1

        is_valid, metrics = validate_image_quality(
            image_path=img_path,
            min_megapixels=min_megapixels,
            min_file_size_kb=min_file_size_kb,
            check_sharpness=check_sharpness,
            min_sharpness=min_sharpness,
        )

        mp = metrics.get("megapixels", 0.0)
        if mp > 0:
            herbarium_stats[inst]["mp_list"].append(mp)

        if is_valid:
            accepted_catalogs.add(cat_num)
            herbarium_stats[inst]["accepted"] += 1
        else:
            herbarium_stats[inst]["rejected"] += 1
            rejected_records.append({
                "catalogNumber": cat_num,
                "institutionCode": inst,
                "path": img_path,
                "reason": metrics.get("reason", "unknown"),
                "megapixels": mp,
                "file_size_kb": metrics.get("file_size_kb", 0.0),
                "sharpness": metrics.get("sharpness"),
            })

    total_accepted = len(accepted_catalogs)
    total_rejected = len(rejected_records)

    # Output detailed report
    logger.info("\n" + "=" * 80)
    logger.info("                  HERBARIUM IMAGE QUALITY AUDIT BREAKDOWN                  ")
    logger.info("=" * 80)
    header = f"{'Institution':<16} | {'Total':>7} | {'Accepted':>8} | {'Rejected':>8} | {'Median MP':>10} | {'Min MP':>8} | {'Max MP':>8}"
    logger.info(header)
    logger.info("-" * 80)

    for inst, s in sorted(herbarium_stats.items(), key=lambda x: x[1]["count"], reverse=True):
        if s["count"] == 0:
            continue
        med_mp = sorted(s["mp_list"])[len(s["mp_list"])//2] if s["mp_list"] else 0.0
        min_mp = min(s["mp_list"]) if s["mp_list"] else 0.0
        max_mp = max(s["mp_list"]) if s["mp_list"] else 0.0
        line = f"{inst:<16} | {s['count']:>7} | {s['accepted']:>8} | {s['rejected']:>8} | {med_mp:>10.2f} | {min_mp:>8.2f} | {max_mp:>8.2f}"
        logger.info(line)

    logger.info("-" * 80)
    logger.info(f"Total Vouchers Audited : {total_images:>7,}")
    logger.info(f"Accepted (High Quality): {total_accepted:>7,} ({(total_accepted/total_images*100.0 if total_images else 0):>5.1f}%)")
    logger.info(f"Rejected (Substandard) : {total_rejected:>7,} ({(total_rejected/total_images*100.0 if total_images else 0):>5.1f}%)")
    logger.info("=" * 80)

    # Perform Pruning / Quarantining Actions if requested
    if mode in {"quarantine", "delete"} and total_rejected > 0:
        if mode == "quarantine":
            quarantine_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"\nMoving {total_rejected} substandard voucher images to quarantine: {quarantine_dir}...")
            for rec in rejected_records:
                src_file = rec["path"]
                dest_file = quarantine_dir / src_file.name
                shutil.move(str(src_file), str(dest_file))
        elif mode == "delete":
            logger.info(f"\nPermanently deleting {total_rejected} substandard voucher images from disk...")
            for rec in rejected_records:
                src_file = rec["path"]
                if src_file.exists():
                    src_file.unlink()

        # Synchronize curated_vouchers.csv table
        if not df.empty and "catalogNumber" in df.columns:
            original_len = len(df)
            df_synced = df[df["catalogNumber"].astype(str).isin(accepted_catalogs)].copy().reset_index(drop=True)
            df_synced.to_csv(csv_path, index=False, encoding="utf-8")
            logger.info(
                f"Synchronized metadata table {csv_path}: "
                f"Updated record count {original_len:,} -> {len(df_synced):,} (Removed {original_len - len(df_synced):,} pruned entries)."
            )

    return {
        "total_images": total_images,
        "total_accepted": total_accepted,
        "total_rejected": total_rejected,
        "herbarium_stats": dict(herbarium_stats),
        "rejected_records": rejected_records,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Audit and prune low-resolution or poor-quality specimen images from Packera dataset."
    )
    parser.add_argument(
        "--raw-dir",
        type=str,
        default=str(DEFAULT_RAW_DIR),
        help=f"Directory containing raw voucher image files (default: {DEFAULT_RAW_DIR})."
    )
    parser.add_argument(
        "--csv-path",
        type=str,
        default=str(DEFAULT_CURATED_CSV),
        help=f"Path to curated metadata CSV (default: {DEFAULT_CURATED_CSV})."
    )
    parser.add_argument(
        "--quarantine-dir",
        type=str,
        default=str(DEFAULT_QUARANTINE_DIR),
        help=f"Destination directory for quarantined images (default: {DEFAULT_QUARANTINE_DIR})."
    )
    parser.add_argument(
        "--min-megapixels",
        type=float,
        default=DEFAULT_MIN_MEGAPIXELS,
        help=f"Minimum resolution threshold in Megapixels (default: {DEFAULT_MIN_MEGAPIXELS})."
    )
    parser.add_argument(
        "--min-file-size-kb",
        type=float,
        default=DEFAULT_MIN_FILE_SIZE_KB,
        help=f"Minimum image file size in KB (default: {DEFAULT_MIN_FILE_SIZE_KB})."
    )
    parser.add_argument(
        "--check-sharpness",
        action="store_true",
        default=False,
        help="Enable Laplacian variance edge sharpness evaluation."
    )
    parser.add_argument(
        "--min-sharpness",
        type=float,
        default=DEFAULT_MIN_SHARPNESS_LAPLACIAN,
        help=f"Minimum Laplacian variance threshold (default: {DEFAULT_MIN_SHARPNESS_LAPLACIAN})."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Perform dry-run audit report without modifying files or CSV."
    )
    parser.add_argument(
        "--quarantine",
        action="store_true",
        default=False,
        help="Move substandard images to quarantine folder and update CSV."
    )
    parser.add_argument(
        "--delete",
        action="store_true",
        default=False,
        help="Permanently delete substandard images from disk and update CSV."
    )

    args = parser.parse_args()

    mode = "dry-run"
    if args.quarantine:
        mode = "quarantine"
    elif args.delete:
        mode = "delete"
    elif args.dry_run:
        mode = "dry-run"

    logger = setup_logger()
    logger.info(f"Running Voucher Image Audit & Pruning Utility (Mode: {mode.upper()})...")

    audit_and_prune_dataset(
        raw_dir=Path(args.raw_dir),
        csv_path=Path(args.csv_path),
        quarantine_dir=Path(args.quarantine_dir),
        min_megapixels=args.min_megapixels,
        min_file_size_kb=args.min_file_size_kb,
        check_sharpness=args.check_sharpness,
        min_sharpness=args.min_sharpness,
        mode=mode,
        logger=logger,
    )


if __name__ == "__main__":
    main()
