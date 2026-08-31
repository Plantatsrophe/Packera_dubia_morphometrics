#!/usr/bin/env python3
"""
===============================================================================
Script: prepare_lm2_dataset.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Dual-purpose Pre-Processor and Dataset Generator for LeafMachine2 (LM2).
    Coordinates SAM 2 binary mask conversion to COCO 1.0 JSON polygon format
    and herbarium sheet image staging / symlink creation.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

from scripts.data_prep.coco_exporter import (
    OMITTED_CLASSES,
    PCD_CLASS_MAPPING,
    PCD_COCO_CATEGORIES,
    convert_masks_to_coco_dataset,
    save_coco_json,
    split_coco_dataset,
)
from scripts.data_prep.staging_utils import (
    setup_lm2_directories,
    stage_voucher_symlinks,
    write_manifest_csv,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("PrepareLM2Dataset")


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for LeafMachine2 dataset preparation."""
    parser = argparse.ArgumentParser(
        description="Prepare LeafMachine2 dataset: Stage voucher images and export COCO annotations."
    )
    # Staging options
    parser.add_argument(
        "--lm2-root",
        type=Path,
        default=Path("LM2_Project"),
        help="Root directory for LeafMachine2 project (default: LM2_Project)",
    )
    parser.add_argument(
        "--input-dirs",
        type=Path,
        nargs="+",
        default=[Path("data/raw_vouchers")],
        help="Input directories containing raw voucher images (default: data/raw_vouchers)",
    )
    parser.add_argument(
        "--skip-symlinks",
        action="store_true",
        help="Skip staging voucher image symlinks",
    )
    parser.add_argument(
        "--manifest-out",
        type=Path,
        default=Path("outputs/reports/staged_images_manifest.csv"),
        help="Destination path for staged image manifest CSV",
    )

    # COCO export options
    parser.add_argument(
        "--skip-coco",
        action="store_true",
        help="Skip SAM 2 mask conversion to COCO JSON format",
    )
    parser.add_argument(
        "--masks-dir",
        type=Path,
        default=Path("data/raw_annotations/masks"),
        help="Directory containing binary PNG instance masks (default: data/raw_annotations/masks)",
    )
    parser.add_argument(
        "--coco-output",
        type=Path,
        default=Path("LM2_Project/Data/annotations/coco_pcd_packera.json"),
        help="Output destination for primary COCO dataset JSON",
    )
    parser.add_argument(
        "--val-split",
        type=float,
        default=0.0,
        help="Fraction of dataset to partition into a validation split (e.g. 0.20 for 80/20)",
    )
    parser.add_argument(
        "--min-area",
        type=float,
        default=50.0,
        help="Minimum contour polygon area threshold in pixels (default: 50.0)",
    )

    return parser.parse_args()


def run_pipeline(args: argparse.Namespace) -> None:
    """Executes the dataset staging and annotation compilation workflow."""
    logger.info("=" * 80)
    logger.info("        LEAFMACHINE2 DATASET PREPARATION & STAGING PIPELINE        ")
    logger.info("=" * 80)

    dirs = setup_lm2_directories(args.lm2_root)

    # Step 1: Stage voucher images via symlinks
    if not args.skip_symlinks:
        logger.info("\n--- Phase 1: Staging Herbarium Voucher Imagery ---")
        created, skipped, records = stage_voucher_symlinks(
            source_dirs=args.input_dirs,
            target_images_dir=dirs["images"],
            use_relative_symlinks=True,
        )
        if records and args.manifest_out:
            write_manifest_csv(records, args.manifest_out)
    else:
        logger.info("Skipping voucher image symlink staging (--skip-symlinks requested).")

    # Step 2: Convert SAM 2 binary masks to COCO 1.0 JSON format
    if not args.skip_coco:
        logger.info("\n--- Phase 2: Converting SAM 2 Binary Masks to COCO 1.0 JSON ---")
        coco_data = convert_masks_to_coco_dataset(
            masks_dir=args.masks_dir,
            min_area_px=args.min_area,
        )

        if args.val_split > 0.0:
            train_doc, val_doc = split_coco_dataset(coco_data, val_ratio=args.val_split)
            stem = args.coco_output.stem
            train_path = args.coco_output.with_name(f"{stem}_train.json")
            val_path = args.coco_output.with_name(f"{stem}_val.json")

            save_coco_json(train_doc, train_path)
            save_coco_json(val_doc, val_path)
            save_coco_json(coco_data, args.coco_output)
            logger.info(f"Generated train split ({len(train_doc['images'])} images) -> {train_path}")
            logger.info(f"Generated val split ({len(val_doc['images'])} images) -> {val_path}")
        else:
            save_coco_json(coco_data, args.coco_output)

        summary_file = Path("outputs/reports/dataset_build_summary.json")
        summary_file.parent.mkdir(parents=True, exist_ok=True)
        summary = {
            "total_images": len(coco_data["images"]),
            "total_annotations": len(coco_data["annotations"]),
            "categories": [c["name"] for c in coco_data["categories"]],
            "masks_source_dir": str(args.masks_dir),
            "coco_output_path": str(args.coco_output),
        }
        import json
        with open(summary_file, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
    else:
        logger.info("Skipping COCO dataset compilation (--skip-coco requested).")

    logger.info("\n" + "=" * 80)
    logger.info("LeafMachine2 dataset staging & preparation completed successfully.")
    logger.info("=" * 80)


if __name__ == "__main__":
    cli_args = parse_args()
    run_pipeline(cli_args)
