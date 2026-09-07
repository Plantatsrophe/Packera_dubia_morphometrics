"""
===============================================================================
Script: prepare_lm2_dataset.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Production Pre-Processor and Dataset Staging Utility for LeafMachine2 (LM2).
    Coordinates herbarium sheet image symlink staging and voucher asset manifest
    generation without interactive labeling dependencies.
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
        description="Prepare LeafMachine2 dataset: Stage voucher images and generate manifest."
    )
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
    return parser.parse_args()


def run_pipeline(args: argparse.Namespace) -> None:
    """Executes the voucher dataset staging workflow."""
    logger.info("=" * 80)
    logger.info("        LEAFMACHINE2 DATASET PREPARATION & STAGING PIPELINE        ")
    logger.info("=" * 80)

    dirs = setup_lm2_directories(args.lm2_root)

    # Stage voucher images via symlinks
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

    logger.info("\n" + "=" * 80)
    logger.info("LeafMachine2 dataset staging completed successfully.")
    logger.info("=" * 80)


if __name__ == "__main__":
    cli_args = parse_args()
    run_pipeline(cli_args)

