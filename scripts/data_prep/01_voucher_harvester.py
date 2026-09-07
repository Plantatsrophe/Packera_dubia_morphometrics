#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    CLI entry point to execute the botanical voucher harvesting pipeline.
    Queries GBIF occurrences, downloads high-resolution herbarium specimen
    imagery, applies 3-tier determiner authority scoring and geographic
    exclusion filters, and atomically generates curated datasets.
===============================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add project root directory to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    DEFAULT_MIN_FILE_SIZE_KB,
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_MIN_SHARPNESS_LAPLACIAN,
    DEFAULT_OUTPUT_CSV,
    DEFAULT_RAW_DIR,
    DEFAULT_SUMMARY_LOG,
    DEFAULT_TARGET_TAXA,
    DEFAULT_WORKSPACE,
)
from scripts.core.harvester import VoucherHarvester, setup_logger


def build_cli_parser() -> argparse.ArgumentParser:
    """Constructs command-line argument parser for the voucher harvester."""
    parser = argparse.ArgumentParser(
        description="Automated GBIF Voucher Harvester & Determiner Authority Scorer for Packera dubia complex."
    )
    parser.add_argument(
        "--taxon",
        "--taxa",
        dest="taxa",
        nargs="+",
        default=DEFAULT_TARGET_TAXA,
        help="One or more scientific binomials to harvest from GBIF (e.g. 'Packera dubia').",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=str(DEFAULT_RAW_DIR),
        help=f"Directory for storing raw downloaded voucher images (default: {DEFAULT_RAW_DIR}).",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=str(DEFAULT_OUTPUT_CSV),
        help=f"Path for curated metadata CSV output (default: {DEFAULT_OUTPUT_CSV}).",
    )
    parser.add_argument(
        "--max-records",
        "--max-records-per-taxon",
        dest="max_records",
        type=int,
        default=5000,
        help="Maximum records to retain per taxon query (default: 5000).",
    )
    parser.add_argument(
        "--max-uncertainty",
        type=float,
        default=5000.0,
        help="Maximum allowed georeferencing uncertainty in meters (default: 5000.0).",
    )
    parser.add_argument(
        "--exclude-western",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude records from western US states (> TX & OK) (default: True).",
    )
    parser.add_argument(
        "--download-images",
        action="store_true",
        default=False,
        help="Flag to enable asynchronous download of high-resolution specimen images.",
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=15,
        help="Maximum concurrent asynchronous image downloads (default: 15).",
    )
    parser.add_argument(
        "--min-megapixels",
        type=float,
        default=DEFAULT_MIN_MEGAPIXELS,
        help=f"Minimum resolution threshold in Megapixels (default: {DEFAULT_MIN_MEGAPIXELS}).",
    )
    parser.add_argument(
        "--min-file-size-kb",
        type=float,
        default=DEFAULT_MIN_FILE_SIZE_KB,
        help=f"Minimum image file size in KB (default: {DEFAULT_MIN_FILE_SIZE_KB}).",
    )
    parser.add_argument(
        "--check-sharpness",
        action="store_true",
        default=False,
        help="Enable Laplacian variance edge sharpness evaluation.",
    )
    parser.add_argument(
        "--min-sharpness",
        type=float,
        default=DEFAULT_MIN_SHARPNESS_LAPLACIAN,
        help=f"Minimum Laplacian variance score for sharpness (default: {DEFAULT_MIN_SHARPNESS_LAPLACIAN}).",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default=str(DEFAULT_SUMMARY_LOG),
        help=f"Path for summary log file (default: {DEFAULT_SUMMARY_LOG}).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose debug logging.",
    )
    return parser


def main() -> None:
    """CLI execution entrypoint."""
    parser = build_cli_parser()
    args = parser.parse_args()

    log_file_path = Path(args.log_file)
    logger = setup_logger(log_file_path=log_file_path, verbose=args.verbose)

    harvester = VoucherHarvester(
        taxa=args.taxa,
        max_uncertainty_meters=args.max_uncertainty,
        max_records_per_taxon=args.max_records,
        exclude_western=args.exclude_western,
        min_megapixels=args.min_megapixels,
        min_file_size_kb=args.min_file_size_kb,
        check_sharpness=args.check_sharpness,
        min_sharpness=args.min_sharpness,
        concurrency=args.concurrency,
        output_csv=Path(args.output_csv),
        raw_dir=Path(args.out_dir),
        workspace_dir=DEFAULT_WORKSPACE,
        logger=logger,
    )

    harvester.run(download_images=args.download_images)


if __name__ == "__main__":
    main()
