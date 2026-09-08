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

from scripts.core.config import PipelineConfig
from scripts.core.harvester import VoucherHarvester, setup_logger


def build_cli_parser() -> argparse.ArgumentParser:
    """Constructs command-line argument parser for the voucher harvester."""
    cfg = PipelineConfig.from_yaml()
    parser = argparse.ArgumentParser(
        description="Automated GBIF Voucher Harvester & Determiner Authority Scorer for Packera dubia complex."
    )
    parser.add_argument(
        "--taxon",
        "--taxa",
        dest="taxa",
        nargs="+",
        default=cfg.taxa.target_species,
        help="One or more scientific binomials to harvest from GBIF (e.g. 'Packera dubia').",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=str(cfg.paths.raw_vouchers_dir),
        help=f"Directory for storing raw downloaded voucher images (default: {cfg.paths.raw_vouchers_dir}).",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=str(cfg.paths.curated_vouchers_csv),
        help=f"Path for curated metadata CSV output (default: {cfg.paths.curated_vouchers_csv}).",
    )
    parser.add_argument(
        "--max-records",
        "--max-records-per-taxon",
        dest="max_records",
        type=int,
        default=cfg.harvesting.max_records_per_taxon,
        help=f"Maximum records to retain per taxon query (default: {cfg.harvesting.max_records_per_taxon}).",
    )
    parser.add_argument(
        "--max-uncertainty",
        type=float,
        default=cfg.harvesting.max_uncertainty_meters,
        help=f"Maximum allowed georeferencing uncertainty in meters (default: {cfg.harvesting.max_uncertainty_meters}).",
    )
    parser.add_argument(
        "--exclude-western",
        action=argparse.BooleanOptionalAction,
        default=cfg.harvesting.exclude_western,
        help=f"Exclude records from western US states (> TX & OK) (default: {cfg.harvesting.exclude_western}).",
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
        default=cfg.harvesting.download_concurrency,
        help=f"Maximum concurrent asynchronous image downloads (default: {cfg.harvesting.download_concurrency}).",
    )
    parser.add_argument(
        "--min-megapixels",
        type=float,
        default=cfg.thresholds.min_megapixels,
        help=f"Minimum resolution threshold in Megapixels (default: {cfg.thresholds.min_megapixels}).",
    )
    parser.add_argument(
        "--min-file-size-kb",
        type=float,
        default=cfg.thresholds.min_file_size_kb,
        help=f"Minimum image file size in KB (default: {cfg.thresholds.min_file_size_kb}).",
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
        default=cfg.thresholds.min_sharpness_laplacian,
        help=f"Minimum Laplacian variance score for sharpness (default: {cfg.thresholds.min_sharpness_laplacian}).",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default=str(cfg.paths.summary_log),
        help=f"Path for summary log file (default: {cfg.paths.summary_log}).",
    )
    parser.add_argument(
        "--force",
        "--overwrite",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-download and overwrite existing cached voucher images.",
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
    cfg = PipelineConfig.from_yaml()
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
        workspace_dir=cfg.paths.workspace_root,
        logger=logger,
        force=args.force,
    )

    harvester.run(download_images=args.download_images)


if __name__ == "__main__":
    main()
