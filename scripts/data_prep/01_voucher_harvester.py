#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    CLI entry point to execute the botanical voucher harvesting pipeline.
    Queries GBIF occurrences (dynamically across all Packera species genus-wide
    by default, or strictly for user-specified --taxon binomials), downloads
    high-resolution herbarium specimen imagery, applies 3-tier determiner
    authority scoring and geographic exclusion filters, and atomically
    generates curated datasets.
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

from scripts.core.config import PipelineConfig, resolve_max_workers
from scripts.core.harvester import VoucherHarvester
from scripts.core.taxonomy_authority import (
    resolve_geographic_scope_settings,
    setup_logger,
)


def build_cli_parser() -> argparse.ArgumentParser:
    """Constructs command-line argument parser for the voucher harvester."""
    cfg = PipelineConfig.from_yaml()
    parser = argparse.ArgumentParser(
        description="Automated GBIF Voucher Harvester & Determiner Authority Scorer for genus Packera."
    )
    parser.add_argument(
        "--taxon",
        dest="taxa",
        nargs="+",
        default=list(cfg.taxa.target_species),
        help=(
            "One or more scientific binomials to harvest strictly from GBIF. "
            "When omitted (empty), dynamically harvests all Packera species across the genus."
        ),
    )
    parser.add_argument(
        "--taxa",
        dest="taxa",
        nargs="+",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--expand-synonyms",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Dynamically unify unaccepted Packera synonyms and historical Senecio/Cineraria "
            "basionyms under accepted Packera species via GBIF Backbone taxonomy (default: True)."
        ),
    )
    parser.add_argument(
        "--geographic-scope",
        choices=["eastern_na", "north_america", "global"],
        default=cfg.harvesting.geographic_scope,
        help=f"Geographic scope for GBIF harvesting and longitudinal filtering (default: {cfg.harvesting.geographic_scope}).",
    )
    parser.add_argument(
        "--include-western",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--countries",
        nargs="+",
        default=None,
        help="ISO-2 country codes to query in GBIF (e.g., --countries US CA MX). Defaults to config countries for eastern_na or US CA MX for north_america.",
    )
    parser.add_argument(
        "--home-herbarium",
        type=str,
        default=cfg.harvesting.home_herbarium,
        help=f"Home herbarium institution code prioritized during duplicate exsiccatae selection (default: {cfg.harvesting.home_herbarium}).",
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
        dest="max_records",
        type=int,
        default=cfg.harvesting.max_records_per_taxon,
        help=f"Maximum records to retain per taxon query (default: {cfg.harvesting.max_records_per_taxon}).",
    )
    parser.add_argument(
        "--max-records-per-taxon",
        dest="max_records",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
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
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--download-images",
        action="store_true",
        default=False,
        help="Flag to enable asynchronous download of high-resolution specimen images.",
    )
    default_conc = resolve_max_workers(default=cfg.harvesting.download_concurrency, cap=32, cfg=cfg)
    parser.add_argument(
        "--concurrency",
        type=int,
        default=default_conc,
        help=f"Maximum concurrent asynchronous image downloads (default: {default_conc}, auto-scaled up to 32 from --threads).",
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
        action=argparse.BooleanOptionalAction,
        default=cfg.harvesting.check_sharpness,
        help=f"Enable Laplacian variance edge sharpness evaluation (default: {cfg.harvesting.check_sharpness}).",
    )
    parser.add_argument(
        "--harvest-full-continental-range",
        action=argparse.BooleanOptionalAction,
        default=cfg.harvesting.harvest_full_continental_range,
        help=(
            "After passing the Eastern NA species residency screen, harvest the full US + CA "
            f"continental range for resident species (default: {cfg.harvesting.harvest_full_continental_range})."
        ),
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
        dest="force",
        action="store_true",
        default=False,
        help="Force re-download and overwrite existing cached voucher images.",
    )
    parser.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--usda-county-filter",
        action=argparse.BooleanOptionalAction,
        default=cfg.harvesting.usda_county_filter,
        help=(
            "Drop non-exempt US vouchers outside USDA PLANTS documented counties + 1-ring/buffer "
            "when enabled; when False, retains all vouchers and records soft usda_county_status in SQLite "
            f"(default: {cfg.harvesting.usda_county_filter})."
        ),
    )
    parser.add_argument(
        "--usda-county-buffer-km",
        type=float,
        default=cfg.harvesting.usda_county_buffer_km,
        help=(
            "Geodesic buffer distance in km around USDA PLANTS documented counties "
            f"(default: {cfg.harvesting.usda_county_buffer_km})."
        ),
    )
    parser.add_argument(
        "--usda-exempt-gold",
        action=argparse.BooleanOptionalAction,
        default=cfg.harvesting.usda_exempt_gold,
        help=(
            "Exempt Tier_1_Gold specialist determinations from USDA county filtering "
            f"(nomenclatural types are always exempt; default: {cfg.harvesting.usda_exempt_gold})."
        ),
    )
    parser.add_argument(
        "--include-ungeoreferenced-state-vouchers",
        action=argparse.BooleanOptionalAction,
        default=cfg.harvesting.include_ungeoreferenced_state_vouchers,
        help=(
            "Rescue ungeoreferenced (hasCoordinate=False) cis-Mississippi / eastern state herbarium "
            f"sheets for morphological and authority curation (default: {cfg.harvesting.include_ungeoreferenced_state_vouchers})."
        ),
    )
    parser.add_argument(
        "--db",
        type=str,
        default=None,
        help=f"Optional SQLite database path override (default: {cfg.paths.db_path}).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose debug logging.",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for 01_voucher_harvester.py (defaults to empty argv when imported under unittest runner)."""
    if argv is None and len(sys.argv) > 0 and Path(sys.argv[0]).name in {"run_audit.py", "unittest", "__main__.py"}:
        argv = []
    return build_cli_parser().parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """CLI execution entrypoint."""
    cfg = PipelineConfig.from_yaml()
    args = parse_args(argv)

    log_file_path = Path(args.log_file)
    logger = setup_logger(log_file_path=log_file_path, verbose=args.verbose)

    target_taxa = list(args.taxa) if args.taxa else []

    explicit_countries = getattr(args, "countries", None)
    arg_scope = getattr(args, "geographic_scope", cfg.harvesting.geographic_scope)
    arg_include_western = bool(getattr(args, "include_western", False))
    if (
        explicit_countries is None
        and arg_scope == cfg.harvesting.geographic_scope
        and not arg_include_western
        and arg_scope == "eastern_na"
    ):
        explicit_countries = cfg.harvesting.countries

    eff_scope, eff_exclude_western, eff_countries = resolve_geographic_scope_settings(
        geographic_scope=arg_scope,
        exclude_western=args.exclude_western,
        include_western=arg_include_western,
        countries=explicit_countries,
    )

    harvester = VoucherHarvester(
        taxa=target_taxa,
        max_uncertainty_meters=args.max_uncertainty,
        max_records_per_taxon=args.max_records,
        exclude_western=eff_exclude_western,
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
        geographic_scope=eff_scope,
        countries=eff_countries,
        include_western=arg_include_western,
        western_longitude_threshold=cfg.thresholds.western_longitude_threshold,
        expand_synonyms=bool(getattr(args, "expand_synonyms", True)),
        cfg_synonyms=cfg.taxa.synonyms,
        home_herbarium=getattr(args, "home_herbarium", cfg.harvesting.home_herbarium),
        db_path=getattr(args, "db", None),
        usda_county_filter=bool(getattr(args, "usda_county_filter", cfg.harvesting.usda_county_filter)),
        usda_county_buffer_km=float(getattr(args, "usda_county_buffer_km", cfg.harvesting.usda_county_buffer_km)),
        usda_exempt_gold=bool(getattr(args, "usda_exempt_gold", cfg.harvesting.usda_exempt_gold)),
        harvest_full_continental_range=bool(
            getattr(args, "harvest_full_continental_range", cfg.harvesting.harvest_full_continental_range)
        ),
        min_eastern_residency_records=int(cfg.harvesting.min_eastern_residency_records),
        cis_mississippi_longitude_threshold=float(cfg.harvesting.cis_mississippi_longitude_threshold),
        min_cis_mississippi_fraction=float(cfg.harvesting.min_cis_mississippi_fraction),
        include_ungeoreferenced_state_vouchers=bool(
            getattr(
                args,
                "include_ungeoreferenced_state_vouchers",
                cfg.harvesting.include_ungeoreferenced_state_vouchers,
            )
        ),
        webp_quality=int(cfg.harvesting.webp_quality),
    )

    harvester.run(download_images=args.download_images)


if __name__ == "__main__":
    main()
