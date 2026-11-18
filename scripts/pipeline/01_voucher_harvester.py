#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py (Pipeline Entrypoint)
Project: Plantmatics Integrative Morphometrics and Phylogenomics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Pipeline CLI entrypoint re-exporting `build_cli_parser`, `parse_args`,
    `main`, `PipelineConfig`, `VoucherHarvester`, `standardize_taxon_name`, and
    `standardize_packera_taxon`, wiring `--portal-sources`, `--symbiota-api-url`,
    and `--infraspecific-policy` into multi-portal voucher ingestion.
===============================================================================
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path
from typing import Any, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core import config as _core_config
from scripts.core.config import PipelineConfig
from scripts.core.harvester import (
    VoucherHarvester,
    standardize_packera_taxon,
    standardize_taxon_name,
)

_data_prep_harvester = importlib.import_module("scripts.data_prep.01_voucher_harvester")


def _resolve_portal_cli_defaults(cfg: PipelineConfig) -> Tuple[List[str], List[str], str]:
    """Resolve portal sources, Symbiota endpoint URLs, and infraspecific policy defaults."""
    raw_portals = getattr(cfg.harvesting, "portal_sources", ["gbif"])
    default_portals = [str(p).strip() for p in (raw_portals or ["gbif"]) if str(p).strip()] or ["gbif"]
    raw_sym = getattr(cfg.harvesting, "symbiota_api_urls", [])
    default_sym = [str(u).strip() for u in (raw_sym or []) if str(u).strip()]
    fn = getattr(_core_config, "resolve_infraspecific_rank_policy", None)
    policy = "collapse_to_species"
    if callable(fn):
        try:
            cand = str(fn(cfg)).strip().lower()
            if cand in {"collapse_to_species", "preserve_accepted_trinomials"}:
                policy = cand
        except Exception:
            pass
    else:
        env_pol = (
            os.environ.get("PLANTMATICS_INFRASPECIFIC_POLICY", "").strip().lower()
            or os.environ.get("PLANTMATICS_INFRASPECIFIC_RANK_POLICY", "").strip().lower()
        )
        if env_pol in {"collapse_to_species", "preserve_accepted_trinomials"}:
            policy = env_pol
    return default_portals, default_sym, policy


def build_cli_parser() -> argparse.ArgumentParser:
    """Construct command-line argument parser including M1B multi-portal and trinomial flags."""
    cfg = PipelineConfig.from_yaml()
    parser = _data_prep_harvester.build_cli_parser()
    existing_opts = {opt for act in parser._actions for opt in act.option_strings}
    default_portals, default_sym, default_policy = _resolve_portal_cli_defaults(cfg)

    if "--portal-sources" not in existing_opts:
        parser.add_argument(
            "--portal-sources",
            nargs="+",
            default=default_portals,
            help=(
                "Occurrence portals to query ('gbif', 'idigbio', 'symbiota'). "
                f"Default: {default_portals}."
            ),
        )
    if "--symbiota-api-url" not in existing_opts:
        parser.add_argument(
            "--symbiota-api-url",
            dest="symbiota_api_urls",
            nargs="+",
            default=default_sym,
            help="One or more Symbiota/Symbiota2 occurrence REST endpoint URLs.",
        )
    if "--symbiota-api-urls" not in existing_opts:
        parser.add_argument(
            "--symbiota-api-urls",
            dest="symbiota_api_urls",
            nargs="+",
            default=argparse.SUPPRESS,
            help=argparse.SUPPRESS,
        )
    if "--infraspecific-policy" not in existing_opts:
        parser.add_argument(
            "--infraspecific-policy",
            choices=["collapse_to_species", "preserve_accepted_trinomials"],
            default=default_policy,
            help=(
                "Taxonomic handling of accepted subsp./var. trinomials "
                f"(default: {default_policy})."
            ),
        )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for 01_voucher_harvester.py."""
    if argv is None and len(sys.argv) > 0 and Path(sys.argv[0]).name in {"run_audit.py", "unittest", "__main__.py"}:
        argv = []
    return build_cli_parser().parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """CLI execution entrypoint forwarding portal and infraspecific policy configuration."""
    cfg = PipelineConfig.from_yaml()
    args = parse_args(argv)

    logger = _data_prep_harvester.setup_logger(log_file_path=Path(args.log_file), verbose=args.verbose)
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

    eff_scope, eff_exclude_western, eff_countries = _data_prep_harvester.resolve_geographic_scope_settings(
        geographic_scope=arg_scope,
        exclude_western=args.exclude_western,
        include_western=arg_include_western,
        countries=explicit_countries,
    )
    default_portals, default_sym, default_policy = _resolve_portal_cli_defaults(cfg)

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
        webp_quality=int(cfg.harvesting.webp_quality),
        portal_sources=getattr(args, "portal_sources", default_portals),
        symbiota_api_urls=getattr(args, "symbiota_api_urls", default_sym),
        infraspecific_policy=getattr(args, "infraspecific_policy", default_policy),
    )
    harvester.run(download_images=args.download_images)


if __name__ == "__main__":
    main()
