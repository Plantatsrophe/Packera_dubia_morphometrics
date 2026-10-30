#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py (Pipeline Alias)
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Pipeline module entry point re-exporting the voucher harvester CLI,
    dynamic genus-wide GBIF Packera species discovery, and dynamic taxon
    standardizer from scripts.data_prep.01_voucher_harvester and
    scripts.core.harvester.
===============================================================================
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig
from scripts.core.harvester import (
    VoucherHarvester,
    discover_genus_packera_species_gbif,
    setup_logger,
    standardize_packera_taxon,
)

_data_prep_harvester = importlib.import_module("scripts.data_prep.01_voucher_harvester")
build_cli_parser = _data_prep_harvester.build_cli_parser
main = _data_prep_harvester.main

if __name__ == "__main__":
    main()
