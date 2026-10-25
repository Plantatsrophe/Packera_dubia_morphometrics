#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py (Pipeline Alias)
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Pipeline module entry point re-exporting the voucher harvester CLI,
    Eastern North American Packera taxon definitions, and Senecio basionym
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
    EASTERN_NA_GBIF_QUERY_TAXA,
    EASTERN_NA_PACKERA_TAXA,
    VoucherHarvester,
    setup_logger,
    standardize_packera_taxon,
)

_data_prep_harvester = importlib.import_module("scripts.data_prep.01_voucher_harvester")
build_cli_parser = _data_prep_harvester.build_cli_parser
main = _data_prep_harvester.main

if __name__ == "__main__":
    main()
