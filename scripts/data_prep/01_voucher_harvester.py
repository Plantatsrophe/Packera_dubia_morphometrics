#!/usr/bin/env python3
"""
===============================================================================
Script: 01_voucher_harvester.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    CLI entry point to execute the botanical voucher harvesting pipeline.
    Queries GBIF and iDigBio repositories, downloads high-resolution herbarium
    specimen imagery, extracts EXIF/DPI metadata, filters by spatial/morphological
    criteria, and constructs curated datasets.
===============================================================================
"""

import sys
from pathlib import Path

# Add project root directory (two levels up from scripts/data_prep) to sys.path
# This ensures absolute module imports like `from scripts.core...` resolve correctly
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import main harvesting workflow function from core module
from scripts.core.harvester import main

if __name__ == "__main__":
    main()
