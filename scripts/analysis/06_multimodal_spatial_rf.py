#!/usr/bin/env python3
"""
===============================================================================
Script: 06_multimodal_spatial_rf.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Standardized CLI entry point for Step 06 (Python):
    Multimodal Spatial Macroecology, SoilGrids 250m Pedology, WorldClim v2.1
    Bioclimatics, Moran's Eigenvector Maps (MEMs), Spatial Random Forests,
    Cross-Modal Consensus Verification, and Warren's Niche Identity Tests.

Usage:
    python scripts/analysis/06_multimodal_spatial_rf.py \
        --vouchers data/tables/curated_vouchers.csv \
        --morphometrics data/tables/morphometrics_misidentification_flags.csv \
        --vision-audit data/tables/label_noise_audit.csv \
        --env-dir data/environmental/ \
        --output-flags data/tables/multimodal_conflict_flags.csv \
        --output-plot outputs/figures/spatial_rf_niche_importance.pdf \
        --output-summary outputs/reports/multimodal_spatial_rf_summary.csv
===============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure repository root is in sys.path
_repo_root = str(Path(__file__).resolve().parents[2])
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from scripts.analysis.run_multimodal_spatial_rf import main

if __name__ == "__main__":
    main()
