#!/usr/bin/env python3
"""
===============================================================================
Script: 07_triage_dashboard_synthesis.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Standardized CLI entry point for Step 07 (Python):
    Multi-Evidence Taxonomic Decision Matrix, Expert Triage Queue Generation,
    Publication 6-Panel Synthesis Plate Rendering, and Formal Taxonomic
    Revision Report Generation.

Usage:
    python scripts/analysis/07_triage_dashboard_synthesis.py \
        --vouchers data/tables/curated_vouchers.csv \
        --morphometrics data/tables/morphometrics_misidentification_flags.csv \
        --vision-audit data/tables/label_noise_audit.csv \
        --multimodal-flags data/tables/multimodal_conflict_flags.csv \
        --gmm-summary outputs/reports/gmm_bayes_factors_summary.csv \
        --niche-summary outputs/reports/multimodal_spatial_rf_summary.csv \
        --output-queue data/tables/triage_queue.csv \
        --output-plot outputs/figures/Figure_Integrative_Packera_dubia_Revision.pdf \
        --output-report outputs/reports/Packera_dubia_Taxonomic_Revision_Summary.md
===============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure repository root is in sys.path
_repo_root = str(Path(__file__).resolve().parents[2])
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from scripts.analysis.run_triage_dashboard_synthesis import main

if __name__ == "__main__":
    main()
