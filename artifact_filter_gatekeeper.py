#!/usr/bin/env python3
"""
===============================================================================
Module: artifact_filter_gatekeeper.py (Root Interface Proxy)
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU) /
             Google Antigravity Image Processing & Pipeline Optimization
Author: Image Processing & Pipeline Optimization Engineer
Date: August 2026

Description:
    Root-level import proxy exposing all classes, dataclasses, functions,
    and the synthetic test runner from scripts/artifact_filter_gatekeeper.py.
===============================================================================
"""

import sys
from pathlib import Path

# Ensure scripts/ directory is in Python module search path
scripts_dir = Path(__file__).resolve().parent / "scripts"
if str(scripts_dir) not in sys.path:
    sys.path.insert(0, str(scripts_dir))

from scripts.artifact_filter_gatekeeper import (
    ArtifactDetection,
    GeometricMetrics,
    SpectralMetrics,
    TextureMetrics,
    FilterResult,
    ArtifactFilterGatekeeper,
    generate_synthetic_leaf,
    generate_synthetic_herbarium_label,
    generate_synthetic_color_chart_swatch,
    generate_synthetic_scale_ruler,
    generate_synthetic_mounting_tape,
    generate_synthetic_clumped_rosette,
    run_synthetic_test_suite,
    TestArtifactFilterGatekeeper,
    main
)

if __name__ == "__main__":
    main()
