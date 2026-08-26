#!/usr/bin/env python3
"""
===============================================================================
Script: train_artifact_robust_yolo.py (Root Launcher Wrapper)
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: Deep Learning Architect & Senior Computer Vision Engineer
Date: August 2026

Description:
    Root launcher wrapper for `scripts/train_artifact_robust_yolo.py`.
    Delegates directly to the robust YOLOv8 training and evaluation pipeline.
===============================================================================
"""

import sys
from pathlib import Path

# Add project root and scripts directory to python path
project_root = Path(__file__).resolve().parent
scripts_dir = project_root / "scripts"
if str(scripts_dir) not in sys.path:
    sys.path.insert(0, str(scripts_dir))

from train_artifact_robust_yolo import main

if __name__ == "__main__":
    main()
