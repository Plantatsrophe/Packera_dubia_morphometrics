#!/usr/bin/env python3
"""
===============================================================================
Script: 02_postprocess_lm2_routing.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Standardized CLI entry point for Step 02:
    LeafMachine2 output ingestion, DBSCAN multi-plant spatial clustering,
    geometric quality gatekeeping, 4-tier morphometric routing, and master
    QC logging.

Usage:
    python scripts/vision/02_postprocess_lm2_routing.py \
        --lm2-dir LM2_Project/Data/output/Packera_dubia_LM2/ \
        --vouchers data/tables/curated_vouchers.csv \
        --raw-images data/raw_vouchers/ \
        --output-dir data/ \
        --min-solidity 0.72 \
        --min-ucs 0.85
===============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure repository root is in sys.path
_repo_root = str(Path(__file__).resolve().parents[2])
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from scripts.vision.postprocess_lm2_routing import parse_args, run_postprocess_routing

if __name__ == "__main__":
    args = parse_args()
    run_postprocess_routing(
        lm2_dir=args.lm2_dir,
        vouchers_path=args.vouchers,
        raw_images_dir=args.raw_images,
        output_dir=args.output_dir,
        min_solidity=args.min_solidity,
        min_ucs=args.min_ucs,
        limit=args.limit
    )
