#!/usr/bin/env python3
"""
scripts/vision/02_hierarchical_leaf_extractor.py
=================================================
Hierarchical Precision Leaf Extractor for botanical morphometrics.
"""

import sys
import argparse
from pathlib import Path

# Add project root directory (two levels up from scripts/vision/) to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.logger import setup_logging
from scripts.core.leaf_extraction import run_pipeline, DEFAULT_MODEL_PATH, DEFAULT_SAM2_CHECKPOINT, DEFAULT_SAM2_CONFIG
from scripts.core.config import DEFAULT_RAW_DIR

logger = setup_logging()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Hierarchical Precision Leaf Extractor & Rosette Disentanglement Pipeline."
    )
    parser.add_argument(
        "--raw-dir", "--input-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help="Path to directory containing raw voucher images (default: data/raw_vouchers/)"
    )
    parser.add_argument(
        "--model-path", "--weights",
        type=Path,
        default=DEFAULT_MODEL_PATH,
        help="Path to YOLOv8-seg weights (default: models/yolov8_leaf_best.pt)"
    )
    parser.add_argument(
        "--conf-threshold",
        type=float,
        default=0.25,
        help="Confidence threshold for organ instance detection (default: 0.25)"
    )
    parser.add_argument(
        "--use-sam2",
        action="store_true",
        default=True,
        help="Enable SAM 2 point prompting for rosette disentanglement (default: True)"
    )
    parser.add_argument(
        "--no-sam2",
        dest="use_sam2",
        action="store_false",
        help="Disable SAM 2 and use marker-controlled watershed segmentation"
    )
    parser.add_argument(
        "--sam2-checkpoint",
        type=Path,
        default=DEFAULT_SAM2_CHECKPOINT,
        help="Path to SAM 2 model weights (.pt)"
    )
    parser.add_argument(
        "--sam2-model-cfg",
        type=str,
        default=DEFAULT_SAM2_CONFIG,
        help="SAM 2 model configuration name/path"
    )
    parser.add_argument(
        "--limit", "--max-vouchers",
        type=int,
        default=None,
        help="Limit number of vouchers to process"
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Purge prior output masks, cropped patches, and QC logs before execution"
    )
    parser.add_argument(
        "--no-overlays",
        action="store_true",
        help="Disable generation of QC visual overlays"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose DEBUG logging"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run_pipeline(
        raw_dir=args.raw_dir,
        model_path=args.model_path,
        conf_threshold=args.conf_threshold,
        use_sam2=args.use_sam2,
        sam2_checkpoint=args.sam2_checkpoint,
        sam2_model_cfg=args.sam2_model_cfg,
        limit=args.limit,
        save_overlays=not args.no_overlays,
        clean=args.clean
    )


if __name__ == "__main__":
    main()
