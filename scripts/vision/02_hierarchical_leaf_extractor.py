#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path
from scripts.core.logger import setup_logging
from scripts.core.leaf_extraction import run_pipeline

logger = setup_logging()

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Hierarchical Precision Leaf Extractor for botanical morphometrics."
    )
    parser.add_argument(
        "--input-dir", type=Path, default=Path("data/raw_vouchers"),
        help="Input directory containing raw voucher images."
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/masks"),
        help="Root output directory for extracted silhouettes and measurements."
    )
    parser.add_argument(
        "--qc-dir", type=Path, default=Path("outputs/extraction_qc"),
        help="Directory to save extraction QC verification figures."
    )
    parser.add_argument(
        "--weights", type=Path, default=Path("models/yolov8_leaf_best.pt"),
        help="Path to YOLOv8x-seg weights for artifact gatekeeping."
    )
    parser.add_argument(
        "--conf-threshold", type=float, default=0.25,
        help="Confidence threshold for YOLO artifact and basal region detection."
    )
    parser.add_argument(
        "--max-vouchers", type=int, default=None,
        help="Limit number of vouchers to process."
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel processing workers."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable debug logging."
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run_pipeline(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        qc_dir=args.qc_dir,
        weights_path=args.weights,
        conf_threshold=args.conf_threshold,
        max_vouchers=args.max_vouchers,
        num_workers=args.workers,
        verbose=args.verbose
    )
