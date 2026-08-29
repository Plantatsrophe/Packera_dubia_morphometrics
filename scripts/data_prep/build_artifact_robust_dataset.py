#!/usr/bin/env python3
import sys
import argparse
import logging
from pathlib import Path

# Add project root directory (two levels up from scripts/data_prep) to sys.path
# This guarantees that absolute imports starting with `scripts.` resolve regardless of CWD or launch method
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.logger import setup_logging
from scripts.core.config import DEFAULT_RAW_DIR, DEFAULT_CURATED_CSV, DEFAULT_OUTPUT_DIR, DEFAULT_CONFIG_PATH, DEFAULT_QC_DIR
from scripts.core.dataset_builder import build_artifact_robust_dataset

logger = setup_logging()

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Construct an artifact-robust multi-class YOLO segmentation dataset for botanical herbarium specimens."
    )
    parser.add_argument(
        "--raw-images-dir", "--raw-dir", dest="raw_dir", type=Path, default=DEFAULT_RAW_DIR,
        help="Path to directory containing raw voucher JPG images."
    )
    parser.add_argument(
        "--annotations-dir", type=Path, default=None,
        help="Optional directory containing pre-computed raw annotations."
    )
    parser.add_argument(
        "--curated-csv", type=Path, default=DEFAULT_CURATED_CSV,
        help="Path to curated_vouchers.csv containing metadata for stratification."
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for the formatted YOLO dataset."
    )
    parser.add_argument(
        "--config-yaml", type=Path, default=DEFAULT_CONFIG_PATH,
        help="Path to export the Ultralytics dataset_config.yaml."
    )
    parser.add_argument(
        "--qc-dir", type=Path, default=DEFAULT_QC_DIR,
        help="Directory to save QC overlay visualization figures."
    )
    parser.add_argument(
        "--negative-ratio", type=float, default=0.09,
        help="Proportion of dataset comprising pure background negative sheets (default: 0.09 / ~9%)."
    )
    parser.add_argument(
        "--augment-prob", type=float, default=0.75,
        help="Probability of applying synthetic copy-paste augmentations to training samples."
    )
    parser.add_argument(
        "--copy-paste-prob", type=float, default=None,
        help="Alias for copy-paste augmentation probability."
    )
    parser.add_argument(
        "--train-ratio", type=float, default=0.70,
        help="Fraction of dataset allocated to training partition (default 0.70)."
    )
    parser.add_argument(
        "--val-ratio", type=float, default=0.15,
        help="Fraction of dataset allocated to validation partition (default 0.15)."
    )
    parser.add_argument(
        "--test-ratio", type=float, default=0.15,
        help="Fraction of dataset allocated to test partition (default 0.15)."
    )
    parser.add_argument(
        "--min-instance-area", type=float, default=150.0,
        help="Minimum bounding/contour area in square pixels for a botanical instance (default 150)."
    )
    parser.add_argument(
        "--num-qc-plots", type=int, default=25,
        help="Number of verification QC overlay images to generate in qc-dir."
    )
    parser.add_argument(
        "--num-workers", type=int, default=4,
        help="Number of processing threads/workers (default 4)."
    )
    parser.add_argument(
        "--limit", type=int, default=1500,
        help="Optional limit on the number of vouchers to process."
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for deterministic stratification and augmentations."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable detailed debug logging."
    )
    return parser.parse_args()

if __name__ == "__main__":
    cli_args = parse_args()
    if cli_args.verbose:
        logger.setLevel(logging.DEBUG)

    build_artifact_robust_dataset(
        raw_vouchers_dir=cli_args.raw_dir,
        annotations_dir=cli_args.annotations_dir,
        curated_csv_path=cli_args.curated_csv,
        output_dir=cli_args.output_dir,
        config_yaml_path=cli_args.config_yaml,
        qc_output_dir=cli_args.qc_dir,
        negative_ratio=cli_args.negative_ratio,
        augment_prob=cli_args.augment_prob,
        copy_paste_prob=cli_args.copy_paste_prob,
        train_ratio=cli_args.train_ratio,
        val_ratio=cli_args.val_ratio,
        test_ratio=cli_args.test_ratio,
        min_instance_area=cli_args.min_instance_area,
        num_qc_plots=cli_args.num_qc_plots,
        num_workers=cli_args.num_workers,
        limit=cli_args.limit,
        seed=cli_args.seed
    )
