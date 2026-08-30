#!/usr/bin/env python3
"""
Wrapper entrypoint for Artifact-Robust Botanical Dataset Builder.
Delegates to scripts.build_artifact_robust_dataset.
"""
import sys
from pathlib import Path

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.dataset_builder import (
    build_artifact_robust_dataset,
    parse_args,
    run_synthetic_test_suite,
    logger
)

if __name__ == "__main__":
    cli_args = parse_args()
    if cli_args.test:
        success = run_synthetic_test_suite()
        sys.exit(0 if success else 1)

    build_artifact_robust_dataset(
        raw_vouchers_dir=cli_args.raw_dir,
        annotations_dir=cli_args.annotations_dir,
        curated_csv_path=cli_args.curated_csv,
        output_dir=cli_args.output_dir,
        config_yaml_path=cli_args.config_yaml,
        qc_output_dir=cli_args.qc_dir,
        negative_ratio=cli_args.negative_ratio,
        train_ratio=cli_args.train_ratio,
        val_ratio=cli_args.val_ratio,
        test_ratio=cli_args.test_ratio,
        min_instance_area=cli_args.min_instance_area,
        num_qc_plots=cli_args.num_qc_plots,
        limit=cli_args.limit,
        seed=cli_args.seed
    )
