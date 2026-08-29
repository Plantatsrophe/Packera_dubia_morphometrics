#!/usr/bin/env python3
import sys
import argparse
from pathlib import Path

# Add project root directory (two levels up from scripts/vision) to sys.path
# This ensures absolute module imports from `scripts.*` resolve correctly
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.logger import setup_logging
from scripts.core.gatekeeper_engine import ArtifactFilterGatekeeper, run_synthetic_test_suite

logger = setup_logging()

def main():
    """Command-line execution interface for artifact filtering and batch validation."""
    parser = argparse.ArgumentParser(
        description="Production Artifact Filter Gatekeeper for Botanical Herbarium Morphometrics."
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Execute synthetic unit test suite verifying 100%% rejection of artifacts and retention of leaves."
    )
    parser.add_argument(
        "--archive-dir",
        type=str,
        default="data/cropped_patches/annotations",
        help="Target archive directory for routed text and label patches."
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=10,
        help="Padding boundary in pixels for Stage 1 hard blanking."
    )
    parser.add_argument(
        "--rect-threshold",
        type=float,
        default=0.86,
        help="Maximum rectangularity threshold for Stage 2a geometric filter."
    )
    parser.add_argument(
        "--solidity-threshold",
        type=float,
        default=0.72,
        help="Minimum solidity threshold for Stage 2c intact leaf filter."
    )

    args = parser.parse_args()

    if args.test or len(sys.argv) == 1:
        success = run_synthetic_test_suite()
        sys.exit(0 if success else 1)
    else:
        logger.info("Artifact Filter Gatekeeper initialized with archive dir: %s", args.archive_dir)

if __name__ == "__main__":
    main()
