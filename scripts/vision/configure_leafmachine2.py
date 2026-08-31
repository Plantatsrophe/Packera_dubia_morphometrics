#!/usr/bin/env python3
"""
===============================================================================
Script: configure_leafmachine2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    CLI entry point and programmatic facade to construct, validate, and update
    LeafMachine2 configuration files for high-performance execution on Packera.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

from scripts.vision.lm2_config_builder import (
    generate_high_performance_config,
    get_default_lm2_template,
    load_config_yaml,
    resolve_pcd_weights,
    save_config_yaml,
    update_main_lm2_yaml,
    validate_lm2_config,
)
from scripts.vision.lm2_config_template import (
    DEFAULT_EFD_ORDER,
    DEFAULT_MIN_LEAF_AREA,
    DEFAULT_NMS_THRESH,
    DEFAULT_PACKERA_CONFIG_TEMPLATE,
    DEFAULT_PCD_CONFIDENCE,
    DEFAULT_PCD_WEIGHTS_NAME,
    DEFAULT_SEG_CONFIDENCE,
    DEFAULT_SEG_MODEL_NAME,
    DEFAULT_SUBDIVISION_STEPS,
)

__all__ = [
    "DEFAULT_PACKERA_CONFIG_TEMPLATE",
    "DEFAULT_PCD_WEIGHTS_NAME",
    "DEFAULT_SEG_MODEL_NAME",
    "DEFAULT_PCD_CONFIDENCE",
    "DEFAULT_NMS_THRESH",
    "DEFAULT_SUBDIVISION_STEPS",
    "DEFAULT_MIN_LEAF_AREA",
    "DEFAULT_SEG_CONFIDENCE",
    "DEFAULT_EFD_ORDER",
    "get_default_lm2_template",
    "resolve_pcd_weights",
    "generate_high_performance_config",
    "validate_lm2_config",
    "save_config_yaml",
    "load_config_yaml",
    "update_main_lm2_yaml",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("LM2_Configurator")


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Programmatically configure LeafMachine2 for Packera dubia morphometrics."
    )
    parser.add_argument("-o", "--output-yaml", type=Path, default=Path("LM2_Project/configs/lm2_packera_highperf.yaml"))
    parser.add_argument("--update-main-config", action="store_true", help="Synchronize LeafMachine2/LeafMachine2.yaml")
    parser.add_argument("--pcd-weights", type=str, default=None)
    parser.add_argument("--pcd-confidence", type=float, default=DEFAULT_PCD_CONFIDENCE)
    parser.add_argument("--nms-thresh", type=float, default=DEFAULT_NMS_THRESH)
    parser.add_argument("--pointrend-subdivision-steps", type=int, default=DEFAULT_SUBDIVISION_STEPS)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--num-workers", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Main execution flow."""
    args = parse_args()
    cfg = generate_high_performance_config(
        pcd_weights=args.pcd_weights,
        pcd_confidence=args.pcd_confidence,
        nms_thresh=args.nms_thresh,
        pointrend_subdivision_steps=args.pointrend_subdivision_steps,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )

    valid, err = validate_lm2_config(cfg)
    if not valid:
        logger.error(f"Configuration validation failed: {err}")
        sys.exit(1)

    if args.dry_run:
        import yaml
        print(yaml.dump(cfg, default_flow_style=False, sort_keys=False))
        logger.info("Dry-run inspection complete. No files written.")
        return

    if args.update_main_config:
        update_main_lm2_yaml(project_yaml_path=args.output_yaml)
    else:
        save_config_yaml(cfg, args.output_yaml)


if __name__ == "__main__":
    main()
