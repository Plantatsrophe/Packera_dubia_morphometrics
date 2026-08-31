#!/usr/bin/env python3
"""
===============================================================================
Script: train_lm2_pointrend.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    PyTorch and Detectron2 Fine-Tuning Pipeline for LeafMachine2's (LM2)
    Plant Component Detector (PCD) using PointRend Instance Segmentation.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

from scripts.vision.pointrend_config import build_pointrend_cfg, freeze_backbone_stages
from scripts.vision.pointrend_export import export_lm2_compatible_checkpoint
from scripts.vision.pointrend_trainer import PointRendPackeraTrainer

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("TrainPointRend")


def register_packera_coco_dataset(
    dataset_name: str,
    coco_json_path: Path,
    images_dir: Path
) -> None:
    """Registers COCO dataset with Detectron2 dataset catalog."""
    try:
        from detectron2.data.datasets import register_coco_instances
        from detectron2.data import MetadataCatalog
        register_coco_instances(dataset_name, {}, str(coco_json_path), str(images_dir))
        MetadataCatalog.get(dataset_name).thing_classes = ["ideal_leaf", "partial_leaf"]
        logger.info(f"Registered COCO dataset: '{dataset_name}' from {coco_json_path}")
    except ImportError:
        logger.warning("Detectron2 not installed; dataset registration skipped.")


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Fine-tune LeafMachine2 PointRend PCD on Packera specimens."
    )
    parser.add_argument("--coco-json", type=Path, default=Path("LM2_Project/Data/annotations/coco_pcd_packera.json"))
    parser.add_argument("--val-coco-json", type=Path, default=None)
    parser.add_argument("--images-dir", type=Path, default=Path("LM2_Project/Data/images"))
    parser.add_argument("--output-dir", type=Path, default=Path("LM2_Project/Data/output/pcd_training"))
    parser.add_argument("--weights", type=str, default=None)
    parser.add_argument("--iterations", type=int, default=2500)
    parser.add_argument("--lr", type=float, default=0.0001)
    parser.add_argument("--freeze-stages", type=int, default=2)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Orchestrates PointRend training flow."""
    args = parse_args()
    logger.info("=" * 80)
    logger.info("       LEAFMACHINE2 POINTREND PLANT COMPONENT DETECTOR TRAINING       ")
    logger.info("=" * 80)

    train_ds_name = "packera_pcd_train"
    val_ds_name = "packera_pcd_val" if args.val_coco_json else None

    register_packera_coco_dataset(train_ds_name, args.coco_json, args.images_dir)
    if val_ds_name:
        register_packera_coco_dataset(val_ds_name, args.val_coco_json, args.images_dir)

    cfg = build_pointrend_cfg(
        train_dataset_name=train_ds_name,
        val_dataset_name=val_ds_name,
        output_dir=args.output_dir,
        base_weights_path=args.weights,
        base_lr=args.lr,
        max_iters=args.iterations,
    )

    if args.dry_run:
        logger.info(f"Dry-run mode active. Config verified for {args.iterations} iterations.")
        return

    trainer = PointRendPackeraTrainer(cfg)
    freeze_backbone_stages(trainer.model, freeze_stages=args.freeze_stages)
    trainer.resume_or_load(resume=False)
    trainer.train()

    # Export final model artifacts
    final_weights = Path(cfg.OUTPUT_DIR) / "model_final.pth"
    if final_weights.exists():
        export_lm2_compatible_checkpoint(final_weights, cfg)


if __name__ == "__main__":
    main()
