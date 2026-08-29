#!/usr/bin/env python3
"""
===============================================================================
Script: train_yolo.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Entrypoint script that fine-tunes a high-capacity YOLOv8 Instance Segmentation
    model (YOLOv8x-seg or YOLOv8m-seg) specifically configured to maximize botanical
    organ recall while rigorously suppressing false-positive detections on herbarium
    sheet mounting artifacts.
    
    Refactored to use the modular `train/` package.
===============================================================================
"""

import sys
import logging
from pathlib import Path
import os
import psutil

# Add project root directory to sys.path to ensure absolute module imports resolve
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.train.config import parse_arguments, resolve_default_paths
from scripts.train.dataset import prepare_tiled_dataset_split
from scripts.train.trainer import RobustYOLOTrainer

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("ArtifactRobustYOLO")

def main():
    """
    Main orchestration entry point.
    """
    args = parse_arguments()
    paths = resolve_default_paths()

    # Determine which dataset YAML configuration to use
    if args.data is not None:
        data_config_path = Path(args.data) if Path(args.data).is_absolute() else paths["root"] / args.data
    else:
        data_config_path = paths["data_config"]

    # Explicitly re-partition tiled dataset if flag is supplied
    if args.split_tiled_dataset:
        logger.info("Executing manual tiled dataset re-partitioning (--split-tiled-dataset requested)...")
        prepare_tiled_dataset_split(
            tiled_dir=paths["tiled_dir"],
            config_yaml_path=paths["tiled_config"],
            train_ratio=args.train_ratio,
            val_ratio=args.val_ratio,
            test_ratio=args.test_ratio,
            force_resplit=True
        )
        data_config_path = paths["tiled_config"]

    if args.dry_run:
        logger.info("Executing DRY RUN: Setting epochs=1, batch=2, imgsz=640 for rapid pipeline verification.")
        args.epochs = 1
        args.batch = 2
        args.imgsz = 640

    cache_val = False if args.cache == "none" else args.cache

    trainer = RobustYOLOTrainer(
        weights=args.weights,
        data_config=data_config_path,
        imgsz=args.imgsz,
        batch=args.batch,
        epochs=args.epochs,
        device=args.device,
        workers=args.workers,
        cache=cache_val,
        resume=args.resume,
        experiment_name=args.name
    )

    if not args.eval_only:
        try:
            # 1. Execute Training Loop
            trainer.train()

            # 2. Export Best Model Checkpoint to models/yolov8_leaf_best.pt
            best_checkpoint_path = trainer.export_best_checkpoint()

            # 3. Evaluate and Export Comprehensive Metrics
            trainer.evaluator.evaluate_and_export_metrics(checkpoint_path=best_checkpoint_path, split="val")
        except KeyboardInterrupt:
            logger.warning("KeyboardInterrupt received! Cleaning up child PyTorch DataLoader workers to prevent zombies...")
            parent = psutil.Process(os.getpid())
            children = parent.children(recursive=True)
            for child in children:
                try:
                    child.kill()
                except psutil.NoSuchProcess:
                    pass
            logger.info(f"Successfully killed {len(children)} background workers. Exiting safely.")
            sys.exit(1)
    else:
        # Evaluation Only Mode
        eval_checkpoint = args.checkpoint if args.checkpoint else trainer.paths["best_model_export"]
        logger.info(f"Running evaluation-only mode using checkpoint: {eval_checkpoint}")
        trainer.evaluator.evaluate_and_export_metrics(checkpoint_path=eval_checkpoint, split="val")

    logger.info("Pipeline execution finished successfully.")


if __name__ == "__main__":
    main()
