#!/usr/bin/env python3
"""
===============================================================================
Script: finetune_lm2_pcd.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Fine-tuning runner for LeafMachine2's Plant Component Detector (PCD) using
    Detectron2 with a PointRend backbone on ground-truth COCO annotations.
    Features automated CUDA single-GPU batch size fallback, backbone layer freezing,
    isolated artifact paths, and checkpoint export to `models/lm2_packera_pcd_finetuned.pth`.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import torch

# Ensure project root is in sys.path
_current = Path(__file__).resolve()
PROJECT_ROOT = _current.parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# LeafMachine2 and Detectron2 vendored paths
LM2_DIR = PROJECT_ROOT / "LeafMachine2" / "leafmachine2"
if LM2_DIR.exists() and str(LM2_DIR) not in sys.path:
    sys.path.insert(0, str(LM2_DIR))

DETECTRON2_PATH = LM2_DIR / "segmentation" / "detectron2"
if DETECTRON2_PATH.exists() and str(DETECTRON2_PATH) not in sys.path:
    sys.path.insert(0, str(DETECTRON2_PATH))

try:
    from detectron2.engine import DefaultTrainer
    from detectron2.evaluation import COCOEvaluator
except ImportError:
    DefaultTrainer = object
    COCOEvaluator = None


class PointRendPackeraTrainer(DefaultTrainer):
    """
    Custom Detectron2 trainer for PointRend fine-tuning on Packera leaf morphology.
    """

    @classmethod
    def build_evaluator(cls, cfg: Any, dataset_name: str, output_folder: Optional[str] = None):
        if output_folder is None:
            output_folder = os.path.join(cfg.OUTPUT_DIR, "validation_eval")
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        if COCOEvaluator is not None:
            return COCOEvaluator(dataset_name, output_dir=output_folder)
        return None

    def build_hooks(self):
        hooks = super().build_hooks() if hasattr(super(), "build_hooks") else []
        return hooks


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("FineTuneLM2PCD")

DEFAULT_THING_CLASSES = ["ideal_leaf", "partial_leaf"]


def freeze_backbone_stages(model: Any, freeze_stages: int = 2) -> None:
    """
    Freezes early stages of ResNet / FPN backbone to prevent catastrophic forgetting.
    """
    if hasattr(model, "backbone") and hasattr(model.backbone, "bottom_up"):
        bottom_up = model.backbone.bottom_up
        if hasattr(bottom_up, "stem"):
            for param in bottom_up.stem.parameters():
                param.requires_grad = False
            logger.info("Froze backbone stem layer.")

        for stage_idx in range(2, freeze_stages + 1):
            stage_name = f"res{stage_idx}"
            if hasattr(bottom_up, stage_name):
                stage = getattr(bottom_up, stage_name)
                for param in stage.parameters():
                    param.requires_grad = False
                logger.info(f"Froze backbone {stage_name} parameters.")


def build_pointrend_cfg(
    train_dataset_name: str,
    val_dataset_name: Optional[str] = None,
    output_dir: Union[str, Path] = "data/annotations/pcd_training_output",
    base_weights_path: Optional[str] = None,
    num_classes: int = 2,
    base_lr: float = 0.0001,
    max_iters: int = 1000,
    batch_size: int = 2,
    num_workers: int = 2,
    device: str = "cuda",
) -> Any:
    """
    Builds a Detectron2 CfgNode configured with the PointRend ResNet-50 FPN architecture.
    """
    try:
        from detectron2.config import get_cfg
        from detectron2.projects import point_rend
    except ImportError as e:
        logger.error(f"Detectron2 / PointRend import failure: {e}")
        return None

    cfg = get_cfg()
    point_rend.add_pointrend_config(cfg)

    pointrend_yaml = (
        DETECTRON2_PATH
        / "projects"
        / "PointRend"
        / "configs"
        / "InstanceSegmentation"
        / "pointrend_rcnn_R_50_FPN_3x_coco.yaml"
    )
    if pointrend_yaml.exists():
        cfg.merge_from_file(str(pointrend_yaml))
        logger.info(f"Loaded PointRend architecture configuration from {pointrend_yaml.name}")

    cfg.DATASETS.TRAIN = (train_dataset_name,)
    cfg.DATASETS.TEST = (val_dataset_name,) if val_dataset_name else ()
    cfg.DATALOADER.NUM_WORKERS = num_workers

    cfg.OUTPUT_DIR = str(output_dir)
    Path(cfg.OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    if base_weights_path and Path(base_weights_path).exists():
        cfg.MODEL.WEIGHTS = str(base_weights_path)
    else:
        cfg.MODEL.WEIGHTS = "detectron2://ImageNetPretrained/MSRA/R-50.pkl"

    cfg.MODEL.DEVICE = device
    cfg.MODEL.ROI_HEADS.NUM_CLASSES = num_classes
    if hasattr(cfg.MODEL, "POINT_HEAD"):
        cfg.MODEL.POINT_HEAD.NUM_CLASSES = num_classes

    # Solver
    cfg.SOLVER.IMS_PER_BATCH = batch_size
    cfg.SOLVER.BASE_LR = base_lr
    cfg.SOLVER.MAX_ITER = max_iters
    cfg.SOLVER.STEPS = (int(max_iters * 0.6), int(max_iters * 0.8))
    cfg.SOLVER.GAMMA = 0.5
    cfg.SOLVER.WARMUP_ITERS = min(100, max(10, max_iters // 10))
    cfg.SOLVER.CHECKPOINT_PERIOD = max(200, max_iters // 2)

    return cfg


def register_coco_dataset(
    dataset_name: str,
    coco_json_path: Path,
    images_dir: Path,
    thing_classes: Optional[List[str]] = None,
) -> bool:
    """Registers COCO dataset with Detectron2 dataset catalog."""
    try:
        from detectron2.data.datasets import register_coco_instances
        from detectron2.data import MetadataCatalog, DatasetCatalog

        if dataset_name in DatasetCatalog.list():
            logger.info(f"Dataset '{dataset_name}' already registered.")
            return True

        register_coco_instances(dataset_name, {}, str(coco_json_path), str(images_dir))
        classes = thing_classes or DEFAULT_THING_CLASSES
        MetadataCatalog.get(dataset_name).thing_classes = classes
        logger.info(f"Registered COCO dataset '{dataset_name}' with classes {classes} from {coco_json_path}")
        return True
    except ImportError as e:
        logger.error(f"Detectron2 not available: {e}")
        return False


def train_with_batch_fallback(
    cfg: Any,
    freeze_stages: int = 2,
    initial_batch_size: int = 2,
) -> bool:
    """
    Executes training loop with automated single-GPU CUDA batch size fallback
    upon encountering OutOfMemoryError.
    """
    batch_size = initial_batch_size

    while batch_size >= 1:
        cfg.SOLVER.IMS_PER_BATCH = batch_size
        logger.info(f"Attempting PointRend training with IMS_PER_BATCH = {batch_size}...")

        try:
            trainer = PointRendPackeraTrainer(cfg)
            freeze_backbone_stages(trainer.model, freeze_stages=freeze_stages)
            trainer.resume_or_load(resume=False)
            trainer.train()
            logger.info(f"Training completed successfully at batch size {batch_size}.")
            return True

        except (torch.cuda.OutOfMemoryError, RuntimeError) as err:
            err_msg = str(err).lower()
            if "out of memory" in err_msg or "cuda oom" in err_msg:
                logger.warning(f"CUDA Out of Memory caught at IMS_PER_BATCH = {batch_size}: {err}")
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

                if batch_size > 1:
                    batch_size -= 1
                    cfg.MODEL.ROI_HEADS.BATCH_SIZE_PER_IMAGE = max(32, cfg.MODEL.ROI_HEADS.BATCH_SIZE_PER_IMAGE // 2)
                    logger.info(f"Retrying training with reduced batch size {batch_size} and smaller ROI proposal batch.")
                    continue
                else:
                    logger.error("CUDA Out of Memory occurred even at minimum batch size (IMS_PER_BATCH = 1).")
                    raise
            else:
                raise

    return False


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for LM2 PCD fine-tuning runner."""
    parser = argparse.ArgumentParser(
        description="Fine-tune LeafMachine2 PointRend PCD on Packera ground-truth annotations."
    )
    parser.add_argument(
        "--coco-json",
        type=Path,
        default=Path("data/annotations/annotations_packera_train.json"),
        help="Path to training COCO dataset JSON",
    )
    parser.add_argument(
        "--val-coco-json",
        type=Path,
        default=None,
        help="Optional path to validation COCO dataset JSON",
    )
    parser.add_argument(
        "--images-dir",
        type=Path,
        default=Path("data/raw_vouchers"),
        help="Directory containing voucher sheets",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/annotations/pcd_training_output"),
        help="Directory for training loss logs and intermediate checkpoints",
    )
    parser.add_argument(
        "--output-weights",
        type=Path,
        default=Path("models/lm2_packera_pcd_finetuned.pth"),
        help="Destination path for final fine-tuned model checkpoint",
    )
    parser.add_argument(
        "--base-weights",
        type=str,
        default=None,
        help="Path to initial weights (defaults to models/Packera_LeafPriority.pth if present)",
    )
    parser.add_argument("--iterations", type=int, default=1000, help="Total training iterations")
    parser.add_argument("--lr", type=float, default=0.0001, help="Base learning rate")
    parser.add_argument("--batch-size", type=int, default=2, help="Initial images per batch")
    parser.add_argument("--freeze-stages", type=int, default=2, help="Number of backbone ResNet stages to freeze")
    parser.add_argument("--device", type=str, default=None, help="Device ('cuda' or 'cpu')")
    parser.add_argument("--dry-run", action="store_true", help="Validate configuration and data loading without training")
    return parser.parse_args()


def main() -> None:
    """Main execution flow for LM2 PCD fine-tuning."""
    args = parse_args()

    logger.info("=" * 80)
    logger.info("       LEAFMACHINE2 PCD / POINTREND FINE-TUNING PIPELINE       ")
    logger.info("=" * 80)

    if not args.coco_json.exists():
        logger.error(f"Training COCO annotation file not found: {args.coco_json}")
        sys.exit(1)

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Target execution device: {device}")

    base_weights = args.base_weights
    if not base_weights:
        default_model = Path("models/Packera_LeafPriority.pth")
        if default_model.exists():
            base_weights = str(default_model)
            logger.info(f"Using existing base checkpoint: {base_weights}")

    train_ds_name = "packera_pcd_finetune_train"
    val_ds_name = "packera_pcd_finetune_val" if args.val_coco_json else None

    registered = register_coco_dataset(
        train_ds_name, args.coco_json, args.images_dir, thing_classes=DEFAULT_THING_CLASSES
    )
    if not registered:
        logger.error("Failed to register COCO dataset with Detectron2.")
        sys.exit(1)

    if val_ds_name:
        register_coco_dataset(
            val_ds_name, args.val_coco_json, args.images_dir, thing_classes=DEFAULT_THING_CLASSES
        )

    cfg = build_pointrend_cfg(
        train_dataset_name=train_ds_name,
        val_dataset_name=val_ds_name,
        output_dir=args.output_dir,
        base_weights_path=base_weights,
        num_classes=len(DEFAULT_THING_CLASSES),
        base_lr=args.lr,
        max_iters=args.iterations,
        batch_size=args.batch_size,
        device=device,
    )

    if cfg is None:
        logger.error("Failed to construct PointRend configuration.")
        sys.exit(1)

    if args.dry_run:
        logger.info("[Dry Run] PointRend fine-tuning configuration validated successfully.")
        logger.info(f"[Dry Run] Max iterations: {args.iterations}, Device: {device}, Initial batch: {args.batch_size}")
        logger.info(f"[Dry Run] Target output weight destination: {args.output_weights}")
        return

    # Execute training with automated single-GPU OOM fallback
    success = train_with_batch_fallback(
        cfg=cfg,
        freeze_stages=args.freeze_stages,
        initial_batch_size=args.batch_size,
    )

    if success:
        final_model_src = Path(cfg.OUTPUT_DIR) / "model_final.pth"
        if final_model_src.exists():
            args.output_weights.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(final_model_src, args.output_weights)
            logger.info(f"Successfully exported fine-tuned checkpoint -> {args.output_weights}")
        else:
            logger.warning(f"Expected model_final.pth not found in {cfg.OUTPUT_DIR}")


if __name__ == "__main__":
    main()
