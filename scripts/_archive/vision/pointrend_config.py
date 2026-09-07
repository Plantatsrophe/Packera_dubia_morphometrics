"""
===============================================================================
Module: pointrend_config.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Detectron2 and PointRend configuration tree construction, backbone layer
    freezing, and hyperparameter configuration for Plant Component Detector training.
===============================================================================
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

logger = logging.getLogger("PointRendConfig")


def freeze_backbone_stages(model: Any, freeze_stages: int = 2) -> None:
    """
    Freezes early stages of ResNet / FPN backbone to avoid catastrophic forgetting.

    Args:
        model: Detectron2 GeneralizedRCNN model instance.
        freeze_stages: Stage index up to which parameters will be frozen (e.g. 2 for stem+res2).
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
    output_dir: Union[str, Path] = "LM2_Project/Data/output/pcd_training",
    base_weights_path: Optional[str] = None,
    num_classes: int = 3,
    base_lr: float = 0.0001,
    max_iters: int = 2500,
    batch_size_per_image: int = 128,
    num_workers: int = 4,
    device: str = "cuda",
    min_rpn_size: float = 32.0,
    anchor_sizes: Optional[List[List[int]]] = None,
) -> Any:
    """
    Constructs a Detectron2 CfgNode configured with PointRend head for Packera leaves.
    """
    try:
        from detectron2.config import get_cfg
        from detectron2.projects import point_rend
    except ImportError:
        logger.warning("Detectron2 not importable in active python environment.")
        return None

    cfg = get_cfg()
    point_rend.add_pointrend_config(cfg)

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
    cfg.MODEL.ROI_HEADS.BATCH_SIZE_PER_IMAGE = batch_size_per_image

    if hasattr(cfg.MODEL, "POINT_HEAD"):
        cfg.MODEL.POINT_HEAD.NUM_CLASSES = num_classes

    # Multi-level FPN anchor sizes: 5 levels (p2-p6) removing small sub-lobar anchors
    if anchor_sizes is not None:
        cfg.MODEL.ANCHOR_GENERATOR.SIZES = anchor_sizes
    else:
        cfg.MODEL.ANCHOR_GENERATOR.SIZES = [[64], [128], [256], [512], [1024]]

    # Minimum proposal size to discard tiny lobe fragments in RPN
    cfg.MODEL.RPN.MIN_SIZE = min_rpn_size

    # Solver / Optimizer hyperparameters
    cfg.SOLVER.IMS_PER_BATCH = 2
    cfg.SOLVER.BASE_LR = base_lr
    cfg.SOLVER.MAX_ITER = max_iters
    cfg.SOLVER.STEPS = (int(max_iters * 0.6), int(max_iters * 0.8))
    cfg.SOLVER.GAMMA = 0.5
    cfg.SOLVER.WARMUP_ITERS = 100
    cfg.SOLVER.CHECKPOINT_PERIOD = 500

    return cfg
