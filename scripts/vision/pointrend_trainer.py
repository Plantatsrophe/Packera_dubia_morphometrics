"""
===============================================================================
Module: pointrend_trainer.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Custom Detectron2 DefaultTrainer implementation with evaluation hooks,
    loss logging, and WandB experiment tracking integration for PointRend.
===============================================================================
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger("PointRendTrainer")

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
