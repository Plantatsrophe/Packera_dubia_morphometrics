"""
===============================================================================
Module: pointrend_export.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Model weight checkpoint packaging, configuration YAML export, and metadata
    synchronization for drop-in deployment with LeafMachine2's Detector_LM2.
===============================================================================
"""

from __future__ import annotations

import json
import logging
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Union

import yaml

logger = logging.getLogger("PointRendExport")


def export_lm2_compatible_checkpoint(
    trained_weights_path: Union[str, Path],
    cfg: Any,
    output_model_dir: Union[str, Path] = "models",
    model_name: str = "Packera_LeafPriority"
) -> Path:
    """
    Exports final model weights and configuration into the models/ directory
    and LeafMachine2 segmentation models tree.
    """
    out_dir = Path(output_model_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dest_pth = out_dir / f"{model_name}.pth"
    src_pth = Path(trained_weights_path)
    if src_pth.exists():
        shutil.copy2(src_pth, dest_pth)
        logger.info(f"Exported production model checkpoint -> {dest_pth}")

    # Write accompanying configuration YAML
    cfg_yaml_path = out_dir / "cfg_output.yaml"
    if cfg is not None and hasattr(cfg, "dump"):
        with open(cfg_yaml_path, "w", encoding="utf-8") as f:
            f.write(cfg.dump())
        logger.info(f"Exported accompanying Detectron2 config -> {cfg_yaml_path}")

    # Write model metadata descriptor
    metadata = {
        "model_name": model_name,
        "classes": ["ideal_leaf", "partial_leaf"],
        "architecture": "PointRend (Mask R-CNN ResNet-50 FPN)",
        "date_trained": datetime.now().isoformat(),
        "checkpoint": str(dest_pth.name),
    }
    meta_path = out_dir / "metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    return dest_pth
