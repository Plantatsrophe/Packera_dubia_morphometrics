"""
===============================================================================
Module: lm2_config_builder.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Programmatically constructs, validates, merges, loads, and saves LeafMachine2
    configuration dictionaries and YAML configuration files.
===============================================================================
"""

from __future__ import annotations

import copy
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import yaml

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

logger = logging.getLogger("LM2_ConfigBuilder")


def get_default_lm2_template() -> Dict[str, Any]:
    """Returns a deep copy of the canonical LeafMachine2 configuration template."""
    return copy.deepcopy(DEFAULT_PACKERA_CONFIG_TEMPLATE)


def resolve_pcd_weights(
    weights_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Path] = None
) -> Tuple[str, Path]:
    """
    Resolves the weights name and path to the Plant Component Detector weights file.

    Returns:
        Tuple[str, Path]: (weights_name, resolved_path)
    """
    root = project_root or Path(__file__).resolve().parents[2]
    target_name = Path(weights_path).name if weights_path else DEFAULT_PCD_WEIGHTS_NAME

    if weights_path:
        p = Path(weights_path)
        if p.exists():
            return p.name, p.resolve()

    candidates = [
        root / "models" / target_name,
        root / "models" / "Packera_LeafPriority.pth",
        root / "models" / "LeafPriority.pth",
        root / "models" / "LeafPriority.pt",
        root / "LeafMachine2" / "leafmachine2" / "component_detector" / "models" / target_name,
    ]
    for c in candidates:
        if c.exists():
            return target_name, c.resolve()

    return target_name, (root / "models" / target_name).resolve()


def generate_high_performance_config(
    images_dir: Union[str, Path] = "LM2_Project/Data/images",
    output_dir: Union[str, Path] = "LM2_Project/Data/output",
    run_name: str = "Packera_dubia_LM2",
    pcd_weights: Optional[str] = None,
    pcd_confidence: float = DEFAULT_PCD_CONFIDENCE,
    nms_thresh: float = DEFAULT_NMS_THRESH,
    pointrend_subdivision_steps: int = DEFAULT_SUBDIVISION_STEPS,
    min_leaf_area: int = DEFAULT_MIN_LEAF_AREA,
    seg_confidence: float = DEFAULT_SEG_CONFIDENCE,
    batch_size: int = 50,
    num_workers: int = 8,
    device: str = "cuda",
) -> Dict[str, Any]:
    """
    Generates a high-performance LeafMachine2 configuration dictionary.
    """
    cfg = get_default_lm2_template()
    proj = cfg["leafmachine"]["project"]
    proj["dir_images_local"] = str(images_dir)
    proj["dir_output"] = str(output_dir)
    proj["run_name"] = run_name

    weights_name = pcd_weights or DEFAULT_PCD_WEIGHTS_NAME

    pcd = cfg["leafmachine"]["plant_component_detector"]
    pcd["detector_weights"] = weights_name
    pcd["minimum_confidence_threshold"] = pcd_confidence
    pcd["PCD_confidence"] = pcd_confidence
    pcd["iou_threshold"] = nms_thresh
    pcd["batch_size"] = batch_size
    pcd["num_workers"] = num_workers
    pcd["device"] = device

    seg = cfg["leafmachine"]["leaf_segmentation"]
    seg["segmentation_model"] = DEFAULT_SEG_MODEL_NAME
    seg["minimum_confidence_threshold"] = seg_confidence
    seg["NMS_thresh"] = nms_thresh
    seg["pointrend_subdivision_steps"] = pointrend_subdivision_steps
    seg["min_leaf_area"] = min_leaf_area
    seg["minimum_leaf_area_px"] = min_leaf_area
    seg["calculate_elliptic_fourier_descriptors"] = True
    seg["elliptic_fourier_descriptor_order"] = DEFAULT_EFD_ORDER
    seg["efd_order"] = DEFAULT_EFD_ORDER
    seg["device"] = device
    seg["num_workers_seg"] = num_workers

    ruler = cfg["leafmachine"]["ruler"]
    ruler["num_workers_ruler"] = num_workers

    return cfg


def validate_lm2_config(config_dict: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """
    Validates structural integrity and data types of an LM2 configuration dictionary.
    """
    if not isinstance(config_dict, dict) or "leafmachine" not in config_dict:
        return False, "Missing top-level 'leafmachine' key"

    lm = config_dict["leafmachine"]
    if "project" not in lm:
        return False, "Missing 'project' subsection"

    pcd = lm.get("plant_component_detector", {})
    if pcd.get("minimum_confidence_threshold", 0.0) <= 0.0 or pcd.get("minimum_confidence_threshold", 1.0) > 1.0:
        return False, "Invalid PCD confidence threshold (must be in (0, 1])"

    return True, None


def save_config_yaml(config_dict: Dict[str, Any], output_path: Union[str, Path]) -> None:
    """Serializes configuration dictionary to a YAML file."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
    logger.info(f"Saved LM2 configuration YAML to {p}")


def load_config_yaml(yaml_path: Union[str, Path]) -> Dict[str, Any]:
    """Loads a YAML configuration file into a Python dictionary."""
    p = Path(yaml_path)
    if not p.exists():
        raise FileNotFoundError(f"Configuration file not found: {p}")
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def update_main_lm2_yaml(
    base_yaml_path: Union[str, Path] = "LeafMachine2/LeafMachine2.yaml",
    project_yaml_path: Union[str, Path] = "LM2_Project/configs/lm2_packera_highperf.yaml",
    **override_kwargs
) -> None:
    """
    Constructs high-performance configuration and synchronizes both base LeafMachine2.yaml
    and the project configuration file.
    """
    cfg = generate_high_performance_config(**override_kwargs)
    valid, err = validate_lm2_config(cfg)
    if not valid:
        raise ValueError(f"Generated configuration is invalid: {err}")

    save_config_yaml(cfg, project_yaml_path)

    base_p = Path(base_yaml_path)
    if base_p.parent.exists():
        save_config_yaml(cfg, base_p)
        logger.info(f"Synchronized main base config at {base_p}")
