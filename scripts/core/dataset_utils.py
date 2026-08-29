
import os
import sys
import shutil
import cv2
import math
import numpy as np
import random
import yaml
import glob
import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from pathlib import Path

# Common imports for the project
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE, DEFAULT_RAW_DIR, DEFAULT_CURATED_CSV, DEFAULT_OUTPUT_DIR, DEFAULT_CONFIG_PATH, DEFAULT_QC_DIR
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

# Initialize shared structured pipeline logger
logger = setup_logging()


def stratify_and_partition_dataset(
    vouchers: List[Dict[str, Any]],
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    rng_seed: int = 42
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Deterministically partitions vouchers into Train, Validation, and Test splits,
    stratified by herbarium institution code and voucher quality tier.
    
    Args:
        vouchers: List of voucher metadata records.
        train_ratio: Proportion allocated to training (default 0.70).
        val_ratio: Proportion allocated to validation (default 0.15).
        test_ratio: Proportion allocated to testing (default 0.15).
        rng_seed: Random seed for deterministic reproducibility.
        
    Returns:
        Dictionary mapping split names ("train", "val", "test") to lists of records.
    """
    rng = random.Random(rng_seed)

    # Group vouchers by stratum: (institutionCode, determiner_tier)
    strata: Dict[str, List[Dict[str, Any]]] = {}
    for v in vouchers:
        inst = str(v.get("institutionCode", "UNKNOWN"))
        tier = str(v.get("determiner_tier", "Tier_Default"))
        stratum_key = f"{inst}_{tier}"
        if stratum_key not in strata:
            strata[stratum_key] = []
        strata[stratum_key].append(v)

    splits: Dict[str, List[Dict[str, Any]]] = {
        "train": [],
        "val": [],
        "test": []
    }

    # Allocate each stratum proportionally
    for stratum_key, group in strata.items():
        rng.shuffle(group)
        n = len(group)
        n_train = int(round(n * train_ratio))
        n_val = int(round(n * val_ratio))
        # Ensure at least 1 sample in train if group >= 1
        if n_train == 0 and n > 0:
            n_train = 1

        train_items = group[:n_train]
        val_items = group[n_train:n_train + n_val]
        test_items = group[n_train + n_val:]

        splits["train"].extend(train_items)
        splits["val"].extend(val_items)
        splits["test"].extend(test_items)

    logger.info(
        f"Stratified Dataset Partition: Train={len(splits['train'])} | "
        f"Val={len(splits['val'])} | Test={len(splits['test'])} "
        f"(Total: {len(vouchers)} vouchers across {len(strata)} strata)"
    )

    return splits


def render_qc_verification_overlay(
    image_bgr: np.ndarray,
    annotations: List[InstanceAnnotation],
    catalog_number: str,
    output_path: Path
) -> None:
    """
    Renders high-visibility multi-class bounding boxes, polygon contours,
    class labels, and synthetic augmentation flags for visual verification.
    """
    overlay = image_bgr.copy()
    h, w = overlay.shape[:2]

    # Create semi-transparent overlay layer for masks
    mask_layer = np.zeros_like(overlay, dtype=np.uint8)

    for ann in annotations:
        color = CLASS_COLORS_BGR.get(ann.class_id, (255, 255, 255))
        # Draw polygon mask if available
        if len(ann.polygon) >= 3:
            pts = ann.polygon.astype(np.int32)
            cv2.fillPoly(mask_layer, [pts], color)
            cv2.polylines(overlay, [pts], isClosed=True, color=color, thickness=2)

        # Draw bounding box
        x1, y1, x2, y2 = [int(v) for v in ann.bbox]
        cv2.rectangle(overlay, (x1, y1), (x2, y2), color, thickness=2)

        # Label tag
        tag_str = f"{ann.class_name}"
        if ann.is_synthetic:
            tag_str += " [SYN_AUG]"

        (tw, th), _ = cv2.getTextSize(tag_str, cv2.FONT_HERSHEY_SIMPLEX, 0.5, 1)
        cv2.rectangle(overlay, (x1, max(0, y1 - th - 6)), (x1 + tw + 6, max(0, y1)), color, -1)
        cv2.putText(
            overlay, tag_str, (x1 + 3, max(th, y1 - 3)),
            cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1, cv2.LINE_AA
        )

    # Blend polygon masks with 30% alpha
    cv2.addWeighted(mask_layer, 0.35, overlay, 0.65, 0, overlay)

    # Header banner with catalogNumber and instance tally
    cv2.rectangle(overlay, (0, 0), (w, 40), (20, 20, 20), -1)
    title = f"QC Verification: {catalog_number} | Instances: {len(annotations)}"
    cv2.putText(overlay, title, (15, 26), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (255, 255, 255), 2, cv2.LINE_AA)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(output_path), overlay)


def export_ultralytics_dataset_yaml(
    dataset_root: Path,
    output_yaml_path: Path
) -> None:
    """
    Exports the dataset YAML configuration file for training Ultralytics YOLOv8 / YOLOv11.
    """
    # Relative or absolute POSIX paths for compatibility across platforms
    dataset_dict = {
        "path": str(dataset_root.resolve()).replace("\\", "/"),
        "train": "images/train",
        "val": "images/val",
        "test": "images/test",
        "names": {idx: name for idx, name in enumerate(CLASS_NAMES)}
    }

    output_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_yaml_path, "w", encoding="utf-8") as f:
        yaml.dump(dataset_dict, f, default_flow_style=False, sort_keys=False)

    logger.info(f"Ultralytics dataset configuration written to: {output_yaml_path}")


