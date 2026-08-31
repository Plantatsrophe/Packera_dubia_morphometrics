"""
===============================================================================
Module: coco_exporter.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Converts binary instance segmentation masks into COCO 1.0 polygon annotations
    aligned with LeafMachine2's Plant Component Detector (PCD) taxonomy.
===============================================================================
"""

from __future__ import annotations

import json
import logging
import random
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import cv2
import numpy as np

logger = logging.getLogger("COCOExporter")

# Strict class mapping from SAM 2 botanical annotations to LeafMachine2 PCD taxonomy
PCD_CLASS_MAPPING: Dict[str, str] = {
    "basal_leaf_whole": "ideal_leaf",
    "basal_leaf_partial": "partial_leaf",
}

# Categories intentionally omitted so Detectron2 treats them as background noise
OMITTED_CLASSES: Set[str] = {
    "cauline_leaf",
    "cauline_stem",
    "root_rhizome",
    "basal_rosette_clump",
    "capitulum",
}

# Standardized COCO 1.0 categories for LeafMachine2 Plant Component Detector
PCD_COCO_CATEGORIES: List[Dict[str, Any]] = [
    {
        "id": 1,
        "name": "ideal_leaf",
        "supercategory": "plant_component",
    },
    {
        "id": 2,
        "name": "partial_leaf",
        "supercategory": "plant_component",
    },
]

CATEGORY_NAME_TO_ID: Dict[str, int] = {
    cat["name"]: cat["id"] for cat in PCD_COCO_CATEGORIES
}


def mask_to_polygons(
    binary_mask: np.ndarray,
    min_area_px: float = 50.0,
    approx_epsilon: float = 1.0,
) -> Tuple[List[List[float]], float, List[float]]:
    """
    Extracts vectorized polygon contours and bounding boxes from a 2D binary mask.

    Args:
        binary_mask: 2D uint8 binary array (0=background, >0=foreground).
        min_area_px: Minimum contour area threshold in pixels.
        approx_epsilon: Polygon approximation epsilon for contour point decimation.

    Returns:
        Tuple: (segmentation_polygons, total_area, [x_min, y_min, width, height])
    """
    if binary_mask is None or binary_mask.size == 0 or np.count_nonzero(binary_mask) == 0:
        return [], 0.0, [0.0, 0.0, 0.0, 0.0]

    contours, hierarchy = cv2.findContours(
        binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
    )

    polygons: List[List[float]] = []
    total_area = 0.0
    all_x: List[float] = []
    all_y: List[float] = []

    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area < min_area_px:
            continue

        if approx_epsilon > 0:
            cnt = cv2.approxPolyDP(cnt, approx_epsilon, True)

        if len(cnt) < 3:
            continue

        flattened = cnt.flatten().astype(float).tolist()
        if len(flattened) >= 6:
            polygons.append(flattened)
            total_area += area
            all_x.extend(flattened[0::2])
            all_y.extend(flattened[1::2])

    if not polygons or not all_x:
        return [], 0.0, [0.0, 0.0, 0.0, 0.0]

    x_min = min(all_x)
    y_min = min(all_y)
    x_max = max(all_x)
    y_max = max(all_y)
    bbox = [round(x_min, 2), round(y_min, 2), round(x_max - x_min, 2), round(y_max - y_min, 2)]

    return polygons, round(total_area, 2), bbox


def parse_mask_filename(mask_filename: str) -> Tuple[str, str, int]:
    """
    Parses catalogNumber, class label, and instance ID from a mask filename.

    Naming pattern: {catalogNumber}_{class_label}_instance_{id}.png
    Example: 'NCU00001234_basal_leaf_whole_instance_1.png'
    """
    stem = Path(mask_filename).stem
    m = re.match(r"^(.*?)(?:_instance_(\d+))?$", stem)
    if m:
        core_name = m.group(1)
        inst_id = int(m.group(2)) if m.group(2) else 1
    else:
        core_name = stem
        inst_id = 1

    matched_class = None
    for known_class in list(PCD_CLASS_MAPPING.keys()) + list(OMITTED_CLASSES):
        if f"_{known_class}" in core_name:
            matched_class = known_class
            cat_num = core_name.replace(f"_{known_class}", "")
            return cat_num, matched_class, inst_id

    return core_name, "unknown", inst_id


def convert_masks_to_coco_dataset(
    masks_dir: Path,
    image_dim_map: Optional[Dict[str, Tuple[int, int]]] = None,
    default_width: int = 4000,
    default_height: int = 6000,
    min_area_px: float = 50.0,
    image_extension: str = ".jpg",
) -> Dict[str, Any]:
    """
    Scans a directory of binary masks and compiles a COCO 1.0 JSON format dictionary.

    Args:
        masks_dir: Path to directory containing binary PNG masks.
        image_dim_map: Dict mapping image file base stems to (width, height).
        default_width: Fallback width if dimensions cannot be inferred.
        default_height: Fallback height if dimensions cannot be inferred.
        min_area_px: Minimum polygon area threshold.
        image_extension: Target herbarium sheet extension (default .jpg).

    Returns:
        Dict[str, Any]: Standardized COCO 1.0 format dataset.
    """
    masks_dir = Path(masks_dir)
    if not masks_dir.exists():
        logger.warning(f"Masks directory not found: {masks_dir}")
        return {"images": [], "annotations": [], "categories": PCD_COCO_CATEGORIES}

    mask_files = sorted([f for f in masks_dir.iterdir() if f.suffix.lower() == ".png"])
    logger.info(f"Found {len(mask_files)} mask files in {masks_dir}")

    images_dict: Dict[str, Dict[str, Any]] = {}
    annotations_list: List[Dict[str, Any]] = []

    image_id_counter = 1
    annotation_id_counter = 1

    for mask_path in mask_files:
        cat_num, class_label, inst_id = parse_mask_filename(mask_path.name)

        if class_label in OMITTED_CLASSES:
            continue

        pcd_category = PCD_CLASS_MAPPING.get(class_label)
        if not pcd_category or pcd_category not in CATEGORY_NAME_TO_ID:
            continue

        category_id = CATEGORY_NAME_TO_ID[pcd_category]

        img_file_name = f"{cat_num}{image_extension}"
        if img_file_name not in images_dict:
            w, h = default_width, default_height
            if image_dim_map and cat_num in image_dim_map:
                w, h = image_dim_map[cat_num]

            images_dict[img_file_name] = {
                "id": image_id_counter,
                "file_name": img_file_name,
                "width": w,
                "height": h,
            }
            curr_image_id = image_id_counter
            image_id_counter += 1
        else:
            curr_image_id = images_dict[img_file_name]["id"]

        mask_img = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
        if mask_img is None:
            continue

        if images_dict[img_file_name]["width"] == default_width:
            h, w = mask_img.shape[:2]
            images_dict[img_file_name]["width"] = w
            images_dict[img_file_name]["height"] = h

        polygons, area, bbox = mask_to_polygons(mask_img, min_area_px=min_area_px)
        if not polygons or area < min_area_px:
            continue

        annotations_list.append({
            "id": annotation_id_counter,
            "image_id": curr_image_id,
            "category_id": category_id,
            "segmentation": polygons,
            "area": area,
            "bbox": bbox,
            "iscrowd": 0,
        })
        annotation_id_counter += 1

    coco_doc = {
        "info": {
            "description": "LeafMachine2 Plant Component Detector - Packera Dataset",
            "version": "1.0",
            "year": 2026,
            "contributor": "NCU Herbarium / UNC Chapel Hill",
        },
        "licenses": [],
        "images": list(images_dict.values()),
        "annotations": annotations_list,
        "categories": PCD_COCO_CATEGORIES,
    }

    logger.info(
        f"Generated COCO dataset: {len(coco_doc['images'])} images, "
        f"{len(coco_doc['annotations'])} annotations across {len(PCD_COCO_CATEGORIES)} categories."
    )
    return coco_doc


def split_coco_dataset(
    coco_data: Dict[str, Any],
    val_ratio: float = 0.20,
    seed: int = 42
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Partitions a COCO dataset into training and validation splits at the image level.
    """
    random.seed(seed)
    images = coco_data.get("images", [])
    annotations = coco_data.get("annotations", [])
    categories = coco_data.get("categories", PCD_COCO_CATEGORIES)
    info = coco_data.get("info", {})

    image_ids = [img["id"] for img in images]
    random.shuffle(image_ids)

    n_val = int(len(image_ids) * val_ratio)
    val_image_ids = set(image_ids[:n_val])
    train_image_ids = set(image_ids[n_val:])

    train_images = [img for img in images if img["id"] in train_image_ids]
    val_images = [img for img in images if img["id"] in val_image_ids]

    train_annotations = [ann for ann in annotations if ann["image_id"] in train_image_ids]
    val_annotations = [ann for ann in annotations if ann["image_id"] in val_image_ids]

    train_doc = {
        "info": {**info, "split": "train"},
        "licenses": coco_data.get("licenses", []),
        "images": train_images,
        "annotations": train_annotations,
        "categories": categories,
    }

    val_doc = {
        "info": {**info, "split": "val"},
        "licenses": coco_data.get("licenses", []),
        "images": val_images,
        "annotations": val_annotations,
        "categories": categories,
    }

    return train_doc, val_doc


def save_coco_json(coco_data: Dict[str, Any], output_path: Path) -> None:
    """Writes COCO dataset dictionary to a JSON file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(coco_data, f, indent=2)
    logger.info(f"Saved COCO dataset to {output_path}")
