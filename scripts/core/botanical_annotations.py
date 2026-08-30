"""
Botanical Annotation Ingestion and Parsing Module
=================================================
Provides clean ingestion, parsing, and validation for 7-class botanical
phenotyping labels across human-annotated YOLO polygon .txt and manual JSON formats.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import cv2
import numpy as np

from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR
from scripts.core.data_structures import InstanceAnnotation
from scripts.core.logger import setup_logging

logger = setup_logging()

BOTANICAL_CLASSES: Dict[int, str] = {i: name for i, name in enumerate(CLASS_NAMES)}
BOTANICAL_CLASS_MAP: Dict[str, int] = CLASS_MAP
BotanicalAnnotation = InstanceAnnotation

# Botanical class aliases for flexible ingestion
BOTANICAL_ALIASES: Dict[str, str] = {
    "basal_leaf_blade": "basal_leaf_blade",
    "basal_leaf": "basal_leaf_blade",
    "leaf_blade": "basal_leaf_blade",
    "blade": "basal_leaf_blade",
    "basal_blade": "basal_leaf_blade",
    "leaf": "basal_leaf_blade",
    "leaf_petiole": "leaf_petiole",
    "petiole": "leaf_petiole",
    "leaf_stalk": "leaf_petiole",
    "stalk": "leaf_petiole",
    "cauline_leaf": "cauline_leaf",
    "cauline_leaves": "cauline_leaf",
    "stem_leaf": "cauline_leaf",
    "stem_leaves": "cauline_leaf",
    "cauline_stem": "cauline_stem",
    "stem": "cauline_stem",
    "flowering_stem": "cauline_stem",
    "peduncle": "cauline_stem",
    "scape": "cauline_stem",
    "root_rhizome": "root_rhizome",
    "root": "root_rhizome",
    "roots": "root_rhizome",
    "rhizome": "root_rhizome",
    "caudex": "root_rhizome",
    "basal_rosette_clump": "basal_rosette_clump",
    "basal_rosette": "basal_rosette_clump",
    "rosette": "basal_rosette_clump",
    "rosette_clump": "basal_rosette_clump",
    "crown": "basal_rosette_clump",
    "capitulum": "capitulum",
    "capitula": "capitulum",
    "flower_head": "capitulum",
    "inflorescence": "capitulum",
    "head": "capitulum",
    "involucre": "capitulum",
}


def normalize_botanical_class(raw_class: Union[int, str]) -> Optional[int]:
    """Maps an integer index or string class name/alias to canonical class_id (0..6)."""
    if isinstance(raw_class, int):
        if 0 <= raw_class < len(CLASS_NAMES):
            return raw_class
        return None

    clean = str(raw_class).strip().lower().replace("-", "_").replace(" ", "_")
    if clean.isdigit():
        idx = int(clean)
        if 0 <= idx < len(CLASS_NAMES):
            return idx
        return None

    canonical = BOTANICAL_ALIASES.get(clean)
    if canonical and canonical in CLASS_MAP:
        return CLASS_MAP[canonical]
    return None


def parse_botanical_yolo_txt(
    txt_path: Path,
    img_w: int,
    img_h: int,
    min_instance_area: float = 50.0
) -> List[InstanceAnnotation]:
    """
    Parses a verified YOLO segmentation polygon .txt file into InstanceAnnotations.
    """
    annotations: List[InstanceAnnotation] = []
    if not txt_path.exists() or txt_path.stat().st_size == 0:
        return annotations

    with open(txt_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for line_num, line in enumerate(lines, start=1):
        parts = line.strip().split()
        if not parts:
            continue

        raw_cls = parts[0]
        class_id = normalize_botanical_class(raw_cls)
        if class_id is None:
            logger.warning(f"Unknown class '{raw_cls}' in {txt_path.name}:{line_num}. Skipping.")
            continue

        coords = [float(v) for v in parts[1:]]
        if len(coords) < 4:
            continue

        if len(coords) == 4:
            # Bounding box format
            xc, yc, bw, bh = coords
            if max(xc, yc, bw, bh) <= 1.05:
                xc, yc, bw, bh = xc * img_w, yc * img_h, bw * img_w, bh * img_h
            x1, y1 = max(0.0, xc - bw / 2.0), max(0.0, yc - bh / 2.0)
            x2, y2 = min(float(img_w), xc + bw / 2.0), min(float(img_h), yc + bh / 2.0)
            poly = np.array([[x1, y1], [x2, y1], [x2, y2], [x1, y2]], dtype=np.float32)
        else:
            # Polygon format
            pts = np.array(coords, dtype=np.float32).reshape(-1, 2)
            if np.max(pts) <= 1.05:
                pts[:, 0] *= img_w
                pts[:, 1] *= img_h
            pts[:, 0] = np.clip(pts[:, 0], 0.0, float(img_w))
            pts[:, 1] = np.clip(pts[:, 1], 0.0, float(img_h))
            poly = pts

        if len(poly) < 3:
            continue

        area = float(cv2.contourArea(poly.astype(np.float32)))
        if area < min_instance_area:
            continue

        xmin, ymin = float(np.min(poly[:, 0])), float(np.min(poly[:, 1]))
        xmax, ymax = float(np.max(poly[:, 0])), float(np.max(poly[:, 1]))

        ann = InstanceAnnotation(
            class_id=class_id,
            polygon=poly,
            bbox=(xmin, ymin, xmax, ymax),
            confidence=1.0,
            tag="verified_human_annotation"
        )
        annotations.append(ann)

    return annotations


def parse_botanical_json(
    json_path: Path,
    img_w: int,
    img_h: int,
    min_instance_area: float = 50.0
) -> List[InstanceAnnotation]:
    """
    Parses a LabelMe, COCO, or generic JSON file into InstanceAnnotations.
    """
    annotations: List[InstanceAnnotation] = []
    if not json_path.exists():
        return annotations

    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        logger.error(f"Failed to read JSON annotation {json_path}: {e}")
        return annotations

    # LabelMe format
    if isinstance(data, dict) and "shapes" in data:
        for shape in data["shapes"]:
            label = shape.get("label", "")
            class_id = normalize_botanical_class(label)
            if class_id is None:
                continue

            raw_pts = shape.get("points", [])
            if len(raw_pts) < 3:
                continue

            poly = np.array(raw_pts, dtype=np.float32)
            if np.max(poly) <= 1.05:
                poly[:, 0] *= img_w
                poly[:, 1] *= img_h
            poly[:, 0] = np.clip(poly[:, 0], 0.0, float(img_w))
            poly[:, 1] = np.clip(poly[:, 1], 0.0, float(img_h))

            area = float(cv2.contourArea(poly.astype(np.float32)))
            if area < min_instance_area:
                continue

            xmin, ymin = float(np.min(poly[:, 0])), float(np.min(poly[:, 1]))
            xmax, ymax = float(np.max(poly[:, 0])), float(np.max(poly[:, 1]))

            annotations.append(InstanceAnnotation(
                class_id=class_id,
                polygon=poly,
                bbox=(xmin, ymin, xmax, ymax),
                confidence=1.0,
                tag="verified_labelme_annotation"
            ))

    # Generic annotations list
    elif isinstance(data, dict) and "annotations" in data:
        for item in data["annotations"]:
            label = item.get("class_name") or item.get("label") or item.get("class_id") or item.get("category_id")
            class_id = normalize_botanical_class(label)
            if class_id is None:
                continue

            raw_poly = item.get("polygon") or item.get("segmentation") or item.get("points")
            if not raw_poly:
                continue

            if isinstance(raw_poly, list) and len(raw_poly) > 0:
                if isinstance(raw_poly[0], list) and isinstance(raw_poly[0][0], (int, float)):
                    if len(raw_poly[0]) >= 6 and not isinstance(raw_poly[0][0], list):
                        poly = np.array([float(v) for v in raw_poly[0]], dtype=np.float32).reshape(-1, 2)
                    else:
                        poly = np.array(raw_poly, dtype=np.float32)
                else:
                    poly = np.array(raw_poly, dtype=np.float32).reshape(-1, 2)
            else:
                continue

            if len(poly) < 3:
                continue

            if np.max(poly) <= 1.05:
                poly[:, 0] *= img_w
                poly[:, 1] *= img_h
            poly[:, 0] = np.clip(poly[:, 0], 0.0, float(img_w))
            poly[:, 1] = np.clip(poly[:, 1], 0.0, float(img_h))

            area = float(cv2.contourArea(poly.astype(np.float32)))
            if area < min_instance_area:
                continue

            xmin, ymin = float(np.min(poly[:, 0])), float(np.min(poly[:, 1]))
            xmax, ymax = float(np.max(poly[:, 0])), float(np.max(poly[:, 1]))

            annotations.append(InstanceAnnotation(
                class_id=class_id,
                polygon=poly,
                bbox=(xmin, ymin, xmax, ymax),
                confidence=1.0,
                tag="verified_json_annotation"
            ))

    return annotations


def extract_botanical_annotations(
    image_bgr: np.ndarray,
    artifact_anns: Optional[List[InstanceAnnotation]] = None,
    annotations_dir: Optional[Path] = None,
    catalog_number: Optional[str] = None,
    min_instance_area: float = 50.0
) -> List[InstanceAnnotation]:
    """
    Ingests verified botanical annotations from disk if available for the given voucher.
    """
    h, w = image_bgr.shape[:2]
    annotations: List[InstanceAnnotation] = []

    if annotations_dir and catalog_number:
        txt_path = annotations_dir / f"{catalog_number}.txt"
        if txt_path.exists():
            annotations.extend(parse_botanical_yolo_txt(txt_path, img_w=w, img_h=h, min_instance_area=min_instance_area))

        json_path = annotations_dir / f"{catalog_number}.json"
        if json_path.exists():
            annotations.extend(parse_botanical_json(json_path, img_w=w, img_h=h, min_instance_area=min_instance_area))

    return annotations


class VerifiedBotanicalLabelParser:
    """
    Parser interface for reading YOLO segmentation .txt and JSON annotation files.
    """
    def __init__(self, min_instance_area: float = 50.0):
        self.min_instance_area = min_instance_area

    def parse_yolo_txt(self, txt_path: Path, img_w: int, img_h: int) -> List[InstanceAnnotation]:
        return parse_botanical_yolo_txt(txt_path, img_w=img_w, img_h=img_h, min_instance_area=self.min_instance_area)

    def parse_json_file(self, json_path: Path, img_w: int, img_h: int) -> List[InstanceAnnotation]:
        return parse_botanical_json(json_path, img_w=img_w, img_h=img_h, min_instance_area=self.min_instance_area)
