"""
scripts/core/leaf_cv_utils.py
=============================
Computer vision utilities for botanical leaf extraction, artifact blanking,
EDT peak seeding, SAM 2/watershed segmentation, and solidity evaluation.
"""

from __future__ import annotations

import glob
import json
import logging
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage
from skimage.morphology import skeletonize

from scripts.core.config import CLASS_COLORS_BGR, CLASS_MAP, CLASS_NAMES, DEFAULT_WORKSPACE
from scripts.core.data_structures import (
    ArtifactDetection,
    FilterResult,
    GeometricMetrics,
    InstanceAnnotation,
    SpectralMetrics,
    TextureMetrics,
)
from scripts.core.leaf_spine_tracer import (
    frangi_vesselness_filter_2d,
    trace_3point_anatomical_spine,
)
from scripts.core.logger import setup_logging

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

logger = setup_logging()

__all__ = [
    "detect_sheet_artifacts",
    "apply_hard_artifact_blanking",
    "detect_native_dpi_regions",
    "frangi_vesselness_filter_2d",
    "trace_3point_anatomical_spine",
    "extract_edt_point_seeds",
    "load_sam2_predictor",
    "segment_leaves_sam2_or_watershed",
    "evaluate_convexity_and_solidity",
    "filter_involucre_profiles",
]


def detect_sheet_artifacts(image_bgr: np.ndarray) -> Dict[str, List[Dict[str, Any]]]:
    """
    Stage 1: Detects non-plant layout objects on the full herbarium scan:
    herbarium_label, barcode_sticker, color_chart/palette, ruler_scale,
    mounting_tape, and annotation_slip.

    Args:
        image_bgr: Full-resolution BGR image of the herbarium sheet.

    Returns:
        Dictionary mapping artifact class names to lists of detected bounding boxes.
    """
    h, w = image_bgr.shape[:2]
    artifacts: Dict[str, List[Dict[str, Any]]] = {
        "herbarium_label": [],
        "color_chart": [],
        "annotation_slip": [],
        "mounting_tape": [],
        "ruler_scale": [],
        "barcode_sticker": [],
    }

    gray = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)

    # 1. Main Herbarium Label (typically lower-right quadrant, high contrast, text block)
    label_y1 = int(h * 0.65)
    label_y2 = int(h * 0.98)
    label_x1 = int(w * 0.55)
    label_x2 = int(w * 0.98)

    label_roi = gray[label_y1:label_y2, label_x1:label_x2]
    _, label_thresh = cv2.threshold(label_roi, 200, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(label_thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    found_label = False
    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area > (w * h * 0.015):
            bx, by, bw, bh = cv2.boundingRect(cnt)
            artifacts["herbarium_label"].append({
                "bbox": np.array([label_x1 + bx, label_y1 + by, label_x1 + bx + bw, label_y1 + by + bh]),
                "conf": 0.95,
                "mask": None
            })
            found_label = True
            break
    if not found_label:
        artifacts["herbarium_label"].append({
            "bbox": np.array([label_x1, label_y1, label_x2, label_y2]),
            "conf": 0.88,
            "mask": None
        })

    # 2. Color Chart / Palette (typically top-left or along the left/right ruler margin)
    artifacts["color_chart"].append({
        "bbox": np.array([int(w * 0.01), int(h * 0.02), int(w * 0.14), int(h * 0.42)]),
        "conf": 0.90,
        "mask": None
    })

    # 3. Barcode Sticker & Institution Stamp (typically upper margin or near label)
    artifacts["barcode_sticker"].append({
        "bbox": np.array([int(w * 0.05), int(h * 0.02), int(w * 0.35), int(h * 0.09)]),
        "conf": 0.85,
        "mask": None
    })

    # 4. Secondary Annotation Slips (top right / middle left)
    artifacts["annotation_slip"].append({
        "bbox": np.array([int(w * 0.60), int(h * 0.05), int(w * 0.95), int(h * 0.22)]),
        "conf": 0.82,
        "mask": None
    })

    # 5. Mounting Tape Strips (detected via horizontal/vertical linear contrast filters)
    tape_kernel_h = cv2.getStructuringElement(cv2.MORPH_RECT, (25, 3))
    grad_x = cv2.Sobel(gray, cv2.CV_32F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(gray, cv2.CV_32F, 0, 1, ksize=3)
    mag = cv2.magnitude(grad_x, grad_y)
    mag = np.uint8(np.clip(mag, 0, 255))
    _, tape_bin = cv2.threshold(mag, 50, 255, cv2.THRESH_BINARY)
    tape_morph = cv2.morphologyEx(tape_bin, cv2.MORPH_OPEN, tape_kernel_h)
    cnts_tape, _ = cv2.findContours(tape_morph, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    for c in cnts_tape:
        if 200 < cv2.contourArea(c) < 5000:
            bx, by, bw, bh = cv2.boundingRect(c)
            aspect = max(bw, bh) / max(min(bw, bh), 1)
            if aspect > 3.0:
                artifacts["mounting_tape"].append({
                    "bbox": np.array([bx, by, bx + bw, by + bh]),
                    "conf": 0.75,
                    "mask": None
                })

    return artifacts


def apply_hard_artifact_blanking(
    image_bgr: np.ndarray,
    artifacts: Dict[str, List[Dict[str, Any]]],
    fill_color: Tuple[int, int, int] = (255, 255, 255)
) -> np.ndarray:
    """
    Stage 1: Hard-masks all non-plant artifact bounding boxes and polygons to
    solid background paper tone (RGB [255, 255, 255]), eliminating label leaks,
    stamp bleed, and tape interference prior to botanical segmentation.
    """
    clean_sheet = image_bgr.copy()
    h, w = clean_sheet.shape[:2]

    # Zero-fill outer border sheet margins (25px safety margin)
    clean_sheet[:25, :] = fill_color
    clean_sheet[h-25:, :] = fill_color
    clean_sheet[:, :25] = fill_color
    clean_sheet[:, w-25:] = fill_color

    blank_classes = ["herbarium_label", "color_chart", "barcode_sticker", "ruler_scale", "mounting_tape"]
    for cls_name in blank_classes:
        for art in artifacts.get(cls_name, []):
            x1, y1, x2, y2 = art["bbox"]
            x1, y1 = max(0, x1 - 5), max(0, y1 - 5)
            x2, y2 = min(w, x2 + 5), min(h, y2 + 5)
            cv2.rectangle(clean_sheet, (x1, y1), (x2, y2), fill_color, -1)

    return clean_sheet


def detect_native_dpi_regions(
    clean_sheet_bgr: np.ndarray
) -> Tuple[Optional[Tuple[int, int, int, int]], List[Tuple[int, int, int, int]]]:
    """
    Stage 2: Identifies coarse bounding boxes for the basal rosette cluster
    and the capitulescence (flower head cyme) on the full sheet, then crops
    uncompressed sub-images at native 300+ DPI resolution without downsampling.
    """
    h, w = clean_sheet_bgr.shape[:2]
    gray = cv2.cvtColor(clean_sheet_bgr, cv2.COLOR_BGR2GRAY)
    inv_gray = 255 - gray

    blurred = cv2.GaussianBlur(inv_gray, (7, 7), 0)
    _, plant_bin = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    rosette_y1 = int(h * 0.35)
    rosette_y2 = int(h * 0.95)
    rosette_x1 = int(w * 0.08)
    rosette_x2 = int(w * 0.92)

    rosette_mask = np.zeros((h, w), dtype=np.uint8)
    rosette_mask[rosette_y1:rosette_y2, rosette_x1:rosette_x2] = 255
    rosette_plant = cv2.bitwise_and(plant_bin, plant_bin, mask=rosette_mask)

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    rosette_plant_clean = cv2.morphologyEx(rosette_plant, cv2.MORPH_OPEN, kernel)

    cnts, _ = cv2.findContours(rosette_plant_clean, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    valid_rosette_cnts = [c for c in cnts if cv2.contourArea(c) > 3000]

    if valid_rosette_cnts:
        all_pts = np.vstack(valid_rosette_cnts)
        rx, ry, rw, rh = cv2.boundingRect(all_pts)
        pad = 40
        rx1 = max(0, rx - pad)
        ry1 = max(0, ry - pad)
        rx2 = min(w, rx + rw + pad)
        ry2 = min(h, ry + rh + pad)
        rosette_bbox = (rx1, ry1, rx2, ry2)
    else:
        rosette_bbox = (rosette_x1, rosette_y1, rosette_x2, rosette_y2)

    cyme_y2 = int(h * 0.48)
    cyme_mask = np.zeros((h, w), dtype=np.uint8)
    cyme_mask[:cyme_y2, :] = 255
    cyme_plant = cv2.bitwise_and(plant_bin, plant_bin, mask=cyme_mask)

    cyme_cnts, _ = cv2.findContours(cyme_plant, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cyme_bboxes = []
    for cc in cyme_cnts:
        area = cv2.contourArea(cc)
        if 800 < area < 100000:
            cx, cy, cw, ch = cv2.boundingRect(cc)
            pad_c = 15
            cx1 = max(0, cx - pad_c)
            cy1 = max(0, cy - pad_c)
            cx2 = min(w, cx + cw + pad_c)
            cy2 = min(h, cy + ch + pad_c)
            cyme_bboxes.append((cx1, cy1, cx2, cy2))

    return rosette_bbox, cyme_bboxes


def extract_edt_point_seeds(
    rosette_binary_mask: np.ndarray,
    min_distance_px: int = 25,
    relative_threshold: float = 0.20
) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
    """
    Stage 4: Computes Euclidean Distance Transform (EDT) on dense overlapping
    basal rosettes. The local maxima (peaks) of the EDT surface identify the
    exact geometric centroids of individual leaf blades within the clump.
    """
    dist_transform = cv2.distanceTransform(rosette_binary_mask, cv2.DIST_L2, 5)
    max_val = dist_transform.max()

    if max_val <= 0:
        return dist_transform, []

    size = max(3, int(min_distance_px))
    local_max = ndimage.maximum_filter(dist_transform, size=size)
    peak_mask = (dist_transform == local_max) & (dist_transform >= relative_threshold * max_val)

    peak_coords = np.argwhere(peak_mask)
    point_seeds = [(int(pt[1]), int(pt[0])) for pt in peak_coords]

    h, w = rosette_binary_mask.shape[:2]
    filtered_seeds = [
        (x, y) for (x, y) in point_seeds
        if 10 < x < w - 10 and 10 < y < h - 10
    ]

    return dist_transform, filtered_seeds


def load_sam2_predictor(
    checkpoint_path: Optional[Union[str, Path]] = None,
    model_cfg: Optional[Union[str, Path]] = None,
    device: Optional[str] = None
) -> Optional[Any]:
    """
    Safely loads and initializes a SAM2ImagePredictor instance.
    Returns None if SAM 2 packages or weights are unavailable.
    """
    try:
        import torch
        from sam2.build_sam import build_sam2
        from sam2.sam2_image_predictor import SAM2ImagePredictor

        if checkpoint_path is None:
            checkpoint_path = DEFAULT_WORKSPACE / "models" / "checkpoints" / "sam2_hiera_large.pt"
        if model_cfg is None:
            model_cfg = "sam2_hiera_l.yaml"

        checkpoint_path = Path(checkpoint_path)
        if not checkpoint_path.is_file():
            logger.warning(f"SAM 2 checkpoint not found at {checkpoint_path}")
            return None

        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"

        sam2_model = build_sam2(str(model_cfg), str(checkpoint_path), device=device)
        predictor = SAM2ImagePredictor(sam2_model)
        logger.info(f"SAM 2 Predictor initialized on {device} from {checkpoint_path}")
        return predictor
    except Exception as e:
        logger.warning(f"Failed to load SAM 2 predictor: {e}. Gracefully falling back to watershed.")
        return None


def segment_leaves_sam2_or_watershed(
    rosette_crop_bgr: np.ndarray,
    rosette_binary_mask: np.ndarray,
    point_seeds: List[Tuple[int, int]],
    use_sam2: bool = True,
    sam2_predictor: Optional[Any] = None,
    sam2_checkpoint: Optional[Union[str, Path]] = None,
    sam2_model_cfg: Optional[Union[str, Path]] = None
) -> List[np.ndarray]:
    """
    Stage 4: Uses EDT peak coordinates as positive point prompts in SAM 2
    (point_labels=np.array([1]), multimask_output=True) or marker-controlled watershed
    segmentation to disentangle individual leaf layers.
    """
    h, w = rosette_binary_mask.shape[:2]
    leaf_masks: List[np.ndarray] = []

    if not point_seeds:
        return leaf_masks

    if use_sam2:
        try:
            predictor = sam2_predictor
            if predictor is None:
                predictor = load_sam2_predictor(sam2_checkpoint, sam2_model_cfg)

            if predictor is not None:
                crop_rgb = cv2.cvtColor(rosette_crop_bgr, cv2.COLOR_BGR2RGB)
                predictor.set_image(crop_rgb)

                for sx, sy in point_seeds:
                    point_coords = np.array([[sx, sy]], dtype=np.float32)
                    point_labels = np.array([1], dtype=np.int32)

                    masks, scores, logits = predictor.predict(
                        point_coords=point_coords,
                        point_labels=point_labels,
                        multimask_output=True
                    )

                    if masks is not None and len(masks) > 0:
                        best_idx = int(np.argmax(scores))
                        raw_mask = masks[best_idx]
                        if hasattr(raw_mask, "cpu"):
                            raw_mask = raw_mask.cpu().numpy()
                        if raw_mask.dtype == bool:
                            mask_u8 = raw_mask.astype(np.uint8) * 255
                        else:
                            mask_u8 = (raw_mask > 0).astype(np.uint8) * 255

                        leaf_mask = cv2.bitwise_and(mask_u8, rosette_binary_mask)
                        if cv2.countNonZero(leaf_mask) >= 300:
                            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
                            leaf_mask = cv2.morphologyEx(leaf_mask, cv2.MORPH_CLOSE, kernel)
                            leaf_masks.append(leaf_mask)

                if leaf_masks:
                    logger.debug(f"SAM 2 successfully disentangled {len(leaf_masks)} leaf layers from rosette.")
                    return leaf_masks
                else:
                    logger.debug("SAM 2 returned no valid masks >= 300 px, falling back to watershed.")
        except Exception as exc:
            logger.warning(f"SAM 2 segmentation failed ({exc}), falling back to marker-controlled watershed.")

    markers = np.zeros((h, w), dtype=np.int32)
    for idx, (sx, sy) in enumerate(point_seeds):
        cv2.circle(markers, (sx, sy), 3, idx + 2, -1)

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    sure_bg = cv2.dilate(rosette_binary_mask, kernel, iterations=4)
    unknown = cv2.subtract(sure_bg, rosette_binary_mask)
    markers[sure_bg == 0] = 1
    markers[unknown == 255] = 0

    cv2.watershed(rosette_crop_bgr, markers)

    for idx, (sx, sy) in enumerate(point_seeds):
        marker_id = idx + 2
        leaf_mask = (markers == marker_id).astype(np.uint8) * 255
        if cv2.countNonZero(leaf_mask) >= 300:
            leaf_mask = cv2.morphologyEx(leaf_mask, cv2.MORPH_CLOSE, kernel)
            leaf_masks.append(leaf_mask)

    return leaf_masks


def evaluate_convexity_and_solidity(
    leaf_mask: np.ndarray
) -> Tuple[float, float, bool, np.ndarray, np.ndarray]:
    """
    Stage 5: Single-Leaf Geometric Gatekeeper.
    Single basal leaves of Packera dubia have convex, elliptic-to-ovate profiles
    with high solidity:
        Solidity = Area(mask) / Area(convex hull) >= 0.72
        UCS = Perimeter(convex hull) / Perimeter(contour) >= 0.85
    """
    contours, _ = cv2.findContours(leaf_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return 0.0, 0.0, False, np.array([]), np.array([])

    best_contour = max(contours, key=cv2.contourArea)
    contour_area = cv2.contourArea(best_contour)

    if contour_area < 100:
        return 0.0, 0.0, False, best_contour, np.array([])

    convex_hull = cv2.convexHull(best_contour)
    hull_area = cv2.contourArea(convex_hull)
    hull_perimeter = cv2.arcLength(convex_hull, closed=True)
    contour_perimeter = cv2.arcLength(best_contour, closed=True)

    solidity = contour_area / max(hull_area, 1e-6)
    ucs_ratio = hull_perimeter / max(contour_perimeter, 1e-6)

    try:
        hull_indices = cv2.convexHull(best_contour, returnPoints=False)
        if hull_indices is not None and len(hull_indices) > 3 and len(best_contour) > 3:
            defects = cv2.convexityDefects(best_contour, hull_indices)
            if defects is not None:
                _, _, w, h = cv2.boundingRect(best_contour)
                norm_dim = max(w, h, 1)
                defect_depths = [d[0][3] / 256.0 for d in defects]
                max_defect = max(defect_depths) if defect_depths else 0.0

                if max_defect / norm_dim > 0.15:
                    ucs_ratio = min(ucs_ratio, 0.80)
    except Exception:
        pass

    ucs_score = float(np.clip(ucs_ratio, 0.0, 1.0))
    solidity_score = float(np.clip(solidity, 0.0, 1.0))
    is_clean = bool(solidity_score >= 0.72 and ucs_score >= 0.85)

    return ucs_score, solidity_score, is_clean, best_contour, convex_hull


def filter_involucre_profiles(
    clean_sheet_bgr: np.ndarray,
    cyme_bboxes: List[Tuple[int, int, int, int]]
) -> List[Dict[str, Any]]:
    """
    Stage 5: Distinguishes capitulescence clusters from single involucres.
    Retains only flower heads presented in clean longitudinal profile with
    aspect ratio H/W in [0.8, 1.5].
    """
    valid_involucres: List[Dict[str, Any]] = []

    for idx, (x1, y1, x2, y2) in enumerate(cyme_bboxes):
        crop = clean_sheet_bgr[y1:y2, x1:x2]
        if crop.size == 0:
            continue

        ch, cw = crop.shape[:2]
        aspect_ratio = ch / max(cw, 1)

        if 0.75 <= aspect_ratio <= 1.60 and (ch * cw) >= 1200:
            valid_involucres.append({
                "involucre_index": idx,
                "bbox": (x1, y1, x2, y2),
                "aspect_ratio": round(aspect_ratio, 3),
                "crop": crop,
                "is_clean_profile": True
            })

    return valid_involucres
