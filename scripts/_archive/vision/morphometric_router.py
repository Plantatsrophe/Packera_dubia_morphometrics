#!/usr/bin/env python3
"""
scripts/vision/morphometric_router.py
=====================================
4-Tier Morphometric Routing implementations:
  - Tier 1: Pristine Closed Silhouettes for 12-harmonic closed EFA
  - Tier 2: Hemi-Blade Bilateral Symmetry Reflection
  - Tier 3: Open Margin Curves & Caliper Measurements (petiole length, blade width, apex angle)
  - Tier 4: Whole-Rosette Dense Patches for DINOv2 embeddings
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import List, Optional, Tuple

import cv2
import numpy as np
import pandas as pd

try:
    from scripts.vision.lm2_data_loader import LeafCandidate
except ImportError:
    from lm2_data_loader import LeafCandidate


def route_tier1_silhouette(
    mask: np.ndarray,
    candidate: LeafCandidate,
    output_dir: Path
) -> str:
    """Tier 1: Direct Pristine Silhouette export."""
    save_path = output_dir / "masks" / "tier1_intact" / f"{candidate.catalog_number}_p{candidate.plant_individual_id}_leaf{candidate.leaf_id}.png"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    bin_mask = (mask > 0).astype(np.uint8) * 255
    cv2.imwrite(str(save_path), bin_mask)
    return str(save_path)


def route_tier2_reflected(
    mask: np.ndarray,
    candidate: LeafCandidate,
    output_dir: Path,
    p_apex: Tuple[int, int],
    p_base: Tuple[int, int]
) -> Tuple[Optional[str], bool]:
    """
    Tier 2: Hemi-Blade Bilateral Symmetry Reflection.
    Aligns mask along midrib, checks upper vs lower half integrity,
    and reflects intact half to synthesize a 100% complete bilateral silhouette.
    """
    h, w = mask.shape[:2]
    if h < 15 or w < 15:
        return None, False

    dx = p_base[0] - p_apex[0]
    dy = p_base[1] - p_apex[1]
    angle = math.degrees(math.atan2(dy, dx))

    cx, cy = w / 2.0, h / 2.0
    rot_mat = cv2.getRotationMatrix2D((cx, cy), -angle, 1.0)
    aligned = cv2.warpAffine(mask, rot_mat, (w, h), flags=cv2.INTER_NEAREST)

    y_indices, _ = np.where(aligned > 0)
    if len(y_indices) < 30:
        return None, False

    y_midrib = int(np.median(y_indices))

    upper_half = np.zeros_like(aligned)
    lower_half = np.zeros_like(aligned)
    upper_half[:y_midrib, :] = aligned[:y_midrib, :]
    lower_half[y_midrib:, :] = aligned[y_midrib:, :]

    upper_profile, lower_profile = [], []
    for x in range(w):
        y_up = np.where(upper_half[:, x] > 0)[0]
        y_low = np.where(lower_half[:, x] > 0)[0]
        upper_profile.append(np.min(y_up) if len(y_up) > 0 else y_midrib)
        lower_profile.append(np.max(y_low) if len(y_low) > 0 else y_midrib)

    upper_diff = np.abs(np.diff(upper_profile))
    lower_diff = np.abs(np.diff(lower_profile))
    upper_defects = np.sum(upper_diff > 8) * 10.0 + float(np.std(upper_diff))
    lower_defects = np.sum(lower_diff > 8) * 10.0 + float(np.std(lower_diff))

    threshold_defect = 30.0
    if upper_defects < threshold_defect and upper_defects <= lower_defects:
        selected_half = "upper"
    elif lower_defects < threshold_defect:
        selected_half = "lower"
    else:
        return None, False

    reflected = np.zeros_like(aligned)
    if selected_half == "upper":
        reflected[:y_midrib, :] = upper_half[:y_midrib, :]
        for y in range(y_midrib):
            target_y = y_midrib + (y_midrib - y)
            if 0 <= target_y < h:
                reflected[target_y, :] = upper_half[y, :]
    else:
        reflected[y_midrib:, :] = lower_half[y_midrib:, :]
        for y in range(y_midrib, h):
            target_y = y_midrib - (y - y_midrib)
            if 0 <= target_y < h:
                reflected[target_y, :] = lower_half[y, :]

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    reflected = cv2.morphologyEx(reflected, cv2.MORPH_CLOSE, kernel)
    reflected = (reflected > 127).astype(np.uint8) * 255

    save_path = output_dir / "masks" / "tier2_reflected" / f"{candidate.catalog_number}_p{candidate.plant_individual_id}_leaf{candidate.leaf_id}_reflected.png"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(save_path), reflected)
    return str(save_path), True


def route_tier3_open_curve(
    mask: np.ndarray,
    candidate: LeafCandidate,
    output_dir: Path,
    scale_mm_per_px: float
) -> Tuple[str, float, float, float]:
    """
    Tier 3: Open Margin Curves & Caliper Measurements.
    Extracts continuous margin coordinate series and measures scalar traits.
    """
    h, w = mask.shape[:2]
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return "", 0.0, 0.0, 0.0

    main_contour = max(contours, key=cv2.contourArea).reshape(-1, 2)
    num_pts = len(main_contour)

    if num_pts < 30:
        resampled = main_contour
    else:
        indices = np.linspace(0, num_pts - 1, 100).astype(int)
        resampled = main_contour[indices]

    min_x, min_y = np.min(resampled, axis=0)
    max_x, max_y = np.max(resampled, axis=0)
    span_x = max(max_x - min_x, 1)
    span_y = max(max_y - min_y, 1)

    curve_df = pd.DataFrame({
        "catalogNumber": candidate.catalog_number,
        "plant_individual_id": candidate.plant_individual_id,
        "leaf_id": candidate.leaf_id,
        "point_index": np.arange(len(resampled)),
        "x_px": resampled[:, 0],
        "y_px": resampled[:, 1],
        "x_norm": (resampled[:, 0] - min_x) / span_x,
        "y_norm": (resampled[:, 1] - min_y) / span_y,
    })

    save_path = output_dir / "masks" / "tier3_open_curves" / f"{candidate.catalog_number}_p{candidate.plant_individual_id}_leaf{candidate.leaf_id}_curve.csv"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    curve_df.to_csv(save_path, index=False)

    col_widths = [np.sum(mask[:, x] > 0) for x in range(w)]
    max_width_px = float(np.max(col_widths)) if col_widths else 0.0
    width_mm = max_width_px * scale_mm_per_px

    petiole_cutoff_px = 0.25 * max_width_px
    petiole_cols = [x for x in range(int(w * 0.5), w) if col_widths[x] <= petiole_cutoff_px]
    petiole_length_mm = float(len(petiole_cols)) * scale_mm_per_px

    apex_angle_deg = 45.0
    probe_x = int(w * 0.10)
    if probe_x < w:
        probe_y = np.where(mask[:, probe_x] > 0)[0]
        if len(probe_y) >= 2:
            y_top, y_bot = np.min(probe_y), np.max(probe_y)
            apex_x, apex_y = 0, h // 2
            v1 = np.array([probe_x - apex_x, y_top - apex_y], dtype=np.float32)
            v2 = np.array([probe_x - apex_x, y_bot - apex_y], dtype=np.float32)
            cos_ang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-6)
            apex_angle_deg = float(math.degrees(math.acos(np.clip(cos_ang, -1.0, 1.0))))

    return str(save_path), petiole_length_mm, width_mm, apex_angle_deg


def crop_dense_rosette_tier4(
    raw_image_path: Path,
    plant_candidates: List[LeafCandidate],
    plant_id: int,
    output_dir: Path
) -> Optional[str]:
    """Tier 4: Dense whole-rosette bounding crop from native-DPI sheet."""
    if not raw_image_path.exists():
        return None

    raw_img = cv2.imread(str(raw_image_path))
    if raw_img is None:
        return None

    img_h, img_w = raw_img.shape[:2]

    all_ymin, all_xmin, all_ymax, all_xmax = [], [], [], []
    for cand in plant_candidates:
        ymin, xmin, ymax, xmax = cand.bbox
        if ymax > ymin and xmax > xmin:
            all_ymin.append(ymin)
            all_xmin.append(xmin)
            all_ymax.append(ymax)
            all_xmax.append(xmax)

    if not all_ymin:
        ymin, ymax = int(img_h * 0.40), int(img_h * 0.85)
        xmin, xmax = int(img_w * 0.20), int(img_w * 0.80)
    else:
        min_y, max_y = min(all_ymin), max(all_ymax)
        min_x, max_x = min(all_xmin), max(all_xmax)
        pad_y = int(0.15 * (max_y - min_y + 10))
        pad_x = int(0.15 * (max_x - min_x + 10))
        ymin = max(0, min_y - pad_y)
        ymax = min(img_h, max_y + pad_y)
        xmin = max(0, min_x - pad_x)
        xmax = min(img_w, max_x + pad_x)

    rosette_crop = raw_img[ymin:ymax, xmin:xmax]
    if rosette_crop.size == 0:
        return None

    cat = plant_candidates[0].catalog_number
    save_path = output_dir / "cropped_patches" / "rosettes_dense" / f"{cat}_p{plant_id}_rosette.jpg"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(save_path), rosette_crop)
    return str(save_path)
