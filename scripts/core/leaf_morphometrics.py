
import os
import sys
import logging
import math
import numpy as np
import cv2
import json
import glob
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
from pathlib import Path
from scipy import ndimage
from skimage.morphology import skeletonize

# Common imports
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable


def detect_leaf_midrib_axis(mask: np.ndarray, contour: np.ndarray) -> Tuple[float, Tuple[int, int], Tuple[int, int], np.ndarray]:
    """
    Detects the primary longitudinal midrib axis from apex to petiole base.
    """
    if len(contour) < 5:
        return 0.0, (0, 0), (0, 0), np.zeros(4)

    line_params = cv2.fitLine(contour, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]

    angle_rad = math.atan2(vy, vx)
    angle_deg = math.degrees(angle_rad)

    pts = contour.reshape(-1, 2).astype(np.float32)
    diff = pts - np.array([x0, y0])
    projections = diff[:, 0] * vx + diff[:, 1] * vy

    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))

    pt_min = (int(pts[min_idx][0]), int(pts[min_idx][1]))
    pt_max = (int(pts[max_idx][0]), int(pts[max_idx][1]))

    M = cv2.moments(mask)
    if M["m00"] > 0:
        cx = M["m10"] / M["m00"]
        cy = M["m01"] / M["m00"]
        dist_min_centroid = math.hypot(pt_min[0] - cx, pt_min[1] - cy)
        dist_max_centroid = math.hypot(pt_max[0] - cx, pt_max[1] - cy)

        if dist_min_centroid > dist_max_centroid:
            base_pt = pt_min
            apex_pt = pt_max
        else:
            base_pt = pt_max
            apex_pt = pt_min
    else:
        apex_pt = pt_min
        base_pt = pt_max

    return angle_deg, apex_pt, base_pt, line_params


def align_mask_horizontally(mask: np.ndarray, apex_pt: Tuple[int, int], base_pt: Tuple[int, int]) -> Tuple[np.ndarray, float, Tuple[int, int]]:
    """
    Rotates and centers the binary leaf mask such that the midrib axis is exactly
    horizontal, with the Apex on the Left (x=0) and Petiole Base on the Right.
    """
    dx = base_pt[0] - apex_pt[0]
    dy = base_pt[1] - apex_pt[1]
    current_angle_deg = math.degrees(math.atan2(dy, dx))
    rotation_needed = -current_angle_deg

    h, w = mask.shape[:2]
    cx, cy = w / 2.0, h / 2.0

    rot_mat = cv2.getRotationMatrix2D((cx, cy), rotation_needed, 1.0)
    cos = np.abs(rot_mat[0, 0])
    sin = np.abs(rot_mat[0, 1])
    new_w = int((h * sin) + (w * cos))
    new_h = int((h * cos) + (w * sin))

    rot_mat[0, 2] += (new_w / 2.0) - cx
    rot_mat[1, 2] += (new_h / 2.0) - cy

    rotated_mask = cv2.warpAffine(
        mask, rot_mat, (new_w, new_h),
        flags=cv2.INTER_NEAREST,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0
    )

    contours, _ = cv2.findContours(rotated_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        x, y, bw, bh = cv2.boundingRect(max(contours, key=cv2.contourArea))
        pad = 5
        x0 = max(0, x - pad)
        y0 = max(0, y - pad)
        x1 = min(new_w, x + bw + pad)
        y1 = min(new_h, y + bh + pad)
        cropped_aligned = rotated_mask[y0:y1, x0:x1]
    else:
        cropped_aligned = rotated_mask

    return cropped_aligned, rotation_needed, (0, cropped_aligned.shape[0] // 2)


def reconstruct_bilateral_symmetry(aligned_mask: np.ndarray) -> Tuple[Optional[np.ndarray], bool, str]:
    """
    Tier 2 Algorithm: Evaluates upper vs lower half-margin completeness and
    performs bilateral symmetry reflection across the longitudinal midrib axis.
    """
    h, w = aligned_mask.shape[:2]
    if h < 20 or w < 20:
        return None, False, "none"

    col_has_leaf = (aligned_mask > 0).any(axis=0)
    if not col_has_leaf.any():
        return None, False, "none"

    y_indices, _ = np.where(aligned_mask > 0)
    y_midrib = int(np.median(y_indices))

    upper_half = np.zeros_like(aligned_mask)
    lower_half = np.zeros_like(aligned_mask)

    upper_half[:y_midrib, :] = aligned_mask[:y_midrib, :]
    lower_half[y_midrib:, :] = aligned_mask[y_midrib:, :]

    upper_profile = []
    lower_profile = []

    for x in range(w):
        y_up = np.where(upper_half[:, x] > 0)[0]
        y_low = np.where(lower_half[:, x] > 0)[0]

        upper_profile.append(np.min(y_up) if len(y_up) > 0 else y_midrib)
        lower_profile.append(np.max(y_low) if len(y_low) > 0 else y_midrib)

    upper_profile = np.array(upper_profile, dtype=np.float32)
    lower_profile = np.array(lower_profile, dtype=np.float32)

    upper_diff = np.abs(np.diff(upper_profile))
    lower_diff = np.abs(np.diff(lower_profile))

    upper_step_cuts = np.sum(upper_diff > 8)
    lower_step_cuts = np.sum(lower_diff > 8)

    upper_defect_score = upper_step_cuts * 10.0 + np.std(upper_diff)
    lower_defect_score = lower_step_cuts * 10.0 + np.std(lower_diff)

    threshold_defect = 25.0
    if upper_defect_score < threshold_defect and upper_defect_score <= lower_defect_score:
        selected_half = "upper"
    elif lower_defect_score < threshold_defect:
        selected_half = "lower"
    else:
        return None, False, "none"

    reflected_canvas = np.zeros_like(aligned_mask)

    if selected_half == "upper":
        reflected_canvas[:y_midrib, :] = upper_half[:y_midrib, :]
        for y in range(y_midrib):
            target_y = y_midrib + (y_midrib - y)
            if 0 <= target_y < h:
                reflected_canvas[target_y, :] = upper_half[y, :]
    else:
        reflected_canvas[y_midrib:, :] = lower_half[y_midrib:, :]
        for y in range(y_midrib, h):
            target_y = y_midrib - (y - y_midrib)
            if 0 <= target_y < h:
                reflected_canvas[target_y, :] = lower_half[y, :]

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    reflected_canvas = cv2.morphologyEx(reflected_canvas, cv2.MORPH_CLOSE, kernel)
    reflected_canvas = (reflected_canvas > 127).astype(np.uint8) * 255

    return reflected_canvas, True, selected_half


def extract_open_margin_curve_and_traits(
    aligned_mask: np.ndarray,
    catalog_number: str
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Tier 3 Algorithm: Extracts the longest continuous unoccluded margin curve
    and measures scalar morphometric caliper traits.
    """
    h, w = aligned_mask.shape[:2]
    contours, _ = cv2.findContours(aligned_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    if not contours:
        return pd.DataFrame(), {"petiole_length_px": 0.0, "max_width_px": 0.0, "apex_angle_deg": 0.0}

    main_contour = max(contours, key=cv2.contourArea).reshape(-1, 2)

    num_pts = len(main_contour)
    if num_pts < 30:
        resampled_pts = main_contour
    else:
        indices = np.linspace(0, num_pts - 1, 100).astype(int)
        resampled_pts = main_contour[indices]

    min_x, min_y = np.min(resampled_pts, axis=0)
    max_x, max_y = np.max(resampled_pts, axis=0)
    span_x = max(max_x - min_x, 1)
    span_y = max(max_y - min_y, 1)

    curve_df = pd.DataFrame({
        "catalogNumber": catalog_number,
        "point_index": np.arange(len(resampled_pts)),
        "x_px": resampled_pts[:, 0],
        "y_px": resampled_pts[:, 1],
        "x_norm": (resampled_pts[:, 0] - min_x) / span_x,
        "y_norm": (resampled_pts[:, 1] - min_y) / span_y,
    })

    col_widths = [np.sum(aligned_mask[:, x] > 0) for x in range(w)]
    max_width_px = float(np.max(col_widths)) if col_widths else 0.0

    petiole_cutoff_width = 0.25 * max_width_px
    petiole_cols = [x for x in range(int(w * 0.5), w) if col_widths[x] <= petiole_cutoff_width]
    petiole_length_px = float(len(petiole_cols))

    apex_x, apex_y = 0, h // 2
    probe_x = int(w * 0.10)
    probe_y_indices = np.where(aligned_mask[:, probe_x] > 0)[0] if probe_x < w else []

    if len(probe_y_indices) >= 2:
        y_top = np.min(probe_y_indices)
        y_bot = np.max(probe_y_indices)
        v1 = np.array([probe_x - apex_x, y_top - apex_y], dtype=np.float32)
        v2 = np.array([probe_x - apex_x, y_bot - apex_y], dtype=np.float32)

        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-6)
        apex_angle_deg = float(math.degrees(math.acos(np.clip(cos_angle, -1.0, 1.0))))
    else:
        apex_angle_deg = 45.0

    traits_dict = {
        "petiole_length_px": petiole_length_px,
        "max_width_px": max_width_px,
        "apex_angle_deg": apex_angle_deg,
    }

    return curve_df, traits_dict


