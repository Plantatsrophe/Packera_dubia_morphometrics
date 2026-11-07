#!/usr/bin/env python3
"""
===============================================================================
Module: lm2_geometry_utils.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Consolidated geometric quality gatekeeping, midrib pose estimation,
    bilateral hemi-blade reflection, adaptive DBSCAN spatial clustering,
    and standardized contour coordinate extraction.

Key Capabilities:
    1. Geometric Quality Gatekeeping:
       - Unoccluded Completeness Score (UCS >= 0.85): Area_mask / (0.70 * Area_min_rect)
       - Solidity (Solidity >= 0.72): Area_mask / Area_convex_hull
       - Longitudinal midrib axis angle and apex/base pose estimation
    2. Vectorized Hemi-Blade Reflection (Tier 2):
       - Primary midrib alignment via affine rotation
       - Vectorized margin profile defect calculation across upper/lower halves
       - Fully vectorized slice reflection in NumPy (zero manual pixel loops)
       - Morphological seam closure
    3. Adaptive DBSCAN Spatial Clustering:
       - Groups detected organs into plant_individual_id clusters to prevent
         trait averaging across multi-plant specimen sheets.
    4. Standardized Contour & Caliper Extraction:
       - Resamples contour boundary to standardized normalized (x, y) coordinates
       - Decoupled asynchronous Hough-transform ruler detection
===============================================================================
"""

from __future__ import annotations

import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN


# =============================================================================
# 1. Geometric Quality Gatekeeping & Pose Estimation
# =============================================================================

def _sample_cross_sectional_widths(
    mask: np.ndarray,
    pt_start: Tuple[int, int],
    pt_end: Tuple[int, int],
    sample_dist: float,
    step_px: int = 2,
    max_span: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized helper to sample perpendicular cross-sectional mask widths along
    the ray from pt_start toward pt_end over distance sample_dist.
    """
    h, w = mask.shape[:2]
    v_dir = np.array([pt_end[0] - pt_start[0], pt_end[1] - pt_start[1]], dtype=np.float64)
    L = float(np.hypot(v_dir[0], v_dir[1]))
    if L < 1e-5:
        return np.array([]), np.array([])
    u = v_dir / L
    n = np.array([-u[1], u[0]], dtype=np.float64)
    s_steps = np.arange(0, sample_dist, max(1, step_px), dtype=np.float64)
    if len(s_steps) == 0:
        return np.array([]), np.array([])
    centers = np.array(pt_start, dtype=np.float64)[None, :] + s_steps[:, None] * u[None, :]
    w_span = max(min(max(h, w), 350), int(0.65 * max(h, w)))
    if max_span is not None and max_span > 0:
        w_span = min(w_span, max_span)
    t_steps = np.arange(-w_span, w_span + 1, 1.0, dtype=np.float64)
    gx = np.round(centers[:, 0, None] + t_steps[None, :] * n[0]).astype(np.int32)
    gy = np.round(centers[:, 1, None] + t_steps[None, :] * n[1]).astype(np.int32)
    valid = (gx >= 0) & (gx < w) & (gy >= 0) & (gy < h)
    grid_mask = np.zeros((len(s_steps), len(t_steps)), dtype=bool)
    grid_mask[valid] = mask[gy[valid], gx[valid]] > 0
    widths = np.zeros(len(s_steps), dtype=np.float64)
    for i in range(len(s_steps)):
        row = grid_mask[i]
        nz = np.where(row)[0]
        widths[i] = float(nz[-1] - nz[0] + 1) if len(nz) > 0 else 0.0
    return s_steps, widths


def _disambiguate_base_apex_endpoints(
    mask: np.ndarray,
    ep1: Tuple[int, int],
    ep2: Tuple[int, int],
    rosette_centroid: Optional[Tuple[float, float]] = None,
    max_petiole_width: float = 22.0,
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """
    Disambiguates which longitudinal axis endpoint is p_base (petiole attachment)
    vs. p_apex (lamina tip):
      1. If rosette_centroid is provided, the endpoint closer to the rosette cluster
         centroid is p_base and the opposite end is p_apex.
      2. Otherwise, measures cross-sectional mask widths across the terminal 20%
         of each end: the narrower, steeper-tapering stalk end is p_base and the
         opposite (broader lamina) end is p_apex.

    Returns
    -------
    Tuple[Tuple[int, int], Tuple[int, int]]
        (p_apex, p_base)
    """
    if rosette_centroid is not None:
        d1 = math.hypot(ep1[0] - rosette_centroid[0], ep1[1] - rosette_centroid[1])
        d2 = math.hypot(ep2[0] - rosette_centroid[0], ep2[1] - rosette_centroid[1])
        if abs(d1 - d2) > 1e-3:
            if d1 < d2:
                return ep2, ep1
            return ep1, ep2

    axis_len = math.hypot(ep2[0] - ep1[0], ep2[1] - ep1[1])
    if axis_len < 10.0:
        return ep1, ep2

    eff_max_petiole_w = (
        max(float(max_petiole_width), 0.085 * axis_len)
        if (max_petiole_width == 22.0 and axis_len > 450.0)
        else float(max_petiole_width)
    )

    probe_dist = max(8.0, 0.20 * axis_len)
    s1, w1 = _sample_cross_sectional_widths(mask, ep1, ep2, probe_dist, step_px=2)
    s2, w2 = _sample_cross_sectional_widths(mask, ep2, ep1, probe_dist, step_px=2)

    def _terminal_profile_stats(s_arr: np.ndarray, w_arr: np.ndarray) -> Tuple[float, float, float]:
        if len(w_arr) == 0:
            return 0.0, float("inf"), 0.0
        # Contiguous narrow stalk length from the tip inward
        stalk_len = 0.0
        for idx, w_val in enumerate(w_arr):
            if w_val <= eff_max_petiole_w:
                stalk_len = float(s_arr[idx])
            elif idx > 1:
                break
        pos_w = w_arr[w_arr > 0]
        mean_w = float(np.mean(pos_w)) if len(pos_w) > 0 else float("inf")
        dw = np.diff(w_arr) / 2.0 if len(w_arr) >= 2 else np.array([0.0])
        max_flare = float(np.max(dw)) if len(dw) > 0 else 0.0
        return stalk_len, mean_w, max_flare

    stalk1, mean_w1, flare1 = _terminal_profile_stats(s1, w1)
    stalk2, mean_w2, flare2 = _terminal_profile_stats(s2, w2)

    # 1. Explicit narrow petiole stalk on one end
    if stalk1 >= 6.0 and stalk1 >= stalk2 + 4.0:
        return ep2, ep1  # ep1 is p_base
    if stalk2 >= 6.0 and stalk2 >= stalk1 + 4.0:
        return ep1, ep2  # ep2 is p_base

    # 2. Narrower terminal 20% cross-sectional width (cuneate/attenuate base vs broader apical lobe)
    if abs(mean_w1 - mean_w2) > 0.75:
        if mean_w1 < mean_w2:
            return ep2, ep1  # ep1 is narrower -> p_base
        return ep1, ep2      # ep2 is narrower -> p_base

    # 3. Steeper tapering flare gradient at the terminal 20%
    if abs(flare1 - flare2) > 0.10:
        if flare1 > flare2:
            return ep2, ep1
        return ep1, ep2

    # Deterministic tie-breaker for perfectly symmetric synthetic shapes:
    # lower point in image (larger y, or larger x if horizontal) is p_base
    if (ep1[1], ep1[0]) > (ep2[1], ep2[0]):
        return ep2, ep1
    return ep1, ep2


def compute_geometric_metrics(
    mask: np.ndarray,
    rosette_centroid: Optional[Tuple[float, float]] = None,
) -> Tuple[float, float, float, Tuple[int, int], Tuple[int, int]]:
    """
    Computes geometric quality and midrib pose metrics for a binary leaf mask:
      - ucs: Unoccluded Completeness Score = Area_mask / (0.70 * Area_min_rect)
      - solidity: Area_mask / Area_convex_hull
      - midrib_angle_deg: Angle of longitudinal midrib axis from apex to base
      - p_apex: Estimated apex coordinate (x, y)
      - p_base: Estimated base/petiole coordinate (x, y)

    Parameters
    ----------
    mask : np.ndarray
        Binary uint8 mask of the detected leaf instance.
    rosette_centroid : Optional[Tuple[float, float]]
        Optional (x, y) centroid of the rosette cluster to disambiguate p_base.

    Returns
    -------
    Tuple[float, float, float, Tuple[int, int], Tuple[int, int]]
        (ucs, solidity, midrib_angle_deg, p_apex, p_base)
    """
    if mask is None or mask.size == 0:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    h_full, w_full = mask_u8.shape[:2]
    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if bw == 0 or bh == 0:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    x0_pad = max(0, bx - 2)
    y0_pad = max(0, by - 2)
    x1_pad = min(w_full, bx + bw + 2)
    y1_pad = min(h_full, by + bh + 2)
    sub_mask = mask_u8[y0_pad:y1_pad, x0_pad:x1_pad]

    area_mask = float(np.count_nonzero(sub_mask))
    if area_mask < 30:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    contours, _ = cv2.findContours(sub_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    main_contour = max(contours, key=cv2.contourArea)

    # 1. Convex Hull & Solidity
    hull = cv2.convexHull(main_contour)
    area_hull = float(cv2.contourArea(hull))
    solidity = (area_mask / area_hull) if area_hull > 0 else 0.0
    solidity = min(max(solidity, 0.0), 1.0)

    # 2. Expected Area via Minimum Area Rect (Oriented Bounding Box)
    rect = cv2.minAreaRect(main_contour)
    rect_w, rect_h = rect[1]
    area_rect = max(rect_w * rect_h, 1.0)
    expected_area = 0.70 * area_rect
    ucs = min(max(area_mask / expected_area, 0.0), 1.0)

    # 3. Longitudinal Midrib Axis (moments principal axis with fitLine fallback)
    m = cv2.moments(sub_mask, binaryImage=True)
    if m["m00"] > 0 and (abs(m["mu20"] - m["mu02"]) > 1e-5 or abs(m["mu11"]) > 1e-5):
        x0 = float(m["m10"] / m["m00"])
        y0 = float(m["m01"] / m["m00"])
        theta_m = 0.5 * math.atan2(2.0 * m["mu11"], m["mu20"] - m["mu02"])
        vx, vy = math.cos(theta_m), math.sin(theta_m)
    else:
        line_params = cv2.fitLine(main_contour, cv2.DIST_L2, 0, 0.01, 0.01)
        vx, vy, x0, y0 = [float(v[0]) for v in line_params]

    # Project contour points along longitudinal midrib axis and select on-axis extreme tips
    pts = main_contour.reshape(-1, 2).astype(np.float64)
    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy
    perp_dist = np.abs(-(pts[:, 0] - x0) * vy + (pts[:, 1] - y0) * vx)

    p_min = float(np.min(projections))
    p_max = float(np.max(projections))
    near_min = np.where(projections <= p_min + 1.5)[0]
    near_max = np.where(projections >= p_max - 1.5)[0]

    min_idx = int(near_min[np.argmin(perp_dist[near_min])]) if len(near_min) > 0 else int(np.argmin(projections))
    max_idx = int(near_max[np.argmin(perp_dist[near_max])]) if len(near_max) > 0 else int(np.argmax(projections))

    ep1_sub = (int(round(pts[min_idx, 0])), int(round(pts[min_idx, 1])))
    ep2_sub = (int(round(pts[max_idx, 0])), int(round(pts[max_idx, 1])))

    rc_sub = (
        (float(rosette_centroid[0]) - x0_pad, float(rosette_centroid[1]) - y0_pad)
        if rosette_centroid is not None
        else None
    )

    # Disambiguate p_base vs p_apex polarity via terminal 20% width profile / rosette centroid
    p_apex_sub, p_base_sub = _disambiguate_base_apex_endpoints(
        sub_mask, ep1_sub, ep2_sub, rosette_centroid=rc_sub
    )
    p_apex = (p_apex_sub[0] + x0_pad, p_apex_sub[1] + y0_pad)
    p_base = (p_base_sub[0] + x0_pad, p_base_sub[1] + y0_pad)
    angle_deg = math.degrees(
        math.atan2(float(p_base[1] - p_apex[1]), float(p_base[0] - p_apex[0]))
    )

    return ucs, solidity, angle_deg, p_apex, p_base


# Alias for backwards compatibility with geometric_gatekeeper.py
compute_geometric_metrics_and_pose = compute_geometric_metrics


def is_tier1_pristine(
    ucs: float,
    solidity: float,
    min_ucs: float = 0.85,
    min_solidity: float = 0.72
) -> bool:
    """
    Evaluates whether a leaf instance satisfies Tier 1 Pristine QC thresholds.
    """
    return bool(ucs >= min_ucs and solidity >= min_solidity)


def trim_petiole_tail(
    mask: np.ndarray,
    p_base: Optional[Tuple[int, int]] = None,
    p_apex: Optional[Tuple[int, int]] = None,
    step_px: int = 2,
    max_petiole_width: float = 16.0,
    gradient_threshold: float = 0.75,
    max_search_fraction: float = 0.45,
    min_petiole_length_px: float = 14.0,
    aspect_ratio_protrusion_thresh: float = 3.0,
    return_metadata: bool = False,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> Union[np.ndarray, Tuple[np.ndarray, Dict[str, Any]]]:
    """
    Width-Inflection Junction Trimmer.

    Detects and perpendicularly trims narrow petiole stalks extending below the leaf blade
    at the petiole-lamina junction. In Packera, petioles are long and slender. If segmentation
    masks include accidental 5-15 mm segments of the petiole stalk below the blade, closed EFA
    is hyper-dominated by the narrow linear "tail", distorting harmonics A1...A3 and masking
    true blade shape variation.

    Parameters
    ----------
    mask : np.ndarray
        Binary uint8 mask of the detected leaf instance (0 and 255/nonzero).
    p_base : Optional[Tuple[int, int]]
        Detected or estimated basal insertion point (x, y). If None, auto-detected
        along primary longitudinal midrib axis.
    p_apex : Optional[Tuple[int, int]]
        Detected or estimated blade apex point (x, y). If None, auto-detected.
    step_px : int
        Step interval in pixels along the medial axis (default: 2 px).
    max_petiole_width : float
        Maximum cross-sectional width in pixels to consider a narrow petiole stalk
        (Packera petiole stalks typically maintain W ~ 5-15 px, default threshold: 22.0 px).
    gradient_threshold : float
        Width gradient inflection threshold (dW/dL > threshold) where the petiole stalk
        rapidly flares into the lamina (default: 0.75).
    max_search_fraction : float
        Fraction of longitudinal axis length to search from the base (default: 0.45).
    min_petiole_length_px : float
        Minimum longitudinal stalk length in pixels required before triggering a junction
        trim (default: 8.0 px), preventing clipping of acute blade tips on clean leaves.
    aspect_ratio_protrusion_thresh : float
        Protrusion aspect ratio (length / width) threshold in the basal 5% segment
        for the fallback morphological width check (default: 3.0).
    return_metadata : bool
        If True, returns a tuple of (trimmed_mask, metadata_dict).
    sheet_shape : Optional[Tuple[int, int]]
        Optional (sheet_h, sheet_w) reference dimensions when operating on a cropped mask.

    Returns
    -------
    Union[np.ndarray, Tuple[np.ndarray, Dict[str, Any]]]
        Trimmed binary mask, or (trimmed_mask, metadata) if return_metadata=True.
    """
    if mask is None or mask.size == 0:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": p_base,
            "p_apex": p_apex,
            "original_area": 0,
            "trimmed_area": 0,
            "stalk_width_mean": 0.0,
            "reason": "EMPTY_OR_TINY_MASK",
        }
        fallback_mask = mask.copy() if mask is not None else np.zeros((10, 10), dtype=np.uint8)
        return (fallback_mask, meta) if return_metadata else fallback_mask

    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    h_full, w_full = mask_u8.shape[:2]
    orig_h, orig_w = sheet_shape if sheet_shape is not None else (h_full, w_full)

    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if bw == 0 or bh == 0:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": p_base,
            "p_apex": p_apex,
            "original_area": 0,
            "trimmed_area": 0,
            "stalk_width_mean": 0.0,
            "reason": "EMPTY_OR_TINY_MASK",
        }
        return (mask.copy(), meta) if return_metadata else mask.copy()

    x0_pad = max(0, bx - 2)
    y0_pad = max(0, by - 2)
    x1_pad = min(w_full, bx + bw + 2)
    y1_pad = min(h_full, by + bh + 2)
    sub_mask = mask_u8[y0_pad:y1_pad, x0_pad:x1_pad].copy()
    init_area = int(np.count_nonzero(sub_mask))
    if init_area < 30:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": p_base,
            "p_apex": p_apex,
            "original_area": init_area,
            "trimmed_area": init_area,
            "stalk_width_mean": 0.0,
            "reason": "EMPTY_OR_TINY_MASK",
        }
        return (mask.copy(), meta) if return_metadata else mask.copy()

    def _to_full(local_arr: np.ndarray) -> np.ndarray:
        if x0_pad == 0 and y0_pad == 0 and x1_pad == w_full and y1_pad == h_full:
            return local_arr.copy()
        out = np.zeros((h_full, w_full), dtype=np.uint8)
        out[y0_pad:y1_pad, x0_pad:x1_pad] = local_arr
        return out

    def _shift_pt(pt: Optional[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        return (int(pt[0]) + x0_pad, int(pt[1]) + y0_pad) if pt is not None else None

    p_base_sub = (int(p_base[0]) - x0_pad, int(p_base[1]) - y0_pad) if p_base is not None else None
    p_apex_sub = (int(p_apex[0]) - x0_pad, int(p_apex[1]) - y0_pad) if p_apex is not None else None

    h, w = sub_mask.shape[:2]
    sub_mask = bridge_dissected_rachis_components(sub_mask)
    contours, _ = cv2.findContours(sub_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": p_base,
            "p_apex": p_apex,
            "original_area": init_area,
            "trimmed_area": init_area,
            "stalk_width_mean": 0.0,
            "reason": "NO_CONTOURS",
        }
        return (mask.copy(), meta) if return_metadata else mask.copy()

    main_cnt = max(contours, key=cv2.contourArea)
    pts = main_cnt.reshape(-1, 2).astype(np.float32)
    area_main = float(cv2.contourArea(pts))
    hull_main = float(cv2.contourArea(cv2.convexHull(pts)))
    solidity_sub = (area_main / hull_main) if hull_main > 0.0 else 1.0
    rect_b = cv2.minAreaRect(pts)
    blade_width_est = float(min(rect_b[1])) if min(rect_b[1]) > 0 else float(min(orig_h, orig_w))
    if blade_width_est > 250.0 and min(orig_h, orig_w) > 500:
        if max_petiole_width == 16.0:
            # On deeply lobed/dissected high-DPI blades, keep petiole width cap tighter so basal pinna pairs
            # are never mistaken for a thick petiole stalk.
            scale_pet = 0.065 if solidity_sub < 0.72 else 0.105
            max_petiole_width = max(16.0, scale_pet * blade_width_est)
        cordate_guard_w = max(float(max_petiole_width) * 1.18, 22.0, 0.165 * blade_width_est)
    else:
        cordate_guard_w = max(float(max_petiole_width), 19.0)

    # 1. Medial Longitudinal Axis via cv2.fitLine
    line_params = cv2.fitLine(main_cnt, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]
    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy

    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))
    ep1 = (int(pts[min_idx, 0]), int(pts[min_idx, 1]))
    ep2 = (int(pts[max_idx, 0]), int(pts[max_idx, 1]))

    max_local_span = int(math.ceil(math.hypot(h, w))) + 2
    base_span = max(min(max(orig_h, orig_w), 350), int(0.65 * max(orig_h, orig_w)))
    w_span_eff = min(base_span, max_local_span)

    # Vectorized helper to sample perpendicular cross-sectional widths
    def _sample_widths(pt_start: Tuple[int, int], pt_end: Tuple[int, int], sample_dist: float) -> Tuple[np.ndarray, np.ndarray]:
        v_dir = np.array([pt_end[0] - pt_start[0], pt_end[1] - pt_start[1]], dtype=np.float64)
        L = np.hypot(v_dir[0], v_dir[1])
        if L < 1e-5:
            return np.array([]), np.array([])
        u = v_dir / L
        n = np.array([-u[1], u[0]], dtype=np.float64)
        s_steps = np.arange(0, sample_dist, step_px, dtype=np.float64)
        if len(s_steps) == 0:
            return np.array([]), np.array([])
        centers = np.array(pt_start, dtype=np.float64)[None, :] + s_steps[:, None] * u[None, :]
        t_steps = np.arange(-w_span_eff, w_span_eff + 1, 1.0, dtype=np.float64)
        gx = np.round(centers[:, 0, None] + t_steps[None, :] * n[0]).astype(np.int32)
        gy = np.round(centers[:, 1, None] + t_steps[None, :] * n[1]).astype(np.int32)
        valid = (gx >= 0) & (gx < w) & (gy >= 0) & (gy < h)
        grid_mask = np.zeros((len(s_steps), len(t_steps)), dtype=bool)
        grid_mask[valid] = sub_mask[gy[valid], gx[valid]] > 0
        widths = np.zeros(len(s_steps), dtype=np.float64)
        for i in range(len(s_steps)):
            row = grid_mask[i]
            nz = np.where(row)[0]
            widths[i] = float(nz[-1] - nz[0] + 1) if len(nz) > 0 else 0.0
        return s_steps, widths

    # Resolve endpoints and orientation
    if p_base_sub is not None and p_apex_sub is not None:
        cand_base, cand_apex = p_base_sub, p_apex_sub
    else:
        # Probe both endpoints along the fitted longitudinal axis
        probe_dist = min(60.0, np.hypot(ep2[0] - ep1[0], ep2[1] - ep1[1]) * 0.35)
        s1, w1 = _sample_widths(ep1, ep2, probe_dist)
        s2, w2 = _sample_widths(ep2, ep1, probe_dist)

        narrow_1 = np.where(w1 <= max_petiole_width)[0] if len(w1) else np.array([])
        narrow_2 = np.where(w2 <= max_petiole_width)[0] if len(w2) else np.array([])
        stalk_len_1 = float(s1[narrow_1[-1]]) if len(narrow_1) else 0.0
        stalk_len_2 = float(s2[narrow_2[-1]]) if len(narrow_2) else 0.0

        if stalk_len_1 >= min_petiole_length_px and stalk_len_1 >= stalk_len_2 + 6.0:
            cand_base, cand_apex = ep1, ep2
        elif stalk_len_2 >= min_petiole_length_px and stalk_len_2 >= stalk_len_1 + 6.0:
            cand_base, cand_apex = ep2, ep1
        else:
            cand_base = p_base_sub if p_base_sub is not None else ep2
            cand_apex = p_apex_sub if p_apex_sub is not None else ep1

    axis_len = float(math.hypot(cand_apex[0] - cand_base[0], cand_apex[1] - cand_base[1]))
    if axis_len < 15.0:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": _shift_pt(cand_base),
            "p_apex": _shift_pt(cand_apex),
            "original_area": init_area,
            "trimmed_area": init_area,
            "stalk_width_mean": 0.0,
            "reason": "AXIS_TOO_SHORT",
        }
        return (mask.copy(), meta) if return_metadata else mask.copy()

    # Step along medial axis from basal insertion toward apex at step_px intervals
    search_dist = min(axis_len * max_search_fraction, axis_len - 10.0)
    s_steps, widths = _sample_widths(cand_base, cand_apex, search_dist)
    if len(s_steps) < 5:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": _shift_pt(cand_base),
            "p_apex": _shift_pt(cand_apex),
            "original_area": init_area,
            "trimmed_area": init_area,
            "stalk_width_mean": 0.0,
            "reason": "INSUFFICIENT_STEPS",
        }
        return (mask.copy(), meta) if return_metadata else mask.copy()

    # 2. Guard against truncate / cordate leaves: broad base is never clipped
    w_init = np.median(widths[:min(3, len(widths))])
    if w_init > cordate_guard_w:
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": _shift_pt(cand_base),
            "p_apex": _shift_pt(cand_apex),
            "original_area": init_area,
            "trimmed_area": init_area,
            "stalk_width_mean": 0.0,
            "reason": "CORDATE_OR_BROAD_BASE_PROTECTED",
        }
        return (mask.copy(), meta) if return_metadata else mask.copy()

    # Find contiguous narrow petiole stalk segment (W <= max_petiole_width)
    # On lobed/dissected leaves (solidity_sub < 0.72), halt the petiole stalk search immediately
    # if the width profile flares into a basal lateral pinna lobe (>= 2.1x initial stalk width),
    # preventing the search from crossing into an inter-lobe rachis internode.
    lobe_flare_ceil = (
        min(float(max_petiole_width), max(18.0, 2.1 * float(max(w_init, 6.0))))
        if solidity_sub < 0.72
        else float(max_petiole_width)
    )
    narrow_indices = []
    for k in range(len(s_steps)):
        if widths[k] <= lobe_flare_ceil:
            narrow_indices.append(k)
        else:
            break

    num_l0, lbl0, st0, _ = cv2.connectedComponentsWithStats((sub_mask > 0).astype(np.uint8), connectivity=8)
    if num_l0 > 2:
        best_l0 = 1 + int(np.argmax(st0[1:, cv2.CC_STAT_AREA]))
        sub_mask = np.where(lbl0 == best_l0, 255, 0).astype(np.uint8)
    cnts0, _ = cv2.findContours(sub_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if cnts0:
        m_fill = np.zeros_like(sub_mask)
        cv2.drawContours(m_fill, [max(cnts0, key=cv2.contourArea)], -1, 255, -1)
        sub_mask = m_fill

    cleaned_area = int(np.count_nonzero(sub_mask))

    if not narrow_indices:
        out_unclipped = _to_full(sub_mask)
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": _shift_pt(cand_base),
            "p_apex": _shift_pt(cand_apex),
            "original_area": cleaned_area,
            "trimmed_area": cleaned_area,
            "stalk_width_mean": 0.0,
            "reason": "NO_NARROW_BASE",
        }
        return (out_unclipped, meta) if return_metadata else out_unclipped

    stalk_len = s_steps[narrow_indices[-1]]
    junction_idx = None
    reason = None

    # Fallback Morphological Width Check if stalk is short
    if stalk_len < min_petiole_length_px:
        idx_5 = np.where(s_steps <= 0.05 * axis_len)[0]
        if len(idx_5) >= 2:
            w_5 = np.mean(widths[idx_5])
            L_5 = s_steps[idx_5[-1]]
            ar = L_5 / max(w_5, 1.0)
            if ar > aspect_ratio_protrusion_thresh and w_5 <= lobe_flare_ceil and stalk_len >= 6.0:
                junction_idx = narrow_indices[-1]
                reason = "FALLBACK_MORPHOLOGICAL_PROTRUSION_TRIMMED"
            else:
                reason = "CLEAN_BLADE_UNCLIPPED"
        else:
            reason = "CLEAN_BLADE_UNCLIPPED"
    else:
        # A narrow petiole stalk exists; search for width inflection junction into lamina
        req_flare_w = lobe_flare_ceil * 1.25 if solidity_sub < 0.72 else max_petiole_width * 1.3
        for k in narrow_indices:
            if s_steps[k] < min_petiole_length_px:
                continue
            span = min(3, len(s_steps) - 1 - k)
            if span > 0:
                dw_dl = (widths[k + span] - widths[k]) / (s_steps[k + span] - s_steps[k])
                if dw_dl >= gradient_threshold:
                    future_w = widths[k: min(k + 8, len(widths))]
                    if np.max(future_w) >= req_flare_w:
                        j_idx = k if solidity_sub < 0.72 else max(k, narrow_indices[-1])
                        if (
                            solidity_sub >= 0.72
                            and j_idx + 2 < len(widths)
                            and widths[j_idx + 1] <= max_petiole_width * 1.25
                            and widths[j_idx + 2] >= max_petiole_width * 1.5
                        ):
                            j_idx += 2
                        junction_idx = j_idx
                        reason = "WIDTH_INFLECTION_JUNCTION_TRIMMED"
                        break
        if junction_idx is None:
            last_narrow = narrow_indices[-1]
            if last_narrow < len(widths) - 2 and widths[last_narrow + 1] > lobe_flare_ceil:
                junction_idx = last_narrow
                reason = "WIDTH_INFLECTION_JUNCTION_TRIMMED"

    if junction_idx is None or junction_idx == 0:
        out_unclipped = _to_full(sub_mask)
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": _shift_pt(cand_base),
            "p_apex": _shift_pt(cand_apex),
            "original_area": cleaned_area,
            "trimmed_area": cleaned_area,
            "stalk_width_mean": 0.0,
            "reason": reason or "CLEAN_BLADE_UNCLIPPED",
        }
        return (out_unclipped, meta) if return_metadata else out_unclipped

    # 3. Clip mask perpendicularly at detected junction
    L_junction = s_steps[junction_idx]
    u_dir = np.array([cand_apex[0] - cand_base[0], cand_apex[1] - cand_base[1]], dtype=np.float64) / axis_len
    p_junction = np.array(cand_base, dtype=np.float64) + L_junction * u_dir
    pj_x, pj_y = p_junction[0], p_junction[1]

    Y, X = np.ogrid[:h, :w]
    proj = (X - pj_x) * u_dir[0] + (Y - pj_y) * u_dir[1]
    clipped = sub_mask.copy()
    clipped[proj < 0] = 0

    if solidity_sub < 0.72:
        clipped = bridge_dissected_rachis_components(clipped)

    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(clipped, connectivity=8)
    if num_labels > 1:
        largest_label = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
        clipped = np.where(labels == largest_label, 255, 0).astype(np.uint8)
    cnts_ext, _ = cv2.findContours(clipped, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if cnts_ext:
        clipped_filled = np.zeros_like(clipped)
        cv2.drawContours(clipped_filled, [max(cnts_ext, key=cv2.contourArea)], -1, 255, -1)
        clipped = clipped_filled

    trimmed_area_val = int(np.count_nonzero(clipped))
    if solidity_sub < 0.72 and trimmed_area_val < int(0.65 * cleaned_area):
        out_unclipped = _to_full(sub_mask)
        meta = {
            "is_trimmed": False,
            "trim_distance_px": 0.0,
            "junction_pt": None,
            "p_base": _shift_pt(cand_base),
            "p_apex": _shift_pt(cand_apex),
            "original_area": cleaned_area,
            "trimmed_area": cleaned_area,
            "stalk_width_mean": 0.0,
            "reason": "DISSECTED_BLADE_PROTECTED",
        }
        return (out_unclipped, meta) if return_metadata else out_unclipped

    j_pt_full = (int(round(pj_x)) + x0_pad, int(round(pj_y)) + y0_pad)
    out_clipped = _to_full(clipped)
    meta = {
        "is_trimmed": True,
        "trim_distance_px": float(L_junction),
        "junction_pt": j_pt_full,
        "p_base": j_pt_full,
        "p_apex": _shift_pt(cand_apex),
        "original_area": cleaned_area,
        "trimmed_area": int(np.count_nonzero(clipped)),
        "stalk_width_mean": float(np.mean(widths[s_steps < L_junction])) if np.any(s_steps < L_junction) else 0.0,
        "reason": reason,
    }
    return (out_clipped, meta) if return_metadata else out_clipped


def homologize_contour_starting_point(
    contour: np.ndarray,
    petiole_attachment_pt: Tuple[int, int]
) -> np.ndarray:
    """
    Homologizes the contour starting point to an anatomical anchor (petiole attachment
    or basal insertion point) and enforces strict clockwise perimeter tracing.

    In Elliptic Fourier Analysis (e.g. Momocs::efourier), if contours start at arbitrary
    perimeter locations (e.g., base vs. apex), Fourier phase angles shift by ~pi (180 deg),
    flipping harmonic signs and artificially splitting identical shapes in PCA morphospace.

    Parameters
    ----------
    contour : np.ndarray
        Closed 2D contour array of shape (N, 2) or (N, 1, 2).
    petiole_attachment_pt : Tuple[int, int]
        (x, y) coordinates of the petiole-blade junction or base-most point along
        the primary longitudinal axis.

    Returns
    -------
    np.ndarray
        Homologized 2D contour array of shape (N, 2) starting at the vertex closest
        to petiole_attachment_pt and tracing strictly clockwise.
    """
    if contour is None or len(contour) == 0:
        return np.empty((0, 2), dtype=np.float32)

    pts_2d = np.asarray(contour).reshape(-1, 2)
    n_pts = len(pts_2d)

    if n_pts < 3 or petiole_attachment_pt is None:
        return pts_2d.copy()

    # 1. Calculate Euclidean distances from all contour vertices to petiole_attachment_pt
    p_anchor = np.array([float(petiole_attachment_pt[0]), float(petiole_attachment_pt[1])], dtype=np.float64)
    diffs = pts_2d.astype(np.float64) - p_anchor
    distances = np.hypot(diffs[:, 0], diffs[:, 1])

    # Find closest vertex index to petiole attachment point
    start_idx = int(np.argmin(distances))

    # Rotate (roll) coordinate array so start_idx becomes Index 0
    homologized_contour = np.roll(pts_2d, -start_idx, axis=0)

    # 2. Enforce Strict Clockwise Tracing Orientation
    # In OpenCV, positive area indicates counter-clockwise orientation; negative indicates clockwise.
    area = cv2.contourArea(homologized_contour.reshape(-1, 1, 2).astype(np.float32), oriented=True)

    # Guard against zero-area or degenerated contours (area == 0.0 remains unchanged)
    if area > 0:
        # Invert vertex order while preserving Index 0
        homologized_contour = np.vstack([homologized_contour[0:1], homologized_contour[1:][::-1]])

    return homologized_contour


def _align_and_normalize_contour(
    x_resamp: np.ndarray,
    y_resamp: np.ndarray,
    p_base: Tuple[int, int],
    p_apex: Tuple[int, int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Rotates resampled contour coordinates (x_resamp, y_resamp) around Index 0 (p_base)
    so the longitudinal midrib vector (p_base -> p_apex) lies along the positive +X axis,
    and normalizes x_norm and y_norm isotropically by dividing BOTH axes by max(span_x, span_y).
    """
    dx_axis = float(p_apex[0] - p_base[0])
    dy_axis = float(p_apex[1] - p_base[1])
    if math.hypot(dx_axis, dy_axis) < 1e-5 and len(x_resamp) > 2:
        dists = np.hypot(x_resamp - x_resamp[0], y_resamp - y_resamp[0])
        far_idx = int(np.argmax(dists))
        dx_axis = float(x_resamp[far_idx] - x_resamp[0])
        dy_axis = float(y_resamp[far_idx] - y_resamp[0])

    if math.hypot(dx_axis, dy_axis) >= 1e-5:
        theta = math.atan2(dy_axis, dx_axis)
        cos_t = math.cos(-theta)
        sin_t = math.sin(-theta)
        bx = float(x_resamp[0])
        by = float(y_resamp[0])
        dx_c = x_resamp - bx
        dy_c = y_resamp - by
        x_rot = bx + dx_c * cos_t - dy_c * sin_t
        y_rot = by + dx_c * sin_t + dy_c * cos_t
    else:
        x_rot = x_resamp.copy()
        y_rot = y_resamp.copy()

    min_x, max_x = float(np.min(x_rot)), float(np.max(x_rot))
    min_y, max_y = float(np.min(y_rot)), float(np.max(y_rot))
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)
    max_span = max(span_x, span_y, 1e-6)

    x_norm = np.round((x_rot - min_x) / max_span, 6)
    y_norm = np.round((y_rot - min_y) / max_span, 6)
    return x_rot, y_rot, x_norm, y_norm


def resample_normalized_contour(
    mask: np.ndarray,
    num_points: int = 120,
    petiole_attachment_pt: Optional[Tuple[int, int]] = None
) -> Optional[np.ndarray]:
    """
    Extracts resampled, normalized 2D coordinates (x_norm, y_norm) of shape (num_points, 2)
    from a binary mask, with consistent clockwise orientation, midrib +X alignment,
    and isotropic normalization by max(span_x, span_y).
    """
    if mask is not None:
        mask = trim_petiole_tail(mask)

    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return None

    largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2)
    if len(largest_cnt) < 10:
        return None

    _, _, _, p_apex_m, p_base_m = compute_geometric_metrics(mask)
    if petiole_attachment_pt is None:
        petiole_attachment_pt = p_base_m
        p_apex_pt = p_apex_m
    else:
        d_to_apex = math.hypot(petiole_attachment_pt[0] - p_apex_m[0], petiole_attachment_pt[1] - p_apex_m[1])
        d_to_base = math.hypot(petiole_attachment_pt[0] - p_base_m[0], petiole_attachment_pt[1] - p_base_m[1])
        p_apex_pt = p_apex_m if d_to_apex >= d_to_base else p_base_m

    largest_cnt = homologize_contour_starting_point(largest_cnt, petiole_attachment_pt)

    pts = largest_cnt.astype(np.float64)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    valid = dt > 1e-7
    if np.sum(valid) < 5:
        return None

    pts = pts[np.concatenate([[True], valid])]
    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    cum_dist = np.concatenate([[0.0], np.cumsum(dt)])
    total_len = cum_dist[-1]
    if total_len <= 0:
        return None

    target_distances = np.linspace(0.0, total_len, num_points, endpoint=False)
    x_resamp = np.interp(target_distances, cum_dist, pts[:, 0])
    y_resamp = np.interp(target_distances, cum_dist, pts[:, 1])

    # Enforce strict clockwise tracing orientation & preserve anchor
    resampled_homologized = homologize_contour_starting_point(
        np.column_stack([x_resamp, y_resamp]), petiole_attachment_pt
    )
    x_resamp = resampled_homologized[:, 0]
    y_resamp = resampled_homologized[:, 1]

    _, _, x_norm, y_norm = _align_and_normalize_contour(
        x_resamp, y_resamp, petiole_attachment_pt, p_apex_pt
    )
    return np.column_stack([x_norm, y_norm])


# =============================================================================
# 2. Vectorized Hemi-Blade Reflection (Tier 2) & Tier 1 Extraction
# =============================================================================

def extract_tier1_pristine(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    output_dir: Path
) -> str:
    """
    Tier 1: Direct Pristine Silhouette extraction.
    Saves binary mask to output_dir/masks/ and output_dir/masks/tier1_intact/.
    """
    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    mask_dir = out_base / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    tier1_dir = mask_dir / "tier1_intact"
    tier1_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{catalog_number}_leaf{leaf_id}.png"
    main_save_path = mask_dir / filename
    tier_save_path = tier1_dir / filename

    bin_mask = (mask > 0).astype(np.uint8) * 255
    ok, png_buf = cv2.imencode(".png", bin_mask)
    if ok:
        png_bytes = png_buf.tobytes()
        main_save_path.write_bytes(png_bytes)
        tier_save_path.write_bytes(png_bytes)
    else:
        cv2.imwrite(str(main_save_path), bin_mask)
        shutil.copyfile(str(main_save_path), str(tier_save_path))
    return str(main_save_path)


# Backward compatibility alias
route_tier1_silhouette = extract_tier1_pristine


def extract_tier2_reflected(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    p_apex: Tuple[int, int],
    p_base: Tuple[int, int],
    output_dir: Optional[Path] = None,
    threshold_defect: float = 35.0,
    min_solidity: float = 0.68,
    min_ucs: float = 0.50
) -> Tuple[Optional[str], Optional[np.ndarray]]:
    """
    Tier 2: Vectorized Hemi-Blade Bilateral Symmetry Reflection.
    Crops around the tight bounding box with 25% padding, aligns mask along the
    midrib axis using +angle (so the midrib rotates to 0 degrees horizontal in
    OpenCV's top-left origin frame), assesses margin defect profile of upper vs
    lower half, reflects the clean half across the midrib axis symmetrically
    using NumPy vectorization (no nested pixel loops), applies morphological
    closing to produce a watertight closed silhouette, and inverse-rotates back
    to the original canvas frame.
    """
    h, w = mask.shape[:2]
    if h < 15 or w < 15:
        return None, None

    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if bw == 0 or bh == 0:
        return None, None

    x_min, x_max = bx, bx + bw - 1
    y_min, y_max = by, by + bh - 1
    pad = max(16, int(0.25 * max(bw, bh)))

    tight_crop = (mask_u8[y_min:y_max + 1, x_min:x_max + 1] > 0).astype(np.uint8) * 255
    if int(np.count_nonzero(tight_crop)) < 30:
        return None, None

    if p_apex == p_base or p_apex == (0, 0):
        _, _, _, p_apex, p_base = compute_geometric_metrics(mask_u8)

    cnts, _ = cv2.findContours(tight_crop, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        return None, None
    in_cnt = max(cnts, key=cv2.contourArea)
    in_cnt_f = in_cnt.astype(np.float32)
    in_rect = cv2.minAreaRect(in_cnt_f)
    in_minor_w = float(min(in_rect[1]))
    in_major_l = float(max(in_rect[1]))
    in_ar = in_minor_w / max(in_major_l, 1.0)

    in_area = float(cv2.contourArea(in_cnt_f))
    in_hull = cv2.convexHull(in_cnt_f)
    in_hull_area = float(cv2.contourArea(in_hull))
    in_solidity = float(in_area / in_hull_area) if in_hull_area > 0 else 0.0

    is_high_dpi_mask = max(bw, bh) > 160
    if is_high_dpi_mask and (in_solidity < 0.68 or in_ar > 0.84):
        return None, None

    crop_mask = cv2.copyMakeBorder(
        tight_crop, pad, pad, pad, pad, cv2.BORDER_CONSTANT, value=0
    )
    ch, cw = crop_mask.shape[:2]

    ax_c = float(p_apex[0] - x_min + pad)
    ay_c = float(p_apex[1] - y_min + pad)
    bx_c = float(p_base[0] - x_min + pad)
    by_c = float(p_base[1] - y_min + pad)

    # 2. Compute midrib angle from p_apex to p_base and pass +angle to rotate to 0 deg horizontal
    dx = float(p_base[0] - p_apex[0])
    dy = float(p_base[1] - p_apex[1])
    angle = math.degrees(math.atan2(dy, dx))

    cx = float(round(0.5 * (ax_c + bx_c)))
    cy = float(round(0.5 * (ay_c + by_c)))
    if not (2.0 < cx < cw - 3.0 and 2.0 < cy < ch - 3.0):
        cx, cy = cw / 2.0, ch / 2.0

    rot_mat = cv2.getRotationMatrix2D((cx, cy), angle, 1.0)
    aligned = cv2.warpAffine(crop_mask, rot_mat, (cw, ch), flags=cv2.INTER_NEAREST)

    y_indices, x_indices = np.where(aligned > 0)
    if len(y_indices) < 30:
        return None, None

    y_midrib = int(round(cy))
    if y_midrib <= 1 or y_midrib >= ch - 2:
        y_midrib = int(np.median(y_indices))
    if y_midrib <= 1 or y_midrib >= ch - 2:
        return None, None

    upper_half = np.zeros_like(aligned)
    lower_half = np.zeros_like(aligned)
    upper_half[:y_midrib, :] = aligned[:y_midrib, :]
    lower_half[y_midrib:, :] = aligned[y_midrib:, :]

    # Fully vectorized margin profile detection across active leaf columns
    upper_has_fg = upper_half > 0
    lower_has_fg = lower_half > 0

    col_has_upper = np.any(upper_has_fg, axis=0)
    col_has_lower = np.any(lower_has_fg, axis=0)
    active_cols = np.where(col_has_upper | col_has_lower)[0]
    if len(active_cols) < 5:
        return None, None
    c_start, c_end = int(active_cols[0]), int(active_cols[-1])

    # Fast column min/max using vectorized argmax
    upper_profile = np.where(col_has_upper, np.argmax(upper_has_fg, axis=0), y_midrib)[c_start:c_end + 1]
    lower_profile = np.where(col_has_lower, ch - 1 - np.argmax(lower_has_fg[::-1, :], axis=0), y_midrib)[c_start:c_end + 1]

    # Assess margin defect score (jumps/discontinuities)
    upper_diff = np.abs(np.diff(upper_profile))
    lower_diff = np.abs(np.diff(lower_profile))
    upper_defects = float(np.sum(upper_diff > 8) * 10.0 + np.std(upper_diff))
    lower_defects = float(np.sum(lower_diff > 8) * 10.0 + np.std(lower_diff))

    upper_area = int(np.count_nonzero(upper_half))
    lower_area = int(np.count_nonzero(lower_half))

    if upper_defects < threshold_defect and lower_defects < threshold_defect:
        if abs(upper_defects - lower_defects) < 5.0 and upper_area != lower_area:
            selected_half = "upper" if upper_area > lower_area else "lower"
        elif upper_defects <= lower_defects:
            selected_half = "upper"
        else:
            selected_half = "lower"
    elif upper_defects < threshold_defect:
        selected_half = "upper"
    elif lower_defects < threshold_defect:
        selected_half = "lower"
    else:
        return None, None

    # Fully vectorized slice reflection in NumPy (zero manual pixel loops)
    reflected = np.zeros_like(aligned)
    max_dy = min(y_midrib, ch - y_midrib)
    if selected_half == "upper":
        reflected[:y_midrib, :] = upper_half[:y_midrib, :]
        if max_dy > 0:
            reflected[y_midrib:y_midrib + max_dy, :] = np.flip(
                upper_half[y_midrib - max_dy:y_midrib, :], axis=0
            )
    else:
        reflected[y_midrib:, :] = lower_half[y_midrib:, :]
        if max_dy > 0:
            reflected[y_midrib - max_dy:y_midrib, :] = np.flip(
                lower_half[y_midrib:y_midrib + max_dy, :], axis=0
            )

    # Morphological closing to seal midrib seam in horizontal aligned space
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    reflected = cv2.morphologyEx(reflected, cv2.MORPH_CLOSE, kernel)

    # Inverse-rotate synthesized silhouette back (-angle) into original pose
    inv_rot_mat = cv2.getRotationMatrix2D((cx, cy), -angle, 1.0)
    unaligned_crop = cv2.warpAffine(reflected, inv_rot_mat, (cw, ch), flags=cv2.INTER_NEAREST)
    unaligned_crop = cv2.morphologyEx(unaligned_crop, cv2.MORPH_CLOSE, kernel)

    # Evaluate and clean in tight ROI coordinates before allocating full (h, w) canvas
    src_x0 = max(0, pad - x_min)
    src_y0 = max(0, pad - y_min)
    dst_x0 = max(0, x_min - pad)
    dst_y0 = max(0, y_min - pad)
    copy_w = min(cw - src_x0, w - dst_x0)
    copy_h = min(ch - src_y0, h - dst_y0)
    if copy_w > 0 and copy_h > 0:
        roi_bin = (
            unaligned_crop[src_y0:src_y0 + copy_h, src_x0:src_x0 + copy_w] > 127
        ).astype(np.uint8) * 255
    else:
        roi_bin = (unaligned_crop > 127).astype(np.uint8) * 255

    num_lr, lbl_r, st_r, _ = cv2.connectedComponentsWithStats((roi_bin > 0).astype(np.uint8), connectivity=8)
    if num_lr > 1:
        best_lr = 1 + int(np.argmax(st_r[1:, cv2.CC_STAT_AREA]))
        roi_bin = np.where(lbl_r == best_lr, 255, 0).astype(np.uint8)
    cnts_r, _ = cv2.findContours(roi_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if cnts_r:
        ref_filled = np.zeros_like(roi_bin)
        cv2.drawContours(ref_filled, [max(cnts_r, key=cv2.contourArea)], -1, 255, -1)
        roi_bin = ref_filled

    # Verify synthesized silhouette quality on tight ROI
    ref_ucs, ref_solidity, _, _, _ = compute_geometric_metrics(roi_bin)
    if ref_solidity < min_solidity or ref_ucs < min_ucs:
        return None, None

    if is_high_dpi_mask:
        ref_cnts, _ = cv2.findContours(roi_bin, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not ref_cnts:
            return None, None
        ref_cnt = max(ref_cnts, key=cv2.contourArea)
        ref_cnt_f = ref_cnt.astype(np.float32)
        ref_rect = cv2.minAreaRect(ref_cnt_f)
        ref_minor_w = float(min(ref_rect[1]))
        ref_major_l = float(max(ref_rect[1]))
        ref_ar = ref_minor_w / max(ref_major_l, 1.0)

        if ref_ar > 0.84 or ref_solidity < 0.78 or ref_ucs < 0.82:
            return None, None

    if copy_w > 0 and copy_h > 0:
        reflected_bin = np.zeros((h, w), dtype=np.uint8)
        reflected_bin[dst_y0:dst_y0 + copy_h, dst_x0:dst_x0 + copy_w] = roi_bin
    else:
        reflected_bin = roi_bin

    main_save_path = None
    if output_dir is not None:
        out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
        mask_dir = out_base / "masks"
        mask_dir.mkdir(parents=True, exist_ok=True)
        tier2_dir = mask_dir / "tier2_reflected"
        tier2_dir.mkdir(parents=True, exist_ok=True)

        filename = f"{catalog_number}_leaf{leaf_id}.png"
        main_save_path = mask_dir / filename
        tier_save_path = tier2_dir / filename

        cv2.imwrite(str(main_save_path), reflected_bin)
        shutil.copyfile(str(main_save_path), str(tier_save_path))

    return str(main_save_path) if main_save_path else None, reflected_bin


def route_tier2_reflected(
    mask: np.ndarray,
    candidate: Any,
    output_dir: Path,
    p_apex: Tuple[int, int],
    p_base: Tuple[int, int]
) -> Tuple[Optional[str], bool]:
    """
    Backward-compatible adapter for route_tier2_reflected signature.
    """
    cat_num = getattr(candidate, "catalog_number", "UNKNOWN")
    leaf_id = getattr(candidate, "leaf_id", 1)
    path, reflected_mask = extract_tier2_reflected(
        mask=mask,
        catalog_number=cat_num,
        leaf_id=leaf_id,
        p_apex=p_apex,
        p_base=p_base,
        output_dir=output_dir
    )
    return path, (path is not None and reflected_mask is not None)


# =============================================================================
# 3. Standardized Coordinate Contour Exporter
# =============================================================================

def export_standardized_contour(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    output_dir: Path,
    num_points: int = 120,
    petiole_attachment_pt: Optional[Tuple[int, int]] = None,
    contour: Optional[np.ndarray] = None,
    p_apex_pt: Optional[Tuple[int, int]] = None,
    return_coords: bool = False,
) -> Any:
    """
    Extracts vectorized 2D (x, y) boundary coordinates from a binary silhouette mask,
    homologizes starting point to the petiole-blade junction anchor, enforces consistent
    clockwise orientation, aligns the midrib vector (p_base -> p_apex) along +X,
    normalizes coordinates (x_norm, y_norm) isotropically by max(span_x, span_y) for
    downstream Elliptic Fourier Analysis, and saves to
    output_dir/contours/{catalogNumber}_leaf{id}.csv.
    """
    empty_ret: Any = (None, None) if return_coords else None
    if contour is not None and len(contour) >= 10:
        largest_cnt = contour.reshape(-1, 2)
        if petiole_attachment_pt is not None and p_apex_pt is not None:
            p_apex_m, p_base_m = p_apex_pt, petiole_attachment_pt
        elif mask is not None and np.count_nonzero(mask) >= 30:
            _, _, _, p_apex_m, p_base_m = compute_geometric_metrics(mask)
        else:
            p_apex_m, p_base_m = (0, 0), (0, 0)
    else:
        if mask is not None:
            mask = trim_petiole_tail(mask)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            return empty_ret
        largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2)
        if len(largest_cnt) < 10:
            return empty_ret
        _, _, _, p_apex_m, p_base_m = compute_geometric_metrics(mask)

    # Resolve petiole attachment anchor point and apex point
    if petiole_attachment_pt is None:
        petiole_attachment_pt = p_base_m
        resolved_apex_pt = p_apex_m if p_apex_pt is None else p_apex_pt
    elif p_apex_pt is not None:
        resolved_apex_pt = p_apex_pt
    else:
        d_to_apex = math.hypot(petiole_attachment_pt[0] - p_apex_m[0], petiole_attachment_pt[1] - p_apex_m[1])
        d_to_base = math.hypot(petiole_attachment_pt[0] - p_base_m[0], petiole_attachment_pt[1] - p_base_m[1])
        resolved_apex_pt = p_apex_m if d_to_apex >= d_to_base else p_base_m
    p_apex_pt = resolved_apex_pt

    # Homologize starting point to petiole base & enforce clockwise orientation
    largest_cnt = homologize_contour_starting_point(largest_cnt, petiole_attachment_pt)

    # Resample smoothly along perimeter to standardized point count
    pts = largest_cnt.astype(np.float64)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    valid = dt > 1e-7
    if np.sum(valid) < 5:
        return empty_ret

    pts = pts[np.concatenate([[True], valid])]
    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    cum_dist = np.concatenate([[0.0], np.cumsum(dt)])
    total_len = cum_dist[-1]
    if total_len <= 0:
        return empty_ret

    target_distances = np.linspace(0.0, total_len, num_points, endpoint=False)
    x_resamp = np.interp(target_distances, cum_dist, pts[:, 0])
    y_resamp = np.interp(target_distances, cum_dist, pts[:, 1])

    # Re-verify homologization on resampled contour points
    resampled_homologized = homologize_contour_starting_point(
        np.column_stack([x_resamp, y_resamp]), petiole_attachment_pt
    )
    x_resamp = resampled_homologized[:, 0]
    y_resamp = resampled_homologized[:, 1]

    # Rotate so vector p_base -> p_apex lies along +X (anchored at p_base) and normalize isotropically
    x_rot, y_rot, x_norm, y_norm = _align_and_normalize_contour(
        x_resamp, y_resamp, petiole_attachment_pt, p_apex_pt
    )

    x_r = np.round(x_rot, 2).tolist()
    y_r = np.round(y_rot, 2).tolist()
    xn_l = np.asarray(x_norm, dtype=np.float64).tolist()
    yn_l = np.asarray(y_norm, dtype=np.float64).tolist()

    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    contour_dir = out_base / "contours"
    contour_dir.mkdir(parents=True, exist_ok=True)
    out_csv = contour_dir / f"{catalog_number}_leaf{leaf_id}.csv"
    prefix = f"{catalog_number},{leaf_id},"
    lines = ["catalogNumber,leaf_id,point_index,x,y,x_norm,y_norm\n"]
    lines.extend(
        f"{prefix}{idx},{xv},{yv},{xnv},{ynv}\n"
        for idx, (xv, yv, xnv, ynv) in enumerate(zip(x_r, y_r, xn_l, yn_l))
    )
    out_csv.write_text("".join(lines), encoding="utf-8")
    if return_coords:
        xy_coords = [[round(float(xv), 4), round(float(yv), 4)] for xv, yv in zip(x_r, y_r)]
        return str(out_csv), xy_coords
    return str(out_csv)


# =============================================================================
# Re-exported Gatekeeper & Organ/Crown Symbols for Backward Compatibility
# =============================================================================

from scripts.vision.leaf_gatekeeper import (  # noqa: E402
    BoundaryContinuityResult,
    CANONICAL_TIER_MAP,
    FoldDetectionResult,
    LeafRoutingResult,
    TIER_1_DIRECT,
    TIER_2_REFLECTED,
    _has_artificial_straight_occlusion_cut,
    bridge_dissected_rachis_components,
    check_boundary_tissue_continuity,
    classify_and_route_leaf,
    compute_classical_shape_descriptors,
    compute_sheet_adaptive_min_blade_area_cm2,
    detect_folded_leaf,
    evaluate_leaf_quality,
    is_botanical_dissection,
    is_cauline_or_stem_contaminant,
    normalize_assigned_tier,
    reflect_folded_mask,
)
from scripts.vision.organ_crown_metrics import (  # noqa: E402
    SAM2PromptBundle,
    _estimate_single_plant_crown,
    _extract_xy_point,
    _parse_basal_leaf_pose,
    approximate_petiole_and_scape_metrics,
    cluster_plant_individuals,
    cluster_voucher_plants_dbscan,
    crop_dense_rosette_tier4,
    detect_ruler_scale_hough,
    estimate_rosette_crown_point,
    generate_sam2_prompts_from_lm2_mask,
    is_topologically_connected_mask,
    route_tier3_open_curve,
)

__all__ = [
    "_sample_cross_sectional_widths",
    "_disambiguate_base_apex_endpoints",
    "compute_geometric_metrics",
    "compute_geometric_metrics_and_pose",
    "is_tier1_pristine",
    "BoundaryContinuityResult",
    "_has_artificial_straight_occlusion_cut",
    "check_boundary_tissue_continuity",
    "compute_classical_shape_descriptors",
    "reflect_folded_mask",
    "FoldDetectionResult",
    "detect_folded_leaf",
    "bridge_dissected_rachis_components",
    "is_botanical_dissection",
    "trim_petiole_tail",
    "homologize_contour_starting_point",
    "_align_and_normalize_contour",
    "resample_normalized_contour",
    "TIER_1_DIRECT",
    "TIER_2_REFLECTED",
    "CANONICAL_TIER_MAP",
    "normalize_assigned_tier",
    "compute_sheet_adaptive_min_blade_area_cm2",
    "is_cauline_or_stem_contaminant",
    "LeafRoutingResult",
    "classify_and_route_leaf",
    "evaluate_leaf_quality",
    "extract_tier1_pristine",
    "route_tier1_silhouette",
    "extract_tier2_reflected",
    "route_tier2_reflected",
    "export_standardized_contour",
    "SAM2PromptBundle",
    "generate_sam2_prompts_from_lm2_mask",
    "is_topologically_connected_mask",
    "cluster_plant_individuals",
    "cluster_voucher_plants_dbscan",
    "detect_ruler_scale_hough",
    "route_tier3_open_curve",
    "crop_dense_rosette_tier4",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
]
