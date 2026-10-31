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


class BoundaryContinuityResult(dict):
    """
    Dictionary-compatible and boolean-convertible result from boundary-continuity
    lobe rejection (`check_boundary_tissue_continuity`).
    """

    def __bool__(self) -> bool:
        return bool(self.get("is_severed_lobe", False))

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, bool):
            return bool(self) == other
        return super().__eq__(other)

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(f"'BoundaryContinuityResult' object has no attribute '{name}'") from exc


def _has_artificial_straight_occlusion_cut(
    mask_u8: np.ndarray,
) -> Tuple[bool, float]:
    """
    Detects razor-straight non-biological truncation edges (e.g., leaf blades partially
    hidden behind a ColorChecker Passport case, herbarium label, barcode sticker, or
    mounting-tape strip) on native-DPI candidate masks.

    Single longitudinal fold creases are NOT flagged here (those have 1 axial straight
    chord and are handled by `detect_folded_leaf`); this check specifically flags:
      - >= 2 distinct non-collinear razor-straight edges totaling >= 26% of perimeter, OR
      - a transverse razor-straight cut (>= 32 deg off the longitudinal axis) spanning
        >= 18% of perimeter (>= 75 px) with orthogonal RMS line-fit error <= 1.40 px.
    """
    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if max(bw, bh) < 140 or min(bw, bh) < 40:
        return False, 0.0

    cnts, _ = cv2.findContours(
        mask_u8[by : by + bh, bx : bx + bw],
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_NONE,
    )
    if not cnts:
        return False, 0.0
    raw_cnt = max(cnts, key=cv2.contourArea)
    n_pts = len(raw_cnt)
    if n_pts < 80:
        return False, 0.0

    area_px = float(cv2.contourArea(raw_cnt))
    hull_area = float(cv2.contourArea(cv2.convexHull(raw_cnt)))
    solidity = (area_px / hull_area) if hull_area > 0.0 else 0.0
    if solidity < 0.62:
        return False, 0.0

    perim = float(cv2.arcLength(raw_cnt, True))
    if perim <= 1.0:
        return False, 0.0

    pts = raw_cnt.reshape(-1, 2).astype(np.float64)
    m = cv2.moments(raw_cnt)
    if abs(m["mu20"] - m["mu02"]) > 1e-5 or abs(m["mu11"]) > 1e-5:
        axis_theta = 0.5 * math.atan2(2.0 * m["mu11"], m["mu20"] - m["mu02"])
    else:
        axis_theta = 0.0
    ax_vx, ax_vy = math.cos(axis_theta), math.sin(axis_theta)

    max_dim = float(max(bw, bh))
    min_dim = float(min(bw, bh))
    # Physical straight edges (ColorChecker rims, paper labels, barcode stickers, canvas clips)
    # exhibit fixed pixel-level linearity regardless of leaf size, whereas curved biological
    # margins accumulate larger chord sagitta on 1000-2500 px leaves.
    eps_poly = 2.4
    poly = cv2.approxPolyDP(raw_cnt, epsilon=eps_poly, closed=True).reshape(-1, 2).astype(np.float64)
    if len(poly) < 3:
        return False, 0.0

    # Map each approxPolyDP vertex back to its exact index in raw_cnt
    poly_indices: List[int] = []
    for pv in poly:
        d2 = (pts[:, 0] - pv[0]) ** 2 + (pts[:, 1] - pv[1]) ** 2
        poly_indices.append(int(np.argmin(d2)))

    min_seg_len = max(55.0, 0.16 * max_dim)
    rms_tol = 1.25
    max_tol = 2.75
    trans_rms_tol = 1.15
    trans_max_tol = 2.50
    straight_angles: List[float] = []
    straight_len_total = 0.0
    max_off_axis_deg = 0.0
    has_frame_aligned_cut = False
    has_transverse_cut = False

    for k in range(len(poly_indices)):
        i0 = poly_indices[k]
        i1 = poly_indices[(k + 1) % len(poly_indices)]
        if i0 == i1:
            continue
        p0 = pts[i0]
        p1 = pts[i1]
        dx = float(p1[0] - p0[0])
        dy = float(p1[1] - p0[1])
        chord_len = math.hypot(dx, dy)
        if chord_len < min_seg_len:
            continue

        lo_i, hi_i = (i0, i1) if i0 < i1 else (i1, i0)
        arc_direct = pts[lo_i : hi_i + 1]
        arc_wrap = np.vstack([pts[hi_i:], pts[: lo_i + 1]])
        arc_pts = arc_direct if len(arc_direct) <= len(arc_wrap) else arc_wrap
        if len(arc_pts) < 12:
            continue
        # Ensure we picked the direct arc along the polygon edge rather than the complementary loop
        arc_steps = np.hypot(np.diff(arc_pts[:, 0]), np.diff(arc_pts[:, 1]))
        arc_len = float(np.sum(arc_steps))
        if arc_len > 1.12 * chord_len:
            continue

        ux, uy = dx / chord_len, dy / chord_len
        perp = (arc_pts[:, 0] - p0[0]) * (-uy) + (arc_pts[:, 1] - p0[1]) * ux
        rms_dev = float(np.sqrt(np.mean(perp * perp)))
        max_dev = float(np.max(np.abs(perp)))

        if rms_dev <= rms_tol and max_dev <= max_tol:
            straight_len_total += chord_len
            seg_angle = math.atan2(dy, dx) % math.pi
            # Count distinct non-collinear straight cut directions (>= 15 deg apart)
            if not any(
                min(abs(seg_angle - prev_a), math.pi - abs(seg_angle - prev_a)) < math.radians(15.0)
                for prev_a in straight_angles
            ):
                straight_angles.append(seg_angle)

            deg_mod_90 = math.degrees(seg_angle) % 90.0
            if min(deg_mod_90, 90.0 - deg_mod_90) <= 6.0:
                has_frame_aligned_cut = True

            cos_align = min(1.0, abs(ux * ax_vx + uy * ax_vy))
            angle_off_axis = math.degrees(math.acos(cos_align))
            if angle_off_axis > max_off_axis_deg:
                max_off_axis_deg = angle_off_axis
            if (
                angle_off_axis >= 40.0
                and chord_len >= max(85.0, 0.34 * min_dim)
                and rms_dev <= trans_rms_tol
                and max_dev <= trans_max_tol
            ):
                has_transverse_cut = True

    straight_ratio = straight_len_total / max(2.0 * (max_dim + min_dim), 1.0)
    is_cut = bool(
        (
            len(straight_angles) >= 2
            and straight_len_total >= 0.44 * max_dim
            and (has_frame_aligned_cut or max_off_axis_deg >= 35.0)
        )
        or has_transverse_cut
    )
    return is_cut, round(min(1.0, straight_ratio), 6)


def check_boundary_tissue_continuity(
    crop_img: Optional[np.ndarray],
    binary_mask: Optional[np.ndarray],
    dilation_px: int = 8,
    max_tissue_border_ratio: float = 0.15,
    exclude_mask: Optional[np.ndarray] = None,
    petiole_base_pt: Optional[Tuple[float, float]] = None,
    inner_gap_px: Optional[int] = None,
) -> BoundaryContinuityResult:
    """
    Verifies whether a candidate leaf segmentation mask is a severed single lateral lobe
    (or internal sub-lamina / stem-chunk / label-clipped fragment) rather than an intact
    whole leaf blade:
      1. Extracts a buffered exterior perimeter ring `(dilated_outer > 0) & (dilated_inner == 0)`
         so 1-3 px PointRend serration margin conservatism and boundary anti-aliasing are not
         misclassified as exterior tissue.
      2. Masks out the trimmed petiole insertion sector around `petiole_base_pt` and any
         separately segmented distinct neighboring organ masks (`exclude_mask`).
      3. Matches exterior pixels against the candidate mask's interior lamina luminance and
         CIE L*a*b* chromaticity envelope (preventing neutral grey paper drop-shadows and
         photometrically distinct background organs from inflating `tissue_border_ratio`).
      4. Also flags multi-edge / transverse razor-straight artificial occlusions (e.g., leaves
         clipped behind ColorChecker rims, barcode stickers, or herbarium labels).
    """
    default_res = BoundaryContinuityResult({
        "is_severed_lobe": False,
        "tissue_border_ratio": 0.0,
        "ring_pixel_count": 0,
        "tissue_ring_pixels": 0,
        "otsu_threshold": 0.0,
    })
    if crop_img is None or binary_mask is None or crop_img.size == 0 or binary_mask.size == 0:
        return default_res

    mask_u8 = (binary_mask > 0).astype(np.uint8) * 255
    if int(np.count_nonzero(mask_u8)) == 0:
        return default_res

    mh, mw = mask_u8.shape[:2]
    if crop_img.shape[:2] != (mh, mw):
        img_work = cv2.resize(crop_img, (mw, mh), interpolation=cv2.INTER_AREA)
    else:
        img_work = crop_img

    bx_m, by_m, bw_m, bh_m = cv2.boundingRect(mask_u8)
    max_dim = max(bw_m, bh_m)
    if inner_gap_px is None:
        gap_px = 1 if max_dim < 280 else max(3, min(14, int(round(0.006 * max_dim))))
    else:
        gap_px = max(0, int(inner_gap_px))

    outer_rad = max(1, int(dilation_px), gap_px + max(5, int(dilation_px) - 2))
    if gap_px > 0:
        k_inner = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * gap_px + 1, 2 * gap_px + 1))
        dilated_inner = cv2.dilate(mask_u8, k_inner, iterations=1)
    else:
        dilated_inner = mask_u8

    k_outer = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * outer_rad + 1, 2 * outer_rad + 1))
    dilated_outer = cv2.dilate(mask_u8, k_outer, iterations=1)
    ring_mask = (dilated_outer > 0) & (dilated_inner == 0)

    if petiole_base_pt is not None and len(petiole_base_pt) >= 2:
        px_b = int(round(float(petiole_base_pt[0])))
        py_b = int(round(float(petiole_base_pt[1])))
        if 0 <= px_b < mw and 0 <= py_b < mh and (px_b, py_b) != (0, 0):
            pet_rad = max(10, int(round(0.08 * max_dim)))
            pet_excl = np.zeros((mh, mw), dtype=np.uint8)
            cv2.circle(pet_excl, (px_b, py_b), pet_rad, 255, -1)
            trimmed_ring = ring_mask & (pet_excl == 0)
            if int(np.count_nonzero(trimmed_ring)) >= int(0.50 * np.count_nonzero(ring_mask)):
                ring_mask = trimmed_ring

    ring_count = int(np.count_nonzero(ring_mask))
    if ring_count == 0:
        return default_res

    excl_bool: Optional[np.ndarray] = None
    if exclude_mask is not None and exclude_mask.size > 0:
        if exclude_mask.shape[:2] != (mh, mw):
            excl_rs = cv2.resize(
                (exclude_mask > 0).astype(np.uint8) * 255,
                (mw, mh),
                interpolation=cv2.INTER_NEAREST,
            )
            excl_bool = (excl_rs > 0) & (mask_u8 == 0)
        else:
            excl_bool = (exclude_mask > 0) & (mask_u8 == 0)

    lab_img: Optional[np.ndarray] = None
    if img_work.ndim == 3 and img_work.shape[2] >= 3:
        bgr_u8 = img_work[:, :, :3].astype(np.uint8)
        lab_img = cv2.cvtColor(bgr_u8, cv2.COLOR_BGR2LAB)
        lab_l = lab_img[:, :, 0]
        gray = cv2.cvtColor(bgr_u8, cv2.COLOR_BGR2GRAY)
        lum = np.minimum(lab_l, gray)
    else:
        lum = img_work.astype(np.uint8)

    # Guard against uniform synthetic canvases (e.g., blank 200-gray test backgrounds)
    dyn_range = float(np.ptp(lum))
    if dyn_range < 10.0:
        return default_res

    blurred = cv2.GaussianBlur(lum, (3, 3), 0)
    int_vals = blurred[mask_u8 > 0]
    far_bg_mask = dilated_outer == 0
    if excl_bool is not None:
        clean_far_bg = far_bg_mask & (~excl_bool)
        if int(np.count_nonzero(clean_far_bg)) >= 64:
            far_bg_mask = clean_far_bg
    far_ext_vals = blurred[far_bg_mask]
    int_med = float(np.median(int_vals)) if int_vals.size > 0 else 128.0
    ext_ref = float(np.percentile(far_ext_vals, 75)) if far_ext_vals.size >= 32 else float(np.percentile(blurred, 85))

    # If neither interior nor perimeter ring differs from background paper, no plant contrast is rendered
    ring_vals = blurred[ring_mask]
    if abs(int_med - ext_ref) < 10.0 and float(np.ptp(ring_vals)) < 10.0:
        return default_res

    if int_med > ext_ref + 25.0:
        # Inverted synthetic fixture (bright leaf foreground on dark background)
        otsu_thresh, tissue_bin = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        tissue_mask = (tissue_bin > 0) & (blurred > (ext_ref + 12.0))
    else:
        # Standard herbarium mount sheet (dark plant tissue on bright mount paper)
        otsu_thresh, tissue_bin = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        paper_floor = max(0.0, ext_ref - 15.0)
        int_p80 = float(np.percentile(int_vals, 80)) if int_vals.size > 0 else int_med
        contrast = max(0.0, ext_ref - int_med)
        lamina_ceil = min(float(otsu_thresh) + 3.0, paper_floor, int_p80 + 0.25 * contrast)
        tissue_mask = (tissue_bin > 0) & (blurred <= lamina_ceil)

        if lab_img is not None and np.any(tissue_mask & ring_mask):
            lab_f = lab_img.astype(np.float32)
            # 1. Neutral paper drop-shadow filter: cast shadows retain background paper (a*, b*) chromaticity
            paper_sample_mask = far_bg_mask & (blurred >= float(np.percentile(far_ext_vals, 55))) if far_ext_vals.size >= 32 else far_bg_mask
            if int(np.count_nonzero(paper_sample_mask)) >= 32:
                ext_a = float(np.median(lab_f[:, :, 1][paper_sample_mask]))
                ext_b = float(np.median(lab_f[:, :, 2][paper_sample_mask]))
                int_a = float(np.median(lab_f[:, :, 1][mask_u8 > 0]))
                int_b = float(np.median(lab_f[:, :, 2][mask_u8 > 0]))
                leaf_paper_chroma = math.hypot(int_a - ext_a, int_b - ext_b)
                if leaf_paper_chroma >= 6.5:
                    px_paper_chroma = np.hypot(lab_f[:, :, 1] - ext_a, lab_f[:, :, 2] - ext_b)
                    is_paper_shadow = (blurred > int_med) & (px_paper_chroma < (0.42 * leaf_paper_chroma))
                    tissue_mask = tissue_mask & (~is_paper_shadow)

            # 2. Same-organ interior margin LAB continuity match (distinguishes foreground leaf from
            #    photometrically distinct background rosette leaves, woody stems, or color targets)
            rim_rad = max(3, min(12, int(round(0.08 * min(bw_m, bh_m)))))
            k_rim = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * rim_rad + 1, 2 * rim_rad + 1))
            eroded_m = cv2.erode(mask_u8, k_rim, iterations=1)
            inner_rim = (mask_u8 > 0) & (eroded_m == 0)
            if int(np.count_nonzero(inner_rim)) < 24:
                inner_rim = mask_u8 > 0
            rim_pixels = lab_f[inner_rim]
            if rim_pixels.shape[0] >= 16:
                med_lab = np.median(rim_pixels, axis=0)
                rim_dists = np.linalg.norm(rim_pixels - med_lab, axis=1)
                rim_p90 = float(np.percentile(rim_dists, 90))
                max_lab_dist = max(24.0, 1.85 * rim_p90)
                ring_lab_dist = np.linalg.norm(lab_f - med_lab, axis=2)
                tissue_mask = tissue_mask & (ring_lab_dist <= max_lab_dist)

    if excl_bool is not None:
        tissue_mask = tissue_mask & (~excl_bool)

    tissue_ring_px = int(np.count_nonzero(ring_mask & tissue_mask))
    border_ratio = float(tissue_ring_px) / float(ring_count)
    is_severed = bool(border_ratio > float(max_tissue_border_ratio))

    straight_cut_flag = False
    straight_cut_ratio = 0.0
    if not is_severed and max_dim >= 140:
        straight_cut_flag, straight_cut_ratio = _has_artificial_straight_occlusion_cut(mask_u8)
        if straight_cut_flag:
            fold_chk = detect_folded_leaf(mask_u8)
            if not fold_chk.get("is_folded", False):
                is_severed = True
                border_ratio = max(border_ratio, straight_cut_ratio)

    return BoundaryContinuityResult({
        "is_severed_lobe": is_severed,
        "tissue_border_ratio": round(border_ratio, 6),
        "ring_pixel_count": ring_count,
        "tissue_ring_pixels": tissue_ring_px,
        "otsu_threshold": float(otsu_thresh),
        "straight_occlusion_cut": straight_cut_flag,
    })


def compute_classical_shape_descriptors(
    mask: np.ndarray,
    bbox: Optional[Tuple[int, int, int, int]] = None,
) -> Dict[str, float]:
    """
    Computes classical dimensionless shape descriptors (Krieger 2017; Russ 2011):
      - Solidity = Area / Convex Hull Area
      - Circularity (Form Factor) = 4 * pi * Area / Perimeter^2
      - Extent = Mask Area / Oriented Bounding Box Area (scaled by external PCD `bbox`
        inflation so 45-degree mounted leaves remain rotation-invariant while tiny
        fragments inside oversized PCD boxes remain penalized)
    """
    if mask is None or mask.size == 0:
        return {
            "area_px": 0.0,
            "perimeter_px": 0.0,
            "convex_hull_area_px": 0.0,
            "bbox_area_px": 0.0,
            "solidity": 0.0,
            "circularity": 0.0,
            "extent": 0.0,
        }

    mask_u8 = (mask > 0).astype(np.uint8) * 255
    area_px = float(np.count_nonzero(mask_u8))
    if area_px <= 0.0:
        return {
            "area_px": 0.0,
            "perimeter_px": 0.0,
            "convex_hull_area_px": 0.0,
            "bbox_area_px": 0.0,
            "solidity": 0.0,
            "circularity": 0.0,
            "extent": 0.0,
        }

    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    sub = mask_u8[by : by + bh, bx : bx + bw]
    contours, _ = cv2.findContours(sub, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return {
            "area_px": area_px,
            "perimeter_px": 0.0,
            "convex_hull_area_px": area_px,
            "bbox_area_px": float(max(bw * bh, 1)),
            "solidity": 1.0,
            "circularity": 0.0,
            "extent": min(1.0, area_px / float(max(bw * bh, 1))),
        }

    main_cnt = max(contours, key=cv2.contourArea)
    cnt_area = max(float(cv2.contourArea(main_cnt)), area_px)
    perim_px = float(cv2.arcLength(main_cnt, True))
    hull = cv2.convexHull(main_cnt)
    hull_area = float(cv2.contourArea(hull))

    solidity = min(max(area_px / hull_area, 0.0), 1.0) if hull_area > 0.0 else 0.0
    circularity = (
        min(max((4.0 * math.pi * cnt_area) / (perim_px * perim_px), 0.0), 1.0)
        if perim_px > 0.0
        else 0.0
    )

    tight_axis_area = float(max(bw * bh, 1))
    rect_w, rect_h = cv2.minAreaRect(main_cnt)[1]
    oriented_area = min(tight_axis_area, float(max(rect_w * rect_h, 1.0)))

    if bbox is not None and len(bbox) == 4:
        ymin, xmin, ymax, xmax = [int(v) for v in bbox]
        box_h = max(0, ymax - ymin)
        box_w = max(0, xmax - xmin)
        ext_bbox_area = float(box_h * box_w) if (box_h > 0 and box_w > 0) else tight_axis_area
        inflation = max(1.0, ext_bbox_area / tight_axis_area)
        bbox_area = oriented_area * inflation
    else:
        bbox_area = oriented_area

    extent = min(max(area_px / max(bbox_area, 1.0), 0.0), 1.0)

    return {
        "area_px": area_px,
        "perimeter_px": perim_px,
        "convex_hull_area_px": hull_area,
        "bbox_area_px": bbox_area,
        "solidity": solidity,
        "circularity": circularity,
        "extent": extent,
    }


def reflect_folded_mask(
    mask: np.ndarray,
    p1: Tuple[int, int],
    p2: Tuple[int, int],
    padding: Optional[int] = None
) -> np.ndarray:
    """
    Symmetrically reflects a folded hemi-blade mask across the straight chord line
    passing through p1 and p2, synthesizing a complete bilateral silhouette.
    Fully vectorized in OpenCV using an affine reflection matrix and bitwise_or.
    Applies morphological closing to seal the midrib seam.

    Parameters
    ----------
    mask : np.ndarray
        Binary uint8 mask of the folded leaf.
    p1 : Tuple[int, int]
        Starting endpoint of the fold chord (apex or base).
    p2 : Tuple[int, int]
        Ending endpoint of the fold chord (base or apex).
    padding : Optional[int]
        Canvas margin to avoid clipping reflected geometry.

    Returns
    -------
    np.ndarray
        Synthesized bilateral silhouette mask (uint8, 0 or 255).
    """
    h, w = mask.shape[:2]
    if h < 5 or w < 5:
        return mask.copy()

    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if bw == 0 or bh == 0:
        return mask.copy()

    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])
    dx = x2 - x1
    dy = y2 - y1
    chord_len = math.hypot(dx, dy)
    if chord_len < 1e-6:
        return mask.copy()

    theta = math.atan2(dy, dx)
    cos2 = math.cos(2.0 * theta)
    sin2 = math.sin(2.0 * theta)
    t_x = x1 * (1.0 - cos2) - y1 * sin2
    t_y = y1 * (1.0 + cos2) - x1 * sin2

    # Compute tight bounding box encompassing both original foreground and its reflection
    corners_x = np.array([bx, bx + bw, bx, bx + bw], dtype=np.float64)
    corners_y = np.array([by, by, by + bh, by + bh], dtype=np.float64)
    ref_x = cos2 * corners_x + sin2 * corners_y + t_x
    ref_y = sin2 * corners_x - cos2 * corners_y + t_y

    roi_x0 = int(math.floor(min(float(bx), float(np.min(ref_x))))) - 8
    roi_y0 = int(math.floor(min(float(by), float(np.min(ref_y))))) - 8
    roi_x1 = int(math.ceil(max(float(bx + bw), float(np.max(ref_x))))) + 8
    roi_y1 = int(math.ceil(max(float(by + bh), float(np.max(ref_y))))) + 8

    roi_h = max(1, roi_y1 - roi_y0)
    roi_w = max(1, roi_x1 - roi_x0)
    local_canvas = np.zeros((roi_h, roi_w), dtype=np.uint8)
    local_canvas[by - roi_y0 : by - roi_y0 + bh, bx - roi_x0 : bx - roi_x0 + bw] = (
        mask_u8[by : by + bh, bx : bx + bw]
    )

    lx1 = x1 - float(roi_x0)
    ly1 = y1 - float(roi_y0)

    # 2D Affine Reflection Matrix across the line passing through (lx1, ly1) at angle theta
    M = np.array([
        [cos2, sin2, lx1 * (1.0 - cos2) - ly1 * sin2],
        [sin2, -cos2, ly1 * (1.0 + cos2) - lx1 * sin2]
    ], dtype=np.float32)

    reflected = cv2.warpAffine(local_canvas, M, (roi_w, roi_h), flags=cv2.INTER_NEAREST)
    combined = cv2.bitwise_or(local_canvas, reflected)

    # Morphological closing to seal the midrib seam
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    sealed = cv2.morphologyEx(combined, cv2.MORPH_CLOSE, kernel)

    ly_indices, lx_indices = np.where(sealed > 0)
    if len(ly_indices) == 0:
        return mask.copy()

    gy_min = int(np.min(ly_indices)) + roi_y0
    gy_max = int(np.max(ly_indices)) + roi_y0
    gx_min = int(np.min(lx_indices)) + roi_x0
    gx_max = int(np.max(lx_indices)) + roi_x0

    # If the synthesized silhouette fits within original bounds, place back into original frame
    if gy_min >= 0 and gy_max < h and gx_min >= 0 and gx_max < w:
        out = np.zeros((h, w), dtype=np.uint8)
        sy0 = max(0, -roi_y0)
        sx0 = max(0, -roi_x0)
        dy0 = max(0, roi_y0)
        dx0 = max(0, roi_x0)
        ch_copy = min(roi_h - sy0, h - dy0)
        cw_copy = min(roi_w - sx0, w - dx0)
        if ch_copy > 0 and cw_copy > 0:
            out[dy0 : dy0 + ch_copy, dx0 : dx0 + cw_copy] = (
                sealed[sy0 : sy0 + ch_copy, sx0 : sx0 + cw_copy] > 127
            ).astype(np.uint8) * 255
        return out

    # Otherwise return tight bounding box with 5px padding clamped to padded canvas bounds
    pad = max(h, w) if padding is None else padding
    ph = h + 2 * pad
    pw = w + 2 * pad
    min_y_p = max(0, gy_min + pad - 5)
    max_y_p = min(ph, gy_max + pad + 6)
    min_x_p = max(0, gx_min + pad - 5)
    max_x_p = min(pw, gx_max + pad + 6)
    out_h = max(1, max_y_p - min_y_p)
    out_w = max(1, max_x_p - min_x_p)
    out_crop = np.zeros((out_h, out_w), dtype=np.uint8)
    p_roi_y0 = roi_y0 + pad
    p_roi_x0 = roi_x0 + pad
    sy0 = max(0, min_y_p - p_roi_y0)
    sx0 = max(0, min_x_p - p_roi_x0)
    dy0 = max(0, p_roi_y0 - min_y_p)
    dx0 = max(0, p_roi_x0 - min_x_p)
    ch_copy = min(roi_h - sy0, out_h - dy0)
    cw_copy = min(roi_w - sx0, out_w - dx0)
    if ch_copy > 0 and cw_copy > 0:
        out_crop[dy0 : dy0 + ch_copy, dx0 : dx0 + cw_copy] = (
            sealed[sy0 : sy0 + ch_copy, sx0 : sx0 + cw_copy] > 127
        ).astype(np.uint8) * 255
    return out_crop


class FoldDetectionResult(dict):
    """
    Subclass of dict holding fold detection results.
    Evaluates as a boolean according to `is_folded` for direct boolean comparisons/assertions,
    while remaining fully compatible with dictionary key indexing.
    """
    def __bool__(self) -> bool:
        return bool(self.get("is_folded", False))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, bool):
            return bool(self.get("is_folded", False)) == other
        return super().__eq__(other)


def detect_folded_leaf(
    mask: np.ndarray,
    contour: Optional[np.ndarray] = None,
    max_dev_ratio: float = 0.035,
    min_r2: float = 0.97,
    max_wl_ratio: float = 0.42,
    min_chord_length_px: float = 5.0,
) -> FoldDetectionResult:
    """
    Detects whether a leaf is folded in half along the midrib.

    Identifies the primary longitudinal axis connecting the basal petiole insertion
    to the blade apex, splits the perimeter into Side A and Side B lateral contours,
    and evaluates linearity and chord straightness.

    Detection Rules:
      - If exactly one lateral boundary is a straight chord (max deviation < 3.5% of
        chord length, R^2 > 0.97 to linear fit) AND W/L <= 0.45:
          Classified as is_folded = True, with the straight chord as the true midrib_axis.
          Symmetrically reflects the curved margin across the fold chord.
      - If folded irregularly or diagonally (neither side is collinear with petiole-to-apex
        axis), classified as is_irregular_fold = True.

    Parameters
    ----------
    mask : np.ndarray
        Binary leaf mask.
    contour : Optional[np.ndarray]
        External contour points. If None, extracted from mask.
    max_dev_ratio : float
        Maximum perpendicular distance to chord as fraction of chord length (default 0.035).
    min_r2 : float
        Minimum linear fit R^2 for straight boundary (default 0.97).
    max_wl_ratio : float
        Maximum width-to-length ratio to be classified as folded (default 0.45).

    Returns
    -------
    FoldDetectionResult
        Dictionary containing:
          - is_folded: bool
          - is_irregular_fold: bool
          - straight_side: Optional[str] ("Side A" or "Side B")
          - midrib_axis: Optional[Tuple[Tuple[int, int], Tuple[int, int]]]
          - wl_ratio: float
          - synthesized_mask: Optional[np.ndarray]
          - side_a_stats: Dict[str, Any]
          - side_b_stats: Dict[str, Any]
    """
    default_res: FoldDetectionResult = FoldDetectionResult({
        "is_folded": False,
        "is_irregular_fold": False,
        "straight_side": None,
        "midrib_axis": None,
        "wl_ratio": 0.0,
        "synthesized_mask": None,
        "side_a_stats": {"dev_ratio": 1.0, "r2": 0.0, "is_straight": False},
        "side_b_stats": {"dev_ratio": 1.0, "r2": 0.0, "is_straight": False},
    })

    if mask is None or np.count_nonzero(mask) < 20:
        return default_res

    if contour is None:
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            return default_res
        main_cnt = max(contours, key=cv2.contourArea)
    else:
        main_cnt = contour

    pts = main_cnt.reshape(-1, 2).astype(np.float32)
    n_pts = len(pts)
    if n_pts < 10:
        return default_res

    # 1. Width-to-length ratio via Oriented Bounding Box
    rect = cv2.minAreaRect(pts)
    rect_w, rect_h = rect[1]
    box_len = max(rect_w, rect_h)
    box_wid = min(rect_w, rect_h)
    wl_ratio = float(box_wid / box_len) if box_len > 1e-6 else 0.0
    default_res["wl_ratio"] = round(wl_ratio, 4)

    # 2. Longitudinal Midrib Axis via cv2.fitLine
    line_params = cv2.fitLine(pts, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]
    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy

    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))
    p_apex = (int(pts[min_idx, 0]), int(pts[min_idx, 1]))
    p_base = (int(pts[max_idx, 0]), int(pts[max_idx, 1]))

    chord_dx = float(p_base[0] - p_apex[0])
    chord_dy = float(p_base[1] - p_apex[1])
    chord_len = math.hypot(chord_dx, chord_dy)
    if chord_len < min_chord_length_px or min_idx == max_idx:
        return default_res

    # 3. Split boundary perimeter into Side A and Side B
    idx_1 = min(min_idx, max_idx)
    idx_2 = max(min_idx, max_idx)
    side_a_pts = pts[idx_1 : idx_2 + 1]
    side_b_pts = np.vstack([pts[idx_2:], pts[: idx_1 + 1]])

    def _calc_side_metrics(side_pts: np.ndarray) -> Tuple[float, float, bool]:
        if len(side_pts) < 3:
            return 1.0, 0.0, False
        p1 = side_pts[0]
        p2 = side_pts[-1]
        L = float(np.hypot(p2[0] - p1[0], p2[1] - p1[1]))
        if L < 1e-6:
            return 1.0, 0.0, False

        dx = float(p2[0] - p1[0])
        dy = float(p2[1] - p1[1])
        # Vectorized perpendicular distance to chord line
        cross = np.abs(dy * side_pts[:, 0] - dx * side_pts[:, 1] + p2[0] * p1[1] - p2[1] * p1[0])
        dists = cross / L
        max_dev = float(np.max(dists))
        dev_ratio = max_dev / L

        # R^2 linear fit along chord
        centroid = np.mean(side_pts, axis=0)
        ss_tot = float(np.sum(np.sum((side_pts - centroid) ** 2, axis=1)))
        ss_res = float(np.sum(dists ** 2))
        r2 = max(0.0, min(1.0, 1.0 - (ss_res / ss_tot))) if ss_tot > 1e-6 else 1.0
        is_straight = (dev_ratio < max_dev_ratio) and (r2 > min_r2)
        return dev_ratio, r2, is_straight

    dev_a, r2_a, straight_a = _calc_side_metrics(side_a_pts)
    dev_b, r2_b, straight_b = _calc_side_metrics(side_b_pts)

    default_res["side_a_stats"] = {"dev_ratio": round(dev_a, 4), "r2": round(r2_a, 4), "is_straight": straight_a}
    default_res["side_b_stats"] = {"dev_ratio": round(dev_b, 4), "r2": round(r2_b, 4), "is_straight": straight_b}

    # 4. Detection Rule Evaluation
    if straight_a ^ straight_b:
        # Exactly one lateral boundary is a straight chord
        if wl_ratio <= max_wl_ratio:
            default_res["is_folded"] = True
            default_res["straight_side"] = "Side A" if straight_a else "Side B"
            default_res["midrib_axis"] = (p_apex, p_base)
            default_res["synthesized_mask"] = reflect_folded_mask(mask, p_apex, p_base)
            return default_res

    # Check for irregular or diagonal fold:
    # If aspect ratio is narrow (W/L <= 0.45) but neither side matches the straight midrib chord
    if wl_ratio <= max_wl_ratio and (not straight_a and not straight_b):
        # Examine if an off-axis straight edge cuts diagonally across the contour
        hull = cv2.convexHull(main_cnt)
        area_hull = float(cv2.contourArea(hull))
        area_mask = float(np.count_nonzero(mask))
        solidity = area_mask / area_hull if area_hull > 0 else 0.0

        approx = cv2.approxPolyDP(pts.reshape(-1, 1, 2), 0.02 * cv2.arcLength(pts.reshape(-1, 1, 2), True), True).reshape(-1, 2)
        has_diagonal_crease = False
        is_sheared_quad_crease = False
        ang_midrib = math.atan2(chord_dy, chord_dx)

        for i in range(len(approx)):
            p_s = approx[i]
            p_e = approx[(i + 1) % len(approx)]
            seg_len = math.hypot(p_e[0] - p_s[0], p_e[1] - p_s[1])
            ang_seg = math.atan2(p_e[1] - p_s[1], p_e[0] - p_s[0])
            diff_ang = abs((ang_seg - ang_midrib + math.pi) % (2.0 * math.pi) - math.pi)
            diff_deg = math.degrees(diff_ang)
            diff_deg = min(diff_deg, 180.0 - diff_deg)
            if seg_len > 0.55 * chord_len and 20.0 < diff_deg < 70.0:
                has_diagonal_crease = True
                break
            if len(approx) == 4 and seg_len > 0.45 * chord_len and 12.0 < diff_deg < 70.0:
                is_sheared_quad_crease = True

        if (has_diagonal_crease and solidity < 0.72) or is_sheared_quad_crease:
            default_res["is_irregular_fold"] = True

    return default_res


def is_botanical_dissection(
    contour_or_mask: Optional[np.ndarray] = None,
    /,
    solidity: Optional[float] = None,
    p_apex: Optional[Tuple[int, int]] = None,
    p_base: Optional[Tuple[int, int]] = None,
    blade_width: Optional[float] = None,
    min_solidity: float = 0.18,
    max_solidity: float = 0.72,
    defect_depth_ratio: float = 0.08,
    min_defect_count: int = 3,
    check_solidity_range: bool = True,
    *,
    contour: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    bipinnatisect_min_solidity: Optional[float] = None,
) -> bool:
    """
    Evaluates whether a leaf contour with lower solidity (0.18 <= Solidity < 0.72)
    represents authentic botanical dissection (bipinnatisect, pinnatifid, or lyrate lobes
    with bilateral sinuses) rather than asymmetric foreign leaf occlusion or damage.

    Detection Rules:
      - Filter out minor crenation noise: retain only significant defects where
        defect depth exceeds 8% of maximum blade width.
      - Evaluate Defect Topology:
        * Foreign leaf occlusion: 1 or 2 large, asymmetric, localized deficits on one side.
        * Finely dissected / bipinnatisect (0.18 <= solidity < 0.48): >= 4 deep bilateral
          convexity defects distributed on both sides of the longitudinal midline.
        * Lyrate / pinnatifid (0.48 <= solidity < 0.72): >= 3 alternating, bilaterally
          distributed sinuses along the margin.
    """
    active_input = contour if contour is not None else (contour_or_mask if contour_or_mask is not None else mask)
    if active_input is None or active_input.size == 0:
        return False
    contour = active_input

    if contour.ndim == 2 and contour.shape[1] > 2:
        cnts, _ = cv2.findContours(
            (contour > 0).astype(np.uint8) * 255,
            cv2.RETR_EXTERNAL,
            cv2.CHAIN_APPROX_SIMPLE,
        )
        if not cnts:
            return False
        cnt = max(cnts, key=cv2.contourArea).reshape(-1, 2).astype(np.int32)
    else:
        cnt = contour.reshape(-1, 2).astype(np.int32)

    if len(cnt) < 4:
        return False

    cnt_f = cnt.astype(np.float32)
    if solidity is None:
        area = float(cv2.contourArea(cnt_f))
        hull = cv2.convexHull(cnt_f)
        hull_area = float(cv2.contourArea(hull))
        solidity = float(area / hull_area) if hull_area > 0 else 0.0

    eff_min_solidity = (
        float(bipinnatisect_min_solidity)
        if bipinnatisect_min_solidity is not None
        else float(min_solidity)
    )
    if check_solidity_range:
        if not (eff_min_solidity <= solidity < max_solidity):
            return False

    rect = cv2.minAreaRect(cnt_f)
    minor_w = float(min(rect[1]))
    major_l = float(max(rect[1]))

    # Reject multi-leaf overlapping clumps masquerading as dissected leaves on high-DPI scans
    is_high_dpi = (
        (contour.ndim == 2 and max(contour.shape[:2]) > 500)
        or cnt_f[:, 0].max() > 500
        or cnt_f[:, 1].max() > 500
        or (solidity >= 0.48 and major_l > 240.0)
        or minor_w > 360.0
    )
    if is_high_dpi and major_l > 220.0 and (minor_w > 360.0 or (solidity >= 0.48 and (minor_w / max(major_l, 1.0)) > 0.78)):
        return False

    if blade_width is None or blade_width <= 0:
        blade_width = minor_w
    if blade_width <= 0:
        return False

    min_defect_depth = defect_depth_ratio * blade_width

    cnt_cv = cnt.reshape(-1, 1, 2)
    hull_indices = cv2.convexHull(cnt_cv, returnPoints=False)
    if hull_indices is None or len(hull_indices) < 3:
        return False

    try:
        defects = cv2.convexityDefects(cnt_cv, hull_indices)
    except cv2.error:
        try:
            hull_sorted = np.unique(hull_indices).reshape(-1, 1)
            if len(hull_sorted) < 3:
                return False
            defects = cv2.convexityDefects(cnt_cv, hull_sorted)
        except cv2.error:
            return False
    if defects is None or len(defects) == 0:
        return False

    # Extract defect depths (OpenCV stores fixpt_depth as distance * 256.0)
    def_arr = defects.reshape(-1, 4)
    depths = def_arr[:, 3].astype(np.float64) / 256.0
    farthest_indices = def_arr[:, 2].astype(np.int32)

    # Filter significant defects exceeding depth threshold
    is_native_scan = (
        (contour.ndim == 2 and max(contour.shape[:2]) > 500)
        or cnt_f[:, 0].max() > 500
        or cnt_f[:, 1].max() > 500
        or major_l > 350.0
    )
    sig_mask = depths >= min_defect_depth
    if solidity < 0.48:
        eff_min_defect_count = max(min_defect_count, 4)
    else:
        eff_min_defect_count = (
            min(min_defect_count, 2)
            if (is_native_scan and solidity >= 0.48 and np.any(sig_mask) and float(np.max(depths[sig_mask])) >= 0.11 * blade_width)
            else min_defect_count
        )
    if np.sum(sig_mask) < eff_min_defect_count:
        return False

    sig_farthest_pts = cnt[farthest_indices[sig_mask]]

    # Longitudinal midline vector from apex to base
    if p_apex is None or p_base is None or solidity < 0.48:
        line_params = cv2.fitLine(cnt_f, cv2.DIST_L2, 0, 0.01, 0.01)
        vx, vy, x0, y0 = [float(v[0]) for v in line_params]
        projections = (cnt_f[:, 0] - x0) * vx + (cnt_f[:, 1] - y0) * vy
        t_min = float(np.min(projections))
        t_max = float(np.max(projections))
        p_apex = (x0 + t_min * vx, y0 + t_min * vy)
        p_base = (x0 + t_max * vx, y0 + t_max * vy)

    dx = float(p_base[0] - p_apex[0])
    dy = float(p_base[1] - p_apex[1])
    line_len = math.hypot(dx, dy)
    if line_len < 1e-6:
        return False

    # Signed cross product to determine left vs right side of midline:
    # cross = dx * (y - y_apex) - dy * (x - x_apex)
    cross_vals = dx * (sig_farthest_pts[:, 1] - p_apex[1]) - dy * (sig_farthest_pts[:, 0] - p_apex[0])
    tol = 0.01 * line_len * blade_width

    side_a_count = int(np.sum(cross_vals > tol))
    side_b_count = int(np.sum(cross_vals < -tol))
    total_sig = side_a_count + side_b_count

    # Botanical dissection: >= eff_min_defect_count defects distributed on both sides (side_a >= 1 and side_b >= 1)
    is_dissected = (total_sig >= eff_min_defect_count) and (side_a_count >= 1) and (side_b_count >= 1)
    return bool(is_dissected)


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
    rect_b = cv2.minAreaRect(pts)
    blade_width_est = float(min(rect_b[1])) if min(rect_b[1]) > 0 else float(min(orig_h, orig_w))
    if blade_width_est > 250.0 and min(orig_h, orig_w) > 500:
        if max_petiole_width == 16.0:
            max_petiole_width = max(16.0, 0.105 * blade_width_est)
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
    narrow_indices = []
    for k in range(len(s_steps)):
        if widths[k] <= max_petiole_width:
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
            if ar > aspect_ratio_protrusion_thresh and w_5 <= max_petiole_width and stalk_len >= 6.0:
                junction_idx = narrow_indices[-1]
                reason = "FALLBACK_MORPHOLOGICAL_PROTRUSION_TRIMMED"
            else:
                reason = "CLEAN_BLADE_UNCLIPPED"
        else:
            reason = "CLEAN_BLADE_UNCLIPPED"
    else:
        # A narrow petiole stalk exists; search for width inflection junction into lamina
        for k in narrow_indices:
            if s_steps[k] < min_petiole_length_px:
                continue
            span = min(3, len(s_steps) - 1 - k)
            if span > 0:
                dw_dl = (widths[k + span] - widths[k]) / (s_steps[k + span] - s_steps[k])
                if dw_dl >= gradient_threshold:
                    future_w = widths[k: min(k + 8, len(widths))]
                    if np.max(future_w) >= max_petiole_width * 1.3:
                        j_idx = max(k, narrow_indices[-1])
                        if (
                            j_idx + 2 < len(widths)
                            and widths[j_idx + 1] <= max_petiole_width * 1.25
                            and widths[j_idx + 2] >= max_petiole_width * 1.5
                        ):
                            j_idx += 2
                        junction_idx = j_idx
                        reason = "WIDTH_INFLECTION_JUNCTION_TRIMMED"
                        break
        if junction_idx is None:
            last_narrow = narrow_indices[-1]
            if last_narrow < len(widths) - 2 and widths[last_narrow + 1] > max_petiole_width:
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

    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(clipped, connectivity=8)
    if num_labels > 1:
        largest_label = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
        clipped = np.where(labels == largest_label, 255, 0).astype(np.uint8)
    cnts_ext, _ = cv2.findContours(clipped, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if cnts_ext:
        clipped_filled = np.zeros_like(clipped)
        cv2.drawContours(clipped_filled, [max(cnts_ext, key=cv2.contourArea)], -1, 255, -1)
        clipped = clipped_filled

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



TIER_1_DIRECT = "Tier_1_Direct"
TIER_2_REFLECTED = "Tier_2_Reflected"

CANONICAL_TIER_MAP: Dict[str, str] = {
    "tier1": TIER_1_DIRECT,
    "tier_1": TIER_1_DIRECT,
    "tier 1": TIER_1_DIRECT,
    "tier 1 (dissected)": TIER_1_DIRECT,
    "tier 1 (bipinnatisect)": TIER_1_DIRECT,
    "tier_1_direct": TIER_1_DIRECT,
    "Tier_1_Direct": TIER_1_DIRECT,
    "tier2": TIER_2_REFLECTED,
    "tier_2": TIER_2_REFLECTED,
    "tier 2": TIER_2_REFLECTED,
    "tier_2_reflected": TIER_2_REFLECTED,
    "Tier_2_Reflected": TIER_2_REFLECTED,
}


def normalize_assigned_tier(tier_value: str) -> str:
    """
    Normalizes routing tier values ('tier1', 'tier2', 'Tier 1', 'Tier 2')
    to the canonical cross-language CSV contract names ('Tier_1_Direct', 'Tier_2_Reflected').
    """
    if not tier_value:
        return "rejected"
    clean = str(tier_value).strip()
    return CANONICAL_TIER_MAP.get(clean, CANONICAL_TIER_MAP.get(clean.lower(), clean))


def compute_sheet_adaptive_min_blade_area_cm2(
    candidate_areas_cm2: Any,
    config_min_blade_area_cm2: float = 0.75,
    absolute_floor_cm2: float = 0.15,
) -> float:
    """
    Computes a cohort/sheet-adaptive physical minimum blade area threshold (cm^2)
    for calibrated herbarium sheets so mature basal leaves of dwarf/alpine Packera
    taxa (e.g., ~0.25-0.45 cm^2) pass while tiny noise fragments on standard
    macro-leaved sheets remain rejected:
        effective_min_area_cm2 = min(config_min_blade_area_cm2, max(0.15, 0.25 * median_sheet_leaf_area_cm2))
    Falls back to `config_min_blade_area_cm2` when `candidate_areas_cm2` is empty
    or contains no positive finite areas.
    """
    cfg_floor = float(config_min_blade_area_cm2)
    abs_floor = float(absolute_floor_cm2)
    if candidate_areas_cm2 is None:
        return cfg_floor
    valid_areas: List[float] = []
    for val in candidate_areas_cm2:
        if val is None:
            continue
        try:
            fval = float(val)
        except (TypeError, ValueError):
            continue
        if np.isfinite(fval) and fval > 0.0:
            valid_areas.append(fval)
    if not valid_areas:
        return cfg_floor
    median_sheet_leaf_area_cm2 = float(np.median(valid_areas))
    return float(min(cfg_floor, max(abs_floor, 0.25 * median_sheet_leaf_area_cm2)))


def is_cauline_or_stem_contaminant(
    *,
    cy: float,
    area: float,
    aspect_ratio: float,
    solidity: float,
    basal_cy_med: float,
    med_area: float,
    img_h: float,
    capitula_count: int = 0,
    cap_cy: Optional[float] = None,
    has_petiole: bool = True,
    med_rosette_radius: float = 0.0,
    radial_dist_from_rosette: Optional[float] = None,
) -> Tuple[bool, Optional[str]]:
    """
    Heterophylly / cauline-bract / stem-segment filter that operates both when
    reproductive capitula are detected (`capitula_count > 0`) and on vegetative
    or non-flowering sheets (`capitula_count == 0`).

    Rejects upper-stem cauline leaves and stem segments before routing to
    `Tier_1_Direct` or `Tier_2_Reflected` if:
      (a) `capitula_count > 0` and the leaf centroid lies in the upper corymb/cauline
          zone (`cy < cap_cy + 0.30 * (basal_cy_med - cap_cy)`) with `area < 0.50 * med_area`, OR
      (b) the candidate is a narrow, sessile/clasping upper cauline bract or stem segment
          (`aspect_ratio > 3.8` and `solidity < 0.78` and `cy < basal_cy_med - 0.18 * img_h`), OR
      (c) the candidate lacks a petiole attachment and is displaced `> 2.2 * med_rosette_radius`
          above the basal rosette centroid with `area < 0.35 * med_area`.
    """
    if (
        capitula_count > 0
        and cap_cy is not None
        and not np.isnan(cap_cy)
        and basal_cy_med > cap_cy
        and cy < (cap_cy + 0.30 * (basal_cy_med - cap_cy))
        and area < (0.50 * med_area)
    ):
        return True, "CAULINE_CORYMB_ZONE"

    if (
        aspect_ratio > 3.8
        and solidity < 0.78
        and cy < (basal_cy_med - 0.18 * float(img_h))
    ):
        return True, "CAULINE_SESSILE_BRACT_OR_STEM"

    upward_disp = float(basal_cy_med - cy)
    eff_disp = (
        max(upward_disp, float(radial_dist_from_rosette))
        if (radial_dist_from_rosette is not None and upward_disp > 0)
        else upward_disp
    )
    if (
        (not has_petiole)
        and med_rosette_radius > 0.0
        and upward_disp > 0.0
        and eff_disp > (2.2 * med_rosette_radius)
        and area < (0.35 * med_area)
    ):
        return True, "DISPLACED_SESSILE_CAULINE"

    return False, None


@dataclass
class LeafRoutingResult:
    """Structured diagnostic and extraction result from geometric quality gatekeeping."""
    tier: str                                       # "Tier 1", "Tier 2", or "Failed QC"
    assigned_tier: str                              # "tier1", "tier2", or "rejected"
    solidity: float
    ucs: float
    is_folded: bool
    is_irregular_fold: bool
    is_dissected: bool
    processed_mask: np.ndarray
    contour_coords: Optional[np.ndarray] = None
    rejection_reason: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def canonical_tier(self) -> str:
        """Returns the canonical cross-language CSV tier ('Tier_1_Direct', 'Tier_2_Reflected', or 'rejected')."""
        return normalize_assigned_tier(self.assigned_tier)

    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        if key in self.metadata:
            return self.metadata[key]
        raise KeyError(key)

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default


def classify_and_route_leaf(
    mask: np.ndarray,
    config: Optional[Any] = None,
    min_ucs: float = 0.85,
    min_solidity: float = 0.72,
    dissection_min_solidity: float = 0.48,
    extract_contours: bool = True,
    num_contour_points: int = 120,
    taxon: Optional[str] = None,
    taxa_with_lyrate_tendency: Optional[List[str]] = None,
    max_chord_deviation_ratio: Optional[float] = None,
    max_folded_aspect_ratio: Optional[float] = None,
    sinus_defect_min_depth_ratio: Optional[float] = None,
    min_bilateral_sinus_count: Optional[int] = None,
    min_chord_length_px: Optional[float] = None,
    fold_detection_enabled: Optional[bool] = None,
    dissection_enabled: Optional[bool] = None,
    rosette_centroid: Optional[Tuple[float, float]] = None,
    already_trimmed: bool = False,
    trim_meta: Optional[Dict[str, Any]] = None,
    sheet_shape: Optional[Tuple[int, int]] = None,
    crop_img: Optional[np.ndarray] = None,
    boundary_dilation_px: int = 8,
    max_tissue_border_ratio: float = 0.15,
    exclude_mask: Optional[np.ndarray] = None,
) -> LeafRoutingResult:
    """
    Comprehensive 100% label-blind geometric quality gatekeeper and routing engine:
      1. Evaluates midrib pose and initial solidity/UCS metrics.
      2. Straight-chord fold detection:
         - If folded along midrib (is_folded = True), symmetrically reflects hemi-blade
           into synthesized bilateral silhouette and routes to Tier 2.
         - If irregular/diagonal fold (is_irregular_fold = True), routes to Failed QC.
      3. Botanical dissection inspection (strictly intrinsic bilateral sinus geometry):
         - If 0.48 <= solidity < min_solidity and bilateral convexity defects indicate
           pinnatifid/lyrate lobing via is_botanical_dissection(), accepts into Tier 1 (Dissected).
         - If 0.18 <= solidity < 0.48 with >= 4 deep bilateral convexity-hull sinuses,
           accepts into Tier 1 (Bipinnatisect).
         - Note: `taxon` and `taxa_with_lyrate_tendency` are retained strictly as deprecated
           no-ops for backward compatibility and never influence routing.
      4. Standard pristine gatekeeping:
         - Leaves with UCS >= min_ucs and Solidity >= min_solidity route to Tier 1.
      5. Fallback hemi-blade reflection (Tier 2):
         - Partially occluded leaves attempt 2-path reflection or route to Failed QC.
    """
    _ = (taxon, taxa_with_lyrate_tendency)  # Deprecated no-ops; routing is 100% label-blind
    if config is not None:
        thresh = getattr(config, "thresholds", config)
        if hasattr(thresh, "min_ucs"):
            min_ucs = getattr(thresh, "min_ucs", min_ucs)
        elif isinstance(thresh, dict):
            min_ucs = thresh.get("min_ucs", min_ucs)

        if hasattr(thresh, "min_solidity"):
            min_solidity = getattr(thresh, "min_solidity", min_solidity)
        elif isinstance(thresh, dict) and "min_solidity" in thresh:
            min_solidity = thresh.get("min_solidity", min_solidity)

        # Dissection configuration (strictly geometric; no species-name overrides)
        diss_cfg = getattr(thresh, "dissection", None) or getattr(thresh, "solidity", None)
        if diss_cfg is not None and not isinstance(diss_cfg, dict):
            if dissection_enabled is None:
                dissection_enabled = getattr(diss_cfg, "enabled", True)
            dissection_min_solidity = getattr(
                diss_cfg, "min_solidity_dissected", getattr(diss_cfg, "min_dissected", dissection_min_solidity)
            )
            if sinus_defect_min_depth_ratio is None:
                sinus_defect_min_depth_ratio = getattr(diss_cfg, "sinus_defect_min_depth_ratio", 0.08)
            if min_bilateral_sinus_count is None:
                min_bilateral_sinus_count = getattr(diss_cfg, "min_bilateral_sinus_count", 3)
        elif isinstance(thresh, dict):
            diss_dict = thresh.get("dissection") or thresh.get("solidity")
            if isinstance(diss_dict, dict):
                if dissection_enabled is None:
                    dissection_enabled = diss_dict.get("enabled", True)
                dissection_min_solidity = diss_dict.get(
                    "min_solidity_dissected", diss_dict.get("min_dissected", dissection_min_solidity)
                )
                if sinus_defect_min_depth_ratio is None:
                    sinus_defect_min_depth_ratio = diss_dict.get("sinus_defect_min_depth_ratio", 0.08)
                if min_bilateral_sinus_count is None:
                    min_bilateral_sinus_count = diss_dict.get("min_bilateral_sinus_count", 3)

        # Fold detection configuration
        fold_cfg = getattr(thresh, "fold_detection", None)
        if fold_cfg is not None and not isinstance(fold_cfg, dict):
            if fold_detection_enabled is None:
                fold_detection_enabled = getattr(fold_cfg, "enabled", True)
            if max_chord_deviation_ratio is None:
                max_chord_deviation_ratio = getattr(fold_cfg, "max_chord_deviation_ratio", None)
            if max_folded_aspect_ratio is None:
                max_folded_aspect_ratio = getattr(fold_cfg, "max_folded_aspect_ratio", None)
            if min_chord_length_px is None:
                min_chord_length_px = getattr(fold_cfg, "min_chord_length_px", None)
        elif isinstance(thresh, dict) and isinstance(thresh.get("fold_detection"), dict):
            fold_dict = thresh["fold_detection"]
            if fold_detection_enabled is None:
                fold_detection_enabled = fold_dict.get("enabled", True)
            if max_chord_deviation_ratio is None:
                max_chord_deviation_ratio = fold_dict.get("max_chord_deviation_ratio")
            if max_folded_aspect_ratio is None:
                max_folded_aspect_ratio = fold_dict.get("max_folded_aspect_ratio")
            if min_chord_length_px is None:
                min_chord_length_px = fold_dict.get("min_chord_length_px")

    if max_chord_deviation_ratio is None:
        max_chord_deviation_ratio = 0.035
    if max_folded_aspect_ratio is None:
        max_folded_aspect_ratio = 0.42
    if min_chord_length_px is None:
        min_chord_length_px = 150.0
    if fold_detection_enabled is None:
        fold_detection_enabled = True
    if dissection_enabled is None:
        dissection_enabled = True
    if sinus_defect_min_depth_ratio is None:
        sinus_defect_min_depth_ratio = 0.08
    if min_bilateral_sinus_count is None:
        min_bilateral_sinus_count = 3

    if mask is None or mask.size == 0 or np.count_nonzero(mask) < 20:
        return LeafRoutingResult(
            tier="Failed QC",
            assigned_tier="rejected",
            solidity=0.0,
            ucs=0.0,
            is_folded=False,
            is_irregular_fold=False,
            is_dissected=False,
            processed_mask=np.zeros_like(mask) if mask is not None else np.zeros((10, 10), dtype=np.uint8),
            rejection_reason="mask_empty_or_too_small"
        )

    # Compute initial geometric metrics and disambiguated (p_apex, p_base), then trim petiole tail
    ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(
        mask, rosette_centroid=rosette_centroid
    )
    if not already_trimmed or rosette_centroid is not None:
        trimmed_mask, trim_meta = trim_petiole_tail(
            mask,
            p_base=p_base,
            p_apex=p_apex,
            return_metadata=True,
            sheet_shape=sheet_shape,
        )
        if trim_meta.get("is_trimmed", False):
            mask = trimmed_mask
            ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(
                mask, rosette_centroid=rosette_centroid
            )
    elif trim_meta is None:
        trim_meta = {"is_trimmed": False, "reason": "PRE_TRIMMED"}
    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    bx_m, by_m, bw_m, bh_m = cv2.boundingRect(mask_u8)
    if bw_m == 0 or bh_m == 0:
        return LeafRoutingResult(
            tier="Failed QC",
            assigned_tier="rejected",
            solidity=solidity,
            ucs=ucs,
            is_folded=False,
            is_irregular_fold=False,
            is_dissected=False,
            processed_mask=mask,
            rejection_reason="no_contour_found"
        )
    x0_c = max(0, bx_m - 2)
    y0_c = max(0, by_m - 2)
    x1_c = min(mask_u8.shape[1], bx_m + bw_m + 2)
    y1_c = min(mask_u8.shape[0], by_m + bh_m + 2)
    contours, _ = cv2.findContours(
        mask_u8[y0_c:y1_c, x0_c:x1_c],
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_NONE,
        offset=(x0_c, y0_c),
    )
    if not contours:
        return LeafRoutingResult(
            tier="Failed QC",
            assigned_tier="rejected",
            solidity=solidity,
            ucs=ucs,
            is_folded=False,
            is_irregular_fold=False,
            is_dissected=False,
            processed_mask=mask,
            rejection_reason="no_contour_found"
        )

    main_cnt = max(contours, key=cv2.contourArea)

    # 1. Straight-chord fold detection
    if fold_detection_enabled:
        fold_info = detect_folded_leaf(
            mask,
            main_cnt,
            max_dev_ratio=max_chord_deviation_ratio,
            max_wl_ratio=max_folded_aspect_ratio,
            min_chord_length_px=min_chord_length_px,
        )
    else:
        fold_info = {
            "is_folded": False,
            "is_irregular_fold": False,
            "straight_side": None,
            "midrib_axis": None,
            "wl_ratio": 0.0,
            "synthesized_mask": None,
            "side_a_stats": {"dev_ratio": 1.0, "r2": 0.0, "is_straight": False},
            "side_b_stats": {"dev_ratio": 1.0, "r2": 0.0, "is_straight": False},
        }
    is_folded = fold_info["is_folded"]
    is_irregular_fold = fold_info["is_irregular_fold"]

    meta = {
        "midrib_angle_deg": angle_deg,
        "p_apex": p_apex,
        "p_base": p_base,
        "fold_info": fold_info,
        "reflection_applied": False,
        "petiole_trim_info": trim_meta,
        "is_severed_lobe": False,
        "tissue_border_ratio": 0.0,
    }

    if crop_img is not None:
        bc_res = check_boundary_tissue_continuity(
            crop_img,
            mask_u8,
            dilation_px=boundary_dilation_px,
            max_tissue_border_ratio=max_tissue_border_ratio,
            exclude_mask=exclude_mask,
            petiole_base_pt=p_base,
        )
        meta["boundary_continuity"] = bc_res
        is_bc_severed = bool(bc_res["is_severed_lobe"])
        if is_bc_severed and bc_res.get("straight_occlusion_cut", False) and (is_folded or is_irregular_fold):
            is_bc_severed = False
        meta["is_severed_lobe"] = is_bc_severed
        meta["tissue_border_ratio"] = float(bc_res["tissue_border_ratio"])
        if is_bc_severed:
            return LeafRoutingResult(
                tier="Failed QC",
                assigned_tier="rejected",
                solidity=solidity,
                ucs=ucs,
                is_folded=False,
                is_irregular_fold=False,
                is_dissected=False,
                processed_mask=mask,
                rejection_reason="SEVERED_LOBE",
                metadata=meta,
            )

    if is_folded:
        # Route folded leaf directly to Tier 2 with synthesized bilateral silhouette
        synth_mask = fold_info["synthesized_mask"]
        if synth_mask is not None and np.count_nonzero(synth_mask) > 0:
            coords = (
                resample_normalized_contour(synth_mask, num_contour_points, petiole_attachment_pt=p_base)
                if extract_contours
                else None
            )
            synth_ucs, synth_sol, _, _, _ = compute_geometric_metrics(
                synth_mask, rosette_centroid=rosette_centroid
            )
            meta["reflection_applied"] = True
            return LeafRoutingResult(
                tier="Tier 2",
                assigned_tier="tier2",
                solidity=synth_sol,
                ucs=synth_ucs,
                is_folded=True,
                is_irregular_fold=False,
                is_dissected=False,
                processed_mask=synth_mask,
                contour_coords=coords,
                metadata=meta
            )

    if is_irregular_fold:
        # Check if finely dissected / bipinnatisect leaf (0.18 <= solidity < 0.48) triggered crease heuristic
        if dissection_enabled and (0.18 <= solidity < 0.48) and ucs >= 0.55 and is_botanical_dissection(
            main_cnt,
            solidity=solidity,
            p_apex=p_apex,
            p_base=p_base,
            min_solidity=0.18,
            max_solidity=min_solidity,
            defect_depth_ratio=sinus_defect_min_depth_ratio,
            min_defect_count=max(min_bilateral_sinus_count, 4),
        ):
            is_irregular_fold = False
        else:
            return LeafRoutingResult(
                tier="Failed QC",
                assigned_tier="rejected",
                solidity=solidity,
                ucs=ucs,
                is_folded=False,
                is_irregular_fold=True,
                is_dissected=False,
                processed_mask=mask,
                rejection_reason="IRREGULAR_FOLD",
                metadata=meta
            )

    # 2. Standard Tier 1 Pristine Check
    if ucs >= min_ucs and solidity >= min_solidity:
        coords = resample_normalized_contour(mask, num_contour_points, petiole_attachment_pt=p_base) if extract_contours else None
        meta["reflection_applied"] = False
        return LeafRoutingResult(
            tier="Tier 1",
            assigned_tier="tier1",
            solidity=solidity,
            ucs=ucs,
            is_folded=False,
            is_irregular_fold=False,
            is_dissected=False,
            processed_mask=mask,
            contour_coords=coords,
            metadata=meta
        )

    # 3. Botanical Dissection Inspection (0.18 <= solidity < min_solidity)
    # Strictly label-blind: relies solely on intrinsic bilateral sinus geometry via is_botanical_dissection().
    ref_shape = sheet_shape if sheet_shape is not None else mask.shape[:2]
    is_high_res_mask = max(ref_shape) > 500
    is_bipinnatisect_window = 0.18 <= solidity < 0.48
    eff_diss_min_sol = min(dissection_min_solidity, 0.18)
    hull_cnt = cv2.convexHull(main_cnt)
    rect_cnt = cv2.minAreaRect(main_cnt)
    hull_ucs = min(
        1.0,
        float(cv2.contourArea(hull_cnt)) / max(0.70 * float(rect_cnt[1][0] * rect_cnt[1][1]), 1.0),
    )
    eff_diss_ucs = (
        max(ucs, hull_ucs)
        if (is_bipinnatisect_window or not is_high_res_mask or (ucs >= 0.52 and hull_ucs >= 0.80))
        else ucs
    )
    meta["hull_ucs"] = round(float(hull_ucs), 4)
    allow_dissection = True
    if is_bipinnatisect_window:
        if eff_diss_ucs < 0.55:
            allow_dissection = False
    elif is_high_res_mask and (eff_diss_ucs < 0.80 or ucs < 0.52 or solidity < max(dissection_min_solidity, 0.48)):
        allow_dissection = False

    if dissection_enabled and (eff_diss_min_sol <= solidity < min_solidity) and allow_dissection:
        req_defects = max(min_bilateral_sinus_count, 4) if is_bipinnatisect_window else min_bilateral_sinus_count
        dissected = is_botanical_dissection(
            main_cnt,
            solidity=solidity,
            p_apex=p_apex,
            p_base=p_base,
            min_solidity=eff_diss_min_sol,
            max_solidity=min_solidity,
            defect_depth_ratio=sinus_defect_min_depth_ratio,
            min_defect_count=req_defects,
        )
        min_req_ucs = 0.55 if is_bipinnatisect_window else (0.80 if is_high_res_mask else 0.70)
        if dissected and (eff_diss_ucs >= min_req_ucs or (not is_high_res_mask and eff_diss_ucs >= 0.55)):
            coords = resample_normalized_contour(mask, num_contour_points, petiole_attachment_pt=p_base) if extract_contours else None
            meta["reflection_applied"] = False
            dissected_tier = "Tier 1 (Bipinnatisect)" if is_bipinnatisect_window else "Tier 1 (Dissected)"
            return LeafRoutingResult(
                tier=dissected_tier,
                assigned_tier="tier1",
                solidity=solidity,
                ucs=eff_diss_ucs,
                is_folded=False,
                is_irregular_fold=False,
                is_dissected=True,
                processed_mask=mask,
                contour_coords=coords,
                metadata=meta
            )

    # 4. Fallback: Tier 2 Hemi-Blade Reflection for Occluded/Damaged Margins
    # Require plausible hemi-blade geometry (solidity >= 0.58 and ucs >= 0.62) to prevent false reflections of stems/clumps
    if solidity < 0.58 or ucs < 0.62:
        return LeafRoutingResult(
            tier="Failed QC",
            assigned_tier="rejected",
            solidity=solidity,
            ucs=ucs,
            is_folded=False,
            is_irregular_fold=False,
            is_dissected=False,
            processed_mask=mask,
            rejection_reason="UNRESOLVED_CLUMP",
            metadata=meta
        )

    path, reflected_mask = extract_tier2_reflected(
        mask=mask,
        catalog_number="TEMP",
        leaf_id=1,
        p_apex=p_apex,
        p_base=p_base,
        output_dir=None
    )

    if reflected_mask is not None:
        ref_ucs, ref_sol, _, _, _ = compute_geometric_metrics(
            reflected_mask, rosette_centroid=rosette_centroid
        )
        coords = (
            resample_normalized_contour(reflected_mask, num_contour_points, petiole_attachment_pt=p_base)
            if extract_contours
            else None
        )
        meta["reflection_applied"] = True
        return LeafRoutingResult(
            tier="Tier 2",
            assigned_tier="tier2",
            solidity=ref_sol,
            ucs=ref_ucs,
            is_folded=False,
            is_irregular_fold=False,
            is_dissected=False,
            processed_mask=reflected_mask,
            contour_coords=coords,
            metadata=meta
        )

    # 5. Failed QC: Determine diagnostic rejection reason
    if is_irregular_fold:
        rejection_reason = "IRREGULAR_FOLD"
    elif is_high_res_mask and (ucs < 0.74 or solidity < max(dissection_min_solidity, 0.65)):
        rejection_reason = "UNRESOLVED_CLUMP"
    elif solidity < dissection_min_solidity or ucs < 0.60:
        rejection_reason = "UNRESOLVED_CLUMP"
    elif solidity < min_solidity:
        rejection_reason = "FAILED_SOLIDITY"
    else:
        rejection_reason = "UNRESOLVED_CLUMP"

    return LeafRoutingResult(
        tier="Failed QC",
        assigned_tier="rejected",
        solidity=solidity,
        ucs=ucs,
        is_folded=False,
        is_irregular_fold=False,
        is_dissected=False,
        processed_mask=mask,
        rejection_reason=rejection_reason,
        metadata=meta
    )


def evaluate_leaf_quality(
    mask_or_contour: np.ndarray,
    config: Optional[Any] = None,
    **kwargs: Any,
) -> LeafRoutingResult:
    """
    Evaluates leaf silhouette quality and botanical dissection geometry, routing intact
    and bipinnatisect/lyrate leaves (solidity >= 0.18 with >= 4 deep bilateral sinuses
    and ucs >= 0.55) to Tier_1_Direct (`is_dissected=True`).
    Accepts either a binary mask `(H, W)` or a contour `(N, 2)` / `(N, 1, 2)`.
    """
    arr = np.asarray(mask_or_contour)
    if arr.ndim >= 2 and (arr.shape[-1] == 2 and (arr.ndim == 3 or arr.shape[0] >= 4 and arr.shape[1] == 2)):
        pts = arr.reshape(-1, 2).astype(np.int32)
        max_x = max(400, int(np.max(pts[:, 0])) + 40) if pts.size > 0 else 400
        max_y = max(400, int(np.max(pts[:, 1])) + 40) if pts.size > 0 else 400
        mask = np.zeros((max_y, max_x), dtype=np.uint8)
        if len(pts) >= 3:
            cv2.fillPoly(mask, [pts], 255)
    else:
        mask = arr
    kwargs.pop("contour", None)
    return classify_and_route_leaf(mask, config=config, **kwargs)


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
# 3B. Automated LM2 -> SAM 2 Prompt Generation & Topological Gatekeeping
# =============================================================================

class SAM2PromptBundle(dict):
    """
    Dictionary-compatible and attribute-accessible bundle of automated SAM 2 prompts
    derived from an initial Stage 2 PointRend binary leaf mask `m_init`.
    """

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(f"'SAM2PromptBundle' object has no attribute '{name}'") from exc


def generate_sam2_prompts_from_lm2_mask(
    m_init: np.ndarray,
    num_lobe_peaks: int = 2,
    erode_kernel_size: int = 5,
    dilate_kernel_size: int = 7,
    min_peak_dist_px: Optional[float] = None,
    max_concavity_points: int = 4,
    logit_scale: float = 6.0,
) -> SAM2PromptBundle:
    """
    Automatically generates deterministic SAM 2 prompt inputs from an initial
    Stage 2 PointRend binary leaf crop mask `m_init` without human intervention:

      1. Positive foreground prompt points (`label = 1`) along the medial skeleton /
         Euclidean distance-transform local maxima (`polygon_interior_point` pole of
         inaccessibility + top `num_lobe_peaks` lobe distance peaks) strictly inside
         the eroded `m_init`.
      2. Negative background prompt points (`label = 0`) in the exterior bounding-box
         corners and deep concavity / sinus pockets (inside the convex hull but
         strictly outside the dilated `m_init`).
      3. Crop bounding box (`[x_min, y_min, x_max, y_max]`) and low-resolution mask
         logits of `m_init` (`(1, 256, 256)` in `[-logit_scale, +logit_scale]`) for
         `SAM2ImagePredictor.predict()`.
    """
    if m_init is None or m_init.size == 0 or np.count_nonzero(m_init) == 0:
        empty_u8 = (
            np.zeros_like(m_init, dtype=np.uint8)
            if (m_init is not None and m_init.ndim == 2)
            else np.zeros((10, 10), dtype=np.uint8)
        )
        return SAM2PromptBundle({
            "point_coords": np.zeros((0, 2), dtype=np.float32),
            "point_labels": np.zeros((0,), dtype=np.int32),
            "positive_points": np.zeros((0, 2), dtype=np.float32),
            "negative_points": np.zeros((0, 2), dtype=np.float32),
            "corner_points": np.zeros((0, 2), dtype=np.float32),
            "concavity_points": np.zeros((0, 2), dtype=np.float32),
            "box": np.zeros((4,), dtype=np.float32),
            "mask_input": np.full((1, 256, 256), -float(logit_scale), dtype=np.float32),
            "eroded_mask": empty_u8,
            "dilated_mask": empty_u8.copy(),
        })

    mask_u8 = (m_init > 0).astype(np.uint8) * 255
    h, w = mask_u8.shape[:2]
    ys, xs = np.where(mask_u8 > 0)

    xmin, xmax = int(np.min(xs)), int(np.max(xs))
    ymin, ymax = int(np.min(ys)), int(np.max(ys))
    span_w = max(1, xmax - xmin)
    span_h = max(1, ymax - ymin)
    pad_box = max(2, int(round(0.04 * max(span_w, span_h))))

    bx0 = max(0, xmin - pad_box)
    by0 = max(0, ymin - pad_box)
    bx1 = min(w - 1, xmax + pad_box)
    by1 = min(h - 1, ymax + pad_box)
    box = np.array([float(bx0), float(by0), float(bx1), float(by1)], dtype=np.float32)

    # 1. Erode m_init and compute positive foreground points via distance transform peaks
    ek = max(3, int(erode_kernel_size) | 1)
    erode_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (ek, ek))
    eroded = cv2.erode(mask_u8, erode_kernel, iterations=1)
    if np.count_nonzero(eroded) < 12 and ek > 3:
        eroded = cv2.erode(mask_u8, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)), iterations=1)
    if np.count_nonzero(eroded) == 0:
        eroded = mask_u8.copy()

    dist_map = cv2.distanceTransform(eroded, cv2.DIST_L2, 5)
    pos_pts: List[Tuple[float, float]] = []

    # Primary interior pole of inaccessibility (matches polygon_interior_point)
    _, max_val, _, max_loc = cv2.minMaxLoc(dist_map)
    if max_val > 0:
        pos_pts.append((float(max_loc[0]), float(max_loc[1])))
    else:
        pos_pts.append((float(xs[len(xs) // 2]), float(ys[len(ys) // 2])))

    # Top `num_lobe_peaks` additional lobe distance peaks via adaptive NMS
    if num_lobe_peaks > 0 and max_val > 0:
        base_min_dist = (
            float(min_peak_dist_px)
            if min_peak_dist_px is not None
            else max(8.0, 0.12 * float(math.hypot(span_w, span_h)))
        )
        work_dist = dist_map.copy()
        p0_x, p0_y = int(round(pos_pts[0][0])), int(round(pos_pts[0][1]))
        r0 = int(round(max(base_min_dist, 1.15 * float(max_val))))
        cv2.circle(work_dist, (p0_x, p0_y), max(1, r0), 0.0, -1)

        # If adaptive radius cleared the entire mask (e.g. small unlobed blade), use compact suppression
        if float(np.max(work_dist)) <= 0.5:
            work_dist = dist_map.copy()
            r0_small = max(3, int(round(min(base_min_dist, 0.25 * max(span_w, span_h)))))
            cv2.circle(work_dist, (p0_x, p0_y), r0_small, 0.0, -1)

        for _ in range(int(num_lobe_peaks)):
            _, pk_val, _, pk_loc = cv2.minMaxLoc(work_dist)
            if pk_val <= 0.0:
                break
            px_i, py_i = int(pk_loc[0]), int(pk_loc[1])
            if eroded[py_i, px_i] > 0:
                pos_pts.append((float(px_i), float(py_i)))
            r_i = int(round(max(base_min_dist, 1.10 * float(pk_val))))
            cv2.circle(work_dist, (px_i, py_i), max(2, r_i), 0.0, -1)

    positive_points = np.array(pos_pts, dtype=np.float32).reshape(-1, 2)

    # 2. Dilate m_init and compute negative background points (corners + concavity pockets)
    dk = max(3, int(dilate_kernel_size) | 1)
    dilate_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (dk, dk))
    dilated = cv2.dilate(mask_u8, dilate_kernel, iterations=1)

    # 2a. Exterior bounding-box corners outside dilated mask
    cx0 = max(1, xmin - dk)
    cy0 = max(1, ymin - dk)
    cx1 = min(w - 2, xmax + dk)
    cy1 = min(h - 2, ymax + dk)
    candidate_corners = [
        (cx0, cy0),
        (cx1, cy0),
        (cx0, cy1),
        (cx1, cy1),
        (1, 1),
        (max(1, w - 2), 1),
        (1, max(1, h - 2)),
        (max(1, w - 2), max(1, h - 2)),
    ]
    corner_pts: List[Tuple[float, float]] = []
    seen_corners = set()
    for qx, qy in candidate_corners:
        qx_c = max(0, min(w - 1, int(qx)))
        qy_c = max(0, min(h - 1, int(qy)))
        if dilated[qy_c, qx_c] == 0 and (qx_c, qy_c) not in seen_corners:
            corner_pts.append((float(qx_c), float(qy_c)))
            seen_corners.add((qx_c, qy_c))
            if len(corner_pts) >= 4:
                break

    # 2b. Deep concavity pockets (lyrate sinuses) inside convex hull but outside dilated m_init
    concavity_pts: List[Tuple[float, float]] = []
    contours, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if contours and max_concavity_points > 0:
        main_cnt = max(contours, key=cv2.contourArea)
        hull = cv2.convexHull(main_cnt)
        hull_mask = np.zeros((h, w), dtype=np.uint8)
        cv2.fillPoly(hull_mask, [hull], 255)

        pocket_mask = np.where((hull_mask > 0) & (dilated == 0), 255, 0).astype(np.uint8)
        if np.count_nonzero(pocket_mask) > 0:
            pocket_dist = cv2.distanceTransform(pocket_mask, cv2.DIST_L2, 5)
            num_pockets, pocket_labels, pocket_stats, _ = cv2.connectedComponentsWithStats(
                pocket_mask, connectivity=8
            )
            pocket_candidates: List[Tuple[float, float, float]] = []
            for lbl in range(1, num_pockets):
                area_p = int(pocket_stats[lbl, cv2.CC_STAT_AREA])
                if area_p < 8:
                    continue
                comp_dist = np.where(pocket_labels == lbl, pocket_dist, 0.0)
                _, p_max_val, _, p_max_loc = cv2.minMaxLoc(comp_dist)
                if p_max_val >= 1.5:
                    qx_p, qy_p = int(p_max_loc[0]), int(p_max_loc[1])
                    if dilated[qy_p, qx_p] == 0:
                        pocket_candidates.append((float(p_max_val), float(qx_p), float(qy_p)))

            pocket_candidates.sort(key=lambda item: item[0], reverse=True)
            for _, qx_p, qy_p in pocket_candidates[: int(max_concavity_points)]:
                concavity_pts.append((qx_p, qy_p))

            # If a single pocket spans multiple sinuses, extract additional peaks via NMS on pocket_dist
            if len(concavity_pts) < int(max_concavity_points):
                work_pdist = pocket_dist.copy()
                supp_p = max(8, int(round(0.10 * max(span_w, span_h))))
                for qx_p, qy_p in concavity_pts:
                    cv2.circle(work_pdist, (int(round(qx_p)), int(round(qy_p))), max(1, supp_p), 0.0, -1)
                for _ in range(int(max_concavity_points) * 3):
                    if len(concavity_pts) >= int(max_concavity_points):
                        break
                    _, p_val, _, p_loc = cv2.minMaxLoc(work_pdist)
                    if p_val < 2.0:
                        break
                    qx_p, qy_p = int(p_loc[0]), int(p_loc[1])
                    if dilated[qy_p, qx_p] == 0:
                        concavity_pts.append((float(qx_p), float(qy_p)))
                    cv2.circle(work_pdist, (qx_p, qy_p), max(1, supp_p), 0.0, -1)

    corner_points = (
        np.array(corner_pts, dtype=np.float32).reshape(-1, 2)
        if corner_pts
        else np.zeros((0, 2), dtype=np.float32)
    )
    concavity_points = (
        np.array(concavity_pts, dtype=np.float32).reshape(-1, 2)
        if concavity_pts
        else np.zeros((0, 2), dtype=np.float32)
    )

    neg_list = corner_pts + concavity_pts
    negative_points = (
        np.array(neg_list, dtype=np.float32).reshape(-1, 2)
        if neg_list
        else np.zeros((0, 2), dtype=np.float32)
    )

    if len(positive_points) > 0 and len(negative_points) > 0:
        point_coords = np.vstack([positive_points, negative_points]).astype(np.float32)
        point_labels = np.concatenate([
            np.ones(len(positive_points), dtype=np.int32),
            np.zeros(len(negative_points), dtype=np.int32),
        ])
    elif len(positive_points) > 0:
        point_coords = positive_points.copy()
        point_labels = np.ones(len(positive_points), dtype=np.int32)
    else:
        point_coords = negative_points.copy()
        point_labels = np.zeros(len(negative_points), dtype=np.int32)

    # 3. Low-resolution mask logits (1, 256, 256) in [-logit_scale, +logit_scale]
    prob_float = (mask_u8 > 0).astype(np.float32)
    prob_smooth = cv2.GaussianBlur(prob_float, (5, 5), 1.0)
    prob_256 = cv2.resize(prob_smooth, (256, 256), interpolation=cv2.INTER_LINEAR)
    logits_256 = ((prob_256 * 2.0) - 1.0) * float(logit_scale)
    mask_input = logits_256[np.newaxis, :, :].astype(np.float32)

    return SAM2PromptBundle({
        "point_coords": point_coords,
        "point_labels": point_labels,
        "positive_points": positive_points,
        "negative_points": negative_points,
        "corner_points": corner_points,
        "concavity_points": concavity_points,
        "box": box,
        "mask_input": mask_input,
        "eroded_mask": eroded,
        "dilated_mask": dilated,
    })


def is_topologically_connected_mask(
    mask: np.ndarray,
    reference_mask: Optional[np.ndarray] = None,
    min_area_px: int = 30,
    min_iou_with_ref: float = 0.35,
    max_hole_area_px: float = 15.0,
) -> bool:
    """
    Validates that a refined binary leaf mask satisfies strict topological connectedness:
      1. Non-empty with foreground pixel count >= `min_area_px`.
      2. Strictly 1 connected foreground component (`num_labels == 2` in 8-connectivity),
         preventing disconnected lobe islands or severed blade fragments.
      3. Simply connected (genus 0; no internal topological holes with area >= `max_hole_area_px`).
      4. When `reference_mask` (`m_init`) is provided, preserves macro-scale alignment
         (IoU >= `min_iou_with_ref` and area ratio within [0.35, 2.50]).
    """
    if mask is None or mask.size == 0:
        return False

    bin_u8 = (mask > 0).astype(np.uint8)
    fg_count = int(np.count_nonzero(bin_u8))
    if fg_count < min_area_px:
        return False

    num_labels, _, _, _ = cv2.connectedComponentsWithStats(bin_u8, connectivity=8)
    if num_labels != 2:
        return False

    # Check for internal holes via two-level contour hierarchy
    contours, hierarchy = cv2.findContours(bin_u8 * 255, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    if hierarchy is not None and len(contours) > 1:
        hier = hierarchy[0]
        for idx, h_node in enumerate(hier):
            parent_idx = int(h_node[3])
            if parent_idx != -1:
                hole_area = float(cv2.contourArea(contours[idx]))
                if hole_area >= max_hole_area_px:
                    return False

    if reference_mask is not None and reference_mask.size == mask.size:
        ref_bin = (reference_mask > 0)
        ref_area = int(np.count_nonzero(ref_bin))
        if ref_area > 0:
            inter = int(np.count_nonzero(bin_u8.astype(bool) & ref_bin))
            union = int(np.count_nonzero(bin_u8.astype(bool) | ref_bin))
            iou = float(inter) / float(max(1, union))
            area_ratio = float(fg_count) / float(ref_area)
            if iou < min_iou_with_ref or not (0.35 <= area_ratio <= 2.50):
                return False

    return True


# =============================================================================
# 4. Adaptive DBSCAN Spatial Clustering
# =============================================================================

def cluster_plant_individuals(
    instances: List[Any],
    sheet_width: int,
    sheet_height: int,
    eps_ratio: float = 0.15
) -> List[Any]:
    """
    Vectorized centroid clustering across detected organs on a specimen sheet.
    Assigns plant_individual_id to avoid trait averaging across distinct plants.

    Parameters
    ----------
    instances : List[Any]
        List of candidate objects having `.bbox` (ymin, xmin, ymax, xmax) and
        `.plant_individual_id`.
    sheet_width : int
        Width of specimen image sheet.
    sheet_height : int
        Height of specimen image sheet.
    eps_ratio : float
        Proportion of max sheet dimension used as DBSCAN neighbourhood radius.

    Returns
    -------
    List[Any]
        Updated instances with plant_individual_id assigned.
    """
    if not instances:
        return instances

    if len(instances) == 1:
        instances[0].plant_individual_id = 0
        return instances

    centroids = []
    for inst in instances:
        ymin, xmin, ymax, xmax = inst.bbox
        cx = (xmin + xmax) / 2.0
        cy = (ymin + ymax) / 2.0
        centroids.append([cx, cy])

    coords = np.array(centroids, dtype=np.float32)
    eps = eps_ratio * float(max(sheet_width, sheet_height))

    clustering = DBSCAN(eps=eps, min_samples=1).fit(coords)
    labels = clustering.labels_

    unique_labels = sorted(list(set(labels)))
    label_map = {lbl: idx for idx, lbl in enumerate(unique_labels)}

    for inst, lbl in zip(instances, labels):
        inst.plant_individual_id = label_map.get(lbl, 0)

    return instances


def cluster_voucher_plants_dbscan(
    candidates: List[Any],
    sheet_width: int = 3000,
    sheet_height: int = 4000
) -> List[Any]:
    """
    Applies adaptive DBSCAN spatial clustering (eps ≈ 0.15 * sheet_width, min_samples=2)
    to group detected organs into distinct plant_individual_id clusters.
    Outlier noise leaves (label == -1) receive unique individual IDs so no leaf is lost.
    """
    if not candidates:
        return candidates

    centroids = []
    for cand in candidates:
        ymin, xmin, ymax, xmax = cand.bbox
        if ymax > ymin and xmax > xmin:
            cx = (xmin + xmax) / 2.0
            cy = (ymin + ymax) / 2.0
        else:
            cx = sheet_width / 2.0
            cy = sheet_height / 2.0
        centroids.append([cx, cy])

    coords = np.array(centroids, dtype=np.float32)

    if len(candidates) <= 1:
        candidates[0].plant_individual_id = 0
        return candidates

    adaptive_eps = 0.15 * float(sheet_width)
    db = DBSCAN(eps=adaptive_eps, min_samples=2).fit(coords)
    labels = db.labels_

    unique_clusters = sorted([l for l in set(labels) if l >= 0])
    cluster_mapping = {c: i for i, c in enumerate(unique_clusters)}

    next_id = len(cluster_mapping)
    for i, cand in enumerate(candidates):
        raw_label = labels[i]
        if raw_label >= 0:
            cand.plant_individual_id = cluster_mapping[raw_label]
        else:
            cand.plant_individual_id = next_id
            next_id += 1

    return candidates


# =============================================================================
# 5. Decoupled Asynchronous Ruler Detector (Hough Transform)
# =============================================================================

def detect_ruler_scale_hough(image: np.ndarray) -> Optional[float]:
    """
    Asynchronous Hough-transform ruler detection.
    Scans margin regions of the herbarium sheet for periodic tick marks.
    Returns pixels_per_mm or None if not confidently detected.
    """
    if image is None or image.size == 0:
        return None

    h, w = image.shape[:2]
    # Rulers are typically along sheet borders (margins)
    margins = [
        image[int(h * 0.75):, :],               # Bottom
        image[:int(h * 0.25), :],               # Top
        image[:, :int(w * 0.20)],               # Left
        image[:, int(w * 0.80):],               # Right
    ]

    candidate_scales: List[float] = []

    for region in margins:
        if region.size == 0:
            continue
        gray = cv2.cvtColor(region, cv2.COLOR_BGR2GRAY) if len(region.shape) == 3 else region

        # Downscale if extremely large for rapid Hough transform
        rh, rw = gray.shape[:2]
        max_dim = max(rh, rw)
        scale_factor = 1.0
        if max_dim > 1500:
            scale_factor = 1500.0 / max_dim
            gray = cv2.resize(gray, (int(rw * scale_factor), int(rh * scale_factor)), interpolation=cv2.INTER_AREA)

        # Bilateral / Gaussian blur + Canny edge detection
        blurred = cv2.GaussianBlur(gray, (3, 3), 0)
        edges = cv2.Canny(blurred, 30, 100)

        # Detect line segments
        min_line_len = max(int(12 * scale_factor), 10)
        lines = cv2.HoughLinesP(edges, rho=1, theta=np.pi / 180, threshold=25,
                                minLineLength=min_line_len, maxLineGap=4)
        if lines is None or len(lines) < 8:
            continue

        # Extract line angles and lengths
        line_data = []
        for line in lines:
            pts = line.reshape(-1)
            if len(pts) < 4:
                continue
            x1, y1, x2, y2 = pts[:4]
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy)
            angle = math.degrees(math.atan2(dy, dx)) % 180.0
            line_data.append((x1, y1, x2, y2, angle, length))

        # Check for dominant parallel angle (horizontal or vertical ticks)
        angles = np.array([ld[4] for ld in line_data])
        hist, bin_edges = np.histogram(angles, bins=18, range=(0, 180))
        dominant_bin = np.argmax(hist)
        if hist[dominant_bin] < 6:
            continue

        dom_angle_low = bin_edges[dominant_bin]
        dom_angle_high = bin_edges[dominant_bin + 1]

        # Filter lines matching dominant tick orientation
        matching_lines = [ld for ld in line_data if dom_angle_low <= ld[4] <= dom_angle_high]
        if len(matching_lines) < 6:
            continue

        # Project lines onto orthogonal axis to determine inter-tick pitch
        dom_angle = (dom_angle_low + dom_angle_high) / 2.0
        ortho_rad = math.radians(dom_angle + 90.0)
        proj_coords = []
        for ld in matching_lines:
            mid_x, mid_y = (ld[0] + ld[2]) / 2.0, (ld[1] + ld[3]) / 2.0
            p = mid_x * math.cos(ortho_rad) + mid_y * math.sin(ortho_rad)
            proj_coords.append(p)

        proj_coords = np.sort(proj_coords)
        tick_positions = []
        for p in proj_coords:
            if not tick_positions or abs(p - tick_positions[-1]) > 5.0:
                tick_positions.append(p)

        if len(tick_positions) >= 5:
            deltas = np.diff(tick_positions)
            med_delta = float(np.median(deltas))
            mad = float(np.median(np.abs(deltas - med_delta)))
            # Consistent periodic tick spacing condition
            if mad < 0.25 * med_delta and (8.0 <= med_delta / scale_factor <= 150.0):
                px_per_mm = med_delta / scale_factor
                candidate_scales.append(px_per_mm)

    if candidate_scales:
        return float(np.median(candidate_scales))
    return None


# =============================================================================
# 6. Additional Caliper & Rosette Utilities (Tier 3 & Tier 4)
# =============================================================================

def route_tier3_open_curve(
    mask: np.ndarray,
    candidate: Any,
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

    cat_num = getattr(candidate, "catalog_number", "UNKNOWN")
    plant_id = getattr(candidate, "plant_individual_id", 0)
    leaf_id = getattr(candidate, "leaf_id", 1)

    curve_df = pd.DataFrame({
        "catalogNumber": cat_num,
        "plant_individual_id": plant_id,
        "leaf_id": leaf_id,
        "point_index": np.arange(len(resampled)),
        "x_px": resampled[:, 0],
        "y_px": resampled[:, 1],
        "x_norm": (resampled[:, 0] - min_x) / span_x,
        "y_norm": (resampled[:, 1] - min_y) / span_y,
    })

    save_path = output_dir / "masks" / "tier3_open_curves" / f"{cat_num}_p{plant_id}_leaf{leaf_id}_curve.csv"
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
    plant_candidates: List[Any],
    plant_id: int,
    output_dir: Path,
    raw_img: Optional[np.ndarray] = None,
    catalog_number: Optional[str] = None,
    precropped_rosette_bgr: Optional[np.ndarray] = None,
) -> Optional[str]:
    """Tier 4: Dense whole-rosette bounding crop from native-DPI sheet."""
    if precropped_rosette_bgr is not None:
        rosette_crop = precropped_rosette_bgr
    else:
        if raw_img is None:
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

    cat = (
        catalog_number
        or (getattr(plant_candidates[0], "catalog_number", None) if plant_candidates else None)
        or (raw_image_path.stem if raw_image_path else None)
        or "UNKNOWN"
    )
    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    save_path = out_base / "cropped_patches" / "rosettes_dense" / f"{cat}_p{plant_id}_rosette.jpg"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(save_path), rosette_crop)
    return str(save_path)


# =============================================================================
# 7. Empirical Basal Petiole Length & Flowering Scape Height Approximations
# =============================================================================

def _extract_xy_point(item: Any, ref_pt: Optional[Tuple[float, float]] = None) -> Optional[Tuple[float, float]]:
    """Convert a 2-tuple point, 4-tuple bbox, or dict into a sheet-space (x, y) coordinate."""
    if item is None:
        return None
    if isinstance(item, dict):
        if "x" in item and "y" in item:
            return float(item["x"]), float(item["y"])
        if "centroid" in item and item["centroid"] is not None:
            c = item["centroid"]
            return float(c[0]), float(c[1])
        if "p_bend" in item and item["p_bend"] is not None:
            b = item["p_bend"]
            return float(b[0]), float(b[1])
        if "bbox" in item and item["bbox"] is not None:
            return _extract_xy_point(item["bbox"], ref_pt=ref_pt)
        return None
    if hasattr(item, "bbox") and not isinstance(item, (tuple, list, np.ndarray)):
        return _extract_xy_point(getattr(item, "bbox"), ref_pt=ref_pt)
    if isinstance(item, (tuple, list, np.ndarray)):
        seq = list(item)
        if len(seq) >= 3 and isinstance(seq[0], (tuple, list, np.ndarray)):
            # Stage 1 detection tuple: (bbox, score, cls_id) or (bbox, mask, score, cls_id)
            cls_id = int(seq[-1]) if isinstance(seq[-1], (int, float, np.integer)) else None
            if cls_id is not None and cls_id not in (0, 1, 3, 4, 5, 6, 7, 8, 9, 10):
                return None
            return _extract_xy_point(seq[0], ref_pt=ref_pt)
        if len(seq) == 2 and all(isinstance(v, (int, float, np.number)) for v in seq):
            return float(seq[0]), float(seq[1])
        if len(seq) == 4 and all(isinstance(v, (int, float, np.number)) for v in seq):
            v0, v1, v2, v3 = float(seq[0]), float(seq[1]), float(seq[2]), float(seq[3])
            # Canonical repository bbox order is (ymin, xmin, ymax, xmax) -> (cx, cy) = (0.5*(v1+v3), 0.5*(v0+v2))
            cand_yx = (0.5 * (v1 + v3), 0.5 * (v0 + v2))
            cand_xy = (0.5 * (v0 + v2), 0.5 * (v1 + v3))
            if ref_pt is not None:
                d_yx = math.hypot(cand_yx[0] - ref_pt[0], cand_yx[1] - ref_pt[1])
                d_xy = math.hypot(cand_xy[0] - ref_pt[0], cand_xy[1] - ref_pt[1])
                if d_xy + 15.0 < d_yx:
                    return cand_xy
            return cand_yx
    return None


def _parse_basal_leaf_pose(item: Any) -> Optional[Dict[str, Any]]:
    """Extract standardized (p_base, p_apex, plant_individual_id, geodesic_petiole_skeleton_px, blade_length_px) from a leaf record."""
    if item is None:
        return None

    pid = 0
    ox, oy = 0.0, 0.0
    p_base: Optional[Tuple[float, float]] = None
    p_apex: Optional[Tuple[float, float]] = None
    geodesic_px = 0.0
    blade_px: Optional[float] = None
    leaf_id: Optional[int] = None

    if isinstance(item, dict):
        pid = int(item.get("plant_individual_id", 0) or 0)
        leaf_id = item.get("leaf_id")
        origin = item.get("mask_origin") or item.get("origin")
        if origin is not None and len(origin) == 2 and not item.get("coords_in_sheet_space", False):
            ox, oy = float(origin[0]), float(origin[1])

        if item.get("p_base_sheet") is not None:
            pb = item["p_base_sheet"]
            p_base = (float(pb[0]), float(pb[1]))
        elif item.get("p_base") is not None:
            pb = item["p_base"]
            p_base = (float(pb[0]) + ox, float(pb[1]) + oy)

        if item.get("p_apex_sheet") is not None:
            pa = item["p_apex_sheet"]
            p_apex = (float(pa[0]), float(pa[1]))
        elif item.get("p_apex") is not None:
            pa = item["p_apex"]
            p_apex = (float(pa[0]) + ox, float(pa[1]) + oy)

        if (p_base is None or p_apex is None) and isinstance(item.get("mask"), np.ndarray):
            _, _, _, pa_m, pb_m = compute_geometric_metrics(item["mask"])
            if p_base is None:
                p_base = (float(pb_m[0]) + ox, float(pb_m[1]) + oy)
            if p_apex is None:
                p_apex = (float(pa_m[0]) + ox, float(pa_m[1]) + oy)

        trim_info = item.get("trim_meta") or item.get("petiole_trim_info") or {}
        geodesic_px = float(
            item.get(
                "geodesic_petiole_skeleton_px",
                item.get("trim_distance_px", trim_info.get("trim_distance_px", 0.0)),
            )
            or 0.0
        )
        if item.get("blade_length_px") is not None and float(item["blade_length_px"]) > 0:
            blade_px = float(item["blade_length_px"])

    elif isinstance(item, np.ndarray) and item.ndim == 2:
        trimmed_m, trim_m = trim_petiole_tail(item, return_metadata=True)
        _, _, _, pa_m, pb_m = compute_geometric_metrics(trimmed_m)
        p_base = (float(pb_m[0]), float(pb_m[1]))
        p_apex = (float(pa_m[0]), float(pa_m[1]))
        geodesic_px = float(trim_m.get("trim_distance_px", 0.0) or 0.0)

    elif isinstance(item, (tuple, list)) and len(item) >= 2:
        first, second = item[0], item[1]
        if (
            isinstance(first, (tuple, list, np.ndarray))
            and isinstance(second, (tuple, list, np.ndarray))
            and len(first) == 2
            and len(second) == 2
        ):
            p_base = (float(first[0]), float(first[1]))
            p_apex = (float(second[0]), float(second[1]))
            if len(item) >= 3 and isinstance(item[2], (int, float, np.number)):
                geodesic_px = float(item[2])
            if len(item) >= 4 and isinstance(item[3], (int, np.integer)):
                pid = int(item[3])

    elif hasattr(item, "p_base") or hasattr(item, "mask"):
        pid = int(getattr(item, "plant_individual_id", 0) or 0)
        leaf_id = getattr(item, "leaf_id", None)
        origin = getattr(item, "mask_origin", (0, 0)) or (0, 0)
        ox, oy = float(origin[0]), float(origin[1])
        pb = getattr(item, "p_base", None)
        pa = getattr(item, "p_apex", None)
        if (pb is None or pa is None) and getattr(item, "mask", None) is not None:
            _, _, _, pa, pb = compute_geometric_metrics(getattr(item, "mask"))
        if pb is not None and pa is not None:
            p_base = (float(pb[0]) + ox, float(pb[1]) + oy)
            p_apex = (float(pa[0]) + ox, float(pa[1]) + oy)
        trim_info = getattr(item, "trim_meta", None) or {}
        geodesic_px = float(
            getattr(item, "geodesic_petiole_skeleton_px", trim_info.get("trim_distance_px", 0.0)) or 0.0
        )

    if p_base is None or p_apex is None:
        return None

    chord_len = math.hypot(p_apex[0] - p_base[0], p_apex[1] - p_base[1])
    if blade_px is None or blade_px <= 0.0:
        blade_px = chord_len

    return {
        "leaf_id": leaf_id,
        "plant_individual_id": pid,
        "p_base": p_base,
        "p_apex": p_apex,
        "geodesic_petiole_skeleton_px": max(0.0, geodesic_px),
        "blade_length_px": max(0.0, float(blade_px)),
    }


def _estimate_single_plant_crown(
    poses: List[Dict[str, Any]],
    root_pts: List[Tuple[float, float]],
) -> Tuple[float, float]:
    """
    Compute the caudex crown point `p_crown = (x_crown, y_crown)` for a single plant individual
    via least-squares intersection of inward basal-leaf midrib rays anchored by Stage 1 Class 9
    root detections and/or the innermost `p_base` coordinates.
    """
    if not poses:
        if root_pts:
            rx = float(np.median([p[0] for p in root_pts]))
            ry = float(np.median([p[1] for p in root_pts]))
            return (rx, ry)
        return (0.0, 0.0)

    bases = np.array([p["p_base"] for p in poses], dtype=np.float64)
    apices = np.array([p["p_apex"] for p in poses], dtype=np.float64)
    vecs = bases - apices  # Inward midrib vector: apex -> base -> caudex crown
    lengths = np.hypot(vecs[:, 0], vecs[:, 1])
    valid_mask = lengths > 1e-4

    mean_base = np.median(bases, axis=0)
    mean_len = float(np.median(lengths[valid_mask])) if np.any(valid_mask) else 50.0

    # Filter root points to those associated with this plant individual
    near_roots: List[Tuple[float, float]] = []
    if root_pts:
        max_root_dist = max(4.0 * mean_len, 450.0)
        for rp in root_pts:
            if math.hypot(rp[0] - float(mean_base[0]), rp[1] - float(mean_base[1])) <= max_root_dist:
                near_roots.append(rp)
        if not near_roots and len(root_pts) == 1:
            near_roots = list(root_pts)

    # Determine innermost p_base coordinate anchor
    if near_roots:
        root_anchor = np.array(
            [np.median([p[0] for p in near_roots]), np.median([p[1] for p in near_roots])],
            dtype=np.float64,
        )
        dists_to_root = np.hypot(bases[:, 0] - root_anchor[0], bases[:, 1] - root_anchor[1])
        inner_idx = int(np.argmin(dists_to_root))
        innermost_base = bases[inner_idx]
        p_anchor = 0.70 * root_anchor + 0.30 * innermost_base
        has_root_anchor = True
    else:
        has_root_anchor = False
        if np.any(valid_mask):
            unit_dirs = vecs[valid_mask] / lengths[valid_mask, None]
            mean_dir = np.sum(unit_dirs, axis=0)
            norm_dir = float(np.linalg.norm(mean_dir))
            if norm_dir > 0.25 and len(bases) >= 2:
                mean_dir = mean_dir / norm_dir
                # Innermost p_base is furthest along the inward midrib direction and close to the lateral median
                proj_inward = (bases - mean_base) @ mean_dir
                lat_vec = np.array([-mean_dir[1], mean_dir[0]], dtype=np.float64)
                lat_dist = np.abs((bases - mean_base) @ lat_vec)
                scores = proj_inward - 0.35 * lat_dist
                k_inner = max(1, int(math.ceil(0.5 * len(bases))))
                top_idx = np.argsort(scores)[-k_inner:]
                innermost_base = np.median(bases[top_idx], axis=0)
            else:
                dists_center = np.hypot(bases[:, 0] - mean_base[0], bases[:, 1] - mean_base[1])
                innermost_base = bases[int(np.argmin(dists_center))]
        else:
            innermost_base = mean_base
        p_anchor = innermost_base.copy()

    if int(np.sum(valid_mask)) < 2:
        if has_root_anchor:
            if int(np.sum(valid_mask)) == 1:
                idx0 = int(np.where(valid_mask)[0][0])
                u0 = vecs[idx0] / lengths[idx0]
                b0 = bases[idx0]
                t_root = float(np.dot(p_anchor - b0, u0))
                if t_root >= 0.0:
                    proj_pt = b0 + min(t_root, 2.5 * mean_len) * u0
                    crown_pt = 0.65 * proj_pt + 0.35 * p_anchor
                    return float(crown_pt[0]), float(crown_pt[1])
            return float(p_anchor[0]), float(p_anchor[1])
        return float(innermost_base[0]), float(innermost_base[1])

    # Least-squares intersection of 2D inward midrib lines:
    # Each line i has normal projector P_i = I - u_i u_i^T passing through b_i.
    A_mat = np.zeros((2, 2), dtype=np.float64)
    b_vec = np.zeros(2, dtype=np.float64)
    valid_indices = np.where(valid_mask)[0]
    unit_list = []
    for idx_i in valid_indices:
        u_i = vecs[idx_i] / lengths[idx_i]
        unit_list.append(u_i)
        P_i = np.eye(2, dtype=np.float64) - np.outer(u_i, u_i)
        A_mat += P_i
        b_vec += P_i @ bases[idx_i]

    eigvals = np.linalg.eigvalsh(A_mat)
    min_eig = float(np.min(eigvals))
    max_reach = max(2.5 * float(np.max(lengths[valid_mask])), 60.0)

    if min_eig > 1e-3:
        # Rays have distinct angles -> solve least-squares intersection with light anchor regularization
        reg_w = 0.35 if has_root_anchor else 1e-4
        p_ls = np.linalg.solve(A_mat + reg_w * np.eye(2, dtype=np.float64), b_vec + reg_w * p_anchor)
        u_arr = np.array(unit_list, dtype=np.float64)
        b_arr = bases[valid_indices]
        t_vals = np.sum((p_ls[None, :] - b_arr) * u_arr, axis=1)
        dist_from_inner = float(np.linalg.norm(p_ls - innermost_base))

        if float(np.median(t_vals)) >= -0.15 * mean_len and dist_from_inner <= max_reach:
            return float(p_ls[0]), float(p_ls[1])
        if dist_from_inner > max_reach and float(np.median(t_vals)) > 0.0:
            # Clamp excessive ray intersection distance to max_reach along the ray direction from innermost_base
            direction = (p_ls - innermost_base) / max(dist_from_inner, 1e-6)
            clamped = innermost_base + max_reach * direction
            if has_root_anchor:
                clamped = 0.5 * clamped + 0.5 * p_anchor
            return float(clamped[0]), float(clamped[1])

    # Fallback when rays are nearly parallel or divergent: anchor by roots / innermost p_base
    return float(p_anchor[0]), float(p_anchor[1])


def estimate_rosette_crown_point(
    basal_leaves: Any = None,
    root_boxes: Optional[Sequence[Any]] = None,
    p_bases: Optional[Sequence[Tuple[float, float]]] = None,
    p_apices: Optional[Sequence[Tuple[float, float]]] = None,
    return_by_plant: bool = False,
    **kwargs: Any,
) -> Union[Tuple[float, float], Dict[int, Tuple[float, float]]]:
    """
    Compute the caudex crown point `p_crown` per `plant_individual_id` via least-squares
    intersection of inward basal-leaf midrib rays (`p_apex -> p_base`) anchored by
    Stage 1 Class 9 (`roots`) detections and/or the innermost `p_base` coordinates.

    Returns:
      - `Tuple[float, float]` `(x_crown, y_crown)` when invoked on a single plant's leaf sequence
        (with `return_by_plant=False`).
      - `Dict[int, Tuple[float, float]]` mapping `plant_individual_id -> (x_crown, y_crown)` when
        `basal_leaves` is a dict keyed by `plant_individual_id` or when `return_by_plant=True`.
    """
    raw_roots = (
        root_boxes
        if root_boxes is not None
        else (
            kwargs.get("root_detections")
            or kwargs.get("roots")
            or kwargs.get("root_pts")
            or kwargs.get("stage1_roots")
        )
    )

    if isinstance(basal_leaves, dict) and not ("p_base" in basal_leaves or "mask" in basal_leaves):
        out_by_pid: Dict[int, Tuple[float, float]] = {}
        for pid_k, leaves_k in basal_leaves.items():
            pt_k = estimate_rosette_crown_point(
                leaves_k,
                root_boxes=raw_roots,
                return_by_plant=False,
            )
            out_by_pid[int(pid_k)] = (float(pt_k[0]), float(pt_k[1]))
        return out_by_pid

    poses: List[Dict[str, Any]] = []
    if p_bases is not None and p_apices is not None:
        for pb, pa in zip(p_bases, p_apices):
            parsed = _parse_basal_leaf_pose((pb, pa))
            if parsed is not None:
                poses.append(parsed)
    elif isinstance(basal_leaves, (dict, np.ndarray)):
        parsed = _parse_basal_leaf_pose(basal_leaves)
        if parsed is not None:
            poses.append(parsed)
    elif basal_leaves is not None:
        for item in basal_leaves:
            parsed = _parse_basal_leaf_pose(item)
            if parsed is not None:
                poses.append(parsed)

    ref_pt: Optional[Tuple[float, float]] = None
    if poses:
        ref_pt = (
            float(np.median([p["p_base"][0] for p in poses])),
            float(np.median([p["p_base"][1] for p in poses])),
        )

    root_pts: List[Tuple[float, float]] = []
    if raw_roots is not None:
        if isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 2 and all(
            isinstance(v, (int, float, np.number)) for v in raw_roots
        ):
            root_pts.append((float(raw_roots[0]), float(raw_roots[1])))
        elif isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 4 and all(
            isinstance(v, (int, float, np.number)) for v in raw_roots
        ):
            pt = _extract_xy_point(raw_roots, ref_pt=ref_pt)
            if pt is not None:
                root_pts.append(pt)
        else:
            for r_item in raw_roots:
                if isinstance(r_item, (tuple, list)) and len(r_item) >= 3 and isinstance(r_item[0], (tuple, list, np.ndarray)):
                    cls_id = int(r_item[-1]) if isinstance(r_item[-1], (int, float, np.integer)) else 9
                    if cls_id != 9:
                        continue
                elif isinstance(r_item, dict) and "class_id" in r_item and int(r_item["class_id"]) != 9:
                    continue
                pt = _extract_xy_point(r_item, ref_pt=ref_pt)
                if pt is not None:
                    root_pts.append(pt)

    if return_by_plant:
        by_pid: Dict[int, List[Dict[str, Any]]] = {}
        for p in poses:
            by_pid.setdefault(int(p["plant_individual_id"]), []).append(p)
        if not by_pid:
            return {0: _estimate_single_plant_crown([], root_pts)}
        return {
            pid_k: _estimate_single_plant_crown(pid_poses, root_pts)
            for pid_k, pid_poses in sorted(by_pid.items())
        }

    return _estimate_single_plant_crown(poses, root_pts)


def approximate_petiole_and_scape_metrics(
    basal_leaves: Any = None,
    p_crown: Optional[Any] = None,
    px_per_cm: Optional[float] = None,
    p_corymb: Optional[Tuple[float, float]] = None,
    p_bend: Optional[Tuple[float, float]] = None,
    upper_cauline_detections: Optional[Sequence[Any]] = None,
    stem_skeleton_mask: Optional[np.ndarray] = None,
    root_boxes: Optional[Sequence[Any]] = None,
    capitula_detections: Optional[Sequence[Any]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    Compute empirical approximations for basal petiole length and flowering plant height:
      1. Per-leaf:
         - `petiole_length_px = max(geodesic_petiole_skeleton_px, ||p_base - p_crown||_2)`
         - `blade_length_px = ||p_apex - p_base||_2`
         - `petiole_to_blade_ratio = petiole_length_px / blade_length_px`
         - `petiole_length_mm` and `blade_length_mm` (using `px_per_cm`, where `px_per_mm = px_per_cm / 10.0`)
      2. Specimen-level outer-rosette aggregation (75th percentile across basal leaves per specimen):
         - `petiole_length_mm`, `blade_length_mm`, `petiole_to_blade_ratio`
      3. Per-voucher flowering scape height:
         - Traces the piecewise path from `p_crown` through any upper stem-bend apex (`p_bend`
           from upper cauline detections / specimen stem skeleton when bent) to the flowering
           capitula centroid (`p_corymb`), exporting `plant_height_cm`, `scape_to_basal_leaf_ratio`,
           and `is_bent_stem`.
    """
    # Resolve scale conversion (px_per_cm and px_per_mm = px_per_cm / 10.0)
    eff_px_per_cm: Optional[float] = None
    if px_per_cm is not None and not np.isnan(px_per_cm) and float(px_per_cm) > 0.0:
        eff_px_per_cm = float(px_per_cm)
    else:
        ppm = kwargs.get("pixels_per_mm")
        if ppm is not None and not np.isnan(ppm) and float(ppm) > 0.0:
            eff_px_per_cm = float(ppm) * 10.0
    px_per_mm = (eff_px_per_cm / 10.0) if (eff_px_per_cm is not None and eff_px_per_cm > 0.0) else None

    # Parse basal leaf poses
    poses: List[Dict[str, Any]] = []
    if isinstance(basal_leaves, (dict, np.ndarray)):
        parsed = _parse_basal_leaf_pose(basal_leaves)
        if parsed is not None:
            poses.append(parsed)
    elif basal_leaves is not None:
        for item in basal_leaves:
            parsed = _parse_basal_leaf_pose(item)
            if parsed is not None:
                poses.append(parsed)
    elif kwargs.get("p_base") is not None and kwargs.get("p_apex") is not None:
        parsed = _parse_basal_leaf_pose({
            "p_base": kwargs["p_base"],
            "p_apex": kwargs["p_apex"],
            "geodesic_petiole_skeleton_px": kwargs.get("geodesic_petiole_skeleton_px", 0.0),
            "blade_length_px": kwargs.get("blade_length_px"),
        })
        if parsed is not None:
            poses.append(parsed)

    # Resolve p_crown per plant_individual_id and primary voucher crown
    crowns_by_pid: Dict[int, Tuple[float, float]] = {}
    if isinstance(p_crown, dict):
        for k, v in p_crown.items():
            if v is not None and len(v) == 2:
                crowns_by_pid[int(k)] = (float(v[0]), float(v[1]))
    elif p_crown is not None and isinstance(p_crown, (tuple, list, np.ndarray)) and len(p_crown) == 2:
        pt_c = (float(p_crown[0]), float(p_crown[1]))
        for p in poses:
            crowns_by_pid[int(p["plant_individual_id"])] = pt_c
        crowns_by_pid.setdefault(0, pt_c)

    if not crowns_by_pid:
        est_crowns = estimate_rosette_crown_point(
            poses,
            root_boxes=root_boxes or kwargs.get("root_detections") or kwargs.get("roots"),
            return_by_plant=True,
        )
        if isinstance(est_crowns, dict):
            crowns_by_pid.update(est_crowns)
        else:
            crowns_by_pid[0] = est_crowns

    primary_crown = crowns_by_pid.get(0) or next(iter(crowns_by_pid.values()), (0.0, 0.0))

    # Compute per-leaf metrics
    leaf_metrics: List[Dict[str, Any]] = []
    for idx_l, p in enumerate(poses):
        pid = int(p["plant_individual_id"])
        c_pt = crowns_by_pid.get(pid, primary_crown)
        pb = p["p_base"]
        pa = p["p_apex"]
        geodesic_px = float(p["geodesic_petiole_skeleton_px"])
        crown_dist_px = float(math.hypot(pb[0] - c_pt[0], pb[1] - c_pt[1]))
        petiole_px = max(geodesic_px, crown_dist_px)
        blade_px = float(p["blade_length_px"])
        if blade_px <= 0.0:
            blade_px = float(math.hypot(pa[0] - pb[0], pa[1] - pb[1]))
        pet_to_blade = (petiole_px / blade_px) if blade_px > 1e-6 else np.nan
        petiole_mm = (petiole_px / px_per_mm) if (px_per_mm is not None and px_per_mm > 0.0) else np.nan
        blade_mm = (blade_px / px_per_mm) if (px_per_mm is not None and px_per_mm > 0.0) else np.nan

        leaf_metrics.append({
            "leaf_id": p["leaf_id"] if p["leaf_id"] is not None else (idx_l + 1),
            "plant_individual_id": pid,
            "p_base": pb,
            "p_apex": pa,
            "p_crown": c_pt,
            "geodesic_petiole_skeleton_px": round(geodesic_px, 4),
            "crown_distance_px": round(crown_dist_px, 4),
            "petiole_length_px": round(petiole_px, 4),
            "blade_length_px": round(blade_px, 4),
            "petiole_to_blade_ratio": round(float(pet_to_blade), 4) if np.isfinite(pet_to_blade) else np.nan,
            "petiole_length_mm": round(float(petiole_mm), 3) if np.isfinite(petiole_mm) else np.nan,
            "blade_length_mm": round(float(blade_mm), 3) if np.isfinite(blade_mm) else np.nan,
        })

    def _agg_p75(records: List[Dict[str, Any]], key: str, digits: int = 4) -> float:
        vals = [float(r[key]) for r in records if r.get(key) is not None and np.isfinite(float(r[key]))]
        if not vals:
            return np.nan
        return round(float(np.percentile(vals, 75)), digits)

    # Aggregate 75th-percentile outer-rosette metrics per specimen (plant_individual_id) and overall
    by_plant_summary: Dict[int, Dict[str, Any]] = {}
    pids = sorted({int(r["plant_individual_id"]) for r in leaf_metrics})
    for pid in pids:
        sub_r = [r for r in leaf_metrics if int(r["plant_individual_id"]) == pid]
        by_plant_summary[pid] = {
            "p_crown": crowns_by_pid.get(pid, primary_crown),
            "petiole_length_px": _agg_p75(sub_r, "petiole_length_px", 4),
            "blade_length_px": _agg_p75(sub_r, "blade_length_px", 4),
            "petiole_length_mm": _agg_p75(sub_r, "petiole_length_mm", 3),
            "blade_length_mm": _agg_p75(sub_r, "blade_length_mm", 3),
            "petiole_to_blade_ratio": _agg_p75(sub_r, "petiole_to_blade_ratio", 4),
        }

    petiole_px_p75 = _agg_p75(leaf_metrics, "petiole_length_px", 4)
    blade_px_p75 = _agg_p75(leaf_metrics, "blade_length_px", 4)
    petiole_mm_p75 = _agg_p75(leaf_metrics, "petiole_length_mm", 3)
    blade_mm_p75 = _agg_p75(leaf_metrics, "blade_length_mm", 3)
    pet_to_blade_p75 = _agg_p75(leaf_metrics, "petiole_to_blade_ratio", 4)

    # Resolve flowering capitula centroid p_corymb
    resolved_corymb: Optional[Tuple[float, float]] = None
    if p_corymb is not None and len(p_corymb) == 2:
        resolved_corymb = (float(p_corymb[0]), float(p_corymb[1]))
    else:
        raw_caps = capitula_detections or kwargs.get("capitula_boxes") or kwargs.get("capitula_pts")
        if raw_caps:
            cap_pts = []
            for c_item in raw_caps:
                pt = _extract_xy_point(c_item, ref_pt=primary_crown)
                if pt is not None:
                    cap_pts.append(pt)
            if cap_pts:
                resolved_corymb = (
                    float(np.median([p[0] for p in cap_pts])),
                    float(np.median([p[1] for p in cap_pts])),
                )

    # Trace piecewise scape path from p_crown through any upper stem-bend apex (p_bend) to p_corymb
    is_bent_stem = False
    resolved_bend: Optional[Tuple[float, float]] = None
    scape_length_px = np.nan
    plant_height_cm = np.nan
    scape_to_basal_leaf_ratio = np.nan

    if resolved_corymb is not None:
        straight_px = float(math.hypot(resolved_corymb[0] - primary_crown[0], resolved_corymb[1] - primary_crown[1]))
        explicit_bend = p_bend or kwargs.get("bend_apex") or kwargs.get("stem_bend_pt")
        cand_bends: List[Tuple[Tuple[float, float], bool]] = []
        if explicit_bend is not None:
            pt_b = _extract_xy_point(explicit_bend, ref_pt=primary_crown)
            if pt_b is not None:
                cand_bends.append((pt_b, True))

        raw_cauline = (
            upper_cauline_detections
            or kwargs.get("cauline_detections")
            or kwargs.get("cauline_boxes")
            or kwargs.get("upper_cauline_pts")
            or kwargs.get("stem_points")
        )
        if raw_cauline:
            for c_item in raw_cauline:
                pt_b = _extract_xy_point(c_item, ref_pt=primary_crown)
                if pt_b is not None:
                    cand_bends.append((pt_b, False))

        skel_mask = stem_skeleton_mask if stem_skeleton_mask is not None else kwargs.get("specimen_mask")
        if isinstance(skel_mask, np.ndarray) and np.count_nonzero(skel_mask) > 0:
            ys_s, xs_s = np.where(skel_mask > 0)
            # Uppermost apex of the stem skeleton
            min_y_idx = int(np.argmin(ys_s))
            top_y = float(ys_s[min_y_idx])
            top_x = float(np.median(xs_s[ys_s <= top_y + 3]))
            cand_bends.append(((top_x, top_y), False))
            # Point of maximum lateral excursion from chord p_crown -> p_corymb
            if straight_px > 1e-4:
                chord_u = np.array(
                    [resolved_corymb[0] - primary_crown[0], resolved_corymb[1] - primary_crown[1]],
                    dtype=np.float64,
                ) / straight_px
                pts_rel = np.column_stack([xs_s - primary_crown[0], ys_s - primary_crown[1]])
                perp_d = np.abs(pts_rel[:, 0] * (-chord_u[1]) + pts_rel[:, 1] * chord_u[0])
                max_p_idx = int(np.argmax(perp_d))
                cand_bends.append(((float(xs_s[max_p_idx]), float(ys_s[max_p_idx])), False))

        best_piecewise_px = straight_px
        for bend_pt, is_explicit in cand_bends:
            d1 = float(math.hypot(bend_pt[0] - primary_crown[0], bend_pt[1] - primary_crown[1]))
            d2 = float(math.hypot(resolved_corymb[0] - bend_pt[0], resolved_corymb[1] - bend_pt[1]))
            pw_len = d1 + d2
            if straight_px > 1e-4:
                ux = (resolved_corymb[0] - primary_crown[0]) / straight_px
                uy = (resolved_corymb[1] - primary_crown[1]) / straight_px
                perp_dist = abs((bend_pt[0] - primary_crown[0]) * (-uy) + (bend_pt[1] - primary_crown[1]) * ux)
            else:
                perp_dist = d1

            # Inverted-V fold (bend apex higher on sheet than corymb) or significant lateral stem deflection
            is_upper_fold = bend_pt[1] < (min(primary_crown[1], resolved_corymb[1]) - 8.0)
            is_lateral_bend = (
                bend_pt[1] < primary_crown[1]
                and perp_dist >= max(15.0, 0.12 * straight_px)
                and pw_len >= 1.04 * straight_px
            )
            if (is_explicit and (pw_len > straight_px + 1.0 or perp_dist > 2.0)) or is_upper_fold or is_lateral_bend:
                if pw_len >= best_piecewise_px:
                    best_piecewise_px = pw_len
                    resolved_bend = bend_pt
                    is_bent_stem = True

        scape_length_px = best_piecewise_px if is_bent_stem else straight_px
        if eff_px_per_cm is not None and eff_px_per_cm > 0.0:
            plant_height_cm = round(float(scape_length_px / eff_px_per_cm), 3)

        # Basal leaf reference length for scape_to_basal_leaf_ratio
        if np.isfinite(petiole_px_p75) and np.isfinite(blade_px_p75) and (petiole_px_p75 + blade_px_p75) > 1e-6:
            basal_ref_px = float(petiole_px_p75 + blade_px_p75)
        elif np.isfinite(blade_px_p75) and blade_px_p75 > 1e-6:
            basal_ref_px = float(blade_px_p75)
        else:
            basal_ref_px = np.nan

        if np.isfinite(basal_ref_px) and basal_ref_px > 1e-6:
            scape_to_basal_leaf_ratio = round(float(scape_length_px / basal_ref_px), 4)

    return {
        "p_crown": primary_crown,
        "crowns_by_plant": crowns_by_pid,
        "p_corymb": resolved_corymb,
        "p_bend": resolved_bend,
        "is_bent_stem": bool(is_bent_stem),
        "scape_length_px": round(float(scape_length_px), 4) if np.isfinite(scape_length_px) else np.nan,
        "plant_height_cm": plant_height_cm,
        "scape_to_basal_leaf_ratio": scape_to_basal_leaf_ratio,
        "petiole_length_px": petiole_px_p75,
        "blade_length_px": blade_px_p75,
        "petiole_length_mm": petiole_mm_p75,
        "blade_length_mm": blade_mm_p75,
        "petiole_to_blade_ratio": pet_to_blade_p75,
        "leaf_metrics": leaf_metrics,
        "by_plant": by_plant_summary,
    }

