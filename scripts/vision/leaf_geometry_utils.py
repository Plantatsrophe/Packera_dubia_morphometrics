#!/usr/bin/env python3
"""
===============================================================================
Module: leaf_geometry_utils.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Consolidated LeafMachine3 (lm3) geometric quality gatekeeping, 31-keypoint
    midrib/petiole pose estimation, 3-class (lamina/petiole/hole) mask routing,
    bilateral hemi-blade reflection, adaptive DBSCAN spatial clustering,
    and standardized contour coordinate extraction.

Key Capabilities:
    1. Geometric Quality Gatekeeping (LeafMachine3 YOLO26x & UNet++ v9):
       - Unoccluded Completeness Score (UCS >= 0.85): Area_mask / (0.70 * Area_min_rect)
       - Solidity (Solidity >= 0.72): Area_mask / Area_convex_hull
       - Longitudinal midrib axis angle and apex/base pose estimation anchored by
         LM3 31-keypoint anatomical landmarks (lamina_base, lamina_tip, midvein_1..15)
    2. Vectorized Hemi-Blade Reflection (Tier 2):
       - Primary midrib alignment via affine rotation or LM3 15-point midvein trace
       - Vectorized margin profile defect calculation across upper/lower halves
       - Fully vectorized slice reflection in NumPy (zero manual pixel loops)
       - Morphological seam closure
    3. Adaptive DBSCAN Spatial Clustering:
       - Groups detected organs into plant_individual_id clusters to prevent
         trait averaging across multi-plant specimen sheets.
    4. Standardized Contour & Caliper Extraction:
       - Resamples contour boundary to standardized normalized (x, y) coordinates
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

def _endpoint_local_dt_radius(
    dt_map: np.ndarray,
    ep_start: Tuple[int, int],
    ep_other: Tuple[int, int],
) -> float:
    """Estimates local inscribed radius near `ep_start` along the terminal 20% of the axis."""
    h, w = dt_map.shape[:2]
    axis_len = math.hypot(ep_other[0] - ep_start[0], ep_other[1] - ep_start[1])
    if axis_len < 4.0:
        ix = max(0, min(w - 1, int(ep_start[0])))
        iy = max(0, min(h - 1, int(ep_start[1])))
        return float(dt_map[iy, ix])
    win_r = max(2, min(20, int(round(0.08 * axis_len))))
    radii: List[float] = []
    for frac in (0.04, 0.08, 0.12, 0.16, 0.20):
        px = int(round(ep_start[0] + frac * (ep_other[0] - ep_start[0])))
        py = int(round(ep_start[1] + frac * (ep_other[1] - ep_start[1])))
        x0 = max(0, px - win_r)
        x1 = min(w, px + win_r + 1)
        y0 = max(0, py - win_r)
        y1 = min(h, py + win_r + 1)
        patch = dt_map[y0:y1, x0:x1]
        if patch.size > 0:
            radii.append(float(np.max(patch)))
    return float(np.mean(radii)) if radii else 0.0


def find_lamina_midrib_endpoints(
    mask: np.ndarray,
    rosette_centroid: Optional[Tuple[float, float]] = None,
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """
    Finds the longitudinal midrib endpoints `(p_apex, p_base)` for a binary leaf mask.
    If `rosette_centroid` is provided, the endpoint closer to `rosette_centroid` is
    `p_base` (and the farther is `p_apex`); otherwise disambiguates using distance-transform
    local radius near the two endpoints (the narrower endpoint is `p_base`, or higher `y`
    in tip-up coordinates if radii are tied).
    """
    _, _, _, p_apex, p_base = compute_geometric_metrics(mask, rosette_centroid=rosette_centroid)
    return p_apex, p_base


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

    # Disambiguate p_base vs p_apex using rosette_centroid or distance-transform local radius
    if rc_sub is not None:
        d1 = math.hypot(ep1_sub[0] - rc_sub[0], ep1_sub[1] - rc_sub[1])
        d2 = math.hypot(ep2_sub[0] - rc_sub[0], ep2_sub[1] - rc_sub[1])
        if abs(d1 - d2) > 1e-3:
            p_apex_sub, p_base_sub = (ep2_sub, ep1_sub) if d1 < d2 else (ep1_sub, ep2_sub)
        else:
            p_apex_sub, p_base_sub = (
                (ep2_sub, ep1_sub) if (ep1_sub[1], ep1_sub[0]) > (ep2_sub[1], ep2_sub[0]) else (ep1_sub, ep2_sub)
            )
    else:
        dt_sub = cv2.distanceTransform((sub_mask > 0).astype(np.uint8), cv2.DIST_L2, 5)
        r1 = _endpoint_local_dt_radius(dt_sub, ep1_sub, ep2_sub)
        r2 = _endpoint_local_dt_radius(dt_sub, ep2_sub, ep1_sub)
        if abs(r1 - r2) > 0.35:
            p_apex_sub, p_base_sub = (ep2_sub, ep1_sub) if r1 < r2 else (ep1_sub, ep2_sub)
        elif (ep1_sub[1], ep1_sub[0]) > (ep2_sub[1], ep2_sub[0]):
            p_apex_sub, p_base_sub = ep2_sub, ep1_sub
        else:
            p_apex_sub, p_base_sub = ep1_sub, ep2_sub

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


def _prepare_lamina_mask(
    mask: np.ndarray,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Applies LM3 `petiole_mask` subtraction and topology-aware `hole_mask` handling:
      - Subtracts `petiole_mask` (`class_id == 1`) from `work_mask` when present.
      - Differentiates true interior lamina herbivory holes from enclosed inter-lobe
        background sinuses on pinnatifid/dissected leaves:
          * True interior holes (compact cavities > 6 px inside the outer contour and
            < 6% of lamina area) are filled so internal insect damage does not distort EFA.
          * Enclosed inter-lobe sinuses (cavities within <= 6 px of the exterior background
            formed when adjacent pinna tips touch, or large background voids >= 6% of lamina
            area) are subtracted and their narrow contact bridge to the exterior background
            is opened so `cv2.RETR_EXTERNAL` traces into the sinus around both pinna lobes.
    """
    work_mask = (mask > 0).astype(np.uint8) * 255
    if (
        petiole_mask is not None
        and isinstance(petiole_mask, np.ndarray)
        and petiole_mask.shape[:2] == work_mask.shape[:2]
        and np.count_nonzero(petiole_mask) > 0
    ):
        cand = np.where((work_mask > 0) & (petiole_mask == 0), 255, 0).astype(np.uint8)
        if np.count_nonzero(cand) >= max(10, int(0.35 * np.count_nonzero(work_mask))):
            work_mask = cand

    lamina_fg = int(np.count_nonzero(work_mask))
    if lamina_fg <= 0:
        return work_mask

    # Identify exterior background from outer contour fill
    ext_cnts, _ = cv2.findContours(work_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not ext_cnts:
        return work_mask
    filled_ext = np.zeros_like(work_mask)
    cv2.drawContours(filled_ext, ext_cnts, -1, 255, thickness=cv2.FILLED)

    # Combine explicit hole_mask with any enclosed cavities inside filled_ext
    internal_voids = ((filled_ext > 0) & (work_mask == 0)).astype(np.uint8) * 255
    if (
        hole_mask is not None
        and isinstance(hole_mask, np.ndarray)
        and hole_mask.shape[:2] == work_mask.shape[:2]
        and np.count_nonzero(hole_mask) > 0
    ):
        internal_voids = np.maximum(internal_voids, ((hole_mask > 0) & (filled_ext > 0)).astype(np.uint8) * 255)

    if int(np.count_nonzero(internal_voids)) == 0:
        return work_mask

    # Distance from each pixel inside filled_ext to the exterior background (filled_ext == 0)
    dist_to_ext = cv2.distanceTransform((filled_ext > 0).astype(np.uint8), cv2.DIST_L2, 3)
    ext_bg_mask = (filled_ext == 0).astype(np.uint8)
    _, ext_labels = cv2.distanceTransformWithLabels(
        (filled_ext > 0).astype(np.uint8), cv2.DIST_L2, 3, labelType=cv2.DIST_LABEL_PIXEL
    )
    ext_coords_y, ext_coords_x = np.where(ext_bg_mask > 0)

    n_voids, void_lbls, void_stats, _ = cv2.connectedComponentsWithStats(
        (internal_voids > 0).astype(np.uint8), connectivity=8
    )
    has_hole_arr = (
        hole_mask is not None
        and isinstance(hole_mask, np.ndarray)
        and hole_mask.shape[:2] == work_mask.shape[:2]
        and np.count_nonzero(hole_mask) > 0
    )
    hole_bool_arr = (hole_mask > 0) if has_hole_arr else None
    max_bridge_px = 24.0
    for v_id in range(1, n_voids):
        v_mask = void_lbls == v_id
        v_area = int(void_stats[v_id, cv2.CC_STAT_AREA])
        if v_area <= 0:
            continue
        min_d_ext = float(np.min(dist_to_ext[v_mask])) if np.any(v_mask) else 999.0
        area_frac = float(v_area) / float(max(lamina_fg, 1))
        hole_overlap = (
            float(np.count_nonzero(v_mask & hole_bool_arr)) / float(max(v_area, 1))
            if hole_bool_arr is not None
            else 0.0
        )
        # Enclosed inter-lobe sinus / background paper pocket: separated from exterior background
        # by a lobe bridge (<= 24 px) with area >= 25 px, or a large enclosed background region.
        is_enclosed_sinus = (
            (v_area >= 25 and 0.0 < min_d_ext <= max_bridge_px and (area_frac >= 0.008 or hole_overlap > 0.15))
            or (min_d_ext <= max_bridge_px and area_frac >= 0.05)
            or (area_frac >= 0.08)
        )
        if is_enclosed_sinus:
            work_mask[v_mask] = 0
            # Carve open the narrow contact bottleneck to exterior background so RETR_EXTERNAL enters the sinus
            v_ys, v_xs = np.where(v_mask)
            if len(v_ys) > 0 and len(ext_coords_y) > 0:
                d_vals = dist_to_ext[v_ys, v_xs]
                best_i = int(np.argmin(d_vals))
                vy0, vx0 = int(v_ys[best_i]), int(v_xs[best_i])
                # Find nearest exterior background pixel within a local window
                win = max(8, int(math.ceil(min_d_ext)) + 4)
                y_lo, y_hi = max(0, vy0 - win), min(work_mask.shape[0], vy0 + win + 1)
                x_lo, x_hi = max(0, vx0 - win), min(work_mask.shape[1], vx0 + win + 1)
                bg_sub_y, bg_sub_x = np.where(filled_ext[y_lo:y_hi, x_lo:x_hi] == 0)
                if len(bg_sub_y) > 0:
                    bg_y = bg_sub_y + y_lo
                    bg_x = bg_sub_x + x_lo
                    nearest_bg_idx = int(np.argmin((bg_y - vy0) ** 2 + (bg_x - vx0) ** 2))
                    ey0, ex0 = int(bg_y[nearest_bg_idx]), int(bg_x[nearest_bg_idx])
                    cv2.line(work_mask, (vx0, vy0), (ex0, ey0), 0, thickness=3)
        else:
            # True interior insect hole -> fill in work_mask
            work_mask[v_mask] = 255

    return work_mask


def _rotate_contour_to_tip_up(
    resamp_pts: np.ndarray,
    anchor_pt: Optional[Tuple[float, float]] = None,
    p_apex_pt: Optional[Tuple[float, float]] = None,
    keypoints_31: Optional[Any] = None,
) -> np.ndarray:
    """
    Rotates a resampled clockwise contour `(N, 2)` around its centroid by LeafMachine3's
    clockwise orientation angle `angle_cw` (`compute_orientation` or `_angle_cw_to_up`)
    so the anatomical `lamina_base -> lamina_tip` axis points strictly tip-up `(0, -1)`
    in image coordinates (`y` down). Leaves already oriented tip-up (`|angle_cw| <= 0.25 deg`)
    untouched.
    """
    if resamp_pts is None or len(resamp_pts) < 4:
        return resamp_pts

    angle_cw: Optional[float] = None
    if keypoints_31 is not None:
        try:
            from scripts.vision.lm3_inference_engines import normalize_lm3_keypoints_dict
            from leafmachine3.core.orientation import compute_orientation

            norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
            pts_2d = {
                k: (float(v[0]), float(v[1]))
                for k, v in norm_kps.items()
                if len(v) >= 3 and float(v[2]) >= 0.25
            }
            ori = compute_orientation(pts_2d)
            if ori.success and ori.angle_cw is not None:
                angle_cw = float(ori.angle_cw)
        except Exception:
            angle_cw = None

    if angle_cw is None:
        base_xy = (
            np.array([float(anchor_pt[0]), float(anchor_pt[1])], dtype=np.float64)
            if anchor_pt is not None
            else np.asarray(resamp_pts[0], dtype=np.float64)
        )
        if p_apex_pt is not None:
            apex_xy = np.array([float(p_apex_pt[0]), float(p_apex_pt[1])], dtype=np.float64)
        else:
            mid_idx = len(resamp_pts) // 2
            apex_xy = np.asarray(resamp_pts[mid_idx], dtype=np.float64)
        dx = float(apex_xy[0] - base_xy[0])
        dy = float(apex_xy[1] - base_xy[1])
        if math.hypot(dx, dy) > 1e-3:
            angle_cw = -math.degrees(math.atan2(dx, -dy))
        else:
            angle_cw = 0.0

    if angle_cw is not None and abs(float(angle_cw)) > 0.25:
        rad = math.radians(float(angle_cw))
        cos_a = math.cos(rad)
        sin_a = math.sin(rad)
        centroid = np.mean(resamp_pts, axis=0)
        centered = resamp_pts - centroid
        x_rot = centered[:, 0] * cos_a - centered[:, 1] * sin_a + centroid[0]
        y_rot = centered[:, 0] * sin_a + centered[:, 1] * cos_a + centroid[1]
        resamp_pts = np.column_stack([x_rot, y_rot])

    return resamp_pts


def resample_normalized_contour(
    mask: Optional[np.ndarray],
    num_points: int = 120,
    petiole_attachment_pt: Optional[Tuple[float, float]] = None,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
    lamina_base: Optional[Tuple[float, float]] = None,
    momocs_outline_xy: Optional[np.ndarray] = None,
    p_apex_pt: Optional[Tuple[float, float]] = None,
    keypoints_31: Optional[Any] = None,
) -> Optional[np.ndarray]:
    """
    Extracts resampled, normalized 2D coordinates (x_norm, y_norm) of shape (num_points, 2)
    from a binary mask or outline coordinate array, anchored at vertex 0 (`petiole_attachment_pt`,
    `lamina_base`, or basal-most point `argmax(y)`), oriented tip-up via LeafMachine3 native
    orientation, with clockwise winding and isotropic normalization by `max(span_x, span_y, 1e-6)`.
    """
    if momocs_outline_xy is not None and len(np.asarray(momocs_outline_xy)) >= 10:
        largest_cnt = np.asarray(momocs_outline_xy, dtype=np.float64).reshape(-1, 2)
    else:
        if mask is None or mask.size == 0:
            return None
        work_mask = _prepare_lamina_mask(mask, petiole_mask=petiole_mask, hole_mask=hole_mask)
        contours, _ = cv2.findContours(work_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            return None
        largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
        if len(largest_cnt) < 10:
            return None

    anchor_pt = petiole_attachment_pt if petiole_attachment_pt is not None else lamina_base
    if anchor_pt is not None:
        diffs = largest_cnt - np.array([float(anchor_pt[0]), float(anchor_pt[1])], dtype=np.float64)
        start_idx = int(np.argmin(np.hypot(diffs[:, 0], diffs[:, 1])))
    else:
        start_idx = int(np.argmax(largest_cnt[:, 1]))

    pts = np.roll(largest_cnt, -start_idx, axis=0)
    signed_area = cv2.contourArea(pts.reshape(-1, 1, 2).astype(np.float32), oriented=True)
    if signed_area > 0:
        pts = np.vstack([pts[0:1], pts[1:][::-1]])

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

    resamp_pts = np.column_stack([x_resamp, y_resamp])
    signed_area_r = cv2.contourArea(resamp_pts.reshape(-1, 1, 2).astype(np.float32), oriented=True)
    if signed_area_r > 0:
        resamp_pts = np.vstack([resamp_pts[0:1], resamp_pts[1:][::-1]])
    resamp_pts = _rotate_contour_to_tip_up(
        resamp_pts,
        anchor_pt=anchor_pt,
        p_apex_pt=p_apex_pt,
        keypoints_31=keypoints_31,
    )
    x_resamp = resamp_pts[:, 0]
    y_resamp = resamp_pts[:, 1]

    min_x, max_x = float(np.min(x_resamp)), float(np.max(x_resamp))
    min_y, max_y = float(np.min(y_resamp)), float(np.max(y_resamp))
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)
    max_span = max(span_x, span_y, 1e-6)

    x_norm = np.round((x_resamp - min_x) / max_span, 6)
    y_norm = np.round((y_resamp - min_y) / max_span, 6)
    return np.column_stack([x_norm, y_norm])


# =============================================================================
# 2. Vectorized Hemi-Blade Reflection (Tier 2) & Tier 1 Extraction
# =============================================================================

def extract_tier1_pristine(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    output_dir: Path,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
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

    bin_mask = _prepare_lamina_mask(mask, petiole_mask=petiole_mask, hole_mask=hole_mask)
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
    min_ucs: float = 0.50,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
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

    mask_u8 = _prepare_lamina_mask(mask, petiole_mask=petiole_mask, hole_mask=hole_mask)
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
    mask: Optional[np.ndarray],
    catalog_number: str,
    leaf_id: int,
    output_dir: Path,
    num_points: int = 120,
    petiole_attachment_pt: Optional[Tuple[float, float]] = None,
    contour: Optional[np.ndarray] = None,
    p_apex_pt: Optional[Tuple[float, float]] = None,
    return_coords: bool = False,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
    lamina_base: Optional[Tuple[float, float]] = None,
    momocs_outline_xy: Optional[np.ndarray] = None,
    keypoints_31: Optional[Any] = None,
) -> Any:
    """
    Extracts vectorized 2D (x, y) boundary coordinates from a binary silhouette mask,
    anchors vertex 0 to `petiole_attachment_pt` / `lamina_base` (or the basal-most
    point `argmax(y)`), enforces clockwise orientation, rotates the contour into
    LeafMachine3 canonical tip-up orientation (`lamina_base -> lamina_tip` along `(0, -1)`),
    normalizes coordinates `(x_norm, y_norm)` isotropically by `max(span_x, span_y, 1e-6)`,
    and saves to `output_dir/contours/{catalogNumber}_leaf{id}.csv`.
    """
    empty_ret: Any = (None, None) if return_coords else None
    work_mask: Optional[np.ndarray] = None
    if momocs_outline_xy is not None and len(np.asarray(momocs_outline_xy)) >= 10:
        largest_cnt = np.asarray(momocs_outline_xy, dtype=np.float64).reshape(-1, 2)
    elif contour is not None and len(contour) >= 10:
        largest_cnt = np.asarray(contour, dtype=np.float64).reshape(-1, 2)
    else:
        if mask is None or mask.size == 0:
            return empty_ret
        work_mask = _prepare_lamina_mask(mask, petiole_mask=petiole_mask, hole_mask=hole_mask)
        contours, _ = cv2.findContours(work_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            return empty_ret
        largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
        if len(largest_cnt) < 10:
            return empty_ret

    anchor_pt = petiole_attachment_pt if petiole_attachment_pt is not None else lamina_base
    if keypoints_31 is not None:
        if isinstance(keypoints_31, dict):
            lb_kp = keypoints_31.get("lamina_base") or keypoints_31.get("base_2") or keypoints_31.get("base_center")
            lt_kp = keypoints_31.get("lamina_tip") or keypoints_31.get("apex_2") or keypoints_31.get("apex_center")
            if anchor_pt is None and lb_kp is not None and len(lb_kp) >= 2:
                c1 = float(lb_kp[2]) if len(lb_kp) >= 3 else 1.0
                if c1 >= 0.25 and (abs(float(lb_kp[0])) + abs(float(lb_kp[1]))) > 0.0:
                    anchor_pt = (float(lb_kp[0]), float(lb_kp[1]))
            if p_apex_pt is None and lt_kp is not None and len(lt_kp) >= 2:
                c0 = float(lt_kp[2]) if len(lt_kp) >= 3 else 1.0
                if c0 >= 0.25 and (abs(float(lt_kp[0])) + abs(float(lt_kp[1]))) > 0.0:
                    p_apex_pt = (float(lt_kp[0]), float(lt_kp[1]))
        else:
            kp_arr = np.asarray(keypoints_31, dtype=np.float64)
            if kp_arr.ndim == 2 and kp_arr.shape[0] >= 2 and kp_arr.shape[1] >= 2:
                c0 = float(kp_arr[0, 2]) if kp_arr.shape[1] >= 3 else 1.0
                c1 = float(kp_arr[1, 2]) if kp_arr.shape[1] >= 3 else 1.0
                if anchor_pt is None and c1 >= 0.25 and (abs(kp_arr[1, 0]) + abs(kp_arr[1, 1])) > 0.0:
                    anchor_pt = (float(kp_arr[1, 0]), float(kp_arr[1, 1]))
                if p_apex_pt is None and c0 >= 0.25 and (abs(kp_arr[0, 0]) + abs(kp_arr[0, 1])) > 0.0:
                    p_apex_pt = (float(kp_arr[0, 0]), float(kp_arr[0, 1]))
    if anchor_pt is None and work_mask is not None:
        try:
            _, _, _, calc_apex, calc_base = compute_geometric_metrics(work_mask)
            anchor_pt = calc_base
            if p_apex_pt is None:
                p_apex_pt = calc_apex
        except Exception:
            pass
    if anchor_pt is not None:
        diffs = largest_cnt - np.array([float(anchor_pt[0]), float(anchor_pt[1])], dtype=np.float64)
        start_idx = int(np.argmin(np.hypot(diffs[:, 0], diffs[:, 1])))
    else:
        start_idx = int(np.argmax(largest_cnt[:, 1]))

    pts = np.roll(largest_cnt, -start_idx, axis=0)
    signed_area = cv2.contourArea(pts.reshape(-1, 1, 2).astype(np.float32), oriented=True)
    if signed_area > 0:
        pts = np.vstack([pts[0:1], pts[1:][::-1]])

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

    resamp_pts = np.column_stack([x_resamp, y_resamp])
    signed_area_r = cv2.contourArea(resamp_pts.reshape(-1, 1, 2).astype(np.float32), oriented=True)
    if signed_area_r > 0:
        resamp_pts = np.vstack([resamp_pts[0:1], resamp_pts[1:][::-1]])
    resamp_pts = _rotate_contour_to_tip_up(
        resamp_pts,
        anchor_pt=anchor_pt,
        p_apex_pt=p_apex_pt,
        keypoints_31=keypoints_31,
    )
    x_resamp = resamp_pts[:, 0]
    y_resamp = resamp_pts[:, 1]

    min_x, max_x = float(np.min(x_resamp)), float(np.max(x_resamp))
    min_y, max_y = float(np.min(y_resamp)), float(np.max(y_resamp))
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)
    max_span = max(span_x, span_y, 1e-6)

    x_norm = np.round((x_resamp - min_x) / max_span, 6)
    y_norm = np.round((y_resamp - min_y) / max_span, 6)

    x_r = np.round(x_resamp, 2).tolist()
    y_r = np.round(y_resamp, 2).tolist()
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
    _extract_lm3_midvein_trace,
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
    _estimate_single_plant_crown,
    _extract_xy_point,
    _parse_basal_leaf_pose,
    approximate_petiole_and_scape_metrics,
    cluster_plant_individuals,
    cluster_voucher_plants_dbscan,
    estimate_rosette_crown_point,
    is_topologically_connected_mask,
    route_tier3_open_curve,
)

__all__ = [
    "find_lamina_midrib_endpoints",
    "compute_geometric_metrics",
    "compute_geometric_metrics_and_pose",
    "is_tier1_pristine",
    "BoundaryContinuityResult",
    "_extract_lm3_midvein_trace",
    "_has_artificial_straight_occlusion_cut",
    "check_boundary_tissue_continuity",
    "compute_classical_shape_descriptors",
    "reflect_folded_mask",
    "FoldDetectionResult",
    "detect_folded_leaf",
    "bridge_dissected_rachis_components",
    "is_botanical_dissection",
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
    "is_topologically_connected_mask",
    "cluster_plant_individuals",
    "cluster_voucher_plants_dbscan",
    "route_tier3_open_curve",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
]
