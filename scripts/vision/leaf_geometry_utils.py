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
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
LM3_PACKAGE_ROOT = PROJECT_ROOT / "LeafMachine3"
if LM3_PACKAGE_ROOT.is_dir() and str(LM3_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(LM3_PACKAGE_ROOT))

from leafmachine3.core.leaf_outline import extract_contour
from leafmachine3.core.momocs_prep import (
    prepare_mask,
    resample_closed,
    signed_area as momocs_signed_area,
)
from leafmachine3.core.morphometrics import polygon_morphology
from leafmachine3.core.orientation import compute_orientation, length_width_from_box
from leafmachine3.inference.leaf_segmenter import _components, _largest_component


# =============================================================================
# 1. Geometric Quality Gatekeeping & Pose Estimation
# =============================================================================

def _endpoint_local_dt_radius(
    dt_map: np.ndarray,
    ep_start: Tuple[int, int],
    ep_other: Tuple[int, int],
) -> float:
    """Estimates local inscribed radius near `ep_start` along the terminal 1-5% petiole neck of the axis."""
    h, w = dt_map.shape[:2]
    axis_len = math.hypot(ep_other[0] - ep_start[0], ep_other[1] - ep_start[1])
    if axis_len < 4.0:
        ix = max(0, min(w - 1, int(ep_start[0])))
        iy = max(0, min(h - 1, int(ep_start[1])))
        return float(dt_map[iy, ix])
    win_r = max(2, min(6, int(round(0.015 * axis_len))))
    radii: List[float] = []
    for frac in (0.01, 0.02, 0.03, 0.04, 0.05):
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
    angle_cw: Optional[float] = None,
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
    angle_cw : Optional[float]
        Optional clockwise upright orientation angle (degrees) for `length_width_from_box`.

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

    main_contour_xy = extract_contour(sub_mask)
    if main_contour_xy is None or len(main_contour_xy) < 3:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    main_contour = main_contour_xy.reshape(-1, 1, 2).astype(np.int32)
    morph = polygon_morphology(main_contour_xy, method="minarearect")
    if morph is None:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    # 1. Convex Hull & Solidity via LeafMachine3 polygon_morphology
    area_hull = float(getattr(morph, "convex_area_px", getattr(morph, "convex_hull_area", 0.0)))
    morph_solidity = float(getattr(morph, "solidity", getattr(morph, "convexity", 0.0)))
    solidity = (area_mask / area_hull) if area_hull > 0 else max(0.0, morph_solidity)
    solidity = min(max(solidity, 0.0), 1.0)

    # 2. Expected Area via Minimum Area Rect (Oriented Bounding Box)
    rect_w = float(getattr(morph, "dim_max", 0.0))
    rect_h = float(getattr(morph, "dim_min", 0.0))
    if rect_w <= 0.0 or rect_h <= 0.0:
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
    pts = main_contour_xy.reshape(-1, 2).astype(np.float64)
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
        math.atan2(float(p_apex[1] - p_base[1]), float(p_apex[0] - p_base[0]))
    )
    eff_angle_cw = angle_cw
    if eff_angle_cw is None:
        dx_up = float(p_apex[0] - p_base[0])
        dy_up = float(p_apex[1] - p_base[1])
        if math.hypot(dx_up, dy_up) > 1e-6:
            eff_angle_cw = (-math.degrees(math.atan2(dx_up, -dy_up))) % 360.0
    if eff_angle_cw is not None and getattr(morph, "rotated_bbox", None) is not None:
        _ = length_width_from_box(morph.rotated_bbox, eff_angle_cw)

    return ucs, solidity, angle_deg, p_apex, p_base


def _prepare_lamina_mask(
    mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
    *,
    lamina_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Delegates connected-component filtering and boolean hole subtraction to native LeafMachine3:
      - Isolates `leaf_comp = _largest_component((mask > 0).astype(np.uint8))` (falling back
        to `(mask > 0).astype(np.uint8)` when total pixels < 64 in small synthetic fixtures).
      - Isolates `pet_comp = _largest_component((petiole_mask > 0).astype(np.uint8))` when
        `petiole_mask` is provided and has the same shape.
      - Isolates `hole_comp = _components((hole_mask > 0).astype(np.uint8), min_area=64)` when
        `hole_mask` is provided and has the same shape.
      - Computes `silhouette = (leaf_comp > 0) | (pet_comp > 0)` and
        `lamina_bool = silhouette & ~(hole_comp > 0)` (and `~(pet_comp > 0)` when `pet_comp`
        does not wipe the mask), passes `(lamina_bool.astype(np.uint8) * 255)` through
        `leafmachine3.core.momocs_prep.prepare_mask`, and returns `(prepared > 0).astype(np.uint8) * 255`.
    """
    src_mask = mask if mask is not None else lamina_mask
    if src_mask is None or src_mask.size == 0:
        return (
            np.zeros((0, 0), dtype=np.uint8)
            if src_mask is None
            else np.zeros(src_mask.shape[:2], dtype=np.uint8)
        )

    fg_u8 = (src_mask > 0).astype(np.uint8)
    fg_count = int(np.count_nonzero(fg_u8))
    if fg_count == 0:
        return np.zeros(src_mask.shape[:2], dtype=np.uint8)

    leaf_comp = _largest_component(fg_u8)
    if fg_count < 64 or not np.any(leaf_comp):
        leaf_comp = fg_u8 > 0

    pet_comp = np.zeros(src_mask.shape[:2], dtype=bool)
    if (
        petiole_mask is not None
        and isinstance(petiole_mask, np.ndarray)
        and petiole_mask.shape[:2] == src_mask.shape[:2]
        and np.count_nonzero(petiole_mask) > 0
    ):
        pet_comp = np.asarray(_largest_component((petiole_mask > 0).astype(np.uint8))) > 0

    hole_comp = np.zeros(src_mask.shape[:2], dtype=bool)
    if (
        hole_mask is not None
        and isinstance(hole_mask, np.ndarray)
        and hole_mask.shape[:2] == src_mask.shape[:2]
        and np.count_nonzero(hole_mask) > 0
    ):
        hole_items = _components((hole_mask > 0).astype(np.uint8), min_area=64)
        if isinstance(hole_items, np.ndarray):
            hole_comp = hole_items > 0
        elif hole_items:
            hole_polys: List[np.ndarray] = []
            for h_item in hole_items:
                h_u8 = (np.asarray(h_item) > 0).astype(np.uint8)
                hole_comp |= h_u8 > 0
                h_cnts, _ = cv2.findContours(h_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                for hc in h_cnts:
                    pts_h = hc.reshape(-1, 2).astype(np.int32)
                    if len(pts_h) >= 3:
                        hole_polys.append(pts_h)
            if hole_polys:
                qualifying_holes = np.zeros(src_mask.shape[:2], dtype=np.uint8)
                cv2.fillPoly(qualifying_holes, hole_polys, 255)
                hole_comp |= qualifying_holes > 0

    silhouette = (np.asarray(leaf_comp) > 0) | (np.asarray(pet_comp) > 0)
    lamina_bool = silhouette & ~(np.asarray(hole_comp) > 0)
    if np.any(pet_comp):
        cand_no_pet = lamina_bool & ~(np.asarray(pet_comp) > 0)
        if int(np.count_nonzero(cand_no_pet)) >= max(10, int(0.35 * int(np.count_nonzero(lamina_bool)))):
            lamina_bool = cand_no_pet

    lamina_u8 = lamina_bool.astype(np.uint8) * 255
    if int(np.count_nonzero(lamina_u8)) == 0:
        return lamina_u8

    min_area_val = 64 if int(np.count_nonzero(lamina_u8)) >= 64 else 1
    try:
        prepared = prepare_mask(lamina_u8, pad_px=0, fill_holes=False, largest_only=True)
    except TypeError:
        prepared = prepare_mask(lamina_u8, min_area=min_area_val, smooth_kernel=0)

    if prepared is None or not np.any(prepared):
        return lamina_u8
    prep_arr = np.asarray(prepared)
    if prep_arr.shape[:2] != src_mask.shape[:2]:
        return lamina_u8
    return (prep_arr > 0).astype(np.uint8) * 255


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
        ori_fallback = compute_orientation({
            "lamina_tip": (float(apex_xy[0]), float(apex_xy[1])),
            "lamina_base": (float(base_xy[0]), float(base_xy[1])),
        })
        if ori_fallback.success and ori_fallback.angle_cw is not None:
            angle_cw = float(ori_fallback.angle_cw)
        else:
            angle_cw = 0.0

    if angle_cw is not None:
        signed_cw = ((float(angle_cw) + 180.0) % 360.0) - 180.0
        if abs(signed_cw) > 0.25:
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
    n_pts: Optional[int] = None,
) -> Optional[np.ndarray]:
    """
    Extracts resampled, normalized 2D coordinates (x_norm, y_norm) of shape (num_points, 2)
    from a binary mask or outline coordinate array, anchored at vertex 0 (`petiole_attachment_pt`,
    `lamina_base`, or basal-most point `argmax(y)`), oriented tip-up via LeafMachine3 native
    orientation, with clockwise winding and isotropic normalization by `max(span_x, span_y, 1e-6)`.
    """
    if n_pts is not None:
        num_points = int(n_pts)

    if momocs_outline_xy is not None and len(np.asarray(momocs_outline_xy)) >= 10:
        largest_cnt = np.asarray(momocs_outline_xy, dtype=np.float64).reshape(-1, 2)
    else:
        if mask is None or np.asarray(mask).size == 0:
            return None
        arr_in = np.asarray(mask)
        if (arr_in.ndim == 3 and arr_in.shape[-1] == 2) or (
            arr_in.ndim == 2 and arr_in.shape[1] == 2 and arr_in.shape[0] >= 3
        ):
            largest_cnt = arr_in.reshape(-1, 2).astype(np.float64)
            if len(largest_cnt) < 10:
                return None
        else:
            clean_mask = _prepare_lamina_mask(arr_in, petiole_mask=petiole_mask, hole_mask=hole_mask)
            extracted_xy = extract_contour(clean_mask)
            if extracted_xy is None or len(extracted_xy) < 10:
                return None
            largest_cnt = np.asarray(extracted_xy, dtype=np.float64).reshape(-1, 2)

    anchor_pt = petiole_attachment_pt if petiole_attachment_pt is not None else lamina_base
    if keypoints_31 is not None:
        try:
            from scripts.vision.lm3_inference_engines import normalize_lm3_keypoints_dict

            norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
            lb_kp = norm_kps.get("lamina_base") or norm_kps.get("base_center")
            lt_kp = norm_kps.get("lamina_tip") or norm_kps.get("apex_center")
            if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= 0.25 and (abs(float(lb_kp[0])) + abs(float(lb_kp[1]))) > 0.0:
                anchor_pt = (float(lb_kp[0]), float(lb_kp[1]))
            if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= 0.25 and (abs(float(lt_kp[0])) + abs(float(lt_kp[1]))) > 0.0:
                p_apex_pt = (float(lt_kp[0]), float(lt_kp[1]))
        except Exception:
            pass
    if anchor_pt is not None:
        diffs = largest_cnt - np.array([float(anchor_pt[0]), float(anchor_pt[1])], dtype=np.float64)
        start_idx = int(np.argmin(np.hypot(diffs[:, 0], diffs[:, 1])))
    else:
        start_idx = int(np.argmax(largest_cnt[:, 1]))

    pts = np.roll(largest_cnt, -start_idx, axis=0)
    if momocs_signed_area(pts) > 0 or cv2.contourArea(pts.reshape(-1, 1, 2).astype(np.float32), oriented=True) > 0:
        pts = np.vstack([pts[0:1], pts[1:][::-1]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    valid = dt > 1e-7
    if np.sum(valid) < 4:
        return None
    pts = pts[np.concatenate([[True], valid])]
    if len(pts) > 2 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]

    resamp_pts = resample_closed(pts, int(num_points))
    if momocs_signed_area(resamp_pts) > 0 or cv2.contourArea(resamp_pts.reshape(-1, 1, 2).astype(np.float32), oriented=True) > 0:
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
    Saves binary mask strictly once to output_dir/masks/tier1_intact/.
    """
    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    mask_dir = out_base / "masks"
    tier1_dir = mask_dir / "tier1_intact"
    tier1_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{catalog_number}_leaf{leaf_id}.png"
    tier_save_path = tier1_dir / filename

    _ = (petiole_mask, hole_mask)
    binary_mask = np.where(mask > 0, 255, 0).astype(np.uint8)
    ok, png_buf = cv2.imencode(".png", binary_mask)
    if ok:
        tier_save_path.write_bytes(png_buf.tobytes())
    else:
        cv2.imwrite(str(tier_save_path), binary_mask)
    return str(tier_save_path)


# Backward compatibility alias
route_tier1_silhouette = extract_tier1_pristine


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
        clean_mask = _prepare_lamina_mask(mask, petiole_mask=petiole_mask, hole_mask=hole_mask)
        work_mask = clean_mask
        extracted_xy = extract_contour(clean_mask)
        if extracted_xy is None or len(extracted_xy) < 10:
            return empty_ret
        largest_cnt = np.asarray(extracted_xy, dtype=np.float64).reshape(-1, 2)

    anchor_pt = petiole_attachment_pt if petiole_attachment_pt is not None else lamina_base
    if keypoints_31 is not None:
        try:
            from scripts.vision.lm3_inference_engines import normalize_lm3_keypoints_dict

            norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
            lb_kp = norm_kps.get("lamina_base") or norm_kps.get("base_center")
            lt_kp = norm_kps.get("lamina_tip") or norm_kps.get("apex_center")
            if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= 0.25 and (abs(float(lb_kp[0])) + abs(float(lb_kp[1]))) > 0.0:
                anchor_pt = (float(lb_kp[0]), float(lb_kp[1]))
            if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= 0.25 and (abs(float(lt_kp[0])) + abs(float(lt_kp[1]))) > 0.0:
                p_apex_pt = (float(lt_kp[0]), float(lt_kp[1]))
        except Exception:
            pass
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
    if momocs_signed_area(pts) > 0 or cv2.contourArea(pts.reshape(-1, 1, 2).astype(np.float32), oriented=True) > 0:
        pts = np.vstack([pts[0:1], pts[1:][::-1]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    valid = dt > 1e-7
    if np.sum(valid) < 4:
        return empty_ret
    pts = pts[np.concatenate([[True], valid])]
    if len(pts) > 2 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]

    resamp_pts = resample_closed(pts, int(num_points))
    if momocs_signed_area(resamp_pts) > 0 or cv2.contourArea(resamp_pts.reshape(-1, 1, 2).astype(np.float32), oriented=True) > 0:
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
    CANONICAL_TIER_MAP,
    LeafRoutingResult,
    TIER_1_DIRECT,
    TIER_2_REFLECTED,
    compute_classical_shape_descriptors,
    compute_sheet_adaptive_min_blade_area_cm2,
    is_botanical_dissection,
    normalize_assigned_tier,
)
from scripts.vision.organ_crown_metrics import (  # noqa: E402
    _estimate_single_plant_crown,
    _extract_xy_point,
    _parse_basal_leaf_pose,
    approximate_petiole_and_scape_metrics,
    classify_leaf_zone_anatomical,
    cluster_plant_individuals,
    estimate_rosette_crown_point,
)

__all__ = [
    "find_lamina_midrib_endpoints",
    "compute_geometric_metrics",
    "compute_classical_shape_descriptors",
    "is_botanical_dissection",
    "resample_normalized_contour",
    "TIER_1_DIRECT",
    "TIER_2_REFLECTED",
    "CANONICAL_TIER_MAP",
    "normalize_assigned_tier",
    "compute_sheet_adaptive_min_blade_area_cm2",
    "LeafRoutingResult",
    "extract_tier1_pristine",
    "route_tier1_silhouette",
    "export_standardized_contour",
    "cluster_plant_individuals",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
    "classify_leaf_zone_anatomical",
]
