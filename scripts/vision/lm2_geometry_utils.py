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

def compute_geometric_metrics(
    mask: np.ndarray
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

    Returns
    -------
    Tuple[float, float, float, Tuple[int, int], Tuple[int, int]]
        (ucs, solidity, midrib_angle_deg, p_apex, p_base)
    """
    area_mask = float(np.count_nonzero(mask))
    if area_mask < 30:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
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

    # 3. Longitudinal Midrib Axis via cv2.fitLine
    line_params = cv2.fitLine(main_contour, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]
    angle_rad = math.atan2(vy, vx)
    angle_deg = math.degrees(angle_rad)

    # Project contour points along longitudinal midrib axis
    pts = main_contour.reshape(-1, 2).astype(np.float32)
    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy

    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))

    p_apex = (int(pts[min_idx, 0]), int(pts[min_idx, 1]))
    p_base = (int(pts[max_idx, 0]), int(pts[max_idx, 1]))

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
    if h < 5 or w < 5 or np.count_nonzero(mask) == 0:
        return mask.copy()

    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])
    dx = x2 - x1
    dy = y2 - y1
    chord_len = math.hypot(dx, dy)
    if chord_len < 1e-6:
        return mask.copy()

    theta = math.atan2(dy, dx)
    pad = max(h, w) if padding is None else padding
    padded = cv2.copyMakeBorder(mask, pad, pad, pad, pad, cv2.BORDER_CONSTANT, value=0)

    px1 = x1 + float(pad)
    py1 = y1 + float(pad)

    cos2 = math.cos(2.0 * theta)
    sin2 = math.sin(2.0 * theta)

    # 2D Affine Reflection Matrix across the line passing through (px1, py1) at angle theta
    M = np.array([
        [cos2, sin2, px1 * (1.0 - cos2) - py1 * sin2],
        [sin2, -cos2, py1 * (1.0 + cos2) - px1 * sin2]
    ], dtype=np.float32)

    ph, pw = padded.shape[:2]
    reflected = cv2.warpAffine(padded, M, (pw, ph), flags=cv2.INTER_NEAREST)
    combined = cv2.bitwise_or(padded, reflected)

    # Morphological closing to seal the midrib seam
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    sealed = cv2.morphologyEx(combined, cv2.MORPH_CLOSE, kernel)

    y_indices, x_indices = np.where(sealed > 0)
    if len(y_indices) == 0:
        return mask.copy()

    # If the synthesized silhouette fits within original bounds, crop back to original frame
    orig_crop = sealed[pad : pad + h, pad : pad + w]
    crop_fg = np.count_nonzero(orig_crop > 0)
    total_fg = len(y_indices)

    if crop_fg == total_fg:
        return (orig_crop > 127).astype(np.uint8) * 255

    # Otherwise return tight bounding box with 5px padding
    min_y = max(0, int(np.min(y_indices)) - 5)
    max_y = min(ph, int(np.max(y_indices)) + 6)
    min_x = max(0, int(np.min(x_indices)) - 5)
    max_x = min(pw, int(np.max(x_indices)) + 6)
    return (sealed[min_y:max_y, min_x:max_x] > 127).astype(np.uint8) * 255


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
        solidity = float(np.count_nonzero(mask)) / area_hull if area_hull > 0 else 0.0

        approx = cv2.approxPolyDP(pts.reshape(-1, 1, 2), 0.02 * cv2.arcLength(pts.reshape(-1, 1, 2), True), True).reshape(-1, 2)
        has_diagonal_crease = False
        ang_midrib = math.atan2(chord_dy, chord_dx)

        for i in range(len(approx)):
            p_s = approx[i]
            p_e = approx[(i + 1) % len(approx)]
            seg_len = math.hypot(p_e[0] - p_s[0], p_e[1] - p_s[1])
            if seg_len > 0.35 * chord_len:
                ang_seg = math.atan2(p_e[1] - p_s[1], p_e[0] - p_s[0])
                diff_ang = abs((ang_seg - ang_midrib + math.pi) % (2.0 * math.pi) - math.pi)
                diff_deg = math.degrees(diff_ang)
                diff_deg = min(diff_deg, 180.0 - diff_deg)
                if diff_deg > 15.0:
                    has_diagonal_crease = True
                    break

        if has_diagonal_crease or solidity >= 0.72:
            default_res["is_irregular_fold"] = True

    return default_res


def is_botanical_dissection(
    contour: np.ndarray,
    solidity: Optional[float] = None,
    p_apex: Optional[Tuple[int, int]] = None,
    p_base: Optional[Tuple[int, int]] = None,
    blade_width: Optional[float] = None,
    min_solidity: float = 0.48,
    max_solidity: float = 0.72,
    defect_depth_ratio: float = 0.08,
    min_defect_count: int = 3,
    check_solidity_range: bool = True
) -> bool:
    """
    Evaluates whether a leaf contour with lower solidity (0.48 <= Solidity < 0.72)
    represents authentic botanical dissection (lyrate/pinnatifid lobes with bilateral sinuses)
    rather than asymmetric foreign leaf occlusion or damage.

    Detection Rules:
      - Filter out minor crenation noise: retain only significant defects where
        defect depth exceeds 8% of maximum blade width.
      - Evaluate Defect Topology:
        * Foreign leaf occlusion: 1 or 2 large, asymmetric, localized deficits on one side.
        * Botanical dissection: Multiple (>= 3), alternating, bilaterally distributed
          sinuses along the margin (present on both sides of the longitudinal midline).
      - If significant defects are >= 3 and distributed on both sides of midline:
        Classify as is_botanical_dissection = True.

    Parameters
    ----------
    contour : np.ndarray
        Contour points of the leaf instance.
    solidity : Optional[float]
        Calculated solidity of the leaf instance. If None, computed from contour.
    p_apex : Optional[Tuple[int, int]]
        Apex coordinate (x, y). If None, estimated via cv2.fitLine along contour.
    p_base : Optional[Tuple[int, int]]
        Base/petiole coordinate (x, y). If None, estimated via cv2.fitLine along contour.
    blade_width : Optional[float]
        Maximum blade width. If None, computed via cv2.minAreaRect.
    min_solidity : float
        Lower bound of solidity range for dissection check (default 0.48).
    max_solidity : float
        Upper bound of solidity range for dissection check (default 0.72).
    defect_depth_ratio : float
        Defect depth threshold relative to blade width (default 0.08).
    min_defect_count : int
        Minimum significant defect count (default 3).
    check_solidity_range : bool
        Whether to enforce min_solidity <= solidity < max_solidity (default True).

    Returns
    -------
    bool
        True if the contour exhibits genuine botanical dissection; False otherwise.
    """
    cnt = contour.reshape(-1, 2).astype(np.int32)
    if len(cnt) < 4:
        return False

    cnt_f = cnt.astype(np.float32)
    if solidity is None:
        area = float(cv2.contourArea(cnt_f))
        hull = cv2.convexHull(cnt_f)
        hull_area = float(cv2.contourArea(hull))
        solidity = float(area / hull_area) if hull_area > 0 else 0.0

    if check_solidity_range:
        if not (min_solidity <= solidity < max_solidity):
            return False

    if blade_width is None or blade_width <= 0:
        rect = cv2.minAreaRect(cnt_f)
        blade_width = float(min(rect[1]))
    if blade_width <= 0:
        return False

    min_defect_depth = defect_depth_ratio * blade_width

    cnt_cv = cnt.reshape(-1, 1, 2)
    hull_indices = cv2.convexHull(cnt_cv, returnPoints=False)
    if hull_indices is None or len(hull_indices) < 3:
        return False

    defects = cv2.convexityDefects(cnt_cv, hull_indices)
    if defects is None or len(defects) == 0:
        return False

    # Extract defect depths (OpenCV stores fixpt_depth as distance * 256.0)
    def_arr = defects.reshape(-1, 4)
    depths = def_arr[:, 3].astype(np.float64) / 256.0
    farthest_indices = def_arr[:, 2].astype(np.int32)

    # Filter significant defects exceeding depth threshold
    sig_mask = depths >= min_defect_depth
    if np.sum(sig_mask) < min_defect_count:
        return False

    sig_farthest_pts = cnt[farthest_indices[sig_mask]]

    # Longitudinal midline vector from apex to base
    if p_apex is None or p_base is None:
        line_params = cv2.fitLine(cnt_f, cv2.DIST_L2, 0, 0.01, 0.01)
        vx, vy, x0, y0 = [float(v[0]) for v in line_params]
        projections = (cnt_f[:, 0] - x0) * vx + (cnt_f[:, 1] - y0) * vy
        min_idx = int(np.argmin(projections))
        max_idx = int(np.argmax(projections))
        p_apex = (int(cnt_f[min_idx, 0]), int(cnt_f[min_idx, 1]))
        p_base = (int(cnt_f[max_idx, 0]), int(cnt_f[max_idx, 1]))

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

    # Botanical dissection: >= min_defect_count defects distributed on both sides (side_a >= 1 and side_b >= 1)
    is_dissected = (total_sig >= min_defect_count) and (side_a_count >= 1) and (side_b_count >= 1)
    return bool(is_dissected)


def resample_normalized_contour(
    mask: np.ndarray,
    num_points: int = 120
) -> Optional[np.ndarray]:
    """
    Extracts resampled, normalized 2D coordinates (x_norm, y_norm) of shape (num_points, 2)
    from a binary mask, with consistent clockwise orientation and min-max normalization.
    """
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return None

    largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2)
    if len(largest_cnt) < 10:
        return None

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

    # Enforce consistent clockwise contour orientation (positive signed area)
    signed_area = 0.5 * float(np.sum(x_resamp * np.roll(y_resamp, -1) - np.roll(x_resamp, -1) * y_resamp))
    if signed_area < 0:
        x_resamp = x_resamp[::-1]
        y_resamp = y_resamp[::-1]

    min_x, max_x = np.min(x_resamp), np.max(x_resamp)
    min_y, max_y = np.min(y_resamp), np.max(y_resamp)
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)

    x_norm = np.round((x_resamp - min_x) / span_x, 6)
    y_norm = np.round((y_resamp - min_y) / span_y, 6)
    return np.column_stack([x_norm, y_norm])


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


def classify_and_route_leaf(
    mask: np.ndarray,
    config: Optional[Any] = None,
    min_ucs: float = 0.85,
    min_solidity: float = 0.72,
    dissection_min_solidity: float = 0.50,
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
) -> LeafRoutingResult:
    """
    Comprehensive geometric quality gatekeeper and routing engine:
      1. Evaluates midrib pose and initial solidity/UCS metrics.
      2. Straight-chord fold detection:
         - If folded along midrib (is_folded = True), symmetrically reflects hemi-blade
           into synthesized bilateral silhouette and routes to Tier 2.
         - If irregular/diagonal fold (is_irregular_fold = True), routes to Failed QC.
      3. Botanical dissection inspection:
         - If min_dissected <= solidity < min_solidity and bilateral convexity defects indicate
           botanical lobing/dissection (or taxon has lyrate tendency), accepts leaf into Tier 1 (Dissected).
      4. Standard pristine gatekeeping:
         - Leaves with UCS >= min_ucs and Solidity >= min_solidity route to Tier 1.
      5. Fallback hemi-blade reflection (Tier 2):
         - Partially occluded leaves attempt 2-path reflection or route to Failed QC.

    Parameters
    ----------
    mask : np.ndarray
        Binary leaf mask.
    config : Optional[Any]
        Optional PipelineConfig or configuration object with threshold parameters.
    min_ucs : float
        UCS threshold for Tier 1 (default 0.85).
    min_solidity : float
        Solidity threshold for Tier 1 (default 0.72).
    dissection_min_solidity : float
        Lower solidity bound for botanical dissection inspection (default 0.50).
    extract_contours : bool
        Whether to extract resampled normalized contour coordinates (default True).
    num_contour_points : int
        Number of standardized contour points (default 120).
    taxon : Optional[str]
        Botanical taxon / scientific name for lyrate dissection tendency check.
    taxa_with_lyrate_tendency : Optional[List[str]]
        List of taxa with known lyrate/pinnatifid basal leaf dissection tendency.
    max_chord_deviation_ratio : Optional[float]
        Fold detection maximum chord deviation ratio (default 0.035).
    max_folded_aspect_ratio : Optional[float]
        Fold detection maximum width-to-length ratio (default 0.42).
    sinus_defect_min_depth_ratio : Optional[float]
        Sinus convexity defect minimum depth ratio (default 0.08).
    min_bilateral_sinus_count : Optional[int]
        Minimum bilateral sinus defect count for botanical dissection (default 3).
    min_chord_length_px : Optional[float]
        Minimum chord length in pixels for fold detection (default 150.0).
    fold_detection_enabled : Optional[bool]
        Whether fold detection is enabled (default True).
    dissection_enabled : Optional[bool]
        Whether dissection inspection is enabled (default True).

    Returns
    -------
    LeafRoutingResult
        Structured metadata and processed mask for downstream extraction.
    """
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

        # Dissection configuration
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
            if taxa_with_lyrate_tendency is None:
                taxa_with_lyrate_tendency = getattr(diss_cfg, "taxa_with_lyrate_tendency", None)
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
                if taxa_with_lyrate_tendency is None:
                    taxa_with_lyrate_tendency = diss_dict.get("taxa_with_lyrate_tendency")

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

    if taxa_with_lyrate_tendency is None:
        taxa_with_lyrate_tendency = [
            "Packera paupercula",
            "Packera plattensis",
            "Packera paupercula var. paupercula",
            "Packera paupercula var. savannarum",
        ]
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

    if mask is None or np.count_nonzero(mask) < 20:
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

    ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
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

    # Taxon-aware lyrate tendency check
    taxon_has_lyrate_tendency = False
    if taxon and taxa_with_lyrate_tendency:
        t_clean = taxon.strip().lower()
        for lyr in taxa_with_lyrate_tendency:
            if lyr.strip().lower() in t_clean or t_clean in lyr.strip().lower():
                taxon_has_lyrate_tendency = True
                break

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
        "taxon_has_lyrate_tendency": taxon_has_lyrate_tendency,
    }

    if is_folded:
        # Route folded leaf directly to Tier 2 with synthesized bilateral silhouette
        synth_mask = fold_info["synthesized_mask"]
        if synth_mask is not None and np.count_nonzero(synth_mask) > 0:
            coords = resample_normalized_contour(synth_mask, num_contour_points) if extract_contours else None
            synth_ucs, synth_sol, _, _, _ = compute_geometric_metrics(synth_mask)
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
        coords = resample_normalized_contour(mask, num_contour_points) if extract_contours else None
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

    # 3. Botanical Dissection Inspection (min_dissected <= solidity < min_solidity)
    if dissection_enabled and (dissection_min_solidity <= solidity < min_solidity):
        dissected = is_botanical_dissection(
            main_cnt,
            solidity=solidity,
            p_apex=p_apex,
            p_base=p_base,
            min_solidity=dissection_min_solidity,
            max_solidity=min_solidity,
            defect_depth_ratio=sinus_defect_min_depth_ratio,
            min_defect_count=min_bilateral_sinus_count,
        )
        if dissected or taxon_has_lyrate_tendency:
            coords = resample_normalized_contour(mask, num_contour_points) if extract_contours else None
            meta["reflection_applied"] = False
            return LeafRoutingResult(
                tier="Tier 1 (Dissected)",
                assigned_tier="tier1",
                solidity=solidity,
                ucs=ucs,
                is_folded=False,
                is_irregular_fold=False,
                is_dissected=True,
                processed_mask=mask,
                contour_coords=coords,
                metadata=meta
            )

    # 4. Fallback: Tier 2 Hemi-Blade Reflection for Occluded/Damaged Margins
    path, reflected_mask = extract_tier2_reflected(
        mask=mask,
        catalog_number="TEMP",
        leaf_id=1,
        p_apex=p_apex,
        p_base=p_base,
        output_dir=None
    )

    if reflected_mask is not None:
        ref_ucs, ref_sol, _, _, _ = compute_geometric_metrics(reflected_mask)
        coords = resample_normalized_contour(reflected_mask, num_contour_points) if extract_contours else None
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
    cv2.imwrite(str(main_save_path), bin_mask)
    cv2.imwrite(str(tier_save_path), bin_mask)
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
    Aligns mask along midrib axis, assesses margin defect profile of upper vs
    lower half, reflects the clean half across the midrib axis symmetrically
    using NumPy vectorization (no nested pixel loops), and applies morphological
    closing to produce a watertight closed silhouette.

    Parameters
    ----------
    mask : np.ndarray
        Original binary leaf mask.
    catalog_number : str
        Specimen voucher catalog identifier.
    leaf_id : int
        Instance index on specimen sheet.
    p_apex : Tuple[int, int]
        Apex coordinate (x, y).
    p_base : Tuple[int, int]
        Base/petiole coordinate (x, y).
    output_dir : Optional[Path]
        Directory to save mask outputs. If None, saves are skipped.
    threshold_defect : float
        Maximum defect score to consider a half-blade intact.
    min_solidity : float
        Minimum solidity required for the synthesized bilateral silhouette.
    min_ucs : float
        Minimum UCS required for the synthesized bilateral silhouette.

    Returns
    -------
    Tuple[Optional[str], Optional[np.ndarray]]
        (saved_path, reflected_mask_binary) or (None, None) if reflection fails QC.
    """
    h, w = mask.shape[:2]
    if h < 15 or w < 15:
        return None, None

    dx = p_base[0] - p_apex[0]
    dy = p_base[1] - p_apex[1]
    angle = math.degrees(math.atan2(dy, dx))

    cx, cy = w / 2.0, h / 2.0
    rot_mat = cv2.getRotationMatrix2D((cx, cy), -angle, 1.0)
    aligned = cv2.warpAffine(mask, rot_mat, (w, h), flags=cv2.INTER_NEAREST)

    y_indices, _ = np.where(aligned > 0)
    if len(y_indices) < 30:
        return None, None

    y_midrib = int(np.median(y_indices))
    if y_midrib <= 1 or y_midrib >= h - 2:
        return None, None

    upper_half = np.zeros_like(aligned)
    lower_half = np.zeros_like(aligned)
    upper_half[:y_midrib, :] = aligned[:y_midrib, :]
    lower_half[y_midrib:, :] = aligned[y_midrib:, :]

    # Fully vectorized margin profile detection across columns
    upper_has_fg = upper_half > 0
    lower_has_fg = lower_half > 0

    col_has_upper = np.any(upper_has_fg, axis=0)
    col_has_lower = np.any(lower_has_fg, axis=0)

    # Fast column min/max using vectorized argmax
    upper_profile = np.where(col_has_upper, np.argmax(upper_has_fg, axis=0), y_midrib)
    lower_profile = np.where(col_has_lower, h - 1 - np.argmax(lower_has_fg[::-1, :], axis=0), y_midrib)

    # Assess margin defect score (jumps/discontinuities)
    upper_diff = np.abs(np.diff(upper_profile))
    lower_diff = np.abs(np.diff(lower_profile))
    upper_defects = float(np.sum(upper_diff > 8) * 10.0 + np.std(upper_diff))
    lower_defects = float(np.sum(lower_diff > 8) * 10.0 + np.std(lower_diff))

    if upper_defects < threshold_defect and upper_defects <= lower_defects:
        selected_half = "upper"
    elif lower_defects < threshold_defect:
        selected_half = "lower"
    else:
        return None, None

    # Fully vectorized slice reflection in NumPy (zero manual pixel loops)
    reflected = np.zeros_like(aligned)
    if selected_half == "upper":
        reflected[:y_midrib, :] = upper_half[:y_midrib, :]
        max_dy = min(y_midrib, h - y_midrib)
        if max_dy > 0:
            reflected[y_midrib:y_midrib + max_dy, :] = np.flip(
                upper_half[y_midrib - max_dy:y_midrib, :], axis=0
            )
    else:
        reflected[y_midrib:, :] = lower_half[y_midrib:, :]
        max_dy = min(y_midrib, h - y_midrib)
        if max_dy > 0:
            reflected[y_midrib - max_dy:y_midrib, :] = np.flip(
                lower_half[y_midrib:y_midrib + max_dy, :], axis=0
            )

    # Morphological closing to seal midrib seam
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    reflected = cv2.morphologyEx(reflected, cv2.MORPH_CLOSE, kernel)
    reflected_bin = (reflected > 127).astype(np.uint8) * 255

    # Verify synthesized silhouette quality
    ref_ucs, ref_solidity, _, _, _ = compute_geometric_metrics(reflected_bin)
    if ref_solidity < min_solidity or ref_ucs < min_ucs:
        return None, None

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
        cv2.imwrite(str(tier_save_path), reflected_bin)

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
    num_points: int = 120
) -> Optional[str]:
    """
    Extracts vectorized 2D (x, y) boundary coordinates from a binary silhouette mask,
    enforces consistent clockwise orientation, normalizes coordinates (x_norm, y_norm)
    for downstream Elliptic Fourier Analysis, and saves to
    output_dir/contours/{catalogNumber}_leaf{id}.csv.
    """
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return None

    largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2)
    if len(largest_cnt) < 10:
        return None

    # Resample smoothly along perimeter to standardized point count
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

    # Enforce consistent clockwise contour orientation (in image space, clockwise has positive signed area)
    signed_area = 0.5 * float(np.sum(x_resamp * np.roll(y_resamp, -1) - np.roll(x_resamp, -1) * y_resamp))
    if signed_area < 0:
        x_resamp = x_resamp[::-1]
        y_resamp = y_resamp[::-1]

    min_x, max_x = np.min(x_resamp), np.max(x_resamp)
    min_y, max_y = np.min(y_resamp), np.max(y_resamp)
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)

    contour_df = pd.DataFrame({
        "catalogNumber": catalog_number,
        "leaf_id": leaf_id,
        "point_index": np.arange(num_points),
        "x": np.round(x_resamp, 2),
        "y": np.round(y_resamp, 2),
        "x_norm": np.round((x_resamp - min_x) / span_x, 6),
        "y_norm": np.round((y_resamp - min_y) / span_y, 6),
    })

    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    contour_dir = out_base / "contours"
    contour_dir.mkdir(parents=True, exist_ok=True)
    out_csv = contour_dir / f"{catalog_number}_leaf{leaf_id}.csv"
    contour_df.to_csv(out_csv, index=False)
    return str(out_csv)


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

    cat = getattr(plant_candidates[0], "catalog_number", "UNKNOWN") if plant_candidates else "UNKNOWN"
    save_path = output_dir / "cropped_patches" / "rosettes_dense" / f"{cat}_p{plant_id}_rosette.jpg"
    save_path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(save_path), rosette_crop)
    return str(save_path)
