#!/usr/bin/env python3
"""
===============================================================================
Module: organ_crown_metrics.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    SAM 2 prompt generation, topological connectivity gatekeeping, multi-plant
    DBSCAN spatial clustering, Hough ruler scale calibration, Tier 4 dense
    rosette patch masking, and rosette crown / petiole / bent-scape metrics.
===============================================================================
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from scripts.vision.lm2_geometry_utils import (
    compute_geometric_metrics,
    trim_petiole_tail,
)

__all__ = [
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
        # Downscale before grayscale conversion if extremely large for rapid Hough transform
        rh, rw = region.shape[:2]
        max_dim = max(rh, rw)
        scale_factor = 1.0
        if max_dim > 1500:
            scale_factor = 1500.0 / max_dim
            region_small = cv2.resize(
                region,
                (int(rw * scale_factor), int(rh * scale_factor)),
                interpolation=cv2.INTER_AREA,
            )
            gray = cv2.cvtColor(region_small, cv2.COLOR_BGR2GRAY) if len(region_small.shape) == 3 else region_small
        else:
            gray = cv2.cvtColor(region, cv2.COLOR_BGR2GRAY) if len(region.shape) == 3 else region

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
            resolved_raw_path = raw_image_path
            if not resolved_raw_path.exists():
                for _ext in (".webp", ".jpg", ".jpeg", ".png"):
                    _cand = resolved_raw_path.with_suffix(_ext)
                    if _cand.exists():
                        resolved_raw_path = _cand
                        break
            if not resolved_raw_path.exists():
                return None
            raw_img = cv2.imread(str(resolved_raw_path))
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
    try:
        cv2.imwrite(str(save_path), rosette_crop, [int(cv2.IMWRITE_JPEG_QUALITY), 85])
    except TypeError:
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
