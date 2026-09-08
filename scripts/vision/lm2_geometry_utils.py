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
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

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
