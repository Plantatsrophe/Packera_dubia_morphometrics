#!/usr/bin/env python3
"""
===============================================================================
Module: organ_crown_metrics.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    LeafMachine3 (YOLO26x / UNet++ v9) topological connectivity gatekeeping,
    multi-plant DBSCAN spatial clustering,
    and 31-keypoint rosette crown / petiole / bent-scape metrics.
===============================================================================
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from scripts.vision.leaf_geometry_utils import (
    compute_geometric_metrics,
    find_lamina_midrib_endpoints,
)

__all__ = [
    "is_topologically_connected_mask",
    "cluster_plant_individuals",
    "cluster_voucher_plants_dbscan",
    "route_tier3_open_curve",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
    "classify_leaf_zone_anatomical",
]


# =============================================================================
# 3B. Topological Gatekeeping
# =============================================================================

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
# 6. Additional Caliper Utilities (Tier 3)
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
    """
    Extract standardized `(p_base, p_apex, plant_individual_id, geodesic_petiole_skeleton_px,
    blade_length_px)` and LeafMachine3 31-keypoint pose attributes (`p_petiole_tip`,
    `pose_petiole_polyline_px`, `ray_dir`, `petiole_width_px`, `petiole_width_mm`) from a leaf record.
    """
    if item is None:
        return None

    pid = 0
    ox, oy = 0.0, 0.0
    p_base: Optional[Tuple[float, float]] = None
    p_apex: Optional[Tuple[float, float]] = None
    p_petiole_tip: Optional[Tuple[float, float]] = None
    ray_dir: Optional[Tuple[float, float]] = None
    has_lm3_petiole_pose = False
    pose_petiole_polyline_px = 0.0
    geodesic_px = 0.0
    blade_px: Optional[float] = None
    petiole_width_px: Optional[float] = None
    petiole_width_mm: Optional[float] = None
    leaf_id: Optional[int] = None
    keypoints_31: Optional[Dict[str, Any]] = None
    kp_in_sheet_space = False

    def _safe_pos_float(val: Any) -> Optional[float]:
        if val is None:
            return None
        try:
            fval = float(val)
        except (TypeError, ValueError):
            return None
        return fval if (np.isfinite(fval) and fval > 0.0) else None

    if isinstance(item, dict):
        pid = int(item.get("plant_individual_id", 0) or 0)
        leaf_id = item.get("leaf_id")
        origin = item.get("mask_origin") or item.get("origin")
        if origin is not None and len(origin) == 2 and not item.get("coords_in_sheet_space", False):
            ox, oy = float(origin[0]), float(origin[1])

        keypoints_31 = item.get("keypoints_31")
        kp_in_sheet_space = bool(item.get("keypoints_in_sheet_space", item.get("coords_in_sheet_space", False)))

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

        if item.get("p_petiole_tip") is not None:
            ptip = item["p_petiole_tip"]
            p_petiole_tip = (float(ptip[0]) + (0.0 if kp_in_sheet_space else ox), float(ptip[1]) + (0.0 if kp_in_sheet_space else oy))
        if item.get("ray_dir") is not None and len(item["ray_dir"]) == 2:
            rd = item["ray_dir"]
            ray_dir = (float(rd[0]), float(rd[1]))
        if item.get("has_lm3_petiole_pose") is not None:
            has_lm3_petiole_pose = bool(item["has_lm3_petiole_pose"])
        if item.get("pose_petiole_polyline_px") is not None:
            pose_petiole_polyline_px = max(0.0, float(item["pose_petiole_polyline_px"] or 0.0))

        petiole_width_mm = _safe_pos_float(item.get("petiole_width_mm"))
        petiole_width_px = _safe_pos_float(item.get("petiole_width_px", item.get("petiole_width")))

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
        pa_m, pb_m = find_lamina_midrib_endpoints(item)
        p_base = (float(pb_m[0]), float(pb_m[1]))
        p_apex = (float(pa_m[0]), float(pa_m[1]))
        geodesic_px = 0.0

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

    elif hasattr(item, "p_base") or hasattr(item, "mask") or hasattr(item, "keypoints_31"):
        pid = int(getattr(item, "plant_individual_id", 0) or 0)
        leaf_id = getattr(item, "leaf_id", None)
        origin = getattr(item, "mask_origin", (0, 0)) or (0, 0)
        ox, oy = float(origin[0]), float(origin[1])
        keypoints_31 = getattr(item, "keypoints_31", None)
        kp_in_sheet_space = bool(getattr(item, "keypoints_in_sheet_space", False))
        petiole_width_mm = _safe_pos_float(getattr(item, "petiole_width_mm", None))
        petiole_width_px = _safe_pos_float(getattr(item, "petiole_width_px", getattr(item, "petiole_width", None)))
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

    # Parse LeafMachine3 31-keypoint pose skeleton when present (supports dict or (31, 3) ndarray)
    if isinstance(keypoints_31, np.ndarray) and keypoints_31.ndim == 2 and keypoints_31.shape[0] >= 8 and keypoints_31.shape[1] >= 2:
        _kp_arr = keypoints_31
        def _row_tup(r: np.ndarray) -> Tuple[float, float, float]:
            return (float(r[0]), float(r[1]), float(r[2]) if len(r) >= 3 else 1.0)
        _kp_dict: Dict[str, Tuple[float, float, float]] = {
            "lamina_tip": _row_tup(_kp_arr[0]),
            "lamina_base": _row_tup(_kp_arr[1]),
        }
        for _pi in range(min(6, _kp_arr.shape[0] - 2)):
            _kp_dict[f"petiole_{_pi + 1}"] = _row_tup(_kp_arr[2 + _pi])
        for _mi in range(min(15, max(0, _kp_arr.shape[0] - 8))):
            _kp_dict[f"midvein_{_mi + 1}"] = _row_tup(_kp_arr[8 + _mi])
        keypoints_31 = _kp_dict

    if isinstance(keypoints_31, dict) and len(keypoints_31) > 0:
        kp_ox = 0.0 if kp_in_sheet_space else ox
        kp_oy = 0.0 if kp_in_sheet_space else oy

        def _extract_kp(name: str, min_conf: float = 0.50) -> Optional[Tuple[float, float, float]]:
            val = keypoints_31.get(name)
            if val is None:
                return None
            if isinstance(val, dict):
                x = float(val.get("x", 0.0))
                y = float(val.get("y", 0.0))
                c = float(val.get("conf", val.get("confidence", 1.0)))
            elif isinstance(val, (tuple, list, np.ndarray)) and len(val) >= 2:
                x = float(val[0])
                y = float(val[1])
                c = float(val[2]) if len(val) >= 3 else 1.0
            else:
                return None
            if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(c)):
                return None
            if c < min_conf or (x == 0.0 and y == 0.0):
                return None
            return (x + kp_ox, y + kp_oy, c)

        lb_kp = _extract_kp("lamina_base", min_conf=0.50)
        lt_kp = _extract_kp("lamina_tip", min_conf=0.50)
        mv_kps = [kp for kp in (_extract_kp(f"midvein_{i}", min_conf=0.50) for i in range(1, 16)) if kp is not None]
        pet_kps = [(i, kp) for i in range(1, 7) for kp in (_extract_kp(f"petiole_{i}", min_conf=0.50),) if kp is not None]

        if lb_kp is not None:
            p_base = (lb_kp[0], lb_kp[1])
        elif pet_kps and p_base is None:
            p_base = (pet_kps[0][1][0], pet_kps[0][1][1])

        if lt_kp is not None:
            p_apex = (lt_kp[0], lt_kp[1])
        elif mv_kps and p_apex is None and p_base is not None:
            far_mv = max(mv_kps, key=lambda kp: math.hypot(kp[0] - p_base[0], kp[1] - p_base[1]))
            p_apex = (far_mv[0], far_mv[1])

        if pet_kps:
            has_lm3_petiole_pose = True
            # Distal crown-ward endpoint of petiole skeleton (petiole_6 or highest valid index)
            p_petiole_tip = (pet_kps[-1][1][0], pet_kps[-1][1][1])

            # Polyline arc length along lamina_base -> petiole_1 -> ... -> petiole_6 when >= 2 valid petiole_1..6 keypoints exist
            if len(pet_kps) >= 2:
                poly_pts: List[Tuple[float, float]] = []
                if lb_kp is not None:
                    poly_pts.append((lb_kp[0], lb_kp[1]))
                for _, pkp in pet_kps:
                    pt_xy = (pkp[0], pkp[1])
                    if not poly_pts or math.hypot(pt_xy[0] - poly_pts[-1][0], pt_xy[1] - poly_pts[-1][1]) > 1e-4:
                        poly_pts.append(pt_xy)
                if len(poly_pts) >= 2:
                    pose_petiole_polyline_px = float(
                        sum(
                            math.hypot(poly_pts[k + 1][0] - poly_pts[k][0], poly_pts[k + 1][1] - poly_pts[k][1])
                            for k in range(len(poly_pts) - 1)
                        )
                    )

            # Align inward ray direction from lamina_tip/midvein through lamina_base and petiole_1..6
            spine_pts: List[Tuple[float, float]] = []
            if lt_kp is not None:
                spine_pts.append((lt_kp[0], lt_kp[1]))
            for mkp in reversed(mv_kps):
                spine_pts.append((mkp[0], mkp[1]))
            if lb_kp is not None:
                spine_pts.append((lb_kp[0], lb_kp[1]))
            for _, pkp in pet_kps:
                spine_pts.append((pkp[0], pkp[1]))

            if len(spine_pts) >= 2 and p_base is not None:
                pet_vec_x = p_petiole_tip[0] - p_base[0]
                pet_vec_y = p_petiole_tip[1] - p_base[1]
                pet_vec_len = math.hypot(pet_vec_x, pet_vec_y)
                spine_arr = np.asarray(spine_pts, dtype=np.float32)
                lp = cv2.fitLine(spine_arr, cv2.DIST_L2, 0, 0.01, 0.01)
                vx_s, vy_s = float(lp[0][0]), float(lp[1][0])
                ref_x = spine_pts[-1][0] - spine_pts[0][0]
                ref_y = spine_pts[-1][1] - spine_pts[0][1]
                if (ref_x * vx_s + ref_y * vy_s) < 0.0:
                    vx_s, vy_s = -vx_s, -vy_s
                if pet_vec_len > 2.0:
                    ux_pet, uy_pet = pet_vec_x / pet_vec_len, pet_vec_y / pet_vec_len
                    bx_dir = 0.65 * ux_pet + 0.35 * vx_s
                    by_dir = 0.65 * uy_pet + 0.35 * vy_s
                    bnorm = math.hypot(bx_dir, by_dir)
                    ray_dir = (bx_dir / bnorm, by_dir / bnorm) if bnorm > 1e-6 else (ux_pet, uy_pet)
                else:
                    ray_dir = (vx_s, vy_s)

        # Extract Module 13 petiole_width or width_1/width_2 fallback when present
        if petiole_width_mm is None:
            petiole_width_mm = _safe_pos_float(keypoints_31.get("petiole_width_mm"))
        if petiole_width_px is None:
            petiole_width_px = _safe_pos_float(
                keypoints_31.get("petiole_width_px", keypoints_31.get("petiole_width"))
            )
        if petiole_width_px is None and petiole_width_mm is None:
            w1_kp = _extract_kp("width_1", min_conf=0.50)
            w2_kp = _extract_kp("width_2", min_conf=0.50)
            if w1_kp is not None and w2_kp is not None:
                w_dist = math.hypot(w1_kp[0] - w2_kp[0], w1_kp[1] - w2_kp[1])
                if w_dist > 0.0:
                    petiole_width_px = float(w_dist)

    if p_base is None or p_apex is None:
        return None

    chord_len = math.hypot(p_apex[0] - p_base[0], p_apex[1] - p_base[1])
    if blade_px is None or blade_px <= 0.0:
        blade_px = chord_len

    if ray_dir is None and chord_len > 1e-6:
        ray_dir = ((p_base[0] - p_apex[0]) / chord_len, (p_base[1] - p_apex[1]) / chord_len)

    return {
        "leaf_id": leaf_id,
        "plant_individual_id": pid,
        "p_base": p_base,
        "p_apex": p_apex,
        "p_petiole_tip": p_petiole_tip,
        "ray_dir": ray_dir,
        "has_lm3_petiole_pose": bool(has_lm3_petiole_pose),
        "pose_petiole_polyline_px": max(0.0, float(pose_petiole_polyline_px)),
        "geodesic_petiole_skeleton_px": max(0.0, max(float(geodesic_px), float(pose_petiole_polyline_px))),
        "blade_length_px": max(0.0, float(blade_px)),
        "petiole_width_px": petiole_width_px,
        "petiole_width_mm": petiole_width_mm,
    }


def _estimate_single_plant_crown(
    poses: List[Dict[str, Any]],
    root_pts: List[Tuple[float, float]],
) -> Tuple[float, float]:
    """
    Compute the caudex crown point `p_crown = (x_crown, y_crown)` for a single plant individual
    via weighted least-squares intersection of inward basal-leaf midrib/petiole rays anchored by
    LeafMachine3 `petiole_1..6` distal endpoints (`petiole_6`), root detections, and/or the
    innermost `p_base` coordinates.
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

    # Incorporate aligned LM3 ray directions and distal petiole_6 endpoints when available
    petiole_tips: List[Tuple[float, float]] = []
    ray_weights = np.ones(len(poses), dtype=np.float64)
    ray_origins = bases.copy()
    for idx_p, p in enumerate(poses):
        if p.get("has_lm3_petiole_pose"):
            ray_weights[idx_p] = 2.5
        ptip = p.get("p_petiole_tip")
        if ptip is not None and len(ptip) == 2:
            petiole_tips.append((float(ptip[0]), float(ptip[1])))
            ray_origins[idx_p] = np.array([float(ptip[0]), float(ptip[1])], dtype=np.float64)
        rdir = p.get("ray_dir")
        if rdir is not None and len(rdir) == 2 and p.get("has_lm3_petiole_pose"):
            rnorm = math.hypot(float(rdir[0]), float(rdir[1]))
            if rnorm > 1e-6:
                eff_len = max(float(lengths[idx_p]), 25.0)
                vecs[idx_p] = np.array([float(rdir[0]) / rnorm, float(rdir[1]) / rnorm], dtype=np.float64) * eff_len
                lengths[idx_p] = eff_len

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

    has_petiole_tip_anchor = len(petiole_tips) > 0
    pet_tip_anchor = (
        np.array(
            [np.median([p[0] for p in petiole_tips]), np.median([p[1] for p in petiole_tips])],
            dtype=np.float64,
        )
        if has_petiole_tip_anchor
        else None
    )

    # Determine innermost p_base / petiole_6 / root coordinate anchor
    if near_roots:
        root_anchor = np.array(
            [np.median([p[0] for p in near_roots]), np.median([p[1] for p in near_roots])],
            dtype=np.float64,
        )
        dists_to_root = np.hypot(bases[:, 0] - root_anchor[0], bases[:, 1] - root_anchor[1])
        inner_idx = int(np.argmin(dists_to_root))
        innermost_base = bases[inner_idx]
        if pet_tip_anchor is not None:
            p_anchor = 0.50 * pet_tip_anchor + 0.35 * root_anchor + 0.15 * innermost_base
        else:
            p_anchor = 0.70 * root_anchor + 0.30 * innermost_base
        has_root_anchor = True
    elif pet_tip_anchor is not None:
        has_root_anchor = True
        innermost_base = pet_tip_anchor.copy()
        p_anchor = pet_tip_anchor.copy()
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
            if has_petiole_tip_anchor:
                return float(p_anchor[0]), float(p_anchor[1])
            if int(np.sum(valid_mask)) == 1:
                idx0 = int(np.where(valid_mask)[0][0])
                u0 = vecs[idx0] / lengths[idx0]
                b0 = ray_origins[idx0]
                t_root = float(np.dot(p_anchor - b0, u0))
                if t_root >= 0.0:
                    proj_pt = b0 + min(t_root, 2.5 * mean_len) * u0
                    crown_pt = 0.65 * proj_pt + 0.35 * p_anchor
                    return float(crown_pt[0]), float(crown_pt[1])
            return float(p_anchor[0]), float(p_anchor[1])
        return float(innermost_base[0]), float(innermost_base[1])

    # Weighted least-squares intersection of 2D inward midrib/petiole lines:
    # Each line i has normal projector P_i = w_i * (I - u_i u_i^T) passing through ray_origins[i].
    A_mat = np.zeros((2, 2), dtype=np.float64)
    b_vec = np.zeros(2, dtype=np.float64)
    valid_indices = np.where(valid_mask)[0]
    unit_list = []
    for idx_i in valid_indices:
        u_i = vecs[idx_i] / lengths[idx_i]
        w_i = float(ray_weights[idx_i])
        unit_list.append(u_i)
        P_i = w_i * (np.eye(2, dtype=np.float64) - np.outer(u_i, u_i))
        A_mat += P_i
        b_vec += P_i @ ray_origins[idx_i]

    eigvals = np.linalg.eigvalsh(A_mat)
    min_eig = float(np.min(eigvals))
    max_reach = max(2.5 * float(np.max(lengths[valid_mask])), 60.0)

    if min_eig > 1e-3:
        # Rays have distinct angles -> solve least-squares intersection with anchor regularization
        reg_w = (0.65 * float(len(petiole_tips))) if has_petiole_tip_anchor else (0.35 if has_root_anchor else 1e-4)
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

    # Fallback when rays are nearly parallel or divergent: anchor by roots / petiole_6 / innermost p_base
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
    Compute the caudex crown point `p_crown` per `plant_individual_id` via weighted least-squares
    intersection of inward basal-leaf midrib/petiole rays (`p_apex -> p_base -> petiole_6`)
    anchored by LeafMachine3 `petiole_1..6` pose skeletons, `roots` detections, and/or the
    innermost `p_base` coordinates.

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

    if isinstance(basal_leaves, dict) and not (
        "p_base" in basal_leaves or "mask" in basal_leaves or "keypoints_31" in basal_leaves
    ):
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
    Compute empirical approximations for basal petiole length, petiole width, and flowering plant height:
      1. Per-leaf:
         - When `keypoints_31` has >= 2 valid `petiole_1..6` keypoints (confidence >= 0.50),
           computes the polyline arc length along `lamina_base -> petiole_1 -> ... -> petiole_6`
           (plus continuation to `p_crown` if `petiole_6` terminates short of `p_crown`), taking
           `petiole_length_px = max(d_pose_polyline, geodesic_petiole_skeleton_px, ||p_base - p_crown||_2)`
         - `blade_length_px = ||p_apex - p_base||_2`
         - `petiole_to_blade_ratio = petiole_length_px / blade_length_px`
         - `petiole_length_mm`, `blade_length_mm`, and `petiole_width_mm` (from LM3 Module 13
           `petiole_width` or `width_1`/`width_2` when available, scaled via `px_per_mm = px_per_cm / 10.0`)
      2. Specimen-level outer-rosette aggregation (75th percentile across basal leaves per specimen):
         - `petiole_length_mm`, `blade_length_mm`, `petiole_to_blade_ratio`, `petiole_width_mm`
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
    if basal_leaves is None and kwargs.get("leaf_records") is not None:
        basal_leaves = kwargs["leaf_records"]
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
            "keypoints_31": kwargs.get("keypoints_31"),
            "petiole_width_mm": kwargs.get("petiole_width_mm"),
            "petiole_width_px": kwargs.get("petiole_width_px"),
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

        blade_px = float(p["blade_length_px"])
        if blade_px <= 0.0:
            blade_px = float(math.hypot(pa[0] - pb[0], pa[1] - pb[1]))

        base_poly_px = float(p.get("pose_petiole_polyline_px", 0.0) or 0.0)
        p_tip = p.get("p_petiole_tip")
        d_pose_polyline = base_poly_px
        if base_poly_px > 0.0 and p_tip is not None:
            tip_to_crown_px = float(math.hypot(p_tip[0] - c_pt[0], p_tip[1] - c_pt[1]))
            # Add continuation to p_crown if petiole_6 terminates short of p_crown along the inward direction
            dot_inward = (p_tip[0] - pb[0]) * (c_pt[0] - p_tip[0]) + (p_tip[1] - pb[1]) * (c_pt[1] - p_tip[1])
            if dot_inward > 0.0 and tip_to_crown_px > 1.0 and tip_to_crown_px <= 1.75 * max(blade_px, 1.0):
                d_pose_polyline = base_poly_px + tip_to_crown_px

        crown_fallback_px = crown_dist_px if (blade_px > 0.0 and crown_dist_px <= 1.75 * blade_px) else 0.0
        petiole_px = max(d_pose_polyline, geodesic_px, crown_fallback_px)
        pet_to_blade = (petiole_px / blade_px) if blade_px > 1e-6 else np.nan
        petiole_mm = (petiole_px / px_per_mm) if (px_per_mm is not None and px_per_mm > 0.0) else np.nan
        blade_mm = (blade_px / px_per_mm) if (px_per_mm is not None and px_per_mm > 0.0) else np.nan

        pw_px = p.get("petiole_width_px")
        pw_mm = p.get("petiole_width_mm")
        if pw_mm is None and pw_px is not None and px_per_mm is not None and px_per_mm > 0.0:
            pw_mm = float(pw_px) / px_per_mm
        elif pw_px is None and pw_mm is not None and px_per_mm is not None and px_per_mm > 0.0:
            pw_px = float(pw_mm) * px_per_mm

        leaf_metrics.append({
            "leaf_id": p["leaf_id"] if p["leaf_id"] is not None else (idx_l + 1),
            "plant_individual_id": pid,
            "p_base": pb,
            "p_apex": pa,
            "p_petiole_tip": p_tip,
            "p_crown": c_pt,
            "pose_petiole_polyline_px": round(float(d_pose_polyline), 4),
            "geodesic_petiole_skeleton_px": round(max(geodesic_px, float(d_pose_polyline)), 4),
            "crown_distance_px": round(crown_dist_px, 4),
            "petiole_length_px": round(petiole_px, 4),
            "blade_length_px": round(blade_px, 4),
            "petiole_to_blade_ratio": round(float(pet_to_blade), 4) if np.isfinite(pet_to_blade) else np.nan,
            "petiole_length_mm": round(float(petiole_mm), 3) if np.isfinite(petiole_mm) else np.nan,
            "blade_length_mm": round(float(blade_mm), 3) if np.isfinite(blade_mm) else np.nan,
            "petiole_width_px": round(float(pw_px), 4) if (pw_px is not None and np.isfinite(float(pw_px))) else np.nan,
            "petiole_width_mm": round(float(pw_mm), 3) if (pw_mm is not None and np.isfinite(float(pw_mm))) else np.nan,
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
            "petiole_width_px": _agg_p75(sub_r, "petiole_width_px", 4),
            "petiole_width_mm": _agg_p75(sub_r, "petiole_width_mm", 3),
        }

    petiole_px_p75 = _agg_p75(leaf_metrics, "petiole_length_px", 4)
    blade_px_p75 = _agg_p75(leaf_metrics, "blade_length_px", 4)
    petiole_mm_p75 = _agg_p75(leaf_metrics, "petiole_length_mm", 3)
    blade_mm_p75 = _agg_p75(leaf_metrics, "blade_length_mm", 3)
    pet_to_blade_p75 = _agg_p75(leaf_metrics, "petiole_to_blade_ratio", 4)
    petiole_width_px_p75 = _agg_p75(leaf_metrics, "petiole_width_px", 4)
    petiole_width_mm_p75 = _agg_p75(leaf_metrics, "petiole_width_mm", 3)

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
        "petiole_width_px": petiole_width_px_p75,
        "petiole_width_mm": petiole_width_mm_p75,
        "leaf_metrics": leaf_metrics,
        "by_plant": by_plant_summary,
    }


def classify_leaf_zone_anatomical(
    leaf_records: Sequence[Mapping[str, Any]],
    crowns_by_pid: Optional[Mapping[int, Tuple[float, float]]] = None,
    rosette_centroids_by_pid: Optional[Mapping[int, Optional[Tuple[float, float]]]] = None,
    basal_stats_by_pid: Optional[Mapping[int, Tuple[float, float, float, float]]] = None,
    capitula_detections: Optional[Sequence[Any]] = None,
    root_boxes: Optional[Sequence[Any]] = None,
    is_bent_stem: bool = False,
    sheet_shape: Optional[Tuple[int, int]] = None,
    **kwargs: Any,
) -> Dict[int, Dict[str, Any]]:
    """
    Post-segmentation anatomical leaf zone classifier (`basal` vs. `cauline`)
    operating per `plant_individual_id`.

    Never drops or rejects any leaf. Returns a dictionary mapping:
        int(leaf_id) -> {
            "leaf_zone": "basal" | "cauline",
            "leaf_zone_confidence": float,
            "dist_from_crown_norm": float,
        }
    """
    _ = sheet_shape
    if not leaf_records:
        return {}

    crowns_map: Dict[int, Tuple[float, float]] = {}
    if crowns_by_pid:
        for k, v in crowns_by_pid.items():
            if v is not None and len(v) >= 2 and np.isfinite(v[0]) and np.isfinite(v[1]):
                crowns_map[int(k)] = (float(v[0]), float(v[1]))
    explicit_crown = kwargs.get("p_crown")
    if explicit_crown is not None and len(explicit_crown) >= 2 and np.isfinite(explicit_crown[0]) and np.isfinite(explicit_crown[1]):
        crowns_map.setdefault(0, (float(explicit_crown[0]), float(explicit_crown[1])))

    rosette_map: Dict[int, Tuple[float, float]] = {}
    if rosette_centroids_by_pid:
        for k, v in rosette_centroids_by_pid.items():
            if v is not None and len(v) >= 2 and np.isfinite(v[0]) and np.isfinite(v[1]):
                rosette_map[int(k)] = (float(v[0]), float(v[1]))

    stats_map: Dict[int, Tuple[float, float, float, float]] = {}
    if basal_stats_by_pid:
        for k, v in basal_stats_by_pid.items():
            if v is not None and len(v) >= 4:
                stats_map[int(k)] = (float(v[0]), float(v[1]), float(v[2]), float(v[3]))

    # Extract all capitula candidate coordinates
    raw_caps = (
        capitula_detections
        or kwargs.get("capitula_boxes")
        or kwargs.get("capitula_pts")
        or []
    )
    cap_pts_all: List[Tuple[float, float]] = []
    for c_item in raw_caps:
        pt = _extract_xy_point(c_item)
        if pt is not None and np.isfinite(pt[0]) and np.isfinite(pt[1]):
            cap_pts_all.append((float(pt[0]), float(pt[1])))
    explicit_corymb = kwargs.get("p_corymb")
    if explicit_corymb is not None and len(explicit_corymb) >= 2:
        if np.isfinite(explicit_corymb[0]) and np.isfinite(explicit_corymb[1]):
            cap_pts_all.append((float(explicit_corymb[0]), float(explicit_corymb[1])))

    # Group leaf records by plant_individual_id
    records_by_pid: Dict[int, List[Tuple[int, Mapping[str, Any]]]] = {}
    for idx, rec in enumerate(leaf_records):
        raw_lid = rec.get("leaf_id")
        lid = int(raw_lid) if raw_lid is not None else (idx + 1)
        raw_pid = rec.get("plant_individual_id", 0)
        pid = int(raw_pid) if raw_pid is not None else 0
        records_by_pid.setdefault(pid, []).append((lid, rec))

    def _resolve_leaf_attachment_pt(rec: Mapping[str, Any]) -> Tuple[float, float]:
        for key in ("p_petiole_tip", "p_base_sheet", "p_base"):
            val = rec.get(key)
            if val is not None and isinstance(val, (tuple, list, np.ndarray)) and len(val) >= 2:
                if np.isfinite(val[0]) and np.isfinite(val[1]):
                    return (float(val[0]), float(val[1]))
        cx_val = rec.get("cx", rec.get("centroid_x"))
        cy_val = rec.get("cy", rec.get("centroid_y"))
        if cx_val is not None and cy_val is not None and np.isfinite(float(cx_val)) and np.isfinite(float(cy_val)):
            return (float(cx_val), float(cy_val))
        bbox = rec.get("bbox")
        if bbox is not None and len(bbox) == 4:
            ymin, xmin, ymax, xmax = [float(v) for v in bbox]
            return (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))
        return (0.0, 0.0)

    def _resolve_leaf_blade_length_px(rec: Mapping[str, Any]) -> float:
        for key in ("blade_length_px", "length_px"):
            val = rec.get(key)
            if val is not None and np.isfinite(float(val)) and float(val) > 0.0:
                return float(val)
        pb = rec.get("p_base_sheet") or rec.get("p_base")
        pa = rec.get("p_apex_sheet") or rec.get("p_apex")
        if (
            pb is not None
            and pa is not None
            and len(pb) >= 2
            and len(pa) >= 2
            and np.isfinite(pb[0])
            and np.isfinite(pa[0])
        ):
            d = float(math.hypot(float(pa[0]) - float(pb[0]), float(pa[1]) - float(pb[1])))
            if d > 0.0:
                return d
        area_px = float(rec.get("area_px", 0.0) or 0.0)
        if area_px > 0.0:
            return float(math.sqrt(area_px) * 1.4)
        return 40.0

    def _resolve_leaf_area(rec: Mapping[str, Any]) -> float:
        for key in ("blade_area_cm2", "area_px", "raw_area"):
            val = rec.get(key)
            if val is not None and np.isfinite(float(val)) and float(val) > 0.0:
                return float(val)
        return 1.0

    def _resolve_intrinsic_petiole_ratio(
        rec: Mapping[str, Any],
        blade_px: float,
        dist_crown_px: float,
    ) -> Tuple[float, bool]:
        """
        Evaluates leaf-intrinsic petiole ratio and petiole presence without
        conflating upper-stem crown distance with true petiole length.
        """
        has_pet_mask = bool(rec.get("has_petiole_mask", False))
        pet_mask = rec.get("petiole_mask")
        if isinstance(pet_mask, np.ndarray) and np.count_nonzero(pet_mask) > 0:
            has_pet_mask = True

        pose_poly = float(rec.get("pose_petiole_polyline_px", 0.0) or 0.0)
        geo_skel = float(rec.get("geodesic_petiole_skeleton_px", 0.0) or 0.0)
        p_tip = rec.get("p_petiole_tip")
        p_base = rec.get("p_base_sheet") or rec.get("p_base")
        tip_base_px = 0.0
        if (
            p_tip is not None
            and p_base is not None
            and len(p_tip) >= 2
            and len(p_base) >= 2
            and np.isfinite(p_tip[0])
            and np.isfinite(p_base[0])
        ):
            tip_base_px = float(math.hypot(float(p_tip[0]) - float(p_base[0]), float(p_tip[1]) - float(p_base[1])))

        direct_pet_px = max(pose_poly, geo_skel, tip_base_px)
        if direct_pet_px > 0.0 and blade_px > 1e-6:
            return float(direct_pet_px / blade_px), (has_pet_mask or direct_pet_px >= 0.12 * blade_px)

        raw_ratio = rec.get("petiole_to_blade_ratio")
        if raw_ratio is not None and np.isfinite(float(raw_ratio)) and float(raw_ratio) >= 0.0:
            r_val = float(raw_ratio)
            # Guard against legacy crown_distance_px / blade_length_px inflation on distant cauline leaves
            if dist_crown_px > 1.75 * max(blade_px, 1.0) and not has_pet_mask and direct_pet_px <= 0.0:
                crown_ratio_est = dist_crown_px / max(blade_px, 1.0)
                if abs(r_val - crown_ratio_est) <= 0.25 * max(crown_ratio_est, 1.0):
                    return 0.0, False
            return r_val, (has_pet_mask or r_val >= 0.12)

        pet_px = rec.get("petiole_length_px")
        if pet_px is not None and np.isfinite(float(pet_px)) and float(pet_px) > 0.0 and blade_px > 1e-6:
            p_val = float(pet_px)
            if dist_crown_px > 1.75 * max(blade_px, 1.0) and not has_pet_mask and abs(p_val - dist_crown_px) <= 5.0:
                return 0.0, False
            r_val = p_val / blade_px
            return r_val, (has_pet_mask or r_val >= 0.12)

        pet_mm = rec.get("petiole_length_mm")
        blade_mm = rec.get("blade_length_mm")
        if (
            pet_mm is not None
            and blade_mm is not None
            and np.isfinite(float(pet_mm))
            and np.isfinite(float(blade_mm))
            and float(blade_mm) > 1e-6
        ):
            r_val = float(pet_mm) / float(blade_mm)
            if dist_crown_px > 1.75 * max(blade_px, 1.0) and not has_pet_mask and r_val > 1.75:
                return 0.0, False
            return r_val, (has_pet_mask or r_val >= 0.12)

        return (0.18 if has_pet_mask else 0.0), has_pet_mask

    results: Dict[int, Dict[str, Any]] = {}

    for pid, pid_items in records_by_pid.items():
        att_pts = [_resolve_leaf_attachment_pt(rec) for _, rec in pid_items]
        blade_lengths = [_resolve_leaf_blade_length_px(rec) for _, rec in pid_items]
        areas = [_resolve_leaf_area(rec) for _, rec in pid_items]
        med_blade_px = float(np.median(blade_lengths)) if blade_lengths else 50.0
        max_area_pid = max(areas) if areas else 1.0

        # 1. Resolve basal crown / rosette origin c_pt
        c_pt: Optional[Tuple[float, float]] = (
            crowns_map.get(pid)
            or rosette_map.get(pid)
            or (stats_map[pid][:2] if pid in stats_map else None)
        )
        if c_pt is None:
            # Estimate from lowest / largest petiolate leaves on this plant
            poses_pid = []
            for idx_p, (lid, rec) in enumerate(pid_items):
                rec_d = dict(rec)
                rec_d.setdefault("leaf_id", lid if lid is not None else (idx_p + 1))
                parsed_p = _parse_basal_leaf_pose(rec_d)
                if parsed_p is not None:
                    poses_pid.append(parsed_p)
            if poses_pid:
                root_pts_est = [
                    pt
                    for r_item in (root_boxes or [])
                    for pt in (_extract_xy_point(r_item),)
                    if pt is not None
                ]
                c_pt = _estimate_single_plant_crown(poses_pid, root_pts_est)
            elif att_pts:
                # Fallback to lowest (largest y in image coords) attachment cluster
                ys = np.array([pt[1] for pt in att_pts], dtype=np.float64)
                xs = np.array([pt[0] for pt in att_pts], dtype=np.float64)
                y_thresh = float(np.percentile(ys, 65))
                low_mask = ys >= y_thresh
                c_pt = (float(np.median(xs[low_mask])), float(np.median(ys[low_mask])))
            else:
                c_pt = (0.0, 0.0)

        # Estimate characteristic basal rosette radius on this plant
        if pid in stats_map and np.isfinite(stats_map[pid][2]) and stats_map[pid][2] > 10.0:
            rosette_rad = float(stats_map[pid][2])
        else:
            dists_all = [float(math.hypot(pt[0] - c_pt[0], pt[1] - c_pt[1])) for pt in att_pts]
            low_dists = sorted(dists_all)[: max(1, int(math.ceil(0.6 * len(dists_all))))]
            est_rad = float(np.median(low_dists)) * 1.35 if low_dists else 0.0
            rosette_rad = max(est_rad, 1.35 * med_blade_px, 45.0)

        # 2. Resolve apical inflorescence corymb pole cap_pt
        cap_pt: Optional[Tuple[float, float]] = None
        if cap_pts_all:
            if len(records_by_pid) == 1:
                cap_pt = (
                    float(np.median([p[0] for p in cap_pts_all])),
                    float(np.median([p[1] for p in cap_pts_all])),
                )
            else:
                plant_xs = [pt[0] for pt in att_pts] + [c_pt[0]]
                med_x = float(np.median(plant_xs))
                span_x = max(float(np.ptp(plant_xs)), 2.5 * rosette_rad, 150.0)
                near_caps = [p for p in cap_pts_all if abs(p[0] - med_x) <= 1.5 * span_x]
                use_caps = near_caps if near_caps else cap_pts_all
                cap_pt = (
                    float(np.median([p[0] for p in use_caps])),
                    float(np.median([p[1] for p in use_caps])),
                )

        dists_from_crown = [
            float(math.hypot(pt[0] - c_pt[0], pt[1] - c_pt[1])) for pt in att_pts
        ]
        max_dist_on_pid = max(dists_from_crown) if dists_from_crown else 0.0

        scape_vec: Optional[Tuple[float, float]] = None
        scape_span = 0.0
        if cap_pt is not None:
            dx_s = float(cap_pt[0] - c_pt[0])
            dy_s = float(cap_pt[1] - c_pt[1])
            scape_span = float(math.hypot(dx_s, dy_s))
            if scape_span > 25.0:
                scape_vec = (dx_s, dy_s)

        for idx_item, (lid, rec) in enumerate(pid_items):
            p_att = att_pts[idx_item]
            blade_px = blade_lengths[idx_item]
            area_i = areas[idx_item]
            dist_crown_px = dists_from_crown[idx_item]

            raw_cy_rel = rec.get("cy_rel_plant_bbox")
            cy_rel: Optional[float] = (
                float(raw_cy_rel)
                if (raw_cy_rel is not None and np.isfinite(float(raw_cy_rel)))
                else None
            )

            # 3. Normalized progression from basal crown (0.0) to upper stem / corymb (1.0)
            if scape_vec is not None and scape_span > 25.0:
                vx = float(p_att[0] - c_pt[0])
                vy = float(p_att[1] - c_pt[1])
                axial_proj = (vx * scape_vec[0] + vy * scape_vec[1]) / (scape_span * scape_span)
                axial_proj = max(0.0, min(1.0, float(axial_proj)))
                radial_norm = min(1.0, dist_crown_px / max(scape_span, 1.0))
                dist_norm = max(0.0, min(1.0, 0.70 * axial_proj + 0.30 * radial_norm))
            elif len(pid_items) >= 2 and max_dist_on_pid > 1.0:
                denom = max(max_dist_on_pid, 2.25 * rosette_rad, 1.0)
                radial_norm = min(1.0, dist_crown_px / denom)
                if not is_bent_stem and cy_rel is not None:
                    vert_norm = max(0.0, min(1.0, 1.0 - cy_rel))
                    dist_norm = max(0.0, min(1.0, 0.72 * radial_norm + 0.28 * vert_norm))
                else:
                    dist_norm = radial_norm
            elif not is_bent_stem and cy_rel is not None and cy_rel < 0.45:
                dist_norm = max(0.0, min(1.0, 1.0 - cy_rel))
            else:
                dist_norm = min(1.0, dist_crown_px / max(2.5 * rosette_rad, 100.0))

            # 4. Continuous evidence scores in [0.0, 1.0]
            s_pos = max(0.0, min(1.0, 1.0 - dist_norm))
            pet_ratio, has_pet = _resolve_intrinsic_petiole_ratio(rec, blade_px, dist_crown_px)
            if has_pet:
                s_pet = min(1.0, max(0.55, pet_ratio / 0.35))
            else:
                s_pet = min(0.45, max(0.0, pet_ratio / 0.25))

            rel_area = max(0.0, min(1.0, area_i / max(max_area_pid, 1e-6)))
            s_area = min(1.0, math.sqrt(rel_area))

            p_basal = 0.50 * s_pos + 0.30 * s_pet + 0.20 * s_area

            # Anatomical cauline bract vs. basal rosette decision boundaries
            is_upper_zone = (dist_norm > 0.45) or (
                not is_bent_stem and cy_rel is not None and cy_rel < 0.40
            )
            is_sessile_or_detached_from_crown = (pet_ratio < 0.12) or (
                dist_crown_px > 1.8 * rosette_rad
            )
            is_upper_reduced_leaf = (dist_norm > 0.55) and (rel_area < 0.55)

            if (is_upper_zone and is_sessile_or_detached_from_crown) or is_upper_reduced_leaf:
                p_basal = min(
                    p_basal,
                    max(0.04, 0.42 * (1.0 - dist_norm) + 0.05 * rel_area),
                )
            elif (
                dist_norm <= 0.35
                or dist_crown_px <= 1.35 * rosette_rad
                or (dist_norm <= 0.45 and has_pet)
            ):
                p_basal = max(
                    p_basal,
                    min(0.98, 0.58 + 0.22 * s_pos + 0.12 * s_pet + 0.08 * s_area),
                )

            p_basal = max(0.0, min(1.0, float(p_basal)))
            zone_label = "basal" if p_basal >= 0.50 else "cauline"
            zone_conf = round(float(max(p_basal, 1.0 - p_basal)), 4)

            results[int(lid)] = {
                "leaf_zone": zone_label,
                "leaf_zone_confidence": zone_conf,
                "dist_from_crown_norm": round(float(dist_norm), 4),
            }

    return results
