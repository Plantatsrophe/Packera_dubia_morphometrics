#!/usr/bin/env python3
"""
===============================================================================
Module: organ_crown_metrics.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Multi-plant DBSCAN spatial clustering, 31-keypoint rosette crown / petiole /
    bent-scape linear trait estimation, and 4-case vascular-topology basal vs.
    cauline leaf zone classification on LeafMachine3 outputs.
===============================================================================
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from scripts.vision.lm3_compat import DBSCAN

def _mask_basal_apical_fallback(mask: np.ndarray) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Minimal contour endpoint fallback `(p_apex, p_base)` when bare 2D masks lack LM3 pose keypoints."""
    if mask is None or mask.size == 0:
        return (0, 0), (0, 0)
    m_u8 = (mask > 0).astype(np.uint8) * 255
    cnts, _ = cv2.findContours(m_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts:
        return (0, 0), (0, 0)
    pts = max(cnts, key=cv2.contourArea).reshape(-1, 2)
    base_idx = int(np.argmax(pts[:, 1]))
    pb = (int(pts[base_idx, 0]), int(pts[base_idx, 1]))
    apex_idx = int(np.argmax(np.hypot(pts[:, 0] - pb[0], pts[:, 1] - pb[1])))
    pa = (int(pts[apex_idx, 0]), int(pts[apex_idx, 1]))
    return pa, pb

__all__ = [
    "cluster_plant_individuals",
    "cluster_voucher_plants_dbscan",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
    "classify_leaf_zone_anatomical",
]


# =============================================================================
# 4. Adaptive DBSCAN Spatial Clustering
# =============================================================================

def cluster_plant_individuals(
    instances: List[Any],
    sheet_width: int,
    sheet_height: int,
    eps_ratio: float = 0.15,
    specimen_mask: Optional[np.ndarray] = None,
    root_boxes: Optional[Sequence[Any]] = None,
) -> List[Any]:
    """
    Stem-aware spatial clustering across detected leaves on a herbarium sheet.
    Starts from DBSCAN centroid neighbourhoods and merges clusters connected by:
      1. Shared Specimen_Segmenter foreground plant silhouette components,
      2. Vertical bolting-scape alignment within a single shoot column (unless
         anchored to distinct horizontally separated root crowns), or
      3. Splayed single-leaf outliers within rosette reach of a multi-leaf plant.
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
    max_dim = float(max(sheet_width, sheet_height, 1))
    sw = float(max(sheet_width, 1))
    sh = float(max(sheet_height, 1))
    eps = eps_ratio * max_dim

    clustering = DBSCAN(eps=eps, min_samples=1).fit(coords)
    raw_labels = [int(lbl) for lbl in clustering.labels_]
    uniq_init = sorted(set(raw_labels))
    init_map = {lbl: idx for idx, lbl in enumerate(uniq_init)}
    inst_clusters = [init_map[lbl] for lbl in raw_labels]
    n_cl = len(uniq_init)

    if n_cl > 1:
        parent = list(range(n_cl))

        def _find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def _union(a: int, b: int) -> None:
            ra, rb = _find(a), _find(b)
            if ra != rb:
                if ra < rb:
                    parent[rb] = ra
                else:
                    parent[ra] = rb

        # 1. Merge leaf clusters connected on the bridged Specimen_Segmenter mask
        if (
            isinstance(specimen_mask, np.ndarray)
            and specimen_mask.ndim == 2
            and specimen_mask.shape[0] == sheet_height
            and specimen_mask.shape[1] == sheet_width
            and np.any(specimen_mask > 0)
        ):
            k_bridge = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
            bridged = cv2.dilate((specimen_mask > 0).astype(np.uint8), k_bridge, iterations=1)
            _, cc_map = cv2.connectedComponents(bridged, connectivity=8)
            comp_owner: Dict[int, int] = {}
            for idx_i, (cx_i, cy_i) in enumerate(centroids):
                xi = int(np.clip(round(cx_i), 0, sheet_width - 1))
                yi = int(np.clip(round(cy_i), 0, sheet_height - 1))
                cid = int(cc_map[yi, xi])
                if cid == 0:
                    b_ymin, b_xmin, b_ymax, b_xmax = instances[idx_i].bbox
                    y0 = max(0, min(yi - 8, int(math.floor(b_ymin))))
                    y1 = min(sheet_height, max(yi + 9, int(math.ceil(b_ymax))))
                    x0 = max(0, min(xi - 8, int(math.floor(b_xmin))))
                    x1 = min(sheet_width, max(xi + 9, int(math.ceil(b_xmax))))
                    nz = cc_map[y0:y1, x0:x1]
                    nz = nz[nz > 0]
                    if nz.size > 0:
                        vals, cnts = np.unique(nz, return_counts=True)
                        cid = int(vals[int(np.argmax(cnts))])
                if cid > 0:
                    if cid in comp_owner:
                        _union(inst_clusters[idx_i], comp_owner[cid])
                    else:
                        comp_owner[cid] = inst_clusters[idx_i]

        # Extract optional root crown points to protect genuine side-by-side rooted plants
        root_pts: List[Tuple[float, float]] = []
        if root_boxes:
            for rb in root_boxes:
                pt_r = _extract_xy_point(rb)
                if pt_r is not None and np.isfinite(pt_r[0]) and np.isfinite(pt_r[1]):
                    root_pts.append((float(pt_r[0]), float(pt_r[1])))

        # 2. Merge vertically aligned bolting-scape clusters in the same shoot column
        cl_members: Dict[int, List[int]] = {}
        for idx_i, cid_i in enumerate(inst_clusters):
            cl_members.setdefault(_find(cid_i), []).append(idx_i)

        cl_keys = sorted(cl_members.keys())
        cl_med_xy: Dict[int, Tuple[float, float]] = {
            ck: (
                float(np.median(coords[idxs, 0])),
                float(np.median(coords[idxs, 1])),
            )
            for ck, idxs in cl_members.items()
        }

        def _nearest_root_idx(xy: Tuple[float, float]) -> Optional[int]:
            if len(root_pts) < 2:
                return None
            dists = [math.hypot(xy[0] - rx, xy[1] - ry) for rx, ry in root_pts]
            best_i = int(np.argmin(dists))
            return best_i if dists[best_i] <= 0.28 * max_dim else None

        cl_root_anchor = {ck: _nearest_root_idx(xy) for ck, xy in cl_med_xy.items()}
        for i_pos, ck_a in enumerate(cl_keys):
            for ck_b in cl_keys[i_pos + 1 :]:
                if _find(ck_a) == _find(ck_b):
                    continue
                ra_idx = cl_root_anchor.get(ck_a)
                rb_idx = cl_root_anchor.get(ck_b)
                if (
                    ra_idx is not None
                    and rb_idx is not None
                    and ra_idx != rb_idx
                    and abs(root_pts[ra_idx][0] - root_pts[rb_idx][0]) > 0.20 * sw
                ):
                    continue
                xa, ya = cl_med_xy[ck_a]
                xb, yb = cl_med_xy[ck_b]
                pts_a = coords[cl_members[ck_a]]
                pts_b = coords[cl_members[ck_b]]
                d_min_ab = float(
                    np.min(np.hypot(pts_a[:, 0:1] - pts_b[:, 0:1].T, pts_a[:, 1:2] - pts_b[:, 1:2].T))
                )
                diag_s = math.hypot(sw, sh)
                has_spec_mask = bool(
                    isinstance(specimen_mask, np.ndarray)
                    and specimen_mask.ndim == 2
                    and specimen_mask.shape[0] == sheet_height
                    and specimen_mask.shape[1] == sheet_width
                    and np.any(specimen_mask > 0)
                )
                has_bridge = False
                if has_spec_mask:
                    ia_x, ia_y = int(np.clip(round(xa), 0, sheet_width - 1)), int(np.clip(round(ya), 0, sheet_height - 1))
                    ib_x, ib_y = int(np.clip(round(xb), 0, sheet_width - 1)), int(np.clip(round(yb), 0, sheet_height - 1))
                    n_chk = 16
                    hits_sb = sum(
                        1
                        for s_i in range(n_chk + 1)
                        if specimen_mask[
                            int(np.clip(round(ia_y + (ib_y - ia_y) * (s_i / n_chk)), 0, sheet_height - 1)),
                            int(np.clip(round(ia_x + (ib_x - ia_x) * (s_i / n_chk)), 0, sheet_width - 1)),
                        ] > 0
                    )
                    has_bridge = hits_sb >= int(0.65 * (n_chk + 1))
                v_gap_close = abs(ya - yb) <= 0.14 * sh or d_min_ab <= 0.22 * diag_s
                if has_spec_mask:
                    if abs(xa - xb) <= 0.20 * sw and (v_gap_close or (has_bridge and abs(ya - yb) <= 0.55 * sh)):
                        _union(ck_a, ck_b)
                elif abs(xa - xb) <= 0.12 * sw and abs(ya - yb) <= 0.60 * sh:
                    _union(ck_a, ck_b)

        # 3. Absorb isolated 1-leaf outliers within rosette/petiole reach of a multi-leaf cluster
        merged_members: Dict[int, List[int]] = {}
        for idx_i, cid_i in enumerate(inst_clusters):
            merged_members.setdefault(_find(cid_i), []).append(idx_i)
        multi_clusters = [ck for ck, idxs in merged_members.items() if len(idxs) >= 2]
        if multi_clusters:
            absorb_rad = 0.34 * max_dim
            for ck_s, idxs_s in list(merged_members.items()):
                if len(idxs_s) != 1:
                    continue
                pt_s = coords[idxs_s[0]]
                best_mc: Optional[int] = None
                best_d = float("inf")
                for ck_m in multi_clusters:
                    pts_m = coords[merged_members[ck_m]]
                    d_m = float(np.min(np.hypot(pts_m[:, 0] - pt_s[0], pts_m[:, 1] - pt_s[1])))
                    if d_m < best_d:
                        best_d = d_m
                        best_mc = ck_m
                if best_mc is not None and best_d <= absorb_rad:
                    _union(ck_s, best_mc)

        ordered_roots: List[int] = []
        root_to_pid: Dict[int, int] = {}
        for idx_i, cid_i in enumerate(inst_clusters):
            r_id = _find(cid_i)
            if r_id not in root_to_pid:
                root_to_pid[r_id] = len(ordered_roots)
                ordered_roots.append(r_id)
            instances[idx_i].plant_individual_id = root_to_pid[r_id]
        return instances

    for inst, cid_i in zip(instances, inst_clusters):
        inst.plant_individual_id = cid_i

    return instances


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
            # Canonical repository bbox order is strictly (ymin, xmin, ymax, xmax) -> (cx, cy) = (0.5*(v1+v3), 0.5*(v0+v2))
            return (0.5 * (v1 + v3), 0.5 * (v0 + v2))
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
    leaf_zone: str = "basal"
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
        raw_pid = item.get("plant_individual_id")
        pid = int(raw_pid) if raw_pid is not None else 0
        leaf_id = item.get("leaf_id")
        leaf_zone = str(item.get("leaf_zone", "basal") or "basal").strip().lower()
        origin = item.get("mask_origin") or item.get("origin")
        if origin is not None and len(origin) == 2 and not item.get("coords_in_sheet_space", False):
            ox, oy = float(origin[0]), float(origin[1])
        local_kp_origin = item.get("local_mask_origin")

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
            pa_m, pb_m = _mask_basal_apical_fallback(item["mask"])
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
        pa_m, pb_m = _mask_basal_apical_fallback(item)
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
        raw_pid = getattr(item, "plant_individual_id", None)
        pid = int(raw_pid) if raw_pid is not None else 0
        leaf_id = getattr(item, "leaf_id", None)
        leaf_zone = str(getattr(item, "leaf_zone", "basal") or "basal").strip().lower()
        origin = getattr(item, "mask_origin", (0, 0)) or (0, 0)
        ox, oy = float(origin[0]), float(origin[1])
        keypoints_31 = getattr(item, "keypoints_31", None)
        kp_in_sheet_space = bool(getattr(item, "keypoints_in_sheet_space", False))
        petiole_width_mm = _safe_pos_float(getattr(item, "petiole_width_mm", None))
        petiole_width_px = _safe_pos_float(getattr(item, "petiole_width_px", getattr(item, "petiole_width", None)))
        pb = getattr(item, "p_base", None)
        pa = getattr(item, "p_apex", None)
        if (pb is None or pa is None) and getattr(item, "mask", None) is not None:
            pa, pb = _mask_basal_apical_fallback(getattr(item, "mask"))
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
        if kp_in_sheet_space:
            kp_ox, kp_oy = 0.0, 0.0
        elif (
            isinstance(item, dict)
            and item.get("local_mask_origin") is not None
            and len(item["local_mask_origin"]) == 2
        ):
            kp_ox = float(item["local_mask_origin"][0])
            kp_oy = float(item["local_mask_origin"][1])
        else:
            kp_ox, kp_oy = ox, oy
        if (
            kp_ox == 0.0
            and kp_oy == 0.0
            and not kp_in_sheet_space
            and isinstance(item, dict)
            and p_base is not None
        ):
            raw_lb = keypoints_31.get("lamina_base")
            if isinstance(raw_lb, (tuple, list, np.ndarray)) and len(raw_lb) >= 2:
                rx, ry = float(raw_lb[0]), float(raw_lb[1])
            elif isinstance(raw_lb, dict):
                rx, ry = float(raw_lb.get("x", 0.0)), float(raw_lb.get("y", 0.0))
            else:
                rx, ry = 0.0, 0.0
            if (
                item.get("p_base_sheet") is not None
                and (rx != 0.0 or ry != 0.0)
                and math.hypot(p_base[0] - rx, p_base[1] - ry) > 64.0
            ):
                kp_ox = p_base[0] - rx
                kp_oy = p_base[1] - ry

        def _extract_kp(name: str, min_conf: float = 0.25) -> Optional[Tuple[float, float, float]]:
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

        lb_kp = _extract_kp("lamina_base", min_conf=0.25)
        lt_kp = _extract_kp("lamina_tip", min_conf=0.25)
        if "midvein_0" in keypoints_31:
            mv_kps = [kp for kp in (_extract_kp(f"midvein_{i}", min_conf=0.25) for i in range(15)) if kp is not None]
        else:
            mv_kps = [kp for kp in (_extract_kp(f"midvein_{i}", min_conf=0.25) for i in range(1, 16)) if kp is not None]

        raw_pet_seq: List[Tuple[int, Tuple[float, float, float]]] = []
        if "petiole_0" in keypoints_31:
            for i in range(5):
                pkp = _extract_kp(f"petiole_{i}", min_conf=0.25)
                if pkp is not None:
                    raw_pet_seq.append((i + 1, pkp))
            ptip_kp = _extract_kp("petiole_tip", min_conf=0.25)
            if ptip_kp is not None:
                raw_pet_seq.append((6, ptip_kp))
        else:
            for i in range(1, 7):
                pkp = _extract_kp(f"petiole_{i}", min_conf=0.25)
                if pkp is not None:
                    raw_pet_seq.append((i, pkp))
            ptip_kp = _extract_kp("petiole_tip", min_conf=0.25)
            if ptip_kp is not None:
                raw_pet_seq.append((7, ptip_kp))

        pet_kps: List[Tuple[int, Tuple[float, float, float]]] = []
        for idx_p, pkp in raw_pet_seq:
            if not pet_kps or math.hypot(pkp[0] - pet_kps[-1][1][0], pkp[1] - pet_kps[-1][1][1]) > 1e-4:
                pet_kps.append((idx_p, pkp))

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
            # Distal crown-ward endpoint of petiole skeleton (petiole_tip / petiole_6 or highest valid index)
            p_petiole_tip = (pet_kps[-1][1][0], pet_kps[-1][1][1])

            # Polyline arc length along lamina_base -> petiole_0..4 -> petiole_tip when valid petiole keypoints exist
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

            # Align inward ray direction from lamina_tip/midvein through lamina_base and petiole keypoints
            spine_pts: List[Tuple[float, float]] = []
            if lt_kp is not None:
                spine_pts.append((lt_kp[0], lt_kp[1]))
            for mkp in mv_kps:
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

        # Extract Module 13 petiole_width or width_left/width_right / width_1/width_2 fallback when present
        if petiole_width_mm is None:
            petiole_width_mm = _safe_pos_float(keypoints_31.get("petiole_width_mm"))
        if petiole_width_px is None:
            petiole_width_px = _safe_pos_float(
                keypoints_31.get("petiole_width_px", keypoints_31.get("petiole_width"))
            )
        if petiole_width_px is None and petiole_width_mm is None:
            w1_kp = _extract_kp("width_left", min_conf=0.25) or _extract_kp("width_1", min_conf=0.25)
            w2_kp = _extract_kp("width_right", min_conf=0.25) or _extract_kp("width_2", min_conf=0.25)
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
        "leaf_zone": leaf_zone if leaf_zone in ("basal", "cauline") else "basal",
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
    other_plant_centers: Optional[Sequence[Tuple[float, float]]] = None,
) -> Tuple[float, float]:
    """
    Compute the caudex crown point `p_crown = (x_crown, y_crown)` for a single plant individual
    via weighted least-squares intersection of inward basal-leaf midrib/petiole rays anchored by
    LeafMachine3 `petiole_0..4` / `petiole_tip` distal endpoints, root detections, and/or the
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

    # Incorporate aligned LM3 ray directions and distal petiole_tip endpoints when available
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

    # Filter root points to those associated with this plant individual (and closer to this plant than to any other plant)
    near_roots: List[Tuple[float, float]] = []
    if root_pts:
        max_root_dist = max(4.0 * mean_len, 450.0)
        for rp in root_pts:
            d_self = float(np.min(np.hypot(bases[:, 0] - rp[0], bases[:, 1] - rp[1])))
            d_center = math.hypot(rp[0] - float(mean_base[0]), rp[1] - float(mean_base[1]))
            if other_plant_centers:
                d_other = min(math.hypot(rp[0] - oc[0], rp[1] - oc[1]) for oc in other_plant_centers)
                if d_center > 1.05 * d_other and d_self > d_other:
                    continue
            if min(d_self, d_center) <= max_root_dist:
                near_roots.append(rp)
        if not near_roots and len(root_pts) == 1 and not other_plant_centers:
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
                    cls_id = int(r_item[-1]) if isinstance(r_item[-1], (int, float, np.integer)) else 7
                    if cls_id not in (7, 8, 9):
                        continue
                elif isinstance(r_item, dict) and "class_id" in r_item and int(r_item["class_id"]) not in (7, 8, 9):
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
        pid_centers: Dict[int, Tuple[float, float]] = {
            pid_k: (
                float(np.median([p["p_base"][0] for p in pid_poses])),
                float(np.median([p["p_base"][1] for p in pid_poses])),
            )
            for pid_k, pid_poses in by_pid.items()
        }
        return {
            pid_k: _estimate_single_plant_crown(
                pid_poses,
                root_pts,
                other_plant_centers=[c for k, c in pid_centers.items() if k != pid_k],
            )
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
         - When `keypoints_31` has valid `petiole_0..4` / `petiole_tip` (or `petiole_1..6`)
           keypoints (`confidence >= 0.25` via `_parse_basal_leaf_pose`), computes the polyline
           arc length along `lamina_base -> petiole_0..4 -> petiole_tip`
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
    # Resolve scale conversion strictly when ruler-calibrated (`is_calibrated_mm` is not False)
    is_cal_arg = kwargs.get("is_calibrated_mm", True)
    eff_px_per_cm: Optional[float] = None
    if is_cal_arg is not False:
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

        is_caul_leaf = str(p.get("leaf_zone", "basal") or "basal").strip().lower() == "cauline"
        base_poly_px = float(p.get("pose_petiole_polyline_px", 0.0) or 0.0)
        p_tip = p.get("p_petiole_tip")
        d_pose_polyline = base_poly_px
        if not is_caul_leaf and base_poly_px > 0.0 and p_tip is not None:
            tip_to_crown_px = float(math.hypot(p_tip[0] - c_pt[0], p_tip[1] - c_pt[1]))
            # Add continuation to p_crown if petiole_6 terminates short of p_crown along the inward direction
            dot_inward = (p_tip[0] - pb[0]) * (c_pt[0] - p_tip[0]) + (p_tip[1] - pb[1]) * (c_pt[1] - p_tip[1])
            if dot_inward > 0.0 and tip_to_crown_px > 1.0 and tip_to_crown_px <= 1.75 * max(blade_px, 1.0):
                d_pose_polyline = base_poly_px + tip_to_crown_px

        crown_fallback_px = (
            crown_dist_px
            if (not is_caul_leaf and blade_px > 0.0 and crown_dist_px <= 1.75 * blade_px)
            else 0.0
        )
        direct_petiole_px = max(float(d_pose_polyline), float(geodesic_px))
        if direct_petiole_px > 0.0:
            petiole_px = direct_petiole_px
        elif is_caul_leaf:
            petiole_px = 0.0
        elif crown_fallback_px > 0.0:
            petiole_px = crown_fallback_px
        else:
            petiole_px = np.nan
        pet_to_blade = (petiole_px / blade_px) if (np.isfinite(petiole_px) and blade_px > 1e-6) else np.nan
        petiole_mm = (petiole_px / px_per_mm) if (np.isfinite(petiole_px) and px_per_mm is not None and px_per_mm > 0.0) else np.nan
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
            "leaf_zone": "cauline" if is_caul_leaf else "basal",
            "p_base": pb,
            "p_apex": pa,
            "p_petiole_tip": p_tip,
            "p_crown": c_pt,
            "pose_petiole_polyline_px": round(float(d_pose_polyline), 4),
            "geodesic_petiole_skeleton_px": round(max(geodesic_px, float(d_pose_polyline)), 4),
            "crown_distance_px": round(crown_dist_px, 4),
            "petiole_length_px": round(float(petiole_px), 4) if np.isfinite(petiole_px) else np.nan,
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
        sub_all = [r for r in leaf_metrics if int(r["plant_individual_id"]) == pid]
        sub_bas = [r for r in sub_all if r.get("leaf_zone") != "cauline"]
        sub_r = sub_bas if sub_bas else sub_all
        pid_cr = crowns_by_pid.get(pid, primary_crown)
        pid_pet_p75 = _agg_p75(sub_r, "petiole_length_px", 4)
        pid_bld_p75 = _agg_p75(sub_r, "blade_length_px", 4)
        pid_scape_px = np.nan
        pid_height_cm = np.nan
        pid_scape_ratio = np.nan
        pid_bent = False
        caps_by_pid = kwargs.get("capitula_by_pid")
        if isinstance(caps_by_pid, dict):
            pid_caps = caps_by_pid.get(pid) or []
            pid_cap_pts = [pt for c_it in pid_caps if (pt := _extract_xy_point(c_it, ref_pt=pid_cr)) is not None]
            if pid_cap_pts:
                pid_cor = (float(np.median([p[0] for p in pid_cap_pts])), float(np.median([p[1] for p in pid_cap_pts])))
                pid_scape_px = float(math.hypot(pid_cor[0] - pid_cr[0], pid_cor[1] - pid_cr[1]))
                if eff_px_per_cm is not None and eff_px_per_cm > 0.0:
                    pid_height_cm = round(float(pid_scape_px / eff_px_per_cm), 3)
                pid_ref = (pid_pet_p75 + pid_bld_p75) if (np.isfinite(pid_pet_p75) and np.isfinite(pid_bld_p75) and (pid_pet_p75 + pid_bld_p75) > 1e-6) else (pid_bld_p75 if (np.isfinite(pid_bld_p75) and pid_bld_p75 > 1e-6) else np.nan)
                if np.isfinite(pid_ref) and pid_ref > 1e-6:
                    pid_scape_ratio = round(float(pid_scape_px / pid_ref), 4)
        by_plant_summary[pid] = {
            "p_crown": pid_cr,
            "petiole_length_px": pid_pet_p75,
            "blade_length_px": pid_bld_p75,
            "petiole_length_mm": _agg_p75(sub_r, "petiole_length_mm", 3),
            "blade_length_mm": _agg_p75(sub_r, "blade_length_mm", 3),
            "petiole_to_blade_ratio": _agg_p75(sub_r, "petiole_to_blade_ratio", 4),
            "petiole_width_px": _agg_p75(sub_r, "petiole_width_px", 4),
            "petiole_width_mm": _agg_p75(sub_r, "petiole_width_mm", 3),
            "scape_length_px": round(float(pid_scape_px), 2) if np.isfinite(pid_scape_px) else np.nan,
            "plant_height_cm": pid_height_cm,
            "scape_to_basal_leaf_ratio": pid_scape_ratio,
            "is_bent_stem": pid_bent,
        }

    bas_all = [r for r in leaf_metrics if r.get("leaf_zone") != "cauline"]
    agg_pool = bas_all if bas_all else leaf_metrics
    petiole_px_p75 = _agg_p75(agg_pool, "petiole_length_px", 4)
    blade_px_p75 = _agg_p75(agg_pool, "blade_length_px", 4)
    petiole_mm_p75 = _agg_p75(agg_pool, "petiole_length_mm", 3)
    blade_mm_p75 = _agg_p75(agg_pool, "blade_length_mm", 3)
    pet_to_blade_p75 = _agg_p75(agg_pool, "petiole_to_blade_ratio", 4)
    petiole_width_px_p75 = _agg_p75(agg_pool, "petiole_width_px", 4)
    petiole_width_mm_p75 = _agg_p75(agg_pool, "petiole_width_mm", 3)

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
    specimen_mask: Optional[np.ndarray] = None,
    p_crown: Optional[Tuple[float, float]] = None,
    p_corymb: Optional[Tuple[float, float]] = None,
    plant_bbox: Optional[Tuple[float, float, float, float]] = None,
    **kwargs: Any,
) -> Dict[int, Dict[str, Any]]:
    """
    Post-segmentation anatomical leaf zone classifier (`basal` vs. `cauline`)
    operating per `plant_individual_id` via vascular attachment topology and
    organ-pole anchoring across all 4 organ-presence states:
      - Case 1 (`Roots` + `Capitula`): Geodesic/piecewise stem path s in [0, 1]
        from `7: Roots` (s = 0.0) to `2..6: Capitula` (s = 1.0), supporting
        upright, diagonal, inverted (`NCU00385199`), and V-bent (`PH00171299`) mounts.
      - Case 2 (`Roots` Only, No `Capitula`): `7: Roots` locks s = 0.0; leaves
        whose inward petiole rays converge to the root/caudex hub (`NCU00438219`)
        receive s_i ~ 0.0 (`basal`), while only leaves attached along an elongated
        bolting stem extending away from the root hub are `cauline`.
      - Case 3 (`Capitula` Only, No `Roots`): `2..6: Capitula` locks s = 1.0 and
        the opposite extremity defines s = 0.0. Distinguishes plants clipped at
        the soil line (basal convergence hub at s = 0.0) from clipped upper
        flowering stems (no basal convergence hub -> all internodal leaves `cauline`).
      - Case 4 (Neither `Roots` Nor `Capitula`): Distinguishes detached sterile
        basal rosettes/tufts (`NCSC00103660 pid=1`, `62226 pid=1`, radial ray
        convergence -> 100% `basal`) from elongated sterile stem segments.

    Never drops or rejects any leaf. Returns a dictionary mapping:
        int(leaf_id) -> {
            "leaf_zone": "basal" | "cauline",
            "leaf_zone_confidence": float,
            "dist_from_crown_norm": float,
        }
    """
    _ = (sheet_shape, plant_bbox)
    if not leaf_records:
        return {}

    crowns_map: Dict[int, Tuple[float, float]] = {}
    if crowns_by_pid:
        for k, v in crowns_by_pid.items():
            if v is not None and len(v) >= 2 and np.isfinite(v[0]) and np.isfinite(v[1]):
                crowns_map[int(k)] = (float(v[0]), float(v[1]))
    explicit_crown = p_crown if p_crown is not None else kwargs.get("p_crown")
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

    # Extract all root candidate coordinates on the sheet
    raw_roots = (
        root_boxes
        if root_boxes is not None
        else (
            kwargs.get("root_detections")
            or kwargs.get("roots")
            or kwargs.get("root_pts")
            or kwargs.get("stage1_roots")
            or []
        )
    )
    root_pts_all: List[Tuple[float, float]] = []
    if isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 2 and all(
        isinstance(v, (int, float, np.number)) for v in raw_roots
    ):
        root_pts_all.append((float(raw_roots[0]), float(raw_roots[1])))
    elif isinstance(raw_roots, (tuple, list, np.ndarray)) and len(raw_roots) == 4 and all(
        isinstance(v, (int, float, np.number)) for v in raw_roots
    ):
        pt = _extract_xy_point(raw_roots)
        if pt is not None and np.isfinite(pt[0]) and np.isfinite(pt[1]):
            root_pts_all.append((float(pt[0]), float(pt[1])))
    elif raw_roots is not None:
        for r_item in raw_roots:
            if isinstance(r_item, (tuple, list)) and len(r_item) >= 3 and isinstance(r_item[0], (tuple, list, np.ndarray)):
                cls_id = int(r_item[-1]) if isinstance(r_item[-1], (int, float, np.integer)) else 7
                if cls_id not in (7, 8, 9):
                    continue
            elif isinstance(r_item, dict) and "class_id" in r_item and int(r_item["class_id"]) not in (7, 8, 9):
                continue
            pt = _extract_xy_point(r_item)
            if pt is not None and np.isfinite(pt[0]) and np.isfinite(pt[1]):
                root_pts_all.append((float(pt[0]), float(pt[1])))

    # Extract all capitula candidate coordinates on the sheet
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
    explicit_corymb = p_corymb if p_corymb is not None else kwargs.get("p_corymb")
    if explicit_corymb is not None and len(explicit_corymb) >= 2:
        if np.isfinite(explicit_corymb[0]) and np.isfinite(explicit_corymb[1]):
            cap_pts_all.append((float(explicit_corymb[0]), float(explicit_corymb[1])))

    # Optional full-sheet Specimen_Segmenter connected components for plant-piece isolation and petiole tracing
    if specimen_mask is None:
        specimen_mask = kwargs.get("specimen_paper_cleaned_mask")
    comp_labels: Optional[np.ndarray] = None
    spec_bin: Optional[np.ndarray] = None
    if isinstance(specimen_mask, np.ndarray) and specimen_mask.ndim == 2 and specimen_mask.size > 0:
        spec_bin = (specimen_mask > 0).astype(np.uint8)
        if np.count_nonzero(spec_bin) > 0:
            k_bridge = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9))
            bridged = cv2.dilate(spec_bin, k_bridge, iterations=1)
            _, comp_labels = cv2.connectedComponents(bridged, connectivity=8)

    def _comp_id_at(pt: Tuple[float, float]) -> int:
        if comp_labels is None:
            return 0
        h_c, w_c = comp_labels.shape[:2]
        ix = int(round(float(pt[0])))
        iy = int(round(float(pt[1])))
        if 0 <= iy < h_c and 0 <= ix < w_c:
            cid = int(comp_labels[iy, ix])
            if cid > 0:
                return cid
        y0, y1 = max(0, iy - 6), min(h_c, iy + 7)
        x0, x1 = max(0, ix - 6), min(w_c, ix + 7)
        patch = comp_labels[y0:y1, x0:x1]
        nz = patch[patch > 0]
        if nz.size > 0:
            vals, counts = np.unique(nz, return_counts=True)
            return int(vals[int(np.argmax(counts))])
        return 0

    # Group leaf records by plant_individual_id
    records_by_pid: Dict[int, List[Tuple[int, Mapping[str, Any]]]] = {}
    for idx, rec in enumerate(leaf_records):
        raw_lid = rec.get("leaf_id")
        lid = int(raw_lid) if raw_lid is not None else (idx + 1)
        raw_pid = rec.get("plant_individual_id")
        pid = int(raw_pid) if raw_pid is not None else 0
        records_by_pid.setdefault(pid, []).append((lid, rec))

    def _resolve_leaf_base_and_apex(
        rec: Mapping[str, Any],
        lid: int,
    ) -> Dict[str, Any]:
        rec_d = dict(rec)
        rec_d.setdefault("leaf_id", lid)
        parsed = _parse_basal_leaf_pose(rec_d)
        p_base: Optional[Tuple[float, float]] = None
        p_apex: Optional[Tuple[float, float]] = None
        p_tip: Optional[Tuple[float, float]] = None
        ray_dir: Optional[Tuple[float, float]] = None
        has_pose_pet = False
        pose_poly_px = 0.0

        if parsed is not None:
            p_base = parsed.get("p_base")
            p_apex = parsed.get("p_apex")
            p_tip = parsed.get("p_petiole_tip")
            ray_dir = parsed.get("ray_dir")
            has_pose_pet = bool(parsed.get("has_lm3_petiole_pose", False))
            pose_poly_px = float(parsed.get("pose_petiole_polyline_px", 0.0) or 0.0)

        for key in ("p_base_sheet", "p_base"):
            val = rec.get(key)
            if p_base is None and val is not None and isinstance(val, (tuple, list, np.ndarray)) and len(val) >= 2:
                if np.isfinite(val[0]) and np.isfinite(val[1]):
                    p_base = (float(val[0]), float(val[1]))
        for key in ("p_apex_sheet", "p_apex"):
            val = rec.get(key)
            if p_apex is None and val is not None and isinstance(val, (tuple, list, np.ndarray)) and len(val) >= 2:
                if np.isfinite(val[0]) and np.isfinite(val[1]):
                    p_apex = (float(val[0]), float(val[1]))
        val_tip = rec.get("p_petiole_tip")
        if p_tip is None and val_tip is not None and isinstance(val_tip, (tuple, list, np.ndarray)) and len(val_tip) >= 2:
            if np.isfinite(val_tip[0]) and np.isfinite(val_tip[1]):
                p_tip = (float(val_tip[0]), float(val_tip[1]))
        val_rd = rec.get("ray_dir")
        if ray_dir is None and val_rd is not None and isinstance(val_rd, (tuple, list, np.ndarray)) and len(val_rd) >= 2:
            if np.isfinite(val_rd[0]) and np.isfinite(val_rd[1]):
                rn = math.hypot(float(val_rd[0]), float(val_rd[1]))
                if rn > 1e-6:
                    ray_dir = (float(val_rd[0]) / rn, float(val_rd[1]) / rn)

        if p_base is None:
            cx_val = rec.get("cx", rec.get("centroid_x"))
            cy_val = rec.get("cy", rec.get("centroid_y"))
            if cx_val is not None and cy_val is not None and np.isfinite(float(cx_val)) and np.isfinite(float(cy_val)):
                p_base = (float(cx_val), float(cy_val))
            else:
                bbox = rec.get("bbox")
                if bbox is not None and len(bbox) == 4:
                    ymin, xmin, ymax, xmax = [float(v) for v in bbox]
                    p_base = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))
                else:
                    p_base = (0.0, 0.0)

        if p_apex is None:
            p_apex = p_base
        if p_tip is None:
            p_tip = p_base

        if ray_dir is None:
            dx_tb = float(p_tip[0] - p_base[0])
            dy_tb = float(p_tip[1] - p_base[1])
            norm_tb = math.hypot(dx_tb, dy_tb)
            if norm_tb > 3.0:
                ray_dir = (dx_tb / norm_tb, dy_tb / norm_tb)
            else:
                dx_ab = float(p_base[0] - p_apex[0])
                dy_ab = float(p_base[1] - p_apex[1])
                norm_ab = math.hypot(dx_ab, dy_ab)
                if norm_ab > 1e-6:
                    ray_dir = (dx_ab / norm_ab, dy_ab / norm_ab)
                else:
                    ray_dir = (0.0, 1.0)

        return {
            "p_base": p_base,
            "p_apex": p_apex,
            "p_tip": p_tip,
            "ray_dir": ray_dir,
            "has_pose_pet": has_pose_pet,
            "pose_poly_px": pose_poly_px,
            "parsed_pose": parsed,
        }

    def _resolve_leaf_blade_length_px(rec: Mapping[str, Any], geom: Mapping[str, Any]) -> float:
        for key in ("blade_length_px", "length_px"):
            val = rec.get(key)
            if val is not None and np.isfinite(float(val)) and float(val) > 0.0:
                return float(val)
        pb = geom["p_base"]
        pa = geom["p_apex"]
        d = float(math.hypot(float(pa[0]) - float(pb[0]), float(pa[1]) - float(pb[1])))
        if d > 0.0:
            return d
        area_px = float(rec.get("area_px", 0.0) or 0.0)
        if area_px > 0.0:
            return float(math.sqrt(area_px) * 1.4)
        return 40.0

    def _resolve_intrinsic_petiole_ratio(
        rec: Mapping[str, Any],
        geom: Mapping[str, Any],
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

        pose_poly = max(
            float(rec.get("pose_petiole_polyline_px", 0.0) or 0.0),
            float(geom.get("pose_poly_px", 0.0) or 0.0),
        )
        geo_skel = float(rec.get("geodesic_petiole_skeleton_px", 0.0) or 0.0)
        p_tip = geom["p_tip"]
        p_base = geom["p_base"]
        tip_base_px = float(math.hypot(float(p_tip[0]) - float(p_base[0]), float(p_tip[1]) - float(p_base[1])))

        direct_pet_px = max(pose_poly, geo_skel, tip_base_px)
        if direct_pet_px > 0.0 and blade_px > 1e-6:
            return float(direct_pet_px / blade_px), (has_pet_mask or direct_pet_px >= 0.12 * blade_px)

        raw_ratio = rec.get("petiole_to_blade_ratio")
        if raw_ratio is not None and np.isfinite(float(raw_ratio)) and float(raw_ratio) >= 0.0:
            r_val = float(raw_ratio)
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

    # Precompute per-PID leaf anchor clouds and connected component IDs for strict multi-plant scoping
    pid_geoms: Dict[int, List[Dict[str, Any]]] = {}
    pid_anchor_pts: Dict[int, List[Tuple[float, float]]] = {}
    pid_comp_sets: Dict[int, set[int]] = {}
    pid_med_blade: Dict[int, float] = {}
    pid_span: Dict[int, float] = {}

    for pid, pid_items in records_by_pid.items():
        g_list = [_resolve_leaf_base_and_apex(rec, lid) for lid, rec in pid_items]
        pid_geoms[pid] = g_list
        b_lens = sorted(_resolve_leaf_blade_length_px(rec, g_list[i]) for i, (_, rec) in enumerate(pid_items))
        if len(b_lens) >= 4:
            med_b = float(np.median(b_lens[len(b_lens) // 2 :]))
        elif b_lens:
            med_b = float(np.median(b_lens))
        else:
            med_b = 50.0
        pid_med_blade[pid] = med_b
        pts_cloud: List[Tuple[float, float]] = []
        comp_set: set[int] = set()
        for g in g_list:
            pts_cloud.append(g["p_tip"])
            pts_cloud.append(g["p_base"])
            for pt_cand in (g["p_tip"], g["p_base"], g["p_apex"]):
                cid = _comp_id_at(pt_cand)
                if cid > 0:
                    comp_set.add(cid)
        pid_anchor_pts[pid] = pts_cloud
        pid_comp_sets[pid] = comp_set
        if pts_cloud:
            xs_c = [p[0] for p in pts_cloud]
            ys_c = [p[1] for p in pts_cloud]
            pid_span[pid] = max(float(np.ptp(xs_c)), float(np.ptp(ys_c)), 2.5 * med_b, 160.0)
        else:
            pid_span[pid] = max(2.5 * med_b, 160.0)

    def _min_dist_to_pid(pt: Tuple[float, float], target_pid: int) -> float:
        pts_t = pid_anchor_pts.get(target_pid, [])
        if not pts_t:
            return float("inf")
        return min(math.hypot(pt[0] - p[0], pt[1] - p[1]) for p in pts_t)

    def _scope_organ_pts_to_pid(
        organ_pts: Sequence[Tuple[float, float]],
        target_pid: int,
        max_mult_span: float,
        max_mult_blade: float,
        min_abs_px: float,
    ) -> List[Tuple[float, float]]:
        if not organ_pts:
            return []
        if len(records_by_pid) == 1:
            return list(organ_pts)

        scoped: List[Tuple[float, float]] = []
        max_allow = max(
            max_mult_span * pid_span.get(target_pid, 200.0),
            max_mult_blade * pid_med_blade.get(target_pid, 50.0),
            min_abs_px,
        )
        t_comps = pid_comp_sets.get(target_pid, set())
        other_pids = [k for k in records_by_pid.keys() if k != target_pid]

        for op in organ_pts:
            d_self = _min_dist_to_pid(op, target_pid)
            d_other = min((_min_dist_to_pid(op, ok) for ok in other_pids), default=float("inf"))
            op_cid = _comp_id_at(op)
            if op_cid > 0 and t_comps:
                in_self_comp = op_cid in t_comps
                in_other_comp = any(op_cid in pid_comp_sets.get(ok, set()) for ok in other_pids)
                if in_other_comp and not in_self_comp:
                    continue
                if in_self_comp and not in_other_comp and d_self <= 1.5 * max_allow:
                    scoped.append(op)
                    continue
            if d_self <= max_allow and d_self <= 1.05 * d_other:
                scoped.append(op)
        return scoped

    results: Dict[int, Dict[str, Any]] = {}

    for pid, pid_items in records_by_pid.items():
        g_list = pid_geoms[pid]
        blade_lengths = [_resolve_leaf_blade_length_px(rec, g_list[i]) for i, (_, rec) in enumerate(pid_items)]
        med_blade_px = pid_med_blade[pid]
        plant_span = pid_span[pid]

        roots_pid = _scope_organ_pts_to_pid(root_pts_all, pid, 1.65, 4.5, 380.0)
        caps_pid = _scope_organ_pts_to_pid(cap_pts_all, pid, 1.85, 5.5, 450.0)
        has_roots = len(roots_pid) > 0
        has_caps = len(caps_pid) > 0

        poses_pid = [g["parsed_pose"] for g in g_list if g.get("parsed_pose") is not None]
        other_centers = [
            (
                float(np.median([p[0] for p in pid_anchor_pts[ok]])),
                float(np.median([p[1] for p in pid_anchor_pts[ok]])),
            )
            for ok in records_by_pid.keys()
            if ok != pid and pid_anchor_pts.get(ok)
        ]

        # Compute ray-convergence hub p_hub from leaf inward rays on this plant piece
        p_hub = (
            _estimate_single_plant_crown(poses_pid, roots_pid, other_plant_centers=other_centers)
            if poses_pid
            else (
                (
                    float(np.median([g["p_tip"][0] for g in g_list])),
                    float(np.median([g["p_tip"][1] for g in g_list])),
                )
                if g_list
                else (0.0, 0.0)
            )
        )

        # Resolve apical corymb pole cap_pt when capitula belong to this plant piece
        cap_pt: Optional[Tuple[float, float]] = None
        if has_caps:
            cap_pt = (
                float(np.median([p[0] for p in caps_pid])),
                float(np.median([p[1] for p in caps_pid])),
            )

        # Resolve basal crown origin c_pt according to the 4 organ-pole presence states
        c_pt: Optional[Tuple[float, float]] = None
        if has_roots:
            # Cases 1 & 2: Roots explicitly anchor the basal pole (never overwrite top-mounted roots!)
            root_med = (
                float(np.median([p[0] for p in roots_pid])),
                float(np.median([p[1] for p in roots_pid])),
            )
            if pid in crowns_map and math.hypot(crowns_map[pid][0] - root_med[0], crowns_map[pid][1] - root_med[1]) <= max(2.0 * med_blade_px, 240.0):
                c_pt = crowns_map[pid]
            else:
                c_pt = _estimate_single_plant_crown(poses_pid, roots_pid, other_plant_centers=other_centers) if poses_pid else root_med
        elif has_caps and cap_pt is not None:
            # Case 3: Capitula lock s = 1.0; basal origin s = 0.0 is at the opposite extremity of the plant piece
            cand_c = crowns_map.get(pid) or rosette_map.get(pid) or (stats_map[pid][:2] if pid in stats_map else None)
            # Find leaves farthest from cap_pt to anchor the proximal base of the stem
            dists_to_cap = [
                math.hypot(g["p_tip"][0] - cap_pt[0], g["p_tip"][1] - cap_pt[1]) for g in g_list
            ]
            max_d_cap = max(dists_to_cap) if dists_to_cap else 0.0
            far_poses = [
                g["parsed_pose"]
                for idx_g, g in enumerate(g_list)
                if g.get("parsed_pose") is not None and dists_to_cap[idx_g] >= 0.55 * max_d_cap
            ]
            far_hub = _estimate_single_plant_crown(far_poses, [], other_plant_centers=other_centers) if far_poses else p_hub
            if cand_c is not None and math.hypot(cand_c[0] - cap_pt[0], cand_c[1] - cap_pt[1]) >= 0.45 * max_d_cap:
                c_pt = (float(cand_c[0]), float(cand_c[1]))
            else:
                c_pt = far_hub
        else:
            # Case 4: Neither Roots nor Capitula on this plant individual
            # Use ray-convergence hub p_hub (or explicit crown if consistent with p_hub)
            cand_c = crowns_map.get(pid)
            if cand_c is not None and math.hypot(cand_c[0] - p_hub[0], cand_c[1] - p_hub[1]) <= max(1.8 * med_blade_px, 200.0):
                c_pt = (float(cand_c[0]), float(cand_c[1]))
            else:
                c_pt = p_hub

        # Characteristic basal hub radius r_hub
        stat_rad = 0.0
        if pid in stats_map:
            st_tup = stats_map[pid]
            if len(st_tup) >= 4 and np.isfinite(st_tup[3]) and float(st_tup[3]) > 10.0:
                stat_rad = float(st_tup[3])
            elif len(st_tup) == 3 and np.isfinite(st_tup[2]) and 10.0 < float(st_tup[2]) <= 2500.0:
                stat_rad = float(st_tup[2])
        dists_tip_c = [float(math.hypot(g["p_tip"][0] - c_pt[0], g["p_tip"][1] - c_pt[1])) for g in g_list]
        low_dists = sorted(dists_tip_c)[: max(1, min(3, int(math.ceil(0.30 * len(dists_tip_c)))))]
        est_rad = float(np.median(low_dists)) * 1.20 if low_dists else 0.0
        if stat_rad > 10.0:
            rosette_rad = max(min(stat_rad, max(est_rad * 1.35, 2.35 * med_blade_px, 125.0)) if has_caps else stat_rad, 1.15 * med_blade_px, 45.0)
        else:
            rosette_rad = max(est_rad, 1.25 * med_blade_px, 45.0)
        if has_caps:
            r_hub = min(max(rosette_rad, 1.25 * med_blade_px, 50.0), max(2.35 * med_blade_px, 125.0))
        else:
            r_hub = max(rosette_rad, 1.35 * med_blade_px, 55.0)

        # Trace each leaf's inward petiole ray R_i(t) = p_tip + t * u_in toward c_pt and/or specimen_mask
        traced_att_pts: List[Tuple[float, float]] = []
        ray_crown_miss_px: List[float] = []
        ray_converges_to_crown: List[bool] = []
        intrinsic_pet_ratios: List[float] = []
        intrinsic_has_pet: List[bool] = []
        effective_pet_ratios: List[float] = []
        effective_has_pet: List[bool] = []
        pid_has_explicit_crown = bool(pid in crowns_map or pid in rosette_map)

        # Preliminary stem vector from c_pt to cap_pt (if cap_pt exists)
        prelim_scape_dx = float(cap_pt[0] - c_pt[0]) if cap_pt is not None else 0.0
        prelim_scape_dy = float(cap_pt[1] - c_pt[1]) if cap_pt is not None else 0.0
        prelim_scape_len = math.hypot(prelim_scape_dx, prelim_scape_dy)

        for idx_i, (_, rec_i) in enumerate(pid_items):
            g_i = g_list[idx_i]
            blade_px_i = blade_lengths[idx_i]
            p_base_i = g_i["p_base"]
            p_tip_i = g_i["p_tip"]
            u_in_i = g_i["ray_dir"]

            raw_d_crown = math.hypot(p_tip_i[0] - c_pt[0], p_tip_i[1] - c_pt[1])
            pet_ratio_i, has_pet_i = _resolve_intrinsic_petiole_ratio(rec_i, g_i, blade_px_i, raw_d_crown)
            intrinsic_pet_ratios.append(pet_ratio_i)
            intrinsic_has_pet.append(has_pet_i)

            # Project inward petiole ray forward from p_tip_i toward c_pt
            max_ray_reach = 3.8 * max(blade_px_i, med_blade_px, 50.0)
            vec_to_c_x = float(c_pt[0] - p_tip_i[0])
            vec_to_c_y = float(c_pt[1] - p_tip_i[1])
            t_proj = vec_to_c_x * u_in_i[0] + vec_to_c_y * u_in_i[1]
            t_clamped = max(0.0, min(max_ray_reach, t_proj))
            ray_closest_pt = (
                p_tip_i[0] + t_clamped * u_in_i[0],
                p_tip_i[1] + t_clamped * u_in_i[1],
            )
            d_perp_crown = math.hypot(ray_closest_pt[0] - c_pt[0], ray_closest_pt[1] - c_pt[1])

            # Check if specimen_mask supports vascular continuity along the ray from p_tip_i toward ray_closest_pt
            mask_ray_ok = True
            if spec_bin is not None and t_clamped > 12.0:
                h_s, w_s = spec_bin.shape[:2]
                n_steps = max(4, min(24, int(round(t_clamped / 12.0))))
                hits = 0
                valid_steps = 0
                for s_idx in range(1, n_steps + 1):
                    frac = float(s_idx) / float(n_steps)
                    sx = int(round(p_tip_i[0] + frac * t_clamped * u_in_i[0]))
                    sy = int(round(p_tip_i[1] + frac * t_clamped * u_in_i[1]))
                    if 0 <= sy < h_s and 0 <= sx < w_s:
                        valid_steps += 1
                        y0, y1 = max(0, sy - 5), min(h_s, sy + 6)
                        x0, x1 = max(0, sx - 5), min(w_s, sx + 6)
                        if np.any(spec_bin[y0:y1, x0:x1] > 0):
                            hits += 1
                if valid_steps >= 4 and (float(hits) / float(valid_steps)) < 0.35:
                    mask_ray_ok = False

            # Guard against aerial cauline leaves along a flowering scape whose downward angle points parallel to the stem
            is_aerial_stem_node = False
            if has_caps and prelim_scape_len > 60.0 and raw_d_crown > 1.35 * r_hub:
                s_prelim = ((p_tip_i[0] - c_pt[0]) * prelim_scape_dx + (p_tip_i[1] - c_pt[1]) * prelim_scape_dy) / (
                    prelim_scape_len * prelim_scape_len
                )
                s_cap_thresh = 0.28 if (has_pet_i and pet_ratio_i >= 0.40) else 0.18
                if s_prelim > s_cap_thresh:
                    is_aerial_stem_node = True

            converges = (
                raw_d_crown <= r_hub
                or (
                    t_proj >= -0.15 * blade_px_i
                    and raw_d_crown <= max_ray_reach
                    and d_perp_crown <= 1.05 * r_hub
                    and mask_ray_ok
                    and not is_aerial_stem_node
                    and (
                        not has_caps
                        or (has_pet_i and pet_ratio_i >= 0.35 and raw_d_crown <= 2.15 * r_hub)
                        or (has_pet_i and raw_d_crown <= 1.55 * r_hub)
                        or raw_d_crown <= 1.25 * r_hub
                    )
                )
            )

            if (has_roots or pid_has_explicit_crown) and converges and t_clamped > 0.0 and d_perp_crown < raw_d_crown:
                p_att_i = ray_closest_pt
                traced_pet_px = math.hypot(p_base_i[0] - p_att_i[0], p_base_i[1] - p_att_i[1])
                if blade_px_i > 1e-6 and traced_pet_px > pet_ratio_i * blade_px_i:
                    pet_ratio_i = float(traced_pet_px / blade_px_i)
                    if pet_ratio_i >= 0.15:
                        has_pet_i = True
            else:
                p_att_i = p_tip_i

            traced_att_pts.append(p_att_i)
            ray_crown_miss_px.append(d_perp_crown)
            ray_converges_to_crown.append(converges)
            effective_pet_ratios.append(pet_ratio_i)
            effective_has_pet.append(has_pet_i)

        n_leaves_pid = len(pid_items)
        # Second pass for dense long-petiolate basal tufts (e.g. NCU00438219) where an overlapping leaf's petiole is occluded
        conv_pet_vals = [
            effective_pet_ratios[k]
            for k in range(n_leaves_pid)
            if ray_converges_to_crown[k] and effective_pet_ratios[k] > 0.0
        ]
        if (
            n_leaves_pid >= 4
            and len(conv_pet_vals) >= max(3, n_leaves_pid - 1)
            and float(np.median(conv_pet_vals)) >= 0.75
        ):
            med_conv_pet = float(np.median(conv_pet_vals))
            for k in range(n_leaves_pid):
                if (
                    not ray_converges_to_crown[k]
                    and blade_lengths[k] >= 0.75 * med_blade_px
                    and (ray_crown_miss_px[k] <= 0.85 * r_hub or med_conv_pet >= 0.85)
                ):
                    ray_converges_to_crown[k] = True
                    traced_att_pts[k] = c_pt
                    traced_pet_k = math.hypot(g_list[k]["p_base"][0] - c_pt[0], g_list[k]["p_base"][1] - c_pt[1])
                    if blade_lengths[k] > 1e-6:
                        effective_pet_ratios[k] = max(effective_pet_ratios[k], float(traced_pet_k / blade_lengths[k]))
                        effective_has_pet[k] = True

        n_converging = sum(1 for c_flag in ray_converges_to_crown if c_flag)
        frac_converging = float(n_converging) / float(max(n_leaves_pid, 1))

        # PCA elongation and span of traced attachment points (for Case 2 & Case 4 stem vs. rosette topology)
        att_span_px = 0.0
        att_elongation = 1.0
        stem_axis_vec: Tuple[float, float] = (0.0, -1.0)
        if n_leaves_pid >= 2:
            att_arr = np.asarray(traced_att_pts, dtype=np.float64)
            centered = att_arr - np.mean(att_arr, axis=0, keepdims=True)
            cov = (centered.T @ centered) / float(max(n_leaves_pid - 1, 1))
            eigvals_pca, eigvecs_pca = np.linalg.eigh(cov)
            l_max = float(max(eigvals_pca[-1], 0.0))
            l_min = float(max(eigvals_pca[0], 1.0))
            att_elongation = math.sqrt(l_max / max(l_min, 1.0))
            axis_u = eigvecs_pca[:, -1]
            projs_pca = centered @ axis_u
            att_span_px = float(np.ptp(projs_pca))
            # Orient stem_axis_vec using outward leaf vectors (axillary stem leaves point toward s = 1.0)
            out_vecs = np.asarray(
                [(-g["ray_dir"][0], -g["ray_dir"][1]) for g in g_list],
                dtype=np.float64,
            )
            mean_out = np.mean(out_vecs, axis=0)
            if float(np.dot(mean_out, axis_u)) < 0.0:
                axis_u = -axis_u
            stem_axis_vec = (float(axis_u[0]), float(axis_u[1]))

        # Detect V-bent stem arch midpoint p_bend when is_bent_stem is set or when leaves arch far off [c_pt, cap_pt]
        p_bend: Optional[Tuple[float, float]] = kwargs.get("p_bend")
        if p_bend is None and cap_pt is not None and n_leaves_pid >= 2:
            chord_dx = float(cap_pt[0] - c_pt[0])
            chord_dy = float(cap_pt[1] - c_pt[1])
            chord_len = math.hypot(chord_dx, chord_dy)
            if is_bent_stem or chord_len < 0.55 * plant_span:
                # Find the attachment point with maximum perpendicular/arch excursion from the [c_pt, cap_pt] segment
                best_arch_d = 0.0
                best_arch_pt: Optional[Tuple[float, float]] = None
                for pt_a in traced_att_pts:
                    if chord_len > 20.0:
                        ux_c, uy_c = chord_dx / chord_len, chord_dy / chord_len
                        proj_c = (pt_a[0] - c_pt[0]) * ux_c + (pt_a[1] - c_pt[1]) * uy_c
                        cl_c = max(0.0, min(chord_len, proj_c))
                        d_arch = math.hypot(
                            pt_a[0] - (c_pt[0] + cl_c * ux_c),
                            pt_a[1] - (c_pt[1] + cl_c * uy_c),
                        )
                    else:
                        d_arch = math.hypot(pt_a[0] - c_pt[0], pt_a[1] - c_pt[1])
                    if d_arch > best_arch_d:
                        best_arch_d = d_arch
                        best_arch_pt = pt_a
                if best_arch_pt is not None and best_arch_d >= max(0.32 * plant_span, 1.5 * r_hub):
                    p_bend = best_arch_pt

        def _project_along_stem_polyline(
            pt: Tuple[float, float],
            p0: Tuple[float, float],
            p1: Tuple[float, float],
            pb_mid: Optional[Tuple[float, float]] = None,
        ) -> float:
            if pb_mid is not None:
                seg1_dx, seg1_dy = float(pb_mid[0] - p0[0]), float(pb_mid[1] - p0[1])
                seg2_dx, seg2_dy = float(p1[0] - pb_mid[0]), float(p1[1] - pb_mid[1])
                l1 = math.hypot(seg1_dx, seg1_dy)
                l2 = math.hypot(seg2_dx, seg2_dy)
                tot_l = l1 + l2
                if tot_l > 25.0 and l1 > 10.0 and l2 > 10.0:
                    t1 = max(0.0, min(1.0, ((pt[0] - p0[0]) * seg1_dx + (pt[1] - p0[1]) * seg1_dy) / (l1 * l1)))
                    near1 = (p0[0] + t1 * seg1_dx, p0[1] + t1 * seg1_dy)
                    d1 = math.hypot(pt[0] - near1[0], pt[1] - near1[1])

                    t2 = max(0.0, min(1.0, ((pt[0] - pb_mid[0]) * seg2_dx + (pt[1] - pb_mid[1]) * seg2_dy) / (l2 * l2)))
                    near2 = (pb_mid[0] + t2 * seg2_dx, pb_mid[1] + t2 * seg2_dy)
                    d2 = math.hypot(pt[0] - near2[0], pt[1] - near2[1])

                    s_poly = (t1 * l1) / tot_l if d1 <= d2 else (l1 + t2 * l2) / tot_l
                    return max(0.0, min(1.0, float(s_poly)))

            dx_s = float(p1[0] - p0[0])
            dy_s = float(p1[1] - p0[1])
            l_s = math.hypot(dx_s, dy_s)
            if l_s <= 20.0:
                return 0.0
            t_s = ((pt[0] - p0[0]) * dx_s + (pt[1] - p0[1]) * dy_s) / (l_s * l_s)
            return max(0.0, min(1.0, float(t_s)))

        dists_att_crown = [
            float(math.hypot(pt_a[0] - c_pt[0], pt_a[1] - c_pt[1])) for pt_a in traced_att_pts
        ]
        max_dist_att = max(dists_att_crown) if dists_att_crown else 0.0

        # Within-individual anatomical heteroblasty descriptors (Barkley 1962; Trock 2006)
        areas_pid = [
            float(rec["area_px"])
            if rec.get("area_px") is not None and np.isfinite(float(rec.get("area_px"))) and float(rec.get("area_px")) > 0.0
            else float("nan")
            for _, rec in pid_items
        ]
        finite_areas = [a for a in areas_pid if np.isfinite(a)]
        max_area_pid = max(finite_areas) if finite_areas else 0.0
        ar_vals_pid = [
            float(rec["aspect_ratio"])
            if rec.get("aspect_ratio") is not None and np.isfinite(float(rec.get("aspect_ratio"))) and float(rec.get("aspect_ratio")) > 0.0
            else float("nan")
            for _, rec in pid_items
        ]
        rel_areas_pid: List[float] = []
        basal_morph_scores: List[float] = []
        has_morph_covs = len(finite_areas) > 0 or any(np.isfinite(v) for v in ar_vals_pid)
        for i_m, (_, rec_m) in enumerate(pid_items):
            ra_i = (areas_pid[i_m] / max_area_pid) if (max_area_pid > 0.0 and np.isfinite(areas_pid[i_m])) else 0.65
            rel_areas_pid.append(ra_i)
            ra_sc = float(np.clip((ra_i - 0.18) / 0.52, 0.0, 1.0)) if (len(finite_areas) >= 3 and max_area_pid > 0.0) else 0.60
            ar_v = ar_vals_pid[i_m]
            ar_sc = float(np.clip(1.0 - (ar_v - 1.35) / 2.25, 0.0, 1.0)) if np.isfinite(ar_v) else 0.55
            pw_v = float(rec_m.get("petiole_width_px") or 0.0)
            bw_v = float(rec_m.get("blade_width_px") or 0.0)
            slender_sc = float(np.clip(1.0 - ((pw_v / max(bw_v, 1.0)) - 0.18) / 0.30, 0.10, 1.0)) if (pw_v > 0.0 and bw_v > 1.0) else 0.55
            pet_cue = min(1.0, 0.38 + 0.65 * effective_pet_ratios[i_m]) if effective_has_pet[i_m] else 0.15
            bms_i = 0.30 * ra_sc + 0.25 * ar_sc + 0.15 * slender_sc + 0.30 * pet_cue
            basal_morph_scores.append(float(np.clip(bms_i, 0.0, 1.0)))

        # Evaluate Case 3B: Clipped upper flowering stem fragment (has_caps, no roots, no genuine petiolate rosette hub)
        is_clipped_upper_flowering_stem = False
        if has_caps and not has_roots:
            n_prox_pet = sum(
                1
                for i in range(n_leaves_pid)
                if (intrinsic_has_pet[i] and intrinsic_pet_ratios[i] >= 0.25 and basal_morph_scores[i] >= 0.48)
                or (intrinsic_pet_ratios[i] >= 0.30 and dists_att_crown[i] <= 1.85 * r_hub)
            )
            if n_prox_pet == 0:
                is_clipped_upper_flowering_stem = True

        # Evaluate Case 4A: Detached sterile basal rosette / basal tuft (no roots, no capitula, radial/fasciculate hub)
        is_detached_sterile_rosette = False
        n_pet_pid = sum(1 for hp in effective_has_pet if hp)
        finite_ar_pid = [v for v in ar_vals_pid if np.isfinite(v)]
        max_ar_pid = max(finite_ar_pid) if finite_ar_pid else 1.8
        mean_bms_pid = float(np.mean(basal_morph_scores)) if basal_morph_scores else 0.60
        if not has_roots and not has_caps:
            is_sessile_stem_fragment = (
                n_pet_pid == 0
                and max(effective_pet_ratios, default=0.0) < 0.16
                and (
                    (n_leaves_pid >= 2 and att_span_px > 1.15 * med_blade_px)
                    or (has_morph_covs and max_ar_pid >= 2.45 and mean_bms_pid < 0.46)
                )
            )
            if not is_sessile_stem_fragment and (
                (n_leaves_pid <= 2 and max_dist_att <= 2.2 * r_hub and (n_pet_pid >= 1 or att_span_px <= 1.15 * med_blade_px))
                or (frac_converging >= 0.85 and max_dist_att <= 1.85 * r_hub and (n_pet_pid >= 1 or max_ar_pid < 2.4))
                or (att_span_px < 1.85 * med_blade_px and (n_pet_pid >= 1 or frac_converging >= 0.65))
                or (att_elongation < 2.4 and max_dist_att <= 1.85 * r_hub and n_pet_pid >= 1)
            ):
                is_detached_sterile_rosette = True

        # Evaluate Case 2 sterile basal tuft (roots present, no capitula, all or nearly all leaves converge to root hub)
        is_rooted_sterile_tuft = False
        if has_roots and not has_caps:
            if frac_converging >= 0.70 or att_span_px < 2.2 * r_hub:
                is_rooted_sterile_tuft = True

        for idx_item, (lid, rec) in enumerate(pid_items):
            p_att = traced_att_pts[idx_item]
            dist_crown_px = dists_att_crown[idx_item]
            converges_i = ray_converges_to_crown[idx_item]
            pet_ratio = effective_pet_ratios[idx_item]
            has_pet = effective_has_pet[idx_item]
            bms_item = basal_morph_scores[idx_item]
            rel_a_item = rel_areas_pid[idx_item]
            ar_item = ar_vals_pid[idx_item]
            morph_adj = 0.10 * (bms_item - 0.50) if has_morph_covs else 0.0

            # Compute normalized progression s_i = dist_norm in [0.0, 1.0] and p_basal in [0.0, 1.0]
            if is_detached_sterile_rosette:
                # Case 4A: Detached sterile rosette / basal tuft -> 100% basal regardless of sheet y-position
                dist_norm = min(0.24, dist_crown_px / max(3.5 * r_hub, 150.0))
                p_basal = max(0.76, min(0.98, 0.92 - 0.35 * dist_norm + (0.06 if has_pet else 0.0)))
            elif is_rooted_sterile_tuft and converges_i:
                # Case 2 (tuft): Erect or spreading basal tuft anchored by Roots (`NCU00438219`)
                dist_norm = min(0.22, dist_crown_px / max(3.5 * r_hub, 150.0))
                p_basal = max(0.78, min(0.98, 0.94 - 0.35 * dist_norm + (0.04 if has_pet else 0.0)))
            elif is_clipped_upper_flowering_stem:
                # Case 3B: Clipped upper flowering stem fragment with Capitula and no basal rosette hub
                if cap_pt is not None:
                    s_stem = _project_along_stem_polyline(p_att, c_pt, cap_pt, p_bend)
                    dist_norm = max(0.55, min(1.0, 0.55 + 0.45 * s_stem))
                else:
                    dist_norm = 0.75
                p_basal = max(0.04, min(0.32, 0.35 * (1.0 - dist_norm)))
            elif has_caps and cap_pt is not None:
                # Cases 1 & 3A: Stem path from basal crown c_pt (s = 0.0) to apical corymb cap_pt (s = 1.0)
                p_tip_raw = g_list[idx_item]["p_tip"]
                raw_d_tip = float(math.hypot(p_tip_raw[0] - c_pt[0], p_tip_raw[1] - c_pt[1]))
                s_stem_att = _project_along_stem_polyline(p_att, c_pt, cap_pt, p_bend)
                s_stem_tip = _project_along_stem_polyline(p_tip_raw, c_pt, cap_pt, p_bend)
                s_stem = max(s_stem_att, 0.85 * s_stem_tip)
                scape_span = max(math.hypot(cap_pt[0] - c_pt[0], cap_pt[1] - c_pt[1]), 1.5 * r_hub, 100.0)
                radial_norm = min(1.0, raw_d_tip / scape_span)
                dist_norm = max(0.0, min(1.0, 0.78 * s_stem_tip + 0.22 * radial_norm))

                is_long_pet_rosette = (
                    has_pet
                    and pet_ratio >= 0.35
                    and s_stem_tip <= 0.32
                    and (converges_i or raw_d_tip <= 2.15 * r_hub)
                    and (not has_morph_covs or (bms_item >= 0.46 and rel_a_item >= 0.32))
                )
                is_sessile_narrow_cauline = (
                    not has_pet
                    and pet_ratio < 0.15
                    and has_morph_covs
                    and np.isfinite(ar_item)
                    and ar_item >= 2.55
                    and rel_a_item < 0.45
                )

                if is_sessile_narrow_cauline and s_stem_tip > 0.08:
                    p_basal = max(0.05, min(0.36, 0.38 * (1.0 - dist_norm) + morph_adj))
                elif converges_i and s_stem_tip <= 0.16 and (raw_d_tip <= 1.15 * r_hub or (has_pet and pet_ratio >= 0.22 and raw_d_tip <= 1.55 * r_hub)):
                    dist_norm = min(dist_norm, 0.18)
                    p_basal = max(0.68, min(0.98, 0.90 - 0.55 * dist_norm + min(0.10, 0.15 * pet_ratio) + 0.5 * morph_adj))
                elif is_long_pet_rosette:
                    dist_norm = min(dist_norm, 0.24)
                    p_basal = max(0.64, min(0.95, 0.82 - 0.35 * dist_norm + min(0.10, 0.12 * pet_ratio) + 0.5 * morph_adj))
                elif s_stem_tip <= 0.12 and raw_d_tip <= 1.05 * r_hub:
                    p_basal = max(0.62, min(0.95, 0.84 - 0.60 * dist_norm + (0.08 if has_pet else 0.0) + 0.5 * morph_adj))
                elif s_stem_tip <= 0.18 and has_pet and pet_ratio >= 0.25 and raw_d_tip <= 1.45 * r_hub:
                    p_basal = max(0.56, min(0.88, 0.78 - 0.55 * dist_norm + min(0.12, 0.15 * pet_ratio) + 0.5 * morph_adj))
                else:
                    p_basal = max(0.04, min(0.46, 0.40 * (1.0 - dist_norm) + (0.04 if (has_pet and s_stem_tip < 0.28) else 0.0) + morph_adj))
            elif has_roots and not has_caps:
                # Case 2 (bolting/decapitated stem with roots): leaves outside root hub along stem axis
                denom = max(max_dist_att, 2.5 * r_hub, 120.0)
                dist_norm = max(0.0, min(1.0, dist_crown_px / denom))
                if converges_i or dist_crown_px <= 1.25 * r_hub or (dist_norm <= 0.42 and has_pet):
                    dist_norm = min(dist_norm, 0.30)
                    p_basal = max(0.65, min(0.98, 0.88 - 0.45 * dist_norm + (0.08 if has_pet else 0.0)))
                else:
                    p_basal = max(0.06, min(0.46, 0.44 * (1.0 - dist_norm)))
            else:
                # Case 4B: Elongated sterile leafy stem segment (no roots, no capitula, elongated internodal axis)
                denom_c4 = max(max_dist_att, 2.5 * r_hub, 120.0)
                radial_c4 = max(0.0, min(1.0, dist_crown_px / denom_c4))
                proj_vals = [
                    (pt_a[0] - c_pt[0]) * stem_axis_vec[0] + (pt_a[1] - c_pt[1]) * stem_axis_vec[1]
                    for pt_a in traced_att_pts
                ]
                p_min = min(proj_vals) if proj_vals else 0.0
                p_ptp = max(max(proj_vals) - p_min, 2.5 * r_hub, 120.0) if proj_vals else 120.0
                s_axis = max(0.0, min(1.0, (proj_vals[idx_item] - p_min) / p_ptp))
                dist_norm = max(radial_c4, s_axis) if n_converging >= 2 else s_axis
                if converges_i or dist_crown_px <= 1.15 * r_hub or (dist_norm <= 0.32 and has_pet):
                    dist_norm = min(dist_norm, 0.28)
                    p_basal = max(0.64, min(0.96, 0.86 - 0.45 * dist_norm + (0.06 if has_pet else 0.0)))
                else:
                    p_basal = max(0.06, min(0.46, 0.44 * (1.0 - dist_norm)))

            p_basal = max(0.0, min(1.0, float(p_basal)))
            zone_label = "basal" if p_basal >= 0.50 else "cauline"
            zone_conf = round(float(max(p_basal, 1.0 - p_basal)), 4)

            results[int(lid)] = {
                "leaf_zone": zone_label,
                "leaf_zone_confidence": zone_conf,
                "dist_from_crown_norm": round(float(dist_norm), 4),
            }

    return results
