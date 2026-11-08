#!/usr/bin/env python3
"""
===============================================================================
Module: leaf_gatekeeper.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Botanical leaf quality gatekeeping for LeafMachine3 (YOLO26x / UNet++ v9),
    3-class (lamina, petiole, hole) mask routing, 31-keypoint midvein trace
    fold detection/reflection, UNet++ v9 whole-specimen rachis bridging,
    boundary tissue continuity verification, cauline bract rejection, and
    2-tier silhouette routing (Tier_1_Direct / Tier_2_Reflected).
===============================================================================
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import cv2
import numpy as np

from scripts.vision.leaf_geometry_utils import (
    compute_geometric_metrics,
    extract_tier2_reflected,
    resample_normalized_contour,
)

__all__ = [
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
    "TIER_1_DIRECT",
    "TIER_2_REFLECTED",
    "CANONICAL_TIER_MAP",
    "normalize_assigned_tier",
    "compute_sheet_adaptive_min_blade_area_cm2",
    "is_cauline_or_stem_contaminant",
    "LeafRoutingResult",
    "classify_and_route_leaf",
    "evaluate_leaf_quality",
]


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
         so 1-3 px LeafMachine3 (YOLO26x) serration margin conservatism and boundary anti-aliasing are not
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


def _normalize_keypoints_31_input(
    keypoints_31: Any,
) -> Optional[Dict[str, Tuple[float, float, float]]]:
    """Normalizes a dict or (31, 3) ndarray of LeafMachine3 keypoints into a canonical dict."""
    if keypoints_31 is None:
        return None
    if isinstance(keypoints_31, dict):
        return keypoints_31 if len(keypoints_31) > 0 else None
    if isinstance(keypoints_31, np.ndarray) and keypoints_31.ndim == 2 and keypoints_31.shape[0] >= 8 and keypoints_31.shape[1] >= 2:
        kp_dict: Dict[str, Tuple[float, float, float]] = {}
        def _row_to_tup(row: np.ndarray) -> Tuple[float, float, float]:
            return (float(row[0]), float(row[1]), float(row[2]) if len(row) >= 3 else 1.0)
        kp_dict["lamina_tip"] = _row_to_tup(keypoints_31[0])
        kp_dict["lamina_base"] = _row_to_tup(keypoints_31[1])
        for i in range(min(6, keypoints_31.shape[0] - 2)):
            kp_dict[f"petiole_{i + 1}"] = _row_to_tup(keypoints_31[2 + i])
        for i in range(min(15, max(0, keypoints_31.shape[0] - 8))):
            kp_dict[f"midvein_{i + 1}"] = _row_to_tup(keypoints_31[8 + i])
        return kp_dict
    return None


def _extract_lm3_midvein_trace(
    keypoints_31: Optional[Any],
    min_points: int = 3,
    min_mean_conf: float = 0.50,
) -> Dict[str, Any]:
    """
    Parses LeafMachine3 (YOLO26x) 31-keypoint midvein trace (`midvein_1`..`midvein_15`,
    plus `lamina_base` and `lamina_tip`) and fits a least-squares anatomical midrib line
    anchored from `lamina_base` -> `lamina_tip` when >= `min_points` valid points exist
    with mean confidence >= `min_mean_conf`.
    """
    default_trace: Dict[str, Any] = {
        "is_valid": False,
        "points": np.zeros((0, 2), dtype=np.float64),
        "confidences": np.zeros((0,), dtype=np.float64),
        "mean_conf": 0.0,
        "vx": 1.0,
        "vy": 0.0,
        "x0": 0.0,
        "y0": 0.0,
        "p_base": None,
        "p_apex": None,
        "lamina_base": None,
        "lamina_tip": None,
    }
    keypoints_31 = _normalize_keypoints_31_input(keypoints_31)
    if keypoints_31 is None:
        return default_trace

    def _parse_kp(val: Any) -> Optional[Tuple[float, float, float]]:
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
        if c <= 0.05 or (x == 0.0 and y == 0.0):
            return None
        return (x, y, c)

    lb_kp = _parse_kp(keypoints_31.get("lamina_base"))
    lt_kp = _parse_kp(keypoints_31.get("lamina_tip"))
    if lb_kp is not None and lb_kp[2] >= 0.50:
        default_trace["lamina_base"] = (lb_kp[0], lb_kp[1])
    if lt_kp is not None and lt_kp[2] >= 0.50:
        default_trace["lamina_tip"] = (lt_kp[0], lt_kp[1])

    trace_pts: List[Tuple[float, float]] = []
    trace_confs: List[float] = []
    if lb_kp is not None and lb_kp[2] >= 0.25:
        trace_pts.append((lb_kp[0], lb_kp[1]))
        trace_confs.append(lb_kp[2])

    for idx_m in range(1, 16):
        mkp = _parse_kp(keypoints_31.get(f"midvein_{idx_m}"))
        if mkp is not None and mkp[2] >= 0.25:
            trace_pts.append((mkp[0], mkp[1]))
            trace_confs.append(mkp[2])

    if lt_kp is not None and lt_kp[2] >= 0.25:
        trace_pts.append((lt_kp[0], lt_kp[1]))
        trace_confs.append(lt_kp[2])

    if len(trace_pts) < min_points:
        return default_trace

    confs_arr = np.asarray(trace_confs, dtype=np.float64)
    mean_conf = float(np.mean(confs_arr))
    default_trace["mean_conf"] = round(mean_conf, 4)
    if mean_conf < min_mean_conf:
        return default_trace

    pts_arr = np.asarray(trace_pts, dtype=np.float64)
    line_params = cv2.fitLine(pts_arr.astype(np.float32), cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]

    # Orient (vx, vy) from lamina_base toward lamina_tip (or first -> last trace point)
    if default_trace["lamina_base"] is not None and default_trace["lamina_tip"] is not None:
        ref_dx = default_trace["lamina_tip"][0] - default_trace["lamina_base"][0]
        ref_dy = default_trace["lamina_tip"][1] - default_trace["lamina_base"][1]
    elif default_trace["lamina_base"] is not None:
        ref_dx = pts_arr[-1, 0] - default_trace["lamina_base"][0]
        ref_dy = pts_arr[-1, 1] - default_trace["lamina_base"][1]
    elif default_trace["lamina_tip"] is not None:
        ref_dx = default_trace["lamina_tip"][0] - pts_arr[0, 0]
        ref_dy = default_trace["lamina_tip"][1] - pts_arr[0, 1]
    else:
        ref_dx = pts_arr[-1, 0] - pts_arr[0, 0]
        ref_dy = pts_arr[-1, 1] - pts_arr[0, 1]

    if (ref_dx * vx + ref_dy * vy) < 0.0:
        vx, vy = -vx, -vy

    proj_t = (pts_arr[:, 0] - x0) * vx + (pts_arr[:, 1] - y0) * vy
    t_base = float(np.min(proj_t))
    t_apex = float(np.max(proj_t))
    if default_trace["lamina_base"] is not None:
        tb_anch = (default_trace["lamina_base"][0] - x0) * vx + (default_trace["lamina_base"][1] - y0) * vy
        t_base = min(t_base, float(tb_anch))
    if default_trace["lamina_tip"] is not None:
        ta_anch = (default_trace["lamina_tip"][0] - x0) * vx + (default_trace["lamina_tip"][1] - y0) * vy
        t_apex = max(t_apex, float(ta_anch))

    if abs(t_apex - t_base) < 1e-4:
        return default_trace

    p_base_fit = (int(round(x0 + t_base * vx)), int(round(y0 + t_base * vy)))
    p_apex_fit = (int(round(x0 + t_apex * vx)), int(round(y0 + t_apex * vy)))

    default_trace.update({
        "is_valid": True,
        "points": pts_arr,
        "confidences": confs_arr,
        "vx": float(vx),
        "vy": float(vy),
        "x0": float(x0),
        "y0": float(y0),
        "p_base": p_base_fit,
        "p_apex": p_apex_fit,
    })
    return default_trace


def reflect_folded_mask(
    mask: np.ndarray,
    p1: Optional[Tuple[int, int]] = None,
    p2: Optional[Tuple[int, int]] = None,
    padding: Optional[int] = None,
    keypoints_31: Optional[Any] = None,
    lm3_symmetry_score: Optional[float] = None,
    p_apex: Optional[Tuple[int, int]] = None,
    p_base: Optional[Tuple[int, int]] = None,
) -> np.ndarray:
    """
    Symmetrically reflects a folded hemi-blade mask across the straight chord line
    passing through p1 and p2 (or across the fitted LeafMachine3 15-point midvein line
    anchored at `lamina_base -> lamina_tip` when `keypoints_31` has >= 3 valid midvein
    keypoints with mean confidence >= 0.50), synthesizing a complete bilateral silhouette.
    Fully vectorized in OpenCV using an affine reflection matrix and bitwise_or.
    Applies morphological closing to seal the midrib seam.

    Parameters
    ----------
    mask : np.ndarray
        Binary uint8 mask of the folded leaf.
    p1 : Optional[Tuple[int, int]]
        Starting endpoint of the fold chord (apex or base). Optional when `keypoints_31`
        is provided or auto-detected.
    p2 : Optional[Tuple[int, int]]
        Ending endpoint of the fold chord (base or apex). Optional when `keypoints_31`
        is provided or auto-detected.
    padding : Optional[int]
        Canvas margin to avoid clipping reflected geometry.
    keypoints_31 : Optional[Dict[str, Tuple[float, float, float]]]
        Optional LeafMachine3 31-keypoint dictionary (`midvein_1..15`, `lamina_base`, `lamina_tip`).
    lm3_symmetry_score : Optional[float]
        Optional LeafMachine3 bilateral symmetry score diagnostic.

    Returns
    -------
    np.ndarray
        Synthesized bilateral silhouette mask (uint8, 0 or 255).
    """
    _ = lm3_symmetry_score
    h, w = mask.shape[:2]
    if h < 5 or w < 5:
        return mask.copy()

    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if bw == 0 or bh == 0:
        return mask.copy()

    if p1 is None and p_apex is not None:
        p1 = p_apex
    if p2 is None and p_base is not None:
        p2 = p_base
    trace = _extract_lm3_midvein_trace(keypoints_31, min_points=3, min_mean_conf=0.50)
    if trace["is_valid"] and trace["p_apex"] is not None and trace["p_base"] is not None:
        p1 = trace["p_apex"]
        p2 = trace["p_base"]
    elif p1 is None or p2 is None:
        _, _, _, pa_m, pb_m = compute_geometric_metrics(mask_u8)
        p1 = pa_m if p1 is None else p1
        p2 = pb_m if p2 is None else p2

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
    keypoints_31: Optional[Any] = None,
    lm3_symmetry_score: Optional[float] = None,
    p_apex: Optional[Tuple[int, int]] = None,
    p_base: Optional[Tuple[int, int]] = None,
) -> FoldDetectionResult:
    """
    Detects whether a leaf is folded in half along the midrib.

    Identifies the primary longitudinal axis connecting the basal petiole insertion
    to the blade apex (using LeafMachine3's 15-point midvein trace `midvein_1..15`
    anchored at `lamina_base -> lamina_tip` when >= 3 valid keypoints have mean
    confidence >= 0.50, or whole-mask principal line otherwise), splits the perimeter
    into Side A and Side B lateral contours, and evaluates linearity, chord straightness,
    and midvein-to-margin lamina width asymmetry.

    Detection Rules:
      - If exactly one lateral boundary is a straight chord (max deviation < 3.5% of
        chord length, R^2 > 0.97 to linear fit) AND W/L <= 0.45:
          Classified as is_folded = True, with the straight chord (or fitted LM3 midvein
          axis) as the true midrib_axis. Symmetrically reflects the curved margin across
          the fold chord.
      - When `keypoints_31` provides a valid midvein trace (>= 3 points, mean conf >= 0.50):
          Computes how close the traced midvein points lie to the outer contour boundary
          and how asymmetric the perpendicular lamina widths are on the left vs. right
          side of the traced midvein axis (incorporating `lm3_symmetry_score` when low),
          fitting the bilateral reflection axis (`p_base -> p_apex`) directly to the
          traced midvein line.
      - If folded irregularly or diagonally (neither side is collinear with petiole-to-apex
        axis), classified as is_irregular_fold = True.
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
        "midvein_trace_used": False,
        "midvein_boundary_ratio": 1.0,
        "lamina_width_asymmetry": 0.0,
        "lm3_symmetry_score": lm3_symmetry_score,
    })

    if mask is None or np.count_nonzero(mask) < 20:
        return default_res

    mask_u8 = mask if mask.dtype == np.uint8 else ((mask > 0).astype(np.uint8) * 255)
    if contour is None:
        contours, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
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

    # 2. Longitudinal Midrib Axis (LM3 15-point midvein trace when available, else cv2.fitLine)
    _ = (p_apex, p_base)
    trace = _extract_lm3_midvein_trace(keypoints_31, min_points=3, min_mean_conf=0.50)
    midvein_folded_trigger = False
    midvein_narrow_side: Optional[str] = None

    if trace["is_valid"] and trace["p_apex"] is not None and trace["p_base"] is not None:
        default_res["midvein_trace_used"] = True
        vx, vy, x0, y0 = trace["vx"], trace["vy"], trace["x0"], trace["y0"]
        projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy
        min_idx = int(np.argmin(projections))
        max_idx = int(np.argmax(projections))
        p_apex = trace["p_apex"]
        p_base = trace["p_base"]

        # Measure proximity of traced midvein points to outer contour boundary
        t_pts = trace["points"]
        d_to_cnt = np.min(
            np.hypot(t_pts[:, None, 0] - pts[None, :, 0], t_pts[:, None, 1] - pts[None, :, 1]),
            axis=1,
        )
        midvein_boundary_ratio = float(np.median(d_to_cnt)) / max(float(box_wid), 1.0)
        default_res["midvein_boundary_ratio"] = round(midvein_boundary_ratio, 4)

        # Measure left vs. right lamina widths perpendicular to the traced midvein axis
        mh, mw = mask_u8.shape[:2]
        nx, ny = -vy, vx
        max_probe = max(12, int(round(max(box_wid, 20.0) * 1.25)))
        t_steps = np.arange(1, max_probe + 1, dtype=np.float64)
        w_pos_list: List[float] = []
        w_neg_list: List[float] = []
        for mp_x, mp_y in t_pts:
            px_p = np.round(mp_x + t_steps * nx).astype(np.int32)
            py_p = np.round(mp_y + t_steps * ny).astype(np.int32)
            vp = (px_p >= 0) & (px_p < mw) & (py_p >= 0) & (py_p < mh)
            fg_p = np.zeros(len(t_steps), dtype=bool)
            fg_p[vp] = mask_u8[py_p[vp], px_p[vp]] > 0
            nz_p = np.where(fg_p)[0]
            w_pos_list.append(float(nz_p[-1] + 1) if len(nz_p) > 0 else 0.0)

            px_n = np.round(mp_x - t_steps * nx).astype(np.int32)
            py_n = np.round(mp_y - t_steps * ny).astype(np.int32)
            vn = (px_n >= 0) & (px_n < mw) & (py_n >= 0) & (py_n < mh)
            fg_n = np.zeros(len(t_steps), dtype=bool)
            fg_n[vn] = mask_u8[py_n[vn], px_n[vn]] > 0
            nz_n = np.where(fg_n)[0]
            w_neg_list.append(float(nz_n[-1] + 1) if len(nz_n) > 0 else 0.0)

        mean_w_pos = float(np.mean(w_pos_list)) if w_pos_list else 0.0
        mean_w_neg = float(np.mean(w_neg_list)) if w_neg_list else 0.0
        lamina_asym = abs(mean_w_pos - mean_w_neg) / max(mean_w_pos + mean_w_neg, 1.0)
        default_res["lamina_width_asymmetry"] = round(lamina_asym, 4)
        midvein_narrow_side = "Side A" if mean_w_pos <= mean_w_neg else "Side B"

        low_sym_bonus = lm3_symmetry_score is not None and float(lm3_symmetry_score) < 0.60
        eff_wl_cap = max(float(max_wl_ratio), 0.55 if low_sym_bonus else 0.48)
        if wl_ratio <= eff_wl_cap and (
            midvein_boundary_ratio <= 0.20
            or lamina_asym >= 0.60
            or (low_sym_bonus and (midvein_boundary_ratio <= 0.28 or lamina_asym >= 0.48))
        ):
            midvein_folded_trigger = True
    else:
        pts_f64 = pts.astype(np.float64)
        mean_pt = np.mean(pts_f64, axis=0)
        centered = pts_f64 - mean_pt
        cov = np.cov(centered, rowvar=False)
        if cov.shape == (2, 2) and np.all(np.isfinite(cov)):
            eigvals, eigvecs = np.linalg.eigh(cov)
            principal = eigvecs[:, int(np.argmax(eigvals))]
            vx, vy = float(principal[0]), float(principal[1])
            x0, y0 = float(mean_pt[0]), float(mean_pt[1])
        else:
            line_params = cv2.fitLine(pts, cv2.DIST_L2, 0, 0.01, 0.01)
            vx, vy, x0, y0 = [float(v[0]) for v in line_params]
        projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy

        min_idx = int(np.argmin(projections))
        max_idx = int(np.argmax(projections))
        p_apex = (int(pts[min_idx, 0]), int(pts[min_idx, 1]))
        p_base = (int(pts[max_idx, 0]), int(pts[max_idx, 1]))
        if trace["lamina_base"] is not None:
            lb_x, lb_y = trace["lamina_base"]
            if math.hypot(p_apex[0] - lb_x, p_apex[1] - lb_y) < math.hypot(p_base[0] - lb_x, p_base[1] - lb_y):
                p_apex, p_base = p_base, p_apex

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
    if (straight_a ^ straight_b and wl_ratio <= max_wl_ratio) or midvein_folded_trigger:
        default_res["is_folded"] = True
        if straight_a ^ straight_b:
            default_res["straight_side"] = "Side A" if straight_a else "Side B"
        else:
            default_res["straight_side"] = midvein_narrow_side or ("Side A" if dev_a <= dev_b else "Side B")
        default_res["midrib_axis"] = (p_apex, p_base)
        default_res["synthesized_mask"] = reflect_folded_mask(
            mask_u8,
            p_apex,
            p_base,
            keypoints_31=keypoints_31,
            lm3_symmetry_score=lm3_symmetry_score,
        )
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


def bridge_dissected_rachis_components(
    mask: np.ndarray,
    petiole_mask: Optional[np.ndarray] = None,
    max_gap_fraction: float = 0.22,
    crop_bgr: Optional[np.ndarray] = None,
    specimen_paper_cleaned_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Reconstructs micro-severed rachis internodes in multi-component dissected leaf masks.

    In deeply pinnatisect/bipinnatisect leaves (e.g. Packera glabella), narrow or winged
    rachis segments between the terminal lobe and lateral pinna pairs are occasionally
    fragmented into 2+ foreground islands by LeafMachine3 (YOLO26x) thresholding (or assigned to
    LeafMachine3 class 1 `petiole` while the pinna lobes are assigned to class 0 `lamina`).

    Invariants:
      - Fast no-op (`0` pixels modified) whenever `mask` already has <= 1 foreground
        connected component (`num_labels <= 2`).
      - Label-blind: relies purely on component area, longitudinal axis collinearity,
        border-clipping exclusion, and anatomical rachis-spine bridging so outer
        inter-lobe sinuses remain open concavity defects.
      - When `specimen_paper_cleaned_mask` (from LeafMachine3 UNet++ v9 `LLR_Hyst_KeepBlobs`
        whole-specimen segmenter) is provided, constrains rachis-bridging corridors to
        non-paper plant pixels: boosts bridging affinity where `specimen_paper_cleaned_mask > 0`
        covers >= 50% of the corridor between adjacent pinnatifid lobes, while rejecting
        corridors that cross > 70% pure paper background (`specimen_paper_cleaned_mask == 0`).
      - When `crop_bgr` is provided, verifies that any synthetic bridge segment > 10 px
        is supported by non-paper botanical tissue (or `petiole_mask`) rather than drawing
        artificial straight lines across bright mounting paper between separate leaves.
    """
    if mask is None or mask.size == 0:
        return mask

    mask_u8 = (mask > 0).astype(np.uint8) * 255
    num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(
        (mask_u8 > 0).astype(np.uint8), connectivity=8
    )
    if num_labels <= 2:
        return mask.copy()

    h, w = mask_u8.shape[:2]
    spec_clean_bool: Optional[np.ndarray] = None
    if specimen_paper_cleaned_mask is not None and specimen_paper_cleaned_mask.size > 0:
        if specimen_paper_cleaned_mask.shape[:2] == (h, w):
            spec_clean_bool = specimen_paper_cleaned_mask > 0
        else:
            spec_rs = cv2.resize(
                (specimen_paper_cleaned_mask > 0).astype(np.uint8) * 255,
                (w, h),
                interpolation=cv2.INTER_NEAREST,
            )
            spec_clean_bool = spec_rs > 0

    areas = stats[1:, cv2.CC_STAT_AREA]
    main_id = 1 + int(np.argmax(areas))
    main_area = int(stats[main_id, cv2.CC_STAT_AREA])
    if main_area < 25:
        return mask.copy()

    min_pinna_area = max(14, int(0.015 * main_area))
    cand_ids: List[int] = [main_id]
    for cid in range(1, num_labels):
        if cid == main_id:
            continue
        c_area = int(stats[cid, cv2.CC_STAT_AREA])
        if c_area < min_pinna_area:
            continue
        cx0 = int(stats[cid, cv2.CC_STAT_LEFT])
        cy0 = int(stats[cid, cv2.CC_STAT_TOP])
        cw0 = int(stats[cid, cv2.CC_STAT_WIDTH])
        ch0 = int(stats[cid, cv2.CC_STAT_HEIGHT])
        touches_border = (cx0 <= 1) or (cy0 <= 1) or (cx0 + cw0 >= w - 1) or (cy0 + ch0 >= h - 1)
        if touches_border:
            comp_bool = labels == cid
            border_px = (
                int(np.count_nonzero(comp_bool[:, :2]))
                + int(np.count_nonzero(comp_bool[:, max(0, w - 2) :]))
                + int(np.count_nonzero(comp_bool[:2, :]))
                + int(np.count_nonzero(comp_bool[max(0, h - 2) :, :]))
            )
            has_pet_overlap = False
            if petiole_mask is not None and petiole_mask.shape[:2] == mask_u8.shape[:2]:
                has_pet_overlap = int(np.count_nonzero(comp_bool & (petiole_mask > 0))) >= max(4, int(0.15 * c_area))
            if border_px >= max(6, int(0.08 * max(cw0, ch0))) and not has_pet_overlap:
                continue
        cand_ids.append(cid)

    if len(cand_ids) <= 1:
        # Strip tiny unbridged noise specks (< 1.5% of main_area) so the primary component is clean
        return np.where(labels == main_id, 255, 0).astype(np.uint8)

    union_mask = np.isin(labels, cand_ids)
    ys_u, xs_u = np.where(union_mask)
    if len(xs_u) < 20:
        return np.where(labels == main_id, 255, 0).astype(np.uint8)

    pts_u = np.column_stack([xs_u, ys_u]).astype(np.float32)
    if petiole_mask is not None and petiole_mask.shape[:2] == mask_u8.shape[:2]:
        ys_p, xs_p = np.where(petiole_mask > 0)
        if len(xs_p) >= 8:
            pts_fit = np.vstack([pts_u, np.column_stack([xs_p, ys_p]).astype(np.float32)])
        else:
            pts_fit = pts_u
    else:
        pts_fit = pts_u

    line_params = cv2.fitLine(pts_fit, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]

    proj_u = (pts_u[:, 0] - x0) * vx + (pts_u[:, 1] - y0) * vy
    norm_u = (pts_u[:, 0] - x0) * (-vy) + (pts_u[:, 1] - y0) * vx
    s_min, s_max = float(np.min(proj_u)), float(np.max(proj_u))
    major_span = max(1.0, s_max - s_min)
    minor_span = max(1.0, float(np.max(norm_u) - np.min(norm_u)))
    max_gap_px = max(16.0, float(max_gap_fraction) * max(major_span, minor_span))
    max_axis_offset = max(20.0, 0.45 * max(major_span, minor_span))

    # Optional grayscale tissue map from crop_bgr to prevent bridging across blank mounting paper
    paper_thresh: Optional[float] = None
    gray_crop: Optional[np.ndarray] = None
    if crop_bgr is not None and crop_bgr.shape[:2] == mask_u8.shape[:2] and crop_bgr.size > 0:
        gray_crop = (
            cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2GRAY)
            if crop_bgr.ndim == 3
            else crop_bgr
        )
        bg_px = gray_crop[mask_u8 == 0]
        fg_px = gray_crop[mask_u8 > 0]
        if bg_px.size >= 32 and fg_px.size >= 32:
            bg_med = float(np.percentile(bg_px, 65))
            fg_p75 = float(np.percentile(fg_px, 75))
            if bg_med > fg_p75 + 10.0:
                paper_thresh = 0.55 * bg_med + 0.45 * fg_p75

    def _bridge_path_has_tissue(pt_a: Tuple[int, int], pt_b: Tuple[int, int]) -> bool:
        seg_len = math.hypot(float(pt_b[0] - pt_a[0]), float(pt_b[1] - pt_a[1]))
        n_samples = max(8, int(round(max(seg_len, 4.0))))
        xs_s = np.clip(np.round(np.linspace(pt_a[0], pt_b[0], n_samples)).astype(int), 0, w - 1)
        ys_s = np.clip(np.round(np.linspace(pt_a[1], pt_b[1], n_samples)).astype(int), 0, h - 1)

        if spec_clean_bool is not None:
            # Evaluate corridor coverage on UNet++ v9 whole-specimen binary mask (with 1px tolerance)
            corridor_hits = np.zeros(n_samples, dtype=bool)
            for dy_c in (-1, 0, 1):
                for dx_c in (-1, 0, 1):
                    yc = np.clip(ys_s + dy_c, 0, h - 1)
                    xc = np.clip(xs_s + dx_c, 0, w - 1)
                    corridor_hits |= spec_clean_bool[yc, xc]
            plant_frac = float(np.mean(corridor_hits))
            if plant_frac < 0.60:
                return False
            return True

        if seg_len <= 10.0:
            return True
        if gray_crop is None or paper_thresh is None:
            return True
        # Check a neighborhood along the bridge segment
        has_pet = False
        if petiole_mask is not None and petiole_mask.shape[:2] == mask_u8.shape[:2]:
            pet_hits = int(np.count_nonzero(petiole_mask[ys_s, xs_s] > 0))
            if float(pet_hits) / float(n_samples) >= 0.25:
                has_pet = True
        if has_pet:
            return True
        if seg_len > 0.14 * max(major_span, minor_span):
            return False
        vals = gray_crop[ys_s, xs_s].astype(np.float32)
        non_mask = mask_u8[ys_s, xs_s] == 0
        if int(np.count_nonzero(non_mask)) <= 3:
            return True
        gap_vals = vals[non_mask]
        tissue_frac = float(np.mean(gap_vals <= paper_thresh))
        return tissue_frac >= 0.45

    # Filter secondary components whose centroids lie within the leaf's lateral envelope
    valid_sec_ids: List[int] = []
    for cid in cand_ids[1:]:
        cx_c, cy_c = float(centroids[cid][0]), float(centroids[cid][1])
        perp_c = abs((cx_c - x0) * (-vy) + (cy_c - y0) * vx)
        if perp_c <= max_axis_offset:
            valid_sec_ids.append(cid)

    if not valid_sec_ids:
        return np.where(labels == main_id, 255, 0).astype(np.uint8)

    bridged = np.where(np.isin(labels, [main_id] + valid_sec_ids), 255, 0).astype(np.uint8)
    bridge_w = max(3, min(10, int(round(0.03 * minor_span))))

    # If native LeafMachine3 (YOLO26x) petiole_mask overlaps the intra-lamina span between lobes, include its rachis pixels
    if petiole_mask is not None and petiole_mask.shape[:2] == mask_u8.shape[:2]:
        ys_pm, xs_pm = np.where(petiole_mask > 0)
        if len(xs_pm) > 0:
            s_pm = (xs_pm - x0) * vx + (ys_pm - y0) * vy
            n_pm = np.abs((xs_pm - x0) * (-vy) + (ys_pm - y0) * vx)
            in_rachis = (s_pm >= s_min - 4.0) & (s_pm <= s_max + 4.0) & (n_pm <= 0.25 * minor_span + 6.0)
            if np.any(in_rachis):
                bridged[ys_pm[in_rachis], xs_pm[in_rachis]] = 255

    # Include verified UNet++ v9 specimen_paper_cleaned_mask rachis tissue along the central midvein corridor
    if spec_clean_bool is not None:
        ys_sp, xs_sp = np.where(spec_clean_bool)
        if len(xs_sp) > 0:
            s_sp = (xs_sp - x0) * vx + (ys_sp - y0) * vy
            n_sp = np.abs((xs_sp - x0) * (-vy) + (ys_sp - y0) * vx)
            in_spec_rachis = (s_sp >= s_min - 2.0) & (s_sp <= s_max + 2.0) & (n_sp <= max(6.0, 0.08 * minor_span))
            if np.any(in_spec_rachis):
                bridged[ys_sp[in_spec_rachis], xs_sp[in_spec_rachis]] = 255

    # Iteratively bridge any remaining disconnected valid components along the rachis spine
    rejected_cids: Set[int] = set()
    for _ in range(len(valid_sec_ids)):
        n_cur, lbl_cur, st_cur, _ = cv2.connectedComponentsWithStats((bridged > 0).astype(np.uint8), connectivity=8)
        if n_cur <= 2:
            break
        cur_main = 1 + int(np.argmax(st_cur[1:, cv2.CC_STAT_AREA]))
        main_cnts, _ = cv2.findContours(
            (lbl_cur == cur_main).astype(np.uint8) * 255, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )
        if not main_cnts:
            break
        pts_main = np.vstack([c.reshape(-1, 2) for c in main_cnts]).astype(np.float64)
        s_main = (pts_main[:, 0] - x0) * vx + (pts_main[:, 1] - y0) * vy
        n_main = np.abs((pts_main[:, 0] - x0) * (-vy) + (pts_main[:, 1] - y0) * vx)

        best_pair = None
        best_dist = float("inf")
        best_cid = -1
        for cid in range(1, n_cur):
            if cid == cur_main or cid in rejected_cids:
                continue
            if int(st_cur[cid, cv2.CC_STAT_AREA]) < min_pinna_area:
                continue
            c_cnts, _ = cv2.findContours(
                (lbl_cur == cid).astype(np.uint8) * 255, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
            )
            if not c_cnts:
                continue
            pts_c = np.vstack([c.reshape(-1, 2) for c in c_cnts]).astype(np.float64)
            # Pairwise distances between pts_c and pts_main (subsampled if huge)
            step_m = max(1, len(pts_main) // 200)
            step_c = max(1, len(pts_c) // 120)
            pm_sub = pts_main[::step_m]
            pc_sub = pts_c[::step_c]
            diff = pc_sub[:, None, :] - pm_sub[None, :, :]
            dmat = np.hypot(diff[:, :, 0], diff[:, :, 1])
            idx_flat = int(np.argmin(dmat))
            ic, im = divmod(idx_flat, len(pm_sub))
            d_min = float(dmat[ic, im])
            if d_min < best_dist:
                best_dist = d_min
                best_cid = cid
                best_pair = (pts_c, pc_sub[ic], pm_sub[im], d_min)

        if best_pair is None or best_dist > max_gap_px:
            break

        pts_c, pt_c_direct, pt_m_direct, d_min = best_pair
        # Determine whether to bridge directly or via the anatomical central rachis line
        n_c_direct = abs((pt_c_direct[0] - x0) * (-vy) + (pt_c_direct[1] - y0) * vx)
        if d_min <= 14.0 or n_c_direct <= 0.16 * minor_span:
            p1 = (int(round(pt_c_direct[0])), int(round(pt_c_direct[1])))
            p2 = (int(round(pt_m_direct[0])), int(round(pt_m_direct[1])))
            if _bridge_path_has_tissue(p1, p2):
                cv2.line(bridged, p1, p2, 255, thickness=bridge_w)
            else:
                bridged[lbl_cur == best_cid] = 0
                rejected_cids.add(best_cid)
        else:
            # Find inner rachis-facing base of the lateral pinna component
            s_c = (pts_c[:, 0] - x0) * vx + (pts_c[:, 1] - y0) * vy
            n_c = np.abs((pts_c[:, 0] - x0) * (-vy) + (pts_c[:, 1] - y0) * vx)
            d_to_m = np.hypot(pts_c[:, 0] - pt_m_direct[0], pts_c[:, 1] - pt_m_direct[1])
            in_idx = int(np.argmin(n_c + 0.20 * d_to_m))
            pt_in = pts_c[in_idx]
            s_in = float(s_c[in_idx])
            pt_rachis = np.array([x0 + s_in * vx, y0 + s_in * vy], dtype=np.float64)
            # Nearest rachis-aligned point on main component
            rachis_score_m = np.abs(s_main - s_in) + 1.5 * n_main
            m_idx = int(np.argmin(rachis_score_m))
            pt_m_rachis = pts_main[m_idx]

            p_in = (int(round(pt_in[0])), int(round(pt_in[1])))
            p_rac = (
                int(np.clip(round(pt_rachis[0]), 0, w - 1)),
                int(np.clip(round(pt_rachis[1]), 0, h - 1)),
            )
            p_mr = (int(round(pt_m_rachis[0])), int(round(pt_m_rachis[1])))
            if _bridge_path_has_tissue(p_in, p_rac) and _bridge_path_has_tissue(p_rac, p_mr):
                cv2.line(bridged, p_in, p_rac, 255, thickness=bridge_w)
                cv2.line(bridged, p_rac, p_mr, 255, thickness=bridge_w)
            else:
                bridged[lbl_cur == best_cid] = 0
                rejected_cids.add(best_cid)

    # Keep largest connected component after bridging
    n_fin, lbl_fin, st_fin, _ = cv2.connectedComponentsWithStats((bridged > 0).astype(np.uint8), connectivity=8)
    if n_fin > 2:
        best_fin = 1 + int(np.argmax(st_fin[1:, cv2.CC_STAT_AREA]))
        bridged = np.where(lbl_fin == best_fin, 255, 0).astype(np.uint8)

    return bridged


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
        bin_in = (contour > 0).astype(np.uint8) * 255
        cnts_raw, _ = cv2.findContours(bin_in, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not cnts_raw:
            return False
        if len(cnts_raw) > 1:
            bin_in = bridge_dissected_rachis_components(bin_in)
            cnts, _ = cv2.findContours(bin_in, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if not cnts:
                return False
        else:
            cnts = cnts_raw
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

    # Reject equidimensional multi-leaf overlapping blobs masquerading as dissected leaves on high-DPI scans,
    # without penalizing physically wide high-DPI dissected blades (where minor_w routinely exceeds 360 px).
    is_high_dpi = (
        (contour.ndim == 2 and max(contour.shape[:2]) > 500)
        or cnt_f[:, 0].max() > 500
        or cnt_f[:, 1].max() > 500
        or (solidity >= 0.48 and major_l > 240.0)
        or minor_w > 360.0
    )
    if is_high_dpi and major_l > 220.0 and solidity >= 0.55 and (minor_w / max(major_l, 1.0)) > 0.88:
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
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None,
    specimen_paper_cleaned_mask: Optional[np.ndarray] = None,
    lm3_symmetry_score: Optional[float] = None,
) -> LeafRoutingResult:
    """
    Comprehensive 100% label-blind geometric quality gatekeeper and routing engine
    for LeafMachine3 (`yolo26x_seg_1024` 3-class masks & `yolo26x_pose_640` 31-keypoints):
      0. Pre-processes LM3 3-class segmentation masks when provided:
         - Keeps the largest 8-connected component of `mask > 0` as `silhouette` without
           subtracting `petiole_mask`, and subtracts 8-connected `hole_mask > 0` components
           with area >= 64 px (`lamina = silhouette & ~hole`).
      1. Evaluates midrib pose (disambiguating `p_base` vs. `p_apex` using `lamina_base`
         from `keypoints_31` when confidence >= 0.50) and initial solidity/UCS metrics.
      2. Straight-chord & 15-point midvein trace fold detection:
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

    _ = specimen_paper_cleaned_mask
    petiole_subtracted = False
    petiole_mask_area_px = 0
    petiole_tip_pt: Optional[Tuple[int, int]] = None
    if (
        petiole_mask is not None
        and isinstance(petiole_mask, np.ndarray)
        and petiole_mask.shape[:2] == mask.shape[:2]
        and np.count_nonzero(petiole_mask) > 0
    ):
        petiole_mask_area_px = int(np.count_nonzero(petiole_mask > 0))
        py_idx, px_idx = np.where(petiole_mask > 0)
        if len(py_idx) > 0:
            petiole_tip_pt = (int(round(float(np.mean(px_idx)))), int(round(float(np.mean(py_idx)))))

    fg_u8 = (mask > 0).astype(np.uint8)
    n_cc, cc_lbls, cc_stats, _ = cv2.connectedComponentsWithStats(fg_u8, connectivity=8)
    if n_cc > 1:
        largest_idx = 1 + int(np.argmax(cc_stats[1:, cv2.CC_STAT_AREA]))
        largest_comp = (cc_lbls == largest_idx).astype(np.uint8)
        ext_cnts, _ = cv2.findContours(largest_comp, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        ext_polys = [c.reshape(-1, 2).astype(np.int32) for c in ext_cnts if len(c.reshape(-1, 2)) >= 3]
        work_mask = np.zeros(mask.shape[:2], dtype=np.uint8)
        if ext_polys:
            cv2.fillPoly(work_mask, ext_polys, 255)
        else:
            work_mask = largest_comp * 255
    else:
        work_mask = (mask > 0).astype(np.uint8) * 255

    hole_filled = False
    hole_mask_area_px = 0
    filled_hole_area_px = 0
    if (
        hole_mask is not None
        and isinstance(hole_mask, np.ndarray)
        and hole_mask.shape[:2] == work_mask.shape[:2]
        and np.count_nonzero(hole_mask) > 0
    ):
        min_area_px = 64
        hole_u8 = (hole_mask > 0).astype(np.uint8)
        hole_mask_area_px = int(np.count_nonzero(hole_u8))
        n_hc, hc_lbls, hc_stats, _ = cv2.connectedComponentsWithStats(hole_u8, connectivity=8)
        qualifying_holes = np.zeros(work_mask.shape[:2], dtype=np.uint8)
        hole_polys: List[np.ndarray] = []
        for hid in range(1, n_hc):
            if int(hc_stats[hid, cv2.CC_STAT_AREA]) >= min_area_px:
                h_comp = (hc_lbls == hid).astype(np.uint8)
                h_cnts, _ = cv2.findContours(h_comp, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                for hc in h_cnts:
                    pts_h = hc.reshape(-1, 2).astype(np.int32)
                    if len(pts_h) >= 3:
                        hole_polys.append(pts_h)
        if hole_polys:
            cv2.fillPoly(qualifying_holes, hole_polys, 255)
            work_mask[qualifying_holes > 0] = 0

    if int(np.count_nonzero(work_mask)) > 0:
        n_out, out_lbls, out_stats, _ = cv2.connectedComponentsWithStats(
            (work_mask > 0).astype(np.uint8), connectivity=8
        )
        if n_out > 2:
            best_out = 1 + int(np.argmax(out_stats[1:, cv2.CC_STAT_AREA]))
            work_mask = np.where(out_lbls == best_out, 255, 0).astype(np.uint8)
        mask = (work_mask > 0).astype(np.uint8) * 255

    trace_kp = _extract_lm3_midvein_trace(keypoints_31, min_points=3, min_mean_conf=0.50)

    def _apply_kp31_base_apex_disambiguation(
        pa_in: Tuple[int, int], pb_in: Tuple[int, int]
    ) -> Tuple[Tuple[int, int], Tuple[int, int], Optional[float]]:
        if trace_kp["lamina_base"] is not None:
            lbx, lby = trace_kp["lamina_base"]
            d_base = math.hypot(pb_in[0] - lbx, pb_in[1] - lby)
            d_apex = math.hypot(pa_in[0] - lbx, pa_in[1] - lby)
            if d_apex + 2.0 < d_base:
                pa_in, pb_in = pb_in, pa_in
            ang_out = math.degrees(math.atan2(float(pb_in[1] - pa_in[1]), float(pb_in[0] - pa_in[0])))
            return pa_in, pb_in, ang_out
        return pa_in, pb_in, None

    # Compute initial geometric metrics and disambiguated (p_apex, p_base)
    ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(
        mask, rosette_centroid=rosette_centroid
    )
    p_apex, p_base, ang_override = _apply_kp31_base_apex_disambiguation(p_apex, p_base)
    if ang_override is not None:
        angle_deg = ang_override

    _ = already_trimmed
    if trim_meta is None:
        trim_meta = {
            "is_trimmed": petiole_subtracted,
            "reason": "LM3_PETIOLE_MASK_SUBTRACTED" if petiole_subtracted else "PRE_TRIMMED",
        }
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

    # 1. Straight-chord & LM3 15-point midvein trace fold detection
    if fold_detection_enabled:
        fold_info = detect_folded_leaf(
            mask,
            main_cnt,
            max_dev_ratio=max_chord_deviation_ratio,
            max_wl_ratio=max_folded_aspect_ratio,
            min_chord_length_px=min_chord_length_px,
            keypoints_31=keypoints_31,
            lm3_symmetry_score=lm3_symmetry_score,
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
        "petiole_mask_subtracted": petiole_subtracted,
        "petiole_mask_area_px": petiole_mask_area_px,
        "petiole_tip_pt": petiole_tip_pt,
        "hole_mask_filled": hole_filled,
        "hole_mask_area_px": hole_mask_area_px,
        "filled_hole_area_px": filled_hole_area_px,
        "keypoints_31_used": bool(keypoints_31),
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
        # Check if finely dissected / bipinnatisect or pinnatisect leaf (0.18 <= solidity < min_solidity) triggered crease heuristic
        if dissection_enabled and (0.18 <= solidity < min_solidity) and is_botanical_dissection(
            main_cnt,
            solidity=solidity,
            p_apex=p_apex,
            p_base=p_base,
            min_solidity=0.18,
            max_solidity=min_solidity,
            defect_depth_ratio=sinus_defect_min_depth_ratio,
            min_defect_count=max(min_bilateral_sinus_count, 4) if solidity < 0.48 else min_bilateral_sinus_count,
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
        if (is_bipinnatisect_window or not is_high_res_mask or (ucs >= 0.35 and hull_ucs >= 0.72))
        else ucs
    )
    meta["hull_ucs"] = round(float(hull_ucs), 4)
    allow_dissection = True
    if is_bipinnatisect_window:
        if eff_diss_ucs < 0.55:
            allow_dissection = False
    elif is_high_res_mask and (eff_diss_ucs < 0.75 or ucs < 0.35):
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
        min_req_ucs = 0.55 if is_bipinnatisect_window else (0.75 if is_high_res_mask else 0.70)
        if dissected and (eff_diss_ucs >= min_req_ucs or (not is_high_res_mask and eff_diss_ucs >= 0.55)):
            coords = resample_normalized_contour(mask, num_contour_points, petiole_attachment_pt=p_base) if extract_contours else None
            meta["reflection_applied"] = False
            dissected_tier = "Tier 1 (Bipinnatisect)" if is_bipinnatisect_window else "Tier 1 (Dissected)"
            return LeafRoutingResult(
                tier=dissected_tier,
                assigned_tier="tier1",
                solidity=solidity,
                ucs=max(eff_diss_ucs, 0.80 if not is_bipinnatisect_window else eff_diss_ucs),
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
