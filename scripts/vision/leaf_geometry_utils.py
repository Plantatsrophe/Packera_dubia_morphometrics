#!/usr/bin/env python3
"""
===============================================================================
Module: leaf_geometry_utils.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Lean LeafMachine3 (LM3) mask persistence, native polygon morphology readout,
    and canonical tip-up 120-vertex contour coordinate exporter for Step 02.
    Delegates contour extraction, polygon morphology, orientation, and closed
    resampling directly to LeafMachine3 without mutating LM3 masks.
===============================================================================
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.lm3_compat import (
    compute_orientation,
    extract_contour,
    length_width_from_box,
    momocs_signed_area,
    polygon_morphology,
    resample_closed,
)
from scripts.vision.organ_crown_metrics import (
    _estimate_single_plant_crown,
    _extract_xy_point,
    _parse_basal_leaf_pose,
    approximate_petiole_and_scape_metrics,
    classify_leaf_zone_anatomical,
    cluster_plant_individuals,
    estimate_rosette_crown_point,
)


def open_hairline_sinus_bridges_and_measure_voids(
    mask: Optional[np.ndarray],
    hole_mask: Optional[np.ndarray] = None,
    max_bridge_px: float = 3.5,
) -> Tuple[np.ndarray, float, int]:
    """
    Inspects interior `0`-pixel voids inside an unfilled lamina `mask`:
      1. Opens hairline `<= max_bridge_px` lobe-tip contact bridges separating non-hole
         dissected sinus pockets from the exterior background (without disconnecting the
         lamina) so single-contour `extract_contour` enters the true sinus.
      2. Measures `enclosed_void_frac` across any remaining broadly bridged (`> max_bridge_px`)
         interior non-hole voids so Step 02/03 can withhold distorted single-loop EFA harmonics
         while retaining the hole-preserved mask for topological ECT and solidity.
    Returns `(contour_ready_mask_u8, enclosed_void_frac, n_opened_bridges)`.
    """
    if mask is None or np.asarray(mask).size == 0:
        return np.zeros((0, 0), dtype=np.uint8), 0.0, 0
    work_u8 = (np.asarray(mask) > 0).astype(np.uint8) * 255
    fg_area = int(np.count_nonzero(work_u8))
    if fg_area < 32:
        return work_u8, 0.0, 0

    h, w = work_u8.shape[:2]
    padded_inv = np.ones((h + 2, w + 2), dtype=np.uint8)
    padded_inv[1 : h + 1, 1 : w + 1] = (work_u8 == 0).astype(np.uint8)
    num_bg, bg_labels, bg_stats, _ = cv2.connectedComponentsWithStats(padded_inv, connectivity=8)
    if num_bg <= 2:
        return work_u8, 0.0, 0

    ext_lbl = int(bg_labels[0, 0])
    ext_bg = (bg_labels[1 : h + 1, 1 : w + 1] == ext_lbl)
    min_void_px = max(16, int(round(0.003 * float(fg_area))))
    hole_bool = (
        (np.asarray(hole_mask) > 0)
        if (hole_mask is not None and isinstance(hole_mask, np.ndarray) and hole_mask.shape[:2] == (h, w))
        else None
    )

    cand_void_lbls = [
        lbl
        for lbl in range(1, num_bg)
        if lbl != ext_lbl and int(bg_stats[lbl, cv2.CC_STAT_AREA]) >= min_void_px
    ]
    if not cand_void_lbls:
        return work_u8, 0.0, 0

    dt_from_ext = cv2.distanceTransform((~ext_bg).astype(np.uint8), cv2.DIST_L2, 5)
    n_fg_init, _, fg_stats_init, _ = cv2.connectedComponentsWithStats((work_u8 > 0).astype(np.uint8), connectivity=8)
    sig_fg_init = sum(1 for c_i in range(1, n_fg_init) if int(fg_stats_init[c_i, cv2.CC_STAT_AREA]) >= 32)

    opened_count = 0
    unopened_void_px = 0
    inner_bg_lbls = bg_labels[1 : h + 1, 1 : w + 1]

    for v_lbl in cand_void_lbls:
        v_mask = inner_bg_lbls == v_lbl
        v_area = int(np.count_nonzero(v_mask))
        if v_area < min_void_px:
            continue
        if hole_bool is not None:
            overlap_h = float(np.count_nonzero(v_mask & hole_bool)) / float(max(v_area, 1))
            if overlap_h > 0.30:
                continue

        v_dt_vals = dt_from_ext[v_mask]
        if v_dt_vals.size == 0:
            continue
        min_dt = float(np.min(v_dt_vals))
        bridge_px = max(0.0, min_dt - 1.0)
        if bridge_px <= float(max_bridge_px):
            v_ys, v_xs = np.where(v_mask)
            best_v_i = int(np.argmin(v_dt_vals))
            vy0, vx0 = int(v_ys[best_v_i]), int(v_xs[best_v_i])
            rad_w = int(math.ceil(max_bridge_px)) + 4
            wy0, wy1 = max(0, vy0 - rad_w), min(h, vy0 + rad_w + 1)
            wx0, wx1 = max(0, vx0 - rad_w), min(w, vx0 + rad_w + 1)
            ext_sub = ext_bg[wy0:wy1, wx0:wx1]
            if np.any(ext_sub):
                ey_s, ex_s = np.where(ext_sub)
                d2_ext = (ey_s + wy0 - vy0) ** 2 + (ex_s + wx0 - vx0) ** 2
                best_e_i = int(np.argmin(d2_ext))
                ey0, ex0 = int(ey_s[best_e_i] + wy0), int(ex_s[best_e_i] + wx0)
                trial_u8 = work_u8.copy()
                cv2.line(trial_u8, (vx0, vy0), (ex0, ey0), 0, thickness=2)
                n_fg_tr, _, fg_stats_tr, _ = cv2.connectedComponentsWithStats((trial_u8 > 0).astype(np.uint8), connectivity=8)
                sig_fg_tr = sum(1 for c_i in range(1, n_fg_tr) if int(fg_stats_tr[c_i, cv2.CC_STAT_AREA]) >= 32)
                if sig_fg_tr <= max(1, sig_fg_init):
                    work_u8 = trial_u8
                    ext_bg = ext_bg | v_mask | (work_u8 == 0)
                    dt_from_ext = cv2.distanceTransform((~ext_bg).astype(np.uint8), cv2.DIST_L2, 5)
                    opened_count += 1
                    continue

        unopened_void_px += v_area

    enclosed_void_frac = float(unopened_void_px) / float(max(fg_area, 1))
    return work_u8, round(enclosed_void_frac, 4), opened_count


def compute_lm3_shape_metrics(
    mask: Optional[np.ndarray],
    angle_cw: Optional[float] = None,
    bbox: Optional[Tuple[int, int, int, int]] = None,
) -> Dict[str, float]:
    """Computes polygon morphology scalars (`area_px`, `solidity`, `circularity`, `extent`, anatomical `aspect_ratio`) via LeafMachine3 `polygon_morphology`, preserving interior sinus voids (`cv2.RETR_CCOMP`) and multi-component pinna islands."""
    empty = {
        "area_px": 0.0,
        "solidity": 0.0,
        "circularity": 0.0,
        "extent": 0.0,
        "aspect_ratio": 0.0,
    }
    if mask is None or mask.size == 0:
        return empty
    mask_u8 = (mask > 0).astype(np.uint8) * 255
    area_px = float(np.count_nonzero(mask_u8))
    if area_px <= 0.0:
        return empty
    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    if bw <= 0 or bh <= 0:
        return empty
    sub = mask_u8[by : by + bh, bx : bx + bw]
    cnts_ccomp, hier_ccomp = cv2.findContours(sub, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    outer_cnts = [
        cnts_ccomp[i]
        for i in range(len(cnts_ccomp))
        if hier_ccomp is not None and int(hier_ccomp[0][i][3]) == -1 and len(cnts_ccomp[i]) >= 3
    ]
    hole_cnts = [
        cnts_ccomp[i]
        for i in range(len(cnts_ccomp))
        if hier_ccomp is not None and int(hier_ccomp[0][i][3]) != -1 and abs(float(cv2.contourArea(cnts_ccomp[i]))) >= 12.0
    ]
    sig_cnts = [c for c in outer_cnts if abs(float(cv2.contourArea(c))) >= 16.0]
    cnt_xy = extract_contour(sub)
    if cnt_xy is None or len(cnt_xy) < 3:
        return {
            "area_px": area_px,
            "solidity": 1.0,
            "circularity": 0.0,
            "extent": min(1.0, area_px / float(max(bw * bh, 1))),
            "aspect_ratio": float(max(bw, bh)) / float(max(min(bw, bh), 1)),
        }
    morph = polygon_morphology(cnt_xy, method="minarearect")
    if len(sig_cnts) > 1 or len(hole_cnts) > 0:
        use_outers = sig_cnts if sig_cnts else outer_cnts
        all_pts_mc = np.vstack([c.reshape(-1, 2) for c in use_outers]).astype(np.float64)
        outer_area_sum = float(sum(abs(cv2.contourArea(c)) for c in use_outers))
        hole_area_sum = float(sum(abs(cv2.contourArea(hc)) for hc in hole_cnts))
        net_poly_area = max(1.0, outer_area_sum - hole_area_sum) if hole_cnts else outer_area_sum
        tot_perim_mc = float(sum(cv2.arcLength(c, True) for c in use_outers)) + float(
            sum(cv2.arcLength(hc, True) for hc in hole_cnts)
        )
        hull_mc = cv2.convexHull(all_pts_mc.astype(np.float32).reshape(-1, 1, 2))
        hull_area_mc = float(abs(cv2.contourArea(hull_mc))) if hull_mc is not None and len(hull_mc) >= 3 else 0.0
        eff_sol_num = min(area_px, hull_area_mc) if hole_cnts else net_poly_area
        solidity = (eff_sol_num / hull_area_mc) if hull_area_mc > 0.0 else 0.0
        circularity = (4.0 * math.pi * net_poly_area / (tot_perim_mc * tot_perim_mc)) if tot_perim_mc > 0.0 else 0.0
        pts_for_span = all_pts_mc
    else:
        solidity = float(getattr(morph, "solidity", getattr(morph, "convexity", 0.0)) or 0.0) if morph is not None else 0.0
        circularity = float(getattr(morph, "circularity", 0.0) or 0.0) if morph is not None else 0.0
        pts_for_span = np.asarray(cnt_xy, dtype=np.float64).reshape(-1, 2)
    aspect_ratio = float(getattr(morph, "aspect_ratio", 0.0) or 0.0) if morph is not None else 0.0
    if angle_cw is not None:
        rad_ar = math.radians(float(angle_cw))
        cos_ar, sin_ar = math.cos(rad_ar), math.sin(rad_ar)
        pts_ar = pts_for_span
        xr_ar = pts_ar[:, 0] * cos_ar - pts_ar[:, 1] * sin_ar
        yr_ar = pts_ar[:, 0] * sin_ar + pts_ar[:, 1] * cos_ar
        l_px = float(np.ptp(yr_ar))
        w_px = float(np.ptp(xr_ar))
        if (l_px <= 0.0 or w_px <= 0.0) and morph is not None and getattr(morph, "rotated_bbox", None) is not None:
            rbox = np.asarray(morph.rotated_bbox, dtype=np.float64)
            if rbox.shape == (4, 2):
                l_px, w_px = length_width_from_box(rbox, float(angle_cw))
        if l_px > 0.0 and w_px > 0.0:
            aspect_ratio = float(l_px) / float(max(w_px, 1e-6))
    rect_w = float(getattr(morph, "dim_max", 0.0) or 0.0) if morph is not None else 0.0
    rect_h = float(getattr(morph, "dim_min", 0.0) or 0.0) if morph is not None else 0.0
    if (rect_w <= 0.0 or rect_h <= 0.0) and morph is not None and getattr(morph, "rotated_bbox", None) is not None:
        rbox0 = np.asarray(morph.rotated_bbox, dtype=np.float64)
        if rbox0.shape == (4, 2):
            d01 = float(np.hypot(rbox0[1, 0] - rbox0[0, 0], rbox0[1, 1] - rbox0[0, 1]))
            d12 = float(np.hypot(rbox0[2, 0] - rbox0[1, 0], rbox0[2, 1] - rbox0[1, 1]))
            rect_w, rect_h = max(d01, d12), min(d01, d12)
    if rect_w <= 0.0 or rect_h <= 0.0:
        (_, (rw_cv, rh_cv), _) = cv2.minAreaRect(cnt_xy.astype(np.float32))
        rect_w = float(max(rw_cv, rh_cv, 1.0))
        rect_h = float(max(min(rw_cv, rh_cv), 1.0))
    if aspect_ratio <= 0.0:
        aspect_ratio = float(max(rect_w, rect_h)) / float(max(min(rect_w, rect_h), 1e-6))
    if bbox is not None and len(bbox) == 4:
        ext_box_area = max(1.0, float(max(0, bbox[2] - bbox[0]) * max(0, bbox[3] - bbox[1])))
        tight_area = float(max(bw * bh, 1))
        oriented_area = ext_box_area if ext_box_area > 1.5 * tight_area else max(1.0, min(ext_box_area, rect_w * rect_h))
    else:
        oriented_area = max(1.0, min(float(max(bw * bh, 1)), rect_w * rect_h))
    extent = min(max(area_px / oriented_area, 0.0), 1.0)
    return {
        "area_px": area_px,
        "solidity": min(max(solidity, 0.0), 1.0),
        "circularity": min(max(circularity, 0.0), 1.0),
        "extent": extent,
        "aspect_ratio": aspect_ratio,
    }


def save_leaf_mask_png(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    output_dir: Path,
) -> str:
    """Saves LeafMachine3 binary leaf mask directly to `output_dir/masks/{catalogNumber}_leaf{id}.png`."""
    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    mask_dir = out_base / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    save_path = mask_dir / f"{catalog_number}_leaf{leaf_id}.png"
    binary_mask = np.where(mask > 0, 255, 0).astype(np.uint8)
    ok, png_buf = cv2.imencode(".png", binary_mask)
    if ok:
        save_path.write_bytes(png_buf.tobytes())
    else:
        cv2.imwrite(str(save_path), binary_mask)
    return str(save_path)


def _rotate_contour_to_tip_up(
    resamp_pts: np.ndarray,
    anchor_pt: Optional[Tuple[float, float]] = None,
    p_apex_pt: Optional[Tuple[float, float]] = None,
    keypoints_31: Optional[Any] = None,
) -> np.ndarray:
    """Rotates resampled clockwise contour `(N, 2)` by LeafMachine3 `compute_orientation` so `lamina_base -> lamina_tip` points tip-up `(0, -1)`."""
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
        apex_xy = (
            np.array([float(p_apex_pt[0]), float(p_apex_pt[1])], dtype=np.float64)
            if p_apex_pt is not None
            else np.asarray(resamp_pts[len(resamp_pts) // 2], dtype=np.float64)
        )
        ori_fb = compute_orientation({
            "lamina_tip": (float(apex_xy[0]), float(apex_xy[1])),
            "lamina_base": (float(base_xy[0]), float(base_xy[1])),
        })
        angle_cw = float(ori_fb.angle_cw) if (ori_fb.success and ori_fb.angle_cw is not None) else 0.0

    signed_cw = ((float(angle_cw) + 180.0) % 360.0) - 180.0
    if abs(signed_cw) > 0.25:
        rad = math.radians(float(angle_cw))
        cos_a, sin_a = math.cos(rad), math.sin(rad)
        centroid = np.mean(resamp_pts, axis=0)
        centered = resamp_pts - centroid
        x_rot = centered[:, 0] * cos_a - centered[:, 1] * sin_a + centroid[0]
        y_rot = centered[:, 0] * sin_a + centered[:, 1] * cos_a + centroid[1]
        resamp_pts = np.column_stack([x_rot, y_rot])
    return resamp_pts


def _extract_resampled_contour_pts(
    mask: Optional[np.ndarray],
    num_points: int = 120,
    petiole_attachment_pt: Optional[Tuple[float, float]] = None,
    hole_mask: Optional[np.ndarray] = None,
    lamina_base: Optional[Tuple[float, float]] = None,
    momocs_outline_xy: Optional[np.ndarray] = None,
    contour: Optional[np.ndarray] = None,
    p_apex_pt: Optional[Tuple[float, float]] = None,
    keypoints_31: Optional[Any] = None,
    already_oriented: bool = False,
) -> Optional[np.ndarray]:
    is_preoriented = bool(already_oriented)
    if momocs_outline_xy is not None and len(np.asarray(momocs_outline_xy)) >= 10:
        largest_cnt = np.asarray(momocs_outline_xy, dtype=np.float64).reshape(-1, 2)
        is_preoriented = True
    elif contour is not None and len(contour) >= 10:
        largest_cnt = np.asarray(contour, dtype=np.float64).reshape(-1, 2)
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
            raw_u8 = (arr_in > 0).astype(np.uint8) * 255
            if not is_preoriented and hole_mask is not None and isinstance(hole_mask, np.ndarray) and hole_mask.shape[:2] == raw_u8.shape[:2]:
                raw_u8 = np.where((raw_u8 > 0) & (hole_mask == 0), 255, 0).astype(np.uint8)
            raw_u8, _, _ = open_hairline_sinus_bridges_and_measure_voids(raw_u8, hole_mask=hole_mask)
            extracted_xy = extract_contour(raw_u8)
            if extracted_xy is None or len(extracted_xy) < 10:
                return None
            largest_cnt = np.asarray(extracted_xy, dtype=np.float64).reshape(-1, 2)

    anchor_pt = petiole_attachment_pt if petiole_attachment_pt is not None else lamina_base
    if not is_preoriented and keypoints_31 is not None:
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

    if momocs_outline_xy is not None and petiole_attachment_pt is None and lamina_base is None:
        start_idx = 0
    elif anchor_pt is not None:
        diffs = largest_cnt - np.array([float(anchor_pt[0]), float(anchor_pt[1])], dtype=np.float64)
        start_idx = int(np.argmin(np.hypot(diffs[:, 0], diffs[:, 1])))
    elif is_preoriented:
        y_max_c = float(np.max(largest_cnt[:, 1]))
        y_span_c = max(float(np.ptp(largest_cnt[:, 1])), 1.0)
        x_med_c = float(np.median(largest_cnt[:, 0]))
        bot_mask = largest_cnt[:, 1] >= (y_max_c - 0.08 * y_span_c)
        bot_idx = np.where(bot_mask)[0]
        if len(bot_idx) > 0:
            sub_b = largest_cnt[bot_idx]
            pick_b = int(np.argmin((y_max_c - sub_b[:, 1]) + 0.65 * np.abs(sub_b[:, 0] - x_med_c)))
            start_idx = int(bot_idx[pick_b])
        else:
            start_idx = int(np.argmax(largest_cnt[:, 1]))
    else:
        start_idx = int(np.argmax(largest_cnt[:, 1]))

    pts = np.roll(largest_cnt, -start_idx, axis=0)
    if momocs_signed_area(pts) > 0 or cv2.contourArea(pts.reshape(-1, 1, 2).astype(np.float32), oriented=True) > 0:
        pts = np.vstack([pts[0:1], pts[1:][::-1]])

    dt = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    valid = dt > 1e-7
    if np.sum(valid) < 4:
        return None
    pts = pts[np.concatenate([[True], valid])]
    if len(pts) > 2 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]

    resamp_pts = resample_closed(pts, int(num_points))
    if momocs_signed_area(resamp_pts) > 0 or cv2.contourArea(resamp_pts.reshape(-1, 1, 2).astype(np.float32), oriented=True) > 0:
        resamp_pts = np.vstack([resamp_pts[0:1], resamp_pts[1:][::-1]])
    if is_preoriented:
        return resamp_pts
    return _rotate_contour_to_tip_up(
        resamp_pts,
        anchor_pt=anchor_pt,
        p_apex_pt=p_apex_pt,
        keypoints_31=keypoints_31,
    )


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
    already_oriented: bool = False,
) -> Optional[np.ndarray]:
    """Extracts resampled, tip-up unit-scaled 2D coordinates via LeafMachine3 helpers."""
    _ = petiole_mask
    if n_pts is not None:
        num_points = int(n_pts)
    resamp_pts = _extract_resampled_contour_pts(
        mask=mask,
        num_points=num_points,
        petiole_attachment_pt=petiole_attachment_pt,
        hole_mask=hole_mask,
        lamina_base=lamina_base,
        momocs_outline_xy=momocs_outline_xy,
        p_apex_pt=p_apex_pt,
        keypoints_31=keypoints_31,
        already_oriented=already_oriented,
    )
    if resamp_pts is None:
        return None
    x_resamp, y_resamp = resamp_pts[:, 0], resamp_pts[:, 1]
    min_x, max_x = float(np.min(x_resamp)), float(np.max(x_resamp))
    min_y, max_y = float(np.min(y_resamp)), float(np.max(y_resamp))
    max_span = max(max_x - min_x, max_y - min_y, 1e-6)
    return np.column_stack([
        np.round((x_resamp - min_x) / max_span, 6),
        np.round((y_resamp - min_y) / max_span, 6),
    ])


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
    already_oriented: bool = False,
) -> Any:
    """Exports LeafMachine3 120-vertex clockwise tip-up contour CSV to `output_dir/contours/{catalogNumber}_leaf{id}.csv`."""
    _ = petiole_mask
    empty_ret: Any = (None, None) if return_coords else None
    resamp_pts = _extract_resampled_contour_pts(
        mask=mask,
        num_points=num_points,
        petiole_attachment_pt=petiole_attachment_pt,
        hole_mask=hole_mask,
        lamina_base=lamina_base,
        momocs_outline_xy=momocs_outline_xy,
        contour=contour,
        p_apex_pt=p_apex_pt,
        keypoints_31=keypoints_31,
        already_oriented=already_oriented,
    )
    if resamp_pts is None:
        return empty_ret

    x_resamp, y_resamp = resamp_pts[:, 0], resamp_pts[:, 1]
    x_r = np.round(x_resamp, 2).tolist()
    y_r = np.round(y_resamp, 2).tolist()

    out_base = (output_dir / "data") if (output_dir.name != "data" and (output_dir / "data").is_dir()) else output_dir
    contour_dir = out_base / "contours"
    contour_dir.mkdir(parents=True, exist_ok=True)
    out_csv = contour_dir / f"{catalog_number}_leaf{leaf_id}.csv"
    prefix = f"{catalog_number},{leaf_id},"
    lines = ["catalogNumber,leaf_id,point_index,x,y\n"]
    lines.extend(
        f"{prefix}{idx},{xv},{yv}\n"
        for idx, (xv, yv) in enumerate(zip(x_r, y_r))
    )
    out_csv.write_text("".join(lines), encoding="utf-8")
    if return_coords:
        xy_coords = [[round(float(xv), 4), round(float(yv), 4)] for xv, yv in zip(x_r, y_r)]
        return str(out_csv), xy_coords
    return str(out_csv)


__all__ = [
    "open_hairline_sinus_bridges_and_measure_voids",
    "compute_lm3_shape_metrics",
    "save_leaf_mask_png",
    "resample_normalized_contour",
    "export_standardized_contour",
    "cluster_plant_individuals",
    "_extract_xy_point",
    "_parse_basal_leaf_pose",
    "_estimate_single_plant_crown",
    "estimate_rosette_crown_point",
    "approximate_petiole_and_scape_metrics",
    "classify_leaf_zone_anatomical",
]
