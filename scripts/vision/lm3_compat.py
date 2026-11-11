#!/usr/bin/env python3
"""
===============================================================================
Script: lm3_compat.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline

Description:
    Optional dependency compatibility layer providing pure NumPy / OpenCV
    fallbacks for `sklearn.cluster.DBSCAN` and `leafmachine3` geometric helpers
    when running outside the dedicated `.venv_LM3` environment (e.g., unit tests,
    Step 01 harvesting, or lightweight CPU QA environments).
    When `leafmachine3` IS installed in the active interpreter, all calls
    delegate 100% unchanged to the canonical `leafmachine3` implementation.
===============================================================================
"""

from __future__ import annotations

import importlib.util
import math
import sys
import types
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import numpy as np

try:
    from sklearn.cluster import DBSCAN as _SklearnDBSCAN  # type: ignore[import-untyped]
except ImportError:
    _SklearnDBSCAN = None


class _NumpyDBSCAN:
    """Exact 2D/ND Euclidean DBSCAN fallback matching `sklearn.cluster.DBSCAN(eps, min_samples).fit(X)`."""

    def __init__(self, eps: float = 0.5, min_samples: int = 5, **_kwargs: Any) -> None:
        self.eps = float(eps)
        self.min_samples = int(min_samples)
        self.labels_: np.ndarray = np.empty(0, dtype=np.int64)

    def fit(self, X: Any) -> "_NumpyDBSCAN":
        arr = np.asarray(X, dtype=np.float64)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        n = int(arr.shape[0])
        if n == 0:
            self.labels_ = np.empty(0, dtype=np.int64)
            return self

        diff = arr[:, None, :] - arr[None, :, :]
        dist_mat = np.sqrt(np.sum(diff * diff, axis=-1))
        neighbors = [np.where(dist_mat[i] <= self.eps)[0] for i in range(n)]
        is_core = np.array([len(nb) >= self.min_samples for nb in neighbors], dtype=bool)

        labels = np.full(n, -1, dtype=np.int64)
        visited = np.zeros(n, dtype=bool)
        cluster_id = 0

        for i in range(n):
            if visited[i]:
                continue
            visited[i] = True
            if not is_core[i]:
                continue
            labels[i] = cluster_id
            queue = list(int(x) for x in neighbors[i] if int(x) != i)
            queued_set = set(queue)
            idx = 0
            while idx < len(queue):
                q = queue[idx]
                idx += 1
                if not visited[q]:
                    visited[q] = True
                    if is_core[q]:
                        for nb_pt in neighbors[q]:
                            nb_int = int(nb_pt)
                            if nb_int not in queued_set and labels[nb_int] == -1:
                                queue.append(nb_int)
                                queued_set.add(nb_int)
                if labels[q] == -1:
                    labels[q] = cluster_id
            cluster_id += 1

        self.labels_ = labels
        return self

    def fit_predict(self, X: Any) -> np.ndarray:
        return self.fit(X).labels_


DBSCAN = _SklearnDBSCAN if _SklearnDBSCAN is not None else _NumpyDBSCAN


def _fb_signed_area(pts: np.ndarray) -> float:
    xy = np.asarray(pts, dtype=np.float64).reshape(-1, 2)
    if len(xy) < 3:
        return 0.0
    x = xy[:, 0]
    y = xy[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _fb_extract_contour(mask: np.ndarray) -> Optional[np.ndarray]:
    if mask is None or mask.size == 0:
        return None
    u8 = (np.asarray(mask) > 0).astype(np.uint8) * 255
    cnts, _ = cv2.findContours(u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not cnts:
        return None
    best = max(cnts, key=cv2.contourArea)
    pts = best.reshape(-1, 2).astype(np.float64)
    return pts if len(pts) >= 3 else None


def _fb_prepare_mask(
    mask: np.ndarray,
    min_area: int = 1,
    smooth_kernel: int = 0,
    pad_px: int = 0,
    fill_holes: bool = True,
    largest_only: bool = True,
    **_kwargs: Any,
) -> np.ndarray:
    u8 = (np.asarray(mask) > 0).astype(np.uint8)
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return u8.astype(bool)
    num_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(u8, connectivity=8)
    if num_lbl <= 1:
        return np.zeros_like(u8, dtype=bool)
    areas = stats[1:, cv2.CC_STAT_AREA]
    valid_ids = [1 + idx for idx, a in enumerate(areas) if int(a) >= int(min_area)]
    if not valid_ids:
        valid_ids = [1 + int(np.argmax(areas))]
    if largest_only:
        best_id = max(valid_ids, key=lambda cid: int(stats[cid, cv2.CC_STAT_AREA]))
        sel_u8 = np.where(labels == best_id, 255, 0).astype(np.uint8)
    else:
        sel_u8 = np.where(np.isin(labels, valid_ids), 255, 0).astype(np.uint8)
    if fill_holes:
        cnts, _ = cv2.findContours(sel_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        filled = np.zeros_like(sel_u8)
        if cnts:
            cv2.drawContours(filled, cnts, -1, 255, thickness=cv2.FILLED)
        sel_u8 = filled
    if smooth_kernel and int(smooth_kernel) > 1:
        k = int(smooth_kernel) | 1
        kern = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
        sel_u8 = cv2.morphologyEx(sel_u8, cv2.MORPH_CLOSE, kern)
    if pad_px and int(pad_px) > 0:
        p = int(pad_px)
        sel_u8 = np.pad(sel_u8, ((p, p), (p, p)), mode="constant", constant_values=0)
    return sel_u8 > 0


def _fb_resample_closed(pts: np.ndarray, n_points: int = 120) -> np.ndarray:
    xy = np.asarray(pts, dtype=np.float64).reshape(-1, 2)
    if len(xy) < 2 or n_points <= 0:
        return xy
    closed = np.vstack([xy, xy[0:1]])
    diffs = np.diff(closed, axis=0)
    seg_lens = np.hypot(diffs[:, 0], diffs[:, 1])
    cum = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total = float(cum[-1])
    if total <= 1e-9:
        return np.repeat(xy[0:1], n_points, axis=0)
    targets = np.linspace(0.0, total, int(n_points), endpoint=False)
    rx = np.interp(targets, cum, closed[:, 0])
    ry = np.interp(targets, cum, closed[:, 1])
    return np.column_stack([rx, ry])


def _fb_polygon_morphology(contour_xy: np.ndarray, method: str = "minarearect") -> Optional[Any]:
    pts = np.asarray(contour_xy, dtype=np.float32).reshape(-1, 1, 2)
    if len(pts) < 3:
        return None
    area = float(abs(cv2.contourArea(pts)))
    perim = float(cv2.arcLength(pts, True))
    hull = cv2.convexHull(pts)
    hull_area = float(abs(cv2.contourArea(hull))) if hull is not None and len(hull) >= 3 else 0.0
    solidity = (area / hull_area) if hull_area > 0.0 else 0.0
    circularity = (4.0 * math.pi * area / (perim * perim)) if perim > 0.0 else 0.0
    rect = cv2.minAreaRect(pts)
    rw, rh = float(rect[1][0]), float(rect[1][1])
    dim_max = max(rw, rh)
    dim_min = min(rw, rh)
    box_pts = cv2.boxPoints(rect).astype(np.float64)
    bx, by, bw, bh = cv2.boundingRect(pts.astype(np.int32))
    extent_rect = (area / max(dim_max * dim_min, 1.0)) if dim_max * dim_min > 0 else 0.0
    return SimpleNamespace(
        area_px=area,
        perimeter_px=perim,
        convex_area_px=hull_area,
        convex_hull_area=hull_area,
        solidity=min(max(solidity, 0.0), 1.0),
        convexity=min(max(solidity, 0.0), 1.0),
        circularity=min(max(circularity, 0.0), 1.0),
        extent=min(max(extent_rect, 0.0), 1.0),
        dim_max=dim_max,
        dim_min=dim_min,
        rotated_bbox=box_pts,
        axis_bbox=(bx, by, bw, bh),
        method=method,
    )


def _fb_components(mask: np.ndarray, min_area: int = 1, **_kwargs: Any) -> List[np.ndarray]:
    u8 = (np.asarray(mask) > 0).astype(np.uint8)
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return []
    num_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(u8, connectivity=8)
    comps: List[np.ndarray] = []
    for cid in range(1, num_lbl):
        if int(stats[cid, cv2.CC_STAT_AREA]) >= int(min_area):
            comps.append(labels == cid)
    return comps


def _fb_largest_component(mask: np.ndarray) -> np.ndarray:
    u8 = (np.asarray(mask) > 0).astype(np.uint8)
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return np.zeros_like(u8, dtype=bool)
    num_lbl, labels, stats, _ = cv2.connectedComponentsWithStats(u8, connectivity=8)
    if num_lbl <= 1:
        return np.zeros_like(u8, dtype=bool)
    best = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
    return labels == best


def _fb_outer_polygons(mask: np.ndarray) -> List[np.ndarray]:
    u8 = (np.asarray(mask) > 0).astype(np.uint8) * 255
    cnts, _ = cv2.findContours(u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    return [c.reshape(-1, 2).astype(np.float64) for c in cnts if len(c) >= 3]


def _fb_polygon_mask(polygons: Sequence[np.ndarray], shape: Tuple[int, int]) -> np.ndarray:
    h, w = int(shape[0]), int(shape[1])
    canvas = np.zeros((h, w), dtype=np.uint8)
    valid = [np.round(np.asarray(p)).astype(np.int32).reshape(-1, 1, 2) for p in polygons if len(p) >= 3]
    if valid:
        cv2.fillPoly(canvas, valid, 1)
    return canvas > 0


def _fb_mask_bbox(mask: np.ndarray) -> Optional[Tuple[int, int, int, int]]:
    if mask is None or mask.size == 0:
        return None
    ys, xs = np.where(np.asarray(mask) > 0)
    if len(xs) == 0:
        return None
    return int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1


def _fb_crop_to_box(
    arr: Optional[np.ndarray],
    bbox: Tuple[int, int, int, int],
    pad: int = 0,
    **_kwargs: Any,
) -> Optional[np.ndarray]:
    if arr is None:
        return None
    x0, y0, x1, y1 = [int(v) for v in bbox]
    p = max(0, int(pad))
    h, w = arr.shape[:2]
    x0 = max(0, min(w, x0 - p))
    x1 = max(x0, min(w, x1 + p))
    y0 = max(0, min(h, y0 - p))
    y1 = max(y0, min(h, y1 + p))
    return arr[y0:y1, x0:x1].copy()


def _fb_oriented_affine(shape_hw: Tuple[int, int], angle_cw: float) -> Tuple[np.ndarray, Tuple[int, int]]:
    h, w = int(shape_hw[0]), int(shape_hw[1])
    ang = float(angle_cw) % 360.0
    if abs(ang) <= 1e-6 or abs(ang - 360.0) <= 1e-6:
        return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64), (h, w)
    cx, cy = (w - 1) * 0.5, (h - 1) * 0.5
    # In image coords (x right, y down), OpenCV positive angle is CCW, so -angle_cw is CW
    mat = cv2.getRotationMatrix2D((cx, cy), -ang, 1.0)
    cos_v = abs(float(mat[0, 0]))
    sin_v = abs(float(mat[0, 1]))
    new_w = max(1, int(math.ceil(w * cos_v + h * sin_v)))
    new_h = max(1, int(math.ceil(w * sin_v + h * cos_v)))
    mat[0, 2] += (new_w - 1) * 0.5 - cx
    mat[1, 2] += (new_h - 1) * 0.5 - cy
    return mat.astype(np.float64), (new_h, new_w)


def _fb_apply_affine(mat: np.ndarray, pts: np.ndarray) -> np.ndarray:
    arr = np.asarray(pts, dtype=np.float64)
    if arr.size == 0:
        return arr.reshape(0, 2)
    m = np.asarray(mat, dtype=np.float64)
    xy = arr.reshape(-1, 2)
    out = (xy @ m[:, :2].T) + m[:, 2]
    return out.reshape(arr.shape)


def _fb_rotate_image(
    img: np.ndarray,
    angle_cw: float,
    bg: Any = 0,
    nearest: bool = True,
) -> np.ndarray:
    h, w = int(img.shape[0]), int(img.shape[1])
    mat, (new_h, new_w) = _fb_oriented_affine((h, w), angle_cw)
    interp = cv2.INTER_NEAREST if nearest else cv2.INTER_LINEAR
    border_val = (int(bg),) * (img.shape[2] if img.ndim == 3 else 1) if np.isscalar(bg) else tuple(bg)
    return cv2.warpAffine(
        img,
        mat,
        (new_w, new_h),
        flags=interp,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=border_val,
    )


def _fb_compute_orientation(
    keypoints_or_mask: Any = None,
    p_base: Optional[Tuple[float, float]] = None,
    p_apex: Optional[Tuple[float, float]] = None,
    **_kwargs: Any,
) -> Any:
    if isinstance(keypoints_or_mask, dict):
        kp = keypoints_or_mask
        if p_base is None:
            for k in ("lamina_base", "base_center", "midvein_14"):
                if k in kp and kp[k] is not None and len(kp[k]) >= 2:
                    p_base = (float(kp[k][0]), float(kp[k][1]))
                    break
        if p_apex is None:
            for k in ("lamina_tip", "apex_center", "midvein_0"):
                if k in kp and kp[k] is not None and len(kp[k]) >= 2:
                    p_apex = (float(kp[k][0]), float(kp[k][1]))
                    break
    if p_base is not None and p_apex is not None:
        dx = float(p_apex[0] - p_base[0])
        dy = float(p_apex[1] - p_base[1])
        if math.hypot(dx, dy) > 1e-6:
            angle_cw = (-math.degrees(math.atan2(dx, -dy))) % 360.0
            return SimpleNamespace(
                angle_cw=angle_cw,
                angle_deg=math.degrees(math.atan2(dy, dx)),
                p_base=p_base,
                p_apex=p_apex,
                success=True,
                method="landmarks",
            )
    if isinstance(keypoints_or_mask, np.ndarray) and keypoints_or_mask.ndim >= 2:
        pts = _fb_extract_contour(keypoints_or_mask)
        if pts is not None and len(pts) >= 3:
            m = cv2.moments((keypoints_or_mask > 0).astype(np.uint8), binaryImage=True)
            if m["m00"] > 0:
                theta = 0.5 * math.atan2(2.0 * m["mu11"], m["mu20"] - m["mu02"])
                vx, vy = math.cos(theta), math.sin(theta)
                angle_cw = (-math.degrees(math.atan2(vx, -vy))) % 360.0
                return SimpleNamespace(
                    angle_cw=angle_cw,
                    angle_deg=math.degrees(theta),
                    p_base=None,
                    p_apex=None,
                    success=True,
                    method="pca_moments",
                )
    return SimpleNamespace(angle_cw=0.0, angle_deg=-90.0, p_base=p_base, p_apex=p_apex, success=False, method="fallback")


def _fb_length_width_from_box(rotated_bbox: np.ndarray, angle_cw: float) -> Tuple[float, float]:
    pts = np.asarray(rotated_bbox, dtype=np.float64).reshape(-1, 2)
    rad = math.radians(float(angle_cw))
    cos_a, sin_a = math.cos(rad), math.sin(rad)
    xr = pts[:, 0] * cos_a - pts[:, 1] * sin_a
    yr = pts[:, 0] * sin_a + pts[:, 1] * cos_a
    length_px = float(np.max(yr) - np.min(yr))
    width_px = float(np.max(xr) - np.min(xr))
    return max(length_px, 0.0), max(width_px, 0.0)


@dataclass
class _FbOrientedLeaf:
    specimen_id: int
    leaf_id: int
    detection_id: int
    stem: str
    silhouette: np.ndarray
    lamina: np.ndarray
    holes: np.ndarray
    petiole: Optional[np.ndarray]
    kpts: Dict[str, Tuple[float, float]]
    kpt_conf: Dict[str, float]
    angle_cw: float
    det_box: Tuple[int, int, int, int]
    crop_shape: Tuple[int, int]
    truncated: bool = False


def _fb_measure_leaf(oleaf: Any, min_score: float = 0.70) -> Optional[Any]:
    lam = np.asarray(oleaf.lamina, dtype=bool)
    if lam.size == 0 or int(np.count_nonzero(lam)) == 0:
        return None
    h, w = lam.shape[:2]
    num_lbl, _, stats, _ = cv2.connectedComponentsWithStats(lam.astype(np.uint8), connectivity=8)
    total_fg = float(np.count_nonzero(lam))
    largest_fg = float(np.max(stats[1:, cv2.CC_STAT_AREA])) if num_lbl > 1 else total_fg
    largest_frac = largest_fg / max(total_fg, 1.0)

    kpts = getattr(oleaf, "kpts", None) or {}
    mid_xs: List[float] = []
    for k in ("lamina_tip", "lamina_base"):
        if k in kpts and kpts[k] is not None:
            mid_xs.append(float(kpts[k][0]))
    for i in range(15):
        mk = f"midvein_{i}"
        if mk in kpts and kpts[mk] is not None:
            mid_xs.append(float(kpts[mk][0]))

    ys_fg, xs_fg = np.where(lam)
    x_axis = float(np.median(mid_xs)) if mid_xs else float(np.mean(xs_fg))

    # Reflect foreground pixels across vertical midrib axis x = x_axis
    xs_refl = np.round(2.0 * x_axis - xs_fg.astype(np.float64)).astype(np.int64)
    in_bounds = (xs_refl >= 0) & (xs_refl < w)
    refl_mask = np.zeros_like(lam, dtype=bool)
    refl_mask[ys_fg[in_bounds], xs_refl[in_bounds]] = True

    inter = float(np.count_nonzero(lam & refl_mask))
    denom = float(np.count_nonzero(lam)) + float(np.count_nonzero(refl_mask))
    dice = (2.0 * inter / denom) if denom > 0.0 else 0.0

    area_left = float(np.count_nonzero(xs_fg < x_axis))
    area_right = float(np.count_nonzero(xs_fg > x_axis))
    lr_sum = area_left + area_right
    si_a = ((area_right - area_left) / lr_sum) if lr_sum > 0.0 else 0.0

    archetype_score = float(
        np.clip(0.65 * dice + 0.35 * (1.0 - min(1.0, abs(si_a))) * largest_frac, 0.0, 1.0)
    )
    gates_pass = bool(dice >= 0.70 and abs(si_a) < 0.35 and largest_frac >= 0.85)
    is_archetypal = bool(gates_pass and archetype_score >= float(min_score))
    return SimpleNamespace(
        dice=dice,
        si_a=si_a,
        largest_frac=largest_frac,
        archetype_score=archetype_score,
        gates_pass=gates_pass,
        is_archetypal=is_archetypal,
    )


def _fb_ect_product_for(include_petiole: bool = False, include_holes: bool = True) -> str:
    if include_petiole and include_holes:
        return "laminapet"
    if include_petiole and not include_holes:
        return "silhouettepet"
    if not include_petiole and include_holes:
        return "lamina"
    return "silhouette"


def _fb_render_leaf_products(
    _rgb: np.ndarray,
    masks_dict: Dict[str, Optional[np.ndarray]],
    angle_cw: Optional[float] = None,
    want: Optional[Sequence[str]] = None,
) -> Dict[str, Optional[np.ndarray]]:
    keys = list(want) if want is not None else list(masks_dict.keys())
    out: Dict[str, Optional[np.ndarray]] = {}
    ang = float(angle_cw) if angle_cw is not None else 0.0
    for k in keys:
        m = masks_dict.get(k)
        if m is None:
            out[k] = None
            continue
        u8 = (np.asarray(m) > 0).astype(np.uint8) * 255
        if abs(ang % 360.0) > 0.25:
            u8 = _fb_rotate_image(u8, ang, bg=0, nearest=True)
        out[k] = u8 > 0
    return out


def _fb_compute_ect(
    mask_u8: np.ndarray,
    num_dirs: int = 32,
    bound_radius: float = 1.0,
    want_simple: bool = False,
    num_thresh: int = 32,
) -> Optional[Any]:
    u8 = (np.asarray(mask_u8) > 0).astype(np.uint8) * 255
    if u8.size == 0 or int(np.count_nonzero(u8)) == 0:
        return None
    cnts_raw, hier = cv2.findContours(u8, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    if not cnts_raw or hier is None:
        return None
    cycles: List[Tuple[float, np.ndarray]] = []
    for idx_c, c_raw in enumerate(cnts_raw):
        c_xy = c_raw.reshape(-1, 2).astype(np.float64)
        if len(c_xy) < 3 or abs(float(cv2.contourArea(c_raw))) < 4.0:
            continue
        sign_c = 1.0 if int(hier[0][idx_c][3]) < 0 else -1.0
        n_res = min(256, max(48, len(c_xy)))
        cycles.append((sign_c, _fb_resample_closed(c_xy, n_points=n_res)))
    if not cycles:
        cnt_xy = _fb_extract_contour(u8)
        if cnt_xy is None or len(cnt_xy) < 3:
            return None
        cycles = [(1.0, _fb_resample_closed(cnt_xy, n_points=min(256, max(64, len(cnt_xy)))))]

    outer_stack = np.vstack([pts_c for s_c, pts_c in cycles if s_c > 0] or [cycles[0][1]])
    centroid = np.mean(outer_stack, axis=0)
    centered_outer = outer_stack - centroid
    max_r = float(np.max(np.hypot(centered_outer[:, 0], centered_outer[:, 1])))
    scale = (float(bound_radius) * 0.95 / max_r) if max_r > 1e-9 else 1.0
    outline_norm = ((cycles[0][1] - centroid) * scale).astype(np.float32)

    n_dirs = max(4, int(num_dirs))
    n_thr = max(4, int(num_thresh))
    thetas = np.linspace(0.0, 2.0 * np.pi, n_dirs, endpoint=False, dtype=np.float32)
    thresholds = np.linspace(-float(bound_radius), float(bound_radius), n_thr, dtype=np.float32)

    # Multi-component polygonal Euler Characteristic Transform: sum_c sign_c * (#V_c(<=t) - #E_c(<=t))
    dirs = np.stack([np.cos(thetas), np.sin(thetas)], axis=1)  # (D, 2)
    ect_acc = np.zeros((n_thr, n_dirs), dtype=np.float32)
    for sign_c, pts_c in cycles:
        c_norm = ((pts_c - centroid) * scale).astype(np.float32)
        heights = c_norm @ dirs.T  # (N_c, D)
        heights_next = np.roll(heights, -1, axis=0)
        edge_max = np.maximum(heights, heights_next)
        v_le = (heights[None, :, :] <= thresholds[:, None, None]).sum(axis=1)
        e_le = (edge_max[None, :, :] <= thresholds[:, None, None]).sum(axis=1)
        ect_acc += np.float32(sign_c) * (v_le - e_le).astype(np.float32)
    ect_mat = ect_acc.astype(np.float32)
    return SimpleNamespace(
        ect_matrix=ect_mat,
        thetas=thetas,
        thresholds=thresholds,
        bound_radius=float(bound_radius),
        outline_norm=outline_norm,
        want_simple=bool(want_simple),
    )


def _fb_tip_up_directions(ect_raw: np.ndarray, thetas: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mat = np.asarray(ect_raw, dtype=np.float32)
    th = np.asarray(thetas, dtype=np.float32)
    if mat.ndim != 2 or len(th) != mat.shape[1] or len(th) == 0:
        return mat, th
    # Tip-up direction in image coordinates is (0, -1) -> angle 3*pi/2 (or -pi/2)
    target = 1.5 * math.pi
    diffs = np.abs(np.angle(np.exp(1j * (th - target))))
    shift = int(np.argmin(diffs))
    return np.roll(mat, -shift, axis=1), np.roll(th, -shift)


def _fb_visual_log(mat: np.ndarray) -> np.ndarray:
    arr = np.asarray(mat, dtype=np.float32)
    return (np.sign(arr) * np.log1p(np.abs(arr))).astype(np.float32)


def _fb_render_radial_ect_overlay(
    ect_vis: np.ndarray,
    thetas: np.ndarray,
    thresholds: np.ndarray,
    bound_radius: float,
    outline_norm: Optional[np.ndarray] = None,
    cmap: str = "viridis",
    size: int = 256,
) -> np.ndarray:
    canvas = np.full((size, size, 3), 248, dtype=np.uint8)
    mat = np.asarray(ect_vis, dtype=np.float32)
    if mat.size > 0:
        yy, xx = np.mgrid[0:size, 0:size].astype(np.float32)
        cx = cy = (size - 1) * 0.5
        rad_max = size * 0.44
        dx = xx - cx
        dy = yy - cy
        rr = np.hypot(dx, dy) / max(rad_max, 1.0)
        inside = rr <= 1.0
        ang = np.mod(np.arctan2(dy, dx), 2.0 * np.pi)
        t_idx = np.clip((rr * (mat.shape[0] - 1)).astype(np.int32), 0, mat.shape[0] - 1)
        d_idx = np.clip((ang / (2.0 * np.pi) * mat.shape[1]).astype(np.int32), 0, mat.shape[1] - 1)
        vals = mat[t_idx, d_idx]
        vmin, vmax = float(np.min(vals)), float(np.max(vals))
        norm_u8 = (
            np.round((vals - vmin) / max(vmax - vmin, 1e-6) * 255.0).astype(np.uint8)
            if vmax > vmin
            else np.full_like(t_idx, 128, dtype=np.uint8)
        )
        cv_cmap = getattr(cv2, f"COLORMAP_{str(cmap).upper()}", cv2.COLORMAP_VIRIDIS)
        colored = cv2.applyColorMap(norm_u8, cv_cmap)
        canvas[inside] = colored[inside]
    if outline_norm is not None and len(outline_norm) >= 3:
        cx = cy = (size - 1) * 0.5
        rad_max = size * 0.44 / max(float(bound_radius), 1e-6)
        pts_px = np.round(np.asarray(outline_norm, dtype=np.float32) * rad_max + np.array([cx, cy])).astype(np.int32)
        cv2.polylines(canvas, [pts_px.reshape(-1, 1, 2)], True, (255, 255, 255), 2, cv2.LINE_AA)
    return canvas


def _fb_measure_petiole(
    pet_u8: np.ndarray,
    lam_u8: np.ndarray,
    pet_chain: Sequence[Tuple[float, float]],
) -> Optional[Any]:
    p_bin = (np.asarray(pet_u8) > 0).astype(np.uint8)
    l_bin = (np.asarray(lam_u8) > 0).astype(np.uint8)
    if p_bin.size == 0 or int(np.count_nonzero(p_bin)) == 0:
        return None
    arc_len = 0.0
    for i in range(len(pet_chain) - 1):
        arc_len += math.hypot(
            float(pet_chain[i + 1][0]) - float(pet_chain[i][0]),
            float(pet_chain[i + 1][1]) - float(pet_chain[i][1]),
        )
    dt = cv2.distanceTransform(p_bin, cv2.DIST_L2, 3)
    h, w = p_bin.shape[:2]
    radii: List[float] = []
    for px, py in pet_chain:
        ix = max(0, min(w - 1, int(round(float(px)))))
        iy = max(0, min(h - 1, int(round(float(py)))))
        if dt[iy, ix] > 0.0:
            radii.append(2.0 * float(dt[iy, ix]))
    if not radii:
        nz = dt[dt > 0]
        width_px = float(2.0 * np.median(nz)) if nz.size > 0 else 1.0
    else:
        width_px = float(np.median(radii))
    touches = bool(
        l_bin.shape == p_bin.shape
        and np.any(cv2.dilate(p_bin, np.ones((3, 3), dtype=np.uint8)) & l_bin)
    )
    return SimpleNamespace(width_px=width_px, length_px=arc_len, touches_leaf=touches)


def _fb_compute_measurements(pts_2d: Dict[str, Tuple[float, float]]) -> Optional[Any]:
    if not pts_2d:
        return None
    lt = pts_2d.get("lamina_tip") or pts_2d.get("midvein_0")
    lb = pts_2d.get("lamina_base") or pts_2d.get("midvein_14")
    tip_base_len = (
        math.hypot(float(lt[0]) - float(lb[0]), float(lt[1]) - float(lb[1]))
        if (lt is not None and lb is not None)
        else None
    )
    mv_pts: List[Tuple[float, float]] = []
    if lt is not None:
        mv_pts.append((float(lt[0]), float(lt[1])))
    for i in range(15):
        mk = f"midvein_{i}"
        if mk in pts_2d and pts_2d[mk] is not None:
            pt = (float(pts_2d[mk][0]), float(pts_2d[mk][1]))
            if not mv_pts or math.hypot(pt[0] - mv_pts[-1][0], pt[1] - mv_pts[-1][1]) > 1e-6:
                mv_pts.append(pt)
    if lb is not None:
        pt_b = (float(lb[0]), float(lb[1]))
        if not mv_pts or math.hypot(pt_b[0] - mv_pts[-1][0], pt_b[1] - mv_pts[-1][1]) > 1e-6:
            mv_pts.append(pt_b)
    lam_trace = (
        sum(
            math.hypot(mv_pts[i + 1][0] - mv_pts[i][0], mv_pts[i + 1][1] - mv_pts[i][1])
            for i in range(len(mv_pts) - 1)
        )
        if len(mv_pts) >= 2
        else tip_base_len
    )
    pet_pts: List[Tuple[float, float]] = []
    if lb is not None:
        pet_pts.append((float(lb[0]), float(lb[1])))
    for k in ("petiole_0", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_tip"):
        if k in pts_2d and pts_2d[k] is not None:
            pt = (float(pts_2d[k][0]), float(pts_2d[k][1]))
            if not pet_pts or math.hypot(pt[0] - pet_pts[-1][0], pt[1] - pet_pts[-1][1]) > 1e-6:
                pet_pts.append(pt)
    pet_trace = (
        sum(
            math.hypot(pet_pts[i + 1][0] - pet_pts[i][0], pet_pts[i + 1][1] - pet_pts[i][1])
            for i in range(len(pet_pts) - 1)
        )
        if len(pet_pts) >= 2
        else None
    )
    w1 = pts_2d.get("width_1") or pts_2d.get("width_left")
    w2 = pts_2d.get("width_2") or pts_2d.get("width_right")
    leaf_w = (
        math.hypot(float(w2[0]) - float(w1[0]), float(w2[1]) - float(w1[1]))
        if (w1 is not None and w2 is not None)
        else None
    )
    curv = (lam_trace / tip_base_len) if (lam_trace and tip_base_len and tip_base_len > 1e-6) else 1.0
    return SimpleNamespace(
        lamina_trace_length_px=lam_trace,
        lamina_tip_base_length_px=tip_base_len,
        petiole_trace_length_px=pet_trace,
        leaf_width_px=leaf_w,
        apex_angle_deg=None,
        base_angle_deg=None,
        lamina_curvature=curv,
    )


def _register_virtual_leafmachine3_modules() -> None:
    """Register fallback `leafmachine3.*` submodules in `sys.modules` ONLY when `leafmachine3` is absent."""
    try:
        if importlib.util.find_spec("leafmachine3") is not None:
            return
    except Exception:
        pass
    if "leafmachine3" in sys.modules:
        return

    def _mk_pkg(name: str) -> types.ModuleType:
        m = types.ModuleType(name)
        m.__path__ = []  # type: ignore[attr-defined]
        sys.modules[name] = m
        return m

    def _mk_mod(parent: types.ModuleType, attr: str, full_name: str, **attrs: Any) -> types.ModuleType:
        m = types.ModuleType(full_name)
        for k, v in attrs.items():
            setattr(m, k, v)
        sys.modules[full_name] = m
        setattr(parent, attr, m)
        return m

    root_pkg = _mk_pkg("leafmachine3")
    core_pkg = _mk_pkg("leafmachine3.core")
    inf_pkg = _mk_pkg("leafmachine3.inference")
    rep_pkg = _mk_pkg("leafmachine3.reporting")
    root_pkg.core = core_pkg  # type: ignore[attr-defined]
    root_pkg.inference = inf_pkg  # type: ignore[attr-defined]
    root_pkg.reporting = rep_pkg  # type: ignore[attr-defined]

    _mk_mod(core_pkg, "leaf_outline", "leafmachine3.core.leaf_outline", extract_contour=_fb_extract_contour)
    _mk_mod(
        core_pkg,
        "momocs_prep",
        "leafmachine3.core.momocs_prep",
        prepare_mask=_fb_prepare_mask,
        resample_closed=_fb_resample_closed,
        signed_area=_fb_signed_area,
    )
    _mk_mod(core_pkg, "morphometrics", "leafmachine3.core.morphometrics", polygon_morphology=_fb_polygon_morphology)
    _mk_mod(
        core_pkg,
        "orientation",
        "leafmachine3.core.orientation",
        compute_orientation=_fb_compute_orientation,
        length_width_from_box=_fb_length_width_from_box,
    )
    _mk_mod(
        core_pkg,
        "imaging",
        "leafmachine3.core.imaging",
        polygon_mask=_fb_polygon_mask,
        mask_bbox=_fb_mask_bbox,
        crop_to_box=_fb_crop_to_box,
        rotate_image=_fb_rotate_image,
    )
    _mk_mod(
        core_pkg,
        "bilateral",
        "leafmachine3.core.bilateral",
        OrientedLeaf=_FbOrientedLeaf,
        oriented_affine=_fb_oriented_affine,
        apply_affine=_fb_apply_affine,
        measure_leaf=_fb_measure_leaf,
    )
    _mk_mod(
        core_pkg,
        "ect_compute",
        "leafmachine3.core.ect_compute",
        ect_product_for=_fb_ect_product_for,
        compute_ect=_fb_compute_ect,
        tip_up_directions=_fb_tip_up_directions,
    )
    _mk_mod(core_pkg, "petiole", "leafmachine3.core.petiole", measure_petiole=_fb_measure_petiole)
    _mk_mod(
        core_pkg,
        "landmark_metrics",
        "leafmachine3.core.landmark_metrics",
        compute_measurements=_fb_compute_measurements,
    )
    _mk_mod(
        inf_pkg,
        "leaf_segmenter",
        "leafmachine3.inference.leaf_segmenter",
        _components=_fb_components,
        _largest_component=_fb_largest_component,
        _outer_polygons=_fb_outer_polygons,
    )
    _mk_mod(
        rep_pkg,
        "leaf_products",
        "leafmachine3.reporting.leaf_products",
        render_leaf_products=_fb_render_leaf_products,
    )
    _mk_mod(
        rep_pkg,
        "ect_viz",
        "leafmachine3.reporting.ect_viz",
        visual_log=_fb_visual_log,
        render_radial_ect_overlay=_fb_render_radial_ect_overlay,
    )


_register_virtual_leafmachine3_modules()

from leafmachine3.core.bilateral import OrientedLeaf, apply_affine, measure_leaf, oriented_affine  # noqa: E402
from leafmachine3.core.ect_compute import compute_ect, ect_product_for, tip_up_directions  # noqa: E402
from leafmachine3.core.imaging import crop_to_box, mask_bbox, polygon_mask, rotate_image  # noqa: E402
from leafmachine3.core.landmark_metrics import compute_measurements  # noqa: E402
from leafmachine3.core.leaf_outline import extract_contour  # noqa: E402
from leafmachine3.core.momocs_prep import (  # noqa: E402
    prepare_mask,
    resample_closed,
    signed_area as momocs_signed_area,
)
from leafmachine3.core.morphometrics import polygon_morphology  # noqa: E402
from leafmachine3.core.orientation import compute_orientation, length_width_from_box  # noqa: E402
from leafmachine3.core.petiole import measure_petiole  # noqa: E402
from leafmachine3.inference.leaf_segmenter import (  # noqa: E402
    _components,
    _largest_component,
    _outer_polygons,
)
from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log  # noqa: E402
from leafmachine3.reporting.leaf_products import render_leaf_products  # noqa: E402
