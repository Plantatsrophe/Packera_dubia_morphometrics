#!/usr/bin/env python3
"""
Step 02 calibrated CIELAB blade photometry, physical-mm sinus-bridge harmonization,
and corymb cluster capitulum deconvolution utilities (`scripts/vision/step02_photometry_corymb_utils.py`).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np


def resolve_effective_px_per_mm(
    pixels_per_mm: Optional[float],
    sheet_shape_hw: Tuple[int, int],
    fallback_px_per_mm: float = 10.0,
) -> Tuple[float, bool]:
    """
    Resolve active linear pixel density `(px/mm, is_empirical)` using verified
    ruler calibration when valid, else LeafMachine3's megapixel prior
    (`2.7112084 * sqrt(MP)` when `MP >= 0.5`) or `fallback_px_per_mm`.
    """
    if pixels_per_mm is not None:
        try:
            ppm = float(pixels_per_mm)
            if np.isfinite(ppm) and ppm > 0.0:
                return (ppm, True)
        except (TypeError, ValueError):
            pass

    if sheet_shape_hw is not None and len(sheet_shape_hw) >= 2:
        h, w = max(0, int(sheet_shape_hw[0])), max(0, int(sheet_shape_hw[1]))
    else:
        h, w = 0, 0
    mp = (float(h) * float(w)) / 1e6
    if mp >= 0.5:
        return (float(2.7112084 * math.sqrt(mp)), False)
    return (float(fallback_px_per_mm), False)


def estimate_sheet_paper_cie_l_ref_and_ab(
    image_bgr: Optional[np.ndarray],
) -> Tuple[float, float, float]:
    """
    Estimate herbarium mounting-paper CIE L*, a*, b* reference from an interior
    border band (3%..15% from edges), filtering for low-chroma paper pixels.
    """
    if (
        image_bgr is None
        or not isinstance(image_bgr, np.ndarray)
        or image_bgr.ndim != 3
        or image_bgr.shape[2] < 3
        or image_bgr.size == 0
    ):
        return (90.0, 0.0, 2.0)

    sample = image_bgr[::8, ::8, :3]
    if sample.size == 0:
        sample = image_bgr[:, :, :3]
    lab = cv2.cvtColor(sample, cv2.COLOR_BGR2LAB)
    l_cie = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    a_c = lab[:, :, 1].astype(np.float32) - 128.0
    b_c = lab[:, :, 2].astype(np.float32) - 128.0
    chroma = np.hypot(a_c, b_c)

    sh, sw = l_cie.shape
    if sh >= 24 and sw >= 24:
        r0, r1 = max(1, int(round(0.03 * sh))), max(3, int(round(0.15 * sh)))
        c0, c1 = max(1, int(round(0.03 * sw))), max(3, int(round(0.15 * sw)))
        yy, xx = np.ogrid[:sh, :sw]
        outer_ok = (yy >= r0) & (yy < sh - r0) & (xx >= c0) & (xx < sw - c0)
        core_mask = (yy >= r1) & (yy < sh - r1) & (xx >= c1) & (xx < sw - c1)
        border_band = outer_ok & ~core_mask
    else:
        border_band = np.ones((sh, sw), dtype=bool)

    paper_mask = border_band & (l_cie >= 70.0) & (l_cie <= 99.0) & (chroma <= 12.0)
    if int(np.count_nonzero(paper_mask)) < 32:
        paper_mask = (l_cie >= 65.0) & (l_cie <= 99.0) & (chroma <= 18.0)

    if int(np.count_nonzero(paper_mask)) >= 16:
        l_ref = float(np.percentile(l_cie[paper_mask], 90.0))
        a_ref = float(np.median(a_c[paper_mask]))
        b_ref = float(np.median(b_c[paper_mask]))
    else:
        l_ref = float(np.percentile(l_cie, 85.0))
        a_ref, b_ref = 0.0, 2.0
    return (float(np.clip(l_ref, 65.0, 99.0)), a_ref, b_ref)


def estimate_sheet_paper_cie_l_ref(image_bgr: Optional[np.ndarray]) -> float:
    """Estimate mounting-paper CIE L* white-point reference (0..100 scale)."""
    return estimate_sheet_paper_cie_l_ref_and_ab(image_bgr)[0]


def _cie_l_to_rel_y(l_val: float) -> float:
    """Convert CIELAB L* lightness (0..100) to relative luminance Y/Y_n."""
    return ((l_val + 16.0) / 116.0) ** 3 if l_val > 8.0 else l_val / 903.3


def compute_calibrated_blade_photometry(
    sheet_bgr: np.ndarray,
    blade_mask: np.ndarray,
    mask_origin: Tuple[int, int] = (0, 0),
    paper_l_ref: Optional[float] = None,
    petiole_mask: Optional[np.ndarray] = None,
    hole_mask: Optional[np.ndarray] = None,
    pixels_per_mm: Optional[float] = None,
) -> Dict[str, float]:
    """
    Compute paper-adapted lamina CIELAB photometry (`blade_cie_l_norm`,
    `blade_cie_l_sd`, `blade_cie_a_delta`, `blade_cie_b_delta`) using physical-mm
    inward margin erosion, mounting-strap exclusion, and linearized Y/Y_n scaling.
    """
    nan_out: Dict[str, float] = {
        "blade_cie_l_norm": np.nan,
        "blade_cie_l_sd": np.nan,
        "blade_cie_a_delta": np.nan,
        "blade_cie_b_delta": np.nan,
    }
    if (
        sheet_bgr is None
        or blade_mask is None
        or not isinstance(sheet_bgr, np.ndarray)
        or not isinstance(blade_mask, np.ndarray)
        or sheet_bgr.ndim != 3
        or sheet_bgr.shape[2] < 3
        or blade_mask.ndim < 2
        or blade_mask.size == 0
    ):
        return nan_out

    sh, sw = sheet_bgr.shape[:2]
    mh, mw = blade_mask.shape[:2]
    ox, oy = int(mask_origin[0]), int(mask_origin[1])
    x0, y0 = max(0, ox), max(0, oy)
    x1, y1 = min(sw, ox + mw), min(sh, oy + mh)
    if x1 <= x0 or y1 <= y0:
        return nan_out

    mx0, my0 = x0 - ox, y0 - oy
    mx1, my1 = mx0 + (x1 - x0), my0 + (y1 - y0)
    crop_bgr = sheet_bgr[y0:y1, x0:x1, :3]
    lam_slice = blade_mask[my0:my1, mx0:mx1]
    lam = (lam_slice[:, :, 0] > 0) if lam_slice.ndim == 3 else (lam_slice > 0)

    for excl_mask in (hole_mask, petiole_mask):
        if excl_mask is not None and isinstance(excl_mask, np.ndarray) and excl_mask.shape[:2] == (mh, mw):
            ex_slice = excl_mask[my0:my1, mx0:mx1]
            ex_bool = (ex_slice[:, :, 0] > 0) if ex_slice.ndim == 3 else (ex_slice > 0)
            lam = lam & ~ex_bool

    if int(np.count_nonzero(lam)) < 9:
        return nan_out

    eff_px_per_mm, _ = resolve_effective_px_per_mm(pixels_per_mm, sheet_bgr.shape[:2])
    dist = cv2.distanceTransform(lam.astype(np.uint8), cv2.DIST_L2, 3)
    max_d = float(dist.max())
    if max_d <= 0.0:
        return nan_out

    raw_px = float(np.clip(0.40 * eff_px_per_mm, 2.0, 12.0))
    margin_px = min(raw_px, max(1.0, 0.35 * max_d))
    interior = dist >= margin_px
    if int(np.count_nonzero(interior)) < 16:
        interior = dist >= 1.0
    if int(np.count_nonzero(interior)) < 9:
        interior = lam

    auto_l, paper_a, paper_b = estimate_sheet_paper_cie_l_ref_and_ab(sheet_bgr)
    eff_paper_l = auto_l
    if paper_l_ref is not None:
        try:
            p_val = float(paper_l_ref)
            if np.isfinite(p_val) and p_val >= 50.0:
                eff_paper_l = p_val
        except (TypeError, ValueError):
            pass

    lab = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2LAB)
    l_cie = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    a_c = lab[:, :, 1].astype(np.float32) - 128.0
    b_c = lab[:, :, 2].astype(np.float32) - 128.0
    chroma = np.hypot(a_c, b_c)

    if eff_paper_l >= 72.0:
        is_strap = (l_cie >= (eff_paper_l - 3.5)) & (chroma <= 6.0)
        non_strap = interior & (~is_strap)
        if int(np.count_nonzero(non_strap)) >= 16:
            interior = non_strap
    if int(np.count_nonzero(interior)) < 9:
        return nan_out

    med_l = float(np.median(l_cie[interior]))
    med_a = float(np.median(a_c[interior]))
    med_b = float(np.median(b_c[interior]))

    y_blade = _cie_l_to_rel_y(med_l)
    y_paper = max(_cie_l_to_rel_y(eff_paper_l), 1e-4)
    y_ratio = max(0.0, y_blade / y_paper)
    if y_ratio > (6.0 / 29.0) ** 3:
        l_norm = float(np.clip(1.16 * (y_ratio ** (1.0 / 3.0)) - 0.16, 0.0, 1.25))
    else:
        l_norm = float(np.clip(9.03296296 * y_ratio, 0.0, 1.25))

    sigma_fine = float(np.clip(0.08 * eff_px_per_mm, 0.6, 4.0))
    sigma_coarse = float(max(2.5 * sigma_fine, np.clip(0.35 * eff_px_per_mm, 2.2, 16.0)))
    blur_fine = cv2.GaussianBlur(l_cie, (0, 0), sigmaX=sigma_fine, sigmaY=sigma_fine)
    blur_coarse = cv2.GaussianBlur(l_cie, (0, 0), sigmaX=sigma_coarse, sigmaY=sigma_coarse)
    bp_vals = (blur_fine - blur_coarse)[interior]
    bp_med = float(np.median(bp_vals))
    mad_sd = 1.4826 * float(np.median(np.abs(bp_vals - bp_med))) * (90.0 / max(eff_paper_l, 50.0))

    return {
        "blade_cie_l_norm": round(l_norm, 4),
        "blade_cie_l_sd": round(max(0.0, mad_sd), 4),
        "blade_cie_a_delta": round(med_a - paper_a, 4),
        "blade_cie_b_delta": round(med_b - paper_b, 4),
    }


def harmonize_upright_mask_sinus_bridges(
    upright_mask_u8: np.ndarray,
    pixels_per_mm: Optional[float] = None,
    default_max_bridge_px: float = 3.5,
    hole_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, float]:
    """
    Open physical-mm hairline sinus bridges (`0.30 * px/mm` clamped to `[2, 8] px`)
    on an upright binary lamina mask while measuring residual enclosed voids.
    """
    from scripts.vision.leaf_geometry_utils import open_hairline_sinus_bridges_and_measure_voids

    max_bridge_px = float(default_max_bridge_px)
    if pixels_per_mm is not None:
        try:
            ppm = float(pixels_per_mm)
            if np.isfinite(ppm) and ppm > 0.0:
                max_bridge_px = float(np.clip(0.30 * ppm, 2.0, 8.0))
        except (TypeError, ValueError):
            pass

    cleaned_u8, evf, _ = open_hairline_sinus_bridges_and_measure_voids(
        upright_mask_u8,
        hole_mask=hole_mask,
        max_bridge_px=max_bridge_px,
    )
    return cleaned_u8, float(evf)


def _extract_item_bbox(item: Dict[str, Any]) -> Optional[Tuple[float, float, float, float]]:
    """Extract `(ymin, xmin, ymax, xmax)` floats from an item dictionary."""
    bb = item.get("bbox")
    if bb is None or len(bb) < 4:
        return None
    return (float(bb[0]), float(bb[1]), float(bb[2]), float(bb[3]))


def estimate_corymb_cluster_head_count(
    cluster_items: List[Dict[str, Any]],
    single_head_items: List[Dict[str, Any]],
    pixels_per_mm: Optional[float] = None,
    sheet_shape_hw: Tuple[int, int] = (0, 0),
) -> Tuple[int, int]:
    """
    Estimate additional capitula inside `Seed_Fruit_MANY` (3) and `Flower_MANY` (5)
    corymb clusters after subtracting overlapping single heads (`2`, `4`, `6`).
    """
    clusters = [
        it for it in (cluster_items or [])
        if isinstance(it, dict) and it.get("class_id") is not None and int(it["class_id"]) in (3, 5)
    ]
    if not clusters:
        return (0, 0)

    singles = [
        it for it in (single_head_items or [])
        if isinstance(it, dict) and it.get("class_id") is not None and int(it["class_id"]) in (2, 4, 6)
    ]
    eff_px_per_mm, _ = resolve_effective_px_per_mm(pixels_per_mm, sheet_shape_hw)

    single_areas: List[float] = []
    single_centers: List[Tuple[float, float]] = []
    for item in singles:
        s_mask = item.get("mask")
        s_bb = _extract_item_bbox(item)
        area_val = 0.0
        if isinstance(s_mask, np.ndarray) and s_mask.size > 0:
            nz = int(np.count_nonzero(s_mask > 0))
            if nz > 0:
                area_val = float(nz)
        if area_val <= 0.0 and s_bb is not None:
            bb_area = max(0.0, s_bb[2] - s_bb[0]) * max(0.0, s_bb[3] - s_bb[1])
            if bb_area > 0.0:
                area_val = 0.55 * bb_area
        if area_val > 0.0:
            single_areas.append(area_val)

        cxy = item.get("center_xy")
        if cxy is not None and len(cxy) >= 2:
            single_centers.append((float(cxy[0]), float(cxy[1])))
        elif s_bb is not None:
            single_centers.append((0.5 * (s_bb[1] + s_bb[3]), 0.5 * (s_bb[0] + s_bb[2])))

    ref_head_area = float(np.median(single_areas)) if single_areas else (28.0 * (eff_px_per_mm ** 2))

    extra_heads = 0
    for cl in clusters:
        cl_bb = _extract_item_bbox(cl)
        if cl_bb is not None:
            ymin, xmin, ymax, xmax = cl_bb
            bbox_area = max(0.0, ymax - ymin) * max(0.0, xmax - xmin)
            n_overlap = sum(1 for cx, cy in single_centers if xmin <= cx <= xmax and ymin <= cy <= ymax)
        else:
            bbox_area, n_overlap = 0.0, 0

        cl_mask, cl_crop = cl.get("mask"), cl.get("crop_bgr")
        fg_area: Optional[float] = None
        if isinstance(cl_mask, np.ndarray) and cl_mask.size > 0:
            nz = int(np.count_nonzero(cl_mask > 0))
            if nz > 0:
                fg_area = float(nz)
        if fg_area is None and isinstance(cl_crop, np.ndarray) and cl_crop.ndim == 3 and cl_crop.shape[2] >= 3 and cl_crop.size > 0:
            crop3 = cl_crop[:, :, :3]
            l_star = cv2.cvtColor(crop3, cv2.COLOR_BGR2LAB)[:, :, 0].astype(np.float32) * (100.0 / 255.0)
            s_norm = cv2.cvtColor(crop3, cv2.COLOR_BGR2HSV)[:, :, 1].astype(np.float32) / 255.0
            fg_area = float(np.count_nonzero((l_star < 78.0) | (s_norm > 0.16)))
        if fg_area is None:
            fg_area = 0.55 * bbox_area

        est_n = int(np.clip(round(fg_area / max(ref_head_area, 16.0)), 2, 35))
        extra_heads += max(0, est_n - n_overlap)

    return (int(extra_heads), int(len(clusters)))


__all__ = [
    "resolve_effective_px_per_mm",
    "estimate_sheet_paper_cie_l_ref_and_ab",
    "estimate_sheet_paper_cie_l_ref",
    "compute_calibrated_blade_photometry",
    "harmonize_upright_mask_sinus_bridges",
    "estimate_corymb_cluster_head_count",
]
