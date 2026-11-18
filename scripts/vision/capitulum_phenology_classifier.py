#!/usr/bin/env python3
"""
===============================================================================
Module: capitulum_phenology_classifier.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Vectorized "Crop-and-Classify" visual biomarker module for Asteraceae capitula
    (LeafMachine3 reproductive Classes 2..6: Seed_Fruit_ONE, Seed_Fruit_MANY,
    Flower_ONE, Flower_MANY, and Bud). Resolves phenological lifecycle transitions
    across three distinct developmental regimes:
      1. Vegetative Bud (HEAD_BUD): Compact, erect, smooth green involucre;
         absence of expanded pappus plume or exposed yellow ligules.
      2. Active Flowering (HEAD_ANTHESIS): Fresh, saturated yellow ligulate and
         disc corollas exposed; erect to campanulate involucre profile.
      3. Fruit / Achene Dispersal (HEAD_FRUIT): Expanded capillary pappus plume
         (high-luminance, low-saturation fibrous bristles on upper dome) with or
         without reflexed/flaring phyllaries.

    Aggregates head-level states into voucher-level phenological classifications
    governed by the Botanical Determinate Cyme Rule:
      - N_anthesis >= 1  -> voucher_phenological_state = "anthesis"
      - N_fruit > 0 and N_anthesis == 0 -> voucher_phenological_state = "fruit"
      - N_total == 0     -> voucher_phenological_state = "sterile"
      - N_bud > 0 and N_anthesis == 0 and N_fruit == 0 -> voucher_phenological_state = "bud"

Execution Constraints:
    - Pure vectorized NumPy and OpenCV calculations (< 0.05s per crop).
    - Robust guards against empty/occluded crops and zero-division.
===============================================================================
"""

from __future__ import annotations

import logging
import math
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np

logger = logging.getLogger(__name__)

# State Constants
HEAD_ANTHESIS = "anthesis"
HEAD_FRUIT = "fruit"
HEAD_BUD = "bud"
HEAD_UNKNOWN = "unknown"

VOUCHER_ANTHESIS = "anthesis"
VOUCHER_FRUIT = "fruit"
VOUCHER_BUD = "bud"
VOUCHER_STERILE = "sterile"
VOUCHER_UNKNOWN = "unknown"


class CapitulumResult(dict):
    """Specialized dictionary returning head metrics while supporting string equality against head_state."""

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.get("head_state") == other
        return super().__eq__(other)

    def __hash__(self) -> int:
        return hash(self.get("head_state"))


HeadPhenologyResult = CapitulumResult


class VoucherPhenologyResult(dict):
    """Specialized dictionary returning voucher metrics while supporting string equality against voucher_phenological_state."""

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, str):
            return self.get("voucher_phenological_state") == other
        return super().__eq__(other)

    def __hash__(self) -> int:
        return hash(self.get("voucher_phenological_state"))


def extract_biomarker_pappus_plume(
    crop_rgb: np.ndarray,
    upper_pct: float = 0.45,
    min_l_star: float = 75.0,
    max_s: float = 0.18,
    min_laplacian_var: float = 40.0,
    min_edge_density: float = 0.03,
    plume_threshold: float = 0.18,
) -> Tuple[bool, float, Dict[str, float]]:
    """
    Biomarker A: Expanded Capillary Pappus Plume (Fruit Indicator).

    Extracts high-luminance, low-saturation pixels in the upper ROI of the capitulum
    dome, confirmed by high-frequency texture (Laplacian variance / Canny edge density)
    to differentiate fine capillary bristles from plain white mounting paper.

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image crop of a detected capitulum.
    upper_pct : float
        Fraction of upper vertical height representing the apical dome (default 0.45).
    min_l_star : float
        Minimum CIELAB lightness (L* >= 75.0) for bright pappus bristles.
    max_s : float
        Maximum HSV saturation (S <= 0.18) for neutral white/cream bristles.
    min_laplacian_var : float
        Threshold for Laplacian variance confirming fibrous bristle micro-texture.
    min_edge_density : float
        Threshold for Canny edge density confirming fine structural fibers.
    plume_threshold : float
        Threshold ratio of white fibrous pixels in upper ROI to assert plume presence (18%).

    Returns
    -------
    Tuple[bool, float, Dict[str, float]]
        - has_pappus_plume (bool)
        - white_fibrous_ratio (float, 0.0 to 1.0)
        - diagnostics (dict containing intermediate ratios and texture metrics)
    """
    h, w = crop_rgb.shape[:2]
    if h < 3 or w < 3:
        return False, 0.0, {
            "white_ratio": 0.0,
            "white_fibrous_ratio": 0.0,
            "laplacian_var": 0.0,
            "edge_density": 0.0,
        }

    upper_h = max(1, int(round(h * upper_pct)))
    upper_roi = crop_rgb[:upper_h, :]
    n_upper = max(1, upper_roi.shape[0] * upper_roi.shape[1])

    # Convert upper dome to LAB and HSV
    upper_lab = cv2.cvtColor(upper_roi, cv2.COLOR_RGB2LAB)
    upper_hsv = cv2.cvtColor(upper_roi, cv2.COLOR_RGB2HSV)
    upper_gray = cv2.cvtColor(upper_roi, cv2.COLOR_RGB2GRAY)

    # L* in OpenCV uint8 is scaled: L* = L_cv2 * (100.0 / 255.0)
    l_star = upper_lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    s_norm = upper_hsv[:, :, 1].astype(np.float32) / 255.0

    white_mask = (l_star >= min_l_star) & (s_norm <= max_s)
    white_pixel_ratio = float(np.count_nonzero(white_mask) / n_upper)

    # High-frequency texture analysis
    laplacian = cv2.Laplacian(upper_gray, cv2.CV_64F)
    laplacian_var = float(laplacian.var())

    edges = cv2.Canny(upper_gray, 50, 150)
    edge_count = np.count_nonzero(edges)
    edge_density = float(edge_count / n_upper)

    has_texture = (laplacian_var >= min_laplacian_var) or (edge_density >= min_edge_density)

    # Differentiate textured white bristles from uniform white mounting paper:
    # Dilate edges slightly (3x3 kernel) to encompass fibrous bristle shafts and gradients
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
    dilated_edges = cv2.dilate(edges, kernel) > 0
    fibrous_white_mask = white_mask & (dilated_edges | (np.abs(laplacian) > 15.0))
    fibrous_white_count = np.count_nonzero(fibrous_white_mask)
    white_fibrous_ratio = float(fibrous_white_count / n_upper)

    has_pappus_plume = bool((white_fibrous_ratio > plume_threshold) and has_texture)

    diagnostics = {
        "white_ratio": round(white_pixel_ratio, 4),
        "white_fibrous_ratio": round(white_fibrous_ratio, 4),
        "laplacian_var": round(laplacian_var, 2),
        "edge_density": round(edge_density, 4),
    }

    return has_pappus_plume, float(round(white_fibrous_ratio, 4)), diagnostics


def extract_biomarker_yellow_corollas(
    crop_rgb: np.ndarray,
    h_min_deg: float = 18.0,
    h_max_deg: float = 45.0,
    s_min: float = 0.32,
    v_min: float = 0.35,
    corolla_threshold: float = 0.12,
) -> Tuple[bool, float, Dict[str, float]]:
    """
    Biomarker B: Active Yellow Corollas / Ligules (Anthesis Indicator).

    Isolates the bright golden-yellow hue band across the capitulum disc and ray margins:
      18° <= Hue <= 45°, Saturation >= 0.32, Value >= 0.35.
    Saturated yellow area > 12% indicates fresh, active anthesis.

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image crop of a detected capitulum.
    h_min_deg : float
        Minimum hue in degrees (default 18.0°).
    h_max_deg : float
        Maximum hue in degrees (default 45.0°).
    s_min : float
        Minimum saturation (default 0.32).
    v_min : float
        Minimum brightness value (default 0.35).
    corolla_threshold : float
        Proportion of saturated yellow pixels to confirm fresh corollas (default 0.12).

    Returns
    -------
    Tuple[bool, float, Dict[str, float]]
        - has_fresh_yellow_corollas (bool)
        - yellow_ratio (float, 0.0 to 1.0)
        - diagnostics (dict containing yellow ratio and thresholds)
    """
    h, w = crop_rgb.shape[:2]
    n_total = max(1, h * w)

    if h < 3 or w < 3:
        return False, 0.0, {"yellow_ratio": 0.0}

    crop_hsv = cv2.cvtColor(crop_rgb, cv2.COLOR_RGB2HSV)

    # In OpenCV uint8, H is in [0, 179], representing 0° to 360° (H_deg = H_cv2 * 2.0)
    h_deg = crop_hsv[:, :, 0].astype(np.float32) * 2.0
    s_norm = crop_hsv[:, :, 1].astype(np.float32) / 255.0
    v_norm = crop_hsv[:, :, 2].astype(np.float32) / 255.0

    yellow_mask = (
        (h_deg >= h_min_deg)
        & (h_deg <= h_max_deg)
        & (s_norm >= s_min)
        & (v_norm >= v_min)
    )

    yellow_pixels = np.count_nonzero(yellow_mask)
    yellow_ratio = float(yellow_pixels / n_total)
    has_fresh_yellow_corollas = bool(yellow_ratio > corolla_threshold)

    diagnostics = {
        "yellow_ratio": round(yellow_ratio, 4),
        "yellow_pixel_count": int(yellow_pixels),
    }

    return has_fresh_yellow_corollas, float(round(yellow_ratio, 4)), diagnostics


def extract_biomarker_involucre_posture(
    crop_rgb: np.ndarray,
    lower_pct: float = 0.45,
    min_aspect_ratio: float = 1.0,
    max_reflexed_ar: float = 0.75,
    roughness_thresh: float = 1.15,
) -> Tuple[bool, bool, Dict[str, float]]:
    """
    Biomarker C: Involucre Posture (Lower Contour Profile).

    Segments the herbaceous lower involucre from background mounting paper and
    characterizes the basal profile:
      - Erect involucre: Smooth, tapered U-shape (H_inv / W_inv >= 1.0, low roughness).
      - Reflexed involucre: Flattened or flared downward (H_inv / W_inv <= 0.75, high boundary roughness).

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image crop of a detected capitulum.
    lower_pct : float
        Vertical split point starting the lower involucre (default 0.45).
    min_aspect_ratio : float
        Minimum H_inv / W_inv for an erect cylindrical/campanulate involucre (1.0).
    max_reflexed_ar : float
        Maximum H_inv / W_inv for reflexed/flaring phyllaries (0.75).
    roughness_thresh : float
        Contour perimeter / convex hull perimeter threshold for jagged reflexed phyllaries (1.15).

    Returns
    -------
    Tuple[bool, bool, Dict[str, float]]
        - erect_involucre (bool)
        - reflexed_involucre (bool)
        - diagnostics (dict with H_inv, W_inv, aspect_ratio, boundary_roughness)
    """
    h, w = crop_rgb.shape[:2]
    if h < 4 or w < 4:
        return True, False, {
            "h_inv_px": 0.0,
            "w_inv_px": 0.0,
            "involucre_aspect_ratio": 1.0,
            "boundary_roughness": 1.0,
        }

    lower_y = int(round(h * lower_pct))
    lower_roi = crop_rgb[lower_y:, :]
    lh, lw = lower_roi.shape[:2]
    if lh < 2 or lw < 2:
        return True, False, {
            "h_inv_px": float(lh),
            "w_inv_px": float(lw),
            "involucre_aspect_ratio": 1.0,
            "boundary_roughness": 1.0,
        }

    lower_lab = cv2.cvtColor(lower_roi, cv2.COLOR_RGB2LAB)
    lower_hsv = cv2.cvtColor(lower_roi, cv2.COLOR_RGB2HSV)

    l_star = lower_lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    s_norm = lower_hsv[:, :, 1].astype(np.float32) / 255.0
    h_deg = lower_hsv[:, :, 0].astype(np.float32) * 2.0

    # Plant tissue segmentation: filter out high-luminance neutral mounting paper
    is_paper = (l_star >= 75.0) & (s_norm <= 0.20)
    is_plant = ~is_paper

    # Clean morphological structure
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    plant_mask = cv2.morphologyEx(is_plant.astype(np.uint8) * 255, cv2.MORPH_CLOSE, kernel)

    contours, _ = cv2.findContours(plant_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    if contours and cv2.contourArea(max(contours, key=cv2.contourArea)) >= 10:
        main_c = max(contours, key=cv2.contourArea)
        bx, by, bw, bh = cv2.boundingRect(main_c)
        h_inv = float(bh)
        w_inv = float(max(1, bw))
        ar_inv = h_inv / w_inv

        peri = float(cv2.arcLength(main_c, True))
        hull = cv2.convexHull(main_c)
        hull_peri = float(cv2.arcLength(hull, True))
        roughness = float(peri / max(1e-6, hull_peri))

        erect_involucre = bool(ar_inv >= min_aspect_ratio or (ar_inv > max_reflexed_ar and roughness < roughness_thresh))
        reflexed_involucre = bool(ar_inv <= max_reflexed_ar and roughness >= roughness_thresh)
    else:
        # Fallback to lower ROI bounding aspect ratio
        h_inv = float(lh)
        w_inv = float(max(1, lw))
        ar_inv = h_inv / w_inv
        roughness = 1.0
        erect_involucre = bool(ar_inv >= 0.85)
        reflexed_involucre = bool(ar_inv <= 0.65)

    diagnostics = {
        "h_inv_px": round(h_inv, 2),
        "w_inv_px": round(w_inv, 2),
        "involucre_aspect_ratio": round(ar_inv, 4),
        "boundary_roughness": round(roughness, 4),
    }

    return erect_involucre, reflexed_involucre, diagnostics


def orient_capitulum_crop_upright(
    crop_rgb: np.ndarray,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
) -> np.ndarray:
    """
    Reorients a capitulum RGB crop so its basal-to-distal axis points upward.

    Uses `peduncle_dir_xy` when supplied and falls back to crop biomarker centroids
    between distal florets or pappus bristles and the herbaceous involucre cup.
    """
    if (
        crop_rgb is None
        or not isinstance(crop_rgb, np.ndarray)
        or crop_rgb.size == 0
        or crop_rgb.ndim < 3
        or crop_rgb.shape[2] < 3
        or min(crop_rgb.shape[:2]) < 5
    ):
        return crop_rgb

    h, w = crop_rgb.shape[:2]
    dx: Optional[float] = None
    dy: Optional[float] = None

    if peduncle_dir_xy is not None and len(peduncle_dir_xy) >= 2:
        try:
            px = float(peduncle_dir_xy[0])
            py = float(peduncle_dir_xy[1])
        except (TypeError, ValueError):
            px, py = float("nan"), float("nan")
        if np.isfinite(px) and np.isfinite(py):
            nrm = math.hypot(px, py)
            if nrm >= 1e-3:
                dx = px / nrm
                dy = py / nrm

    if dx is None or dy is None:
        rgb3 = crop_rgb[:, :, :3]
        lab = cv2.cvtColor(rgb3, cv2.COLOR_RGB2LAB)
        hsv = cv2.cvtColor(rgb3, cv2.COLOR_RGB2HSV)
        gray = cv2.cvtColor(rgb3, cv2.COLOR_RGB2GRAY)

        l_star = lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
        b_star = lab[:, :, 2].astype(np.float32) - 128.0
        h_deg = hsv[:, :, 0].astype(np.float32) * 2.0
        s_norm = hsv[:, :, 1].astype(np.float32) / 255.0
        v_norm = hsv[:, :, 2].astype(np.float32) / 255.0

        non_paper = ~((l_star >= 75.0) & (s_norm <= 0.20))
        floret_m = (
            (h_deg >= 10.0)
            & (h_deg <= 50.0)
            & (s_norm >= 0.18)
            & (v_norm >= 0.28)
            & (b_star >= 8.0)
        )
        edges = cv2.Canny(gray, 50, 150)
        dil_edges = cv2.dilate(edges, cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))) > 0
        pappus_m = (l_star >= 75.0) & (s_norm <= 0.18) & dil_edges
        distal_m = floret_m | pappus_m
        cup_m = non_paper & (~floret_m)

        if np.count_nonzero(distal_m) >= 8 and np.count_nonzero(cup_m) >= 8:
            ys_d, xs_d = np.where(distal_m)
            ys_c, xs_c = np.where(cup_m)
            cy_d, cx_d = float(np.mean(ys_d)), float(np.mean(xs_d))
            cy_c, cx_c = float(np.mean(ys_c)), float(np.mean(xs_c))
            vx = cx_d - cx_c
            vy = cy_d - cy_c
            if math.hypot(vx, vy) >= 3.0:
                dx = vx
                dy = vy

    if dx is None or dy is None:
        return crop_rgb

    theta_vec = math.degrees(math.atan2(dy, dx))
    rot_deg = ((theta_vec + 90.0 + 180.0) % 360.0) - 180.0
    if abs(rot_deg) <= 25.0:
        return crop_rgb

    rad = math.radians(rot_deg)
    abs_cos = abs(math.cos(rad))
    abs_sin = abs(math.sin(rad))
    new_w = max(1, int(round(h * abs_sin + w * abs_cos)))
    new_h = max(1, int(round(h * abs_cos + w * abs_sin)))
    rot_mat = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), rot_deg, 1.0)
    rot_mat[0, 2] += (new_w / 2.0) - (w / 2.0)
    rot_mat[1, 2] += (new_h / 2.0) - (h / 2.0)
    return cv2.warpAffine(
        crop_rgb,
        rot_mat,
        (new_w, new_h),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=(240, 240, 240),
    )


def extract_biomarker_aged_ochre_corollas(
    crop_rgb: np.ndarray,
    h_min_deg: float = 10.0,
    h_max_deg: float = 48.0,
    s_min: float = 0.18,
    s_max: float = 0.32,
    v_min: float = 0.28,
    min_b_star_excess: float = 8.0,
    ochre_threshold: float = 0.10,
) -> Tuple[bool, float, Dict[str, float]]:
    """
    Biomarker B2: Aged Ochre and Amber Corollas (Herbarium Anthesis Indicator).

    Detects oxidized ochre, amber, and browned ligules or disc corollas in the
    upper capitulum dome by thresholding HSV hue/saturation/value and CIE b*
    excess over border mounting paper while excluding bright achromatic pappus plumes.
    """
    if (
        crop_rgb is None
        or not isinstance(crop_rgb, np.ndarray)
        or crop_rgb.size == 0
        or crop_rgb.ndim < 3
        or crop_rgb.shape[2] < 3
    ):
        return False, 0.0, {
            "ochre_ratio": 0.0,
            "ochre_pixel_count": 0,
            "paper_b_star": 0.0,
        }

    h, w = crop_rgb.shape[:2]
    if h < 3 or w < 3:
        return False, 0.0, {
            "ochre_ratio": 0.0,
            "ochre_pixel_count": 0,
            "paper_b_star": 0.0,
        }

    rgb3 = crop_rgb[:, :, :3]
    crop_lab = cv2.cvtColor(rgb3, cv2.COLOR_RGB2LAB)
    crop_hsv = cv2.cvtColor(rgb3, cv2.COLOR_RGB2HSV)

    l_full = crop_lab[:, :, 0].astype(np.float32) * (100.0 / 255.0)
    b_full = crop_lab[:, :, 2].astype(np.float32) - 128.0
    h_full = crop_hsv[:, :, 0].astype(np.float32) * 2.0
    s_full = crop_hsv[:, :, 1].astype(np.float32) / 255.0
    v_full = crop_hsv[:, :, 2].astype(np.float32) / 255.0

    # Estimate border mounting paper CIE b* baseline from outermost 3 px border
    border_mask = np.zeros((h, w), dtype=bool)
    bw = min(3, max(1, h // 4), max(1, w // 4))
    border_mask[:bw, :] = True
    border_mask[-bw:, :] = True
    border_mask[:, :bw] = True
    border_mask[:, -bw:] = True
    paper_border = border_mask & (l_full >= 65.0) & (s_full <= 0.22)
    b_paper = float(np.median(b_full[paper_border])) if np.any(paper_border) else 0.0

    upper_h = max(1, int(round(h * 0.58)))
    l_star = l_full[:upper_h, :]
    b_star = b_full[:upper_h, :]
    h_deg = h_full[:upper_h, :]
    s_norm = s_full[:upper_h, :]
    v_norm = v_full[:upper_h, :]

    is_achromatic_pappus = (l_star >= 75.0) & (s_norm <= 0.18)
    ochre_mask = (
        (h_deg >= h_min_deg)
        & (h_deg <= h_max_deg)
        & (s_norm >= s_min)
        & (s_norm < max(s_max, s_min + 1e-4))
        & (v_norm >= v_min)
        & ((b_star - b_paper) >= min_b_star_excess)
        & (~is_achromatic_pappus)
    )

    ochre_pixels = int(np.count_nonzero(ochre_mask))
    n_upper = max(1, upper_h * w)
    ochre_ratio = float(ochre_pixels / n_upper)
    has_aged_ochre = bool(ochre_ratio >= ochre_threshold)

    diagnostics = {
        "ochre_ratio": round(ochre_ratio, 4),
        "ochre_pixel_count": ochre_pixels,
        "paper_b_star": round(b_paper, 2),
    }

    return has_aged_ochre, float(round(ochre_ratio, 4)), diagnostics


def classify_single_capitulum(
    crop_rgb: np.ndarray,
    min_crop_dim: int = 5,
    min_pappus_luminance: float = 75.0,
    max_pappus_saturation: float = 0.18,
    min_pappus_area_ratio: float = 0.18,
    min_floret_yellow_saturation: float = 0.32,
    min_floret_yellow_ratio: float = 0.12,
    class_id: Optional[int] = None,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
) -> CapitulumResult:
    """
    Performs head-level phenological biomarker classification on a single cropped capitulum.

    Decision Logic:
      1. If has_pappus_plume and not has_fresh_yellow_corollas -> HEAD_FRUIT.
      2. If has_fresh_yellow_corollas and erect_involucre -> HEAD_ANTHESIS.
      3. If both pappus and fresh yellow are absent and aged ochre corollas pass
         class/posture gating -> HEAD_ANTHESIS.
      4. If erect_involucre and low yellow and low white -> HEAD_BUD.
      5. Fallback arbitration:
         - has_pappus_plume -> HEAD_FRUIT
         - has_fresh_yellow_corollas -> HEAD_ANTHESIS
         - reflexed_involucre -> HEAD_FRUIT
         - else -> HEAD_BUD

    Parameters
    ----------
    crop_rgb : np.ndarray
        (H, W, 3) RGB image patch of a single detected capitulum.
    min_crop_dim : int
        Minimum height and width required to process crop (default 5 px).
    class_id : Optional[int]
        Optional LeafMachine3 reproductive class identifier (2..6).
    peduncle_dir_xy : Optional[Tuple[float, float]]
        Optional basal-to-distal peduncle direction vector in image (x, y) coordinates.

    Returns
    -------
    CapitulumResult
        Dictionary with classification state, booleans, and extracted biomarkers.
    """
    if (
        crop_rgb is None
        or crop_rgb.size == 0
        or crop_rgb.ndim < 3
        or crop_rgb.shape[2] < 3
        or crop_rgb.shape[0] < min_crop_dim
        or crop_rgb.shape[1] < min_crop_dim
    ):
        return CapitulumResult({
            "head_state": HEAD_UNKNOWN,
            "has_pappus_plume": False,
            "has_fresh_yellow_corollas": False,
            "has_aged_ochre_corollas": False,
            "erect_involucre": False,
            "reflexed_involucre": False,
            "white_fibrous_ratio": 0.0,
            "yellow_ratio": 0.0,
            "ochre_ratio": 0.0,
            "involucre_aspect_ratio": 0.0,
            "boundary_roughness": 1.0,
            "diagnostics": {},
        })

    oriented_rgb = orient_capitulum_crop_upright(crop_rgb, peduncle_dir_xy=peduncle_dir_xy)

    # Extract Biomarker A: Pappus Plume
    has_pappus_plume, white_fibrous_ratio, diag_a = extract_biomarker_pappus_plume(
        oriented_rgb,
        min_l_star=float(min_pappus_luminance),
        max_s=float(max_pappus_saturation),
        plume_threshold=float(min_pappus_area_ratio),
    )

    # Extract Biomarker B: Active Yellow Corollas
    has_fresh_yellow_corollas, yellow_ratio, diag_b = extract_biomarker_yellow_corollas(
        oriented_rgb,
        s_min=float(min_floret_yellow_saturation),
        corolla_threshold=float(min_floret_yellow_ratio),
    )

    # Extract Biomarker B2: Aged Ochre / Amber Corollas
    has_aged_ochre, ochre_ratio, diag_ochre = extract_biomarker_aged_ochre_corollas(
        oriented_rgb,
        s_max=float(min_floret_yellow_saturation),
    )

    # Extract Biomarker C: Involucre Posture
    erect_involucre, reflexed_involucre, diag_c = extract_biomarker_involucre_posture(oriented_rgb)
    inv_ar = float(diag_c.get("involucre_aspect_ratio", 1.0))

    aged_anthesis_ok = bool(
        (not has_pappus_plume)
        and (not has_fresh_yellow_corollas)
        and has_aged_ochre
        and (class_id != 6)
        and (
            (class_id in (4, 5))
            or (erect_involucre and inv_ar >= 0.95 and (not reflexed_involucre))
        )
    )

    # Head Decision Logic
    if has_pappus_plume and (not has_fresh_yellow_corollas):
        head_state = HEAD_FRUIT
    elif has_fresh_yellow_corollas and erect_involucre:
        head_state = HEAD_ANTHESIS
    elif aged_anthesis_ok:
        head_state = HEAD_ANTHESIS
    elif erect_involucre and (not has_fresh_yellow_corollas) and (not has_pappus_plume):
        head_state = HEAD_BUD
    else:
        # Resolving intermediate/ambiguous boundary cases
        if has_pappus_plume:
            head_state = HEAD_FRUIT
        elif has_fresh_yellow_corollas:
            head_state = HEAD_ANTHESIS
        elif reflexed_involucre:
            head_state = HEAD_FRUIT
        else:
            head_state = HEAD_BUD

    return CapitulumResult({
        "head_state": head_state,
        "has_pappus_plume": has_pappus_plume,
        "has_fresh_yellow_corollas": has_fresh_yellow_corollas,
        "has_aged_ochre_corollas": bool(aged_anthesis_ok),
        "erect_involucre": erect_involucre,
        "reflexed_involucre": reflexed_involucre,
        "white_fibrous_ratio": white_fibrous_ratio,
        "yellow_ratio": yellow_ratio,
        "ochre_ratio": ochre_ratio,
        "involucre_aspect_ratio": diag_c.get("involucre_aspect_ratio", 1.0),
        "boundary_roughness": diag_c.get("boundary_roughness", 1.0),
        "diagnostics": {
            "pappus": diag_a,
            "corolla": diag_b,
            "aged_ochre": diag_ochre,
            "involucre": diag_c,
        },
    })


def classify_voucher_phenology(
    catalog_number: str,
    head_classifications: List[Dict[str, Any]],
    min_anthesis_heads_for_flowering: int = 1,
) -> Dict[str, Any]:
    """
    Aggregates head-level classifications for a voucher into a single phenological state
    governed by the Botanical Determinate Cyme Rule.

    Biological Cyme Rule:
      - N_anthesis >= 1: Scored as "anthesis" (active flowering head present; valid for Delta DOY).
      - N_fruit > 0 and N_anthesis == 0: Scored as "fruit" (entirely post-anthesis/dispersal).
      - N_total == 0: Scored as "sterile" (vegetative rosette only).
      - N_bud > 0 and N_anthesis == 0 and N_fruit == 0: Scored as "bud" (pre-anthesis).

    Parameters
    ----------
    catalog_number : str
        Herbarium voucher accession / catalog identifier.
    head_classifications : List[Dict[str, Any]]
        List of head-level classification dictionaries from `classify_single_capitulum`.

    Returns
    -------
    Dict[str, Any]
        Voucher phenology record formatted for manifest export.
    """
    total_capitula = len(head_classifications)
    n_anthesis = sum(1 for h in head_classifications if h.get("head_state") == HEAD_ANTHESIS)
    n_fruiting = sum(1 for h in head_classifications if h.get("head_state") == HEAD_FRUIT)
    n_bud = sum(1 for h in head_classifications if h.get("head_state") == HEAD_BUD)

    if n_anthesis >= max(1, int(min_anthesis_heads_for_flowering)):
        voucher_state = VOUCHER_ANTHESIS
    elif n_fruiting > 0 and n_anthesis == 0:
        voucher_state = VOUCHER_FRUIT
    elif total_capitula == 0:
        voucher_state = VOUCHER_STERILE
    elif n_bud > 0 and n_anthesis == 0 and n_fruiting == 0:
        voucher_state = VOUCHER_BUD
    else:
        voucher_state = VOUCHER_UNKNOWN

    pappus_ratios = [
        float(h.get("white_fibrous_ratio", 0.0))
        for h in head_classifications
        if "white_fibrous_ratio" in h
    ]
    dominant_pappus_ratio = round(max(pappus_ratios), 4) if pappus_ratios else 0.0

    return VoucherPhenologyResult({
        "catalogNumber": str(catalog_number).strip(),
        "total_capitula": int(total_capitula),
        "n_anthesis": int(n_anthesis),
        "n_fruiting": int(n_fruiting),
        "n_bud": int(n_bud),
        "voucher_phenological_state": voucher_state,
        "dominant_pappus_ratio": float(dominant_pappus_ratio),
    })


def _extract_item_crop_rgb(
    item: Any,
    sheet_rgb: Optional[np.ndarray],
    bbox: Optional[Tuple[int, int, int, int]],
) -> Optional[np.ndarray]:
    """Resolve an RGB crop array from an item dict/object/array or by slicing `sheet_rgb`."""
    if isinstance(item, np.ndarray):
        return item
    crop_rgb = item.get("crop_rgb") if isinstance(item, dict) else getattr(item, "crop_rgb", None)
    if crop_rgb is None:
        crop_rgb = item.get("crop") if isinstance(item, dict) else getattr(item, "crop", None)
    if crop_rgb is None:
        crop_rgb = item.get("image") if isinstance(item, dict) else getattr(item, "image", None)
    if isinstance(crop_rgb, np.ndarray):
        return crop_rgb
    crop_bgr = item.get("crop_bgr") if isinstance(item, dict) else getattr(item, "crop_bgr", None)
    if isinstance(crop_bgr, np.ndarray) and crop_bgr.size > 0 and crop_bgr.ndim >= 3 and crop_bgr.shape[2] >= 3:
        return cv2.cvtColor(crop_bgr[:, :, :3], cv2.COLOR_BGR2RGB)
    if (
        isinstance(sheet_rgb, np.ndarray)
        and sheet_rgb.ndim >= 3
        and bbox is not None
        and len(bbox) >= 4
    ):
        sh, sw = sheet_rgb.shape[:2]
        y0 = max(0, min(sh, int(round(float(bbox[0])))))
        x0 = max(0, min(sw, int(round(float(bbox[1])))))
        y1 = max(0, min(sh, int(round(float(bbox[2])))))
        x1 = max(0, min(sw, int(round(float(bbox[3])))))
        if y1 > y0 and x1 > x0:
            return sheet_rgb[y0:y1, x0:x1, :3]
    return None


def classify_asteraceae_voucher_crops(
    items: List[Any],
    sheet_rgb: Optional[np.ndarray] = None,
    catalog_number: str = "",
    min_anthesis_heads_for_flowering: int = 1,
    capitulum_kwargs: Optional[Dict[str, Any]] = None,
) -> VoucherPhenologyResult:
    """
    Delegate helper classifying Asteraceae capitulum crops/boxes into voucher phenology
    using `classify_single_capitulum` and `classify_voucher_phenology`.
    """
    cap_kw = dict(capitulum_kwargs or {})
    head_results: List[Dict[str, Any]] = []
    confs: List[float] = []

    for item in items:
        if isinstance(item, dict):
            hp = item.get("head_pheno")
            if isinstance(hp, dict) and "head_state" in hp:
                head_results.append(hp)
                sc = item.get("confidence", item.get("score"))
                if sc is not None and np.isfinite(float(sc)):
                    confs.append(float(np.clip(float(sc), 0.0, 1.0)))
                continue
            if "head_state" in item:
                head_results.append(item)
                sc = item.get("confidence", item.get("score"))
                if sc is not None and np.isfinite(float(sc)):
                    confs.append(float(np.clip(float(sc), 0.0, 1.0)))
                continue

        bbox: Optional[Tuple[int, int, int, int]] = None
        cls_id: Optional[int] = None
        ped_dir: Optional[Tuple[float, float]] = None
        score_val: Optional[float] = None
        lbl_str: Optional[str] = None

        if isinstance(item, dict):
            bbox = item.get("bbox")
            raw_cid = item.get("class_id", item.get("cls_id", item.get("cls")))
            if raw_cid is not None:
                try:
                    cls_id = int(raw_cid)
                except (TypeError, ValueError):
                    cls_id = None
            ped_dir = item.get("peduncle_dir_xy")
            score_val = item.get("confidence", item.get("score", item.get("conf")))
            lbl_str = item.get("label", item.get("class_name", item.get("name")))
        elif isinstance(item, (tuple, list)) and len(item) >= 2:
            if isinstance(item[0], np.ndarray):
                if isinstance(item[1], (int, np.integer)):
                    cls_id = int(item[1])
                elif isinstance(item[1], str):
                    lbl_str = item[1]
            elif len(item) == 4 and all(isinstance(v, (int, float, np.number)) for v in item):
                bbox = (int(item[0]), int(item[1]), int(item[2]), int(item[3]))
            elif len(item) == 3 and isinstance(item[0], (tuple, list, np.ndarray)):
                bbox = tuple(int(v) for v in item[0][:4])  # type: ignore[assignment]
                score_val = float(item[1]) if isinstance(item[1], (int, float, np.number)) else None
                if isinstance(item[2], (int, np.integer)):
                    cls_id = int(item[2])
                elif isinstance(item[2], str):
                    lbl_str = item[2]
            elif len(item) == 4 and isinstance(item[0], (tuple, list, np.ndarray)):
                bbox = tuple(int(v) for v in item[0][:4])  # type: ignore[assignment]
                score_val = float(item[2]) if isinstance(item[2], (int, float, np.number)) else None
                if isinstance(item[3], (int, np.integer)):
                    cls_id = int(item[3])
                elif isinstance(item[3], str):
                    lbl_str = item[3]
        elif not isinstance(item, np.ndarray):
            bbox = getattr(item, "bbox", None)
            raw_cid = getattr(item, "class_id", getattr(item, "cls_id", None))
            if raw_cid is not None:
                try:
                    cls_id = int(raw_cid)
                except (TypeError, ValueError):
                    cls_id = None
            ped_dir = getattr(item, "peduncle_dir_xy", None)
            score_val = getattr(item, "confidence", getattr(item, "score", None))
            lbl_str = getattr(item, "label", getattr(item, "class_name", None))

        if cls_id is not None and cls_id not in (2, 3, 4, 5, 6):
            continue

        crop_rgb = _extract_item_crop_rgb(item, sheet_rgb=sheet_rgb, bbox=bbox)
        if crop_rgb is not None:
            head_res = classify_single_capitulum(
                crop_rgb,
                class_id=cls_id,
                peduncle_dir_xy=ped_dir,
                **cap_kw,
            )
            head_results.append(head_res)
        elif cls_id in (2, 3, 4, 5, 6) or (lbl_str and str(lbl_str).strip()):
            low_l = str(lbl_str or "").strip().lower()
            if cls_id == 6 or "bud" in low_l:
                fb_state = HEAD_BUD
            elif cls_id in (4, 5) or "flower" in low_l or "anthesis" in low_l:
                fb_state = HEAD_ANTHESIS
            elif cls_id in (2, 3) or "fruit" in low_l or "seed" in low_l:
                fb_state = HEAD_FRUIT
            else:
                continue
            head_results.append(CapitulumResult({"head_state": fb_state, "white_fibrous_ratio": 0.0}))
        else:
            continue

        if score_val is not None and np.isfinite(float(score_val)):
            confs.append(float(np.clip(float(score_val), 0.0, 1.0)))

    base = classify_voucher_phenology(
        catalog_number,
        head_results,
        min_anthesis_heads_for_flowering=min_anthesis_heads_for_flowering,
    )
    out = VoucherPhenologyResult(dict(base))
    out["total_reproductive_organs"] = int(base["total_capitula"])
    out["n_senescent"] = 0
    out["is_verified_anthesis"] = bool(
        int(base["n_anthesis"]) >= max(1, int(min_anthesis_heads_for_flowering))
    )
    out["reproductive_stage_confidence"] = (
        round(float(np.mean(confs)), 4)
        if confs
        else (1.0 if int(base["total_capitula"]) > 0 else 0.0)
    )
    out["phenology_plugin"] = "asteraceae_radiate"
    return out
