#!/usr/bin/env python3
"""
===============================================================================
Module: capitulum_phenology_classifier.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Vectorized "Crop-and-Classify" visual biomarker module for Asteraceae capitula
    (Class 6 detections). Resolves phenological lifecycle transitions across three
    distinct developmental regimes:
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
from typing import Any, Dict, List, Optional, Tuple, Union

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


def classify_single_capitulum(
    crop_rgb: np.ndarray,
    min_crop_dim: int = 5,
) -> Dict[str, Any]:
    """
    Performs head-level phenological biomarker classification on a single cropped capitulum.

    Decision Logic:
      1. If has_pappus_plume and not has_fresh_yellow_corollas -> HEAD_FRUIT.
      2. If has_fresh_yellow_corollas and erect_involucre -> HEAD_ANTHESIS.
      3. If erect_involucre and low yellow and low white -> HEAD_BUD.
      4. Fallback arbitration:
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

    Returns
    -------
    Dict[str, Any]
        Dictionary with classification state, booleans, and extracted biomarkers.
    """
    if (
        crop_rgb is None
        or crop_rgb.size == 0
        or crop_rgb.shape[0] < min_crop_dim
        or crop_rgb.shape[1] < min_crop_dim
    ):
        return {
            "head_state": HEAD_UNKNOWN,
            "has_pappus_plume": False,
            "has_fresh_yellow_corollas": False,
            "erect_involucre": False,
            "reflexed_involucre": False,
            "white_fibrous_ratio": 0.0,
            "yellow_ratio": 0.0,
            "involucre_aspect_ratio": 0.0,
            "boundary_roughness": 1.0,
            "diagnostics": {},
        }

    # Extract Biomarker A: Pappus Plume
    has_pappus_plume, white_fibrous_ratio, diag_a = extract_biomarker_pappus_plume(crop_rgb)

    # Extract Biomarker B: Active Yellow Corollas
    has_fresh_yellow_corollas, yellow_ratio, diag_b = extract_biomarker_yellow_corollas(crop_rgb)

    # Extract Biomarker C: Involucre Posture
    erect_involucre, reflexed_involucre, diag_c = extract_biomarker_involucre_posture(crop_rgb)

    # Head Decision Logic
    if has_pappus_plume and (not has_fresh_yellow_corollas):
        head_state = HEAD_FRUIT
    elif has_fresh_yellow_corollas and erect_involucre:
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

    return {
        "head_state": head_state,
        "has_pappus_plume": has_pappus_plume,
        "has_fresh_yellow_corollas": has_fresh_yellow_corollas,
        "erect_involucre": erect_involucre,
        "reflexed_involucre": reflexed_involucre,
        "white_fibrous_ratio": white_fibrous_ratio,
        "yellow_ratio": yellow_ratio,
        "involucre_aspect_ratio": diag_c.get("involucre_aspect_ratio", 1.0),
        "boundary_roughness": diag_c.get("boundary_roughness", 1.0),
        "diagnostics": {
            "pappus": diag_a,
            "corolla": diag_b,
            "involucre": diag_c,
        },
    }


def classify_voucher_phenology(
    catalog_number: str,
    head_classifications: List[Dict[str, Any]],
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
    n_fruit = sum(1 for h in head_classifications if h.get("head_state") == HEAD_FRUIT)
    n_bud = sum(1 for h in head_classifications if h.get("head_state") == HEAD_BUD)

    if n_anthesis >= 1:
        voucher_state = VOUCHER_ANTHESIS
    elif n_fruit > 0 and n_anthesis == 0:
        voucher_state = VOUCHER_FRUIT
    elif total_capitula == 0:
        voucher_state = VOUCHER_STERILE
    elif n_bud > 0 and n_anthesis == 0 and n_fruit == 0:
        voucher_state = VOUCHER_BUD
    else:
        voucher_state = VOUCHER_UNKNOWN

    pappus_ratios = [
        float(h.get("white_fibrous_ratio", 0.0))
        for h in head_classifications
        if "white_fibrous_ratio" in h
    ]
    dominant_pappus_ratio = round(max(pappus_ratios), 4) if pappus_ratios else 0.0

    return {
        "catalogNumber": str(catalog_number).strip(),
        "total_capitula": int(total_capitula),
        "n_anthesis": int(n_anthesis),
        "n_fruit": int(n_fruit),
        "n_bud": int(n_bud),
        "voucher_phenological_state": voucher_state,
        "dominant_pappus_ratio": float(dominant_pappus_ratio),
    }
