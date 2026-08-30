"""
scripts/core/gatekeeper_metrics.py
==================================
Deterministic metric computation functions for botanical artifact gatekeeper.
Calculates geometric morphology, spectral saturation, and typographic edge metrics.
"""

from __future__ import annotations

import math
from typing import List, Tuple, Optional, Any, Dict

import cv2
import numpy as np

from scripts.core.data_structures import (
    GeometricMetrics,
    SpectralMetrics,
    TextureMetrics,
)


def compute_geometric_metrics(
    candidate_mask: np.ndarray,
    max_rectangularity_threshold: float = 0.86,
    min_solidity_threshold: float = 0.72,
    dp_epsilon_ratio: float = 0.02,
    orthogonal_angle_range: Tuple[float, float] = (80.0, 100.0)
) -> GeometricMetrics:
    """
    Compute deterministic geometric morphology metrics from a binary leaf candidate mask.
    
    Calculates:
      1. Rectangularity: Ratio of mask area to minimum bounding oriented rectangle area.
      2. Douglas-Peucker polygon approximation and 4-corner internal angle analysis.
      3. Solidity: Ratio of mask area to convex hull area.
      
    Args:
        candidate_mask: Binary single-channel uint8 mask (foreground=255/1, background=0).
        max_rectangularity_threshold: Threshold for rectangular rejection (0.86).
        min_solidity_threshold: Minimum allowable solidity (0.72).
        dp_epsilon_ratio: Epsilon factor for approxPolyDP (0.02).
        orthogonal_angle_range: Internal angle tolerance in degrees (80.0, 100.0).
        
    Returns:
        GeometricMetrics dataclass containing all morphological parameters.
    """
    bin_mask = (candidate_mask > 0).astype(np.uint8) * 255
    mask_area = float(cv2.countNonZero(bin_mask))

    if mask_area < 1.0:
        return GeometricMetrics(
            mask_area=0.0,
            min_area_rect_area=0.0,
            rectangularity=0.0,
            convex_hull_area=0.0,
            solidity=0.0,
            num_approx_vertices=0,
            corner_angles_deg=[],
            is_rectangular=False,
            is_orthogonal_quad=False,
            is_valid_solidity=False
        )

    contours, _ = cv2.findContours(bin_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return GeometricMetrics(
            mask_area=mask_area,
            min_area_rect_area=0.0,
            rectangularity=0.0,
            convex_hull_area=0.0,
            solidity=0.0,
            num_approx_vertices=0,
            corner_angles_deg=[],
            is_rectangular=False,
            is_orthogonal_quad=False,
            is_valid_solidity=False
        )

    cnt = max(contours, key=cv2.contourArea)
    cnt_area = float(cv2.contourArea(cnt))
    if cnt_area < 1.0:
        cnt_area = mask_area

    # 1. Minimum Bounding Oriented Rectangle & Rectangularity
    min_rect = cv2.minAreaRect(cnt)
    rect_w, rect_h = min_rect[1]
    min_area_rect_area = float(rect_w * rect_h)

    if min_area_rect_area > 0.0:
        rectangularity = float(mask_area / min_area_rect_area)
    else:
        rectangularity = 0.0

    is_rectangular = (rectangularity > max_rectangularity_threshold)

    # 2. Douglas-Peucker Polygon Approximation & Corner Angle Analysis
    perimeter = cv2.arcLength(cnt, True)
    epsilon = dp_epsilon_ratio * perimeter
    approx = cv2.approxPolyDP(cnt, epsilon, True)
    num_vertices = len(approx)

    corner_angles: List[float] = []
    is_orthogonal_quad = False
    ortho_min, ortho_max = orthogonal_angle_range

    if num_vertices == 4:
        pts = approx.reshape((4, 2)).astype(np.float64)
        for i in range(4):
            p_prev = pts[(i - 1) % 4]
            p_curr = pts[i]
            p_next = pts[(i + 1) % 4]

            v1 = p_prev - p_curr
            v2 = p_next - p_curr

            norm1 = np.linalg.norm(v1)
            norm2 = np.linalg.norm(v2)

            if norm1 > 1e-6 and norm2 > 1e-6:
                dot_prod = np.dot(v1, v2)
                cos_theta = dot_prod / (norm1 * norm2)
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                angle_deg = math.degrees(math.acos(cos_theta))
            else:
                angle_deg = 0.0

            corner_angles.append(float(angle_deg))

        if len(corner_angles) == 4 and all(ortho_min <= ang <= ortho_max for ang in corner_angles):
            is_orthogonal_quad = True

    # 3. Convex Hull & Solidity Filter
    hull = cv2.convexHull(cnt)
    convex_hull_area = float(cv2.contourArea(hull))

    if convex_hull_area > 0.0:
        solidity = float(mask_area / convex_hull_area)
    else:
        solidity = 0.0

    is_valid_solidity = (solidity >= min_solidity_threshold)

    return GeometricMetrics(
        mask_area=mask_area,
        min_area_rect_area=min_area_rect_area,
        rectangularity=rectangularity,
        convex_hull_area=convex_hull_area,
        solidity=solidity,
        num_approx_vertices=num_vertices,
        corner_angles_deg=corner_angles,
        is_rectangular=is_rectangular,
        is_orthogonal_quad=is_orthogonal_quad,
        is_valid_solidity=is_valid_solidity
    )


def compute_spectral_metrics(
    candidate_patch: np.ndarray,
    candidate_mask: np.ndarray,
    high_sat_threshold: float = 0.45,
    max_color_swatch_ratio: float = 0.15,
    is_rgb: bool = False
) -> SpectralMetrics:
    """
    Analyze the colorimetric and HSV saturation profile of foreground mask pixels.
    
    Args:
        candidate_patch: Multi-channel image patch (H, W, 3).
        candidate_mask: Binary single-channel uint8 mask (H, W).
        high_sat_threshold: HSV saturation threshold (0.45).
        max_color_swatch_ratio: Max fraction of high-sat pixels (0.15).
        is_rgb: Set True if candidate_patch is RGB, False if BGR.
        
    Returns:
        SpectralMetrics dataclass containing saturation and color distribution statistics.
    """
    if candidate_patch is None or candidate_patch.size == 0 or candidate_mask is None:
        return SpectralMetrics(
            mean_hue=0.0,
            mean_saturation=0.0,
            mean_value=0.0,
            high_saturation_pixel_count=0,
            high_saturation_ratio=0.0,
            is_color_swatch=False
        )

    if is_rgb:
        bgr_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_RGB2BGR)
    else:
        bgr_patch = candidate_patch

    hsv_patch = cv2.cvtColor(bgr_patch, cv2.COLOR_BGR2HSV)
    fg_indices = np.where(candidate_mask > 0)
    total_fg_pixels = len(fg_indices[0])

    if total_fg_pixels == 0:
        return SpectralMetrics(
            mean_hue=0.0,
            mean_saturation=0.0,
            mean_value=0.0,
            high_saturation_pixel_count=0,
            high_saturation_ratio=0.0,
            is_color_swatch=False
        )

    fg_hsv = hsv_patch[fg_indices]
    h_vals = fg_hsv[:, 0].astype(np.float64) * 2.0
    s_vals = fg_hsv[:, 1].astype(np.float64) / 255.0
    v_vals = fg_hsv[:, 2].astype(np.float64) / 255.0

    mean_hue = float(np.mean(h_vals))
    mean_saturation = float(np.mean(s_vals))
    mean_value = float(np.mean(v_vals))

    high_sat_count = int(np.sum(s_vals > high_sat_threshold))
    high_sat_ratio = float(high_sat_count / total_fg_pixels)
    is_color_swatch = (high_sat_ratio > max_color_swatch_ratio)

    return SpectralMetrics(
        mean_hue=mean_hue,
        mean_saturation=mean_saturation,
        mean_value=mean_value,
        high_saturation_pixel_count=high_sat_count,
        high_saturation_ratio=high_sat_ratio,
        is_color_swatch=is_color_swatch
    )


def compute_texture_metrics(
    candidate_patch: np.ndarray,
    candidate_mask: np.ndarray,
    paper_mean_val_threshold: float = 205.0,
    paper_max_sat_threshold: float = 35.0,
    laplacian_text_var_threshold: float = 450.0,
    canny_text_edge_threshold: float = 0.15,
    is_rgb: bool = False
) -> TextureMetrics:
    """
    Quantify interior high-frequency gradient variance, typographic edge stroke density,
    and paper substrate gating.
    
    Args:
        candidate_patch: Multi-channel or grayscale image patch (H, W, 3) or (H, W).
        candidate_mask: Binary single-channel uint8 mask (H, W).
        paper_mean_val_threshold: Minimum HSV brightness V for paper substrate (205.0).
        paper_max_sat_threshold: Maximum HSV saturation S for paper substrate (35.0).
        laplacian_text_var_threshold: Minimum Laplacian variance for text (450.0).
        canny_text_edge_threshold: Minimum Canny edge density for text (0.15).
        is_rgb: Set True if candidate_patch is RGB, False if BGR.
        
    Returns:
        TextureMetrics dataclass summarizing gradient variance, edge density, and text indicators.
    """
    if candidate_patch is None or candidate_patch.size == 0 or candidate_mask is None:
        return TextureMetrics(
            laplacian_variance=0.0,
            canny_edge_density=0.0,
            horizontal_stroke_density=0.0,
            vertical_stroke_density=0.0,
            is_printed_text=False,
            mean_val=0.0,
            mean_sat=0.0,
            is_paper_substrate=False
        )

    if candidate_patch.ndim == 3:
        if is_rgb:
            bgr_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_RGB2BGR)
            gray_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_RGB2GRAY)
        else:
            bgr_patch = candidate_patch
            gray_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_BGR2GRAY)
        hsv_patch = cv2.cvtColor(bgr_patch, cv2.COLOR_BGR2HSV)
    else:
        gray_patch = candidate_patch.copy()
        bgr_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_GRAY2BGR)
        hsv_patch = cv2.cvtColor(bgr_patch, cv2.COLOR_BGR2HSV)

    erode_kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5, 5))
    interior_mask = cv2.erode((candidate_mask > 0).astype(np.uint8) * 255, erode_kernel)
    interior_fg_pixels = int(cv2.countNonZero(interior_mask))

    if interior_fg_pixels < 25:
        interior_mask = (candidate_mask > 0).astype(np.uint8) * 255
        interior_fg_pixels = int(cv2.countNonZero(interior_mask))

    if interior_fg_pixels == 0:
        return TextureMetrics(
            laplacian_variance=0.0,
            canny_edge_density=0.0,
            horizontal_stroke_density=0.0,
            vertical_stroke_density=0.0,
            is_printed_text=False,
            mean_val=0.0,
            mean_sat=0.0,
            is_paper_substrate=False
        )

    interior_pts = (interior_mask > 0)
    s_channel = hsv_patch[:, :, 1].astype(np.float64)
    v_channel = hsv_patch[:, :, 2].astype(np.float64)
    mean_sat = float(np.mean(s_channel[interior_pts]))
    mean_val = float(np.mean(v_channel[interior_pts]))

    is_paper_substrate = bool(
        mean_val > paper_mean_val_threshold and
        mean_sat < paper_max_sat_threshold
    )

    laplacian = cv2.Laplacian(gray_patch, cv2.CV_64F)
    interior_laplacian = laplacian[interior_pts]
    laplacian_variance = float(np.var(interior_laplacian)) if len(interior_laplacian) > 0 else 0.0

    canny_edges = cv2.Canny(gray_patch, 50, 150)
    interior_canny_count = int(np.sum((canny_edges > 0) & interior_pts))
    canny_edge_density = float(interior_canny_count / interior_fg_pixels)

    sobel_x = cv2.Sobel(gray_patch, cv2.CV_64F, 1, 0, ksize=3)
    sobel_y = cv2.Sobel(gray_patch, cv2.CV_64F, 0, 1, ksize=3)

    interior_sobel_x = np.abs(sobel_x[interior_pts])
    interior_sobel_y = np.abs(sobel_y[interior_pts])

    horizontal_stroke_density = float(np.mean(interior_sobel_y)) if len(interior_sobel_y) > 0 else 0.0
    vertical_stroke_density = float(np.mean(interior_sobel_x)) if len(interior_sobel_x) > 0 else 0.0

    is_printed_text = bool(
        is_paper_substrate and
        laplacian_variance > laplacian_text_var_threshold and
        canny_edge_density > canny_text_edge_threshold
    )

    return TextureMetrics(
        laplacian_variance=laplacian_variance,
        canny_edge_density=canny_edge_density,
        horizontal_stroke_density=horizontal_stroke_density,
        vertical_stroke_density=vertical_stroke_density,
        is_printed_text=is_printed_text,
        mean_val=mean_val,
        mean_sat=mean_sat,
        is_paper_substrate=is_paper_substrate
    )
