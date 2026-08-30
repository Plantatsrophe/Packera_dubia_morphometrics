"""
scripts/core/leaf_spine_tracer.py
=================================
Botanical leaf anatomical spine tracing and Frangi vesselness filtering
following the LeafMachine2 3-point keypoint protocol.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple

import cv2
import numpy as np
from scipy import ndimage


def frangi_vesselness_filter_2d(
    image_gray: np.ndarray,
    sigmas: List[float] = [1.0, 2.0, 3.0]
) -> np.ndarray:
    """
    Computes 2D multi-scale Hessian ridge filter (Frangi vesselness)
    to trace thin, faint petioles obscured by arachnoid tomentum or mounting tape.
    
    Args:
        image_gray: Grayscale uint8 image.
        sigmas: List of Gaussian smoothing scales.
        
    Returns:
        Normalized vesselness response array (float32, [0, 1]).
    """
    max_response = np.zeros_like(image_gray, dtype=np.float32)
    inv_img = 255.0 - image_gray.astype(np.float32)

    for sigma in sigmas:
        smoothed = ndimage.gaussian_filter(inv_img, sigma=sigma)
        dy, dx = np.gradient(smoothed)
        dyy, dyx = np.gradient(dy)
        dxy, dxx = np.gradient(dx)

        trace = dxx + dyy
        det = dxx * dyy - dxy * dyx
        discriminant = np.sqrt(np.maximum(trace**2 - 4 * det, 0))

        lambda1 = 0.5 * (trace + discriminant)
        lambda2 = 0.5 * (trace - discriminant)

        c = 15.0
        beta = 0.5

        rb = np.abs(lambda1) / (np.abs(lambda2) + 1e-6)
        s2 = lambda1**2 + lambda2**2

        vesselness = np.exp(-(rb**2) / (2 * beta**2)) * (1.0 - np.exp(-s2 / (2 * c**2)))
        vesselness[lambda2 > 0] = 0.0

        max_response = np.maximum(max_response, vesselness)

    norm = max_response.max()
    if norm > 0:
        max_response = max_response / norm

    return max_response


def trace_3point_anatomical_spine(
    leaf_mask: np.ndarray,
    leaf_gray: np.ndarray
) -> Dict[str, Any]:
    """
    Traces the 3-point anatomical spine (LeafMachine2 protocol)
    modeling each basal leaf with three keypoints:
        - p_apex: Lamina apex (distal tip)
        - p_transition: Lamina-to-petiole junction (blade base expansion)
        - p_caudex: Petiole insertion point at the caudex / rootstock
        
    Args:
        leaf_mask: Binary leaf silhouette mask (uint8).
        leaf_gray: Grayscale leaf crop image (uint8).
        
    Returns:
        Dictionary with keypoint coordinates, sub-segment lengths, and total length.
    """
    h, w = leaf_mask.shape[:2]
    contours, _ = cv2.findContours(leaf_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return {
            "p_apex": (0, 0),
            "p_transition": (0, 0),
            "p_caudex": (0, 0),
            "petiole_length_px": 0.0,
            "lamina_length_px": 0.0,
            "total_spine_length_px": 0.0,
            "spine_points": []
        }

    main_cnt = max(contours, key=cv2.contourArea)
    pts = main_cnt.reshape(-1, 2).astype(np.float32)

    line_params = cv2.fitLine(main_cnt, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]

    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy
    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))

    pt_a = (int(pts[min_idx][0]), int(pts[min_idx][1]))
    pt_b = (int(pts[max_idx][0]), int(pts[max_idx][1]))

    M = cv2.moments(leaf_mask)
    if M["m00"] > 0:
        cx = M["m10"] / M["m00"]
        cy = M["m01"] / M["m00"]
        dist_a = math.hypot(pt_a[0] - cx, pt_a[0] - cy)
        dist_b = math.hypot(pt_b[0] - cx, pt_b[0] - cy)
        if dist_a > dist_b:
            p_caudex, p_apex = pt_a, pt_b
        else:
            p_caudex, p_apex = pt_b, pt_a
    else:
        p_apex, p_caudex = pt_a, pt_b

    spine_vec = np.array([p_caudex[0] - p_apex[0], p_caudex[1] - p_apex[1]], dtype=np.float32)
    spine_len = float(np.linalg.norm(spine_vec))
    if spine_len < 5:
        return {
            "p_apex": p_apex,
            "p_transition": p_apex,
            "p_caudex": p_caudex,
            "petiole_length_px": 0.0,
            "lamina_length_px": spine_len,
            "total_spine_length_px": spine_len,
            "spine_points": [p_apex, p_caudex]
        }

    spine_unit = spine_vec / spine_len
    perp_unit = np.array([-spine_unit[1], spine_unit[0]], dtype=np.float32)

    num_samples = 50
    widths = []
    sample_points = []
    for t in np.linspace(0.1, 0.95, num_samples):
        center_pt = np.array(p_apex, dtype=np.float32) + t * spine_vec
        sample_points.append(center_pt)

        span = 0
        for s in range(-int(w / 2), int(w / 2)):
            probe = center_pt + s * perp_unit
            px, py = int(round(probe[0])), int(round(probe[1]))
            if 0 <= px < w and 0 <= py < h and leaf_mask[py, px] > 0:
                span += 1
        widths.append(span)

    widths = np.array(widths, dtype=np.float32)
    max_w = np.max(widths) if len(widths) > 0 else 1.0

    narrow_indices = np.where(widths <= 0.35 * max_w)[0]
    if len(narrow_indices) > 0:
        trans_idx = narrow_indices[0]
        trans_pt = sample_points[trans_idx]
        p_transition = (int(round(trans_pt[0])), int(round(trans_pt[1])))
    else:
        mid_pt = np.array(p_apex, dtype=np.float32) + 0.60 * spine_vec
        p_transition = (int(round(mid_pt[0])), int(round(mid_pt[1])))

    vesselness = frangi_vesselness_filter_2d(leaf_gray, sigmas=[1.0, 2.0, 3.0])
    petiole_len = float(math.hypot(p_caudex[0] - p_transition[0], p_caudex[1] - p_transition[1]))
    lamina_len = float(math.hypot(p_transition[0] - p_apex[0], p_transition[1] - p_apex[1]))

    return {
        "p_apex": p_apex,
        "p_transition": p_transition,
        "p_caudex": p_caudex,
        "petiole_length_px": petiole_len,
        "lamina_length_px": lamina_len,
        "total_spine_length_px": spine_len,
        "vesselness_peak": float(vesselness.max()),
        "spine_points": [p_apex, p_transition, p_caudex]
    }
