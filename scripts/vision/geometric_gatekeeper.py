#!/usr/bin/env python3
"""
scripts/vision/geometric_gatekeeper.py
======================================
Geometric quality gatekeeping and midrib pose estimation module.
Calculates Unoccluded Completeness Score (UCS), Solidity, and Longitudinal Midrib Axis.
"""

from __future__ import annotations

import math
from typing import Tuple

import cv2
import numpy as np


def compute_geometric_metrics_and_pose(
    mask: np.ndarray
) -> Tuple[float, float, float, Tuple[int, int], Tuple[int, int]]:
    """
    Computes:
      - UCS (Unoccluded Completeness Score): Area_mask / Area_expected
      - Solidity: Area_mask / Area_convex_hull
      - Midrib Angle (deg): Principal longitudinal inertia / line fit axis
      - Apex Point (p_apex)
      - Transition / Base Point (p_transition)
    """
    area_mask = float(np.count_nonzero(mask))
    if area_mask < 50:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    main_contour = max(contours, key=cv2.contourArea)

    # 1. Convex Hull & Solidity
    hull = cv2.convexHull(main_contour)
    area_hull = float(cv2.contourArea(hull))
    solidity = (area_mask / area_hull) if area_hull > 0 else 0.0
    solidity = min(max(solidity, 0.0), 1.0)

    # 2. Expected Area (Oriented Bounding Box / Minimum Area Rect)
    rect = cv2.minAreaRect(main_contour)
    rect_w, rect_h = rect[1]
    area_rect = max(rect_w * rect_h, 1.0)
    expected_area = 0.70 * area_rect
    ucs = min(max(area_mask / expected_area, 0.0), 1.0)

    # 3. Longitudinal Midrib Axis Fit
    line_params = cv2.fitLine(main_contour, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]
    angle_rad = math.atan2(vy, vx)
    angle_deg = math.degrees(angle_rad)

    # Find apex and base points along fitted line
    pts = main_contour.reshape(-1, 2).astype(np.float32)
    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy

    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))

    p_apex = (int(pts[min_idx, 0]), int(pts[min_idx, 1]))
    p_transition = (int(pts[max_idx, 0]), int(pts[max_idx, 1]))

    return ucs, solidity, angle_deg, p_apex, p_transition
