"""
===============================================================================
Module: sam2_geometry.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Geometric and contour processing routines for interactive botanical annotation,
    including knife-line mask severing, polygon contour simplification, freehand
    lasso rasterization, and normalized YOLO bounding box formatting.
===============================================================================
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import cv2
import numpy as np


def clip_box_to_image(
    box: Tuple[int, int, int, int],
    img_w: int,
    img_h: int
) -> Tuple[int, int, int, int]:
    """
    Clips bounding box coordinates (x0, y0, x1, y1) to valid image boundaries.
    """
    x0, y0, x1, y1 = box
    x_min = max(0, min(x0, x1))
    y_min = max(0, min(y0, y1))
    x_max = min(img_w, max(x0, x1))
    y_max = min(img_h, max(y0, y1))
    return x_min, y_min, x_max, y_max


def split_mask_with_knife_line(
    binary_mask: np.ndarray,
    line_start: Tuple[int, int],
    line_end: Tuple[int, int],
    line_thickness: int = 3
) -> np.ndarray:
    """
    Sever a binary mask using a knife cut line (e.g. to isolate blade from petiole/roots).

    Args:
        binary_mask: 2D uint8 binary mask array.
        line_start: (x, y) start coordinate.
        line_end: (x, y) end coordinate.
        line_thickness: Cut line stroke width in pixels.

    Returns:
        np.ndarray: Mask with knife cut line zeroed out.
    """
    cut_mask = binary_mask.copy()
    cv2.line(cut_mask, line_start, line_end, 0, thickness=line_thickness)
    return cut_mask


def rasterize_lasso_polygon(
    lasso_points: List[Tuple[int, int]],
    img_h: int,
    img_w: int
) -> np.ndarray:
    """
    Rasterizes a sequence of freehand lasso boundary points into a filled binary mask.

    Args:
        lasso_points: List of (x, y) polygon vertices.
        img_h: Image height in pixels.
        img_w: Image width in pixels.

    Returns:
        np.ndarray: 2D uint8 binary mask (0 background, 255 filled polygon).
    """
    mask = np.zeros((img_h, img_w), dtype=np.uint8)
    if len(lasso_points) < 3:
        return mask
    pts = np.array(lasso_points, dtype=np.int32).reshape((-1, 1, 2))
    cv2.fillPoly(mask, [pts], 255)
    return mask


def mask_to_yolo_bbox(
    binary_mask: np.ndarray,
    class_id: int
) -> Optional[Tuple[int, float, float, float, float]]:
    """
    Computes normalized YOLO format bounding box (class_id, x_center, y_center, width, height).

    Args:
        binary_mask: 2D uint8 binary mask.
        class_id: Integer taxonomic class index.

    Returns:
        Optional[Tuple[int, float, float, float, float]]: Normalized bounding box or None if empty.
    """
    y_indices, x_indices = np.where(binary_mask > 0)
    if len(x_indices) == 0 or len(y_indices) == 0:
        return None

    h, w = binary_mask.shape[:2]
    x_min, x_max = float(np.min(x_indices)), float(np.max(x_indices))
    y_min, y_max = float(np.min(y_indices)), float(np.max(y_indices))

    box_w = (x_max - x_min + 1.0) / w
    box_h = (y_max - y_min + 1.0) / h
    x_center = (x_min + (x_max - x_min) / 2.0) / w
    y_center = (y_min + (y_max - y_min) / 2.0) / h

    return (
        class_id,
        round(x_center, 6),
        round(y_center, 6),
        round(box_w, 6),
        round(box_h, 6),
    )


def mask_to_normalized_polygon(
    binary_mask: np.ndarray,
    class_id: int,
    approx_epsilon: float = 1.0
) -> Optional[str]:
    """
    Converts binary mask to normalized YOLO segmentation string: 'class_id x1 y1 x2 y2 ...'
    """
    contours, _ = cv2.findContours(binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None

    largest_cnt = max(contours, key=cv2.contourArea)
    if approx_epsilon > 0:
        largest_cnt = cv2.approxPolyDP(largest_cnt, approx_epsilon, True)

    if len(largest_cnt) < 3:
        return None

    h, w = binary_mask.shape[:2]
    pts = largest_cnt.reshape(-1, 2)
    norm_pts = []
    for x, y in pts:
        norm_pts.append(f"{round(x / w, 6):.6f}")
        norm_pts.append(f"{round(y / h, 6):.6f}")

    return f"{class_id} " + " ".join(norm_pts)


def polygon_to_bounding_box(
    polygon_points: List[Tuple[int, int]]
) -> Optional[Tuple[int, int, int, int]]:
    """
    Computes the enclosing bounding box (min_x, min_y, max_x, max_y) from a sequence of polygon vertices.

    Args:
        polygon_points: List of (x, y) integer pixel coordinates.

    Returns:
        Optional[Tuple[int, int, int, int]]: (min_x, min_y, max_x, max_y) or None if empty.
    """
    if not polygon_points:
        return None
    xs = [p[0] for p in polygon_points]
    ys = [p[1] for p in polygon_points]
    return min(xs), min(ys), max(xs), max(ys)


def polygon_interior_point(
    polygon_points: List[Tuple[int, int]],
    img_h: int,
    img_w: int
) -> Optional[Tuple[float, float]]:
    """
    Finds the optimal interior point (pole of inaccessibility) inside an arbitrary polygon
    using Euclidean distance transform.

    Args:
        polygon_points: List of (x, y) vertices defining the selection polygon.
        img_h: Image height in pixels.
        img_w: Image width in pixels.

    Returns:
        Optional[Tuple[float, float]]: (x, y) interior coordinate to guide SAM 2 prompt.
    """
    if len(polygon_points) < 3:
        return None
    mask = rasterize_lasso_polygon(polygon_points, img_h, img_w)
    if np.count_nonzero(mask) == 0:
        return None
    dist = cv2.distanceTransform(mask, cv2.DIST_L2, 5)
    _, max_val, _, max_loc = cv2.minMaxLoc(dist)
    if max_val > 0:
        return float(max_loc[0]), float(max_loc[1])
    return float(polygon_points[0][0]), float(polygon_points[0][1])

