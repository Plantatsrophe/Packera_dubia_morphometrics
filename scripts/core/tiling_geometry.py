import os
import sys
import logging
import math
import random
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
import numpy as np
import cv2
import json
import glob
from shapely.geometry import Polygon, MultiPolygon, box
from shapely.validation import make_valid

# Common imports
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

logger = setup_logging()


class HerbariumAnnotation:
    """
    Represents a single botanical or artifact instance on a full herbarium sheet.
    Stores class ID, original polygon contour, bounding box, and original area.
    """

    def __init__(
        self,
        class_id: int,
        polygon: Polygon,
        confidence: float = 1.0,
        original_line: str = ""
    ):
        """
        Initialize an annotation with geometric attributes.

        Args:
            class_id: Integer class index (0-8).
            polygon: Shapely Polygon object in absolute pixel coordinates.
            confidence: Detection score (1.0 for ground truth).
            original_line: Raw line string from source YOLO format file.
        """
        self.class_id = class_id
        # Ensure polygon geometry is geometrically valid (repair self-intersections)
        if not polygon.is_valid:
            polygon = make_valid(polygon)
        self.polygon = polygon
        self.confidence = confidence
        self.original_line = original_line
        self.original_area = polygon.area if hasattr(polygon, "area") else 0.0

    @classmethod
    def from_yolo_line(
        cls,
        line: str,
        img_width: int,
        img_height: int
    ) -> Optional["HerbariumAnnotation"]:
        """
        Parses a YOLO-format line (polygon segmentation or bounding box) into
        absolute full-sheet pixel coordinates.

        Args:
            line: Text line with space-separated values.
            img_width: Full sheet pixel width.
            img_height: Full sheet pixel height.

        Returns:
            HerbariumAnnotation instance or None if parsing fails.
        """
        parts = line.strip().split()
        if not parts:
            return None

        try:
            class_id = int(parts[0])
            coords = [float(p) for p in parts[1:]]

            # Case A: Standard YOLO Bounding Box (class_id, x_center, y_center, w, h)
            if len(coords) == 4:
                cx_norm, cy_norm, w_norm, h_norm = coords
                cx = cx_norm * img_width
                cy = cy_norm * img_height
                w = w_norm * img_width
                h = h_norm * img_height
                x1 = cx - w / 2.0
                y1 = cy - h / 2.0
                x2 = cx + w / 2.0
                y2 = cy + h / 2.0
                poly = box(x1, y1, x2, y2)
                return cls(class_id=class_id, polygon=poly, original_line=line)

            # Case B: YOLOv8 Segmentation Polygon (class_id, x1, y1, x2, y2, ..., xn, yn)
            elif len(coords) >= 6 and len(coords) % 2 == 0:
                pts = []
                for i in range(0, len(coords), 2):
                    px = coords[i] * img_width
                    py = coords[i + 1] * img_height
                    pts.append((px, py))

                if len(pts) < 3:
                    return None

                poly = Polygon(pts)
                if not poly.is_valid:
                    poly = make_valid(poly)
                return cls(class_id=class_id, polygon=poly, original_line=line)

            else:
                logger.warning("Unsupported YOLO line format with %d tokens: %s", len(coords), line)
                return None

        except Exception as err:
            logger.debug("Failed parsing YOLO annotation line: %s | Error: %s", line, err)
            return None


class DynamicGeometricReprojector:
    """
    Performs geometric intersection and re-projection of full-sheet annotation
    polygons into local tile coordinate systems $[0, W_{\\text{tile}}]$ and $[0, H_{\\text{tile}}]$.
    Filters partial instances where the visible tile area is $< 15\\%$ of original.
    """

    def __init__(self, min_area_ratio: float = 0.15):
        """
        Initialize geometric re-projector.

        Args:
            min_area_ratio: Minimum ratio of visible area within tile to original
                            area (default 0.15 = 15%).
        """
        self.min_area_ratio = float(min_area_ratio)

    def reproject_annotations_to_tile(
        self,
        annotations: List[HerbariumAnnotation],
        window: Tuple[int, int, int, int],
        tile_width: int,
        tile_height: int
    ) -> Tuple[List[Tuple[int, List[Tuple[float, float]]]], bool]:
        """
        Clips and transforms full-sheet annotations to local normalized tile coordinates.

        Args:
            annotations: List of HerbariumAnnotation objects on full sheet.
            window: Tile bounding box in full sheet space (x1, y1, x2, y2).
            tile_width: Width of the tile in pixels.
            tile_height: Height of the tile in pixels.

        Returns:
            Tuple containing:
              - List of (class_id, normalized_polygon_points) for surviving instances.
              - Boolean flag True if any polygon fragment was dropped due to min_area_ratio.
        """
        x1, y1, x2, y2 = window
        win_box = box(x1, y1, x2, y2)
        surviving_instances: List[Tuple[int, List[Tuple[float, float]]]] = []
        dropped_any_fragment = False

        for ann in annotations:
            # Check bounding box envelope overlap before expensive boolean operation
            if not ann.polygon.intersects(win_box):
                continue

            try:
                # Perform 2D geometric boolean intersection
                clipped = ann.polygon.intersection(win_box)
                if clipped.is_empty:
                    continue

                # Make valid if clipping produced non-standard geometry
                if not clipped.is_valid:
                    clipped = make_valid(clipped)

                # Extract individual polygons from MultiPolygon or GeometryCollection
                poly_parts: List[Polygon] = []
                if isinstance(clipped, Polygon):
                    poly_parts = [clipped]
                elif isinstance(clipped, MultiPolygon):
                    poly_parts = list(clipped.geoms)
                elif hasattr(clipped, "geoms"):
                    # GeometryCollection: extract polygon sub-elements
                    for g in clipped.geoms:
                        if isinstance(g, Polygon):
                            poly_parts.append(g)

                for part in poly_parts:
                    part_area = part.area
                    # Filter out degenerate/zero-area polygons
                    if part_area <= 1e-4:
                        continue

                    # Filter out partial instances where visible tile area is < 15% of original
                    visible_ratio = part_area / max(1e-6, ann.original_area)
                    if visible_ratio < self.min_area_ratio:
                        logger.debug(
                            "Filtered partial class %d: visible ratio %.3f < threshold %.2f",
                            ann.class_id, visible_ratio, self.min_area_ratio
                        )
                        dropped_any_fragment = True
                        continue

                    # Extract exterior boundary coordinates
                    ext_coords = list(part.exterior.coords)
                    if len(ext_coords) < 3:
                        continue

                    # Map coordinates from full sheet to local tile space [0, 1]
                    norm_pts: List[Tuple[float, float]] = []
                    for px, py in ext_coords:
                        # Translate by tile window origin
                        lx = px - x1
                        ly = py - y1
                        # Normalize to [0.0, 1.0] within local tile
                        nx = max(0.0, min(1.0, lx / float(tile_width)))
                        ny = max(0.0, min(1.0, ly / float(tile_height)))
                        norm_pts.append((round(nx, 6), round(ny, 6)))

                    # Deduplicate consecutive identical vertices
                    clean_pts: List[Tuple[float, float]] = []
                    for pt in norm_pts:
                        if not clean_pts or pt != clean_pts[-1]:
                            clean_pts.append(pt)

                    if len(clean_pts) >= 3:
                        surviving_instances.append((ann.class_id, clean_pts))

            except Exception as err:
                logger.debug("Error clipping polygon for class %d: %s", ann.class_id, err)
                continue

        return surviving_instances, dropped_any_fragment


class BackgroundPaperFilter:
    """
    Sub-samples pure empty paper tiles to avoid background imbalance during
    YOLO training. Retains only 5% of pure white background paper patches.
    """

    def __init__(self, keep_prob: float = 0.05, seed: int = 42):
        """
        Initialize background paper filter.

        Args:
            keep_prob: Retention probability for empty tiles (default 0.05 = 5%).
            seed: Random seed for reproducible sub-sampling.
        """
        self.keep_prob = float(keep_prob)
        self.rng = random.Random(seed)

    def should_keep_empty_tile(self, tile_image: np.ndarray) -> bool:
        """
        Evaluates an empty tile (0 annotations) and decides whether to keep it as
        a hard negative background sample.
        Applies a max local standard deviation filter (64x64 grid, threshold 15.0)
        to strictly discard tiles containing artifacts before randomly sub-sampling.

        Args:
            tile_image: RGB or BGR numpy image tile array.

        Returns:
            True if tile should be retained, False to discard.
        """
        # 1. Max Local Standard Deviation Gating
        tile_gray = cv2.cvtColor(tile_image, cv2.COLOR_BGR2GRAY)
        h, w = tile_gray.shape
        grid_size = 64
        
        # Ensure dimensions are divisible by grid_size (e.g. 1024 / 64 = 16)
        if h % grid_size == 0 and w % grid_size == 0:
            grid = tile_gray.reshape(h // grid_size, grid_size, w // grid_size, grid_size)
            local_stds = grid.std(axis=(1, 3))
            max_local_std = local_stds.max()
        else:
            # Fallback if tile was not padded correctly (should not happen in our pipeline)
            max_local_std = tile_gray.std()
            
        if max_local_std > 15.0:
            # Reject tile: contains too much structural variance (roots, tape, shadows)
            return False

        # 2. Random Sub-sampling for pure paper
        if self.rng.random() <= self.keep_prob:
            return True
        return False
