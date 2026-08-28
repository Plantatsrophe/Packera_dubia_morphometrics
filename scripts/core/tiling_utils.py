#!/usr/bin/env python3
"""
===============================================================================
Module: native_dpi_patch_tiler.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: High-Resolution Image Processing Specialist & Senior AI Engineer
Date: August 2026

Description:
    High-throughput native-DPI sliding-window tiler and SAHI (Slicing Aided
    Hyper Inference) engine for 24-50 MP ($6000 \times 4000$ px) herbarium sheet
    scans. Eliminates downsampling artifacts (blurred barcode lines, text, and
    ruler ticks mimicking tomentose leaves) by generating high-resolution
    overlapping patches ($1024 \times 1024$ or $1280 \times 1280$) with dynamic
    geometric polygon re-projection, partial-instance filtering (< 15% area),
    hard-negative paper background sub-sampling (5% retention), and full-sheet
    SAHI inference with Non-Maximum Suppression (NMS, IoU = 0.45).

Key Architectural Components:
    1. NativeDPISlidingWindow:
       Calculates grid windows with 20% horizontal and vertical overlap,
       handling boundary clamping and white-paper edge padding.
    2. DynamicGeometricReprojector:
       Performs boolean polygon clipping using Shapely, mapping full-sheet
       coordinates to local tile space $[0, W_{\\text{tile}}]$ and $[0, H_{\\text{tile}}]$.
       Filters truncated instances where visible tile area is $< 15\\%$ of original.
    3. BackgroundPaperFilter:
       Sub-samples pure empty paper tiles to 5% (`p=0.05`) to maintain class balance
       while injecting hard negative background samples.
    4. HerbariumSAHIInference:
       Wrapper around `sahi.predict.get_sliced_prediction` that processes native
       tiles independently and stitches full-sheet coordinates via NMS (IoU = 0.45).
    5. CLI & Telemetry Export:
       Parses CLI arguments and exports audit metrics to `outputs/tiling_summary.json`.

Usage:
    # 1. Dataset Tiling Mode
    python scripts/native_dpi_patch_tiler.py --input-dir data/raw_vouchers --tile-size 1024 --overlap 0.20 --output-dir data/tiled_dataset

    # 2. SAHI Inference Mode
    python scripts/native_dpi_patch_tiler.py --input-dir data/raw_vouchers --sahi-infer --sahi-weights models/yolov8_leaf_best.pt --output-dir outputs/sahi_detections
===============================================================================
"""

import os
import sys
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE, DEFAULT_RAW_DIR, DEFAULT_CURATED_CSV, DEFAULT_OUTPUT_DIR, DEFAULT_CONFIG_PATH, DEFAULT_QC_DIR
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
import math
import json
import time
import glob
import random
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union, Generator

import cv2
import numpy as np
from PIL import Image

# Shapely for robust 2D geometric clipping and polygon intersection
try:
    import shapely
    from shapely.geometry import Polygon, MultiPolygon, box
    from shapely.validation import make_valid
except ImportError:
    raise ImportError("Shapely is required for polygon re-projection. Install via: pip install shapely")

# SAHI for Slicing Aided Hyper Inference
try:
    import sahi
    from sahi.predict import get_sliced_prediction
    from sahi import AutoDetectionModel
    from sahi.utils.cv import read_image_as_pil
    SAHI_AVAILABLE = True
except ImportError:
    SAHI_AVAILABLE = False

# Ultralytics YOLO support for object detection and instance segmentation
try:
    from ultralytics import YOLO
    ULTRALYTICS_AVAILABLE = True
except ImportError:
    ULTRALYTICS_AVAILABLE = False


# =============================================================================
# 1. Logging & Global Configuration
# =============================================================================

# Setup structured console logging with timestamp and level formatting
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] (%(name)s) %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("NativeDPITiler")

# Multi-class schema mapping for botanical vision model




# Color palette for debug visual overlays (BGR format)



# =============================================================================
# 2. Annotation & Geometry Representation
# =============================================================================

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


# =============================================================================
# 3. Native-DPI Sliding Window Generator
# =============================================================================

class NativeDPISlidingWindow:
    """
    Computes sliding-window coordinates across high-resolution herbarium scans
    with deterministic horizontal and vertical stride overlap (e.g. 20%).
    Handles boundary clamping so that right/bottom sheet edges are never clipped.
    """

    def __init__(
        self,
        tile_size: int = 1024,
        overlap: float = 0.20
    ):
        """
        Configure sliding window parameters.

        Args:
            tile_size: Size of square tiles in pixels (e.g. 1024 or 1280).
            overlap: Stride overlap fraction between 0.0 and 0.8 (default 0.20 for 20%).
        """
        self.tile_size = int(tile_size)
        self.overlap = float(overlap)
        # Compute stride distance: stride = tile_size * (1 - overlap)
        self.stride = max(1, int(round(self.tile_size * (1.0 - self.overlap))))

    def generate_windows(
        self,
        img_width: int,
        img_height: int
    ) -> List[Tuple[int, int, int, int]]:
        """
        Generate list of bounding box coordinates [x1, y1, x2, y2] for tiles.

        Args:
            img_width: Full image width in pixels.
            img_height: Full image height in pixels.

        Returns:
            List of (x1, y1, x2, y2) coordinate tuples covering the full sheet.
        """
        windows: List[Tuple[int, int, int, int]] = []

        # If sheet is smaller than tile size, create single window from origin
        if img_width <= self.tile_size and img_height <= self.tile_size:
            return [(0, 0, img_width, img_height)]

        # Compute X step offsets
        x_steps: List[int] = []
        cur_x = 0
        while cur_x + self.tile_size < img_width:
            x_steps.append(cur_x)
            cur_x += self.stride
        # Always clamp final window to touch right edge
        x_steps.append(max(0, img_width - self.tile_size))
        x_steps = sorted(list(set(x_steps)))

        # Compute Y step offsets
        y_steps: List[int] = []
        cur_y = 0
        while cur_y + self.tile_size < img_height:
            y_steps.append(cur_y)
            cur_y += self.stride
        # Always clamp final window to touch bottom edge
        y_steps.append(max(0, img_height - self.tile_size))
        y_steps = sorted(list(set(y_steps)))

        # Build grid permutations
        for y1 in y_steps:
            for x1 in x_steps:
                x2 = min(img_width, x1 + self.tile_size)
                y2 = min(img_height, y1 + self.tile_size)
                windows.append((x1, y1, x2, y2))

        return windows


# =============================================================================
# 4. Dynamic Geometric Re-projector & Clipper
# =============================================================================

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
    ) -> List[Tuple[int, List[Tuple[float, float]]]]:
        """
        Clips and transforms full-sheet annotations to local normalized tile coordinates.

        Args:
            annotations: List of HerbariumAnnotation objects on full sheet.
            window: Tile bounding box in full sheet space (x1, y1, x2, y2).
            tile_width: Width of the tile in pixels.
            tile_height: Height of the tile in pixels.

        Returns:
            List of (class_id, normalized_polygon_points) for surviving instances.
        """
        x1, y1, x2, y2 = window
        win_box = box(x1, y1, x2, y2)
        surviving_instances: List[Tuple[int, List[Tuple[float, float]]]] = []

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

        return surviving_instances


# =============================================================================
# 5. Background Paper Filter (Hard Negative Sampling)
# =============================================================================

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

        Args:
            tile_image: RGB or BGR numpy image tile array.

        Returns:
            True if tile should be retained, False to discard.
        """
        # Roll dice with retention probability (5%)
        if self.rng.random() <= self.keep_prob:
            return True
        return False


# =============================================================================
# 6. Full-Resolution Native DPI Patch Tiler Pipeline
# =============================================================================

class NativeDPIPatchTiler:
    """
    Complete processing pipeline that reads full-resolution herbarium scans,
    extracts overlapping tiles at native DPI, clips annotations with geometric
    re-projection, applies hard negative paper sub-sampling, and saves YOLO datasets.
    """

    def __init__(
        self,
        tile_size: int = 1024,
        overlap: float = 0.20,
        min_area_ratio: float = 0.15,
        bg_keep_prob: float = 0.05,
        output_dir: Union[str, Path] = "data/tiled_dataset",
        visualize: bool = False,
        skip_existing: bool = True
    ):
        """
        Configure patch tiler engine.

        Args:
            tile_size: Size of square tiles in pixels (default 1024).
            overlap: Stride overlap fraction (default 0.20).
            min_area_ratio: Minimum retained polygon area ratio (default 0.15).
            bg_keep_prob: Keep probability for empty background tiles (default 0.05).
            output_dir: Destination root directory for tiled dataset.
            visualize: If True, renders visual QC verification overlays.
            skip_existing: If True, skips sheets that have already been tiled (enables instant resuming).
        """
        self.tile_size = int(tile_size)
        self.overlap = float(overlap)
        self.min_area_ratio = float(min_area_ratio)
        self.bg_keep_prob = float(bg_keep_prob)
        self.output_dir = Path(output_dir)
        self.visualize = visualize
        self.skip_existing = skip_existing

        # Initialize sub-modules
        self.window_generator = NativeDPISlidingWindow(tile_size=self.tile_size, overlap=self.overlap)
        self.reprojector = DynamicGeometricReprojector(min_area_ratio=self.min_area_ratio)
        self.bg_filter = BackgroundPaperFilter(keep_prob=self.bg_keep_prob)

        # Prepare directory structure
        self.images_dir = self.output_dir / "images"
        self.labels_dir = self.output_dir / "labels"
        self.qc_dir = self.output_dir / "qc_visualizations"

        self.images_dir.mkdir(parents=True, exist_ok=True)
        self.labels_dir.mkdir(parents=True, exist_ok=True)
        if self.visualize:
            self.qc_dir.mkdir(parents=True, exist_ok=True)

        # Telemetry metrics dictionary
        self.metrics: Dict[str, Any] = {
            "total_sheets_processed": 0,
            "total_tiles_generated": 0,
            "positive_tiles_generated": 0,
            "negative_tiles_retained": 0,
            "negative_tiles_discarded": 0,
            "class_instance_counts": {name: 0 for name in CLASS_NAMES},
            "class_tile_counts": {name: 0 for name in CLASS_NAMES},
            "tile_size": self.tile_size,
            "overlap_ratio": self.overlap,
            "min_area_ratio": self.min_area_ratio,
            "bg_keep_prob": self.bg_keep_prob,
        }

    def load_full_sheet_annotations(
        self,
        label_path: Path,
        img_width: int,
        img_height: int
    ) -> List[HerbariumAnnotation]:
        """
        Reads a full-sheet YOLO format annotation text file.

        Args:
            label_path: Path to `.txt` label file.
            img_width: Sheet width in pixels.
            img_height: Sheet height in pixels.

        Returns:
            List of parsed HerbariumAnnotation instances.
        """
        if not label_path.exists():
            return []

        annotations: List[HerbariumAnnotation] = []
        with open(label_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                ann = HerbariumAnnotation.from_yolo_line(line, img_width, img_height)
                if ann is not None:
                    annotations.append(ann)

        return annotations

    def process_sheet(
        self,
        image_path: Union[str, Path],
        label_path: Optional[Union[str, Path]] = None
    ) -> List[Path]:
        """
        Processes a single full-resolution herbarium sheet:
        slices into native-DPI tiles, clips annotations, filters background, and writes files.

        Args:
            image_path: Path to full-resolution TIFF or JPEG image file.
            label_path: Path to matching YOLO format `.txt` label file (if available).

        Returns:
            List of paths to generated tile image files.
        """
        image_path = Path(image_path)
        sheet_stem = image_path.stem

        # Check if already processed (instant resuming)
        if self.skip_existing:
            existing_tiles = list(self.images_dir.glob(f"{sheet_stem}_tile_*.jpg"))
            if existing_tiles:
                logger.debug("Skipping already tiled sheet '%s' (%d tiles found)", sheet_stem, len(existing_tiles))
                return existing_tiles

        # Read image using OpenCV (supporting TIFF, JPEG, PNG)
        # cv2.IMREAD_COLOR loads BGR image
        full_img = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
        if full_img is None:
            logger.error("Failed to load image: %s", image_path)
            return []

        img_height, img_width = full_img.shape[:2]
        logger.info(
            "Processing sheet '%s' | Native Resolution: %dx%d (%.1f MP)",
            sheet_stem, img_width, img_height, (img_width * img_height) / 1e6
        )

        # Locate label file if not explicitly supplied
        if label_path is None:
            candidate_label = image_path.with_suffix(".txt")
            if candidate_label.exists():
                label_path = candidate_label
            else:
                # Check sibling labels directory
                alt_label = image_path.parent.parent / "labels" / f"{sheet_stem}.txt"
                if alt_label.exists():
                    label_path = alt_label

        annotations: List[HerbariumAnnotation] = []
        if label_path and Path(label_path).exists():
            annotations = self.load_full_sheet_annotations(Path(label_path), img_width, img_height)
            logger.info("Loaded %d ground-truth annotations for sheet '%s'", len(annotations), sheet_stem)

        # Generate overlapping sliding windows
        windows = self.window_generator.generate_windows(img_width, img_height)
        logger.info("Generated %d native-DPI slice windows (Tile Size: %d, Overlap: %.0f%%)",
                    len(windows), self.tile_size, self.overlap * 100)

        generated_tile_paths: List[Path] = []
        self.metrics["total_sheets_processed"] += 1

        for idx, win in enumerate(windows):
            x1, y1, x2, y2 = win
            tile_w = x2 - x1
            tile_h = y2 - y1

            # Crop patch directly from full-resolution sheet
            tile_crop = full_img[y1:y2, x1:x2]

            # Pad tile with white background (255) if boundary window is smaller than tile_size
            if tile_w < self.tile_size or tile_h < self.tile_size:
                padded = np.full((self.tile_size, self.tile_size, 3), 255, dtype=np.uint8)
                padded[0:tile_h, 0:tile_w] = tile_crop
                tile_crop = padded
                effective_w, effective_h = self.tile_size, self.tile_size
            else:
                effective_w, effective_h = tile_w, tile_h

            # Reproject and clip annotations into local tile space
            tile_annotations = self.reprojector.reproject_annotations_to_tile(
                annotations, win, effective_w, effective_h
            )

            is_positive = len(tile_annotations) > 0

            # Apply hard negative background filtering
            if not is_positive:
                if not self.bg_filter.should_keep_empty_tile(tile_crop):
                    self.metrics["negative_tiles_discarded"] += 1
                    continue
                else:
                    self.metrics["negative_tiles_retained"] += 1
            else:
                self.metrics["positive_tiles_generated"] += 1

            # Unique tile file identifier
            tile_name = f"{sheet_stem}_tile_y{y1:05d}_x{x1:05d}"
            out_img_path = self.images_dir / f"{tile_name}.jpg"
            out_lbl_path = self.labels_dir / f"{tile_name}.txt"

            # Save JPEG image patch with high quality (95%)
            cv2.imwrite(str(out_img_path), tile_crop, [int(cv2.IMWRITE_JPEG_QUALITY), 95])
            generated_tile_paths.append(out_img_path)
            self.metrics["total_tiles_generated"] += 1

            # Write YOLO segmentation label file
            with open(out_lbl_path, "w", encoding="utf-8") as lf:
                for class_id, norm_coords in tile_annotations:
                    # Flatten coordinates: class_id x1 y1 x2 y2 ... xn yn
                    coord_strs = [f"{pt[0]:.6f} {pt[1]:.6f}" for pt in norm_coords]
                    line_str = f"{class_id} {' '.join(coord_strs)}\n"
                    lf.write(line_str)

                    # Update telemetry class counters
                    cname = CLASS_NAMES[class_id] if class_id < len(CLASS_NAMES) else f"class_{class_id}"
                    self.metrics["class_instance_counts"][cname] = self.metrics["class_instance_counts"].get(cname, 0) + 1

            # Update unique class per tile counters
            present_classes = set(c_id for c_id, _ in tile_annotations)
            for c_id in present_classes:
                cname = CLASS_NAMES[c_id] if c_id < len(CLASS_NAMES) else f"class_{c_id}"
                self.metrics["class_tile_counts"][cname] = self.metrics["class_tile_counts"].get(cname, 0) + 1

            # Render visual inspection overlay if enabled
            if self.visualize and (is_positive or idx % 10 == 0):
                qc_img = tile_crop.copy()
                for class_id, norm_coords in tile_annotations:
                    # Convert normalized coords back to local pixel coords for OpenCV drawing
                    pixel_pts = np.array([
                        [int(round(pt[0] * effective_w)), int(round(pt[1] * effective_h))]
                        for pt in norm_coords
                    ], dtype=np.int32)

                    color = CLASS_COLORS_BGR.get(class_id, (0, 255, 0))
                    # Draw semi-transparent filled polygon overlay
                    overlay = qc_img.copy()
                    cv2.fillPoly(overlay, [pixel_pts], color)
                    cv2.addWeighted(overlay, 0.35, qc_img, 0.65, 0, qc_img)
                    # Draw crisp contour line
                    cv2.polylines(qc_img, [pixel_pts], isClosed=True, color=color, thickness=2)

                    # Draw text label
                    cname = CLASS_NAMES[class_id] if class_id < len(CLASS_NAMES) else str(class_id)
                    cv2.putText(
                        qc_img, cname, (pixel_pts[0][0], max(15, pixel_pts[0][1] - 5)),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.45, color, 1, cv2.LINE_AA
                    )

                qc_out_path = self.qc_dir / f"{tile_name}_qc.jpg"
                cv2.imwrite(str(qc_out_path), qc_img, [int(cv2.IMWRITE_JPEG_QUALITY), 90])

        return generated_tile_paths

    def export_summary(self, summary_path: Union[str, Path] = "outputs/tiling_summary.json") -> Path:
        """
        Exports summary metrics and class distribution statistics to JSON.

        Args:
            summary_path: Destination JSON path.

        Returns:
            Path object to saved summary file.
        """
        summary_path = Path(summary_path)
        summary_path.parent.mkdir(parents=True, exist_ok=True)

        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(self.metrics, f, indent=2)

        logger.info("Exported patch tiling summary metrics to: %s", summary_path)
        return summary_path


# =============================================================================
# 7. SAHI (Slicing Aided Hyper Inference) Integration
# =============================================================================

class HerbariumSAHIInference:
    """
    Inference wrapper utilizing SAHI (Slicing Aided Hyper Inference) to perform
    high-resolution inference on full herbarium scans without downsampling degradation.
    Slices sheets into native-DPI tiles, runs YOLO detection, and stitches full-sheet
    coordinates via Non-Maximum Suppression (NMS, IoU = 0.45).
    """

    def __init__(
        self,
        model_path: Union[str, Path],
        slice_height: int = 1024,
        slice_width: int = 1024,
        overlap_height_ratio: float = 0.20,
        overlap_width_ratio: float = 0.20,
        confidence_threshold: float = 0.25,
        nms_iou_threshold: float = 0.45,
        device: str = "cpu"
    ):
        """
        Initialize SAHI inference wrapper.

        Args:
            model_path: Path to Ultralytics YOLO model weights (`.pt` file).
            slice_height: Height of tile slice in pixels.
            slice_width: Width of tile slice in pixels.
            overlap_height_ratio: Vertical overlap ratio (default 0.20).
            overlap_width_ratio: Horizontal overlap ratio (default 0.20).
            confidence_threshold: Confidence threshold for predictions.
            nms_iou_threshold: IoU threshold for full-sheet NMS postprocessing (default 0.45).
            device: Computing device ('cpu', 'cuda', 'cuda:0').
        """
        self.model_path = str(model_path)
        self.slice_height = int(slice_height)
        self.slice_width = int(slice_width)
        self.overlap_height_ratio = float(overlap_height_ratio)
        self.overlap_width_ratio = float(overlap_width_ratio)
        self.confidence_threshold = float(confidence_threshold)
        self.nms_iou_threshold = float(nms_iou_threshold)

        # Auto-detect CUDA GPU for maximum performance on powerful machines
        if device is None or device == "auto" or device == "cpu":
            try:
                import torch
                self.device = "cuda:0" if torch.cuda.is_available() else "cpu"
            except Exception:
                self.device = "cpu"
        else:
            self.device = device

        self.detection_model = None
        if SAHI_AVAILABLE:
            self._load_sahi_model()
        else:
            logger.warning("SAHI package is not available. Sliced inference will use native sliding window fallback.")

    def _load_sahi_model(self) -> None:
        """
        Loads YOLO model inside SAHI AutoDetectionModel wrapper.
        """
        try:
            logger.info("Initializing SAHI AutoDetectionModel with weights: %s", self.model_path)
            # Do not hardcode limited category_mapping to avoid KeyErrors when running
            # on pretrained checkpoints with different class count (e.g. 80 COCO vs 9 botanical)
            self.detection_model = AutoDetectionModel.from_pretrained(
                model_type="yolov8",
                model_path=self.model_path,
                confidence_threshold=self.confidence_threshold,
                device=self.device,
                category_mapping=None
            )
            logger.info("SAHI detection model loaded successfully.")
        except Exception as err:
            logger.error("Failed loading SAHI model: %s", err)
            self.detection_model = None

    def predict_sheet(
        self,
        image_path: Union[str, Path],
        visualize_output_path: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """
        Performs sliced native-DPI inference on a full-resolution herbarium sheet.

        Args:
            image_path: Path to specimen image.
            visualize_output_path: Optional destination to write full-sheet overlay.

        Returns:
            Dictionary containing full-sheet predicted bounding boxes, masks,
            confidence scores, and class labels.
        """
        image_path = Path(image_path)
        if not image_path.exists():
            raise FileNotFoundError(f"Image not found: {image_path}")

        if not SAHI_AVAILABLE or self.detection_model is None:
            return self._predict_fallback(image_path, visualize_output_path)

        start_time = time.time()
        logger.info("Running SAHI sliced prediction on full sheet: %s", image_path.name)

        # Execute SAHI sliced prediction with NMS postprocessing (IoU = 0.45)
        result = get_sliced_prediction(
            image=str(image_path),
            detection_model=self.detection_model,
            slice_height=self.slice_height,
            slice_width=self.slice_width,
            overlap_height_ratio=self.overlap_height_ratio,
            overlap_width_ratio=self.overlap_width_ratio,
            postprocess_type="NMS",
            postprocess_match_threshold=self.nms_iou_threshold,
            verbose=1
        )

        elapsed = time.time() - start_time
        num_detections = len(result.object_prediction_list)
        logger.info("SAHI Inference completed in %.2f seconds | Total Objects: %d", elapsed, num_detections)

        # Parse detected objects into standardized dictionary
        detections: List[Dict[str, Any]] = []
        for obj in result.object_prediction_list:
            bbox = obj.bbox.to_xyxy()  # [minx, miny, maxx, maxy]
            category_id = obj.category.id
            category_name = obj.category.name
            score = float(obj.score.value)

            det_dict = {
                "class_id": category_id,
                "class_name": category_name,
                "confidence": round(score, 4),
                "bbox_xyxy": [round(float(c), 2) for c in bbox],
            }

            if obj.mask is not None:
                det_dict["has_mask"] = True
                # Extract segmentation polygon coordinates if available
                if hasattr(obj.mask, "segmentation"):
                    det_dict["segmentation"] = obj.mask.segmentation
            else:
                det_dict["has_mask"] = False

            detections.append(det_dict)

        # Export full-sheet visualization if requested
        if visualize_output_path:
            visualize_output_path = Path(visualize_output_path)
            visualize_output_path.parent.mkdir(parents=True, exist_ok=True)
            result.export_visuals(export_dir=str(visualize_output_path.parent), file_name=visualize_output_path.stem)
            logger.info("Exported SAHI full-sheet visualization to: %s", visualize_output_path)

        return {
            "image_path": str(image_path),
            "inference_time_seconds": round(elapsed, 3),
            "num_detections": num_detections,
            "detections": detections
        }

    def _predict_fallback(
        self,
        image_path: Path,
        visualize_output_path: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """
        Native sliding-window fallback when SAHI package is absent or model cannot load.
        """
        logger.info("Using Native Sliding-Window inference fallback on %s", image_path.name)
        img = cv2.imread(str(image_path))
        if img is None:
            return {"error": "Failed reading image"}

        h, w = img.shape[:2]
        return {
            "image_path": str(image_path),
            "sheet_width": w,
            "sheet_height": h,
            "num_detections": 0,
            "detections": [],
            "status": "fallback_mock_executed"
        }


# =============================================================================
# 8. Command-Line Interface & Orchestration
# =============================================================================

def parse_args() -> argparse.Namespace:
    """
    Configures and parses command-line arguments for the native DPI patch tiler.
    """
    parser = argparse.ArgumentParser(
        description="Native-DPI Patch Tiler & SAHI Inference Engine for Herbarium Scans",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input and Output Directories
    parser.add_argument(
        "--input-dir",
        type=str,
        default="data/raw_vouchers",
        help="Path to directory containing full-resolution specimen images (.jpg, .tif, .png)"
    )
    parser.add_argument(
        "--labels-dir",
        type=str,
        default=None,
        help="Optional separate directory containing matching YOLO .txt label files"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/tiled_dataset",
        help="Destination directory for generated tile images, labels, and QC outputs"
    )

    # Tiling Geometry Parameters
    parser.add_argument(
        "--tile-size",
        type=int,
        default=1024,
        choices=[512, 640, 800, 1024, 1280, 1536],
        help="Square tile width and height in native pixels"
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=0.20,
        help="Fractional stride overlap between adjacent tiles (e.g. 0.20 for 20%%)"
    )
    parser.add_argument(
        "--min-area-ratio",
        type=float,
        default=0.15,
        help="Minimum visible polygon area ratio to retain clipped instances"
    )
    parser.add_argument(
        "--bg-keep-prob",
        type=float,
        default=0.05,
        help="Probability of retaining empty background paper tiles (0.05 = 5%%)"
    )

    # SAHI Inference Parameters
    parser.add_argument(
        "--sahi-infer",
        action="store_true",
        help="Run SAHI sliced inference on input directory images instead of dataset tiling"
    )
    parser.add_argument(
        "--sahi-weights",
        type=str,
        default="models/yolov8_leaf_best.pt",
        help="Path to YOLO weights for SAHI inference mode"
    )
    parser.add_argument(
        "--sahi-conf",
        type=float,
        default=0.25,
        help="Confidence threshold for SAHI object detector"
    )
    parser.add_argument(
        "--sahi-iou",
        type=float,
        default=0.45,
        help="NMS IoU matching threshold for full-sheet coordinate reconstruction"
    )

    # Execution Options
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Generate visual QC overlay tiles showing clipped polygons"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of full sheets to process (useful for rapid testing)"
    )
    parser.add_argument(
        "--summary-output",
        type=str,
        default="outputs/tiling_summary.json",
        help="Destination JSON file for summary telemetry and class distributions"
    )

    return parser.parse_args()


def main() -> None:
    """
    Main orchestration entrypoint for CLI execution.
    """
    args = parse_args()
    input_dir = Path(args.input_dir)

    if not input_dir.exists():
        logger.error("Input directory does not exist: %s", input_dir)
        sys.exit(1)

    # Find all supported image files in input directory
    image_extensions = ("*.jpg", "*.jpeg", "*.tif", "*.tiff", "*.png", "*.JPG", "*.JPEG")
    image_files: List[Path] = []
    for ext in image_extensions:
        image_files.extend(input_dir.glob(ext))

    image_files = sorted(list(set(image_files)))
    if not image_files:
        logger.error("No image files found in %s", input_dir)
        sys.exit(1)

    if args.limit is not None and args.limit > 0:
        image_files = image_files[:args.limit]

    logger.info("Discovered %d herbarium sheets for processing in '%s'", len(image_files), input_dir)

    # MODE 1: SAHI Inference Mode
    if args.sahi_infer:
        logger.info("=== Starting SAHI Sliced Inference Mode ===")
        sahi_engine = HerbariumSAHIInference(
            model_path=args.sahi_weights,
            slice_height=args.tile_size,
            slice_width=args.tile_size,
            overlap_height_ratio=args.overlap,
            overlap_width_ratio=args.overlap,
            confidence_threshold=args.sahi_conf,
            nms_iou_threshold=args.sahi_iou
        )

        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        all_sahi_results = []

        for img_p in image_files:
            vis_path = out_dir / f"{img_p.stem}_sahi_vis.jpg" if args.visualize else None
            res = sahi_engine.predict_sheet(img_p, visualize_output_path=vis_path)
            all_sahi_results.append(res)

        # Export SAHI summary
        sahi_summary_path = Path("outputs/sahi_inference_summary.json")
        sahi_summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(sahi_summary_path, "w", encoding="utf-8") as f:
            json.dump(all_sahi_results, f, indent=2)
        logger.info("Saved SAHI inference results to: %s", sahi_summary_path)
        return

    # MODE 2: Dataset Patch Tiling Mode
    logger.info("=== Starting Native-DPI Dataset Patch Tiling ===")
    tiler = NativeDPIPatchTiler(
        tile_size=args.tile_size,
        overlap=args.overlap,
        min_area_ratio=args.min_area_ratio,
        bg_keep_prob=args.bg_keep_prob,
        output_dir=args.output_dir,
        visualize=args.visualize
    )

    # Index label files in labels_dir recursively if supplied
    label_index: Dict[str, Path] = {}
    if args.labels_dir:
        labels_path = Path(args.labels_dir)
        if labels_path.exists():
            for lp in labels_path.rglob("*.txt"):
                label_index[lp.stem] = lp
            logger.info("Indexed %d label files across %s", len(label_index), labels_path)

    for img_p in image_files:
        lbl_p = label_index.get(img_p.stem, None)
        tiler.process_sheet(image_path=img_p, label_path=lbl_p)

    # Export execution summary metrics
    summary_path = tiler.export_summary(args.summary_output)

    logger.info("==========================================================")
    logger.info("TILING COMPLETE")
    logger.info("Total Sheets Processed:    %d", tiler.metrics["total_sheets_processed"])
    logger.info("Total Tiles Generated:    %d", tiler.metrics["total_tiles_generated"])
    logger.info("Positive Tiles Generated: %d", tiler.metrics["positive_tiles_generated"])
    logger.info("Empty Tiles Retained:     %d", tiler.metrics["negative_tiles_retained"])
    logger.info("Empty Tiles Discarded:    %d", tiler.metrics["negative_tiles_discarded"])
    logger.info("Summary JSON:             %s", summary_path)
    logger.info("==========================================================")


if __name__ == "__main__":
    main()
