"""
scripts/core/tiling_utils.py
============================
Core sliding-window and native-DPI patch tiling engine for high-resolution
herbarium voucher scans.
"""

from __future__ import annotations

import json
import logging
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np

from scripts.core.config import CLASS_COLORS_BGR, CLASS_MAP, CLASS_NAMES, DEFAULT_WORKSPACE
from scripts.core.data_structures import (
    ArtifactDetection,
    FilterResult,
    GeometricMetrics,
    InstanceAnnotation,
    SpectralMetrics,
    TextureMetrics,
)
from scripts.core.logger import setup_logging
from scripts.core.tiling_geometry import (
    BackgroundPaperFilter,
    DynamicGeometricReprojector,
    HerbariumAnnotation,
)

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

logger = setup_logging()

__all__ = [
    "NativeDPISlidingWindow",
    "NativeDPIPatchTiler",
    "HerbariumAnnotation",
]


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

        if img_width <= self.tile_size and img_height <= self.tile_size:
            return [(0, 0, img_width, img_height)]

        x_steps: List[int] = []
        cur_x = 0
        while cur_x + self.tile_size < img_width:
            x_steps.append(cur_x)
            cur_x += self.stride
        x_steps.append(max(0, img_width - self.tile_size))
        x_steps = sorted(list(set(x_steps)))

        y_steps: List[int] = []
        cur_y = 0
        while cur_y + self.tile_size < img_height:
            y_steps.append(cur_y)
            cur_y += self.stride
        y_steps.append(max(0, img_height - self.tile_size))
        y_steps = sorted(list(set(y_steps)))

        for y1 in y_steps:
            for x1 in x_steps:
                x2 = min(img_width, x1 + self.tile_size)
                y2 = min(img_height, y1 + self.tile_size)
                windows.append((x1, y1, x2, y2))

        return windows


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

        self.window_generator = NativeDPISlidingWindow(tile_size=self.tile_size, overlap=self.overlap)
        self.reprojector = DynamicGeometricReprojector(min_area_ratio=self.min_area_ratio)
        self.bg_filter = BackgroundPaperFilter(keep_prob=self.bg_keep_prob)

        self.images_dir = self.output_dir / "images"
        self.labels_dir = self.output_dir / "labels"
        self.qc_dir = self.output_dir / "qc_visualizations"

        self.images_dir.mkdir(parents=True, exist_ok=True)
        self.labels_dir.mkdir(parents=True, exist_ok=True)
        if self.visualize:
            self.qc_dir.mkdir(parents=True, exist_ok=True)

        self.metrics: Dict[str, Any] = {
            "total_sheets_processed": 0,
            "total_tiles_generated": 0,
            "positive_tiles_generated": 0,
            "negative_tiles_retained": 0,
            "negative_tiles_discarded": 0,
            "negative_tiles_discarded_due_to_fragments": 0,
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
        label_path: Optional[Union[str, Path]] = None,
        split_name: Optional[str] = None
    ) -> List[Path]:
        """
        Processes a single full-resolution herbarium sheet:
        slices into native-DPI tiles, clips annotations, filters background, and writes files.
        """
        image_path = Path(image_path)
        sheet_stem = image_path.stem

        if self.skip_existing:
            existing_tiles = list(self.images_dir.glob(f"{sheet_stem}_tile_*.jpg"))
            if existing_tiles:
                logger.debug("Skipping already tiled sheet '%s' (%d tiles found)", sheet_stem, len(existing_tiles))
                return existing_tiles

        full_img = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
        if full_img is None:
            logger.error("Failed to load image: %s", image_path)
            return []

        img_height, img_width = full_img.shape[:2]
        logger.info(
            "Processing sheet '%s' | Native Resolution: %dx%d (%.1f MP)",
            sheet_stem, img_width, img_height, (img_width * img_height) / 1e6
        )

        if label_path is None:
            candidate_label = image_path.with_suffix(".txt")
            if candidate_label.exists():
                label_path = candidate_label
            else:
                alt_label = image_path.parent.parent / "labels" / f"{sheet_stem}.txt"
                if alt_label.exists():
                    label_path = alt_label

        annotations: List[HerbariumAnnotation] = []
        if label_path and Path(label_path).exists():
            annotations = self.load_full_sheet_annotations(Path(label_path), img_width, img_height)
            logger.info("Loaded %d ground-truth annotations for sheet '%s'", len(annotations), sheet_stem)

        windows = self.window_generator.generate_windows(img_width, img_height)
        logger.info("Generated %d native-DPI slice windows (Tile Size: %d, Overlap: %.0f%%)",
                    len(windows), self.tile_size, self.overlap * 100)

        generated_tile_paths: List[Path] = []
        self.metrics["total_sheets_processed"] += 1

        for idx, win in enumerate(windows):
            x1, y1, x2, y2 = win
            tile_w = x2 - x1
            tile_h = y2 - y1

            tile_crop = full_img[y1:y2, x1:x2]

            if tile_w < self.tile_size or tile_h < self.tile_size:
                padded = np.full((self.tile_size, self.tile_size, 3), 255, dtype=np.uint8)
                padded[0:tile_h, 0:tile_w] = tile_crop
                tile_crop = padded
                effective_w, effective_h = self.tile_size, self.tile_size
            else:
                effective_w, effective_h = tile_w, tile_h

            tile_annotations, dropped_any_fragment = self.reprojector.reproject_annotations_to_tile(
                annotations, win, effective_w, effective_h
            )

            is_positive = len(tile_annotations) > 0

            if not is_positive:
                if dropped_any_fragment:
                    self.metrics["negative_tiles_discarded_due_to_fragments"] += 1
                    continue
                elif not self.bg_filter.should_keep_empty_tile(tile_crop):
                    self.metrics["negative_tiles_discarded"] += 1
                    continue
                else:
                    self.metrics["negative_tiles_retained"] += 1
            else:
                self.metrics["positive_tiles_generated"] += 1

            tile_name = f"{sheet_stem}_tile_y{y1:05d}_x{x1:05d}"
            out_img_path = self.images_dir / f"{tile_name}.jpg"
            out_lbl_path = self.labels_dir / f"{tile_name}.txt"
            if split_name:
                out_img_dir = self.images_dir / split_name
                out_lbl_dir = self.labels_dir / split_name
                out_img_dir.mkdir(parents=True, exist_ok=True)
                out_lbl_dir.mkdir(parents=True, exist_ok=True)
                out_img_path = out_img_dir / f"{tile_name}.jpg"
                out_lbl_path = out_lbl_dir / f"{tile_name}.txt"

            cv2.imwrite(str(out_img_path), tile_crop, [int(cv2.IMWRITE_JPEG_QUALITY), 95])
            generated_tile_paths.append(out_img_path)
            self.metrics["total_tiles_generated"] += 1

            with open(out_lbl_path, "w", encoding="utf-8") as lf:
                for class_id, norm_coords in tile_annotations:
                    coord_strs = [f"{pt[0]:.6f} {pt[1]:.6f}" for pt in norm_coords]
                    line_str = f"{class_id} {' '.join(coord_strs)}\n"
                    lf.write(line_str)

                    cname = CLASS_NAMES[class_id] if class_id < len(CLASS_NAMES) else f"class_{class_id}"
                    self.metrics["class_instance_counts"][cname] = self.metrics["class_instance_counts"].get(cname, 0) + 1

            present_classes = set(c_id for c_id, _ in tile_annotations)
            for c_id in present_classes:
                cname = CLASS_NAMES[c_id] if c_id < len(CLASS_NAMES) else f"class_{c_id}"
                self.metrics["class_tile_counts"][cname] = self.metrics["class_tile_counts"].get(cname, 0) + 1

            if self.visualize and (is_positive or idx % 10 == 0):
                qc_img = tile_crop.copy()
                for class_id, norm_coords in tile_annotations:
                    pixel_pts = np.array([
                        [int(round(pt[0] * effective_w)), int(round(pt[1] * effective_h))]
                        for pt in norm_coords
                    ], dtype=np.int32)

                    color = CLASS_COLORS_BGR.get(class_id, (0, 255, 0))
                    overlay = qc_img.copy()
                    cv2.fillPoly(overlay, [pixel_pts], color)
                    cv2.addWeighted(overlay, 0.35, qc_img, 0.65, 0, qc_img)
                    cv2.polylines(qc_img, [pixel_pts], isClosed=True, color=color, thickness=2)

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
        """
        summary_path = Path(summary_path)
        summary_path.parent.mkdir(parents=True, exist_ok=True)

        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(self.metrics, f, indent=2)

        logger.info("Exported patch tiling summary metrics to: %s", summary_path)
        return summary_path
