#!/usr/bin/env python3
"""
===============================================================================
Test Suite: test_native_dpi_patch_tiler.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Author: High-Resolution Image Processing Specialist & Senior AI Engineer
Date: August 2026

Description:
    Unit test and validation suite for native_dpi_patch_tiler.py.
    Tests sliding window coverage, geometric polygon clipping, visible area ratio
    filtering, background sub-sampling, and SAHI inference integration.
===============================================================================
"""

import os
import sys
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import cv2
import numpy as np
from shapely.geometry import Polygon, box

# Add project root directory (two levels up from scripts/tests) to python path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(PROJECT_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from scripts.core.tiling_utils import (
    CLASS_NAMES,
    CLASS_MAP,
    HerbariumAnnotation,
    NativeDPISlidingWindow,
    DynamicGeometricReprojector,
    BackgroundPaperFilter,
    NativeDPIPatchTiler
)
from scripts.core.sahi_inference import HerbariumSAHIInference
from scripts.vision.run_dpi_tiler import run_dpi_tiling
from scripts.vision.run_sahi_inference import run_sahi_inference


class TestNativeDPISlidingWindow(unittest.TestCase):
    """
    Test sliding window calculation across various image dimensions.
    """

    def test_window_generation_6000x4000(self):
        """
        Verify grid windows on typical 24 MP herbarium scan (6000x4000 px).
        """
        tiler = NativeDPISlidingWindow(tile_size=1024, overlap=0.20)
        windows = tiler.generate_windows(img_width=6000, img_height=4000)

        self.assertGreater(len(windows), 0)
        # Verify stride calculation: 1024 * 0.8 = 819
        self.assertEqual(tiler.stride, 819)

        # Verify all windows stay within image boundaries
        for x1, y1, x2, y2 in windows:
            self.assertGreaterEqual(x1, 0)
            self.assertGreaterEqual(y1, 0)
            self.assertLessEqual(x2, 6000)
            self.assertLessEqual(y2, 4000)
            self.assertLessEqual(x2 - x1, 1024)
            self.assertLessEqual(y2 - y1, 1024)

        # Verify the rightmost and bottommost windows touch image boundaries exactly
        max_x2 = max(w[2] for w in windows)
        max_y2 = max(w[3] for w in windows)
        self.assertEqual(max_x2, 6000)
        self.assertEqual(max_y2, 4000)

    def test_small_image_single_window(self):
        """
        Verify that images smaller than tile size produce a single window.
        """
        tiler = NativeDPISlidingWindow(tile_size=1024, overlap=0.20)
        windows = tiler.generate_windows(img_width=800, img_height=600)
        self.assertEqual(len(windows), 1)
        self.assertEqual(windows[0], (0, 0, 800, 600))


class TestDynamicGeometricReprojector(unittest.TestCase):
    """
    Test polygon clipping, re-projection, and visible area ratio filtering.
    """

    def setUp(self):
        self.reprojector = DynamicGeometricReprojector(min_area_ratio=0.15)

    def test_fully_contained_polygon(self):
        """
        A polygon fully inside a tile should be preserved with 100% visible area.
        """
        # Create a polygon inside window (100, 100, 1124, 1124)
        poly = Polygon([(200, 200), (300, 200), (300, 300), (200, 300)])
        ann = HerbariumAnnotation(class_id=0, polygon=poly)

        window = (100, 100, 1124, 1124)
        results, dropped = self.reprojector.reproject_annotations_to_tile(
            annotations=[ann], window=window, tile_width=1024, tile_height=1024
        )

        self.assertEqual(len(results), 1)
        self.assertFalse(dropped)
        class_id, norm_pts = results[0]
        self.assertEqual(class_id, 0)
        # Verify local coordinate re-projection:
        # All x values in {100/1024, 200/1024} and y values in {100/1024, 200/1024}
        x_vals = [pt[0] for pt in norm_pts]
        y_vals = [pt[1] for pt in norm_pts]
        self.assertAlmostEqual(min(x_vals), 100.0 / 1024.0, places=4)
        self.assertAlmostEqual(max(x_vals), 200.0 / 1024.0, places=4)
        self.assertAlmostEqual(min(y_vals), 100.0 / 1024.0, places=4)
        self.assertAlmostEqual(max(y_vals), 200.0 / 1024.0, places=4)

    def test_partial_polygon_filtered_below_15_percent(self):
        """
        A polygon with < 15% visible area in tile should be filtered out.
        """
        # Total polygon: 100x100 = 10,000 px^2 spanning (50, 50) to (150, 150)
        poly = Polygon([(50, 50), (150, 50), (150, 150), (50, 150)])
        ann = HerbariumAnnotation(class_id=0, polygon=poly)

        # Window intersecting only from (140, 50) to (150, 150) -> width 10, height 100 -> area 1000 = 10%
        window = (140, 0, 1164, 1024)
        results, dropped = self.reprojector.reproject_annotations_to_tile(
            annotations=[ann], window=window, tile_width=1024, tile_height=1024
        )

        # 10% < 15% threshold -> should be filtered out
        self.assertEqual(len(results), 0)
        self.assertTrue(dropped)

    def test_partial_polygon_retained_above_15_percent(self):
        """
        A polygon with >= 15% visible area in tile should be retained.
        """
        # Total polygon: 100x100 = 10,000 px^2 spanning (50, 50) to (150, 150)
        poly = Polygon([(50, 50), (150, 50), (150, 150), (50, 150)])
        ann = HerbariumAnnotation(class_id=0, polygon=poly)

        # Window intersecting from (120, 50) to (150, 150) -> width 30, height 100 -> area 3000 = 30%
        window = (120, 0, 1144, 1024)
        results, dropped = self.reprojector.reproject_annotations_to_tile(
            annotations=[ann], window=window, tile_width=1024, tile_height=1024
        )

        # 30% >= 15% threshold -> should be preserved
        self.assertEqual(len(results), 1)
        self.assertFalse(dropped)
        class_id, norm_pts = results[0]
        self.assertEqual(class_id, 0)
        self.assertGreaterEqual(len(norm_pts), 3)


class TestBackgroundPaperFilter(unittest.TestCase):
    """
    Test background hard negative paper retention rate.
    """

    def test_sub_sampling_rate(self):
        """
        Verify that 5% keep probability retains approximately 5% over 1000 trials.
        """
        bg_filter = BackgroundPaperFilter(keep_prob=0.05, seed=42)
        dummy_tile = np.full((1024, 1024, 3), 240, dtype=np.uint8)

        trials = 1000
        kept = sum(1 for _ in range(trials) if bg_filter.should_keep_empty_tile(dummy_tile))
        rate = kept / trials

        # Expected rate ~ 0.05 (tolerance [0.03, 0.08])
        self.assertGreaterEqual(rate, 0.03)
        self.assertLessEqual(rate, 0.08)


class TestNativeDPIPatchTilerPipeline(unittest.TestCase):
    """
    Integration test for end-to-end tiling on synthetic high-resolution herbarium sheet.
    """

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="packera_tiler_test_")
        self.temp_path = Path(self.temp_dir)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_end_to_end_tiling(self):
        """
        Create a synthetic 3000x2000 image with 2 annotated polygons and tile it.
        """
        img_w, img_h = 3000, 2000
        synthetic_img = np.full((img_h, img_w, 3), 245, dtype=np.uint8)

        # Draw a synthetic leaf (class 0) and a synthetic ruler (class 6)
        cv2.circle(synthetic_img, (500, 500), 100, (30, 120, 30), -1)
        cv2.rectangle(synthetic_img, (2000, 1200), (2200, 1500), (0, 140, 255), -1)

        img_file = self.temp_path / "test_specimen.jpg"
        cv2.imwrite(str(img_file), synthetic_img)

        # Create YOLO format annotations
        # Leaf: circle approximated polygon (normalized)
        # Ruler: normalized bounding box [class cx cy w h]
        lbl_file = self.temp_path / "test_specimen.txt"
        with open(lbl_file, "w", encoding="utf-8") as f:
            # Class 0: polygon points around (500, 500)
            f.write(f"0 {400/img_w:.6f} {400/img_h:.6f} {600/img_w:.6f} {400/img_h:.6f} {600/img_w:.6f} {600/img_h:.6f} {400/img_w:.6f} {600/img_h:.6f}\n")
            # Class 6: ruler [class 6, cx=(2000+100)/3000, cy=(1200+150)/2000, w=200/3000, h=300/2000]
            f.write(f"6 {2100/img_w:.6f} {1350/img_h:.6f} {200/img_w:.6f} {300/img_h:.6f}\n")

        # Run tiler
        out_dir = self.temp_path / "tiled_output"
        tiler = NativeDPIPatchTiler(
            tile_size=1024,
            overlap=0.20,
            min_area_ratio=0.15,
            bg_keep_prob=0.50,  # higher probability for test visibility
            output_dir=out_dir,
            visualize=True
        )

        generated_tiles = tiler.process_sheet(image_path=img_file, label_path=lbl_file)
        self.assertGreater(len(generated_tiles), 0)

        # Export and verify summary metrics
        summary_path = out_dir / "tiling_summary.json"
        tiler.export_summary(summary_path)

        self.assertTrue(summary_path.exists())
        with open(summary_path, "r", encoding="utf-8") as sf:
            metrics = json.load(sf)

        self.assertEqual(metrics["total_sheets_processed"], 1)
        self.assertGreater(metrics["total_tiles_generated"], 0)
        self.assertGreater(metrics["class_instance_counts"]["basal_leaf_blade"], 0)
        self.assertGreater(metrics["class_instance_counts"]["capitulum"], 0)


class TestSplitPipelineComponents(unittest.TestCase):
    """
    Test the split pipeline runners: run_dpi_tiler, run_sahi_inference, and orchestrator.
    """

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="packera_split_test_")
        self.temp_path = Path(self.temp_dir)

        # Create mock directories
        self.raw_dir = self.temp_path / "raw_vouchers"
        self.labels_dir = self.temp_path / "labels"
        self.tiled_dir = self.temp_path / "tiled_dataset"
        self.outputs_dir = self.temp_path / "outputs"

        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.labels_dir.mkdir(parents=True, exist_ok=True)

        # Create a synthetic image and label
        img_w, img_h = 2000, 2000
        synthetic_img = np.full((img_h, img_w, 3), 240, dtype=np.uint8)
        cv2.circle(synthetic_img, (500, 500), 80, (20, 100, 20), -1)

        img_file = self.raw_dir / "sheet_001.jpg"
        cv2.imwrite(str(img_file), synthetic_img)

        lbl_file = self.labels_dir / "sheet_001.txt"
        with open(lbl_file, "w", encoding="utf-8") as f:
            f.write(f"0 {420/img_w:.6f} {420/img_h:.6f} {580/img_w:.6f} {420/img_h:.6f} {580/img_w:.6f} {580/img_h:.6f} {420/img_w:.6f} {580/img_h:.6f}\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_run_dpi_tiling_execution(self):
        """
        Verify standalone run_dpi_tiling correctly processes vouchers and writes summary.
        """
        summary_file = self.outputs_dir / "tiling_summary.json"
        metrics = run_dpi_tiling(
            input_dir=self.raw_dir,
            labels_dir=self.labels_dir,
            output_dir=self.tiled_dir,
            summary_output=summary_file,
            tile_size=1024,
            overlap=0.20,
            num_workers=2,
            force=True
        )

        self.assertEqual(metrics["total_sheets_processed"], 1)
        self.assertGreater(metrics["total_tiles_generated"], 0)
        self.assertTrue(summary_file.exists())

    def test_run_sahi_inference_empty_dir(self):
        """
        Verify run_sahi_inference handles an empty directory gracefully without crashing.
        """
        empty_dir = self.temp_path / "empty_dir"
        empty_dir.mkdir(parents=True, exist_ok=True)
        results = run_sahi_inference(
            weights="dummy_weights.pt",
            input_dir=empty_dir,
            output_dir=self.outputs_dir,
            summary_output=self.outputs_dir / "sahi_summary.json"
        )
        self.assertEqual(results, [])


if __name__ == "__main__":
    unittest.main()

