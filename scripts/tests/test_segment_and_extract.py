#!/usr/bin/env python3
"""
===============================================================================
Script: test_segment_and_extract.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 02:
    Direct PointRend segmentation, 2-path botanical leaf extraction (Tier 1
    Pristine and Tier 2 Bilateral Reflection), asynchronous Hough ruler detection,
    centroid spatial clustering, and standardized contour CSV generation.
===============================================================================
"""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

import importlib

PROJECT_ROOT = Path(__file__).resolve().parents[2]

step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
DetectedInstance = step02.DetectedInstance
cluster_plant_individuals = step02.cluster_plant_individuals
compute_geometric_metrics = step02.compute_geometric_metrics
detect_ruler_scale_hough = step02.detect_ruler_scale_hough
export_standardized_contour = step02.export_standardized_contour
extract_tier1_pristine = step02.extract_tier1_pristine
extract_tier2_reflected = step02.extract_tier2_reflected


class TestSegmentAndExtract(unittest.TestCase):
    """Test suite for Step 02 direct segmentation & 2-path extraction pipeline."""

    def setUp(self) -> None:
        self.temp_dir = Path(tempfile.mkdtemp(prefix="test_step02_"))

    def tearDown(self) -> None:
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_geometric_metrics_pristine_ellipse(self) -> None:
        """Verify that a solid elliptic silhouette meets Tier 1 thresholds."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 150), (40, 100), 0, 0, 360, 255, -1)

        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)
        self.assertGreaterEqual(solidity, 0.72, "Solid ellipse must pass Tier 1 solidity threshold (>=0.72).")
        self.assertGreaterEqual(ucs, 0.85, "Solid ellipse must pass Tier 1 UCS threshold (>=0.85).")
        self.assertNotEqual(p_apex, p_base, "Apex and base points must be distinct.")

    def test_tier1_extraction_export(self) -> None:
        """Verify that Tier 1 pristine silhouette is exported cleanly as binary mask."""
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.ellipse(mask, (100, 100), (30, 70), 0, 0, 360, 255, -1)

        saved_path = extract_tier1_pristine(mask, "TEST001", 1, self.temp_dir)
        self.assertTrue(Path(saved_path).exists(), "Exported Tier 1 mask file must exist.")

        read_mask = cv2.imread(saved_path, cv2.IMREAD_GRAYSCALE)
        self.assertIsNotNone(read_mask)
        self.assertEqual(np.count_nonzero(read_mask > 128), np.count_nonzero(mask > 0))

    def test_tier2_bilateral_reflection(self) -> None:
        """Verify that a partially occluded leaf with clean half-blade synthesizes a closed silhouette."""
        # Create an aligned leaf ellipse where lower half is intact, upper half is partially eroded/notched
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 150), (100, 40), 0, 0, 360, 255, -1)  # horizontal ellipse

        # Cut out a chunk in the upper half to simulate occlusion
        mask[:130, 100:180] = 0

        # Prior to reflection, mask has lower solidity/ucs
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)

        saved_path, reflected_mask = extract_tier2_reflected(
            mask, "TEST002", 1, p_apex, p_base, self.temp_dir
        )
        self.assertIsNotNone(saved_path, "Tier 2 bilateral reflection should succeed for intact half.")
        self.assertIsNotNone(reflected_mask)

        # Reflected silhouette should have high solidity and be symmetric
        ref_ucs, ref_solidity, _, _, _ = compute_geometric_metrics(reflected_mask)
        self.assertGreaterEqual(ref_solidity, 0.70, "Synthesized Tier 2 silhouette must have high solidity.")

    def test_hough_ruler_detection(self) -> None:
        """Verify that periodic ruler tick marks yield pixels_per_mm estimation."""
        # Synthesize a herbarium sheet margin containing a 20-tick ruler
        sheet = np.full((1000, 800, 3), 240, dtype=np.uint8)
        ruler_x_start = 50
        tick_pitch = 25  # 25 pixels per mm

        # Draw ruler ticks in the bottom margin
        for i in range(25):
            x = ruler_x_start + i * tick_pitch
            cv2.line(sheet, (x, 850), (x, 900), (20, 20, 20), 2)
        # Add ruler baseline
        cv2.line(sheet, (ruler_x_start, 900), (ruler_x_start + 24 * tick_pitch, 900), (20, 20, 20), 3)

        detected_scale = detect_ruler_scale_hough(sheet)
        self.assertIsNotNone(detected_scale, "Ruler detector should identify tick pitch.")
        self.assertAlmostEqual(detected_scale, float(tick_pitch), delta=4.0)

    def test_hough_ruler_detection_fallback_on_blank(self) -> None:
        """Verify that absent ruler safely returns None without exception."""
        blank_sheet = np.full((500, 500, 3), 240, dtype=np.uint8)
        scale = detect_ruler_scale_hough(blank_sheet)
        self.assertIsNone(scale, "Blank sheet without ticks should gracefully return None.")

    def test_standardized_contour_export(self) -> None:
        """Verify that 2D coordinates CSV is exported with normalized columns."""
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(mask, (100, 100), 50, 255, -1)

        csv_path = export_standardized_contour(mask, "TEST003", 2, self.temp_dir, num_points=120)
        self.assertIsNotNone(csv_path)
        self.assertTrue(Path(csv_path).exists())

        df = pd.read_csv(csv_path)
        self.assertEqual(len(df), 120, "Resampled contour must contain exactly requested points.")
        expected_cols = {"catalogNumber", "leaf_id", "point_index", "x", "y", "x_norm", "y_norm"}
        self.assertTrue(expected_cols.issubset(set(df.columns)))
        self.assertGreaterEqual(df["x_norm"].min(), 0.0)
        self.assertLessEqual(df["x_norm"].max(), 1.0)

    def test_in_sheet_centroid_clustering(self) -> None:
        """Verify that distant leaf candidates receive different plant_individual_id values."""
        dummy_mask = np.ones((10, 10), dtype=np.uint8)
        inst1 = DetectedInstance(
            catalog_number="TEST004", leaf_id=1, bbox=(50, 50, 150, 150),
            mask=dummy_mask, score=0.9, class_id=0
        )
        inst2 = DetectedInstance(
            catalog_number="TEST004", leaf_id=2, bbox=(60, 60, 160, 160),
            mask=dummy_mask, score=0.88, class_id=0
        )
        inst3 = DetectedInstance(
            catalog_number="TEST004", leaf_id=3, bbox=(1800, 1800, 1900, 1900),
            mask=dummy_mask, score=0.85, class_id=0
        )

        clustered = cluster_plant_individuals([inst1, inst2, inst3], sheet_width=2000, sheet_height=2000)
        self.assertEqual(clustered[0].plant_individual_id, clustered[1].plant_individual_id)
        self.assertNotEqual(clustered[0].plant_individual_id, clustered[2].plant_individual_id)


if __name__ == "__main__":
    unittest.main()
