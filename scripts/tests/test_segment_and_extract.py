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
SegmentAndExtractPipeline = step02.SegmentAndExtractPipeline
parse_args = step02.parse_args
cluster_plant_individuals = step02.cluster_plant_individuals
compute_geometric_metrics = step02.compute_geometric_metrics
detect_ruler_scale_hough = step02.detect_ruler_scale_hough
export_standardized_contour = step02.export_standardized_contour
extract_tier1_pristine = step02.extract_tier1_pristine
extract_tier2_reflected = step02.extract_tier2_reflected
compute_capitulum_metrics = step02.compute_capitulum_metrics
homologize_contour_starting_point = step02.homologize_contour_starting_point


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

    def test_compute_capitulum_metrics_valid(self) -> None:
        """Verify valid cylindrical capitulum metric extraction and mm scaling."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        # Draw a vertical ellipse representing an upright cylindrical capitulum (H=140, W=100)
        cv2.ellipse(mask, (150, 150), (50, 70), 0, 0, 360, 255, -1)

        metrics = compute_capitulum_metrics(mask, (80, 100, 220, 200), pixels_per_mm=10.0)
        self.assertIsNotNone(metrics)
        self.assertGreaterEqual(metrics["capitulum_aspect_ratio"], 0.7)
        self.assertLessEqual(metrics["capitulum_aspect_ratio"], 2.0)
        self.assertAlmostEqual(metrics["capitulum_aspect_ratio"], 1.4, delta=0.15)
        self.assertAlmostEqual(metrics["involucre_height_mm"], 14.0, delta=1.5)
        self.assertAlmostEqual(metrics["involucre_width_mm"], 10.0, delta=1.5)

    def test_compute_capitulum_metrics_filter_extremes(self) -> None:
        """Verify aberrant non-cylindrical shapes are rejected by 0.7 <= H/W <= 2.0 gating."""
        # 1. Extreme flattened horizontal disc (AR < 0.7)
        mask_flat = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(mask_flat, (20, 90), (180, 110), 255, -1)  # W=160, H=20, AR=0.125
        metrics_flat = compute_capitulum_metrics(mask_flat, (90, 20, 110, 180))
        self.assertIsNone(metrics_flat, "Excessively flat outline must be rejected (AR < 0.7).")

        # 2. Extreme vertical stem fragment (AR > 2.0)
        mask_stem = np.zeros((300, 100), dtype=np.uint8)
        cv2.rectangle(mask_stem, (40, 20), (60, 280), 255, -1)  # W=20, H=260, AR=13.0
        metrics_stem = compute_capitulum_metrics(mask_stem, (20, 40, 280, 60))
        self.assertIsNone(metrics_stem, "Excessively elongated stem fragment must be rejected (AR > 2.0).")

    def test_homologize_contour_starting_point_anchoring_and_orientation(self) -> None:
        """Verify contour start point rolls to closest anchor vertex and enforces clockwise orientation."""
        # 1. Counter-clockwise oriented polygon (positive area in OpenCV)
        # Coordinates: bottom (100, 170), left (30, 100), top (100, 30), right (170, 100)
        pts_ccw = np.array([[100, 170], [30, 100], [100, 30], [170, 100]], dtype=np.float32)
        self.assertGreater(cv2.contourArea(pts_ccw, oriented=True), 0)

        # Petiole anchor near top (100, 25)
        homo_top = homologize_contour_starting_point(pts_ccw, (100, 25))
        self.assertEqual(len(homo_top), len(pts_ccw), "Vertex count must be preserved exactly.")
        np.testing.assert_array_almost_equal(homo_top[0], [100.0, 30.0], err_msg="Index 0 must be closest to anchor.")
        self.assertLessEqual(cv2.contourArea(homo_top, oriented=True), 0, "Oriented area must be negative (CW).")

        # Petiole anchor near bottom (100, 180)
        homo_bottom = homologize_contour_starting_point(pts_ccw, (100, 180))
        self.assertEqual(len(homo_bottom), len(pts_ccw))
        np.testing.assert_array_almost_equal(homo_bottom[0], [100.0, 170.0])
        self.assertLessEqual(cv2.contourArea(homo_bottom, oriented=True), 0)

        # 2. Already clockwise oriented polygon (negative area in OpenCV)
        pts_cw = np.array([[100, 170], [170, 100], [100, 30], [30, 100]], dtype=np.float32)
        self.assertLess(cv2.contourArea(pts_cw, oriented=True), 0)

        homo_cw = homologize_contour_starting_point(pts_cw, (25, 100))
        self.assertEqual(len(homo_cw), len(pts_cw))
        np.testing.assert_array_almost_equal(homo_cw[0], [30.0, 100.0])
        self.assertLessEqual(cv2.contourArea(homo_cw, oriented=True), 0)

    def test_homologize_contour_starting_point_guards(self) -> None:
        """Verify homologize_contour_starting_point gracefully handles empty, degenerate, and edge cases."""
        # Empty array
        empty = homologize_contour_starting_point(np.empty((0, 2)), (100, 100))
        self.assertEqual(len(empty), 0)

        # Single point
        single = homologize_contour_starting_point(np.array([[50, 50]]), (100, 100))
        self.assertEqual(len(single), 1)

        # Collinear points (zero area)
        collinear = np.array([[10, 10], [20, 20], [30, 30]], dtype=np.float32)
        homo_collinear = homologize_contour_starting_point(collinear, (28, 28))
        self.assertEqual(len(homo_collinear), 3)
        np.testing.assert_array_almost_equal(homo_collinear[0], [30.0, 30.0])

        # None anchor
        unchanged = homologize_contour_starting_point(collinear, None)
        np.testing.assert_array_almost_equal(unchanged, collinear)

    def test_export_standardized_contour_homologized_start_point(self) -> None:
        """Verify exported CSV contour begins at the petiole anchor and has clockwise orientation."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        # Vertical ellipse: center (150, 150), axes (40, 100) -> apex ~(150, 50), base ~(150, 250)
        cv2.ellipse(mask, (150, 150), (40, 100), 0, 0, 360, 255, -1)

        p_base = (150, 250)
        csv_path = export_standardized_contour(
            mask, "TEST_HOMO", 1, self.temp_dir, num_points=120, petiole_attachment_pt=p_base
        )
        self.assertIsNotNone(csv_path)
        df = pd.read_csv(csv_path)

        # First vertex should be closest to p_base
        dist_pt0 = np.hypot(df.iloc[0]["x"] - p_base[0], df.iloc[0]["y"] - p_base[1])
        self.assertLessEqual(dist_pt0, 5.0, "Index 0 must be located at the petiole base.")

        # Resampled contour must trace clockwise: oriented area should be negative
        pts = df[["x", "y"]].to_numpy(dtype=np.float32)
        area = cv2.contourArea(pts, oriented=True)
        self.assertLess(area, 0, "Resampled standardized contour must be strictly clockwise.")

    """Test suite for Step 02 inference resumption, artifact skipping, and CLI flags."""

    def setUp(self) -> None:
        self.temp_dir = Path(tempfile.mkdtemp(prefix="test_resumption_"))
        self.vouchers_csv = self.temp_dir / "test_curated_vouchers.csv"
        self.output_dir = self.temp_dir / "pipeline_data" / "data"
        self.contours_dir = self.output_dir / "contours"
        self.masks_dir = self.output_dir / "masks"
        self.tables_dir = self.output_dir / "tables"

        for d in [self.contours_dir, self.masks_dir, self.tables_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Create 2 dummy images and vouchers CSV
        self.img1 = self.temp_dir / "NCU001.jpg"
        self.img2 = self.temp_dir / "NCU002.jpg"
        self.img1.write_bytes(b"\xff\xd8\xff" + b"1" * 500)
        self.img2.write_bytes(b"\xff\xd8\xff" + b"2" * 500)

        df = pd.DataFrame([
            {"catalogNumber": "NCU001", "image_path": str(self.img1)},
            {"catalogNumber": "NCU002", "image_path": str(self.img2)},
        ])
        df.to_csv(self.vouchers_csv, index=False)

    def tearDown(self) -> None:
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_is_voucher_completed_identifies_artifacts(self) -> None:
        """Verify is_voucher_completed identifies existing contour CSV or mask files."""
        pipeline = SegmentAndExtractPipeline(
            vouchers_csv=self.vouchers_csv,
            model_weights=Path("models/dummy.pth"),
            output_dir=self.output_dir,
            device="cpu",
            force=False,
        )
        self.assertFalse(pipeline.is_voucher_completed("NCU001"))

        # Create contour file
        contour_file = self.contours_dir / "NCU001_leaf1.csv"
        contour_file.write_text("catalogNumber,leaf_id,point_index,x,y\nNCU001,1,0,10,20\n")
        self.assertTrue(pipeline.is_voucher_completed("NCU001"))
        self.assertFalse(pipeline.is_voucher_completed("NCU002"))

        # Create mask file for NCU002
        mask_file = self.masks_dir / "NCU002_leaf1.png"
        mask_file.write_bytes(b"dummy_png_bytes")
        self.assertTrue(pipeline.is_voucher_completed("NCU002"))
        pipeline.shutdown()

    def test_run_resumes_and_skips_completed_specimens(self) -> None:
        """Verify run skips inference on completed vouchers and retains prior manifest entries."""
        # Setup pre-existing contour for NCU001
        contour_file = self.contours_dir / "NCU001_leaf1.csv"
        contour_file.write_text("catalogNumber,leaf_id\nNCU001,1\n")

        # Setup pre-existing manifest with entry for NCU001
        manifest_file = self.tables_dir / "extracted_leaves_manifest.csv"
        pd.DataFrame([{
            "catalogNumber": "NCU001",
            "plant_individual_id": 1,
            "leaf_id": 1,
            "assigned_tier": "tier1",
            "ucs_score": 0.90,
            "solidity": 0.85,
            "midrib_angle_deg": 12.0,
            "pixels_per_mm": 10.5,
            "mask_path": "data/masks/NCU001_leaf1.png",
            "contour_path": str(contour_file),
        }]).to_csv(manifest_file, index=False)

        pipeline = SegmentAndExtractPipeline(
            vouchers_csv=self.vouchers_csv,
            model_weights=Path("models/dummy.pth"),
            output_dir=self.output_dir,
            device="cpu",
            force=False,
        )

        mock_extracted = [{
            "catalogNumber": "NCU002",
            "plant_individual_id": 1,
            "leaf_id": 1,
            "assigned_tier": "tier1",
            "ucs_score": 0.88,
            "solidity": 0.82,
            "midrib_angle_deg": 5.0,
            "pixels_per_mm": 10.0,
            "mask_path": "data/masks/NCU002_leaf1.png",
            "contour_path": "data/contours/NCU002_leaf1.csv",
        }]

        from unittest.mock import patch
        with patch.object(pipeline, "process_voucher", return_value=(mock_extracted, None)) as mock_process:
            extracted_df, failed_df = pipeline.run()

            # process_voucher must be called ONLY once (for NCU002, skipping NCU001)
            mock_process.assert_called_once()
            self.assertEqual(mock_process.call_args[0][0], "NCU002")

            # Final manifest must retain records for BOTH NCU001 and NCU002
            self.assertEqual(len(extracted_df), 2)
            self.assertListEqual(list(extracted_df["catalogNumber"]), ["NCU001", "NCU002"])

    def test_run_force_flag_reprocesses_all(self) -> None:
        """Verify that force=True re-runs inference even if artifacts exist."""
        contour_file = self.contours_dir / "NCU001_leaf1.csv"
        contour_file.write_text("catalogNumber,leaf_id\nNCU001,1\n")

        pipeline = SegmentAndExtractPipeline(
            vouchers_csv=self.vouchers_csv,
            model_weights=Path("models/dummy.pth"),
            output_dir=self.output_dir,
            device="cpu",
            force=True,
        )

        from unittest.mock import patch
        with patch.object(pipeline, "process_voucher", return_value=([], None)) as mock_process:
            pipeline.run()
            # process_voucher must be called for BOTH NCU001 and NCU002
            self.assertEqual(mock_process.call_count, 2)

    def test_parse_args_force_and_overwrite(self) -> None:
        """Verify parse_args supports --force and --overwrite flags."""
        import sys
        from unittest.mock import patch
        with patch.object(sys, "argv", ["02_segment_and_extract.py", "--force"]):
            args = parse_args()
            self.assertTrue(args.force)

        with patch.object(sys, "argv", ["02_segment_and_extract.py", "--overwrite"]):
            args = parse_args()
            self.assertTrue(args.force)


class TestRoutingAndManifestSchema(unittest.TestCase):
    """Test suite for the 5-step routing decision tree and diagnostic manifest schema."""

    def setUp(self) -> None:
        self.temp_dir = Path(tempfile.mkdtemp(prefix="test_routing_"))
        self.vouchers_csv = self.temp_dir / "vouchers.csv"
        self.output_dir = self.temp_dir / "data"
        self.img_path = self.temp_dir / "NCU123.jpg"
        self.img_path.write_bytes(b"\xff\xd8\xff" + b"0" * 200)

        pd.DataFrame([
            {"catalogNumber": "NCU123", "image_path": str(self.img_path), "scientificName": "Packera paupercula"}
        ]).to_csv(self.vouchers_csv, index=False)

        self.pipeline = SegmentAndExtractPipeline(
            vouchers_csv=self.vouchers_csv,
            model_weights=Path("models/dummy.pth"),
            output_dir=self.output_dir,
            device="cpu",
            force=True,
        )

    def tearDown(self) -> None:
        self.pipeline.shutdown()
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_manifest_schema_columns(self) -> None:
        """Verify extracted leaf manifest contains required diagnostic classification columns."""
        # Create a mock extracted record
        mock_records = [{
            "catalogNumber": "NCU123",
            "plant_individual_id": 1,
            "leaf_id": 1,
            "tier": "Tier 1",
            "assigned_tier": "tier1",
            "solidity": 0.88,
            "ucs": 0.92,
            "ucs_score": 0.92,
            "is_folded": False,
            "is_dissected": False,
            "reflection_applied": False,
            "midrib_angle_deg": 10.0,
            "pixels_per_mm": 12.0,
            "mask_path": "data/masks/NCU123_leaf1.png",
            "contour_path": "data/contours/NCU123_leaf1.csv",
        }]

        from unittest.mock import patch
        with patch.object(self.pipeline, "process_voucher", return_value=(mock_records, None)):
            extracted_df, _ = self.pipeline.run()

        expected_prefix = [
            "catalogNumber", "leaf_id", "tier", "solidity", "ucs",
            "is_folded", "is_dissected", "reflection_applied", "contour_path"
        ]
        self.assertEqual(list(extracted_df.columns[:9]), expected_prefix)
        self.assertEqual(extracted_df.iloc[0]["tier"], "Tier 1")
        self.assertFalse(bool(extracted_df.iloc[0]["reflection_applied"]))
        self.assertFalse(bool(extracted_df.iloc[0]["is_folded"]))
        self.assertFalse(bool(extracted_df.iloc[0]["is_dissected"]))

    def test_routing_pristine_whole_leaf(self) -> None:
        """Verify pristine whole leaf routes to Tier 1 with reflection_applied=False."""
        # Pristine ellipse (solidity > 0.72, ucs > 0.85, W/L > 0.45)
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 150), (65, 120), 0, 0, 360, 255, -1)

        from unittest.mock import patch
        with patch.object(self.pipeline.model_engine, "predict", return_value=[((30, 85, 270, 215), mask, 0.95, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.pipeline.process_voucher("NCU123", self.img_path)

        self.assertIsNone(failure)
        self.assertEqual(len(records), 1)
        rec = records[0]
        self.assertEqual(rec["tier"], "Tier 1")
        self.assertEqual(rec["assigned_tier"], "tier1")
        self.assertFalse(rec["reflection_applied"])
        self.assertFalse(rec["is_folded"])
        self.assertFalse(rec["is_dissected"])

    def test_routing_folded_leaf_to_tier2(self) -> None:
        """Verify folded leaf along straight chord routes to Tier 2 with reflection_applied=True."""
        # Half-ellipse with straight edge along midrib
        mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(mask, (200, 200), (35, 130), 0, -90, 90, 255, -1)

        from unittest.mock import patch
        with patch.object(self.pipeline.model_engine, "predict", return_value=[((70, 160, 330, 240), mask, 0.92, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.pipeline.process_voucher("NCU123", self.img_path)

        self.assertIsNone(failure)
        self.assertEqual(len(records), 1)
        rec = records[0]
        self.assertEqual(rec["tier"], "Tier 2")
        self.assertEqual(rec["assigned_tier"], "tier2")
        self.assertTrue(rec["reflection_applied"])
        self.assertTrue(rec["is_folded"])
        self.assertFalse(rec["is_dissected"])

    def test_routing_botanical_dissection_to_tier1_dissected(self) -> None:
        """Verify dissected leaf with lyrate taxon routes to Tier 1 (Dissected)."""
        # Create a synthetic pinnatifid/lyrate leaf
        mask = np.zeros((400, 400), dtype=np.uint8)
        cx, cy = 200, 200
        # Terminal lobe
        cv2.ellipse(mask, (cx, cy - 70), (35, 45), 0, 0, 360, 255, -1)
        # Midrib
        cv2.line(mask, (cx, cy - 100), (cx, cy + 100), 255, thickness=16)
        # Bilateral lobes
        for y_off in [-20, 20, 60]:
            cv2.ellipse(mask, (cx - 35, cy + y_off), (35, 12), -15, 0, 360, 255, -1)
            cv2.ellipse(mask, (cx + 35, cy + y_off), (35, 12), 15, 0, 360, 255, -1)

        from unittest.mock import patch
        with patch.object(self.pipeline.model_engine, "predict", return_value=[((100, 160, 300, 240), mask, 0.90, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.pipeline.process_voucher("NCU123", self.img_path, taxon="Packera paupercula")

        self.assertIsNone(failure)
        self.assertEqual(len(records), 1)
        rec = records[0]
        self.assertEqual(rec["tier"], "Tier 1 (Dissected)")
        self.assertEqual(rec["assigned_tier"], "tier1")
        self.assertTrue(rec["is_dissected"])
        self.assertFalse(rec["reflection_applied"])
        self.assertFalse(rec["is_folded"])

    def test_routing_rejection_reasons(self) -> None:
        """Verify rejection reasons FAILED_SOLIDITY, IRREGULAR_FOLD, and UNRESOLVED_CLUMP."""
        # Diagonal fold -> IRREGULAR_FOLD
        diag_mask = np.zeros((400, 400), dtype=np.uint8)
        pts = np.array([[190, 60], [230, 180], [200, 320], [160, 200]], dtype=np.int32)
        cv2.fillPoly(diag_mask, [pts], 255)

        from unittest.mock import patch
        with patch.object(self.pipeline.model_engine, "predict", return_value=[((60, 160, 320, 230), diag_mask, 0.90, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.pipeline.process_voucher("NCU123", self.img_path)

        self.assertEqual(len(records), 0)
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_reason"], "IRREGULAR_FOLD")


if __name__ == "__main__":
    unittest.main()
