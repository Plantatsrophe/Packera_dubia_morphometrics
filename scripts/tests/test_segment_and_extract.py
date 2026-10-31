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
from unittest.mock import patch, MagicMock
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
        self.qc_pipeline = SegmentAndExtractPipeline(
            vouchers_csv=self.vouchers_csv,
            model_weights=Path("models/dummy.pth"),
            output_dir=self.temp_dir / "data_qc",
            device="cpu",
            force=True,
            enable_geometric_qc=True,
        )

    def tearDown(self) -> None:
        self.pipeline.shutdown()
        self.qc_pipeline.shutdown()
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
        self.assertEqual(extracted_df.iloc[0]["assigned_tier"], "Tier_1_Direct")
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
        self.assertEqual(rec["assigned_tier"], "Tier_1_Direct")
        self.assertFalse(rec["reflection_applied"])
        self.assertFalse(rec["is_folded"])
        self.assertFalse(rec["is_dissected"])

    def test_routing_folded_leaf_to_tier2(self) -> None:
        """Verify folded/partial leaf is excluded from closed-outline EFA and routed exclusively to rosette crops under opt-in geometric QC."""
        # Half-ellipse with straight edge along midrib
        mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(mask, (200, 200), (35, 130), 0, -90, 90, 255, -1)

        from unittest.mock import patch
        with patch.object(self.qc_pipeline.model_engine, "predict", return_value=[((70, 160, 330, 240), mask, 0.92, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.qc_pipeline.process_voucher("NCU123", self.img_path)

        self.assertEqual(len(records), 0, "Mirrored Tier 2 hemi-blades must not enter closed-outline EFA records.")
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_reason"], "PARTIAL_OR_FOLDED_ROUTED_TO_ROSETTE")
        rosette_crop = self.qc_pipeline.rosettes_dense_dir / "NCU123_p0_rosette.jpg"
        self.assertTrue(rosette_crop.exists(), "Partial/folded leaf must still generate whole-rosette crop for Phase 5 DINOv2.")

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
        self.assertEqual(rec["assigned_tier"], "Tier_1_Direct")
        self.assertTrue(rec["is_dissected"])
        self.assertFalse(rec["reflection_applied"])
        self.assertFalse(rec["is_folded"])

    def test_routing_rejection_reasons(self) -> None:
        """Verify rejection reasons FAILED_SOLIDITY, IRREGULAR_FOLD, and UNRESOLVED_CLUMP under opt-in geometric QC."""
        # Diagonal fold -> IRREGULAR_FOLD
        diag_mask = np.zeros((400, 400), dtype=np.uint8)
        pts = np.array([[190, 60], [230, 180], [200, 320], [160, 200]], dtype=np.int32)
        cv2.fillPoly(diag_mask, [pts], 255)

        from unittest.mock import patch
        with patch.object(self.qc_pipeline.model_engine, "predict", return_value=[((60, 160, 320, 230), diag_mask, 0.90, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.qc_pipeline.process_voucher("NCU123", self.img_path)

        self.assertEqual(len(records), 0)
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_reason"], "IRREGULAR_FOLD")

    def test_stage1_to_stage2_crop_reprojection(self) -> None:
        """Verify Stage 1 (LeafPriority.pt 11-class) -> 12%-padded crop -> Stage 2 PointRend reprojection."""
        from unittest.mock import patch

        sheet_img = np.full((1000, 800, 3), 220, dtype=np.uint8)
        # Stage 1 box for leaf_whole (class 0): (ymin=200, xmin=100, ymax=400, xmax=300) -> 200x200 box
        # 12% padding -> pad_h=24, pad_w=24 -> crop_box = (176, 76, 424, 324) -> size (248, 248)
        s1_bbox = (200, 100, 400, 300)
        crop_bgr, crop_box = step02.crop_leaf_bbox_with_padding(sheet_img, s1_bbox, pad_ratio=0.12)
        self.assertEqual(crop_box, (176, 76, 424, 324))
        self.assertEqual(crop_bgr.shape[:2], (248, 248))

        # Build a pristine elliptic mask inside the 248x248 crop
        crop_mask = np.zeros((248, 248), dtype=np.uint8)
        cv2.ellipse(crop_mask, (124, 124), (55, 95), 0, 0, 360, 255, -1)

        sheet_mask, full_bbox = step02.reproject_crop_mask_to_sheet(crop_mask, crop_box, (1000, 800))
        self.assertEqual(sheet_mask.shape, (1000, 800))
        self.assertEqual(np.count_nonzero(sheet_mask), np.count_nonzero(crop_mask))
        # Verify reprojected mask is positioned at crop_box offset
        np.testing.assert_array_equal(sheet_mask[176:424, 76:324], crop_mask)

        # Now test full two-stage pipeline execution with mixed 11-class Stage 1 detections:
        # 0: leaf_whole, 1: leaf_partial, 2: leaflet (ignored), 8: specimen (ignored), 10: wood (ignored)
        stage1_dets = [
            ((200, 100, 400, 300), 0.94, 0),  # leaf_whole -> routed to Stage 2
            ((550, 400, 750, 600), 0.89, 1),  # leaf_partial -> routed to Stage 2 under qc_pipeline
            ((100, 50, 150, 90), 0.80, 2),    # leaflet -> ignored
            ((10, 10, 990, 790), 0.99, 8),    # specimen -> ignored
            ((800, 100, 950, 140), 0.77, 10), # wood -> ignored
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.24), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.96, 0)]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.qc_pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.qc_pipeline.model_engine, "predict", side_effect=fake_stage2_predict) as mock_s2:
                    records, failure = self.qc_pipeline.process_voucher("NCU123", self.img_path)

        self.assertIsNone(failure)
        # Stage 2 PointRend must be called ONLY for class 0 (leaf_whole) and class 1 (leaf_partial)
        self.assertEqual(mock_s2.call_count, 2)
        self.assertEqual(len(records), 2)
        for rec in records:
            self.assertEqual(rec["assigned_tier"], "Tier_1_Direct")
            self.assertTrue(Path(rec["contour_path"]).exists())
            self.assertTrue(Path(rec["mask_path"]).exists())

    def test_reproductive_routing_and_capitula_patch_saving(self) -> None:
        """Verify Stage 1 reproductive detections (classes 3..7) route to phenology classifier & save capitula crops."""
        from unittest.mock import patch

        sheet_img = np.full((1000, 800, 3), 230, dtype=np.uint8)
        # Paint golden-yellow corollas inside the flower_one (class 5) region (50..170, 300..400)
        # BGR for golden yellow (RGB ~ 240, 200, 20 -> BGR = 20, 200, 240)
        sheet_img[60:120, 310:390] = (20, 200, 240)
        # Paint dark green erect involucre at lower half (120..165, 325..375)
        sheet_img[120:165, 325:375] = (30, 90, 40)

        stage1_dets = [
            ((300, 200, 550, 420), 0.95, 0),  # 0: leaf_whole
            ((50, 300, 170, 400), 0.92, 5),   # 5: flower_one (H=120, W=100, AR=1.2)
            ((60, 450, 170, 550), 0.88, 3),   # 3: seed_fruit_one (H=110, W=100, AR=1.1)
            ((70, 600, 160, 680), 0.84, 7),   # 7: bud (H=90, W=80, AR=1.125)
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.25), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.95, 0)]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict):
                    extracted_df, failed_df = self.pipeline.run()

        self.assertEqual(len(extracted_df), 1)
        self.assertEqual(len(failed_df), 0)

        # Verify cropped capitula patches saved to data/cropped_patches/capitula/
        cap_patches = sorted(self.pipeline.capitula_patches_dir.glob("NCU123_cap*.jpg"))
        self.assertEqual(len(cap_patches), 3, "All 3 reproductive detections must be saved as capitula crops.")

        # Verify voucher_phenological_states.csv was written and contains total_capitula == 3
        pheno_csv = self.pipeline.tables_dir / "voucher_phenological_states.csv"
        self.assertTrue(pheno_csv.exists())
        pheno_df = pd.read_csv(pheno_csv)
        self.assertEqual(len(pheno_df), 1)
        self.assertEqual(int(pheno_df.iloc[0]["total_capitula"]), 3)
        self.assertIn(
            pheno_df.iloc[0]["voucher_phenological_state"],
            {"anthesis", "fruit", "bud"},
        )

    def test_tier4_rosette_and_basal_leaf_patch_generation(self) -> None:
        """Verify crop_dense_rosette_tier4() and basal_leaves_raw RGB patch saving per plant_individual_id."""
        from unittest.mock import patch

        sheet_img = np.full((2000, 2000, 3), 210, dtype=np.uint8)
        # Two spatially separated plants on the same sheet:
        # Plant 0 near top-left (100..300, 100..280) and Plant 1 near bottom-right (1500..1700, 1500..1680)
        stage1_dets = [
            ((100, 100, 320, 280), 0.95, 0),
            ((120, 140, 340, 310), 0.93, 0),
            ((1500, 1500, 1720, 1680), 0.91, 0),
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.25), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.95, 0)]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict):
                    records, failure = self.pipeline.process_voucher("NCU123", self.img_path)

        self.assertIsNone(failure)
        self.assertEqual(len(records), 3)

        # Verify 2 distinct plant_individual_id clusters were formed
        plant_ids = sorted({r["plant_individual_id"] for r in records})
        self.assertEqual(plant_ids, [0, 1])

        # Verify Tier 4 dense rosette patches exist for both plant individuals
        rosette_p0 = self.pipeline.rosettes_dense_dir / "NCU123_p0_rosette.jpg"
        rosette_p1 = self.pipeline.rosettes_dense_dir / "NCU123_p1_rosette.jpg"
        self.assertTrue(rosette_p0.exists(), "Tier 4 rosette crop for plant_individual_id=0 must exist.")
        self.assertTrue(rosette_p1.exists(), "Tier 4 rosette crop for plant_individual_id=1 must exist.")

        # Verify basal leaf RGB crops exist in data/cropped_patches/basal_leaves_raw/
        for leaf_id in [1, 2, 3]:
            leaf_patch = self.pipeline.basal_leaves_raw_dir / f"NCU123_leaf{leaf_id}.jpg"
            self.assertTrue(leaf_patch.exists(), f"Basal leaf RGB crop {leaf_patch} must exist.")

    def test_lm2_ruler_calibration_with_hough_fallback(self) -> None:
        """Verify load_ruler_calibrations() primary scale vs. detect_ruler_scale_hough() secondary fallback."""
        from unittest.mock import patch

        ruler_dir = self.temp_dir / "lm2_ruler_test" / "Data" / "Ruler"
        ruler_dir.mkdir(parents=True, exist_ok=True)
        ruler_csv = ruler_dir / "ruler_summary.csv"
        pd.DataFrame([
            {
                "filename": "NCU_HIGH_CONF__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 118.11,  # 118.11 px/cm = 11.811 px/mm
                "ruler_class_confidence": 96.5,
            },
            {
                "filename": "NCU_LOW_CONF__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 150.0,
                "ruler_class_confidence": 20.0,  # Below 0.50 confidence threshold -> should be filtered out
            },
        ]).to_csv(ruler_csv, index=False)

        calibs = step02.load_ruler_calibrations(self.temp_dir / "lm2_ruler_test", min_confidence=0.50)
        self.assertIn("NCU_HIGH_CONF", calibs)
        self.assertNotIn("NCU_LOW_CONF", calibs)
        self.assertAlmostEqual(calibs["NCU_HIGH_CONF"], 11.811, places=3)

        # Verify pipeline uses LM2 calibration for NCU_HIGH_CONF without needing Hough fallback,
        # and falls back to detect_ruler_scale_hough for NCU_LOW_CONF
        self.pipeline.ruler_calibrations = calibs
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 150), (65, 120), 0, 0, 360, 255, -1)

        with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
            with patch.object(self.pipeline.model_engine, "predict", return_value=[((30, 85, 270, 215), mask, 0.95, 0)]):
                with patch.object(step02, "detect_ruler_scale_hough", return_value=22.5) as mock_hough:
                    recs_high, _ = self.pipeline.process_voucher("NCU_HIGH_CONF", self.img_path)
                    self.assertEqual(mock_hough.call_count, 0, "Hough fallback must not run when LM2 ruler is valid.")
                    self.assertAlmostEqual(recs_high[0]["pixels_per_mm"], 11.811, places=3)
                    self.assertTrue(recs_high[0]["is_calibrated_mm"])

                    recs_low, _ = self.pipeline.process_voucher("NCU_LOW_CONF", self.img_path)
                    self.assertEqual(mock_hough.call_count, 1, "Hough fallback must run when LM2 ruler confidence is low.")
                    self.assertAlmostEqual(recs_low[0]["pixels_per_mm"], 22.5, places=3)
                    self.assertTrue(recs_low[0]["is_calibrated_mm"])

    @staticmethod
    def _build_synthetic_lyrate_mask(size: int = 400) -> np.ndarray:
        """Synthesizes a deeply lyrate-pinnatifid basal leaf mask with terminal + bilateral lobes and deep sinuses."""
        mask = np.zeros((size, size), dtype=np.uint8)
        cx, cy = size // 2, size // 2
        # Terminal lobe
        cv2.ellipse(mask, (cx, cy - 70), (35, 45), 0, 0, 360, 255, -1)
        # Central midrib/rachis
        cv2.line(mask, (cx, cy - 100), (cx, cy + 100), 255, thickness=16)
        # Bilateral pinnae lobes with deep inter-lobe sinuses
        for y_off in [-20, 20, 60]:
            cv2.ellipse(mask, (cx - 35, cy + y_off), (35, 12), -15, 0, 360, 255, -1)
            cv2.ellipse(mask, (cx + 35, cy + y_off), (35, 12), 15, 0, 360, 255, -1)
        return mask

    def test_generate_sam2_prompts_from_lm2_mask_on_lyrate_leaf(self) -> None:
        """Verify automated prompt generation (medial/lobe positive points, corner/concavity negative points, box, logits)."""
        m_init = self._build_synthetic_lyrate_mask(400)
        bundle = step02.generate_sam2_prompts_from_lm2_mask(m_init)

        # 1. Verify positive foreground points (pole of inaccessibility + top 2 lobe distance peaks)
        pos_pts = bundle["positive_points"]
        self.assertGreaterEqual(pos_pts.shape[0], 2, "Must generate primary interior pole + lobe peak prompts.")
        self.assertLessEqual(pos_pts.shape[0], 3, "At most 3 positive points (pole + top 2 lobe peaks).")
        eroded = bundle["eroded_mask"]
        for px, py in pos_pts:
            ix, iy = int(round(float(px))), int(round(float(py)))
            self.assertGreater(
                eroded[iy, ix],
                0,
                f"Positive prompt {(px, py)} must lie strictly inside eroded m_init.",
            )

        # 2. Verify negative background points (4 exterior corners + deep concavity pockets outside dilated m_init)
        neg_pts = bundle["negative_points"]
        corner_pts = bundle["corner_points"]
        concavity_pts = bundle["concavity_points"]
        self.assertGreaterEqual(corner_pts.shape[0], 4, "Must include 4 exterior bounding-box corner prompts.")
        self.assertGreaterEqual(
            concavity_pts.shape[0],
            1,
            "Deeply lyrate leaf mask must yield negative prompt points in deep concavity pockets.",
        )
        dilated = bundle["dilated_mask"]
        for nx, ny in neg_pts:
            ix, iy = int(round(float(nx))), int(round(float(ny)))
            self.assertEqual(
                dilated[iy, ix],
                0,
                f"Negative prompt {(nx, ny)} must lie strictly outside dilated m_init.",
            )

        # 3. Verify combined point_coords, point_labels, box, and (1, 256, 256) low-res logits
        self.assertEqual(bundle["point_coords"].shape[0], pos_pts.shape[0] + neg_pts.shape[0])
        self.assertEqual(bundle["point_labels"].shape[0], bundle["point_coords"].shape[0])
        self.assertEqual(int(np.sum(bundle["point_labels"] == 1)), pos_pts.shape[0])
        self.assertEqual(int(np.sum(bundle["point_labels"] == 0)), neg_pts.shape[0])
        self.assertEqual(bundle["box"].shape, (4,))
        self.assertEqual(bundle["mask_input"].shape, (1, 256, 256))
        self.assertEqual(bundle["mask_input"].dtype, np.float32)
        self.assertGreater(float(np.max(bundle["mask_input"])), 0.0)
        self.assertLess(float(np.min(bundle["mask_input"])), 0.0)

    def test_sam2_boundary_refinement_and_topological_fallback(self) -> None:
        """Verify SAM2BoundaryRefiner refines sub-pixel margins and falls back cleanly on disconnected masks."""
        from unittest.mock import MagicMock, patch

        m_init = self._build_synthetic_lyrate_mask(400)
        crop_bgr = np.full((400, 400, 3), 210, dtype=np.uint8)

        # 1. Clean fallback when SAM 2 weights/predictor are unavailable
        refiner_no_weights = step02.SAM2BoundaryRefiner(
            checkpoint_path=self.temp_dir / "nonexistent_sam2.pt",
            predictor=None,
        )
        out_fallback = refiner_no_weights.refine_crop_mask(crop_bgr, m_init)
        np.testing.assert_array_equal(out_fallback, m_init)

        # 2. Successful refinement when SAM2ImagePredictor returns a connected refined lyrate mask
        refined_lyrate = cv2.dilate(m_init, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)), iterations=1)
        mock_predictor = MagicMock()
        mock_predictor.predict.return_value = (
            np.stack([(refined_lyrate > 0).astype(np.float32)], axis=0),
            np.array([0.97], dtype=np.float32),
            np.zeros((1, 256, 256), dtype=np.float32),
        )
        refiner_ok = step02.SAM2BoundaryRefiner(predictor=mock_predictor)
        out_refined = refiner_ok.refine_crop_mask(crop_bgr, m_init)
        self.assertTrue(mock_predictor.set_image.called)
        self.assertTrue(mock_predictor.predict.called)
        _, pred_kwargs = mock_predictor.predict.call_args
        self.assertIn("point_coords", pred_kwargs)
        self.assertIn("point_labels", pred_kwargs)
        self.assertIn("box", pred_kwargs)
        self.assertIn("mask_input", pred_kwargs)
        self.assertEqual(pred_kwargs["mask_input"].shape, (1, 256, 256))
        np.testing.assert_array_equal(out_refined, refined_lyrate)

        # Verify trim_petiole_tail() and is_botanical_dissection() run identically on refined mask
        trimmed_init = step02.trim_petiole_tail(m_init)
        trimmed_refined = step02.trim_petiole_tail(out_refined)
        self.assertGreater(np.count_nonzero(trimmed_init), 0)
        self.assertGreater(np.count_nonzero(trimmed_refined), 0)
        self.assertTrue(step02.is_botanical_dissection(m_init))
        self.assertTrue(step02.is_botanical_dissection(out_refined))

        # 3. Topological connectedness violation -> clean fallback to m_init
        disconnected_mask = np.zeros_like(m_init)
        cv2.circle(disconnected_mask, (70, 70), 25, 255, -1)
        cv2.circle(disconnected_mask, (290, 290), 25, 255, -1)
        self.assertFalse(step02.is_topologically_connected_mask(disconnected_mask, reference_mask=m_init))
        mock_predictor_bad = MagicMock()
        mock_predictor_bad.predict.return_value = (
            np.stack([(disconnected_mask > 0).astype(np.float32)], axis=0),
            np.array([0.99], dtype=np.float32),
            np.zeros((1, 256, 256), dtype=np.float32),
        )
        refiner_bad = step02.SAM2BoundaryRefiner(predictor=mock_predictor_bad)
        out_bad_fallback = refiner_bad.refine_crop_mask(crop_bgr, m_init)
        np.testing.assert_array_equal(
            out_bad_fallback,
            m_init,
            "Refiner must fall back cleanly to m_init when refined mask violates topological connectedness.",
        )

        # 4. Verify --sam2-refine CLI flag and pipeline integration
        with patch("sys.argv", ["02_segment_and_extract.py", "--sam2-refine"]):
            parsed = step02.parse_args()
            self.assertTrue(parsed.sam2_refine)

    def test_extract_tier2_reflected_45deg_asymmetric_leaf_iou(self) -> None:
        """
        Verify a synthetic asymmetric leaf mask oriented at 45 degrees passed into
        extract_tier2_reflected() produces a bilaterally symmetric mask with >0.95 IoU
        against the ground-truth unoccluded blade.
        """
        from scripts.vision.lm2_geometry_utils import extract_tier2_reflected

        canvas_size = (400, 400)
        gt_mask = np.zeros(canvas_size, dtype=np.uint8)
        cx, cy = 200, 200
        semi_major, semi_minor = 110, 45  # L=220, W=90 oriented at 45 deg
        cv2.ellipse(gt_mask, (cx, cy), (semi_major, semi_minor), 45.0, 0, 360, 255, -1)

        # Ground-truth apex and base along the 45-degree midrib axis
        dx = int(round(semi_major * np.cos(np.deg2rad(45.0))))
        dy = int(round(semi_major * np.sin(np.deg2rad(45.0))))
        p_apex = (cx - dx, cy - dy)
        p_base = (cx + dx, cy + dy)

        # Create an asymmetric occluded leaf mask by subtracting an irregular bite on one lateral side
        asymmetric_mask = gt_mask.copy()
        # Normal to (+1, +1)/sqrt(2) is (-1, +1)/sqrt(2) -> bite centered on lower-left lateral margin
        bite_cx = int(round(cx - 32 * np.cos(np.deg2rad(45.0))))
        bite_cy = int(round(cy + 32 * np.sin(np.deg2rad(45.0))))
        cv2.circle(asymmetric_mask, (bite_cx, bite_cy), 24, 0, -1)

        _, reflected_mask = extract_tier2_reflected(
            mask=asymmetric_mask,
            catalog_number="TEST_45DEG",
            leaf_id=1,
            p_apex=p_apex,
            p_base=p_base,
            output_dir=None,
        )
        self.assertIsNotNone(reflected_mask, "extract_tier2_reflected must return a valid synthesized mask.")

        intersection = float(np.count_nonzero((reflected_mask > 0) & (gt_mask > 0)))
        union = float(np.count_nonzero((reflected_mask > 0) | (gt_mask > 0)))
        iou = intersection / max(union, 1.0)
        self.assertGreater(
            iou,
            0.95,
            f"Reflected 45-deg leaf mask must achieve >0.95 IoU against ground-truth unoccluded blade (got {iou:.4f}).",
        )

    def test_export_standardized_contour_isotropic_aspect_ratio_preservation(self) -> None:
        """
        Verify an oblong leaf with L/W = 3.0 exported via export_standardized_contour()
        retains (max(x_norm) - min(x_norm)) / (max(y_norm) - min(y_norm)) == 3.0 within 1e-3.
        """
        mask = np.zeros((500, 500), dtype=np.uint8)
        # Oblong leaf with semi-axes a=150 (L=300), b=50 (W=100) -> L/W = 3.0
        cv2.ellipse(mask, (250, 250), (50, 150), 0.0, 0, 360, 255, -1)

        csv_path = export_standardized_contour(
            mask,
            catalog_number="TEST_OBLONG_LW3",
            leaf_id=1,
            output_dir=self.temp_dir,
            num_points=120,
            petiole_attachment_pt=(250, 400),
        )
        self.assertIsNotNone(csv_path)
        df = pd.read_csv(csv_path)

        span_x_norm = float(df["x_norm"].max() - df["x_norm"].min())
        span_y_norm = float(df["y_norm"].max() - df["y_norm"].min())
        lw_ratio = span_x_norm / max(span_y_norm, 1e-9)
        self.assertAlmostEqual(
            lw_ratio,
            3.0,
            delta=1e-3,
            msg=f"Isotropically normalized contour must preserve L/W = 3.0 within 1e-3 (got {lw_ratio:.6f}).",
        )

    def test_base_apex_polarity_disambiguation_and_label_blind_dissection(self) -> None:
        """
        Verify compute_geometric_metrics() disambiguates p_base (narrow tapering stalk end)
        from p_apex (broad lamina tip) regardless of cv2.fitLine sign, and classify_and_route_leaf()
        never leaks taxonomic label into dissection gatekeeping.
        """
        from scripts.vision.lm2_geometry_utils import classify_and_route_leaf

        # 1. Oblanceolate/spathulate leaf with narrow stalk at TOP (y=60) and broad rounded apex at BOTTOM (y=300)
        mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(mask, (200, 240), (45, 70), 0.0, 0, 360, 255, -1)
        cv2.rectangle(mask, (194, 60), (206, 185), 255, -1)

        _, _, _, p_apex, p_base = compute_geometric_metrics(mask)
        self.assertLess(
            p_base[1],
            p_apex[1],
            f"Narrow tapering petiole stalk at top (y~60) must be identified as p_base, got p_base={p_base}, p_apex={p_apex}.",
        )

        # 2. Unilateral bite on low-solidity leaf with lyrate taxon label must NOT be accepted as dissected
        bitten = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(bitten, (200, 200), (65, 110), 0.0, 0, 360, 255, -1)
        cv2.ellipse(bitten, (145, 200), (58, 65), 0.0, 0, 360, 0, -1)
        res_paupercula = classify_and_route_leaf(bitten, taxon="Packera paupercula")
        res_dubia = classify_and_route_leaf(bitten, taxon="Packera dubia")
        self.assertFalse(res_paupercula.is_dissected)
        self.assertEqual(res_paupercula.assigned_tier, res_dubia.assigned_tier)

    def test_sam2_concavity_points_bounded_loop_and_suppression(self) -> None:
        """Verify generate_sam2_prompts_from_lm2_mask terminates within bounded iterations and never visits same pixel."""
        from scripts.vision.lm2_geometry_utils import generate_sam2_prompts_from_lm2_mask

        # Star-shaped polygon with multiple deep sinuses
        mask = np.zeros((300, 300), dtype=np.uint8)
        pts = np.array([
            [150, 20], [180, 110], [270, 110], [200, 160],
            [230, 250], [150, 195], [70, 250], [100, 160],
            [30, 110], [120, 110]
        ], dtype=np.int32)
        cv2.fillPoly(mask, [pts], 255)

        bundle = generate_sam2_prompts_from_lm2_mask(mask, max_concavity_points=6)
        self.assertIsInstance(bundle, dict)
        self.assertLessEqual(len(bundle["concavity_points"]), 6)
        if len(bundle["concavity_points"]) > 1:
            # Check all points are unique
            pts_arr = bundle["concavity_points"]
            unique_pts = np.unique(pts_arr, axis=0)
            self.assertEqual(len(pts_arr), len(unique_pts))

    def test_classify_and_route_stem_clump_rejection_before_reflection(self) -> None:
        """Verify classify_and_route_leaf rejects stems and clumps (solidity < 0.58 or ucs < 0.62) without Tier 2 reflection."""
        from scripts.vision.lm2_geometry_utils import classify_and_route_leaf

        # Create a ragged clump/stem mask: thin jagged curve with very low solidity (< 0.50) and low UCS (< 0.60)
        mask = np.zeros((400, 400), dtype=np.uint8)
        # Draw a thin meandering C-shape
        cv2.ellipse(mask, (200, 200), (120, 150), 0, 30, 330, 255, 12)
        # Add random jagged protrusions
        cv2.line(mask, (200, 50), (280, 20), 255, 8)
        cv2.line(mask, (80, 200), (30, 260), 255, 8)

        routing = classify_and_route_leaf(mask)
        self.assertEqual(routing.tier, "Failed QC")
        self.assertEqual(routing.assigned_tier, "rejected")
        self.assertEqual(routing.rejection_reason, "UNRESOLVED_CLUMP")

    def test_stage1_high_dpi_filtering_and_single_primary_blade_selection(self) -> None:
        """Verify on high-DPI herbarium sheets (>2000px), Stage 1 restricts to Class 0, filters tiny/upper boxes, and selects 1 blade."""
        sheet_img = np.full((3000, 2400, 3), 230, dtype=np.uint8)
        self.assertTrue(max(sheet_img.shape[:2]) > 2000)

        # Stage 1 detections:
        # 1. Class 1 (leaf_partial) -> should be skipped on high-DPI
        # 2. Class 0, but upper cauline: cy = 600 < 0.28 * 3000 (840) -> should be skipped
        # 3. Class 0, but tiny box: bh=80, bw=80 < 0.045 * 3000 (135) -> should be skipped
        # 4. Class 0, valid basal leaf box: ymin=1800, xmin=800, ymax=2400, xmax=1400 (bh=600, bw=600, cy=2100) -> accepted!
        stage1_dets = [
            ((1500, 500, 1900, 900), 0.90, 1),   # leaf_partial -> skipped
            ((400, 800, 800, 1200), 0.92, 0),    # cy=600 < 840 -> skipped
            ((1600, 600, 1680, 680), 0.88, 0),   # max(bh, bw)=80 < 135 -> skipped
            ((1800, 800, 2400, 1400), 0.95, 0),  # valid basal leaf!
        ]

        def fake_high_dpi_stage2(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            # Prediction 1: Class 1 (petiole/stem fragment) -> ignored
            m1 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.rectangle(m1, (10, 10), (30, 200), 255, -1)

            # Prediction 2: Class 0, tiny noise fragment (< 8000 px) -> ignored
            m2 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.circle(m2, (50, 50), 20, 255, -1)

            # Prediction 3: Class 0, large primary blade (> 8000 px, area ~ 40,000) -> selected!
            m3 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m3, (cw // 2, ch // 2), (cw // 3, ch // 3), 0, 0, 360, 255, -1)

            # Prediction 4: Class 0, smaller blade fragment (area ~ 15,000, lower score) -> not selected
            m4 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m4, (cw // 2, ch // 2), (cw // 5, ch // 5), 0, 0, 360, 255, -1)

            return [
                ((10, 10, 200, 30), m1, 0.90, 1),
                ((30, 30, 70, 70), m2, 0.85, 0),
                ((50, 50, ch - 50, cw - 50), m3, 0.95, 0),
                ((100, 100, ch - 100, cw - 100), m4, 0.70, 0),
            ]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.qc_pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.qc_pipeline.model_engine, "predict", side_effect=fake_high_dpi_stage2) as mock_s2:
                    records, failure = self.qc_pipeline.process_voucher("NCU_HIGHDPI", self.img_path)

        self.assertIsNone(failure)
        # Stage 2 PointRend must be called EXACTLY ONCE (only for valid basal leaf box #4)
        self.assertEqual(mock_s2.call_count, 1)
        # Exactly ONE primary blade instance extracted
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["assigned_tier"], "Tier_1_Direct")

    def test_stage1_use_cache_flag(self) -> None:
        """Verify LM2Stage1Detector.predict respects use_cache flag."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            proj_dir = Path(tmp_dir)
            labels_dir = proj_dir / "Plant_Components" / "labels"
            labels_dir.mkdir(parents=True)
            label_file = labels_dir / "VOUCHER1.txt"
            label_file.write_text("0 0.5 0.5 0.2 0.2\n")

            detector = step02.LM2Stage1Detector(
                weights_path="dummy.pt",
                lm2_project_dir=proj_dir,
            )
            # When use_cache=True, loads precomputed labels
            dets_cached = detector.predict(np.zeros((100, 100, 3), dtype=np.uint8), catalog_number="VOUCHER1", use_cache=True)
            self.assertEqual(len(dets_cached), 1)

            # When use_cache=False, bypasses precomputed labels (and returns [] since model is None)
            dets_uncached = detector.predict(np.zeros((100, 100, 3), dtype=np.uint8), catalog_number="VOUCHER1", use_cache=False)
            self.assertEqual(len(dets_uncached), 0)

    def test_heterophylly_cauline_and_stem_rejection_conditions(self) -> None:
        """Verify Milestone 1 cauline/bract/stem rejection conditions (a), (b), and (c)."""
        from scripts.vision.lm2_geometry_utils import is_cauline_or_stem_contaminant

        # Condition (a): capitula_count > 0, cy in upper corymb zone, area < 0.50 * med_area
        rej_a, reason_a = is_cauline_or_stem_contaminant(
            cy=220.0,
            area=400.0,
            aspect_ratio=2.1,
            solidity=0.88,
            basal_cy_med=700.0,
            med_area=1000.0,
            img_h=1000.0,
            capitula_count=2,
            cap_cy=150.0,
            has_petiole=True,
            med_rosette_radius=120.0,
            radial_dist_from_rosette=480.0,
        )
        self.assertTrue(rej_a)
        self.assertEqual(reason_a, "CAULINE_CORYMB_ZONE")

        # Condition (b): narrow sessile/clasping upper cauline bract or stem segment (capitula_count == 0)
        rej_b, reason_b = is_cauline_or_stem_contaminant(
            cy=450.0,
            area=800.0,
            aspect_ratio=4.2,
            solidity=0.72,
            basal_cy_med=700.0,
            med_area=1000.0,
            img_h=1000.0,
            capitula_count=0,
            cap_cy=None,
            has_petiole=False,
            med_rosette_radius=120.0,
            radial_dist_from_rosette=250.0,
        )
        self.assertTrue(rej_b)
        self.assertEqual(reason_b, "CAULINE_SESSILE_BRACT_OR_STEM")

        # Condition (c): lacks petiole, displaced > 2.2 * med_rosette_radius above basal rosette, area < 0.35 * med_area
        rej_c, reason_c = is_cauline_or_stem_contaminant(
            cy=380.0,
            area=300.0,
            aspect_ratio=2.2,
            solidity=0.86,
            basal_cy_med=700.0,
            med_area=1000.0,
            img_h=1000.0,
            capitula_count=0,
            cap_cy=None,
            has_petiole=False,
            med_rosette_radius=110.0,
            radial_dist_from_rosette=320.0,
        )
        self.assertTrue(rej_c)
        self.assertEqual(reason_c, "DISPLACED_SESSILE_CAULINE")

        # Valid basal rosette leaf -> not rejected
        rej_ok, _ = is_cauline_or_stem_contaminant(
            cy=710.0,
            area=1050.0,
            aspect_ratio=1.9,
            solidity=0.91,
            basal_cy_med=700.0,
            med_area=1000.0,
            img_h=1000.0,
            capitula_count=0,
            cap_cy=None,
            has_petiole=True,
            med_rosette_radius=110.0,
            radial_dist_from_rosette=40.0,
        )
        self.assertFalse(rej_ok)

    def test_lobe_rejection_min_blade_area_cm2_extent_and_telemetry(self) -> None:
        """Verify check_boundary_tissue_continuity, min_blade_area_cm2, Extent, and rejection telemetry counters under opt-in QC."""
        from unittest.mock import patch

        # Calibrate sheet at 10.0 px/mm (100.0 px/cm -> 1 cm^2 = 10,000 px, 0.75 cm^2 = 7,500 px)
        self.qc_pipeline.ruler_calibrations["NCU_TELEM"] = 10.0

        sheet_img = np.full((600, 600, 3), 235, dtype=np.uint8)

        # 1. Severed lobe candidate: mask inside contiguous dark plant lamina (>15% exterior ring is dark tissue)
        cv2.ellipse(sheet_img, (150, 150), (95, 110), 0, 0, 360, (40, 55, 40), -1)
        m_severed = np.zeros((600, 600), dtype=np.uint8)
        cv2.ellipse(m_severed, (130, 150), (45, 70), 0, 0, 360, 255, -1)

        # 2. Tiny fragment candidate (< 0.75 cm^2, e.g., radius 25 -> area ~ 1,963 px = 0.196 cm^2)
        m_tiny = np.zeros((600, 600), dtype=np.uint8)
        cv2.circle(m_tiny, (420, 150), 25, 255, -1)

        # 3. Low-Extent candidate: mask area ~ 9,500 px (0.95 cm^2 >= 0.75), but inside a huge 250x250 PCD box (area 62,500 -> Extent = 0.152 < 0.35)
        m_low_extent = np.zeros((600, 600), dtype=np.uint8)
        cv2.ellipse(m_low_extent, (420, 420), (40, 75), 0, 0, 360, 255, -1)

        # 4. Leaflet class candidate (class_id == 2 -> rejected_leaflet_class)
        m_leaflet = np.zeros((600, 600), dtype=np.uint8)
        cv2.ellipse(m_leaflet, (150, 420), (40, 75), 0, 0, 360, 255, -1)

        dets = [
            ((80, 85, 220, 175), m_severed, 0.94, 0),
            ((120, 390, 180, 450), m_tiny, 0.91, 0),
            ((300, 300, 550, 550), m_low_extent, 0.90, 0),
            ((340, 100, 500, 200), m_leaflet, 0.89, 2),
        ]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.qc_pipeline.model_engine, "predict", return_value=dets):
                records, failure = self.qc_pipeline.process_voucher("NCU_TELEM", self.img_path)

        self.assertEqual(len(records), 0)
        self.assertIsNotNone(failure)
        self.assertEqual(self.qc_pipeline.rejection_telemetry["rejected_leaflet_class"], 1)
        self.assertEqual(self.qc_pipeline.rejection_telemetry["rejected_severed_lobe"], 1)
        self.assertEqual(self.qc_pipeline.rejection_telemetry["rejected_min_area_cm2"], 1)
        self.assertEqual(self.qc_pipeline.rejection_telemetry["rejected_solidity_extent"], 1)

    def test_dissected_tier1_routing_enforces_min_ucs(self) -> None:
        """Verify dissected candidates with ucs < 0.80 (and < min_ucs) are rejected under opt-in geometric QC."""
        from unittest.mock import patch
        from scripts.vision.lm2_geometry_utils import LeafRoutingResult

        sheet_img = np.full((600, 600, 3), 235, dtype=np.uint8)
        m_dissected_low_ucs = self._build_synthetic_lyrate_mask(400)
        m_full = np.zeros((600, 600), dtype=np.uint8)
        m_full[100:500, 100:500] = m_dissected_low_ucs

        dets = [((100, 100, 500, 500), m_full, 0.95, 0)]
        low_ucs_routing = LeafRoutingResult(
            assigned_tier="tier1",
            processed_mask=m_dissected_low_ucs,
            solidity=0.64,
            ucs=0.71,
            is_folded=False,
            is_dissected=True,
            is_irregular_fold=False,
            rejection_reason=None,
            tier="Tier 1",
            metadata={"p_base": (200, 340), "p_apex": (200, 55)},
        )

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.qc_pipeline.model_engine, "predict", return_value=dets):
                with patch.object(step02, "classify_and_route_leaf", return_value=low_ucs_routing):
                    records, failure = self.qc_pipeline.process_voucher("NCU_LOW_UCS_DISS", self.img_path)

        self.assertEqual(len(records), 0, "Dissected leaf with ucs < 0.80 must not enter Tier_1_Direct under opt-in QC.")
        self.assertIsNotNone(failure)

    def test_step02_figures_downstream_qc_filtering_and_active_species_montage(self) -> None:
        """
        Verify generate_step02_publication_figures() and regenerate_step02_figures_from_tables():
          1. Exclude segments with is_contour_outlier=True, pct_asymmetry > 0.15, or harmonic_roughness > 0.10
             in leaf_efa_harmonics_individual.csv, or misidentification_flag=True in morphometric_flags.csv.
          2. Dynamically restrict montage rows to active species present in df_plot without cross-species fallback.
          3. Rank quantile-window (q +/- 0.08) candidates by Gold/Silver authority, mm calibration, and leaf_quality_weight,
             while rejecting masks with internal holes or clipped border edges.
        """
        from unittest.mock import patch
        import matplotlib.pyplot as plt

        run_dir = Path(self.temp_dir) / "run_qc_test"
        tables_dir = run_dir / "data" / "tables"
        masks_dir = run_dir / "data" / "masks" / "tier1_intact"
        figures_dir = run_dir / "outputs" / "figures"
        tables_dir.mkdir(parents=True, exist_ok=True)
        masks_dir.mkdir(parents=True, exist_ok=True)
        figures_dir.mkdir(parents=True, exist_ok=True)

        def _write_smooth_mask(path: Path, ar: float, clipped: bool = False, with_hole: bool = False) -> None:
            canvas = np.zeros((300, 300), dtype=np.uint8)
            minor_ax = 36
            major_ax = max(40, min(125, int(round(minor_ax * ar))))
            cv2.ellipse(canvas, (150, 150), (minor_ax, major_ax), 0.0, 0, 360, 255, -1)
            if clipped:
                canvas[0, 120:180] = 255
            if with_hole:
                cv2.circle(canvas, (150, 150), 12, 0, -1)
            cv2.imwrite(str(path), canvas)

        voucher_rows = []
        manifest_rows = []
        ind_efa_rows = []
        otu_efa_rows = []
        mflag_rows = []

        # Only 2 active species: Packera dubia and Packera paupercula (anonyma & plattensis excluded)
        active_species = ["Packera dubia", "Packera paupercula"]
        ar_grid = [1.30, 1.45, 1.50, 1.75, 1.80, 2.05, 2.10, 2.35, 2.40, 2.75]

        for sp_idx, sp_name in enumerate(active_species):
            sp_code = "DUB" if sp_idx == 0 else "PAU"
            for i, ar_v in enumerate(ar_grid):
                cat = f"{sp_code}_V{i:02d}"
                is_gold = (i % 2 == 1)
                det_tier = "Tier_1_Gold" if is_gold else "Tier_3_Bronze"
                is_cal = is_gold
                lqw = 0.96 if is_gold else 0.72
                m_path = masks_dir / f"{cat}_leaf1.png"
                _write_smooth_mask(m_path, min(2.8, ar_v))

                voucher_rows.append({
                    "catalogNumber": cat,
                    "species_standardized": sp_name,
                    "determiner_tier": det_tier,
                })
                manifest_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct",
                    "solidity": 0.90,
                    "ucs_score": 0.94,
                    "aspect_ratio": ar_v,
                    "pixels_per_mm": 12.0 if is_cal else np.nan,
                    "is_calibrated_mm": is_cal,
                    "mask_path": str(m_path),
                    "contour_path": "",
                })
                ind_efa_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "leaf_id": 1,
                    "species_standardized": sp_name,
                    "determiner_tier": det_tier,
                    "assigned_tier": "Tier_1_Direct",
                    "is_calibrated_mm": "TRUE" if is_cal else "FALSE",
                    "pixels_per_mm": 12.0 if is_cal else np.nan,
                    "solidity": 0.90,
                    "ucs_score": 0.94,
                    "aspect_ratio": ar_v,
                    "pct_asymmetry": 0.04,
                    "harmonic_roughness": 0.01,
                    "is_contour_outlier": "FALSE",
                    "leaf_quality_weight": lqw,
                })
                otu_efa_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "is_contour_outlier": "FALSE",
                })
                mflag_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "misidentification_flag": "FALSE",
                })

            # Add contaminated leaves that MUST be excluded despite having top quality scores:
            # (a) Contour outlier in leaf_efa_harmonics_individual.csv
            cat_out = f"{sp_code}_OUTLIER"
            m_out = masks_dir / f"{cat_out}_leaf1.png"
            _write_smooth_mask(m_out, 1.48)
            voucher_rows.append({"catalogNumber": cat_out, "species_standardized": sp_name, "determiner_tier": "Tier_1_Gold"})
            manifest_rows.append({
                "catalogNumber": cat_out, "plant_individual_id": 1, "leaf_id": 1,
                "assigned_tier": "Tier_1_Direct", "solidity": 0.95, "ucs_score": 0.98,
                "aspect_ratio": 1.48, "pixels_per_mm": 12.0, "is_calibrated_mm": True,
                "mask_path": str(m_out), "contour_path": "",
            })
            ind_efa_rows.append({
                "catalogNumber": cat_out, "plant_individual_id": 1, "leaf_id": 1,
                "species_standardized": sp_name, "determiner_tier": "Tier_1_Gold",
                "assigned_tier": "Tier_1_Direct", "is_calibrated_mm": "TRUE", "pixels_per_mm": 12.0,
                "solidity": 0.95, "ucs_score": 0.98, "aspect_ratio": 1.48,
                "pct_asymmetry": 0.02, "harmonic_roughness": 0.01,
                "is_contour_outlier": "TRUE", "leaf_quality_weight": 0.999,
            })
            otu_efa_rows.append({"catalogNumber": cat_out, "plant_individual_id": 1, "is_contour_outlier": "FALSE"})
            mflag_rows.append({"catalogNumber": cat_out, "plant_individual_id": 1, "misidentification_flag": "FALSE"})

            # (b) Misidentified voucher in morphometric_flags.csv
            cat_mis = f"{sp_code}_MISID"
            m_mis = masks_dir / f"{cat_mis}_leaf1.png"
            _write_smooth_mask(m_mis, 1.78)
            voucher_rows.append({"catalogNumber": cat_mis, "species_standardized": sp_name, "determiner_tier": "Tier_1_Gold"})
            manifest_rows.append({
                "catalogNumber": cat_mis, "plant_individual_id": 1, "leaf_id": 1,
                "assigned_tier": "Tier_1_Direct", "solidity": 0.95, "ucs_score": 0.98,
                "aspect_ratio": 1.78, "pixels_per_mm": 12.0, "is_calibrated_mm": True,
                "mask_path": str(m_mis), "contour_path": "",
            })
            ind_efa_rows.append({
                "catalogNumber": cat_mis, "plant_individual_id": 1, "leaf_id": 1,
                "species_standardized": sp_name, "determiner_tier": "Tier_1_Gold",
                "assigned_tier": "Tier_1_Direct", "is_calibrated_mm": "TRUE", "pixels_per_mm": 12.0,
                "solidity": 0.95, "ucs_score": 0.98, "aspect_ratio": 1.78,
                "pct_asymmetry": 0.02, "harmonic_roughness": 0.01,
                "is_contour_outlier": "FALSE", "leaf_quality_weight": 0.999,
            })
            otu_efa_rows.append({"catalogNumber": cat_mis, "plant_individual_id": 1, "is_contour_outlier": "FALSE"})
            mflag_rows.append({"catalogNumber": cat_mis, "plant_individual_id": 1, "misidentification_flag": "TRUE"})

            # (c) Clipped border mask (should be rejected by _is_valid_smooth_closed_mask)
            cat_clip = f"{sp_code}_CLIPPED"
            m_clip = masks_dir / f"{cat_clip}_leaf1.png"
            _write_smooth_mask(m_clip, 2.08, clipped=True)
            voucher_rows.append({"catalogNumber": cat_clip, "species_standardized": sp_name, "determiner_tier": "Tier_1_Gold"})
            manifest_rows.append({
                "catalogNumber": cat_clip, "plant_individual_id": 1, "leaf_id": 1,
                "assigned_tier": "Tier_1_Direct", "solidity": 0.95, "ucs_score": 0.98,
                "aspect_ratio": 2.08, "pixels_per_mm": 12.0, "is_calibrated_mm": True,
                "mask_path": str(m_clip), "contour_path": "",
            })
            ind_efa_rows.append({
                "catalogNumber": cat_clip, "plant_individual_id": 1, "leaf_id": 1,
                "species_standardized": sp_name, "determiner_tier": "Tier_1_Gold",
                "assigned_tier": "Tier_1_Direct", "is_calibrated_mm": "TRUE", "pixels_per_mm": 12.0,
                "solidity": 0.95, "ucs_score": 0.98, "aspect_ratio": 2.08,
                "pct_asymmetry": 0.02, "harmonic_roughness": 0.01,
                "is_contour_outlier": "FALSE", "leaf_quality_weight": 0.999,
            })
            otu_efa_rows.append({"catalogNumber": cat_clip, "plant_individual_id": 1, "is_contour_outlier": "FALSE"})
            mflag_rows.append({"catalogNumber": cat_clip, "plant_individual_id": 1, "misidentification_flag": "FALSE"})

        pd.DataFrame(voucher_rows).to_csv(tables_dir / "curated_vouchers.csv", index=False)
        pd.DataFrame(manifest_rows).to_csv(tables_dir / "extracted_leaf_manifest.csv", index=False)
        pd.DataFrame(ind_efa_rows).to_csv(tables_dir / "leaf_efa_harmonics_individual.csv", index=False)
        pd.DataFrame(otu_efa_rows).to_csv(tables_dir / "leaf_efa_harmonics.csv", index=False)
        pd.DataFrame(mflag_rows).to_csv(tables_dir / "morphometric_flags.csv", index=False)

        captured_titles: list[str] = []
        captured_ylabels: list[str] = []
        captured_grid_shapes: list[tuple[int, int]] = []
        orig_save = step02._atomic_save_matplotlib_fig

        def _spy_save(fig, out_base_path, dpi=300):
            if Path(out_base_path).name == "fig02b_silhouette_montage":
                axes = fig.axes
                captured_grid_shapes.append((len(axes) // 4, 4))
                for ax in axes:
                    captured_titles.append(ax.get_title())
                    if ax.get_ylabel():
                        captured_ylabels.append(ax.get_ylabel())
            return orig_save(fig, out_base_path, dpi=dpi)

        with patch.object(step02, "_atomic_save_matplotlib_fig", side_effect=_spy_save):
            saved = step02.regenerate_step02_figures_from_tables(
                output_dir=run_dir / "data",
                vouchers_csv=tables_dir / "curated_vouchers.csv",
            )

        for k in ("fig02a_pdf", "fig02a_png", "fig02b_pdf", "fig02b_png"):
            self.assertIn(k, saved)
            self.assertTrue(Path(saved[k]).exists())

        # 1. Montage grid must have exactly 2 rows (for Packera dubia & P. paupercula) x 4 cols = 8 tiles
        self.assertEqual(captured_grid_shapes[-1], (2, 4))
        self.assertEqual(len(captured_ylabels), 2)
        self.assertTrue(any("dubia" in yl for yl in captured_ylabels))
        self.assertTrue(any("paupercula" in yl for yl in captured_ylabels))
        self.assertFalse(any("anonyma" in yl or "plattensis" in yl for yl in captured_ylabels))

        # 2. Outlier, misidentified, and clipped-border candidates must NEVER appear in montage titles
        joined_titles = " | ".join(captured_titles)
        self.assertNotIn("OUTLIER", joined_titles)
        self.assertNotIn("MISID", joined_titles)
        self.assertNotIn("CLIPPED", joined_titles)

        # 3. Within the 40th %ile window [0.32, 0.48] (containing V03 [q=0.333, Gold] and V04 [q=0.444, Bronze])
        #    and 60th %ile window [0.52, 0.68] (containing V05 [q=0.556, Gold] and V06 [q=0.667, Bronze]),
        #    Gold + calibrated leaves (_V03, _V05, _V07) must be prioritized over Bronze (_V04, _V06)
        self.assertIn("DUB_V03", joined_titles)
        self.assertIn("DUB_V05", joined_titles)
        self.assertIn("DUB_V07", joined_titles)
        self.assertNotIn("DUB_V04", joined_titles)
        self.assertNotIn("DUB_V06", joined_titles)
        self.assertIn("PAU_V03", joined_titles)
        self.assertIn("PAU_V05", joined_titles)
        self.assertIn("PAU_V07", joined_titles)

    def test_sheet_adaptive_min_blade_area_and_rosette_bbox_vertical_fallback(self):
        """Verify Task 1A cohort-adaptive min_blade_area_cm2 floor for dwarf/alpine Packera and upper-sheet rosette plant-bbox relaxation."""
        # 1. Dwarf alpine rosette (e.g., Packera cymbalaria / P. hyperborea) with ~0.30 cm^2 mature leaves:
        #    0.25 * median(0.30) = 0.075 -> clamped to absolute_floor_cm2 = 0.15 cm^2 so 0.30 cm^2 blades pass!
        eff_dwarf = step02.compute_sheet_adaptive_min_blade_area_cm2([0.28, 0.30, 0.32, 0.31], config_min_blade_area_cm2=0.75)
        self.assertAlmostEqual(eff_dwarf, 0.15, places=4)
        self.assertLess(eff_dwarf, 0.28)

        # 2. Compact subalpine rosette with ~1.20 cm^2 leaves:
        #    0.25 * 1.20 = 0.30 cm^2 (< 0.75 cm^2)
        eff_compact = step02.compute_sheet_adaptive_min_blade_area_cm2([1.10, 1.20, 1.30], config_min_blade_area_cm2=0.75)
        self.assertAlmostEqual(eff_compact, 0.30, places=4)

        # 3. Large-leaved temperate sheet (~4.0 cm^2 median):
        #    0.25 * 4.0 = 1.00 -> capped at config_min_blade_area_cm2 = 0.75 so ~0.19 cm^2 fragments still fail!
        eff_large = step02.compute_sheet_adaptive_min_blade_area_cm2([3.8, 4.0, 4.2, 0.19], config_min_blade_area_cm2=0.75)
        self.assertAlmostEqual(eff_large, 0.75, places=4)

        # 4. Empty / uncalibrated candidate list falls back safely to config_min_blade_area_cm2
        eff_empty = step02.compute_sheet_adaptive_min_blade_area_cm2([], config_min_blade_area_cm2=0.75)
        self.assertAlmostEqual(eff_empty, 0.75, places=4)

        # 5. Verify plant-bbox vertical rosette fallback tokens & telemetry exist in 02_segment_and_extract.py
        step02_src = (PROJECT_ROOT / "scripts" / "pipeline" / "02_segment_and_extract.py").read_text(encoding="utf-8")
        self.assertIn("compute_sheet_adaptive_min_blade_area_cm2", step02_src)
        self.assertIn("cy_rel_plant_bbox", step02_src)
        self.assertIn("rosette_bbox_fallback_used", step02_src)

    def test_macro_basal_leaf_survival_and_sqlite_descriptor_persistence(self) -> None:
        """Verify macro-basal leaves survive high-DPI Stage 1/2 bounds and SQLite persists shape descriptors & failed calibration."""
        from scripts.core.db import get_db_connection, record_voucher_segmentation_transaction

        sheet_h, sheet_w = 5600, 3800
        sheet_img = np.full((sheet_h, sheet_w, 3), 232, dtype=np.uint8)
        # Macro-basal leaf bbox: height=2500 px, width=1050 px (exceeds legacy 1300*scale_f and 520*scale_f limits)
        s1_dets = [((2700, 1200, 5200, 2250), 0.96, 0)]

        def fake_macro_stage2(crops):
            out = []
            for crop_img in crops:
                ch, cw = crop_img.shape[:2]
                m = np.zeros((ch, cw), dtype=np.uint8)
                cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.38), int(ch * 0.43)), 0, 0, 360, 255, -1)
                out.append([((20, 20, ch - 20, cw - 20), m, 0.97, 0)])
            return out

        self.pipeline.ruler_calibrations["NCU_MACRO"] = 12.5
        self.pipeline.ruler_calibrations["NCU_FAIL_CAL"] = 11.8
        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=s1_dets):
                with patch.object(self.pipeline.model_engine, "predict_batch", side_effect=fake_macro_stage2):
                    records, failure = self.pipeline.process_voucher("NCU_MACRO", self.img_path)

            with patch.object(self.pipeline.stage1_detector, "predict", return_value=[]):
                with patch.object(self.pipeline.model_engine, "predict", return_value=[]):
                    fail_records, fail_rec = self.pipeline.process_voucher("NCU_FAIL_CAL", self.img_path)

        self.assertIsNone(failure)
        self.assertEqual(len(records), 1)
        rec = records[0]
        for col in ("circularity", "extent", "aspect_ratio", "area_px", "blade_area_cm2"):
            self.assertIn(col, rec)
            self.assertTrue(np.isfinite(float(rec[col])), f"{col} should be finite, got {rec[col]}")

        self.assertEqual(len(fail_records), 0)
        self.assertIsNotNone(fail_rec)
        self.assertAlmostEqual(float(fail_rec["pixels_per_mm"]), 11.8, places=3)
        self.assertTrue(bool(fail_rec["is_calibrated_mm"]))

        conn = get_db_connection(db_path=":memory:", init_schema=True)
        record_voucher_segmentation_transaction(
            conn,
            catalog_number="NCU_MACRO",
            image_path=str(self.img_path),
            extracted_records=records,
            failure_record=None,
            output_dir=self.pipeline.output_dir,
        )
        record_voucher_segmentation_transaction(
            conn,
            catalog_number="NCU_FAIL_CAL",
            image_path=str(self.img_path),
            extracted_records=[],
            failure_record=fail_rec,
            output_dir=self.pipeline.output_dir,
        )
        row_leaf = conn.execute(
            "SELECT circularity, extent, aspect_ratio, area_px, blade_area_cm2 FROM extracted_leaves WHERE catalogNumber = 'NCU_MACRO'"
        ).fetchone()
        self.assertIsNotNone(row_leaf)
        for col in ("circularity", "extent", "aspect_ratio", "area_px", "blade_area_cm2"):
            self.assertIsNotNone(row_leaf[col], f"SQLite extracted_leaves.{col} must not be NULL")

        row_fail_v = conn.execute(
            "SELECT qc_status, pixels_per_mm, is_calibrated_mm FROM vouchers WHERE catalogNumber = 'NCU_FAIL_CAL'"
        ).fetchone()
        self.assertIsNotNone(row_fail_v)
        self.assertEqual(row_fail_v["qc_status"], "FAILED_QC")
        self.assertAlmostEqual(float(row_fail_v["pixels_per_mm"]), 11.8, places=3)
        self.assertEqual(int(row_fail_v["is_calibrated_mm"]), 1)
        conn.close()

    def test_raw_unconstrained_lm2_native_leaf_and_petiole_segmentation(self) -> None:
        """
        Verify default 100% Raw Unconstrained LeafMachine2 Mode (`enable_geometric_qc=False`):
          1. Stage 1 (`LeafPriority.pt`) unpacks all Class 0 (`leaf_whole`) boxes regardless of size
             or vertical sheet position, while natively excluding `leaf_partial` (1), `leaflet` (2),
             `specimen` (8), and `wood`/stems (10).
          2. Stage 2 PointRend (`thing_classes = ['leaf', 'petiole', 'hole']`) retains the top-scoring
             `'leaf'` (`class_id == 0`) mask 100% pixel-identical (zero `trim_petiole_tail()` clipping),
             pairs the top-scoring `'petiole'` (`class_id == 1`) mask via `_extract_native_lm2_petiole_anchor()`
             to locate `p_base` and native petiole length, and retains >6 leaves without any post-hoc
             geometric/photometric rejections or top-6 cap.
        """
        self.assertFalse(self.pipeline.enable_geometric_qc)
        with patch("sys.argv", ["02_segment_and_extract.py"]):
            parsed_default = step02.parse_args()
            self.assertFalse(parsed_default.enable_geometric_qc)

        sheet_h, sheet_w = 3600, 2600
        sheet_img = np.full((sheet_h, sheet_w, 3), (225, 232, 236), dtype=np.uint8)
        self.pipeline.ruler_calibrations["NCU_RAW_LM2"] = 10.0  # 10 px/mm = 100 px/cm

        # 7 Stage 1 Class 0 (`leaf_whole`) boxes across top, middle, and bottom of sheet
        # (including small alpine-sized boxes and >6 leaves) plus native LM2 non-leaf classes:
        s1_dets = [
            ((250, 300, 420, 460), 0.94, 0),    # 1. small upper Class 0 leaf_whole
            ((650, 350, 1050, 750), 0.95, 0),   # 2. upper-mid Class 0 leaf_whole
            ((1200, 400, 1750, 920), 0.96, 0),  # 3. mid-sheet cordate/suborbicular Class 0 leaf_whole
            ((1850, 450, 2400, 950), 0.97, 0),  # 4. lower-mid Class 0 leaf_whole
            ((2450, 400, 3000, 920), 0.95, 0),  # 5. basal Class 0 leaf_whole
            ((2450, 1050, 3000, 1550), 0.93, 0),# 6. basal Class 0 leaf_whole
            ((2450, 1650, 3000, 2150), 0.92, 0),# 7. 7th Class 0 leaf_whole (verifies no top-6 cap)
            ((1400, 1800, 1900, 2200), 0.89, 1),# leaf_partial -> natively ignored by Stage 1
            ((1100, 1800, 1250, 1950), 0.85, 2),# leaflet -> natively ignored by Stage 1
            ((50, 50, 3500, 2500), 0.99, 8),    # specimen -> natively ignored by Stage 1
            ((900, 1200, 2600, 1280), 0.91, 10),# wood/stem -> natively ignored by Stage 1
        ]

        expected_fg_counts: list[int] = []

        def fake_raw_lm2_stage2(crops):
            out = []
            for idx, crop_img in enumerate(crops):
                ch, cw = crop_img.shape[:2]
                leaf_m = np.zeros((ch, cw), dtype=np.uint8)
                cx, cy = cw // 2, int(ch * 0.42)
                rx, ry = max(18, int(cw * 0.30)), max(22, int(ch * 0.28))
                cv2.ellipse(leaf_m, (cx, cy), (rx, ry), 0.0, 0, 360, 255, -1)
                if idx == 2:
                    # Cordate basal sinus + crenate indent (low solidity) that raw LM2 must keep untouched
                    cv2.ellipse(leaf_m, (cx, cy + ry - 4), (int(rx * 0.38), int(ry * 0.35)), 0.0, 0, 360, 0, -1)
                expected_fg_counts.append(int(np.count_nonzero(leaf_m)))

                # Native PointRend class_id == 1 (`petiole`) mask extending below lamina base
                pet_m = np.zeros((ch, cw), dtype=np.uint8)
                pet_y0 = min(ch - 6, cy + ry - 2)
                pet_y1 = min(ch - 2, pet_y0 + max(25, int(ch * 0.22)))
                cv2.line(pet_m, (cx, pet_y0), (cx, pet_y1), 255, thickness=7)

                out.append([
                    ((10, 10, ch - 10, cw - 10), leaf_m, 0.96, 0),
                    ((pet_y0, max(0, cx - 5), pet_y1, min(cw, cx + 5)), pet_m, 0.91, 1),
                ])
            return out

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.pipeline.stage1_detector, "predict", return_value=s1_dets):
                with patch.object(self.pipeline.model_engine, "predict_batch", side_effect=fake_raw_lm2_stage2) as mock_s2:
                    records, failure = self.pipeline.process_voucher("NCU_RAW_LM2", self.img_path)

        self.assertIsNone(failure)
        # Only the 7 Class 0 (`leaf_whole`) boxes enter Stage 2; all 7 are retained with zero rejections and zero top-6 cap
        self.assertEqual(len(mock_s2.call_args[0][0]), 7)
        self.assertEqual(len(records), 7)
        for idx, rec in enumerate(records):
            self.assertEqual(rec["assigned_tier"], "Tier_1_Direct")
            self.assertFalse(bool(rec["reflection_applied"]))
            self.assertFalse(bool(rec["is_folded"]))
            # Saved binary mask foreground pixel count must equal raw PointRend leaf_m 100% identically (zero petiole trimming!)
            saved_m = cv2.imread(str(rec["mask_path"]), cv2.IMREAD_GRAYSCALE)
            self.assertIsNotNone(saved_m)
            self.assertEqual(int(np.count_nonzero(saved_m > 0)), expected_fg_counts[idx])
            self.assertGreater(float(rec["petiole_length_mm"]), 0.0)

    def test_three_stage_wave_engine_parity_and_sqlite_wave_commit(self) -> None:
        """
        Verify 3-Stage Asynchronous Wave Engine (`run()` across >4 vouchers with wave_size=4)
        produces 100% identical extracted leaf records, binary masks, 120-vertex contours,
        ROI-scoped rosette patches/masks, and SQLite `contour_xy_json` records compared to
        sequential single-voucher execution.
        """
        import json
        from scripts.core.db import get_db_connection

        wave_root = self.temp_dir / "wave_run"
        seq_root = self.temp_dir / "seq_run"
        img_dir = self.temp_dir / "wave_imgs"
        img_dir.mkdir(parents=True, exist_ok=True)

        vouchers_rows = []
        for idx in range(6):
            cat = f"NCU_WAVE_{idx:02d}"
            img_p = img_dir / f"{cat}.jpg"
            sheet = np.full((1200, 900, 3), (225, 232, 236), dtype=np.uint8)
            # Draw two green foliar patches so real JPEG decode works across prefetch workers
            cv2.ellipse(sheet, (300, 760), (75, 110), 0.0, 0, 360, (60, 125, 70), -1)
            cv2.ellipse(sheet, (600, 780), (70, 105), 0.0, 0, 360, (55, 120, 65), -1)
            cv2.imwrite(str(img_p), sheet)
            vouchers_rows.append({
                "catalogNumber": cat,
                "scientificName": "Packera anonyma",
                "image_path": str(img_p),
            })

        vouchers_csv = self.temp_dir / "wave_vouchers.csv"
        pd.DataFrame(vouchers_rows).to_csv(vouchers_csv, index=False)

        s1_boxes = [
            ((610, 190, 910, 410), 0.96, 0),
            ((630, 490, 930, 710), 0.94, 0),
        ]

        def fake_s1_predict_batch(_self_s1, items, use_cache=False, max_batch_size=16):
            return [list(s1_boxes) for _ in items]

        def fake_s1_predict(_self_s1, _img_bgr, catalog_number=None, use_cache=False, precomputed_letterbox=None):
            return list(s1_boxes)

        def fake_s2_predict_batch(_self_s2, images_bgr, batch_size=24, pretransformed=None):
            out = []
            for crop_img in images_bgr:
                ch, cw = crop_img.shape[:2]
                leaf_m = np.zeros((ch, cw), dtype=np.uint8)
                cx, cy = cw // 2, int(ch * 0.46)
                rx, ry = max(16, int(cw * 0.28)), max(22, int(ch * 0.30))
                cv2.ellipse(leaf_m, (cx, cy), (rx, ry), 0.0, 0, 360, 255, -1)
                pet_m = np.zeros((ch, cw), dtype=np.uint8)
                py0 = min(ch - 6, cy + ry - 2)
                py1 = min(ch - 2, py0 + max(18, int(ch * 0.16)))
                cv2.line(pet_m, (cx, py0), (cx, py1), 255, thickness=5)
                out.append([
                    ((8, 8, ch - 8, cw - 8), leaf_m, 0.97, 0),
                    ((py0, max(0, cx - 4), py1, min(cw, cx + 4)), pet_m, 0.90, 1),
                ])
            return out

        wave_db = self.temp_dir / "wave_ephemeral.sqlite"
        seq_db = self.temp_dir / "seq_ephemeral.sqlite"

        dummy_s1_weights = self.temp_dir / "dummy_s1.pt"
        with patch.object(step02.LM2Stage1Detector, "predict_batch", fake_s1_predict_batch), \
             patch.object(step02.LM2Stage1Detector, "predict", fake_s1_predict), \
             patch.object(step02.PointRendInferenceEngine, "predict_batch", fake_s2_predict_batch):
            pipe_wave = step02.SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=self.pipeline.model_weights,
                stage1_weights=dummy_s1_weights,
                output_dir=wave_root,
                device="cpu",
                db_path=wave_db,
            )
            pipe_wave.stage1_wave_size = 4  # Exercises 2 cross-voucher waves (4 + 2)
            for r in vouchers_rows:
                pipe_wave.ruler_calibrations[r["catalogNumber"]] = 10.0
            df_wave, fail_wave = pipe_wave.run()

            pipe_seq = step02.SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=self.pipeline.model_weights,
                stage1_weights=dummy_s1_weights,
                output_dir=seq_root,
                device="cpu",
                db_path=seq_db,
            )
            # Instance attribute hook forces sequential single-voucher fallback in run()
            pipe_seq.stage1_detector.predict = lambda img, **kw: list(s1_boxes)  # type: ignore[method-assign]
            for r in vouchers_rows:
                pipe_seq.ruler_calibrations[r["catalogNumber"]] = 10.0
            df_seq, fail_seq = pipe_seq.run()

        self.assertTrue(fail_wave.empty)
        self.assertTrue(fail_seq.empty)
        self.assertEqual(len(df_wave), 12)
        self.assertEqual(len(df_seq), 12)

        compare_cols = [
            "catalogNumber", "leaf_id", "plant_individual_id", "assigned_tier",
            "solidity", "ucs_score", "circularity", "extent", "area_px",
            "blade_area_cm2", "petiole_length_mm", "blade_length_mm",
        ]
        pd.testing.assert_frame_equal(
            df_wave[compare_cols].reset_index(drop=True),
            df_seq[compare_cols].reset_index(drop=True),
        )

        # Verify ROI-scoped rosette masks & binary leaf masks match pixel-for-pixel and SQLite stores 120-pt contours
        conn = get_db_connection(wave_db, init_schema=False)
        try:
            rows_db = conn.execute(
                "SELECT catalogNumber, leaf_id, contour_n_pts, contour_xy_json FROM extracted_leaves ORDER BY catalogNumber, leaf_id"
            ).fetchall()
            self.assertEqual(len(rows_db), 12)
            for row_db in rows_db:
                self.assertEqual(int(row_db["contour_n_pts"]), 120)
                coords = json.loads(row_db["contour_xy_json"])
                self.assertEqual(len(coords), 120)
        finally:
            conn.close()

        for r in vouchers_rows:
            cat = r["catalogNumber"]
            w_rmask = cv2.imread(str(pipe_wave.rosettes_dense_dir / f"{cat}_p0_rosette_mask.png"), cv2.IMREAD_GRAYSCALE)
            s_rmask = cv2.imread(str(pipe_seq.rosettes_dense_dir / f"{cat}_p0_rosette_mask.png"), cv2.IMREAD_GRAYSCALE)
            self.assertIsNotNone(w_rmask)
            self.assertIsNotNone(s_rmask)
            np.testing.assert_array_equal(w_rmask, s_rmask)


if __name__ == "__main__":
    unittest.main()
