#!/usr/bin/env python3
"""
===============================================================================
Script: test_segment_and_extract.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 02:
    LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime) segmentation, 2-path botanical
    leaf extraction (Tier 1 Pristine and Tier 2 Bilateral Reflection), LM3 ruler
    scale calibration, centroid spatial clustering, and standardized contour CSV generation.
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
export_standardized_contour = step02.export_standardized_contour
extract_tier1_pristine = step02.extract_tier1_pristine
extract_tier2_reflected = step02.extract_tier2_reflected
compute_capitulum_metrics = step02.compute_capitulum_metrics


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

    def test_lm3_ruler_lattice_and_mp_prior_scale_resolution(self) -> None:
        """Verify measured_from_ruler, measured_from_fieldprism, and predicted_from_megapixels scale resolution in lm3_data_loader."""
        from scripts.vision import lm3_data_loader

        ruler_dir = self.temp_dir / "lm3_scale_res" / "Data" / "Ruler"
        ruler_dir.mkdir(parents=True, exist_ok=True)
        ruler_csv = ruler_dir / "ruler_summary.csv"
        pd.DataFrame([
            {
                "filename": "NCU_RULER__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 118.11,
                "ruler_class_confidence": 96.5,
                "conversion_determination": "measured_from_ruler",
            },
            {
                "filename": "NCU_FIELDPRISM__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 125.0,
                "ruler_class_confidence": 92.0,
                "conversion_determination": "measured_from_fieldprism",
            },
            {
                "filename": "NCU_MP_PRIOR__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 105.0,
                "ruler_class_confidence": 85.0,
                "conversion_determination": "predicted_from_megapixels",
            },
            {
                "filename": "NCU_UNRELIABLE__R__10-10-100-500",
                "ruler_success": False,
                "conversion_mean": 0.0,
                "ruler_class_confidence": 12.0,
            },
        ]).to_csv(ruler_csv, index=False)

        calibs = lm3_data_loader.load_ruler_calibrations(self.temp_dir / "lm3_scale_res", min_confidence=0.50)
        self.assertAlmostEqual(calibs["NCU_RULER"], 11.811, places=3)
        self.assertAlmostEqual(calibs["NCU_FIELDPRISM"], 12.500, places=3)
        self.assertAlmostEqual(calibs["NCU_MP_PRIOR"], 10.500, places=3)
        self.assertNotIn("NCU_UNRELIABLE", calibs)

        empty_calibs = lm3_data_loader.load_ruler_calibrations(self.temp_dir / "nonexistent_lm3_dir")
        self.assertEqual(empty_calibs, {})

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

    def test_lm3_tip_up_contour_export_preserves_lamina_base_anchor(self) -> None:
        """Verify 120-vertex clockwise contour export anchored at lamina_base (petiole_attachment_pt) via resample_normalized_contour and export_standardized_contour."""
        from scripts.vision.leaf_geometry_utils import resample_normalized_contour

        mask = np.zeros((300, 300), dtype=np.uint8)
        # Vertical ellipse: center (150, 150), axes (40, 100) -> apex ~(150, 50), lamina_base ~(150, 250)
        cv2.ellipse(mask, (150, 150), (40, 100), 0, 0, 360, 255, -1)

        p_base = (150, 250)
        resampled = resample_normalized_contour(mask, num_points=120, petiole_attachment_pt=p_base)
        self.assertEqual(resampled.shape, (120, 2))
        self.assertLess(cv2.contourArea(resampled.astype(np.float32), oriented=True), 0)

        csv_path = export_standardized_contour(
            mask, "TEST_HOMO", 1, self.temp_dir, num_points=120, petiole_attachment_pt=p_base
        )
        self.assertIsNotNone(csv_path)
        df = pd.read_csv(csv_path)
        self.assertEqual(len(df), 120)

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
        """Verify folded/partial leaf is excluded from closed-outline EFA under opt-in geometric QC."""
        # Half-ellipse with straight edge along midrib
        mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(mask, (200, 200), (35, 130), 0, -90, 90, 255, -1)

        from unittest.mock import patch
        with patch.object(self.qc_pipeline.model_engine, "predict", return_value=[((70, 160, 330, 240), mask, 0.92, 0)]):
            with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
                records, failure = self.qc_pipeline.process_voucher("NCU123", self.img_path)

        self.assertEqual(len(records), 0, "Mirrored Tier 2 hemi-blades must not enter closed-outline EFA records.")
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_reason"], "PARTIAL_OR_FOLDED_EXCLUDED_FROM_EFA")

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
        """Verify Stage 1 (LM3PlantDetector 11-class) -> 12%-padded crop -> Stage 2 LM3LeafSegmenterEngine reprojection."""
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
        # Stage 2 LM3LeafSegmenterEngine must be called ONLY for class 0 (leaf_whole) and class 1 (leaf_partial)
        self.assertEqual(mock_s2.call_count, 2)
        self.assertEqual(len(records), 2)
        for rec in records:
            self.assertEqual(rec["assigned_tier"], "Tier_1_Direct")
            self.assertTrue(Path(rec["contour_path"]).exists())
            self.assertTrue(Path(rec["mask_path"]).exists())

    def test_reproductive_routing_and_capitula_patch_saving(self) -> None:
        """Verify Stage 1 reproductive detections (classes 2..6) route to phenology classifier & save capitula crops."""
        from unittest.mock import patch

        sheet_img = np.full((1000, 800, 3), 230, dtype=np.uint8)
        # Paint golden-yellow corollas inside the Flower_ONE (class 4) region (50..170, 300..400)
        # BGR for golden yellow (RGB ~ 240, 200, 20 -> BGR = 20, 200, 240)
        sheet_img[60:120, 310:390] = (20, 200, 240)
        # Paint dark green erect involucre at lower half (120..165, 325..375)
        sheet_img[120:165, 325:375] = (30, 90, 40)

        stage1_dets = [
            ((300, 200, 550, 420), 0.95, 0),  # 0: Leaf_WHOLE
            ((50, 300, 170, 400), 0.92, 4),   # 4: Flower_ONE (H=120, W=100, AR=1.2)
            ((60, 450, 170, 550), 0.88, 2),   # 2: Seed_Fruit_ONE (H=110, W=100, AR=1.1)
            ((70, 600, 160, 680), 0.84, 6),   # 6: Bud (H=90, W=80, AR=1.125)
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

    def test_multi_plant_clustering_and_tier1_mask_generation(self) -> None:
        """Verify multi-plant DBSCAN clustering and Tier 1 binary mask saving per plant_individual_id."""
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

        # Verify Tier 1 binary leaf masks exist in data/masks/tier1_intact/
        for leaf_id in [1, 2, 3]:
            leaf_mask = self.pipeline.masks_dir / "tier1_intact" / f"NCU123_leaf{leaf_id}.png"
            self.assertTrue(leaf_mask.exists(), f"Tier 1 binary leaf mask {leaf_mask} must exist.")

    def test_lm3_ruler_calibration_resolution(self) -> None:
        """Verify load_ruler_calibrations() resolves measured_from_ruler and predicted_from_megapixels scales."""
        from unittest.mock import patch

        ruler_dir = self.temp_dir / "lm3_ruler_test" / "Data" / "Ruler"
        ruler_dir.mkdir(parents=True, exist_ok=True)
        ruler_csv = ruler_dir / "ruler_summary.csv"
        pd.DataFrame([
            {
                "filename": "NCU_HIGH_CONF__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 118.11,  # 118.11 px/cm = 11.811 px/mm
                "ruler_class_confidence": 96.5,
                "conversion_determination": "measured_from_ruler",
            },
            {
                "filename": "NCU_MP_PRIOR__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 105.0,  # 105.0 px/cm = 10.5 px/mm
                "ruler_class_confidence": 85.0,
                "conversion_determination": "predicted_from_megapixels",
            },
            {
                "filename": "NCU_LOW_CONF__R__10-10-100-500",
                "ruler_success": True,
                "conversion_mean": 150.0,
                "ruler_class_confidence": 20.0,  # Below 0.50 confidence threshold -> should be filtered out
            },
        ]).to_csv(ruler_csv, index=False)

        calibs = step02.load_ruler_calibrations(self.temp_dir / "lm3_ruler_test", min_confidence=0.50)
        self.assertIn("NCU_HIGH_CONF", calibs)
        self.assertIn("NCU_MP_PRIOR", calibs)
        self.assertNotIn("NCU_LOW_CONF", calibs)
        self.assertAlmostEqual(calibs["NCU_HIGH_CONF"], 11.811, places=3)
        self.assertAlmostEqual(calibs["NCU_MP_PRIOR"], 10.5, places=3)

        self.pipeline.ruler_calibrations = calibs
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 150), (65, 120), 0, 0, 360, 255, -1)

        with patch("cv2.imread", return_value=np.full((400, 400, 3), 200, dtype=np.uint8)):
            with patch.object(self.pipeline.model_engine, "predict", return_value=[((30, 85, 270, 215), mask, 0.95, 0)]):
                recs_high, _ = self.pipeline.process_voucher("NCU_HIGH_CONF", self.img_path)
                self.assertAlmostEqual(recs_high[0]["pixels_per_mm"], 11.811, places=3)
                self.assertTrue(recs_high[0]["is_calibrated_mm"])

                recs_mp, _ = self.pipeline.process_voucher("NCU_MP_PRIOR", self.img_path)
                self.assertAlmostEqual(recs_mp[0]["pixels_per_mm"], 10.5, places=3)
                self.assertTrue(recs_mp[0]["is_calibrated_mm"])

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

    def test_extract_tier2_reflected_45deg_asymmetric_leaf_iou(self) -> None:
        """
        Verify a synthetic asymmetric leaf mask oriented at 45 degrees passed into
        extract_tier2_reflected() produces a bilaterally symmetric mask with >0.95 IoU
        against the ground-truth unoccluded blade.
        """
        from scripts.vision.leaf_geometry_utils import extract_tier2_reflected

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
        Verify an oblong tip-up leaf with L/W = 3.0 exported via export_standardized_contour()
        retains (max(y_norm) - min(y_norm)) / (max(x_norm) - min(x_norm)) == 3.0 within 1e-3.
        """
        mask = np.zeros((500, 500), dtype=np.uint8)
        # Oblong tip-up leaf with semi-axes a=150 (L=300 along y), b=50 (W=100 along x) -> L/W = 3.0
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
        lw_ratio = max(span_x_norm, span_y_norm) / max(min(span_x_norm, span_y_norm), 1e-9)
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
        from scripts.vision.leaf_geometry_utils import classify_and_route_leaf

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

    def test_classify_and_route_stem_clump_rejection_before_reflection(self) -> None:
        """Verify classify_and_route_leaf rejects stems and clumps (solidity < 0.58 or ucs < 0.62) without Tier 2 reflection."""
        from scripts.vision.leaf_geometry_utils import classify_and_route_leaf

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
        """Verify on high-DPI herbarium sheets (>2000px), Stage 1 forwards all Class 0 (`Leaf_WHOLE`) boxes to Stage 2 while skipping Class 1 (`Leaf_PARTIAL`), and post-segmentation filters tiny noise."""
        sheet_img = np.full((3000, 2400, 3), 230, dtype=np.uint8)
        self.assertTrue(max(sheet_img.shape[:2]) > 2000)

        # Stage 1 detections:
        # 1. Class 1 (leaf_partial) -> skipped when Class 0 is present and segment_partial_leaves=False
        # 2. Class 0, upper-sheet plant leaf: cy = 600 -> forwarded to Stage 2 and retained!
        # 3. Class 0, tiny 80x80 box near lower rosette -> forwarded to Stage 2 and filtered post-segmentation
        # 4. Class 0, lower-sheet basal leaf box: ymin=1800, xmin=800, ymax=2400, xmax=1400 -> forwarded and retained!
        stage1_dets = [
            ((1500, 500, 1900, 900), 0.90, 1),   # leaf_partial -> skipped
            ((400, 800, 800, 1200), 0.92, 0),    # upper-sheet Class 0 leaf_whole -> forwarded
            ((1600, 600, 1680, 680), 0.88, 0),   # 80x80 Class 0 box -> forwarded to Stage 2, rejected post-seg
            ((1800, 800, 2400, 1400), 0.95, 0),  # lower-sheet Class 0 leaf_whole -> forwarded
        ]

        def fake_high_dpi_stage2(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            # Prediction 1: Class 1 (petiole/stem fragment)
            m1 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.rectangle(m1, (10, 10), (min(cw - 1, 30), min(ch - 1, 200)), 255, -1)

            # Prediction 2: Class 0, tiny noise fragment
            m2 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.circle(m2, ( min(cw - 1, 50), min(ch - 1, 50)), 5, 255, -1)

            # Prediction 3: Class 0, primary blade in crop
            m3 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m3, (cw // 2, ch // 2), (max(8, cw // 3), max(8, ch // 3)), 0, 0, 360, 255, -1)

            # Prediction 4: Class 0, smaller blade fragment (lower score)
            m4 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m4, (cw // 2, ch // 2), (max(5, cw // 5), max(5, ch // 5)), 0, 0, 360, 255, -1)

            return [
                ((10, 10, min(ch, 200), min(cw, 30)), m1, 0.90, 1),
                ((10, 10, min(ch, 30), min(cw, 30)), m2, 0.85, 0),
                ((10, 10, max(20, ch - 10), max(20, cw - 10)), m3, 0.95, 0),
                ((20, 20, max(30, ch - 20), max(30, cw - 20)), m4, 0.70, 0),
            ]

        with patch("cv2.imread", return_value=sheet_img):
            with patch.object(self.qc_pipeline.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(self.qc_pipeline.model_engine, "predict", side_effect=fake_high_dpi_stage2) as mock_s2:
                    records, failure = self.qc_pipeline.process_voucher("NCU_HIGHDPI", self.img_path)

        self.assertIsNone(failure)
        # Stage 2 LM3LeafSegmenterEngine is called for all 3 Class 0 (`Leaf_WHOLE`) boxes (Class 1 skipped)
        self.assertEqual(mock_s2.call_count, 3)
        # Upper-sheet and lower-sheet whole leaves are retained; tiny 80x80 candidate is rejected post-segmentation
        self.assertEqual(len(records), 2)
        self.assertTrue(all(r["assigned_tier"] == "Tier_1_Direct" for r in records))

    def test_stage1_use_cache_flag(self) -> None:
        """Verify LM3PlantDetector.predict respects use_cache flag."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            proj_dir = Path(tmp_dir)
            labels_dir = proj_dir / "Plant_Components" / "labels"
            labels_dir.mkdir(parents=True)
            label_file = labels_dir / "VOUCHER1.txt"
            label_file.write_text("0 0.5 0.5 0.2 0.2\n")

            detector = step02.LM3PlantDetector(
                weights_path="dummy.pt",
                lm3_project_dir=proj_dir,
            )
            # When use_cache=True, loads precomputed labels
            dets_cached = detector.predict(np.zeros((100, 100, 3), dtype=np.uint8), catalog_number="VOUCHER1", use_cache=True)
            self.assertEqual(len(dets_cached), 1)

            # When use_cache=False, bypasses precomputed labels (and returns [] since model is None)
            dets_uncached = detector.predict(np.zeros((100, 100, 3), dtype=np.uint8), catalog_number="VOUCHER1", use_cache=False)
            self.assertEqual(len(dets_uncached), 0)

    def test_heterophylly_cauline_and_stem_rejection_conditions(self) -> None:
        """Verify Milestone 1 cauline/bract/stem rejection conditions (a), (b), and (c)."""
        from scripts.vision.leaf_geometry_utils import is_cauline_or_stem_contaminant

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
        from scripts.vision.leaf_geometry_utils import LeafRoutingResult

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

    def test_raw_unconstrained_lm3_native_leaf_and_petiole_segmentation(self) -> None:
        """
        Verify default 100% Raw Unconstrained LeafMachine3 Mode (`enable_geometric_qc=False`):
          1. Stage 1 (`LM3PlantDetector`) unpacks all Class 0 (`leaf_whole`) boxes regardless of size
             or vertical sheet position, while natively excluding `leaf_partial` (1), `leaflet` (2),
             `specimen` (8), and `wood`/stems (10).
          2. Stage 2 `LM3LeafSegmenterEngine` (`thing_classes = ['leaf', 'petiole', 'hole']`) retains the top-scoring
             `'leaf'` (`class_id == 0`) mask 100% pixel-identical (zero heuristic petiole clipping),
             pairs the top-scoring `'petiole'` (`class_id == 1`) mask to locate `p_base` and native petiole length,
             and retains >6 leaves without any post-hoc geometric/photometric rejections or top-6 cap.
        """
        self.assertFalse(self.pipeline.enable_geometric_qc)
        with patch("sys.argv", ["02_segment_and_extract.py"]):
            parsed_default = step02.parse_args()
            self.assertFalse(parsed_default.enable_geometric_qc)

        sheet_h, sheet_w = 3600, 2600
        sheet_img = np.full((sheet_h, sheet_w, 3), (225, 232, 236), dtype=np.uint8)
        self.pipeline.ruler_calibrations["NCU_RAW_LM3"] = 10.0  # 10 px/mm = 100 px/cm

        # 7 Stage 1 Class 0 (`leaf_whole`) boxes across top, middle, and bottom of sheet
        # (including small alpine-sized boxes and >6 leaves) plus native LM3 non-leaf classes:
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

        def fake_raw_lm3_stage2(crops):
            out = []
            for idx, crop_img in enumerate(crops):
                ch, cw = crop_img.shape[:2]
                leaf_m = np.zeros((ch, cw), dtype=np.uint8)
                cx, cy = cw // 2, int(ch * 0.42)
                rx, ry = max(18, int(cw * 0.30)), max(22, int(ch * 0.28))
                cv2.ellipse(leaf_m, (cx, cy), (rx, ry), 0.0, 0, 360, 255, -1)
                if idx == 2:
                    # Cordate basal sinus + crenate indent (low solidity) that raw LM3 must keep untouched
                    cv2.ellipse(leaf_m, (cx, cy + ry - 4), (int(rx * 0.38), int(ry * 0.35)), 0.0, 0, 360, 0, -1)
                expected_fg_counts.append(int(np.count_nonzero(leaf_m)))

                # Native LM3LeafSegmenterEngine class_id == 1 (`petiole`) mask extending below lamina base
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
                with patch.object(self.pipeline.model_engine, "predict_batch", side_effect=fake_raw_lm3_stage2) as mock_s2:
                    records, failure = self.pipeline.process_voucher("NCU_RAW_LM3", self.img_path)

        self.assertIsNone(failure)
        # Only the 7 Class 0 (`leaf_whole`) boxes enter Stage 2; all 7 are retained with zero rejections and zero top-6 cap
        self.assertEqual(len(mock_s2.call_args[0][0]), 7)
        self.assertEqual(len(records), 7)
        for idx, rec in enumerate(records):
            self.assertEqual(rec["assigned_tier"], "Tier_1_Direct")
            self.assertFalse(bool(rec["reflection_applied"]))
            self.assertFalse(bool(rec["is_folded"]))
            # Saved binary mask foreground pixel count must equal raw LM3 leaf_m 100% identically (zero petiole trimming!)
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
        with patch.object(step02.LM3PlantDetector, "predict_batch", fake_s1_predict_batch), \
             patch.object(step02.LM3PlantDetector, "predict", fake_s1_predict), \
             patch.object(step02.LM3LeafSegmenterEngine, "predict_batch", fake_s2_predict_batch):
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

        # Verify SQLite stores 120-pt contours
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


class TestDissectedLeafSegmentationEnhancements(unittest.TestCase):
    """Tests for high-DPI dissected leaf geometry, rachis bridging, and petiole trimming."""

    @staticmethod
    def _build_high_dpi_pinnatisect_mask(sever_rachis: bool = False, include_petiole: bool = False) -> np.ndarray:
        """Synthesize a high-DPI (1200x900, minor_w > 480 px) lyrate-pinnatisect leaf mask."""
        mask = np.zeros((1200, 900), dtype=np.uint8)
        cx = 450
        # Central rachis spine (from y=180 terminal lobe down to y=880 lowest pinna pair)
        if not sever_rachis:
            cv2.rectangle(mask, (cx - 10, 200), (cx + 10, 880), 255, -1)

        # Terminal rounded-lobate segment
        cv2.ellipse(mask, (cx, 210), (135, 95), 0, 0, 360, 255, -1)
        cv2.ellipse(mask, (cx - 95, 220), (55, 40), -15, 0, 360, 255, -1)
        cv2.ellipse(mask, (cx + 95, 220), (55, 40), 15, 0, 360, 255, -1)

        # 3 pairs of opposing lateral pinnae spanning width ~520 px (> 360 px high-DPI threshold)
        pinna_specs = [
            (410, 250, 48, 12),
            (610, 235, 44, 10),
            (800, 195, 38, 8),
        ]
        for y_c, span_x, lobe_h, angle in pinna_specs:
            # Left pinna
            cv2.ellipse(mask, (cx - span_x // 2, y_c), (span_x // 2, lobe_h), -angle, 0, 360, 255, -1)
            # Right pinna
            cv2.ellipse(mask, (cx + span_x // 2, y_c), (span_x // 2, lobe_h), angle, 0, 360, 255, -1)
            if sever_rachis:
                # Carve a narrow gap across the central rachis internode above this pinna pair
                cv2.rectangle(mask, (cx - 35, y_c - lobe_h - 32), (cx + 35, y_c - lobe_h - 4), 0, -1)

        if include_petiole:
            # Narrow basal petiole stalk below the lowest pinna pair (y=840..1090, width=14 px)
            cv2.rectangle(mask, (cx - 7, 835), (cx + 7, 1090), 255, -1)

        return mask

    def test_high_dpi_dissected_leaf_passes_is_botanical_dissection_and_tier1_routing(self) -> None:
        """High-DPI dissected leaves (minor_w > 360 px) must pass is_botanical_dissection and route to Tier 1."""
        from scripts.vision.leaf_geometry_utils import (
            classify_and_route_leaf,
            is_botanical_dissection,
        )

        mask = self._build_high_dpi_pinnatisect_mask(sever_rachis=False, include_petiole=False)
        cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        main_cnt = max(cnts, key=cv2.contourArea)
        (_, (w_rect, h_rect), _) = cv2.minAreaRect(main_cnt)
        minor_w = float(min(w_rect, h_rect))
        self.assertGreater(minor_w, 400.0, "Fixture must exercise native high-DPI minor_w > 360 px regime.")

        self.assertTrue(
            is_botanical_dissection(mask),
            "High-DPI pinnatisect leaf with minor_w > 360 px must be recognized as botanical dissection.",
        )

        res = classify_and_route_leaf(mask)
        self.assertTrue(res.is_dissected)
        self.assertIn("Tier 1", res.tier)
        self.assertGreaterEqual(res.ucs, 0.80)

    def test_bridge_dissected_rachis_components_reconnects_severed_pinnae_and_noop_on_single_component(self) -> None:
        """Severed pinna islands must be reconnected along the central rachis while single-component masks stay untouched."""
        from scripts.vision.leaf_geometry_utils import bridge_dissected_rachis_components

        # 1. Single-component mask must be a strict 0-pixel no-op
        intact = self._build_high_dpi_pinnatisect_mask(sever_rachis=False)
        bridged_intact = bridge_dissected_rachis_components(intact)
        np.testing.assert_array_equal(intact, bridged_intact)

        # 2. Severed rachis mask with disconnected pinna pairs
        severed = self._build_high_dpi_pinnatisect_mask(sever_rachis=True)
        n_before, _ = cv2.connectedComponents((severed > 0).astype(np.uint8))
        self.assertGreater(n_before, 2, "Severed fixture must contain multiple disconnected pinna components.")

        bridged = bridge_dissected_rachis_components(severed)
        n_after, _ = cv2.connectedComponents((bridged > 0).astype(np.uint8))
        self.assertEqual(n_after, 2, "All valid pinna components must be unified into 1 foreground component.")

        cnts_after, _ = cv2.findContours(bridged, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        main_after = max(cnts_after, key=cv2.contourArea)
        # Ensure total contour area retains all pinna lobes rather than just the terminal lobe
        self.assertGreater(cv2.contourArea(main_after), 0.90 * float(np.count_nonzero(severed)))

        # 3. Also verify _merge_and_bridge_dissected_crop_mask uses Class 1 intra-lamina rachis pixels
        petiole_rachis = np.zeros_like(severed)
        cv2.rectangle(petiole_rachis, (442, 260), (458, 820), 255, -1)
        merged = step02._merge_and_bridge_dissected_crop_mask(severed > 0, petiole_rachis > 0)
        n_merged, _ = cv2.connectedComponents(merged)
        self.assertEqual(n_merged, 2)

    def test_lm3_explicit_petiole_and_hole_mask_routing_on_dissected_leaf(self) -> None:
        """classify_and_route_leaf strips petiole_mask (class_id == 1) and fills hole_mask (class_id == 2) while preserving basal pinnae on dissected Packera leaves."""
        from scripts.vision.leaf_geometry_utils import classify_and_route_leaf

        mask_no_stalk = self._build_high_dpi_pinnatisect_mask(sever_rachis=False, include_petiole=False)
        mask_with_stalk = self._build_high_dpi_pinnatisect_mask(sever_rachis=False, include_petiole=True)
        lamina_area = float(np.count_nonzero(mask_no_stalk))

        petiole_mask = np.zeros_like(mask_with_stalk)
        cv2.rectangle(petiole_mask, (450 - 7, 835), (450 + 7, 1090), 255, -1)

        hole_mask = np.zeros_like(mask_with_stalk)
        cv2.circle(hole_mask, (450, 210), 18, 255, -1)
        raw_dissected = mask_with_stalk.copy()
        raw_dissected[hole_mask > 0] = 0

        routed = classify_and_route_leaf(
            raw_dissected,
            petiole_mask=petiole_mask,
            hole_mask=hole_mask,
        )
        self.assertIsNotNone(routed.processed_mask)
        routed_area = float(np.count_nonzero(routed.processed_mask))

        # Processed mask must retain >= 95% of the true lamina area (all 3 lateral pinna pairs + terminal lobe)
        self.assertGreaterEqual(
            routed_area / lamina_area,
            0.95,
            f"Explicit petiole/hole routing lost dissected pinnae: retained {routed_area / lamina_area:.3f} of lamina.",
        )
        # Internal hole in terminal lobe must be filled and basal petiole stalk (y > 950) stripped
        self.assertGreater(int(routed.processed_mask[210, 450]), 0)
        self.assertEqual(int(np.count_nonzero(routed.processed_mask[950:, :])), 0)
        self.assertEqual(routed.assigned_tier, "tier1")
        self.assertTrue(routed.is_dissected)

    def test_stage1_rescues_leaf_partial_when_leaf_whole_absent(self) -> None:
        """When Stage 1 finds 0 Class 0 (leaf_whole) boxes, Class 1 (leaf_partial) boxes are rescued for LM3LeafSegmenterEngine."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            dummy_weights = tmp_p / "dummy.pt"
            dummy_weights.write_bytes(b"0")
            v_csv = tmp_p / "v.csv"
            pd.DataFrame([{"catalogNumber": "TEST_GLABELLA_1", "scientificName": "Packera glabella"}]).to_csv(
                v_csv, index=False
            )
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(
                step02, "LM3PlantDetector"
            ):
                pipe = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_weights,
                    stage1_weights=dummy_weights,
                    output_dir=tmp_p / "out",
                    device="cpu",
                    db_path=tmp_p / "test.sqlite",
                )
            sheet_img = np.zeros((2000, 1500, 3), dtype=np.uint8)
            # Case A: Only Class 1 (leaf_partial) detected on sheet -> rescued into candidate_s1_crops
            s1_only_partial = [
                ((100, 200, 600, 550), 0.72, 1),
            ]
            unpacked_a = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_only_partial, 10.0)
            candidate_crops_a = unpacked_a[0]
            self.assertEqual(len(candidate_crops_a), 1)
            self.assertEqual(candidate_crops_a[0][2], 1)

            # Case B: Both Class 0 and Class 1 detected on sheet (with default segment_partial_leaves=False)
            s1_mixed = [
                ((100, 200, 600, 550), 0.85, 0),
                ((250, 400, 550, 600), 0.65, 1),
            ]
            unpacked_b = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_mixed, 10.0)
            candidate_crops_b = unpacked_b[0]
            self.assertEqual(len(candidate_crops_b), 1)
            self.assertEqual(candidate_crops_b[0][2], 0)

    def test_lm3_working_frame_coordinate_unscaling(self) -> None:
        """
        Verify LeafMachine3 working-frame coordinate unscaling (`compute_lm3_work_scale`,
        `unscale_working_bbox_to_sheet`, and `unscale_working_points_to_sheet`) converts
        from a 3,200 px working frame back to native sheet dimensions (3753x5634).
        """
        from scripts.vision.lm3_data_loader import (
            compute_lm3_work_scale,
            unscale_working_bbox_to_sheet,
            unscale_working_points_to_sheet,
        )

        sheet_h, sheet_w = 5634, 3753
        work_scale = compute_lm3_work_scale(sheet_h, sheet_w, max_working_dim=3200)
        expected_scale = 3200.0 / 5634.0
        self.assertAlmostEqual(work_scale, expected_scale, places=6)

        # Sheets already <= 3200 px on longest side must retain work_scale == 1.0
        self.assertAlmostEqual(compute_lm3_work_scale(2400, 1800, max_working_dim=3200), 1.0, places=6)

        # Working-frame bbox -> native sheet bbox (ymin, xmin, ymax, xmax)
        work_bbox = (1600, 1000, 2400, 1500)
        sheet_bbox = unscale_working_bbox_to_sheet(
            work_bbox,
            work_scale,
            sheet_shape=(sheet_h, sheet_w),
        )
        self.assertEqual(len(sheet_bbox), 4)
        ymin, xmin, ymax, xmax = sheet_bbox
        self.assertAlmostEqual(ymin, int(round(1600 / expected_scale)), delta=1)
        self.assertAlmostEqual(xmin, int(round(1000 / expected_scale)), delta=1)
        self.assertAlmostEqual(ymax, int(round(2400 / expected_scale)), delta=1)
        self.assertAlmostEqual(xmax, int(round(1500 / expected_scale)), delta=1)
        self.assertGreaterEqual(ymin, 0)
        self.assertGreaterEqual(xmin, 0)
        self.assertLessEqual(ymax, sheet_h)
        self.assertLessEqual(xmax, sheet_w)

        # Working-frame 2D points / 31-keypoint array -> native sheet coordinates
        work_pts = np.array(
            [
                [1000.0 * expected_scale, 2000.0 * expected_scale, 0.95],
                [1500.0 * expected_scale, 3500.0 * expected_scale, 0.91],
            ],
            dtype=np.float32,
        )
        unscaled_pts = unscale_working_points_to_sheet(
            work_pts,
            work_scale,
            sheet_shape=(sheet_h, sheet_w),
        )
        self.assertEqual(unscaled_pts.shape, (2, 3))
        self.assertAlmostEqual(float(unscaled_pts[0, 0]), 1000.0, delta=1.5)
        self.assertAlmostEqual(float(unscaled_pts[0, 1]), 2000.0, delta=1.5)
        self.assertAlmostEqual(float(unscaled_pts[0, 2]), 0.95, places=4)
        self.assertAlmostEqual(float(unscaled_pts[1, 0]), 1500.0, delta=1.5)
        self.assertAlmostEqual(float(unscaled_pts[1, 1]), 3500.0, delta=1.5)

    def test_lm3_three_class_leaf_petiole_hole_routing_and_31kp_reflection(self) -> None:
        """
        Verify LeafMachine3 3-class `Leaf` / `Petiole` / `Hole` mask routing in
        `classify_and_route_leaf` (petiole_mask subtraction and internal hole_mask filling)
        and 31-keypoint midvein fold reflection (`detect_folded_leaf` and `reflect_folded_mask`).
        """
        from scripts.vision.leaf_geometry_utils import (
            classify_and_route_leaf,
            detect_folded_leaf,
            reflect_folded_mask,
        )

        # 1. 3-class Leaf / Petiole / Hole routing
        canvas_shape = (400, 400)
        pure_lamina = np.zeros(canvas_shape, dtype=np.uint8)
        cv2.ellipse(pure_lamina, (200, 180), (65, 110), 0.0, 0, 360, 255, -1)

        petiole_mask = np.zeros(canvas_shape, dtype=np.uint8)
        cv2.rectangle(petiole_mask, (194, 285), (206, 375), 255, -1)

        hole_mask = np.zeros(canvas_shape, dtype=np.uint8)
        cv2.circle(hole_mask, (200, 175), 18, 255, -1)

        # Combined raw leaf mask with attached petiole stalk and internal herbivory hole
        raw_leaf_mask = cv2.bitwise_or(pure_lamina, petiole_mask)
        raw_leaf_mask[hole_mask > 0] = 0

        routed = classify_and_route_leaf(
            raw_leaf_mask,
            petiole_mask=petiole_mask,
            hole_mask=hole_mask,
        )
        self.assertEqual(routed.assigned_tier, "tier1")
        self.assertIsNotNone(routed.processed_mask)
        # Internal hole at (200, 175) must be filled in processed_mask
        self.assertGreater(int(routed.processed_mask[175, 200]), 0)
        # External petiole stalk at y=350 must be subtracted from processed_mask
        self.assertEqual(int(np.count_nonzero(routed.processed_mask[330:, :])), 0)

        # 2. 31-keypoint midvein fold detection & bilateral reflection
        folded_hemi = np.zeros(canvas_shape, dtype=np.uint8)
        # Right half-blade only (x >= 200), folded directly along vertical midvein x = 200
        cv2.ellipse(folded_hemi, (200, 200), (60, 120), 0.0, -90, 90, 255, -1)

        keypoints_31 = np.zeros((31, 3), dtype=np.float32)
        keypoints_31[0] = [200.0, 80.0, 0.98]   # 0: lamina_tip
        keypoints_31[1] = [200.0, 320.0, 0.98]  # 1: lamina_base
        for i in range(6):
            # 2..7: petiole_1..6
            keypoints_31[2 + i] = [200.0, 325.0 + i * 8.0, 0.92]
        for i in range(23):
            # 8..30: midvein_1..23 along the fold axis x = 200
            y_mv = 320.0 - (i / 22.0) * 240.0
            keypoints_31[8 + i] = [200.0, y_mv, 0.96]

        is_folded = detect_folded_leaf(
            folded_hemi,
            p_apex=(200, 80),
            p_base=(200, 320),
            keypoints_31=keypoints_31,
        )
        self.assertTrue(is_folded, "Hemi-blade with unilateral lamina along 31-kp midvein must be detected as folded.")

        reflected = reflect_folded_mask(
            folded_hemi,
            p_apex=(200, 80),
            p_base=(200, 320),
            keypoints_31=keypoints_31,
        )
        self.assertIsNotNone(reflected)
        hemi_area = float(np.count_nonzero(folded_hemi))
        refl_area = float(np.count_nonzero(reflected))
        self.assertGreater(refl_area, 1.75 * hemi_area)
        # Both left (x < 200) and right (x > 200) lamina halves must now be populated symmetrically
        left_area = float(np.count_nonzero(reflected[:, :200]))
        right_area = float(np.count_nonzero(reflected[:, 200:]))
        self.assertAlmostEqual(left_area / max(right_area, 1.0), 1.0, delta=0.15)

    def test_lm3_unetpp_specimen_mask_rachis_bridging_and_6pt_petiole_pose(self) -> None:
        """
        Verify UNet++ v9 `specimen_paper_cleaned_mask` guidance in `bridge_dissected_rachis_components`
        and LM3 6-point petiole pose skeleton (`petiole_1..6`, `lamina_base`) handling in
        `estimate_rosette_crown_point` and `approximate_petiole_and_scape_metrics`.
        """
        from scripts.vision.leaf_geometry_utils import bridge_dissected_rachis_components

        # 1. UNet++ v9 specimen_paper_cleaned_mask rachis bridging
        severed = self._build_high_dpi_pinnatisect_mask(sever_rachis=True)
        unet_specimen_mask = severed.copy()
        # UNet++ v9 specimen_segmenter preserves the central rachis tissue across the severed gaps
        cv2.rectangle(unet_specimen_mask, (442, 195), (458, 880), 255, -1)

        bridged = bridge_dissected_rachis_components(
            severed,
            specimen_paper_cleaned_mask=unet_specimen_mask,
        )
        n_comp, _ = cv2.connectedComponents((bridged > 0).astype(np.uint8))
        self.assertEqual(
            n_comp,
            2,
            "UNet++ specimen_paper_cleaned_mask rachis corridor must bridge severed pinna islands into 1 component.",
        )

        # 2. LM3 6-point petiole pose skeleton (`petiole_1..6`, `lamina_base`) in crown & petiole metrics
        sheet_h, sheet_w = 1000, 800
        sheet_img = np.full((sheet_h, sheet_w, 3), 230, dtype=np.uint8)
        true_crown_xy = (400.0, 820.0)

        leaf_records = []
        for idx, angle_deg in enumerate([-25.0, 0.0, 25.0], start=1):
            rad = np.deg2rad(angle_deg)
            dx = np.sin(rad)
            dy = -np.cos(rad)
            # Lamina base 90 px outward from crown; tip 220 px outward from crown
            bx = true_crown_xy[0] + 90.0 * dx
            by = true_crown_xy[1] + 90.0 * dy
            tx = true_crown_xy[0] + 220.0 * dx
            ty = true_crown_xy[1] + 220.0 * dy

            kp31 = np.zeros((31, 3), dtype=np.float32)
            kp31[0] = [tx, ty, 0.96]  # lamina_tip
            kp31[1] = [bx, by, 0.96]  # lamina_base
            for p_i in range(6):
                # petiole_1 (near lamina_base) -> petiole_6 (near rosette crown)
                frac = (p_i + 1) / 6.0
                px = bx + frac * (true_crown_xy[0] - bx)
                py = by + frac * (true_crown_xy[1] - by)
                kp31[2 + p_i] = [px, py, 0.94]

            m_sheet = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cx_i, cy_i = int(round((bx + tx) * 0.5)), int(round((by + ty) * 0.5))
            cv2.circle(m_sheet, (cx_i, cy_i), 35, 255, -1)

            leaf_records.append({
                "catalogNumber": "NCU_POSE31",
                "plant_individual_id": 0,
                "leaf_id": idx,
                "bbox": (cy_i - 40, cx_i - 40, cy_i + 40, cx_i + 40),
                "mask": m_sheet,
                "p_base": (int(round(bx)), int(round(by))),
                "p_apex": (int(round(tx)), int(round(ty))),
                "keypoints_31": kp31,
                "pixels_per_mm": 10.0,
                "is_calibrated_mm": True,
            })

        crown_pt = step02.estimate_rosette_crown_point(
            leaf_records,
            sheet_shape=(sheet_h, sheet_w),
        )
        self.assertIsNotNone(crown_pt)
        dist_crown = float(np.hypot(crown_pt[0] - true_crown_xy[0], crown_pt[1] - true_crown_xy[1]))
        self.assertLessEqual(
            dist_crown,
            35.0,
            f"Estimated crown {crown_pt} with 6-pt petiole pose skeleton should converge near {true_crown_xy}.",
        )

        pet_metrics = step02.approximate_petiole_and_scape_metrics(
            leaf_records=leaf_records,
            reproductive_boxes=[],
            sheet_bgr=sheet_img,
            pixels_per_mm=10.0,
        )
        self.assertIn("petiole_length_px", pet_metrics)
        self.assertGreater(float(pet_metrics["petiole_length_px"]), 50.0)
        self.assertGreater(float(pet_metrics["petiole_length_mm"]), 5.0)

    def test_dissected_sinus_topology_opening_and_canonical_tip_up_ect_export(self) -> None:
        """
        Verify topology-aware inter-lobe sinus opening vs interior hole filling,
        canonical tip-up contour rotation in export_standardized_contour,
        LM3ScaleCalibrator megapixel prior, orient_and_score_leaf_lm3, and
        compute_leaf_ect_descriptor.
        """
        from scripts.vision.leaf_geometry_utils import (
            classify_and_route_leaf,
            export_standardized_contour,
        )
        from scripts.vision.lm3_inference_engines import (
            LM3ScaleCalibrator,
            compute_leaf_ect_descriptor,
            orient_and_score_leaf_lm3,
        )

        # 1. Dissected leaf with an enclosed inter-lobe sinus (touching pinna tips)
        #    and a true interior insect hole deep inside the terminal lobe.
        mask = self._build_high_dpi_pinnatisect_mask(sever_rachis=False, include_petiole=False)
        # Bridge the outer tips of the upper right pinna (y=410) and middle right pinna (y=610)
        # with a thin 2-px touch at x=565 so the sinus between y=450..570 becomes an enclosed hole_mask region.
        cv2.line(mask, (565, 410), (565, 610), 255, thickness=2)
        hole_mask = np.zeros_like(mask)
        # Mark the enclosed inter-lobe sinus (x=475..562, y=455..565) in hole_mask
        cv2.rectangle(hole_mask, (475, 455), (562, 565), 255, -1)
        mask[hole_mask > 0] = 0
        # Also carve a small true interior insect hole deep inside the terminal lobe at (450, 210)
        cv2.circle(hole_mask, (450, 210), 14, 255, -1)
        cv2.circle(mask, (450, 210), 14, 0, -1)

        routed = classify_and_route_leaf(mask, hole_mask=hole_mask)
        self.assertIsNotNone(routed.processed_mask)
        # True interior insect hole at (450, 210) is filled (> 0)
        self.assertGreater(int(routed.processed_mask[210, 450]), 0)
        # Enclosed inter-lobe sinus at (515, 510) remains open (== 0) and connected to exterior background
        self.assertEqual(int(routed.processed_mask[510, 515]), 0)
        ext_cnts, _ = cv2.findContours(routed.processed_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        main_ext_area = float(cv2.contourArea(max(ext_cnts, key=cv2.contourArea)))
        # Because the touching-pinna bottleneck was opened, RETR_EXTERNAL area matches foreground pixel area closely
        fg_area = float(np.count_nonzero(routed.processed_mask))
        self.assertLess(abs(main_ext_area - fg_area) / max(fg_area, 1.0), 0.05)

        # 2. Canonical tip-up contour rotation on an oblique 135-degree leaf
        oblique_mask = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(oblique_mask, (200, 200), (95, 38), 45.0, 0, 360, 255, -1)
        kp31 = np.zeros((31, 3), dtype=np.float32)
        # Base at top-left (133, 133), tip at bottom-right (267, 267) -> pointing down-right (+135 deg from tip-up)
        kp31[0] = [267.0, 267.0, 0.98]  # lamina_tip
        kp31[1] = [133.0, 133.0, 0.98]  # lamina_base
        for i in range(23):
            frac = i / 22.0
            kp31[8 + i] = [133.0 + frac * 134.0, 133.0 + frac * 134.0, 0.95]

        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            c_path = export_standardized_contour(
                oblique_mask,
                catalog_number="TEST_OBLIQUE_135",
                leaf_id=1,
                output_dir=tmp_p,
                num_points=120,
                keypoints_31=kp31,
            )
            self.assertIsNotNone(c_path)
            cdf = pd.read_csv(c_path)
            # Point 0 (lamina_base) must sit at the bottom (high y_norm), and point 60 (lamina_tip) at the top (low y_norm)
            self.assertGreater(float(cdf.loc[0, "y_norm"]), 0.65)
            self.assertLess(float(cdf.loc[60, "y_norm"]), 0.15)
            span_y = float(cdf["y_norm"].max() - cdf["y_norm"].min())
            span_x = float(cdf["x_norm"].max() - cdf["x_norm"].min())
            self.assertGreater(span_y, 1.8 * span_x)

            # 3. Module 18 ECT descriptor export
            ect_file = tmp_p / "TEST_OBLIQUE_135_leaf1_ect.npy"
            ect_arr = compute_leaf_ect_descriptor(
                oblique_mask,
                output_path=ect_file,
                num_dirs=32,
                num_thresh=32,
            )
            self.assertIsNotNone(ect_arr)
            self.assertTrue(ect_file.exists())
            self.assertEqual(ect_arr.shape, (32, 32))
            self.assertEqual(ect_arr.dtype, np.float32)

        # 4. LM3 native leaf orientation and bilateral archetype scoring
        orient_res = orient_and_score_leaf_lm3(oblique_mask, keypoints_31=kp31)
        self.assertTrue(orient_res["success"])
        self.assertGreater(abs(float(orient_res["angle_cw"])), 100.0)
        self.assertGreaterEqual(float(orient_res["archetype_score"]), 0.70)

        # 5. LM3ScaleCalibrator megapixel prior fallback when no ruler is present
        calibrator = LM3ScaleCalibrator(device="cpu")
        calibrator.archival_session = None
        sheet_dummy = np.zeros((3000, 2000, 3), dtype=np.uint8)  # 6.0 MP
        cal_res = calibrator.calibrate_sheet(sheet_dummy, catalog_number="NCU_MP_FALLBACK")
        expected_px_per_mm = 27.112084 * np.sqrt(6.0) / 10.0
        self.assertAlmostEqual(float(cal_res["pixels_per_mm"]), expected_px_per_mm, places=3)
        self.assertEqual(cal_res["cf_source"], "predicted_from_megapixels")

    def test_native_lm3_stride32_letterbox_tight_crops_and_cordate_survival(self) -> None:
        """
        Verify Phases 1-4 native LeafMachine3 alignment:
          1. Dynamic stride-32 letterboxing (`auto=True, stride=32`) matches `ultra_rep.letterbox`
             instead of forcing a square canvas.
          2. Default `crop_leaf_bbox_with_padding` uses tight `pad_ratio=0.0` matching `detector_common.py`.
          3. `LM3LeafSegmenterEngine.segment_leaf_crop` unions multiple `Petiole` and `Hole` instances (`>= 64 px`).
          4. Stage 1 `_deduplicate_sheet_boxes` uses `overlap=0.95` (`inter / min_area >= 0.95`) so adjacent
             rosette leaves with moderate overlap (`IoU ~ 0.55`) both survive.
          5. Broad cordate/orbicular leaves (`W/L > 0.78, solidity < 0.84`) survive high-DPI QC.
        """
        from scripts.vision.lm3_inference_engines import (
            LM3LeafSegmenterEngine,
            LM3PlantDetector,
            _letterbox_bgr,
            _letterbox_gain_pad_hw,
            crop_leaf_bbox_with_padding,
        )

        # 1. Dynamic stride-32 letterboxing on a 600x200 crop at img_size=1024
        dummy_crop = np.full((600, 200, 3), 180, dtype=np.uint8)
        lb_bgr, gain, (pad_x, pad_y) = _letterbox_bgr(dummy_crop, 1024)
        self.assertEqual(lb_bgr.shape[0], 1024)
        self.assertLess(lb_bgr.shape[1], 1024)
        self.assertEqual(lb_bgr.shape[1] % 32, 0)
        gain_hw, (pad_hw_x, pad_hw_y) = _letterbox_gain_pad_hw(lb_bgr.shape[:2], dummy_crop.shape[:2])
        self.assertAlmostEqual(gain_hw, gain, places=6)
        self.assertAlmostEqual(pad_hw_x, pad_x, places=6)
        self.assertAlmostEqual(pad_hw_y, pad_y, places=6)

        # 2. Default tight crop (`pad_ratio=0.0`)
        sheet_img = np.full((1000, 800, 3), 235, dtype=np.uint8)
        crop_tight, box_tight = crop_leaf_bbox_with_padding(sheet_img, (200, 150, 500, 450))
        self.assertEqual(box_tight, (200, 150, 500, 450))
        self.assertEqual(crop_tight.shape[:2], (300, 300))

        # 3. Multi-instance Petiole and Hole unioning in LM3LeafSegmenterEngine.segment_leaf_crop
        seg_eng = LM3LeafSegmenterEngine(device="cpu")
        lam_m = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(lam_m, (100, 90), 55, 255, -1)
        pet_1 = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(pet_1, (96, 145), (104, 165), 255, -1)  # 9 * 21 = 189 px >= 64
        pet_2 = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(pet_2, (96, 165), (104, 190), 255, -1)  # 9 * 26 = 234 px >= 64
        hole_1 = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(hole_1, (85, 85), 6, 255, -1)  # ~113 px >= 64
        hole_2 = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(hole_2, (115, 95), 6, 255, -1)  # ~113 px >= 64
        fake_preds = [
            ((35, 45, 145, 155), lam_m, 0.94, 0),
            ((145, 96, 165, 104), pet_1, 0.85, 1),
            ((165, 96, 190, 104), pet_2, 0.80, 1),
            ((79, 79, 91, 91), hole_1, 0.82, 2),
            ((89, 109, 101, 121), hole_2, 0.79, 2),
        ]
        with patch.object(seg_eng, "predict", return_value=fake_preds):
            lam_out, pet_out, hole_out, sc_out = seg_eng.segment_leaf_crop(np.zeros((200, 200, 3), dtype=np.uint8))
        self.assertIsNotNone(lam_out)
        self.assertIsNotNone(pet_out)
        self.assertIsNotNone(hole_out)
        self.assertGreater(int(np.count_nonzero(pet_out)), int(np.count_nonzero(pet_1)))
        self.assertGreater(int(np.count_nonzero(hole_out)), int(np.count_nonzero(hole_1)))
        self.assertAlmostEqual(sc_out, 0.94, places=3)

        # 4. Stage 1 overlap=0.95 deduplication keeps moderately overlapping rosette leaves
        boxes_s1 = [
            ((500, 500, 900, 900), 0.95, 0),
            ((550, 550, 950, 950), 0.92, 0),  # ~62% containment (< 0.95) -> both kept!
            ((510, 510, 890, 890), 0.88, 0),  # 100% contained (>= 0.95) -> suppressed!
        ]
        det_eng = LM3PlantDetector(device="cpu")
        deduped = det_eng._deduplicate_sheet_boxes(boxes_s1)
        self.assertEqual(len(deduped), 2)

    def test_step02_cli_compute_flags_and_back_to_back_minibatch_halving(self) -> None:
        """
        Verify CLI flags (--threads/--workers, --gpu-batch-size/--batch-size, --wave-size,
        --prefetch-workers, --prep-workers, --post-workers) and _run_onnx_minibatches_back_to_back
        automatic mini-batch halving on transient memory errors.
        """
        import main as main_cli
        from scripts.vision.lm3_inference_engines import _run_onnx_minibatches_back_to_back

        with patch(
            "sys.argv",
            [
                "02_segment_and_extract.py",
                "--workers",
                "32",
                "--batch-size",
                "16",
                "--wave-size",
                "24",
                "--prefetch-workers",
                "16",
                "--prep-workers",
                "24",
                "--post-workers",
                "32",
            ],
        ):
            parsed_s2 = step02.parse_args()
            self.assertEqual(parsed_s2.threads, 32)
            self.assertEqual(parsed_s2.gpu_batch_size, 16)
            self.assertEqual(parsed_s2.wave_size, 24)
            self.assertEqual(parsed_s2.prefetch_workers, 16)
            self.assertEqual(parsed_s2.prep_workers, 24)
            self.assertEqual(parsed_s2.post_workers, 32)

        main_parser = main_cli.build_parser()
        parsed_main = main_parser.parse_args(
            [
                "segment",
                "--workers",
                "64",
                "--gpu-batch-size",
                "16",
                "--wave-size",
                "32",
                "--prefetch-workers",
                "24",
                "--prep-workers",
                "32",
                "--post-workers",
                "48",
            ]
        )
        self.assertEqual(parsed_main.threads, 64)
        self.assertEqual(parsed_main.gpu_batch_size, 16)
        self.assertEqual(parsed_main.wave_size, 32)
        self.assertEqual(parsed_main.prefetch_workers, 24)
        self.assertEqual(parsed_main.prep_workers, 32)
        self.assertEqual(parsed_main.post_workers, 48)

        # Verify _run_onnx_minibatches_back_to_back automatically halves sub_bs when B>2 raises an OOM error
        class _FakeBFCExhaustedSession:
            def __init__(self) -> None:
                self.batch_sizes_seen: list[int] = []

            def run(self, _out_names, feed_dict):
                x = feed_dict["images"]
                b = int(x.shape[0])
                self.batch_sizes_seen.append(b)
                if b > 2:
                    raise RuntimeError("BFCArena::AllocateRawInternal Failed to allocate memory")
                return [np.zeros((b, 8, 10), dtype=np.float32)]

        fake_sess = _FakeBFCExhaustedSession()
        tensors = [np.zeros((3, 8, 8), dtype=np.float32) for _ in range(6)]
        records = _run_onnx_minibatches_back_to_back(fake_sess, "images", tensors, initial_sub_bs=8)
        # Initial attempt at B=6 fails -> halves to 3 (fails) -> halves to 1 (succeeds for all 6 items)
        total_covered = sum(end_i - start_i for start_i, end_i, _ in records)
        self.assertEqual(total_covered, 6)
        self.assertGreaterEqual(len(fake_sess.batch_sizes_seen), 3)

    def test_step02_figures_eastern_na_refinements(self) -> None:
        """
        Verify Phase 2 figure refinements in step02_figures.py:
        1. orientation="horizontal" is used for Matplotlib 3.11 boxplot compatibility.
        2. bbox_to_anchor=(0.5, 1.005) places the Panel B2 legend above the horizontal bars.
        3. "20th %ile (Compact)" is used for non-Type Column 0 headers.
        4. Small-N species (n=3 leaves) do not duplicate leaf IDs across a montage row.
        """
        from pathlib import Path
        from unittest.mock import patch
        import scripts.vision.step02_figures as s02_figs

        src_text = Path(s02_figs.__file__).read_text(encoding="utf-8")
        self.assertIn('orientation="horizontal"', src_text)
        self.assertIn("bbox_to_anchor=(0.5, 1.005)", src_text)
        self.assertIn('"20th %ile (Compact)"', src_text)

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_dir = Path(tmp_d)
            mask_paths = []
            for idx_l in (1, 2, 3):
                m_img = np.zeros((160, 120), dtype=np.uint8)
                cv2.ellipse(m_img, (60, 80), (22 + idx_l * 3, 52), 25, 0, 360, 255, -1)
                m_file = tmp_dir / f"NCU000{idx_l}_leaf{idx_l}.png"
                cv2.imwrite(str(m_file), m_img)
                mask_paths.append(str(m_file))

            df_small = pd.DataFrame(
                {
                    "catalogNumber": ["NCU0001", "NCU0002", "NCU0003"],
                    "leaf_id": [1, 2, 3],
                    "species_standardized": ["Packera millefolia"] * 3,
                    "assigned_tier": ["Tier_1_Direct"] * 3,
                    "determiner_tier": ["Tier_1_Gold"] * 3,
                    "is_type_specimen": [False, False, False],
                    "solidity": [0.86, 0.88, 0.90],
                    "aspect_ratio": [1.65, 2.05, 2.55],
                    "ucs_score": [0.93, 0.94, 0.95],
                    "pixels_per_mm": [10.0, 10.0, 10.0],
                    "is_calibrated_mm": [True, True, True],
                    "mask_path": mask_paths,
                }
            )

            captured_titles: list[str] = []
            orig_save = s02_figs._atomic_save_matplotlib_fig

            def _spy_save(fig, out_base_path, dpi=300):
                if Path(out_base_path).name == "fig02b_silhouette_montage":
                    captured_titles.extend(ax.get_title() for ax in fig.axes)
                return orig_save(fig, out_base_path, dpi=dpi)

            with patch.object(s02_figs, "_atomic_save_matplotlib_fig", side_effect=_spy_save), patch.object(
                step02, "_atomic_save_matplotlib_fig", side_effect=_spy_save
            ):
                out_paths = s02_figs.generate_step02_publication_figures(
                    extracted_df=df_small,
                    failed_df=pd.DataFrame(),
                    repro_df=pd.DataFrame(),
                    pheno_df=pd.DataFrame(),
                    figures_dir=tmp_dir / "figures",
                )

            self.assertIn("fig02b_png", out_paths)
            self.assertEqual(len(captured_titles), 4)
            self.assertTrue(captured_titles[0].startswith("20th %ile (Compact)\n"))
            leaf_titles = [t for t in captured_titles if "NCU000" in t]
            self.assertEqual(len(leaf_titles), 3, "Exactly 3 tiles should render leaf silhouettes for n=3.")
            seen_vouchers = [t.split("\n")[-1].split(" ")[0] for t in leaf_titles]
            self.assertEqual(len(seen_vouchers), len(set(seen_vouchers)), "No duplicate leaf IDs across a row.")


if __name__ == "__main__":
    unittest.main()
