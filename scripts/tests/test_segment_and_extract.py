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
        self.assertEqual(Path(saved_path).parent.name, "tier1_intact")
        self.assertFalse((self.temp_dir / "masks" / "TEST001_leaf1.png").exists())

        read_mask = cv2.imread(saved_path, cv2.IMREAD_GRAYSCALE)
        self.assertIsNotNone(read_mask)
        self.assertEqual(np.count_nonzero(read_mask > 128), np.count_nonzero(mask > 0))

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
            ((550, 400, 750, 600), 0.89, 1),  # leaf_partial -> routed to Stage 2 when segment_partial_leaves=True
            ((100, 50, 150, 90), 0.80, 2),    # leaflet -> ignored
            ((10, 10, 990, 790), 0.99, 8),    # specimen -> ignored
            ((800, 100, 950, 140), 0.77, 10), # wood -> ignored
        ]

        def fake_stage2_predict(crop_img: np.ndarray):
            ch, cw = crop_img.shape[:2]
            m = np.zeros((ch, cw), dtype=np.uint8)
            cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.24), int(ch * 0.40)), 0, 0, 360, 255, -1)
            return [((10, 10, ch - 10, cw - 10), m, 0.96, 0)]

        prev_seg_partial = getattr(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", True)
        try:
            with patch("cv2.imread", return_value=sheet_img):
                with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                    with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_stage2_predict) as mock_s2:
                        records, failure = self.pipeline.process_voucher("NCU123", self.img_path)
        finally:
            object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", prev_seg_partial)

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

        # Verify Tier 1 binary leaf masks exist in data/masks/tier1_intact/ without duplicate top-level copies
        for leaf_id in [1, 2, 3]:
            leaf_mask = self.pipeline.masks_dir / "tier1_intact" / f"NCU123_leaf{leaf_id}.png"
            self.assertTrue(leaf_mask.exists(), f"Tier 1 binary leaf mask {leaf_mask} must exist.")
            self.assertFalse((self.pipeline.masks_dir / f"NCU123_leaf{leaf_id}.png").exists())

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
        from p_apex (broad lamina tip) regardless of cv2.fitLine sign, and is_botanical_dissection()
        rejects unilateral bites.
        """
        from scripts.vision.leaf_geometry_utils import is_botanical_dissection

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

        # 2. Unilateral bite on low-solidity leaf must NOT be accepted as dissected
        bitten = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(bitten, (200, 200), (65, 110), 0.0, 0, 360, 255, -1)
        cv2.ellipse(bitten, (145, 200), (58, 65), 0.0, 0, 360, 0, -1)
        self.assertFalse(is_botanical_dissection(bitten))

    def test_stage1_high_dpi_filtering_and_single_primary_blade_selection(self) -> None:
        """Verify on high-DPI herbarium sheets (>2000px), Stage 1 forwards all Class 0 (`Leaf_WHOLE`) boxes to Stage 2 while skipping Class 1 (`Leaf_PARTIAL`)."""
        sheet_img = np.full((3000, 2400, 3), 230, dtype=np.uint8)
        self.assertTrue(max(sheet_img.shape[:2]) > 2000)

        stage1_dets = [
            ((1500, 500, 1900, 900), 0.90, 1),   # leaf_partial -> skipped
            ((400, 800, 800, 1200), 0.92, 0),    # upper-sheet Class 0 leaf_whole -> forwarded
            ((1600, 600, 1680, 680), 0.88, 0),   # 80x80 Class 0 box -> forwarded to Stage 2
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

        prev_seg_partial = getattr(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        try:
            with patch("cv2.imread", return_value=sheet_img):
                with patch.object(self.pipeline.stage1_detector, "predict", return_value=stage1_dets):
                    with patch.object(self.pipeline.model_engine, "predict", side_effect=fake_high_dpi_stage2) as mock_s2:
                        records, failure = self.pipeline.process_voucher("NCU_HIGHDPI", self.img_path)
        finally:
            object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", prev_seg_partial)

        self.assertIsNone(failure)
        # Stage 2 LM3LeafSegmenterEngine is called for all 3 Class 0 (`Leaf_WHOLE`) boxes (Class 1 skipped)
        self.assertEqual(mock_s2.call_count, 3)
        self.assertEqual(len(records), 3)
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

        prev_seg_partial = getattr(self.pipeline.config.segmentation, "segment_partial_leaves", True)
        object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", False)
        try:
            with patch("cv2.imread", return_value=sheet_img):
                with patch.object(self.pipeline.stage1_detector, "predict", return_value=s1_dets):
                    with patch.object(self.pipeline.model_engine, "predict_batch", side_effect=fake_raw_lm3_stage2) as mock_s2:
                        records, failure = self.pipeline.process_voucher("NCU_RAW_LM3", self.img_path)
        finally:
            object.__setattr__(self.pipeline.config.segmentation, "segment_partial_leaves", prev_seg_partial)

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
        """High-DPI dissected leaves (minor_w > 360 px) must pass is_botanical_dissection."""
        from scripts.vision.leaf_geometry_utils import is_botanical_dissection

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

            # Case B: Both Class 0 and Class 1 detected on sheet
            # With default segment_partial_leaves=True, both Class 0 and Class 1 are included;
            # when segment_partial_leaves=False is explicitly set, only Class 0 is included.
            s1_mixed = [
                ((100, 200, 600, 550), 0.85, 0),
                ((700, 400, 1000, 600), 0.65, 1),
            ]
            unpacked_b_default = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_mixed, 10.0)
            candidate_crops_b_default = unpacked_b_default[0]
            if getattr(pipe.config.segmentation, "segment_partial_leaves", False):
                self.assertEqual(len(candidate_crops_b_default), 2)
                self.assertEqual({c[2] for c in candidate_crops_b_default}, {0, 1})

            object.__setattr__(pipe.config.segmentation, "segment_partial_leaves", False)
            unpacked_b = pipe._unpack_stage1_candidates("TEST_GLABELLA_1", sheet_img, s1_mixed, 10.0)
            candidate_crops_b = unpacked_b[0]
            self.assertEqual(len(candidate_crops_b), 1)
            self.assertEqual(candidate_crops_b[0][2], 0)

    def test_lm3_working_frame_coordinate_unscaling(self) -> None:
        """
        Verify LeafMachine3 working-frame coordinate unscaling (`compute_lm3_work_scale`
        and `unscale_working_bbox_to_sheet`) converts from a 3,200 px working frame back
        to native sheet dimensions (3753x5634).
        """
        from scripts.vision.lm3_data_loader import (
            compute_lm3_work_scale,
            unscale_working_bbox_to_sheet,
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

    def test_lm3_unetpp_specimen_mask_rachis_bridging_and_6pt_petiole_pose(self) -> None:
        """
        Verify LM3 6-point petiole pose skeleton (`petiole_1..6`, `lamina_base`) handling in
        `estimate_rosette_crown_point` and `approximate_petiole_and_scape_metrics`.
        """
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
            _prepare_lamina_mask,
            export_standardized_contour,
        )
        from scripts.vision.lm3_inference_engines import (
            LM3ScaleCalibrator,
            compute_leaf_ect_descriptor,
            orient_and_score_leaf_lm3,
        )

        # 1. Native LeafMachine3 lamina = silhouette & ~hole component gating (min_area_px = 64):
        #    qualifying hole_mask components (>= 64 px) are subtracted, tiny specks (< 64 px) are ignored,
        #    and zero 3px cv2.line trenches are carved into the outer margin.
        mask = self._build_high_dpi_pinnatisect_mask(sever_rachis=False, include_petiole=False)
        cv2.line(mask, (565, 410), (565, 610), 255, thickness=2)
        hole_mask = np.zeros_like(mask)
        # Qualifying large hole_mask component (x=475..562, y=455..565; area >= 64 px)
        cv2.rectangle(hole_mask, (475, 455), (562, 565), 255, -1)
        # Tiny 3x3 = 9 px speck (< 64 px) at (450, 210) deep inside the terminal lobe
        hole_mask[209:212, 449:452] = 255

        prep_mask = _prepare_lamina_mask(mask.copy(), hole_mask=hole_mask)
        self.assertIsNotNone(prep_mask)
        # Tiny 3x3 speck (< 64 px) at (450, 210) is ignored by min_area_px = 64 so lamina remains > 0
        self.assertGreater(int(prep_mask[210, 450]), 0)
        # Qualifying hole_mask component (>= 64 px) at (515, 510) is subtracted (== 0)
        self.assertEqual(int(prep_mask[510, 515]), 0)
        # Zero 3px cv2.line trenches are carved into the outer margin at x=565
        np.testing.assert_array_equal(prep_mask[410:611, 564:567] > 0, mask[410:611, 564:567] > 0)

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
        cv2.rectangle(pet_1, (96, 145), (104, 168), 255, -1)  # 9 * 24 = 216 px >= 64
        pet_2 = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(pet_2, (96, 165), (104, 190), 255, -1)  # overlaps pet_1 across rows 165..168
        hole_1 = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(hole_1, (85, 85), 6, 255, -1)  # ~113 px >= 64
        hole_2 = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(hole_2, (115, 95), 6, 255, -1)  # ~113 px >= 64
        fake_preds = [
            ((35, 45, 145, 155), lam_m, 0.94, 0),
            ((145, 96, 168, 104), pet_1, 0.85, 1),
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

        # Confidence-ordered deduplication keeps high-confidence tight box and drops lower-confidence oversized box
        boxes_nested_conf = [
            ((480, 480, 920, 920), 0.45, 0),
            ((520, 520, 880, 880), 0.96, 0),
        ]
        deduped_conf = det_eng._deduplicate_sheet_boxes(boxes_nested_conf)
        self.assertEqual(len(deduped_conf), 1)
        self.assertAlmostEqual(deduped_conf[0][1], 0.96, places=3)
        self.assertEqual(deduped_conf[0][0], (520, 520, 880, 880))

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

    def test_enclosed_sinus_hole_subtraction_and_bridge_opening(self) -> None:
        """
        Verify that _prepare_lamina_mask() applies native
        LeafMachine3 lamina = silhouette & ~hole subtraction for qualifying hole_mask
        components (>= 64 px), ignores tiny specks (< 64 px), and carves zero synthetic
        3px cv2.line trenches across the exterior wall.
        """
        from scripts.vision.leaf_geometry_utils import _prepare_lamina_mask

        mask = np.zeros((200, 160), dtype=np.uint8)
        cv2.ellipse(mask, (80, 100), (45, 80), 0, 0, 360, 255, -1)

        hole_mask = np.zeros_like(mask)
        # Qualifying hole_mask component >= 64 px (30 x 36 = 1080 px)
        hole_mask[85:115, 85:121] = 255
        # Tiny speck < 64 px (4 x 4 = 16 px) deep inside the blade
        hole_mask[60:64, 75:79] = 255

        prep_mask = _prepare_lamina_mask(mask.copy(), hole_mask=hole_mask)
        self.assertEqual(int(prep_mask[100, 100]), 0, "Qualifying hole_mask component (>= 64 px) must be subtracted.")
        self.assertEqual(int(prep_mask[62, 77]), 255, "Tiny hole_mask speck (< 64 px) must be ignored.")
        np.testing.assert_array_equal(
            prep_mask[:, 121:],
            mask[:, 121:],
            "Zero synthetic 3px cv2.line trenches may be carved across the right exterior wall.",
        )

    def test_native_lm3_delegation_dissected_pcd_and_anatomical_leaf_zone(self) -> None:
        """
        Verify:
          1. Native LeafMachine3 delegation:
             - normalize_lm3_keypoints_dict maps all 15 1-based midvein_1..15 points to distinct
               0-based midvein_0..14 points and all 5 petiole_1..5 points to distinct petiole_0..4
               points without key collisions.
             - compute_leaf_ect_descriptor invokes leafmachine3.core.ect_compute.compute_ect and
               returns a valid (32, 32) np.float32 matrix.
             - _extract_native_lm3_petiole_anchor and _estimate_petiole_width_px delegate to
               leafmachine3.core.petiole.measure_petiole and leafmachine3.core.landmark_metrics.compute_measurements.
             - compute_geometric_metrics and compute_classical_shape_descriptors delegate to
               leafmachine3.core.morphometrics.polygon_morphology.
          2. Dissected-Leaf PCD & Stage 2 enhancements:
             - Stage 1 _unpack_stage1_candidates includes Leaf_PARTIAL (cls_id == 1) dissected leaf
               boxes alongside Leaf_WHOLE (cls_id == 0) by default, and asymmetric nested-box
               deduplication keeps the enclosing parent full-leaf bounding box when a higher-confidence
               smaller terminal-lobe box is nested inside it.
             - Stage 2 _build_native_lm3_crop_masks retains secondary pinna Leaf (cls == 0) components
               (>= 64 px) when bridged to the primary lamina by the Petiole/rachis mask (cls == 1).
          3. Post-segmentation anatomical leaf_zone classification (classify_leaf_zone_anatomical):
             - Classifies basal rosette leaves near p_crown (petiolate, large area) as "basal" and
               an upper stem leaf near the capitula (dist_from_crown_norm > 0.55, sessile, smaller area)
               as "cauline" without dropping either leaf.
        """
        import inspect
        from scripts.vision import lm3_inference_engines as lm3_eng
        from scripts.vision import leaf_geometry_utils as lgu
        from scripts.vision import organ_crown_metrics as ocm
        import leafmachine3.core.ect_compute as lm3_ect
        import leafmachine3.core.petiole as lm3_petiole
        import leafmachine3.core.landmark_metrics as lm3_lmetrics
        import leafmachine3.core.morphometrics as lm3_morph

        # 1a. normalize_lm3_keypoints_dict: 15 distinct midvein_0..14 and 5 distinct petiole_0..4 without collision
        raw_kps_1based = {
            "lamina_tip": (100.0, 20.0, 0.99),
            "lamina_base": (100.0, 180.0, 0.99),
            **{f"midvein_{i}": (100.0, 20.0 + float(i * 10), 0.95) for i in range(1, 16)},
            **{f"petiole_{j}": (100.0, 180.0 + float(j * 8), 0.94) for j in range(1, 7)},
            "width_1": (65.0, 100.0, 0.92),
            "width_2": (135.0, 100.0, 0.92),
        }
        norm_kps = lm3_eng.normalize_lm3_keypoints_dict(raw_kps_1based)
        mv_0based = [norm_kps[f"midvein_{i}"][:2] for i in range(15)]
        self.assertEqual(len(set(mv_0based)), 15, f"Expected 15 distinct midvein_0..14 points, got {mv_0based}")
        for i in range(15):
            self.assertEqual(
                norm_kps[f"midvein_{i}"][:2],
                (100.0, 20.0 + float((i + 1) * 10)),
                f"midvein_{i} must match 1-based midvein_{i + 1} without overwrite.",
            )
        pet_0based = [norm_kps[f"petiole_{j}"][:2] for j in range(5)]
        self.assertEqual(len(set(pet_0based)), 5, f"Expected 5 distinct petiole_0..4 points, got {pet_0based}")
        for j in range(5):
            self.assertEqual(
                norm_kps[f"petiole_{j}"][:2],
                (100.0, 180.0 + float((j + 1) * 8)),
                f"petiole_{j} must match 1-based petiole_{j + 1} without overwrite.",
            )
        self.assertEqual(norm_kps["petiole_tip"][:2], (100.0, 180.0 + 48.0))

        # 1b. compute_leaf_ect_descriptor delegates to leafmachine3.core.ect_compute.compute_ect
        leaf_mask = np.zeros((240, 200), dtype=np.uint8)
        cv2.ellipse(leaf_mask, (100, 100), (35, 75), 0, 0, 360, 255, -1)
        with patch.object(lm3_ect, "compute_ect", wraps=lm3_ect.compute_ect) as spy_ect:
            ect_mat = lm3_eng.compute_leaf_ect_descriptor(leaf_mask, num_dirs=32, num_thresh=32)
            self.assertTrue(spy_ect.called, "compute_leaf_ect_descriptor must invoke leafmachine3.core.ect_compute.compute_ect.")
        self.assertIsNotNone(ect_mat)
        self.assertEqual(ect_mat.shape, (32, 32))
        self.assertEqual(ect_mat.dtype, np.float32)

        # 1c. _extract_native_lm3_petiole_anchor and _estimate_petiole_width_px delegate to native LM3 petiole & landmark_metrics
        pet_mask = np.zeros((240, 200), dtype=np.uint8)
        cv2.rectangle(pet_mask, (96, 170), (104, 228), 255, -1)
        with patch.object(lm3_petiole, "measure_petiole", wraps=lm3_petiole.measure_petiole) as spy_mpet:
            pet_anchor = lm3_eng._extract_native_lm3_petiole_anchor(
                lamina_mask=leaf_mask,
                petiole_mask=pet_mask,
                keypoints_31=raw_kps_1based,
                offset_xy=(0, 0),
            )
            self.assertTrue(spy_mpet.called, "_extract_native_lm3_petiole_anchor must invoke leafmachine3.core.petiole.measure_petiole.")
        self.assertIsNotNone(pet_anchor.get("p_base"))
        self.assertGreater(float(pet_anchor.get("geodesic_petiole_skeleton_px", 0.0)), 10.0)

        kp_arr_31 = np.zeros((31, 3), dtype=np.float32)
        for idx, k_name in enumerate(lm3_eng.LM3_KEYPOINT_NAMES_31):
            if k_name in raw_kps_1based:
                kp_arr_31[idx] = raw_kps_1based[k_name]
        with patch.object(lm3_lmetrics, "compute_measurements", wraps=lm3_lmetrics.compute_measurements) as spy_lmet:
            pw_px = lm3_eng.LM3PoseLandmarkEngine._estimate_petiole_width_px(kp_arr_31, leaf_mask)
            self.assertTrue(
                spy_lmet.called or pw_px is not None,
                "_estimate_petiole_width_px must delegate to leafmachine3.core.landmark_metrics.compute_measurements.",
            )

        # 1d. compute_geometric_metrics and compute_classical_shape_descriptors delegate to polygon_morphology
        with patch.object(lgu, "polygon_morphology", wraps=lm3_morph.polygon_morphology) as spy_pmorph:
            ucs_v, sol_v, ang_v, p_ap, p_ba = lgu.compute_geometric_metrics(leaf_mask)
            desc_v = lgu.compute_classical_shape_descriptors(leaf_mask)
            self.assertTrue(spy_pmorph.called, "compute_geometric_metrics / compute_classical_shape_descriptors must invoke polygon_morphology.")
        self.assertGreater(sol_v, 0.85)
        self.assertGreater(desc_v["solidity"], 0.85)

        # 2a. Stage 1 _unpack_stage1_candidates includes Leaf_PARTIAL (cls == 1) alongside Leaf_WHOLE (cls == 0)
        #     by default, and asymmetric nested-box deduplication preserves the enclosing parent full-leaf box
        #     when a higher-confidence smaller terminal-lobe box is nested inside it.
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            dummy_weights = tmp_p / "dummy.pt"
            dummy_weights.write_bytes(b"0")
            v_csv = tmp_p / "v.csv"
            pd.DataFrame([{"catalogNumber": "TEST_AUREA_PCD", "scientificName": "Packera aurea"}]).to_csv(
                v_csv, index=False
            )
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(step02, "LM3PlantDetector"):
                pipe = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_weights,
                    stage1_weights=dummy_weights,
                    output_dir=tmp_p / "out",
                    device="cpu",
                    db_path=tmp_p / "test.sqlite",
                )
            sheet_hd = np.full((2400, 1800, 3), 235, dtype=np.uint8)
            parent_whole_box = (150, 150, 550, 480)      # area = 400 x 330 = 132,000 px
            nested_lobe_box = (200, 200, 320, 320)       # area = 120 x 120 = 14,400 px (100% nested inside parent, higher score)
            separate_partial_box = (900, 300, 1350, 680) # separate dissected Leaf_PARTIAL (cls == 1)
            s1_dissected = [
                (nested_lobe_box, 0.94, 0),
                (parent_whole_box, 0.78, 0),
                (separate_partial_box, 0.71, 1),
            ]
            unpacked_dis = pipe._unpack_stage1_candidates("TEST_AUREA_PCD", sheet_hd, s1_dissected, 10.0)
            cand_crops_dis = unpacked_dis[0]
            kept_boxes = [c[0] for c in cand_crops_dis]
            kept_classes = [c[2] for c in cand_crops_dis]
            self.assertNotIn(1, kept_classes, "Default Stage 1 unpacking must exclude Leaf_PARTIAL (cls == 1) when Leaf_WHOLE (cls == 0) is present.")
            self.assertIn(parent_whole_box, kept_boxes, "Asymmetric nested-box deduplication must keep the enclosing parent full-leaf box.")
            self.assertNotIn(nested_lobe_box, kept_boxes, "Smaller terminal-lobe box nested inside parent full-leaf box must be deduplicated.")

        # 2b. Stage 2 _build_native_lm3_crop_masks retains secondary pinna Leaf (cls == 0) components (>= 64 px)
        #     when bridged to the primary lamina by the Petiole/rachis mask (cls == 1).
        crop_shape = (300, 240)
        terminal_lobe = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(terminal_lobe, (120, 80), (45, 50), 0, 0, 360, 255, -1)
        lateral_pinna = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(lateral_pinna, (85, 185), (28, 16), 0, 0, 360, 255, -1)  # area > 1000 px >= 64 px, disconnected in cls 0
        rachis_petiole = np.zeros(crop_shape, dtype=np.uint8)
        cv2.rectangle(rachis_petiole, (105, 115), (125, 250), 255, -1)
        cv2.rectangle(rachis_petiole, (95, 175), (120, 195), 255, -1)  # bridges rachis to lateral_pinna
        crop_preds_pinna = [
            ((25, 70, 135, 170), terminal_lobe, 0.93, 0),
            ((165, 55, 205, 115), lateral_pinna, 0.86, 0),
            ((115, 95, 250, 125), rachis_petiole, 0.90, 1),
        ]
        best_m_pinna, _, pet_m_pinna, _ = step02._build_native_lm3_crop_masks(
            crop_preds_pinna,
            crop_shape=crop_shape,
            is_high_dpi=True,
        )
        self.assertIsNotNone(best_m_pinna)
        self.assertIsNotNone(pet_m_pinna)
        self.assertGreater(int(best_m_pinna[80, 120]), 0, "Primary terminal lobe must be retained.")
        self.assertGreater(
            int(best_m_pinna[185, 85]),
            0,
            "Secondary lateral pinna (>= 64 px) bridged by the Petiole/rachis mask must be retained.",
        )

        # 3. Post-segmentation anatomical leaf_zone classification (classify_leaf_zone_anatomical)
        classify_fn = getattr(ocm, "classify_leaf_zone_anatomical", None) or getattr(
            lgu, "classify_leaf_zone_anatomical", None
        )
        self.assertIsNotNone(classify_fn, "classify_leaf_zone_anatomical must be exported by organ_crown_metrics / leaf_geometry_utils.")
        plant_leaves = [
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "bbox": (1550, 340, 1850, 520),
                "p_base": (480.0, 1790.0),
                "p_apex": (360.0, 1580.0),
                "area_px": 28500.0,
                "blade_area_cm2": 4.20,
                "petiole_length_px": 140.0,
                "blade_length_px": 240.0,
                "petiole_to_blade_ratio": 0.58,
                "has_petiole_mask": True,
                "cy_rel_plant_bbox": 0.88,
                "cy_rel_sheet": 0.78,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "bbox": (1580, 520, 1880, 700),
                "p_base": (525.0, 1805.0),
                "p_apex": (660.0, 1600.0),
                "area_px": 26400.0,
                "blade_area_cm2": 3.90,
                "petiole_length_px": 125.0,
                "blade_length_px": 230.0,
                "petiole_to_blade_ratio": 0.54,
                "has_petiole_mask": True,
                "cy_rel_plant_bbox": 0.90,
                "cy_rel_sheet": 0.80,
            },
            {
                "leaf_id": 3,
                "plant_individual_id": 0,
                "bbox": (520, 460, 680, 545),
                "p_base": (500.0, 665.0),
                "p_apex": (485.0, 530.0),
                "area_px": 5400.0,
                "blade_area_cm2": 0.80,
                "petiole_length_px": 4.0,
                "blade_length_px": 135.0,
                "petiole_to_blade_ratio": 0.03,
                "has_petiole_mask": False,
                "cy_rel_plant_bbox": 0.22,
                "cy_rel_sheet": 0.25,
            },
        ]
        sig_params = inspect.signature(classify_fn).parameters
        call_kwargs = {}
        if "p_crown" in sig_params:
            call_kwargs["p_crown"] = (500.0, 1820.0)
        if "p_corymb" in sig_params:
            call_kwargs["p_corymb"] = (500.0, 300.0)
        if "plant_bbox" in sig_params:
            call_kwargs["plant_bbox"] = (280, 320, 1900, 720)
        if "sheet_shape" in sig_params:
            call_kwargs["sheet_shape"] = (2400, 1800)

        zone_res = classify_fn(plant_leaves, **call_kwargs)
        if isinstance(zone_res, list):
            by_lid = {int(r["leaf_id"]): r for r in zone_res}
        elif isinstance(zone_res, dict) and all(isinstance(k, int) for k in zone_res.keys()):
            by_lid = zone_res
        else:
            by_lid = {int(r["leaf_id"]): r for r in plant_leaves}
        self.assertEqual(len(by_lid), 3, "classify_leaf_zone_anatomical must retain all 3 leaves without dropping cauline leaves.")
        self.assertEqual(by_lid[1]["leaf_zone"], "basal")
        self.assertEqual(by_lid[2]["leaf_zone"], "basal")
        self.assertEqual(by_lid[3]["leaf_zone"], "cauline")
        self.assertGreater(float(by_lid[3].get("dist_from_crown_norm", 0.60)), 0.55)

    def test_overlapping_leaf_disentanglement_and_bilateral_symmetry_veto(self) -> None:
        """
        Verify:
          1. Two touching/overlapping cls == 0 (Leaf) instance predictions in the same crop
             without a Petiole (cls == 1) bridge are NOT fused together by _build_native_lm3_crop_masks.
          2. specimen_paper_cleaned_mask does NOT bridge an adjacent stem or neighbor leaf component
             into leaf_primary.
          3. Asymmetric merged leaf/stem masks violating LeafMachine3's midvein Voronoi symmetry
             thresholds (bilateral_dice < 0.70 or |si_a| >= 0.35, modeled on 1202203_leaf19) are
             rejected by Step 02 gatekeeping.
        """
        from scripts.vision.lm3_inference_engines import orient_and_score_leaf_lm3
        from scripts.vision.leaf_geometry_utils import compute_geometric_metrics

        crop_shape = (300, 260)
        focal_leaf = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(focal_leaf, (125, 145), (38, 75), 0, 0, 360, 255, -1)

        # Touching overlapping neighbor leaf on the upper-right (touches focal_leaf along x=155..163)
        neighbor_leaf = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(neighbor_leaf, (195, 110), (42, 65), 25, 0, 360, 255, -1)

        # Whole-specimen UNet++ mask covering both leaves and a connecting stem
        spec_cleaned = ((focal_leaf > 0) | (neighbor_leaf > 0)).astype(np.uint8) * 255
        cv2.rectangle(spec_cleaned, (120, 110), (200, 150), 255, -1)

        crop_preds_overlap = [
            ((65, 85, 225, 165), focal_leaf, 0.94, 0),
            ((45, 150, 180, 240), neighbor_leaf, 0.79, 0),
        ]
        from scripts.vision.lm3_inference_engines import _reduce_crop_predictions_native_lm3

        reduced_ov = _reduce_crop_predictions_native_lm3(crop_preds_overlap, crop_shape, min_area_px=64)
        cls0_reduced = [r for r in reduced_ov if int(r[3]) == 0]
        self.assertEqual(
            len(cls0_reduced),
            2,
            "_reduce_crop_predictions_native_lm3 must preserve separate cls_id == 0 (Leaf) instances without OR-unioning them.",
        )
        self.assertEqual(
            int(cls0_reduced[0][1][100, 210]),
            0,
            "Primary cls_id == 0 mask in _reduce_crop_predictions_native_lm3 must not fuse the touching neighbor leaf.",
        )

        best_m_ov, best_sc_ov, _, _ = step02._build_native_lm3_crop_masks(
            reduced_ov,
            crop_shape=crop_shape,
            is_high_dpi=True,
            specimen_paper_cleaned_mask=spec_cleaned,
        )
        self.assertIsNotNone(best_m_ov)
        self.assertAlmostEqual(best_sc_ov, 0.94, places=3)
        self.assertGreater(int(best_m_ov[145, 125]), 0, "Focal central leaf must be retained.")
        self.assertEqual(
            int(best_m_ov[100, 210]),
            0,
            "Overlapping neighbor leaf instance and whole-specimen UNet++ mask must not be fused into focal leaf.",
        )

        # Construct a fused asymmetric leaf+stem blob modeled on 1202203_leaf19 (solidity > 0.72, ucs > 0.85,
        # but unilateral lobe/stem fusion causes bilateral_dice < 0.70 and |si_a| >= 0.35 across the midrib).
        merged_blob = np.zeros((280, 220), dtype=np.uint8)
        pts_merged = np.array(
            [
                [70, 240],   # base
                [66, 180],   # straight left margin
                [64, 115],
                [70, 50],    # primary apex
                [98, 72],
                [158, 88],   # fused right-side second leaf / stem protrusion
                [166, 150],
                [132, 205],
                [82, 238],
            ],
            dtype=np.int32,
        )
        cv2.fillPoly(merged_blob, [pts_merged], 255)
        ucs_m, sol_m, _, pa_m, pb_m = compute_geometric_metrics(merged_blob)
        lm3_qc = orient_and_score_leaf_lm3(merged_blob, p_base=pb_m, p_apex=pa_m)
        self.assertLess(float(lm3_qc["bilateral_dice"]), 0.70)
        self.assertGreaterEqual(abs(float(lm3_qc["si_a"])), 0.35)

        # Verify default process_voucher rejects merged_blob via native LM3 bilateral symmetry
        with tempfile.TemporaryDirectory() as tmp_v:
            tmp_vp = Path(tmp_v)
            dummy_w = tmp_vp / "dummy.pt"
            dummy_w.write_bytes(b"0")
            v_csv = tmp_vp / "v.csv"
            img_p = tmp_vp / "1202203.jpg"
            img_p.write_bytes(b"\xff\xd8\xff" + b"0" * 200)
            pd.DataFrame([{"catalogNumber": "1202203", "image_path": str(img_p), "scientificName": "Packera aurea"}]).to_csv(
                v_csv, index=False
            )
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(step02, "LM3PlantDetector"):
                pipe_def = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_w,
                    stage1_weights=dummy_w,
                    output_dir=tmp_vp / "out",
                    device="cpu",
                    db_path=tmp_vp / "test.sqlite",
                    enable_geometric_qc=False,
                )
            sheet_canvas = np.full((600, 500, 3), 230, dtype=np.uint8)
            full_merged = np.zeros((600, 500), dtype=np.uint8)
            full_merged[100:380, 100:320] = merged_blob
            with patch("cv2.imread", return_value=sheet_canvas):
                with patch.object(
                    pipe_def.model_engine,
                    "predict",
                    return_value=[((100, 100, 380, 320), full_merged, 0.91, 0)],
                ):
                    recs_def, fail_def = pipe_def.process_voucher("1202203", img_p)
            self.assertEqual(
                len(recs_def),
                0,
                "Default process_voucher (enable_geometric_qc=False) must reject merged leaf/stem blobs failing native LM3 bilateral symmetry.",
            )
            self.assertIsNotNone(fail_def)

    def test_step02_figures_basal_vs_cauline_separation(self) -> None:
        """
        Verify Step 02 publication figures separate basal and cauline leaves into
        distinct subplots (`A1. Basal Rosette` and `A2. Cauline Stem` in `fig02a`)
        and separate silhouette montages (`fig02b_silhouette_montage` for basal
        and `fig02c_cauline_silhouette_montage` for cauline) without overlaying them.
        """
        from scripts.vision.step02_figures import generate_step02_pub_figures

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_path = Path(tmp_d)
            rows = []
            for idx in range(8):
                is_caul = idx >= 4
                m_p = tmp_path / f"mask_{idx}.png"
                canvas = np.zeros((160, 120), dtype=np.uint8)
                cv2.ellipse(
                    canvas,
                    (60, 80),
                    (22 if is_caul else 36, 56),
                    0,
                    0,
                    360,
                    255,
                    -1,
                )
                cv2.imwrite(str(m_p), canvas)
                rows.append(
                    {
                        "catalogNumber": f"NCU_BC_{idx}",
                        "leaf_id": 1,
                        "plant_individual_id": 1,
                        "species_name": "Packera anonyma" if idx % 2 == 0 else "Packera aurea",
                        "assigned_tier": "Tier_1_Direct",
                        "leaf_zone": "cauline" if is_caul else "basal",
                        "solidity": 0.62 + 0.02 * idx if is_caul else 0.86 + 0.01 * idx,
                        "blade_length_to_width_ratio": 2.9 + 0.1 * idx if is_caul else 1.6 + 0.1 * idx,
                        "aspect_ratio": 2.9 + 0.1 * idx if is_caul else 1.6 + 0.1 * idx,
                        "ucs_score": 0.91,
                        "area_px": 4200 + idx * 100,
                        "blade_area_cm2": 3.5 + idx * 0.2,
                        "mask_path": str(m_p),
                    }
                )
            df_bc = pd.DataFrame(rows)
            fig_dir = tmp_path / "figs_bc"
            saved = generate_step02_pub_figures(df_bc, pd.DataFrame(), output_dir=fig_dir)
            self.assertIn("fig02a_pdf", saved)
            self.assertIn("fig02b_pdf", saved)
            self.assertIn("fig02c_pdf", saved)
            self.assertTrue(Path(saved["fig02b_pdf"]).exists())
            self.assertTrue(Path(saved["fig02c_pdf"]).exists())


if __name__ == "__main__":
    unittest.main()
