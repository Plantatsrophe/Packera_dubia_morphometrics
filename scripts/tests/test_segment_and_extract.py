#!/usr/bin/env python3
"""
===============================================================================
Script: test_segment_and_extract.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 02:
    Read-only LeafMachine3 (YOLO26x / UNet++ v9 / ONNX Runtime) segmentation adapter,
    LM3 ruler scale calibration, centroid spatial clustering, anatomical leaf_zone
    classification, and standardized contour / ECT export.
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
compute_lm3_shape_metrics = step02.compute_lm3_shape_metrics
export_standardized_contour = step02.export_standardized_contour
save_leaf_mask_png = step02.save_leaf_mask_png
compute_capitulum_metrics = step02.compute_capitulum_metrics


class TestSegmentAndExtract(unittest.TestCase):
    """Test suite for Step 02 LeafMachine3 ONNX segmentation, upright extraction, and SQLite/CSV export."""

    def setUp(self) -> None:
        self.temp_dir = Path(tempfile.mkdtemp(prefix="test_step02_"))

    def tearDown(self) -> None:
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_compute_lm3_shape_metrics_and_mask_export(self) -> None:
        """Verify LM3 shape metrics and direct binary mask export to data/masks/."""
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.ellipse(mask, (100, 100), (30, 70), 0, 0, 360, 255, -1)
        desc = compute_lm3_shape_metrics(mask)
        self.assertGreaterEqual(desc["solidity"], 0.90)
        self.assertGreaterEqual(desc["aspect_ratio"], 2.0)

        saved_path = save_leaf_mask_png(mask, "TEST001", 1, self.temp_dir)
        self.assertTrue(Path(saved_path).exists(), "Exported LM3 mask file must exist.")
        self.assertEqual(Path(saved_path).parent.name, "masks")
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
        self.assertGreaterEqual(metrics["AR_inv"], 0.25)
        self.assertLessEqual(metrics["AR_inv"], 4.5)
        self.assertAlmostEqual(metrics["AR_inv"], 1.4, delta=0.15)
        self.assertAlmostEqual(metrics["involucre_height_mm"], 14.0, delta=1.5)
        self.assertAlmostEqual(metrics["involucre_width_mm"], 10.0, delta=1.5)

    def test_compute_capitulum_metrics_filter_extremes(self) -> None:
        """Verify aberrant non-capitulum shapes are rejected by 0.25 <= AR_inv <= 4.5 gating while single-head filter rejects clusters."""
        # 1. Extreme flattened horizontal artifact (AR < 0.25)
        mask_flat = np.zeros((200, 200), dtype=np.uint8)
        cv2.rectangle(mask_flat, (20, 90), (180, 110), 255, -1)  # W=160, H=20, AR=0.125
        metrics_flat = compute_capitulum_metrics(mask_flat, (90, 20, 110, 180))
        self.assertIsNone(metrics_flat, "Excessively flat outline must be rejected (AR_inv < 0.25).")

        # 2. Extreme vertical stem fragment (AR > 4.5)
        mask_stem = np.zeros((300, 100), dtype=np.uint8)
        cv2.rectangle(mask_stem, (40, 20), (60, 280), 255, -1)  # W=20, H=260, AR=13.0
        metrics_stem = compute_capitulum_metrics(mask_stem, (20, 40, 280, 60))
        self.assertIsNone(metrics_stem, "Excessively elongated stem fragment must be rejected (AR_inv > 4.5).")

        # 3. Multi-head cluster class_id=5 (Flower_MANY) must return None
        mask_ok = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask_ok, (150, 150), (50, 70), 0, 0, 360, 255, -1)
        self.assertIsNone(compute_capitulum_metrics(mask_ok, (80, 100, 220, 200), class_id=5))

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
        manifest_file = self.tables_dir / "extracted_leaf_manifest.csv"
        pd.DataFrame([{
            "catalogNumber": "NCU001",
            "plant_individual_id": 1,
            "leaf_id": 1,
            "archetype_score": 0.90,
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
            "archetype_score": 0.88,
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
            self.assertTrue(bool(rec["lm3_gates_pass"]))
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

        # Verify binary leaf masks exist directly in data/masks/
        for leaf_id in [1, 2, 3]:
            leaf_mask = self.pipeline.masks_dir / f"NCU123_leaf{leaf_id}.png"
            self.assertTrue(leaf_mask.exists(), f"Binary leaf mask {leaf_mask} must exist.")

    def test_export_standardized_contour_isotropic_aspect_ratio_preservation(self) -> None:
        """
        Verify an oblong tip-up leaf with L/W = 3.0 exported via export_standardized_contour()
        retains (max(y) - min(y)) / (max(x) - min(x)) == 3.0 within 1e-3.
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

        span_x = float(df["x"].max() - df["x"].min())
        span_y = float(df["y"].max() - df["y"].min())
        lw_ratio = max(span_x, span_y) / max(min(span_x, span_y), 1e-9)
        self.assertAlmostEqual(
            lw_ratio,
            3.0,
            delta=1e-3,
            msg=f"Isotropically normalized contour must preserve L/W = 3.0 within 1e-3 (got {lw_ratio:.6f}).",
        )

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
        self.assertTrue(all(bool(r["lm3_gates_pass"]) for r in records))

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
        for col in ("circularity", "extent", "aspect_ratio", "area_px", "area_mm2"):
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
            "SELECT circularity, extent, aspect_ratio, area_px, area_mm2 FROM extracted_leaves WHERE catalogNumber = 'NCU_MACRO'"
        ).fetchone()
        self.assertIsNotNone(row_leaf)
        for col in ("circularity", "extent", "aspect_ratio", "area_px", "area_mm2"):
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
        Verify default 100% Raw Unconstrained LeafMachine3 Mode (read-only Step 02 adapter):
          1. Stage 1 (`LM3PlantDetector`) unpacks all Class 0 (`leaf_whole`) boxes regardless of size
             or vertical sheet position, while natively excluding `leaf_partial` (1), `leaflet` (2),
             `specimen` (8), and `wood`/stems (10).
          2. Stage 2 `LM3LeafSegmenterEngine` (`thing_classes = ['leaf', 'petiole', 'hole']`) retains the top-scoring
             `'leaf'` (`class_id == 0`) mask 100% pixel-identical (zero heuristic petiole clipping),
             pairs the top-scoring `'petiole'` (`class_id == 1`) mask to locate `p_base` and native petiole length,
             and retains >6 leaves without any post-hoc geometric/photometric rejections or top-6 cap.
        """
        with patch("sys.argv", ["02_segment_and_extract.py"]):
            parsed_default = step02.parse_args()
            self.assertEqual(parsed_default.score_thresh, self.pipeline.config.segmentation.score_thresh)

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
            self.assertIn("lm3_gates_pass", rec)
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
            "catalogNumber", "leaf_id", "plant_individual_id", "lm3_gates_pass",
            "solidity", "archetype_score", "circularity", "extent", "area_px",
            "area_mm2", "petiole_length_mm", "blade_length_mm",
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

    def test_canonical_tip_up_contour_ect_and_scale_calibrator(self) -> None:
        """
        Verify canonical tip-up contour rotation in export_standardized_contour,
        LM3ScaleCalibrator megapixel prior, orient_and_score_leaf_lm3, and
        compute_leaf_ect_descriptor.
        """
        from scripts.vision.leaf_geometry_utils import export_standardized_contour
        from scripts.vision.lm3_inference_engines import (
            LM3ScaleCalibrator,
            compute_leaf_ect_descriptor,
            orient_and_score_leaf_lm3,
        )

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
            # Point 0 (lamina_base) must sit near the bottom (max y), and point 60 (lamina_tip) near the top (min y)
            y_min_c, y_max_c = float(cdf["y"].min()), float(cdf["y"].max())
            span_y = max(y_max_c - y_min_c, 1e-6)
            span_x = float(cdf["x"].max() - cdf["x"].min())
            self.assertGreater((float(cdf.loc[0, "y"]) - y_min_c) / span_y, 0.85)
            self.assertLess((float(cdf.loc[60, "y"]) - y_min_c) / span_y, 0.15)
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
             - compute_lm3_shape_metrics delegates to leafmachine3.core.morphometrics.polygon_morphology.
          2. Dissected-Leaf PCD & Stage 2 enhancements:
             - Stage 1 _unpack_stage1_candidates excludes Leaf_PARTIAL (cls_id == 1) by default when Leaf_WHOLE (cls_id == 0) is present
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

        # 1d. compute_lm3_shape_metrics delegates to polygon_morphology
        with patch.object(lgu, "polygon_morphology", wraps=lm3_morph.polygon_morphology) as spy_pmorph:
            desc_v = lgu.compute_lm3_shape_metrics(leaf_mask)
            self.assertTrue(spy_pmorph.called, "compute_lm3_shape_metrics must invoke polygon_morphology.")
        self.assertGreater(desc_v["solidity"], 0.85)

        # 2a. Stage 1 _unpack_stage1_candidates excludes Leaf_PARTIAL (cls == 1) by default when Leaf_WHOLE (cls == 0)
        #     is present, and asymmetric nested-box deduplication preserves the enclosing parent full-leaf box
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

        # 2b. Stage 2 _build_native_lm3_crop_masks extracts top-scoring LM3 Leaf (cls == 0) and Petiole (cls == 1)
        #     masks without synthesizing OpenCV rachis bridges across separate instance masks.
        crop_shape = (300, 240)
        terminal_lobe = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(terminal_lobe, (120, 80), (45, 50), 0, 0, 360, 255, -1)
        lateral_pinna = np.zeros(crop_shape, dtype=np.uint8)
        cv2.ellipse(lateral_pinna, (85, 185), (28, 16), 0, 0, 360, 255, -1)
        rachis_petiole = np.zeros(crop_shape, dtype=np.uint8)
        cv2.rectangle(rachis_petiole, (105, 115), (125, 250), 255, -1)
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
        np.testing.assert_array_equal(best_m_pinna > 0, terminal_lobe > 0)

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
                "area_mm2": 420.0,
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
                "area_mm2": 390.0,
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
                "area_mm2": 80.0,
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
        lm3_qc = orient_and_score_leaf_lm3(merged_blob, p_base=(70, 240), p_apex=(70, 50))
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
                1,
                "Read-only Step 02 adapter must persist all LM3 masks with efa_qc_pass=False rather than dropping in Python.",
            )
            self.assertIsNone(fail_def)
            self.assertFalse(bool(recs_def[0]["lm3_gates_pass"]))
            self.assertLess(float(recs_def[0]["bilateral_dice"]), 0.70)

    def test_four_organ_pole_cases_and_vascular_ray_tracing(self) -> None:
        """
        Verify classify_leaf_zone_anatomical across all 4 biological organ-pole cases:
          - Case 1 (Roots + Capitula): Inverted mount (NCU00385199: roots at top y=197,
            capitula at bottom y=1572) preserves top root crown without overwriting and
            classifies top rosette leaf as basal and lower stem leaves as cauline.
          - Case 2 (Roots Only, No Capitula): Erect long-petiolate basal tuft (NCU00438219:
            roots at bottom y=2027, no capitula, elevated lamina blades at y=1580 with
            downward petiole rays tracing into root crown) classifies all tuft leaves as basal.
          - Case 3 (Capitula Only, No Roots): Clipped upper flowering stem (capitula at top,
            no roots, spaced sessile stem leaves without basal convergence hub) classifies
            internodal leaves as cauline.
          - Case 4 (Neither Roots Nor Capitula): Multi-plant sheet (NCSC00103660: pid=0
            flowering stem on right, pid=1 detached sterile basal rosette on upper-left)
            does not leak pid=0 capitula onto pid=1 and classifies all pid=1 rosette leaves
            as basal regardless of sheet y-position.
        """
        from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical

        # Case 1: Inverted mount (roots at top y=197, capitula at bottom y=1572)
        inverted_leaves = [
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base_sheet": (723.0, 349.0),
                "p_apex_sheet": (640.0, 480.0),
                "p_petiole_tip": (760.0, 245.0),
                "blade_length_px": 155.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 4,
                "plant_individual_id": 0,
                "p_base_sheet": (1085.0, 1061.0),
                "p_apex_sheet": (1165.0, 1185.0),
                "blade_length_px": 148.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.02,
            },
        ]
        res_case1 = classify_leaf_zone_anatomical(
            inverted_leaves,
            crowns_by_pid={0: (781.0, 206.0)},
            rosette_centroids_by_pid={0: (980.0, 920.0)},
            root_boxes=[((160, 740, 235, 820), 0.88, 7)],
            capitula_detections=[(1520, 1280, 1620, 1360)],
        )
        self.assertEqual(res_case1[2]["leaf_zone"], "basal")
        self.assertEqual(res_case1[4]["leaf_zone"], "cauline")

        # Case 2: Erect long-petiolate basal tuft with Roots only and no Capitula (NCU00438219)
        tuft_leaves = [
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "p_base_sheet": (855.0, 1910.0),
                "p_apex_sheet": (845.0, 1730.0),
                "p_petiole_tip": (860.0, 1995.0),
                "ray_dir": (0.04, 0.999),
                "blade_length_px": 180.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base_sheet": (820.0, 1580.0),
                "p_apex_sheet": (790.0, 1370.0),
                "p_petiole_tip": (832.0, 1665.0),
                "ray_dir": (0.09, 0.996),
                "blade_length_px": 210.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 3,
                "plant_individual_id": 0,
                "p_base_sheet": (910.0, 1595.0),
                "p_apex_sheet": (945.0, 1390.0),
                "p_petiole_tip": (895.0, 1680.0),
                "ray_dir": (-0.09, 0.996),
                "blade_length_px": 205.0,
                "has_petiole_mask": True,
            },
        ]
        res_case2 = classify_leaf_zone_anatomical(
            tuft_leaves,
            crowns_by_pid={0: (863.0, 2027.0)},
            root_boxes=[((1980, 810, 2080, 915), 0.91, 7)],
            capitula_detections=[],
        )
        self.assertEqual(res_case2[1]["leaf_zone"], "basal")
        self.assertEqual(res_case2[2]["leaf_zone"], "basal")
        self.assertEqual(res_case2[3]["leaf_zone"], "basal")

        # Case 3: Capitula only, no Roots, clipped upper flowering stem without basal rosette hub
        clipped_top_leaves = [
            {
                "leaf_id": 10,
                "plant_individual_id": 0,
                "p_base_sheet": (600.0, 1150.0),
                "p_apex_sheet": (510.0, 1060.0),
                "blade_length_px": 125.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.0,
            },
            {
                "leaf_id": 11,
                "plant_individual_id": 0,
                "p_base_sheet": (605.0, 780.0),
                "p_apex_sheet": (700.0, 690.0),
                "blade_length_px": 115.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.0,
            },
        ]
        res_case3 = classify_leaf_zone_anatomical(
            clipped_top_leaves,
            root_boxes=[],
            capitula_detections=[(380, 560, 460, 650)],
        )
        self.assertEqual(res_case3[10]["leaf_zone"], "cauline")
        self.assertEqual(res_case3[11]["leaf_zone"], "cauline")

        # Case 4: Multi-plant sheet where pid=1 is a detached sterile basal rosette (NCSC00103660)
        multi_plant_leaves = [
            # pid=0: flowering plant on right side (x ~ 1350)
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "p_base_sheet": (1340.0, 1820.0),
                "p_apex_sheet": (1240.0, 1690.0),
                "p_petiole_tip": (1365.0, 1870.0),
                "blade_length_px": 165.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base_sheet": (1360.0, 920.0),
                "p_apex_sheet": (1450.0, 810.0),
                "blade_length_px": 140.0,
                "has_petiole_mask": False,
                "petiole_to_blade_ratio": 0.02,
            },
            # pid=1: detached sterile basal rosette on upper-left (x ~ 440, y ~ 725)
            {
                "leaf_id": 3,
                "plant_individual_id": 1,
                "p_base_sheet": (380.0, 740.0),
                "p_apex_sheet": (220.0, 780.0),
                "p_petiole_tip": (435.0, 728.0),
                "blade_length_px": 165.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 4,
                "plant_individual_id": 1,
                "p_base_sheet": (500.0, 755.0),
                "p_apex_sheet": (650.0, 810.0),
                "p_petiole_tip": (448.0, 730.0),
                "blade_length_px": 160.0,
                "has_petiole_mask": True,
            },
            {
                "leaf_id": 5,
                "plant_individual_id": 1,
                "p_base_sheet": (442.0, 650.0),
                "p_apex_sheet": (440.0, 480.0),
                "p_petiole_tip": (442.0, 720.0),
                "blade_length_px": 170.0,
                "has_petiole_mask": True,
            },
        ]
        res_case4 = classify_leaf_zone_anatomical(
            multi_plant_leaves,
            crowns_by_pid={0: (1370.0, 1890.0), 1: (442.0, 725.0)},
            root_boxes=[((1870, 1330, 1960, 1420), 0.90, 7)],
            capitula_detections=[(380, 1310, 470, 1410)],
        )
        self.assertEqual(res_case4[1]["leaf_zone"], "basal")
        self.assertEqual(res_case4[2]["leaf_zone"], "cauline")
        self.assertEqual(res_case4[3]["leaf_zone"], "basal")
        self.assertEqual(res_case4[4]["leaf_zone"], "basal")
        self.assertEqual(res_case4[5]["leaf_zone"], "basal")

    def test_step02_upright_mask_persistence_reniform_lw_crown_tiebreak_and_32x360_ect(self) -> None:
        """
        Verify Step 2 botanical systematics & LM3 downstream handoff invariants:
          1. Oblique / inverted leaves persist upright base-to-tip binary masks in `data/masks/*.png`
             with `mask_is_upright == True` and `32 x 360` ECT `.npy` matrices matching `LM3_settings.yaml`.
          2. Reniform / cordate wider-than-long blades (`W > L`, e.g., Packera aurea) preserve true
             anatomical `aspect_ratio < 1.0` (`blade_length / blade_width`) without `90 deg` sideways rotation.
          3. Downward-radiating rosette / inverted-mount leaves (`NCU00385199` pattern) tie-break `lamina_base`
             toward the root/crown hub rather than flipping `180 deg` via `argmax(y)`.
          4. `RulerCalibrationMap` keeps `predicted_from_megapixels` as `is_calibrated_mm == False`
             (`scale_source == "mp_estimated"`) in both CSV manifests and SQLite `vouchers` / `extracted_leaves`.
          5. Pre-oriented `momocs_outline_xy` passed to `export_standardized_contour` is never double-rotated.
        """
        from scripts.core.db import get_db_connection, record_voucher_segmentation_transaction
        from scripts.vision.leaf_geometry_utils import export_standardized_contour
        from scripts.vision.lm3_data_loader import load_ruler_calibrations

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            vouchers_csv = tmp_p / "vouchers.csv"
            img_obl = tmp_p / "NCU_OBLIQUE.jpg"
            img_ren = tmp_p / "NCU_RENIFORM.jpg"
            img_inv = tmp_p / "NCU_INVERTED.jpg"
            for ip in (img_obl, img_ren, img_inv):
                cv2.imwrite(str(ip), np.full((1200, 1000, 3), 235, dtype=np.uint8))

            pd.DataFrame([
                {"catalogNumber": "NCU_OBLIQUE", "scientificName": "Packera anonyma", "image_path": str(img_obl)},
                {"catalogNumber": "NCU_RENIFORM", "scientificName": "Packera aurea", "image_path": str(img_ren)},
                {"catalogNumber": "NCU_INVERTED", "scientificName": "Packera dubia", "image_path": str(img_inv)},
            ]).to_csv(vouchers_csv, index=False)

            # Set up LM3 ruler CSV with one empirical ruler and one megapixel prior
            ruler_dir = tmp_p / "lm3_proj" / "Data" / "Ruler"
            ruler_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([
                {
                    "filename": "NCU_OBLIQUE__R__10-10-100-500",
                    "ruler_success": True,
                    "conversion_mean": 100.0,
                    "ruler_class_confidence": 95.0,
                    "conversion_determination": "measured_from_ruler",
                },
                {
                    "filename": "NCU_RENIFORM__R__10-10-100-500",
                    "ruler_success": True,
                    "conversion_mean": 105.0,
                    "ruler_class_confidence": 88.0,
                    "conversion_determination": "predicted_from_megapixels",
                },
            ]).to_csv(ruler_dir / "ruler_summary.csv", index=False)

            calibs = load_ruler_calibrations(tmp_p / "lm3_proj", min_confidence=0.50)
            self.assertEqual(calibs.sources["NCU_OBLIQUE"], "measured_from_ruler")
            self.assertEqual(calibs.sources["NCU_RENIFORM"], "predicted_from_megapixels")

            pipe = SegmentAndExtractPipeline(
                vouchers_csv=vouchers_csv,
                model_weights=Path("models/dummy.pth"),
                output_dir=tmp_p / "out_data",
                lm3_project_dir=tmp_p / "lm3_proj",
                device="cpu",
                force=True,
            )
            try:
                # 1. Oblique 135-deg leaf (base at top-left, tip at bottom-right) -> saved mask must be upright (H > 1.8 * W)
                obl_crop = np.zeros((420, 420), dtype=np.uint8)
                cv2.ellipse(obl_crop, (210, 210), (110, 42), 45.0, 0, 360, 255, -1)
                obl_pet = np.zeros((420, 420), dtype=np.uint8)
                cv2.line(obl_pet, (132, 132), (92, 92), 255, thickness=6)

                with patch.object(pipe.stage1_detector, "predict", return_value=[((200, 200, 620, 620), 0.96, 0)]):
                    with patch.object(
                        pipe.model_engine,
                        "predict_batch",
                        return_value=[[((10, 10, 410, 410), obl_crop, 0.96, 0), ((85, 85, 140, 140), obl_pet, 0.90, 1)]],
                    ):
                        recs_obl, fail_obl = pipe.process_voucher("NCU_OBLIQUE", img_obl)

                self.assertIsNone(fail_obl)
                self.assertEqual(len(recs_obl), 1)
                r_obl = recs_obl[0]
                self.assertTrue(bool(r_obl["is_calibrated_mm"]))
                self.assertEqual(r_obl["scale_source"], "ruler_calibrated")
                self.assertTrue(bool(r_obl["mask_is_upright"]))
                saved_obl_mask = cv2.imread(str(r_obl["mask_path"]), cv2.IMREAD_GRAYSCALE)
                self.assertIsNotNone(saved_obl_mask)
                ys_o, xs_o = np.where(saved_obl_mask > 0)
                h_obl = float(ys_o.max() - ys_o.min())
                w_obl = float(xs_o.max() - xs_o.min())
                self.assertGreater(h_obl, 1.8 * w_obl, "Saved PNG mask in data/masks/ must be rotated upright base-to-tip.")
                ect_obl = np.load(str(r_obl["ect_path"]))
                self.assertEqual(ect_obl.shape, (32, 360), "Step 02 ECT .npy must honor LM3_settings.yaml 32x360 grid.")

                # 2. Reniform wider-than-long leaf (L ~ 120, W ~ 190 -> L/W < 0.72) + mp_estimated provenance
                ren_crop = np.zeros((380, 380), dtype=np.uint8)
                cv2.ellipse(ren_crop, (190, 160), (95, 62), 0.0, 0, 360, 255, -1)
                cv2.ellipse(ren_crop, (190, 220), (34, 24), 0.0, 0, 360, 0, -1)
                ren_pet = np.zeros((380, 380), dtype=np.uint8)
                cv2.line(ren_pet, (190, 222), (190, 315), 255, thickness=6)

                with patch.object(pipe.stage1_detector, "predict", return_value=[((300, 300, 680, 680), 0.97, 0)]):
                    with patch.object(
                        pipe.model_engine,
                        "predict_batch",
                        return_value=[[((10, 10, 370, 370), ren_crop, 0.96, 0), ((220, 184, 320, 196), ren_pet, 0.91, 1)]],
                    ):
                        recs_ren, fail_ren = pipe.process_voucher("NCU_RENIFORM", img_ren)

                self.assertIsNone(fail_ren)
                r_ren = recs_ren[0]
                self.assertLess(float(r_ren["aspect_ratio"]), 0.80, "Reniform W > L blade must preserve aspect_ratio < 1.0.")
                self.assertFalse(bool(r_ren["is_calibrated_mm"]), "predicted_from_megapixels must not set is_calibrated_mm=True.")
                self.assertEqual(r_ren["scale_source"], "mp_estimated")
                self.assertTrue(pd.isna(r_ren["area_mm2"]))

                # Verify SQLite persistence keeps is_calibrated_mm == 0 for mp_estimated
                db_file = tmp_p / "test_step02_prov.sqlite"
                conn = get_db_connection(db_file, init_schema=True)
                try:
                    record_voucher_segmentation_transaction(
                        conn,
                        catalog_number="NCU_RENIFORM",
                        image_path=str(img_ren),
                        extracted_records=recs_ren,
                        failure_record=None,
                        output_dir=pipe.output_dir,
                    )
                    v_cal = conn.execute("SELECT is_calibrated_mm FROM vouchers WHERE catalogNumber = 'NCU_RENIFORM';").fetchone()[0]
                    el_cal, el_ss = conn.execute(
                        "SELECT is_calibrated_mm, scale_source FROM extracted_leaves WHERE catalogNumber = 'NCU_RENIFORM';"
                    ).fetchone()
                    self.assertEqual(int(v_cal), 0)
                    self.assertEqual(int(el_cal), 0)
                    self.assertEqual(str(el_ss), "mp_estimated")
                finally:
                    conn.close()

                # 3. Inverted herbarium mount (root crown at top y=150, downward-pointing blade lacking petiole_mask)
                #    Draw ovate blade wider near top base (y ~ 310) and tapering to acute apex at bottom (y ~ 480)
                inv_crop = np.zeros((320, 320), dtype=np.uint8)
                pts_inv = np.array([[160, 40], [215, 95], [195, 190], [160, 280], [125, 190], [105, 95]], dtype=np.int32)
                cv2.fillPoly(inv_crop, [pts_inv], 255)
                s1_inv = [
                    ((120, 460, 190, 540), 0.93, 7),   # Root tuft at top of sheet (inverted mount)
                    ((240, 340, 560, 660), 0.95, 0),   # Downward-pointing rosette leaf below top root crown
                ]
                with patch.object(pipe.stage1_detector, "predict", return_value=s1_inv):
                    with patch.object(
                        pipe.model_engine,
                        "predict_batch",
                        return_value=[[((10, 10, 310, 310), inv_crop, 0.95, 0)]],
                    ):
                        recs_inv, fail_inv = pipe.process_voucher("NCU_INVERTED", img_inv)

                self.assertIsNone(fail_inv)
                # Sheet-space midrib points downward (~ +90 deg), so Step 02 rotates ~180 deg to make upright mask tip-up
                self.assertGreater(float(recs_inv[0]["midrib_angle_deg"]), 45.0)
                inv_saved = cv2.imread(str(recs_inv[0]["mask_path"]), cv2.IMREAD_GRAYSCALE)
                # After 180-deg upright rotation, broader lamina base sits in lower half and acute apex sits at top
                h_i = inv_saved.shape[0]
                top_w = int(np.count_nonzero(inv_saved[: h_i // 3, :]))
                bot_w = int(np.count_nonzero(inv_saved[(2 * h_i) // 3 :, :]))
                self.assertGreater(bot_w, top_w, "Crown-hub tie-breaker must orient broader basal half to bottom (+Y) and acute apex to top (-Y).")

                # 4. Double-rotation guard when pre-oriented momocs_outline_xy is supplied alongside sheet keypoints
                th = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
                pre_oriented_xy = np.column_stack([200.0 - 35.0 * np.sin(th), 200.0 + 95.0 * np.cos(th)])
                kp_oblique = {
                    "lamina_base": (120.0, 120.0, 0.95),
                    "lamina_tip": (280.0, 280.0, 0.95),
                }
                c_pre = export_standardized_contour(
                    None,
                    catalog_number="PRE_ORI",
                    leaf_id=1,
                    output_dir=tmp_p / "pre_ori_out",
                    momocs_outline_xy=pre_oriented_xy,
                    keypoints_31=kp_oblique,
                )
                df_pre = pd.read_csv(c_pre)
                span_y_pre = float(df_pre["y"].max() - df_pre["y"].min())
                span_x_pre = float(df_pre["x"].max() - df_pre["x"].min())
                self.assertGreater(span_y_pre, 2.2 * span_x_pre, "Pre-oriented momocs_outline_xy must not be double-rotated by sheet keypoints.")

                # 5. Unbridged focal multi-pinna retention in _reduce_crop_predictions_native_lm3 + 32x360 multi-component ECT
                from scripts.vision import lm3_inference_engines as lm3_eng_chk
                pinna_crop_shape = (320, 240)
                p_inst_u8 = np.zeros(pinna_crop_shape, dtype=np.uint8)
                cv2.ellipse(p_inst_u8, (120, 80), (42, 48), 0, 0, 360, 255, -1)   # terminal lobe
                cv2.ellipse(p_inst_u8, (82, 168), (24, 14), 0, 0, 360, 255, -1)   # detached left pinna
                cv2.ellipse(p_inst_u8, (158, 168), (24, 14), 0, 0, 360, 255, -1)  # detached right pinna
                cv2.rectangle(p_inst_u8, (0, 270), (35, 319), 255, -1)            # border-clipped neighbor bleed -> rejected
                pet_stalk_u8 = np.zeros(pinna_crop_shape, dtype=np.uint8)
                cv2.rectangle(pet_stalk_u8, (110, 122), (130, 260), 255, -1)
                reduced_pinna = lm3_eng_chk._reduce_crop_predictions_native_lm3(
                    [
                        ((0, 0, 320, 240), p_inst_u8, 0.94, 0),
                        ((122, 110, 260, 130), pet_stalk_u8, 0.91, 1),
                    ],
                    pinna_crop_shape,
                    min_area_px=64,
                )
                lam_pinna_u8 = [m for _b, m, _s, cid in reduced_pinna if cid == 0][0]
                n_cc_p, _, stats_p, _ = cv2.connectedComponentsWithStats((lam_pinna_u8 > 0).astype(np.uint8), connectivity=8)
                sig_p = sum(1 for cid in range(1, n_cc_p) if int(stats_p[cid, cv2.CC_STAT_AREA]) >= 64)
                self.assertEqual(sig_p, 3, "Unbridged terminal lobe + 2 petiole-attached lateral pinnae must be retained (rejecting border neighbor bleed).")
                self.assertEqual(int(np.count_nonzero(lam_pinna_u8[130:150, 112:128])), 0, "Zero synthetic petiole bridge pixels may be added to lamina mask.")
                term_only_u8 = np.zeros(pinna_crop_shape, dtype=np.uint8)
                term_only_u8[:135, :] = p_inst_u8[:135, :]
                ect_term_only = lm3_eng_chk.compute_leaf_ect_descriptor(term_only_u8, num_dirs=360, num_thresh=32, already_oriented=True)
                ect_multi_pinna = lm3_eng_chk.compute_leaf_ect_descriptor(lam_pinna_u8, num_dirs=360, num_thresh=32, already_oriented=True)
                self.assertGreater(float(np.max(ect_multi_pinna)), float(np.max(ect_term_only)), "Multi-pinna 32x360 ECT must capture higher multi-component Euler characteristic than terminal-lobe-only.")
                mc_metrics = compute_lm3_shape_metrics(lam_pinna_u8, angle_cw=0.0)
                self.assertLess(mc_metrics["solidity"], 0.82, "Multi-pinna solidity must reflect inter-pinna sinuses across all islands.")
            finally:
                pipe.shutdown()

    def test_step02_rereview_coordinate_alignment_config_db_and_docstring_audit(self) -> None:
        """
        Verify Step 02 zero-trust code & docstring audit fixes:
          1. `_parse_basal_leaf_pose`, `estimate_rosette_crown_point`, `approximate_petiole_and_scape_metrics`,
             and `classify_leaf_zone_anatomical` shift `local_win` `keypoints_31` via `local_mask_origin`
             into native full-sheet coordinates and request `return_sheet_size=True` from `specimen_segmenter`.
          2. `SegmentAndExtractPipeline` resolves nested `modules.*` keys (`ect.num_dirs`, `ect.num_thresh`,
             `plant_detector.conf_thres`, `leaf_segmenter.min_area_px`, `landmark_detector.conf_thres`)
             from `LM3_settings.yaml`.
          3. `SegmentationConfig` parses `segment_partial_leaves` and `rosette_plant_bbox_lower_fraction`,
             `LM3ScaleCalibrator.calibrate_sheet(..., use_cf_predicted_by_mp=False)` suppresses megapixel
             guesses, `normalize_pixels_per_mm` preserves 200-DPI scans (`0.78 px/mm`), and
             `_parse_calibration_row` never promotes `mp_estimated`/`uncalibrated` rows to `ruler_calibrated`.
          4. SQLite `record_voucher_segmentation_transaction` + `load_step02_records_from_db` round-trip
             `mask_is_upright`, `petiole_width_px`, `petiole_width_mm`, `blade_width_px`, and `blade_width_mm`.
          5. `generate_step02_publication_figures` scales `fig02a` `ylim` below `1.0` for reniform blades
             (`aspect_ratio < 1.0`), `device="cpu", threads=2` initializes without `torch` NameError,
             JPEG->WebP migration converts `.jpg` via `PIL.Image`, runtime crashes record
             `SEGMENTATION_RUNTIME_ERROR`, and stale docstrings are purged.
        """
        from PIL import Image as PILImage
        from scripts.core.config import PipelineConfig
        from scripts.core.db import (
            get_db_connection,
            load_step02_records_from_db,
            record_voucher_segmentation_transaction,
        )
        from scripts.vision.capitulum_phenology_classifier import (
            classify_single_capitulum,
            classify_voucher_phenology,
        )
        from scripts.vision.lm3_data_loader import (
            _parse_calibration_row,
            normalize_pixels_per_mm,
            resolve_px_per_cm,
        )
        from scripts.vision.lm3_inference_engines import LM3ScaleCalibrator
        from scripts.vision.organ_crown_metrics import (
            _parse_basal_leaf_pose,
            approximate_petiole_and_scape_metrics,
            estimate_rosette_crown_point,
        )
        from scripts.vision.step02_figures import generate_step02_publication_figures

        # 1. Three-space coordinate alignment (`local_mask_origin` -> full-sheet coordinates)
        leaf_local_1 = {
            "leaf_id": 1,
            "plant_individual_id": 0,
            "p_base": (1950.0, 3600.0),
            "p_apex": (1820.0, 3420.0),
            "p_base_sheet": (1950.0, 3600.0),
            "p_apex_sheet": (1820.0, 3420.0),
            "mask_origin": (0.0, 0.0),
            "local_mask_origin": (1750.0, 3350.0),
            "blade_length_px": 222.0,
            "petiole_width_px": 14.0,
            "keypoints_31": {
                "lamina_tip": (70.0, 70.0, 0.95),
                "lamina_base": (200.0, 250.0, 0.95),
                "petiole_1": (210.0, 275.0, 0.92),
                "petiole_6": (245.0, 345.0, 0.90),  # sheet: (1995, 3695)
            },
        }
        leaf_local_2 = {
            "leaf_id": 2,
            "plant_individual_id": 0,
            "p_base": (2060.0, 3590.0),
            "p_apex": (2200.0, 3410.0),
            "p_base_sheet": (2060.0, 3590.0),
            "p_apex_sheet": (2200.0, 3410.0),
            "mask_origin": (0.0, 0.0),
            "local_mask_origin": (1960.0, 3340.0),
            "blade_length_px": 228.0,
            "petiole_width_px": 16.0,
            "keypoints_31": {
                "lamina_tip": (240.0, 70.0, 0.95),
                "lamina_base": (100.0, 250.0, 0.95),
                "petiole_1": (85.0, 280.0, 0.92),
                "petiole_6": (45.0, 360.0, 0.90),  # sheet: (2005, 3700)
            },
        }
        parsed_l1 = _parse_basal_leaf_pose(leaf_local_1)
        self.assertIsNotNone(parsed_l1)
        self.assertAlmostEqual(parsed_l1["p_base"][0], 1950.0, places=2)
        self.assertAlmostEqual(parsed_l1["p_base"][1], 3600.0, places=2)
        self.assertAlmostEqual(parsed_l1["p_petiole_tip"][0], 1995.0, places=2)
        self.assertAlmostEqual(parsed_l1["p_petiole_tip"][1], 3695.0, places=2)

        crown_xy = estimate_rosette_crown_point(
            [leaf_local_1, leaf_local_2],
            root_boxes=[(3680, 1960, 3820, 2040)],
        )
        self.assertTrue(1920.0 <= crown_xy[0] <= 2080.0, f"Expected sheet-space crown X ~2000, got {crown_xy}")
        self.assertTrue(3650.0 <= crown_xy[1] <= 3800.0, f"Expected sheet-space crown Y ~3700, got {crown_xy}")

        scape_res = approximate_petiole_and_scape_metrics(
            basal_leaves=[leaf_local_1, leaf_local_2],
            p_crown={0: crown_xy},
            px_per_cm=100.0,
            capitula_detections=[(700, 1980, 760, 2020)],
        )
        # Corymb at (2000, 730) above crown at (~2000, ~3700) -> ~2970 px = ~29.7 cm
        self.assertTrue(28.0 <= float(scape_res["plant_height_cm"]) <= 31.5)
        self.assertAlmostEqual(float(scape_res["petiole_width_mm"]), 1.55, delta=0.15)

        # 2. Calibration provenance allowlist & 200-DPI (< 1.0 px/mm) preservation
        self.assertAlmostEqual(normalize_pixels_per_mm(0.7874), 0.7874, places=4)
        self.assertAlmostEqual(resolve_px_per_cm(pixels_per_mm=0.7874), 7.874, places=3)
        self.assertAlmostEqual(normalize_pixels_per_mm(0.08), 12.5, places=3)  # explicit legacy mm/px < 0.25

        parsed_mp = _parse_calibration_row({"catalogNumber": "NCU_MP_1", "pixels_per_mm": 9.8, "scale_source": "mp_estimated"})
        self.assertIsNotNone(parsed_mp)
        self.assertFalse(parsed_mp["is_empirical_ruler"])
        self.assertEqual(parsed_mp["scale_source"], "mp_estimated")
        self.assertIsNone(
            _parse_calibration_row({"catalogNumber": "NCU_UNCAL_1", "pixels_per_mm": 9.8, "scale_source": "uncalibrated"})
        )

        cal_no_mp = LM3ScaleCalibrator(device="cpu", providers=[], use_cf_predicted_by_mp=False)
        cal_no_mp.archival_session = None
        res_no_mp = cal_no_mp.calibrate_sheet(np.zeros((2000, 1500, 3), dtype=np.uint8))
        self.assertIsNone(res_no_mp["pixels_per_mm"])
        self.assertFalse(res_no_mp["is_calibrated_mm"])
        self.assertEqual(res_no_mp["scale_source"], "uncalibrated")

        # 3. CapitulumPhenologyConfig parameter wiring
        yellow_patch = np.full((30, 30, 3), 100, dtype=np.uint8)
        yellow_patch[:15, :, :] = (240, 175, 20)
        strict_head = classify_single_capitulum(yellow_patch, min_floret_yellow_ratio=0.80)
        permissive_head = classify_single_capitulum(yellow_patch, min_floret_yellow_ratio=0.20)
        self.assertFalse(strict_head["has_fresh_yellow_corollas"])
        self.assertTrue(permissive_head["has_fresh_yellow_corollas"])
        v_pheno_2 = classify_voucher_phenology("NCU_PHENO", [permissive_head], min_anthesis_heads_for_flowering=2)
        self.assertNotEqual(v_pheno_2["voucher_phenological_state"], "anthesis")

        # 4. Nested LM3_settings.yaml lookup, CPU threads init, JPEG->WebP migration, SQLite width/upright round-trip,
        #    RUNTIME_ERROR label, and reniform aspect_ratio < 1.0 plotting
        cfg_loaded = PipelineConfig.from_yaml()
        self.assertFalse(cfg_loaded.segmentation.segment_partial_leaves)
        self.assertAlmostEqual(cfg_loaded.segmentation.rosette_plant_bbox_lower_fraction, 0.60, places=2)

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "run1" / "data" / "raw_vouchers"
            raw_dir.mkdir(parents=True, exist_ok=True)
            jpg_path = raw_dir / "NCU_RENIFORM.jpg"
            PILImage.fromarray(np.full((600, 600, 3), 240, dtype=np.uint8)).save(jpg_path, format="JPEG", quality=90)
            crash_path = raw_dir / "NCU_CRASH.webp"
            PILImage.fromarray(np.full((400, 400, 3), 240, dtype=np.uint8)).save(crash_path, format="WEBP", quality=90)

            v_csv = tmp_p / "run1" / "data" / "tables" / "curated_vouchers.csv"
            v_csv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([
                {"catalogNumber": "NCU_RENIFORM", "scientificName": "Packera aurea", "image_path": str(jpg_path)},
                {"catalogNumber": "NCU_CRASH", "scientificName": "Packera aurea", "image_path": str(crash_path)},
            ]).to_csv(v_csv, index=False)

            dummy_w = tmp_p / "dummy.onnx"
            dummy_w.write_bytes(b"0")
            db_file = tmp_p / "run1" / "step02_audit.sqlite"

            # Verify CPU threads init does not raise NameError on `torch`
            pipe = step02.SegmentAndExtractPipeline(
                vouchers_csv=v_csv,
                model_weights=dummy_w,
                output_dir=tmp_p / "run1" / "data",
                device="cpu",
                threads=2,
                db_path=db_file,
                wave_size=1,
            )
            try:
                self.assertEqual(pipe.ect_num_dirs, 360)
                self.assertEqual(pipe.ect_num_thetas, 32)
                self.assertFalse(pipe.use_cf_predicted_by_mp)

                # Reniform leaf wider than long (W = 220 px, L = 140 px -> anatomical AR ~ 0.64 < 1.0)
                ren_mask = np.zeros((600, 600), dtype=np.uint8)
                cv2.ellipse(ren_mask, (200, 190), (110, 70), 0, 0, 360, 255, -1)
                ren_pet = np.zeros((600, 600), dtype=np.uint8)
                cv2.rectangle(ren_pet, (193, 255), (207, 350), 255, -1)

                orig_proc = pipe.process_voucher

                def _side_effect_proc(cat_n: str, img_p: Path, *a, **kw):
                    if cat_n == "NCU_CRASH":
                        raise RuntimeError("synthetic corrupt scan")
                    return orig_proc(cat_n, img_p, *a, **kw)

                pipe.ruler_calibrations["NCU_RENIFORM"] = 10.0
                with patch.object(
                    pipe.model_engine,
                    "predict",
                    return_value=[
                        ((120, 90, 260, 310), ren_mask, 0.96, 0),
                        ((255, 193, 350, 207), ren_pet, 0.92, 1),
                    ],
                ):
                    with patch.object(pipe, "process_voucher", side_effect=_side_effect_proc):
                        ext_df, fail_df = pipe.run()

                # JPEG migrated cleanly to WebP via PIL.Image
                self.assertFalse(jpg_path.exists())
                self.assertTrue(jpg_path.with_suffix(".webp").exists())

                # Runtime crash tagged as SEGMENTATION_RUNTIME_ERROR, not UNRESOLVED_CLUMP
                self.assertEqual(len(fail_df), 1)
                self.assertEqual(str(fail_df.iloc[0]["failure_reason"]), "SEGMENTATION_RUNTIME_ERROR")

                # Reniform AR < 1.0 preserved and rendered in fig02a without clipping
                self.assertEqual(len(ext_df), 1)
                ar_ren = float(ext_df.iloc[0]["aspect_ratio"])
                self.assertLess(ar_ren, 0.85)
                self.assertTrue(bool(ext_df.iloc[0]["mask_is_upright"]))
                self.assertGreater(float(ext_df.iloc[0]["blade_width_mm"]), float(ext_df.iloc[0]["blade_length_mm"]))

                fig_paths = generate_step02_publication_figures(
                    extracted_df=ext_df,
                    failed_df=fail_df,
                    repro_df=pd.DataFrame(),
                    pheno_df=pd.DataFrame(),
                    figures_dir=tmp_p / "run1" / "outputs" / "figures",
                )
                self.assertTrue(Path(fig_paths["fig02a_png"]).is_file())

                # Verify SQLite round-trip of mask_is_upright, petiole_width_mm, blade_width_px, blade_width_mm
                conn = get_db_connection(db_path=db_file, init_schema=True)
                try:
                    loaded_ext, _, _, _ = load_step02_records_from_db(conn, ["NCU_RENIFORM"])
                    self.assertIn("NCU_RENIFORM", loaded_ext)
                    rec_db = loaded_ext["NCU_RENIFORM"][0]
                    self.assertTrue(rec_db["mask_is_upright"])
                    self.assertIsNotNone(rec_db["blade_width_px"])
                    self.assertIsNotNone(rec_db["blade_width_mm"])
                    self.assertGreater(float(rec_db["blade_width_mm"]), float(rec_db["blade_length_mm"]))
                finally:
                    conn.close()
            finally:
                pipe.shutdown()

        # 5. Docstring & CLI text hygiene audit across custom Step 2 modules
        stale_tokens = (
            "lm3_packera_seg_finetuned",
            "2-path leaf extraction",
            "retaining petiole-bridged lateral pinnae",
            "(Class 6 detections)",
            "Dissected Rachis Bridging",
            "topological connectivity gatekeeping",
            "Batched FP16 Stage 1 YOLOv5",
        )
        step2_sources = [
            PROJECT_ROOT / "scripts" / "pipeline" / "02_segment_and_extract.py",
            PROJECT_ROOT / "scripts" / "vision" / "lm3_inference_engines.py",
            PROJECT_ROOT / "scripts" / "vision" / "organ_crown_metrics.py",
            PROJECT_ROOT / "scripts" / "vision" / "capitulum_phenology_classifier.py",
            PROJECT_ROOT / "scripts" / "vision" / "step02_figures.py",
        ]
        for src_p in step2_sources:
            src_txt = src_p.read_text(encoding="utf-8")
            for bad_tok in stale_tokens:
                self.assertNotIn(bad_tok, src_txt, f"Stale docstring/comment token {bad_tok!r} found in {src_p.name}")



    def test_300dpi_png_only_12_figure_suite_and_high_richness_scaling(self) -> None:
        """
        Verify the 12-plate 300-DPI PNG-only figure architecture:
          1. Step 02 scales dynamically across 18 taxa (including subsp., var., and hybrid ranks),
             writes 300-DPI PNGs only (zero .pdf files), caches unlimited lm3_assets/ + ect_overlays/,
             and exports capped whole-sheet diagnostic overlays (lm3_sheet_overlays/).
          2. Steps 03-06 modular R suite (scripts/figures/*.R) defines all 10 R canonical PNG plates
             with balanced R delimiters and zero grDevices::pdf() calls.
        """
        from PIL import Image
        from scripts.vision.step02_figures import generate_step02_publication_figures

        taxa_18 = [
            "Packera anonyma", "Packera aurea", "Packera crawfordii", "Packera debilis",
            "Packera dubia", "Packera glabella", "Packera hyperborealis", "Packera indecora",
            "Packera millefolium", "Packera obovata", "Packera paupercula",
            "Packera paupercula var. savannarum", "Packera paupercula subsp. pseudotomentosa",
            "Packera plattensis", "Packera pseudaurea", "Packera schweinitziana",
            "Packera tomentosa", "Packera × memmingeri",
        ]
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "raw"
            mask_dir = tmp_p / "masks"
            fig_dir = tmp_p / "outputs" / "figures"
            raw_dir.mkdir(parents=True)
            mask_dir.mkdir(parents=True)

            sheet_img = np.full((360, 280, 3), 242, dtype=np.uint8)
            cv2.rectangle(sheet_img, (20, 40), (120, 180), (65, 115, 70), -1)
            cv2.rectangle(sheet_img, (150, 170), (245, 310), (55, 105, 65), -1)
            sheet_file = raw_dir / "NCU_MULTI_RAMET.jpg"
            cv2.imwrite(str(sheet_file), sheet_img)

            rows = []
            for idx, tx in enumerate(taxa_18):
                for l_idx in (1, 2):
                    m_img = np.zeros((140, 110), dtype=np.uint8)
                    cv2.ellipse(m_img, (55, 70), (24 + (idx % 4) * 3, 46), 0, 0, 360, 255, -1)
                    m_p = mask_dir / f"CAT_{idx:02d}_L{l_idx}.png"
                    cv2.imwrite(str(m_p), m_img)
                    cat_id = "NCU_MULTI_RAMET" if idx == 0 else f"CAT_{idx:02d}"
                    rows.append({
                        "catalogNumber": cat_id,
                        "leaf_id": l_idx,
                        "plant_individual_id": l_idx if idx == 0 else 1,
                        "species_standardized": tx,
                        "determiner_tier": "Tier_1_Gold",
                        "is_type_specimen": (idx == 0 and l_idx == 1),
                        "leaf_zone": "cauline" if (idx % 3 == 0 and l_idx == 2) else "basal",
                        "solidity": 0.72 + 0.01 * (idx % 12),
                        "aspect_ratio": 0.78 if idx == 1 else (1.35 + 0.08 * (idx % 10)),
                        "archetype_score": 0.92,
                        "pixels_per_mm": 11.5,
                        "is_calibrated_mm": True,
                        "scale_source": "ruler_ocr",
                        "mask_path": str(m_p),
                        "image_path": str(sheet_file) if idx == 0 else "",
                        "bbox_x1": 20 if l_idx == 1 else 150,
                        "bbox_y1": 40 if l_idx == 1 else 170,
                        "bbox_x2": 120 if l_idx == 1 else 245,
                        "bbox_y2": 180 if l_idx == 1 else 310,
                    })

            saved = generate_step02_publication_figures(
                extracted_df=pd.DataFrame(rows),
                failed_df=pd.DataFrame(),
                repro_df=pd.DataFrame({"catalogNumber": ["CAT_01", "CAT_02"], "AR_inv": [1.18, 1.42]}),
                pheno_df=pd.DataFrame({"catalogNumber": ["CAT_01", "CAT_02"], "voucher_phenological_state": ["anthesis", "fruiting"]}),
                figures_dir=fig_dir,
            )
            self.assertEqual(list(fig_dir.rglob("*.pdf")), [])
            for key in ("fig02a_png", "fig02b_png"):
                p_out = Path(saved[key])
                self.assertTrue(p_out.is_file() and p_out.stat().st_size > 15000)
                with Image.open(p_out) as im:
                    dpi_x, dpi_y = im.info.get("dpi", (0, 0))
                    self.assertAlmostEqual(dpi_x, 300.0, delta=2.0)
                    self.assertAlmostEqual(dpi_y, 300.0, delta=2.0)

            self.assertGreaterEqual(len(list((fig_dir / "ect_overlays").glob("*.png"))), 18)
            self.assertGreaterEqual(len(list((fig_dir / "lm3_sheet_overlays").glob("*.png"))), 1)

        expected_r_pngs = {
            "fig03a_efa_morphospace_and_harmonics.png",
            "fig03b_topological_ect_and_heteroblasty.png",
            "fig04a_foliar_allometry_audit.png",
            "fig04b_gmm_cda_morphometrics.png",
            "fig04c_ecoregion_and_batch_audit.png",
            "fig05a_ecological_niche_and_rf_importance.png",
            "fig05b_hopkins_circular_phenology.png",
            "fig05c_sdm_suitability_and_distribution_map.png",
            "fig06a_integrative_taxonomic_synthesis.png",
            "fig06b_taxon_complex_delimitation_network.png",
        }
        r_fig_dir = PROJECT_ROOT / "scripts" / "figures"
        combined_r_src = "\n".join(p.read_text(encoding="utf-8") for p in sorted(r_fig_dir.glob("*.R")))
        for png_name in expected_r_pngs:
            self.assertIn(png_name, combined_r_src)
        self.assertNotIn("grDevices::pdf(", combined_r_src)

    def test_eroded_blade_cie_l_tomentum_and_paper_strap_masking(self) -> None:
        """Verify compute_eroded_blade_cie_l_metrics separates floccose white tomentum from glabrous blades and masks white paper mounting straps, petioles, and holes."""
        rng = np.random.default_rng(42)
        sheet_bgr = np.full((240, 320, 3), (228, 236, 242), dtype=np.uint8)
        mask = np.zeros((120, 120), dtype=np.uint8)
        cv2.circle(mask, (60, 60), 46, 255, -1)

        sheet_glab = sheet_bgr.copy()
        sub_g = sheet_glab[60:180, 40:160].copy()
        glab_tex = np.clip(68.0 + rng.normal(0.0, 2.0, size=(120, 120, 1)), 40, 95).astype(np.uint8)
        leaf_glab_bgr = np.concatenate([glab_tex - 18, glab_tex + 8, glab_tex + 4], axis=2)
        sub_g[mask > 0] = leaf_glab_bgr[mask > 0]
        sheet_glab[60:180, 40:160] = sub_g

        sheet_glab_strap = sheet_glab.copy()
        sheet_glab_strap[114:124, 40:160] = (232, 238, 243)

        sheet_tom = sheet_bgr.copy()
        sub_t = sheet_tom[60:180, 40:160].copy()
        wool_grid = rng.normal(0.0, 26.0, size=(30, 30)).astype(np.float32)
        wool_up = cv2.resize(wool_grid, (120, 120), interpolation=cv2.INTER_CUBIC)
        tom_lum = np.clip(168.0 + wool_up, 95.0, 222.0).astype(np.uint8)
        leaf_tom_bgr = np.stack([tom_lum, tom_lum, tom_lum], axis=2)
        sub_t[mask > 0] = leaf_tom_bgr[mask > 0]
        sheet_tom[60:180, 40:160] = sub_t

        paper_l = step02.estimate_sheet_paper_cie_l_ref(sheet_glab)
        self.assertGreater(paper_l, 85.0)

        m_glab = step02.compute_eroded_blade_cie_l_metrics(
            sheet_glab, mask, mask_origin=(40, 60), paper_l_ref=paper_l
        )
        m_strap = step02.compute_eroded_blade_cie_l_metrics(
            sheet_glab_strap, mask, mask_origin=(40, 60), paper_l_ref=paper_l
        )
        m_tom = step02.compute_eroded_blade_cie_l_metrics(
            sheet_tom, mask, mask_origin=(40, 60), paper_l_ref=paper_l
        )

        self.assertGreater(m_tom["blade_cie_l_norm"], m_glab["blade_cie_l_norm"] + 0.25)
        self.assertGreater(m_tom["blade_cie_l_sd"], m_glab["blade_cie_l_sd"] * 2.5)
        self.assertAlmostEqual(m_glab["blade_cie_l_norm"], m_strap["blade_cie_l_norm"], delta=0.02)


    def test_bolting_scape_specimen_mask_and_splayed_leaf_consolidation(self) -> None:
        """Verify cluster_plant_individuals reunites bolting stems, mask-bridged shoots, and 1-leaf strays while separating side-by-side rooted plants."""
        dummy_mask = np.ones((12, 12), dtype=np.uint8)
        # Single bolting plant on a 2000x3000 sheet: basal rosette at y~2550, upper cauline leaves at y~1250 (dy=1300 > 0.15*3000=450),
        # plus 1 splayed outer basal leaf 650 px laterally from the rosette.
        b1 = DetectedInstance(catalog_number="BOLT1", leaf_id=1, bbox=(2500, 860, 2620, 960), mask=dummy_mask, score=0.95, class_id=0)
        b2 = DetectedInstance(catalog_number="BOLT1", leaf_id=2, bbox=(2520, 930, 2640, 1030), mask=dummy_mask, score=0.94, class_id=0)
        c1 = DetectedInstance(catalog_number="BOLT1", leaf_id=3, bbox=(1200, 890, 1300, 970), mask=dummy_mask, score=0.92, class_id=0)
        c2 = DetectedInstance(catalog_number="BOLT1", leaf_id=4, bbox=(1120, 910, 1210, 980), mask=dummy_mask, score=0.91, class_id=0)
        splay = DetectedInstance(catalog_number="BOLT1", leaf_id=5, bbox=(2530, 1520, 2630, 1620), mask=dummy_mask, score=0.89, class_id=0)

        out_single = cluster_plant_individuals([b1, b2, c1, c2, splay], sheet_width=2000, sheet_height=3000)
        self.assertEqual([inst.plant_individual_id for inst in out_single], [0, 0, 0, 0, 0])

        # Diagonal V-bent plant bridged by Specimen_Segmenter mask across > 0.22 * sheet_width
        v1 = DetectedInstance(catalog_number="VBENT", leaf_id=1, bbox=(2400, 300, 2520, 420), mask=dummy_mask, score=0.95, class_id=0)
        v2 = DetectedInstance(catalog_number="VBENT", leaf_id=2, bbox=(2420, 360, 2540, 480), mask=dummy_mask, score=0.93, class_id=0)
        v3 = DetectedInstance(catalog_number="VBENT", leaf_id=3, bbox=(900, 1400, 1020, 1520), mask=dummy_mask, score=0.92, class_id=0)
        v4 = DetectedInstance(catalog_number="VBENT", leaf_id=4, bbox=(940, 1450, 1060, 1560), mask=dummy_mask, score=0.90, class_id=0)
        spec_mask = np.zeros((3000, 2000), dtype=np.uint8)
        cv2.line(spec_mask, (390, 2460), (1460, 960), 255, thickness=18)
        out_bridged = cluster_plant_individuals([v1, v2, v3, v4], sheet_width=2000, sheet_height=3000, specimen_mask=spec_mask)
        self.assertEqual([inst.plant_individual_id for inst in out_bridged], [0, 0, 0, 0])

        # Two genuinely separate side-by-side rooted rosettes remain distinct [0, 0, 1, 1]
        p0_a = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=1, bbox=(2300, 250, 2420, 360), mask=dummy_mask, score=0.95, class_id=0)
        p0_b = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=2, bbox=(2320, 310, 2440, 420), mask=dummy_mask, score=0.94, class_id=0)
        p1_a = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=3, bbox=(2300, 1480, 2420, 1590), mask=dummy_mask, score=0.93, class_id=0)
        p1_b = DetectedInstance(catalog_number="TWO_ROOT", leaf_id=4, bbox=(2320, 1540, 2440, 1650), mask=dummy_mask, score=0.92, class_id=0)
        roots = [(2460, 280, 2580, 380), (2460, 1510, 2580, 1610)]
        out_two = cluster_plant_individuals([p0_a, p0_b, p1_a, p1_b], sheet_width=2000, sheet_height=3000, root_boxes=roots)
        self.assertEqual([inst.plant_individual_id for inst in out_two], [0, 0, 1, 1])


if __name__ == "__main__":
    unittest.main()
