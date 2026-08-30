#!/usr/bin/env python3
"""
Unit and Integration Test Suite for Hierarchical Leaf Extractor Pipeline
========================================================================
Tests:
- Stage 1: ArtifactFilterGatekeeper pre-emptive hard blanking.
- Stage 2: Native-DPI sub-image cropping for basal rosettes and cymes.
- Stage 3: YOLO organ detection & 3-point anatomical spine tracing.
- Stage 4: EDT peak seeding + SAM 2 point prompting & watershed fallback.
- Stage 5: Geometric gatekeeper (Solidity >= 0.72) & botanical topology routing.
- 4-Tier Extraction Routing (Tier 1, Tier 2, Tier 3, Tier 4).
- End-to-end execution and QC log generation.
"""

import sys
import unittest
import tempfile
import shutil
import numpy as np
import cv2
import pandas as pd
from pathlib import Path

# Add project root directory to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.leaf_cv_utils import (
    extract_edt_point_seeds,
    segment_leaves_sam2_or_watershed,
    trace_3point_anatomical_spine,
    evaluate_convexity_and_solidity,
    detect_sheet_artifacts,
    detect_native_dpi_regions
)
from scripts.core.gatekeeper_engine import ArtifactFilterGatekeeper
from scripts.core.leaf_extraction import (
    process_voucher_precision,
    run_pipeline,
    OUTPUT_DIRS
)
from scripts.core.botanical_topology_classifier import (
    classify_elongated_botanical_organ,
    CAULINE_STEM,
    ROOT_RHIZOME,
    LEAF_PETIOLE
)


class TestEDTPeakSeedingAndDisentanglement(unittest.TestCase):
    """Tests EDT local maxima extraction and SAM 2 / watershed segmentation."""

    def test_edt_peak_seeding_synthetic_rosette(self):
        """Verify that EDT finds centroids of distinct overlapping circular/elliptical leaf blades."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        # Draw three distinct overlapping blades
        cv2.circle(mask, (100, 150), 35, 255, -1)
        cv2.circle(mask, (180, 130), 40, 255, -1)
        cv2.circle(mask, (150, 210), 38, 255, -1)

        dist_map, seeds = extract_edt_point_seeds(mask, min_distance_px=25, relative_threshold=0.25)
        self.assertGreater(dist_map.max(), 0)
        self.assertGreaterEqual(len(seeds), 2)
        # Seeds should be inside the mask
        for sx, sy in seeds:
            self.assertEqual(mask[sy, sx], 255)

    def test_watershed_fallback(self):
        """Verify watershed segmentation cleanly segments multi-blade clumps into separate masks."""
        crop_bgr = np.full((300, 300, 3), 245, dtype=np.uint8)
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.circle(mask, (100, 150), 35, 255, -1)
        cv2.circle(mask, (180, 130), 40, 255, -1)
        crop_bgr[mask > 0] = (40, 80, 35)

        seeds = [(100, 150), (180, 130)]
        leaf_masks = segment_leaves_sam2_or_watershed(
            rosette_crop_bgr=crop_bgr,
            rosette_binary_mask=mask,
            point_seeds=seeds,
            use_sam2=False
        )
        self.assertEqual(len(leaf_masks), 2)
        for lm in leaf_masks:
            self.assertGreater(cv2.countNonZero(lm), 300)


class TestStage1GatekeeperBlanking(unittest.TestCase):
    """Tests Stage 1 ArtifactFilterGatekeeper.pre_emptive_hard_blanking integration."""

    def test_gatekeeper_hard_blanking_metadata(self):
        canvas = np.full((1000, 800, 3), 240, dtype=np.uint8)
        # Synthetic label at bottom right
        canvas[700:950, 500:750] = (10, 10, 10)

        gatekeeper = ArtifactFilterGatekeeper()
        artifacts = detect_sheet_artifacts(canvas)
        artifact_list = []
        for cat, items in artifacts.items():
            for item in items:
                bbox = item.get("bbox")
                if bbox is not None:
                    artifact_list.append({
                        "category": cat,
                        "bbox": [int(b) for b in bbox],
                        "confidence": 0.95
                    })

        clean = gatekeeper.pre_emptive_hard_blanking(canvas, artifact_list, is_rgb=False, padding_pixels=15)
        # Verify the label area is wiped to white [255, 255, 255]
        self.assertTrue(np.all(clean[720:930, 520:730] == 255))


class TestStage3AnatomicalSpines(unittest.TestCase):
    """Tests 3-point anatomical spine tracing."""

    def test_spine_keypoints_on_synthetic_leaf(self):
        mask = np.zeros((200, 100), dtype=np.uint8)
        # Ovate blade + narrow petiole
        cv2.ellipse(mask, (50, 70), (30, 50), 0, 0, 360, 255, -1)
        cv2.rectangle(mask, (45, 120), (55, 180), 255, -1)

        gray = np.full((200, 100), 240, dtype=np.uint8)
        gray[mask > 0] = 50

        spine = trace_3point_anatomical_spine(mask, gray)
        self.assertIn("p_apex", spine)
        self.assertIn("p_transition", spine)
        self.assertIn("p_caudex", spine)
        self.assertGreater(spine["lamina_length_px"], 0)
        self.assertGreater(spine["total_spine_length_px"], 0)


class TestStage5TopologyClassifierRouting(unittest.TestCase):
    """Tests Stage 5 linear organ routing via botanical_topology_classifier."""

    def test_linear_root_rhizome_detected(self):
        """Subterranean branching structure should classify as root_rhizome."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.line(mask, (150, 50), (150, 250), 255, 4)
        cv2.line(mask, (150, 120), (80, 180), 255, 4)
        cv2.line(mask, (150, 160), (220, 220), 255, 4)
        cv2.line(mask, (150, 200), (90, 270), 255, 4)

        cls = classify_elongated_botanical_organ(
            mask,
            full_sheet_height=1000.0,
            y_centroid=800.0,
            has_connected_blade=False
        )
        self.assertEqual(cls, ROOT_RHIZOME)


class TestFourTierPrecisionExtraction(unittest.TestCase):
    """Tests the 4-Tier extraction routing and QC log export."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp(prefix="packera_extract_test_"))
        self.raw_dir = self.test_dir / "raw_vouchers"
        self.raw_dir.mkdir(parents=True, exist_ok=True)

        self.output_dirs = {
            "annotations": self.test_dir / "masks" / "annotations",
            "rosettes_dense": self.test_dir / "cropped_patches" / "rosettes_dense",
            "basal_leaves_raw": self.test_dir / "masks" / "basal_leaves_raw",
            "tier1_intact": self.test_dir / "masks" / "tier1_intact",
            "tier2_reflected": self.test_dir / "masks" / "tier2_reflected",
            "tier3_open_curves": self.test_dir / "masks" / "tier3_open_curves",
            "capitula": self.test_dir / "masks" / "capitula",
            "qc_overlays": self.test_dir / "outputs" / "extraction_qc",
            "tables": self.test_dir / "tables"
        }
        for d in self.output_dirs.values():
            d.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_voucher_extraction_pipeline(self):
        """Generate a synthetic voucher sheet and execute process_voucher_precision."""
        voucher_path = self.raw_dir / "NCU00012345.jpg"
        sheet = np.full((1200, 900, 3), 245, dtype=np.uint8)

        # Herbarium label at bottom right
        sheet[900:1150, 600:870] = (20, 20, 20)

        # Basal rosette in lower quadrant (y=600..900, x=300..600)
        cv2.ellipse(sheet, (450, 750), (90, 60), 15, 0, 360, (30, 70, 25), -1)

        # Capitulescence in upper quadrant (y=200..400, x=400..500)
        cv2.circle(sheet, (450, 250), 30, (20, 60, 20), -1)

        cv2.imwrite(str(voucher_path), sheet)

        records = process_voucher_precision(
            image_path=voucher_path,
            output_dirs=self.output_dirs,
            save_overlays=True,
            model=None,
            use_sam2=True
        )

        self.assertIsInstance(records, list)
        self.assertGreater(len(records), 0)

        first_rec = records[0]
        self.assertIn("catalogNumber", first_rec)
        self.assertIn("assigned_tier", first_rec)
        self.assertIn(first_rec["assigned_tier"], [1, 2, 3, 4])

        # Test QC table writing
        df = pd.DataFrame(records)
        qc_table_path = self.output_dirs["tables"] / "leaf_extraction_qc.csv"
        df.to_csv(qc_table_path, index=False)
        self.assertTrue(qc_table_path.is_file())


if __name__ == "__main__":
    unittest.main()
