#!/usr/bin/env python3
"""
===============================================================================
Script: test_configure_leafmachine2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit tests for scripts/vision/configure_leafmachine2.py:
      - Validates default template structure
      - Validates Packera_LeafPriority.pth weights assignment
      - Validates PCD confidence threshold (0.20)
      - Validates NMS overlap threshold (0.75)
      - Validates PointRend subdivision steps (5)
      - Validates minimum leaf area (500)
===============================================================================
"""

import unittest
from pathlib import Path
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts.vision.configure_leafmachine2 import (
    get_default_lm2_template,
    generate_high_performance_config,
    resolve_pcd_weights,
    DEFAULT_PCD_WEIGHTS_NAME,
    DEFAULT_PCD_CONFIDENCE,
    DEFAULT_NMS_THRESH,
    DEFAULT_SUBDIVISION_STEPS,
    DEFAULT_MIN_LEAF_AREA,
)


class TestConfigureLeafMachine2(unittest.TestCase):
    """Test suite for LeafMachine2 configuration generator."""

    def test_default_template_structure(self):
        """Ensure canonical LeafMachine2 template contains required top-level sections."""
        tpl = get_default_lm2_template()
        self.assertIn("leafmachine", tpl)
        lm = tpl["leafmachine"]
        self.assertIn("project", lm)
        self.assertIn("plant_component_detector", lm)
        self.assertIn("leaf_segmentation", lm)
        self.assertIn("overlay", lm)

    def test_packera_rosette_thresholds(self):
        """Verify that generate_high_performance_config configures Packera-specific thresholds."""
        cfg = generate_high_performance_config(
            pcd_confidence=0.20,
            nms_thresh=0.75,
            pointrend_subdivision_steps=5,
            min_leaf_area=500,
        )
        lm = cfg["leafmachine"]
        pcd = lm["plant_component_detector"]
        seg = lm["leaf_segmentation"]

        # Plant Component Detector assertions
        self.assertEqual(pcd["detector_weights"], "LeafPriority.pt")
        self.assertEqual(pcd["minimum_confidence_threshold"], 0.20)
        self.assertEqual(pcd["PCD_confidence"], 0.20)

        # Leaf Segmentation / PointRend assertions
        self.assertEqual(seg["segmentation_model"], "Packera_LeafPriority")
        self.assertEqual(seg["NMS_thresh"], 0.75)
        self.assertEqual(seg["pointrend_subdivision_steps"], 5)
        self.assertEqual(seg["min_leaf_area"], 500)
        self.assertEqual(seg["minimum_confidence_threshold"], 0.70)
        self.assertTrue(seg["calculate_elliptic_fourier_descriptors"])
        self.assertEqual(seg["elliptic_fourier_descriptor_order"], 40)

    def test_pcd_weights_resolver(self):
        """Verify that LeafPriority.pt resolves correctly."""
        weights_name, resolved_path = resolve_pcd_weights(DEFAULT_PCD_WEIGHTS_NAME)
        self.assertEqual(weights_name, "LeafPriority.pt")
        self.assertIsNotNone(resolved_path)
        if resolved_path:
            self.assertTrue(resolved_path.exists())

    def test_pointrend_segmentation_artifacts(self):
        """Verify that Packera_LeafPriority segmentation directory contains required PointRend configs."""
        seg_dir = PROJECT_ROOT / "LeafMachine2" / "leafmachine2" / "segmentation" / "models" / "Packera_LeafPriority"
        self.assertTrue(seg_dir.exists())
        self.assertTrue((seg_dir / "cfg_output.yaml").exists())
        self.assertTrue((seg_dir / "metadata.json").exists())
        self.assertTrue((seg_dir / "model_final.pth").exists())


if __name__ == "__main__":
    unittest.main()
