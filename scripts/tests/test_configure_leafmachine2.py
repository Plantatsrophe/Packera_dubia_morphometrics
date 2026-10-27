#!/usr/bin/env python3
"""
===============================================================================
Script: test_configure_leafmachine2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit test suite for scripts/vision/configure_leafmachine2.py:
      - Validates canonical LeafMachine2 template structure and schema
      - Tests dynamic parameter injection from PipelineConfig (pcd_conf, weights)
      - Tests PointRend, batch size, and worker configuration
      - Verifies weights path resolution and output serialization
===============================================================================
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
import yaml

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig
from scripts.vision.configure_leafmachine2 import (
    DEFAULT_BATCH_SIZE,
    DEFAULT_NUM_WORKERS,
    DEFAULT_PCD_CONFIDENCE,
    DEFAULT_SEG_MODEL_NAME,
    DEFAULT_SUBDIVISION_STEPS,
    generate_high_performance_config,
    get_default_lm2_template,
    load_config_yaml,
    resolve_pcd_weights,
    save_config_yaml,
    validate_lm2_config,
)


class TestConfigureLeafMachine2(unittest.TestCase):
    """Test suite for consolidated LeafMachine2 configuration generator."""

    def test_default_template_structure(self):
        """Ensure canonical LeafMachine2 template contains all required top-level sections."""
        tpl = get_default_lm2_template()
        self.assertIn("leafmachine", tpl)
        lm = tpl["leafmachine"]
        self.assertIn("project", lm)
        self.assertIn("plant_component_detector", lm)
        self.assertIn("leaf_segmentation", lm)
        self.assertIn("overlay", lm)
        self.assertIn("ruler_detection", lm)
        self.assertIn("do", lm)
        self.assertIn("data", lm)
        self.assertIn("modules", lm)

    def test_dynamic_pipeline_config_wiring(self):
        """Verify generate_high_performance_config dynamically pulls from PipelineConfig."""
        pipe_cfg = PipelineConfig.from_yaml()
        cfg = generate_high_performance_config(pipeline_cfg=pipe_cfg)

        valid, err = validate_lm2_config(cfg)
        self.assertTrue(valid, f"Generated config should be valid: {err}")

        lm = cfg["leafmachine"]
        pcd = lm["plant_component_detector"]
        seg = lm["leaf_segmentation"]
        proj = lm["project"]

        # Dynamic thresholds wired from config.yaml (Weaver & Smith 2023/2024 empirical range 0.50-0.70, default 0.65)
        expected_pcd_conf = float(getattr(pipe_cfg.segmentation, "minimum_confidence_threshold", DEFAULT_PCD_CONFIDENCE))
        self.assertGreaterEqual(pcd["minimum_confidence_threshold"], 0.50)
        self.assertLessEqual(pcd["minimum_confidence_threshold"], 0.70)
        self.assertAlmostEqual(pcd["minimum_confidence_threshold"], expected_pcd_conf)
        self.assertAlmostEqual(pcd["PCD_confidence"], expected_pcd_conf)
        self.assertEqual(pcd["detector_type"], "Plant_Detector")
        self.assertEqual(pcd["detector_version"], "PLANT_LeafPriority")
        self.assertEqual(pcd["detector_weights"], "LeafPriority.pt")

        # Empirical LeafMachine2 v2.2+ lobe/leaflet & raw PointRend contour guardrails
        self.assertIs(seg["treat_leaflet_as_leaf_whole"], False)
        self.assertIs(proj["treat_leaflet_as_leaf_whole"], False)
        self.assertIs(seg["use_efds_for_png_masks"], False)
        self.assertIs(seg["do_apply_landmarks"], False)
        self.assertIs(lm["do"]["do_apply_landmarks"], False)
        self.assertIs(lm["landmark_detector"]["do_apply_landmarks"], False)
        self.assertIs(lm["landmark_detector"]["landmark_whole_leaves"], False)

        # Batch size and workers
        self.assertEqual(proj["batch_size"], DEFAULT_BATCH_SIZE)
        self.assertEqual(pcd["batch_size"], DEFAULT_BATCH_SIZE)
        self.assertEqual(pcd["num_workers"], DEFAULT_NUM_WORKERS)
        self.assertEqual(proj["num_workers_seg"], DEFAULT_NUM_WORKERS)
        self.assertEqual(proj["num_workers_ruler"], DEFAULT_NUM_WORKERS)

        # PointRend configuration
        self.assertEqual(seg["segmentation_model"], DEFAULT_SEG_MODEL_NAME)
        self.assertEqual(seg["segmentation_type"], "Detectron2_PointRend")
        self.assertTrue(seg["use_pointrend"])
        self.assertEqual(seg["pointrend_subdivision_steps"], DEFAULT_SUBDIVISION_STEPS)
        self.assertTrue(seg["calculate_elliptic_fourier_descriptors"])
        self.assertEqual(seg["elliptic_fourier_descriptor_order"], 40)

        # Path resolution
        self.assertTrue(Path(proj["dir_images_local"]).is_absolute())
        self.assertTrue(Path(proj["dir_output"]).is_absolute())
        self.assertIn("LM2_Project", proj["dir_images_local"])
        self.assertIn("LM2_Project", proj["dir_output"])

    def test_custom_overrides(self):
        """Verify custom arguments override defaults properly."""
        cfg = generate_high_performance_config(
            pcd_confidence=0.88,
            batch_size=25,
            num_workers=4,
            pointrend_subdivision_steps=7,
        )
        lm = cfg["leafmachine"]
        self.assertAlmostEqual(lm["plant_component_detector"]["minimum_confidence_threshold"], 0.88)
        self.assertEqual(lm["project"]["batch_size"], 25)
        self.assertEqual(lm["plant_component_detector"]["num_workers"], 4)
        self.assertEqual(lm["leaf_segmentation"]["pointrend_subdivision_steps"], 7)

    def test_pcd_weights_resolver(self):
        """Verify resolve_pcd_weights locates target or candidate weights file."""
        weights_name, resolved_path = resolve_pcd_weights()
        self.assertIsNotNone(weights_name)
        self.assertIsNotNone(resolved_path)
        self.assertTrue(resolved_path.is_absolute())

    def test_save_load_and_validate(self):
        """Verify serialization and validation to/from disk."""
        cfg = generate_high_performance_config()
        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = Path(tmpdir) / "test_lm2_config.yaml"
            save_config_yaml(cfg, out_file)
            self.assertTrue(out_file.exists())

            loaded = load_config_yaml(out_file)
            self.assertEqual(loaded["leafmachine"]["project"]["run_name"], "Packera_dubia_LM2")
            valid, err = validate_lm2_config(loaded)
            self.assertTrue(valid, err)


if __name__ == "__main__":
    unittest.main()
