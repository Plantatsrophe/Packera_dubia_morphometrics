#!/usr/bin/env python3
"""
Unit tests for dynamic YAML configuration loader and PipelineConfig dataclass.
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import yaml

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    DEFAULT_CONFIG_PATH,
    PROJECT_ROOT as CONFIG_ROOT,
    PipelineConfig,
    PathsConfig,
    TaxaConfig,
    ThresholdsConfig,
    DissectionConfig,
    FoldDetectionConfig,
    MorphometricsConfig,
    HarvestingConfig,
    SegmentationConfig,
)


class TestPipelineConfig(unittest.TestCase):
    """Test suite for PipelineConfig dataclass and config loader."""

    def test_pipeline_config_from_yaml_default(self):
        """Verify that default YAML loads into structured PipelineConfig dataclass."""
        cfg = PipelineConfig.from_yaml()
        self.assertIsInstance(cfg, PipelineConfig)
        self.assertIsInstance(cfg.paths, PathsConfig)
        self.assertIsInstance(cfg.taxa, TaxaConfig)
        self.assertIsInstance(cfg.thresholds, ThresholdsConfig)
        self.assertIsInstance(cfg.morphometrics, MorphometricsConfig)
        self.assertIsInstance(cfg.harvesting, HarvestingConfig)
        self.assertIsInstance(cfg.segmentation, SegmentationConfig)

    def test_relative_path_resolution(self):
        """Verify all paths in PathsConfig are absolute and anchored to project root."""
        cfg = PipelineConfig.from_yaml()
        for field_name in cfg.paths.__dataclass_fields__:
            path_val = getattr(cfg.paths, field_name)
            self.assertIsInstance(path_val, Path, f"{field_name} should be a Path object")
            self.assertTrue(path_val.is_absolute(), f"{field_name} ({path_val}) must be an absolute path")
            self.assertIn(str(CONFIG_ROOT), str(path_val), f"{field_name} should be anchored to {CONFIG_ROOT}")

    def test_thresholds_values(self):
        """Verify key pipeline quality and segmentation thresholds match expected parameters."""
        cfg = PipelineConfig.from_yaml()
        self.assertAlmostEqual(cfg.thresholds.pcd_conf, 0.72)
        self.assertAlmostEqual(cfg.thresholds.min_solidity, 0.72)
        self.assertAlmostEqual(cfg.thresholds.min_ucs, 0.85)
        self.assertAlmostEqual(cfg.thresholds.cleanlab_cutoff, 0.85)
        self.assertAlmostEqual(cfg.thresholds.min_megapixels, 8.0)
        self.assertAlmostEqual(cfg.thresholds.min_file_size_kb, 500.0)
        self.assertAlmostEqual(cfg.thresholds.min_sharpness_laplacian, 80.0)

        # Dissection parameters (calibrated)
        self.assertIsInstance(cfg.thresholds.dissection, DissectionConfig)
        self.assertTrue(cfg.thresholds.dissection.enabled)
        self.assertAlmostEqual(cfg.thresholds.dissection.min_solidity_dissected, 0.50)
        self.assertAlmostEqual(cfg.thresholds.dissection.sinus_defect_min_depth_ratio, 0.08)
        self.assertEqual(cfg.thresholds.dissection.min_bilateral_sinus_count, 3)
        self.assertIn("Packera paupercula", cfg.thresholds.dissection.taxa_with_lyrate_tendency)
        self.assertIn("Packera plattensis", cfg.thresholds.dissection.taxa_with_lyrate_tendency)
        self.assertIn("Packera paupercula var. paupercula", cfg.thresholds.dissection.taxa_with_lyrate_tendency)
        self.assertIn("Packera paupercula var. savannarum", cfg.thresholds.dissection.taxa_with_lyrate_tendency)

        # Backward compatibility for solidity property
        self.assertAlmostEqual(cfg.thresholds.solidity.default, 0.72)
        self.assertAlmostEqual(cfg.thresholds.solidity.min_dissected, 0.50)
        self.assertEqual(cfg.thresholds.solidity.taxa_with_lyrate_tendency, cfg.thresholds.dissection.taxa_with_lyrate_tendency)

        # Fold detection parameters (calibrated)
        self.assertIsInstance(cfg.thresholds.fold_detection, FoldDetectionConfig)
        self.assertTrue(cfg.thresholds.fold_detection.enabled)
        self.assertAlmostEqual(cfg.thresholds.fold_detection.max_chord_deviation_ratio, 0.035)
        self.assertAlmostEqual(cfg.thresholds.fold_detection.max_folded_aspect_ratio, 0.42)
        self.assertEqual(cfg.thresholds.fold_detection.min_chord_length_px, 150)

    def test_thresholds_safe_fallbacks_when_subkeys_missing(self):
        """Verify PipelineConfig.from_yaml safely falls back when threshold sub-keys are absent."""
        minimal_yaml = (
            "paths:\n"
            "  workspace_root: \".\"\n"
            "taxa:\n"
            "  target_species:\n"
            "    - \"Packera dubia\"\n"
            "thresholds:\n"
            "  pcd_conf: 0.75\n"
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            cfg_file = Path(tmp_dir) / "minimal.yaml"
            cfg_file.write_text(minimal_yaml, encoding="utf-8")
            cfg = PipelineConfig.from_yaml(cfg_file)

            # Standard baseline defaults
            self.assertAlmostEqual(cfg.thresholds.pcd_conf, 0.75)
            self.assertAlmostEqual(cfg.thresholds.min_solidity, 0.72)
            self.assertAlmostEqual(cfg.thresholds.min_ucs, 0.85)

            # Dissection safe defaults
            self.assertTrue(cfg.thresholds.dissection.enabled)
            self.assertAlmostEqual(cfg.thresholds.dissection.min_solidity_dissected, 0.50)
            self.assertAlmostEqual(cfg.thresholds.dissection.sinus_defect_min_depth_ratio, 0.08)
            self.assertEqual(cfg.thresholds.dissection.min_bilateral_sinus_count, 3)
            self.assertIn("Packera paupercula", cfg.thresholds.dissection.taxa_with_lyrate_tendency)

            # Fold detection safe defaults
            self.assertTrue(cfg.thresholds.fold_detection.enabled)
            self.assertAlmostEqual(cfg.thresholds.fold_detection.max_folded_aspect_ratio, 0.42)
            self.assertAlmostEqual(cfg.thresholds.fold_detection.max_chord_deviation_ratio, 0.035)
            self.assertEqual(cfg.thresholds.fold_detection.min_chord_length_px, 150)

    def test_taxa_structure_and_synonyms(self):
        """Verify target species, synonyms dictionary, and outgroups."""
        cfg = PipelineConfig.from_yaml()
        self.assertIn("Packera dubia", cfg.taxa.target_species)
        self.assertIn("Packera tomentosa", cfg.taxa.target_species)
        self.assertIn("Packera anonyma", cfg.taxa.target_species)
        self.assertIn("Packera plattensis", cfg.taxa.target_species)
        self.assertIn("Packera paupercula", cfg.taxa.target_species)

        # Backward-compatible property alias
        self.assertEqual(cfg.taxa.target_taxa, cfg.taxa.target_species)

        # Synonyms mapping
        self.assertIsInstance(cfg.taxa.synonyms, dict)
        self.assertIn("Packera dubia", cfg.taxa.synonyms)
        self.assertIn("Senecio tomentosus", cfg.taxa.synonyms["Packera dubia"])
        self.assertIn("Packera anonyma", cfg.taxa.synonyms)
        self.assertIn("Senecio smallii", cfg.taxa.synonyms["Packera anonyma"])

        # Outgroups
        self.assertIsInstance(cfg.taxa.outgroups, list)
        self.assertGreater(len(cfg.taxa.outgroups), 0)

    def test_morphometrics_parameters(self):
        """Verify Fourier harmonics count and normalization invariants."""
        cfg = PipelineConfig.from_yaml()
        self.assertEqual(cfg.morphometrics.nb_harmonics, 12)
        self.assertEqual(cfg.morphometrics.harmonics, 12)
        self.assertEqual(cfg.morphometrics.num_pcs, 5)
        self.assertEqual(cfg.morphometrics.max_k, 8)
        self.assertEqual(cfg.morphometrics.random_seed, 42)

        norm = cfg.morphometrics.normalization
        self.assertTrue(norm.align_major_axis)
        self.assertTrue(norm.scale_invariant)
        self.assertTrue(norm.rotation_invariant)
        self.assertTrue(norm.start_point_invariant)

    def test_missing_file_raises_file_not_found(self):
        """Verify missing config file raises informative FileNotFoundError."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            missing_file = Path(tmp_dir) / "non_existent_config.yaml"
            with self.assertRaises(FileNotFoundError) as exc_info:
                PipelineConfig.from_yaml(missing_file)
            self.assertIn("Configuration file not found", str(exc_info.exception))
            self.assertIn(str(missing_file), str(exc_info.exception))

    def test_invalid_yaml_syntax_raises_value_error(self):
        """Verify malformed YAML syntax raises ValueError pointing to file path."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_yaml = Path(tmp_dir) / "malformed.yaml"
            bad_yaml.write_text("paths:\n  workspace_root: [unclosed_list", encoding="utf-8")
            with self.assertRaises(ValueError) as exc_info:
                PipelineConfig.from_yaml(bad_yaml)
            self.assertIn("Invalid YAML syntax", str(exc_info.exception))
            self.assertIn(str(bad_yaml), str(exc_info.exception))

    def test_empty_yaml_raises_value_error(self):
        """Verify empty YAML raises informative ValueError."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            empty_yaml = Path(tmp_dir) / "empty.yaml"
            empty_yaml.write_text("", encoding="utf-8")
            with self.assertRaises(ValueError) as exc_info:
                PipelineConfig.from_yaml(empty_yaml)
            self.assertIn("empty or not a valid dictionary mapping", str(exc_info.exception))

    def test_cross_module_imports(self):
        """Verify main.py, 02_segment_and_extract.py, and 01_voucher_harvester.py import PipelineConfig."""
        import importlib
        import main
        seg_module = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        harv_module = importlib.import_module("scripts.data_prep.01_voucher_harvester")

        self.assertTrue(hasattr(main, "PipelineConfig"))
        self.assertIs(main.PipelineConfig, PipelineConfig)

        self.assertTrue(hasattr(seg_module, "PipelineConfig"))
        self.assertIs(seg_module.PipelineConfig, PipelineConfig)

        self.assertTrue(hasattr(harv_module, "PipelineConfig"))
        self.assertIs(harv_module.PipelineConfig, PipelineConfig)

    def test_cli_parsers_reflect_pipeline_config(self):
        """Verify CLI argument parsers have defaults aligned with PipelineConfig."""
        import importlib

        cfg = PipelineConfig.from_yaml()
        seg_module = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        harv_module = importlib.import_module("scripts.data_prep.01_voucher_harvester")

        # Segmentation parser
        seg_parse_args = seg_module.parse_args
        with patch.object(sys, "argv", ["02_segment_and_extract.py"]):
            seg_args = seg_parse_args()
        self.assertEqual(seg_args.min_solidity, cfg.thresholds.min_solidity)
        self.assertEqual(seg_args.min_ucs, cfg.thresholds.min_ucs)
        self.assertEqual(seg_args.score_thresh, cfg.segmentation.score_thresh)

        # Harvester parser
        harv_parser = harv_module.build_cli_parser()
        harv_args = harv_parser.parse_args([])
        self.assertEqual(harv_args.max_records, cfg.harvesting.max_records_per_taxon)
        self.assertEqual(harv_args.min_megapixels, cfg.thresholds.min_megapixels)
        self.assertEqual(harv_args.taxa, cfg.taxa.target_species)


if __name__ == "__main__":
    unittest.main()
