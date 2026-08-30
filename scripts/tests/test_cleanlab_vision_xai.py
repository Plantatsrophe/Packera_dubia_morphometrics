#!/usr/bin/env python3
"""
scripts/tests/test_cleanlab_vision_xai.py
========================================
Unit and regression tests for DINOv2 Confident Learning & Grad-CAM XAI.
"""

import importlib
import unittest
from pathlib import Path
import pandas as pd
import numpy as np

# Dynamically import module with numeric prefix
mod_xai = importlib.import_module("scripts.analysis.05_cleanlab_vision_xai")
standardize_packera_taxon = mod_xai.standardize_packera_taxon
TARGET_TAXA = mod_xai.TARGET_TAXA


class TestCleanlabVisionXAI(unittest.TestCase):
    """Test suite for DINOv2 vision embedding extraction & Cleanlab XAI."""

    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[2]
        cls.audit_csv = cls.repo_root / "data" / "tables" / "label_noise_audit.csv"
        cls.figure_png = cls.repo_root / "outputs" / "figures" / "GradCAM_audit_panel.png"

    def test_taxonomic_standardization(self):
        """Test species synonymy standardizer for the 4 core complex taxa."""
        self.assertEqual(standardize_packera_taxon("Senecio smallii Britton"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Packera anonyma (Alph.Wood) W.A.Weber"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Senecio tomentosus Michx."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera tomentosa C.Jeffrey"), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. savannarum"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. pseudotomentosa"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Senecio plattensis Nutt."), "Packera plattensis")

    def test_output_artifacts_exist(self):
        """Verify that label noise audit CSV and Grad-CAM figure were generated."""
        self.assertTrue(self.audit_csv.exists(), f"Missing audit CSV at {self.audit_csv}")
        self.assertTrue(self.figure_png.exists(), f"Missing Grad-CAM figure at {self.figure_png}")
        self.assertGreater(self.figure_png.stat().st_size, 50000, "Figure file size unexpectedly small.")

    def test_audit_table_schema_and_values(self):
        """Verify data integrity of the Cleanlab label noise audit table."""
        df = pd.read_csv(self.audit_csv)
        self.assertGreater(len(df), 0, "Audit table is empty.")

        required_cols = [
            "catalogNumber",
            "patch_path",
            "species_raw",
            "species_standardized",
            "label_idx",
            "determiner_tier",
            "given_label",
            "predicted_label",
            "confidence_given_class",
            "confidence_predicted_class",
            "label_quality_score",
            "c_error",
            "is_cleanlab_issue",
            "is_label_corrupted",
            "triage_action",
            "discordance_reason",
        ]
        for col in required_cols:
            self.assertIn(col, df.columns, f"Missing required column '{col}' in audit table.")

        # Check values
        self.assertTrue(df["c_error"].between(0.0, 1.0).all(), "c_error contains values outside [0, 1].")
        self.assertTrue(df["label_quality_score"].between(0.0, 1.0).all(), "label_quality_score outside [0, 1].")

        # Check consistency of corrupted flag (threshold 0.85)
        expected_corrupted = df["c_error"] > 0.85
        np.testing.assert_array_equal(
            df["is_label_corrupted"].values,
            expected_corrupted.values,
            err_msg="is_label_corrupted flag does not match (c_error > 0.85).",
        )

        # Check triage actions
        self.assertTrue(
            df.loc[df["is_label_corrupted"], "triage_action"]
            .str.contains("Prune & Queue", case=False)
            .all(),
            "Corrupted vouchers not properly assigned to Prune & Queue.",
        )


if __name__ == "__main__":
    unittest.main()
