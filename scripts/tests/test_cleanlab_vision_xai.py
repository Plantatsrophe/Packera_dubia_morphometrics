#!/usr/bin/env python3
"""
scripts/tests/test_cleanlab_vision_xai.py
========================================
Unit and regression tests for DINOv2 feature extraction, Confident Learning thresholds,
taxonomic standardization, and Grad-CAM XAI.
"""

import importlib
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

# Dynamically import module with numeric prefix
mod_xai = importlib.import_module("scripts.analysis.05_cleanlab_vision_xai")
standardize_packera_taxon = mod_xai.standardize_packera_taxon
extract_dinov2_embeddings = mod_xai.extract_dinov2_embeddings
run_confident_learning_audit = mod_xai.run_confident_learning_audit
compute_out_of_fold_probabilities = mod_xai.compute_out_of_fold_probabilities
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
        self.assertEqual(standardize_packera_taxon(None), "Unknown")

    def test_dinov2_feature_extraction_synthetic(self):
        """Test DINOv2 feature extraction returns (N, 768) embeddings for synthetic records."""
        synthetic_records = [
            {
                "catalogNumber": f"SYNTH_{i:03d}",
                "patch_path": f"/nonexistent/path/patch_{i}.jpg",
                "label_idx": i % 4,
                "taxon": TARGET_TAXA[i % 4],
            }
            for i in range(12)
        ]

        features, labels, cat_nums = extract_dinov2_embeddings(
            records=synthetic_records,
            model_name="dinov2_vitb14",
            device="cpu",
            batch_size=4,
        )

        self.assertEqual(features.shape, (12, 768))
        self.assertEqual(len(labels), 12)
        self.assertEqual(len(cat_nums), 12)
        self.assertEqual(cat_nums[0], "SYNTH_000")
        self.assertEqual(features.dtype, np.float32)

    def test_confident_learning_thresholds_and_triage(self):
        """Test Confident Learning audit flags label noise based on error thresholds."""
        # 4 records:
        # Record 0: Given 0, model predicts 1 with high confidence (0.95), c_error = 0.95 -> corrupted
        # Record 1: Given 1, model predicts 1 with high confidence (0.90), c_error = 0.10 -> clean
        # Record 2: Given 2, model predicts 2 with moderate conf (0.60), c_error = 0.40 -> clean at 0.85, issue at 0.50
        # Record 3: Given 3, model predicts 0 with conf (0.88), c_error = 0.88 -> corrupted
        pred_probs = np.array([
            [0.05, 0.95, 0.00, 0.00],
            [0.05, 0.90, 0.03, 0.02],
            [0.20, 0.20, 0.60, 0.00],
            [0.88, 0.00, 0.00, 0.12],
        ])
        labels = np.array([0, 1, 2, 3])

        records_df = pd.DataFrame({
            "catalogNumber": ["CAT001", "CAT002", "CAT003", "CAT004"],
            "patch_path": ["/p1", "/p2", "/p3", "/p4"],
            "species_raw": ["Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"],
            "label_idx": [0, 1, 2, 3],
            "determiner_tier": ["Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze", "Tier_1_Gold"],
        })

        # Test with default error_threshold = 0.85
        audit_85 = run_confident_learning_audit(
            pred_probs=pred_probs,
            labels=labels,
            records_df=records_df,
            class_names=TARGET_TAXA,
            error_threshold=0.85,
        )

        self.assertEqual(len(audit_85), 4)
        self.assertTrue(audit_85.loc[0, "is_label_corrupted"])
        self.assertFalse(audit_85.loc[1, "is_label_corrupted"])
        self.assertFalse(audit_85.loc[2, "is_label_corrupted"])
        self.assertTrue(audit_85.loc[3, "is_label_corrupted"])

        # Verify triage actions
        self.assertIn("Prune & Queue", audit_85.loc[0, "triage_action"])
        self.assertIn("Prune & Queue", audit_85.loc[3, "triage_action"])
        self.assertIn("Retain", audit_85.loc[1, "triage_action"])

        # Test with tighter error_threshold = 0.35 -> Record 2 should now be flagged
        audit_35 = run_confident_learning_audit(
            pred_probs=pred_probs,
            labels=labels,
            records_df=records_df,
            class_names=TARGET_TAXA,
            error_threshold=0.35,
        )
        self.assertTrue(audit_35.loc[2, "is_label_corrupted"])

    def test_compute_out_of_fold_probabilities(self):
        """Test Stratified K-Fold out-of-fold probability computation on synthetic features."""
        np.random.seed(42)
        n_samples = 20
        n_features = 768
        n_classes = 4

        # Generate well-separated synthetic features
        labels = np.array([i % n_classes for i in range(n_samples)])
        features = np.random.randn(n_samples, n_features).astype(np.float32)
        for i in range(n_samples):
            features[i, labels[i] * 100 : (labels[i] + 1) * 100] += 5.0

        pred_probs, acc, f1 = compute_out_of_fold_probabilities(
            features=features,
            labels=labels,
            n_splits=2,
            random_state=42,
        )

        self.assertEqual(pred_probs.shape, (n_samples, n_classes))
        self.assertTrue(np.allclose(pred_probs.sum(axis=1), 1.0, atol=1e-4))
        self.assertGreaterEqual(acc, 0.5)

    def test_output_artifacts_exist(self):
        """Verify that label noise audit CSV and Grad-CAM figure exist if produced."""
        if not self.audit_csv.exists() or not self.figure_png.exists():
            self.skipTest("Production audit CSV or Grad-CAM figure not present on this machine.")
        self.assertTrue(self.audit_csv.exists(), f"Missing audit CSV at {self.audit_csv}")
        self.assertTrue(self.figure_png.exists(), f"Missing Grad-CAM figure at {self.figure_png}")
        self.assertGreater(self.figure_png.stat().st_size, 50000, "Figure file size unexpectedly small.")

    def test_audit_table_schema_and_values(self):
        """Verify data integrity of the Cleanlab label noise audit table if present."""
        if not self.audit_csv.exists():
            self.skipTest("Production label_noise_audit.csv not present.")
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
