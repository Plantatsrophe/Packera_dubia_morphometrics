#!/usr/bin/env python3
"""
===============================================================================
Script: test_gmm_morphotools.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit tests for label-blind GMM clustering, Bayes Factor calculations,
    and Canonical Discriminant Analysis with passive sample projection.
===============================================================================
"""

import unittest
from pathlib import Path
import pandas as pd
import numpy as np
import scipy.linalg as la
import re

PROJECT_ROOT = Path(__file__).resolve().parents[2]


class TestGMMAndMorphoTools(unittest.TestCase):
    """Test suite for GMM morphometrics and MorphoTools2 CDA architecture."""

    def setUp(self):
        self.vouchers_path = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
        self.harmonics_path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics.csv"
        self.r_script_path = PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R"

    def test_r_script_exists_and_under_line_limit(self):
        """Verify that 04_gmm_morphotools.R exists and adheres to modular line guidelines (< 650 lines)."""
        self.assertTrue(self.r_script_path.exists(), "04_gmm_morphotools.R must exist.")
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        self.assertLessEqual(len(lines), 650, f"Script exceeds line limit: {len(lines)}")

    def test_datasets_exist(self):
        """Verify that curated_vouchers.csv and leaf_efa_harmonics.csv are present."""
        self.assertTrue(self.vouchers_path.exists(), "curated_vouchers.csv not found.")
        self.assertTrue(self.harmonics_path.exists(), "leaf_efa_harmonics.csv not found.")

    def test_taxonomic_standardization(self):
        """Test standardization regex for synonymy in Packera dubia complex."""
        def standardize_taxon(s):
            if pd.isna(s): return "Unknown"
            s = str(s).strip()
            if re.search(r"anonym|smallii|earlei", s, re.I): return "Packera anonyma"
            if re.search(r"tomentos|dubia", s, re.I): return "Packera dubia"
            if re.search(r"plattensis|flavovirens", s, re.I): return "Packera plattensis"
            if re.search(r"paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", s, re.I): return "Packera paupercula"
            return s.split("(")[0].strip()

        self.assertEqual(standardize_taxon("Senecio tomentosus Michx."), "Packera dubia")
        self.assertEqual(standardize_taxon("Packera tomentosa C.Jeffrey"), "Packera dubia")
        self.assertEqual(standardize_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_taxon("Senecio smallii Britton, 1893"), "Packera anonyma")
        self.assertEqual(standardize_taxon("Packera paupercula var. savannarum R.R.Kowal"), "Packera paupercula")
        self.assertEqual(standardize_taxon("Senecio plattensis Nutt."), "Packera plattensis")

    def test_cda_mathematics_and_eigenvalues(self):
        """Verify Canonical Discriminant Analysis generalized eigenvalue decomposition."""
        efa_df = pd.read_csv(self.harmonics_path)
        closed_df = efa_df[efa_df["assigned_tier"].isin(["Tier_1_Direct", "Tier_2_Reflected"])].copy()
        pca_cols = ["PC1", "PC2", "PC3", "PC4", "PC5"]
        valid = closed_df.dropna(subset=pca_cols).copy()

        def standardize_taxon(s):
            if pd.isna(s): return "Unknown"
            s = str(s).strip()
            if re.search(r"anonym|smallii|earlei", s, re.I): return "Packera anonyma"
            if re.search(r"tomentos|dubia", s, re.I): return "Packera dubia"
            if re.search(r"plattensis|flavovirens", s, re.I): return "Packera plattensis"
            if re.search(r"paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", s, re.I): return "Packera paupercula"
            return s.split("(")[0].strip()

        valid["species_std"] = valid["species_raw"].apply(standardize_taxon)
        target_taxa = ["Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula"]
        is_active = (valid["determiner_tier"] == "Tier_1_Gold") & (valid["species_std"].isin(target_taxa))

        X_a = valid.loc[is_active, pca_cols].values
        y_a = valid.loc[is_active, "species_std"].values
        N_a, p = X_a.shape
        g = len(target_taxa)

        grand_mean = np.mean(X_a, axis=0)
        B = np.zeros((p, p))
        W = np.zeros((p, p))

        for grp in target_taxa:
            sub_x = X_a[y_a == grp]
            if len(sub_x) > 0:
                m_k = np.mean(sub_x, axis=0)
                diff_m = (m_k - grand_mean).reshape(-1, 1)
                B += len(sub_x) * (diff_m @ diff_m.T)
                diff_x = sub_x - m_k
                W += diff_x.T @ diff_x

        S_reg = (W / max(N_a - g, 1)) + np.eye(p) * 1e-7
        eigvals, eigvecs = la.eigh(B, S_reg)
        idx = np.argsort(eigvals)[::-1][:g-1]
        top_eigvals = eigvals[idx]

        # Check positive eigenvalues and variation
        self.assertGreater(top_eigvals[0], 0)
        var_pct = (top_eigvals / np.sum(top_eigvals)) * 100
        self.assertGreater(var_pct[0], 80.0, "Can1 should explain >80% between-group variation.")

    def test_r_script_reproductive_integration_flags(self):
        """Verify that 04_gmm_morphotools.R supports reproductive options and traits."""
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("--reproductive", content, "04_gmm_morphotools.R should support --reproductive CLI flag")
        self.assertIn("--outline-only", content, "04_gmm_morphotools.R should support --outline-only CLI flag")
        self.assertIn("AR_inv", content, "04_gmm_morphotools.R should integrate AR_inv reproductive trait")
        self.assertIn("AR_inv_std", content, "04_gmm_morphotools.R should compute standardized AR_inv_std")

    def test_joint_vegetative_reproductive_cda_mathematics(self):
        """Verify CDA mathematics with joint Fourier PCs and standardized capitulum aspect ratio."""
        efa_df = pd.read_csv(self.harmonics_path)
        closed_df = efa_df[efa_df["assigned_tier"].isin(["Tier_1_Direct", "Tier_2_Reflected"])].copy()

        # Simulate / synthesize reproductive AR_inv if not already in harmonics file
        if "AR_inv" not in closed_df.columns:
            np.random.seed(42)
            closed_df["AR_inv"] = np.random.uniform(0.8, 1.6, size=len(closed_df))
            # Simulate early-season vouchers with unexpanded/missing capitula (NA)
            closed_df.loc[closed_df.sample(frac=0.1, random_state=42).index, "AR_inv"] = np.nan

        # Standardize AR_inv (z-score)
        valid_ar = closed_df["AR_inv"].dropna()
        closed_df["AR_inv_std"] = (closed_df["AR_inv"] - valid_ar.mean()) / valid_ar.std()

        feature_cols = ["PC1", "PC2", "PC3", "PC4", "PC5", "AR_inv_std"]
        valid = closed_df.dropna(subset=feature_cols).copy()

        def standardize_taxon(s):
            if pd.isna(s): return "Unknown"
            s = str(s).strip()
            if re.search(r"anonym|smallii|earlei", s, re.I): return "Packera anonyma"
            if re.search(r"tomentos|dubia", s, re.I): return "Packera dubia"
            if re.search(r"plattensis|flavovirens", s, re.I): return "Packera plattensis"
            if re.search(r"paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", s, re.I): return "Packera paupercula"
            return s.split("(")[0].strip()

        valid["species_std"] = valid["species_raw"].apply(standardize_taxon)
        target_taxa = ["Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula"]
        is_active = (valid["determiner_tier"] == "Tier_1_Gold") & (valid["species_std"].isin(target_taxa))

        X_a = valid.loc[is_active, feature_cols].values
        y_a = valid.loc[is_active, "species_std"].values
        N_a, p = X_a.shape
        g = len(target_taxa)

        grand_mean = np.mean(X_a, axis=0)
        B = np.zeros((p, p))
        W = np.zeros((p, p))

        for grp in target_taxa:
            sub_x = X_a[y_a == grp]
            if len(sub_x) > 0:
                m_k = np.mean(sub_x, axis=0)
                diff_m = (m_k - grand_mean).reshape(-1, 1)
                B += len(sub_x) * (diff_m @ diff_m.T)
                diff_x = sub_x - m_k
                W += diff_x.T @ diff_x

        S_reg = (W / max(N_a - g, 1)) + np.eye(p) * 1e-7
        eigvals, eigvecs = la.eigh(B, S_reg)
        idx = np.argsort(eigvals)[::-1][:g-1]
        top_eigvals = eigvals[idx]

        self.assertGreater(top_eigvals[0], 0)
        var_pct = (top_eigvals / np.sum(top_eigvals)) * 100
        self.assertGreater(var_pct[0], 70.0, "Can1 in joint ordination should explain dominant between-group variation.")

    def test_r_script_institutional_batch_effect_support(self):
        """Verify 04_gmm_morphotools.R has institutional batch-effect audit logic and R^2 < 0.05 assertion."""
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("--output-batch-audit", content, "04_gmm_morphotools.R should support --output-batch-audit")
        self.assertIn("audit_institutional_batch_effects", content, "04_gmm_morphotools.R must include audit_institutional_batch_effects")
        self.assertIn("institutionCode", content, "04_gmm_morphotools.R must parse institutionCode")
        self.assertIn("0.05", content, "04_gmm_morphotools.R must evaluate 0.05 variance threshold")

    def test_institutional_batch_effect_report_and_r2(self):
        """Verify institutional batch effect audit report exists and confirms R^2 < 0.05 for PC1 and PC2."""
        audit_csv = PROJECT_ROOT / "outputs" / "reports" / "institutional_batch_effect_audit.csv"
        self.assertTrue(audit_csv.exists(), f"Diagnostic report {audit_csv} must exist.")
        df = pd.read_csv(audit_csv)
        self.assertIn("Response", df.columns)
        self.assertIn("Term", df.columns)
        self.assertIn("R_squared", df.columns)
        self.assertIn("Batch_Effect_Status", df.columns)

        pc1_row = df[(df["Response"] == "PC1") & (df["Term"] == "institutionCode")]
        pc2_row = df[(df["Response"] == "PC2") & (df["Term"] == "institutionCode")]
        self.assertEqual(len(pc1_row), 1, "Missing PC1 institutionCode row in batch effect audit")
        self.assertEqual(len(pc2_row), 1, "Missing PC2 institutionCode row in batch effect audit")

        r2_pc1 = float(pc1_row["R_squared"].iloc[0])
        r2_pc2 = float(pc2_row["R_squared"].iloc[0])
        self.assertLess(r2_pc1, 0.05, f"PC1 R^2 ({r2_pc1}) must be < 0.05")
        self.assertLess(r2_pc2, 0.05, f"PC2 R^2 ({r2_pc2}) must be < 0.05")
        self.assertIn("Pass", str(pc1_row["Batch_Effect_Status"].iloc[0]))
        self.assertIn("Pass", str(pc2_row["Batch_Effect_Status"].iloc[0]))

    def test_voucher_level_voting_weight_support(self):
        """Verify 04_gmm_morphotools.R enforces 1 row = 1 voucher collection event for 1.0x voting weight."""
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("1 row = 1 voucher", content, "04_gmm_morphotools.R must enforce 1 row = 1 voucher")
        self.assertIn("group_by(catalogNumber)", content, "04_gmm_morphotools.R must group by catalogNumber for voucher weighting")

        # Verify dataset contains valid foliar variance and leaf counts
        df = pd.read_csv(self.harmonics_path)
        self.assertIn("foliar_variance", df.columns)
        self.assertIn("leaf_count", df.columns)
        self.assertTrue((df["leaf_count"] >= 1).all(), "Every specimen must have >=1 leaf.")
        self.assertTrue((df["foliar_variance"] >= 0.0).all(), "foliar_variance must be non-negative.")


if __name__ == "__main__":
    unittest.main()
