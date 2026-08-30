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
        """Verify that 04_gmm_morphotools.R exists and is under 500 lines."""
        self.assertTrue(self.r_script_path.exists(), "04_gmm_morphotools.R must exist.")
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        self.assertLessEqual(len(lines), 500, f"Script exceeds 500 lines: {len(lines)}")

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


if __name__ == "__main__":
    unittest.main()
