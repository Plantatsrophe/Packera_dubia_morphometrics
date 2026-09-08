#!/usr/bin/env python3
"""
===============================================================================
Script: test_symmetric_fourier_validation.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit test suite validating Symmetric Fourier Decomposition and PERMANOVA
    empirical tier equivalence in scripts/morphometrics/03_fourier_extractor.R
    and its exported diagnostic reports.
===============================================================================
"""

import unittest
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]


class TestSymmetricFourierValidation(unittest.TestCase):
    """Test suite for Symmetric EFA and PERMANOVA empirical tier validation."""

    def setUp(self):
        self.r_script_path = PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R"
        self.harmonics_path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics.csv"
        self.report_path = PROJECT_ROOT / "outputs" / "reports" / "tier_symmetry_validation.csv"
        self.plot_path = PROJECT_ROOT / "outputs" / "figures" / "tier1_vs_tier2_density_overlay.pdf"

    def test_r_script_exists_and_within_line_limit(self):
        """Verify 03_fourier_extractor.R exists and adheres to modular line limits (<=800 lines)."""
        self.assertTrue(self.r_script_path.exists(), "03_fourier_extractor.R must exist.")
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        self.assertLessEqual(len(lines), 800, f"Script exceeds 800 lines: {len(lines)}")
        self.assertGreaterEqual(len(lines), 350, f"Script is unexpectedly short: {len(lines)}")

    def test_r_script_symmetric_decomposition_and_options(self):
        """Verify 03_fourier_extractor.R includes symmetric decomposition and defensive fallbacks."""
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Check CLI arguments
        self.assertIn("--report-out", content)
        self.assertIn("--plot-out", content)
        self.assertIn("--permutations", content)
        self.assertIn("--seed", content)

        # Check symmetric harmonic extraction (An and Dn)
        self.assertIn("sym_harmonics", content)
        self.assertIn("sym_names", content)

        # Check defensive PERMANOVA handling (vegan and manova fallback)
        self.assertIn("vegan::adonis2", content)
        self.assertIn("stats::manova", content)
        self.assertIn("reconstruction_tier", content)

        # Check random seed reproducibility
        self.assertIn("opts$seed", content)

    def test_harmonics_table_symmetric_pca_columns(self):
        """Verify leaf_efa_harmonics.csv contains symmetric PC scores and metadata."""
        self.assertTrue(self.harmonics_path.exists(), "leaf_efa_harmonics.csv must exist.")
        df = pd.read_csv(self.harmonics_path)

        # Verify key metadata columns
        required_cols = [
            "catalogNumber", "assigned_tier", "reconstruction_tier",
            "scientificName", "determiner_tier", "PC1", "PC2", "PC3", "PC4", "PC5"
        ]
        for col in required_cols:
            self.assertIn(col, df.columns, f"Missing column {col} in harmonics table.")

        # Verify closed outlines have valid PC scores
        closed = df[df["reconstruction_tier"].isin(["Tier 1", "Tier 2"])]
        self.assertGreater(len(closed), 1000, "Expected >1000 closed leaf outlines.")
        for p in range(1, 6):
            self.assertFalse(closed[f"PC{p}"].isna().any(), f"Found NaNs in PC{p} for closed outlines.")

        # Verify symmetric harmonic columns are present
        for i in range(1, 13):
            self.assertIn(f"A{i}", df.columns)
            self.assertIn(f"D{i}", df.columns)

    def test_permanova_validation_report(self):
        """Verify PERMANOVA tier validation report schema and statistical properties."""
        self.assertTrue(self.report_path.exists(), "tier_symmetry_validation.csv must exist.")
        df = pd.read_csv(self.report_path)

        # Check column schema
        expected_cols = ["Term", "Df", "SumOfSqs", "R2", "F", "p_value"]
        for col in expected_cols:
            self.assertIn(col, df.columns, f"Missing column {col} in validation report.")

        terms = df["Term"].tolist()
        self.assertIn("scientificName", terms)
        self.assertIn("determiner_tier", terms)
        self.assertIn("reconstruction_tier", terms)
        self.assertIn("Residual", terms)
        self.assertIn("Total", terms)

        # Check scientificName accounts for primary significant variance
        sp_row = df[df["Term"] == "scientificName"].iloc[0]
        self.assertLess(sp_row["p_value"], 0.05, "scientificName should be statistically significant.")
        self.assertGreater(sp_row["SumOfSqs"], 1.0, "scientificName SumOfSqs should be substantial.")

        # Check reconstruction_tier variance explained
        tier_row = df[df["Term"] == "reconstruction_tier"].iloc[0]
        self.assertLessEqual(tier_row["R2"], 0.02, "reconstruction_tier should account for <=2% of variance.")

    def test_density_plot_artifact(self):
        """Verify comparative density plot exists and is a valid non-empty PDF."""
        self.assertTrue(self.plot_path.exists(), "tier1_vs_tier2_density_overlay.pdf must exist.")
        self.assertGreater(self.plot_path.stat().st_size, 1000, "Density plot PDF must be > 1KB.")
        with open(self.plot_path, "rb") as f:
            header = f.read(5)
        self.assertEqual(header, b"%PDF-", "Density plot file must have valid PDF magic bytes.")


if __name__ == "__main__":
    unittest.main()
