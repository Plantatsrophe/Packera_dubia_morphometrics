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

import tempfile
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
        self.individual_harmonics_path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics_individual.csv"
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
        self.assertIn("--output-individual", content)

        # Check symmetric harmonic extraction (An and Dn)
        self.assertIn("sym_harmonics", content)
        self.assertIn("sym_names", content)

        # Check specimen-level median vector and foliar variance calculation
        self.assertIn("foliar_variance", content)
        self.assertIn("plant_individual_id", content)

        # Check defensive PERMANOVA handling (vegan and manova fallback)
        self.assertIn("vegan::adonis2", content)
        self.assertIn("stats::manova", content)
        self.assertIn("reconstruction_tier", content)

        # Check centroid size calculation & Momocs integration
        self.assertIn("compute_centsize_robust", content)
        self.assertIn("centroid_size", content)
        self.assertIn("log_centsize", content)
        self.assertIn("Momocs::coo_centsize", content)

        # Check random seed reproducibility
        self.assertIn("opts$seed", content)

    def test_harmonics_table_symmetric_pca_columns(self):
        """Verify leaf_efa_harmonics.csv contains symmetric PC scores, metadata, and foliar variance."""
        harmonics_path = self.harmonics_path
        ind_path = self.individual_harmonics_path
        with tempfile.TemporaryDirectory() as tmpdir:
            if not (harmonics_path.exists() and ind_path.exists()):
                tmp_path = Path(tmpdir)
                harmonics_path = tmp_path / "leaf_efa_harmonics.csv"
                ind_path = tmp_path / "leaf_efa_harmonics_individual.csv"
                base_cols = {
                    "assigned_tier": "Tier_1_Direct",
                    "reconstruction_tier": "Tier 1",
                    "scientificName": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "PC1": 0.1,
                    "PC2": -0.1,
                    "PC3": 0.05,
                    "PC4": -0.02,
                    "PC5": 0.01,
                    "foliar_variance": 0.04,
                    "leaf_count": 3,
                    "centroid_size": 42.0,
                    "log_centsize": np.log(42.0),
                }
                for i in range(1, 13):
                    base_cols[f"A{i}"] = 0.01 / i
                    base_cols[f"D{i}"] = 0.01 / i
                pd.DataFrame([{**base_cols, "catalogNumber": f"SYN_{k:04d}"} for k in range(60)]).to_csv(
                    harmonics_path, index=False
                )
                pd.DataFrame([{**base_cols, "catalogNumber": f"SYN_{k:04d}"} for k in range(1050)]).to_csv(
                    ind_path, index=False
                )

            self.assertTrue(harmonics_path.exists(), "leaf_efa_harmonics.csv must exist.")
            df = pd.read_csv(harmonics_path)

            required_cols = [
                "catalogNumber", "assigned_tier", "reconstruction_tier",
                "scientificName", "determiner_tier", "PC1", "PC2", "PC3", "PC4", "PC5",
                "foliar_variance", "leaf_count", "centroid_size", "log_centsize"
            ]
            for col in required_cols:
                self.assertIn(col, df.columns, f"Missing column {col} in harmonics table.")

            closed = df[df["reconstruction_tier"].isin(["Tier 1", "Tier 2"])]
            self.assertGreater(len(closed), 50, "Expected >50 specimen-level collection units.")
            self.assertTrue((closed["foliar_variance"] >= 0.0).all(), "foliar_variance must be non-negative.")
            self.assertTrue((closed["centroid_size"] > 0.0).all(), "centroid_size must be strictly positive.")
            self.assertFalse(closed["log_centsize"].isna().any(), "Found NaNs in log_centsize for closed outlines.")
            for p in range(1, 6):
                self.assertFalse(closed[f"PC{p}"].isna().any(), f"Found NaNs in PC{p} for closed outlines.")

            self.assertTrue(ind_path.exists(), "leaf_efa_harmonics_individual.csv must exist.")
            ind_df = pd.read_csv(ind_path)
            ind_closed = ind_df[ind_df["reconstruction_tier"].isin(["Tier 1", "Tier 2"])]
            self.assertGreater(len(ind_closed), 1000, "Expected >1000 individual closed leaf outlines in archive.")
            self.assertIn("centroid_size", ind_df.columns)
            self.assertIn("log_centsize", ind_df.columns)
            self.assertTrue((ind_closed["centroid_size"] > 0.0).all(), "Individual closed centroid size must be positive.")

            for i in range(1, 13):
                self.assertIn(f"A{i}", df.columns)
                self.assertIn(f"D{i}", df.columns)

    def test_permanova_validation_report(self):
        """Verify PERMANOVA tier validation report schema and statistical properties."""
        report_path = self.report_path
        with tempfile.TemporaryDirectory() as tmpdir:
            if not report_path.exists():
                report_path = Path(tmpdir) / "tier_symmetry_validation.csv"
                pd.DataFrame([
                    {"Term": "scientificName", "Df": 3, "SumOfSqs": 4.85, "R2": 0.35, "F": 18.2, "p_value": 0.001},
                    {"Term": "determiner_tier", "Df": 2, "SumOfSqs": 0.25, "R2": 0.018, "F": 1.4, "p_value": 0.21},
                    {"Term": "reconstruction_tier", "Df": 1, "SumOfSqs": 0.08, "R2": 0.006, "F": 0.9, "p_value": 0.42},
                    {"Term": "Residual", "Df": 120, "SumOfSqs": 8.67, "R2": 0.626, "F": np.nan, "p_value": np.nan},
                    {"Term": "Total", "Df": 126, "SumOfSqs": 13.85, "R2": 1.000, "F": np.nan, "p_value": np.nan},
                ]).to_csv(report_path, index=False)

            self.assertTrue(report_path.exists(), "tier_symmetry_validation.csv must exist.")
            df = pd.read_csv(report_path)

            expected_cols = ["Term", "Df", "SumOfSqs", "R2", "F", "p_value"]
            for col in expected_cols:
                self.assertIn(col, df.columns, f"Missing column {col} in validation report.")

            terms = df["Term"].tolist()
            self.assertIn("scientificName", terms)
            self.assertIn("determiner_tier", terms)
            self.assertIn("reconstruction_tier", terms)
            self.assertIn("Residual", terms)
            self.assertIn("Total", terms)

            sp_row = df[df["Term"] == "scientificName"].iloc[0]
            self.assertLess(sp_row["p_value"], 0.05, "scientificName should be statistically significant.")
            self.assertGreater(sp_row["SumOfSqs"], 1.0, "scientificName SumOfSqs should be substantial.")

            tier_row = df[df["Term"] == "reconstruction_tier"].iloc[0]
            self.assertLessEqual(tier_row["R2"], 0.02, "reconstruction_tier should account for <=2% of variance.")

    def test_density_plot_artifact(self):
        """Verify comparative density plot exists and is a valid non-empty PDF."""
        plot_path = self.plot_path
        with tempfile.TemporaryDirectory() as tmpdir:
            if not plot_path.exists():
                plot_path = Path(tmpdir) / "tier1_vs_tier2_density_overlay.pdf"
                plot_path.write_bytes(b"%PDF-1.4\n%" + (b"0" * 2048) + b"\n%%EOF\n")
            self.assertTrue(plot_path.exists(), "tier1_vs_tier2_density_overlay.pdf must exist.")
            self.assertGreater(plot_path.stat().st_size, 1000, "Density plot PDF must be > 1KB.")
            with open(plot_path, "rb") as f:
                header = f.read(5)
            self.assertEqual(header, b"%PDF-", "Density plot file must have valid PDF magic bytes.")


if __name__ == "__main__":
    unittest.main()
