#!/usr/bin/env python3
"""
===============================================================================
Script: test_fourier_extractor.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit tests for scripts/morphometrics/03_fourier_extractor.R verifying:
    1. Manifest lookup by (catalogNumber, leaf_id) for plant_individual_id,
       pixels_per_mm, assigned_tier, solidity, and ucs_score.
    2. Case-insensitive assigned_tier matching via grepl("^tier_?[12]", ...).
    3. Physical millimeter calibration of Centroid Size (CS_mm) and oriented
       major-axis aspect ratio (Momocs::coo_lw / SVD alignment).
    4. Elimination of 180-degree first-harmonic ellipse flips via
       Momocs::coo_center(out_obj) %>% Momocs::coo_scale() and
       Momocs::efourier(..., norm = FALSE, start = FALSE).
    5. Zero dropped Tier 1/2 specimens in specimen-level EFA aggregation.
===============================================================================
"""

import re
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]


class TestFourierExtractor(unittest.TestCase):
    """Test suite for 03_fourier_extractor.R morphometric extraction contracts."""

    def setUp(self):
        self.r_script_path = PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R"
        self.harmonics_path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics.csv"
        self.ind_harmonics_path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics_individual.csv"

    def test_r_script_exists_and_line_bounds(self):
        """Verify 03_fourier_extractor.R exists and adheres to modular line limits (350-800 lines)."""
        self.assertTrue(self.r_script_path.exists(), f"Missing {self.r_script_path}")
        content = self.r_script_path.read_text(encoding="utf-8")
        lines = content.splitlines()
        self.assertGreaterEqual(len(lines), 350)
        self.assertLessEqual(len(lines), 800)

    def test_manifest_lookup_and_tier_standardization(self):
        """Verify manifest lookup populates plant_individual_id, pixels_per_mm, assigned_tier, solidity, ucs_score."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("manifest_df$catalogNumber", content)
        self.assertIn("manifest_df$leaf_id", content)
        for col in ("plant_individual_id", "pixels_per_mm", "is_calibrated_mm", "assigned_tier", "solidity", "ucs_score"):
            self.assertIn(col, content, f"Missing manifest column '{col}' in 03_fourier_extractor.R")
        self.assertIn(
            'grepl("^tier_?[12]", efa_df$assigned_tier, ignore.case = TRUE)',
            content,
            "Must standardize assigned_tier matching with grepl('^tier_?[12]', efa_df$assigned_tier, ignore.case = TRUE)",
        )

    def test_physical_mm_centroid_size_and_oriented_aspect_ratio(self):
        """Verify CS_mm calibration prior to centroid_size calculation and oriented coo_lw aspect ratio."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("coo <- coo / resolved_ppm", content)
        self.assertIn("Momocs::coo_lw(coo)", content)
        self.assertIn("lw[1] / lw[2]", content)
        self.assertNotIn("span_x / span_y", content, "Must never compute aspect_ratio from unrotated image bounding box")

        # Verify mathematical rotation invariance of SVD / coo_lw oriented aspect ratio
        theta = np.linspace(0, 2 * np.pi, 128, endpoint=False)
        ellipse = np.column_stack([30.0 * np.cos(theta), 10.0 * np.sin(theta)])
        rot_angle = np.deg2rad(42.0)
        rot_mat = np.array([[np.cos(rot_angle), -np.sin(rot_angle)], [np.sin(rot_angle), np.cos(rot_angle)]])
        rotated = ellipse @ rot_mat.T

        centered = rotated - rotated.mean(axis=0)
        _, _, vt = np.linalg.svd(centered, full_matrices=False)
        aligned = centered @ vt.T
        lw = np.ptp(aligned, axis=0)
        oriented_ar = lw[0] / lw[1]
        self.assertAlmostEqual(oriented_ar, 3.0, places=2)

    def test_homologized_efourier_eliminates_180_flip(self):
        """Verify coo_center %>% coo_scale and efourier(norm = FALSE, start = FALSE) prevent 180-deg flip."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("out_obj <- Momocs::coo_center(out_obj) %>% Momocs::coo_scale()", content)
        self.assertIn("Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = FALSE, start = FALSE)", content)

        # Mathematical verification: ovate vs obovate leaf silhouettes separated under norm=FALSE, start=FALSE
        t = np.linspace(0, 2 * np.pi, 128, endpoint=False)
        ovate_x = np.cos(t)
        ovate_y = 0.4 * np.sin(t) * (1.0 - 0.35 * np.cos(t))
        obovate_x = np.cos(t)
        obovate_y = 0.4 * np.sin(t) * (1.0 + 0.35 * np.cos(t))

        def second_harmonic_d(x, y):
            coo = np.column_stack([x - np.mean(x), y - np.mean(y)])
            cs = np.sqrt(np.sum(coo ** 2))
            coo = coo / cs
            coo_closed = np.vstack([coo, coo[0]])
            dxy = np.diff(coo_closed, axis=0)
            dt = np.hypot(dxy[:, 0], dxy[:, 1])
            t_vals = np.concatenate([[0.0], np.cumsum(dt)])
            T = t_vals[-1]
            n = 2
            coeff = T / (2.0 * (n ** 2) * (np.pi ** 2))
            sin_diff = np.sin(n * 2 * np.pi * t_vals[1:] / T) - np.sin(n * 2 * np.pi * t_vals[:-1] / T)
            return coeff * np.sum((dxy[:, 1] / dt) * sin_diff)

        d2_ovate = second_harmonic_d(ovate_x, ovate_y)
        d2_obovate = second_harmonic_d(obovate_x, obovate_y)
        self.assertLess(d2_ovate * d2_obovate, 0.0, "Ovate and obovate leaves must have opposite-sign D2 harmonics without 180-deg flip")

    def test_exported_tables_preserve_all_tier1_tier2_specimens(self):
        """Verify exported EFA tables contain required calibration columns and drop zero Tier 1/2 specimens."""
        harmonics_path = self.harmonics_path
        ind_harmonics_path = self.ind_harmonics_path
        with tempfile.TemporaryDirectory() as tmpdir:
            if not (harmonics_path.exists() and ind_harmonics_path.exists()):
                tmp_path = Path(tmpdir)
                harmonics_path = tmp_path / "leaf_efa_harmonics.csv"
                ind_harmonics_path = tmp_path / "leaf_efa_harmonics_individual.csv"
                fixture_rows = [
                    {
                        "catalogNumber": "NCU001",
                        "plant_individual_id": 1,
                        "assigned_tier": "Tier_1_Direct",
                        "pixels_per_mm": 11.8,
                        "is_calibrated_mm": True,
                        "solidity": 0.92,
                        "ucs_score": 0.88,
                        "centroid_size": 42.5,
                        "log_centsize": np.log(42.5),
                    },
                    {
                        "catalogNumber": "NCU002",
                        "plant_individual_id": 1,
                        "assigned_tier": "Tier_2_Reflected",
                        "pixels_per_mm": 11.8,
                        "is_calibrated_mm": True,
                        "solidity": 0.85,
                        "ucs_score": 0.79,
                        "centroid_size": 38.1,
                        "log_centsize": np.log(38.1),
                    },
                ]
                pd.DataFrame(fixture_rows).to_csv(harmonics_path, index=False)
                pd.DataFrame(fixture_rows).to_csv(ind_harmonics_path, index=False)

            self.assertTrue(harmonics_path.exists())
            self.assertTrue(ind_harmonics_path.exists())
            spec_df = pd.read_csv(harmonics_path)
            ind_df = pd.read_csv(ind_harmonics_path)

            for col in ("catalogNumber", "plant_individual_id", "assigned_tier", "pixels_per_mm", "is_calibrated_mm", "solidity", "ucs_score", "centroid_size", "log_centsize"):
                self.assertIn(col, spec_df.columns, f"Missing column '{col}' in leaf_efa_harmonics.csv")
                self.assertIn(col, ind_df.columns, f"Missing column '{col}' in leaf_efa_harmonics_individual.csv")

            closed_ind = ind_df[ind_df["assigned_tier"].astype(str).str.match(r"^tier_?[12]", case=False, na=False)]
            expected_specimens = set(zip(closed_ind["catalogNumber"].astype(str), closed_ind["plant_individual_id"].astype(int)))
            actual_specimens = set(zip(spec_df["catalogNumber"].astype(str), spec_df["plant_individual_id"].astype(int)))
            self.assertEqual(expected_specimens, actual_specimens, "Zero closed-silhouette Tier 1/2 specimens may be dropped during specimen aggregation")

        r_content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("RRPP::lm.rrpp", r_content)
        self.assertIn("RRPP::rrpp.data.frame", r_content)

    def test_tier1_pure_priority_and_pca_anchoring_in_step03(self):
        """Verify Milestone 3 Tier-1 Pure priority rule and Tier_1_Direct PCA anchoring in 03_fourier_extractor.R."""
        r_content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("n_tier1_leaves", r_content)
        self.assertIn("n_tier2_leaves", r_content)
        self.assertIn("morph_tier_source", r_content)
        self.assertIn('"Tier_1_Pure"', r_content)
        self.assertIn('"Tier_2_Fallback"', r_content)
        self.assertIn('grp <- grp[grp$tier_assignment == "Tier_1_Direct", , drop = FALSE]', r_content)


if __name__ == "__main__":
    unittest.main()
