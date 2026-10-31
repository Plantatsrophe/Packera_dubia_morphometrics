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
    5. Strict Tier_1_Direct EFA filtering, unweighted OTU median aggregation,
       adaptive anchor promotion for 0-Gold/Silver taxa, and Step 01/02 metadata handoff.
===============================================================================
"""

import shutil
import subprocess
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
        self.assertLess(
            d2_ovate * d2_obovate,
            0.0,
            "Ovate and obovate leaves must have opposite-sign D2 harmonics without 180-deg flip",
        )

    def test_tier1_pure_priority_and_pca_anchoring_in_step03(self):
        """Verify Milestone 3 Tier-1 Pure priority rule and Tier_1_Direct PCA anchoring in 03_fourier_extractor.R."""
        r_content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("n_tier1_leaves", r_content)
        self.assertIn("n_tier2_leaves", r_content)
        self.assertIn("morph_tier_source", r_content)
        self.assertIn('"Tier_1_Pure"', r_content)
        self.assertIn('"Tier_2_Fallback"', r_content)
        self.assertIn('grp <- grp[grp$tier_assignment == "Tier_1_Direct", , drop = FALSE]', r_content)

    def test_empirical_efa_alignment_harmonic_power_and_otu_aggregation(self):
        """Verify Phase 3 empirical EFA protocol: coo_alignxax, calibrate_harmonicpower_efourier, RRPP, and OTU-level PCA."""
        r_content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("align_contour_basal_apical", r_content)
        self.assertIn("Momocs::coo_alignxax", r_content)
        self.assertIn("circularity > 0.75", r_content)
        self.assertIn("verify_harmonic_power_calibration", r_content)
        self.assertIn("calibrate_harmonicpower_efourier", r_content)
        self.assertIn("heteroblasty_index", r_content)
        self.assertIn("otu_t1_idx", r_content)
        self.assertIn("symmetrize_efa_harmonics", r_content)
        self.assertIn("compute_authority_precision_weights", r_content)
        self.assertIn("RRPP::lm.rrpp", r_content)
        self.assertIn("RRPP::rrpp.data.frame", r_content)

    def test_tier1_only_efa_unweighted_median_and_tier3_passive_pca_invariance(self):
        """Verify only Tier 1 intact contours enter EFA, OTU aggregation uses unweighted median, and Tier_3_Bronze has zero influence on reference PCA."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        def write_contour_csv(path: Path, cat: str, lid: int, a: float, b: float, d2_mod: float = 0.0):
            x = a * np.cos(t)
            y = b * np.sin(t) * (1.0 + d2_mod * np.cos(t))
            pd.DataFrame({
                "catalogNumber": cat,
                "leaf_id": lid,
                "point_index": np.arange(120),
                "x": x + 200.0,
                "y": y + 200.0,
            }).to_csv(path, index=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = [
                {"catalogNumber": "G1", "species_raw": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "G2", "species_raw": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "G3", "species_raw": "Packera anonyma", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "S1", "species_raw": "Packera anonyma", "determiner_tier": "Tier_2_Silver"},
                {"catalogNumber": "S2", "species_raw": "Packera paupercula", "determiner_tier": "Tier_2_Silver"},
                {"catalogNumber": "B1", "species_raw": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
                {"catalogNumber": "T2_ONLY", "species_raw": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            write_contour_csv(contours_dir / "G1_leaf1.csv", "G1", 1, 90.0, 35.0, -0.15)
            write_contour_csv(contours_dir / "G1_leaf2.csv", "G1", 2, 90.0, 45.0, 0.00)
            write_contour_csv(contours_dir / "G1_leaf3.csv", "G1", 3, 90.0, 60.0, 0.20)
            write_contour_csv(contours_dir / "G1_leaf4.csv", "G1", 4, 90.0, 85.0, 0.45)

            write_contour_csv(contours_dir / "G2_leaf1.csv", "G2", 1, 85.0, 42.0, -0.10)
            write_contour_csv(contours_dir / "G3_leaf1.csv", "G3", 1, 95.0, 30.0, 0.18)
            write_contour_csv(contours_dir / "S1_leaf1.csv", "S1", 1, 92.0, 32.0, 0.12)
            write_contour_csv(contours_dir / "S2_leaf1.csv", "S2", 1, 100.0, 28.0, -0.22)
            write_contour_csv(contours_dir / "B1_leaf1.csv", "B1", 1, 88.0, 40.0, 0.05)
            write_contour_csv(contours_dir / "T2_ONLY_leaf1.csv", "T2_ONLY", 1, 80.0, 50.0, 0.10)
            write_contour_csv(contours_dir / "EXCL_PLAT_leaf1.csv", "EXCL_PLAT", 1, 85.0, 45.0, 0.08)

            manifest_rows = [
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.95, "ucs_score": 0.95, "contour_path": str(contours_dir / "G1_leaf1.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 2, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.90, "ucs_score": 0.90, "contour_path": str(contours_dir / "G1_leaf2.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 3, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.85, "ucs_score": 0.85, "contour_path": str(contours_dir / "G1_leaf3.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 4, "assigned_tier": "Tier_2_Reflected", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.92, "ucs_score": 0.92, "contour_path": str(contours_dir / "G1_leaf4.csv")},
                {"catalogNumber": "G2", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.93, "ucs_score": 0.93, "contour_path": str(contours_dir / "G2_leaf1.csv")},
                {"catalogNumber": "G3", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.91, "ucs_score": 0.91, "contour_path": str(contours_dir / "G3_leaf1.csv")},
                {"catalogNumber": "S1", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.92, "ucs_score": 0.92, "contour_path": str(contours_dir / "S1_leaf1.csv")},
                {"catalogNumber": "S2", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.90, "ucs_score": 0.90, "contour_path": str(contours_dir / "S2_leaf1.csv")},
                {"catalogNumber": "B1", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.91, "ucs_score": 0.91, "contour_path": str(contours_dir / "B1_leaf1.csv")},
                {"catalogNumber": "T2_ONLY", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_2_Reflected", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.89, "ucs_score": 0.89, "contour_path": str(contours_dir / "T2_ONLY_leaf1.csv")},
                {"catalogNumber": "EXCL_PLAT", "plant_individual_id": 0, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.94, "ucs_score": 0.94, "contour_path": str(contours_dir / "EXCL_PLAT_leaf1.csv")},
            ]
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            def run_step03(tag: str):
                out_spec = tmp_p / f"leaf_efa_harmonics_{tag}.csv"
                out_ind = tmp_p / f"leaf_efa_harmonics_ind_{tag}.csv"
                out_rep = tmp_p / f"tier_report_{tag}.csv"
                out_fig = tmp_p / f"fig_{tag}.pdf"
                proc = subprocess.run(
                    [
                        "Rscript",
                        str(self.r_script_path),
                        "--contours-dir", str(contours_dir),
                        "--masks-dir", str(tmp_p / "masks"),
                        "--vouchers", str(vouchers_csv),
                        "--manifest", str(manifest_csv),
                        "--output", str(out_spec),
                        "--output-individual", str(out_ind),
                        "--report-out", str(out_rep),
                        "--plot-out", str(out_fig),
                        "--permutations", "19",
                    ],
                    cwd=str(PROJECT_ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, f"Step 03 Rscript failed ({tag}):\n{proc.stderr}\n{proc.stdout}")
                return pd.read_csv(out_spec), pd.read_csv(out_ind)

            spec1, ind1 = run_step03("run1")

            # 1. Only intact Tier_1_Direct contours present in vouchers_df enter leaf_efa_harmonics.csv
            self.assertTrue((ind1["assigned_tier"] == "Tier_1_Direct").all(), "Only Tier_1_Direct leaves may enter individual EFA table.")
            self.assertTrue((spec1["assigned_tier"] == "Tier_1_Direct").all(), "Only Tier_1_Direct leaves may enter specimen EFA table.")
            self.assertNotIn("T2_ONLY", set(spec1["catalogNumber"]), "Tier_2_Reflected-only specimen must not enter EFA table.")
            self.assertNotIn("EXCL_PLAT", set(ind1["catalogNumber"]), "Excluded taxon absent from vouchers_df must not be ingested into individual EFA table.")
            self.assertNotIn("EXCL_PLAT", set(spec1["catalogNumber"]), "Excluded taxon absent from vouchers_df must not enter specimen EFA table.")
            pub_utils_content = (PROJECT_ROOT / "scripts" / "morphometrics" / "pub_figures_utils.R").read_text(encoding="utf-8")
            self.assertNotIn("if (length(disp_tx) < 4) disp_tx <- unique(c(disp_tx, focal_4))[1:4]", pub_utils_content)
            self.assertEqual(int(spec1.loc[spec1["catalogNumber"] == "G1", "n_tier2_leaves"].iloc[0]), 0)

            for col in ("catalogNumber", "plant_individual_id", "assigned_tier", "pixels_per_mm", "is_calibrated_mm", "solidity", "ucs_score", "centroid_size", "log_centsize"):
                self.assertIn(col, spec1.columns)
                self.assertIn(col, ind1.columns)

            # 2. Unweighted coordinate-wise median per harmonic coefficient for G1
            g1_ind = ind1[ind1["catalogNumber"] == "G1"]
            g1_spec = spec1[spec1["catalogNumber"] == "G1"].iloc[0]
            for h_col in ["A1", "A2", "A3", "D1", "D2", "D3"]:
                expected_med = float(np.median(g1_ind[h_col].to_numpy()))
                self.assertAlmostEqual(float(g1_spec[h_col]), expected_med, places=5)

            # 3. Tier_3_Bronze has authority_weight == 0.0 and zero influence on reference PCA
            b1_row = spec1[spec1["catalogNumber"] == "B1"].iloc[0]
            self.assertEqual(float(b1_row["authority_weight"]), 0.0)

            # Perturb B1's contour drastically and re-run Step 03
            write_contour_csv(contours_dir / "B1_leaf1.csv", "B1", 1, 130.0, 22.0, -0.38)
            spec2, ind2 = run_step03("run2")

            ref1 = spec1[spec1["catalogNumber"] != "B1"].sort_values("catalogNumber")
            ref2 = spec2[spec2["catalogNumber"] != "B1"].sort_values("catalogNumber")
            for pc_col in ["PC1", "PC2", "PC3", "Sym_PC1", "Sym_PC2"]:
                np.testing.assert_allclose(
                    ref1[pc_col].to_numpy(),
                    ref2[pc_col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter reference {pc_col} scores.",
                )
            ref_ind1 = ind1[ind1["catalogNumber"] != "B1"].sort_values(["catalogNumber", "leaf_id"])
            ref_ind2 = ind2[ind2["catalogNumber"] != "B1"].sort_values(["catalogNumber", "leaf_id"])
            for pc_col in ["PC1", "PC2", "PC3", "Sym_PC1", "Sym_PC2"]:
                np.testing.assert_allclose(
                    ref_ind1[pc_col].to_numpy(),
                    ref_ind2[pc_col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter individual reference {pc_col} scores.",
                )

    def test_zero_gold_silver_adaptive_anchor_promotion_and_step01_step02_handoff(self):
        """Verify Step 03 promotes clean Bronze anchors for 0-Gold/Silver species and preserves Step 01/02 metadata and petiole/dissection traits."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        def write_contour_csv(path: Path, cat: str, lid: int, a: float, b: float, d2_mod: float = 0.0):
            x = a * np.cos(t)
            y = b * np.sin(t) * (1.0 + d2_mod * np.cos(t))
            pd.DataFrame({
                "catalogNumber": cat,
                "leaf_id": lid,
                "point_index": np.arange(120),
                "x": x + 200.0,
                "y": y + 200.0,
            }).to_csv(path, index=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            voucher_defs = [
                ("D_G1", "Packera dubia", "Tier_1_Gold", "Holotype"),
                ("D_G2", "Packera dubia", "Tier_1_Gold", ""),
                ("D_B1", "Packera dubia", "Tier_3_Bronze", ""),
                ("A_S1", "Packera anonyma", "Tier_2_Silver", ""),
                ("A_S2", "Packera anonyma", "Tier_2_Silver", ""),
                # Rare endemic Packera serpenticola with ZERO Gold and ZERO Silver vouchers (only Bronze)
                ("R_B1", "Packera serpenticola", "Tier_3_Bronze", ""),
                ("R_B2", "Packera serpenticola", "Tier_3_Bronze", ""),
                ("R_B3", "Packera serpenticola", "Tier_3_Bronze", ""),
            ]
            vouchers = [
                {
                    "catalogNumber": cat,
                    "taxon": sp,
                    "scientificName": sp,
                    "species_raw": sp,
                    "institutionCode": "NCU",
                    "determiner_tier": tier,
                    "authority_tier": tier,
                    "type_status": t_stat,
                    "stateProvince": "North Carolina",
                    "county": "Macon",
                    "locality": "Whiteside Mountain granitic dome",
                    "decimalLatitude": 35.081,
                    "decimalLongitude": -83.138,
                    "coordinateUncertainty": 120.0,
                    "year": 2021,
                    "month": 5,
                    "day": 14,
                    "eventDate": "2021-05-14",
                }
                for cat, sp, tier, t_stat in voucher_defs
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            specs = [
                ("D_G1", 92.0, 38.0, -0.14, False),
                ("D_G2", 90.0, 40.0, -0.12, False),
                ("D_B1", 89.0, 41.0, -0.10, False),
                ("A_S1", 96.0, 28.0, 0.16, True),
                ("A_S2", 94.0, 30.0, 0.14, True),
                ("R_B1", 72.0, 58.0, 0.28, False),
                ("R_B2", 74.0, 56.0, 0.26, False),
                ("R_B3", 73.0, 57.0, 0.27, False),
            ]
            manifest_rows = []
            for idx, (cat, a_v, b_v, d2_v, is_diss) in enumerate(specs):
                cpath = contours_dir / f"{cat}_leaf1.csv"
                write_contour_csv(cpath, cat, 1, a_v, b_v, d2_v)
                manifest_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.42 if is_diss else 0.93,
                    "ucs_score": 0.94,
                    "petiole_length_mm": 18.5 + idx,
                    "blade_length_mm": 24.0 + idx,
                    "petiole_to_blade_ratio": round((18.5 + idx) / (24.0 + idx), 4),
                    "is_dissected": is_diss,
                    "contour_path": str(cpath),
                })
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            out_spec = tmp_p / "leaf_efa_harmonics.csv"
            out_ind = tmp_p / "leaf_efa_harmonics_ind.csv"
            proc = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--contours-dir", str(contours_dir),
                    "--masks-dir", str(tmp_p / "masks"),
                    "--vouchers", str(vouchers_csv),
                    "--manifest", str(manifest_csv),
                    "--output", str(out_spec),
                    "--output-individual", str(out_ind),
                    "--report-out", str(tmp_p / "tier_rep.csv"),
                    "--plot-out", str(tmp_p / "tier_fig.pdf"),
                    "--permutations", "19",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Step 03 failed on 0-Gold/Silver & handoff run:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            spec_df = pd.read_csv(out_spec)
            self.assertEqual(len(spec_df), 8)

            # 1. Adaptive Anchor Promotion & Empirical Bayes active weight ordering: Gold > Silver > Promoted Bronze > 0
            self.assertEqual(float(spec_df.loc[spec_df["catalogNumber"] == "D_B1", "authority_weight"].iloc[0]), 0.0)
            rare_rows = spec_df[spec_df["catalogNumber"].isin(["R_B1", "R_B2", "R_B3"])]
            self.assertGreaterEqual(int((rare_rows["authority_weight"] > 0.0).sum()), 2)

            w_gold = float(spec_df.loc[spec_df["catalogNumber"] == "D_G1", "authority_weight"].iloc[0])
            w_silver = float(spec_df.loc[spec_df["catalogNumber"] == "A_S1", "authority_weight"].iloc[0])
            w_promo = float(rare_rows.loc[rare_rows["authority_weight"] > 0.0, "authority_weight"].iloc[0])
            self.assertGreater(w_gold, w_silver)
            self.assertGreater(w_silver, w_promo)
            self.assertGreater(w_promo, 0.0)
            active_w = spec_df.loc[spec_df["authority_weight"] > 0.0, "authority_weight"]
            self.assertAlmostEqual(float(active_w.mean()), 1.0, places=4)

            # 2. Step 01 locality/date metadata & Step 02 petiole/dissection trait handoff
            for df_chk in (ind_df, spec_df):
                for req_col in (
                    "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "is_dissected",
                    "locality", "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
                    "year", "month", "day", "eventDate", "type_status",
                ):
                    self.assertIn(req_col, df_chk.columns)
            row0 = spec_df[spec_df["catalogNumber"] == "D_G1"].iloc[0]
            self.assertAlmostEqual(float(row0["petiole_length_mm"]), 18.5, places=3)
            self.assertAlmostEqual(float(row0["blade_length_mm"]), 24.0, places=3)
            self.assertEqual(str(row0["locality"]), "Whiteside Mountain granitic dome")
            self.assertEqual(str(row0["eventDate"]), "2021-05-14")
            row_diss = spec_df[spec_df["catalogNumber"] == "A_S1"].iloc[0]
            self.assertTrue(bool(row_diss["is_dissected"]))

    def test_dynamic_harmonic_calibration_and_dynamic_pc_export(self):
        """Verify Task 1B dynamic H_opt in [6, 24] harmonic calibration and dynamic PC1..PC{k_export} / Sym_PC1..Sym_PC{k_export} export (>5 PCs when k_95 > 5)."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("Momocs::efourier(out_obj, nb.h = h_opt, norm = FALSE, start = FALSE)", content)
        self.assertIn("efd_harmonic_opt", content)
        self.assertIn("efd_k95_pcs", content)

        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)
        rng = np.random.default_rng(20260430)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = []
            manifest_rows = []
            # Create 18 outlines with 8 orthogonal symmetric harmonic variation modes so k_95 >= 6
            freqs = [2, 3, 4, 5, 6, 7, 8, 9]
            for i in range(18):
                cat = f"DYN_H_{i:02d}"
                sp = "Packera millefolium" if i < 9 else "Packera schweinitziana"
                coefs = rng.normal(0.0, 0.055, size=len(freqs))
                x = 85.0 * np.cos(t) + sum(85.0 * coefs[m] * np.cos(freqs[m] * t) for m in range(len(freqs)))
                y = 42.0 * np.sin(t) + sum(42.0 * coefs[m] * np.sin(freqs[m] * t) for m in range(len(freqs)))
                cpath = contours_dir / f"{cat}_leaf1.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": 1,
                    "point_index": np.arange(120),
                    "x": x + 250.0,
                    "y": y + 250.0,
                }).to_csv(cpath, index=False)
                vouchers.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": "Tier_1_Gold",
                })
                manifest_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.62,
                    "ucs_score": 0.94,
                    "is_dissected": True,
                    "contour_path": str(cpath),
                })

            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)
            pd.DataFrame(manifest_rows).to_csv(manifest_csv, index=False)

            out_spec = tmp_p / "leaf_efa_harmonics.csv"
            out_ind = tmp_p / "leaf_efa_harmonics_ind.csv"
            proc = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--contours-dir", str(contours_dir),
                    "--masks-dir", str(tmp_p / "masks"),
                    "--vouchers", str(vouchers_csv),
                    "--manifest", str(manifest_csv),
                    "--output", str(out_spec),
                    "--output-individual", str(out_ind),
                    "--report-out", str(tmp_p / "tier_rep.csv"),
                    "--plot-out", str(tmp_p / "tier_fig.pdf"),
                    "--permutations", "19",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Dynamic H_opt Step 03 failed:\n{proc.stderr}\n{proc.stdout}")
            spec_df = pd.read_csv(out_spec)
            self.assertIn("efd_harmonic_opt", spec_df.columns)
            self.assertIn("efd_k95_pcs", spec_df.columns)
            h_opt_val = int(spec_df["efd_harmonic_opt"].iloc[0])
            k95_val = int(spec_df["efd_k95_pcs"].iloc[0])
            self.assertGreaterEqual(h_opt_val, 6)
            self.assertLessEqual(h_opt_val, 24)
            self.assertGreaterEqual(k95_val, 6, f"Expected k_95 >= 6 on 8-mode outlines, got {k95_val}")
            self.assertIn("PC6", spec_df.columns)
            self.assertIn("Sym_PC6", spec_df.columns)


if __name__ == "__main__":
    unittest.main()
