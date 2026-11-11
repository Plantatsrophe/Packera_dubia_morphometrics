#!/usr/bin/env python3
"""
===============================================================================
Script: test_fourier_extractor.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit tests for scripts/morphometrics/03_fourier_extractor.R verifying:
    1. Manifest lookup by (catalogNumber, leaf_id) for plant_individual_id,
       pixels_per_mm, scale_source, solidity, and archetype_score.
    2. Dissection-aware LeafMachine3 EFA QC gating via evaluate_step03_lm3_efa_qc().
    3. Physical millimeter calibration of Centroid Size (CS_mm) and oriented
       major-axis aspect ratio (Momocs::coo_lw / SVD alignment).
    4. Elimination of 180-degree first-harmonic ellipse flips via
       Momocs::coo_center(out_obj) %>% Momocs::coo_scale() and
       Momocs::efourier(..., norm = FALSE, start = FALSE).
    5. Strict Step 03 LM3 EFA QC filtering, OTU median aggregation,
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

    def test_manifest_lookup_and_lm3_efa_qc_gating(self):
        """Verify manifest lookup populates plant_individual_id, pixels_per_mm, scale_source, solidity, archetype_score, and calls evaluate_step03_lm3_efa_qc."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("manifest_df$catalogNumber", content)
        self.assertIn("manifest_df$leaf_id", content)
        for col in ("plant_individual_id", "pixels_per_mm", "is_calibrated_mm", "scale_source", "solidity", "archetype_score"):
            self.assertIn(col, content, f"Missing manifest column '{col}' in 03_fourier_extractor.R")
        self.assertIn("evaluate_step03_lm3_efa_qc", content)

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

    def test_basal_rosette_priority_and_pca_anchoring_in_step03(self):
        """Verify Basal_Rosette priority rule and Cauline_Only_Fallback separation in 03_fourier_extractor.R."""
        r_content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("basal_leaf_count", r_content)
        self.assertIn("cauline_leaf_count", r_content)
        self.assertIn("morph_tier_source", r_content)
        self.assertIn('"Basal_Rosette"', r_content)
        self.assertIn('"Cauline_Only_Fallback"', r_content)

    def test_empirical_efa_alignment_harmonic_power_and_otu_aggregation(self):
        """Verify Phase 3 empirical EFA protocol: coo_center, calibrate_harmonicpower_efourier, RRPP, and OTU-level PCA."""
        r_content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("Momocs::efourier", r_content)
        self.assertIn("norm = FALSE, start = FALSE", r_content)
        self.assertIn("Momocs::coo_center", r_content)
        self.assertIn("Momocs::coo_scale", r_content)
        self.assertIn("verify_harmonic_power_calibration", r_content)
        self.assertIn("calibrate_harmonicpower_efourier", r_content)
        self.assertIn("heteroblasty_index", r_content)
        self.assertIn("otu_t1_idx", r_content)
        self.assertIn("symmetrize_efa_harmonics", r_content)
        self.assertIn("compute_authority_precision_weights", r_content)

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
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.95, "archetype_score": 0.95, "contour_path": str(contours_dir / "G1_leaf1.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 2, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.90, "archetype_score": 0.90, "contour_path": str(contours_dir / "G1_leaf2.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 3, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.85, "archetype_score": 0.85, "contour_path": str(contours_dir / "G1_leaf3.csv")},
                {"catalogNumber": "G1", "plant_individual_id": 0, "leaf_id": 4, "efa_qc_pass": False, "qc_flag_reason": "UNRESOLVED_CLUMP", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.92, "archetype_score": 0.40, "contour_path": str(contours_dir / "G1_leaf4.csv")},
                {"catalogNumber": "G2", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.93, "archetype_score": 0.93, "contour_path": str(contours_dir / "G2_leaf1.csv")},
                {"catalogNumber": "G3", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.91, "archetype_score": 0.91, "contour_path": str(contours_dir / "G3_leaf1.csv")},
                {"catalogNumber": "S1", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.92, "archetype_score": 0.92, "contour_path": str(contours_dir / "S1_leaf1.csv")},
                {"catalogNumber": "S2", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.90, "archetype_score": 0.90, "contour_path": str(contours_dir / "S2_leaf1.csv")},
                {"catalogNumber": "B1", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.91, "archetype_score": 0.91, "contour_path": str(contours_dir / "B1_leaf1.csv")},
                {"catalogNumber": "T2_ONLY", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": False, "qc_flag_reason": "UNRESOLVED_CLUMP", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.89, "archetype_score": 0.35, "contour_path": str(contours_dir / "T2_ONLY_leaf1.csv")},
                {"catalogNumber": "EXCL_PLAT", "plant_individual_id": 0, "leaf_id": 1, "efa_qc_pass": True, "qc_flag_reason": "PASS", "pixels_per_mm": 10.0, "is_calibrated_mm": True, "solidity": 0.94, "archetype_score": 0.94, "contour_path": str(contours_dir / "EXCL_PLAT_leaf1.csv")},
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

            # 1. Only QC-passing contours present in vouchers_df enter leaf_efa_harmonics.csv
            self.assertEqual(len(ind1[ind1["catalogNumber"] == "G1"]), 3, "UNRESOLVED_CLUMP leaf 4 must be excluded by Step 03 EFA QC.")
            self.assertNotIn("T2_ONLY", set(spec1["catalogNumber"]), "UNRESOLVED_CLUMP-only specimen must not enter EFA table.")
            self.assertNotIn("EXCL_PLAT", set(ind1["catalogNumber"]), "Excluded taxon absent from vouchers_df must not be ingested into individual EFA table.")
            self.assertNotIn("EXCL_PLAT", set(spec1["catalogNumber"]), "Excluded taxon absent from vouchers_df must not enter specimen EFA table.")
            pub_utils_content = (PROJECT_ROOT / "scripts" / "morphometrics" / "pub_figures_utils.R").read_text(encoding="utf-8")
            self.assertNotIn("if (length(disp_tx) < 4) disp_tx <- unique(c(disp_tx, focal_4))[1:4]", pub_utils_content)

            for col in ("catalogNumber", "plant_individual_id", "pixels_per_mm", "is_calibrated_mm", "solidity", "archetype_score", "centroid_size", "log_centsize"):
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
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.42 if is_diss else 0.93,
                    "archetype_score": 0.94,
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
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.62,
                    "archetype_score": 0.94,
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

    def test_basal_positional_gate_lowest_lamina_fallback_ribbon_filter_and_otu_medians(self):
        """Verify basal positional homology filter, lowest-lamina fallback, and OTU median/q25 dissection summaries."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        def write_xy_csv(path: Path, cat: str, lid: int, xy: np.ndarray):
            pd.DataFrame({
                "catalogNumber": cat,
                "leaf_id": lid,
                "point_index": np.arange(len(xy)),
                "x": xy[:, 0] + 250.0,
                "y": xy[:, 1] + 250.0,
            }).to_csv(path, index=False)

        entire_xy = np.column_stack([85.0 * np.cos(t), 52.0 * np.sin(t)])
        u_t = np.cos(t)
        lobe_env = (0.25 + 0.75 * ((u_t + 1.0) / 2.0) ** 0.6) * (0.35 + 0.65 * np.abs(np.sin(4.5 * np.pi * (u_t + 1.0) / 2.0)))
        pinn_xy = np.column_stack([100.0 * u_t, 42.0 * np.sin(t) * lobe_env])

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = [
                {"catalogNumber": "V_MIX", "taxon": "Packera paupercula", "scientificName": "Packera paupercula", "species_raw": "Packera paupercula", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_DET", "taxon": "Packera paupercula", "scientificName": "Packera paupercula", "species_raw": "Packera paupercula", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF1", "taxon": "Packera obovata", "scientificName": "Packera obovata", "species_raw": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF2", "taxon": "Packera obovata", "scientificName": "Packera obovata", "species_raw": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            leaf_specs = [
                # V_MIX: leaf 1 is upper cauline (cy=0.18 -> excluded); leaf 2 entire basal; leaves 3-4 pinnatifid basal
                ("V_MIX", 1, entire_xy, 0.18, 0.94, False),
                ("V_MIX", 2, entire_xy, 0.72, 0.92, False),
                ("V_MIX", 3, pinn_xy, 0.78, 0.44, True),
                ("V_MIX", 4, pinn_xy * 1.03, 0.82, 0.46, True),
                # V_DET: top-mounted detached rosette (all cy < 0.40, rescued by max - 0.15 fallback)
                ("V_DET", 1, entire_xy * 0.98, 0.22, 0.91, False),
                ("V_DET", 2, pinn_xy * 0.97, 0.31, 0.45, True),
                # Reference vouchers
                ("V_REF1", 1, entire_xy * 1.02, 0.76, 0.93, False),
                ("V_REF2", 1, entire_xy * 0.96, 0.79, 0.92, False),
            ]
            manifest_rows = []
            for cat, lid, xy_arr, cy_p, sol_v, is_diss in leaf_specs:
                cpath = contours_dir / f"{cat}_leaf{lid}.csv"
                write_xy_csv(cpath, cat, lid, xy_arr)
                manifest_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "leaf_id": lid,
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": sol_v,
                    "archetype_score": 0.93,
                    "is_dissected": is_diss,
                    "cy_rel_plant_bbox": cy_p,
                    "cy_rel_sheet": cy_p,
                    "blade_area_cm2": 2.0,
                    "effective_min_blade_area_cm2": 0.4,
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
            self.assertEqual(proc.returncode, 0, f"Step 03 basal/dissection test failed:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            spec_df = pd.read_csv(out_spec)

            vmix_rows = ind_df.loc[ind_df["catalogNumber"] == "V_MIX"]
            if "leaf_zone" in vmix_rows.columns and len(vmix_rows) == 4:
                vmix_basal_lids = sorted(
                    vmix_rows.loc[vmix_rows["leaf_zone"].astype(str).str.lower() == "basal", "leaf_id"]
                    .astype(int)
                    .tolist()
                )
                vmix_caul_lids = sorted(
                    vmix_rows.loc[vmix_rows["leaf_zone"].astype(str).str.lower() == "cauline", "leaf_id"]
                    .astype(int)
                    .tolist()
                )
                self.assertEqual(vmix_basal_lids, [2, 3, 4], f"Expected basal leaves [2, 3, 4], got {vmix_basal_lids}")
                self.assertEqual(vmix_caul_lids, [1], f"Expected cauline leaf [1], got {vmix_caul_lids}")
            else:
                vmix_lids = sorted(vmix_rows["leaf_id"].astype(int).tolist())
                self.assertEqual(vmix_lids, [2, 3, 4], f"Expected basal leaves [2, 3, 4], got {vmix_lids}")
            vdet_lids = sorted(ind_df.loc[ind_df["catalogNumber"] == "V_DET", "leaf_id"].astype(int).tolist())
            self.assertEqual(vdet_lids, [1, 2], f"Expected detached rosette leaves rescued by fallback, got {vdet_lids}")

            vmix_otu = spec_df.loc[spec_df["catalogNumber"] == "V_MIX"].iloc[0]
            self.assertEqual(int(vmix_otu["basal_leaf_count"]), 3)
            # True median of [0.92, 0.44, 0.46] is 0.46 (NOT first-leaf 0.92)
            self.assertAlmostEqual(float(vmix_otu["solidity"]), 0.46, places=2)
            self.assertLess(float(vmix_otu["q25_solidity"]), 0.50)
            self.assertAlmostEqual(float(vmix_otu["dissected_leaf_fraction"]), 2.0 / 3.0, places=3)
            self.assertGreater(float(vmix_otu["dissection_index"]), 1.35)

    def test_dynamic_sqlite_synonym_resolution_and_tier_b_cohort_hygiene(self):
        """Verify dynamic SQLite taxon_synonyms lookup in R and Tier B orphan/synonym refresh in populate_run_active_cohort."""
        import os
        import sqlite3
        import sys
        if str(PROJECT_ROOT) not in sys.path:
            sys.path.insert(0, str(PROJECT_ROOT))
        from scripts.core.db import get_db_connection, init_db_schema, populate_run_active_cohort

        with tempfile.TemporaryDirectory() as tmp:
            eph_db = Path(tmp) / "ephemeral_synonym_test.sqlite"
            conn = get_db_connection(eph_db)
            init_db_schema(conn)
            with conn:
                conn.executemany(
                    "INSERT INTO taxon_synonyms (raw_name_lower, canonical_packera) VALUES (?, ?);",
                    [
                        ("packera synalpha", "Packera canonbeta"),
                        ("packera synalpha c.auth", "Packera canonbeta"),
                        ("senecio synalpha", "Packera canonbeta"),
                    ],
                )
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, qc_status) VALUES (?, ?, ?, ?, 'PASS');",
                    [
                        ("CAT_KEEP", "Packera synalpha C.Auth", "Packera canonbeta", "Tier_1_Gold"),
                        ("CAT_DROP", "Packera gammataxon", "Packera gammataxon", "Tier_2_Silver"),
                    ],
                )
                conn.executemany(
                    "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, ?);",
                    [("CAT_KEEP", 1), ("CAT_DROP", 1)],
                )
            populate_run_active_cohort(
                conn,
                run_id="test_run",
                voucher_rows=[
                    {"catalogNumber": "CAT_KEEP", "species_standardized": "Packera synalpha", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "CAT_DROP", "species_standardized": "Packera gammataxon", "determiner_tier": "Tier_2_Silver"},
                ],
            )
            with conn:
                conn.executemany(
                    "INSERT INTO otu_morphometrics (run_id, catalogNumber, plant_individual_id, species_standardized) VALUES (?, ?, ?, ?);",
                    [
                        ("test_run", "CAT_KEEP", 1, "Packera synalpha"),
                        ("test_run", "CAT_DROP", 1, "Packera gammataxon"),
                    ],
                )
            # Re-populate narrowed cohort with canonical standardized species
            populate_run_active_cohort(
                conn,
                run_id="test_run",
                voucher_rows=[
                    {"catalogNumber": "CAT_KEEP", "species_standardized": "Packera canonbeta", "determiner_tier": "Tier_1_Gold"},
                ],
            )
            rows_om = conn.execute(
                "SELECT catalogNumber, plant_individual_id, species_standardized FROM otu_morphometrics WHERE run_id = 'test_run';"
            ).fetchall()
            conn.close()
            self.assertEqual([(r[0], int(r[1]), r[2]) for r in rows_om], [("CAT_KEEP", 1, "Packera canonbeta")])

            if shutil.which("Rscript"):
                env = dict(os.environ)
                env["PACKERA_DB_PATH"] = str(eph_db)
                r_expr = (
                    "source('scripts/morphometrics/pub_figures_utils.R'); "
                    "res <- standardize_taxon_name(c('Packera synalpha', 'Packera synalpha C.Auth', 'P. synalpha', 'Senecio synalpha', 'Packera canonbeta var. minor')); "
                    "cat(paste(res, collapse='|'))"
                )
                proc = subprocess.run(
                    ["Rscript", "-e", r_expr],
                    cwd=str(PROJECT_ROOT),
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, f"R synonym test failed:\n{proc.stderr}\n{proc.stdout}")
                self.assertEqual(
                    proc.stdout.strip().splitlines()[-1],
                    "|".join(["Packera canonbeta"] * 5),
                )

    def test_vertical_tip_up_contour_horizontal_xax_alignment_before_rm_asym(self) -> None:
        """
        Verify align_contour_horizontal_xax in fourier_parallel_utils.R and 03_fourier_extractor.R
        rotates a vertically oriented tip-up symmetric leaf contour onto the horizontal X-axis
        before Momocs::rm_asym() so pct_asymmetry < 0.05 and is_contour_outlier == False.
        """
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")

        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = []
            manifest_rows = []
            for idx in range(6):
                cat = f"VERT_TIPUP_{idx:02d}"
                sp = "Packera dubia" if idx < 3 else "Packera anonyma"
                # Vertically oriented tip-up bilaterally symmetric leaf:
                # point 0 at lamina_base (x=250, y=340), point 60 at lamina_tip (x=250, y=160)
                # lateral width along X-axis symmetric around x=250
                a_len = 90.0 + idx * 1.5
                b_wid = 36.0 + (idx % 2) * 3.0
                y = 250.0 + a_len * np.cos(t)
                x = 250.0 + b_wid * np.sin(t) * (1.0 - 0.18 * np.cos(t))
                cpath = contours_dir / f"{cat}_leaf1.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": 1,
                    "point_index": np.arange(120),
                    "x": x,
                    "y": y,
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
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.94,
                    "archetype_score": 0.96,
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
            self.assertEqual(proc.returncode, 0, f"Step 03 vertical tip-up test failed:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            self.assertEqual(len(ind_df), 6)
            self.assertTrue(
                (ind_df["pct_asymmetry"].astype(float) < 0.05).all(),
                f"Expected pct_asymmetry < 0.05 on vertically oriented symmetric leaves, got {ind_df['pct_asymmetry'].tolist()}",
            )
            outlier_flags = ind_df["is_contour_outlier"].astype(str).str.upper()
            self.assertTrue(
                (outlier_flags == "FALSE").all(),
                f"Expected zero contour outliers on vertically oriented symmetric leaves, got {outlier_flags.tolist()}",
            )

    def test_step03_preserves_basal_and_cauline_tracks_with_pure_basal_otu_medians(self) -> None:
        """
        Verify that when a specimen has both leaf_zone = "basal" and leaf_zone = "cauline" leaves in Step 03:
          - Both basal and cauline leaves are preserved in leaf_efa_harmonics_ind.csv with their leaf_zone tags.
          - The specimen's primary OTU harmonics (A1..D_H, solidity) in specimen_efa_medians.csv match the
            pure basal leaf outlines (w_cauline = 0.0, uncontaminated by the cauline leaf), while the separate
            cauline track columns (cauline_leaf_count, cauline_solidity_med, cauline_dissection_index_med,
            basal_vs_cauline_length_ratio, heteroblasty_efa_distance) are populated.
        """
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")

        t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            contours_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_mixed_csv = tmp_p / "manifest_mixed.csv"
            manifest_basal_only_csv = tmp_p / "manifest_basal_only.csv"

            def write_contour(cat: str, lid: int, r_vec: np.ndarray) -> Path:
                x = 250.0 + r_vec * np.cos(t)
                y = 250.0 + r_vec * np.sin(t)
                p = contours_dir / f"{cat}_leaf{lid}.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": lid,
                    "point_index": np.arange(120),
                    "x": x,
                    "y": y,
                }).to_csv(p, index=False)
                return p

            # Basal rosette outline (broad ovate, solidity ~ 0.94) vs dissected cauline outline (deep lobes, solidity ~ 0.48)
            r_basal_1 = 80.0 * (1.0 + 0.12 * np.cos(t))
            r_basal_2 = 82.0 * (1.0 + 0.14 * np.cos(t))
            r_cauline = 55.0 * (1.0 + 0.46 * np.cos(6 * t))

            vouchers = [
                {"catalogNumber": "V_HETERO_1", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF_2", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF_3", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
                {"catalogNumber": "V_REF_4", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
            ]
            pd.DataFrame(vouchers).to_csv(vouchers_csv, index=False)

            p_h1_b1 = write_contour("V_HETERO_1", 1, r_basal_1)
            p_h1_b2 = write_contour("V_HETERO_1", 2, r_basal_2)
            p_h1_c3 = write_contour("V_HETERO_1", 3, r_cauline)

            basal_rows = [
                {
                    "catalogNumber": "V_HETERO_1",
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "leaf_zone": "basal",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.94,
                    "dissection_index": 1.08,
                    "blade_length_mm": 42.0,
                    "petiole_length_mm": 28.0,
                    "petiole_to_blade_ratio": 0.67,
                    "cy_rel_plant_bbox": 0.86,
                    "cy_rel_sheet": 0.80,
                    "blade_area_cm2": 4.2,
                    "effective_min_blade_area_cm2": 0.50,
                    "archetype_score": 0.95,
                    "contour_path": str(p_h1_b1),
                },
                {
                    "catalogNumber": "V_HETERO_1",
                    "plant_individual_id": 0,
                    "leaf_id": 2,
                    "leaf_zone": "basal",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.92,
                    "dissection_index": 1.10,
                    "blade_length_mm": 40.0,
                    "petiole_length_mm": 25.0,
                    "petiole_to_blade_ratio": 0.625,
                    "cy_rel_plant_bbox": 0.84,
                    "cy_rel_sheet": 0.78,
                    "blade_area_cm2": 4.0,
                    "effective_min_blade_area_cm2": 0.50,
                    "archetype_score": 0.95,
                    "contour_path": str(p_h1_b2),
                },
            ]
            for idx_ref, cat_ref in enumerate(["V_REF_2", "V_REF_3", "V_REF_4"], start=2):
                p_ref = write_contour(cat_ref, 1, 78.0 * (1.0 + (0.10 + 0.03 * idx_ref) * np.cos(2 * t)))
                basal_rows.append({
                    "catalogNumber": cat_ref,
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "leaf_zone": "basal",
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                    "solidity": 0.91,
                    "dissection_index": 1.12,
                    "blade_length_mm": 38.0,
                    "petiole_length_mm": 20.0,
                    "petiole_to_blade_ratio": 0.53,
                    "cy_rel_plant_bbox": 0.82,
                    "cy_rel_sheet": 0.76,
                    "blade_area_cm2": 3.6,
                    "effective_min_blade_area_cm2": 0.50,
                    "archetype_score": 0.94,
                    "contour_path": str(p_ref),
                })

            cauline_row = {
                "catalogNumber": "V_HETERO_1",
                "plant_individual_id": 0,
                "leaf_id": 3,
                "leaf_zone": "cauline",
                "pixels_per_mm": 10.0,
                "is_calibrated_mm": True,
                "solidity": 0.48,
                "dissection_index": 2.35,
                "blade_length_mm": 20.5,
                "petiole_length_mm": 0.8,
                "petiole_to_blade_ratio": 0.039,
                "cy_rel_plant_bbox": 0.20,
                "cy_rel_sheet": 0.22,
                "blade_area_cm2": 1.1,
                "effective_min_blade_area_cm2": 0.50,
                "archetype_score": 0.90,
                "contour_path": str(p_h1_c3),
            }

            pd.DataFrame(basal_rows).to_csv(manifest_basal_only_csv, index=False)
            pd.DataFrame(basal_rows + [cauline_row]).to_csv(manifest_mixed_csv, index=False)

            figs_b_dir = tmp_p / "figs_b"
            figs_m_dir = tmp_p / "figs_m"
            figs_b_dir.mkdir(parents=True, exist_ok=True)
            figs_m_dir.mkdir(parents=True, exist_ok=True)

            out_spec_basal = tmp_p / "spec_basal_only.csv"
            out_ind_basal = tmp_p / "ind_basal_only.csv"
            proc_b = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--contours-dir", str(contours_dir),
                    "--masks-dir", str(tmp_p / "masks"),
                    "--vouchers", str(vouchers_csv),
                    "--manifest", str(manifest_basal_only_csv),
                    "--output", str(out_spec_basal),
                    "--output-individual", str(out_ind_basal),
                    "--report-out", str(tmp_p / "rep_b.csv"),
                    "--plot-out", str(figs_b_dir / "fig_b.pdf"),
                    "--figures-dir", str(figs_b_dir),
                    "--permutations", "9",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_b.returncode, 0, f"Basal-only Step 03 failed:\n{proc_b.stderr}\n{proc_b.stdout}")

            out_spec_mix = tmp_p / "spec_mixed.csv"
            out_ind_mix = tmp_p / "ind_mixed.csv"
            proc_m = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--contours-dir", str(contours_dir),
                    "--masks-dir", str(tmp_p / "masks"),
                    "--vouchers", str(vouchers_csv),
                    "--manifest", str(manifest_mixed_csv),
                    "--output", str(out_spec_mix),
                    "--output-individual", str(out_ind_mix),
                    "--report-out", str(tmp_p / "rep_m.csv"),
                    "--plot-out", str(figs_m_dir / "fig_m.pdf"),
                    "--figures-dir", str(figs_m_dir),
                    "--permutations", "9",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_m.returncode, 0, f"Mixed basal+cauline Step 03 failed:\n{proc_m.stderr}\n{proc_m.stdout}")

            # Verify Step 03 publication figures separate basal and cauline morphospaces (EFA & ECT)
            self.assertTrue((figs_b_dir / "fig03_efa_morphospace_reconstruction.pdf").exists())
            self.assertFalse((figs_b_dir / "fig03c_cauline_efa_morphospace.pdf").exists())
            self.assertTrue((figs_b_dir / "fig03d_basal_ect_morphospace.pdf").exists())
            self.assertFalse((figs_b_dir / "fig03e_cauline_ect_morphospace.pdf").exists())
            self.assertTrue((figs_m_dir / "fig03_efa_morphospace_reconstruction.pdf").exists())
            self.assertTrue((figs_m_dir / "fig03c_cauline_efa_morphospace.pdf").exists())
            self.assertTrue((figs_m_dir / "fig03d_basal_ect_morphospace.pdf").exists())
            self.assertTrue((figs_m_dir / "fig03e_cauline_ect_morphospace.pdf").exists())

            ind_mix = pd.read_csv(out_ind_mix)
            spec_basal = pd.read_csv(out_spec_basal)
            spec_mix = pd.read_csv(out_spec_mix)

            # 1. Both basal and cauline leaves are preserved in leaf_efa_harmonics_ind.csv with leaf_zone tags
            h1_ind = ind_mix.loc[ind_mix["catalogNumber"] == "V_HETERO_1"].sort_values("leaf_id")
            self.assertEqual(h1_ind["leaf_id"].astype(int).tolist(), [1, 2, 3])
            self.assertEqual(
                h1_ind["leaf_zone"].astype(str).str.lower().tolist(),
                ["basal", "basal", "cauline"],
            )

            # 2. Primary OTU harmonics (A1..D6, solidity) in spec_mix match pure basal leaf outlines (w_cauline = 0.0)
            row_b = spec_basal.loc[spec_basal["catalogNumber"] == "V_HETERO_1"].iloc[0]
            row_m = spec_mix.loc[spec_mix["catalogNumber"] == "V_HETERO_1"].iloc[0]
            self.assertEqual(int(row_m["basal_leaf_count"]), 2)
            self.assertAlmostEqual(float(row_m["solidity"]), float(row_b["solidity"]), places=4)
            for h_col in ["A1", "A2", "A3", "D1", "D2", "D3"]:
                if h_col in row_m and h_col in row_b:
                    self.assertAlmostEqual(
                        float(row_m[h_col]),
                        float(row_b[h_col]),
                        places=5,
                        msg=f"Primary OTU harmonic {h_col} must match pure basal leaves without cauline contamination.",
                    )

            # 3. Separate cauline track columns (EFA & ECT) are populated on the mixed specimen
            self.assertEqual(int(row_m["cauline_leaf_count"]), 1)
            self.assertAlmostEqual(float(row_m["cauline_solidity_med"]), 0.48, places=2)
            self.assertGreater(float(row_m["cauline_dissection_index_med"]), 1.8)
            self.assertGreater(float(row_m["basal_vs_cauline_length_ratio"]), 1.5)
            self.assertGreater(float(row_m["heteroblasty_efa_distance"]), 0.01)
            self.assertIn("ECT_PC1", ind_mix.columns)
            self.assertIn("ect_lobe_complexity", ind_mix.columns)
            self.assertIn("ECT_PC1", spec_mix.columns)
            self.assertIn("Cauline_ECT_PC1", spec_mix.columns)
            self.assertIn("ect_lobe_complexity", spec_mix.columns)
            self.assertIn("cauline_ect_lobe_complexity", spec_mix.columns)
            self.assertIn("heteroblasty_ect_distance", spec_mix.columns)
            self.assertTrue(np.isfinite(float(row_m["ECT_PC1"])))
            self.assertTrue(np.isfinite(float(row_m["ect_lobe_complexity"])))
            self.assertTrue(np.isfinite(float(row_m["cauline_ect_lobe_complexity"])))
            self.assertGreater(float(row_m["heteroblasty_ect_distance"]), 0.0)

    def test_step03_rigid_tipup_dissected_homology_reniform_lw_and_32x360_ect(self) -> None:
        """
        Verify Step 02 -> Step 03 rigid (-dy, -dx) contour handoff, reniform L/W < 1.0 preservation,
        mp_estimated centroid_size gating, and dynamic 32x360 .npy ECT ingestion.
        """
        # 1. Verify fourier_parallel_utils.R and dissected_leaf_utils.R static contracts
        fpu_content = (PROJECT_ROOT / "scripts" / "morphometrics" / "fourier_parallel_utils.R").read_text(encoding="utf-8")
        dlu_content = (PROJECT_ROOT / "scripts" / "morphometrics" / "dissected_leaf_utils.R").read_text(encoding="utf-8")
        dbu_content = (PROJECT_ROOT / "scripts" / "core" / "db_utils.R").read_text(encoding="utf-8")
        self.assertIn("is_step02_tipup", fpu_content)
        self.assertIn("cbind(-dy + cx, -dx + cy)", fpu_content)
        self.assertNotIn("floor(nrow(c_xy) / 2L) + 1L", fpu_content)
        self.assertIn("524288L", dlu_content)
        self.assertIn(".resample_ect_directions_r", dlu_content)
        self.assertIn("mp_estimated", dbu_content)

        # Mathematical verification: asymmetrically pinnatifid Step 02 tip-up leaf (extra lobes on left margin)
        # preserves exact 0.000-deg midrib alignment under (-dy + cx, -dx + cy) whereas vertex-61 drifts off-axis
        t = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
        # Clockwise in image space (x_right, y_down): t=0 base (250, 340), t=pi/2 left (x<250), t=pi tip (250, 160), t=3pi/2 right (x>250)
        y_img = 250.0 + 90.0 * np.cos(t)
        left_mask = np.sin(t) > 0.0
        lobe_asym = np.where(left_mask, 1.0 + 0.42 * np.sin(5.0 * t) ** 2, 1.0)
        x_img = 250.0 - 38.0 * np.sin(t) * lobe_asym
        # Resample by equal perimeter arc-length (just like Step 02 resample_closed) so asymmetric left perimeter shifts index 60 off tip!
        xy_raw = np.column_stack([x_img, y_img])
        closed = np.vstack([xy_raw, xy_raw[0:1]])
        seg = np.hypot(np.diff(closed[:, 0]), np.diff(closed[:, 1]))
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        targets = np.linspace(0.0, cum[-1], 120, endpoint=False)
        coo_step02 = np.column_stack([np.interp(targets, cum, closed[:, 0]), np.interp(targets, cum, closed[:, 1])])

        # Old vertex-60 heuristic drifts significantly off the vertical midrib x=250:
        old_drift_deg = abs(np.degrees(np.arctan2(coo_step02[60, 0] - coo_step02[0, 0], coo_step02[0, 1] - coo_step02[60, 1])))
        self.assertGreater(old_drift_deg, 3.0)

        # New Step 03 rigid transform (-dy + cx, -dx + cy) maps true midrib (x_img == 250) onto constant Y == cy_m (0.000-deg drift)
        cx_m, cy_m = float(coo_step02[:, 0].mean()), float(coo_step02[:, 1].mean())
        dx_m, dy_m = coo_step02[:, 0] - cx_m, coo_step02[:, 1] - cy_m
        momocs_xy = np.column_stack([-dy_m + cx_m, -dx_m + cy_m])
        # Clockwise in Cartesian Momocs space (shoelace signed area < 0)
        shoelace = 0.5 * np.sum(momocs_xy[:, 0] * np.roll(momocs_xy[:, 1], -1) - np.roll(momocs_xy[:, 0], -1) * momocs_xy[:, 1])
        self.assertLess(shoelace, 0.0)
        s_03_chk = (PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R").read_text(encoding="utf-8")
        s_db_chk = (PROJECT_ROOT / "scripts" / "core" / "db_utils.R").read_text(encoding="utf-8")
        self.assertIn("allow_multi_component_ect = TRUE", s_03_chk)
        self.assertIn("allow_multi_component_ect = TRUE", s_db_chk)

        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            contours_dir = tmp_p / "contours"
            ect_dir = tmp_p / "ect"
            contours_dir.mkdir(parents=True, exist_ok=True)
            ect_dir.mkdir(parents=True, exist_ok=True)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            manifest_csv = tmp_p / "extracted_leaf_manifest.csv"

            vouchers = []
            manifest_rows = []
            for idx in range(6):
                cat = f"STEP23_HO_{idx:02d}"
                is_reniform = idx >= 3
                sp = "Packera aurea" if is_reniform else "Packera paupercula"
                if is_reniform:
                    # Reniform wider-than-long tip-up blade: L = 110 (y span), W = 160 (x span) -> L/W = 0.6875 < 1.0
                    y_c = 250.0 + 55.0 * np.cos(t)
                    x_c = 250.0 - 80.0 * np.sin(t)
                    lw_true = 110.0 / 160.0
                    sol_v = 0.93
                else:
                    y_c = coo_step02[:, 1]
                    x_c = coo_step02[:, 0]
                    lw_true = 2.15
                    sol_v = 0.64

                cpath = contours_dir / f"{cat}_leaf1.csv"
                pd.DataFrame({
                    "catalogNumber": cat,
                    "leaf_id": 1,
                    "point_index": np.arange(120),
                    "x": x_c,
                    "y": y_c,
                }).to_csv(cpath, index=False)

                epath = ect_dir / f"{cat}_leaf1_ect.npy"
                theta_grid = np.linspace(0.0, 2.0 * np.pi, 360, endpoint=False, dtype=np.float32)
                rad_grid = np.linspace(-1.0, 1.0, 32, dtype=np.float32)[:, None]
                ect_32x360 = (
                    np.exp(-(rad_grid ** 2) * (2.0 + (0.0 if is_reniform else 1.5 * np.cos(4.0 * theta_grid)[None, :])))
                ).astype(np.float32)
                np.save(str(epath), ect_32x360)

                is_mp = idx == 5
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
                    "pixels_per_mm": 10.5 if is_mp else 10.0,
                    "is_calibrated_mm": (not is_mp),
                    "scale_source": "mp_estimated" if is_mp else "ruler_calibrated",
                    "aspect_ratio": round(lw_true, 4),
                    "blade_length_to_width_ratio": round(lw_true, 4),
                    "solidity": sol_v,
                    "archetype_score": 0.94,
                    "contour_path": str(cpath),
                    "ect_path": str(epath),
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
                    "--output-ind", str(out_ind),
                    "--figures-dir", str(tmp_p / "figs"),
                    "--harmonics", "12",
                    "--num-pcs", "5",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, f"Step 03 handoff Rscript test failed:\n{proc.stderr}\n{proc.stdout}")
            ind_df = pd.read_csv(out_ind)
            ren_rows = ind_df[ind_df["catalogNumber"].isin(["STEP23_HO_03", "STEP23_HO_04", "STEP23_HO_05"])]
            self.assertTrue((ren_rows["aspect_ratio"].astype(float) < 0.80).all(), "Reniform Packera aurea leaves must retain aspect_ratio < 1.0 in Step 03.")
            mp_row = ind_df[ind_df["catalogNumber"] == "STEP23_HO_05"].iloc[0]
            self.assertTrue(pd.isna(mp_row["centroid_size"]), "mp_estimated leaf must keep centroid_size == NA in Step 03.")
            self.assertTrue(ind_df["ECT_PC1"].notna().all(), "Dynamic 32x360 .npy ECT reader must populate ECT_PC1 for all leaves.")


if __name__ == "__main__":
    unittest.main()
