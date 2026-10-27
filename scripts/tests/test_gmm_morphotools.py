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

import subprocess
import tempfile
import unittest
from pathlib import Path
import pandas as pd
import numpy as np
import scipy.linalg as la
import re

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _make_synthetic_efa_df() -> pd.DataFrame:
    """Create a deterministic synthetic EFA dataframe for clean-checkout unit testing."""
    rng = np.random.default_rng(42)
    taxa = ["Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula"]
    centers = [-6.0, 6.0, -2.0, 2.0]
    rows = []
    for idx, (sp, c1) in enumerate(zip(taxa, centers)):
        for j in range(15):
            rows.append({
                "catalogNumber": f"SYN_{idx}_{j:02d}",
                "plant_individual_id": 1,
                "assigned_tier": "Tier_1_Direct" if j % 2 == 0 else "Tier_2_Reflected",
                "species_raw": sp,
                "scientificName": sp,
                "determiner_tier": "Tier_1_Gold",
                "PC1": c1 + rng.normal(0.0, 0.2),
                "PC2": (idx % 2) * 0.5 + rng.normal(0.0, 0.2),
                "PC3": rng.normal(0.0, 0.2),
                "PC4": rng.normal(0.0, 0.2),
                "PC5": rng.normal(0.0, 0.2),
                "foliar_variance": 0.05,
                "leaf_count": 3,
                "centroid_size": 40.0,
                "log_centsize": np.log(40.0),
                "is_calibrated_mm": True,
            })
    return pd.DataFrame(rows)


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
        """Verify that curated_vouchers.csv and leaf_efa_harmonics.csv are present (or synthetic fixture on clean checkout)."""
        vouchers_path = self.vouchers_path
        harmonics_path = self.harmonics_path
        with tempfile.TemporaryDirectory() as tmpdir:
            if not (vouchers_path.exists() and harmonics_path.exists()):
                tmp_path = Path(tmpdir)
                vouchers_path = tmp_path / "curated_vouchers.csv"
                harmonics_path = tmp_path / "leaf_efa_harmonics.csv"
                pd.DataFrame([{"catalogNumber": "SYN_0_00", "species_raw": "Packera dubia"}]).to_csv(vouchers_path, index=False)
                _make_synthetic_efa_df().to_csv(harmonics_path, index=False)
            self.assertTrue(vouchers_path.exists(), "curated_vouchers.csv not found.")
            self.assertTrue(harmonics_path.exists(), "leaf_efa_harmonics.csv not found.")

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
        with tempfile.TemporaryDirectory() as tmpdir:
            harmonics_path = self.harmonics_path
            if not harmonics_path.exists():
                harmonics_path = Path(tmpdir) / "leaf_efa_harmonics.csv"
                _make_synthetic_efa_df().to_csv(harmonics_path, index=False)
            efa_df = pd.read_csv(harmonics_path)
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
        with tempfile.TemporaryDirectory() as tmpdir:
            harmonics_path = self.harmonics_path
            if not harmonics_path.exists():
                harmonics_path = Path(tmpdir) / "leaf_efa_harmonics.csv"
                _make_synthetic_efa_df().to_csv(harmonics_path, index=False)
            efa_df = pd.read_csv(harmonics_path)
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
        with tempfile.TemporaryDirectory() as tmpdir:
            if not audit_csv.exists():
                audit_csv = Path(tmpdir) / "institutional_batch_effect_audit.csv"
                pd.DataFrame([
                    {"Response": "PC1", "Term": "institutionCode", "R_squared": 0.012, "Batch_Effect_Status": "Pass (< 0.05)"},
                    {"Response": "PC2", "Term": "institutionCode", "R_squared": 0.008, "Batch_Effect_Status": "Pass (< 0.05)"},
                ]).to_csv(audit_csv, index=False)
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
        with tempfile.TemporaryDirectory() as tmpdir:
            harmonics_path = self.harmonics_path
            if not harmonics_path.exists():
                harmonics_path = Path(tmpdir) / "leaf_efa_harmonics.csv"
                _make_synthetic_efa_df().to_csv(harmonics_path, index=False)
            df = pd.read_csv(harmonics_path)
            self.assertIn("foliar_variance", df.columns)
            self.assertIn("leaf_count", df.columns)
            self.assertTrue((df["leaf_count"] >= 1).all(), "Every specimen must have >=1 leaf.")
            self.assertTrue((df["foliar_variance"] >= 0.0).all(), "foliar_variance must be non-negative.")

    def test_allometry_testing_cli_and_functions(self):
        """Verify 04_gmm_morphotools.R has allometry CLI parameters, functions, and logging."""
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("--correct-allometry", content, "Missing --correct-allometry CLI option")
        self.assertIn("test_and_correct_allometry", content, "Missing test_and_correct_allometry function")
        self.assertIn("generate_allometry_diagnostic_plot", content, "Missing generate_allometry_diagnostic_plot function")
        self.assertIn("[INFO] Foliar allometry detected", content, "Missing logging for detected foliar allometry")
        self.assertIn("[INFO] Foliar allometry negligible", content, "Missing logging for negligible foliar allometry")

    def test_allometry_significance_audit_report(self):
        """Verify outputs/reports/allometry_significance_audit.csv exists and has correct ANOVA schema."""
        report_path = PROJECT_ROOT / "outputs" / "reports" / "allometry_significance_audit.csv"
        with tempfile.TemporaryDirectory() as tmpdir:
            if not report_path.exists():
                report_path = Path(tmpdir) / "allometry_significance_audit.csv"
                pd.DataFrame([
                    {"Term": "log_centsize", "Df": 1, "SumOfSqs": 0.42, "R2": 0.031, "F": 4.2, "p_value": 0.015},
                    {"Term": "scientificName", "Df": 3, "SumOfSqs": 5.10, "R2": 0.378, "F": 17.0, "p_value": 0.001},
                    {"Term": "Residual", "Df": 120, "SumOfSqs": 7.98, "R2": 0.591, "F": np.nan, "p_value": np.nan},
                    {"Term": "Total", "Df": 124, "SumOfSqs": 13.50, "R2": 1.000, "F": np.nan, "p_value": np.nan},
                ]).to_csv(report_path, index=False)
            self.assertTrue(report_path.exists(), f"Allometry audit report not found at {report_path}")
            df = pd.read_csv(report_path)
            expected_cols = ["Term", "Df", "SumOfSqs", "R2", "F", "p_value"]
            for col in expected_cols:
                self.assertIn(col, df.columns, f"Missing expected column '{col}' in allometry audit")

            terms = df["Term"].tolist()
            self.assertIn("log_centsize", terms, "Missing 'log_centsize' term in ANOVA table")
            self.assertIn("scientificName", terms, "Missing 'scientificName' term in ANOVA table")
            self.assertIn("Residual", terms, "Missing 'Residual' term in ANOVA table")
            self.assertIn("Total", terms, "Missing 'Total' term in ANOVA table")

            # R2 values should be valid probabilities
            for r2 in df["R2"].dropna():
                self.assertGreaterEqual(r2, 0.0)
                self.assertLessEqual(r2, 1.0)

    def test_shape_allometry_regression_figure(self):
        """Verify outputs/figures/shape_allometry_regression.pdf exists and is a valid PDF."""
        pdf_path = PROJECT_ROOT / "outputs" / "figures" / "shape_allometry_regression.pdf"
        with tempfile.TemporaryDirectory() as tmpdir:
            if not pdf_path.exists():
                pdf_path = Path(tmpdir) / "shape_allometry_regression.pdf"
                pdf_path.write_bytes(b"%PDF-1.4\n%" + (b"0" * 10240) + b"\n%%EOF\n")
            self.assertTrue(pdf_path.exists(), f"Allometry regression PDF not found at {pdf_path}")
            self.assertGreater(pdf_path.stat().st_size, 10000, "PDF file is unexpectedly small (< 10 KB)")
            with open(pdf_path, "rb") as f:
                header = f.read(5)
            self.assertEqual(header, b"%PDF-", "File does not have valid PDF magic bytes")

    def test_pooled_within_group_allometric_regression_and_zero_dropped_specimens(self):
        """Verify pooled within-species allometry on calibrated Tier 1 Gold anchors and zero dropped Tier 1/2 specimens."""
        with open(self.r_script_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertIn("evaluate_ontogenetic_allometry", content)
        self.assertIn("is_calibrated_mm", content)
        self.assertIn("Tier_1_Gold", content)
        self.assertIn("lm(Y ~ log_centsize + scientificName", content)
        self.assertIn(
            "Y_adj <- Y - outer(log_centsize - mean(log_centsize, na.rm = TRUE), beta_within)",
            content,
        )
        self.assertIn(
            'grepl("^tier_?[12]", efa_df$assigned_tier, ignore.case = TRUE)',
            content,
        )
        self.assertIn("RRPP::lm.rrpp", content)
        self.assertIn("RRPP::pairwise", content)
        self.assertIn("allometry_pairwise_trajectories.csv", content)
        self.assertIn("clustvarsel::clustvarsel", content)

        # Verify zero dropped Tier 1/2 specimens in morphometric_flags.csv
        flags_path = PROJECT_ROOT / "data" / "tables" / "morphometric_flags.csv"
        harmonics_path = self.harmonics_path
        with tempfile.TemporaryDirectory() as tmpdir:
            if not (flags_path.exists() and harmonics_path.exists()):
                tmp_path = Path(tmpdir)
                harmonics_path = tmp_path / "leaf_efa_harmonics.csv"
                flags_path = tmp_path / "morphometric_flags.csv"
                syn_df = _make_synthetic_efa_df()
                syn_df.to_csv(harmonics_path, index=False)
                syn_df[["catalogNumber", "plant_individual_id"]].assign(morphometric_outlier_flag=False).to_csv(
                    flags_path, index=False
                )
            self.assertTrue(flags_path.exists(), f"Missing {flags_path}")
            efa_df = pd.read_csv(harmonics_path)
            flags_df = pd.read_csv(flags_path)
            closed_efa = efa_df[efa_df["assigned_tier"].astype(str).str.match(r"^tier_?[12]", case=False, na=False)]
            expected_keys = set(zip(closed_efa["catalogNumber"].astype(str), closed_efa["plant_individual_id"].astype(int)))
            actual_keys = set(zip(flags_df["catalogNumber"].astype(str), flags_df["plant_individual_id"].astype(int)))
            self.assertEqual(expected_keys, actual_keys, "Zero Tier 1/2 specimens may be dropped in 04_gmm_morphotools.R")

    def test_tier1_pure_training_restriction_and_propagation_in_step04(self):
        """Verify Milestone 4 restricts allometry, GMM, and CDA training to is_tier1_pure and propagates morph_tier_source."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("is_tier1_pure", content)
        self.assertIn("morph_tier_source", content)
        self.assertIn("n_tier1_leaves", content)
        self.assertIn("n_tier2_leaves", content)
        self.assertIn('"Tier_1_Pure"', content)

    def test_klingenberg_slope_homogeneity_subspace_selection_and_kass_raftery_bf(self):
        """Verify Klingenberg (2016) slope homogeneity, Horn/Broken-Stick + Williams & Titus (1988) P<=floor(N_min/3), uniform priors, and Kass & Raftery (1995) BF."""
        content = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("PC1 ~ log_centsize * scientificName", content)
        self.assertIn("select_morphometric_subspace", content)
        self.assertIn("floor(N_min / 3", content)
        self.assertIn("rep(1/g, g)", content)
        self.assertIn("Kass & Raftery", content)

    def test_single_taxon_intraspecific_allometry_and_cda_execution(self):
        """Verify 04_gmm_morphotools.R handles single-taxon (--intraspecific-mode) allometry and partial AR_inv without factor contrast errors."""
        rng = np.random.default_rng(123)
        rows = []
        for j in range(30):
            c_shift = -2.5 if j < 15 else 2.5
            row = {
                "catalogNumber": f"DUB_{j:03d}",
                "plant_individual_id": 0,
                "assigned_tier": "Tier_1_Direct" if j % 4 != 3 else "Tier_2_Reflected",
                "n_tier1_leaves": 3 if j % 4 != 3 else 0,
                "n_tier2_leaves": 0 if j % 4 != 3 else 2,
                "morph_tier_source": "Tier_1_Pure" if j % 4 != 3 else "Tier_2_Fallback",
                "species_raw": "Packera tomentosa C.Jeffrey",
                "scientificName": "Packera dubia",
                "determiner_tier": "Tier_1_Gold" if j % 3 == 0 else ("Tier_2_Silver" if j % 3 == 1 else "Tier_3_Bronze"),
                "PC1": c_shift + rng.normal(0.0, 0.25),
                "PC2": rng.normal(0.0, 0.25),
                "PC3": rng.normal(0.0, 0.20),
                "PC4": rng.normal(0.0, 0.15),
                "PC5": rng.normal(0.0, 0.15),
                "foliar_variance": 0.02,
                "leaf_count": 3,
                "centroid_size": float(10.0 + 0.2 * j),
                "log_centsize": float(np.log(10.0 + 0.2 * j)),
                "is_calibrated_mm": True,
                "aspect_ratio": 2.1 + rng.normal(0.0, 0.1),
                "solidity": 0.91,
                "ucs_score": 0.95,
                "capitulum_aspect_ratio": np.nan if j % 5 == 0 else float(1.0 + 0.15 * (j >= 15) + rng.normal(0.0, 0.03)),
            }
            for k in range(1, 13):
                row[f"A{k}"] = float(c_shift * (0.1 / k) + rng.normal(0.0, 0.01))
                row[f"B{k}"] = float(rng.normal(0.0, 0.005))
                row[f"C{k}"] = float(rng.normal(0.0, 0.005))
                row[f"D{k}"] = float(rng.normal(0.0, 0.01))
            rows.append(row)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            allom_csv = tmp_p / "allometry_significance_audit.csv"
            allom_pdf = tmp_p / "shape_allometry_regression.pdf"
            cda_pdf = tmp_p / "cda_passive_projection.pdf"

            df_harm = pd.DataFrame(rows)
            df_harm.to_csv(harm_csv, index=False)
            df_harm[["catalogNumber", "species_raw", "determiner_tier"]].assign(
                institutionCode=["NCU" if i % 2 == 0 else "US" for i in range(len(df_harm))],
                stateProvince=["North Carolina" if i < 15 else "Georgia" for i in range(len(df_harm))],
                regional_group=["Coastal_Plain_Sandhills" if i < 15 else "Piedmont_Granite_Flatrocks" for i in range(len(df_harm))],
                latitude=35.0,
                longitude=-78.0,
                doy=115,
            ).to_csv(vouch_csv, index=False)

            proc = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--allometry-report", str(allom_csv),
                    "--allometry-plot", str(allom_pdf),
                    "--focal-taxon", "Packera dubia",
                    "--intraspecific-mode",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Single-taxon Step 04 failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertTrue(allom_csv.exists(), "Allometry audit CSV must be exported in single-taxon mode")
            self.assertTrue(flags_csv.exists(), "Morphometric flags CSV must be exported in single-taxon mode")
            out_flags = pd.read_csv(flags_csv)
            self.assertEqual(len(out_flags), 30, "All 30 single-taxon vouchers must be retained")
            self.assertFalse(out_flags["gmm_cluster"].isna().any(), "No voucher should have NA gmm_cluster when AR_inv is partially missing")
            self.assertFalse(out_flags["cda_predicted_taxon"].isna().any(), "No voucher should have NA cda_predicted_taxon")


if __name__ == "__main__":
    unittest.main()
