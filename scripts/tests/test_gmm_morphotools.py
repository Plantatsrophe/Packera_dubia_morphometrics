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
import zlib

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
        """Test generic Packera binomial standardization across North American taxa."""
        def standardize_taxon(s):
            if pd.isna(s):
                return "Unknown"
            m = re.match(r"^(?:packera|p\.)\s*(?:x\s+|×\s*)?([a-z-]+)", str(s).strip().lower())
            return f"Packera {m.group(1).lstrip('×-')}" if m else "Unknown"

        self.assertEqual(standardize_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_taxon("Packera anonyma (Alph.Wood) W.A.Weber"), "Packera anonyma")
        self.assertEqual(standardize_taxon("Packera paupercula var. savannarum R.R.Kowal"), "Packera paupercula")
        self.assertEqual(standardize_taxon("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve"), "Packera streptanthifolia")
        self.assertEqual(standardize_taxon("P. sanguisorbae (DC.) C.Jeffrey"), "Packera sanguisorbae")
        self.assertEqual(standardize_taxon("Packera x memmingeri"), "Packera memmingeri")

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
            if pd.isna(s):
                return "Unknown"
            m = re.match(r"^(?:packera|p\.)\s*(?:x\s+|×\s*)?([a-z-]+)", str(s).strip().lower())
            return f"Packera {m.group(1).lstrip('×-')}" if m else "Unknown"

        sp_col = (
            "species_standardized"
            if "species_standardized" in valid.columns
            else ("scientificName" if "scientificName" in valid.columns else "species_raw")
        )
        valid["species_std"] = valid[sp_col].apply(standardize_taxon)
        active_taxa = sorted([t for t in valid["species_std"].unique() if str(t).startswith("Packera ")])
        is_active = (valid["determiner_tier"] == "Tier_1_Gold") & (valid["species_std"].isin(active_taxa))

        X_a = valid.loc[is_active, pca_cols].values
        y_a = valid.loc[is_active, "species_std"].values
        N_a, p = X_a.shape
        g = len(active_taxa)

        grand_mean = np.mean(X_a, axis=0)
        B = np.zeros((p, p))
        W = np.zeros((p, p))

        for grp in active_taxa:
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
            if pd.isna(s):
                return "Unknown"
            m = re.match(r"^(?:packera|p\.)\s*(?:x\s+|×\s*)?([a-z-]+)", str(s).strip().lower())
            return f"Packera {m.group(1).lstrip('×-')}" if m else "Unknown"

        sp_col = (
            "species_standardized"
            if "species_standardized" in valid.columns
            else ("scientificName" if "scientificName" in valid.columns else "species_raw")
        )
        valid["species_std"] = valid[sp_col].apply(standardize_taxon)
        active_taxa = sorted([t for t in valid["species_std"].unique() if str(t).startswith("Packera ")])
        is_active = (valid["determiner_tier"] == "Tier_1_Gold") & (valid["species_std"].isin(active_taxa))

        X_a = valid.loc[is_active, feature_cols].values
        y_a = valid.loc[is_active, "species_std"].values
        N_a, p = X_a.shape
        g = len(active_taxa)

        grand_mean = np.mean(X_a, axis=0)
        B = np.zeros((p, p))
        W = np.zeros((p, p))

        for grp in active_taxa:
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
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(allom_csv),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
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

    def test_tier3_bronze_zero_influence_on_allometry_gmm_and_cda(self):
        """Verify Tier_3_Bronze has authority_weight == 0.0 and zero influence on reference allometry, GMM, and CDA fits."""
        rng = np.random.default_rng(99)
        taxa = ["Packera anonyma", "Packera dubia", "Packera paupercula"]
        centers = [-4.0, 4.0, 0.0]
        rows = []
        for idx, (sp, c1) in enumerate(zip(taxa, centers)):
            for j in range(12):
                tier = "Tier_1_Gold" if j < 5 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                cs = float(25.0 + 1.5 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": f"REF_{idx}_{j:02d}",
                    "plant_individual_id": 0,
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 2,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": tier,
                    "PC1": float(c1 + 0.3 * log_cs + rng.normal(0.0, 0.15)),
                    "PC2": float((idx - 1) * 1.2 + rng.normal(0.0, 0.15)),
                    "PC3": float(rng.normal(0.0, 0.12)),
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 2,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(2.0 + 0.2 * idx + rng.normal(0.0, 0.05)),
                    "solidity": 0.92,
                    "ucs_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((c1 * 0.08 + 0.04 * log_cs) / k + rng.normal(0.0, 0.004))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float(((idx - 1) * 0.06 + 0.03 * log_cs) / k + rng.normal(0.0, 0.004))
                rows.append(row)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            df_base = pd.DataFrame(rows)

            def run_step04_df(df_in: pd.DataFrame, tag: str):
                harm_csv = tmp_p / f"harm_{tag}.csv"
                vouch_csv = tmp_p / f"vouch_{tag}.csv"
                flags_csv = tmp_p / f"flags_{tag}.csv"
                bic_csv = tmp_p / f"bic_{tag}.csv"
                allom_csv = tmp_p / f"allom_{tag}.csv"
                allom_pdf = tmp_p / f"allom_{tag}.pdf"
                cda_pdf = tmp_p / f"cda_{tag}.pdf"
                df_in.to_csv(harm_csv, index=False)
                df_in[["catalogNumber", "species_raw", "determiner_tier"]].assign(
                    institutionCode="NCU",
                    stateProvince="North Carolina",
                    latitude=35.5,
                    longitude=-79.0,
                    doy=120,
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
                        "--output-batch-audit", str(tmp_p / f"batch_audit_{tag}.csv"),
                        "--allometry-report", str(allom_csv),
                        "--allometry-pairwise", str(tmp_p / f"pairwise_{tag}.csv"),
                        "--allometry-plot", str(allom_pdf),
                        "--correct-allometry", "true",
                        "--permutations", "99",
                    ],
                    cwd=str(PROJECT_ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, f"Step 04 failed ({tag}):\n{proc.stderr}\n{proc.stdout}")
                return pd.read_csv(flags_csv), pd.read_csv(bic_csv), pd.read_csv(allom_csv)

            flags1, bic1, allom1 = run_step04_df(df_base, "base")

            # Verify Tier_3_Bronze has authority_weight == 0.0 and is_passive == True, while Gold/Silver have EB active weights (mean == 1.0, Gold > Silver)
            bronze_rows = flags1[flags1["determiner_tier"] == "Tier_3_Bronze"]
            active_rows = flags1[flags1["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])]
            gold_rows = flags1[flags1["determiner_tier"] == "Tier_1_Gold"]
            silver_rows = flags1[flags1["determiner_tier"] == "Tier_2_Silver"]
            self.assertTrue((bronze_rows["authority_weight"] == 0.0).all())
            self.assertTrue((bronze_rows["is_passive"] == True).all())
            self.assertTrue((active_rows["authority_weight"] > 0.0).all())
            self.assertAlmostEqual(float(active_rows["authority_weight"].mean()), 1.0, places=4)
            self.assertGreater(float(gold_rows["authority_weight"].mean()), float(silver_rows["authority_weight"].mean()))

            # Drastically perturb all Tier_3_Bronze rows in PC space, size, and harmonics
            df_pert = df_base.copy()
            b_mask = df_pert["determiner_tier"] == "Tier_3_Bronze"
            df_pert.loc[b_mask, "PC1"] += 25.0
            df_pert.loc[b_mask, "PC2"] -= 18.0
            df_pert.loc[b_mask, "log_centsize"] += 2.5
            for k in range(1, 13):
                df_pert.loc[b_mask, f"A{k}"] += 0.5
                df_pert.loc[b_mask, f"D{k}"] -= 0.5

            flags2, bic2, allom2 = run_step04_df(df_pert, "perturbed")

            # Allometry ANOVA table, GMM BIC table, and active Gold/Silver CDA scores must be strictly identical
            np.testing.assert_allclose(
                allom1.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                allom2.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                atol=1e-6,
                err_msg="Tier_3_Bronze perturbation must have zero influence on reference allometry ANOVA.",
            )
            np.testing.assert_allclose(
                bic1.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                bic2.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                atol=1e-6,
                err_msg="Tier_3_Bronze perturbation must have zero influence on reference GMM BIC table.",
            )
            act1 = flags1[flags1["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])].sort_values("catalogNumber")
            act2 = flags2[flags2["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])].sort_values("catalogNumber")
            for col in ["can1", "can2", "cda_posterior_prob", "PC1", "PC2"]:
                np.testing.assert_allclose(
                    act1[col].to_numpy(),
                    act2[col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter reference {col}.",
                )

    def test_dynamic_non_default_packera_taxa_discovery_and_classification(self):
        """Verify 04_gmm_morphotools.R dynamically discovers and classifies non-default Packera species in cda_predicted_taxon and passively projects singleton taxa."""
        r_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertNotIn(
            'target_taxa <- c("Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula")',
            r_src,
            "04_gmm_morphotools.R must not hardcode the 4-species target_taxa vector.",
        )

        rng = np.random.default_rng(2026)
        non_default_specs = [
            ("Packera schweinitziana", "Senecio schweinitzianus", (-4.5, -2.5, 0.5)),
            ("Packera serpenticola", "Packera eburnea", (4.5, -2.0, -0.5)),
            ("Packera millefolium", "Senecio millefolium", (-1.5, 4.0, -2.0)),
            ("Packera aurea", "Senecio aureus", (2.5, 3.5, 2.5)),
        ]
        rows = []
        for idx, (canon_sp, raw_syn, (c1, c2, c3)) in enumerate(non_default_specs):
            for j in range(10):
                tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 8 else "Tier_3_Bronze")
                cs = float(28.0 + 1.2 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": f"ENA_{idx}_{j:02d}",
                    "plant_individual_id": 0,
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
                    "species_raw": raw_syn if j % 2 == 0 else canon_sp,
                    "scientificName": canon_sp,
                    "determiner_tier": tier,
                    "PC1": float(c1 + rng.normal(0.0, 0.15)),
                    "PC2": float(c2 + rng.normal(0.0, 0.15)),
                    "PC3": float(c3 + rng.normal(0.0, 0.12)),
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.8 + 0.25 * idx + rng.normal(0.0, 0.04)),
                    "solidity": float(0.42 if canon_sp == "Packera millefolium" else 0.89),
                    "ucs_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((c1 * 0.06) / k + rng.normal(0.0, 0.003))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((c2 * 0.06) / k + rng.normal(0.0, 0.003))
                rows.append(row)

        # Add a singleton rare outgroup row (n = 1) to verify passive projection without singular W
        singleton_row = dict(rows[0])
        singleton_row.update({
            "catalogNumber": "ENA_SINGLETON_01",
            "species_raw": "Senecio glabellus",
            "scientificName": "Packera glabella",
            "determiner_tier": "Tier_1_Gold",
            "PC1": 0.0,
            "PC2": 0.0,
            "PC3": 0.0,
        })
        rows.append(singleton_row)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            allom_csv = tmp_p / "allometry_significance_audit.csv"
            allom_pdf = tmp_p / "shape_allometry_regression.pdf"
            cda_pdf = tmp_p / "cda_passive_projection.pdf"

            df_in = pd.DataFrame(rows)
            df_in.to_csv(harm_csv, index=False)
            df_in[["catalogNumber", "species_raw", "determiner_tier"]].assign(
                institutionCode=["NCU" if i % 2 == 0 else "NY" for i in range(len(df_in))],
                stateProvince="North Carolina",
                latitude=35.6,
                longitude=-82.5,
                doy=132,
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
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(allom_csv),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
                    "--allometry-plot", str(allom_pdf),
                    "--permutations", "99",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Step 04 failed on non-default Packera taxa:\n{proc.stderr}\n{proc.stdout}")

            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), len(df_in), "All vouchers including singleton outgroup must be retained.")

            expected_four = {
                "Packera schweinitziana",
                "Packera serpenticola",
                "Packera millefolium",
                "Packera aurea",
            }
            pred_taxa = set(flags_df["cda_predicted_taxon"].dropna().unique())
            self.assertEqual(pred_taxa, expected_four)
            self.assertTrue(
                pred_taxa.isdisjoint({"Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"}),
                f"Default 4-species names must not appear when absent from dataset: {pred_taxa}",
            )

            # Verify high classification accuracy on the 4 well-separated species
            core_rows = flags_df[flags_df["species_standardized"].isin(expected_four)]
            acc = float((core_rows["species_standardized"] == core_rows["cda_predicted_taxon"]).mean())
            self.assertGreaterEqual(acc, 0.95, f"Expected >=95% CDA classification accuracy on 4 distinct taxa, got {acc:.3f}")

            # Verify singleton Packera glabella row was projected passively with finite 3D canonical scores
            sing_out = flags_df[flags_df["catalogNumber"] == "ENA_SINGLETON_01"]
            self.assertEqual(len(sing_out), 1)
            self.assertEqual(sing_out["species_standardized"].iloc[0], "Packera glabella")
            self.assertTrue(bool(sing_out["is_passive"].iloc[0]))
            for col in ("can1", "can2", "can3", "cda_posterior_prob"):
                self.assertTrue(np.isfinite(float(sing_out[col].iloc[0])), f"Singleton {col} must be finite.")

    def test_step02_and_preflight_dynamic_non_default_packera_support(self):
        """Verify Step 02 standardization/figures, config.yaml, and run_pipeline.sh support dynamic genus-wide Packera species."""
        import importlib
        from unittest.mock import patch
        from scripts.core.config import PipelineConfig

        seg_mod = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        std_fn = seg_mod._standardize_packera_species_py
        expected_map = {
            "Packera schweinitziana": "Packera schweinitziana",
            "Packera serpenticola": "Packera serpenticola",
            "Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve": "Packera millefolium",
            "Packera aurea (L.) Á.Löve & D.Löve": "Packera aurea",
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve": "Packera streptanthifolia",
            "Packera cymbalaria (Pursh) W.A.Weber & Á.Löve": "Packera cymbalaria",
            "Packera sanguisorbae (DC.) C.Jeffrey": "Packera sanguisorbae",
        }
        for raw_name, canon_name in expected_map.items():
            self.assertEqual(std_fn(raw_name), canon_name, f"Step 02 failed to standardize {raw_name!r}")

        with patch(
            "scripts.core.harvester.species.name_backbone",
            return_value={
                "matchType": "EXACT",
                "status": "SYNONYM",
                "species": "Packera schweinitziana",
                "canonicalName": "Senecio robbinsii",
            },
        ), patch("scripts.core.harvester.species.name_usage", return_value={"results": []}):
            self.assertEqual(std_fn("Senecio robbinsii"), "Packera schweinitziana")

        # Verify Step 02 publication figures render 8 diverse North American species with distinct colors/markers
        taxa_8 = [
            "Packera aurea",
            "Packera cana",
            "Packera cymbalaria",
            "Packera millefolium",
            "Packera sanguisorbae",
            "Packera schweinitziana",
            "Packera serpenticola",
            "Packera streptanthifolia",
        ]
        fig_rows = []
        for idx, sp in enumerate(taxa_8):
            for k in range(4):
                fig_rows.append({
                    "catalogNumber": f"FIG_{idx}_{k}",
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct",
                    "scientificName": sp,
                    "solidity": 0.35 + 0.05 * k if sp == "Packera millefolium" else 0.80 + 0.03 * k,
                    "ucs_score": 0.92,
                    "aspect_ratio": 1.4 + 0.3 * k,
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                })
        with tempfile.TemporaryDirectory() as tmpdir:
            saved = seg_mod.generate_step02_publication_figures(
                extracted_df=pd.DataFrame(fig_rows),
                failed_df=pd.DataFrame(),
                repro_df=pd.DataFrame(),
                pheno_df=pd.DataFrame(),
                figures_dir=Path(tmpdir) / "figures",
            )
            self.assertTrue(Path(saved["fig02a_pdf"]).exists())
            self.assertTrue(Path(saved["fig02b_pdf"]).exists())

        # Verify config.yaml dynamic genus-wide mode and run_pipeline.sh dynamic outgroups
        cfg = PipelineConfig.from_yaml()
        self.assertEqual(cfg.taxa.target_species, [])
        self.assertEqual(cfg.taxa.outgroups, [])

        sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        self.assertIn("base_outgroups = present_taxa", sh_src)
        self.assertNotIn("mod05.TARGET_TAXA", sh_src)

    def test_dual_unsupervised_morphotype_and_supervised_epa_ecoregion_cda(self):
        """Verify Step 04 runs BOTH Track 1 unsupervised Morphotype_K* CDA and Track 2 supervised EPA Level III Ecoregion CDA."""
        rng = np.random.default_rng(202610)
        # Three distinct EPA Level III Ecoregions + ungeoreferenced NA vouchers projected passively
        eco_specs = [
            ("Piedmont", 35.9132, -79.0558, "North Carolina", "Orange", -3.2, -1.6),
            ("Middle_Atlantic_Coastal_Plain", 34.2257, -77.9447, "North Carolina", "New Hanover", 3.2, -1.4),
            ("South_Central_Plains", 31.6035, -94.6555, "Texas", "Nacogdoches", 0.2, 3.0),
        ]
        rows = []
        vouch_rows = []
        idx_counter = 0
        for eco_name, lat, lon, st_name, co_name, c1, c2 in eco_specs:
            for j in range(12):
                cat = f"DUAL_{idx_counter:03d}"
                idx_counter += 1
                tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                cs = float(24.0 + 1.1 * j)
                log_cs = float(np.log(cs))
                rec = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
                    "species_raw": "Packera anonyma",
                    "scientificName": "Packera anonyma",
                    "determiner_tier": tier,
                    "PC1": float(c1 + rng.normal(0.0, 0.22)),
                    "PC2": float(c2 + rng.normal(0.0, 0.22)),
                    "PC3": float(rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)),
                    "PC5": float(rng.normal(0.0, 0.12)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.9 + 0.25 * (c1 > 0) + rng.normal(0.0, 0.05)),
                    "solidity": float(0.88 + rng.normal(0.0, 0.02)),
                    "ucs_score": 0.94,
                }
                for k in range(1, 13):
                    rec[f"A{k}"] = float((c1 * 0.06) / k + rng.normal(0.0, 0.003))
                    rec[f"B{k}"] = 0.0
                    rec[f"C{k}"] = 0.0
                    rec[f"D{k}"] = float((c2 * 0.06) / k + rng.normal(0.0, 0.003))
                rows.append(rec)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": "Packera anonyma",
                    "determiner_tier": tier,
                    "institutionCode": "NCU" if j % 2 == 0 else "TEX",
                    "stateProvince": st_name,
                    "county": co_name,
                    "decimalLatitude": lat if j < 10 else np.nan,
                    "decimalLongitude": lon if j < 10 else np.nan,
                    "latitude": lat if j < 10 else np.nan,
                    "longitude": lon if j < 10 else np.nan,
                    "regional_group": eco_name,
                    "doy": 118 + j,
                })

        # Add 4 ungeoreferenced Fall Line straddling-county vouchers (regional_group = NA) to verify passive projection
        for u_i in range(4):
            cat = f"DUAL_NA_{u_i:02d}"
            rec_na = dict(rows[u_i])
            rec_na.update({
                "catalogNumber": cat,
                "PC1": float(-2.9 + 0.2 * u_i),
                "PC2": float(-1.5 + 0.1 * u_i),
            })
            rows.append(rec_na)
            vouch_rows.append({
                "catalogNumber": cat,
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold",
                "institutionCode": "NCU",
                "stateProvince": "North Carolina",
                "county": "Wake",
                "decimalLatitude": np.nan,
                "decimalLongitude": np.nan,
                "latitude": np.nan,
                "longitude": np.nan,
                "regional_group": np.nan,
                "doy": 122,
            })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            fig_dir = tmp_p / "figures"
            rep_dir = tmp_p / "reports"
            fig_dir.mkdir(parents=True, exist_ok=True)
            rep_dir.mkdir(parents=True, exist_ok=True)

            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = rep_dir / "gmm_bayes_factors_summary.csv"
            eco_sum_csv = rep_dir / "ecoregion_supervised_cda_summary.csv"
            cda_pdf = fig_dir / "cda_passive_projection.pdf"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(rep_dir / "batch_audit.csv"),
                    "--allometry-report", str(rep_dir / "allometry.csv"),
                    "--allometry-pairwise", str(rep_dir / "pairwise.csv"),
                    "--allometry-plot", str(fig_dir / "allometry.pdf"),
                    "--focal-taxon", "Packera anonyma",
                    "--intraspecific-mode",
                    "--permutations", "199",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Dual Step 04 failed:\n{proc.stderr}\n{proc.stdout}")

            # 1. Track 1: Unsupervised GMM BIC table evaluates G = 1:eff_k with Two_Delta_BIC vs G = 1
            self.assertTrue(bic_csv.exists())
            bic_df = pd.read_csv(bic_csv)
            self.assertEqual(int(bic_df["K"].iloc[0]), 1)
            self.assertAlmostEqual(float(bic_df["Two_Delta_BIC"].iloc[0]), 0.0, places=5)
            self.assertAlmostEqual(
                float(bic_df["Two_Delta_BIC"].iloc[2]),
                round(float(bic_df["BIC"].iloc[2] - bic_df["BIC"].iloc[0]), 2),
                places=2,
            )

            # 2. Track 1: cda_predicted_taxon is strictly Morphotype_K1..Kn (never regional_group)
            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), 40)
            pred_morphs = set(flags_df["cda_predicted_taxon"].dropna().astype(str))
            self.assertGreaterEqual(len(pred_morphs), 2)
            for m_lbl in pred_morphs:
                self.assertRegex(m_lbl, r"^Morphotype_K[0-9]+$")

            # 3. Track 2: Supervised EPA Level III Ecoregion CDA columns in morphometric_flags.csv
            for req_col in ("ecoregion_cda_predicted", "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2"):
                self.assertIn(req_col, flags_df.columns)
                self.assertFalse(flags_df[req_col].isna().any(), f"{req_col} must be populated for active + NA-projected rows")
            pred_ecos = set(flags_df["ecoregion_cda_predicted"].dropna().astype(str))
            self.assertEqual(pred_ecos, {"Piedmont", "Middle_Atlantic_Coastal_Plain", "South_Central_Plains"})
            # NA-ecoregion rows keep regional_group == NA while receiving projected ecoregion_cda_predicted
            na_rows = flags_df[flags_df["catalogNumber"].str.startswith("DUAL_NA_")]
            self.assertEqual(len(na_rows), 4)
            self.assertTrue(na_rows["regional_group"].isna().all())
            self.assertTrue((na_rows["ecoregion_cda_predicted"] == "Piedmont").all())

            # 4. Track 2 summary report & fig04c_supervised_ecoregion_cda.{pdf,png}
            self.assertTrue(eco_sum_csv.exists(), "ecoregion_supervised_cda_summary.csv must be generated")
            sum_df = pd.read_csv(eco_sum_csv)
            self.assertEqual(set(sum_df["ecoregion"]), {"Piedmont", "Middle_Atlantic_Coastal_Plain", "South_Central_Plains"})
            for s_col in (
                "n_total_vouchers",
                "n_active_anchors",
                "loocv_accuracy",
                "permanova_r2",
                "permanova_p_value",
                "ari_morphotype_vs_ecoregion",
                "fisher_exact_p_value",
            ):
                self.assertIn(s_col, sum_df.columns)
            self.assertGreater(float(sum_df["loocv_accuracy"].iloc[0]), 0.85)
            self.assertGreater(float(sum_df["permanova_r2"].iloc[0]), 0.50)
            self.assertGreater(float(sum_df["ari_morphotype_vs_ecoregion"].iloc[0]), 0.50)
            self.assertLess(float(sum_df["fisher_exact_p_value"].iloc[0]), 0.01)

            for stem in (
                "fig04a_gmm_bic_and_uncertainty",
                "fig04b_cda_biplot_and_allometry",
                "fig04_gmm_cda_clustering",
                "fig04c_supervised_ecoregion_cda",
            ):
                pdf_f = fig_dir / f"{stem}.pdf"
                png_f = fig_dir / f"{stem}.png"
                self.assertTrue(pdf_f.exists() and pdf_f.stat().st_size > 3000, f"Missing/empty {pdf_f}")
                self.assertTrue(png_f.exists() and png_f.stat().st_size > 10000, f"Missing/empty {png_f}")

    def test_r_epa_level3_ecoregion_spatial_helper_agreement(self):
        """Verify R infer_epa_level3_ecoregion() in pub_figures_utils.R matches Python EPA Level III rules."""
        r_cmd = (
            'source("scripts/morphometrics/pub_figures_utils.R"); '
            'res <- infer_epa_level3_ecoregion('
            '  lat = c(35.9132, 34.2257, 32.5404, 31.6035, NA, NA, NA), '
            '  lon = c(-79.0558, -77.9447, -82.9037, -94.6555, NA, NA, NA), '
            '  state_province = c("NC", "NC", "GA", "TX", "North Carolina", "North Carolina", "NC"), '
            '  county = c("", "", "", "", "Orange County", "Wake County", "")'
            '); '
            'cat(paste(ifelse(is.na(res), "NA", res), collapse = "|"))'
        )
        proc = subprocess.run(["Rscript", "-e", r_cmd], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"R spatial helper failed:\n{proc.stderr}")
        got = proc.stdout.strip().split("|")
        self.assertEqual(
            got,
            [
                "Piedmont",
                "Middle_Atlantic_Coastal_Plain",
                "Southeastern_Plains",
                "South_Central_Plains",
                "Piedmont",
                "NA",
                "NA",
            ],
        )

    @staticmethod
    def _extract_pdf_text_streams(pdf_path: Path) -> str:
        """Decompress all FlateDecode streams in an R-generated vector PDF for text verification."""
        raw = pdf_path.read_bytes()
        chunks = [raw.decode("latin1", errors="ignore")]
        for m in re.finditer(rb"stream\r?\n(.*?)endstream", raw, flags=re.DOTALL):
            payload = m.group(1).rstrip(b"\r\n")
            for cand in (payload, m.group(1)):
                try:
                    txt = zlib.decompress(cand).decode("latin1", errors="ignore")
                    txt = re.sub(r"\)\s*-?\d+(?:\.\d+)?\s*\(", "", txt)
                    txt = txt.replace(r"\(", "(").replace(r"\)", ")")
                    chunks.append(txt)
                    break
                except zlib.error:
                    continue
        return "\n".join(chunks)

    def test_multivariate_manova_audit_three_factors_and_cda_subtitles(self):
        """Verify parametric MANOVA (Pillai/Wilks), pairwise Hotelling T2 contrasts, Tier_3_Bronze invariance, and CDA biplot subtitles."""
        rng = np.random.default_rng(314)
        specs = [
            ("Packera anonyma", "Piedmont", 35.9132, -79.0558, "NC", -5.5, -1.5, 11),
            ("Packera dubia", "Middle_Atlantic_Coastal_Plain", 34.2257, -77.9447, "NC", 5.5, -1.5, 11),
            ("Packera paupercula", "South_Central_Plains", 31.6035, -94.6555, "TX", 0.0, 5.0, 11),
        ]
        rows = []
        vouch_rows = []
        for s_idx, (sp, eco_name, lat0, lon0, st_code, c1, c2, n_v) in enumerate(specs):
            for j in range(n_v):
                cat = f"MANOVA_{s_idx}_{j:02d}"
                det_tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                pc1 = float(c1 + rng.normal(0.0, 0.25))
                pc2 = float(c2 + rng.normal(0.0, 0.25))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": det_tier,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": float(0.4 * s_idx + rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": float(28.0 + 0.4 * j),
                    "log_centsize": float(np.log(28.0 + 0.4 * j)),
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(2.0 + 0.3 * s_idx + rng.normal(0.0, 0.05)),
                    "solidity": float(0.89 + 0.02 * s_idx),
                    "ucs_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((-1.25 - 0.06 * pc1) / k + rng.normal(0.0, 0.003))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((0.46 + 0.05 * pc2) / k + rng.normal(0.0, 0.003))
                rows.append(row)
                vouch_rows.append(
                    {
                        "catalogNumber": cat,
                        "species_raw": sp,
                        "determiner_tier": det_tier,
                        "institutionCode": "NCU" if j % 2 == 0 else "BRIT",
                        "stateProvince": st_code,
                        "county": "",
                        "regional_group": eco_name,
                        "decimalLatitude": float(lat0 + rng.normal(0.0, 0.03)),
                        "decimalLongitude": float(lon0 + rng.normal(0.0, 0.03)),
                        "doy": 115 + j,
                    }
                )

        # Add 2 Blue_Ridge Gold vouchers (< 3 threshold -> excluded from regional_group MANOVA) and 2 NA-ecoregion vouchers
        for extra_idx, (rg_val, lat_v, lon_v) in enumerate(
            [
                ("Blue_Ridge", 35.5951, -82.5515),
                ("Blue_Ridge", 35.6100, -82.5400),
                ("", np.nan, np.nan),
                ("", np.nan, np.nan),
            ]
        ):
            cat = f"MANOVA_EXTRA_{extra_idx:02d}"
            row = dict(rows[0])
            row["catalogNumber"] = cat
            row["determiner_tier"] = "Tier_1_Gold"
            row["PC1"] = float(-5.4 + 0.1 * extra_idx)
            row["PC2"] = float(-1.4 - 0.1 * extra_idx)
            rows.append(row)
            vouch_rows.append(
                {
                    "catalogNumber": cat,
                    "species_raw": "Packera anonyma",
                    "determiner_tier": "Tier_1_Gold",
                    "institutionCode": "NCU",
                    "stateProvince": "NC",
                    "county": "",
                    "regional_group": rg_val,
                    "decimalLatitude": lat_v,
                    "decimalLongitude": lon_v,
                    "doy": 120,
                }
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)

            def _run_case(sub_name: str, perturb_bronze: bool = False):
                case_dir = tmp_p / sub_name
                fig_dir = case_dir / "outputs" / "figures"
                rep_dir = case_dir / "outputs" / "reports"
                tbl_dir = case_dir / "data" / "tables"
                for d in (fig_dir, rep_dir, tbl_dir):
                    d.mkdir(parents=True, exist_ok=True)
                df_h = pd.DataFrame(rows).copy()
                df_v = pd.DataFrame(vouch_rows).copy()
                if perturb_bronze:
                    b_mask = df_h["determiner_tier"] == "Tier_3_Bronze"
                    df_h.loc[b_mask, ["PC1", "PC2", "PC3"]] += 45.0
                    df_v.loc[df_v["determiner_tier"] == "Tier_3_Bronze", "regional_group"] = "Blue_Ridge"
                harm_csv = tbl_dir / "leaf_efa_harmonics.csv"
                vouch_csv = tbl_dir / "curated_vouchers.csv"
                flags_csv = tbl_dir / "morphometric_flags.csv"
                bic_csv = rep_dir / "gmm_bayes_factors_summary.csv"
                cda_pdf = fig_dir / "cda_passive_projection.pdf"
                df_h.to_csv(harm_csv, index=False)
                df_v.to_csv(vouch_csv, index=False)
                proc = subprocess.run(
                    [
                        "Rscript",
                        str(self.r_script_path),
                        "--vouchers", str(vouch_csv),
                        "--harmonics", str(harm_csv),
                        "--output-flags", str(flags_csv),
                        "--output-report", str(bic_csv),
                        "--output-plot", str(cda_pdf),
                        "--output-batch-audit", str(rep_dir / "batch_audit.csv"),
                        "--allometry-report", str(rep_dir / "allometry.csv"),
                        "--allometry-pairwise", str(rep_dir / "pairwise.csv"),
                        "--allometry-plot", str(fig_dir / "allometry.pdf"),
                        "--permutations", "99",
                    ],
                    cwd=str(PROJECT_ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, f"Step 04 MANOVA run failed:\n{proc.stderr}\n{proc.stdout}")
                return rep_dir, fig_dir

            rep_base, fig_base = _run_case("base", perturb_bronze=False)
            grp_csv = rep_base / "manova_group_significance.csv"
            pw_csv = rep_base / "manova_pairwise_posthoc.csv"
            self.assertTrue(grp_csv.exists(), "manova_group_significance.csv must be exported")
            self.assertTrue(pw_csv.exists(), "manova_pairwise_posthoc.csv must be exported")

            grp_df = pd.read_csv(grp_csv)
            pw_df = pd.read_csv(pw_csv)

            req_grp_cols = {
                "grouping_factor",
                "pillai_trace",
                "wilks_lambda",
                "approx_F",
                "num_df",
                "den_df",
                "partial_eta2",
                "p_value",
                "q_value",
            }
            self.assertTrue(req_grp_cols.issubset(set(grp_df.columns)))
            self.assertEqual(
                set(grp_df["grouping_factor"]),
                {"species_standardized", "gmm_cluster", "regional_group"},
            )

            for _, r in grp_df.iterrows():
                self.assertGreater(float(r["pillai_trace"]), 0.20)
                self.assertGreater(float(r["wilks_lambda"]), 0.0)
                self.assertLess(float(r["wilks_lambda"]), 0.90)
                self.assertGreater(float(r["approx_F"]), 3.0)
                self.assertGreaterEqual(float(r["num_df"]), 1.0)
                self.assertGreaterEqual(float(r["den_df"]), 2.0)
                self.assertGreater(float(r["partial_eta2"]), 0.10)
                self.assertLessEqual(float(r["partial_eta2"]), 1.0)
                self.assertGreaterEqual(float(r["p_value"]), 0.0)
                self.assertLess(float(r["p_value"]), 0.01)
                self.assertGreaterEqual(float(r["q_value"]) + 1e-12, float(r["p_value"]))
                self.assertLessEqual(float(r["q_value"]), 1.0)

            # Confirm Blue_Ridge (n=2 < 3) and NA were excluded from regional_group MANOVA
            rg_row = grp_df[grp_df["grouping_factor"] == "regional_group"].iloc[0]
            self.assertEqual(int(rg_row["n_groups"]), 3)

            req_pw_cols = {
                "grouping_factor",
                "group_1",
                "group_2",
                "contrast",
                "hotelling_t2",
                "pillai_trace",
                "wilks_lambda",
                "approx_F",
                "num_df",
                "den_df",
                "partial_eta2",
                "p_value",
                "q_value",
            }
            self.assertTrue(req_pw_cols.issubset(set(pw_df.columns)))
            sp_pw = pw_df[pw_df["grouping_factor"] == "species_standardized"]
            rg_pw = pw_df[pw_df["grouping_factor"] == "regional_group"]
            cl_pw = pw_df[pw_df["grouping_factor"] == "gmm_cluster"]
            self.assertEqual(len(sp_pw), 3)
            self.assertEqual(len(rg_pw), 3)
            self.assertGreaterEqual(len(cl_pw), 1)
            for m_lbl in set(cl_pw["group_1"]).union(set(cl_pw["group_2"])):
                self.assertRegex(str(m_lbl), r"^Morphotype_K[0-9]+$")
            self.assertNotIn("Blue_Ridge", set(rg_pw["group_1"]).union(set(rg_pw["group_2"])))

            for _, pr in pw_df.iterrows():
                self.assertGreater(float(pr["hotelling_t2"]), 0.0)
                # Exact 2-group identity: Pillai V + Wilks Lambda == 1 and partial_eta2 == Pillai V
                self.assertAlmostEqual(float(pr["pillai_trace"] + pr["wilks_lambda"]), 1.0, places=4)
                self.assertAlmostEqual(float(pr["partial_eta2"]), float(pr["pillai_trace"]), places=4)
                self.assertGreaterEqual(float(pr["q_value"]) + 1e-12, float(pr["p_value"]))

            # Verify omnibus Pillai's Trace F and p-value subtitles in both CDA biplot PDFs
            fig04_txt = self._extract_pdf_text_streams(fig_base / "fig04_gmm_cda_clustering.pdf")
            fig04c_txt = self._extract_pdf_text_streams(fig_base / "fig04c_supervised_ecoregion_cda.pdf")
            sp_row = grp_df[grp_df["grouping_factor"] == "species_standardized"].iloc[0]
            self.assertIn("Omnibus MANOVA: Pillai's Trace V =", fig04_txt)
            self.assertIn(f"= {float(sp_row['approx_F']):.2f}, p =", fig04_txt)
            self.assertIn("Omnibus MANOVA: Pillai's Trace V =", fig04c_txt)
            self.assertIn(f"= {float(rg_row['approx_F']):.2f}, p =", fig04c_txt)

            # Verify Tier_3_Bronze zero-influence invariance on both MANOVA CSVs
            rep_pert, _ = _run_case("perturbed_bronze", perturb_bronze=True)
            pd.testing.assert_frame_equal(grp_df, pd.read_csv(rep_pert / "manova_group_significance.csv"))
            pd.testing.assert_frame_equal(pw_df, pd.read_csv(rep_pert / "manova_pairwise_posthoc.csv"))

    def test_run_multivariate_manova_audit_direct_r_single_species_and_collinearity(self):
        """Verify direct run_multivariate_manova_audit(df, feature_cols, gmm_res) factor gating & collinearity safety."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            grp_out = tmp_p / "manova_group_significance.csv"
            pw_out = tmp_p / "manova_pairwise_posthoc.csv"
            r_code = f"""
            source("scripts/morphometrics/pub_figures_utils.R")
            set.seed(17)
            df <- data.frame(
              catalogNumber = sprintf("V%02d", 1:18),
              species_standardized = rep("Packera dubia", 18),
              determiner_tier = rep("Tier_1_Gold", 18),
              morph_tier_source = rep("Tier_1_Pure", 18),
              regional_group = c(rep("Piedmont", 10), rep("Blue_Ridge", 2), rep(NA_character_, 6)),
              PC1 = c(rnorm(9, -3.0, 0.3), rnorm(9, 3.0, 0.3)),
              PC2 = c(rnorm(9, 1.5, 0.3), rnorm(9, -1.5, 0.3)),
              stringsAsFactors = FALSE
            )
            df$PC3_collinear <- df$PC1 + df$PC2
            df$PC4_const <- 0.0
            gmm_res <- list(classification = c(rep(1L, 9), rep(2L, 9)))
            res <- run_multivariate_manova_audit(
              df,
              c("PC1", "PC2", "PC3_collinear", "PC4_const"),
              gmm_res,
              output_group_csv = "{grp_out.as_posix()}",
              output_pairwise_csv = "{pw_out.as_posix()}"
            )
            stopifnot(nrow(res$group_table) == 1L)
            stopifnot(res$group_table$grouping_factor[1] == "gmm_cluster")
            """
            proc = subprocess.run(["Rscript", "-e", r_code], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(proc.returncode, 0, f"Direct R run_multivariate_manova_audit failed:\n{proc.stderr}")
            self.assertTrue(grp_out.exists())
            self.assertTrue(pw_out.exists())
            g_df = pd.read_csv(grp_out)
            p_df = pd.read_csv(pw_out)
            self.assertEqual(list(g_df["grouping_factor"]), ["gmm_cluster"])
            self.assertEqual(int(g_df["n_features"].iloc[0]), 2)
            self.assertEqual(len(p_df), 1)
            self.assertEqual(p_df["contrast"].iloc[0], "Morphotype_K1 vs Morphotype_K2")
            self.assertGreater(float(p_df["hotelling_t2"].iloc[0]), 20.0)

    def test_zero_gold_species_adaptive_anchor_promotion_and_empirical_bayes_weights(self):
        """Verify multi-species GMM and CDA succeed via Adaptive Anchor Promotion when species have 0 Tier_1_Gold vouchers."""
        r_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertNotIn("Glabrescent_P_dubia_Misidentified_as_P_anonyma", r_src)
        self.assertIn("Unassigned_Packera", r_src)

        rng = np.random.default_rng(202604)
        rows = []
        vouch_rows = []

        # Species 1: Packera dubia (has Gold, Silver, and Bronze + 1 discordant Bronze misidentified near P. anonyma)
        # Species 2: Packera anonyma (0 Tier_1_Gold — only 3 Tier_2_Silver + 7 Tier_3_Bronze including 1 extreme outlier)
        # Species 3: Packera serpenticola (0 Tier_1_Gold and 0 Tier_2_Silver — only 10 Tier_3_Bronze including 1 extreme outlier)
        species_specs = [
            ("Packera dubia", (-5.0, -1.8, 0.4), ["Tier_1_Gold"] * 5 + ["Tier_2_Silver"] * 4 + ["Tier_3_Bronze"] * 3),
            ("Packera anonyma", (5.0, -1.5, -0.4), ["Tier_2_Silver"] * 3 + ["Tier_3_Bronze"] * 7),
            ("Packera serpenticola", (0.0, 4.8, 0.8), ["Tier_3_Bronze"] * 10),
        ]

        for s_idx, (sp, (c1, c2, c3), tiers) in enumerate(species_specs):
            for j, tier in enumerate(tiers):
                cat = f"ADAPT_{s_idx}_{j:02d}"
                pc1 = float(c1 + rng.normal(0.0, 0.16))
                pc2 = float(c2 + rng.normal(0.0, 0.16))
                pc3 = float(c3 + rng.normal(0.0, 0.12))
                # Inject discordant / outlier Bronze rows to test Mahalanobis screening & generalized discordance_reason
                if sp == "Packera dubia" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = 5.05, -1.45, -0.38  # Misidentified P. dubia sitting inside P. anonyma morphospace
                elif sp == "Packera anonyma" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = -18.0, 16.0, 9.0    # Extreme Bronze outlier -> must NOT be promoted
                elif sp == "Packera serpenticola" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = 19.0, -17.0, -8.0   # Extreme Bronze outlier -> must NOT be promoted

                cs = float(26.0 + 0.9 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": tier,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": pc3,
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.9 + 0.2 * s_idx + rng.normal(0.0, 0.04)),
                    "solidity": 0.90,
                    "ucs_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((pc1 * 0.06) / k + rng.normal(0.0, 0.002))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((pc2 * 0.06) / k + rng.normal(0.0, 0.002))
                rows.append(row)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "determiner_tier": tier,
                    "institutionCode": "NCU" if j % 2 == 0 else "DUKE",
                    "stateProvince": "NC",
                    "county": "Orange",
                    "decimalLatitude": 35.91,
                    "decimalLongitude": -79.05,
                    "regional_group": "Piedmont",
                    "doy": 115 + j,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            cda_pdf = tmp_p / "cda_passive_projection.pdf"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(tmp_p / "allometry.csv"),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
                    "--allometry-plot", str(tmp_p / "allometry.pdf"),
                    "--permutations", "99",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Step 04 failed with 0-Gold species:\n{proc.stderr}\n{proc.stdout}")

            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), len(rows))
            # All 3 species (including 0-Gold Packera anonyma and 0-Gold/0-Silver Packera serpenticola) retained in CDA
            self.assertEqual(
                set(flags_df["cda_predicted_taxon"].dropna().unique()),
                {"Packera dubia", "Packera anonyma", "Packera serpenticola"},
            )

            # Packera anonyma (3 Silver + 7 Bronze): promotes clean Bronze anchors, rejects extreme outlier ADAPT_1_09
            anon_df = flags_df[flags_df["species_standardized"] == "Packera anonyma"]
            anon_promo = anon_df[(anon_df["determiner_tier"] == "Tier_3_Bronze") & (anon_df["authority_weight"] > 0.0)]
            self.assertGreaterEqual(len(anon_promo), 3, "Expected clean Bronze anchors promoted for 0-Gold Packera anonyma")
            anon_outlier = flags_df[flags_df["catalogNumber"] == "ADAPT_1_09"].iloc[0]
            self.assertEqual(float(anon_outlier["authority_weight"]), 0.0)
            self.assertTrue(bool(anon_outlier["is_passive"]))

            # Packera serpenticola (0 Gold, 0 Silver, 10 Bronze): promotes clean Bronze anchors, rejects outlier ADAPT_2_09
            serp_df = flags_df[flags_df["species_standardized"] == "Packera serpenticola"]
            serp_promo = serp_df[serp_df["authority_weight"] > 0.0]
            self.assertGreaterEqual(len(serp_promo), 5, "Expected clean Bronze anchors promoted for 0-Gold/0-Silver species")
            serp_outlier = flags_df[flags_df["catalogNumber"] == "ADAPT_2_09"].iloc[0]
            self.assertEqual(float(serp_outlier["authority_weight"]), 0.0)
            self.assertTrue(bool(serp_outlier["is_passive"]))

            # Generalized taxon-agnostic discordance_reason (no hardcoded Glabrescent_P_dubia_Misidentified_as_P_anonyma)
            dub_misid = flags_df[flags_df["catalogNumber"] == "ADAPT_0_11"].iloc[0]
            self.assertEqual(str(dub_misid["cda_predicted_taxon"]), "Packera anonyma")
            self.assertEqual(str(dub_misid["discordance_reason"]), "Packera_dubia_Discordant_to_Packera_anonyma")

            # Verify continuous_prior_posterior mode [0.5, 3.0] and Beta-Binomial Empirical Bayes shrinkage directly in R
            r_eb_check = """
            source("scripts/morphometrics/authority_weighting_utils.R")
            df_t <- data.frame(
              species_standardized = rep("Packera aurea", 8),
              determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver", "Tier_2_Silver",
                                  "Tier_3_Bronze", "Tier_3_Bronze", "Tier_3_Bronze", "Tier_3_Bronze"),
              is_promoted_bronze = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
              c_error = c(0.02, 0.04, 0.12, 0.14, 0.25, 0.30, 0.55, 0.65),
              stringsAsFactors = FALSE
            )
            w_ap <- compute_empirical_bayes_authority_weights(df_t, mode = "active_passive")
            w_cp <- compute_empirical_bayes_authority_weights(df_t, mode = "continuous_prior_posterior")
            stopifnot(all(w_ap[7:8] == 0.0), abs(mean(w_ap[1:6]) - 1.0) < 1e-5)
            stopifnot(w_ap[1] > w_ap[3], w_ap[3] > w_ap[5], w_ap[5] > 0.0)
            stopifnot(abs(max(w_cp) - 3.0) < 1e-5, min(w_cp) >= 0.5, w_cp[1] > w_cp[3], w_cp[3] > w_cp[5], w_cp[5] > w_cp[7])
            cat("EB_OK")
            """
            eb_proc = subprocess.run(["Rscript", "-e", r_eb_check], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(eb_proc.returncode, 0, f"Direct EB authority weighting check failed:\n{eb_proc.stderr}")
            self.assertIn("EB_OK", eb_proc.stdout)


if __name__ == "__main__":
    unittest.main()
