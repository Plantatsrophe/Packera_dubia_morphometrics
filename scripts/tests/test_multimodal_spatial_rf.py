#!/usr/bin/env python3
"""
===============================================================================
Script: test_multimodal_spatial_rf.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and regression tests for Step 06:
    Multimodal Spatial Random Forests, SoilGrids 250m / WorldClim feature
    extraction, Moran's Eigenvector Maps (MEMs), and cross-modal consensus.
===============================================================================
"""

import unittest
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts._archive.analysis.run_multimodal_spatial_rf import (
    standardize_packera_taxon,
    extract_environmental_layers,
    execute_crossmodal_consensus,
    classify_ssurgo_outcrop,
    query_ssurgo_sda,
    synthesize_multiscale_edaphics,
    TARGET_TAXA,
)


class TestMultimodalSpatialRF(unittest.TestCase):
    """Test suite for Step 06 Spatial Random Forests, SSURGO Pedology, and Niche Modeling."""

    def test_taxonomic_standardization(self):
        """Verify standardization of raw species names to core taxa."""
        self.assertEqual(standardize_packera_taxon("Senecio tomentosus Michx."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera tomentosa (Michx.) C.Jeffrey"), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Senecio smallii Britton"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Packera anonyma (Alph.Wood) W.A.Weber"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. savannarum"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Senecio plattensis Nutt."), "Packera plattensis")

    def test_environmental_layer_extraction(self):
        """Verify extraction and imputation of edaphic and bioclimatic layers."""
        df = pd.DataFrame({
            "catalogNumber": ["TEST001", "TEST002"],
            "latitude": [35.5, 38.0],
            "longitude": [-79.0, -85.0],
            "regional_group": ["Southeastern_Coastal_Plain", "Interior_Prairie_Midwest"],
            "species_standardized": ["Packera dubia", "Packera anonyma"]
        })
        extracted_df = extract_environmental_layers(df)
        expected_cols = [
            "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density",
            "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality"
        ]
        for col in expected_cols:
            self.assertIn(col, extracted_df.columns, f"Missing environmental feature column: {col}")
            self.assertFalse(extracted_df[col].isna().any(), f"NaN values found in {col}")

    def test_classify_ssurgo_outcrop(self):
        """Verify micro-edaphic outcrop classification across pedologic categories."""
        # Granitic flatrock
        is_out, cls = classify_ssurgo_outcrop("Rock outcrop-Wake complex", "Rock outcrop", "Udepts", "Lithic Dystrudepts", 35.5, -80.0)
        self.assertTrue(is_out)
        self.assertEqual(cls, "Granitic flatrock")

        # Limestone glade
        is_out, cls = classify_ssurgo_outcrop("Gladeville-Rock outcrop complex", "Gladeville", "Udolls", "Lithic Hapludolls", 36.0, -86.5)
        self.assertTrue(is_out)
        self.assertEqual(cls, "Limestone glade")

        # Sandstone/Siliceous
        is_out, cls = classify_ssurgo_outcrop("Ramsey-Rock outcrop complex", "Ramsey", "Udepts", "Lithic Dystrudepts")
        self.assertTrue(is_out)
        self.assertEqual(cls, "Sandstone/Siliceous")

        # Serpentine
        is_out, cls = classify_ssurgo_outcrop("Chrome-Rock outcrop complex", "Chrome", "Udepts", "Lithic Eutrudepts")
        self.assertTrue(is_out)
        self.assertEqual(cls, "Serpentine")

        # Non-outcrop matrix
        is_out, cls = classify_ssurgo_outcrop("Cecil sandy clay loam", "Cecil", "Udults", "Typic Hapludults")
        self.assertFalse(is_out)
        self.assertEqual(cls, "Non-outcrop matrix")

    def test_ssurgo_query_and_caching(self):
        """Verify chunked SDA query and disk caching behavior."""
        df = pd.DataFrame({
            "catalogNumber": ["TEST_NCU_01", "TEST_NCU_02", "TEST_CAN_01"],
            "latitude": [35.5, 36.0, 55.0], # Two CONUS, one non-CONUS
            "longitude": [-80.5, -86.5, -110.0],
            "species_standardized": ["Packera anonyma", "Packera paupercula", "Packera paupercula"],
            "regional_group": ["Piedmont_Granite_Flatrocks", "Appalachian_Highlands", "Canada"]
        })
        cache_dir = PROJECT_ROOT / "data" / "environmental" / "test_ssurgo_cache"
        try:
            queried = query_ssurgo_sda(df, cache_dir=cache_dir, batch_size=2)
            self.assertIn("musym", queried.columns)
            self.assertIn("muname", queried.columns)
            self.assertIn("taxsuborder", queried.columns)
            self.assertIn("taxgreatgroup", queried.columns)
            self.assertIn("is_rock_outcrop", queried.columns)
            self.assertIn("outcrop_class", queried.columns)

            # CONUS records must have non-null pedology
            self.assertTrue(pd.notna(queried.loc[0, "musym"]))
            self.assertTrue(pd.notna(queried.loc[1, "musym"]))
            # Cache files should have been created
            cached_files = list(cache_dir.glob("*.csv"))
            self.assertGreater(len(cached_files), 0)
        finally:
            # Cleanup test cache
            if cache_dir.exists():
                for f in cache_dir.glob("*"):
                    f.unlink()
                cache_dir.rmdir()

    def test_synthesize_multiscale_edaphics(self):
        """Verify contingency table generation and Fisher's exact test calculation."""
        df = pd.DataFrame({
            "catalogNumber": [f"T{i}" for i in range(20)],
            "latitude": [35.0] * 20,
            "longitude": [-80.0] * 20,
            "species_standardized": ["Packera anonyma"] * 10 + ["Packera paupercula"] * 10,
            "regional_group": ["Piedmont_Granite_Flatrocks"] * 10 + ["Interior_Prairie_Midwest"] * 10,
            "musym": ["RoB"] * 8 + ["ApB"] * 2 + ["RoB"] * 1 + ["TaA"] * 9,
            "muname": ["Rock outcrop-Wake complex"] * 8 + ["Appling"] * 2 + ["Rock outcrop-Wake"] * 1 + ["Tama"] * 9,
            "taxsuborder": ["Udepts"] * 8 + ["Udults"] * 2 + ["Udepts"] * 1 + ["Udolls"] * 9,
            "taxgreatgroup": ["Lithic Dystrudepts"] * 8 + ["Typic Kanhapludults"] * 2 + ["Lithic Dystrudepts"] * 1 + ["Typic Argiudolls"] * 9,
            "is_rock_outcrop": [True] * 8 + [False] * 2 + [True] * 1 + [False] * 9,
            "outcrop_class": ["Granitic flatrock"] * 8 + ["Non-outcrop matrix"] * 2 + ["Granitic flatrock"] * 1 + ["Non-outcrop matrix"] * 9,
            "soil_ph": [5.2] * 10 + [7.1] * 10,
            "soil_cec": [12.0] * 10 + [25.0] * 10,
            "soil_sand": [55.0] * 10 + [22.0] * 10,
            "soil_bulk_density": [1.32] * 10 + [1.25] * 10,
        })
        out_table = PROJECT_ROOT / "data" / "tables" / "test_ssurgo_val.csv"
        out_rep = PROJECT_ROOT / "outputs" / "reports" / "test_contingency.csv"
        try:
            res = synthesize_multiscale_edaphics(df, ssurgo_out=out_table, contingency_out=out_rep)
            self.assertIn("contingency", res)
            self.assertIn("odds_ratio", res)
            self.assertGreater(res["odds_ratio"], 1.0)
            self.assertTrue(out_table.exists())
            self.assertTrue(out_rep.exists())
        finally:
            if out_table.exists():
                out_table.unlink()
            if out_rep.exists():
                out_rep.unlink()

    def test_crossmodal_consensus_flags(self):
        """Verify cross-modal consensus flag computation."""
        df = pd.DataFrame({
            "catalogNumber": ["T1", "T2"],
            "species_standardized": ["Packera dubia", "Packera dubia"],
            "determiner_tier": ["Tier 1 (Specialist)", "Tier 3 (Unverified)"],
            "cda_predicted_taxon": ["Packera dubia", "Packera anonyma"],
            "vision_predicted_label": ["Packera dubia", "Packera anonyma"],
            "soil_ph": [4.8, 7.8],
            "soil_sand": [65.0, 15.0],
            "soil_cec": [8.5, 24.0],
            "soil_bulk_density": [1.4, 1.2],
            "doy": [120, 120],
            "doy_sin": [0.5, 0.5],
            "doy_cos": [0.8, 0.8]
        })
        flagged_df = execute_crossmodal_consensus(df)
        self.assertIn("multimodal_conflict_flag", flagged_df.columns)
        self.assertIn("triage_category", flagged_df.columns)
        self.assertFalse(flagged_df.loc[0, "multimodal_conflict_flag"])
        self.assertEqual(flagged_df.loc[0, "triage_category"], "Clean_MultiModal_Consensus")

    def test_entrypoint_scripts_exist(self):
        """Verify that Step 06 R canonical script and archived Python entrypoint exist."""
        r_script = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
        py_archived = PROJECT_ROOT / "scripts" / "_archive" / "analysis" / "06_multimodal_spatial_rf.py"
        self.assertTrue(r_script.exists(), "06_multimodal_spatial_rf.R must exist as canonical R script.")
        self.assertTrue(py_archived.exists(), "06_multimodal_spatial_rf.py must exist in _archive/analysis.")

    def test_production_diagnostic_tables(self):
        """Verify production SSURGO validation and contingency report tables."""
        val_path = PROJECT_ROOT / "data" / "tables" / "ssurgo_edaphic_validation.csv"
        con_path = PROJECT_ROOT / "outputs" / "reports" / "micro_edaphic_outcrop_contingency.csv"
        self.assertTrue(val_path.exists(), "ssurgo_edaphic_validation.csv must exist.")
        self.assertTrue(con_path.exists(), "micro_edaphic_outcrop_contingency.csv must exist.")

        con_df = pd.read_csv(con_path)
        contrast_row = con_df[con_df["Taxon"].str.contains("Contrast")].iloc[0]
        # Odds ratio must reflect significant outcrop fidelity (OR > 1.0, p < 0.001)
        self.assertGreater(float(contrast_row["Odds_Ratio"]), 1.0)
        self.assertLess(float(contrast_row["P_Value"]), 0.001)

    def test_phenological_leap_year_doy_extraction(self):
        """Verify day-of-year calculation accurately accounts for leap years (366 days)."""
        from scripts.analysis.generate_phenological_artifacts import extract_doy
        test_df = pd.DataFrame({
            "eventDate": ["2020-02-29", "2020-03-01", "2021-02-28", "2021-03-01"],
            "year": [2020, 2020, 2021, 2021],
            "month": [2, 3, 2, 3],
            "day": [29, 1, 28, 1]
        })
        doys = extract_doy(test_df).tolist()
        self.assertEqual(doys[0], 60)  # Feb 29 in leap year
        self.assertEqual(doys[1], 61)  # Mar 01 in leap year
        self.assertEqual(doys[2], 59)  # Feb 28 in common year
        self.assertEqual(doys[3], 60)  # Mar 01 in common year

    def test_phenological_baseline_cline_and_hopkins_law(self):
        """Verify empirical baseline slope matches Hopkins' Bioclimatic Law (3.5 - 4.5 days/deg)."""
        from scripts.analysis.generate_phenological_artifacts import model_latitudinal_spring_baseline, extract_doy
        vouchers_path = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
        df = pd.read_csv(vouchers_path, low_memory=False)
        df["doy"] = extract_doy(df)
        df["decimalLatitude"] = pd.to_numeric(df["decimalLatitude"].combine_first(df["latitude"]), errors="coerce")
        baseline = model_latitudinal_spring_baseline(df)

        self.assertGreaterEqual(baseline["slope"], 3.5, "Baseline slope must be >= 3.5 days/degree")
        self.assertLessEqual(baseline["slope"], 4.5, "Baseline slope must be <= 4.5 days/degree")
        self.assertGreater(baseline["r_squared"], 0.40, "R2 must reflect strong latitudinal cline (> 0.40)")
        self.assertLess(baseline["p_value"], 1e-15, "Latitudinal cline must be highly significant")

    def test_phenological_anomalies_and_allochronic_divergence(self):
        """Verify phenological anomalies table and ANOVA / Tukey HSD allochronic separation."""
        anom_path = PROJECT_ROOT / "data" / "tables" / "phenological_anomalies.csv"
        sum_path = PROJECT_ROOT / "outputs" / "reports" / "phenological_anomaly_summary.csv"
        self.assertTrue(anom_path.exists(), "phenological_anomalies.csv must exist.")
        self.assertTrue(sum_path.exists(), "phenological_anomaly_summary.csv must exist.")

        anom_df = pd.read_csv(anom_path)
        self.assertIn("expected_doy", anom_df.columns)
        self.assertIn("delta_doy", anom_df.columns)
        self.assertIn("pheno_timing_category", anom_df.columns)
        valid_anom = anom_df.dropna(subset=["delta_doy"])
        self.assertGreater(len(valid_anom), 2000, "Must have > 2000 valid delta_doy vouchers.")

        sum_df = pd.read_csv(sum_path)
        anova_row = sum_df[sum_df["Analysis_Section"] == "One_Way_ANOVA"].iloc[0]
        f_stat = float(anova_row["Mean_Delta_DOY"])
        p_val = float(anova_row["P_Value"])
        self.assertGreater(f_stat, 100.0, "ANOVA F-statistic must indicate decisive divergence (F > 100)")
        self.assertLess(p_val, 1e-50, "ANOVA p-value must be < 1e-50")

        # Tukey HSD: P. anonyma vs P. dubia
        tukey_dub_ano = sum_df[sum_df["Taxon_or_Comparison"].str.contains("Packera anonyma vs Packera dubia")].iloc[0]
        gap_days = abs(float(tukey_dub_ano["Mean_Delta_DOY"]))
        self.assertGreater(gap_days, 15.0, "Allochronic separation between P. dubia and P. anonyma must be > 15 days")

    def test_phenological_diagnostic_figure_artifact(self):
        """Verify 2-panel publication figure PDF exists and has valid PDF format."""
        fig_path = PROJECT_ROOT / "outputs" / "figures" / "phenological_latitudinal_anomaly.pdf"
        self.assertTrue(fig_path.exists(), "phenological_latitudinal_anomaly.pdf must exist.")
        self.assertGreater(fig_path.stat().st_size, 10000, "Figure PDF must be > 10 KB.")
        with open(fig_path, "rb") as f:
            header = f.read(5)
        self.assertEqual(header, b"%PDF-", "File must have valid PDF magic bytes.")

    def test_r_script_phenological_integration(self):
        """Verify 06_multimodal_spatial_rf.R integrates phenological baseline and delta_doy."""
        r_script = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
        with open(r_script, "r", encoding="utf-8") as f:
            r_content = f.read()
        self.assertIn("--pheno-anomalies", r_content)
        self.assertIn("--pheno-summary", r_content)
        self.assertIn("--pheno-plot", r_content)
        self.assertIn("model_latitudinal_spring_baseline", r_content)
        self.assertIn("compute_phenological_anomalies", r_content)
        self.assertIn("export_phenological_latitudinal_figures", r_content)
        self.assertIn("delta_doy", r_content)


if __name__ == "__main__":
    unittest.main()

