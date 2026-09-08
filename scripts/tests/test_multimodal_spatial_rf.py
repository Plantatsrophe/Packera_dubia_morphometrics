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


if __name__ == "__main__":
    unittest.main()
