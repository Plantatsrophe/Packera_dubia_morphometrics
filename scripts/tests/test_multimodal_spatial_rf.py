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

from scripts.analysis.run_multimodal_spatial_rf import (
    standardize_packera_taxon,
    extract_environmental_layers,
    execute_crossmodal_consensus,
    TARGET_TAXA,
)


class TestMultimodalSpatialRF(unittest.TestCase):
    """Test suite for Step 06 Spatial Random Forests and Niche Modeling."""

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
        """Verify that Step 06 R and Python entrypoints exist."""
        r_script = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
        py_script = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.py"
        self.assertTrue(r_script.exists(), "06_multimodal_spatial_rf.R must exist.")
        self.assertTrue(py_script.exists(), "06_multimodal_spatial_rf.py must exist.")


if __name__ == "__main__":
    unittest.main()
