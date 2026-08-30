#!/usr/bin/env python3
"""
===============================================================================
Script: test_triage_dashboard_synthesis.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and regression tests for Step 07:
    Multi-Evidence Taxonomic Decision Matrix, Triage Queue Generation,
    Priority Categorization, and Taxonomic Revision Report Synthesis.
===============================================================================
"""

import unittest
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts.analysis.run_triage_dashboard_synthesis import (
    standardize_packera_taxon,
    apply_taxonomic_decision_matrix,
    TARGET_TAXA,
)


class TestTriageDashboardSynthesis(unittest.TestCase):
    """Test suite for Step 07 Multi-Evidence Taxonomic Triage Synthesis."""

    def test_decision_matrix_concordant(self):
        """Verify that fully concordant vouchers receive RESOLVED priority."""
        df = pd.DataFrame([{
            "catalogNumber": "NCU001",
            "species_standardized": "Packera dubia",
            "species_raw": "Packera dubia",
            "determiner_tier": "Tier_1_Gold",
            "doy": 120,
            "cda_predicted_taxon": "Packera dubia",
            "cda_posterior_prob": 0.98,
            "vision_predicted_label": "Packera dubia",
            "c_error": 0.05,
            "is_label_corrupted": False,
            "misidentification_flag": False,
            "triage_category": "Clean_MultiModal_Consensus",
            "edaphic_best_fit_taxon": "Packera dubia",
            "soil_ph": 5.0,
            "soil_sand": 70.0,
            "gmm_uncertainty": 0.02,
        }])
        res_df = apply_taxonomic_decision_matrix(df)
        self.assertEqual(res_df.loc[0, "synthesis_triage_priority"], "RESOLVED")
        self.assertEqual(res_df.loc[0, "taxonomic_status_call"], "Species_Confirmed")
        self.assertIn("Accept_Current_Determination", res_df.loc[0, "synthesis_triage_action"])

    def test_decision_matrix_critical_conflict(self):
        """Verify that Tier 1 specialist vouchers with severe multi-stream discordance get Critical priority."""
        df = pd.DataFrame([{
            "catalogNumber": "NCU002",
            "species_standardized": "Packera dubia",
            "species_raw": "Packera dubia",
            "determiner_tier": "Tier_1_Gold",
            "doy": 120,
            "cda_predicted_taxon": "Packera anonyma",
            "cda_posterior_prob": 0.95,
            "vision_predicted_label": "Packera anonyma",
            "c_error": 0.95,
            "is_label_corrupted": True,
            "misidentification_flag": True,
            "triage_category": "Severe_Triple_Stream_Conflict",
            "edaphic_best_fit_taxon": "Packera anonyma",
            "soil_ph": 6.8,
            "soil_sand": 25.0,
            "gmm_uncertainty": 0.02,
        }])
        res_df = apply_taxonomic_decision_matrix(df)
        self.assertEqual(res_df.loc[0, "synthesis_triage_priority"], "CRITICAL")
        self.assertEqual(res_df.loc[0, "taxonomic_status_call"], "Misidentification_Severe")
        self.assertEqual(res_df.loc[0, "recommended_determination"], "Packera anonyma")

    def test_output_artifacts_exist(self):
        """Verify that synthesis report and figure artifacts exist in outputs/."""
        report_path = PROJECT_ROOT / "outputs" / "reports" / "Packera_dubia_Taxonomic_Revision_Summary.md"
        figure_pdf = PROJECT_ROOT / "outputs" / "figures" / "Figure_Integrative_Packera_dubia_Revision.pdf"
        triage_csv = PROJECT_ROOT / "data" / "tables" / "triage_queue.csv"

        self.assertTrue(report_path.exists(), f"Missing report at {report_path}")
        self.assertTrue(figure_pdf.exists(), f"Missing figure at {figure_pdf}")
        self.assertTrue(triage_csv.exists(), f"Missing triage queue at {triage_csv}")

    def test_entrypoint_scripts_exist(self):
        """Verify that Step 07 R and Python entrypoints exist."""
        r_script = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"
        py_script = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.py"
        self.assertTrue(r_script.exists(), "07_triage_dashboard_synthesis.R must exist.")
        self.assertTrue(py_script.exists(), "07_triage_dashboard_synthesis.py must exist.")


if __name__ == "__main__":
    unittest.main()
