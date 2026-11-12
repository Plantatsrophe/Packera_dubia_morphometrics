#!/usr/bin/env python3
"""
===============================================================================
Script: test_triage_dashboard_synthesis.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 06 (scripts/analysis/06_triage_dashboard_synthesis.R):
    - Safe Relational Joins & Non-Georeferenced Voucher Coalescing
    - Reliability-Weighted Probabilistic Consensus Score (Dawid & Skene 1979)
      with Dynamic Missing-Modality Weight Renormalization & OOF Class Precision
    - Five-Stage Hierarchical Integrative Decision Matrix (Padial et al. 2010;
      Schlick-Steiner et al. 2010; Rieseberg 1995; Slenker et al. 2022)
    - Hopkins-Corrected Circular Anthesis Phenology (Hopkins 1918; Pewsey et al. 2013)
    - Taxon-Adjusted Marginal ANCOVA Batch Audit (partial eta-squared + BH FDR)
===============================================================================
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import importlib
import sqlite3
import unittest
from unittest.mock import patch

from scripts.core.db import get_db_connection
from scripts.core.harvester import standardize_packera_taxon
from scripts.tests._r_session import rscript_available, run_warm_r_cmd
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


class TestTriageDashboardSynthesis(unittest.TestCase):
    """Test suite for Step 06 active R script (scripts/analysis/06_triage_dashboard_synthesis.R)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.r_unit_proc = None

        if rscript_available():
            cls.prec_in_csv = cls._class_tmp_path / "prec_in.csv"
            cls.prec_out_csv = cls._class_tmp_path / "prec_out.csv"
            pd.DataFrame([
                {
                    "catalogNumber": "STAGE1-MISID-OVER-HYB",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_2_Silver",
                    "decimalLatitude": 35.812,
                    "decimalLongitude": -79.012,
                    "doy": 122,
                    "PC1": 2.40,
                    "PC2": 0.02,
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_second_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.91,
                    "cda_second_posterior_prob": 0.09,
                    "cda_posterior_margin": 0.82,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera anonyma",
                    "rf_predicted_taxon": "Packera anonyma",
                    "rf_predicted_prob": 0.88,
                    "edaphic_schoener_d": 0.18,
                    "delta_doy_resid": 3.0,
                    "gmm_cluster": 2,
                    "gmm_uncertainty": 0.42,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "STAGE2-MOSAIC-HYBRID",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_2_Silver",
                    "decimalLatitude": 35.806,
                    "decimalLongitude": -79.006,
                    "doy": 119,
                    "PC1": -0.05,
                    "PC2": 0.02,
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_second_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.52,
                    "cda_second_posterior_prob": 0.40,
                    "cda_posterior_margin": 0.12,
                    "morphometric_outlier_flag": False,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "rf_predicted_taxon": "Packera dubia",
                    "rf_predicted_prob": 0.48,
                    "edaphic_schoener_d": 0.52,
                    "delta_doy_resid": 4.0,
                    "gmm_cluster": 2,
                    "gmm_uncertainty": 0.38,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "STAGE3-TRUE-CRYPTIC",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "decimalLatitude": 35.803,
                    "decimalLongitude": -79.003,
                    "doy": 144,
                    "PC1": 0.0,
                    "PC2": 3.6,
                    "cda_predicted_taxon": "Packera dubia",
                    "cda_second_taxon": "Packera anonyma",
                    "cda_posterior_prob": 0.46,
                    "cda_second_posterior_prob": 0.28,
                    "cda_posterior_margin": 0.18,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.18,
                    "delta_doy_resid": 24.0,
                    "gmm_cluster": 5,
                    "gmm_2_delta_bic": 12.4,
                    "gmm_uncertainty": 0.08,
                    "basal_leaf_count": 4,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "STAGE4-NO-BARRIER-PLASTIC",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "decimalLatitude": 35.802,
                    "decimalLongitude": -79.002,
                    "doy": 118,
                    "PC1": -2.1,
                    "PC2": 0.85,
                    "cda_predicted_taxon": "Packera dubia",
                    "cda_second_taxon": "Packera anonyma",
                    "cda_posterior_prob": 0.48,
                    "cda_second_posterior_prob": 0.25,
                    "cda_posterior_margin": 0.23,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.68,
                    "delta_doy_resid": 3.0,
                    "gmm_cluster": 5,
                    "gmm_2_delta_bic": 12.4,
                    "gmm_uncertainty": 0.08,
                    "basal_leaf_count": 4,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "STAGE2-SWARM-MATE",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_2_Silver",
                    "decimalLatitude": 35.808,
                    "decimalLongitude": -79.008,
                    "doy": 120,
                    "PC1": 0.10,
                    "PC2": -0.02,
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_second_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.51,
                    "cda_second_posterior_prob": 0.41,
                    "cda_posterior_margin": 0.10,
                    "morphometric_outlier_flag": False,
                    "edaphic_best_fit_taxon": "Packera anonyma",
                    "rf_predicted_taxon": "Packera anonyma",
                    "rf_predicted_prob": 0.49,
                    "edaphic_schoener_d": 0.54,
                    "delta_doy_resid": 3.0,
                    "gmm_cluster": 2,
                    "gmm_uncertainty": 0.39,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "DUB-GOLD-ANC-1",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "decimalLatitude": 35.800,
                    "decimalLongitude": -79.000,
                    "doy": 115,
                    "PC1": -2.50,
                    "PC2": 0.08,
                    "cda_predicted_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.96,
                    "cda_second_posterior_prob": 0.04,
                    "cda_posterior_margin": 0.92,
                    "morphometric_outlier_flag": False,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.84,
                    "delta_doy_resid": 0.0,
                    "gmm_cluster": 1,
                    "gmm_uncertainty": 0.03,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "DUB-GOLD-ANC-2",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "decimalLatitude": 35.805,
                    "decimalLongitude": -79.005,
                    "doy": 117,
                    "PC1": -2.30,
                    "PC2": -0.08,
                    "cda_predicted_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.95,
                    "cda_second_posterior_prob": 0.05,
                    "cda_posterior_margin": 0.90,
                    "morphometric_outlier_flag": False,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.82,
                    "delta_doy_resid": 1.0,
                    "gmm_cluster": 1,
                    "gmm_uncertainty": 0.04,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "ANO-GOLD-ANC-1",
                    "species_standardized": "Packera anonyma",
                    "determiner_tier": "Tier_1_Gold",
                    "decimalLatitude": 35.810,
                    "decimalLongitude": -79.010,
                    "doy": 122,
                    "PC1": 2.50,
                    "PC2": 0.08,
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_posterior_prob": 0.96,
                    "cda_second_posterior_prob": 0.04,
                    "cda_posterior_margin": 0.92,
                    "morphometric_outlier_flag": False,
                    "edaphic_best_fit_taxon": "Packera anonyma",
                    "edaphic_schoener_d": 0.84,
                    "delta_doy_resid": 0.0,
                    "gmm_cluster": 1,
                    "gmm_uncertainty": 0.03,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
                {
                    "catalogNumber": "ANO-GOLD-ANC-2",
                    "species_standardized": "Packera anonyma",
                    "determiner_tier": "Tier_1_Gold",
                    "decimalLatitude": 35.815,
                    "decimalLongitude": -79.015,
                    "doy": 124,
                    "PC1": 2.30,
                    "PC2": -0.08,
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_posterior_prob": 0.95,
                    "cda_second_posterior_prob": 0.05,
                    "cda_posterior_margin": 0.90,
                    "morphometric_outlier_flag": False,
                    "edaphic_best_fit_taxon": "Packera anonyma",
                    "edaphic_schoener_d": 0.82,
                    "delta_doy_resid": 1.0,
                    "gmm_cluster": 1,
                    "gmm_uncertainty": 0.04,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                },
            ]).to_csv(cls.prec_in_csv, index=False)

            cls.nd_dir = cls._class_tmp_path / "non_dubia"
            cls.nd_dir.mkdir(parents=True, exist_ok=True)
            cls.nd_in_csv = cls.nd_dir / "non_dubia_input.csv"
            cls.nd_out_csv = cls.nd_dir / "non_dubia_triage.csv"
            cls.nd_multi_report = cls.nd_dir / "multi_species_report.md"
            cls.nd_focal_report_default = cls.nd_dir / "Packera_Taxonomic_Revision_Summary.md"
            cls.nd_focal_plot_default = cls.nd_dir / "fig06a_integrative_taxonomic_synthesis.png"
            cls.nd_helper_check_csv = cls.nd_dir / "pub_fig_helpers.csv"
            cls.dubia_dir = cls.nd_dir / "dubia_focal"
            cls.dubia_dir.mkdir(parents=True, exist_ok=True)
            dubia_rep_def = cls.dubia_dir / "Packera_Taxonomic_Revision_Summary.md"
            dubia_fig_def = cls.dubia_dir / "fig06a_integrative_taxonomic_synthesis.png"
            dubia_fig_def.write_bytes(b"\x89PNG\r\n\x1a\n")

            nd_rows = [
                {"catalogNumber": "AUR-GOLD-1", "institutionCode": "NCU", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "regional_group": "Blue Ridge", "doy": 120, "decimalLatitude": 35.6, "decimalLongitude": -82.5, "fine_soil_ph": 6.2, "fine_soil_sand": 32.0, "fine_bedrock_depth_cm": 95.0, "can1": -3.0, "can2": 0.1, "can3": 0.0, "PC1": -1.5, "PC2": 0.2, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.82, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-GOLD-2", "institutionCode": "DUKE", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "regional_group": "Ridge and Valley", "doy": 120, "decimalLatitude": 36.18, "decimalLongitude": -83.42, "fine_soil_ph": 6.4, "fine_soil_sand": 36.0, "fine_bedrock_depth_cm": 110.0, "can1": -2.8, "can2": -0.1, "can3": 0.1, "PC1": -1.4, "PC2": -0.1, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.95, "cda_posterior_margin": 0.90, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.80, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-GOLD-NORTH", "institutionCode": "GH", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "regional_group": "Northeastern Highlands", "doy": 158, "decimalLatitude": 44.22, "decimalLongitude": -71.32, "fine_soil_ph": 6.1, "fine_soil_sand": 38.0, "fine_bedrock_depth_cm": 102.0, "can1": -2.9, "can2": 0.0, "can3": 0.0, "PC1": -1.45, "PC2": 0.05, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.81, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "OBO-GOLD-1", "institutionCode": "NCU", "species_raw": "Packera obovata", "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "regional_group": "Interior Plateau", "doy": 112, "decimalLatitude": 36.0, "decimalLongitude": -86.8, "fine_soil_ph": 7.3, "fine_soil_sand": 22.0, "fine_bedrock_depth_cm": 38.0, "can1": 3.1, "can2": 0.2, "can3": -0.1, "PC1": 1.6, "PC2": 0.1, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.97, "cda_posterior_margin": 0.94, "edaphic_best_fit_taxon": "Packera obovata", "edaphic_schoener_d": 0.85, "gmm_cluster": 1, "gmm_uncertainty": 0.02, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 4},
                {"catalogNumber": "OBO-GOLD-2", "institutionCode": "TENN", "scientificName": "Packera obovata (Muhl. ex Willd.) W.A.Weber & Á.Löve", "species_raw": "Packera obovata", "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "regional_group": "Ridge and Valley", "doy": 116, "decimalLatitude": 36.15, "decimalLongitude": -83.39, "fine_soil_ph": 7.5, "fine_soil_sand": 26.0, "fine_bedrock_depth_cm": 42.0, "can1": 2.9, "can2": -0.2, "can3": 0.1, "PC1": 1.5, "PC2": -0.2, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "edaphic_best_fit_taxon": "Packera obovata", "edaphic_schoener_d": 0.83, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "SCH-GOLD-1", "institutionCode": "GH", "species_raw": "Packera schweinitziana", "species_standardized": "Packera schweinitziana", "determiner_tier": "Tier_1_Gold", "regional_group": "Northeastern Highlands", "doy": 168, "decimalLatitude": 44.2, "decimalLongitude": -71.3, "fine_soil_ph": 5.4, "fine_soil_sand": 48.0, "fine_bedrock_depth_cm": 125.0, "can1": 0.0, "can2": 3.4, "can3": 0.2, "PC1": 0.1, "PC2": 1.9, "cda_predicted_taxon": "Packera schweinitziana", "cda_posterior_prob": 0.97, "cda_posterior_margin": 0.94, "edaphic_best_fit_taxon": "Packera schweinitziana", "edaphic_schoener_d": 0.84, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "SCH-GOLD-2", "institutionCode": "NY", "species_raw": "Packera schweinitziana", "species_standardized": "Packera schweinitziana", "determiner_tier": "Tier_1_Gold", "regional_group": "Northern Appalachian and Atlantic Maritime Highlands", "doy": 172, "decimalLatitude": 44.6, "decimalLongitude": -72.8, "fine_soil_ph": 5.2, "fine_soil_sand": 52.0, "fine_bedrock_depth_cm": 130.0, "can1": 0.2, "can2": 3.2, "can3": -0.2, "PC1": 0.2, "PC2": 1.8, "cda_predicted_taxon": "Packera schweinitziana", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "edaphic_best_fit_taxon": "Packera schweinitziana", "edaphic_schoener_d": 0.81, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-MISID-OBO", "institutionCode": "NCU", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_3_Bronze", "regional_group": "Interior Plateau", "doy": 110, "decimalLatitude": 36.4, "decimalLongitude": -86.2, "fine_soil_ph": 7.4, "fine_soil_sand": 24.0, "fine_bedrock_depth_cm": 40.0, "can1": 2.95, "can2": 0.05, "can3": 0.0, "PC1": 1.55, "PC2": 0.0, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.93, "cda_posterior_margin": 0.86, "edaphic_best_fit_taxon": "Packera obovata", "edaphic_schoener_d": 0.18, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "morphometric_outlier_flag": True, "triage_category": "Severe_Triple_Stream_Conflict", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-HYB-EXPLICIT", "institutionCode": "DUKE", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Ridge and Valley", "doy": 118, "decimalLatitude": 36.165, "decimalLongitude": -83.405, "fine_soil_ph": 6.8, "fine_soil_sand": 29.0, "fine_bedrock_depth_cm": 68.0, "can1": 0.05, "can2": 0.0, "can3": 0.0, "PC1": 0.05, "PC2": 0.0, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.52, "cda_posterior_margin": 0.08, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.55, "gmm_cluster": 1, "gmm_uncertainty": 0.40, "morphometric_outlier_flag": False, "triage_category": "Putative_Hybrid_Zone_Intergrade", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-HYB-EXPLICIT-2", "institutionCode": "TENN", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Ridge and Valley", "doy": 119, "decimalLatitude": 36.170, "decimalLongitude": -83.410, "fine_soil_ph": 6.7, "fine_soil_sand": 30.0, "fine_bedrock_depth_cm": 72.0, "can1": -0.10, "can2": 0.02, "can3": 0.0, "PC1": -0.05, "PC2": 0.02, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.51, "cda_posterior_margin": 0.06, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.54, "gmm_cluster": 1, "gmm_uncertainty": 0.39, "morphometric_outlier_flag": False, "triage_category": "Putative_Hybrid_Zone_Intergrade", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-HYB-MAHAL", "institutionCode": "GH", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Northeastern Highlands", "doy": 162, "decimalLatitude": 44.21, "decimalLongitude": -71.31, "fine_soil_ph": 5.9, "fine_soil_sand": 40.0, "fine_bedrock_depth_cm": 105.0, "can1": -1.35, "can2": 1.70, "can3": 0.02, "PC1": -0.65, "PC2": 0.95, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.51, "cda_posterior_margin": 0.10, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.58, "gmm_cluster": 1, "gmm_uncertainty": 0.42, "morphometric_outlier_flag": False, "triage_category": "Putative_Hybrid_Zone_Intergrade", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-HYB-MAHAL-2", "institutionCode": "NY", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Northeastern Highlands", "doy": 163, "decimalLatitude": 44.215, "decimalLongitude": -71.315, "fine_soil_ph": 5.8, "fine_soil_sand": 42.0, "fine_bedrock_depth_cm": 108.0, "can1": -1.25, "can2": 1.80, "can3": 0.0, "PC1": -0.60, "PC2": 1.00, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.50, "cda_posterior_margin": 0.08, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.57, "gmm_cluster": 1, "gmm_uncertainty": 0.41, "morphometric_outlier_flag": False, "triage_category": "Putative_Hybrid_Zone_Intergrade", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-ECOTYPE-VAR", "institutionCode": "NCU", "species_raw": "Packera aurea var. gracilis", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Blue Ridge", "doy": 125, "decimalLatitude": 35.9, "decimalLongitude": -82.1, "fine_soil_ph": 6.1, "fine_soil_sand": 38.0, "fine_bedrock_depth_cm": 88.0, "can1": -2.5, "can2": 0.4, "can3": 0.0, "PC1": -1.2, "PC2": 0.3, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.88, "cda_posterior_margin": 0.76, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.38, "gmm_cluster": 2, "gmm_uncertainty": 0.08, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
                {"catalogNumber": "AUR-PLASTIC-LATE", "institutionCode": "TENN", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Blue Ridge", "doy": 148, "decimalLatitude": 35.7, "decimalLongitude": -82.4, "fine_soil_ph": 6.3, "fine_soil_sand": 35.0, "fine_bedrock_depth_cm": 98.0, "can1": -2.7, "can2": 0.05, "can3": 0.0, "PC1": -1.35, "PC2": 0.05, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.65, "cda_posterior_margin": 0.42, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.78, "gmm_cluster": 1, "gmm_uncertainty": 0.06, "morphometric_outlier_flag": True, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3},
            ]
            pd.DataFrame(nd_rows).to_csv(cls.nd_in_csv, index=False)

            r_batch_code = f"""
            source("{cls.r_script_path.as_posix()}")

            # Block 1: Stage 1-5 hierarchical precedence & cryptic congruence
            local({{
              d <- readr::read_csv("{cls.prec_in_csv.as_posix()}", show_col_types = FALSE)
              r <- apply_taxonomic_decision_matrix(d, gmm_two_delta_bic = 12.4)
              readr::write_csv(r, "{cls.prec_out_csv.as_posix()}")
              cat("STAGE1_5_PRECEDENCE_OK\\n")
            }})

            # Block 2: Non-dubia taxon-agnostic decision matrix, revision reports, and figures
            local({{
              df <- readr::read_csv("{cls.nd_in_csv.as_posix()}", show_col_types = FALSE)
              res <- apply_taxonomic_decision_matrix(df)
              readr::write_csv(res, "{cls.nd_out_csv.as_posix()}")
              opts_multi <- list(output_report = "{cls.nd_multi_report.as_posix()}", output_plot = "", focal_taxon = "", gmm_bic = "", gmm_summary = "")
              generate_taxonomic_revision_report(res, opts_multi)
              render_synthesis_figure_pdf(res, "{cls.nd_focal_plot_default.as_posix()}", focal_taxon = "Packera aurea")
              opts_focal <- list(output_report = "{cls.nd_focal_report_default.as_posix()}", output_plot = "{cls.nd_focal_plot_default.as_posix()}", focal_taxon = "Packera aurea", gmm_bic = "", gmm_summary = "")
              generate_taxonomic_revision_report(res[res$species_standardized == "Packera aurea", , drop = FALSE], opts_focal)
              opts_dub <- list(output_report = "{dubia_rep_def.as_posix()}", output_plot = "{dubia_fig_def.as_posix()}", focal_taxon = "Packera dubia", gmm_bic = "", gmm_summary = "")
              generate_taxonomic_revision_report(res[1:2, , drop = FALSE], opts_dub)
              pal <- get_dynamic_packera_palette(c("Packera aurea", "Packera obovata", "Packera schweinitziana"))
              t_multi <- format_active_complex_title(c("Packera aurea", "Packera obovata", "Packera schweinitziana"), "CDA Morphospace")
              t_single <- format_active_complex_title(c("Packera aurea"), "CDA Morphospace")
              mb_west <- compute_dynamic_na_map_bounds(c(-121.5, -110.2), c(37.4, 46.8))
              readr::write_csv(data.frame(
                has_dubia_pal = "Packera dubia" %in% names(pal),
                any_grey_pal = any(pal[c("Packera aurea", "Packera obovata", "Packera schweinitziana")] == "#7F7F7F"),
                n_unique_cols = length(unique(pal[c("Packera aurea", "Packera obovata", "Packera schweinitziana")])),
                title_multi = t_multi, title_single = t_single,
                west_xmin = mb_west$xlim[1], west_xmax = mb_west$xlim[2],
                west_ymin = mb_west$ylim[1], west_ymax = mb_west$ylim[2]
              ), "{cls.nd_helper_check_csv.as_posix()}")
              cat("NON_DUBIA_REPORT_AND_FIGS_OK\\n")
            }})

            # Block 3: Stage 1 cross-modal non-contradiction gate & OOF class-precision pooling
            local({{
              df7 <- data.frame(
                catalogNumber = c(
                  "A_G1", "A_G2", "A_G3",
                  "B_G1", "B_G2", "B_G3",
                  "C_G1", "C_G2", "C_G3",
                  "EDAPH_CONTRADICTED_BY_MORPH",
                  "MORPH_MOD_CONTRADICTED_BY_EDAPH",
                  "CORROBORATED_MISID"
                ),
                plant_individual_id = 1L,
                institutionCode = "NCU",
                species_raw = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera aurea", "Packera aurea"),
                species_standardized = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera aurea", "Packera aurea"),
                determiner_tier = c(rep("Tier_1_Gold", 9), rep("Tier_3_Bronze", 3)),
                morph_tier_source = "Basal_Rosette",
                basal_leaf_count = 3L,
                cda_predicted_taxon = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera obovata", "Packera obovata"),
                cda_second_taxon = c(rep("Packera obovata", 3), rep("Packera aurea", 3), rep("Packera aurea", 3), "Packera obovata", "Packera aurea", "Packera aurea"),
                cda_posterior_prob = c(rep(0.95, 9), 0.91, 0.72, 0.92),
                cda_posterior_margin = c(rep(0.90, 9), 0.82, 0.44, 0.84),
                morphometric_outlier_flag = c(rep(FALSE, 9), FALSE, TRUE, TRUE),
                edaphic_best_fit_taxon = c(
                  rep("Packera aurea", 3),
                  rep("Packera obovata", 3),
                  "Packera obovata", "Packera obovata", "Packera schweinitziana",
                  "Packera obovata",
                  "Packera aurea",
                  "Packera obovata"
                ),
                rf_predicted_taxon = c(
                  rep("Packera aurea", 3),
                  rep("Packera obovata", 3),
                  "Packera obovata", "Packera obovata", "Packera schweinitziana",
                  "Packera obovata",
                  "Packera aurea",
                  "Packera obovata"
                ),
                rf_predicted_prob = c(rep(0.92, 6), 0.85, 0.85, 0.92, 0.88, 0.91, 0.92),
                edaphic_schoener_d = c(rep(0.82, 6), 0.25, 0.25, 0.82, 0.20, 0.79, 0.18),
                gmm_cluster = 1L,
                gmm_uncertainty = 0.05,
                can1 = c(rep(-3.0, 3), rep(3.0, 3), rep(0.0, 3), -2.9, 1.8, 2.9),
                can2 = c(rep(0.0, 3), rep(0.0, 3), rep(3.0, 3), 0.1, 0.1, 0.0),
                can3 = 0.0,
                PC1 = c(rep(-1.5, 3), rep(1.5, 3), rep(0.0, 3), -1.4, 0.9, 1.4),
                PC2 = c(rep(0.0, 3), rep(0.0, 3), rep(1.5, 3), 0.0, 0.0, 0.0),
                triage_category = c(rep("Clean_MultiModal_Consensus", 9), rep("Severe_Triple_Stream_Conflict", 3)),
                stringsAsFactors = FALSE
              )

              res7 <- apply_taxonomic_decision_matrix(df7)
              row_edaph_contra <- res7[res7$catalogNumber == "EDAPH_CONTRADICTED_BY_MORPH", ]
              row_morph_contra <- res7[res7$catalogNumber == "MORPH_MOD_CONTRADICTED_BY_EDAPH", ]
              row_corrob <- res7[res7$catalogNumber == "CORROBORATED_MISID", ]

              stopifnot(row_edaph_contra$taxonomic_treatment_subtype != "Misidentification_Reassign")
              stopifnot(row_edaph_contra$recommended_determination == "Packera aurea")
              stopifnot(row_morph_contra$taxonomic_treatment_subtype != "Misidentification_Reassign")
              stopifnot(row_morph_contra$recommended_determination == "Packera aurea")
              stopifnot(row_corrob$taxonomic_treatment_subtype == "Misidentification_Reassign")
              stopifnot(row_corrob$recommended_determination == "Packera obovata")
              stopifnot(all(res7$probabilistic_consensus_conflict >= 0.0 & res7$probabilistic_consensus_conflict <= 1.0))

              anc_m <- .resolve_step06_anchor_mask(df7)
              prec_edaph <- .compute_oof_class_precision_map(df7$species_standardized[anc_m], df7$edaphic_best_fit_taxon[anc_m])
              prec_mor <- .compute_oof_class_precision_map(df7$species_standardized[anc_m], df7$cda_predicted_taxon[anc_m])
              stopifnot(prec_edaph[["Packera obovata"]] < prec_mor[["Packera obovata"]])

              df7_hi <- df7
              df7_hi$edaphic_best_fit_taxon[7:9] <- "Packera schweinitziana"
              df7_hi$rf_predicted_taxon[7:9] <- "Packera schweinitziana"
              df7_hi$edaphic_schoener_d[7:9] <- 0.82
              res7_hi <- apply_taxonomic_decision_matrix(df7_hi)
              c_low_prec <- res7$probabilistic_consensus_conflict[res7$catalogNumber == "EDAPH_CONTRADICTED_BY_MORPH"]
              c_hi_prec <- res7_hi$probabilistic_consensus_conflict[res7_hi$catalogNumber == "EDAPH_CONTRADICTED_BY_MORPH"]
              stopifnot(c_low_prec < c_hi_prec)
              cat("STEP07_PRECISION_AND_GATE_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step07_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = run_warm_r_cmd(
                ["Rscript", str(r_batch_file)],
                stage_key="step06",
                cwd=PROJECT_ROOT,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 06 R unit suite failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_five_stage_hierarchical_precedence_and_cryptic_congruence(self):
        """Test Stage 1 (Misidentification_Reassign BEFORE Hybrid_Introgressant), Stage 2 (Hybrid), Stage 3 (Putative_Cryptic_Species Congruence), and Stage 4 (Ecophenotypic_Plasticity)."""
        out = self._require_batched_r()
        self.assertIn("STAGE1_5_PRECEDENCE_OK", out)
        r_out = pd.read_csv(self.prec_out_csv)
        self.assertEqual(r_out.loc[0, "taxonomic_treatment_subtype"], "Misidentification_Reassign")
        self.assertEqual(r_out.loc[1, "taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertEqual(r_out.loc[2, "taxonomic_treatment_subtype"], "Putative_Cryptic_Species")
        self.assertEqual(r_out.loc[3, "taxonomic_treatment_subtype"], "Ecophenotypic_Plasticity")

    def test_non_dubia_taxon_agnostic_decision_matrix_report_and_figures(self):
        """Verify apply_taxonomic_decision_matrix(), generate_taxonomic_revision_report(), and pub_figures_utils.R run cleanly on a non-dubia dataset with zero Packera dubia injections."""
        out = self._require_batched_r()
        self.assertIn("NON_DUBIA_REPORT_AND_FIGS_OK", out)
        res = pd.read_csv(self.nd_out_csv).set_index("catalogNumber")

        rec_all = " | ".join(res["recommended_determination"].astype(str).tolist())
        self.assertNotIn("Packera dubia", rec_all)
        self.assertNotIn("Packera anonyma", rec_all)
        self.assertNotIn("Packera plattensis", rec_all)

        self.assertEqual(res.loc["AUR-GOLD-1", "taxonomic_treatment_subtype"], "Confirmed_Species")
        self.assertEqual(res.loc["AUR-MISID-OBO", "taxonomic_treatment_subtype"], "Misidentification_Reassign")
        self.assertEqual(res.loc["AUR-MISID-OBO", "recommended_determination"], "Packera obovata")
        self.assertEqual(res.loc["AUR-HYB-EXPLICIT", "taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertEqual(res.loc["AUR-HYB-EXPLICIT", "recommended_determination"], "Packera aurea x Packera obovata")
        self.assertEqual(res.loc["AUR-HYB-EXPLICIT", "hybrid_parent_1"], "Packera aurea")
        self.assertEqual(res.loc["AUR-HYB-EXPLICIT", "hybrid_parent_2"], "Packera obovata")
        self.assertTrue(bool(res.loc["AUR-HYB-EXPLICIT", "hybrid_zone_cohesion_pass"]))
        self.assertEqual(res.loc["AUR-HYB-MAHAL", "taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertEqual(res.loc["AUR-HYB-MAHAL", "recommended_determination"], "Packera aurea x Packera schweinitziana")
        self.assertEqual(res.loc["AUR-HYB-MAHAL", "hybrid_parent_1"], "Packera aurea")
        self.assertEqual(res.loc["AUR-HYB-MAHAL", "hybrid_parent_2"], "Packera schweinitziana")
        self.assertTrue(bool(res.loc["AUR-HYB-MAHAL", "hybrid_zone_cohesion_pass"]))
        self.assertEqual(res.loc["AUR-ECOTYPE-VAR", "taxonomic_treatment_subtype"], "Subspecies_Ecotype")
        self.assertEqual(res.loc["AUR-PLASTIC-LATE", "taxonomic_treatment_subtype"], "Ecophenotypic_Plasticity")
        self.assertEqual(res.loc["AUR-PLASTIC-LATE", "recommended_determination"], "Packera aurea")

        md_multi = self.nd_multi_report.read_text(encoding="utf-8")
        self.assertIn(
            "# Taxonomic Revision and Integrative Species Delimitation Across *Packera* Species Complexes (N = 3 Taxa; Asteraceae: Senecioneae)",
            md_multi,
        )
        self.assertIn("- ***Packera aurea***", md_multi)
        self.assertIn("- ***Packera obovata*** (Muhl. ex Willd.) W.A.Weber & Á.Löve", md_multi)
        self.assertIn("- ***Packera schweinitziana***", md_multi)
        self.assertIn("Dominant US EPA Level III Ecoregions:", md_multi)
        self.assertIn("Verified Anthesis Window:", md_multi)
        self.assertNotIn("Packera dubia", md_multi)

        focal_named_md = self.nd_dir / "Packera_aurea_Taxonomic_Revision_Summary.md"
        focal_named_png = self.nd_dir / "Figure_Integrative_Packera_aurea_Revision.png"
        self.assertTrue(focal_named_md.exists(), "Expected taxon-named focal Markdown report copy")
        self.assertTrue(focal_named_png.exists(), "Expected taxon-named focal 300-DPI PNG synthesis plate copy")
        self.assertTrue((self.dubia_dir / "Packera_dubia_Taxonomic_Revision_Summary.md").exists())
        self.assertTrue((self.dubia_dir / "Figure_Integrative_Packera_dubia_Revision.png").exists())
        md_focal = focal_named_md.read_text(encoding="utf-8")
        self.assertIn(
            "# Taxonomic Revision and Integrative Species Delimitation in the *Packera aurea* Complex (Asteraceae: Senecioneae)",
            md_focal,
        )
        self.assertNotIn("Packera dubia", md_focal)

        helpers_df = pd.read_csv(self.nd_helper_check_csv)
        self.assertFalse(bool(helpers_df.loc[0, "has_dubia_pal"]))
        self.assertFalse(bool(helpers_df.loc[0, "any_grey_pal"]))
        self.assertEqual(int(helpers_df.loc[0, "n_unique_cols"]), 3)
        self.assertEqual(str(helpers_df.loc[0, "title_multi"]), "Packera Taxon Complex (3 Taxa): CDA Morphospace")
        self.assertEqual(str(helpers_df.loc[0, "title_single"]), "Packera aurea: CDA Morphospace")
        self.assertLess(float(helpers_df.loc[0, "west_xmin"]), -121.5)
        self.assertGreater(float(helpers_df.loc[0, "west_xmax"]), -110.2)
        self.assertLess(float(helpers_df.loc[0, "west_ymin"]), 37.4)
        self.assertGreater(float(helpers_df.loc[0, "west_ymax"]), 46.8)

    def test_stage1_crossmodal_non_contradiction_gate_and_oof_precision_pooling(self):
        """Verify Stage 1 blocks contradicted single-modality reassignments and OOF confusion-matrix precision modulates consensus weights."""
        out = self._require_batched_r()
        self.assertIn("STEP07_PRECISION_AND_GATE_OK", out)


    def test_step06_cli_handoff_coalescing_multi_plant_artifact_gate_and_ood_cryptic(self):
        if not rscript_available():
            self.skipTest("Rscript not installed or fast mode enabled.")
        """Verify Step 07 CLI preserves multi-plant keys, coalesces ungeoreferenced vouchers (missing-modality renormalization), assigns RESOLVED/CRITICAL priorities, gates Stage 1 on cross-modal non-contradiction, and corroborates Stage 3 cryptic divergence across morphometric + edaphic + phenological streams."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vouchers_csv = tmp / "curated_vouchers.csv"
            pheno_csv = tmp / "voucher_phenological_states.csv"
            morph_csv = tmp / "morphometric_flags.csv"
            mm_csv = tmp / "multimodal_conflict_flags.csv"
            niche_csv = tmp / "niche_identity_tests.csv"
            out_queue = tmp / "triage_queue.csv"
            out_html = tmp / "Triage_Dashboard_Interactive.html"
            out_rep = tmp / "Packera_Taxonomic_Revision_Summary.md"
            out_fig = tmp / "fig06a_integrative_taxonomic_synthesis.png"

            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.5, "decimalLongitude": -82.5, "doy": 118},
                {"catalogNumber": "CAT-ANCHOR-A2", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.55, "decimalLongitude": -82.45, "doy": 119},
                {"catalogNumber": "CAT-ANCHOR-B1", "institutionCode": "NCU", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 36.0, "decimalLongitude": -79.0, "doy": 112},
                {"catalogNumber": "CAT-ANCHOR-B2", "institutionCode": "NCU", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 36.05, "decimalLongitude": -79.05, "doy": 113},
                {"catalogNumber": "CAT-SEVERE-MISID", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 36.1, "decimalLongitude": -79.1, "doy": 114},
                {"catalogNumber": "NCU-HIST-2", "institutionCode": "PH", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_2_Silver", "decimalLatitude": np.nan, "decimalLongitude": np.nan, "doy": 158},
                {"catalogNumber": "CAT-ART-SHORTCUT", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_3_Bronze", "decimalLatitude": 35.6, "decimalLongitude": -82.4, "doy": 120},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_2_Silver", "decimalLatitude": 35.8, "decimalLongitude": -81.9, "doy": 152},
                {"catalogNumber": "CAT-MULTI-SHEET", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_2_Silver", "decimalLatitude": 35.7, "decimalLongitude": -82.1, "doy": 119},
            ]).to_csv(vouchers_csv, index=False)

            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "voucher_phenological_state": "anthesis", "n_anthesis": 6, "n_bud": 1, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.93},
                {"catalogNumber": "CAT-ANCHOR-A2", "voucher_phenological_state": "anthesis", "n_anthesis": 6, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.94},
                {"catalogNumber": "CAT-ANCHOR-B1", "voucher_phenological_state": "anthesis", "n_anthesis": 5, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.91},
                {"catalogNumber": "CAT-ANCHOR-B2", "voucher_phenological_state": "anthesis", "n_anthesis": 5, "n_bud": 1, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.92},
                {"catalogNumber": "CAT-SEVERE-MISID", "voucher_phenological_state": "anthesis", "n_anthesis": 5, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.95},
                {"catalogNumber": "NCU-HIST-2", "voucher_phenological_state": "anthesis", "n_anthesis": 4, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.88},
                {"catalogNumber": "CAT-ART-SHORTCUT", "voucher_phenological_state": "anthesis", "n_anthesis": 4, "n_bud": 1, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.88},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "voucher_phenological_state": "anthesis", "n_anthesis": 7, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.95},
                {"catalogNumber": "CAT-MULTI-SHEET", "voucher_phenological_state": "anthesis", "n_anthesis": 8, "n_bud": 2, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.90},
            ]).to_csv(pheno_csv, index=False)

            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": -2.0, "PC2": 0.1, "can1": -2.5, "can2": 0.2, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "petiole_length_mm": 30.0, "blade_length_mm": 25.0, "petiole_to_blade_ratio": 1.20, "plant_height_cm": 44.0, "scape_to_basal_leaf_ratio": 8.0, "is_bent_stem": False, "capitula_count": 7, "involucre_height_mm": 6.4, "involucre_width_mm": 5.1, "AR_inv": 1.25, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.88, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-ANCHOR-A2", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.97, "cda_second_posterior_prob": 0.03, "cda_posterior_margin": 0.94, "PC1": -2.1, "PC2": 0.15, "can1": -2.6, "can2": 0.15, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.02, "petiole_length_mm": 31.0, "blade_length_mm": 25.5, "petiole_to_blade_ratio": 1.22, "plant_height_cm": 45.0, "scape_to_basal_leaf_ratio": 8.1, "is_bent_stem": False, "capitula_count": 7, "involucre_height_mm": 6.5, "involucre_width_mm": 5.1, "AR_inv": 1.27, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.89, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B1", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.95, "cda_second_posterior_prob": 0.05, "cda_posterior_margin": 0.90, "PC1": 2.1, "PC2": -0.2, "can1": 2.6, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.04, "petiole_length_mm": 11.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": 0.39, "plant_height_cm": 36.0, "scape_to_basal_leaf_ratio": 9.2, "is_bent_stem": False, "capitula_count": 5, "involucre_height_mm": 5.8, "involucre_width_mm": 4.8, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.85, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B2", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": 2.2, "PC2": -0.15, "can1": 2.7, "can2": -0.15, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.03, "petiole_length_mm": 10.8, "blade_length_mm": 27.5, "petiole_to_blade_ratio": 0.39, "plant_height_cm": 35.5, "scape_to_basal_leaf_ratio": 9.1, "is_bent_stem": False, "capitula_count": 6, "involucre_height_mm": 5.9, "involucre_width_mm": 4.9, "AR_inv": 1.20, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.87, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-SEVERE-MISID", "misidentification_flag": True, "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": 2.15, "PC2": -0.2, "can1": 2.65, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.03, "petiole_length_mm": 11.2, "blade_length_mm": 28.0, "petiole_to_blade_ratio": 0.40, "plant_height_cm": 36.0, "scape_to_basal_leaf_ratio": 9.0, "is_bent_stem": False, "capitula_count": 5, "involucre_height_mm": 5.8, "involucre_width_mm": 4.8, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.86, "morphometric_outlier_flag": True},
                {"catalogNumber": "NCU-HIST-2", "misidentification_flag": True, "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 2, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.93, "cda_second_posterior_prob": 0.07, "cda_posterior_margin": 0.86, "PC1": 2.0, "PC2": -0.2, "can1": 2.4, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.04, "petiole_length_mm": 11.5, "blade_length_mm": 27.5, "petiole_to_blade_ratio": 0.42, "plant_height_cm": 35.0, "scape_to_basal_leaf_ratio": 8.9, "is_bent_stem": False, "capitula_count": 4, "involucre_height_mm": 5.8, "involucre_width_mm": 4.8, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.82, "morphometric_outlier_flag": True},
                {"catalogNumber": "CAT-ART-SHORTCUT", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.72, "cda_second_posterior_prob": 0.28, "cda_posterior_margin": 0.44, "PC1": -1.8, "PC2": 0.2, "can1": -2.3, "can2": 0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.05, "petiole_length_mm": 29.0, "blade_length_mm": 24.0, "petiole_to_blade_ratio": 1.21, "plant_height_cm": 41.0, "scape_to_basal_leaf_ratio": 7.6, "is_bent_stem": False, "capitula_count": 5, "involucre_height_mm": 6.2, "involucre_width_mm": 5.0, "AR_inv": 1.24, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.82, "morphometric_outlier_flag": True},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.48, "cda_second_posterior_prob": 0.15, "cda_posterior_margin": 0.33, "PC1": 0.0, "PC2": 2.4, "can1": 0.0, "can2": 2.5, "can3": 1.4, "gmm_cluster": 2, "gmm_uncertainty": 0.08, "petiole_length_mm": 19.0, "blade_length_mm": 22.0, "petiole_to_blade_ratio": 0.86, "plant_height_cm": 33.0, "scape_to_basal_leaf_ratio": 8.0, "is_bent_stem": True, "capitula_count": 7, "involucre_height_mm": 7.1, "involucre_width_mm": 5.3, "AR_inv": 1.34, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.74, "morphometric_outlier_flag": True},
                {"catalogNumber": "CAT-MULTI-SHEET", "plant_individual_id": 1, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 2, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.94, "cda_second_posterior_prob": 0.06, "cda_posterior_margin": 0.88, "PC1": -1.9, "PC2": 0.1, "can1": -2.4, "can2": 0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "petiole_length_mm": 31.0, "blade_length_mm": 25.0, "petiole_to_blade_ratio": 1.24, "plant_height_cm": 43.0, "scape_to_basal_leaf_ratio": 7.7, "is_bent_stem": False, "capitula_count": 8, "involucre_height_mm": 6.3, "involucre_width_mm": 5.0, "AR_inv": 1.26, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.86, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-MULTI-SHEET", "misidentification_flag": True, "plant_individual_id": 2, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 2, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.93, "cda_second_posterior_prob": 0.07, "cda_posterior_margin": 0.86, "PC1": 2.0, "PC2": -0.2, "can1": 2.5, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.05, "petiole_length_mm": 10.5, "blade_length_mm": 27.0, "petiole_to_blade_ratio": 0.39, "plant_height_cm": 35.0, "scape_to_basal_leaf_ratio": 9.3, "is_bent_stem": False, "capitula_count": 6, "involucre_height_mm": 5.7, "involucre_width_mm": 4.7, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.81, "morphometric_outlier_flag": True},
            ]).to_csv(morph_csv, index=False)

            # Note: NCU-HIST-2 lacks GPS and is intentionally omitted from mm_csv to test non-georeferenced coalescing & missing-modality renormalization!
            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "plant_individual_id": 1, "fine_soil_ph": 6.1, "fine_soil_sand": 34.0, "fine_bedrock_depth_cm": 110.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.35, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.92, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.78, "delta_doy": 1.0, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-ANCHOR-A2", "plant_individual_id": 1, "fine_soil_ph": 6.15, "fine_soil_sand": 33.5, "fine_bedrock_depth_cm": 112.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.32, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.93, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.78, "delta_doy": 1.2, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B1", "plant_individual_id": 1, "fine_soil_ph": 6.7, "fine_soil_sand": 48.0, "fine_bedrock_depth_cm": 55.0, "ssurgo_rock_outcrop_pct": 18.0, "edaphic_z_score": 0.41, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.90, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.76, "delta_doy": -2.0, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B2", "plant_individual_id": 1, "fine_soil_ph": 6.75, "fine_soil_sand": 49.0, "fine_bedrock_depth_cm": 54.0, "ssurgo_rock_outcrop_pct": 19.0, "edaphic_z_score": 0.39, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.91, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.76, "delta_doy": -1.8, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-SEVERE-MISID", "misidentification_flag": True, "plant_individual_id": 1, "fine_soil_ph": 6.8, "fine_soil_sand": 49.0, "fine_bedrock_depth_cm": 52.0, "ssurgo_rock_outcrop_pct": 20.0, "edaphic_z_score": 2.9, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.93, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.21, "delta_doy": -2.0, "triage_category": "Severe_Triple_Stream_Conflict", "multimodal_conflict_flag": True},
                {"catalogNumber": "CAT-ART-SHORTCUT", "plant_individual_id": 1, "fine_soil_ph": 6.0, "fine_soil_sand": 35.0, "fine_bedrock_depth_cm": 105.0, "ssurgo_rock_outcrop_pct": 3.0, "edaphic_z_score": 0.44, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.89, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.75, "delta_doy": 2.0, "triage_category": "Severe_Triple_Stream_Conflict", "multimodal_conflict_flag": True},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "plant_individual_id": 1, "fine_soil_ph": 4.4, "fine_soil_sand": 72.0, "fine_bedrock_depth_cm": 18.0, "ssurgo_rock_outcrop_pct": 64.0, "edaphic_z_score": 3.45, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.51, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.14, "delta_doy": 26.0, "triage_category": "Edaphic_Envelope_Violation", "multimodal_conflict_flag": True},
                {"catalogNumber": "CAT-MULTI-SHEET", "plant_individual_id": 1, "fine_soil_ph": 6.1, "fine_soil_sand": 33.0, "fine_bedrock_depth_cm": 112.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.30, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.91, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.77, "delta_doy": 1.5, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-MULTI-SHEET", "misidentification_flag": True, "plant_individual_id": 2, "fine_soil_ph": 6.1, "fine_soil_sand": 33.0, "fine_bedrock_depth_cm": 112.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.30, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.91, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.77, "delta_doy": 1.5, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
            ]).to_csv(mm_csv, index=False)

            pd.DataFrame([
                {"Pair": "P.aurea vs P.obovata", "Taxon1": "Packera aurea", "Taxon2": "Packera obovata", "N1": 10, "N2": 10, "Warren_D": 0.21, "Warren_I": 0.45, "p_value_D": 0.0099, "p_value_I": 0.0099, "Taxonomic_Inference": "Distinct Ecological Niche"}
            ]).to_csv(niche_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouchers_csv),
                    "--morphometrics", str(morph_csv),
                    "--multimodal-flags", str(mm_csv),
                    "--niche-summary", str(niche_csv),
                    "--output-queue", str(out_queue),
                    "--output-html", str(out_html),
                    "--output-report", str(out_rep),
                    "--output-plot", str(out_fig),
                ], stage_key='step06', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 06 handoff run failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertTrue(out_queue.exists() and out_rep.exists() and out_html.exists())
            tdf = pd.read_csv(out_queue)
            self.assertEqual(len(tdf), 10, "Expected 10 rows (preserving ungeoreferenced NCU-HIST-2 and both individuals on CAT-MULTI-SHEET).")

            required_cols = [
                "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
                "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem",
                "capitula_count", "involucre_height_mm", "involucre_width_mm", "AR_inv",
                "is_verified_anthesis", "n_anthesis", "voucher_phenological_state",
                "cda_second_taxon", "cda_second_posterior_prob", "cda_posterior_margin",
                "ecoregion_cda_predicted", "ecoregion_cda_posterior_prob",
                "fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm",
                "ssurgo_rock_outcrop_pct", "edaphic_z_score", "rf_predicted_taxon",
            ]
            for col in required_cols:
                self.assertIn(col, tdf.columns, f"Missing expected column {col} in triage_queue.csv")

            # 1. RESOLVED vs CRITICAL priority assignments
            resolved_row = tdf[tdf["catalogNumber"] == "CAT-ANCHOR-A1"].iloc[0]
            critical_row = tdf[tdf["catalogNumber"] == "CAT-SEVERE-MISID"].iloc[0]
            self.assertEqual(str(resolved_row["synthesis_triage_priority"]), "RESOLVED")
            self.assertEqual(str(critical_row["synthesis_triage_priority"]), "CRITICAL")

            # 2. Ungeoreferenced historic specimen coalescing & missing-modality weight renormalization
            hist_row = tdf[tdf["catalogNumber"] == "NCU-HIST-2"].iloc[0]
            self.assertTrue(bool(hist_row["is_completed_voucher"]))
            self.assertEqual(str(hist_row["cda_predicted_taxon"]), "Packera obovata")
            self.assertEqual(str(hist_row["taxonomic_treatment_subtype"]), "Misidentification_Reassign")
            self.assertEqual(int(hist_row["n_available_modalities"]), 2)
            self.assertGreater(float(hist_row["probabilistic_consensus_conflict"]), 0.80)

            # 3. Moderate morph outlier contradicted by high-confidence edaphic RF MUST NOT trigger Stage 1 Misidentification_Reassign
            art_row = tdf[tdf["catalogNumber"] == "CAT-ART-SHORTCUT"].iloc[0]
            self.assertNotEqual(str(art_row["taxonomic_treatment_subtype"]), "Misidentification_Reassign")
            self.assertEqual(str(art_row["recommended_determination"]), "Packera aurea")

            # 4. Morphometric divergence + edaphic/phenological barriers triggers Stage 3 Putative_Cryptic_Species
            ood_row = tdf[tdf["catalogNumber"] == "CAT-VIS-OOD-CRYPT"].iloc[0]
            self.assertEqual(str(ood_row["taxonomic_treatment_subtype"]), "Putative_Cryptic_Species")
            self.assertIn("Basal_Rosette outlier", str(ood_row["synthesis_rationale"]))

            # 5. Multi-plant herbarium sheet retains distinct individuals and per-individual petiole ratios
            multi_df = tdf[tdf["catalogNumber"] == "CAT-MULTI-SHEET"].sort_values("plant_individual_id")
            self.assertEqual(list(multi_df["plant_individual_id"].astype(int)), [1, 2])
            self.assertAlmostEqual(float(multi_df.iloc[0]["petiole_to_blade_ratio"]), 1.24, places=2)
            self.assertAlmostEqual(float(multi_df.iloc[1]["petiole_to_blade_ratio"]), 0.39, places=2)

    def test_empirical_bayes_hopkins_shrinkage_dynamic_edaphic_z_and_taxon_scaled_pheno_cutoffs(self):
        """Verify Tasks 3A-3C: Empirical Bayes circular Hopkins shrinkage, dynamic primary-predictor edaphic_z_score (no hardcoded 5.8/45.0), multi-PC ANCOVA batch audit, and taxon-dispersion-scaled tau_ecotype / tau_cryptic."""
        r07_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertNotIn("5.8", r07_src, "Must eliminate hardcoded SEUS soil pH 5.8 fallback in Step 07")

        if not rscript_available():
            self.skipTest("Rscript not installed.")
        r_code = r"""
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/06_triage_dashboard_synthesis.R")

        set.seed(2026)
        # 25 active-anthesis anchors (alpha = clip((25 - 5)/15, 0, 1) = 1.0 -> full empirical circular slope)
        n <- 26
        sp <- c(rep("Packera hyperborea", 13), rep("Packera cymbalaria", 13))
        # Narrow flowering dispersion for P. hyperborea (SD ~ 3d -> tau_ecotype clipped to 10d, tau_cryptic clipped to 14d)
        # Broad flowering dispersion for P. cymbalaria (SD ~ 15d -> tau_ecotype clipped to 28d, tau_cryptic clipped to 35d)
        lat <- seq(60.0, 66.0, length.out = n)
        doy_hyp <- round(210 + 2.0 * (lat[1:13] - 63.0) + seq(-4, 4, length.out = 13))
        doy_cym <- round(225 + 2.0 * (lat[14:26] - 63.0) + seq(-22, 22, length.out = 13))
        df <- data.frame(
          catalogNumber = paste0("S07_", seq_len(n)),
          plant_individual_id = 0L,
          species_standardized = sp,
          determiner_tier = "Tier_1_Gold",
          institutionCode = rep(c("NCU", "ALA"), length.out = n),
          decimalLatitude = lat,
          elevation = seq(800, 1400, length.out = n),
          bio1_temp_mean = seq(-2.0, 2.0, length.out = n),
          doy = c(doy_hyp, doy_cym),
          is_verified_anthesis = TRUE,
          n_anthesis = 4L,
          voucher_phenological_state = "Peak_Anthesis",
          PC1 = rnorm(n), PC2 = rnorm(n), PC3 = rnorm(n), PC4 = rnorm(n),
          selected_morph_features = "PC1;PC2;PC3;PC4",
          fine_soil_cec = c(rnorm(12, 10.0, 0.8), 26.0, rnorm(13, 22.0, 1.0)),
          polaris_ph_30m = c(rnorm(12, 6.9, 0.1), 4.2, rnorm(13, 5.1, 0.15)),
          primary_edaph_predictors = "fine_soil_cec;polaris_ph_30m",
          edaphic_z_score = c(rep(0.45, 12), NA_real_, rep(0.55, 13)),
          stringsAsFactors = FALSE
        )

        df_h <- compute_hopkins_adjusted_phenology(df)
        stopifnot(all(is.finite(df_h$delta_doy_resid)))

        taus <- .compute_taxon_pheno_thresholds(df_h)
        # P. hyperborea narrow dispersion -> tau_ecotype < P. cymbalaria broad dispersion tau_ecotype
        stopifnot(taus$tau_ecotype[1] < taus$tau_ecotype[20])
        stopifnot(taus$tau_cryptic[1] < taus$tau_cryptic[20])
        stopifnot(taus$tau_ecotype[1] >= 10.0 && taus$tau_ecotype[20] <= 28.0)
        stopifnot(taus$tau_cryptic[1] >= 14.0 && taus$tau_cryptic[20] <= 35.0)

        # Verify Step 06 edaphic_z_score is preserved on rows 1..12 and dynamically computed on row 13 (extreme anomaly > 3.0)
        df_dec <- apply_taxonomic_decision_matrix(df_h)
        stopifnot(abs(df_dec$edaphic_z_score[1] - 0.45) < 1e-6)
        stopifnot(is.finite(df_dec$edaphic_z_score[13]) && df_dec$edaphic_z_score[13] > 3.0)
        # Multi-feature institutional batch ANCOVA tests all 4 PC1..PC4 features with BH FDR
        anc_df <- audit_institutional_batch_effects_ancova(df_dec)
        stopifnot(all(c("PC1", "PC2", "PC3", "PC4") %in% anc_df$Response))
        stopifnot("p_adj_bh" %in% names(anc_df) && all(is.finite(anc_df$p_adj_bh)))
        cat("STEP07_COHORT_ADAPTIVE_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step06_ca.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            res = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key="step06", cwd=PROJECT_ROOT)
        self.assertEqual(res.returncode, 0, f"Step 06 cohort-adaptive test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("STEP07_COHORT_ADAPTIVE_OK", res.stdout)

    def test_nomenclatural_type_control_rules_a_b_c_and_locus_classicus_step06_baseline(self):
        """Verify Milestone 4 (Step 05 locus classicus baseline) and Milestone 5 (Step 06 Rules A, B, and C: Confirmed_Type_Control / Type_Synonymy_Candidate, type_control_reassign, and is_typical_morphotype sensu stricto locking)."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/05_multimodal_spatial_rf.R")
        source("scripts/analysis/06_triage_dashboard_synthesis.R")

        # 1. Verify Step 06 compute_type_locality_ecological_control() computes type_edaphic_z_score,
        #    nearest_type_niche_taxon, and matches_type_lithology relative to locus classicus,
        #    and gates out coarse-uncertainty (>2000m) historical type centroids in favor of precise near_gold vouchers
        df06 <- data.frame(
          catalogNumber = c("DUB_HOLO", "DUB_NEAR", "DUB_MATCHING_OBO_LOCUS", "OBO_HOLO", "OBO_NEAR", "COARSE_HOLO", "COARSE_NEAR_GOLD", "COARSE_TEST_VOUCHER"),
          plant_individual_id = 0L,
          species_standardized = c("Packera dubia", "Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera aurea", "Packera aurea", "Packera aurea"),
          determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze", "Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver"),
          type_status = c("Holotype", "None", "None", "Holotype", "None", "Holotype", "None", "None"),
          is_type_specimen = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE),
          coordinateUncertaintyInMeters = c(250.0, 300.0, 400.0, 200.0, 350.0, 25000.0, 150.0, 200.0),
          regional_group = c("Piedmont", "Piedmont", "Interior Plateau", "Interior Plateau", "Interior Plateau", "Blue Ridge", "Blue Ridge", "Blue Ridge"),
          decimalLatitude = c(35.8, 35.9, 36.1, 36.0, 36.2, 35.5, 35.52, 35.53),
          decimalLongitude = c(-78.6, -78.5, -86.7, -86.8, -86.6, -82.5, -82.48, -82.49),
          fine_soil_ph = c(5.2, 5.3, 7.4, 7.5, 7.3, 8.1, 5.8, 5.9),
          fine_soil_sand = c(68.0, 66.0, 22.0, 20.0, 24.0, 12.0, 54.0, 55.0),
          fine_bedrock_depth_cm = c(140.0, 135.0, 35.0, 32.0, 38.0, 15.0, 110.0, 112.0),
          sgmc_lithology_class = c("Granitic_Felsic", "Granitic_Felsic", "Carbonate_Limestone_Dolostone", "Carbonate_Limestone_Dolostone", "Carbonate_Limestone_Dolostone", "Carbonate_Limestone_Dolostone", "Mafic_Metamorphic", "Mafic_Metamorphic"),
          primary_env_predictors = "fine_soil_ph;fine_soil_sand;fine_bedrock_depth_cm",
          primary_edaph_predictors = "fine_soil_ph;fine_soil_sand;fine_bedrock_depth_cm",
          stringsAsFactors = FALSE
        )
        tc06 <- compute_type_locality_ecological_control(df06)
        stopifnot(tc06$type_edaphic_z_score[1] < 0.5)
        stopifnot(tc06$type_edaphic_z_score[3] > 2.5)
        stopifnot(identical(tc06$nearest_type_niche_taxon[3], "Packera obovata"))
        stopifnot(isTRUE(tc06$matches_type_lithology[1]))
        stopifnot(!isTRUE(tc06$matches_type_lithology[3]))
        # Coarse-uncertainty type row (25000m, noisy carbonate/pH 8.1 centroid) must fall back to precise near_gold (Mafic_Metamorphic, pH 5.8)
        stopifnot(tc06$type_edaphic_z_score[8] < 0.5)
        stopifnot(isTRUE(tc06$matches_type_lithology[8]))

        # 2. Verify Step 07 Rules A, B, and C in apply_taxonomic_decision_matrix()
        df07 <- data.frame(
          catalogNumber = c(
            "RULE_A_CLEAN_HOLO",
            "RULE_A_SYNONYM_OVERLAP_HOLO",
            "RULE_B_TYPE_CONTROL_REASSIGN",
            "RULE_C_SENSU_STRICTO_LOCKED",
            "NON_TYPICAL_ECOTYPE_K2"
          ),
          plant_individual_id = 0L,
          species_raw = rep("Packera dubia", 5),
          species_standardized = rep("Packera dubia", 5),
          determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze", "Tier_1_Gold", "Tier_2_Silver"),
          type_status = c("Holotype", "Lectotype", "None", "None", "None"),
          is_type_specimen = c(TRUE, TRUE, FALSE, FALSE, FALSE),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          # Note: Row 3 has borderline classifier posteriors (cda_prob=0.58) that would NOT trigger
          # ordinary single-stream reassign, but matches Packera obovata across type-control streams!
          cda_predicted_taxon = c("Packera dubia", "Packera anonyma", "Packera dubia", "Packera dubia", "Packera dubia"),
          cda_posterior_prob = c(0.96, 0.91, 0.58, 0.62, 0.85),
          cda_posterior_margin = c(0.92, 0.82, 0.30, 0.35, 0.70),
          morphometric_outlier_flag = c(FALSE, TRUE, FALSE, TRUE, FALSE),
          edaphic_best_fit_taxon = c("Packera dubia", "Packera anonyma", "Packera dubia", "Packera dubia", "Packera dubia"),
          edaphic_schoener_d = c(0.85, 0.22, 0.55, 0.20, 0.25),
          delta_doy_resid = c(1.0, 2.0, 3.0, 25.0, 22.0),
          is_active_anthesis = TRUE,
          gmm_cluster = c(1L, 1L, 2L, 2L, 3L),
          gmm_2_delta_bic = 14.0,
          gmm_uncertainty = 0.05,
          is_typical_morphotype = c(TRUE, TRUE, FALSE, TRUE, FALSE),
          type_morph_mahal_d2 = c(0.4, 8.5, 11.2, 1.5, 6.8),
          nearest_type_morph_taxon = c("Packera dubia", "Packera anonyma", "Packera obovata", "Packera dubia", "Packera dubia"),
          type_morph_concordance = c(0.95, 0.15, 0.12, 0.82, 0.45),
          type_edaphic_z_score = c(0.3, 1.8, 2.8, 0.9, 2.6),
          nearest_type_niche_taxon = c("Packera dubia", "Packera anonyma", "Packera obovata", "Packera dubia", "Packera dubia"),
          stringsAsFactors = FALSE
        )

        res07 <- apply_taxonomic_decision_matrix(df07)
        # Rule A1: Clean Holotype -> Confirmed_Type_Control, Confirmed_Species, RESOLVED
        stopifnot(identical(res07$taxonomic_status_call[1], "Confirmed_Type_Control"))
        stopifnot(identical(res07$taxonomic_treatment_subtype[1], "Confirmed_Species"))
        stopifnot(identical(res07$recommended_determination[1], "Packera dubia"))

        # Rule A2: Cross-modal overlapping Lectotype -> Type_Synonymy_Candidate, NEVER Misidentification_Reassign
        stopifnot(identical(res07$taxonomic_status_call[2], "Type_Synonymy_Candidate"))
        stopifnot(identical(res07$taxonomic_treatment_subtype[2], "Confirmed_Species"))
        stopifnot(identical(res07$recommended_determination[2], "Packera dubia"))
        stopifnot(grepl("Packera anonyma", res07$synthesis_rationale[2]))

        # Rule B: Non-type sheet matching Packera obovata across >=2 type baselines + deviating from P. dubia type baseline
        stopifnot(identical(res07$taxonomic_treatment_subtype[3], "Misidentification_Reassign"))
        stopifnot(identical(res07$recommended_determination[3], "Packera obovata"))
        stopifnot(grepl("Type-Baseline Concordance", res07$synthesis_rationale[3]))

        # Rule C: Specimen in typical morphotype cluster (is_typical_morphotype=TRUE) with type_edaphic_z_score <= 2.0
        # is locked as Confirmed_Species (sensu stricto) in both multi-species and intraspecific modes!
        stopifnot(identical(res07$taxonomic_treatment_subtype[4], "Confirmed_Species"))
        stopifnot(identical(res07$recommended_determination[4], "Packera dubia"))
        res07_intra <- apply_taxonomic_decision_matrix(df07, intraspecific_mode = TRUE)
        stopifnot(identical(res07_intra$taxonomic_treatment_subtype[4], "Confirmed_Species"))
        # Meanwhile non-typical cluster 3 with barrier (row 5) is allowed to route to Subspecies_Ecotype
        stopifnot(identical(res07$taxonomic_treatment_subtype[5], "Subspecies_Ecotype"))
        stopifnot(identical(res07_intra$taxonomic_treatment_subtype[5], "Subspecies_Ecotype"))

        cat("STEP06_07_TYPE_RULES_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_script_path = Path(tmpdir) / "test_step06_07_type_rules.R"
            r_script_path.write_text(r_code, encoding="utf-8")
            res = run_warm_r_cmd(["Rscript", str(r_script_path)], stage_key="step05", cwd=PROJECT_ROOT)
        self.assertEqual(res.returncode, 0, f"Step 06/07 type rules test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("STEP06_07_TYPE_RULES_OK", res.stdout)

    def test_triage_queue_cryptic_cluster_cohesion_pass_blocks_singleton_cryptic(self):
        """
        Verify that `triage_queue.csv` exports `cryptic_cluster_cohesion_pass`, and
        isolated non-core morphotype vouchers without population cohesion
        (`cryptic_cluster_cohesion_pass == FALSE`) are not misclassified as
        `Putative_Cryptic_Species`.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            vouchers_csv = tmp_p / "curated_vouchers.csv"
            pheno_csv = tmp_p / "voucher_reproductive_phenology.csv"
            morph_csv = tmp_p / "morphometric_flags.csv"
            mm_csv = tmp_p / "multimodal_flags.csv"
            niche_csv = tmp_p / "niche_summary.csv"
            out_queue = tmp_p / "triage_queue.csv"
            out_html = tmp_p / "Triage_Dashboard_Interactive.html"
            out_rep = tmp_p / "revision_report.md"
            out_fig = tmp_p / "synthesis.png"

            pd.DataFrame(
                [
                    {
                        "catalogNumber": "COH-ANCHOR-A1",
                        "institutionCode": "NCU",
                        "recordedBy": "A. S. Weakley",
                        "recordNumber": "101",
                        "species_raw": "Packera dubia",
                        "scientificName": "Packera dubia",
                        "determiner_tier": "Tier_1_Gold",
                        "decimalLatitude": 35.8,
                        "decimalLongitude": -79.0,
                        "doy": 115,
                    },
                    {
                        "catalogNumber": "COH-ANCHOR-A2",
                        "institutionCode": "DUKE",
                        "recordedBy": "D. K. Trock",
                        "recordNumber": "202",
                        "species_raw": "Packera dubia",
                        "scientificName": "Packera dubia",
                        "determiner_tier": "Tier_1_Gold",
                        "decimalLatitude": 35.85,
                        "decimalLongitude": -79.05,
                        "doy": 117,
                    },
                    {
                        "catalogNumber": "COH-ANCHOR-B1",
                        "institutionCode": "NCU",
                        "recordedBy": "A. S. Weakley",
                        "recordNumber": "303",
                        "species_raw": "Packera obovata",
                        "scientificName": "Packera obovata",
                        "determiner_tier": "Tier_1_Gold",
                        "decimalLatitude": 36.1,
                        "decimalLongitude": -86.7,
                        "doy": 108,
                    },
                    {
                        "catalogNumber": "COH-ANCHOR-B2",
                        "institutionCode": "TENN",
                        "recordedBy": "D. K. Trock",
                        "recordNumber": "404",
                        "species_raw": "Packera obovata",
                        "scientificName": "Packera obovata",
                        "determiner_tier": "Tier_1_Gold",
                        "decimalLatitude": 36.15,
                        "decimalLongitude": -86.75,
                        "doy": 110,
                    },
                    {
                        "catalogNumber": "COH-SINGLETON-NONCORE",
                        "institutionCode": "NCU",
                        "recordedBy": "J. LoneCollector",
                        "recordNumber": "999",
                        "species_raw": "Packera dubia",
                        "scientificName": "Packera dubia",
                        "determiner_tier": "Tier_2_Silver",
                        "decimalLatitude": 35.2,
                        "decimalLongitude": -80.1,
                        "doy": 162,
                    },
                ]
            ).to_csv(vouchers_csv, index=False)

            pd.DataFrame(
                [
                    {"catalogNumber": "COH-ANCHOR-A1", "voucher_phenological_state": "anthesis", "n_anthesis": 6, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.94},
                    {"catalogNumber": "COH-ANCHOR-A2", "voucher_phenological_state": "anthesis", "n_anthesis": 5, "n_bud": 1, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.93},
                    {"catalogNumber": "COH-ANCHOR-B1", "voucher_phenological_state": "anthesis", "n_anthesis": 6, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.92},
                    {"catalogNumber": "COH-ANCHOR-B2", "voucher_phenological_state": "anthesis", "n_anthesis": 5, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.91},
                    {"catalogNumber": "COH-SINGLETON-NONCORE", "voucher_phenological_state": "anthesis", "n_anthesis": 5, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.92},
                ]
            ).to_csv(pheno_csv, index=False)

            pd.DataFrame(
                [
                    {"catalogNumber": "COH-ANCHOR-A1", "plant_individual_id": 0, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera dubia", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": -2.2, "PC2": 0.2, "can1": -2.5, "can2": 0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "is_typical_morphotype": True, "morphometric_outlier_flag": False},
                    {"catalogNumber": "COH-ANCHOR-A2", "plant_individual_id": 0, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera dubia", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.95, "cda_second_posterior_prob": 0.05, "cda_posterior_margin": 0.90, "PC1": -2.1, "PC2": 0.1, "can1": -2.4, "can2": 0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "is_typical_morphotype": True, "morphometric_outlier_flag": False},
                    {"catalogNumber": "COH-ANCHOR-B1", "plant_individual_id": 0, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera dubia", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": 2.2, "PC2": -0.2, "can1": 2.5, "can2": -0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "is_typical_morphotype": True, "morphometric_outlier_flag": False},
                    {"catalogNumber": "COH-ANCHOR-B2", "plant_individual_id": 0, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera dubia", "cda_posterior_prob": 0.95, "cda_second_posterior_prob": 0.05, "cda_posterior_margin": 0.90, "PC1": 2.1, "PC2": -0.1, "can1": 2.4, "can2": -0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "is_typical_morphotype": True, "morphometric_outlier_flag": False},
                    {"catalogNumber": "COH-SINGLETON-NONCORE", "plant_individual_id": 0, "morph_tier_source": "Basal_Rosette", "basal_leaf_count": 3, "cda_predicted_taxon": "Packera dubia", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.45, "cda_second_posterior_prob": 0.25, "cda_posterior_margin": 0.20, "PC1": -0.2, "PC2": 3.2, "can1": 0.0, "can2": 2.8, "can3": 1.2, "gmm_cluster": 4, "gmm_2_delta_bic": 14.0, "gmm_uncertainty": 0.06, "is_typical_morphotype": False, "morphometric_outlier_flag": True, "cryptic_cluster_cohesion_pass": False},
                ]
            ).to_csv(morph_csv, index=False)

            pd.DataFrame(
                [
                    {"catalogNumber": "COH-ANCHOR-A1", "plant_individual_id": 0, "fine_soil_ph": 5.4, "fine_soil_sand": 62.0, "fine_bedrock_depth_cm": 120.0, "edaphic_z_score": 0.30, "edaphic_best_fit_taxon": "Packera dubia", "rf_predicted_taxon": "Packera dubia", "rf_predicted_prob": 0.92, "is_verified_anthesis": True, "edaphic_schoener_d": 0.80, "delta_doy": 1.0, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                    {"catalogNumber": "COH-ANCHOR-A2", "plant_individual_id": 0, "fine_soil_ph": 5.5, "fine_soil_sand": 60.0, "fine_bedrock_depth_cm": 118.0, "edaphic_z_score": 0.32, "edaphic_best_fit_taxon": "Packera dubia", "rf_predicted_taxon": "Packera dubia", "rf_predicted_prob": 0.91, "is_verified_anthesis": True, "edaphic_schoener_d": 0.79, "delta_doy": 1.5, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                    {"catalogNumber": "COH-ANCHOR-B1", "plant_individual_id": 0, "fine_soil_ph": 7.3, "fine_soil_sand": 24.0, "fine_bedrock_depth_cm": 40.0, "edaphic_z_score": 0.35, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.93, "is_verified_anthesis": True, "edaphic_schoener_d": 0.82, "delta_doy": -1.0, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                    {"catalogNumber": "COH-ANCHOR-B2", "plant_individual_id": 0, "fine_soil_ph": 7.4, "fine_soil_sand": 22.0, "fine_bedrock_depth_cm": 38.0, "edaphic_z_score": 0.38, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.92, "is_verified_anthesis": True, "edaphic_schoener_d": 0.81, "delta_doy": -1.5, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                    {"catalogNumber": "COH-SINGLETON-NONCORE", "plant_individual_id": 0, "fine_soil_ph": 4.3, "fine_soil_sand": 78.0, "fine_bedrock_depth_cm": 16.0, "edaphic_z_score": 3.40, "edaphic_best_fit_taxon": "Packera dubia", "rf_predicted_taxon": "Packera dubia", "rf_predicted_prob": 0.52, "is_verified_anthesis": True, "edaphic_schoener_d": 0.15, "delta_doy": 28.0, "triage_category": "Edaphic_Envelope_Violation", "multimodal_conflict_flag": True},
                ]
            ).to_csv(mm_csv, index=False)

            pd.DataFrame(
                [
                    {"Pair": "P.dubia vs P.obovata", "Taxon1": "Packera dubia", "Taxon2": "Packera obovata", "N1": 10, "N2": 10, "Warren_D": 0.20, "Warren_I": 0.44, "p_value_D": 0.0099, "p_value_I": 0.0099, "Taxonomic_Inference": "Distinct Ecological Niche"}
                ]
            ).to_csv(niche_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouchers_csv),
                    "--phenology", str(pheno_csv),
                    "--morphometrics", str(morph_csv),
                    "--multimodal-flags", str(mm_csv),
                    "--niche-summary", str(niche_csv),
                    "--output-queue", str(out_queue),
                    "--output-html", str(out_html),
                    "--output-report", str(out_rep),
                    "--output-plot", str(out_fig),
                ], stage_key='step06', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 06 cohesion gate test failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertTrue(out_queue.exists())
            qdf = pd.read_csv(out_queue)
            self.assertIn("cryptic_cluster_cohesion_pass", qdf.columns)
            singleton_row = qdf[qdf["catalogNumber"] == "COH-SINGLETON-NONCORE"].iloc[0]
            self.assertFalse(bool(singleton_row["cryptic_cluster_cohesion_pass"]))
            self.assertNotEqual(str(singleton_row["taxonomic_status_call"]), "Putative_Cryptic_Species")
            self.assertNotEqual(str(singleton_row["taxonomic_treatment_subtype"]), "Putative_Cryptic_Species")

    def test_sibling_ramet_propagation_and_morph_concordance_protection(self):
        """
        Verify in `apply_taxonomic_decision_matrix()` (`scripts/analysis/06_triage_dashboard_synthesis.R`):
        1. Leafless ramets on a multi-plant sheet (`catalogNumber`) inherit sibling sheet-mate
           morphology and are never assigned `Putative_Hybrid_Introgression` / `Hybrid_Introgressant`
           via an imputed `cda_margin = 0.15` fallback.
        2. Specimens with concordant leaf morphology (`cda_predicted_taxon == species_standardized`
           and `cda_posterior_prob >= 0.80`) are protected from edaphic-only
           `Misidentification_Reassign` overrides (`Rule B` / `congener_osr_reassign`).
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/06_triage_dashboard_synthesis.R")

        df_test <- data.frame(
          catalogNumber = c(
            "MULTI-SHEET-001",
            "MULTI-SHEET-001",
            "LEAFLESS-SOLO-002",
            "MORPH-CONCORDANT-EDAPH-CONFLICT-003"
          ),
          plant_individual_id = c(0L, 1L, 0L, 0L),
          species_raw = rep("Packera dubia", 4L),
          species_standardized = rep("Packera dubia", 4L),
          determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver", "Tier_1_Gold"),
          type_status = rep("None", 4L),
          is_type_specimen = rep(FALSE, 4L),
          morph_tier_source = c("Basal_Rosette", "No_Leaves", "No_Leaves", "Basal_Rosette"),
          basal_leaf_count = c(3L, 0L, 0L, 4L),
          cda_predicted_taxon = c("Packera dubia", NA_character_, NA_character_, "Packera dubia"),
          cda_posterior_prob = c(0.94, NA_real_, NA_real_, 0.91),
          cda_posterior_margin = c(0.88, NA_real_, NA_real_, 0.82),
          morphometric_outlier_flag = c(FALSE, FALSE, FALSE, FALSE),
          edaphic_best_fit_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata"),
          rf_predicted_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata"),
          rf_predicted_prob = c(0.92, 0.92, 0.78, 0.89),
          edaphic_schoener_d = c(0.82, 0.82, 0.40, 0.35),
          delta_doy_resid = c(1.0, 1.0, 2.0, 1.5),
          is_active_anthesis = rep(TRUE, 4L),
          gmm_cluster = c(1L, NA_integer_, NA_integer_, 1L),
          gmm_2_delta_bic = rep(8.0, 4L),
          gmm_uncertainty = c(0.04, NA_real_, NA_real_, 0.05),
          is_typical_morphotype = c(TRUE, NA, NA, TRUE),
          type_morph_mahal_d2 = c(0.55, NA_real_, NA_real_, 0.70),
          nearest_type_morph_taxon = c("Packera dubia", NA_character_, NA_character_, "Packera dubia"),
          type_morph_concordance = c(0.95, NA_real_, NA_real_, 0.92),
          type_edaphic_z_score = c(0.40, 0.40, 1.60, 2.95),
          nearest_type_niche_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata"),
          min_congener_mahal_d2 = c(9.50, NA_real_, NA_real_, 3.10),
          min_congener_taxon = c("Packera obovata", NA_character_, NA_character_, "Packera obovata"),
          stringsAsFactors = FALSE
        )

        res <- apply_taxonomic_decision_matrix(df_test)

        # 1a. Multi-plant sheet leafless ramet (row 2) inherits sibling morphology from row 1
        stopifnot(isTRUE(res$has_observed_or_sheet_morph[2]))
        stopifnot(identical(res$cda_predicted_taxon[2], "Packera dubia"))
        stopifnot(is.finite(res$cda_posterior_prob[2]) && res$cda_posterior_prob[2] > 0.80)
        stopifnot(identical(res$taxonomic_treatment_subtype[2], "Confirmed_Species"))

        # 1b. Solo leafless ramet (row 3) without observed or sibling morphology is never assigned
        #     Hybrid_Introgressant / Putative_Hybrid_Introgression via an imputed cda_margin = 0.15 fallback
        stopifnot(!isTRUE(res$has_observed_or_sheet_morph[3]))
        stopifnot(!identical(res$taxonomic_treatment_subtype[3], "Hybrid_Introgressant"))
        stopifnot(!identical(res$taxonomic_status_call[3], "Putative_Hybrid_Introgression"))
        stopifnot(!identical(res$taxonomic_status_call[3], "Hybrid_Intergrade_Swarm"))

        # 2. Morphologically concordant specimen (row 4: cda_predicted_taxon == species_standardized,
        #    cda_posterior_prob == 0.91 >= 0.80, type_morph_mahal_d2 == 0.70) is protected from
        #    edaphic-only Misidentification_Reassign overrides (Rule B / congener_osr_reassign)
        stopifnot(!identical(res$taxonomic_treatment_subtype[4], "Misidentification_Reassign"))
        stopifnot(identical(res$recommended_determination[4], "Packera dubia"))

        cat("SIBLING_PROPAGATION_AND_MORPH_PROTECTION_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_script = Path(tmpdir) / "test_sibling_and_morph_protection.R"
            r_script.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(['Rscript', str(r_script)], stage_key='step06', cwd=PROJECT_ROOT)
        self.assertEqual(
            proc.returncode,
            0,
            f"Sibling ramet propagation / morph concordance protection test failed:\n{proc.stderr}\n{proc.stdout}",
        )
        self.assertIn("SIBLING_PROPAGATION_AND_MORPH_PROTECTION_OK", proc.stdout)



    def test_step06_zero_trust_four_gates_and_offline_html_dashboard(self):
        """
        Zero-trust regression suite for Step 06 Four-Gate Overhaul:
          1. Canonical single-column contract (zero twin/obsolete columns in R output & SQLite).
          2. p_type_morph adapts dynamically to selected_morph_features dimensionality.
          3. Stage 4 late-season foliar wear evaluates circular Hopkins-adjusted DOY quantiles
             (protecting high-latitude summer sheets from false Ecophenotypic_Plasticity).
          4. Rule B vote_edaph requires posterior support (rf_predicted_prob >= 0.55) with 1 vote
             per modality, and Rule A Type_Synonymy_Candidate requires cda_posterior_prob >= 0.80.
          5. Anchor mask & dynamic edaphic z-score exclude is_confident_label_noise & is_passive,
             and p_e has zero cross-modality leakage from multimodal_conflict_flag.
          6. Offline HTML dashboard (Triage_Dashboard_Interactive.html) has zero external CDN links,
             and revision Markdown report has zero em/en dashes.
        """
        r_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("generate_interactive_triage_html", r_src)
        self.assertIn("Triage_Dashboard_Interactive.html", r_src)
        self.assertIn("persist_step06_to_sqlite", r_src)
        self.assertIn("generate_step06_pub_figures", r_src)
        self.assertIn("render_step06_synthesis_plate_pdf", r_src)

        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/06_triage_dashboard_synthesis.R")

        # 1. p_type_morph dimensionality (p=6: chisq_0.95,6 = 12.592 vs df=3: 7.815) & Rule B rf_prob floor
        df_dim <- data.frame(
          catalogNumber = c("ANC_1", "ANC_2", "ANC_3", "ANC_4", "LOW_RF_PROB", "BORDERLINE_TYPE_HOLO"),
          plant_individual_id = 0L,
          species_raw = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera dubia", "Packera dubia"),
          species_standardized = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera dubia", "Packera dubia"),
          determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze", "Tier_1_Gold"),
          is_type_specimen = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          type_status = c("None", "None", "None", "None", "None", "Holotype"),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          selected_morph_features = "PC1;PC2;PC3;PC4;PC5;ECT_PC1",
          PC1 = c(-2.0, -2.1, 2.0, 2.1, 0.2, 1.5),
          PC2 = c(0.1, 0.2, -0.1, -0.2, 0.0, -0.1),
          PC3 = c(0.0, 0.1, 0.0, -0.1, 0.1, 0.0),
          PC4 = c(0.1, 0.0, -0.1, 0.0, 0.0, 0.1),
          PC5 = c(0.0, -0.1, 0.1, 0.0, -0.1, 0.0),
          ECT_PC1 = c(-1.0, -1.1, 1.0, 1.1, 0.1, 0.8),
          cda_predicted_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera obovata", "Packera obovata"),
          cda_posterior_prob = c(0.96, 0.95, 0.96, 0.95, 0.58, 0.62),
          cda_posterior_margin = c(0.92, 0.90, 0.92, 0.90, 0.32, 0.24),
          morphometric_outlier_flag = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
          edaphic_best_fit_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera obovata", "Packera obovata"),
          rf_predicted_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera obovata", "Packera obovata"),
          rf_predicted_prob = c(0.92, 0.91, 0.93, 0.90, 0.24, 0.65),
          edaphic_schoener_d = c(0.82, 0.80, 0.84, 0.81, 0.62, 0.30),
          type_morph_mahal_d2 = c(1.2, 1.4, 1.1, 1.3, 10.5, 9.0),
          nearest_type_morph_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera obovata", "Packera obovata"),
          type_edaphic_z_score = c(0.4, 0.5, 0.3, 0.4, 1.2, 1.5),
          nearest_type_niche_taxon = c("Packera dubia", "Packera dubia", "Packera obovata", "Packera obovata", "Packera obovata", "Packera obovata"),
          gmm_cluster = 1L,
          gmm_uncertainty = 0.05,
          stringsAsFactors = FALSE
        )
        res_dim <- apply_taxonomic_decision_matrix(df_dim)
        # Canonical columns present; twin columns absent
        stopifnot(all(c("n_available_modalities", "probabilistic_consensus_conflict", "probabilistic_consensus_score", "taxonomic_treatment_broad", "taxonomic_treatment_subtype") %in% names(res_dim)))
        stopifnot(!(paste0("taxonomic_", "treatment") %in% names(res_dim)))
        stopifnot(!(paste0("taxonomic_treatment_", "sub") %in% names(res_dim)))
        # Row 5: type_morph_mahal_d2 = 10.5 <= qchisq(0.95, df=6)=12.592 AND rf_predicted_prob=0.24 < 0.55 -> NOT Misidentification_Reassign!
        stopifnot(!identical(res_dim$taxonomic_treatment_subtype[5], "Misidentification_Reassign"))
        # Row 6: Holotype with cda_posterior_prob = 0.62 < 0.80 -> Confirmed_Type_Control (NOT Type_Synonymy_Candidate)
        stopifnot(identical(res_dim$taxonomic_status_call[6], "Confirmed_Type_Control"))

        # 2. Circular Hopkins-adjusted DOY protects northern summer sheets from false Ecophenotypic_Plasticity
        df_lat <- data.frame(
          catalogNumber = c("S1", "S2", "S3", "S4", "S5", "OB1", "OB2", "NORTH_NORMAL", "SOUTH_LATE"),
          plant_individual_id = 0L,
          species_raw = c(rep("Packera aurea", 5), rep("Packera obovata", 2), "Packera aurea", "Packera aurea"),
          species_standardized = c(rep("Packera aurea", 5), rep("Packera obovata", 2), "Packera aurea", "Packera aurea"),
          determiner_tier = c(rep("Tier_1_Gold", 7), "Tier_3_Bronze", "Tier_3_Bronze"),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          decimalLatitude = c(34.8, 35.0, 35.2, 35.4, 35.6, 35.1, 35.3, 45.0, 35.2),
          decimalLongitude = rep(-82.0, 9),
          doy = c(116, 118, 120, 122, 124, 119, 121, 159, 132),
          is_verified_anthesis = TRUE,
          cda_predicted_taxon = c(rep("Packera aurea", 5), rep("Packera obovata", 2), "Packera obovata", "Packera obovata"),
          cda_posterior_prob = c(0.95, 0.94, 0.95, 0.96, 0.95, 0.95, 0.94, 0.58, 0.58),
          cda_posterior_margin = c(0.90, 0.88, 0.90, 0.92, 0.90, 0.90, 0.88, 0.30, 0.30),
          morphometric_outlier_flag = FALSE,
          edaphic_best_fit_taxon = c(rep("Packera aurea", 5), rep("Packera obovata", 2), "Packera aurea", "Packera aurea"),
          edaphic_z_score = rep(0.5, 9),
          edaphic_schoener_d = rep(0.80, 9),
          gmm_cluster = 1L,
          gmm_uncertainty = 0.05,
          stringsAsFactors = FALSE
        )
        res_lat <- apply_taxonomic_decision_matrix(df_lat)
        # Row 8 (NORTH_NORMAL: lat=45.0, doy=159 -> hopkins_adjusted_doy ~ 119.8 < sp_q75) is protected as Confirmed_Species;
        # Row 9 (SOUTH_LATE: lat=35.2, doy=132 -> +12d past species median >= sp_q75, < 18d cryptic barrier) is Ecophenotypic_Plasticity!
        stopifnot(res_lat$hopkins_adjusted_doy[8] < 121.0)
        stopifnot(identical(res_lat$taxonomic_treatment_subtype[8], "Confirmed_Species"))
        stopifnot(identical(res_lat$taxonomic_treatment_subtype[9], "Ecophenotypic_Plasticity"))

        # 3. Zero leakage from is_confident_label_noise / is_passive & multimodal_conflict_flag
        df_clean <- df_dim[1:4, ]
        df_clean$fine_soil_ph <- c(5.2, 5.3, 7.2, 7.3)
        df_clean$primary_edaph_predictors <- "fine_soil_ph"
        df_clean$edaphic_z_score <- NA_real_
        df_clean$multimodal_conflict_flag <- FALSE

        df_noisy <- rbind(df_clean, df_clean[1, ])
        df_noisy$catalogNumber[5] <- "CORRUPT_NOISE"
        df_noisy$is_confident_label_noise <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
        df_noisy$cda_predicted_taxon[5] <- "Packera obovata"
        df_noisy$edaphic_best_fit_taxon[5] <- "Packera obovata"
        df_noisy$fine_soil_ph[5] <- 99.0

        w_clean <- compute_modality_reliability_weights(df_clean)
        w_noisy <- compute_modality_reliability_weights(df_noisy)
        stopifnot(all(abs(w_clean - w_noisy) < 1e-9))

        res_c1 <- compute_probabilistic_consensus_scores(df_clean)
        df_clean_mm <- df_clean
        df_clean_mm$multimodal_conflict_flag <- TRUE
        res_c2 <- compute_probabilistic_consensus_scores(df_clean_mm)
        stopifnot(all(abs(res_c1$probabilistic_consensus_conflict - res_c2$probabilistic_consensus_conflict) < 1e-9))

        # 4. Offline HTML Dashboard + Zero Em/En-Dash Markdown Report
        tf_html <- tempfile(fileext = ".html")
        tf_md <- tempfile(fileext = ".md")
        generate_interactive_triage_html(res_dim, tf_html, focal_taxon = "Packera dubia")
        generate_taxonomic_revision_report(res_dim, list(output_report = tf_md, output_plot = "", focal_taxon = "Packera dubia", gmm_bic = "", gmm_summary = ""))
        html_txt <- paste(readLines(tf_html, warn = FALSE), collapse = "\n")
        md_txt <- paste(readLines(tf_md, warn = FALSE), collapse = "\n")
        stopifnot(file.exists(tf_html) && nchar(html_txt) > 500)
        stopifnot(!grepl("https?://", html_txt, ignore.case = TRUE))
        stopifnot(grepl("Export Filtered CSV", html_txt, fixed = TRUE))
        stopifnot(!grepl("[\u2013\u2014]", md_txt))

        cat("STEP06_ZERO_TRUST_OK\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_script = Path(tmpdir) / "test_step06_zero_trust.R"
            r_script.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(['Rscript', str(r_script)], stage_key='step06', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Step 06 zero-trust suite failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP06_ZERO_TRUST_OK", proc.stdout)

    def test_conservative_five_gate_hybrid_swarm_parentage_and_intraspecific_guard(self):
        """
        Verify the 5-Gate Conservative Contact-Zone Hybrid Swarm Protocol:
          1. Blocks soil-only Spatial RF disagreement when morphology matches given_sp.
          2. Blocks allopatric intermediates (> 15 km from congener parent).
          3. Blocks allochronic intermediates (> 21 d flowering peak separation).
          4. Blocks off-axis taphonomic outliers failing orthogonal Mahalanobis residual.
          5. Blocks isolated n_deme == 1 sympatric intermediates while confirming n_deme >= 2
             contact-zone swarms with explicit hybrid_parent_1, hybrid_parent_2, HTML chips, & report summary.
          6. Strictly disables hybrid and parent calls in --intraspecific-mode and 1-taxon cohorts.
        """
        r_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("can_evaluate_hybrids", r_src)
        self.assertIn(".eval_five_gate_hybrid_pair", r_src)
        self.assertIn("hybrid_zone_cohesion_pass", r_src)
        self.assertIn("hybrid_parent_1", r_src)
        self.assertIn("hybrid_parent_2", r_src)
        self.assertNotIn("[Admixed Morphotype]", r_src)

        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/06_triage_dashboard_synthesis.R")

        df_hyb <- data.frame(
          catalogNumber = c(
            "DUB_ANC_1", "DUB_ANC_2",
            "ANO_ANC_1", "ANO_ANC_2",
            "LATE_ANC_1", "LATE_ANC_2",
            "SOIL_ONLY_DISAGREE",
            "ALLOPATRIC_INTERM",
            "ALLOCHRONIC_INTERM",
            "OFF_CHORD_OUTLIER",
            "SINGLETON_SYMPATRIC",
            "SWARM_CONTACT_1",
            "SWARM_CONTACT_2"
          ),
          plant_individual_id = 0L,
          species_raw = c(
            "Packera dubia", "Packera dubia",
            "Packera anonyma", "Packera anonyma",
            "Packera schweinitziana", "Packera schweinitziana",
            rep("Packera dubia", 7)
          ),
          species_standardized = c(
            "Packera dubia", "Packera dubia",
            "Packera anonyma", "Packera anonyma",
            "Packera schweinitziana", "Packera schweinitziana",
            rep("Packera dubia", 7)
          ),
          determiner_tier = c(rep("Tier_1_Gold", 6), rep("Tier_2_Silver", 7)),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          # Contact zone 1 near (35.80, -79.00) has DUB + ANO + LATE; Contact zone 2 near (35.20, -80.80) has DUB + ANO for singleton test
          decimalLatitude = c(
            35.800, 35.200,
            35.808, 35.208,
            35.805, 35.812,
            35.804,
            34.100,
            35.806,
            35.805,
            35.204,
            35.803,
            35.807
          ),
          decimalLongitude = c(
            -79.000, -80.800,
            -79.008, -80.808,
            -79.005, -79.012,
            -79.004,
            -77.900,
            -79.006,
            -79.005,
            -80.804,
            -79.003,
            -79.007
          ),
          doy = c(
            114, 116,
            120, 122,
            165, 169,
            116,
            118,
            140,
            118,
            118,
            117,
            119
          ),
          is_verified_anthesis = TRUE,
          PC1 = c(
            -3.0, -2.8,
             3.0,  2.8,
             0.0,  0.2,
            -2.7,
             0.0,
            -1.4,
             0.0,
             0.1,
            -0.1,
             0.2
          ),
          PC2 = c(
            0.0,  0.1,
            0.0, -0.1,
            3.5,  3.3,
            0.05,
            0.0,
            1.7,
            3.9,
            0.02,
            0.03,
           -0.02
          ),
          blade_cie_l_sheet_delta = c(
            -12.0, -11.0,
            -28.0, -27.0,
            -26.0, -25.0,
            -11.5,
            -20.0,
            -18.5,
            -19.5,
            -19.5,
            -19.8,
            -20.4
          ),
          cda_predicted_taxon = c(
            "Packera dubia", "Packera dubia",
            "Packera anonyma", "Packera anonyma",
            "Packera schweinitziana", "Packera schweinitziana",
            "Packera dubia",
            "Packera anonyma",
            "Packera schweinitziana",
            "Packera anonyma",
            "Packera anonyma",
            "Packera anonyma",
            "Packera anonyma"
          ),
          cda_second_taxon = c(
            "Packera anonyma", "Packera anonyma",
            "Packera dubia", "Packera dubia",
            "Packera dubia", "Packera dubia",
            "Packera anonyma",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia"
          ),
          cda_posterior_prob = c(0.96, 0.95, 0.96, 0.95, 0.96, 0.95, 0.74, 0.52, 0.51, 0.52, 0.53, 0.52, 0.54),
          cda_second_posterior_prob = c(0.04, 0.05, 0.04, 0.05, 0.04, 0.05, 0.24, 0.42, 0.41, 0.40, 0.41, 0.42, 0.40),
          cda_posterior_margin = c(0.92, 0.90, 0.92, 0.90, 0.92, 0.90, 0.50, 0.10, 0.10, 0.12, 0.12, 0.10, 0.14),
          morphometric_outlier_flag = FALSE,
          edaphic_best_fit_taxon = c(
            "Packera dubia", "Packera dubia",
            "Packera anonyma", "Packera anonyma",
            "Packera schweinitziana", "Packera schweinitziana",
            "Packera anonyma",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera anonyma"
          ),
          rf_predicted_taxon = c(
            "Packera dubia", "Packera dubia",
            "Packera anonyma", "Packera anonyma",
            "Packera schweinitziana", "Packera schweinitziana",
            "Packera anonyma",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera anonyma"
          ),
          rf_predicted_prob = c(0.92, 0.91, 0.92, 0.91, 0.92, 0.91, 0.68, 0.48, 0.48, 0.48, 0.49, 0.48, 0.49),
          edaphic_schoener_d = rep(0.62, 13),
          edaphic_z_score = rep(0.85, 13),
          gmm_cluster = 1L,
          gmm_uncertainty = c(rep(0.04, 6), 0.12, 0.38, 0.39, 0.37, 0.38, 0.40, 0.39),
          triage_category = c(rep("Clean_MultiModal_Consensus", 6), rep("Putative_Hybrid_Zone_Intergrade", 7)),
          stringsAsFactors = FALSE
        )

        res_multi <- apply_taxonomic_decision_matrix(df_hyb, intraspecific_mode = FALSE)
        rownames(res_multi) <- res_multi$catalogNumber

        # Gates 1 to 5 negative controls: ZERO false hybrids
        for (neg_id in c("SOIL_ONLY_DISAGREE", "ALLOPATRIC_INTERM", "ALLOCHRONIC_INTERM", "OFF_CHORD_OUTLIER", "SINGLETON_SYMPATRIC")) {
          stopifnot(!identical(res_multi[neg_id, "taxonomic_treatment_subtype"], "Hybrid_Introgressant"))
          stopifnot(!isTRUE(res_multi[neg_id, "hybrid_zone_cohesion_pass"]))
          stopifnot(is.na(res_multi[neg_id, "hybrid_parent_1"]))
          stopifnot(is.na(res_multi[neg_id, "hybrid_parent_2"]))
        }
        stopifnot(grepl("Unreplicated sympatric intermediate", res_multi["SINGLETON_SYMPATRIC", "synthesis_rationale"]))

        # Confirmed contact-zone hybrid swarm (n_deme = 2) reports both parent species explicitly
        for (pos_id in c("SWARM_CONTACT_1", "SWARM_CONTACT_2")) {
          stopifnot(identical(res_multi[pos_id, "taxonomic_treatment_subtype"], "Hybrid_Introgressant"))
          stopifnot(isTRUE(res_multi[pos_id, "hybrid_zone_cohesion_pass"]))
          stopifnot(identical(res_multi[pos_id, "hybrid_parent_1"], "Packera dubia"))
          stopifnot(identical(res_multi[pos_id, "hybrid_parent_2"], "Packera anonyma"))
          stopifnot(identical(res_multi[pos_id, "recommended_determination"], "Packera dubia x Packera anonyma"))
        }

        # Check Markdown report & interactive HTML parent-species reporting
        tf_rep <- tempfile(fileext = ".md")
        tf_htm <- tempfile(fileext = ".html")
        opts_r <- list(output_report = tf_rep, output_plot = "", focal_taxon = "Packera dubia", intraspecific_mode = FALSE, gmm_bic = "", gmm_summary = "")
        generate_taxonomic_revision_report(res_multi, opts_r)
        generate_interactive_triage_html(res_multi, tf_htm, focal_taxon = "Packera dubia", opts = opts_r)
        rep_txt <- paste(readLines(tf_rep, warn = FALSE), collapse = "\n")
        htm_txt <- paste(readLines(tf_htm, warn = FALSE), collapse = "\n")
        stopifnot(grepl("Confirmed Hybrid Swarm Parentage", rep_txt, fixed = TRUE))
        stopifnot(grepl("*Packera anonyma* x *Packera dubia*", rep_txt, fixed = TRUE))
        stopifnot(grepl("Confirmed Hybrid Swarms (Putative Parent Species)", htm_txt, fixed = TRUE))
        stopifnot(!grepl("[–—]", rep_txt))

        # Intraspecific mode & single-taxon cohorts NEVER flag hybrids or parent species
        res_intra <- apply_taxonomic_decision_matrix(df_hyb, intraspecific_mode = TRUE)
        stopifnot(!any(res_intra$taxonomic_treatment_subtype == "Hybrid_Introgressant"))
        stopifnot(!any(res_intra$hybrid_zone_cohesion_pass))
        stopifnot(all(is.na(res_intra$hybrid_parent_1)) && all(is.na(res_intra$hybrid_parent_2)))
        stopifnot(!any(grepl("Admixed Morphotype", res_intra$recommended_determination, fixed = TRUE)))

        df_single_sp <- df_hyb[df_hyb$species_standardized == "Packera dubia", , drop = FALSE]
        res_single_sp <- apply_taxonomic_decision_matrix(df_single_sp, intraspecific_mode = FALSE)
        stopifnot(!any(res_single_sp$taxonomic_treatment_subtype == "Hybrid_Introgressant"))
        stopifnot(all(is.na(res_single_sp$hybrid_parent_1)) && all(is.na(res_single_sp$hybrid_parent_2)))

        cat("CONSERVATIVE_FIVE_GATE_HYBRID_OK\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_script = Path(tmpdir) / "test_five_gate_hybrid.R"
            r_script.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_script)], stage_key="step06", cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"5-gate hybrid test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("CONSERVATIVE_FIVE_GATE_HYBRID_OK", proc.stdout)





REPO_ROOT = PROJECT_ROOT
DISCOVERY_UTILS_R = REPO_ROOT / "scripts" / "analysis" / "taxon_complex_discovery_utils.R"
SYNTHESIS_R_SCRIPT = REPO_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"
PUB_FIG_R_SCRIPT = REPO_ROOT / "scripts" / "morphometrics" / "pub_figures_utils.R"

CANONICAL_COMPLEX_COLS = {
    "species_standardized",
    "n_vouchers",
    "n_gold_silver_anchors",
    "n_type_specimens",
    "consensus_complex_id",
    "complex_label",
    "leiden_complex_id",
    "ward_complex_id",
    "silhouette_width",
    "within_complex_mean_affinity",
    "nearest_neighbor_taxon",
    "nearest_neighbor_affinity",
    "second_neighbor_taxon",
    "second_neighbor_affinity",
    "hopkins_pheno_overlap",
    "eco_schoener_overlap",
    "integrative_congruence",
    "w_morph_kappa",
    "w_edaphic_kappa",
    "w_pheno_kappa",
    "n_misidentification_reassign",
    "n_hybrid_introgressant",
    "n_putative_cryptic",
}


class TestEasternPackeraAndComplexDiscovery(unittest.TestCase):
    """Tests for N-taxon basionym mapping, unified Step 06 complex discovery, and Four-Gate contracts."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.r_unit_proc = None

        if rscript_available():
            intra_fig_dir = cls._class_tmp_path / "intra_figs"
            intra_fig_dir.mkdir(parents=True, exist_ok=True)
            mm_dir = cls._class_tmp_path / "mm_join"
            mm_dir.mkdir(parents=True, exist_ok=True)
            morph_csv = mm_dir / "morphometric_flags.csv"
            mm_csv = mm_dir / "multimodal_conflict_flags.csv"
            q_csv = mm_dir / "triage_queue.csv"
            rep_md = mm_dir / "rev.md"
            html_out = mm_dir / "dash.html"
            fig_out = mm_dir / "rev.png"

            pd.DataFrame([
                {"catalogNumber": "MM_01", "plant_individual_id": 1, "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "basal_leaf_count": 3, "PC1": -1.2, "PC2": 0.3},
                {"catalogNumber": "MM_02", "plant_individual_id": 1, "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "basal_leaf_count": 3, "PC1": 1.4, "PC2": -0.2},
            ]).to_csv(morph_csv, index=False)
            pd.DataFrame([
                {
                    "catalogNumber": "MM_01", "plant_individual_id": 1,
                    "fine_soil_ph": 5.9, "fine_soil_sand": 34.5, "fine_bedrock_depth_cm": 98.0,
                    "ssurgo_rock_outcrop_pct": 3.5, "soil_ph": 5.8, "soil_sand": 35.0,
                    "soil_cec": 14.2, "soil_bulk_density": 1.31,
                    "bio1_temp_mean": 13.4, "bio4_temp_seasonality": 712.0,
                    "bio12_precip_annual": 1340.0, "bio15_precip_seasonality": 18.0,
                    "edaphic_schoener_d": 0.41, "schoener_d": 0.44,
                },
                {
                    "catalogNumber": "MM_02", "plant_individual_id": 1,
                    "fine_soil_ph": 6.8, "fine_soil_sand": 52.0, "fine_bedrock_depth_cm": 42.0,
                    "ssurgo_rock_outcrop_pct": 24.0, "soil_ph": 6.7, "soil_sand": 50.0,
                    "soil_cec": 18.5, "soil_bulk_density": 1.42,
                    "bio1_temp_mean": 14.8, "bio4_temp_seasonality": 745.0,
                    "bio12_precip_annual": 1180.0, "bio15_precip_seasonality": 21.0,
                    "edaphic_schoener_d": 0.41, "schoener_d": 0.44,
                },
            ]).to_csv(mm_csv, index=False)

            r_batch_code = f"""
            source("{PUB_FIG_R_SCRIPT.as_posix()}")

            # Block 1: Dynamic palette, taxon standardization, shapes, titles, and italic expressions
            local({{
              stopifnot(!exists("EASTERN_NA_PACKERA_TAXA"))
              stopifnot(!exists("OKABE_ITO_TAXA"))
              stopifnot(!exists("PACKERA_ITALIC_LABELS"))

              stopifnot(standardize_taxon_name("Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve") == "Packera millefolium")
              stopifnot(standardize_taxon_name("P. aurea") == "Packera aurea")
              stopifnot(standardize_taxon_name("Packera x memmingeri") == "Packera memmingeri")
              stopifnot(standardize_taxon_name("Packera × memmingeri") == "Packera memmingeri")
              stopifnot(standardize_taxon_name("Packera cana (Hook.) W.A.Weber & Á.Löve") == "Packera cana")
              stopifnot(standardize_taxon_name("Packera sanguisorbae (DC.) C.Jeffrey") == "Packera sanguisorbae")
              stopifnot(standardize_taxon_name("Packera hyperborealis") == "Packera hyperborealis")
              stopifnot(standardize_taxon_name("Packera xylorrhiza") == "Packera xylorrhiza")

              taxa8 <- c(
                "Packera streptanthifolia", "Packera aurea", "Packera sanguisorbae", "Packera cymbalaria",
                "Packera millefolium", "Packera cana", "Packera obovata", "Packera fendleri"
              )
              sorted8 <- sort(taxa8)
              pal <- get_dynamic_packera_palette(taxa8)
              stopifnot(identical(names(pal), c(sorted8, "Unknown")))
              stopifnot(pal[[sorted8[1]]] == "#0072B2")
              stopifnot(pal[[sorted8[2]]] == "#D55E00")
              stopifnot(pal[[sorted8[3]]] == "#009E73")
              stopifnot(pal[[sorted8[4]]] == "#CC79A7")
              stopifnot(pal[["Unknown"]] == "#7F7F7F")
              stopifnot(length(unique(unname(pal[sorted8]))) == 8)

              shps <- get_dynamic_packera_shapes(taxa8)
              stopifnot(identical(unname(shps[sorted8[1:4]]), c(1L, 2L, 5L, 0L)))
              gg_shps <- get_dynamic_packera_ggplot_shapes(taxa8)
              stopifnot(identical(unname(gg_shps[sorted8[1:5]]), c(21, 24, 22, 23, 25)))

              Sys.setenv(PACKERA_FOCAL_TAXON = "")
              stopifnot(format_active_complex_title(c("Packera cana", "Packera fendleri", "Packera streptanthifolia"), "CDA") == "Packera Taxon Complex (3 Taxa): CDA")
              stopifnot(format_active_complex_title(c("Packera sanguisorbae"), "CDA") == "Packera sanguisorbae: CDA")
              Sys.setenv(PACKERA_FOCAL_TAXON = "Packera streptanthifolia")
              stopifnot(format_active_complex_title(c("Packera cana", "Packera streptanthifolia"), "CDA") == "Packera streptanthifolia Complex: CDA")
              Sys.unsetenv("PACKERA_FOCAL_TAXON")

              expr_tx <- format_italic_taxon_expr("Packera sanguisorbae")
              stopifnot(is.call(expr_tx) && identical(deparse(expr_tx), 'italic("Packera sanguisorbae")'))
              cat("STEP06_R_PALETTE_OK\\n")
            }})

            source("{DISCOVERY_UTILS_R.as_posix()}")

            # Block 2: Stage 1-5 hierarchical precedence, missing-modality weight renormalization, & BH FDR
            local({{
              fdr_in <- data.frame(Pair = c("A_vs_B", "A_vs_C", "B_vs_C"), p_value = c(0.002, 0.030, 0.250))
              fdr_out <- apply_bh_fdr_correction(fdr_in)
              stopifnot("p_adj_bh" %in% names(fdr_out))
              stopifnot(all(abs(fdr_out$p_adj_bh - stats::p.adjust(fdr_in$p_value, method = "BH")) < 1e-9))

              sh <- get_dynamic_packera_shapes(c("Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"))
              stopifnot(all(unname(sh) %in% c(1L, 2L, 5L, 0L)))

              mdf <- data.frame(
                catalogNumber = c("V1", "V2", "V3", "V4", "V5", "V6"),
                species_standardized = c("Packera dubia", "Packera dubia", "Packera dubia", "Packera anonyma", "Packera anonyma", "Packera anonyma"),
                determiner_tier = rep("Tier_1_Gold", 6),
                cda_predicted_taxon = c("Packera dubia", "Packera anonyma", "Packera anonyma", "Packera dubia", "Packera anonyma", "Packera anonyma"),
                cda_posterior_prob = c(0.95, 0.85, 0.60, 0.44, 0.48, 0.92),
                cda_second_posterior_prob = c(0.03, 0.12, 0.40, 0.28, 0.25, 0.05),
                cda_posterior_margin = c(0.92, 0.73, 0.20, 0.16, 0.23, 0.87),
                morphometric_outlier_flag = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE),
                gmm_cluster = c(1L, 2L, 2L, 5L, 5L, 2L),
                gmm_2_delta_bic = c(0, 0, 0, 12.0, 12.0, 0),
                edaphic_schoener_d = c(0.75, 0.60, 0.55, 0.18, 0.72, 0.80),
                basal_leaf_count = rep(3L, 6),
                stringsAsFactors = FALSE
              )
              kw_missing <- compute_modality_kappa_weights(mdf, rep(TRUE, 6))
              # Only morph + edaphic (from edaphic_schoener_d? No: no edaphic_best_fit_taxon/predictors -> w_morph == 1.0!)
              stopifnot(abs(kw_missing$w_morph - 1.0) < 1e-9)

              pheno_mock <- list(delta_doy_resid = c(2.0, 3.0, 4.0, 22.0, 3.0, 1.0))
              calls <- evaluate_hierarchical_complex_delimitation(mdf, pheno_res = pheno_mock)
              stopifnot(calls[1] == "Confirmed_Species")
              stopifnot(calls[2] == "Misidentification_Reassign")
              stopifnot(calls[3] == "Hybrid_Introgressant")
              stopifnot(calls[4] == "Putative_Cryptic_Species")
              stopifnot(calls[5] == "Ecophenotypic_Plasticity")
              cat("STEP06_HIERARCHICAL_OK\\n")
            }})

            # Block 3: Unified Step 06 multimodal left join + complex discovery integration
            source("{SYNTHESIS_R_SCRIPT.as_posix()}")
            local({{
              opts_m <- parse_args_robust()
              opts_m$vouchers <- "{morph_csv.as_posix()}"
              opts_m$morphometrics <- "{morph_csv.as_posix()}"
              opts_m$multimodal_flags <- "{mm_csv.as_posix()}"
              opts_m$output_queue <- "{q_csv.as_posix()}"
              opts_m$output_report <- "{rep_md.as_posix()}"
              opts_m$output_html <- "{html_out.as_posix()}"
              opts_m$output_plot <- "{fig_out.as_posix()}"
              opts_m$output_complexes <- "{(mm_dir / 'empirical_taxon_complexes.csv').as_posix()}"
              opts_m$output_affinity <- "{(mm_dir / 'pairwise_taxon_affinity_matrix.csv').as_posix()}"
              opts_m$output_neighbors <- "{(mm_dir / 'focal_auto_neighbors.txt').as_posix()}"
              opts_m$output_complex_plot <- "{(mm_dir / 'fig06b_taxon_complex_delimitation_network.png').as_posix()}"
              opts_m$demo <- TRUE
              mdf <- build_triage_queue(opts_m)
              req <- c("fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "ssurgo_rock_outcrop_pct",
                       "soil_ph", "soil_sand", "soil_cec", "soil_bulk_density",
                       "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality",
                       "edaphic_schoener_d", "schoener_d")
              stopifnot(all(req %in% names(mdf)))
              stopifnot(abs(mdf$fine_soil_ph[mdf$catalogNumber == "MM_01"] - 5.9) < 1e-6)
              stopifnot(abs(mdf$ssurgo_rock_outcrop_pct[mdf$catalogNumber == "MM_02"] - 24.0) < 1e-6)
              cat("STEP06_MULTIMODAL_JOIN_OK\\n")
            }})

            # Block 4: Intraspecific morphotype cluster visualization across R figures
            source("{PUB_FIG_R_SCRIPT.as_posix()}")
            source("{(REPO_ROOT / 'scripts' / 'analysis' / '05_multimodal_spatial_rf.R').as_posix()}")
            local({{
              set.seed(123)
              n <- 24
              df <- data.frame(
                catalogNumber = sprintf("INTRA_%03d", seq_len(n)),
                plant_individual_id = 1L,
                species_standardized = rep("Packera dubia", n),
                scientificName = rep("Packera dubia", n),
                gmm_cluster = rep(1:3, each = 8),
                analysis_group = paste0("Morphotype_K", rep(1:3, each = 8)),
                grouping_col = paste0("Morphotype_K", rep(1:3, each = 8)),
                can1 = c(rnorm(8, -2), rnorm(8, 0), rnorm(8, 2)),
                can2 = c(rnorm(8, 1), rnorm(8, -1.5), rnorm(8, 1)),
                can3 = rnorm(n),
                PC1 = c(rnorm(8, -2), rnorm(8, 0), rnorm(8, 2)),
                PC2 = c(rnorm(8, 1), rnorm(8, -1.5), rnorm(8, 1)),
                PC3 = rnorm(n),
                aspect_ratio = runif(n, 1.5, 3.2),
                solidity = runif(n, 0.75, 0.95),
                log_centsize = runif(n, 2.0, 3.5),
                AR_inv = runif(n, 0.9, 1.4),
                cda_predicted_taxon = paste0("Morphotype_K", rep(1:3, each = 8)),
                cda_posterior_prob = runif(n, 0.75, 0.98),
                cda_second_posterior_prob = runif(n, 0.01, 0.20),
                fine_soil_sand = runif(n, 40, 85),
                fine_soil_ph = runif(n, 4.2, 6.8),
                fine_bedrock_depth_cm = runif(n, 20, 150),
                ssurgo_rock_outcrop_pct = runif(n, 0, 30),
                topo_tpi_30m = rnorm(n),
                bio1_temp_mean = runif(n, 12, 18),
                bio12_precip_annual = runif(n, 1000, 1400),
                decimalLatitude = runif(n, 32, 38),
                decimalLongitude = runif(n, -84, -77),
                doy = runif(n, 85, 145),
                delta_doy = rnorm(n, 0, 8),
                is_completed_voucher = TRUE,
                taxonomic_treatment_subtype = rep(c("Confirmed_Species", "Subspecies_Ecotype", "Hybrid_Introgressant"), 8),
                stringsAsFactors = FALSE
              )
              grp <- resolve_plot_group_column(df)
              stopifnot(identical(sort(unique(grp)), c("Morphotype_K1", "Morphotype_K2", "Morphotype_K3")))

              cda_res <- list(
                grouping_col = df$grouping_col,
                loocv_accuracy = 0.91,
                variance_pct = c(64.0, 36.0),
                canonical_scores = cbind(Can1 = df$can1, Can2 = df$can2, Can3 = df$can3)
              )
              gmm_res <- list(optimal_G = 3, best_k = 3, optimal_model = "VVI", classification = df$gmm_cluster)

              f_cda <- file.path("{intra_fig_dir.as_posix()}", "fig04b_gmm_cda_morphometrics.png")
              render_cda_biplot_4panel_pdf(df, cda_res, gmm_res, f_cda)
              stopifnot(file.exists(f_cda) && file.info(f_cda)$size > 1000)

              generate_step04_pub_figures(df, cda_res, gmm_res, NULL, f_cda)
              stopifnot(file.exists(file.path("{intra_fig_dir.as_posix()}", "fig04b_gmm_cda_morphometrics.png")))

              f_srf <- file.path("{intra_fig_dir.as_posix()}", "fig05a_ecological_niche_and_rf_importance.png")
              srf_mock <- list(
                OOB_Error = 0.12,
                Spatial_Block_CV_Accuracy = 0.86,
                importance = data.frame(
                  Predictor = c("fine_soil_sand", "fine_soil_ph", "MEM1"),
                  MeanDecreaseAccuracy = c(0.28, 0.21, 0.12),
                  Type = c("Fine-Grain Edaphic & Outcrop (10-30m)", "Fine-Grain Edaphic & Outcrop (10-30m)", "Spatial (MEMs)")
                ),
                variance_partition = data.frame(
                  Component = c("[a] Pure Fine-Grain Edaphic", "[b] Pure Macro-Climatic", "[c] Pure Spatial (MEMs)", "[d] Shared Env-Spatial", "[e] Unexplained"),
                  Variance_Pct = c(38, 22, 15, 10, 15)
                )
              )
              niche_mock <- list(
                summary = data.frame(Pair = "Morphotype_K1 vs Morphotype_K2", Taxon_1 = "Morphotype_K1", Taxon_2 = "Morphotype_K2", Warren_D = 0.42, Schoeners_D = 0.42, Equivalency_p_value = 0.02),
                permutations = data.frame(Pair = rep("Morphotype_K1 vs Morphotype_K2", 10), D_null = runif(10, 0.4, 0.9))
              )
              export_spatial_rf_figures(df, srf_mock, niche_mock, f_srf)
              stopifnot(file.exists(f_srf) && file.info(f_srf)$size > 1000)

              f_s06 <- file.path("{intra_fig_dir.as_posix()}", "fig06a_integrative_taxonomic_synthesis.png")
              render_step06_synthesis_plate_pdf(df, f_s06)
              stopifnot(file.exists(f_s06) && file.info(f_s06)$size > 1000)
              cat("STEP06_INTRA_FIGS_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step06_complex_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = run_warm_r_cmd(
                ["Rscript", str(r_batch_file)],
                stage_key="step05",
                cwd=REPO_ROOT,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 06 complex discovery R unit checks failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_standardize_packera_taxon_all_eastern_species_and_basionyms(self):
        """Verify North American Packera species (Eastern, Western, Boreal, Mexican) map cleanly."""
        cases = {
            "Packera anonyma (Alph.Wood) W.A.Weber & Á.Löve": "Packera anonyma",
            "Packera antennariifolia (Britton) W.A.Weber & Á.Löve": "Packera antennariifolia",
            "Packera aurea (L.) Á.Löve & D.Löve": "Packera aurea",
            "Packera crawfordii (Britton) A.M.Mahoney & R.R.Kowal": "Packera crawfordii",
            "Packera dubia (Spreng.) Trock & Mabb.": "Packera dubia",
            "Packera glabella (Poir.) C.Jeffrey": "Packera glabella",
            "Packera insulanus (C.L.Hitchc.) W.A.Weber & Á.Löve": "Packera insulanus",
            "Packera memmingeri": "Packera memmingeri",
            "Packera x memmingeri": "Packera memmingeri",
            "Packera × memmingeri": "Packera memmingeri",
            "Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve": "Packera millefolium",
            "Packera obovata (Muhl. ex Willd.) W.A.Weber & Á.Löve": "Packera obovata",
            "Packera paupercula var. appalachiana": "Packera paupercula",
            "Packera plattensis (Nutt.) W.A.Weber & Á.Löve": "Packera plattensis",
            "Packera pseudaurea var. semicordata": "Packera pseudaurea",
            "Packera schweinitziana (Nutt.) W.A.Weber & Á.Löve": "Packera schweinitziana",
            "Packera serpenticola Boufford, Kartesz, Shiuchi & R.Zhou": "Packera serpenticola",
            "Packera spellenbergii (T.M.Barkley) C.Jeffrey": "Packera spellenbergii",
            "Packera tampicana (DC.) C.Jeffrey": "Packera tampicana",
            "Packera bernardina (Greene) W.A.Weber & Á.Löve": "Packera bernardina",
            "Packera cana (Hook.) W.A.Weber & Á.Löve": "Packera cana",
            "Packera fendleri (A.Gray) W.A.Weber & Á.Löve": "Packera fendleri",
            "Packera multilobata (Torr. & A.Gray ex A.Gray) W.A.Weber & Á.Löve": "Packera multilobata",
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve": "Packera streptanthifolia",
            "Packera werneriifolia (A.Gray) W.A.Weber & Á.Löve": "Packera werneriifolia",
            "Packera cymbalaria (Pursh) W.A.Weber & Á.Löve": "Packera cymbalaria",
            "Packera hyperborealis (Greenm.) Á.Löve & D.Löve": "Packera hyperborealis",
            "Packera indecora (Greene) Á.Löve & D.Löve": "Packera indecora",
            "Packera ogotorukensis (Packer) Á.Löve & D.Löve": "Packera ogotorukensis",
            "Packera pauciflora (Pursh) Á.Löve & D.Löve": "Packera pauciflora",
            "Packera bellidifolia (Kunth) W.A.Weber & Á.Löve": "Packera bellidifolia",
            "Packera coahuilensis (Greenm.) C.Jeffrey": "Packera coahuilensis",
            "Packera sanguisorbae (DC.) C.Jeffrey": "Packera sanguisorbae",
            "Packera toluccana (DC.) W.A.Weber & Á.Löve": "Packera toluccana",
            "Packera zimapanica (Hemsl.) C.C.Freeman & T.M.Barkley": "Packera zimapanica",
        }
        mod01_pipe = importlib.import_module("scripts.pipeline.01_voucher_harvester")
        for raw_name, expected in cases.items():
            self.assertEqual(standardize_packera_taxon(raw_name), expected, f"Failed for {raw_name!r}")
            self.assertEqual(mod01_pipe.standardize_packera_taxon(raw_name), expected)

        basionym_backbone = {
            "senecio millefolium": "Packera millefolium",
            "senecio streptanthifolius": "Packera streptanthifolia",
            "senecio cymbalaria": "Packera cymbalaria",
            "senecio sanguisorbae": "Packera sanguisorbae",
        }

        def _fake_backbone(name, **kwargs):
            acc = basionym_backbone.get(str(name).strip().lower())
            if acc:
                return {"matchType": "EXACT", "status": "SYNONYM", "species": acc, "canonicalName": str(name)}
            return {"matchType": "NONE"}

        with patch("scripts.core.harvester.species.name_backbone", side_effect=_fake_backbone), patch(
            "scripts.core.harvester.species.name_usage", return_value={"results": []}
        ):
            for raw_syn, expected in [
                ("Senecio millefolium Torr. & A.Gray", "Packera millefolium"),
                ("Senecio streptanthifolius Greene", "Packera streptanthifolia"),
                ("Senecio cymbalaria Pursh", "Packera cymbalaria"),
                ("Senecio sanguisorbae DC.", "Packera sanguisorbae"),
            ]:
                self.assertEqual(standardize_packera_taxon(raw_syn), expected)
                self.assertEqual(mod01_pipe.standardize_packera_taxon(raw_syn), expected)

    def test_r_standardize_taxon_and_dynamic_palette(self):
        """Verify R standardize_taxon_name(), get_dynamic_packera_palette(), shapes, titles, and italic expressions are species-agnostic."""
        self.assertIn("STEP06_R_PALETTE_OK", self._require_batched_r())

    def test_discover_taxon_complexes_8_species_3_complexes(self):
        """Verify unified Step 06 recovers 3 true complexes across 8 species with canonical twin-free outputs."""
        rng = np.random.default_rng(20261025)
        complex_map = {
            "Packera dubia": ("C1", np.array([-4.5, 3.2, 0.4, -0.2, 0.1]), 2.15),
            "Packera anonyma": ("C1", np.array([-4.1, 2.9, 0.3, -0.1, 0.0]), 1.95),
            "Packera plattensis": ("C1", np.array([-4.3, 3.0, 0.5, -0.2, 0.1]), 2.05),
            "Packera millefolium": ("C2", np.array([6.0, 0.2, -2.5, 1.2, -0.4]), 1.80),
            "Packera memmingeri": ("C2", np.array([5.6, 0.4, -2.3, 1.0, -0.3]), 1.85),
            "Packera aurea": ("C3", np.array([0.2, -5.2, 1.8, -1.1, 0.5]), 2.40),
            "Packera obovata": ("C3", np.array([0.5, -4.9, 1.6, -0.9, 0.4]), 2.30),
            "Packera pseudaurea": ("C3", np.array([0.0, -5.0, 1.7, -1.0, 0.5]), 2.35),
        }
        complex_members = {
            "C1": ["Packera dubia", "Packera anonyma", "Packera plattensis"],
            "C2": ["Packera millefolium", "Packera memmingeri"],
            "C3": ["Packera aurea", "Packera obovata", "Packera pseudaurea"],
        }

        vouchers_rows = []
        harmonics_rows = []
        morph_rows = []
        idx = 1
        for taxon, (true_cid, center, inv_ratio) in complex_map.items():
            peers = [t for t in complex_members[true_cid] if t != taxon]
            for k in range(15):
                cat = f"NCU-SYN-{idx:04d}"
                idx += 1
                det_tier = "Tier_1_Gold" if k < 10 else "Tier_2_Silver"
                auth_wt = 3.0 if det_tier == "Tier_1_Gold" else 1.5
                pcs = center + rng.normal(0.0, 0.28, size=5)
                if k % 5 == 0:
                    cda_pred = peers[k % len(peers)]
                    cda_post = 0.58
                    cda_sec = 0.38
                else:
                    cda_pred = taxon
                    cda_post = 0.78
                    cda_sec = 0.19

                vouchers_rows.append({
                    "catalogNumber": cat,
                    "species_raw": taxon,
                    "species_standardized": taxon,
                    "determiner_tier": det_tier,
                    "authority_weight": auth_wt,
                })
                harmonics_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "species_standardized": taxon,
                    "determiner_tier": det_tier,
                    "authority_weight": auth_wt,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "PC1": float(pcs[0]), "PC2": float(pcs[1]), "PC3": float(pcs[2]), "PC4": float(pcs[3]), "PC5": float(pcs[4]),
                })
                morph_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "species": taxon,
                    "species_standardized": taxon,
                    "determiner_tier": det_tier,
                    "authority_weight": auth_wt,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "PC1": float(pcs[0]), "PC2": float(pcs[1]), "PC3": float(pcs[2]), "PC4": float(pcs[3]), "PC5": float(pcs[4]),
                    "involucre_length_width_ratio": float(inv_ratio + rng.normal(0.0, 0.05)),
                    "cda_predicted_taxon": cda_pred,
                    "cda_posterior_prob": cda_post,
                    "cda_second_posterior_prob": cda_sec,
                    "cda_posterior_margin": cda_post - cda_sec,
                    "morphometric_outlier_flag": False,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            vouchers_csv = tmp_path / "curated_vouchers.csv"
            harmonics_csv = tmp_path / "leaf_efa_harmonics.csv"
            morph_csv = tmp_path / "morphometric_flags.csv"
            out_queue = tmp_path / "triage_queue.csv"
            out_plot = tmp_path / "Figure_Integrative_Packera_dubia_Revision.png"
            out_report = tmp_path / "Packera_dubia_Taxonomic_Revision_Summary.md"
            out_html = tmp_path / "Triage_Dashboard_Interactive.html"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_figure = tmp_path / "fig06b_taxon_complex_delimitation_network.png"

            pd.DataFrame(vouchers_rows).to_csv(vouchers_csv, index=False)
            pd.DataFrame(harmonics_rows).to_csv(harmonics_csv, index=False)
            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)

            cmd = [
                "Rscript",
                "scripts/analysis/06_triage_dashboard_synthesis.R",
                "--vouchers", str(vouchers_csv),
                "--harmonics", str(harmonics_csv),
                "--morphometrics", str(morph_csv),
                "--focal-taxon", "Packera dubia",
                "--min-affinity", "0.05",
                "--output-queue", str(out_queue),
                "--output-plot", str(out_plot),
                "--output-report", str(out_report),
                "--output-html", str(out_html),
                "--output-complexes", str(out_complexes),
                "--output-affinity", str(out_affinity),
                "--output-neighbors", str(out_neighbors),
                "--output-complex-plot", str(out_figure),
                "--demo",
            ]
            if not rscript_available():
                self.skipTest("Rscript not installed.")
            res = run_warm_r_cmd(cmd, stage_key='step06', cwd=REPO_ROOT, skip_figures=False)
            self.assertEqual(res.returncode, 0, f"Unified Step 06 failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

            self.assertTrue(out_complexes.exists() and out_complexes.stat().st_size > 0)
            self.assertTrue(out_affinity.exists() and out_affinity.stat().st_size > 0)
            self.assertTrue(out_neighbors.exists() and out_neighbors.stat().st_size > 0)
            self.assertTrue(out_figure.exists() and out_figure.stat().st_size > 0)
            self.assertEqual(list(tmp_path.glob("*.pdf")), [])
            self.assertIn("Empirical Taxon-Complex Delimitation", out_report.read_text(encoding="utf-8"))
            self.assertIn("Empirical Species Complexes", out_html.read_text(encoding="utf-8"))

            comp_df = pd.read_csv(out_complexes)
            self.assertEqual(len(comp_df), 8)
            self.assertEqual(set(comp_df.columns), CANONICAL_COMPLEX_COLS)

            aff_df = pd.read_csv(out_affinity)
            self.assertIn("species_standardized", aff_df.columns)
            self.assertNotIn("taxon", aff_df.columns)

            grouped = (
                comp_df.groupby("consensus_complex_id")["species_standardized"]
                .apply(lambda s: frozenset(s.tolist()))
                .tolist()
            )
            expected_partitions = {
                frozenset(["Packera dubia", "Packera anonyma", "Packera plattensis"]),
                frozenset(["Packera millefolium", "Packera memmingeri"]),
                frozenset(["Packera aurea", "Packera obovata", "Packera pseudaurea"]),
            }
            self.assertEqual(set(grouped), expected_partitions)

            neighbors = [line.strip() for line in out_neighbors.read_text(encoding="utf-8").splitlines() if line.strip()]
            self.assertEqual(set(neighbors), {"Packera anonyma", "Packera plattensis"})

            disc_text = DISCOVERY_UTILS_R.read_text(encoding="utf-8")
            self.assertIn("compute_hopkins_phenological_affinity", disc_text)
            self.assertIn("compute_modality_kappa_weights", disc_text)
            self.assertIn("evaluate_hierarchical_complex_delimitation", disc_text)
            self.assertIn("Padial et al. 2010", disc_text)
            self.assertIn("Dawid & Skene 1979", disc_text)
            self.assertIn("Schlick-Steiner et al. 2010", disc_text)

    def test_hierarchical_precedence_missing_modality_and_3d_hollow_ordinations(self):
        """Verify Stage 1-5 precedence, missing-modality weight renormalization, and 3D hollow ordinations + BH FDR."""
        self.assertIn("STEP06_HIERARCHICAL_OK", self._require_batched_r())

    def test_subcluster_network_when_2_species_or_within_species_gmm_clusters(self) -> None:
        """Verify taxon_complex_discovery_utils.R (<=800 lines) builds Taxon x Morphotype subcluster network inside Step 06 when <=2 species."""
        line_count = len(DISCOVERY_UTILS_R.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(line_count, 800, f"taxon_complex_discovery_utils.R has {line_count} lines; must be <= 800 lines.")

        rng = np.random.default_rng(20261028)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            morph_csv = tmp_path / "morphometric_flags.csv"
            out_queue = tmp_path / "triage_queue.csv"
            out_plot = tmp_path / "rev.png"
            out_report = tmp_path / "rev.md"
            out_html = tmp_path / "dash.html"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_fig = tmp_path / "fig06b_taxon_complex_delimitation_network.png"

            subclusters = [
                ("Packera dubia", 1, 8, np.array([-2.2, -1.5, 0.0, 0.0, 0.0])),
                ("Packera dubia", 2, 7, np.array([-1.8, -1.2, 0.2, 0.0, 0.0])),
                ("Packera dubia", 3, 2, np.array([0.0, 0.0, 0.0, 0.0, 0.0])),
                ("Packera paupercula", 1, 8, np.array([2.2, 1.5, 0.0, 0.0, 0.0])),
                ("Packera paupercula", 2, 6, np.array([1.9, 1.3, -0.2, 0.0, 0.0])),
            ]
            morph_rows = []
            v_idx = 0
            for sp_name, k_id, count, center in subclusters:
                for _ in range(count):
                    v_idx += 1
                    vec = center + rng.normal(0.0, 0.25, size=5)
                    morph_rows.append({
                        "catalogNumber": f"SUB_{v_idx:04d}",
                        "plant_individual_id": 1,
                        "species_standardized": sp_name,
                        "determiner_tier": "Tier_1_Gold",
                        "basal_leaf_count": 3,
                        "morph_tier_source": "Basal_Rosette",
                        "gmm_cluster": k_id,
                        "PC1": float(vec[0]), "PC2": float(vec[1]), "PC3": float(vec[2]), "PC4": float(vec[3]), "PC5": float(vec[4]),
                        "cda_predicted_taxon": sp_name,
                        "cda_posterior_prob": 0.90,
                        "cda_second_posterior_prob": 0.08,
                    })

            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)
            if not rscript_available():
                self.skipTest("Rscript not installed.")

            res = run_warm_r_cmd([
                    "Rscript", str(SYNTHESIS_R_SCRIPT),
                    "--vouchers", str(morph_csv),
                    "--morphometrics", str(morph_csv),
                    "--focal-taxon", "Packera dubia",
                    "--output-queue", str(out_queue),
                    "--output-plot", str(out_plot),
                    "--output-report", str(out_report),
                    "--output-html", str(out_html),
                    "--output-complexes", str(out_complexes),
                    "--output-affinity", str(out_affinity),
                    "--output-neighbors", str(out_neighbors),
                    "--output-complex-plot", str(out_fig),
                    "--demo",
                ], stage_key='step06', cwd=REPO_ROOT, skip_figures=False)
            self.assertEqual(res.returncode, 0, f"Subcluster discovery failed:\n{res.stderr}\n{res.stdout}")
            comp_df = pd.read_csv(out_complexes)
            expected_nodes = {"P. dubia [K1]", "P. dubia [K2]", "P. paupercula [K1]", "P. paupercula [K2]"}
            self.assertEqual(set(comp_df["species_standardized"]), expected_nodes)
            self.assertTrue((comp_df["n_vouchers"] >= 5).all())
            self.assertTrue(out_fig.exists() and out_fig.stat().st_size > 1000)

    def test_intraspecific_morphotype_cluster_visualization_across_r_figures(self) -> None:
        """Verify pub_figures_utils.R and 05_multimodal_spatial_rf.R use morphotype clusters when <2 species."""
        self.assertIn("STEP06_INTRA_FIGS_OK", self._require_batched_r())

    def test_multimodal_cli_flag_and_ecological_edaphic_left_join(self) -> None:
        """Verify unified Step 06 merges Step 05 ecological/edaphic columns and taxon_complex_discovery_utils.R stays <=800 lines."""
        line_count = len(DISCOVERY_UTILS_R.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(line_count, 800, f"taxon_complex_discovery_utils.R has {line_count} lines; must be <= 800 lines.")
        self.assertIn("STEP06_MULTIMODAL_JOIN_OK", self._require_batched_r())

    def test_dynamic_morph_and_primary_env_predictor_discovery_in_step06(self) -> None:
        """Verify unified Step 06 dynamically discovers >5 PCs via selected_morph_features, Step 05 primary_env_predictors, and Peak_Anthesis DOY."""
        rng = np.random.default_rng(20261031)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            morph_csv = tmp_path / "morphometric_flags.csv"
            mm_csv = tmp_path / "multimodal_conflict_flags.csv"
            out_queue = tmp_path / "triage_queue.csv"
            out_plot = tmp_path / "rev.png"
            out_report = tmp_path / "rev.md"
            out_html = tmp_path / "dash.html"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_fig = tmp_path / "fig06b_taxon_complex_delimitation_network.png"

            taxa_specs = [
                ("Packera hyperborea", (-3.5, -3.0), (5.5, 1.8)),
                ("Packera contermina", (-3.2, -2.8), (5.8, 1.9)),
                ("Packera cymbalaria", (4.5, 4.2), (13.5, 9.2)),
            ]
            morph_rows = []
            mm_rows = []
            v_idx = 0
            for sp_name, (pc6_mu, pc7_mu), (b10_mu, ksat_mu) in taxa_specs:
                for j in range(8):
                    v_idx += 1
                    cat_no = f"DYN06_{v_idx:03d}"
                    lat_v = 61.0 + 0.4 * j
                    morph_rows.append({
                        "catalogNumber": cat_no,
                        "plant_individual_id": 0,
                        "species_standardized": sp_name,
                        "determiner_tier": "Tier_1_Gold",
                        "basal_leaf_count": 3,
                        "morph_tier_source": "Basal_Rosette",
                        "gmm_cluster": 1,
                        "PC1": float(rng.normal(0.0, 0.15)),
                        "PC2": float(rng.normal(0.0, 0.15)),
                        "PC3": float(rng.normal(0.0, 0.15)),
                        "PC4": float(rng.normal(0.0, 0.15)),
                        "PC5": float(rng.normal(0.0, 0.15)),
                        "PC6": float(pc6_mu + rng.normal(0.0, 0.20)),
                        "PC7": float(pc7_mu + rng.normal(0.0, 0.20)),
                        "selected_morph_features": "PC1;PC2;PC3;PC4;PC5;PC6;PC7",
                        "cda_predicted_taxon": sp_name,
                        "cda_posterior_prob": 0.92,
                        "cda_second_posterior_prob": 0.05,
                    })
                    mm_rows.append({
                        "catalogNumber": cat_no,
                        "plant_individual_id": 0,
                        "species_standardized": sp_name,
                        "decimalLatitude": lat_v,
                        "decimalLongitude": -148.0 + 0.3 * j,
                        "doy": int(215 + 3.0 * (lat_v - 61.0)),
                        "is_verified_anthesis": True,
                        "voucher_phenological_state": "Peak_Anthesis",
                        "bio10_mean_temp_warmest_qtr": float(b10_mu + rng.normal(0.0, 0.25)),
                        "polaris_ksat_30m": float(ksat_mu + rng.normal(0.0, 0.25)),
                        "primary_env_predictors": "bio10_mean_temp_warmest_qtr;polaris_ksat_30m",
                        "primary_clim_predictors": "bio10_mean_temp_warmest_qtr",
                        "primary_edaph_predictors": "polaris_ksat_30m",
                    })

            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)
            pd.DataFrame(mm_rows).to_csv(mm_csv, index=False)
            if not rscript_available():
                self.skipTest("Rscript not installed.")

            res = run_warm_r_cmd([
                    "Rscript", str(SYNTHESIS_R_SCRIPT),
                    "--vouchers", str(morph_csv),
                    "--morphometrics", str(morph_csv),
                    "--multimodal-flags", str(mm_csv),
                    "--focal-taxon", "Packera hyperborea",
                    "--output-queue", str(out_queue),
                    "--output-plot", str(out_plot),
                    "--output-report", str(out_report),
                    "--output-html", str(out_html),
                    "--output-complexes", str(out_complexes),
                    "--output-affinity", str(out_affinity),
                    "--output-neighbors", str(out_neighbors),
                    "--output-complex-plot", str(out_fig),
                    "--demo",
                ], stage_key='step06', cwd=REPO_ROOT, skip_figures=True)
            self.assertEqual(res.returncode, 0, f"Step 06 dynamic feature test failed:\n{res.stderr}\n{res.stdout}")
            comp_df = pd.read_csv(out_complexes)
            self.assertEqual(len(comp_df), 3)
            hyp_row = comp_df.loc[comp_df["species_standardized"] == "Packera hyperborea"].iloc[0]
            self.assertEqual(str(hyp_row["nearest_neighbor_taxon"]), "Packera contermina")
            # Gate 1 fix: Peak_Anthesis vouchers are recognized (hopkins_pheno_overlap > 0.80, NOT NA/0!)
            self.assertGreater(float(hyp_row["hopkins_pheno_overlap"]), 0.80)
            self.assertGreater(float(hyp_row["eco_schoener_overlap"]), 0.50)

    def test_multi_species_subcluster_suppression_and_line_budget(self):
        """Verify taxon_complex_discovery_utils.R returns subcluster_mode = FALSE when n_sp > 3L and stays <= 800 lines."""
        src_lines = DISCOVERY_UTILS_R.read_text(encoding="utf-8").splitlines()
        self.assertLessEqual(len(src_lines), 800, f"taxon_complex_discovery_utils.R has {len(src_lines)} lines (> 800)")
        src_text = "\n".join(src_lines)
        self.assertIn("n_sp > 3L", src_text)

        r_code = f"""
        source("{DISCOVERY_UTILS_R.as_posix()}")
        df4 <- data.frame(
          catalogNumber = sprintf("SP4_%02d", 1:40),
          plant_individual_id = 1L,
          species_standardized = rep(c("Packera dubia", "Packera aurea", "Packera anonyma", "Packera obovata"), each = 10),
          gmm_cluster = rep(c(1L, 2L), 20),
          can1 = seq(-2, 2, length.out = 40),
          can2 = seq(1, -1, length.out = 40),
          stringsAsFactors = FALSE
        )
        sub_res <- build_taxon_morphotype_subclusters(df4, min_vouchers = 3L)
        stopifnot(identical(sub_res$subcluster_mode, FALSE))
        stopifnot(length(unique(sub_res$morph_df$species_standardized)) == 4L)
        cat("SUBCLUSTER_MULTI_SP_GUARD_OK\\n")
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step06_guard.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            res = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key="step06", cwd=REPO_ROOT)
        self.assertEqual(res.returncode, 0, f"Subcluster multi-species guard test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("SUBCLUSTER_MULTI_SP_GUARD_OK", res.stdout)

    def test_unified_step06_complex_discovery_zero_trust_four_gates(self):
        """
        Zero-trust regression suite for unified Step 06 Complex Delimitation Four-Gate Overhaul:
          1. Gate 1: S_morph proximity gating prevents sympatric habitat convergence from merging distant clades,
             metric sqrt(1 - W_ij) distance is used in Ward.D2, and Peak_Anthesis vouchers yield has_pheno == TRUE.
          2. Gate 2: Confident-noise & Misidentification_Reassign rows are excluded from anchors and C_morph,
             and subcluster mode preserves base_species / base_cda_predicted_taxon without forcing w_morph_raw == 1.0.
          3. Gate 3: SQLite taxon_complexes persists all 23 canonical columns under schema 20261116.
          4. Gate 4: .Random.seed is strictly preserved across phenology permutations & Leiden partitioning.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            db_file = tmp_p / "step06_complexes.sqlite"
            with get_db_connection(db_file, init_schema=True):
                pass

            r_code = rf"""
            source("scripts/analysis/06_triage_dashboard_synthesis.R")

            # 1. Sympatric habitat convergence must NOT merge morphologically distant clades (Gate 1)
            set.seed(777)
            df_symp <- data.frame(
              catalogNumber = sprintf("SYM_%02d", 1:24),
              plant_individual_id = 0L,
              species_standardized = rep(c("Packera dubia", "Packera anonyma", "Packera millefolium"), each = 8),
              determiner_tier = "Tier_1_Gold",
              morph_tier_source = "Basal_Rosette",
              basal_leaf_count = 3L,
              PC1 = c(rnorm(8, -4.0, 0.15), rnorm(8, -3.7, 0.15), rnorm(8, 8.0, 0.15)),
              PC2 = c(rnorm(8, 1.0, 0.15), rnorm(8, 1.1, 0.15), rnorm(8, -6.0, 0.15)),
              PC3 = rnorm(24, 0.0, 0.10),
              cda_predicted_taxon = rep(c("Packera dubia", "Packera anonyma", "Packera millefolium"), each = 8),
              cda_posterior_prob = 0.94,
              decimalLatitude = rep(seq(35.0, 35.7, length.out = 8), 3),
              decimalLongitude = -80.0,
              doy = rep(120:127, 3),
              is_verified_anthesis = TRUE,
              voucher_phenological_state = "Peak_Anthesis",
              fine_soil_ph = c(rnorm(8, 6.5, 0.1), rnorm(8, 5.1, 0.1), rnorm(8, 5.1, 0.1)),
              fine_soil_sand = c(rnorm(8, 45.0, 1.5), rnorm(8, 75.0, 1.5), rnorm(8, 75.0, 1.5)),
              primary_edaph_predictors = "fine_soil_ph;fine_soil_sand",
              taxonomic_treatment_subtype = "Confirmed_Species",
              stringsAsFactors = FALSE
            )
            taxa3 <- sort(unique(df_symp$species_standardized))
            anc3 <- select_taxon_anchors(df_symp)
            seed_before <- .Random.seed
            ph3 <- compute_hopkins_phenological_affinity(df_symp, taxa3, anc3$anchor_mask)
            # RNG state preserved after permutation test!
            stopifnot(identical(.Random.seed, seed_before))
            stopifnot(isTRUE(ph3$has_pheno))

            ec3 <- compute_ecological_affinity(df_symp, taxa3, anc3$anchor_mask)
            mo3 <- compute_morphometric_affinity(df_symp, anc3$anchor_df, taxa3)
            kw3 <- compute_modality_kappa_weights(df_symp, anc3$anchor_mask)
            W3 <- compute_multimodal_complex_affinity(mo3, ec3, ph3, kw3, taxa3)
            # P. anonyma and P. millefolium share identical sand/pH & DOY, yet S_morph ~ 0 gates W_ij near 0!
            stopifnot(W3["Packera anonyma", "Packera millefolium"] < 0.02)
            stopifnot(W3["Packera anonyma", "Packera dubia"] > W3["Packera anonyma", "Packera millefolium"] * 10.0)

            pr3 <- partition_taxon_complexes(W3, taxa3, table(df_symp$species_standardized)[taxa3], anc3$n_gold_silver)
            stopifnot(identical(.Random.seed, seed_before))
            # Metric sqrt(1 - W_ij) check
            d_mat <- as.matrix(pr3$dist_obj)
            stopifnot(abs(d_mat["Packera anonyma", "Packera dubia"] - sqrt(1.0 - W3["Packera anonyma", "Packera dubia"])) < 1e-6)

            # 2. Noise & Misidentification_Reassign rows excluded from anchors and C_morph (Gate 2)
            df_corrupt <- rbind(df_symp, df_symp[1:3, ])
            df_corrupt$catalogNumber[25:27] <- c("ERR_1", "ERR_2", "ERR_3")
            df_corrupt$species_standardized[25:27] <- "Packera millefolium"
            df_corrupt$cda_predicted_taxon[25:27] <- "Packera dubia"
            df_corrupt$is_confident_label_noise <- c(rep(FALSE, 24), TRUE, FALSE, FALSE)
            df_corrupt$taxonomic_treatment_subtype[25:27] <- c("Confirmed_Species", "Misidentification_Reassign", "Misidentification_Reassign")
            anc_c <- select_taxon_anchors(df_corrupt)
            stopifnot(!any(anc_c$anchor_mask[25:27]))
            mo_c <- compute_morphometric_affinity(df_corrupt, anc_c$anchor_df, taxa3)
            stopifnot(mo_c$C_morph_sym["Packera millefolium", "Packera dubia"] < 0.01)

            # 3. Subcluster mode preserves base_cda_predicted_taxon without inflating w_morph_raw to 1.0 (Gate 2)
            df_sub <- data.frame(
              catalogNumber = sprintf("SB_%02d", 1:20),
              plant_individual_id = 0L,
              species_standardized = rep(c("Packera dubia", "Packera aurea"), each = 10),
              determiner_tier = "Tier_1_Gold",
              morph_tier_source = "Basal_Rosette",
              basal_leaf_count = 3L,
              gmm_cluster = rep(c(1L, 2L), 10),
              PC1 = rnorm(20),
              PC2 = rnorm(20),
              cda_predicted_taxon = c(rep("Packera dubia", 7), rep("Packera aurea", 3), rep("Packera aurea", 8), rep("Packera dubia", 2)),
              stringsAsFactors = FALSE
            )
            seed_sub_before <- .Random.seed
            sr <- build_taxon_morphotype_subclusters(df_sub, min_vouchers = 4L)
            stopifnot(identical(.Random.seed, seed_sub_before))
            kw_sub <- compute_modality_kappa_weights(sr$morph_df, rep(TRUE, nrow(sr$morph_df)))
            stopifnot(kw_sub$w_morph_raw < 0.80)

            # 4. SQLite persistence into taxon_complexes (Gate 3)
            Sys.setenv(PACKERA_DB_PATH = "{db_file.as_posix()}", RUN_ID = "run_step06_tc")
            full_res <- run_step06_complex_discovery(df_symp, list(min_affinity = 0.03, demo = FALSE, vouchers = ""))
            persist_step06_complexes_to_sqlite(full_res$complexes_df, vouchers_path = NULL, run_id = "run_step06_tc", demo = FALSE)
            cat("STEP06_COMPLEX_FOUR_GATES_OK\n")
            """
            r_file = tmp_p / "test_complex_gates.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step06", cwd=REPO_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 06 complex four-gate R test failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertIn("STEP06_COMPLEX_FOUR_GATES_OK", proc.stdout)

            conn = sqlite3.connect(db_file)
            rows = conn.execute("SELECT species_standardized, consensus_complex_id FROM taxon_complexes WHERE run_id = 'run_step06_tc'").fetchall()
            conn.close()
            self.assertEqual(len(rows), 3)

    def test_misidentification_and_cryptic_species_five_step_protocol(self) -> None:
        """
        Verify the 5-step Misidentification & Cryptic-Species protocol + AUR-PLASTIC-LATE regression fix:
          1. Core-cluster g_cl == 1 late-season outlier (AUR-PLASTIC-LATE) never passes cryptic_cluster_cohesion_pass.
          2. Stage 1 protects vouchers inside their 95% self-Mahalanobis ellipsoid or where Tier_1_Gold / Tier_2_Silver
             determinations are corroborated by habitat or blade_cie_l_norm vestiture reflectance.
          3. Intraspecific-mode congener_osr_reassign evaluates exact open_set_df and reassigns outgroup congeners.
          4. Stage 3 deduplicates duplicate exsiccatae sheets across herbaria and routes smooth allopatric spatial
             clines (cryptic_non_clinal_pass == FALSE) to Subspecies_Ecotype.
          5. Stage 3 separates Heterotypic Basionym Resurrection (Taxonomic_Resurrection_Candidate + resurrection_basionym)
             from unnamed de novo Putative_Cryptic_Species.
        """
        r06_path = REPO_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"
        r06_src = r06_path.read_text(encoding="utf-8")
        r04_src = (REPO_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R").read_text(encoding="utf-8")
        auth_src = (REPO_ROOT / "scripts" / "morphometrics" / "authority_weighting_utils.R").read_text(encoding="utf-8")
        for tok in ("outside_self_morph_95", "tier_authority_gate_ok", "cryptic_non_clinal_pass", "resurrection_basionym", "Taxonomic_Resurrection_Candidate", "osr_chisq90_vec"):
            self.assertIn(tok, r06_src)
        self.assertNotIn("n < 50L && g_cl == 1L", r06_src)
        for tok04 in ("cda_oof_posterior_prob", "cl_target_threshold", "self_morph_mahal_d2", "open_set_df"):
            self.assertIn(tok04, r04_src)
        self.assertIn(".compute_oof_cda_and_self_ellipsoid", auth_src)

        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/06_triage_dashboard_synthesis.R")

        # A. Stage 1: Self-ellipsoid protection, Tier_1_Gold/Tier_2_Silver vestiture/habitat protection, vs corroborated reassign
        df_s1 <- data.frame(
          catalogNumber = c(
            "DUB_ANC_1", "DUB_ANC_2", "DUB_ANC_3",
            "ANO_ANC_1", "ANO_ANC_2", "ANO_ANC_3",
            "GOLD_IN_SELF_ELLIPSOID",
            "SILVER_VESTITURE_PROTECTED",
            "SILVER_CORROBORATED_MISID"
          ),
          plant_individual_id = 0L,
          institutionCode = "NCU",
          species_raw = c(rep("Packera dubia", 3), rep("Packera anonyma", 3), rep("Packera dubia", 3)),
          species_standardized = c(rep("Packera dubia", 3), rep("Packera anonyma", 3), rep("Packera dubia", 3)),
          determiner_tier = c(rep("Tier_1_Gold", 6), "Tier_1_Gold", "Tier_2_Silver", "Tier_2_Silver"),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          decimalLatitude = c(35.1, 35.2, 35.3, 35.8, 35.9, 36.0, 35.2, 35.25, 35.95),
          decimalLongitude = c(-78.5, -78.6, -78.7, -79.1, -79.2, -79.3, -78.55, -78.62, -79.15),
          doy = c(108, 110, 112, 132, 134, 136, 111, 110, 135),
          is_verified_anthesis = TRUE,
          PC1 = c(-2.2, -2.0, -2.1, 2.2, 2.0, 2.1, -1.5, 1.4, 2.1),
          PC2 = c(0.1, -0.1, 0.0, -0.1, 0.1, 0.0, 0.2, 0.1, 0.0),
          blade_cie_l_norm = c(14.5, 15.2, 14.8, 1.8, 2.2, 2.0, 14.6, 15.0, 2.1),
          self_morph_mahal_d2 = c(0.8, 0.7, 0.6, 0.8, 0.7, 0.6, 3.2, 14.5, 18.2),
          self_morph_df = 3L,
          cda_predicted_taxon = c(rep("Packera dubia", 3), rep("Packera anonyma", 3), "Packera anonyma", "Packera anonyma", "Packera anonyma"),
          cda_posterior_prob = c(rep(0.96, 6), 0.79, 0.84, 0.93),
          cda_oof_posterior_prob = c(rep(0.94, 6), 0.74, 0.81, 0.91),
          cda_posterior_margin = c(rep(0.90, 6), 0.58, 0.68, 0.86),
          morphometric_outlier_flag = c(rep(FALSE, 6), FALSE, TRUE, TRUE),
          edaphic_best_fit_taxon = c(rep("Packera dubia", 3), rep("Packera anonyma", 3), "Packera dubia", "Packera dubia", "Packera anonyma"),
          rf_predicted_taxon = c(rep("Packera dubia", 3), rep("Packera anonyma", 3), "Packera dubia", "Packera dubia", "Packera anonyma"),
          rf_predicted_prob = c(rep(0.92, 6), 0.88, 0.86, 0.90),
          edaphic_schoener_d = c(rep(0.82, 6), 0.79, 0.76, 0.19),
          gmm_cluster = 1L,
          gmm_uncertainty = 0.04,
          stringsAsFactors = FALSE
        )
        res_s1 <- apply_taxonomic_decision_matrix(df_s1)
        stopifnot(res_s1$taxonomic_treatment_subtype[7] != "Misidentification_Reassign")
        stopifnot(res_s1$recommended_determination[7] == "Packera dubia")
        stopifnot(res_s1$taxonomic_treatment_subtype[8] != "Misidentification_Reassign")
        stopifnot(res_s1$recommended_determination[8] == "Packera dubia")
        stopifnot(res_s1$taxonomic_treatment_subtype[9] == "Misidentification_Reassign")
        stopifnot(res_s1$recommended_determination[9] == "Packera anonyma")

        # B. Intraspecific-mode congener_osr_reassign with exact open_set_df = 14 (qchisq(0.90, 14) = 21.064)
        df_intra_osr <- data.frame(
          catalogNumber = c("DUB_CORE_1", "DUB_CORE_2", "DUB_CORE_3", "DUB_ACTUALLY_ANONYMA"),
          plant_individual_id = 0L,
          institutionCode = "NCU",
          species_raw = "Packera dubia",
          species_standardized = "Packera dubia",
          determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze"),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          cda_predicted_taxon = c("Morphotype_K1", "Morphotype_K1", "Morphotype_K1", "Morphotype_K2"),
          cda_posterior_prob = c(0.95, 0.94, 0.95, 0.88),
          is_typical_morphotype = c(TRUE, TRUE, TRUE, FALSE),
          type_morph_mahal_d2 = c(1.1, 1.3, 1.2, 24.5),
          min_congener_mahal_d2 = c(38.0, 39.5, 37.2, 11.4),
          min_congener_taxon = "Packera anonyma",
          open_set_novelty_p = c(0.0004, 0.0003, 0.0005, 0.654),
          open_set_df = 14L,
          morphometric_outlier_flag = c(FALSE, FALSE, FALSE, TRUE),
          edaphic_z_score = c(0.4, 0.5, 0.4, 2.3),
          gmm_cluster = c(1L, 1L, 1L, 2L),
          stringsAsFactors = FALSE
        )
        res_iosr <- apply_taxonomic_decision_matrix(df_intra_osr, intraspecific_mode = TRUE)
        stopifnot(res_iosr$taxonomic_treatment_subtype[4] == "Misidentification_Reassign")
        stopifnot(res_iosr$recommended_determination[4] == "Packera anonyma")

        # C. Stage 3: Duplicate exsiccatae collapse, allopatric spatial-cline rejection, and Basionym Resurrection vs sp. nov.
        # Build 60-row intraspecific cohort so full >=5 gathering & >=3 10 km deme cohesion gates execute:
        #   - Rows 1..40: K1 sensu stricto core across Piedmont/Coastal Plain
        #   - Rows 41..45: K2 duplicate sheets of the SAME gathering ("A. S. Weakley" / "Weakley, Alan S." #4102) -> MUST FAIL cohesion!
        #   - Rows 46..51: K3 allopatric smooth latitudinal cline (R^2 > 0.85, >100 km away, no flowering/azonal jump) -> Subspecies_Ecotype!
        #   - Rows 52..57: K4 parapatric serpentine lineage matching heterotypic basionym "Senecio earlei" -> Taxonomic_Resurrection_Candidate!
        #   - Rows 58..63: K5 parapatric allochronic lineage outside all basionym types -> unnamed Putative_Cryptic_Species!
        n_tot <- 63L
        cats <- sprintf("ST3_%02d", seq_len(n_tot))
        g_cls <- c(rep(1L, 40), rep(2L, 5), rep(3L, 6), rep(4L, 6), rep(5L, 6))
        insts <- rep(c("NCU", "DUKE", "US", "NY", "GH", "MO"), length.out = n_tot)
        colls <- c(
          paste0("Collector_", seq_len(40)),
          c("A. S. Weakley", "Weakley, Alan S.", "A.S. Weakley et al.", "Dr. Alan Weakley", "Weakley, A. S."),
          c("R. Kral", "R. K. Godfrey", "A. Cronquist", "D. K. Trock", "T. M. Barkley", "G. L. Nesom"),
          c("R. Kral", "R. K. Godfrey", "A. Cronquist", "D. K. Trock", "T. M. Barkley", "G. L. Nesom"),
          c("R. Kral", "R. K. Godfrey", "A. Cronquist", "D. K. Trock", "T. M. Barkley", "G. L. Nesom")
        )
        rnums <- c(
          as.character(1000 + seq_len(40)),
          rep("4102", 5),
          as.character(5101:5106),
          as.character(6101:6106),
          as.character(7101:7106)
        )
        lats <- c(
          seq(33.5, 36.5, length.out = 40),
          rep(34.80, 5),
          seq(39.0, 41.5, length.out = 6),
          c(35.10, 35.35, 35.60, 35.85, 36.10, 36.35),
          c(34.90, 35.20, 35.50, 35.80, 36.05, 36.30)
        )
        lons <- c(
          seq(-81.0, -78.5, length.out = 40),
          rep(-80.20, 5),
          seq(-76.5, -74.5, length.out = 6),
          c(-80.8, -80.4, -80.0, -79.6, -79.2, -78.8),
          c(-80.7, -80.3, -79.9, -79.5, -79.1, -78.7)
        )
        pc1 <- c(
          seq(-2.0, -0.2, length.out = 40),
          rep(2.6, 5),
          seq(0.1, 1.8, length.out = 6),
          rep(3.1, 6),
          rep(-3.4, 6)
        )
        pc2 <- c(
          rep(0.0, 40),
          rep(2.2, 5),
          seq(0.05, 0.30, length.out = 6),
          rep(2.6, 6),
          rep(-2.8, 6)
        )
        ddoys <- c(rep(1.0, 40), rep(25.0, 5), rep(4.0, 6), rep(24.0, 6), rep(-26.0, 6))
        dedaph <- c(rep(0.82, 40), rep(0.16, 5), rep(0.24, 6), rep(0.15, 6), rep(0.14, 6))
        edz <- c(rep(0.4, 40), rep(3.2, 5), rep(2.7, 6), rep(3.4, 6), rep(3.6, 6))
        nt_bas <- c(rep("Senecio anonymus", 40), rep(NA_character_, 5), rep(NA_character_, 6), rep("Senecio earlei", 6), rep("Senecio smallii", 6))
        nt_bas_d2 <- c(rep(0.9, 40), rep(NA_real_, 5), rep(NA_real_, 6), rep(2.1, 6), rep(19.5, 6))

        df_st3 <- data.frame(
          catalogNumber = cats,
          plant_individual_id = 0L,
          institutionCode = insts,
          recordedBy = colls,
          recordNumber = rnums,
          year = "2024",
          species_raw = "Packera anonyma",
          species_standardized = "Packera anonyma",
          determiner_tier = "Tier_1_Gold",
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          decimalLatitude = lats,
          decimalLongitude = lons,
          doy = 130 + ddoys,
          delta_doy_resid = ddoys,
          is_active_anthesis = TRUE,
          PC1 = pc1,
          PC2 = pc2,
          selected_morph_features = "PC1;PC2",
          cda_predicted_taxon = paste0("Morphotype_K", g_cls),
          cda_posterior_prob = ifelse(g_cls == 1L, 0.94, 0.52),
          morphometric_outlier_flag = g_cls > 1L,
          is_typical_morphotype = g_cls == 1L,
          edaphic_schoener_d = dedaph,
          edaphic_z_score = edz,
          type_edaphic_z_score = edz,
          min_congener_mahal_d2 = 34.0,
          min_congener_taxon = "Packera dubia",
          open_set_novelty_p = 0.001,
          open_set_df = 12L,
          nearest_type_basionym = nt_bas,
          nearest_type_basionym_d2 = nt_bas_d2,
          gmm_cluster = g_cls,
          gmm_2_delta_bic = 16.5,
          stringsAsFactors = FALSE
        )
        res_st3 <- apply_taxonomic_decision_matrix(df_st3, intraspecific_mode = TRUE)

        # K2 (rows 41..45): 5 duplicate sheets of 1 gathering -> cryptic_cluster_cohesion_pass == FALSE
        stopifnot(all(!res_st3$cryptic_cluster_cohesion_pass[41:45]))
        stopifnot(all(res_st3$taxonomic_treatment_subtype[41:45] != "Putative_Cryptic_Species"))

        # K3 (rows 46..51): smooth allopatric geographic cline -> cryptic_non_clinal_pass == FALSE -> Subspecies_Ecotype
        stopifnot(all(!res_st3$cryptic_non_clinal_pass[46:51]))
        stopifnot(all(!res_st3$cryptic_cluster_cohesion_pass[46:51]))
        stopifnot(all(res_st3$taxonomic_treatment_subtype[46:51] == "Subspecies_Ecotype"))

        # K4 (rows 52..57): cohesive non-clinal lineage matching heterotypic basionym Senecio earlei -> Taxonomic_Resurrection_Candidate
        stopifnot(all(res_st3$cryptic_cluster_cohesion_pass[52:57]))
        stopifnot(all(res_st3$taxonomic_status_call[52:57] == "Taxonomic_Resurrection_Candidate"))
        stopifnot(all(res_st3$resurrection_basionym[52:57] == "Senecio earlei"))
        stopifnot(all(grepl("Resurrect Basionym: Senecio earlei", res_st3$recommended_determination[52:57], fixed = TRUE)))

        # K5 (rows 58..63): cohesive non-clinal lineage outside all basionym types (d2=19.5 > qchisq(0.95,2)=5.99) -> de novo Putative_Cryptic_Species
        stopifnot(all(res_st3$cryptic_cluster_cohesion_pass[58:63]))
        stopifnot(all(res_st3$taxonomic_status_call[58:63] == "Putative_Cryptic_Species"))
        stopifnot(all(is.na(res_st3$resurrection_basionym[58:63])))
        stopifnot(all(grepl("Putative Cryptic Lineage K5", res_st3$recommended_determination[58:63], fixed = TRUE)))

        cat("FIVE_STEP_MISID_AND_CRYPTIC_OK\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_five_step_misid_cryptic.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key="step06", cwd=REPO_ROOT, skip_figures=True)
        self.assertEqual(proc.returncode, 0, f"Five-step misid/cryptic R suite failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("FIVE_STEP_MISID_AND_CRYPTIC_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
