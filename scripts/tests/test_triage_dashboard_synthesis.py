#!/usr/bin/env python3
"""
===============================================================================
Script: test_triage_dashboard_synthesis.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 07 (scripts/analysis/07_triage_dashboard_synthesis.R):
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
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


class TestTriageDashboardSynthesis(unittest.TestCase):
    """Test suite for Step 07 active R script (scripts/analysis/07_triage_dashboard_synthesis.R)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.r_unit_proc = None

        if shutil.which("Rscript"):
            cls.prec_in_csv = cls._class_tmp_path / "prec_in.csv"
            cls.prec_out_csv = cls._class_tmp_path / "prec_out.csv"
            pd.DataFrame([
                {
                    "catalogNumber": "STAGE1-MISID-OVER-HYB",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_2_Silver",
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_posterior_prob": 0.58,
                    "cda_second_posterior_prob": 0.40,
                    "cda_posterior_margin": 0.18,
                    "vision_predicted_label": "Packera anonyma",
                    "c_error": 0.82,
                    "is_label_corrupted": True,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.55,
                    "delta_doy_resid": 3.0,
                    "gmm_cluster": 2,
                    "gmm_uncertainty": 0.42,
                    "n_tier1_leaves": 3,
                    "morph_tier_source": "Tier_1_Pure",
                },
                {
                    "catalogNumber": "STAGE2-MOSAIC-HYBRID",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_2_Silver",
                    "cda_predicted_taxon": "Packera anonyma",
                    "cda_posterior_prob": 0.65,
                    "cda_second_posterior_prob": 0.30,
                    "cda_posterior_margin": 0.35,
                    "vision_predicted_label": "Packera paupercula",
                    "c_error": 0.72,
                    "is_label_corrupted": True,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.52,
                    "delta_doy_resid": 4.0,
                    "gmm_cluster": 2,
                    "gmm_uncertainty": 0.38,
                    "n_tier1_leaves": 3,
                    "morph_tier_source": "Tier_1_Pure",
                },
                {
                    "catalogNumber": "STAGE3-TRUE-CRYPTIC",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "cda_predicted_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.46,
                    "cda_second_posterior_prob": 0.28,
                    "cda_posterior_margin": 0.18,
                    "vision_predicted_label": "Packera dubia",
                    "c_error": 0.32,
                    "is_label_corrupted": False,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.18,
                    "delta_doy_resid": 24.0,
                    "gmm_cluster": 5,
                    "gmm_2_delta_bic": 12.4,
                    "gmm_uncertainty": 0.08,
                    "n_tier1_leaves": 4,
                    "morph_tier_source": "Tier_1_Pure",
                },
                {
                    "catalogNumber": "STAGE4-NO-BARRIER-PLASTIC",
                    "species_standardized": "Packera dubia",
                    "determiner_tier": "Tier_1_Gold",
                    "cda_predicted_taxon": "Packera dubia",
                    "cda_posterior_prob": 0.48,
                    "cda_second_posterior_prob": 0.25,
                    "cda_posterior_margin": 0.23,
                    "vision_predicted_label": "Packera dubia",
                    "c_error": 0.25,
                    "is_label_corrupted": False,
                    "morphometric_outlier_flag": True,
                    "edaphic_best_fit_taxon": "Packera dubia",
                    "edaphic_schoener_d": 0.68,
                    "delta_doy_resid": 3.0,
                    "gmm_cluster": 5,
                    "gmm_2_delta_bic": 12.4,
                    "gmm_uncertainty": 0.08,
                    "n_tier1_leaves": 4,
                    "morph_tier_source": "Tier_1_Pure",
                },
            ]).to_csv(cls.prec_in_csv, index=False)

            cls.nd_dir = cls._class_tmp_path / "non_dubia"
            cls.nd_dir.mkdir(parents=True, exist_ok=True)
            cls.nd_in_csv = cls.nd_dir / "non_dubia_input.csv"
            cls.nd_out_csv = cls.nd_dir / "non_dubia_triage.csv"
            cls.nd_multi_report = cls.nd_dir / "multi_species_report.md"
            cls.nd_focal_report_default = cls.nd_dir / "Packera_Taxonomic_Revision_Summary.md"
            cls.nd_focal_plot_default = cls.nd_dir / "Figure_Integrative_Packera_Revision.pdf"
            cls.nd_helper_check_csv = cls.nd_dir / "pub_fig_helpers.csv"
            cls.dubia_dir = cls.nd_dir / "dubia_focal"
            cls.dubia_dir.mkdir(parents=True, exist_ok=True)
            dubia_rep_def = cls.dubia_dir / "Packera_Taxonomic_Revision_Summary.md"
            dubia_fig_def = cls.dubia_dir / "Figure_Integrative_Packera_Revision.pdf"
            dubia_fig_def.write_bytes(b"%PDF-1.4\n%%EOF\n")

            nd_rows = [
                {"catalogNumber": "AUR-GOLD-1", "institutionCode": "NCU", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "regional_group": "Blue Ridge", "doy": 120, "decimalLatitude": 35.6, "decimalLongitude": -82.5, "fine_soil_ph": 6.2, "fine_soil_sand": 32.0, "fine_bedrock_depth_cm": 95.0, "can1": -3.0, "can2": 0.1, "can3": 0.0, "PC1": -1.5, "PC2": 0.2, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "vision_predicted_label": "Packera aurea", "c_error": 0.04, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.82, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "AUR-GOLD-2", "institutionCode": "DUKE", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "regional_group": "Ridge and Valley", "doy": 132, "decimalLatitude": 36.2, "decimalLongitude": -81.7, "fine_soil_ph": 6.4, "fine_soil_sand": 36.0, "fine_bedrock_depth_cm": 110.0, "can1": -2.8, "can2": -0.1, "can3": 0.1, "PC1": -1.4, "PC2": -0.1, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.95, "cda_posterior_margin": 0.90, "vision_predicted_label": "Packera aurea", "c_error": 0.05, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.80, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "OBO-GOLD-1", "institutionCode": "NCU", "species_raw": "Packera obovata", "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "regional_group": "Interior Plateau", "doy": 105, "decimalLatitude": 36.0, "decimalLongitude": -86.8, "fine_soil_ph": 7.3, "fine_soil_sand": 22.0, "fine_bedrock_depth_cm": 38.0, "can1": 3.1, "can2": 0.2, "can3": -0.1, "PC1": 1.6, "PC2": 0.1, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.97, "cda_posterior_margin": 0.94, "vision_predicted_label": "Packera obovata", "c_error": 0.03, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera obovata", "edaphic_schoener_d": 0.85, "gmm_cluster": 1, "gmm_uncertainty": 0.02, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 4},
                {"catalogNumber": "OBO-GOLD-2", "institutionCode": "TENN", "scientificName": "Packera obovata (Muhl. ex Willd.) W.A.Weber & Á.Löve", "species_raw": "Packera obovata", "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "regional_group": "Interior Plateau", "doy": 115, "decimalLatitude": 37.1, "decimalLongitude": -85.9, "fine_soil_ph": 7.5, "fine_soil_sand": 26.0, "fine_bedrock_depth_cm": 42.0, "can1": 2.9, "can2": -0.2, "can3": 0.1, "PC1": 1.5, "PC2": -0.2, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "vision_predicted_label": "Packera obovata", "c_error": 0.04, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera obovata", "edaphic_schoener_d": 0.83, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "SCH-GOLD-1", "institutionCode": "GH", "species_raw": "Packera schweinitziana", "species_standardized": "Packera schweinitziana", "determiner_tier": "Tier_1_Gold", "regional_group": "Northeastern Highlands", "doy": 168, "decimalLatitude": 44.2, "decimalLongitude": -71.3, "fine_soil_ph": 5.4, "fine_soil_sand": 48.0, "fine_bedrock_depth_cm": 125.0, "can1": 0.0, "can2": 3.4, "can3": 0.2, "PC1": 0.1, "PC2": 1.9, "cda_predicted_taxon": "Packera schweinitziana", "cda_posterior_prob": 0.97, "cda_posterior_margin": 0.94, "vision_predicted_label": "Packera schweinitziana", "c_error": 0.03, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera schweinitziana", "edaphic_schoener_d": 0.84, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "SCH-GOLD-2", "institutionCode": "NY", "species_raw": "Packera schweinitziana", "species_standardized": "Packera schweinitziana", "determiner_tier": "Tier_1_Gold", "regional_group": "Northern Appalachian and Atlantic Maritime Highlands", "doy": 176, "decimalLatitude": 44.6, "decimalLongitude": -72.8, "fine_soil_ph": 5.2, "fine_soil_sand": 52.0, "fine_bedrock_depth_cm": 130.0, "can1": 0.2, "can2": 3.2, "can3": -0.2, "PC1": 0.2, "PC2": 1.8, "cda_predicted_taxon": "Packera schweinitziana", "cda_posterior_prob": 0.96, "cda_posterior_margin": 0.92, "vision_predicted_label": "Packera schweinitziana", "c_error": 0.04, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera schweinitziana", "edaphic_schoener_d": 0.81, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "AUR-MISID-OBO", "institutionCode": "NCU", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_3_Bronze", "regional_group": "Interior Plateau", "doy": 110, "decimalLatitude": 36.4, "decimalLongitude": -86.2, "fine_soil_ph": 7.4, "fine_soil_sand": 24.0, "fine_bedrock_depth_cm": 40.0, "can1": 2.95, "can2": 0.05, "can3": 0.0, "PC1": 1.55, "PC2": 0.0, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.93, "cda_posterior_margin": 0.86, "vision_predicted_label": "Packera obovata", "c_error": 0.91, "is_label_corrupted": True, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera obovata", "edaphic_schoener_d": 0.18, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "morphometric_outlier_flag": True, "triage_category": "Severe_Triple_Stream_Conflict", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "AUR-HYB-EXPLICIT", "institutionCode": "DUKE", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Ridge and Valley", "doy": 118, "decimalLatitude": 36.1, "decimalLongitude": -83.4, "fine_soil_ph": 6.8, "fine_soil_sand": 29.0, "fine_bedrock_depth_cm": 68.0, "can1": 0.1, "can2": 0.0, "can3": 0.0, "PC1": 0.05, "PC2": 0.0, "cda_predicted_taxon": "Packera obovata", "cda_posterior_prob": 0.52, "cda_posterior_margin": 0.08, "vision_predicted_label": "Packera aurea", "c_error": 0.35, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.55, "gmm_cluster": 1, "gmm_uncertainty": 0.40, "morphometric_outlier_flag": False, "triage_category": "Putative_Hybrid_Zone_Intergrade", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "AUR-HYB-MAHAL", "institutionCode": "GH", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Northeastern Highlands", "doy": 122, "decimalLatitude": 42.8, "decimalLongitude": -73.2, "fine_soil_ph": 5.9, "fine_soil_sand": 40.0, "fine_bedrock_depth_cm": 105.0, "can1": -0.9, "can2": 2.4, "can3": 0.1, "PC1": -0.4, "PC2": 1.3, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.51, "cda_posterior_margin": 0.10, "vision_predicted_label": "Packera aurea", "c_error": 0.30, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.58, "gmm_cluster": 1, "gmm_uncertainty": 0.42, "morphometric_outlier_flag": False, "triage_category": "Putative_Hybrid_Zone_Intergrade", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "AUR-ECOTYPE-VAR", "institutionCode": "NCU", "species_raw": "Packera aurea var. gracilis", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Blue Ridge", "doy": 125, "decimalLatitude": 35.9, "decimalLongitude": -82.1, "fine_soil_ph": 6.1, "fine_soil_sand": 38.0, "fine_bedrock_depth_cm": 88.0, "can1": -2.5, "can2": 0.4, "can3": 0.0, "PC1": -1.2, "PC2": 0.3, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.88, "cda_posterior_margin": 0.76, "vision_predicted_label": "Packera aurea", "c_error": 0.12, "is_label_corrupted": False, "tomentum_taxon_discordance": False, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.38, "gmm_cluster": 2, "gmm_uncertainty": 0.08, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
                {"catalogNumber": "AUR-PLASTIC-LATE", "institutionCode": "TENN", "species_raw": "Packera aurea", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver", "regional_group": "Blue Ridge", "doy": 148, "decimalLatitude": 35.7, "decimalLongitude": -82.4, "fine_soil_ph": 6.3, "fine_soil_sand": 35.0, "fine_bedrock_depth_cm": 98.0, "can1": -2.7, "can2": 0.05, "can3": 0.0, "PC1": -1.35, "PC2": 0.05, "cda_predicted_taxon": "Packera aurea", "cda_posterior_prob": 0.86, "cda_posterior_margin": 0.72, "vision_predicted_label": "Packera obovata", "c_error": 0.58, "is_label_corrupted": False, "tomentum_taxon_discordance": True, "edaphic_best_fit_taxon": "Packera aurea", "edaphic_schoener_d": 0.78, "gmm_cluster": 1, "gmm_uncertainty": 0.06, "morphometric_outlier_flag": False, "triage_category": "Clean_MultiModal_Consensus", "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3},
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
                  "VIS_CONTRADICTED_BY_MORPH",
                  "MORPH_MOD_CONTRADICTED_BY_VIS",
                  "CORROBORATED_MISID"
                ),
                plant_individual_id = 1L,
                institutionCode = "NCU",
                species_raw = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera aurea", "Packera aurea"),
                species_standardized = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera aurea", "Packera aurea"),
                determiner_tier = c(rep("Tier_1_Gold", 9), rep("Tier_3_Bronze", 3)),
                morph_tier_source = "Tier_1_Pure",
                n_tier1_leaves = 3L,
                cda_predicted_taxon = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera obovata", "Packera obovata"),
                cda_second_taxon = c(rep("Packera obovata", 3), rep("Packera aurea", 3), rep("Packera aurea", 3), "Packera obovata", "Packera aurea", "Packera aurea"),
                cda_posterior_prob = c(rep(0.95, 9), 0.91, 0.76, 0.92),
                cda_posterior_margin = c(rep(0.90, 9), 0.82, 0.52, 0.84),
                morphometric_outlier_flag = c(rep(FALSE, 9), FALSE, TRUE, TRUE),
                vision_predicted_label = c(
                  rep("Packera aurea", 3),
                  rep("Packera obovata", 3),
                  "Packera obovata", "Packera obovata", "Packera schweinitziana",
                  "Packera obovata",
                  "Packera aurea",
                  "Packera obovata"
                ),
                c_error = c(rep(0.05, 6), 0.85, 0.85, 0.05, 0.88, 0.06, 0.92),
                is_label_corrupted = c(rep(FALSE, 6), TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
                is_artifact_driven = FALSE,
                tomentum_taxon_discordance = FALSE,
                vision_is_ood = FALSE,
                vision_ood_score = 0.20,
                edaphic_best_fit_taxon = c(rep("Packera aurea", 3), rep("Packera obovata", 3), rep("Packera schweinitziana", 3), "Packera aurea", "Packera aurea", "Packera obovata"),
                edaphic_schoener_d = c(rep(0.82, 9), 0.80, 0.79, 0.18),
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
              row_vis_contra <- res7[res7$catalogNumber == "VIS_CONTRADICTED_BY_MORPH", ]
              row_morph_contra <- res7[res7$catalogNumber == "MORPH_MOD_CONTRADICTED_BY_VIS", ]
              row_corrob <- res7[res7$catalogNumber == "CORROBORATED_MISID", ]

              stopifnot(row_vis_contra$taxonomic_treatment_subtype != "Misidentification_Reassign")
              stopifnot(row_vis_contra$recommended_determination == "Packera aurea")
              stopifnot(row_morph_contra$taxonomic_treatment_subtype != "Misidentification_Reassign")
              stopifnot(row_morph_contra$recommended_determination == "Packera aurea")
              stopifnot(row_corrob$taxonomic_treatment_subtype == "Misidentification_Reassign")
              stopifnot(row_corrob$recommended_determination == "Packera obovata")
              stopifnot(all(res7$probabilistic_consensus_conflict >= 0.0 & res7$probabilistic_consensus_conflict <= 1.0))

              anc_m <- .resolve_step07_anchor_mask(df7)
              prec_vis <- .compute_oof_class_precision_map(df7$species_standardized[anc_m], df7$vision_predicted_label[anc_m])
              prec_mor <- .compute_oof_class_precision_map(df7$species_standardized[anc_m], df7$cda_predicted_taxon[anc_m])
              stopifnot(prec_vis[["Packera obovata"]] < prec_mor[["Packera obovata"]])

              df7_hi <- df7
              df7_hi$vision_predicted_label[7:9] <- "Packera schweinitziana"
              df7_hi$c_error[7:9] <- 0.05
              df7_hi$is_label_corrupted[7:9] <- FALSE
              res7_hi <- apply_taxonomic_decision_matrix(df7_hi)
              c_low_prec <- res7$probabilistic_consensus_conflict[res7$catalogNumber == "VIS_CONTRADICTED_BY_MORPH"]
              c_hi_prec <- res7_hi$probabilistic_consensus_conflict[res7_hi$catalogNumber == "VIS_CONTRADICTED_BY_MORPH"]
              stopifnot(c_low_prec < c_hi_prec)
              cat("STEP07_PRECISION_AND_GATE_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step07_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = subprocess.run(
                ["Rscript", str(r_batch_file)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                timeout=90,
                check=False,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 07 R unit suite failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_step07_static_contracts_phenology_and_tier1_pure_gating(self):
        """Verify 07_triage_dashboard_synthesis.R retains safe coalescing, Hopkins circular phenology, Tier_1_Pure gating, and dynamic species-agnostic contracts."""
        self.assertTrue(self.r_script_path.exists())
        r_text = self.r_script_path.read_text(encoding="utf-8")
        for token in (
            'unique(c("catalogNumber"',
            "safe_left_join_and_coalesce",
            "morphometric_outlier_flag",
            "gmm_bic_model_selection.csv",
            "--focal-taxon",
            "--intraspecific-mode",
            "outputs/figures/Figure_Integrative_Packera_Revision.pdf",
            "outputs/reports/Packera_Taxonomic_Revision_Summary.md",
            "morph_tier_source",
            "n_tier1_leaves",
            "is_tier1_pure_morph",
            "is_tier1_pure_i",
            "compute_hopkins_adjusted_phenology",
            "hopkins_adjusted_doy",
            "sin_doy",
            "cos_doy",
            "atan2",
            "n_cryptic_modalities >= 2L",
            "Padial et al. 2010",
            "Dawid & Skene 1979",
            "Schlick-Steiner et al. 2010",
            "Pewsey et al. 2013",
            "audit_institutional_batch_effects_ancova",
            "PC1 ~ species_standardized + institutionCode",
            "partial_eta_sq",
        ):
            self.assertIn(token, r_text, f"Missing expected Step 07 token: {token}")

        for forbidden in (
            "TARGET_TAXA",
            "TAXON_COLORS",
            "nomen_citations",
            "Packera_dubia_Taxonomic_Revision_Summary.md",
            "Figure_Integrative_Packera_dubia_Revision.pdf",
            "Canonical Axis 1 (74.2%)",
            "Canonical Axis 2 (18.6%)",
            "10,815.5",
        ):
            self.assertNotIn(forbidden, r_text)

        r06_text = (PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R").read_text(encoding="utf-8")
        self.assertIn("morph_tier_source", r06_text)
        self.assertIn("n_tier1_leaves", r06_text)

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
        self.assertEqual(res.loc["AUR-HYB-MAHAL", "taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertEqual(res.loc["AUR-HYB-MAHAL", "recommended_determination"], "Packera aurea x Packera schweinitziana")
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
        focal_named_pdf = self.nd_dir / "Figure_Integrative_Packera_aurea_Revision.pdf"
        self.assertTrue(focal_named_md.exists(), "Expected taxon-named focal Markdown report copy")
        self.assertTrue(focal_named_pdf.exists(), "Expected taxon-named focal PDF synthesis plate copy")
        self.assertTrue((self.dubia_dir / "Packera_dubia_Taxonomic_Revision_Summary.md").exists())
        self.assertTrue((self.dubia_dir / "Figure_Integrative_Packera_dubia_Revision.pdf").exists())
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

    @unittest.skipUnless(shutil.which("Rscript") is not None, "Rscript required for full cross-stage handoff test")
    def test_step07_cli_handoff_coalescing_multi_plant_artifact_gate_and_ood_cryptic(self):
        """Verify Step 07 CLI preserves multi-plant keys, coalesces ungeoreferenced vouchers (missing-modality renormalization), assigns RESOLVED/CRITICAL priorities, gates Stage 1 on !is_artifact_driven, and corroborates Stage 3 via vision_is_ood."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vouchers_csv = tmp / "curated_vouchers.csv"
            pheno_csv = tmp / "voucher_phenological_states.csv"
            morph_csv = tmp / "morphometric_flags.csv"
            vision_csv = tmp / "label_noise_audit.csv"
            mm_csv = tmp / "multimodal_conflict_flags.csv"
            niche_csv = tmp / "niche_identity_tests.csv"
            out_queue = tmp / "triage_queue.csv"
            out_triage = tmp / "final_curation_triage.csv"
            out_rep = tmp / "Packera_Taxonomic_Revision_Summary.md"
            out_fig = tmp / "Figure_Integrative_Packera_Revision.pdf"

            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.5, "decimalLongitude": -82.5, "doy": 118},
                {"catalogNumber": "CAT-ANCHOR-A2", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.55, "decimalLongitude": -82.45, "doy": 119},
                {"catalogNumber": "CAT-ANCHOR-B1", "institutionCode": "NCU", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 36.0, "decimalLongitude": -79.0, "doy": 112},
                {"catalogNumber": "CAT-ANCHOR-B2", "institutionCode": "NCU", "species_raw": "Packera obovata", "scientificName": "Packera obovata", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 36.05, "decimalLongitude": -79.05, "doy": 113},
                {"catalogNumber": "CAT-SEVERE-MISID", "institutionCode": "NCU", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 36.1, "decimalLongitude": -79.1, "doy": 114},
                {"catalogNumber": "NCU-HIST-2", "institutionCode": "PH", "species_raw": "Packera aurea", "scientificName": "Packera aurea", "determiner_tier": "Tier_2_Silver", "decimalLatitude": np.nan, "decimalLongitude": np.nan, "doy": np.nan},
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
                {"catalogNumber": "NCU-HIST-2", "voucher_phenological_state": "sterile", "n_anthesis": 0, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.80},
                {"catalogNumber": "CAT-ART-SHORTCUT", "voucher_phenological_state": "anthesis", "n_anthesis": 4, "n_bud": 1, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.88},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "voucher_phenological_state": "anthesis", "n_anthesis": 7, "n_bud": 0, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.95},
                {"catalogNumber": "CAT-MULTI-SHEET", "voucher_phenological_state": "anthesis", "n_anthesis": 8, "n_bud": 2, "n_fruiting": 0, "n_senescent": 0, "reproductive_stage_confidence": 0.90},
            ]).to_csv(pheno_csv, index=False)

            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": -2.0, "PC2": 0.1, "Sym_PC1": -1.9, "Sym_PC2": 0.1, "can1": -2.5, "can2": 0.2, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.03, "petiole_length_mm": 30.0, "blade_length_mm": 25.0, "petiole_to_blade_ratio": 1.20, "plant_height_cm": 44.0, "scape_to_basal_leaf_ratio": 8.0, "is_bent_stem": False, "capitula_count": 7, "involucre_height_mm": 6.4, "involucre_width_mm": 5.1, "AR_inv": 1.25, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.88, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-ANCHOR-A2", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.97, "cda_second_posterior_prob": 0.03, "cda_posterior_margin": 0.94, "PC1": -2.1, "PC2": 0.15, "Sym_PC1": -2.0, "Sym_PC2": 0.15, "can1": -2.6, "can2": 0.15, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.02, "petiole_length_mm": 31.0, "blade_length_mm": 25.5, "petiole_to_blade_ratio": 1.22, "plant_height_cm": 45.0, "scape_to_basal_leaf_ratio": 8.1, "is_bent_stem": False, "capitula_count": 7, "involucre_height_mm": 6.5, "involucre_width_mm": 5.1, "AR_inv": 1.27, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.89, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B1", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.95, "cda_second_posterior_prob": 0.05, "cda_posterior_margin": 0.90, "PC1": 2.1, "PC2": -0.2, "Sym_PC1": 2.0, "Sym_PC2": -0.2, "can1": 2.6, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.04, "petiole_length_mm": 11.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": 0.39, "plant_height_cm": 36.0, "scape_to_basal_leaf_ratio": 9.2, "is_bent_stem": False, "capitula_count": 5, "involucre_height_mm": 5.8, "involucre_width_mm": 4.8, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.85, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B2", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": 2.2, "PC2": -0.15, "Sym_PC1": 2.1, "Sym_PC2": -0.15, "can1": 2.7, "can2": -0.15, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.03, "petiole_length_mm": 10.8, "blade_length_mm": 27.5, "petiole_to_blade_ratio": 0.39, "plant_height_cm": 35.5, "scape_to_basal_leaf_ratio": 9.1, "is_bent_stem": False, "capitula_count": 6, "involucre_height_mm": 5.9, "involucre_width_mm": 4.9, "AR_inv": 1.20, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.87, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-SEVERE-MISID", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.96, "cda_second_posterior_prob": 0.04, "cda_posterior_margin": 0.92, "PC1": 2.15, "PC2": -0.2, "Sym_PC1": 2.05, "Sym_PC2": -0.2, "can1": 2.65, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.03, "petiole_length_mm": 11.2, "blade_length_mm": 28.0, "petiole_to_blade_ratio": 0.40, "plant_height_cm": 36.0, "scape_to_basal_leaf_ratio": 9.0, "is_bent_stem": False, "capitula_count": 5, "involucre_height_mm": 5.8, "involucre_width_mm": 4.8, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.86, "morphometric_outlier_flag": True},
                {"catalogNumber": "NCU-HIST-2", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 2, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.93, "cda_second_posterior_prob": 0.07, "cda_posterior_margin": 0.86, "PC1": 2.0, "PC2": -0.2, "Sym_PC1": 1.9, "Sym_PC2": -0.2, "can1": 2.4, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.04, "petiole_length_mm": 11.5, "blade_length_mm": 27.5, "petiole_to_blade_ratio": 0.42, "plant_height_cm": 35.0, "scape_to_basal_leaf_ratio": 8.9, "is_bent_stem": False, "capitula_count": 0, "involucre_height_mm": np.nan, "involucre_width_mm": np.nan, "AR_inv": np.nan, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.82, "morphometric_outlier_flag": True},
                {"catalogNumber": "CAT-ART-SHORTCUT", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.91, "cda_second_posterior_prob": 0.09, "cda_posterior_margin": 0.82, "PC1": -1.8, "PC2": 0.2, "Sym_PC1": -1.7, "Sym_PC2": 0.2, "can1": -2.3, "can2": 0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.05, "petiole_length_mm": 29.0, "blade_length_mm": 24.0, "petiole_to_blade_ratio": 1.21, "plant_height_cm": 41.0, "scape_to_basal_leaf_ratio": 7.6, "is_bent_stem": False, "capitula_count": 5, "involucre_height_mm": 6.2, "involucre_width_mm": 5.0, "AR_inv": 1.24, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.82, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 3, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.48, "cda_second_posterior_prob": 0.15, "cda_posterior_margin": 0.33, "PC1": 0.0, "PC2": 2.4, "Sym_PC1": 0.1, "Sym_PC2": 2.3, "can1": 0.0, "can2": 2.5, "can3": 1.4, "gmm_cluster": 2, "gmm_uncertainty": 0.08, "petiole_length_mm": 19.0, "blade_length_mm": 22.0, "petiole_to_blade_ratio": 0.86, "plant_height_cm": 33.0, "scape_to_basal_leaf_ratio": 8.0, "is_bent_stem": True, "capitula_count": 7, "involucre_height_mm": 7.1, "involucre_width_mm": 5.3, "AR_inv": 1.34, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.74, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-MULTI-SHEET", "plant_individual_id": 1, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 2, "cda_predicted_taxon": "Packera aurea", "cda_second_taxon": "Packera obovata", "cda_posterior_prob": 0.94, "cda_second_posterior_prob": 0.06, "cda_posterior_margin": 0.88, "PC1": -1.9, "PC2": 0.1, "Sym_PC1": -1.8, "Sym_PC2": 0.1, "can1": -2.4, "can2": 0.1, "can3": 0.0, "gmm_cluster": 1, "gmm_uncertainty": 0.04, "petiole_length_mm": 31.0, "blade_length_mm": 25.0, "petiole_to_blade_ratio": 1.24, "plant_height_cm": 43.0, "scape_to_basal_leaf_ratio": 7.7, "is_bent_stem": False, "capitula_count": 8, "involucre_height_mm": 6.3, "involucre_width_mm": 5.0, "AR_inv": 1.26, "ecoregion_cda_predicted": "Blue_Ridge", "ecoregion_cda_posterior_prob": 0.86, "morphometric_outlier_flag": False},
                {"catalogNumber": "CAT-MULTI-SHEET", "plant_individual_id": 2, "morph_tier_source": "Tier_1_Pure", "n_tier1_leaves": 2, "cda_predicted_taxon": "Packera obovata", "cda_second_taxon": "Packera aurea", "cda_posterior_prob": 0.93, "cda_second_posterior_prob": 0.07, "cda_posterior_margin": 0.86, "PC1": 2.0, "PC2": -0.2, "Sym_PC1": 1.9, "Sym_PC2": -0.2, "can1": 2.5, "can2": -0.2, "can3": 0.0, "gmm_cluster": 2, "gmm_uncertainty": 0.05, "petiole_length_mm": 10.5, "blade_length_mm": 27.0, "petiole_to_blade_ratio": 0.39, "plant_height_cm": 35.0, "scape_to_basal_leaf_ratio": 9.3, "is_bent_stem": False, "capitula_count": 6, "involucre_height_mm": 5.7, "involucre_width_mm": 4.7, "AR_inv": 1.21, "ecoregion_cda_predicted": "Piedmont", "ecoregion_cda_posterior_prob": 0.81, "morphometric_outlier_flag": True},
            ]).to_csv(morph_csv, index=False)

            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "vision_predicted_label": "Packera aurea", "confidence_predicted_class": 0.95, "vision_second_prob": 0.05, "vision_margin": 0.90, "c_error": 0.04, "is_label_corrupted": False, "vision_ood_score": 0.18, "vision_is_ood": False, "foreground_attention_ratio": 0.84, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 52.0, "mean_saturation": 0.31, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-ANCHOR-A2", "vision_predicted_label": "Packera aurea", "confidence_predicted_class": 0.96, "vision_second_prob": 0.04, "vision_margin": 0.92, "c_error": 0.03, "is_label_corrupted": False, "vision_ood_score": 0.17, "vision_is_ood": False, "foreground_attention_ratio": 0.85, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 51.0, "mean_saturation": 0.32, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-ANCHOR-B1", "vision_predicted_label": "Packera obovata", "confidence_predicted_class": 0.94, "vision_second_prob": 0.06, "vision_margin": 0.88, "c_error": 0.05, "is_label_corrupted": False, "vision_ood_score": 0.20, "vision_is_ood": False, "foreground_attention_ratio": 0.81, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 55.0, "mean_saturation": 0.29, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-ANCHOR-B2", "vision_predicted_label": "Packera obovata", "confidence_predicted_class": 0.95, "vision_second_prob": 0.05, "vision_margin": 0.90, "c_error": 0.04, "is_label_corrupted": False, "vision_ood_score": 0.19, "vision_is_ood": False, "foreground_attention_ratio": 0.83, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 54.0, "mean_saturation": 0.30, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-SEVERE-MISID", "vision_predicted_label": "Packera obovata", "confidence_predicted_class": 0.96, "vision_second_prob": 0.04, "vision_margin": 0.92, "c_error": 0.95, "is_label_corrupted": True, "vision_ood_score": 0.21, "vision_is_ood": False, "foreground_attention_ratio": 0.86, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 55.0, "mean_saturation": 0.29, "tomentum_taxon_discordance": False},
                {"catalogNumber": "NCU-HIST-2", "vision_predicted_label": "Packera obovata", "confidence_predicted_class": 0.92, "vision_second_prob": 0.08, "vision_margin": 0.84, "c_error": 0.91, "is_label_corrupted": True, "vision_ood_score": 0.22, "vision_is_ood": False, "foreground_attention_ratio": 0.80, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 56.0, "mean_saturation": 0.28, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-ART-SHORTCUT", "vision_predicted_label": "Packera obovata", "confidence_predicted_class": 0.92, "vision_second_prob": 0.08, "vision_margin": 0.84, "c_error": 0.91, "is_label_corrupted": True, "vision_ood_score": 0.25, "vision_is_ood": False, "foreground_attention_ratio": 0.42, "is_artifact_driven": True, "surface_orientation": "adaxial", "mean_luminance": 58.0, "mean_saturation": 0.28, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "vision_predicted_label": "Packera obovata", "confidence_predicted_class": 0.54, "vision_second_prob": 0.46, "vision_margin": 0.08, "c_error": 0.82, "is_label_corrupted": True, "vision_ood_score": 0.78, "vision_is_ood": True, "foreground_attention_ratio": 0.87, "is_artifact_driven": False, "surface_orientation": "abaxial", "mean_luminance": 74.0, "mean_saturation": 0.12, "tomentum_taxon_discordance": False},
                {"catalogNumber": "CAT-MULTI-SHEET", "vision_predicted_label": "Packera aurea", "confidence_predicted_class": 0.89, "vision_second_prob": 0.11, "vision_margin": 0.78, "c_error": 0.10, "is_label_corrupted": False, "vision_ood_score": 0.22, "vision_is_ood": False, "foreground_attention_ratio": 0.79, "is_artifact_driven": False, "surface_orientation": "adaxial", "mean_luminance": 53.0, "mean_saturation": 0.30, "tomentum_taxon_discordance": False},
            ]).to_csv(vision_csv, index=False)

            # Note: NCU-HIST-2 lacks GPS and is intentionally omitted from mm_csv to test non-georeferenced coalescing & missing-modality renormalization!
            pd.DataFrame([
                {"catalogNumber": "CAT-ANCHOR-A1", "plant_individual_id": 1, "fine_soil_ph": 6.1, "fine_soil_sand": 34.0, "fine_bedrock_depth_cm": 110.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.35, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.92, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.78, "delta_doy": 1.0, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-ANCHOR-A2", "plant_individual_id": 1, "fine_soil_ph": 6.15, "fine_soil_sand": 33.5, "fine_bedrock_depth_cm": 112.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.32, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.93, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.78, "delta_doy": 1.2, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B1", "plant_individual_id": 1, "fine_soil_ph": 6.7, "fine_soil_sand": 48.0, "fine_bedrock_depth_cm": 55.0, "ssurgo_rock_outcrop_pct": 18.0, "edaphic_z_score": 0.41, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.90, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.76, "delta_doy": -2.0, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-ANCHOR-B2", "plant_individual_id": 1, "fine_soil_ph": 6.75, "fine_soil_sand": 49.0, "fine_bedrock_depth_cm": 54.0, "ssurgo_rock_outcrop_pct": 19.0, "edaphic_z_score": 0.39, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.91, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.76, "delta_doy": -1.8, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-SEVERE-MISID", "plant_individual_id": 1, "fine_soil_ph": 6.8, "fine_soil_sand": 49.0, "fine_bedrock_depth_cm": 52.0, "ssurgo_rock_outcrop_pct": 20.0, "edaphic_z_score": 2.9, "edaphic_best_fit_taxon": "Packera obovata", "rf_predicted_taxon": "Packera obovata", "rf_predicted_prob": 0.93, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.21, "delta_doy": -2.0, "triage_category": "Severe_Triple_Stream_Conflict", "multimodal_conflict_flag": True},
                {"catalogNumber": "CAT-ART-SHORTCUT", "plant_individual_id": 1, "fine_soil_ph": 6.0, "fine_soil_sand": 35.0, "fine_bedrock_depth_cm": 105.0, "ssurgo_rock_outcrop_pct": 3.0, "edaphic_z_score": 0.44, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.89, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.75, "delta_doy": 2.0, "triage_category": "Severe_Triple_Stream_Conflict", "multimodal_conflict_flag": True},
                {"catalogNumber": "CAT-VIS-OOD-CRYPT", "plant_individual_id": 1, "fine_soil_ph": 4.4, "fine_soil_sand": 72.0, "fine_bedrock_depth_cm": 18.0, "ssurgo_rock_outcrop_pct": 64.0, "edaphic_z_score": 3.45, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.51, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.14, "delta_doy": 26.0, "triage_category": "Edaphic_Envelope_Violation", "multimodal_conflict_flag": True},
                {"catalogNumber": "CAT-MULTI-SHEET", "plant_individual_id": 1, "fine_soil_ph": 6.1, "fine_soil_sand": 33.0, "fine_bedrock_depth_cm": 112.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.30, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.91, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.77, "delta_doy": 1.5, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
                {"catalogNumber": "CAT-MULTI-SHEET", "plant_individual_id": 2, "fine_soil_ph": 6.1, "fine_soil_sand": 33.0, "fine_bedrock_depth_cm": 112.0, "ssurgo_rock_outcrop_pct": 2.0, "edaphic_z_score": 0.30, "edaphic_best_fit_taxon": "Packera aurea", "rf_predicted_taxon": "Packera aurea", "rf_predicted_prob": 0.91, "spatial_block_cv_acc": 0.88, "is_verified_anthesis": True, "edaphic_schoener_d": 0.77, "delta_doy": 1.5, "triage_category": "Clean_MultiModal_Consensus", "multimodal_conflict_flag": False},
            ]).to_csv(mm_csv, index=False)

            pd.DataFrame([
                {"Pair": "P.aurea vs P.obovata", "Taxon1": "Packera aurea", "Taxon2": "Packera obovata", "N1": 10, "N2": 10, "Warren_D": 0.21, "Warren_I": 0.45, "p_value_D": 0.0099, "p_value_I": 0.0099, "Taxonomic_Inference": "Distinct Ecological Niche"}
            ]).to_csv(niche_csv, index=False)

            proc = subprocess.run(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouchers_csv),
                    "--morphometrics", str(morph_csv),
                    "--vision-audit", str(vision_csv),
                    "--multimodal-flags", str(mm_csv),
                    "--niche-summary", str(niche_csv),
                    "--output-queue", str(out_queue),
                    "--output-triage", str(out_triage),
                    "--output-report", str(out_rep),
                    "--output-plot", str(out_fig),
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Step 07 handoff run failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertTrue(out_triage.exists() and out_rep.exists())
            tdf = pd.read_csv(out_triage)
            self.assertEqual(len(tdf), 10, "Expected 10 rows (preserving ungeoreferenced NCU-HIST-2 and both individuals on CAT-MULTI-SHEET).")

            required_cols = [
                "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
                "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem",
                "capitula_count", "involucre_height_mm", "involucre_width_mm", "AR_inv",
                "is_verified_anthesis", "n_anthesis", "voucher_phenological_state",
                "cda_second_taxon", "cda_second_posterior_prob", "cda_posterior_margin",
                "Sym_PC1", "Sym_PC2", "ecoregion_cda_predicted", "ecoregion_cda_posterior_prob",
                "vision_ood_score", "vision_is_ood", "foreground_attention_ratio",
                "is_artifact_driven", "surface_orientation", "tomentum_taxon_discordance",
                "fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm",
                "ssurgo_rock_outcrop_pct", "edaphic_z_score", "rf_predicted_taxon",
            ]
            for col in required_cols:
                self.assertIn(col, tdf.columns, f"Missing expected column {col} in final_curation_triage.csv")

            # 1. RESOLVED vs CRITICAL priority assignments
            resolved_row = tdf[tdf["catalogNumber"] == "CAT-ANCHOR-A1"].iloc[0]
            critical_row = tdf[tdf["catalogNumber"] == "CAT-SEVERE-MISID"].iloc[0]
            self.assertEqual(str(resolved_row["synthesis_triage_priority"]), "RESOLVED")
            self.assertEqual(str(critical_row["synthesis_triage_priority"]), "CRITICAL")

            # 2. Ungeoreferenced historic specimen coalescing & missing-modality weight renormalization
            hist_row = tdf[tdf["catalogNumber"] == "NCU-HIST-2"].iloc[0]
            self.assertTrue(bool(hist_row["is_completed_voucher"]))
            self.assertEqual(str(hist_row["cda_predicted_taxon"]), "Packera obovata")
            self.assertEqual(str(hist_row["vision_predicted_label"]), "Packera obovata")
            self.assertEqual(str(hist_row["taxonomic_treatment_subtype"]), "Misidentification_Reassign")
            self.assertEqual(int(hist_row["n_available_modalities"]), 2)
            self.assertGreater(float(hist_row["composite_conflict_score"]), 0.80)

            # 3. Background-shortcut vision prediction MUST NOT trigger Stage 1 Misidentification_Reassign
            art_row = tdf[tdf["catalogNumber"] == "CAT-ART-SHORTCUT"].iloc[0]
            self.assertNotEqual(str(art_row["taxonomic_treatment_subtype"]), "Misidentification_Reassign")
            self.assertEqual(str(art_row["recommended_determination"]), "Packera aurea")

            # 4. High vision_ood_score (>=0.65 / vision_is_ood=TRUE) + barrier triggers Stage 3 Putative_Cryptic_Species
            ood_row = tdf[tdf["catalogNumber"] == "CAT-VIS-OOD-CRYPT"].iloc[0]
            self.assertEqual(str(ood_row["taxonomic_treatment_subtype"]), "Putative_Cryptic_Species")
            self.assertIn("vision_ood", str(ood_row["synthesis_rationale"]))

            # 5. Multi-plant herbarium sheet retains distinct individuals and per-individual petiole ratios
            multi_df = tdf[tdf["catalogNumber"] == "CAT-MULTI-SHEET"].sort_values("plant_individual_id")
            self.assertEqual(list(multi_df["plant_individual_id"].astype(int)), [1, 2])
            self.assertAlmostEqual(float(multi_df.iloc[0]["petiole_to_blade_ratio"]), 1.24, places=2)
            self.assertAlmostEqual(float(multi_df.iloc[1]["petiole_to_blade_ratio"]), 0.39, places=2)

    def test_empirical_bayes_hopkins_shrinkage_dynamic_edaphic_z_and_taxon_scaled_pheno_cutoffs(self):
        """Verify Tasks 3A-3C: Empirical Bayes circular Hopkins shrinkage, dynamic primary-predictor edaphic_z_score (no hardcoded 5.8/45.0), multi-PC ANCOVA batch audit, and taxon-dispersion-scaled tau_ecotype / tau_cryptic."""
        r07_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertNotIn("5.8", r07_src, "Must eliminate hardcoded SEUS soil pH 5.8 fallback in Step 07")

        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/analysis/07_triage_dashboard_synthesis.R")

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
        res = subprocess.run(["Rscript", "-e", r_code], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
        self.assertEqual(res.returncode, 0, f"Step 07 cohort-adaptive test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("STEP07_COHORT_ADAPTIVE_OK", res.stdout)


if __name__ == "__main__":
    unittest.main()
