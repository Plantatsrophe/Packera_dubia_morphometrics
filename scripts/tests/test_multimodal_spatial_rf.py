#!/usr/bin/env python3
"""
===============================================================================
Script: test_multimodal_spatial_rf.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 06 (scripts/analysis/06_multimodal_spatial_rf.R):
    - 0-15 cm Trapezoidal Rhizosphere Integration & Non-Blended Coalesce Hierarchy
    - 6-Class USGS SGMC Bedrock Lithology & Azonal Contact Proximity
    - Coordinate-Uncertainty Buffered Focal Extraction (>10000m Training Exclusion)
    - Block-Wise Collinearity Pruning (|r| < 0.70, VIF < 5.0) & Partitioned Niche Tests
    - Multi-Class SGMC/SSURGO Substrate Contingency & Organic PCA-env Ecotypes
    - Authority-Anchored Spatial Block CV, OOF Posterior Pooling, & Step 07 Handshake
===============================================================================
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon


class TestMultimodalSpatialRF(unittest.TestCase):
    """Test suite for Step 06 active R script (scripts/analysis/06_multimodal_spatial_rf.R)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
        cls.util_path = PROJECT_ROOT / "scripts" / "analysis" / "edaphic_fine_grain_utils.R"
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.con_out = cls._class_tmp_path / "contingency.csv"
        cls.pca_env_dir = cls._class_tmp_path / "pca_env_plots"
        cls.pca_env_dir.mkdir(parents=True, exist_ok=True)
        cls.pca_out_pdf = cls.pca_env_dir / "spatial_rf_niche_importance.pdf"
        cls.r_unit_proc = None

        if shutil.which("Rscript"):
            ssurgo_out = cls._class_tmp_path / "ssurgo_val.csv"
            r_batch_code = f"""
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/06_multimodal_spatial_rf.R")

            # 1. Taxonomic standardization
            local({{
              stopifnot(standardize_packera_taxon("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve") == "Packera streptanthifolia")
              stopifnot(standardize_packera_taxon("P. sanguisorbae (DC.) C.Jeffrey") == "Packera sanguisorbae")
              stopifnot(standardize_packera_taxon("Packera x memmingeri") == "Packera memmingeri")
              cat("STD_TAXON_OK\\n")
            }})

            # 2. 0-15 cm trapezoidal rhizosphere integration, log-ratio texture, & strict non-blended coalesce
            local({{
              v_0_5 <- c(4.5, 6.0, NA, 5.2)
              v_5_15 <- c(5.1, 6.6, 5.8, NA)
              rhizo <- compute_rhizosphere_0_15cm(v_0_5, v_5_15)
              stopifnot(abs(rhizo[1] - ((5 * 4.5 + 10 * 5.1) / 15)) < 1e-8)
              stopifnot(abs(rhizo[2] - ((5 * 6.0 + 10 * 6.6) / 15)) < 1e-8)
              stopifnot(abs(rhizo[3] - 5.8) < 1e-8)
              stopifnot(abs(rhizo[4] - 5.2) < 1e-8)

              tex <- compute_soil_texture_log_ratios(sand_pct = c(70.0, 30.0), clay_pct = c(10.0, 35.0))
              stopifnot(all(c("soil_silt_pct", "alr_sand_clay", "alr_silt_clay") %in% names(tex)))
              stopifnot(tex$alr_sand_clay[1] > 0 && tex$alr_sand_clay[2] < 0)

              df_in <- data.frame(
                catalogNumber = c("V1", "V2", "V3"),
                latitude = c(35.8, 36.1, 50.5),
                longitude = c(-79.1, -86.4, -85.0),
                is_rock_outcrop = c(TRUE, FALSE, FALSE),
                outcrop_class = c("Granitic flatrock", "Non-outcrop matrix", "Non-outcrop matrix"),
                ssurgo_rock_outcrop_pct = c(55.0, 0.0, NA),
                ssurgo_bedrock_depth_cm = c(18.0, 180.0, NA),
                ssurgo_sand_pct = c(74.0, 32.0, NA),
                ssurgo_clay_pct = c(8.0, 28.0, NA),
                ssurgo_ph = c(4.7, 6.5, NA),
                ssurgo_cec = c(4.2, 19.5, NA),
                ssurgo_aws_0_50cm = c(2.2, 9.0, NA),
                soil_ph = c(6.9, 4.1, 5.8),
                soil_sand = c(20.0, 85.0, 62.0),
                soil_cec = c(30.0, 3.0, 14.5),
                stringsAsFactors = FALSE
              )
              df_out <- extract_polaris_and_solus_predictors(df_in, env_dir = tempdir(), demo = TRUE)
              req_cols <- c("bdricm", "resdept_r", "cfvo", "awc", "fine_soil_cec", "alr_sand_clay", "alr_silt_clay")
              stopifnot(all(req_cols %in% names(df_out)))
              stopifnot(df_out$bdricm[1] < df_out$bdricm[2])
              stopifnot(df_out$cfvo[1] > df_out$cfvo[2])
              stopifnot(abs(df_out$fine_soil_ph[1] - 4.7) < 1e-6 && abs(df_out$fine_soil_ph[2] - 6.5) < 1e-6 && abs(df_out$fine_soil_ph[3] - 5.8) < 1e-6)
              stopifnot(abs(df_out$fine_soil_sand[1] - 74.0) < 1e-6 && abs(df_out$fine_soil_sand[2] - 32.0) < 1e-6 && abs(df_out$fine_soil_sand[3] - 62.0) < 1e-6)
              stopifnot(abs(df_out$fine_soil_cec[1] - 4.2) < 1e-6 && abs(df_out$fine_soil_cec[2] - 19.5) < 1e-6 && abs(df_out$fine_soil_cec[3] - 14.5) < 1e-6)

              sub_pts <- data.frame(gid = 1L, catalogNumber = "V1", latitude = 35.8, longitude = -79.1, stringsAsFactors = FALSE)
              tab_res <- data.frame(
                gid = 1L, mukey = "101", musym = "WkB", muname = "Wake-Rock outcrop complex",
                brockdepmin = 18.0, aws050wta = 2.4, cokey = "101_1", compname = "Wake",
                comppct_r = 65.0, compkind = "Series", majcompflag = "Yes",
                taxsuborder = "Udipsamments", taxgreatgroup = "Quartzipsamments",
                hzdept_r = 0.0, hzdepb_r = 15.0, sandtotal_r = 76.0, claytotal_r = 8.0,
                ph1to1h2o_r = 4.8, om_r = 1.7, cec7_r = 4.1, awc_r = 0.05,
                resdept_r = 18.0, reskind = "Lithic bedrock", stringsAsFactors = FALSE
              )
              agg1 <- aggregate_ssurgo_tabular_to_points(sub_pts, tab_res)
              stopifnot(nrow(agg1) == 1 && isTRUE(agg1$is_rock_outcrop[1]) && agg1$outcrop_class[1] == "Granitic flatrock")
              cat("RHIZOSPHERE_PEDOLOGY_OK\\n")
            }})

            # 3. 6-class USGS SGMC bedrock lithology & azonal contact proximity
            local({{
              stopifnot(length(SGMC_LITHOLOGY_CLASSES) == 6L)
              expected_classes <- c(
                "Ultramafic_Serpentine",
                "Mafic_Igneous_Metamorphic",
                "Carbonate_Calcareous",
                "Felsic_Granitic_Gneiss",
                "Siliciclastic_Consolidated",
                "Unconsolidated_Sedimentary"
              )
              stopifnot(identical(SGMC_LITHOLOGY_CLASSES, expected_classes))
              mapped <- standardize_sgmc_lithology(
                rocktype1 = c("Serpentinite and peridotite", "Diabase and gabbro", "Limestone and dolostone", "Biotite granite gneiss", "Quartzite and sandstone", "Unconsolidated coastal plain sand"),
                muname = c("Chrome silt loam", "Iredell loam", "Gladeville-Rock outcrop", "Wake-Rock outcrop complex", "Ramsey sandy loam", "Wagram loamy sand")
              )
              stopifnot(identical(as.character(mapped), expected_classes))

              df_geo <- data.frame(
                catalogNumber = paste0("G", 1:6),
                latitude = c(39.7, 35.9, 36.1, 35.8, 35.5, 34.2),
                longitude = c(-76.1, -78.9, -86.4, -80.2, -82.5, -78.1),
                sgmc_rocktype1 = c("Serpentinite", "Diabase dike", "Ordovician limestone", "Granitic gneiss", "Sandstone", "Alluvial sand"),
                is_rock_outcrop = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
                ssurgo_rock_outcrop_pct = c(45, 35, 50, 65, 5, 0),
                fine_bedrock_depth_cm = c(25, 35, 18, 12, 95, 200),
                stringsAsFactors = FALSE
              )
              df_geo <- extract_sgmc_lithology(df_geo, env_dir = tempdir(), demo = TRUE)
              stopifnot(all(c("sgmc_lithology_class", "is_azonal_lithology", "azonal_contact_dist_m", "azonal_proximity_index") %in% names(df_geo)))
              stopifnot(all(df_geo$is_azonal_lithology[1:4]))
              stopifnot(!df_geo$is_azonal_lithology[6])
              stopifnot(df_geo$azonal_contact_dist_m[4] < df_geo$azonal_contact_dist_m[6])
              stopifnot(df_geo$azonal_proximity_index[4] > df_geo$azonal_proximity_index[6])
              cat("SGMC_LITHOLOGY_OK\\n")
            }})

            # 4. Coordinate-uncertainty buffered extraction & >10000m training exclusion
            local({{
              radii <- compute_uncertainty_buffer_radius(c(NA, 25, 250, 600, 3500, 25000))
              stopifnot(identical(as.numeric(radii), c(250, 100, 250, 600, 1000, 1000)))

              df <- data.frame(
                catalogNumber = c("U1", "U2", "U3", "U4"),
                latitude = c(35.8, 35.9, 36.0, 34.5),
                longitude = c(-80.1, -80.2, -79.8, -78.5),
                coordinateUncertaintyInMeters = c(NA, 45, 750, 15000),
                fine_soil_ph = c(4.9, 5.1, 6.2, 5.4),
                fine_bedrock_depth_cm = c(22, 35, 110, 160),
                cfvo = c(32, 24, 12, 6),
                ssurgo_rock_outcrop_pct = c(55, 30, 5, 0),
                is_rock_outcrop = c(TRUE, TRUE, FALSE, FALSE),
                stringsAsFactors = FALSE
              )
              df_buf <- extract_uncertainty_buffered_predictors(df, env_dir = tempdir(), demo = TRUE)
              req <- c("buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
                       "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer")
              stopifnot(all(req %in% names(df_buf)))
              stopifnot(identical(df_buf$high_coord_uncertainty_flag, c(FALSE, FALSE, FALSE, TRUE)))
              stopifnot(identical(df_buf$is_training_eligible_coord, c(TRUE, TRUE, TRUE, FALSE)))
              stopifnot(all(df_buf$ph_q10_buffer <= df_buf$fine_soil_ph))
              stopifnot(all(df_buf$ph_q90_buffer >= df_buf$fine_soil_ph))
              stopifnot(all(df_buf$min_bedrock_depth_buffer_cm <= df_buf$fine_bedrock_depth_cm))
              stopifnot(all(df_buf$max_cfvo_buffer_pct >= df_buf$cfvo))
              stopifnot(df_buf$azonal_lithology_presence_buffer[1] == TRUE)
              cat("UNCERTAINTY_BUFFER_OK\\n")
            }})

            # 5. Block-wise collinearity pruning & partitioned D_clim vs D_edaph PCA-env niche tests
            local({{
              set.seed(42)
              n <- 40
              sp <- rep(c("Packera anonyma", "Packera dubia"), each = n / 2)
              is_dub <- sp == "Packera dubia"
              df <- data.frame(
                catalogNumber = sprintf("V_%03d", seq_len(n)),
                species_standardized = sp,
                analysis_group = sp,
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"), length.out = n),
                latitude = ifelse(is_dub, runif(n, 33.5, 36.0), runif(n, 33.8, 36.2)),
                longitude = ifelse(is_dub, runif(n, -81.0, -78.0), runif(n, -81.2, -78.2)),
                coordinateUncertaintyInMeters = c(rep(200, n - 2), 18000, 25000),
                bio1_temp_mean = rnorm(n, 16.5, 0.8),
                bio4_temp_seasonality = rnorm(n, 620, 25),
                bio12_precip_annual = rnorm(n, 1250, 60),
                bio15_precip_seasonality = rnorm(n, 18, 2.0),
                delta_doy = ifelse(is_dub, rnorm(n, -10, 4), rnorm(n, 8, 4)),
                fine_soil_ph = ifelse(is_dub, rnorm(n, 4.6, 0.25), rnorm(n, 5.9, 0.30)),
                fine_soil_sand = ifelse(is_dub, rnorm(n, 78, 5), rnorm(n, 42, 6)),
                fine_bedrock_depth_cm = ifelse(is_dub, rnorm(n, 22, 6), rnorm(n, 135, 20)),
                cfvo = ifelse(is_dub, rnorm(n, 38, 6), rnorm(n, 9, 3)),
                awc = ifelse(is_dub, rnorm(n, 0.06, 0.015), rnorm(n, 0.16, 0.02)),
                alr_sand_clay = ifelse(is_dub, rnorm(n, 2.1, 0.25), rnorm(n, 0.4, 0.25)),
                azonal_proximity_index = ifelse(is_dub, rnorm(n, 0.82, 0.08), rnorm(n, 0.15, 0.06)),
                stringsAsFactors = FALSE
              )
              df$bio1_collinear_dup <- df$bio1_temp_mean * 1.01 + rnorm(n, 0, 0.005)
              env_vars <- c("bio1_temp_mean", "bio1_collinear_dup", "bio4_temp_seasonality", "bio12_precip_annual",
                            "bio15_precip_seasonality", "delta_doy", "fine_soil_ph", "fine_soil_sand",
                            "fine_bedrock_depth_cm", "cfvo", "awc", "alr_sand_clay", "azonal_proximity_index")
              bs <- screen_collinear_predictors_by_block(df[, env_vars], cor.threshold = 0.70, vif.threshold = 5.0)
              stopifnot(length(bs$climatic_variables) >= 1)
              stopifnot(length(bs$edaphic_variables) >= 1)
              stopifnot(!(all(c("bio1_temp_mean", "bio1_collinear_dup") %in% bs$climatic_variables)))

              niche_out <- run_warren_niche_identity_tests(df, env_vars, n_perm = 15)
              sum_df <- niche_out$summary
              req_cols <- c("Pair", "Taxon1", "Taxon2", "N1", "N2",
                            "D_clim", "I_clim", "p_equiv_D_clim", "p_bg_D_clim",
                            "D_edaph", "I_edaph", "p_equiv_D_edaph", "p_bg_D_edaph",
                            "D_joint", "I_joint", "Schoener_D", "Warren_D", "Warren_I",
                            "p_value_D", "p_value_I", "p_bg_D_joint", "Niche_Divergence_Axis", "Taxonomic_Inference")
              stopifnot(all(req_cols %in% names(sum_df)))
              stopifnot((sum_df$N1[1] + sum_df$N2[1]) == (n - 2))
              stopifnot(sum_df$D_edaph[1] < sum_df$D_clim[1])
              cat("PARTITIONED_NICHE_OK\\n")
            }})

            # 6. Multi-class lithology & substrate contingency on non-dubia species
            local({{
              df <- data.frame(
                catalogNumber = sprintf("ND_%03d", 1:45),
                latitude = c(rep(36.15, 15), rep(35.80, 15), rep(42.20, 15)),
                longitude = c(rep(-81.65, 15), rep(-79.10, 15), rep(-72.80, 15)),
                species_standardized = c(
                  rep("Packera serpenticola", 15),
                  rep("Packera aurea", 15),
                  rep("Packera schweinitziana", 15)
                ),
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver", "Tier_1_Gold"), 15),
                coordinateUncertaintyInMeters = rep(150.0, 45),
                is_rock_outcrop = c(rep(TRUE, 14), FALSE, rep(FALSE, 15), rep(FALSE, 14), TRUE),
                outcrop_class = c(
                  rep("Serpentine", 14), "Non-outcrop matrix",
                  rep("Non-outcrop matrix", 15),
                  rep("Non-outcrop matrix", 14), "Sandstone/Siliceous"
                ),
                taxsuborder = c(
                  rep("Udepts", 15),
                  rep("Aquults", 14), "Udults",
                  rep("Orthents", 15)
                ),
                sgmc_lithology_class = c(
                  rep("Ultramafic_Serpentine", 14), "Mafic_Igneous_Metamorphic",
                  rep("Unconsolidated_Sedimentary", 15),
                  rep("Siliciclastic_Consolidated", 14), "Carbonate_Calcareous"
                ),
                is_azonal_lithology = c(rep(TRUE, 15), rep(FALSE, 15), rep(FALSE, 14), TRUE),
                ssurgo_rock_outcrop_pct = c(rep(55.0, 14), 5.0, rep(0.0, 15), rep(2.0, 14), 22.0),
                fine_bedrock_depth_cm = c(rep(24.0, 14), 85.0, rep(180.0, 15), rep(135.0, 14), 45.0),
                fine_soil_sand = c(rep(42.0, 15), rep(28.0, 15), rep(51.0, 15)),
                fine_soil_ph = c(rep(6.8, 15), rep(5.4, 15), rep(4.9, 15)),
                ssurgo_aws_0_50cm = c(rep(3.1, 15), rep(9.4, 15), rep(6.8, 15)),
                topo_tpi_30m = c(rep(3.6, 15), rep(-1.8, 15), rep(0.4, 15)),
                topo_twi_30m = c(rep(4.1, 15), rep(9.8, 15), rep(6.2, 15)),
                doy = rep(132, 45),
                stringsAsFactors = FALSE
              )
              res <- synthesize_multiscale_edaphics(
                df,
                ssurgo_out = "{ssurgo_out.as_posix()}",
                contingency_out = "{cls.con_out.as_posix()}",
                demo = FALSE
              )
              con_df <- res$contingency
              mc_df <- res$multiclass_contingency
              stopifnot(nrow(con_df) == 4)
              stopifnot(con_df$Taxon[1] == "Omnibus_Substrate_Contingency_Across_Active_Groups")
              stopifnot(!any(grepl("dubia|anonyma|paupercula|plattensis", con_df$Taxon, ignore.case = TRUE)))

              req_cols <- c("Taxon", "Outcrop_Count", "Matrix_Count", "Total_Vouchers",
                            "Outcrop_Affinity_Pct", "Odds_Ratio", "CI_Lower", "CI_Upper",
                            "P_Value", "P_Value_BH", "Dominant_Lithology_Class",
                            "Dominant_Substrate_Class", "Top_Enriched_Lithology", "Top_Substrate_Odds_Ratio")
              stopifnot(all(req_cols %in% names(con_df)))

              serp_row <- con_df[con_df$Taxon == "Packera serpenticola", ]
              stopifnot(serp_row$Dominant_Lithology_Class == "Ultramafic_Serpentine")
              stopifnot(as.numeric(serp_row$Top_Substrate_Odds_Ratio) > 10.0)
              stopifnot(as.numeric(serp_row$P_Value_BH) < 0.01)

              aurea_row <- con_df[con_df$Taxon == "Packera aurea", ]
              stopifnot(aurea_row$Dominant_Lithology_Class == "Unconsolidated_Sedimentary")
              stopifnot(aurea_row$Dominant_Substrate_Class %in% c("Unconsolidated_Sedimentary", "Aquults"))

              stopifnot(all(c("SGMC_Lithology", "SSURGO_Outcrop", "SSURGO_Drainage_Suborder") %in% unique(mc_df$Substrate_Domain)))
              stopifnot("P_Value_BH" %in% names(mc_df))

              cons_df <- execute_crossmodal_consensus(df, intraspecific_mode = FALSE)
              stopifnot("substrate_specialist_flag" %in% names(cons_df))
              stopifnot("piedmont_outcrop_specialist_flag" %in% names(cons_df))
              stopifnot(sum(cons_df$substrate_specialist_flag[cons_df$species_standardized == "Packera serpenticola"]) >= 14)
              stopifnot(sum(cons_df$substrate_specialist_flag[cons_df$species_standardized == "Packera aurea"]) == 0)
              cat("NON_DUBIA_CONTINGENCY_OK\\n")
            }})

            # 7. Organic PCA-env intraspecific ecotypes & dynamic Spatial RF plots
            local({{
              set.seed(123)
              n <- 28
              df <- data.frame(
                catalogNumber = sprintf("INTRA_%03d", seq_len(n)),
                species_standardized = rep("Packera aurea", n),
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver"), length.out = n),
                latitude = c(seq(34.0, 35.4, length.out = 14), seq(40.0, 41.8, length.out = 14)),
                longitude = c(seq(-83.0, -81.6, length.out = 14), seq(-76.0, -74.2, length.out = 14)),
                coordinateUncertaintyInMeters = rep(200.0, n),
                fine_bedrock_depth_cm = c(rnorm(14, mean = 25, sd = 3), rnorm(14, mean = 175, sd = 5)),
                bio1_temp_mean = c(rnorm(14, mean = 16.5, sd = 0.4), rnorm(14, mean = 8.2, sd = 0.4)),
                fine_soil_ph = rnorm(n, mean = 5.8, sd = 0.05),
                fine_soil_sand = rnorm(n, mean = 45.0, sd = 0.5),
                topo_twi_30m = c(rnorm(14, mean = 4.5, sd = 0.3), rnorm(14, mean = 9.5, sd = 0.3)),
                doy = rep(130, n),
                delta_doy = rnorm(n, mean = 0, sd = 2),
                stringsAsFactors = FALSE
              )

              resolved <- resolve_analysis_group(df, intraspecific_mode = TRUE)
              u_grps <- sort(unique(resolved$analysis_group))
              stopifnot(length(u_grps) == 2)
              stopifnot(!any(grepl("Acidic_Sand_Coastal|Upland_Piedmont_Outcrop", u_grps)))
              stopifnot(all(grepl("^Ecotype_[12]_[A-Za-z0-9_]+_(Low|High)$", u_grps)))

              df_flat <- data.frame(
                catalogNumber = sprintf("FLAT_%02d", 1:8),
                species_standardized = rep("Packera aurea", 8),
                fine_soil_ph = rep(6.0, 8),
                stringsAsFactors = FALSE
              )
              res_flat <- resolve_analysis_group(df_flat, intraspecific_mode = TRUE)
              stopifnot(identical(sort(unique(res_flat$analysis_group)), c("Ecotype_Cluster_1", "Ecotype_Cluster_2")))

              env_vars <- c("fine_bedrock_depth_cm", "bio1_temp_mean", "fine_soil_ph", "fine_soil_sand", "topo_twi_30m")
              srf_res <- compute_spatial_rf_mems(resolved, env_vars, n_trees = 50, intraspecific_mode = TRUE)
              top_sel <- select_top_rf_predictors(srf_res$importance, resolved)
              stopifnot(top_sel$top_env_1 %in% c("fine_bedrock_depth_cm", "bio1_temp_mean", "topo_twi_30m"))
              stopifnot(top_sel$top_edaph_1 %in% c("fine_bedrock_depth_cm", "topo_twi_30m"))

              niche_res <- run_warren_niche_identity_tests(resolved, env_vars, n_perm = 10, intraspecific_mode = TRUE)
              export_spatial_rf_figures(resolved, srf_res, niche_res, "{cls.pca_out_pdf.as_posix()}")
              cat("ORGANIC_PCA_ENV_AND_PLOTS_OK\\n")
            }})

            # 8. Composite key (catalogNumber, plant_individual_id) multi-plant join & expanded morph/vision columns
            local({{
              vouchers_df <- data.frame(
                catalogNumber = c("NCU_MULTI", "NCU_SINGLE"),
                species_standardized = c("Packera aurea", "Packera obovata"),
                latitude = c(35.6, 36.1),
                longitude = c(-82.5, -79.1),
                stringsAsFactors = FALSE
              )
              morph_df <- data.frame(
                catalogNumber = c("NCU_MULTI", "NCU_MULTI", "NCU_SINGLE"),
                plant_individual_id = c(1L, 2L, 1L),
                morph_tier_source = rep("Tier_1_Pure", 3),
                n_tier1_leaves = c(3L, 2L, 4L),
                gmm_cluster = c(1L, 2L, 1L),
                gmm_uncertainty = c(0.04, 0.18, 0.03),
                cda_predicted_taxon = c("Packera aurea", "Packera obovata", "Packera obovata"),
                cda_posterior_prob = c(0.94, 0.88, 0.96),
                cda_second_taxon = c("Packera obovata", "Packera aurea", "Packera aurea"),
                cda_second_posterior_prob = c(0.06, 0.12, 0.04),
                cda_posterior_margin = c(0.88, 0.76, 0.92),
                Sym_PC1 = c(-1.4, 1.3, 1.5),
                Sym_PC2 = c(0.2, -0.3, -0.1),
                petiole_length_mm = c(28.0, 11.5, 10.0),
                blade_length_mm = c(25.0, 29.0, 30.0),
                petiole_to_blade_ratio = c(1.12, 0.3966, 0.3333),
                plant_height_cm = c(42.0, 35.0, 38.0),
                scape_to_basal_leaf_ratio = c(7.9, 8.6, 9.5),
                is_bent_stem = c(FALSE, TRUE, FALSE),
                capitula_count = c(12L, 8L, 10L),
                involucre_height_mm = c(6.5, 5.8, 6.0),
                involucre_width_mm = c(5.2, 4.9, 5.0),
                AR_inv = c(1.25, 1.18, 1.20),
                ecoregion_cda_predicted = c("Blue_Ridge", "Piedmont", "Piedmont"),
                ecoregion_cda_posterior_prob = c(0.84, 0.79, 0.91),
                stringsAsFactors = FALSE
              )
              vis_df <- data.frame(
                catalogNumber = c("NCU_MULTI", "NCU_SINGLE"),
                plant_individual_id = c(0L, 0L),
                vision_predicted_label = c("Packera aurea", "Packera obovata"),
                confidence_predicted_class = c(0.91, 0.95),
                vision_second_prob = c(0.09, 0.05),
                vision_margin = c(0.82, 0.90),
                c_error = c(0.08, 0.04),
                is_label_corrupted = c(FALSE, FALSE),
                vision_ood_score = c(0.21, 0.72),
                vision_is_ood = c(FALSE, TRUE),
                foreground_attention_ratio = c(0.83, 0.51),
                is_artifact_driven = c(FALSE, TRUE),
                surface_orientation = c("adaxial", "abaxial"),
                mean_luminance = c(54.2, 71.8),
                mean_saturation = c(0.32, 0.14),
                tomentum_taxon_discordance = c(FALSE, TRUE),
                stringsAsFactors = FALSE
              )

              m_joined <- left_join_by_voucher_individual(vouchers_df, morph_df, names(morph_df))
              mv_joined <- left_join_by_voucher_individual(m_joined, vis_df, names(vis_df))
              stopifnot(nrow(mv_joined) == 3L)
              multi_sub <- mv_joined[mv_joined$catalogNumber == "NCU_MULTI", ]
              stopifnot(identical(sort(as.integer(multi_sub$plant_individual_id)), c(1L, 2L)))
              stopifnot(all(c("cda_second_taxon", "Sym_PC1", "petiole_to_blade_ratio", "plant_height_cm",
                              "AR_inv", "ecoregion_cda_predicted", "vision_ood_score", "vision_is_ood",
                              "foreground_attention_ratio", "is_artifact_driven", "surface_orientation",
                              "tomentum_taxon_discordance") %in% names(mv_joined)))
              stopifnot(all(multi_sub$foreground_attention_ratio == 0.83))
              cat("STEP06_COMPOSITE_HANDOFF_OK\\n")
            }})

            # 9. Spatial Block-CV OOF RF posterior integration in edaphic_best_fit_taxon
            local({{
              set.seed(106)
              n6 <- 24
              sp6 <- rep(c("Packera aurea", "Packera obovata"), each = n6 / 2)
              df6 <- data.frame(
                catalogNumber = sprintf("S6_%02d", seq_len(n6)),
                plant_individual_id = 1L,
                species_standardized = sp6,
                analysis_group = sp6,
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver"), length.out = n6),
                latitude = c(seq(35.0, 36.1, length.out = 12), seq(36.2, 37.3, length.out = 12)),
                longitude = c(seq(-82.5, -81.4, length.out = 12), seq(-80.0, -78.9, length.out = 12)),
                coordinateUncertaintyInMeters = rep(150.0, n6),
                fine_soil_ph = c(rnorm(12, 5.2, 0.15), rnorm(12, 7.1, 0.15)),
                fine_soil_sand = c(rnorm(12, 65.0, 3.0), rnorm(12, 28.0, 3.0)),
                fine_bedrock_depth_cm = c(rnorm(12, 120.0, 8.0), rnorm(12, 35.0, 6.0)),
                ssurgo_rock_outcrop_pct = c(rep(2.0, 12), rep(28.0, 12)),
                topo_twi_30m = c(rnorm(12, 9.2, 0.3), rnorm(12, 4.8, 0.3)),
                doy = rep(120, n6),
                stringsAsFactors = FALSE
              )
              df6$fine_soil_ph[1] <- 7.2
              df6$fine_soil_sand[1] <- 26.0
              df6$fine_bedrock_depth_cm[1] <- 33.0
              df6$topo_twi_30m[1] <- 4.7
              df6$determiner_tier[1] <- "Tier_3_Bronze"

              srf6 <- compute_spatial_rf_mems(
                df6,
                c("fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "topo_twi_30m"),
                n_trees = 50
              )
              stopifnot(!is.null(srf6$oof_prob_matrix))
              stopifnot(nrow(srf6$oof_prob_matrix) == nrow(df6))

              cons_mahal_only <- execute_crossmodal_consensus(df6, intraspecific_mode = FALSE, rf_res = NULL)
              cons_pooled <- execute_crossmodal_consensus(df6, intraspecific_mode = FALSE, rf_res = srf6)
              stopifnot("edaphic_best_fit_taxon" %in% names(cons_pooled))
              stopifnot("edaphic_z_score" %in% names(cons_pooled))
              stopifnot(cons_mahal_only$edaphic_best_fit_taxon[1] == "Packera obovata")
              stopifnot(cons_pooled$edaphic_best_fit_taxon[1] == "Packera obovata")

              df_mid <- df6
              df_mid$fine_soil_ph[2] <- mean(df6$fine_soil_ph)
              df_mid$fine_soil_sand[2] <- mean(df6$fine_soil_sand)
              df_mid$fine_bedrock_depth_cm[2] <- mean(df6$fine_bedrock_depth_cm)
              df_mid$ssurgo_rock_outcrop_pct[2] <- mean(df6$ssurgo_rock_outcrop_pct)
              df_mid$topo_twi_30m[2] <- mean(df6$topo_twi_30m)
              synth_rf <- list(
                oof_prob_matrix = matrix(
                  c(rep(c(0.90, 0.10), 1), c(0.03, 0.97), rep(c(0.90, 0.10), 10), rep(c(0.10, 0.90), 12)),
                  nrow = n6,
                  ncol = 2,
                  byrow = TRUE,
                  dimnames = list(NULL, c("Packera aurea", "Packera obovata"))
                ),
                spatial_block_cv_acc = 0.92
              )
              cons_rf_tiebreak <- execute_crossmodal_consensus(df_mid, intraspecific_mode = FALSE, rf_res = synth_rf)
              stopifnot(cons_rf_tiebreak$edaphic_best_fit_taxon[2] == "Packera obovata")
              cat("STEP06_RF_POSTERIOR_POOLING_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step06_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = subprocess.run(
                ["Rscript", str(r_batch_file)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"
        self.util_path = PROJECT_ROOT / "scripts" / "analysis" / "edaphic_fine_grain_utils.R"

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 06 R unit suite failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_r_script_phenological_and_fine_grain_edaphic_contracts(self):
        """Verify 06_multimodal_spatial_rf.R and edaphic_fine_grain_utils.R retain all phenological, ENM, SSURGO, POLARIS, 3DEP, and Spatial RF contracts."""
        self.assertTrue(self.r_script_path.exists(), "06_multimodal_spatial_rf.R must exist.")
        self.assertTrue(self.util_path.exists(), "edaphic_fine_grain_utils.R must exist.")

        r_content = self.r_script_path.read_text(encoding="utf-8")
        for token in (
            "--pheno-anomalies",
            "--pheno-summary",
            "--pheno-plot",
            "--demo",
            "model_latitudinal_spring_baseline",
            "compute_phenological_anomalies",
            "export_phenological_latitudinal_figures",
            "delta_doy",
            'dplyr::left_join(vouchers_df, pheno_states_df, by = "catalogNumber")',
            "vegan::varpart",
            "geodata::worldclim_global",
            "geodata::soil_world",
            "adespatial::dbmem",
            "ranger::ranger",
            "blockCV::cv_spatial",
            "Spatial_Block_CV_Accuracy",
            "ecospat::ecospat.grid.clim.dyn",
            "ecospat::ecospat.niche.overlap",
            "ecospat::ecospat.niche.equivalency.test",
        ):
            self.assertIn(token, r_content)
        self.assertNotIn("* 100 - 10", r_content)
        self.assertNotIn("* 100 - 5", r_content)

        combined_r = r_content + "\n" + self.util_path.read_text(encoding="utf-8")
        for token in (
            "ST_Buffer(ST_Point(",
            "0.00045",
            "ssurgo_rock_outcrop_pct",
            "ssurgo_bedrock_depth_cm",
            "ssurgo_sand_pct",
            "ssurgo_clay_pct",
            "ssurgo_ph",
            "ssurgo_om_pct",
            "ssurgo_cec",
            "ssurgo_aws_0_50cm",
            "extract_polaris_and_solus_predictors",
            "polaris_ph_30m",
            "polaris_sand_30m",
            "polaris_clay_30m",
            "polaris_om_30m",
            "polaris_bd_30m",
            "polaris_ksat_30m",
            "fine_soil_ph",
            "fine_soil_sand",
            "fine_soil_clay",
            "fine_bedrock_depth_cm",
            "fine_rock_outcrop_prob",
            "topo_tpi_30m",
            "topo_roughness_30m",
            "topo_slope_deg",
            "topo_twi_30m",
            "fine_grain_edaphic_cache.csv",
            "[A] Fine-Grain Edaphic & Outcrop Lithology",
            "[B] Macro-Climate",
            "[C] Spatial MEMs",
            "edaphic_z_score",
            "fill = NA",
            "spatialRF::auto_cor",
            "spatialRF::auto_vif",
            "screen_collinear_predictors",
            "Schoener_D",
            "Warren_I",
            "residual_moran_i",
            "hopkins_adjusted_doy",
        ):
            self.assertIn(token, combined_r, f"Expected fine-grain token '{token}' in Step 06 R code.")

    def test_taxonomic_standardization(self):
        """Verify generic species-agnostic standardization across North American Packera binomials in Python and R."""
        self.assertEqual(standardize_packera_taxon("Packera dubia (Spreng.) Trock & Mabb."), "Packera dubia")
        self.assertEqual(standardize_packera_taxon("Packera anonyma (Alph.Wood) W.A.Weber"), "Packera anonyma")
        self.assertEqual(standardize_packera_taxon("Packera paupercula var. savannarum"), "Packera paupercula")
        self.assertEqual(standardize_packera_taxon("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve"), "Packera streptanthifolia")
        self.assertEqual(standardize_packera_taxon("Packera sanguisorbae (DC.) C.Jeffrey"), "Packera sanguisorbae")
        self.assertEqual(standardize_packera_taxon("Packera x memmingeri"), "Packera memmingeri")
        self.assertIn("STD_TAXON_OK", self._require_batched_r())

    def test_rhizosphere_0_15cm_sgmc_uncertainty_and_collinearity_niche_tests(self):
        """Verify 0-15 cm trapezoidal rhizosphere pedology, 6-class SGMC lithology, uncertainty buffering, and partitioned climatic vs edaphic niche tests."""
        out = self._require_batched_r()
        for marker in (
            "RHIZOSPHERE_PEDOLOGY_OK",
            "SGMC_LITHOLOGY_OK",
            "UNCERTAINTY_BUFFER_OK",
            "PARTITIONED_NICHE_OK",
        ):
            self.assertIn(marker, out)

    def test_multiclass_contingency_organic_pca_env_composite_join_and_rf_pooling(self):
        """Verify multi-class SGMC/SSURGO contingency, organic PCA-env ecotypes, multi-plant composite join, and Spatial Block-CV OOF RF posterior pooling."""
        out = self._require_batched_r()
        for marker in (
            "NON_DUBIA_CONTINGENCY_OK",
            "ORGANIC_PCA_ENV_AND_PLOTS_OK",
            "STEP06_COMPOSITE_HANDOFF_OK",
            "STEP06_RF_POSTERIOR_POOLING_OK",
        ):
            self.assertIn(marker, out)

        con_csv_df = pd.read_csv(self.con_out)
        self.assertEqual(con_csv_df.loc[0, "Taxon"], "Omnibus_Substrate_Contingency_Across_Active_Groups")
        self.assertIn("Dominant_Lithology_Class", con_csv_df.columns)
        self.assertIn("P_Value_BH", con_csv_df.columns)

        self.assertTrue(self.pca_out_pdf.exists(), "spatial_rf_niche_importance.pdf must be generated.")
        self.assertTrue((self.pca_env_dir / "fig06_ecological_niche_phenology.pdf").exists())
        self.assertTrue((self.pca_env_dir / "fig06_ecological_niche_phenology.png").exists())

    def test_authority_anchored_spatial_block_cv_zero_synthetic_cache_and_step07_contract(self):
        """Verify authority-anchored Spatial RF training, zero synthetic disk caching under --demo, passive Bronze/high-uncertainty projection, and Step 07 compatibility."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vouchers_csv = tmp / "curated_vouchers.csv"
            morph_csv = tmp / "morphometric_flags.csv"
            vision_csv = tmp / "label_noise_audit.csv"
            pheno_states_csv = tmp / "voucher_phenological_states.csv"
            out_flags = tmp / "multimodal_conflict_flags.csv"
            out_summary = tmp / "multimodal_spatial_rf_summary.csv"
            out_ssurgo = tmp / "ssurgo_edaphic_validation.csv"
            out_contingency = tmp / "micro_edaphic_outcrop_contingency.csv"
            out_pheno_anom = tmp / "phenological_anomalies.csv"
            out_pheno_sum = tmp / "phenological_anomaly_summary.csv"
            out_plot = tmp / "spatial_rf_niche_summary.pdf"
            out_pheno_plot = tmp / "phenological_latitudinal_anomaly.pdf"

            fixture_taxa = ["Packera aurea", "Packera millefolium", "Packera obovata", "Packera schweinitziana"]
            rows = []
            for idx in range(32):
                sp = fixture_taxa[idx % len(fixture_taxa)]
                lat = 33.0 + (idx % 8) * 1.05
                lon = -84.0 + (idx % 6) * 0.85
                doy = int(round(-30.0 + 3.9 * lat))
                tier = "Tier_1_Gold" if idx % 4 in (0, 1) else ("Tier_2_Silver" if idx % 4 == 2 else "Tier_3_Bronze")
                coord_unc = 16000.0 if idx == 31 else float(120 + (idx % 5) * 150)
                rows.append({
                    "catalogNumber": f"ANC_{idx:03d}",
                    "institutionCode": "NCU",
                    "scientificName": sp,
                    "species_raw": sp,
                    "determiner_tier": tier,
                    "latitude": lat,
                    "longitude": lon,
                    "decimalLatitude": lat,
                    "decimalLongitude": lon,
                    "coordinateUncertaintyInMeters": coord_unc,
                    "eventDate": f"2022-04-{(idx % 20) + 1:02d}",
                    "year": 2022,
                    "month": 4,
                    "day": (idx % 20) + 1,
                    "doy": doy,
                    "reproductiveCondition": "fl",
                    "regional_group": "Piedmont_Granite_Flatrocks" if sp == "Packera anonyma" else "Southeastern_Coastal_Plain",
                    "is_primary_duplicate": True,
                })
            pd.DataFrame(rows).to_csv(vouchers_csv, index=False)

            pd.DataFrame([
                {
                    "catalogNumber": r["catalogNumber"],
                    "plant_individual_id": 1,
                    "morph_tier_source": "Tier_1_Pure",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "cda_predicted_taxon": r["species_raw"],
                    "cda_posterior_prob": 0.92,
                    "can1": 1.1,
                    "can2": -0.3,
                    "can3": 0.2,
                    "gmm_cluster": (i % 4) + 1,
                    "gmm_uncertainty": 0.05,
                    "morphometric_outlier_flag": False,
                    "misidentification_flag": False,
                }
                for i, r in enumerate(rows)
            ]).to_csv(morph_csv, index=False)

            pd.DataFrame([
                {
                    "catalogNumber": r["catalogNumber"],
                    "predicted_label": r["species_raw"],
                    "confidence_predicted_class": 0.91,
                    "c_error": 0.08,
                    "is_label_corrupted": False,
                }
                for r in rows
            ]).to_csv(vision_csv, index=False)

            pd.DataFrame([
                {
                    "catalogNumber": r["catalogNumber"],
                    "total_capitula": 4,
                    "n_anthesis": 3,
                    "n_fruit": 1,
                    "n_bud": 0,
                    "voucher_phenological_state": "anthesis",
                    "dominant_pappus_ratio": 0.10,
                }
                for r in rows
            ]).to_csv(pheno_states_csv, index=False)

            cmd = [
                "Rscript",
                str(self.r_script_path),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(morph_csv),
                "--vision-audit", str(vision_csv),
                "--pheno-states", str(pheno_states_csv),
                "--env-dir", str(tmp / "env"),
                "--ssurgo-cache", str(tmp / "ssurgo_cache"),
                "--ssurgo-table", str(out_ssurgo),
                "--contingency-table", str(out_contingency),
                "--pheno-anomalies", str(out_pheno_anom),
                "--pheno-summary", str(out_pheno_sum),
                "--pheno-plot", str(out_pheno_plot),
                "--output-flags", str(out_flags),
                "--output-summary", str(out_summary),
                "--output-plot", str(out_plot),
                "--permutations", "10",
                "--n-trees", "50",
                "--demo",
            ]
            res = subprocess.run(cmd, cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(res.returncode, 0, f"Step 06 Rscript failed: {res.stderr}\n{res.stdout}")

            # Verify synthetic SSURGO records are never cached into ssurgo_cache under --demo
            cache_dir = tmp / "ssurgo_cache"
            if cache_dir.exists():
                cached_files = list(cache_dir.glob("*.rds")) + list(cache_dir.glob("*.csv"))
                self.assertEqual(cached_files, [], "Synthetic SSURGO records must never be cached via saveRDS/write_csv.")

            flags_df = pd.read_csv(out_flags)
            self.assertEqual(len(flags_df), 32)
            for col in [
                "multimodal_conflict_flag",
                "triage_category",
                "edaphic_best_fit_taxon",
                "rf_predicted_taxon",
                "rf_predicted_prob",
                "rf_given_taxon_prob",
                "rf_prediction_source",
                "spatial_block_cv_acc",
                "sgmc_lithology_class",
                "bdricm",
                "cfvo",
                "awc",
                "alr_sand_clay",
                "ph_q10_buffer",
                "ph_q90_buffer",
                "min_bedrock_depth_buffer_cm",
                "max_cfvo_buffer_pct",
                "azonal_lithology_presence_buffer",
            ]:
                self.assertIn(col, flags_df.columns, f"Missing column '{col}' in multimodal_conflict_flags.csv")

            gold_silver_rows = flags_df[
                flags_df["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])
                & (~flags_df["high_coord_uncertainty_flag"])
            ]
            bronze_rows = flags_df[
                (flags_df["determiner_tier"] == "Tier_3_Bronze")
                & (~flags_df["high_coord_uncertainty_flag"])
            ]
            high_unc_rows = flags_df[flags_df["catalogNumber"] == "ANC_031"]
            self.assertTrue((gold_silver_rows["rf_prediction_source"] == "Spatial_Block_OOF_Anchor").all())
            self.assertTrue((bronze_rows["rf_prediction_source"] == "Passive_Projection_Bronze").all())
            self.assertTrue((high_unc_rows["rf_prediction_source"] == "Passive_Projection_High_Uncertainty").all())

            sum_df = pd.read_csv(out_summary)
            for col in ["D_clim", "I_clim", "D_edaph", "I_edaph", "D_joint", "I_joint", "Schoener_D", "Warren_D", "Warren_I"]:
                self.assertIn(col, sum_df.columns, f"Missing column '{col}' in multimodal_spatial_rf_summary.csv")

            # Verify Step 07 07_triage_dashboard_synthesis.R accepts the Step 06 outputs and completes without error
            step07_script = PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"
            out_queue = tmp / "triage_queue.csv"
            out_report = tmp / "revision_summary.md"
            out_triage_plot = tmp / "triage_plot.pdf"
            cmd07 = [
                "Rscript",
                str(step07_script),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(morph_csv),
                "--vision-audit", str(vision_csv),
                "--multimodal-flags", str(out_flags),
                "--niche-summary", str(out_summary),
                "--output-queue", str(out_queue),
                "--output-report", str(out_report),
                "--output-plot", str(out_triage_plot),
            ]
            res07 = subprocess.run(cmd07, cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(res07.returncode, 0, f"Step 07 execution failed: {res07.stderr}\n{res07.stdout}")
            self.assertTrue(out_queue.exists(), "Step 07 triage_queue.csv must be generated from Step 06 outputs.")
            queue_df = pd.read_csv(out_queue)
            self.assertEqual(len(queue_df), 32)

    def test_cohort_adaptive_19_bioclim_rf_primary_gating_arctic_phenology_and_coord_free_lithology(self):
        """Verify Tasks 2A-2E: all 19 WorldClim bio1..bio19, anchor-ranked 3-block VIF + RF primary predictor gating, RF-weighted edaphic/climatic z-scores, Arctic late-summer circular phenology, and coordinate-free non-demo lithology/outcrop classification."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        r_code = """
        suppressPackageStartupMessages({ library(readr); library(dplyr) })
        source("scripts/core/db_utils.R")
        source("scripts/morphometrics/authority_weighting_utils.R")
        source("scripts/analysis/edaphic_fine_grain_utils.R")
        source("scripts/analysis/06_multimodal_spatial_rf.R")

        # 2A: All 19 WorldClim bioclimatic variables present and aligned with SQLite schema
        stopifnot(length(WC_BIO_LOOKUP) == 19L)
        stopifnot(identical(unname(WC_BIO_LOOKUP), PACKERA_DB_WC_BIO_COLS))

        # 2E: Coordinate-free non-demo SGMC & SSURGO classification (Carolina/Georgia coords must NOT trigger outcrop/granitic without attributes)
        sg_nd <- standardize_sgmc_lithology(c("", ""), lat = c(34.7, 35.9), lon = c(-80.5, -79.1), is_demo = FALSE)
        stopifnot(all(sg_nd == "Unknown_Lithology"))
        ss_nd <- classify_ssurgo_outcrop(muname = "", compname = "", lat = 34.7, lon = -80.5, is_demo = FALSE)
        stopifnot(ss_nd$outcrop_class == "Non-outcrop matrix" && !isTRUE(ss_nd$is_rock_outcrop))
        ss_ak <- classify_ssurgo_outcrop(muname = "Rock outcrop", compname = "Wake", soil_sand = 82, soil_ph = 5.1, bedrock_depth_cm = 22, rock_outcrop_pct = 28, lat = 64.5, lon = -148.0, is_demo = FALSE)
        stopifnot(ss_ak$outcrop_class == "Granitic flatrock" && isTRUE(ss_ak$is_rock_outcrop))

        # 2D: Cohort-adaptive circular anthesis gate for Arctic/alpine late-summer cohort (DOY 228..255, Lat 62..67.5)
        arctic_df <- data.frame(
          catalogNumber = paste0("ARC_", 1:12),
          species_standardized = rep(c("Packera hyperborea", "Packera contermina"), each = 6),
          determiner_tier = "Tier_1_Gold",
          latitude = seq(62.0, 67.5, length.out = 12),
          longitude = seq(-149.0, -138.0, length.out = 12),
          doy = c(228, 232, 235, 238, 241, 245, 236, 240, 244, 248, 252, 255),
          n_anthesis = c(4, 3, 5, 2, 4, 3, 3, 5, 4, 2, 0, 0),
          total_capitula = rep(5, 12),
          voucher_phenological_state = c(rep("Peak_Anthesis", 10), "Unknown", "Unknown"),
          stringsAsFactors = FALSE
        )
        ph_res <- evaluate_circular_phenology(arctic_df, demo = FALSE)
        stopifnot(sum(ph_res$df$is_verified_anthesis, na.rm = TRUE) >= 10L)
        stopifnot(all(is.finite(ph_res$df$delta_doy)))

        # 2B & 2C: Anchor-separated tie-breaking in screen_collinear_predictors_by_block + RF-importance-weighted z-scores
        set.seed(42)
        n <- 24
        sp <- rep(c("Packera hyperborea", "Packera contermina"), each = 12)
        # bio10 separates taxa strongly, while collinear bio1 is mostly shared + noise -> anchor tie-breaking must keep bio10 over bio1!
        base_t <- seq(-4, 4, length.out = n)
        df_tb <- data.frame(
          catalogNumber = paste0("TB_", seq_len(n)),
          species_standardized = sp,
          determiner_tier = "Tier_1_Gold",
          latitude = seq(61.0, 66.0, length.out = n),
          longitude = seq(-150.0, -140.0, length.out = n),
          bio10_mean_temp_warmest_qtr = base_t + ifelse(sp == "Packera hyperborea", -2.5, 2.5),
          bio1_temp_mean = base_t + ifelse(sp == "Packera hyperborea", -0.2, 0.2) + rnorm(n, 0, 0.05),
          fine_soil_cec = ifelse(sp == "Packera hyperborea", 8.0, 24.0) + rnorm(n, 0, 0.8),
          topo_tpi_30m = rnorm(n, 0, 1.0),
          stringsAsFactors = FALSE
        )
        # Ensure high collinearity between bio10 and bio1 (|r| > 0.75) while bio10 has much higher anchor F/KW statistic
        df_tb$bio1_temp_mean <- 0.92 * df_tb$bio10_mean_temp_warmest_qtr + rnorm(n, 0, 0.9)
        anc_g <- resolve_active_anchor_groups(df_tb)
        cand_tb <- c("bio1_temp_mean", "bio10_mean_temp_warmest_qtr", "fine_soil_cec", "topo_tpi_30m")
        blk_tb <- screen_collinear_predictors_by_block(df_tb[, cand_tb], cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = anc_g)
        stopifnot("bio10_mean_temp_warmest_qtr" %in% blk_tb$climatic_variables)

        rf_tb <- compute_spatial_rf_mems(df_tb, blk_tb$kept_vars, n_trees = 60, vif_res = blk_tb)
        prim_tb <- select_primary_predictors_from_rf(rf_tb, blk_tb, df_tb)
        stopifnot(length(prim_tb$primary_env_vars) >= 1L)
        cons_tb <- execute_crossmodal_consensus(
          df_tb, rf_res = rf_tb,
          primary_env_vars = prim_tb$primary_env_vars,
          primary_clim_vars = prim_tb$primary_clim_vars,
          primary_edaph_vars = prim_tb$primary_edaph_vars,
          rf_var_weights = prim_tb$rf_var_weights
        )
        stopifnot(all(c("edaphic_z_score", "climatic_z_score", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors") %in% names(cons_tb)))
        stopifnot(all(is.finite(cons_tb$edaphic_z_score)))
        cat("STEP06_COHORT_ADAPTIVE_OK\\n")
        """
        res = subprocess.run(["Rscript", "-e", r_code], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
        self.assertEqual(res.returncode, 0, f"Step 06 cohort-adaptive test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("STEP06_COHORT_ADAPTIVE_OK", res.stdout)


if __name__ == "__main__":
    unittest.main()
