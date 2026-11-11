#!/usr/bin/env python3
"""
===============================================================================
Script: test_multimodal_spatial_rf.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 05 (scripts/analysis/05_multimodal_spatial_rf.R):
    - 0-15 cm Trapezoidal Rhizosphere Integration & Non-Blended Coalesce Hierarchy
    - 6-Class USGS SGMC Bedrock Lithology & Azonal Contact Proximity
    - Coordinate-Uncertainty Buffered Focal Extraction (>10000m Training Exclusion)
    - Block-Wise Collinearity Pruning (|r| < 0.70, VIF < 5.0) & Partitioned Niche Tests
    - Multi-Class SGMC/SSURGO Substrate Contingency & Organic PCA-env Ecotypes
    - Authority-Anchored Spatial Block CV, OOF Posterior Pooling, & Step 06 Handshake
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
    """Test suite for Step 05 active R script (scripts/analysis/05_multimodal_spatial_rf.R)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"
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
            source("scripts/analysis/05_multimodal_spatial_rf.R")

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

              tex <- compute_soil_texture_log_ratios(sand_pct = c(70.0, 30.0, NA), clay_pct = c(10.0, 35.0, 20.0))
              stopifnot(all(c("soil_silt_pct", "ilr_texture_1", "ilr_texture_2") %in% names(tex)))
              stopifnot(!any(c("alr_sand_clay", "alr_silt_clay") %in% names(tex)))
              stopifnot(tex$ilr_texture_1[1] > 0 && tex$ilr_texture_1[2] < 0)
              stopifnot(is.na(tex$ilr_texture_1[3]) && is.na(tex$ilr_texture_2[3]))

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
              req_cols <- c("fine_bedrock_depth_cm", "fine_cfvo_pct", "fine_soil_cec", "ilr_texture_1", "ilr_texture_2")
              stopifnot(all(req_cols %in% names(df_out)))
              stopifnot(df_out$fine_bedrock_depth_cm[1] < df_out$fine_bedrock_depth_cm[2])
              stopifnot(df_out$fine_cfvo_pct[1] > df_out$fine_cfvo_pct[2])
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
              stopifnot(all(c("sgmc_lithology_class", "is_azonal_lithology") %in% names(df_geo)))
              stopifnot(!any(c("azonal_contact_dist_m", "azonal_proximity_index") %in% names(df_geo)))
              stopifnot(all(df_geo$is_azonal_lithology[1:4]))
              stopifnot(!df_geo$is_azonal_lithology[6])
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
                fine_cfvo_pct = c(32, 24, 12, 6),
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
              stopifnot(all(df_buf$max_cfvo_buffer_pct >= df_buf$fine_cfvo_pct))
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
                ilr_texture_1 = ifelse(is_dub, rnorm(n, 1.8, 0.22), rnorm(n, 0.2, 0.22)),
                ilr_texture_2 = ifelse(is_dub, rnorm(n, -0.3, 0.18), rnorm(n, 0.4, 0.18)),
                fine_bedrock_depth_cm = ifelse(is_dub, rnorm(n, 22, 6), rnorm(n, 135, 20)),
                fine_cfvo_pct = ifelse(is_dub, rnorm(n, 38, 6), rnorm(n, 9, 3)),
                topo_hli_30m = ifelse(is_dub, rnorm(n, 0.96, 0.04), rnorm(n, 0.72, 0.04)),
                stringsAsFactors = FALSE
              )
              df$bio1_collinear_dup <- df$bio1_temp_mean * 1.01 + rnorm(n, 0, 0.005)
              env_vars <- c("bio1_temp_mean", "bio1_collinear_dup", "bio4_temp_seasonality", "bio12_precip_annual",
                            "bio15_precip_seasonality", "fine_soil_ph", "ilr_texture_1", "ilr_texture_2",
                            "fine_bedrock_depth_cm", "fine_cfvo_pct", "topo_hli_30m")
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
              exp_active_n <- sum(df$determiner_tier %in% c("Tier_1_Gold", "Tier_2_Silver") & df$coordinateUncertaintyInMeters <= 10000)
              stopifnot((sum_df$N1[1] + sum_df$N2[1]) == exp_active_n)
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
              stopifnot(!("piedmont_outcrop_specialist_flag" %in% names(cons_df)))
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
              stopifnot(identical(unique(res_flat$analysis_group), "Packera aurea"))

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
                morph_tier_source = rep("Basal_Rosette", 3),
                basal_leaf_count = c(3L, 2L, 4L),
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
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                check=False,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"
        self.util_path = PROJECT_ROOT / "scripts" / "analysis" / "edaphic_fine_grain_utils.R"

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 05 R unit suite failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_r_script_phenological_and_fine_grain_edaphic_contracts(self):
        """Verify 05_multimodal_spatial_rf.R and edaphic_fine_grain_utils.R retain all phenological, ENM, SSURGO, POLARIS, 3DEP, and Spatial RF contracts."""
        self.assertTrue(self.r_script_path.exists(), "05_multimodal_spatial_rf.R must exist.")
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
            "fine_cfvo_pct",
            "ilr_texture_1",
            "ilr_texture_2",
            "topo_tpi_30m",
            "topo_hli_30m",
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
            self.assertIn(token, combined_r, f"Expected fine-grain token '{token}' in Step 05 R code.")

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
        self.assertTrue((self.pca_env_dir / "fig05_ecological_niche_phenology.pdf").exists())
        self.assertTrue((self.pca_env_dir / "fig05_ecological_niche_phenology.png").exists())

    def test_authority_anchored_spatial_block_cv_zero_synthetic_cache_and_step07_contract(self):
        """Verify authority-anchored Spatial RF training, zero synthetic disk caching under --demo, passive Bronze/high-uncertainty projection, and Step 07 compatibility."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            vouchers_csv = tmp / "curated_vouchers.csv"
            morph_csv = tmp / "morphometric_flags.csv"
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
                    "morph_tier_source": "Basal_Rosette",
                    "basal_leaf_count": 3,
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
            res = subprocess.run(cmd, cwd=str(PROJECT_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
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
                "fine_bedrock_depth_cm",
                "fine_cfvo_pct",
                "ilr_texture_1",
                "ilr_texture_2",
                "topo_hli_30m",
                "n_observed_modalities",
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

            # Verify Step 06 06_triage_dashboard_synthesis.R accepts the Step 05 outputs and completes without error
            step06_script = PROJECT_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"
            out_queue = tmp / "triage_queue.csv"
            out_report = tmp / "revision_summary.md"
            out_triage_plot = tmp / "triage_plot.pdf"
            cmd06 = [
                "Rscript",
                str(step06_script),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(morph_csv),
                "--multimodal-flags", str(out_flags),
                "--niche-summary", str(out_summary),
                "--output-queue", str(out_queue),
                "--output-report", str(out_report),
                "--output-plot", str(out_triage_plot),
            ]
            res07 = subprocess.run(cmd06, cwd=str(PROJECT_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
            self.assertEqual(res07.returncode, 0, f"Step 06 execution failed: {res07.stderr}\n{res07.stdout}")
            self.assertTrue(out_queue.exists(), "Step 06 triage_queue.csv must be generated from Step 05 outputs.")
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
        source("scripts/analysis/05_multimodal_spatial_rf.R")

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
        ph_res <- compute_phenological_anomalies(arctic_df, demo = FALSE)
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
        with tempfile.TemporaryDirectory() as tmp_ca:
            r_ca_file = Path(tmp_ca) / "test_cohort_adaptive.R"
            r_ca_file.write_text(r_code, encoding="utf-8")
            res = subprocess.run(["Rscript", str(r_ca_file)], cwd=str(PROJECT_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
        self.assertEqual(res.returncode, 0, f"Step 05 cohort-adaptive test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("STEP06_COHORT_ADAPTIVE_OK", res.stdout)

    def test_na_raster_cropping_deflate_and_maxent_spatial_rf_sdm_maps(self):
        """Verify automatic North American raster cropping/DEFLATE compression and default MaxEnt + Spatial RF SDM GeoTIFF & 4-panel figure generation."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmp_d:
            r_code = f"""
            suppressPackageStartupMessages({{ library(readr); library(dplyr); library(terra) }})
            source("scripts/core/db_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/05_multimodal_spatial_rf.R")

            tmp_root <- {repr(str(tmp_d))}
            env_dir <- file.path(tmp_root, "env")
            out_dir <- file.path(tmp_root, "outputs")
            dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

            # 1. Create a global uncropped GeoTIFF (-180..180, -60..90) and verify automatic cropping to NA (-130..-50, 14..66)
            r_glob <- terra::rast(xmin = -180, xmax = 180, ymin = -60, ymax = 90, res = 1.0, crs = "EPSG:4326")
            terra::values(r_glob) <- seq_len(terra::ncell(r_glob)) %% 50
            names(r_glob) <- "bio1_temp_mean"
            glob_tif <- file.path(env_dir, "wc2.1_bio_1.tif")
            terra::writeRaster(r_glob, glob_tif, overwrite = TRUE)

            crop_and_compress_na_rasters(env_dir = env_dir)
            r_after <- terra::rast(glob_tif)
            ex_after <- terra::ext(r_after)
            stopifnot(terra::xmin(ex_after) >= -130.1 && terra::xmax(ex_after) <= -49.9)
            stopifnot(terra::ymin(ex_after) >= 13.9 && terra::ymax(ex_after) <= 66.1)
            stopifnot(length(list.files(env_dir, pattern = "\\\\.crop_tmp\\\\.tif$")) == 0L)

            # 2. Fit MaxEnt + Spatial RF SDM and verify GeoTIFFs, sdm_model_performance.csv, and fig05b PDF/PNG
            set.seed(123)
            n <- 30
            sp <- rep(c("Packera dubia", "Packera anonyma"), each = 15)
            df_sdm <- data.frame(
              catalogNumber = paste0("SDM_", seq_len(n)),
              plant_individual_id = 0L,
              species_standardized = sp,
              analysis_group = sp,
              determiner_tier = c("Tier_1_Gold", rep("Tier_2_Silver", n - 1)),
              type_status = c("Holotype", rep("None", n - 1)),
              is_type_specimen = c(TRUE, rep(FALSE, n - 1)),
              Combined_Prior_Weight = c(1.0, rep(0.85, n - 1)),
              coordinateUncertaintyInMeters = rep(250, n),
              latitude = ifelse(sp == "Packera dubia", seq(32.5, 36.2, length.out = 15), seq(34.0, 38.5, length.out = 15)) + rnorm(n, 0, 0.15),
              longitude = ifelse(sp == "Packera dubia", seq(-81.5, -78.2, length.out = 15), seq(-83.5, -77.5, length.out = 15)) + rnorm(n, 0, 0.15),
              bio1_temp_mean = ifelse(sp == "Packera dubia", 16.5, 13.2) + rnorm(n, 0, 0.8),
              bio12_precip_annual = ifelse(sp == "Packera dubia", 1180, 1320) + rnorm(n, 0, 45),
              fine_soil_sand = ifelse(sp == "Packera dubia", 82.0, 44.0) + rnorm(n, 0, 4.0),
              fine_soil_ph = ifelse(sp == "Packera dubia", 4.9, 6.1) + rnorm(n, 0, 0.25),
              stringsAsFactors = FALSE
            )
            kept_v <- c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_sand", "fine_soil_ph")
            srf_res <- compute_spatial_rf_mems(df_sdm, kept_v, n_trees = 50)
            df_sdm$spatial_block_id <- srf_res$spatial_block_id
            df_sdm$rf_predicted_taxon <- srf_res$rf_predicted_taxon
            df_sdm$rf_predicted_prob <- srf_res$rf_predicted_prob
            df_sdm$rf_given_taxon_prob <- srf_res$rf_given_taxon_prob

            # 2a. Verify empty_env with demo=FALSE NEVER fabricates IDW GeoTIFF rasters, while still reporting tabular CV performance
            out_empty <- file.path(tmp_root, "outputs_empty")
            suppressWarnings(fit_and_project_maxent_sdm_maps(
              df = df_sdm,
              srf_res = srf_res,
              env_dir = file.path(tmp_root, "empty_env"),
              out_dir = out_empty,
              kept_vars = kept_v,
              primary_clim_vars = c("bio1_temp_mean", "bio12_precip_annual"),
              primary_edaph_vars = c("fine_soil_sand", "fine_soil_ph"),
              focal_taxon = "Packera dubia",
              intraspecific_mode = FALSE,
              demo = FALSE
            ))
            stopifnot(file.exists(file.path(out_empty, "reports", "sdm_model_performance.csv")))
            stopifnot(!file.exists(file.path(out_empty, "rasters", "sdm_packera_dubia_maxent_cloglog.tif")))

            # 2b. Provide real continuous GeoTIFF rasters in env_dir and verify full MaxEnt + Spatial RF projection (fig05b + fig05c)
            r_tmpl <- terra::rast(xmin = -90, xmax = -72, ymin = 28, ymax = 42, res = 0.25, crs = "EPSG:4326")
            xy_g <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
            r_b1 <- r_tmpl; terra::values(r_b1) <- 22.0 - 0.55 * (xy_g[, 2] - 28.0); names(r_b1) <- "bio1_temp_mean"
            r_b12 <- r_tmpl; terra::values(r_b12) <- 1100.0 + 18.0 * (xy_g[, 1] + 85.0); names(r_b12) <- "bio12_precip_annual"
            r_snd <- r_tmpl; terra::values(r_snd) <- pmin(92.0, pmax(15.0, 65.0 + 3.5 * (xy_g[, 1] + 80.0))); names(r_snd) <- "fine_soil_sand"
            r_ph <- r_tmpl; terra::values(r_ph) <- pmin(7.5, pmax(4.2, 5.4 - 0.12 * (xy_g[, 1] + 80.0))); names(r_ph) <- "fine_soil_ph"
            terra::writeRaster(c(r_b1, r_b12, r_snd, r_ph), file.path(env_dir, "regional_predictors.tif"), overwrite = TRUE)

            sdm_out <- fit_and_project_maxent_sdm_maps(
              df = df_sdm,
              srf_res = srf_res,
              env_dir = env_dir,
              out_dir = out_dir,
              kept_vars = kept_v,
              primary_clim_vars = c("bio1_temp_mean", "bio12_precip_annual"),
              primary_edaph_vars = c("fine_soil_sand", "fine_soil_ph"),
              focal_taxon = "Packera dubia",
              intraspecific_mode = FALSE,
              demo = FALSE
            )

            perf_csv <- file.path(out_dir, "reports", "sdm_model_performance.csv")
            fig_pdf <- file.path(out_dir, "figures", "fig05b_sdm_habitat_suitability_map.pdf")
            fig_png <- file.path(out_dir, "figures", "fig05b_sdm_habitat_suitability_map.png")
            grid_pdf <- file.path(out_dir, "figures", "fig05c_multispecies_maxent_sdm_grid.pdf")
            tif_maxent <- file.path(out_dir, "rasters", "sdm_packera_dubia_maxent_cloglog.tif")
            tif_ens <- file.path(out_dir, "rasters", "sdm_packera_dubia_ensemble.tif")

            stopifnot(file.exists(perf_csv) && file.info(perf_csv)$size > 0)
            stopifnot(file.exists(fig_pdf) && file.info(fig_pdf)$size > 0)
            stopifnot(file.exists(fig_png) && file.info(fig_png)$size > 0)
            stopifnot(file.exists(grid_pdf) && file.info(grid_pdf)$size > 0)
            stopifnot(file.exists(tif_maxent) && file.info(tif_maxent)$size > 0)
            stopifnot(file.exists(tif_ens) && file.info(tif_ens)$size > 0)

            perf_df <- readr::read_csv(perf_csv, show_col_types = FALSE)
            stopifnot("Packera dubia" %in% perf_df$Taxon_or_Morphotype)
            stopifnot(all(c("AUC_Block_MaxEnt", "CBI_MaxEnt", "AUC_Block_SpatialRF", "CBI_SpatialRF", "AUC_Block_Ensemble", "CBI_Ensemble") %in% names(perf_df)))
            cat("SDM_AND_CROP_OK\\n")
            """
            r_sdm_file = Path(tmp_d) / "test_na_raster_sdm.R"
            r_sdm_file.write_text(r_code, encoding="utf-8")
            res = subprocess.run(["Rscript", str(r_sdm_file)], cwd=str(PROJECT_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
            self.assertEqual(res.returncode, 0, f"SDM & raster cropping test failed:\n{res.stderr}\n{res.stdout}")
            self.assertIn("SDM_AND_CROP_OK", res.stdout)

    def test_noncircular_k1_ibe_background_similarity_and_empty_raster_hygiene(self):
        """Verify Step 05 non-circular K=1 ecotype gate, IBE partial RDA, accessible background M similarity, and zero synthetic POLARIS/3DEP fabrication when demo=FALSE."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmp_d:
            r_code = f"""
            suppressPackageStartupMessages({{ library(readr); library(dplyr) }})
            source("scripts/core/db_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/morphometrics/pub_figures_utils.R")
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/05_multimodal_spatial_rf.R")

            tmp_root <- {repr(str(tmp_d))}
            empty_env <- file.path(tmp_root, "empty_env")
            dir.create(file.path(empty_env, "polaris_30m"), recursive = TRUE, showWarnings = FALSE)
            dir.create(file.path(empty_env, "solus_100m"), recursive = TRUE, showWarnings = FALSE)
            dir.create(file.path(empty_env, "dem_3dep_30m"), recursive = TRUE, showWarnings = FALSE)

            # 1. Verify zero synthetic POLARIS / SOLUS / 3DEP / buffer fabrication when demo = FALSE and raster folders are empty
            df_real <- data.frame(
              catalogNumber = c("R1", "R2", "R3", "R4"),
              latitude = c(35.6, 35.7, 34.2, 33.9),
              longitude = c(-79.1, -80.2, -78.4, -79.5),
              coordinateUncertaintyInMeters = c(120, 450, 2500, 14000),
              muname = c("Rion sandy loam, 8 to 15 percent slopes", "Mecklenburg clay loam, 2 to 8 percent slopes", "Wake-Rock outcrop complex, 2 to 15 percent slopes", "Wagram loamy sand, 0 to 6 percent slopes"),
              compname = c("Rion", "Mecklenburg", "Wake", "Wagram"),
              ssurgo_ph = c(5.3, 6.2, 4.8, 5.1),
              ssurgo_sand_pct = c(66.0, 38.0, 76.0, 84.0),
              ssurgo_clay_pct = c(15.0, 32.0, 8.0, 6.0),
              ssurgo_cec = c(5.8, 19.5, 4.1, 3.6),
              ssurgo_aws_0_50cm = c(6.4, 7.2, 2.1, 4.0),
              ssurgo_bedrock_depth_cm = c(160.0, 110.0, 18.0, 200.0),
              ssurgo_rock_outcrop_pct = c(0.0, 0.0, 45.0, 0.0),
              stringsAsFactors = FALSE
            )
            cl_rion <- classify_ssurgo_outcrop(muname = df_real$muname[1], compname = "Rion", bedrock_depth_cm = 160.0, rock_outcrop_pct = 0.0, is_demo = FALSE)
            cl_meck <- classify_ssurgo_outcrop(muname = df_real$muname[2], compname = "Mecklenburg", bedrock_depth_cm = 110.0, rock_outcrop_pct = 0.0, is_demo = FALSE)
            cl_wake <- classify_ssurgo_outcrop(muname = df_real$muname[3], compname = "Wake", bedrock_depth_cm = 18.0, rock_outcrop_pct = 45.0, is_demo = FALSE)
            stopifnot(!isTRUE(cl_rion$is_rock_outcrop) && cl_rion$outcrop_class == "Non-outcrop matrix")
            stopifnot(!isTRUE(cl_meck$is_rock_outcrop) && cl_meck$outcrop_class == "Non-outcrop matrix")
            stopifnot(isTRUE(cl_wake$is_rock_outcrop) && cl_wake$outcrop_class == "Granitic flatrock")

            df_real$is_rock_outcrop <- c(FALSE, FALSE, TRUE, FALSE)
            df_real$outcrop_class <- c("Non-outcrop matrix", "Non-outcrop matrix", "Granitic flatrock", "Non-outcrop matrix")
            df_ps <- extract_polaris_and_solus_predictors(df_real, env_dir = empty_env, demo = FALSE)
            df_ps <- extract_micro_topography_3dep(df_ps, env_dir = empty_env, demo = FALSE)
            df_ps <- extract_sgmc_lithology(df_ps, env_dir = empty_env, demo = FALSE)
            df_ps <- extract_uncertainty_buffered_predictors(df_ps, env_dir = empty_env, demo = FALSE)

            stopifnot(all(is.na(df_ps$polaris_ph_30m)))
            stopifnot(all(is.na(df_ps$solus_fragvol_pct)))
            stopifnot(all(is.na(df_ps$topo_tpi_30m)))
            stopifnot(all(is.na(df_ps$topo_twi_30m)))
            stopifnot(all(is.finite(df_ps$topo_slope_deg)))
            stopifnot(all(is.na(df_ps$ph_q10_buffer)))
            stopifnot(all(is.na(df_ps$azonal_lithology_presence_buffer)))
            stopifnot(identical(df_ps$is_edaphic_eligible_coord, c(TRUE, TRUE, FALSE, FALSE)))
            stopifnot(identical(df_ps$is_clim_eligible_coord, c(TRUE, TRUE, TRUE, FALSE)))
            stopifnot(!df_ps$is_azonal_lithology[1] && !df_ps$is_azonal_lithology[2] && df_ps$is_azonal_lithology[3])

            # 2. Verify sympatrically intermixed Packera dubia with G=1 unimodal morphology resolves to K=1 (not false Ecotype_1 / Ecotype_2!)
            set.seed(2026)
            n_dub <- 60
            # Sympatric coordinates spanning Coastal Plain & Piedmont (intermixed substrates in same localities)
            base_lat <- rep(seq(32.0, 37.5, length.out = 30), each = 2) + rnorm(n_dub, 0, 0.03)
            base_lon <- rep(seq(-84.0, -77.0, length.out = 30), each = 2) + rnorm(n_dub, 0, 0.03)
            is_flatrock_microsite <- rep(c(FALSE, TRUE), 30)
            df_sym <- data.frame(
              catalogNumber = sprintf("DUB_%03d", seq_len(n_dub)),
              plant_individual_id = 1L,
              species_standardized = rep("Packera dubia", n_dub),
              determiner_tier = c(rep(c("Tier_1_Gold", "Tier_2_Silver"), 28), rep("Tier_3_Bronze", 3), "Tier_1_Gold"),
              coordinateUncertaintyInMeters = c(rep(180.0, n_dub - 1), 16000.0),
              latitude = base_lat,
              longitude = base_lon,
              gmm_cluster = rep(1L, n_dub),
              Sym_PC1 = rnorm(n_dub, 0, 1.0),
              Sym_PC2 = rnorm(n_dub, 0, 0.8),
              blade_length_to_width_ratio = rnorm(n_dub, 2.2, 0.25),
              doy = round(95 + 4.0 * (base_lat - 32.0) + rnorm(n_dub, 0, 6)),
              delta_doy = rnorm(n_dub, 0, 5.5),
              fine_soil_ph = rnorm(n_dub, 5.3, 0.35),
              fine_soil_sand = ifelse(is_flatrock_microsite, rnorm(n_dub, 74, 5), rnorm(n_dub, 78, 6)),
              fine_soil_clay = rnorm(n_dub, 9.0, 2.0),
              fine_soil_cec = rnorm(n_dub, 4.8, 1.1),
              fine_bedrock_depth_cm = ifelse(is_flatrock_microsite, rnorm(n_dub, 28, 8), rnorm(n_dub, 165, 20)),
              ssurgo_rock_outcrop_pct = ifelse(is_flatrock_microsite, runif(n_dub, 20, 55), rep(0.0, n_dub)),
              ilr_texture_1 = rnorm(n_dub, 1.6, 0.25),
              ilr_texture_2 = rnorm(n_dub, 0.1, 0.20),
              is_azonal_lithology = is_flatrock_microsite,
              bio10_mean_temp_warmest_qtr = 27.2 - 0.55 * (base_lat - 32.0) + rnorm(n_dub, 0, 0.35),
              bio12_precip_annual = rnorm(n_dub, 1220, 75),
              bio15_precip_seasonality = rnorm(n_dub, 16, 2.2),
              stringsAsFactors = FALSE
            )
            cand_v <- discover_candidate_eco_predictors(df_sym, include_pheno = FALSE)
            stopifnot(!("is_azonal_lithology" %in% cand_v))
            df_res <- resolve_analysis_group(df_sym, intraspecific_mode = TRUE)
            stopifnot(identical(unique(df_res$analysis_group), "Packera dubia"))
            stopifnot(identical(unique(resolve_plot_group_column(df_res)), "Packera dubia"))

            srf_k1 <- compute_spatial_rf_mems(df_res, cand_v, n_trees = 60, intraspecific_mode = TRUE)
            stopifnot(isTRUE(srf_k1$single_lineage_mode))
            stopifnot(all(srf_k1$rf_predicted_taxon == "Packera dubia"))
            stopifnot(srf_k1$rf_prediction_source[n_dub] == "Passive_Projection_High_Uncertainty")
            stopifnot(srf_k1$rf_prediction_source[n_dub - 1] == "Passive_Projection_Bronze")
            stopifnot(!is.null(srf_k1$ibe_summary) && nrow(srf_k1$ibe_summary) >= 1L)
            stopifnot(all(srf_k1$ibe_summary$Partial_RDA_p > 0.01))

            niche_k1 <- run_warren_niche_identity_tests(df_res, cand_v, n_perm = 15, intraspecific_mode = TRUE)
            stopifnot(nrow(niche_k1$summary) == 1L)
            stopifnot(niche_k1$summary$Schoener_D[1] == 1.0)
            fig_k1_pdf <- file.path(tmp_root, "k1_fig", "spatial_rf_niche_summary.pdf")
            export_spatial_rf_figures(df_res, srf_k1, niche_k1, fig_k1_pdf)
            stopifnot(file.exists(fig_k1_pdf))
            stopifnot(file.exists(file.path(tmp_root, "k1_fig", "fig05_ecological_niche_phenology.png")))
            cat("NONCIRCULAR_K1_IBE_AND_HYGIENE_OK\n")
            """
            r_script_file = Path(tmp_d) / "test_noncircular_k1_hygiene.R"
            r_script_file.write_text(r_code, encoding="utf-8")
            res = subprocess.run(["Rscript", str(r_script_file)], cwd=str(PROJECT_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
            self.assertEqual(res.returncode, 0, f"Non-circular K=1 & edaphic hygiene test failed:\n{res.stderr}\n{res.stdout}")
            self.assertIn("NONCIRCULAR_K1_IBE_AND_HYGIENE_OK", res.stdout)

    def test_step05_eastern_na_figure_refinements(self):
        sdm_src = (PROJECT_ROOT / "scripts" / "analysis" / "step05_sdm_maxent.R").read_text(encoding="utf-8")
        pheno_src = (PROJECT_ROOT / "scripts" / "analysis" / "step05_phenology_ecotypes.R").read_text(encoding="utf-8")

        self.assertIn("keep_pairs <- head(unique(perm_df$pair[order(!grepl(\"dubia\", perm_df$pair, ignore.case = TRUE), seq_len(nrow(perm_df)))]), 6L)", sdm_src)
        self.assertIn("focal_artifacts$X_comb[focal_artifacts$p_vec == 1L, , drop = FALSE]", sdm_src)
        self.assertIn("med_row <- apply(pres_mat, 2, stats::median, na.rm = TRUE)", sdm_src)
        self.assertIn("decimalLatitude >= 22.0", pheno_src)
        self.assertIn("doy = pmin(365.0, pmax(1.0, b0_use + b1_use * lat_seq))", pheno_src)
        self.assertIn("hopkins_expected = pmin(365.0, pmax(1.0, pred_min_lat + 4.0 * (lat_seq - min_lat)))", pheno_src)


    def test_step05_zero_trust_scientific_and_anti_leakage_suite(self):
        """
        Zero-trust scientific, mathematical, and anti-leakage regression suite for Step 05:
          1. Zero target-label leakage when northern/non-CONUS rows have NA edaphic predictors.
          2. Strict 23-predictor physiological ENM whitelist (excludes bio8/9/18/19, bio3/7,
             elevation, delta_doy, uncertainty-buffer order statistics, and raw duplicate soil layers).
          3. Hopkins' Bioclimatic Law exactness (+4.0 d / 400 ft = +0.0328084 d/m, NOT 0.0125 d/m)
             plus application of b_lon and b_year in compute_phenological_anomalies().
          4. Exact zero-influence invariance to Tier_3_Bronze and is_confident_label_noise corruption.
          5. Sterile/fruiting/undated specimen protection (never flagged as Phenological_Anomaly_Outlier)
             and passive retention of non-georeferenced vouchers (Passive_Ungeoreferenced).
          6. Real 30 m DEM McCune & Keon (2002) HLI (SW > NE) and D8 hydrological TWI (valley > ridge).
        """
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmp_d:
            r_code = f"""
            suppressPackageStartupMessages({{ library(readr); library(dplyr); library(terra) }})
            source("scripts/core/db_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/morphometrics/pub_figures_utils.R")
            source("scripts/analysis/edaphic_fine_grain_utils.R")
            source("scripts/analysis/05_multimodal_spatial_rf.R")

            tmp_root <- {repr(str(tmp_d))}

            # 1. Strict 23-predictor whitelist check
            df_w <- data.frame(
              catalogNumber = paste0("W", 1:12),
              fine_soil_ph = seq(4.5, 6.7, length.out = 12),
              ilr_texture_1 = seq(-1.2, 1.5, length.out = 12),
              ilr_texture_2 = seq(-0.6, 0.8, length.out = 12),
              fine_soil_sand = seq(25, 80, length.out = 12),
              fine_soil_clay = seq(35, 8, length.out = 12),
              soil_silt_pct = seq(40, 12, length.out = 12),
              fine_soil_cec = seq(4, 22, length.out = 12),
              fine_cfvo_pct = seq(2, 42, length.out = 12),
              fine_bedrock_depth_cm = seq(20, 180, length.out = 12),
              polaris_ph_30m = seq(4.5, 6.7, length.out = 12),
              ssurgo_ph = seq(4.6, 6.8, length.out = 12),
              ph_q10_buffer = seq(4.1, 6.2, length.out = 12),
              min_bedrock_depth_buffer_cm = seq(10, 140, length.out = 12),
              elevation_m = seq(100, 1200, length.out = 12),
              topo_slope_deg = seq(1.5, 24.0, length.out = 12),
              topo_tpi_30m = seq(-3.0, 4.0, length.out = 12),
              topo_hli_30m = seq(0.62, 1.02, length.out = 12),
              topo_twi_30m = seq(4.5, 11.5, length.out = 12),
              bio1_temp_mean = seq(7.0, 18.0, length.out = 12),
              bio3_isothermality = seq(28, 39, length.out = 12),
              bio7_temp_annual_range = seq(26, 37, length.out = 12),
              bio8_mean_temp_wettest_qtr = seq(15, 24, length.out = 12),
              bio9_mean_temp_driest_qtr = seq(2, 14, length.out = 12),
              bio12_precip_annual = seq(850, 1450, length.out = 12),
              bio18_precip_warmest_qtr = seq(240, 420, length.out = 12),
              bio19_precip_coldest_qtr = seq(180, 340, length.out = 12),
              delta_doy = seq(-15, 15, length.out = 12),
              doy = seq(95, 150, length.out = 12),
              stringsAsFactors = FALSE
            )
            cand_w <- discover_candidate_eco_predictors(df_w)
            forbidden_w <- c(
              "bio3_isothermality", "bio7_temp_annual_range",
              "bio8_mean_temp_wettest_qtr", "bio9_mean_temp_driest_qtr",
              "bio18_precip_warmest_qtr", "bio19_precip_coldest_qtr",
              "delta_doy", "doy", "elevation_m", "ph_q10_buffer", "min_bedrock_depth_buffer_cm",
              "polaris_ph_30m", "ssurgo_ph", "soil_silt_pct", "fine_soil_sand", "fine_soil_clay"
            )
            stopifnot(!any(forbidden_w %in% cand_w))
            expected_w <- c(
              "fine_soil_ph", "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm",
              "topo_slope_deg", "topo_tpi_30m", "topo_hli_30m", "topo_twi_30m",
              "bio1_temp_mean", "bio12_precip_annual"
            )
            stopifnot(all(expected_w %in% cand_w))

            # 2. Hopkins' Bioclimatic Law (+4.0 d / 400 ft = +0.0328084 d/m) + b_lon + b_year exactness
            set.seed(77)
            n_h <- 24
            lat_h <- seq(34.0, 38.0, length.out = n_h)
            lon_h <- rep(c(-83.0, -78.0), length.out = n_h) # 5 deg lon span -> +4.0 d eastward
            elev_h <- c(rep(200.0, n_h - 2), 1200.0, 1400.0) # +1000..1200 m Blue Ridge sheets
            yr_h <- rep(seq(1960, 2020, length.out = 6), 4)
            # Generate exact Hopkins synthetic truth: -0.10 d/yr secular shift
            true_doy <- round(115.0 + 4.0 * (lat_h - 34.0) + 0.80 * (lon_h - median(lon_h)) + (4.0 / 121.92) * (elev_h - median(elev_h)) - 0.10 * (yr_h - median(yr_h)))
            df_hop <- data.frame(
              catalogNumber = paste0("HOP_", seq_len(n_h)),
              species_standardized = "Packera dubia",
              determiner_tier = "Tier_1_Gold",
              decimalLatitude = lat_h,
              decimalLongitude = lon_h,
              elevation_m = elev_h,
              year = yr_h,
              doy = true_doy,
              reproductiveCondition = "fl",
              voucher_phenological_state = "anthesis",
              total_capitula = 5,
              n_anthesis = 4,
              stringsAsFactors = FALSE
            )
            hop_res <- compute_phenological_anomalies(df_hop, demo = FALSE)
            stopifnot(abs(hop_res$baseline$b_elev - (4.0 / 121.92)) < 1e-6)
            stopifnot(abs(hop_res$baseline$b_lon - 0.80) < 1e-6)
            stopifnot(hop_res$baseline$b_year < -0.04)
            # High-elevation Blue Ridge specimens (rows 23-24) must have near-zero residual (< 3 days), NOT +20 days!
            stopifnot(all(abs(hop_res$df$delta_doy[23:24]) < 3.0))

            # 3. Zero target-label leakage on trans-border NA soil values (<50% sparse pruning + >=50% fold-median non-leakage)
            set.seed(101)
            n_lk <- 80
            sp_lk <- rep(c("Packera dubia", "Packera paupercula"), each = 40)
            lat_lk <- rep(seq(36.0, 47.0, length.out = 40), 2) + rnorm(n_lk, 0, 0.05)
            lon_lk <- rep(seq(-84.0, -76.0, length.out = 40), 2) + rnorm(n_lk, 0, 0.05)
            # 3a. Sparse (<50% finite) CONUS-only predictor (6/80 = 7.5% finite) must be pruned rather than median-spiked
            bed_sparse <- rep(NA_real_, n_lk)
            bed_sparse[c(1, 2, 3)] <- c(30.0, 32.0, 34.0)
            bed_sparse[c(41, 42, 43)] <- c(150.0, 155.0, 160.0)
            df_lk_sparse <- data.frame(
              catalogNumber = paste0("LK_", seq_len(n_lk)),
              species_standardized = sp_lk,
              analysis_group = sp_lk,
              determiner_tier = "Tier_1_Gold",
              coordinateUncertaintyInMeters = 200.0,
              latitude = lat_lk,
              longitude = lon_lk,
              bio1_temp_mean = rnorm(n_lk, 11.0, 1.0),
              bio12_precip_annual = rnorm(n_lk, 1050.0, 40.0),
              fine_bedrock_depth_cm = bed_sparse,
              stringsAsFactors = FALSE
            )
            rf_noleak <- compute_spatial_rf_mems(df_lk_sparse, c("bio1_temp_mean", "bio12_precip_annual", "fine_bedrock_depth_cm"), n_trees = 100)
            stopifnot(!("fine_bedrock_depth_cm" %in% rf_noleak$importance$Predictor))
            stopifnot(rf_noleak$Spatial_Block_CV_Accuracy < 0.72)

            # 3b. Retained predictor (44/80 = 55% finite, 45% NA) with bimodal symmetric overlap (Bayes ceiling = 52.5%):
            # Leaky species-median imputation would assign 39.0 to all sp1 NAs and 161.0 to all sp2 NAs (>75% acc);
            # taxon-blind fold-median imputation keeps accuracy near 50% (< 0.72)!
            bed_55 <- rep(NA_real_, n_lk)
            bed_55[1:22] <- c(rep(39.0, 12), rep(160.0, 10))
            bed_55[41:62] <- c(rep(40.0, 10), rep(161.0, 12))
            df_lk_55 <- df_lk_sparse
            df_lk_55$fine_bedrock_depth_cm <- bed_55
            rf_noleak_55 <- compute_spatial_rf_mems(df_lk_55, c("bio1_temp_mean", "bio12_precip_annual", "fine_bedrock_depth_cm"), n_trees = 100)
            stopifnot("fine_bedrock_depth_cm" %in% rf_noleak_55$importance$Predictor)
            stopifnot(rf_noleak_55$Spatial_Block_CV_Accuracy < 0.72)

            # 4. Tier_3_Bronze & is_confident_label_noise zero-influence invariance
            df_clean <- data.frame(
              catalogNumber = paste0("INV_", 1:24),
              species_standardized = rep(c("Packera dubia", "Packera anonyma"), each = 12),
              analysis_group = rep(c("Packera dubia", "Packera anonyma"), each = 12),
              determiner_tier = rep(c(rep("Tier_1_Gold", 5), rep("Tier_2_Silver", 4), rep("Tier_3_Bronze", 2), "Tier_1_Gold"), 2),
              is_confident_label_noise = rep(c(rep(FALSE, 11), TRUE), 2),
              coordinateUncertaintyInMeters = 180.0,
              latitude = rep(seq(34.0, 36.5, length.out = 12), 2),
              longitude = rep(seq(-81.0, -78.5, length.out = 12), 2),
              fine_soil_ph = c(rnorm(12, 4.9, 0.2), rnorm(12, 6.3, 0.2)),
              ilr_texture_1 = c(rnorm(12, 1.6, 0.2), rnorm(12, -0.3, 0.2)),
              bio1_temp_mean = c(rnorm(12, 16.0, 0.4), rnorm(12, 13.5, 0.4)),
              is_verified_anthesis = c(rep(TRUE, 10), FALSE, FALSE, rep(TRUE, 10), FALSE, FALSE),
              doy = c(rep(108, 10), NA, 260, rep(132, 10), NA, 25),
              cda_predicted_taxon = rep(c("Packera dubia", "Packera anonyma"), each = 12),
              cda_posterior_prob = 0.94,
              stringsAsFactors = FALSE
            )
            df_corrupt <- df_clean
            bad_rows <- which(df_corrupt$determiner_tier == "Tier_3_Bronze" | df_corrupt$is_confident_label_noise)
            df_corrupt$fine_soil_ph[bad_rows] <- c(9.5, 2.5, 9.8, 2.2, 9.2, 2.8)
            df_corrupt$ilr_texture_1[bad_rows] <- c(-5.0, 5.0, -4.8, 4.8, -5.2, 5.2)
            df_corrupt$bio1_temp_mean[bad_rows] <- c(-10.0, 35.0, -8.0, 38.0, -12.0, 40.0)

            rf_c1 <- compute_spatial_rf_mems(df_clean, c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), n_trees = 80)
            rf_c2 <- compute_spatial_rf_mems(df_corrupt, c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), n_trees = 80)
            act_rows <- setdiff(seq_len(24), bad_rows)
            stopifnot(max(abs(rf_c1$oof_prob_matrix[act_rows, ] - rf_c2$oof_prob_matrix[act_rows, ])) < 1e-9)

            cons_c1 <- execute_crossmodal_consensus(df_clean, rf_res = rf_c1, primary_env_vars = c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), primary_edaph_vars = c("fine_soil_ph", "ilr_texture_1"), primary_clim_vars = "bio1_temp_mean")
            cons_c2 <- execute_crossmodal_consensus(df_corrupt, rf_res = rf_c2, primary_env_vars = c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean"), primary_edaph_vars = c("fine_soil_ph", "ilr_texture_1"), primary_clim_vars = "bio1_temp_mean")
            stopifnot(max(abs(cons_c1$edaphic_z_score[act_rows] - cons_c2$edaphic_z_score[act_rows])) < 1e-9)
            # Sterile / non-anthesis rows (11, 12, 23, 24) with extreme DOY 260 / 25 or NA must NEVER be flagged as Phenological_Anomaly_Outlier!
            stopifnot(!any(cons_c1$triage_category[c(11, 12, 23, 24)] == "Phenological_Anomaly_Outlier"))
            stopifnot(all(cons_c1$n_observed_modalities[c(11, 12, 23, 24)] == 2L))

            # 5. Real 30 m DEM McCune & Keon (2002) HLI & D8 Hydrological TWI
            dem_dir <- file.path(tmp_root, "dem_env", "dem_3dep_30m")
            dir.create(dem_dir, recursive = TRUE, showWarnings = FALSE)
            r_dem <- terra::rast(xmin = -79.20, xmax = -78.80, ymin = 35.40, ymax = 35.80, nrows = 160, ncols = 160, crs = "EPSG:4326")
            xy_d <- terra::xyFromCell(r_dem, seq_len(terra::ncell(r_dem)))
            # Create a V-shaped valley draining south along lon = -79.00, plus an E-W ridge at lat = 35.70 for SSW vs NNE aspect check
            v_valley <- 200.0 + 1800.0 * abs(xy_d[, 1] - (-79.00)) + 400.0 * (xy_d[, 2] - 35.40)
            terra::values(r_dem) <- v_valley
            terra::writeRaster(r_dem, file.path(dem_dir, "N35W080.tif"), overwrite = TRUE)

            pts_dem <- data.frame(
              catalogNumber = c("VALLEY_FLOOR", "UPPER_FLANK_WEST_FACING", "UPPER_FLANK_EAST_FACING"),
              latitude = c(35.55, 35.55, 35.55),
              longitude = c(-79.00, -78.95, -79.05),
              stringsAsFactors = FALSE
            )
            dem_out <- extract_micro_topography_3dep(pts_dem, env_dir = file.path(tmp_root, "dem_env"), demo = FALSE)
            stopifnot(all(is.finite(dem_out$topo_hli_30m)))
            stopifnot(all(is.finite(dem_out$topo_twi_30m)))
            stopifnot(all(is.finite(dem_out$elevation_m)))
            # Valley floor (-79.00) receives D8 upslope convergence -> higher TWI & negative TPI vs upper flanks!
            stopifnot(dem_out$topo_twi_30m[1] > dem_out$topo_twi_30m[2])
            stopifnot(dem_out$topo_tpi_30m[1] < 0.0)
            # East side of N-S valley (-78.95) faces West/SW (warmer HLI) vs West side (-79.05) facing East/NE (cooler HLI)!
            stopifnot(dem_out$topo_hli_30m[2] > dem_out$topo_hli_30m[3])

            # 6. Verify ilr_texture_1/2 edaphic group routing & fast IPP MaxEnt convergence with glmnet.control hygiene
            imp_ilr <- rf_c1$importance
            ilr_grp <- imp_ilr$Group[imp_ilr$Predictor == "ilr_texture_1"]
            stopifnot(identical(as.character(ilr_grp), "[A] Fine-Grain Edaphic & Outcrop Lithology"))
            top_ilr <- select_top_rf_predictors(imp_ilr, df_clean)
            stopifnot(top_ilr$top_edaph_1 %in% c("fine_soil_ph", "ilr_texture_1"))
            t0_mx <- proc.time()[["elapsed"]]
            mx_sep <- .fit_single_maxent_model(
              p_vec = c(rep(1L, 15), rep(0L, 15)),
              X_df = data.frame(
                bio1_temp_mean = c(rnorm(15, 16.5, 0.3), rnorm(15, 11.0, 0.3)),
                fine_soil_ph = c(rnorm(15, 4.8, 0.15), rnorm(15, 6.6, 0.15))
              ),
              w_vec = c(1.0, rep(0.85, 14), rep(1.0, 15))
            )
            dt_mx <- proc.time()[["elapsed"]] - t0_mx
            stopifnot(dt_mx < 10.0)
            stopifnot(all(is.finite(mx_sep$predict_cloglog(data.frame(bio1_temp_mean = c(16.5, 11.0), fine_soil_ph = c(4.8, 6.6))))))
            if (requireNamespace("glmnet", quietly = TRUE)) stopifnot(glmnet::glmnet.control()$fdev > 0)

            cat("STEP05_ZERO_TRUST_SUITE_OK\n")
            """
            r_file = Path(tmp_d) / "test_step05_zero_trust.R"
            r_file.write_text(r_code, encoding="utf-8")
            res = subprocess.run(["Rscript", str(r_file)], cwd=str(PROJECT_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
            self.assertEqual(res.returncode, 0, f"Step 05 zero-trust suite failed:\n{res.stderr}\n{res.stdout}")
            self.assertIn("STEP05_ZERO_TRUST_SUITE_OK", res.stdout)


if __name__ == "__main__":
    unittest.main()
