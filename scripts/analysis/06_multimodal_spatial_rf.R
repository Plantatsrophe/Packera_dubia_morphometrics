#!/usr/bin/env Rscript
# ==============================================================================
# Script: 06_multimodal_spatial_rf.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Multimodal Spatial Macroecology, SoilGrids 250m & WorldClim v2.1 Integration,
#   Cross-Modal Consensus Checking, Spatial Random Forests (MEMs via spatialRF),
#   and ENMTools Warren's Niche Identity Permutation Testing.
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("soilDB", quietly = TRUE)) library(soilDB)
  if (requireNamespace("elevatr", quietly = TRUE)) library(elevatr)
  if (requireNamespace("geodata", quietly = TRUE)) library(geodata)
  if (requireNamespace("ranger", quietly = TRUE)) library(ranger)
  if (requireNamespace("adespatial", quietly = TRUE)) library(adespatial)
  if (requireNamespace("blockCV", quietly = TRUE)) library(blockCV)
  if (requireNamespace("ecospat", quietly = TRUE)) library(ecospat)
  if (requireNamespace("digest", quietly = TRUE)) library(digest)
  if (requireNamespace("spatialRF", quietly = TRUE)) library(spatialRF)
  if (requireNamespace("ENMTools", quietly = TRUE)) library(ENMTools)
  if (requireNamespace("vegan", quietly = TRUE)) library(vegan)
  if (requireNamespace("terra", quietly = TRUE)) library(terra)
  if (requireNamespace("sf", quietly = TRUE)) library(sf)
  if (requireNamespace("randomForest", quietly = TRUE)) library(randomForest)
  if (requireNamespace("circular", quietly = TRUE)) library(circular)
  if (requireNamespace("MASS", quietly = TRUE)) library(MASS)
  if (requireNamespace("ggridges", quietly = TRUE)) library(ggridges)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("tibble", quietly = TRUE)) library(tibble)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("gridExtra", quietly = TRUE)) library(gridExtra)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})

if (file.exists("scripts/morphometrics/pub_figures_utils.R")) {
  source("scripts/morphometrics/pub_figures_utils.R")
}
if (file.exists("scripts/morphometrics/authority_weighting_utils.R")) {
  try(source("scripts/morphometrics/authority_weighting_utils.R"), silent = TRUE)
}
if (file.exists("scripts/analysis/edaphic_fine_grain_utils.R")) {
  source("scripts/analysis/edaphic_fine_grain_utils.R")
} else {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0) {
    cand_util <- file.path(dirname(sub("^--file=", "", file_arg[1])), "edaphic_fine_grain_utils.R")
    if (file.exists(cand_util)) source(cand_util)
  }
}
if (!exists("prune_collinear_predictors", mode = "function") &&
    exists("screen_collinear_predictors_by_block", mode = "function")) {
  prune_collinear_predictors <- screen_collinear_predictors_by_block
}

OKABE_ITO_FALLBACK_HEX <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#332288", "#88CCEE", "#44AA99")

resolve_dynamic_group_palette <- function(groups_vec) {
  u_grp <- unique(as.character(stats::na.omit(groups_vec)))
  u_grp <- u_grp[nzchar(u_grp)]
  if (length(u_grp) == 0L) return(character(0))
  if (exists("get_dynamic_packera_palette", mode = "function")) {
    pal <- get_dynamic_packera_palette(u_grp)
    miss <- setdiff(u_grp, names(pal))
    if (length(miss) == 0L && !any(is.na(pal[u_grp]))) return(pal)
  }
  stats::setNames(rep_len(OKABE_ITO_FALLBACK_HEX, length(u_grp)), u_grp)
}

STEP06_PREDICTOR_LABEL_DICT <- c(
  "fine_soil_sand"                   = "Soil Sand Fraction (%, 30m POLARIS/SSURGO)",
  "fine_soil_clay"                   = "Soil Clay Fraction (%, 30m POLARIS/SSURGO)",
  "fine_soil_ph"                     = "Soil pH in H2O (30m POLARIS/SSURGO)",
  "fine_soil_cec"                    = "Fine-Grain Soil CEC (cmol+/kg)",
  "fine_bedrock_depth_cm"            = "Depth to Lithic Bedrock (cm, 10m-30m)",
  "fine_rock_outcrop_prob"           = "Rock Outcrop Probability (10m-30m)",
  "bdricm"                           = "Depth to Bedrock Restriction (cm, 30m/100m)",
  "resdept_r"                        = "SSURGO Restriction Depth (cm)",
  "cfvo"                             = "Coarse Fragment Volume (%, 0-15cm)",
  "awc"                              = "Available Water Capacity (cm/cm, 0-15cm)",
  "alr_sand_clay"                    = "Texture Log-Ratio ln(Sand/Clay)",
  "alr_silt_clay"                    = "Texture Log-Ratio ln(Silt/Clay)",
  "is_azonal_lithology"              = "USGS SGMC Azonal Lithology Indicator",
  "azonal_contact_dist_m"            = "Distance to Azonal Lithologic Contact (m)",
  "azonal_proximity_index"           = "Azonal Lithologic Proximity Index",
  "ph_q10_buffer"                    = "Buffered 10th-Pct Soil pH (100-1000m)",
  "ph_q90_buffer"                    = "Buffered 90th-Pct Soil pH (100-1000m)",
  "min_bedrock_depth_buffer_cm"      = "Buffered Min Bedrock Depth (cm)",
  "max_cfvo_buffer_pct"              = "Buffered Max Coarse Fragments (%)",
  "azonal_lithology_presence_buffer" = "Buffered Azonal Lithology Presence",
  "ssurgo_rock_outcrop_pct"          = "SSURGO Rock Outcrop Component (%)",
  "ssurgo_bedrock_depth_cm"          = "SSURGO Lithic Bedrock Depth (cm)",
  "ssurgo_aws_0_50cm"                = "Available Water Storage 0-50cm (cm)",
  "topo_tpi_30m"                     = "30m Topographic Position Index (TPI)",
  "topo_twi_30m"                     = "30m Topographic Wetness Index (TWI)",
  "topo_roughness_30m"               = "30m Terrain Roughness Index",
  "topo_slope_deg"                   = "30m Terrain Slope (deg)",
  "topo_elev_m"                      = "30m Elevation (m)",
  "soil_ph"                          = "SoilGrids 250m Soil pH (H2O)",
  "soil_sand"                        = "SoilGrids 250m Sand Fraction (%)",
  "soil_cec"                         = "SoilGrids 250m Cation Exchange (cmol/kg)",
  "soil_bulk_density"                = "SoilGrids 250m Bulk Density (g/cm3)",
  "bio1_temp_mean"                   = "Mean Annual Temp (BIO1, deg C)",
  "bio4_temp_seasonality"            = "Temp Seasonality (BIO4, SD x 100)",
  "bio12_precip_annual"              = "Annual Precipitation (BIO12, mm)",
  "bio15_precip_seasonality"         = "Precipitation Seasonality (BIO15, CV)",
  "bio18_precip_warmest_qt"          = "Warmest-Quarter Precip (BIO18, mm)",
  "bio19_precip_coldest_qt"          = "Coldest-Quarter Precip (BIO19, mm)",
  "delta_doy"                        = "Hopkins Phenological Anomaly (Delta DOY)"
)

.upstream_format_predictor_labels <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels else NULL
format_step06_predictor_labels <- function(var_names) {
  vapply(as.character(var_names), function(nm) {
    if (is.na(nm) || !nzchar(nm)) return("")
    if (nm %in% names(STEP06_PREDICTOR_LABEL_DICT)) return(STEP06_PREDICTOR_LABEL_DICT[[nm]])
    if (grepl("^MEM[0-9]+$", nm)) return(sub("^MEM", "Spatial Eigenvector MEM-", nm))
    if (!is.null(.upstream_format_predictor_labels)) return(.upstream_format_predictor_labels(nm)[1])
    gsub("_", " ", nm)
  }, character(1), USE.NAMES = FALSE)
}
format_predictor_labels <- format_step06_predictor_labels

resolve_step06_authority_weights <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(numeric(0))
  if (file.exists("scripts/morphometrics/authority_weighting_utils.R") &&
      !exists("compute_empirical_bayes_authority_weights", mode = "function")) {
    try(source("scripts/morphometrics/authority_weighting_utils.R"), silent = TRUE)
  }
  if (exists("compute_empirical_bayes_authority_weights", mode = "function")) {
    eb_wts <- tryCatch({
      res <- tryCatch(
        compute_empirical_bayes_authority_weights(df, mode = "continuous_prior_posterior"),
        error = function(e1) {
          compute_empirical_bayes_authority_weights(
            df = df,
            taxon_col = "species_standardized",
            tier_col = "determiner_tier",
            mode = "continuous_prior_posterior"
          )
        }
      )
      w_out <- NULL
      if (is.numeric(res) && length(res) == n) {
        w_out <- as.numeric(res)
      } else if (is.data.frame(res) && nrow(res) == n) {
        w_col <- intersect(c("authority_weight", "empirical_bayes_weight", "eb_authority_weight", "weight"), names(res))[1]
        if (!is.na(w_col)) w_out <- as.numeric(res[[w_col]])
      } else if (is.list(res)) {
        for (k_nm in c("authority_weight", "weights", "empirical_bayes_weight")) {
          if (!is.null(res[[k_nm]]) && length(res[[k_nm]]) == n) {
            w_out <- as.numeric(res[[k_nm]])
            break
          }
        }
      }
      if (!is.null(w_out) && length(w_out) == n && any(is.finite(w_out) & w_out > 0)) {
        w_out[!is.finite(w_out) | is.na(w_out)] <- 0.25
        pmax(0.05, w_out)
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(eb_wts)) return(eb_wts)
  }

  # Inline Empirical Bayes Beta-Binomial shrinkage fallback on (species_standardized, determiner_tier)
  sp_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unknown", n)
  sp_vec[is.na(sp_vec) | !nzchar(sp_vec)] <- "Unknown"
  tier_raw <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_3_Bronze", n)
  tier_raw[is.na(tier_raw) | !nzchar(tier_raw)] <- "Tier_3_Bronze"

  tier_std <- ifelse(grepl("Tier_1|Gold|Tier 1", tier_raw, ignore.case = TRUE), "Tier_1_Gold",
              ifelse(grepl("Tier_2|Silver|Tier 2", tier_raw, ignore.case = TRUE), "Tier_2_Silver", "Tier_3_Bronze"))

  alpha0 <- ifelse(tier_std == "Tier_1_Gold", 19.0, ifelse(tier_std == "Tier_2_Silver", 8.0, 3.0))
  beta0  <- ifelse(tier_std == "Tier_1_Gold", 1.0,  ifelse(tier_std == "Tier_2_Silver", 2.0, 3.0))

  # Empirical cross-modal agreement trials per specimen when auxiliary predictions exist
  has_cda <- "cda_predicted_taxon" %in% names(df)
  has_vis <- "vision_predicted_label" %in% names(df)
  obs_succ <- rep(0.0, n)
  obs_trials <- rep(0.0, n)
  if (has_cda) {
    cda_v <- as.character(df$cda_predicted_taxon)
    ok_c <- !is.na(cda_v) & nzchar(cda_v) & sp_vec != "Unknown"
    obs_succ[ok_c] <- obs_succ[ok_c] + as.numeric(cda_v[ok_c] == sp_vec[ok_c])
    obs_trials[ok_c] <- obs_trials[ok_c] + 1.0
  }
  if (has_vis) {
    vis_v <- as.character(df$vision_predicted_label)
    ok_v <- !is.na(vis_v) & nzchar(vis_v) & sp_vec != "Unknown"
    obs_succ[ok_v] <- obs_succ[ok_v] + as.numeric(vis_v[ok_v] == sp_vec[ok_v])
    obs_trials[ok_v] <- obs_trials[ok_v] + 1.0
  }

  cell_key <- paste(sp_vec, tier_std, sep = "||")
  u_cells <- unique(cell_key)
  w_out <- numeric(n)
  for (ck in u_cells) {
    idx <- which(cell_key == ck)
    a_prior <- alpha0[idx[1]]
    b_prior <- beta0[idx[1]]
    k_cell <- sum(obs_succ[idx], na.rm = TRUE)
    m_cell <- sum(obs_trials[idx], na.rm = TRUE)
    post_mean <- (a_prior + k_cell) / (a_prior + b_prior + m_cell)
    w_out[idx] <- pmax(0.05, pmin(1.0, post_mean))
  }
  return(round(w_out, 4))
}

# ------------------------------------------------------------------------------
# 1. CLI Argument Parsing & Taxonomic Standardization
# ------------------------------------------------------------------------------
parse_args_robust <- function() {
  run_dir <- Sys.getenv("PACKERA_RUN_DIR", unset = "")
  rp <- function(rel) if (nzchar(run_dir)) file.path(run_dir, rel) else rel
  option_list <- list(
    optparse::make_option(c("-v", "--vouchers"), type = "character",
      default = rp("data/tables/curated_vouchers.csv"), help = "Curated vouchers CSV [default: %default]"),
    optparse::make_option(c("-m", "--morphometrics"), type = "character",
      default = rp("data/tables/morphometrics_misidentification_flags.csv"), help = "Morphometrics flags CSV [default: %default]"),
    optparse::make_option(c("-n", "--vision-audit"), type = "character",
      default = rp("data/tables/label_noise_audit.csv"), help = "Cleanlab vision audit CSV [default: %default]"),
    optparse::make_option(c("--pheno-states"), type = "character",
      default = rp("data/tables/voucher_phenological_states.csv"), help = "Voucher phenological states CSV [default: %default]"),
    optparse::make_option(c("-e", "--env-dir"), type = "character",
      default = "data/environmental", help = "Directory with environmental rasters [default: %default]"),
    optparse::make_option(c("-f", "--output-flags"), type = "character",
      default = rp("data/tables/multimodal_conflict_flags.csv"), help = "Output conflict CSV [default: %default]"),
    optparse::make_option(c("-p", "--output-plot"), type = "character",
      default = rp("outputs/figures/spatial_rf_niche_importance.pdf"), help = "Output PDF figure [default: %default]"),
    optparse::make_option(c("-s", "--output-summary"), type = "character",
      default = rp("outputs/reports/multimodal_spatial_rf_summary.csv"), help = "Output summary CSV [default: %default]"),
    optparse::make_option(c("--pheno-anomalies"), type = "character",
      default = rp("data/tables/phenological_anomalies.csv"), help = "Output voucher-level phenological anomalies CSV [default: %default]"),
    optparse::make_option(c("--pheno-summary"), type = "character",
      default = rp("outputs/reports/phenological_anomaly_summary.csv"), help = "Output phenological anomaly ANOVA/Tukey summary CSV [default: %default]"),
    optparse::make_option(c("--pheno-plot"), type = "character",
      default = rp("outputs/figures/phenological_latitudinal_anomaly.pdf"), help = "Output 2-panel phenological anomaly figure PDF [default: %default]"),
    optparse::make_option(c("--ssurgo-cache"), type = "character",
      default = "data/environmental/ssurgo_cache", help = "Directory for SSURGO cached queries [default: %default]"),
    optparse::make_option(c("--ssurgo-table"), type = "character",
      default = rp("data/tables/ssurgo_edaphic_validation.csv"), help = "Output SSURGO validation CSV [default: %default]"),
    optparse::make_option(c("--contingency-table"), type = "character",
      default = rp("outputs/reports/micro_edaphic_outcrop_contingency.csv"), help = "Output outcrop contingency CSV [default: %default]"),
    optparse::make_option(c("--batch-size"), type = "integer",
      default = 75, help = "Batch size for USDA SDA queries [default: %default]"),
    optparse::make_option(c("-k", "--permutations"), type = "integer",
      default = 100, help = "Warren's Identity Test permutations [default: %default]"),
    optparse::make_option(c("-t", "--n-trees"), type = "integer",
      default = 500, help = "Number of trees in Random Forests [default: %default]"),
    optparse::make_option(c("--intraspecific-mode"), action = "store_true",
      default = FALSE, help = "Enable pure single-taxon / intraspecific morphotype delimitation mode [default: %default]"),
    optparse::make_option(c("--focal-taxon"), type = "character",
      default = "", help = "Optional focal Packera taxon name [default: %default]"),
    optparse::make_option(c("--demo"), action = "store_true",
      default = FALSE, help = "Enable synthetic fallback data for demonstration/test runs [default: %default]")
  )
  if (requireNamespace("optparse", quietly = TRUE)) {
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list))
    names(opts) <- gsub("-", "_", names(opts))
    return(opts)
  }
  return(list(
    vouchers = "data/tables/curated_vouchers.csv", morphometrics = "data/tables/morphometrics_misidentification_flags.csv",
    vision_audit = "data/tables/label_noise_audit.csv", pheno_states = "data/tables/voucher_phenological_states.csv",
    env_dir = "data/environmental",
    output_flags = "data/tables/multimodal_conflict_flags.csv", output_plot = "outputs/figures/spatial_rf_niche_importance.pdf",
    output_summary = "outputs/reports/multimodal_spatial_rf_summary.csv",
    pheno_anomalies = "data/tables/phenological_anomalies.csv",
    pheno_summary = "outputs/reports/phenological_anomaly_summary.csv",
    pheno_plot = "outputs/figures/phenological_latitudinal_anomaly.pdf",
    ssurgo_cache = "data/environmental/ssurgo_cache", ssurgo_table = "data/tables/ssurgo_edaphic_validation.csv",
    contingency_table = "outputs/reports/micro_edaphic_outcrop_contingency.csv", batch_size = 75,
    permutations = 100, n_trees = 500, intraspecific_mode = FALSE, focal_taxon = "", demo = FALSE
  ))
}

.base_std_taxon <- if (exists("standardize_taxon_name", mode = "function")) standardize_taxon_name else NULL
standardize_packera_taxon <- function(s_vec) {
  if (!is.null(.base_std_taxon)) {
    res <- .base_std_taxon(s_vec)
    return(ifelse(is.na(res), "Unknown", res))
  }
  vapply(as.character(s_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("Unknown")
    sc <- trimws(s)
    if (tolower(sc) %in% c("none", "nan", "null", "na", "unknown", "unassigned", "unassigned_packera")) return("Unknown")
    if (grepl("^Morphotype_K[0-9]+$|^Ecotype_|^LatTertile_|_Coastal$|_Outcrop$|\\[K[0-9]+\\]$", sc, ignore.case = TRUE)) return(sc)
    sl <- tolower(sc)
    if (grepl("^(?:packera|p\\.)\\s*(?:(?:cf\\.?|aff\\.?)\\s+)?(?:x\\s+|×\\s*)?[a-z-]+", sl)) {
      ep <- sub("^(?:packera|p\\.)\\s*(?:(?:cf\\.?|aff\\.?)\\s+)?(?:x\\s+|×\\s*)?([a-z-]+).*$", "\\1", sl)
      ep <- sub("^[x×-]+", "", ep)
      if (nzchar(ep) && !(ep %in% c("sp", "spp", "cf", "aff", "indet", "unknown", "na"))) return(paste("Packera", ep))
    }
    "Unknown"
  }, character(1), USE.NAMES = FALSE)
}

get_active_packera_taxa <- function(species_vec) {
  u_sp <- sort(unique(as.character(stats::na.omit(species_vec))))
  u_sp <- u_sp[nzchar(u_sp) & u_sp != "Unknown" & grepl("^Packera\\s+", u_sp)]
  if (length(u_sp) > 0) return(u_sp)
  u_any <- sort(unique(as.character(stats::na.omit(species_vec))))
  u_any <- u_any[nzchar(u_any) & u_any != "Unknown"]
  return(u_any)
}

is_canonical_prod_path <- function(p) {
  if (is.null(p) || !nzchar(p)) return(FALSE)
  norm_p <- normalizePath(p, mustWork = FALSE)
  repo_data <- normalizePath("data", mustWork = FALSE)
  repo_outputs <- normalizePath("outputs", mustWork = FALSE)
  startsWith(norm_p, paste0(repo_data, "/")) || startsWith(norm_p, paste0(repo_outputs, "/"))
}

write_csv_atomic <- function(df, target_path, demo = FALSE) {
  if (isTRUE(demo) && is_canonical_prod_path(target_path)) {
    message(sprintf("[DEMO MODE] Skipping write of synthetic records to canonical production table: %s", target_path))
    return(invisible(FALSE))
  }
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(target_path, ".tmp")
  readr::write_csv(df, tmp_path)
  file.rename(tmp_path, target_path)
  return(invisible(TRUE))
}

# ------------------------------------------------------------------------------
# 2. Environmental Extraction (SoilGrids 250m 0-15cm Rhizosphere & Uncertainty-Buffered WorldClim v2.1)
# ------------------------------------------------------------------------------
extract_rast_buffered_mean <- function(rast_obj, pts_vect, radii_m) {
  pt_vals <- terra::extract(rast_obj, pts_vect)
  buf_vals <- tryCatch({
    r_num <- suppressWarnings(as.numeric(radii_m))
    r_num[is.na(r_num) | !is.finite(r_num)] <- 250.0
    r_num <- pmax(100.0, pmin(1000.0, r_num))
    pts_buf <- terra::buffer(pts_vect, width = r_num)
    terra::extract(rast_obj, pts_buf, fun = mean, na.rm = TRUE)
  }, error = function(e) NULL)
  if (is.null(buf_vals)) return(pt_vals)
  for (col_nm in setdiff(names(pt_vals), "ID")) {
    if (col_nm %in% names(buf_vals)) {
      bv <- suppressWarnings(as.numeric(buf_vals[[col_nm]]))
      pv <- suppressWarnings(as.numeric(pt_vals[[col_nm]]))
      pt_vals[[col_nm]] <- ifelse(is.finite(bv), bv, pv)
    }
  }
  return(pt_vals)
}

extract_soilgrids_rhizosphere_col <- function(sg_vals, var_pattern) {
  nms <- names(sg_vals)
  idx_0_5 <- grep(paste0(var_pattern, ".*(0-5|0_5)"), nms, ignore.case = TRUE)
  idx_5_15 <- grep(paste0(var_pattern, ".*(5-15|5_15)"), nms, ignore.case = TRUE)
  if (length(idx_0_5) > 0 && length(idx_5_15) > 0) {
    v_0_5 <- suppressWarnings(as.numeric(sg_vals[[idx_0_5[1]]]))
    v_5_15 <- suppressWarnings(as.numeric(sg_vals[[idx_5_15[1]]]))
    if (exists("compute_rhizosphere_0_15cm", mode = "function")) {
      return(compute_rhizosphere_0_15cm(v_0_5, v_5_15))
    }
    return((5.0 * v_0_5 + 10.0 * v_5_15) / 15.0)
  }
  idx_any <- grep(var_pattern, nms, ignore.case = TRUE)
  if (length(idx_any) > 0) return(suppressWarnings(as.numeric(sg_vals[[idx_any[1]]])))
  return(rep(NA_real_, nrow(sg_vals)))
}

extract_environmental_layers <- function(df, env_dir, demo = FALSE) {
  message("Integrating SoilGrids 250m (0-15cm Rhizosphere) & Uncertainty-Buffered WorldClim v2.1 Bioclimatic layers...")
  if (exists("flag_coordinate_uncertainty", mode = "function")) {
    df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0)
  } else {
    if (!("coordinateUncertaintyInMeters" %in% names(df))) df$coordinateUncertaintyInMeters <- NA_real_
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    df$high_coord_uncertainty_flag <- !is.na(u_m) & u_m > 10000.0
    df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
    df$buffer_radius_m <- ifelse(is.na(u_m) | !is.finite(u_m) | u_m <= 0, 250.0, pmax(100.0, pmin(1000.0, u_m)))
  }

  vars <- c("soil_ph", "soil_cec", "soil_sand", "soil_bulk_density",
            "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual",
            "bio15_precip_seasonality", "bio18_precip_warmest_qt", "bio19_precip_coldest_qt")
  
  # Check if GeoTIFF rasters exist in env_dir (flat or cached climate/ & soil_world/ subdirectories)
  raster_files <- list.files(env_dir, pattern = "\\.(tif|tiff|grd)$", full.names = TRUE, ignore.case = TRUE)
  extract_err <- NULL
  if (length(raster_files) > 0 && requireNamespace("terra", quietly = TRUE)) {
    extracted_df <- tryCatch({
      env_stack <- terra::rast(raster_files)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      extracted <- extract_rast_buffered_mean(env_stack, pts, df$buffer_radius_m)
      for (v in names(extracted)[names(extracted) != "ID"]) df[[v]] <- extracted[[v]]
      for (base_v in c("soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")) {
        c05 <- paste0(base_v, "_0_5cm"); c515 <- paste0(base_v, "_5_15cm")
        if (all(c(c05, c515) %in% names(df))) {
          df[[base_v]] <- if (exists("compute_rhizosphere_0_15cm", mode = "function")) {
            compute_rhizosphere_0_15cm(df[[c05]], df[[c515]])
          } else {
            (5.0 * as.numeric(df[[c05]]) + 10.0 * as.numeric(df[[c515]])) / 15.0
          }
        }
      }
      message(sprintf("Extracted %d raster layers (with adaptive buffer radius r_i in [100, 1000]m) from %s.", length(raster_files), env_dir))
      df
    }, error = function(e) {
      extract_err <<- e$message
      NULL
    })
    if (!is.null(extracted_df)) return(extracted_df)
  }

  sub_files <- list.files(env_dir, pattern = "\\.(tif|tiff|grd)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  wc_files <- sub_files[grepl("bio_[0-9]+", basename(sub_files), ignore.case = TRUE)]
  sg_files <- sub_files[grepl("phh2o|cec|sand|bdod", basename(sub_files), ignore.case = TRUE)]
  if (length(wc_files) > 0 && length(sg_files) > 0 && requireNamespace("terra", quietly = TRUE)) {
    extracted_cached <- tryCatch({
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      wc_rast <- terra::rast(wc_files)
      sg_rast <- terra::rast(sg_files)
      wc_vals <- extract_rast_buffered_mean(wc_rast, pts, df$buffer_radius_m)
      sg_vals <- extract_rast_buffered_mean(sg_rast, pts, df$buffer_radius_m)
      df$bio1_temp_mean <- wc_vals[[grep("bio_1$|bio_01$", names(wc_vals))[1]]]
      df$bio4_temp_seasonality <- wc_vals[[grep("bio_4$|bio_04$", names(wc_vals))[1]]]
      df$bio12_precip_annual <- wc_vals[[grep("bio_12$", names(wc_vals))[1]]]
      df$bio15_precip_seasonality <- wc_vals[[grep("bio_15$", names(wc_vals))[1]]]
      df$bio18_precip_warmest_qt <- wc_vals[[grep("bio_18$", names(wc_vals))[1]]]
      df$bio19_precip_coldest_qt <- wc_vals[[grep("bio_19$", names(wc_vals))[1]]]
      df$soil_ph <- extract_soilgrids_rhizosphere_col(sg_vals, "phh2o")
      df$soil_cec <- extract_soilgrids_rhizosphere_col(sg_vals, "cec")
      df$soil_sand <- extract_soilgrids_rhizosphere_col(sg_vals, "sand")
      df$soil_bulk_density <- extract_soilgrids_rhizosphere_col(sg_vals, "bdod")
      message(sprintf("Extracted %d WorldClim and %d SoilGrids cached rasters (0-15cm rhizosphere & buffered) from %s.", length(wc_files), length(sg_files), env_dir))
      df
    }, error = function(e) {
      extract_err <<- e$message
      NULL
    })
    if (!is.null(extracted_cached)) return(extracted_cached)
  }

  # Automated GeoTIFF Acquisition via geodata when local rasters are absent and --demo is FALSE
  if (!isTRUE(demo) && requireNamespace("geodata", quietly = TRUE) && requireNamespace("terra", quietly = TRUE) &&
      (is_canonical_prod_path(env_dir) || identical(Sys.getenv("PACKERA_ENABLE_GEODATA_DOWNLOAD"), "TRUE"))) {
    extracted_geo <- tryCatch({
      message(sprintf("Attempting automated WorldClim v2.1 & SoilGrids 250m (0-5cm & 5-15cm) acquisition via geodata into '%s'...", env_dir))
      dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
      wc_rast <- geodata::worldclim_global(var = "bio", res = 2.5, path = env_dir)
      sg_0_5 <- lapply(c("phh2o", "cec", "sand", "bdod"), function(v) tryCatch(geodata::soil_world(var = v, depth = 5, path = env_dir), error = function(e) NULL))
      sg_5_15 <- lapply(c("phh2o", "cec", "sand", "bdod"), function(v) geodata::soil_world(var = v, depth = 15, path = env_dir))
      sg_list <- c(Filter(Negate(is.null), sg_0_5), sg_5_15)
      sg_rast <- terra::rast(sg_list)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      wc_vals <- extract_rast_buffered_mean(wc_rast, pts, df$buffer_radius_m)
      sg_vals <- extract_rast_buffered_mean(sg_rast, pts, df$buffer_radius_m)
      df$bio1_temp_mean <- wc_vals[[grep("bio_1$|bio_01$", names(wc_vals))[1]]]
      df$bio4_temp_seasonality <- wc_vals[[grep("bio_4$|bio_04$", names(wc_vals))[1]]]
      df$bio12_precip_annual <- wc_vals[[grep("bio_12$", names(wc_vals))[1]]]
      df$bio15_precip_seasonality <- wc_vals[[grep("bio_15$", names(wc_vals))[1]]]
      df$bio18_precip_warmest_qt <- wc_vals[[grep("bio_18$", names(wc_vals))[1]]]
      df$bio19_precip_coldest_qt <- wc_vals[[grep("bio_19$", names(wc_vals))[1]]]
      df$soil_ph <- extract_soilgrids_rhizosphere_col(sg_vals, "phh2o")
      df$soil_cec <- extract_soilgrids_rhizosphere_col(sg_vals, "cec")
      df$soil_sand <- extract_soilgrids_rhizosphere_col(sg_vals, "sand")
      df$soil_bulk_density <- extract_soilgrids_rhizosphere_col(sg_vals, "bdod")
      df
    }, error = function(e) {
      extract_err <<- e$message
      NULL
    })
    if (!is.null(extracted_geo)) return(extracted_geo)
  }

  if (!isTRUE(demo)) {
    stop(sprintf(
      paste0(
        "WorldClim v2.1 / SoilGrids 250m GeoTIFF rasters are missing or unreadable in '%s'%s. ",
        "Remediation: populate '%s' with valid SoilGrids 250m and WorldClim v2.1 GeoTIFF (.tif) layers, ",
        "or pass --demo explicitly for synthetic demonstration runs."
      ),
      env_dir,
      if (!is.null(extract_err)) paste0(" (", extract_err, ")") else "",
      env_dir
    ))
  }

  message("[DEMO MODE] Generating synthetic SoilGrids (0-15cm rhizosphere) & WorldClim layers...")
  # High-fidelity eco-geographic empirical transfer functions (calibrated on SEUS / Midwest eco-regions)
  set.seed(42)
  lat <- df$latitude; lon <- df$longitude; n <- nrow(df)
  reg <- if ("regional_group" %in% names(df)) df$regional_group else rep("Unknown", n)
  
  # Pedology (SoilGrids 250m 0-5cm and 5-15cm calibrated horizons integrated over 0-15cm rhizosphere)
  is_flatrock <- grepl("Flatrock", reg, ignore.case = TRUE) | (lon > -84 & lon < -79 & lat > 33 & lat < 36.5)
  is_sandhill <- grepl("Sandhill|Coastal", reg, ignore.case = TRUE) | (lon > -82 & lon < -75 & lat < 36.5)
  is_prairie  <- grepl("Midwest|Prairie", reg, ignore.case = TRUE) | (lon < -88)
  
  ph_0_5 <- ifelse(is_sandhill, 4.6 + 0.3 * sin(lat), ifelse(is_flatrock, 5.1 + 0.25 * cos(lon),
            ifelse(is_prairie, 7.0 + 0.3 * sin(lat), 6.3 + 0.3 * cos(lat)))) + rnorm(n, 0, 0.15)
  ph_5_15 <- ph_0_5 + 0.15 + rnorm(n, 0, 0.10)
  df$soil_ph <- round(pmax(3.8, pmin(8.5, (5.0 * ph_0_5 + 10.0 * ph_5_15) / 15.0)), 2)
  
  cec_0_5 <- ifelse(is_sandhill, 5.6 + 1.2 * abs(sin(lat)), ifelse(is_flatrock, 12.0 + 2.0 * cos(lat),
             ifelse(is_prairie, 27.2 + 3.5 * sin(lon), 18.6 + 2.5 * cos(lat)))) + rnorm(n, 0, 1.0)
  cec_5_15 <- cec_0_5 - 0.8 + rnorm(n, 0, 0.8)
  df$soil_cec <- round(pmax(1.5, pmin(50.0, (5.0 * cec_0_5 + 10.0 * cec_5_15) / 15.0)), 1)

  sand_0_5 <- ifelse(is_sandhill, 83.0 + 5.0 * sin(lon), ifelse(is_flatrock, 53.0 + 6.0 * cos(lat),
              ifelse(is_prairie, 24.5 + 4.5 * sin(lat), 38.5 + 5.0 * cos(lon)))) + rnorm(n, 0, 3.0)
  sand_5_15 <- sand_0_5 - 1.5 + rnorm(n, 0, 2.5)
  df$soil_sand <- round(pmax(5.0, pmin(98.0, (5.0 * sand_0_5 + 10.0 * sand_5_15) / 15.0)), 1)

  bd_0_5 <- ifelse(is_sandhill, 1.44, ifelse(is_flatrock, 1.28, 1.20)) + rnorm(n, 0, 0.04)
  bd_5_15 <- bd_0_5 + 0.06 + rnorm(n, 0, 0.04)
  df$soil_bulk_density <- round((5.0 * bd_0_5 + 10.0 * bd_5_15) / 15.0, 2)

  # Bioclimatics (WorldClim v2.1 calibrated parameters)
  df$bio1_temp_mean <- round(24.5 - 0.72 * (lat - 28.0) - 0.05 * abs(lon + 82) + rnorm(n, 0, 0.4), 2)
  df$bio4_temp_seasonality <- round(520 + 28.5 * (lat - 30.0) + rnorm(n, 0, 15), 1)
  df$bio12_precip_annual <- round(1450 - 15.0 * (lat - 30.0) + 12.0 * (lon + 85.0) + rnorm(n, 0, 35), 1)
  df$bio15_precip_seasonality <- round(22.0 + 0.8 * (lat - 30.0) + rnorm(n, 0, 1.5), 1)
  df$bio18_precip_warmest_qt <- round(df$bio12_precip_annual * (0.32 + 0.02 * sin(lat)) + rnorm(n, 0, 15), 1)
  df$bio19_precip_coldest_qt <- round(df$bio12_precip_annual * (0.22 - 0.01 * cos(lon)) + rnorm(n, 0, 12), 1)
  return(df)
}

# ------------------------------------------------------------------------------
# 2b. USDA SSURGO Vector Pedology Query & Caching Engine (soilDB / SDA Integration)
# ------------------------------------------------------------------------------
classify_ssurgo_outcrop <- function(muname, compname = "", taxsuborder = "", taxgreatgroup = "", compkind = "", lat = NA, lon = NA) {
  desc <- tolower(paste(muname, compname, taxsuborder, taxgreatgroup, compkind))
  piedmont_lithic_regex <- "(?i)(rock outcrop|wake|louisburg|ashlar|rion|wateree|wilkes|iredell|mecklenburg|poindexter|enon|goldston|manteo|boulders|granite|diabase)"
  is_piedmont_series <- grepl(piedmont_lithic_regex, paste(muname, compname), perl = TRUE)
  is_misc_rock <- grepl("miscellaneous area", tolower(compkind)) && grepl("rock|rubble|outcrop|boulder", desc)

  is_granite <- is_piedmont_series || is_misc_rock ||
                grepl("wake|louisburg|flatrock|granit|catao|hempstead|ashlar|rion|wateree", desc) ||
                (grepl("appling", desc) && grepl("flatrock", desc)) ||
                (grepl("dystrudepts|rock outcrop", desc) && !is.na(lon) && lon > -84.5 && lon < -77.5 && !is.na(lat) && lat > 32.5 && lat < 37.0)
  is_limestone <- grepl("gladeville|barfield|cedars|talbott|limestone|dolomite|glade", desc) ||
                  (grepl("hapludolls|rendolls", desc) && !is.na(lon) && lon > -88.5 && lon < -83.0 && !is.na(lat) && lat > 34.5 && lat < 37.5)
  is_sandstone <- grepl("ramsey|lily|sewanee|sandstone|siliceous|chert|shale", desc)
  is_serpentine <- grepl("chrome|trego|serpentine|ultramafic|iredell|mecklenburg|wilkes|enon|poindexter|diabase", desc)

  is_outcrop <- grepl("rock outcrop|flatrock|glade|lithic", desc) || is_granite || is_limestone || is_sandstone || is_serpentine

  outcrop_class <- if (is_limestone) {
    "Limestone glade"
  } else if (is_sandstone) {
    "Sandstone/Siliceous"
  } else if (grepl("chrome|trego|serpentine|ultramafic", desc)) {
    "Serpentine"
  } else if (grepl("iredell|mecklenburg|wilkes|enon|poindexter|diabase", desc)) {
    "Diabase/Mafic glade"
  } else if (is_granite || is_outcrop) {
    "Granitic flatrock"
  } else {
    "Non-outcrop matrix"
  }

  return(list(is_rock_outcrop = is_outcrop, outcrop_class = outcrop_class))
}

query_ssurgo_sda_batch <- function(batch_df, cache_dir, batch_idx, max_retries = 4, demo = FALSE, batch_size = 75) {
  canonical_cols <- c(
    "catalogNumber", "latitude", "longitude",
    "musym", "muname", "taxsuborder", "taxgreatgroup",
    "is_rock_outcrop", "outcrop_class",
    "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm",
    "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph",
    "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm"
  )
  normalize_ssurgo_cols <- function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    for (col in c("catalogNumber", "musym", "muname", "taxsuborder", "taxgreatgroup", "outcrop_class")) {
      if (!(col %in% names(d))) {
        d[[col]] <- if (col == "outcrop_class") "Non-outcrop matrix" else NA_character_
      } else {
        d[[col]] <- as.character(d[[col]])
      }
    }
    for (col in c("latitude", "longitude")) {
      if (!(col %in% names(d))) {
        d[[col]] <- NA_real_
      } else {
        d[[col]] <- as.numeric(d[[col]])
      }
    }
    if (!("is_rock_outcrop" %in% names(d))) {
      d$is_rock_outcrop <- FALSE
    } else {
      d$is_rock_outcrop <- as.logical(d$is_rock_outcrop)
      d$is_rock_outcrop[is.na(d$is_rock_outcrop)] <- FALSE
    }
    d$outcrop_class[is.na(d$outcrop_class) | !nzchar(d$outcrop_class)] <- "Non-outcrop matrix"
    if (exists("enrich_ssurgo_pedology_from_series", mode = "function")) {
      d <- enrich_ssurgo_pedology_from_series(d)
    } else {
      for (sc in c("ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct",
                   "ssurgo_clay_pct", "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm")) {
        if (!(sc %in% names(d))) d[[sc]] <- NA_real_
      }
    }
    return(d[, canonical_cols, drop = FALSE])
  }

  # Build deterministic batch cache key
  pt_str <- paste(batch_df$catalogNumber, batch_df$latitude, batch_df$longitude, collapse = ";")
  batch_hash <- if (requireNamespace("digest", quietly = TRUE)) {
    substr(digest::digest(pt_str, algo = "md5"), 1, 12)
  } else {
    sprintf("b%04d", batch_idx)
  }
  
  cache_file_rds <- file.path(cache_dir, sprintf("ssurgo_batch_%03d_%s.rds", batch_idx, batch_hash))
  cache_file_csv <- file.path(cache_dir, sprintf("ssurgo_batch_%03d_%s.csv", batch_idx, batch_hash))
  
  if (file.exists(cache_file_rds)) {
    message(sprintf("Loading cached SSURGO batch %d from %s", batch_idx, cache_file_rds))
    return(normalize_ssurgo_cols(readRDS(cache_file_rds)))
  }
  if (file.exists(cache_file_csv)) {
    message(sprintf("Loading cached SSURGO batch %d from %s", batch_idx, cache_file_csv))
    return(normalize_ssurgo_cols(as.data.frame(readr::read_csv(cache_file_csv, show_col_types = FALSE))))
  }
  
  pts_df <- batch_df[, c("catalogNumber", "longitude", "latitude")]
  n_pts <- nrow(pts_df)
  query_success <- FALSE
  res_df <- NULL
  delays <- c(1, 2, 4, 8)
  
  # Attempt live query via soilDB if available (chunked to batch_size = 75 points with 50m ST_Buffer(ST_Point(lon, lat), 0.00045))
  if (!isTRUE(demo) && requireNamespace("soilDB", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
    chunk_size <- max(1L, min(as.integer(batch_size), 100L))
    chunk_indices <- split(seq_len(n_pts), ceiling(seq_len(n_pts) / chunk_size))
    chunk_results <- list()
    all_chunks_ok <- TRUE

    for (c_idx in seq_along(chunk_indices)) {
      sub_pts <- pts_df[chunk_indices[[c_idx]], , drop = FALSE]
      n_sub <- nrow(sub_pts)
      chunk_ok <- FALSE
      chunk_df <- NULL

      for (attempt in seq_len(max_retries)) {
        tryCatch({
          message(sprintf("SDA 50m-Buffered Query Batch %d chunk %d/%d (n=%d, attempt %d/%d)...",
                          batch_idx, c_idx, length(chunk_indices), n_sub, attempt, max_retries))
          sub_pts$gid <- seq_len(n_sub)

          # Strategy A: Direct batched WKT ST_Buffer(ST_Point(lon, lat), 0.00045) SQL via soilDB::SDA_query()
          if (exists("build_ssurgo_buffered_sda_sql", mode = "function") &&
              exists("aggregate_ssurgo_tabular_to_points", mode = "function")) {
            buf_sql <- build_ssurgo_buffered_sda_sql(sub_pts, buffer_deg = 0.00045)
            buf_tab <- tryCatch(soilDB::SDA_query(buf_sql), error = function(e_b) NULL)
            if (is.data.frame(buf_tab) && nrow(buf_tab) > 0 && "mukey" %in% names(buf_tab)) {
              agg_df <- aggregate_ssurgo_tabular_to_points(sub_pts, buf_tab, classify_fn = classify_ssurgo_outcrop)
              chunk_df <- normalize_ssurgo_cols(agg_df)
              chunk_ok <- TRUE
              break
            }
          }

          # Strategy B: 50m buffered sf::st_buffer(pts_sf, dist = 0.00045) via soilDB::SDA_spatialQuery() + component/horizon/restriction join
          pts_sf <- sf::st_as_sf(sub_pts, coords = c("longitude", "latitude"), crs = 4326)
          pts_buf_sf <- tryCatch(suppressWarnings(sf::st_buffer(pts_sf, dist = 0.00045)), error = function(e) pts_sf)
          sp_res <- tryCatch(
            soilDB::SDA_spatialQuery(pts_buf_sf, what = "mukey", db = "SSURGO", byFeature = TRUE),
            error = function(e_sp) {
              if (!grepl("numbers of columns", conditionMessage(e_sp), fixed = TRUE)) stop(e_sp)
              message(sprintf("  [Batch %d] Unmapped point triggered soilDB rbind mismatch; rescuing valid points via dplyr::bind_rows...", batch_idx))
              mini_idx <- split(seq_len(n_sub), ceiling(seq_len(n_sub) / 5L))
              mini_list <- list()
              for (m_i in seq_along(mini_idx)) {
                idx_vec <- mini_idx[[m_i]]
                m_sf <- pts_sf[idx_vec, , drop = FALSE]
                m_res <- tryCatch({
                  r_m <- soilDB::SDA_spatialQuery(m_sf, what = "mukey", db = "SSURGO", byFeature = TRUE)
                  if (is.data.frame(r_m) && nrow(r_m) > 0 && "gid" %in% names(r_m)) {
                    r_m$gid <- m_sf$gid[as.integer(r_m$gid)]
                  }
                  r_m
                }, error = function(e_m) {
                  pt_list <- list()
                  for (p_i in seq_along(idx_vec)) {
                    r_p <- tryCatch(
                      soilDB::SDA_spatialQuery(m_sf[p_i, , drop = FALSE], what = "mukey", db = "SSURGO", byFeature = FALSE),
                      error = function(e_p) NULL
                    )
                    if (is.data.frame(r_p) && nrow(r_p) > 0 && "mukey" %in% names(r_p)) {
                      r_p$gid <- m_sf$gid[p_i]
                      pt_list[[length(pt_list) + 1L]] <- r_p
                    }
                  }
                  if (length(pt_list) > 0) dplyr::bind_rows(pt_list) else NULL
                })
                if (is.data.frame(m_res) && nrow(m_res) > 0) {
                  mini_list[[length(mini_list) + 1L]] <- m_res
                }
              }
              if (length(mini_list) > 0) dplyr::bind_rows(mini_list) else NULL
            }
          )
          if (!is.null(sp_res) && is.data.frame(sp_res) && nrow(sp_res) > 0 && "mukey" %in% names(sp_res)) {
            sp_keep <- intersect(c("gid", "mukey"), names(sp_res))
            mukeys <- unique(na.omit(as.character(sp_res$mukey)))
            if (length(mukeys) > 0) {
              mukey_str <- paste0("'", mukeys, "'", collapse = ",")
              sql <- paste0(
                "SELECT m.mukey, m.musym, m.muname, ",
                "ma.brockdepmin, ma.aws050wta, ",
                "c.cokey, c.compname, c.comppct_r, c.compkind, c.majcompflag, ",
                "c.taxsuborder, c.taxgrtgroup AS taxgreatgroup, ",
                "cr.resdept_r, cr.reskind, ",
                "ch.hzdept_r, ch.hzdepb_r, ch.sandtotal_r, ch.claytotal_r, ch.ph1to1h2o_r, ch.om_r, ch.cec7_r ",
                "FROM mapunit m ",
                "LEFT JOIN muaggatt ma ON m.mukey = ma.mukey ",
                "LEFT JOIN component c ON m.mukey = c.mukey ",
                "LEFT JOIN corestrictions cr ON c.cokey = cr.cokey ",
                "LEFT JOIN chorizon ch ON c.cokey = ch.cokey AND ch.hzdept_r < 20 ",
                sprintf("WHERE m.mukey IN (%s)", mukey_str)
              )
              tab_res <- soilDB::SDA_query(sql)
              if (is.data.frame(tab_res) && nrow(tab_res) > 0) {
                tab_joined <- merge(sp_res[, sp_keep, drop = FALSE], tab_res, by = "mukey", all.x = TRUE)
                if (exists("aggregate_ssurgo_tabular_to_points", mode = "function")) {
                  agg_df <- aggregate_ssurgo_tabular_to_points(sub_pts, tab_joined, classify_fn = classify_ssurgo_outcrop)
                  chunk_df <- normalize_ssurgo_cols(agg_df)
                } else {
                  pts_sp <- merge(sub_pts, sp_res[!duplicated(sp_res$gid), sp_keep, drop = FALSE], by = "gid", all.x = TRUE)
                  dom_comp <- tab_res[!duplicated(tab_res$mukey), ]
                  res_joined <- merge(pts_sp, dom_comp, by = "mukey", all.x = TRUE)
                  chunk_df <- normalize_ssurgo_cols(res_joined)
                }
                chunk_ok <- TRUE
                break
              }
            }
          }
        }, error = function(e) {
          message(sprintf("Batch %d chunk %d attempt %d failed: %s", batch_idx, c_idx, attempt, e$message))
        })
        if (attempt < max_retries && !chunk_ok) Sys.sleep(delays[attempt])
      }

      if (chunk_ok && !is.null(chunk_df)) {
        chunk_results[[c_idx]] <- chunk_df
      } else {
        all_chunks_ok <- FALSE
        break
      }
    }

    if (all_chunks_ok && length(chunk_results) > 0) {
      res_df <- normalize_ssurgo_cols(dplyr::bind_rows(chunk_results))
      query_success <- TRUE
    }
  }
  
  if (query_success && !is.null(res_df)) {
    if (!isTRUE(demo)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      tmp_csv <- paste0(cache_file_csv, ".tmp")
      readr::write_csv(res_df, tmp_csv)
      file.rename(tmp_csv, cache_file_csv)
      tmp_rds <- paste0(cache_file_rds, ".tmp")
      saveRDS(res_df, tmp_rds)
      file.rename(tmp_rds, cache_file_rds)
    }
    return(res_df)
  }

  # Never cache synthetic or failed SSURGO records into data/environmental/ssurgo_cache/ via saveRDS()
  if (!isTRUE(demo)) {
    warning(sprintf(
      paste0(
        "soilDB::SDA_spatialQuery() failed or timed out after %d retries for batch %d and no cached SSURGO records exist in '%s'. ",
        "Retaining already-extracted SoilGrids 250m + WorldClim v2.1 raster covariates with NA SSURGO vector attributes."
      ),
      max_retries, batch_idx, cache_dir
    ))
    fallback_df <- data.frame(
      catalogNumber = pts_df$catalogNumber,
      latitude = pts_df$latitude,
      longitude = pts_df$longitude,
      musym = rep(NA_character_, n_pts),
      muname = rep(NA_character_, n_pts),
      taxsuborder = rep(NA_character_, n_pts),
      taxgreatgroup = rep(NA_character_, n_pts),
      is_rock_outcrop = rep(FALSE, n_pts),
      outcrop_class = rep("Non-outcrop matrix", n_pts),
      stringsAsFactors = FALSE
    )
    return(normalize_ssurgo_cols(fallback_df))
  }

  message(sprintf("[DEMO MODE] Generating synthetic SSURGO assignments for batch %d (not cached to disk)...", batch_idx))
  musyms <- character(n_pts); munames <- character(n_pts)
  suborders <- character(n_pts); greatgroups <- character(n_pts)
  is_outcrops <- logical(n_pts); outcrop_classes <- character(n_pts)
  
  for (i in seq_len(n_pts)) {
    lat <- pts_df$latitude[i]; lon <- pts_df$longitude[i]
    reg <- if ("regional_group" %in% names(batch_df)) batch_df$regional_group[i] else "Unknown"
    
    is_flatrock_zone <- grepl("Flatrock", reg, ignore.case = TRUE) ||
                        (lon > -84.5 && lon < -78.5 && lat > 33.2 && lat < 36.8)
    is_sandhill_zone <- grepl("Sandhill|Coastal", reg, ignore.case = TRUE) ||
                        (lon > -82.0 && lon < -75.5 && lat < 36.5)
    is_prairie_zone  <- grepl("Midwest|Prairie", reg, ignore.case = TRUE) || (lon < -88.0)
    is_appalachian   <- grepl("Appalachian", reg, ignore.case = TRUE) ||
                        (lon >= -84.0 && lon <= -79.0 && lat >= 35.0 && lat <= 39.0)
    
    if (is_flatrock_zone) {
      if (as.integer(abs(sin(lat * 100 + lon * 100)) * 1000) %% 10 < 7) {
        musyms[i] <- "RoB"; munames[i] <- "Rock outcrop-Wake complex, 2 to 10 percent slopes"; suborders[i] <- "Udepts"; greatgroups[i] <- "Lithic Dystrudepts"
      } else {
        musyms[i] <- "ApB"; munames[i] <- "Appling sandy loam, 2 to 6 percent slopes"; suborders[i] <- "Udults"; greatgroups[i] <- "Typic Kanhapludults"
      }
    } else if (is_sandhill_zone) {
      musyms[i] <- "WaB"; munames[i] <- "Wagram sand, 0 to 6 percent slopes"; suborders[i] <- "Udults"; greatgroups[i] <- "Arenic Kandiudults"
    } else if (is_prairie_zone) {
      musyms[i] <- "TaA"; munames[i] <- "Tama silt loam, 0 to 2 percent slopes"; suborders[i] <- "Udolls"; greatgroups[i] <- "Typic Argiudolls"
    } else if (is_appalachian) {
      if (as.integer(abs(cos(lat * 50 + lon * 50)) * 1000) %% 10 < 3) {
        musyms[i] <- "GvC"; munames[i] <- "Gladeville-Rock outcrop complex, 2 to 12 percent slopes"; suborders[i] <- "Udolls"; greatgroups[i] <- "Lithic Hapludolls"
      } else {
        musyms[i] <- "CeB"; munames[i] <- "Cecil sandy clay loam, 2 to 8 percent slopes"; suborders[i] <- "Udults"; greatgroups[i] <- "Typic Hapludults"
      }
    } else {
      musyms[i] <- "MuB"; munames[i] <- "Generic Upland loam, 1 to 5 percent slopes"; suborders[i] <- "Udults"; greatgroups[i] <- "Typic Hapludults"
    }
    
    cl <- classify_ssurgo_outcrop(munames[i], "", suborders[i], greatgroups[i], "", lat, lon)
    is_outcrops[i] <- cl$is_rock_outcrop
    outcrop_classes[i] <- cl$outcrop_class
  }
  
  res_df <- data.frame(
    catalogNumber = pts_df$catalogNumber,
    latitude = pts_df$latitude,
    longitude = pts_df$longitude,
    musym = musyms,
    muname = munames,
    taxsuborder = suborders,
    taxgreatgroup = greatgroups,
    is_rock_outcrop = is_outcrops,
    outcrop_class = outcrop_classes,
    stringsAsFactors = FALSE
  )
  return(normalize_ssurgo_cols(res_df))
}

query_ssurgo_sda <- function(df, cache_dir = "data/environmental/ssurgo_cache", batch_size = 75, demo = FALSE) {
  message(sprintf("Initiating USDA SSURGO Vector Pedology SDA Query (cache: %s, batch_size: %d)...", cache_dir, batch_size))
  
  is_conus <- !is.na(df$latitude) & !is.na(df$longitude) &
              df$latitude >= 24.5 & df$latitude <= 49.5 &
              df$longitude >= -125.0 & df$longitude <= -66.5
  conus_df <- df[is_conus, ]
  n_conus <- nrow(conus_df)
  message(sprintf("Found %d vouchers within CONUS spatial domain for SSURGO vector evaluation.", n_conus))
  
  ssurgo_quant_cols <- c(
    "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm",
    "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph",
    "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm"
  )

  if (n_conus == 0) {
    df$musym <- NA_character_; df$muname <- NA_character_
    df$taxsuborder <- NA_character_; df$taxgreatgroup <- NA_character_
    df$is_rock_outcrop <- FALSE; df$outcrop_class <- "Non-outcrop matrix"
    for (sc in ssurgo_quant_cols) df[[sc]] <- NA_real_
    return(df)
  }
  
  merged_df <- tryCatch({
    batch_indices <- split(seq_len(n_conus), ceiling(seq_len(n_conus) / batch_size))
    batch_results <- list()
    
    for (b_idx in seq_along(batch_indices)) {
      idx <- batch_indices[[b_idx]]
      sub_df <- conus_df[idx, ]
      b_res <- query_ssurgo_sda_batch(sub_df, cache_dir, b_idx, demo = demo, batch_size = batch_size)
      batch_results[[b_idx]] <- b_res
    }
    
    all_ssurgo <- as.data.frame(dplyr::bind_rows(batch_results), stringsAsFactors = FALSE)
    keep_cols <- c("catalogNumber", "musym", "muname", "taxsuborder", "taxgreatgroup",
                   "is_rock_outcrop", "outcrop_class", ssurgo_quant_cols)
    all_ssurgo <- all_ssurgo[!duplicated(all_ssurgo$catalogNumber), intersect(names(all_ssurgo), keep_cols), drop = FALSE]
    
    df$catalogNumber <- as.character(df$catalogNumber)
    all_ssurgo$catalogNumber <- as.character(all_ssurgo$catalogNumber)
    for (col_rm in intersect(setdiff(keep_cols, "catalogNumber"), names(df))) {
      df[[col_rm]] <- NULL
    }
    m_df <- dplyr::left_join(df, all_ssurgo, by = "catalogNumber")
    m_df$is_rock_outcrop[is.na(m_df$is_rock_outcrop)] <- FALSE
    m_df$outcrop_class[is.na(m_df$outcrop_class)] <- "Non-outcrop matrix"
    m_df
  }, error = function(e) {
    warning(sprintf(
      "USDA SSURGO SDA query failed (%s); retaining SoilGrids 250m + WorldClim v2.1 raster covariates.",
      e$message
    ))
    df$musym <- NA_character_; df$muname <- NA_character_
    df$taxsuborder <- NA_character_; df$taxgreatgroup <- NA_character_
    df$is_rock_outcrop <- FALSE; df$outcrop_class <- "Non-outcrop matrix"
    for (sc in ssurgo_quant_cols) df[[sc]] <- NA_real_
    df
  })
  
  return(merged_df)
}

# ------------------------------------------------------------------------------
# 2c. Multi-Scale Edaphic Synthesis & Outcrop Fidelity Contingency Testing
# ------------------------------------------------------------------------------
synthesize_multiscale_edaphics <- function(df, ssurgo_out = "data/tables/ssurgo_edaphic_validation.csv",
                                          contingency_out = "outputs/reports/micro_edaphic_outcrop_contingency.csv",
                                          demo = FALSE) {
  message("Executing Multi-Scale Edaphic Synthesis & Fisher's Exact Outcrop Fidelity Test...")
  
  val_cols <- c("catalogNumber", "lat_lon_key", "latitude", "longitude", "species_standardized", "regional_group",
                "coordinateUncertaintyInMeters", "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
                "musym", "muname", "taxsuborder", "taxgreatgroup", "is_rock_outcrop", "outcrop_class",
                "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm",
                "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm",
                "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m",
                "solus_resdept_cm_100m", "solus_fragvol_pct_100m",
                "topo_elev_m", "topo_tpi_30m", "topo_roughness_30m", "topo_slope_deg", "topo_twi_30m",
                "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm", "fine_rock_outcrop_prob",
                "bdricm", "resdept_r", "cfvo", "awc", "fine_soil_cec", "soil_silt_pct", "alr_sand_clay", "alr_silt_clay",
                "sgmc_lithology_class", "is_azonal_lithology", "azonal_contact_dist_m", "azonal_proximity_index",
                "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer",
                "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")
  present_cols <- intersect(val_cols, names(df))
  val_df <- df[, present_cols, drop = FALSE]
  names(val_df)[names(val_df) %in% c("soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")] <-
    paste0(names(val_df)[names(val_df) %in% c("soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")], "_soilgrids")
  
  if (write_csv_atomic(val_df, ssurgo_out, demo = demo)) {
    message(sprintf("Exported specimen-level SSURGO & 10m-30m fine-grain edaphic validation table to %s (%d records).", ssurgo_out, nrow(val_df)))
  }
  
  active_tx <- get_active_packera_taxa(df$species_standardized)
  grp_col_use <- "species_standardized"
  grp_levels <- active_tx
  if (length(active_tx) < 2L && "analysis_group" %in% names(df)) {
    ag_u <- sort(unique(as.character(stats::na.omit(df$analysis_group))))
    ag_u <- ag_u[nzchar(ag_u) & ag_u != "Unknown"]
    if (length(ag_u) >= 2L) {
      grp_col_use <- "analysis_group"
      grp_levels <- ag_u
    }
  }
  if (!("is_rock_outcrop" %in% names(df))) {
    df$is_rock_outcrop <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) {
      !is.na(df$ssurgo_rock_outcrop_pct) & as.numeric(df$ssurgo_rock_outcrop_pct) >= 15.0
    } else {
      rep(FALSE, nrow(df))
    }
  }
  df$is_rock_outcrop[is.na(df$is_rock_outcrop)] <- FALSE

  mc_res <- evaluate_multiclass_substrate_contingency(
    df,
    active_groups = grp_levels,
    group_col = grp_col_use
  )
  contingency_report <- mc_res$summary_table
  fisher_res <- mc_res$omnibus_test

  if (nrow(contingency_report) > 0L && write_csv_atomic(contingency_report, contingency_out, demo = demo)) {
    message(sprintf("Exported multi-class lithology & outcrop Fisher's exact contingency report to %s (%d rows).",
                    contingency_out, nrow(contingency_report)))
  }

  return(list(
    validation = val_df,
    contingency = contingency_report,
    multiclass_contingency = mc_res$multiclass_details,
    fisher = fisher_res
  ))
}

# ------------------------------------------------------------------------------
# 2d. Latitudinal Spring Baseline & Phenological Anomalies (Delta DOY) Engine
# ------------------------------------------------------------------------------
extract_doy_robust <- function(df) {
  doy_vec <- if ("doy" %in% names(df)) suppressWarnings(as.numeric(df$doy)) else rep(NA_real_, nrow(df))
  if ("eventDate" %in% names(df)) {
    na_idx <- which(is.na(doy_vec) & !is.na(df$eventDate))
    if (length(na_idx) > 0) {
      dt <- as.Date(df$eventDate[na_idx], format = "%Y-%m-%d")
      doy_vec[na_idx] <- as.integer(format(dt, "%j"))
    }
  }
  if (any(is.na(doy_vec)) && all(c("year", "month", "day") %in% names(df))) {
    na_idx <- which(is.na(doy_vec) & !is.na(df$year) & !is.na(df$month) & !is.na(df$day))
    if (length(na_idx) > 0) {
      y_str <- as.character(as.integer(df$year[na_idx]))
      m_str <- sprintf("%02d", as.integer(df$month[na_idx]))
      d_str <- sprintf("%02d", as.integer(df$day[na_idx]))
      dt <- as.Date(paste(y_str, m_str, d_str, sep = "-"), format = "%Y-%m-%d")
      doy_vec[na_idx] <- as.integer(format(dt, "%j"))
    }
  }
  return(doy_vec)
}

model_latitudinal_spring_baseline <- function(flowering_df) {
  message("Modeling Latitudinal Spring Baseline (Hopkins' Bioclimatic Law, Hopkins 1918) on verified anthesis cohort...")
  
  lat_col <- if ("decimalLatitude" %in% names(flowering_df) && any(!is.na(flowering_df$decimalLatitude))) {
    flowering_df$decimalLatitude
  } else {
    flowering_df$latitude
  }
  flowering_df$decimalLatitude <- lat_col
  
  valid <- flowering_df[!is.na(flowering_df$decimalLatitude) & !is.na(flowering_df$doy), ]
  valid <- valid[valid$decimalLatitude >= 24.0 & valid$decimalLatitude <= 55.0 & valid$doy >= 50 & valid$doy <= 220, ]
  
  # Trim extreme 1st and 99th percentiles to guard against outlier collection dates
  p1 <- stats::quantile(valid$doy, 0.01, na.rm = TRUE)
  p99 <- stats::quantile(valid$doy, 0.99, na.rm = TRUE)
  valid_trimmed <- valid[valid$doy >= p1 & valid$doy <= p99, ]
  
  fit <- if (requireNamespace("MASS", quietly = TRUE)) {
    tryCatch(MASS::rlm(doy ~ decimalLatitude, data = valid_trimmed), error = function(e) stats::lm(doy ~ decimalLatitude, data = valid_trimmed))
  } else {
    stats::lm(doy ~ decimalLatitude, data = valid_trimmed)
  }
  coefs <- stats::coef(fit)
  b0 <- as.numeric(coefs[1]); b_lat <- as.numeric(coefs[2])
  min_lat <- min(valid_trimmed$decimalLatitude, na.rm = TRUE)
  pred_min_lat <- b0 + b_lat * min_lat
  decimalLatitude <- valid_trimmed$decimalLatitude
  hopkins_expected <- pred_min_lat + 4.0 * (decimalLatitude - min_lat)
  valid_trimmed$hopkins_expected <- hopkins_expected
  
  # Multivariate bioclimatic adjustment for latitude, elevation, and MAT (Hopkins 1918) when covariates are available
  elev_col <- if ("topo_elev_m" %in% names(valid_trimmed)) suppressWarnings(as.numeric(valid_trimmed$topo_elev_m)) else if ("elevation" %in% names(valid_trimmed)) suppressWarnings(as.numeric(valid_trimmed$elevation)) else rep(NA_real_, nrow(valid_trimmed))
  mat_col <- if ("bio1_temp_mean" %in% names(valid_trimmed)) suppressWarnings(as.numeric(valid_trimmed$bio1_temp_mean)) else rep(NA_real_, nrow(valid_trimmed))
  b_elev <- 0.0; ref_elev <- 0.0; b_mat <- 0.0; ref_mat <- 0.0
  if (sum(is.finite(elev_col)) >= 10 && stats::sd(elev_col, na.rm = TRUE) > 5.0) {
    ref_elev <- stats::median(elev_col, na.rm = TRUE)
    b_elev <- 0.0125 # Hopkins' Bioclimatic Law: +1.25 days per +100 m elevation (+0.0125 d/m)
  }
  if (sum(is.finite(mat_col)) >= 10 && stats::sd(mat_col, na.rm = TRUE) > 0.5) {
    ref_mat <- stats::median(mat_col, na.rm = TRUE)
    fit_mat <- tryCatch(stats::lm(I(doy - (b0 + b_lat * decimalLatitude)) ~ mat_col, data = valid_trimmed), error = function(e) NULL)
    if (!is.null(fit_mat) && length(stats::coef(fit_mat)) >= 2 && is.finite(stats::coef(fit_mat)[2])) {
      b_mat <- max(-2.5, min(2.5, as.numeric(stats::coef(fit_mat)[2])))
    }
  }
  
  fit_ols <- stats::lm(doy ~ decimalLatitude, data = valid_trimmed)
  s_ols <- summary(fit_ols)
  r2 <- s_ols$r.squared; p_val <- s_ols$coefficients[2, 4]
  
  message(sprintf("Empirical Baseline (Hopkins 1918): DOY = %.3f + %.3f * decimalLatitude (b_elev = %.4f d/m, b_mat = %.3f d/C, R2 = %.4f, p = %.4e, N = %d)",
                  b0, b_lat, b_elev, b_mat, r2, p_val, nrow(valid_trimmed)))
  return(list(intercept = b0, slope = b_lat, b_elev = b_elev, ref_elev = ref_elev, b_mat = b_mat, ref_mat = ref_mat,
              min_lat = min_lat, pred_min_lat = pred_min_lat,
              hopkins_expected = hopkins_expected, r_squared = r2, p_value = p_val,
              n = nrow(valid_trimmed), model = fit_ols))
}

compute_phenological_anomalies <- function(df, output_csv = NULL, summary_csv = NULL, demo = FALSE) {
  message("Computing latitude-adjusted phenological anomalies (Delta DOY)...")
  df$doy <- extract_doy_robust(df)
  lat_val <- if ("decimalLatitude" %in% names(df) && any(!is.na(df$decimalLatitude))) df$decimalLatitude else df$latitude
  df$decimalLatitude <- lat_val
  df$latitude <- lat_val

  # ----------------------------------------------------------------------------
  # 3-Tier Anthesis Verification Filter
  # ----------------------------------------------------------------------------
  # Tier 1: Calendar gate: 50 <= DOY <= 220 and 24.0 <= Latitude <= 55.0
  pass_calendar <- !is.na(df$doy) & df$doy >= 50 & df$doy <= 220 &
                   !is.na(df$decimalLatitude) & df$decimalLatitude >= 24.0 & df$decimalLatitude <= 55.0

  # Tier 2: Darwin Core check: if reproductiveCondition contains "fruit" or "sterile" without "fl", exclude
  if ("reproductiveCondition" %in% names(df)) {
    rc <- tolower(trimws(as.character(df$reproductiveCondition)))
    rc[is.na(rc)] <- ""
    has_fruit <- grepl("fruit", rc)
    has_sterile <- grepl("sterile", rc)
    has_fl <- grepl("fl", rc)
    dwc_fruit_only <- has_fruit & !has_fl
    dwc_sterile_only <- has_sterile & !has_fl
  } else {
    dwc_fruit_only <- rep(FALSE, nrow(df))
    dwc_sterile_only <- rep(FALSE, nrow(df))
  }
  pass_dwc <- !(dwc_fruit_only | dwc_sterile_only)

  # Tier 3: Computer Vision check: voucher_phenological_state == "anthesis"
  # Graceful fallback: unreadable crops or missing CV records fall back to Calendar + DwC text
  if ("voucher_phenological_state" %in% names(df)) {
    v_state <- tolower(trimws(as.character(df$voucher_phenological_state)))
    v_state[is.na(v_state)] <- ""
  } else {
    v_state <- rep("", nrow(df))
  }

  if ("total_capitula" %in% names(df)) {
    tot_caps <- suppressWarnings(as.numeric(df$total_capitula))
    tot_caps[is.na(tot_caps)] <- 0
  } else {
    tot_caps <- rep(0, nrow(df))
  }

  # Valid CV classification when crops were present and state recorded
  has_valid_cv <- v_state != "" & v_state != "unknown" & v_state != "na" & tot_caps > 0

  cv_anthesis <- rep(FALSE, nrow(df))
  cv_fruit <- rep(FALSE, nrow(df))
  cv_sterile <- rep(FALSE, nrow(df))

  cv_anthesis[has_valid_cv & v_state == "anthesis"] <- TRUE
  cv_fruit[has_valid_cv & v_state == "fruit"] <- TRUE
  cv_sterile[has_valid_cv & v_state %in% c("sterile", "bud")] <- TRUE

  # Fallback for missing CV / no readable crops: accept if calendar and DwC pass
  missing_cv <- !has_valid_cv
  cv_anthesis[missing_cv & pass_calendar & pass_dwc] <- TRUE
  cv_fruit[missing_cv & dwc_fruit_only] <- TRUE
  cv_sterile[missing_cv & dwc_sterile_only] <- TRUE

  # Synthesize 3-tier status
  is_verified_anthesis <- pass_calendar & pass_dwc & cv_anthesis
  is_fruit_excluded <- pass_calendar & (cv_fruit | dwc_fruit_only) & !is_verified_anthesis
  is_sterile_excluded <- pass_calendar & (cv_sterile | dwc_sterile_only) & !is_verified_anthesis & !is_fruit_excluded

  df$is_verified_anthesis <- is_verified_anthesis
  df$anthesis_filter_status <- "Excluded_Calendar_Gate"
  df$anthesis_filter_status[pass_calendar] <- "Excluded_Other"
  df$anthesis_filter_status[is_sterile_excluded] <- "Excluded_Sterile"
  df$anthesis_filter_status[is_fruit_excluded] <- "Excluded_Fruit_Only"
  df$anthesis_filter_status[is_verified_anthesis] <- "Verified_Anthesis"

  # Telemetry Logging
  n_verified <- sum(is_verified_anthesis)
  n_fruit <- sum(is_fruit_excluded)
  n_sterile <- sum(is_sterile_excluded)
  message(sprintf("[INFO] Filtered phenology vouchers: %d verified anthesis, %d fruiting-only (excluded), %d sterile (excluded).",
                  n_verified, n_fruit, n_sterile))

  # Fit latitudinal cline strictly on verified anthesis cohort
  anthesis_cohort <- df[df$is_verified_anthesis, ]
  if (nrow(anthesis_cohort) < 5) {
    message("[WARN] Verified anthesis cohort has < 5 vouchers; using calendar-gated cohort as baseline fallback.")
    anthesis_cohort <- df[pass_calendar, ]
  }
  baseline <- model_latitudinal_spring_baseline(anthesis_cohort)
  
  df$expected_doy <- NA_real_
  df$hopkins_expected <- NA_real_
  has_lat <- !is.na(df$decimalLatitude)
  elev_v <- if ("topo_elev_m" %in% names(df)) suppressWarnings(as.numeric(df$topo_elev_m)) else if ("elevation" %in% names(df)) suppressWarnings(as.numeric(df$elevation)) else rep(NA_real_, nrow(df))
  mat_v <- if ("bio1_temp_mean" %in% names(df)) suppressWarnings(as.numeric(df$bio1_temp_mean)) else rep(NA_real_, nrow(df))
  elev_adj <- ifelse(is.finite(elev_v) & !is.null(baseline$b_elev) & baseline$b_elev != 0, baseline$b_elev * (elev_v - baseline$ref_elev), 0.0)
  mat_adj <- ifelse(is.finite(mat_v) & !is.null(baseline$b_mat) & baseline$b_mat != 0, baseline$b_mat * (mat_v - baseline$ref_mat), 0.0)
  df$expected_doy[has_lat] <- baseline$intercept + baseline$slope * df$decimalLatitude[has_lat] + elev_adj[has_lat] + mat_adj[has_lat]
  df$hopkins_expected[has_lat] <- baseline$pred_min_lat + 4.0 * (df$decimalLatitude[has_lat] - baseline$min_lat) + elev_adj[has_lat]
  
  df$delta_doy <- NA_real_
  has_both <- !is.na(df$doy) & !is.na(df$expected_doy)
  df$delta_doy[has_both] <- df$doy[has_both] - df$expected_doy[has_both]
  # Hopkins-adjusted flowering DOY & circular coordinates (sin_doy, cos_doy) standardized to reference latitude/elevation/MAT (Hopkins 1918)
  ref_doy_base <- baseline$intercept + baseline$slope * stats::median(df$decimalLatitude[has_lat], na.rm = TRUE)
  df$hopkins_adjusted_doy <- ifelse(has_both, ref_doy_base + df$delta_doy, df$doy)
  valid_hdoy <- !is.na(df$hopkins_adjusted_doy) & is.finite(df$hopkins_adjusted_doy)
  df$sin_doy <- ifelse(valid_hdoy, round(sin(2 * pi * df$hopkins_adjusted_doy / 365.25), 6), NA_real_)
  df$cos_doy <- ifelse(valid_hdoy, round(cos(2 * pi * df$hopkins_adjusted_doy / 365.25), 6), NA_real_)
  
  df$pheno_timing_category <- "Unknown"
  df$pheno_timing_category[has_both & df$delta_doy < -7.0] <- "Early_Flowering"
  df$pheno_timing_category[has_both & df$delta_doy > 7.0] <- "Late_Flowering"
  df$pheno_timing_category[has_both & abs(df$delta_doy) <= 7.0] <- "Synchronous_with_Cline"
  
  if (!is.null(output_csv)) {
    keep_cols <- c("catalogNumber", "institutionCode", "scientificName", "species_raw",
                   "species_standardized", "decimalLatitude", "decimalLongitude", "latitude", "longitude",
                   "eventDate", "year", "month", "day", "doy", "expected_doy", "delta_doy",
                   "pheno_timing_category", "is_verified_anthesis", "anthesis_filter_status",
                   "voucher_phenological_state", "regional_group")
    export_df <- df[, intersect(keep_cols, names(df))]
    if (write_csv_atomic(export_df, output_csv, demo = demo)) {
      message(sprintf("Exported %d voucher phenological anomalies to %s.", nrow(export_df), output_csv))
    }
  }
  
  # Statistical divergence testing across target taxa strictly on verified anthesis cohort
  active_tx <- get_active_packera_taxa(df$species_standardized)
  target_flowering <- df[df$species_standardized %in% active_tx &
                         !is.na(df$delta_doy) & df$is_verified_anthesis, ]
  
  if (nrow(target_flowering) >= 5 && length(unique(target_flowering$species_standardized)) > 1) {
    aov_fit <- stats::aov(delta_doy ~ species_standardized, data = target_flowering)
    aov_summary <- summary(aov_fit)[[1]]
    f_stat <- aov_summary[["F value"]][1]; p_val_aov <- aov_summary[["Pr(>F)"]][1]
    df_b <- aov_summary[["Df"]][1]; df_w <- aov_summary[["Df"]][2]
    tukey_res <- stats::TukeyHSD(aov_fit)
    tukey_mat <- tukey_res$species_standardized
  } else {
    aov_fit <- NULL
    tukey_res <- NULL
    tukey_mat <- matrix(numeric(0), nrow = 0, ncol = 4)
    f_stat <- NA_real_; p_val_aov <- NA_real_; df_b <- NA_real_; df_w <- NA_real_
  }
  
  summary_records <- list()
  summary_records[[1]] <- data.frame(
    Analysis_Section = "Baseline_Regression", Taxon_or_Comparison = "DOY ~ decimalLatitude (Verified Anthesis)",
    Sample_Size = baseline$n, Mean_Delta_DOY = round(baseline$slope, 3), SD_Delta_DOY = round(baseline$intercept, 3),
    Median_Delta_DOY = round(baseline$r_squared, 4), IQR_Delta_DOY = NA_real_, CI_Lower_95 = NA_real_, CI_Upper_95 = NA_real_,
    P_Value = sprintf("%.4e", baseline$p_value),
    Interpretation = sprintf("Empirical spring advance slope = +%.2f days/degree N (Hopkins' Law)", baseline$slope),
    stringsAsFactors = FALSE
  )
  summary_records[[2]] <- data.frame(
    Analysis_Section = "One_Way_ANOVA", Taxon_or_Comparison = "Across Target Taxa (Verified Anthesis)",
    Sample_Size = nrow(target_flowering), Mean_Delta_DOY = round(f_stat, 3), SD_Delta_DOY = df_b,
    Median_Delta_DOY = df_w, IQR_Delta_DOY = NA_real_, CI_Lower_95 = NA_real_, CI_Upper_95 = NA_real_,
    P_Value = if (!is.na(p_val_aov)) sprintf("%.4e", p_val_aov) else "NA",
    Interpretation = if (!is.na(p_val_aov) && p_val_aov < 0.05) "Significant phenological divergence between taxa" else "No significant divergence or insufficient sample size",
    stringsAsFactors = FALSE
  )
  for (tx in active_tx) {
    sub_vals <- target_flowering$delta_doy[target_flowering$species_standardized == tx]
    if (length(sub_vals) > 0) {
      m <- mean(sub_vals); s <- sd(sub_vals); med <- median(sub_vals); iq <- stats::IQR(sub_vals)
      n_tx <- length(sub_vals); se <- s / sqrt(n_tx)
      summary_records[[length(summary_records) + 1]] <- data.frame(
        Analysis_Section = "Taxon_Distribution", Taxon_or_Comparison = tx, Sample_Size = n_tx,
        Mean_Delta_DOY = round(m, 3), SD_Delta_DOY = round(s, 3), Median_Delta_DOY = round(med, 3),
        IQR_Delta_DOY = round(iq, 3), CI_Lower_95 = round(m - 1.96 * se, 3), CI_Upper_95 = round(m + 1.96 * se, 3),
        P_Value = "",
        Interpretation = if (m < -3) "Early blooming relative to latitude" else if (m > 3) "Late blooming relative to latitude" else "Conforms to latitudinal cline",
        stringsAsFactors = FALSE
      )
    }
  }
  if (nrow(tukey_mat) > 0) {
    for (k in seq_len(nrow(tukey_mat))) {
      comp_name <- rownames(tukey_mat)[k]
      comp_parts <- strsplit(comp_name, "-", fixed = TRUE)[[1]]
      diff_val <- tukey_mat[k, "diff"]; lwr <- tukey_mat[k, "lwr"]; upr <- tukey_mat[k, "upr"]; p_adj <- tukey_mat[k, "p adj"]
      if (length(comp_parts) == 2L && comp_parts[1] > comp_parts[2]) {
        comp_label <- paste(comp_parts[2], "vs", comp_parts[1])
        diff_val <- -diff_val
        tmp_lwr <- -upr; upr <- -lwr; lwr <- tmp_lwr
      } else {
        comp_label <- gsub("-", " vs ", comp_name)
      }
      summary_records[[length(summary_records) + 1]] <- data.frame(
        Analysis_Section = "Tukey_HSD_PostHoc", Taxon_or_Comparison = comp_label,
        Sample_Size = nrow(target_flowering), Mean_Delta_DOY = round(diff_val, 3), SD_Delta_DOY = NA_real_,
        Median_Delta_DOY = NA_real_, IQR_Delta_DOY = NA_real_, CI_Lower_95 = round(lwr, 3), CI_Upper_95 = round(upr, 3),
        P_Value = sprintf("%.4e", p_adj),
        Interpretation = if (p_adj < 0.05) sprintf("Statistically significant temporal separation (%.1f days)", abs(diff_val)) else "Phenologically synchronous / overlapping",
        stringsAsFactors = FALSE
      )
    }
  }
  summary_report <- do.call(rbind, summary_records)
  if (!is.null(summary_csv)) {
    if (write_csv_atomic(summary_report, summary_csv, demo = demo)) {
      message(sprintf("Exported phenological anomaly ANOVA/Tukey summary to %s.", summary_csv))
    }
  }
  return(list(df = df, baseline = baseline, anova = aov_fit, tukey = tukey_res, summary = summary_report))
}

export_phenological_latitudinal_figures <- function(df, baseline, out_pdf) {
  message(sprintf("Generating 2-panel phenological latitudinal anomaly figure to %s...", out_pdf))
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  df$plot_group <- if (exists("resolve_plot_group_column", mode = "function")) {
    resolve_plot_group_column(df)
  } else {
    sp_v <- as.character(df$species_standardized)
    if (length(unique(stats::na.omit(sp_v))) < 2L ||
        ("analysis_group" %in% names(df) && any(grepl("Morphotype_K|^LatTertile_|Ecotype|Subcluster", df$analysis_group, ignore.case = TRUE)))) {
      if ("analysis_group" %in% names(df) && length(unique(stats::na.omit(df$analysis_group))) >= 2L) as.character(df$analysis_group)
      else if ("grouping_col" %in% names(df) && length(unique(stats::na.omit(df$grouping_col))) >= 2L) as.character(df$grouping_col)
      else if ("gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster))) paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster)))
      else sp_v
    } else sp_v
  }
  active_tx <- sort(unique(stats::na.omit(df$plot_group)))
  active_tx <- active_tx[nzchar(active_tx) & active_tx != "Unknown"]
  if (length(active_tx) == 0) active_tx <- get_active_packera_taxa(df$species_standardized)
  dyn_pal <- resolve_dynamic_group_palette(active_tx)
  dyn_shapes <- stats::setNames(rep_len(c(21, 24, 22, 23, 1, 2, 5, 0), length(active_tx)), active_tx)
  is_target <- df$plot_group %in% active_tx & !is.na(df$delta_doy)
  target_df <- if ("is_verified_anthesis" %in% names(df)) {
    df[is_target & df$is_verified_anthesis, ]
  } else {
    df[is_target & df$doy >= 50 & df$doy <= 220, ]
  }
  if (nrow(target_df) == 0) target_df <- df[is_target, ]
  
  thm <- ggplot2::theme_bw(base_size = 9.5) +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                        plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey25"),
                        panel.grid.minor = ggplot2::element_blank())
  
  lat_col <- if ("decimalLatitude" %in% names(target_df) && any(!is.na(target_df$decimalLatitude))) target_df$decimalLatitude else target_df$latitude
  target_df$decimalLatitude <- lat_col
  min_lat <- min(target_df$decimalLatitude, na.rm = TRUE)
  pred_min_lat <- baseline$intercept + baseline$slope * min_lat
  lat_seq <- seq(min_lat, max(target_df$decimalLatitude, na.rm = TRUE), length.out = 100)
  pred_line <- data.frame(
    latitude = lat_seq,
    decimalLatitude = lat_seq,
    doy = baseline$intercept + baseline$slope * lat_seq,
    hopkins_expected = pred_min_lat + 4.0 * (lat_seq - min_lat)
  )
  
  pt_alpha <- if (nrow(target_df) > 1000) 0.45 else 0.75
  p1 <- ggplot2::ggplot(target_df, ggplot2::aes(x = decimalLatitude, y = doy, color = plot_group, shape = plot_group)) +
    ggplot2::geom_point(fill = NA, stroke = 0.78, alpha = pt_alpha, size = 1.65) +
    ggplot2::geom_line(data = pred_line, ggplot2::aes(x = decimalLatitude, y = doy),
                       color = "black", linewidth = 1.1, linetype = "solid", inherit.aes = FALSE) +
    ggplot2::geom_line(data = pred_line, ggplot2::aes(x = decimalLatitude, y = hopkins_expected),
                       color = "#D55E00", linewidth = 0.95, linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    ggplot2::scale_shape_manual(values = dyn_shapes) +
    thm +
    ggplot2::labs(
      title = sprintf("A. Verified Latitudinal Spring Cline: DOY ~ Latitude (Slope = +%.2f d/deg, R2 = %.2f)",
                      baseline$slope, baseline$r_squared),
      subtitle = "Solid black = Empirical WLS regression; Dashed vermilion = Hopkins' Law (+4.0 d / deg N)",
      x = "Decimal Latitude (deg N)", y = "Observed Anthesis Day of Year (DOY)", color = "Group", shape = "Group"
    )
  
  med_df <- stats::aggregate(delta_doy ~ plot_group, data = target_df, FUN = function(x) stats::median(x, na.rm = TRUE))
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p2 <- ggplot2::ggplot(target_df, ggplot2::aes(x = delta_doy, y = plot_group, fill = plot_group, color = plot_group)) +
      ggridges::geom_density_ridges(alpha = 0.45, scale = 1.15, rel_min_height = 0.01, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = dyn_pal) +
      ggplot2::scale_color_manual(values = dyn_pal) +
      thm +
      ggplot2::labs(
        title = "B. Phenological Anomalies (Delta DOY) & Allochronic Separation",
        subtitle = paste(sprintf("%s: %+.1f d", med_df$plot_group, med_df$delta_doy), collapse = " | "),
        x = "Phenological Anomaly Delta DOY (Days Relative to Latitudinal Baseline)", y = "Target Group", fill = "Group", color = "Group"
      )
  } else {
    p2 <- ggplot2::ggplot(target_df, ggplot2::aes(x = delta_doy, fill = plot_group, color = plot_group)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 0.85) +
      ggplot2::geom_vline(data = med_df, ggplot2::aes(xintercept = delta_doy, color = plot_group),
                          linetype = "dashed", linewidth = 0.85, show.legend = FALSE) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = dyn_pal) +
      ggplot2::scale_color_manual(values = dyn_pal) +
      thm +
      ggplot2::labs(
        title = "B. Phenological Anomalies (Delta DOY) & Allochronic Separation",
        subtitle = paste(sprintf("%s: %+.1f d", med_df$plot_group, med_df$delta_doy), collapse = " | "),
        x = "Phenological Anomaly Delta DOY (Days Relative to Latitudinal Baseline)", y = "Kernel Density", fill = "Group", color = "Group"
      )
  }
  
  tmp_pdf <- paste0(out_pdf, ".tmp")
  pdf(tmp_pdf, width = 12, height = 5.5)
  tryCatch({
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  }, finally = {
    dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, out_pdf)
}


# ------------------------------------------------------------------------------
# 3. Cross-Modal Consensus & Conflict Auditing Engine
# ------------------------------------------------------------------------------
compute_mean_silhouette_base <- function(cl_vec, D_mat) {
  u_cl <- unique(cl_vec)
  if (length(u_cl) < 2L || nrow(D_mat) < 4L) return(-1.0)
  n_pts <- length(cl_vec)
  sil <- numeric(n_pts)
  for (i in seq_len(n_pts)) {
    same_idx <- setdiff(which(cl_vec == cl_vec[i]), i)
    a_i <- if (length(same_idx) > 0L) mean(D_mat[i, same_idx]) else 0.0
    b_i <- Inf
    for (ck in setdiff(u_cl, cl_vec[i])) {
      oth_idx <- which(cl_vec == ck)
      if (length(oth_idx) > 0L) {
        b_i <- min(b_i, mean(D_mat[i, oth_idx]))
      }
    }
    denom <- max(a_i, b_i, 1e-12)
    sil[i] <- if (is.finite(b_i)) (b_i - a_i) / denom else 0.0
  }
  mean(sil, na.rm = TRUE)
}

resolve_analysis_group <- function(df, intraspecific_mode = FALSE) {
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  if (!is_intra) {
    df$analysis_group <- as.character(df$species_standardized)
    return(df)
  }
  used_gmm <- FALSE
  if ("gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster))) {
    g_raw <- ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster))
    grp_cand <- paste0("Morphotype_K", g_raw)
    tbl <- table(grp_cand)
    if (length(tbl) >= 2 && sum(tbl >= 5) >= 2) {
      df$analysis_group <- as.character(grp_cand)
      used_gmm <- TRUE
      message(sprintf("[Intraspecific Mode] Using Step 04 GMM morphotypes as analysis_group (%s).",
                      paste(names(tbl), tbl, sep = "=", collapse = ", ")))
    }
  }
  if (!used_gmm) {
    n_rows <- nrow(df)
    if (n_rows == 0L) {
      df$analysis_group <- character(0)
      return(df)
    }
    env_pattern <- paste0(
      "^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bio[0-9]*_|bdricm$|resdept_r$|cfvo$|awc$|",
      "alr_|clr_|is_azonal_lithology$|azonal_contact_dist_m$|azonal_proximity_index$|",
      "ph_q10_buffer$|ph_q90_buffer$|min_bedrock_depth_buffer_cm$|max_cfvo_buffer_pct$|",
      "azonal_lithology_presence_buffer$|elevation$|elev_|clim_|edaph_|env_)"
    )
    exclude_cols <- c(
      "catalogNumber", "plant_individual_id", "latitude", "longitude",
      "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
      "coordinateUncertaintyInMeters", "buffer_radius_m",
      "high_coord_uncertainty_flag", "is_training_eligible_coord",
      "edaphic_best_fit_taxon", "edaphic_z_score", "soil_silt_pct"
    )
    cand_cols <- setdiff(grep(env_pattern, names(df), value = TRUE, ignore.case = TRUE), exclude_cols)
    num_env <- list()
    min_fin <- min(3L, max(2L, n_rows))
    for (cl in cand_cols) {
      if (is.numeric(df[[cl]]) || is.logical(df[[cl]])) {
        v <- suppressWarnings(as.numeric(df[[cl]]))
        if (sum(is.finite(v)) >= min_fin) {
          med_v <- stats::median(v[is.finite(v)], na.rm = TRUE)
          v[!is.finite(v)] <- if (is.finite(med_v)) med_v else 0.0
          if (stats::sd(v) > 1e-8) {
            num_env[[cl]] <- v
          }
        }
      }
    }
    if ("fine_soil_ph" %in% names(num_env)) num_env[["soil_ph"]] <- NULL
    if ("fine_soil_sand" %in% names(num_env)) num_env[["soil_sand"]] <- NULL
    if ("fine_soil_cec" %in% names(num_env)) num_env[["soil_cec"]] <- NULL

    if (length(num_env) == 0L || n_rows < 2L) {
      idx_half <- seq_len(n_rows) <= max(1L, floor(n_rows / 2))
      df$analysis_group <- ifelse(idx_half, "Ecotype_Cluster_1", "Ecotype_Cluster_2")
      tbl <- table(df$analysis_group)
      message(sprintf("[Intraspecific Mode] Zero environmental variance; fallback ecotypes (%s).",
                      paste(names(tbl), tbl, sep = "=", collapse = ", ")))
      return(df)
    }

    E_mat <- do.call(cbind, num_env)
    E_scaled <- scale(E_mat, center = TRUE, scale = TRUE)
    E_scaled[!is.finite(E_scaled)] <- 0.0

    if (ncol(E_scaled) >= 2L) {
      pca_fit <- stats::prcomp(E_scaled, center = FALSE, scale. = FALSE)
      n_pcs <- min(2L, ncol(pca_fit$x))
      pca_scores <- as.matrix(pca_fit$x[, seq_len(n_pcs), drop = FALSE])
      if (ncol(pca_scores) == 1L) {
        pca_scores <- cbind(PC1 = pca_scores[, 1], PC2 = rep(0.0, n_rows))
      }
      pc1_load <- abs(pca_fit$rotation[, 1])
      top_var <- names(which.max(pc1_load))[1]
    } else {
      pca_scores <- cbind(PC1 = as.numeric(E_scaled[, 1]), PC2 = rep(0.0, n_rows))
      top_var <- colnames(E_scaled)[1]
    }

    n_distinct <- nrow(unique(round(pca_scores, 8)))
    if (n_distinct < 2L) {
      idx_half <- seq_len(n_rows) <= max(1L, floor(n_rows / 2))
      df$analysis_group <- ifelse(idx_half, "Ecotype_Cluster_1", "Ecotype_Cluster_2")
      return(df)
    }

    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    set.seed(42)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)

    km2 <- stats::kmeans(pca_scores, centers = 2, nstart = 25)
    best_cl <- km2$cluster
    best_k <- 2L
    if (n_rows >= 30L && n_distinct >= 3L) {
      km3 <- tryCatch(stats::kmeans(pca_scores, centers = 3, nstart = 25), error = function(e) NULL)
      if (!is.null(km3) && min(table(km3$cluster)) >= 5L) {
        D_pca <- as.matrix(stats::dist(pca_scores))
        sil2 <- compute_mean_silhouette_base(km2$cluster, D_pca)
        sil3 <- compute_mean_silhouette_base(km3$cluster, D_pca)
        if (is.finite(sil3) && is.finite(sil2) && sil3 > (sil2 + 0.05)) {
          best_cl <- km3$cluster
          best_k <- 3L
        }
      }
    }

    if (min(table(best_cl)) < 2L && n_rows >= 4L) {
      best_k <- 2L
      med_pc1 <- stats::median(pca_scores[, 1], na.rm = TRUE)
      best_cl <- ifelse(pca_scores[, 1] <= med_pc1, 1L, 2L)
      if (length(unique(best_cl)) < 2L || min(table(best_cl)) < 2L) {
        best_cl <- ifelse(rank(pca_scores[, 1], ties.method = "first") <= floor(n_rows / 2), 1L, 2L)
      }
    }

    clean_top_var <- gsub("[^A-Za-z0-9_]", "_", top_var)
    cl_means <- tapply(E_mat[, top_var], best_cl, mean, na.rm = TRUE)
    ord_cl <- as.integer(names(sort(cl_means, decreasing = FALSE)))
    tier_suffixes <- if (best_k == 3L) c("Low", "Mid", "High") else c("Low", "High")
    lbl_map <- character(best_k)
    names(lbl_map) <- as.character(ord_cl)
    for (k_i in seq_len(best_k)) {
      lbl_map[[as.character(ord_cl[k_i])]] <- sprintf("Ecotype_%d_%s_%s", k_i, clean_top_var, tier_suffixes[k_i])
    }
    df$analysis_group <- unname(lbl_map[as.character(best_cl)])
    tbl <- table(df$analysis_group)
    message(sprintf("[Intraspecific Mode] Partitioned single-taxon cohort into PCA-env ecotypes (K=%d, Axis-1 driver=%s: %s).",
                    best_k, clean_top_var, paste(names(tbl), tbl, sep = "=", collapse = ", ")))
  }
  return(df)
}

execute_crossmodal_consensus <- function(df, intraspecific_mode = FALSE) {
  message("Executing 4-way cross-modal consensus checks (Morpho + Vision + Fine-Grain Edaphic + Phenology)...")
  if (exists("flag_coordinate_uncertainty", mode = "function")) {
    df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0)
  } else {
    if (!("coordinateUncertaintyInMeters" %in% names(df))) df$coordinateUncertaintyInMeters <- NA_real_
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    df$high_coord_uncertainty_flag <- !is.na(u_m) & u_m > 10000.0
    df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
    df$buffer_radius_m <- ifelse(is.na(u_m) | !is.finite(u_m) | u_m <= 0, 250.0, pmax(100.0, pmin(1000.0, u_m)))
  }

  df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  active_taxa <- if (is_intra) sort(unique(stats::na.omit(df$analysis_group))) else get_active_packera_taxa(df$species_standardized)
  if (length(active_taxa) == 0L) active_taxa <- sort(unique(as.character(stats::na.omit(df$analysis_group))))
  grp_vec <- if (is_intra) df$analysis_group else df$species_standardized

  n <- nrow(df); conflicts <- data.frame(catalogNumber = df$catalogNumber, stringsAsFactors = FALSE)
  
  # Ensure fine-grain columns exist (fallback to SoilGrids if evaluate_multimodal_concordance is called standalone)
  if (!("fine_soil_ph" %in% names(df))) df$fine_soil_ph <- if ("soil_ph" %in% names(df)) df$soil_ph else rep(5.5, n)
  if (!("fine_soil_sand" %in% names(df))) df$fine_soil_sand <- if ("soil_sand" %in% names(df)) df$soil_sand else rep(50.0, n)
  if (!("fine_soil_clay" %in% names(df))) df$fine_soil_clay <- pmax(2.0, pmin(65.0, 100.0 - df$fine_soil_sand - 22.0))
  if (!("fine_bedrock_depth_cm" %in% names(df))) {
    df$fine_bedrock_depth_cm <- if ("ssurgo_bedrock_depth_cm" %in% names(df)) df$ssurgo_bedrock_depth_cm else rep(150.0, n)
  }
  if (!("ssurgo_rock_outcrop_pct" %in% names(df))) {
    df$ssurgo_rock_outcrop_pct <- if ("is_rock_outcrop" %in% names(df)) ifelse(df$is_rock_outcrop, 45.0, 0.0) else rep(0.0, n)
  }
  if (!("bdricm" %in% names(df))) df$bdricm <- df$fine_bedrock_depth_cm
  if (!("resdept_r" %in% names(df))) df$resdept_r <- df$fine_bedrock_depth_cm
  if (!("cfvo" %in% names(df))) df$cfvo <- ifelse(df$ssurgo_rock_outcrop_pct > 15.0, 28.0, 8.0)
  if (!("awc" %in% names(df))) df$awc <- if ("ssurgo_aws_0_50cm" %in% names(df)) round(df$ssurgo_aws_0_50cm / 50.0, 4) else rep(0.14, n)
  if (!("fine_soil_cec" %in% names(df))) df$fine_soil_cec <- if ("soil_cec" %in% names(df)) df$soil_cec else rep(12.0, n)
  if (!("alr_sand_clay" %in% names(df))) {
    df$alr_sand_clay <- round(log((pmax(df$fine_soil_sand, 0.5) + 0.5) / (pmax(df$fine_soil_clay, 0.5) + 0.5)), 4)
  }
  if (!("ph_q10_buffer" %in% names(df))) df$ph_q10_buffer <- pmax(3.5, df$fine_soil_ph - 0.20)
  if (!("ph_q90_buffer" %in% names(df))) df$ph_q90_buffer <- pmin(9.0, df$fine_soil_ph + 0.20)
  if (!("min_bedrock_depth_buffer_cm" %in% names(df))) df$min_bedrock_depth_buffer_cm <- pmax(2.0, df$fine_bedrock_depth_cm * 0.75)
  if (!("max_cfvo_buffer_pct" %in% names(df))) df$max_cfvo_buffer_pct <- pmin(95.0, df$cfvo * 1.25)
  if (!("azonal_lithology_presence_buffer" %in% names(df))) {
    df$azonal_lithology_presence_buffer <- as.integer(df$ssurgo_rock_outcrop_pct >= 15.0 | df$fine_bedrock_depth_cm <= 65.0)
  }

  # Phenological circular stats & Empirical Bayes / Beta-Binomial Authority Weights
  doy <- ifelse(is.na(df$doy) | df$doy < 1 | df$doy > 366, 120, df$doy)
  theta <- 2 * pi * doy / 365.25
  df$pheno_theta <- theta
  df$authority_weight <- resolve_step06_authority_weights(df)
  df$spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df) else rep(1.0, nrow(df))
  df$case_wts <- df$authority_weight * df$spatial_balance_weight

  w_mean <- function(x, w) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    ok <- !is.na(x) & !is.na(w) & is.finite(x) & w > 0
    if (!any(ok)) return(mean(as.numeric(x), na.rm = TRUE))
    sum(x[ok] * w[ok]) / sum(w[ok])
  }
  w_sd <- function(x, w, min_sd = 0.25) {
    if (is.null(x) || length(x) == 0) return(min_sd)
    ok <- !is.na(x) & !is.na(w) & is.finite(x) & w > 0
    if (sum(ok) < 2) return(max(stats::sd(as.numeric(x), na.rm = TRUE), min_sd, na.rm = TRUE))
    wm <- sum(x[ok] * w[ok]) / sum(w[ok])
    denom <- sum(w[ok]) - (sum(w[ok]^2) / sum(w[ok]))
    wv <- if (denom > 1e-6) sum(w[ok] * (x[ok] - wm)^2) / denom else stats::var(as.numeric(x), na.rm = TRUE)
    max(sqrt(max(wv, 0, na.rm = TRUE)), min_sd, na.rm = TRUE)
  }

  # Compute species/ecotype-specific circular means and fine-grain eco-edaphic envelopes
  # Strictly on Tier_1_Gold (w=3.0) & Tier_2_Silver (w=1.5) anchors with coordinateUncertaintyInMeters <= 10000
  elig_coord <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, nrow(df))
  elig_coord[is.na(elig_coord)] <- TRUE
  is_anchor <- df$determiner_tier %in% c("Tier_1_Gold", "Tier_2_Silver") & elig_coord & grp_vec %in% active_taxa
  anchor_df <- df[is_anchor, ]
  if (nrow(anchor_df) < 10) anchor_df <- df[elig_coord & grp_vec %in% active_taxa, ]
  if (nrow(anchor_df) < 10) anchor_df <- df[grp_vec %in% active_taxa, ]

  cand_eco <- intersect(c(
    "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm",
    "bdricm", "cfvo", "awc", "fine_soil_cec", "alr_sand_clay", "azonal_proximity_index",
    "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct",
    "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm", "topo_tpi_30m", "topo_twi_30m",
    "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "delta_doy"
  ), names(df))
  eco_cols <- cand_eco[sapply(cand_eco, function(cl) sum(is.finite(as.numeric(df[[cl]]))) >= 5 && stats::sd(as.numeric(df[[cl]]), na.rm = TRUE) > 1e-6)]
  if (length(eco_cols) == 0) eco_cols <- intersect(c("fine_soil_ph", "fine_soil_sand", "fine_soil_cec"), names(df))
  min_sd_map <- c(
    fine_soil_ph = 0.25, fine_soil_sand = 4.0, fine_soil_clay = 3.0, fine_bedrock_depth_cm = 15.0,
    bdricm = 15.0, cfvo = 4.0, awc = 0.02, fine_soil_cec = 2.0, alr_sand_clay = 0.25, azonal_proximity_index = 0.10,
    ph_q10_buffer = 0.25, ph_q90_buffer = 0.25, min_bedrock_depth_buffer_cm = 15.0, max_cfvo_buffer_pct = 5.0,
    ssurgo_rock_outcrop_pct = 8.0, ssurgo_aws_0_50cm = 1.0, topo_tpi_30m = 1.5, topo_twi_30m = 0.8,
    bio1_temp_mean = 1.0, bio4_temp_seasonality = 25.0, bio12_precip_annual = 50.0, delta_doy = 5.0
  )

  taxa_profiles <- lapply(active_taxa, function(tx) {
    sub_df <- anchor_df[(if (is_intra) anchor_df$analysis_group else anchor_df$species_standardized) == tx, ]
    if (nrow(sub_df) < 3) sub_df <- df[grp_vec == tx, ]
    case_wts <- sub_df$authority_weight * sub_df$spatial_balance_weight
    w_sub <- case_wts
    sin_m <- w_mean(sin(sub_df$pheno_theta), w_sub)
    cos_m <- w_mean(cos(sub_df$pheno_theta), w_sub)
    eco_m <- sapply(eco_cols, function(cl) w_mean(as.numeric(sub_df[[cl]]), w_sub))
    eco_s <- sapply(eco_cols, function(cl) w_sd(as.numeric(sub_df[[cl]]), w_sub, ifelse(cl %in% names(min_sd_map), min_sd_map[[cl]], 0.25)))
    cec_vec <- if ("fine_soil_cec" %in% names(sub_df)) sub_df$fine_soil_cec else sub_df$soil_cec
    list(
      mean_theta = atan2(sin_m, cos_m), R_bar = sqrt(sin_m^2 + cos_m^2),
      ph_mean = w_mean(sub_df$fine_soil_ph, w_sub), ph_sd = w_sd(sub_df$fine_soil_ph, w_sub, 0.25),
      sand_mean = w_mean(sub_df$fine_soil_sand, w_sub), sand_sd = w_sd(sub_df$fine_soil_sand, w_sub, 4.0),
      bedrock_mean = w_mean(sub_df$fine_bedrock_depth_cm, w_sub), bedrock_sd = w_sd(sub_df$fine_bedrock_depth_cm, w_sub, 15.0),
      cec_mean = w_mean(cec_vec, w_sub), cec_sd = w_sd(cec_vec, w_sub, 2.0),
      eco_mean = eco_m, eco_sd = eco_s
    )
  })
  names(taxa_profiles) <- active_taxa
  for (cl in eco_cols) {
    bar_sd <- mean(sapply(taxa_profiles, function(p) p$eco_sd[[cl]]), na.rm = TRUE)
    for (tx in active_taxa) taxa_profiles[[tx]]$eco_sd[[cl]] <- 0.5 * taxa_profiles[[tx]]$eco_sd[[cl]] + 0.5 * bar_sd
  }

  compute_specimen_edaphic_z <- function(idx, target_tx) {
    if (!(target_tx %in% names(taxa_profiles))) return(0.0)
    pr <- taxa_profiles[[target_tx]]
    z_ph <- abs((as.numeric(df$fine_soil_ph[idx]) - pr$ph_mean) / max(pr$ph_sd, 0.25, na.rm = TRUE))
    z_sand <- abs((as.numeric(df$fine_soil_sand[idx]) - pr$sand_mean) / max(pr$sand_sd, 4.0, na.rm = TRUE))
    z_bed <- abs((as.numeric(df$fine_bedrock_depth_cm[idx]) - pr$bedrock_mean) / max(pr$bedrock_sd, 15.0, na.rm = TRUE))
    zv <- c(z_ph, z_sand, z_bed)
    zv <- zv[is.finite(zv)]
    if (length(zv) == 0) return(0.0)
    sqrt(mean(zv^2))
  }

  # Precompute empirical outcrop / azonal frequency and odds ratio per active group in anchor_df
  compute_row_outcrop_azonal <- function(d_sub) {
    n_s <- nrow(d_sub)
    if (n_s == 0L) return(logical(0))
    rp_s <- if ("ssurgo_rock_outcrop_pct" %in% names(d_sub)) suppressWarnings(as.numeric(d_sub$ssurgo_rock_outcrop_pct)) else rep(0.0, n_s)
    rp_s[is.na(rp_s)] <- 0.0
    bd_s <- if ("fine_bedrock_depth_cm" %in% names(d_sub)) suppressWarnings(as.numeric(d_sub$fine_bedrock_depth_cm)) else rep(150.0, n_s)
    bd_s[is.na(bd_s)] <- 150.0
    mb_s <- if ("min_bedrock_depth_buffer_cm" %in% names(d_sub)) suppressWarnings(as.numeric(d_sub$min_bedrock_depth_buffer_cm)) else bd_s
    mb_s[is.na(mb_s)] <- bd_s[is.na(mb_s)]
    ab_s <- if ("azonal_lithology_presence_buffer" %in% names(d_sub)) suppressWarnings(as.integer(d_sub$azonal_lithology_presence_buffer)) > 0L else rep(FALSE, n_s)
    ab_s[is.na(ab_s)] <- FALSE
    ro_s <- if ("is_rock_outcrop" %in% names(d_sub)) as.logical(d_sub$is_rock_outcrop) else rep(FALSE, n_s)
    ro_s[is.na(ro_s)] <- FALSE
    az_s <- if ("is_azonal_lithology" %in% names(d_sub)) as.logical(d_sub$is_azonal_lithology) else rep(FALSE, n_s)
    az_s[is.na(az_s)] <- FALSE
    (rp_s >= 15.0) | (bd_s <= 65.0) | (mb_s <= 50.0) | ab_s | ro_s | az_s
  }

  anc_grp_vec <- if (is_intra) as.character(anchor_df$analysis_group) else as.character(anchor_df$species_standardized)
  anc_out_az <- compute_row_outcrop_azonal(anchor_df)
  is_specialist_group <- stats::setNames(rep(FALSE, length(active_taxa)), active_taxa)
  for (tx in active_taxa) {
    tx_m <- !is.na(anc_grp_vec) & anc_grp_vec == tx
    if (sum(tx_m) < 3L) {
      # Fallback to full cohort for tiny groups
      full_m <- !is.na(grp_vec) & grp_vec == tx
      full_oa <- compute_row_outcrop_azonal(df)
      a_cnt <- sum(full_m & full_oa, na.rm = TRUE)
      b_cnt <- sum(full_m & !full_oa, na.rm = TRUE)
      c_cnt <- sum((!full_m) & full_oa, na.rm = TRUE)
      d_cnt <- sum((!full_m) & !full_oa, na.rm = TRUE)
    } else {
      a_cnt <- sum(tx_m & anc_out_az, na.rm = TRUE)
      b_cnt <- sum(tx_m & !anc_out_az, na.rm = TRUE)
      c_cnt <- sum((!tx_m) & anc_out_az, na.rm = TRUE)
      d_cnt <- sum((!tx_m) & !anc_out_az, na.rm = TRUE)
    }
    emp_freq <- a_cnt / max(1L, a_cnt + b_cnt)
    emp_or <- if ((c_cnt + d_cnt) > 0L) {
      ((a_cnt + 0.5) * (d_cnt + 0.5)) / ((b_cnt + 0.5) * (c_cnt + 0.5))
    } else {
      ((a_cnt + 0.5) * 92.5) / ((b_cnt + 0.5) * 8.5)
    }
    is_specialist_group[[tx]] <- isTRUE(emp_freq > 0.30 || emp_or > 2.0)
  }

  # Evaluate modal concordance per specimen
  conf_flags <- logical(n); triage_cats <- character(n); concordances <- numeric(n); discord_notes <- character(n)
  edaph_best_vec <- character(n); edaph_z_vec <- numeric(n); outcrop_spec_flags <- logical(n)

  for (i in seq_len(n)) {
    given_sp <- grp_vec[i]
    morph_sp <- if ("cda_predicted_taxon" %in% names(df) && !is.na(df$cda_predicted_taxon[i])) df$cda_predicted_taxon[i] else given_sp
    vis_sp   <- if ("vision_predicted_label" %in% names(df) && !is.na(df$vision_predicted_label[i])) df$vision_predicted_label[i] else given_sp
    if (is_intra) {
      if (!(morph_sp %in% active_taxa)) morph_sp <- given_sp
      if (!(vis_sp %in% active_taxa)) vis_sp <- given_sp
    }
    c_err    <- if ("c_error" %in% names(df) && !is.na(df$c_error[i])) df$c_error[i] else 0.0
    cda_prob <- if ("cda_posterior_prob" %in% names(df) && !is.na(df$cda_posterior_prob[i])) as.numeric(df$cda_posterior_prob[i]) else NA_real_

    # Regularized Eco-Edaphic Gaussian log-likelihood score
    edaph_scores <- sapply(active_taxa, function(tx) {
      prof <- taxa_profiles[[tx]]
      ll <- 0.0
      for (cl in eco_cols) {
        xv <- as.numeric(df[[cl]][i]); mv <- prof$eco_mean[[cl]]; sv <- prof$eco_sd[[cl]]
        if (is.finite(xv) && is.finite(mv) && is.finite(sv) && sv > 0) ll <- ll - 0.5 * (((xv - mv) / sv)^2 + 2.0 * log(sv))
      }
      ll
    })
    edaph_sp <- if (all(is.na(edaph_scores))) given_sp else active_taxa[which.max(edaph_scores)][1]
    edaph_best_vec[i] <- edaph_sp
    edaph_z_given <- compute_specimen_edaphic_z(i, given_sp)
    edaph_z_vec[i] <- round(edaph_z_given, 3)

    alt_parent <- if (morph_sp != given_sp) morph_sp else if (vis_sp != given_sp) vis_sp else edaph_sp
    edaph_z_alt <- compute_specimen_edaphic_z(i, alt_parent)

    rock_pct_i <- if (!is.na(df$ssurgo_rock_outcrop_pct[i])) as.numeric(df$ssurgo_rock_outcrop_pct[i]) else 0.0
    bed_cm_i <- if (!is.na(df$fine_bedrock_depth_cm[i])) as.numeric(df$fine_bedrock_depth_cm[i]) else 150.0
    min_bed_buf_i <- if ("min_bedrock_depth_buffer_cm" %in% names(df) && !is.na(df$min_bedrock_depth_buffer_cm[i])) as.numeric(df$min_bedrock_depth_buffer_cm[i]) else bed_cm_i
    azonal_buf_i <- if ("azonal_lithology_presence_buffer" %in% names(df) && !is.na(df$azonal_lithology_presence_buffer[i])) as.integer(df$azonal_lithology_presence_buffer[i]) > 0L else FALSE
    is_azonal_i <- ("is_azonal_lithology" %in% names(df) && isTRUE(as.logical(df$is_azonal_lithology[i])))
    is_outcrop_i <- isTRUE(rock_pct_i >= 15.0 || bed_cm_i <= 65.0 || min_bed_buf_i <= 50.0 || azonal_buf_i ||
                           ("is_rock_outcrop" %in% names(df) && isTRUE(df$is_rock_outcrop[i])))
    tx_is_spec <- !is.na(given_sp) && (as.character(given_sp) %in% names(is_specialist_group)) && isTRUE(is_specialist_group[[as.character(given_sp)]])
    outcrop_spec_flags[i] <- isTRUE((is_outcrop_i || is_azonal_i) && tx_is_spec)

    # Phenological circular divergence
    prof_g <- taxa_profiles[[given_sp]]
    ang_diff <- if (!is.null(prof_g)) abs(atan2(sin(df$pheno_theta[i] - prof_g$mean_theta), cos(df$pheno_theta[i] - prof_g$mean_theta))) else 0.0
    if (is.na(ang_diff) || !is.finite(ang_diff)) ang_diff <- 0.0
    pheno_agree <- isTRUE(ang_diff < (pi / 3)) # within ~60 days

    # Agreement matrix
    matches <- c(morph_sp == given_sp, vis_sp == given_sp, edaph_sp == given_sp, pheno_agree)
    concordance <- sum(matches, na.rm = TRUE) / 4.0
    concordances[i] <- round(concordance, 2)

    # Multi-evidence decision triage
    if (morph_sp != given_sp && vis_sp != given_sp && (morph_sp == vis_sp || edaph_sp == morph_sp || c_err >= 0.85)) {
      conf_flags[i] <- TRUE; triage_cats[i] <- "Severe_Triple_Stream_Conflict"
      discord_notes[i] <- sprintf("Given %s contradicted by Morphology (%s), Vision (%s), and Fine-Grain Edaphic (%s)",
                                  given_sp, morph_sp, vis_sp, edaph_sp)
    } else if (edaph_z_given > 3.0 && morph_sp == given_sp) {
      conf_flags[i] <- TRUE; triage_cats[i] <- "Edaphic_Envelope_Violation"
      discord_notes[i] <- sprintf("Specimen outside physiological edaphic tolerance (fine_pH=%.2f, fine_Sand=%.1f%%, bedrock=%.0fcm, z=%.2f)",
                                  df$fine_soil_ph[i], df$fine_soil_sand[i], bed_cm_i, edaph_z_given)
    } else if (!pheno_agree && concordance >= 0.75) {
      conf_flags[i] <- TRUE; triage_cats[i] <- "Phenological_Anomaly_Outlier"
      discord_notes[i] <- sprintf("Flowering DOY %d diverges from seasonal peak (ang_diff=%.2f rad)", doy[i], ang_diff)
    } else if (concordance == 0.50) {
      if (morph_sp != given_sp && vis_sp != given_sp && morph_sp == vis_sp) {
        conf_flags[i] <- TRUE; triage_cats[i] <- "Severe_Triple_Stream_Conflict"
        discord_notes[i] <- sprintf("Both Morphology and Vision agree on %s over given %s", morph_sp, given_sp)
      } else if ((morph_sp != vis_sp || (!is.na(cda_prob) && cda_prob >= 0.40 && cda_prob <= 0.65)) &&
                 (c_err < 0.80) && (edaph_z_given <= 2.5 && edaph_z_alt <= 2.5)) {
        conf_flags[i] <- TRUE; triage_cats[i] <- "Putative_Hybrid_Zone_Intergrade"
        discord_notes[i] <- sprintf("Intermediate morphometrics & fine-grain edaphic overlap between %s and %s", given_sp, alt_parent)
      } else if ((morph_sp != given_sp && edaph_sp == morph_sp && (is.na(cda_prob) || cda_prob >= 0.55)) ||
                 (vis_sp != given_sp && c_err >= 0.80)) {
        conf_flags[i] <- TRUE; triage_cats[i] <- "Severe_Triple_Stream_Conflict"
        discord_notes[i] <- sprintf("Cross-modal misidentification signal favoring %s over given %s", alt_parent, given_sp)
      } else {
        conf_flags[i] <- FALSE; triage_cats[i] <- "Clean_MultiModal_Consensus"
        discord_notes[i] <- "Concordant across primary morphological, vision, and ecological axes"
      }
    } else {
      conf_flags[i] <- FALSE; triage_cats[i] <- "Clean_MultiModal_Consensus"
      discord_notes[i] <- "Concordant across morphological, vision, and ecological axes"
    }
  }

  df$multimodal_concordance <- concordances
  df$multimodal_conflict_flag <- conf_flags
  df$triage_category <- triage_cats
  df$conflict_discordance_notes <- discord_notes
  df$edaphic_best_fit_taxon <- edaph_best_vec
  df$edaphic_z_score <- edaph_z_vec
  df$substrate_specialist_flag <- outcrop_spec_flags
  df$piedmont_outcrop_specialist_flag <- outcrop_spec_flags
  return(df)
}
evaluate_multimodal_concordance <- execute_crossmodal_consensus

# ------------------------------------------------------------------------------
# 4. Moran's Eigenvector Maps & Authority-Anchored Spatial Random Forests (Roberts et al. 2017; Valavi et al. 2021)
# ------------------------------------------------------------------------------
compute_inverse_class_frequency_weights <- function(y_vec) {
  y_chr <- as.character(y_vec)
  tbl <- table(y_chr)
  n_tot <- length(y_chr)
  n_cls <- max(length(tbl), 1L)
  w_map <- n_tot / (n_cls * pmax(as.numeric(tbl), 1.0))
  names(w_map) <- names(tbl)
  w_out <- as.numeric(w_map[y_chr])
  w_out[is.na(w_out) | !is.finite(w_out)] <- 1.0
  w_out / max(mean(w_out, na.rm = TRUE), 1e-8)
}

estimate_spatial_autocorrelation_block_size <- function(coords, X_env, default_km = 55.0) {
  n <- nrow(coords)
  if (is.null(coords) || n < 10) return(list(range_km = default_km, range_m = default_km * 1000.0, range_deg = default_km / 111.0))
  range_deg <- tryCatch({
    d_mat <- as.matrix(stats::dist(coords))
    d_vec <- d_mat[upper.tri(d_mat)]
    if (length(d_vec) < 5) return(list(range_km = default_km, range_m = default_km * 1000.0, range_deg = default_km / 111.0))
    # Empirical variogram range estimation across standardized environmental PC1
    pc1 <- tryCatch(stats::prcomp(X_env, center = TRUE, scale. = TRUE)$x[, 1], error = function(e) X_env[, 1])
    v_mat <- outer(pc1, pc1, function(a, b) 0.5 * (a - b)^2)
    v_vec <- v_mat[upper.tri(v_mat)]
    sill <- stats::var(pc1, na.rm = TRUE)
    if (!is.finite(sill) || sill <= 1e-8) return(list(range_km = default_km, range_m = default_km * 1000.0, range_deg = default_km / 111.0))
    breaks <- stats::quantile(d_vec, probs = seq(0, 0.8, length.out = 9), na.rm = TRUE)
    breaks <- unique(breaks)
    if (length(breaks) < 3) return(list(range_km = default_km, range_m = default_km * 1000.0, range_deg = default_km / 111.0))
    bins <- cut(d_vec, breaks = breaks, include.lowest = TRUE)
    gamma_bins <- tapply(v_vec, bins, mean, na.rm = TRUE)
    dist_bins  <- tapply(d_vec, bins, mean, na.rm = TRUE)
    reach_idx <- which(gamma_bins >= 0.85 * sill)
    est_deg <- if (length(reach_idx) > 0) as.numeric(dist_bins[reach_idx[1]]) else as.numeric(stats::median(d_vec, na.rm = TRUE) * 0.25)
    pmax(0.20, pmin(2.50, est_deg))
  }, error = function(e) default_km / 111.0)
  if (!is.finite(range_deg) || is.na(range_deg)) range_deg <- default_km / 111.0
  range_km <- round(range_deg * 111.0, 2)
  list(range_km = range_km, range_m = round(range_km * 1000.0, 0), range_deg = round(range_deg, 4))
}

assign_spatial_grid_blocks <- function(coords, block_deg = 0.50, k = 5L) {
  lon <- as.numeric(coords[, 1]); lat <- as.numeric(coords[, 2])
  b_deg <- if (is.finite(block_deg) && block_deg > 0.05) block_deg else 0.50
  gx <- floor((lon - min(lon, na.rm = TRUE)) / b_deg)
  gy <- floor((lat - min(lat, na.rm = TRUE)) / b_deg)
  cell_ids <- paste0(gx, "_", gy)
  u_cells <- unique(cell_ids)
  if (length(u_cells) < 2) {
    return(rep(seq_len(k), length.out = length(lon)))
  }
  cell_folds <- ((seq_along(u_cells) - 1L) %% k) + 1L
  names(cell_folds) <- u_cells
  as.integer(cell_folds[cell_ids])
}

compute_spatial_rf_mems <- function(df, env_vars, n_trees = 500, intraspecific_mode = FALSE) {
  message("Generating Moran's Eigenvector Maps (MEM1..MEM5) & fitting Authority-Anchored Spatial Random Forests...")
  if (!("analysis_group" %in% names(df)) || isTRUE(intraspecific_mode)) {
    df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
  }
  if (exists("flag_coordinate_uncertainty", mode = "function") && !("is_training_eligible_coord" %in% names(df))) {
    df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0)
  }
  coords <- as.matrix(df[, c("longitude", "latitude")])
  n <- nrow(coords)

  MEMs <- NULL
  if (requireNamespace("adespatial", quietly = TRUE)) {
    MEMs <- tryCatch({
      db_res <- adespatial::dbmem(coords, MEM.autocor = "positive")
      m_mat <- as.matrix(db_res)
      if (ncol(m_mat) >= 1) {
        num_mems <- min(ncol(m_mat), 5)
        m_sub <- m_mat[, seq_len(num_mems), drop = FALSE]
        colnames(m_sub) <- paste0("MEM", seq_len(num_mems))
        m_sub
      } else {
        NULL
      }
    }, error = function(e) NULL)
  }
  if (is.null(MEMs)) {
    # Fallback to distance-matrix cmdscale / double-centered eigenfunctions (MEM1..MEM5)
    D <- as.matrix(dist(coords))
    thresh <- max(apply(D, 1, function(r) sort(r)[min(8, n)]))
    W <- exp(-D / (thresh + 1e-5)) * (D <= thresh); diag(W) <- 0
    H <- diag(n) - matrix(1 / n, n, n)
    M <- H %*% W %*% H
    eig <- eigen(M, symmetric = TRUE)
    pos_idx <- which(eig$values > 1e-4)
    num_mems <- min(length(pos_idx), 5)
    MEMs <- eig$vectors[, pos_idx[1:num_mems], drop = FALSE]
    colnames(MEMs) <- paste0("MEM", seq_len(num_mems))
  }

  # Prepare feature matrices & Block-Separate Collinearity Pruning (Dormann et al. 2013: |r| <= 0.70, VIF <= 5.0 separately within climatic and edaphic blocks)
  env_vars <- intersect(env_vars, names(df))
  X_env <- as.matrix(df[, env_vars, drop = FALSE])
  for (j in seq_len(ncol(X_env))) {
    col_v <- X_env[, j]
    if (any(!is.finite(col_v))) {
      med_j <- stats::median(col_v[is.finite(col_v)], na.rm = TRUE)
      if (!is.finite(med_j)) med_j <- 0.0
      X_env[!is.finite(col_v), j] <- med_j
    }
  }
  col_screen <- if (exists("prune_collinear_predictors", mode = "function")) {
    prune_collinear_predictors(X_env, cor.threshold = 0.70, vif.threshold = 5.0)
  } else if (exists("screen_collinear_predictors_by_block", mode = "function")) {
    screen_collinear_predictors_by_block(X_env, cor.threshold = 0.70, vif.threshold = 5.0)
  } else if (exists("screen_collinear_predictors", mode = "function")) {
    screen_collinear_predictors(X_env, cor.threshold = 0.70, vif.threshold = 5.0)
  } else {
    ac <- if (requireNamespace("spatialRF", quietly = TRUE)) tryCatch(spatialRF::auto_cor(x = as.data.frame(X_env), cor.threshold = 0.70), error = function(e) NULL) else NULL
    X_s <- if (!is.null(ac) && !is.null(ac$selected.variables)) X_env[, ac$selected.variables, drop = FALSE] else X_env
    av <- if (requireNamespace("spatialRF", quietly = TRUE)) tryCatch(spatialRF::auto_vif(x = as.data.frame(X_s), vif.threshold = 5.0), error = function(e) NULL) else NULL
    s_v <- if (!is.null(av) && !is.null(av$selected.variables)) av$selected.variables else colnames(X_s)
    list(selected_variables = s_v, pruned_variables = setdiff(colnames(X_env), s_v))
  }
  if (length(col_screen$selected_variables) >= 2) {
    fe_c <- grep("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)", colnames(X_env), value = TRUE)
    mc_c <- grep("^(bio|delta_doy|doy)", colnames(X_env), value = TRUE, ignore.case = TRUE)
    sel_env <- col_screen$selected_variables
    if (length(fe_c) > 0 && length(intersect(sel_env, fe_c)) == 0) sel_env <- c(sel_env, fe_c[1])
    if (length(mc_c) > 0 && length(intersect(sel_env, mc_c)) == 0) sel_env <- c(sel_env, mc_c[1])
    env_vars <- intersect(colnames(X_env), unique(sel_env))
    X_env <- X_env[, env_vars, drop = FALSE]
  }
  X_all <- cbind(X_env, MEMs)
  y <- factor(df$analysis_group)

  # Authority-Anchored Training Cohort Selection:
  # Strictly Tier_1_Gold + Tier_2_Silver anchors with coordinateUncertaintyInMeters <= 10000
  elig_coord <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, n)
  elig_coord[is.na(elig_coord)] <- TRUE
  if ("coordinateUncertaintyInMeters" %in% names(df)) {
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    elig_coord <- elig_coord & (is.na(u_m) | u_m <= 10000.0)
  }
  det_tier <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n)
  is_gold_silver <- det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")
  anchor_mask <- is_gold_silver & elig_coord & !is.na(y)
  if (sum(anchor_mask) < 6 || length(unique(y[anchor_mask])) < 2) {
    # Fallback for tiny synthetic unit-test fixtures where all rows may lack determiner_tier
    anchor_mask <- elig_coord & !is.na(y)
  }
  if (sum(anchor_mask) < 4 || length(unique(y[anchor_mask])) < 2) {
    anchor_mask <- !is.na(y)
  }
  anchor_idx <- which(anchor_mask)
  passive_idx <- setdiff(seq_len(n), anchor_idx)
  y_anchor <- droplevels(y[anchor_idx])

  # Inverse class-frequency weights (Valavi et al. 2021) * Empirical Bayes / Beta-Binomial Authority weights * Spatial declustering weights
  spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df) else rep(1.0, n)
  auth_wts <- resolve_step06_authority_weights(df)
  inv_freq_wts <- rep(1.0, n)
  inv_freq_wts[anchor_idx] <- compute_inverse_class_frequency_weights(y_anchor)
  case_wts <- auth_wts * inv_freq_wts * spatial_balance_weight
  case_wts <- case_wts / max(mean(case_wts[anchor_idx], na.rm = TRUE), 1e-8)

  X_env_anc <- X_env[anchor_idx, , drop = FALSE]
  X_all_anc <- X_all[anchor_idx, , drop = FALSE]
  wts_anc   <- case_wts[anchor_idx]

  if (requireNamespace("ranger", quietly = TRUE)) {
    df_env_anc <- data.frame(y = y_anchor, as.data.frame(X_env_anc), check.names = FALSE)
    df_all_anc <- data.frame(y = y_anchor, as.data.frame(X_all_anc), check.names = FALSE)
    rf_nonspatial <- ranger::ranger(y ~ ., data = df_env_anc, num.trees = n_trees, importance = "permutation", probability = FALSE, case.weights = wts_anc, seed = 42)
    oob_acc_nonspatial <- 1 - rf_nonspatial$prediction.error
    rf_spatial <- ranger::ranger(y ~ ., data = df_all_anc, num.trees = n_trees, importance = "permutation", probability = FALSE, case.weights = wts_anc, seed = 42)
    oob_acc_spatial <- 1 - rf_spatial$prediction.error
    rf_prob_anchor <- ranger::ranger(y ~ ., data = df_env_anc, num.trees = n_trees, probability = TRUE, case.weights = wts_anc, seed = 42)
    imp_vec <- ranger::importance(rf_spatial)
    imp_mat <- cbind(MeanDecreaseAccuracy = imp_vec)
  } else {
    rf_nonspatial <- randomForest::randomForest(x = X_env_anc, y = y_anchor, ntree = n_trees, importance = TRUE)
    oob_acc_nonspatial <- 1 - rf_nonspatial$err.rate[n_trees, "OOB"]
    rf_spatial <- randomForest::randomForest(x = X_all_anc, y = y_anchor, ntree = n_trees, importance = TRUE)
    oob_acc_spatial <- 1 - rf_spatial$err.rate[n_trees, "OOB"]
    rf_prob_anchor <- rf_nonspatial
    imp_mat <- randomForest::importance(rf_spatial, type = 1)
  }

  # Estimate predictor spatial autocorrelation range for Spatial Block Cross-Validation (Roberts et al. 2017)
  autocorr_info <- estimate_spatial_autocorrelation_block_size(coords[anchor_idx, , drop = FALSE], X_env_anc)
  n_anc <- length(anchor_idx)
  fold_ids_anc <- NULL
  if (requireNamespace("blockCV", quietly = TRUE) && requireNamespace("sf", quietly = TRUE) && n_anc >= 10) {
    fold_ids_anc <- tryCatch({
      pts_sf <- sf::st_as_sf(data.frame(longitude = coords[anchor_idx, 1], latitude = coords[anchor_idx, 2], y = y_anchor),
                             coords = c("longitude", "latitude"), crs = 4326)
      sb <- tryCatch(
        suppressMessages(suppressWarnings(
          blockCV::cv_spatial(x = pts_sf, column = "y", size = autocorr_info$range_m, k = 5, selection = "random",
                              iteration = 10, progress = FALSE, report = FALSE, plot = FALSE)
        )),
        error = function(e) {
          suppressMessages(suppressWarnings(
            blockCV::cv_spatial(x = pts_sf, column = "y", k = 5, selection = "random",
                                iteration = 10, progress = FALSE, report = FALSE, plot = FALSE)
          ))
        }
      )
      sb$folds_ids
    }, error = function(e) NULL)
  }
  if (is.null(fold_ids_anc) || length(unique(stats::na.omit(fold_ids_anc))) < 2) {
    fold_ids_anc <- assign_spatial_grid_blocks(coords[anchor_idx, , drop = FALSE], block_deg = autocorr_info$range_deg, k = min(5L, max(2L, floor(n_anc / 2))))
  }

  # Generate Out-of-Fold (OOF) Ecological Probabilities for Gold/Silver Anchors via Spatial Block CV
  cls_levels <- levels(y_anchor)
  oof_prob_mat <- matrix(NA_real_, nrow = n_anc, ncol = length(cls_levels), dimnames = list(NULL, cls_levels))
  oof_preds_anc <- rep(NA_character_, n_anc)

  for (f_id in sort(unique(stats::na.omit(fold_ids_anc)))) {
    tr_sub <- which(fold_ids_anc != f_id)
    te_sub <- which(fold_ids_anc == f_id)
    if (length(tr_sub) >= 4 && length(te_sub) >= 1 && length(unique(y_anchor[tr_sub])) >= 2) {
      if (requireNamespace("ranger", quietly = TRUE)) {
        tr_df <- data.frame(y = droplevels(y_anchor[tr_sub]), as.data.frame(X_env_anc[tr_sub, , drop = FALSE]), check.names = FALSE)
        te_df <- as.data.frame(X_env_anc[te_sub, , drop = FALSE])
        fit_cv <- ranger::ranger(y ~ ., data = tr_df, num.trees = max(100, floor(n_trees / 2)),
                                 probability = TRUE, case.weights = wts_anc[tr_sub], seed = 42 + as.integer(f_id))
        pr_mat <- stats::predict(fit_cv, data = te_df)$predictions
        if (!is.matrix(pr_mat)) pr_mat <- matrix(pr_mat, nrow = length(te_sub), dimnames = list(NULL, names(pr_mat)))
        for (cl_nm in intersect(colnames(pr_mat), cls_levels)) {
          oof_prob_mat[te_sub, cl_nm] <- pr_mat[, cl_nm]
        }
        oof_preds_anc[te_sub] <- colnames(pr_mat)[max.col(pr_mat, ties.method = "first")]
      } else {
        fit_cv <- randomForest::randomForest(x = X_env_anc[tr_sub, , drop = FALSE], y = droplevels(y_anchor[tr_sub]), ntree = max(100, floor(n_trees / 2)))
        pr_mat <- stats::predict(fit_cv, newdata = X_env_anc[te_sub, , drop = FALSE], type = "prob")
        for (cl_nm in intersect(colnames(pr_mat), cls_levels)) {
          oof_prob_mat[te_sub, cl_nm] <- pr_mat[, cl_nm]
        }
        oof_preds_anc[te_sub] <- as.character(stats::predict(fit_cv, newdata = X_env_anc[te_sub, , drop = FALSE]))
      }
    }
  }

  # Fill any unassigned OOF rows (e.g. if a tiny fold had only 1 class) from OOB predictions of rf_prob_anchor
  missing_oof <- which(is.na(oof_preds_anc))
  if (length(missing_oof) > 0 && requireNamespace("ranger", quietly = TRUE)) {
    oob_p <- rf_prob_anchor$predictions
    if (is.matrix(oob_p)) {
      for (cl_nm in intersect(colnames(oob_p), cls_levels)) {
        oof_prob_mat[missing_oof, cl_nm] <- oob_p[missing_oof, cl_nm]
      }
      oof_preds_anc[missing_oof] <- colnames(oob_p)[max.col(oob_p[missing_oof, , drop = FALSE], ties.method = "first")]
    }
  }

  ok_cv <- !is.na(oof_preds_anc) & !is.na(y_anchor)
  spatial_block_cv_acc <- if (any(ok_cv)) {
    sum(wts_anc[ok_cv] * (oof_preds_anc[ok_cv] == as.character(y_anchor[ok_cv]))) / sum(wts_anc[ok_cv])
  } else {
    oob_acc_spatial
  }
  if (is.na(spatial_block_cv_acc) || !is.finite(spatial_block_cv_acc)) spatial_block_cv_acc <- oob_acc_spatial

  # Assemble full specimen-level ecological probabilities:
  # - Out-of-fold (OOF) spatial block CV probabilities for Tier_1_Gold + Tier_2_Silver anchors
  # - Passive projection via predict(rf_prob_anchor) for Tier_3_Bronze & high-coordinate-uncertainty records
  full_prob_mat <- matrix(NA_real_, nrow = n, ncol = length(cls_levels), dimnames = list(NULL, cls_levels))
  full_prob_mat[anchor_idx, ] <- oof_prob_mat
  full_preds <- rep(NA_character_, n)
  full_preds[anchor_idx] <- oof_preds_anc
  full_block_id <- rep(NA_integer_, n)
  full_block_id[anchor_idx] <- as.integer(fold_ids_anc)

  if (length(passive_idx) > 0) {
    if (requireNamespace("ranger", quietly = TRUE)) {
      pass_df <- as.data.frame(X_env[passive_idx, , drop = FALSE])
      pass_pr <- stats::predict(rf_prob_anchor, data = pass_df)$predictions
      if (!is.matrix(pass_pr)) pass_pr <- matrix(pass_pr, nrow = length(passive_idx), dimnames = list(NULL, names(pass_pr)))
      for (cl_nm in intersect(colnames(pass_pr), cls_levels)) {
        full_prob_mat[passive_idx, cl_nm] <- pass_pr[, cl_nm]
      }
      full_preds[passive_idx] <- colnames(pass_pr)[max.col(pass_pr, ties.method = "first")]
    } else {
      pass_pr <- stats::predict(rf_prob_anchor, newdata = X_env[passive_idx, , drop = FALSE], type = "prob")
      for (cl_nm in intersect(colnames(pass_pr), cls_levels)) {
        full_prob_mat[passive_idx, cl_nm] <- pass_pr[, cl_nm]
      }
      full_preds[passive_idx] <- as.character(stats::predict(rf_prob_anchor, newdata = X_env[passive_idx, , drop = FALSE]))
    }
  }

  rf_pred_prob <- apply(full_prob_mat, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
  rf_given_prob <- sapply(seq_len(n), function(i) {
    g_i <- as.character(y[i])
    if (!is.na(g_i) && g_i %in% colnames(full_prob_mat)) full_prob_mat[i, g_i] else NA_real_
  })
  rf_pred_source <- ifelse(!elig_coord, "Passive_Projection_High_Uncertainty",
                    ifelse(seq_len(n) %in% anchor_idx, "Spatial_Block_OOF_Anchor", "Passive_Projection_Bronze"))

  # Spatial vs Environmental Variance Partitioning (vegan::varpart across [A] Fine-Grain Edaphic & Outcrop Lithology, [B] Macro-Climate, [C] Spatial MEMs)
  pred_names <- rownames(imp_mat)
  is_mem_pred <- grepl("^MEM", pred_names)
  is_fine_edaphic_pred <- grepl("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)", pred_names)
  is_climate_pred <- grepl("^(bio|delta_doy|doy)", pred_names, ignore.case = TRUE)

  imp_df <- data.frame(
    Predictor = pred_names,
    MeanDecreaseAccuracy = imp_mat[, 1],
    PositiveImportance = pmax(0, imp_mat[, 1]),
    Type = ifelse(is_mem_pred, "Spatial (MEMs)",
           ifelse(is_fine_edaphic_pred, "Fine-Grain Edaphic & Outcrop (10-30m)", "Environmental")),
    Group = ifelse(
      is_mem_pred, "[C] Spatial MEMs",
      ifelse(
        is_fine_edaphic_pred, "[A] Fine-Grain Edaphic & Outcrop Lithology",
        ifelse(is_climate_pred, "[B] Macro-Climate", "[A] Fine-Grain Edaphic & Outcrop Lithology")
      )
    ),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]

  fine_edaphic_vars <- intersect(env_vars, grep("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)", env_vars, value = TRUE))
  macro_climate_vars <- intersect(env_vars, grep("^(bio|delta_doy|doy)", env_vars, value = TRUE, ignore.case = TRUE))

  varpart_obj <- NULL
  if (requireNamespace("vegan", quietly = TRUE) && length(fine_edaphic_vars) >= 1 && length(macro_climate_vars) >= 1) {
    varpart_obj <- tryCatch({
      Y_mat <- stats::model.matrix(~ y_anchor - 1)
      vegan::varpart(
        Y_mat,
        as.data.frame(X_env_anc[, fine_edaphic_vars, drop = FALSE]),
        as.data.frame(X_env_anc[, macro_climate_vars, drop = FALSE]),
        as.data.frame(MEMs[anchor_idx, , drop = FALSE])
      )
    }, error = function(e) NULL)
  }

  grp_names <- c("[A] Fine-Grain Edaphic & Outcrop Lithology", "[B] Macro-Climate", "[C] Spatial MEMs")
  grp_sums <- sapply(grp_names, function(g) {
    sum(imp_df$PositiveImportance[imp_df$Group == g], na.rm = TRUE)
  })
  tot_pos <- sum(grp_sums, na.rm = TRUE)

  if (!is.null(varpart_obj) && !is.null(varpart_obj$part$indfract)) {
    adj_col <- intersect(c("Adj.R.square", "Adj.R.sq"), colnames(varpart_obj$part$indfract))[1]
    if (!is.na(adj_col)) {
      vp_pos <- pmax(0, as.numeric(varpart_obj$part$indfract[1:3, adj_col]))
      if (sum(vp_pos, na.rm = TRUE) > 0) {
        grp_sums <- vp_pos
        names(grp_sums) <- grp_names
        tot_pos <- sum(grp_sums, na.rm = TRUE)
      }
    }
  }

  if (tot_pos <= 0) {
    var_pct <- round(rep(100 / length(grp_names), length(grp_names)), 2)
  } else {
    var_pct <- round((as.numeric(grp_sums) / tot_pos) * 100, 2)
    var_pct[which.max(var_pct)] <- round(var_pct[which.max(var_pct)] + (100.0 - sum(var_pct)), 2)
  }

  var_part <- data.frame(
    Component = grp_names,
    Variance_Pct = var_pct,
    stringsAsFactors = FALSE
  )

  # Moran's I spatial autocorrelation check on model residuals (Dormann et al. 2013; Benito 2021)
  oob_pred_sp <- if (requireNamespace("ranger", quietly = TRUE)) as.character(rf_spatial$predictions) else as.character(rf_spatial$predicted)
  res_vec <- as.numeric(oob_pred_sp != as.character(y_anchor)) - mean(oob_pred_sp != as.character(y_anchor), na.rm = TRUE)
  res_vec[!is.finite(res_vec)] <- 0.0
  moran_i <- 0.0; moran_p <- 1.0
  if (n_anc >= 5 && sum(res_vec^2) > 1e-10) {
    D_sp <- as.matrix(stats::dist(coords[anchor_idx, , drop = FALSE])); W_sp <- 1.0 / (D_sp + 0.05); diag(W_sp) <- 0.0
    w_sum <- sum(W_sp)
    if (w_sum > 1e-8) {
      moran_i <- as.numeric((n_anc / w_sum) * (t(res_vec) %*% W_sp %*% res_vec) / sum(res_vec^2))
      e_i <- -1.0 / (n_anc - 1)
      z_m <- (moran_i - e_i) / max(sqrt(2.0 / max(n_anc, 4)), 1e-6)
      moran_p <- 2.0 * stats::pnorm(-abs(z_m))
    }
  }
  message(sprintf(
    "Authority-Anchored Spatial RF (%d Gold/Silver anchors, %d passive projections | Block Range = %.1f km): OOB_Acc = %.2f%% | Spatial_Block_CV_Acc = %.2f%% | Residual Moran's I = %.4f (p = %.4f)",
    n_anc, length(passive_idx), autocorr_info$range_km, oob_acc_spatial * 100, spatial_block_cv_acc * 100, moran_i, moran_p
  ))
  return(list(
    rf_spatial = rf_spatial, rf_nonspatial = rf_nonspatial, rf_prob_anchor = rf_prob_anchor,
    importance = imp_df, variance_partition = var_part, varpart_model = varpart_obj, mems = MEMs,
    OOB_Accuracy = oob_acc_spatial, OOB_Error = 1 - oob_acc_spatial,
    NonSpatial_OOB_Accuracy = oob_acc_nonspatial,
    Spatial_Block_CV_Accuracy = spatial_block_cv_acc,
    spatial_autocorr_range_km = autocorr_info$range_km,
    rf_predicted_taxon = full_preds,
    rf_predicted_prob = round(rf_pred_prob, 4),
    rf_given_taxon_prob = round(rf_given_prob, 4),
    rf_prediction_source = rf_pred_source,
    spatial_block_id = full_block_id,
    collinearity_screening = col_screen,
    residual_moran_i = round(moran_i, 4), residual_moran_p = round(moran_p, 4)
  ))
}
run_spatial_random_forests <- compute_spatial_rf_mems
run_spatial_rf_models <- compute_spatial_rf_mems

# ------------------------------------------------------------------------------
# 5. Partitioned PCA-env Niche Overlap: Climatic vs. Edaphic-Geologic vs. Joint
#    (Broennimann et al. 2012; Warren et al. 2008; Dormann et al. 2013)
# ------------------------------------------------------------------------------
compute_pca_env_partition_overlap <- function(bg_mat, grp_col, sp1, sp2, n_perm = 100) {
  if (is.null(bg_mat) || ncol(bg_mat) == 0 || nrow(bg_mat) < 6) {
    return(list(D_emp = NA_real_, I_emp = NA_real_, p_equiv_D = NA_real_, p_equiv_I = NA_real_,
                p_bg_D = NA_real_, p_bg_I = NA_real_, D_null = rep(NA_real_, n_perm), I_null = rep(NA_real_, n_perm)))
  }
  for (j in seq_len(ncol(bg_mat))) {
    col_v <- bg_mat[, j]
    if (any(!is.finite(col_v))) {
      med_j <- stats::median(col_v[is.finite(col_v)], na.rm = TRUE)
      if (!is.finite(med_j)) med_j <- 0.0
      bg_mat[!is.finite(col_v), j] <- med_j
    }
  }
  active_cols <- which(apply(bg_mat, 2, stats::sd, na.rm = TRUE) > 1e-8)
  if (length(active_cols) == 0) {
    return(list(D_emp = NA_real_, I_emp = NA_real_, p_equiv_D = NA_real_, p_equiv_I = NA_real_,
                p_bg_D = NA_real_, p_bg_I = NA_real_, D_null = rep(NA_real_, n_perm), I_null = rep(NA_real_, n_perm)))
  }
  bg_sub <- bg_mat[, active_cols, drop = FALSE]
  scores_glob <- tryCatch({
    if (ncol(bg_sub) >= 2) {
      as.data.frame(stats::prcomp(bg_sub, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE])
    } else {
      z1 <- as.numeric(scale(bg_sub[, 1]))
      data.frame(PC1 = z1, PC2 = 0.05 * seq_along(z1) / length(z1))
    }
  }, error = function(e) NULL)

  idx1 <- which(grp_col == sp1)
  idx2 <- which(grp_col == sp2)
  n1 <- length(idx1); n2 <- length(idx2)

  ecospat_ok <- FALSE
  if (requireNamespace("ecospat", quietly = TRUE) && !is.null(scores_glob) && n1 >= 5 && n2 >= 5) {
    eco_fit <- tryCatch({
      scores_sp1 <- scores_glob[idx1, 1:2, drop = FALSE]
      scores_sp2 <- scores_glob[idx2, 1:2, drop = FALSE]
      grid_R <- if (n_perm < 99) 50 else 100
      z1 <- ecospat::ecospat.grid.clim.dyn(glob = scores_glob, glob1 = scores_glob, sp = scores_sp1, R = grid_R)
      z2 <- ecospat::ecospat.grid.clim.dyn(glob = scores_glob, glob1 = scores_glob, sp = scores_sp2, R = grid_R)
      ov <- ecospat::ecospat.niche.overlap(z1, z2, cor = TRUE)
      if (n_perm >= 99) {
        eq <- tryCatch(
          suppressMessages(suppressWarnings(ecospat::ecospat.niche.equivalency.test(z1, z2, rep = n_perm, overlap.alternative = "lower", ncores = 1))),
          error = function(e) suppressMessages(suppressWarnings(ecospat::ecospat.niche.equivalency.test(z1, z2, rep = n_perm, ncores = 1)))
        )
        bg_sim <- tryCatch(
          suppressMessages(suppressWarnings(ecospat::ecospat.niche.similarity.test(z1, z2, rep = n_perm, overlap.alternative = "lower", rand.type = 2, ncores = 1))),
          error = function(e) tryCatch(
            suppressMessages(suppressWarnings(ecospat::ecospat.niche.similarity.test(z1, z2, rep = n_perm, rand.type = 2, ncores = 1))),
            error = function(e2) NULL
          )
        )
        d_sim <- as.numeric(eq$sim$D); i_sim <- as.numeric(eq$sim$I)
        p_bg_d <- if (!is.null(bg_sim) && !is.null(bg_sim$sim$D)) {
          d_bg <- as.numeric(bg_sim$sim$D)
          (sum(d_bg <= as.numeric(ov$D), na.rm = TRUE) + 1) / (length(d_bg) + 1)
        } else {
          NA_real_
        }
        p_bg_i <- if (!is.null(bg_sim) && !is.null(bg_sim$sim$I)) {
          i_bg <- as.numeric(bg_sim$sim$I)
          (sum(i_bg <= as.numeric(ov$I), na.rm = TRUE) + 1) / (length(i_bg) + 1)
        } else {
          NA_real_
        }
        list(
          D_emp = as.numeric(ov$D),
          I_emp = as.numeric(ov$I),
          p_equiv_D = (sum(d_sim <= as.numeric(ov$D), na.rm = TRUE) + 1) / (length(d_sim) + 1),
          p_equiv_I = (sum(i_sim <= as.numeric(ov$I), na.rm = TRUE) + 1) / (length(i_sim) + 1),
          p_bg_D = p_bg_d,
          p_bg_I = p_bg_i,
          D_null = d_sim,
          I_null = i_sim,
          has_perm = TRUE
        )
      } else {
        list(
          D_emp = as.numeric(ov$D),
          I_emp = as.numeric(ov$I),
          has_perm = FALSE
        )
      }
    }, error = function(e) NULL)
    if (!is.null(eco_fit) && is.finite(eco_fit$D_emp) && is.finite(eco_fit$I_emp)) {
      D_emp <- eco_fit$D_emp; I_emp <- eco_fit$I_emp
      if (isTRUE(eco_fit$has_perm)) {
        p_equiv_D <- eco_fit$p_equiv_D; p_equiv_I <- eco_fit$p_equiv_I
        p_bg_D <- eco_fit$p_bg_D; p_bg_I <- eco_fit$p_bg_I
        D_null <- eco_fit$D_null; I_null <- eco_fit$I_null
        ecospat_ok <- TRUE
      }
    }
  }

  # Calibrated 2D PCA-env Kernel Density & Equivalency / Background Similarity Tests (Broennimann et al. 2012; Warren et al. 2008)
  ord_bg <- if (!is.null(scores_glob)) as.matrix(scores_glob) else as.matrix(bg_sub)
  ord_df1 <- ord_bg[idx1, , drop = FALSE]
  ord_df2 <- ord_bg[idx2, , drop = FALSE]
  m1 <- colMeans(ord_df1); s1 <- apply(ord_df1, 2, stats::sd) + 1e-4
  m2 <- colMeans(ord_df2); s2 <- apply(ord_df2, 2, stats::sd) + 1e-4
  p1 <- exp(-0.5 * rowSums(t((t(ord_bg) - m1) / s1)^2)); p1 <- p1 / max(sum(p1), 1e-12)
  p2 <- exp(-0.5 * rowSums(t((t(ord_bg) - m2) / s2)^2)); p2 <- p2 / max(sum(p2), 1e-12)

  if (!ecospat_ok) {
    if (!exists("D_emp", inherits = FALSE) || !is.finite(D_emp)) {
      D_emp <- 1 - 0.5 * sum(abs(p1 - p2))
      I_emp <- 1 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(p2))^2))
    }

    pooled <- rbind(ord_df1, ord_df2); n_tot <- n1 + n2
    D_null <- numeric(n_perm); I_null <- numeric(n_perm)
    for (k in seq_len(n_perm)) {
      idx_perm <- sample.int(n_tot, n1)
      p_df1 <- pooled[idx_perm, , drop = FALSE]; p_df2 <- pooled[-idx_perm, , drop = FALSE]
      pm1 <- colMeans(p_df1); ps1 <- apply(p_df1, 2, stats::sd) + 1e-4
      pm2 <- colMeans(p_df2); ps2 <- apply(p_df2, 2, stats::sd) + 1e-4
      pp1 <- exp(-0.5 * rowSums(t((t(ord_bg) - pm1) / ps1)^2)); pp1 <- pp1 / max(sum(pp1), 1e-12)
      pp2 <- exp(-0.5 * rowSums(t((t(ord_bg) - pm2) / ps2)^2)); pp2 <- pp2 / max(sum(pp2), 1e-12)
      D_null[k] <- 1 - 0.5 * sum(abs(pp1 - pp2))
      I_null[k] <- 1 - 0.5 * sqrt(sum((sqrt(pp1) - sqrt(pp2))^2))
    }
    p_equiv_D <- (sum(D_null <= D_emp) + 1) / (n_perm + 1)
    p_equiv_I <- (sum(I_null <= I_emp) + 1) / (n_perm + 1)
    p_bg_D <- NA_real_; p_bg_I <- NA_real_
  }

  # Background Similarity Permutation Test (Warren et al. 2008: random background niche shifts across regional PCA-env)
  if (is.na(p_bg_D) || !is.finite(p_bg_D)) {
    n_bg_tot <- nrow(ord_bg)
    D_bg_null <- numeric(n_perm); I_bg_null <- numeric(n_perm)
    for (k in seq_len(n_perm)) {
      bg_idx2 <- sample.int(n_bg_tot, size = min(n2, n_bg_tot), replace = TRUE)
      bg_df2 <- ord_bg[bg_idx2, , drop = FALSE]
      bm2 <- colMeans(bg_df2); bs2 <- apply(bg_df2, 2, stats::sd) + 1e-4
      bp2 <- exp(-0.5 * rowSums(t((t(ord_bg) - bm2) / bs2)^2)); bp2 <- bp2 / max(sum(bp2), 1e-12)
      D_bg_null[k] <- 1 - 0.5 * sum(abs(p1 - bp2))
      I_bg_null[k] <- 1 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(bp2))^2))
    }
    p_bg_D <- (sum(D_bg_null <= D_emp) + 1) / (n_perm + 1)
    p_bg_I <- (sum(I_bg_null <= I_emp) + 1) / (n_perm + 1)
  }

  return(list(
    D_emp = D_emp, I_emp = I_emp,
    p_equiv_D = p_equiv_D, p_equiv_I = p_equiv_I,
    p_bg_D = p_bg_D, p_bg_I = p_bg_I,
    D_null = D_null, I_null = I_null
  ))
}

run_warren_niche_identity_tests <- function(df, env_vars, n_perm = 100, intraspecific_mode = FALSE) {
  message(sprintf("Executing Partitioned PCA-env Niche Equivalency & Background Similarity Tests (%d permutations per pair)...", n_perm))
  if (!("analysis_group" %in% names(df)) || isTRUE(intraspecific_mode)) {
    df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
  }
  # Exclude records with coordinateUncertaintyInMeters > 10000 from calibrating PCA-env niche envelopes
  elig_mask <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, nrow(df))
  elig_mask[is.na(elig_mask)] <- TRUE
  if ("coordinateUncertaintyInMeters" %in% names(df)) {
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    elig_mask <- elig_mask & (is.na(u_m) | u_m <= 10000.0)
  }
  df_calib <- if (sum(elig_mask) >= 10) df[elig_mask, , drop = FALSE] else df

  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df_calib$species_standardized))) < 2)
  grp_col <- if (is_intra) df_calib$analysis_group else df_calib$species_standardized
  cand_groups <- if (is_intra) sort(unique(stats::na.omit(grp_col))) else get_active_packera_taxa(grp_col)
  if (length(cand_groups) < 2) cand_groups <- sort(unique(stats::na.omit(df_calib$analysis_group)))
  pairs <- if (length(cand_groups) >= 2) combn(cand_groups, 2, simplify = FALSE) else list()
  results_list <- list(); perm_dist_list <- list()

  env_vars <- intersect(env_vars, names(df_calib))
  X_raw <- as.matrix(df_calib[, env_vars, drop = FALSE])
  block_screen <- if (exists("screen_collinear_predictors_by_block", mode = "function")) {
    screen_collinear_predictors_by_block(X_raw, cor.threshold = 0.70, vif.threshold = 5.0)
  } else {
    clim_c <- grep("^(bio|delta_doy|doy)", env_vars, value = TRUE, ignore.case = TRUE)
    edaph_c <- setdiff(env_vars, clim_c)
    list(climatic_variables = clim_c, edaphic_variables = edaph_c, selected_variables = env_vars)
  }

  clim_vars <- intersect(block_screen$climatic_variables, colnames(X_raw))
  if (length(clim_vars) == 0) clim_vars <- intersect(grep("^(bio|delta_doy|doy)", colnames(X_raw), value = TRUE, ignore.case = TRUE), colnames(X_raw))
  edaph_vars <- intersect(block_screen$edaphic_variables, colnames(X_raw))
  if (length(edaph_vars) == 0) edaph_vars <- setdiff(colnames(X_raw), clim_vars)
  joint_vars <- unique(c(clim_vars, edaph_vars))

  for (pr in pairs) {
    sp1 <- pr[1]; sp2 <- pr[2]; pair_name <- paste0(gsub("Packera ", "P.", sp1), " vs ", gsub("Packera ", "P.", sp2))
    n1 <- sum(grp_col == sp1, na.rm = TRUE)
    n2 <- sum(grp_col == sp2, na.rm = TRUE)
    if (n1 < 5 || n2 < 5) next

    # (a) Climatic Niche Partition (D_clim, I_clim)
    clim_res <- compute_pca_env_partition_overlap(
      if (length(clim_vars) > 0) X_raw[, clim_vars, drop = FALSE] else X_raw,
      grp_col, sp1, sp2, n_perm = n_perm
    )
    # (b) Edaphic-Geologic Niche Partition (D_edaph, I_edaph)
    edaph_res <- compute_pca_env_partition_overlap(
      if (length(edaph_vars) > 0) X_raw[, edaph_vars, drop = FALSE] else X_raw,
      grp_col, sp1, sp2, n_perm = n_perm
    )
    # (c) Joint Climato-Edaphic Niche Partition
    joint_res <- compute_pca_env_partition_overlap(
      X_raw[, joint_vars, drop = FALSE],
      grp_col, sp1, sp2, n_perm = n_perm
    )

    D_emp <- joint_res$D_emp; I_emp <- joint_res$I_emp
    p_val_D <- joint_res$p_equiv_D; p_val_I <- joint_res$p_equiv_I
    D_clim <- clim_res$D_emp; I_clim <- clim_res$I_emp
    D_edaph <- edaph_res$D_emp; I_edaph <- edaph_res$I_emp

    div_axis <- ifelse(!is.na(D_edaph) && !is.na(D_clim) && (D_edaph + 0.08 < D_clim), "Edaphic_Lithological_Divergence",
                ifelse(!is.na(D_clim) && !is.na(D_edaph) && (D_clim + 0.08 < D_edaph), "Macroclimatic_Divergence",
                       "Joint_Climato_Edaphic"))

    min_p <- min(c(p_val_D, edaph_res$p_equiv_D, clim_res$p_equiv_D), na.rm = TRUE)
    outcome <- ifelse(min_p < 0.01, "Distinct Ecological Niche (Species Boundary Confirmed)",
               ifelse(min_p < 0.05, "Significant Divergence (Subspecific Ecological Race)",
                      "Conserved / Overlapping Niche (Ecophenotypic Variant)"))

    results_list[[pair_name]] <- data.frame(
      Pair = pair_name, Taxon1 = sp1, Taxon2 = sp2, N1 = n1, N2 = n2,
      D_clim = round(D_clim, 4), I_clim = round(I_clim, 4),
      p_equiv_D_clim = round(clim_res$p_equiv_D, 4), p_bg_D_clim = round(clim_res$p_bg_D, 4),
      D_edaph = round(D_edaph, 4), I_edaph = round(I_edaph, 4),
      p_equiv_D_edaph = round(edaph_res$p_equiv_D, 4), p_bg_D_edaph = round(edaph_res$p_bg_D, 4),
      D_joint = round(D_emp, 4), I_joint = round(I_emp, 4),
      Schoener_D = round(D_emp, 4), Warren_D = round(D_emp, 4), Warren_I = round(I_emp, 4),
      p_value_D = round(p_val_D, 4), p_value_I = round(p_val_I, 4),
      p_bg_D_joint = round(joint_res$p_bg_D, 4), p_bg_I_joint = round(joint_res$p_bg_I, 4),
      Niche_Divergence_Axis = div_axis,
      Taxonomic_Inference = outcome, stringsAsFactors = FALSE
    )
    perm_dist_list[[pair_name]] <- data.frame(
      Pair = pair_name,
      D_null = joint_res$D_null, I_null = joint_res$I_null,
      D_clim_null = clim_res$D_null, D_edaph_null = edaph_res$D_null
    )
  }

  summary_table <- do.call(rbind, results_list)
  perm_table <- do.call(rbind, perm_dist_list)
  return(list(summary = summary_table, permutations = perm_table, collinearity_blocks = block_screen))
}
run_niche_identity_tests <- run_warren_niche_identity_tests

# ------------------------------------------------------------------------------
# 6. Multi-Panel Publication Figure & Summary Export
# ------------------------------------------------------------------------------
select_top_rf_predictors <- function(imp_df, plot_df) {
  avail_cols <- names(plot_df)
  edaph_regex <- "^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)"

  default_1 <- if ("fine_soil_sand" %in% avail_cols) "fine_soil_sand" else if ("soil_sand" %in% avail_cols) "soil_sand" else avail_cols[1]
  default_2 <- if ("fine_soil_ph" %in% avail_cols) "fine_soil_ph" else if ("soil_ph" %in% avail_cols) "soil_ph" else avail_cols[min(2L, length(avail_cols))]

  if (is.null(imp_df) || !is.data.frame(imp_df) || nrow(imp_df) == 0L || !("Predictor" %in% names(imp_df))) {
    return(list(top_env_1 = default_1, top_env_2 = default_2, top_edaph_1 = default_1))
  }

  mda <- if ("MeanDecreaseAccuracy" %in% names(imp_df)) {
    suppressWarnings(as.numeric(imp_df$MeanDecreaseAccuracy))
  } else if ("Importance" %in% names(imp_df)) {
    suppressWarnings(as.numeric(imp_df$Importance))
  } else {
    rep(0.0, nrow(imp_df))
  }
  mda[!is.finite(mda)] <- -Inf

  preds <- as.character(imp_df$Predictor)
  is_non_mem <- !grepl("^MEM", preds)
  if ("Group" %in% names(imp_df)) {
    is_non_mem <- is_non_mem & (imp_df$Group %in% c("[A] Fine-Grain Edaphic & Outcrop Lithology", "[B] Macro-Climate"))
  }
  in_df <- preds %in% avail_cols & vapply(preds, function(p) {
    if (!(p %in% avail_cols)) return(FALSE)
    v <- suppressWarnings(as.numeric(plot_df[[p]]))
    any(is.finite(v))
  }, logical(1))

  keep_idx <- which(is_non_mem & in_df)
  if (length(keep_idx) == 0L) {
    return(list(top_env_1 = default_1, top_env_2 = default_2, top_edaph_1 = default_1))
  }

  ord_env <- keep_idx[order(mda[keep_idx], decreasing = TRUE)]
  top_env_1 <- preds[ord_env[1]]
  if (length(ord_env) >= 2L) {
    top_env_2 <- preds[ord_env[2]]
  } else {
    cand_2 <- setdiff(c(default_2, default_1, intersect(c("bio1_temp_mean", "fine_bedrock_depth_cm"), avail_cols)), top_env_1)
    top_env_2 <- if (length(cand_2) > 0L) cand_2[1] else top_env_1
  }

  edaph_mask <- if ("Group" %in% names(imp_df)) {
    imp_df$Group[ord_env] == "[A] Fine-Grain Edaphic & Outcrop Lithology"
  } else {
    grepl(edaph_regex, preds[ord_env])
  }
  edaph_ranked <- preds[ord_env[edaph_mask]]
  top_edaph_1 <- if (length(edaph_ranked) > 0L) edaph_ranked[1] else top_env_1

  list(top_env_1 = top_env_1, top_env_2 = top_env_2, top_edaph_1 = top_edaph_1)
}

export_dynamic_fig06_ecological_niche_phenology <- function(df, srf_res = NULL, niche_res = NULL, out_plot_pdf = NULL) {
  fdir <- if (exists("resolve_active_figures_dir", mode = "function")) {
    resolve_active_figures_dir(reference_path = out_plot_pdf)
  } else if (!is.null(out_plot_pdf) && nzchar(out_plot_pdf)) {
    dirname(out_plot_pdf)
  } else {
    "outputs/figures"
  }
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  tx_clean <- if (exists("resolve_plot_group_column", mode = "function")) {
    resolve_plot_group_column(df)
  } else if ("analysis_group" %in% names(df) && length(unique(stats::na.omit(df$analysis_group))) >= 2L) {
    as.character(df$analysis_group)
  } else {
    as.character(df$species_standardized)
  }
  present_tx <- sort(unique(tx_clean[!is.na(tx_clean) & nzchar(tx_clean) & tx_clean != "Unknown"]))
  if (length(present_tx) == 0L) present_tx <- get_active_packera_taxa(df$species_standardized)
  if (length(present_tx) == 0L) return(invisible(NULL))

  dyn_pal <- resolve_dynamic_group_palette(present_tx)
  hollow_pchs <- c(1L, 2L, 5L, 0L, 6L, 4L, 3L, 8L)
  dyn_pch <- if (exists("get_dynamic_packera_shapes", mode = "function")) {
    sh <- get_dynamic_packera_shapes(present_tx)
    if (any(is.na(sh[present_tx]))) stats::setNames(rep_len(hollow_pchs, length(present_tx)), present_tx) else sh
  } else {
    stats::setNames(rep_len(hollow_pchs, length(present_tx)), present_tx)
  }
  alpha_val <- if (nrow(df) > 1000L) 0.45 else 0.78
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = alpha_val)
  pch_tx  <- dyn_pch[tx_clean]; pch_tx[is.na(pch_tx)] <- 1L

  if (!("fine_soil_sand" %in% names(df)) && "soil_sand" %in% names(df)) df$fine_soil_sand <- df$soil_sand
  if (!("fine_soil_ph" %in% names(df)) && "soil_ph" %in% names(df)) df$fine_soil_ph <- df$soil_ph

  imp_df <- if (!is.null(srf_res) && !is.null(srf_res$importance)) srf_res$importance else NULL
  top_sel <- select_top_rf_predictors(imp_df, df)
  top_env_1 <- top_sel$top_env_1
  top_env_2 <- top_sel$top_env_2
  top_edaph_1 <- top_sel$top_edaph_1
  top_env_lbls <- format_step06_predictor_labels(c(top_env_1, top_env_2))
  top_edaph_lbl <- format_step06_predictor_labels(top_edaph_1)[1]

  render_plate <- function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.5, 11.0, 3.0, 1.2), oma = c(0.5, 0.5, 1.8, 0.5))

    # Panel A: Spatial RF Predictor Importance Ranking
    pred_raw <- c(top_env_2, top_env_1)
    imp_vals <- c(0.15, 0.25)
    if (!is.null(imp_df) && nrow(imp_df) > 0L && all(c("Predictor", "MeanDecreaseAccuracy") %in% names(imp_df))) {
      ord <- utils::tail(order(as.numeric(imp_df$MeanDecreaseAccuracy), decreasing = FALSE), 12)
      pred_raw <- as.character(imp_df$Predictor[ord])
      imp_vals <- pmax(0, as.numeric(imp_df$MeanDecreaseAccuracy[ord]))
    }
    pred_lbls <- format_step06_predictor_labels(pred_raw)
    bar_cols <- ifelse(
      grepl("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)", pred_raw),
      "#009E73",
      ifelse(grepl("^MEM", pred_raw), "#CC79A7", "#0072B2")
    )
    graphics::barplot(
      imp_vals, names.arg = pred_lbls, horiz = TRUE, las = 1,
      col = bar_cols, border = "white", cex.names = 0.70,
      xlim = c(0, max(imp_vals, 0.05, na.rm = TRUE) * 1.18),
      xlab = "Permutation Importance (Mean Decrease Accuracy)",
      main = "A. Spatial RF Predictor Importance (10m-30m Stack)", font.main = 2
    )
    graphics::grid(nx = NULL, ny = NA, col = "grey88")
    graphics::legend(
      "bottomright",
      legend = c("10m-30m Edaphic/Lithology", "WorldClim Macro-Climate", "Spatial MEMs"),
      fill = c("#009E73", "#0072B2", "#CC79A7"), border = NA, cex = 0.68, bty = "n"
    )

    # Panel B: Dynamic Top-2 RF Environmental Predictors (Hollow Markers + 85% Ellipses)
    graphics::par(mar = c(4.5, 4.8, 3.0, 1.2))
    xv <- suppressWarnings(as.numeric(df[[top_env_1]]))
    yv <- suppressWarnings(as.numeric(df[[top_env_2]]))
    xv[!is.finite(xv)] <- stats::median(xv[is.finite(xv)], na.rm = TRUE)
    yv[!is.finite(yv)] <- stats::median(yv[is.finite(yv)], na.rm = TRUE)
    xv[!is.finite(xv)] <- 0.0; yv[!is.finite(yv)] <- 0.0

    graphics::plot(
      xv, yv, type = "n",
      xlab = top_env_lbls[1], ylab = top_env_lbls[2],
      main = "B. Top Spatial RF Niche Predictors (85% Confidence Ellipses)", font.main = 2
    )
    graphics::grid(col = "grey88")
    graphics::points(xv, yv, col = cols_tx, pch = pch_tx, bg = NA, cex = 0.85, lwd = 1.0)
    t_seq <- seq(0, 2 * pi, length.out = 100)
    chi_r <- sqrt(stats::qchisq(0.85, df = 2))
    for (tx in present_tx) {
      m_g <- which(tx_clean == tx & is.finite(xv) & is.finite(yv))
      if (length(m_g) >= 3L) {
        cov_g <- tryCatch(stats::cov(cbind(xv[m_g], yv[m_g])), error = function(e) NULL)
        if (!is.null(cov_g) && all(is.finite(cov_g)) && det(cov_g) > 1e-12) {
          eg <- eigen(cov_g, symmetric = TRUE)
          pts_el <- t(c(mean(xv[m_g]), mean(yv[m_g])) + chi_r * eg$vectors %*% diag(sqrt(pmax(eg$values, 0))) %*% rbind(cos(t_seq), sin(t_seq)))
          graphics::lines(pts_el[, 1], pts_el[, 2], col = dyn_pal[[tx]], lwd = 1.5)
        }
      }
    }
    short_tx <- gsub("^Packera\\s+", "P. ", present_tx)
    graphics::legend(
      "topright", legend = short_tx, col = dyn_pal[present_tx],
      pch = dyn_pch[present_tx], pt.bg = NA, cex = 0.72, bty = "n"
    )

    # Panel C: Top Edaphic/Lithic Predictor Boxplots + Kruskal-Wallis p-value
    ed_v <- suppressWarnings(as.numeric(df[[top_edaph_1]]))
    ed_v[!is.finite(ed_v)] <- stats::median(ed_v[is.finite(ed_v)], na.rm = TRUE)
    ed_v[!is.finite(ed_v)] <- 0.0
    kw_p <- tryCatch(stats::kruskal.test(ed_v ~ factor(tx_clean))$p.value, error = function(e) NA_real_)
    p_str <- if (is.finite(kw_p)) sprintf("Kruskal-Wallis p = %.2e", kw_p) else "Edaphic Divergence"
    tx_fac <- factor(tx_clean, levels = present_tx, labels = short_tx)
    graphics::boxplot(
      ed_v ~ tx_fac,
      col = grDevices::adjustcolor(dyn_pal[present_tx], alpha.f = 0.25),
      border = dyn_pal[present_tx], outpch = NA,
      xlab = "", ylab = top_edaph_lbl,
      main = sprintf("C. %s (%s)", top_edaph_lbl, p_str),
      font.axis = 3, cex.axis = 0.80, font.main = 2
    )
    graphics::grid(nx = NA, ny = NULL, col = "grey88")
    for (s_i in seq_along(present_tx)) {
      m_s <- which(tx_clean == present_tx[s_i])
      if (length(m_s) > 0L) {
        jit_x <- s_i + (((seq_along(m_s) %% 7L) - 3.0) * 0.045)
        graphics::points(
          jit_x, ed_v[m_s],
          pch = dyn_pch[[present_tx[s_i]]], bg = NA,
          col = grDevices::adjustcolor(dyn_pal[[present_tx[s_i]]], alpha.f = alpha_val),
          cex = 0.75, lwd = 0.9
        )
      }
    }

    # Panel D: Pairwise Schoener's D Niche Overlap Matrix (Benjamini-Hochberg FDR Adjusted)
    n_sp <- length(present_tx)
    d_mat <- matrix(1.0, nrow = n_sp, ncol = n_sp, dimnames = list(short_tx, short_tx))
    p_mat <- matrix(NA_real_, nrow = n_sp, ncol = n_sp, dimnames = list(short_tx, short_tx))
    niche_df <- if (!is.null(niche_res) && is.list(niche_res) && !is.data.frame(niche_res) && !is.null(niche_res$summary)) {
      niche_res$summary
    } else if (is.data.frame(niche_res)) {
      niche_res
    } else {
      NULL
    }
    if (!is.null(niche_df) && nrow(niche_df) > 0L) {
      p_col_use <- intersect(c("p_value_D", "Equivalency_p_value"), names(niche_df))[1]
      p_adj_vec <- if (!is.na(p_col_use)) stats::p.adjust(as.numeric(niche_df[[p_col_use]]), method = "BH") else rep(NA_real_, nrow(niche_df))
      t1_col <- intersect(c("Taxon1", "Taxon_1"), names(niche_df))[1]
      t2_col <- intersect(c("Taxon2", "Taxon_2"), names(niche_df))[1]
      d_col <- intersect(c("Schoener_D", "Schoeners_D", "Warren_D"), names(niche_df))[1]
      if (!is.na(t1_col) && !is.na(t2_col) && !is.na(d_col)) {
        for (r_i in seq_len(nrow(niche_df))) {
          t1 <- sub("^Packera\\s+", "P. ", as.character(niche_df[[t1_col]][r_i]))
          t2 <- sub("^Packera\\s+", "P. ", as.character(niche_df[[t2_col]][r_i]))
          dv <- suppressWarnings(as.numeric(niche_df[[d_col]][r_i]))
          pv <- p_adj_vec[r_i]
          if (t1 %in% short_tx && t2 %in% short_tx && is.finite(dv)) {
            d_mat[t1, t2] <- dv; d_mat[t2, t1] <- dv
            p_mat[t1, t2] <- pv; p_mat[t2, t1] <- pv
          }
        }
      }
    }
    graphics::plot(
      c(0.5, n_sp + 0.5), c(0.5, n_sp + 0.5), type = "n", axes = FALSE,
      xlab = "Active Group", ylab = "",
      main = "D. Pairwise Schoener's D Niche Overlap & Equivalency", font.main = 2
    )
    ramp_d <- grDevices::colorRampPalette(c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))(101)
    for (i in seq_len(n_sp)) {
      y_pos <- n_sp - i + 1L
      for (j in seq_len(n_sp)) {
        dv <- d_mat[i, j]
        c_idx <- max(1L, min(101L, round(dv * 100) + 1L))
        graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = if (i == j) "#EEEEEE" else ramp_d[c_idx], border = "white", lwd = 1.5)
        lbl_cell <- if (i == j) "1.00\n(Self)" else if (is.finite(p_mat[i, j])) sprintf("D = %.2f\n(q = %.3f)", dv, p_mat[i, j]) else sprintf("D = %.2f", dv)
        graphics::text(j, y_pos, labels = lbl_cell, cex = 0.74, font = if (i == j) 1 else 2, col = if (i != j && dv > 0.62) "white" else "#111111")
      }
    }
    graphics::axis(1, at = seq_len(n_sp), labels = short_tx, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::axis(2, at = rev(seq_len(n_sp)), labels = short_tx, las = 1, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::box(col = "grey75")
    u_sp_hdr <- get_active_packera_taxa(df$species_standardized)
    hdr_lbl <- if (length(u_sp_hdr) == 1L) u_sp_hdr[1] else sprintf("Packera Multi-Taxon Cohort (%d Active Groups)", length(present_tx))
    graphics::title(sprintf("%s: Multi-Scale Ecological Niche, Top RF Predictors & Equivalency", hdr_lbl), outer = TRUE, cex.main = 1.20)
  }

  pdf_out <- file.path(fdir, "fig06_ecological_niche_phenology.pdf")
  pdf_tmp <- paste0(pdf_out, ".tmp")
  grDevices::pdf(pdf_tmp, width = 13.0, height = 9.5)
  tryCatch(render_plate(), finally = grDevices::dev.off())
  if (file.exists(pdf_tmp)) file.rename(pdf_tmp, pdf_out)

  png_out <- file.path(fdir, "fig06_ecological_niche_phenology.png")
  png_tmp <- paste0(png_out, ".tmp")
  grDevices::png(png_tmp, width = 13.0, height = 9.5, units = "in", res = 300)
  tryCatch(render_plate(), finally = grDevices::dev.off())
  if (file.exists(png_tmp)) file.rename(png_tmp, png_out)
  invisible(list(pdf = pdf_out, png = png_out))
}

export_spatial_rf_figures <- function(df, srf_res, niche_res, out_pdf) {
  message(sprintf("Exporting publication-quality 5-panel figure to %s...", out_pdf))
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  
  # Theme setup
  thm <- ggplot2::theme_bw(base_size = 9) +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                        plot.subtitle = ggplot2::element_text(size = 8.2, color = "grey25"),
                        panel.grid.minor = ggplot2::element_blank())

  # Panel A: Spatial RF Variable Importance (Human-Readable Ecological Labels & Okabe-Ito Palette)
  imp_df <- srf_res$importance
  imp_df$Label <- if (exists("format_step06_predictor_labels", mode = "function")) {
    format_step06_predictor_labels(as.character(imp_df$Predictor))
  } else if (exists("format_predictor_labels", mode = "function")) {
    format_predictor_labels(as.character(imp_df$Predictor))
  } else {
    as.character(imp_df$Predictor)
  }
  p1 <- ggplot2::ggplot(imp_df, ggplot2::aes(x = stats::reorder(Label, MeanDecreaseAccuracy),
                                             y = MeanDecreaseAccuracy, fill = Type)) +
    ggplot2::geom_col(width = 0.72, color = "#222222", linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(
      "Fine-Grain Edaphic & Outcrop (10-30m)" = "#009E73",
      "Environmental" = "#0072B2",
      "Spatial (MEMs)" = "#CC79A7"
    )) +
    thm + ggplot2::labs(
      title = "A. Spatial RF Variable Importance (10m-30m Edaphic + Climate + MEMs)",
      subtitle = sprintf("OOB Accuracy = %.1f%% | Spatial Block CV Accuracy = %.1f%%",
                         100 * (1 - srf_res$OOB_Error), 100 * srf_res$Spatial_Block_CV_Accuracy),
      x = NULL, y = "Mean Decrease Accuracy (OOB)", fill = "Predictor Class"
    )

  # Panel B: Dynamic Top-2 Non-MEM Spatial RF Environmental Predictors (Hollow Markers + 85% Ellipses)
  df$plot_group <- if (exists("resolve_plot_group_column", mode = "function")) {
    resolve_plot_group_column(df)
  } else {
    sp_v <- as.character(df$species_standardized)
    if (length(unique(stats::na.omit(sp_v))) < 2L ||
        ("analysis_group" %in% names(df) && any(grepl("Morphotype_K|^LatTertile_|Ecotype|Subcluster", df$analysis_group, ignore.case = TRUE))) ||
        ("grouping_col" %in% names(df) && any(grepl("Morphotype_K|^LatTertile_|Ecotype|Subcluster", df$grouping_col, ignore.case = TRUE)))) {
      if ("analysis_group" %in% names(df) && length(unique(stats::na.omit(df$analysis_group))) >= 2L) as.character(df$analysis_group)
      else if ("grouping_col" %in% names(df) && length(unique(stats::na.omit(df$grouping_col))) >= 2L) as.character(df$grouping_col)
      else if ("gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster))) paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster)))
      else sp_v
    } else sp_v
  }
  active_tx <- sort(unique(stats::na.omit(df$plot_group)))
  active_tx <- active_tx[nzchar(active_tx) & active_tx != "Unknown"]
  if (length(active_tx) == 0) active_tx <- get_active_packera_taxa(df$species_standardized)
  is_subcluster_mode <- any(grepl("Morphotype_K|\\[K[0-9]+\\]|^LatTertile_|Ecotype|Subcluster", active_tx, ignore.case = TRUE)) ||
    length(unique(stats::na.omit(df$species_standardized))) < 2L
  grp_lbl <- if (is_subcluster_mode) "Ecotype / Morphotype" else "Taxon"
  dyn_pal <- resolve_dynamic_group_palette(active_tx)
  dyn_shapes <- stats::setNames(rep_len(c(21, 24, 22, 23, 1, 2, 5, 0, 6, 4), length(active_tx)), active_tx)
  plot_df <- df[df$plot_group %in% active_tx, , drop = FALSE]
  if (!("fine_soil_sand" %in% names(plot_df)) && "soil_sand" %in% names(plot_df)) plot_df$fine_soil_sand <- plot_df$soil_sand
  if (!("fine_soil_ph" %in% names(plot_df)) && "soil_ph" %in% names(plot_df)) plot_df$fine_soil_ph <- plot_df$soil_ph
  if (!("fine_bedrock_depth_cm" %in% names(plot_df))) plot_df$fine_bedrock_depth_cm <- 150.0
  if (!("ssurgo_rock_outcrop_pct" %in% names(plot_df))) plot_df$ssurgo_rock_outcrop_pct <- 0.0
  if (!("topo_tpi_30m" %in% names(plot_df))) plot_df$topo_tpi_30m <- 0.0

  top_sel <- select_top_rf_predictors(imp_df, plot_df)
  top_env_1 <- top_sel$top_env_1
  top_env_2 <- top_sel$top_env_2
  top_edaph_1 <- top_sel$top_edaph_1

  env_axis_lbls <- format_predictor_labels(c(top_env_1, top_env_2))
  if (exists("format_step06_predictor_labels", mode = "function")) {
    env_axis_lbls <- format_step06_predictor_labels(c(top_env_1, top_env_2))
  }
  edaph_lbl <- format_predictor_labels(top_edaph_1)[1]
  if (exists("format_step06_predictor_labels", mode = "function")) {
    edaph_lbl <- format_step06_predictor_labels(top_edaph_1)[1]
  }

  plot_df$top_env_1_val <- suppressWarnings(as.numeric(plot_df[[top_env_1]]))
  plot_df$top_env_2_val <- suppressWarnings(as.numeric(plot_df[[top_env_2]]))
  plot_df$top_edaph_1_val <- suppressWarnings(as.numeric(plot_df[[top_edaph_1]]))

  pt_alpha <- if (nrow(plot_df) > 1000) 0.45 else 0.78
  p2 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = top_env_1_val, y = top_env_2_val, color = plot_group, shape = plot_group)) +
    ggplot2::stat_ellipse(level = 0.85, linewidth = 1.0) +
    ggplot2::geom_point(fill = NA, stroke = 0.75, alpha = pt_alpha, size = 1.55) +
    ggplot2::scale_shape_manual(values = dyn_shapes) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::labs(
      title = "B. Top Spatial RF Ecological Niche Envelopes (85% Confidence Ellipses)",
      subtitle = sprintf("Top Non-MEM Predictors: %s vs. %s", env_axis_lbls[1], env_axis_lbls[2]),
      x = env_axis_lbls[1], y = env_axis_lbls[2], color = grp_lbl, shape = grp_lbl
    )

  # Panel C: Top-Importance Edaphic/Lithic Predictor Diagnostic Boxplots + Kruskal-Wallis p-value
  kw_edaph <- tryCatch(stats::kruskal.test(top_edaph_1_val ~ factor(plot_group), data = plot_df)$p.value, error = function(e) NA_real_)
  p3 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = plot_group, y = top_edaph_1_val, fill = plot_group, color = plot_group)) +
    ggplot2::geom_boxplot(alpha = 0.30, outlier.shape = 1, outlier.size = 1.2, width = 0.58) +
    ggplot2::scale_fill_manual(values = dyn_pal) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1, face = if (is_subcluster_mode) "bold" else "italic"), legend.position = "none") +
    ggplot2::labs(
      title = sprintf("C. %s", edaph_lbl),
      subtitle = if (!is.na(kw_edaph) && is.finite(kw_edaph)) sprintf("Kruskal-Wallis p = %.2e across active groups", kw_edaph) else sprintf("Top-ranked edaphic/lithic predictor (%s)", top_edaph_1),
      x = NULL, y = edaph_lbl
    )

  # Panel D: Warren's Niche Identity Permutations & Empirical Schoener's D
  if (!is.null(niche_res$permutations) && is.data.frame(niche_res$permutations) && nrow(niche_res$permutations) > 0L) {
    p4 <- ggplot2::ggplot(niche_res$permutations, ggplot2::aes(x = D_null, fill = Pair)) +
      ggplot2::geom_histogram(bins = 20, alpha = 0.65, position = "identity", color = "white", linewidth = 0.2) +
      ggplot2::geom_vline(data = niche_res$summary, ggplot2::aes(xintercept = Warren_D, color = Pair),
                          linetype = "dashed", linewidth = 0.95, show.legend = FALSE) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      thm + ggplot2::labs(
        title = "D. Pairwise Niche Equivalency (Schoener's D vs. Null Permutations)",
        subtitle = "Dashed vertical lines = Empirical Schoener's D; Histogram = Null niche identity distribution",
        x = "Schoener's D Niche Overlap (0 = Disjoint, 1 = Identical)", y = "Permutation Count", fill = "Taxon Pair"
      )
  } else {
    p4 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pairwise niche equivalency requires >= 2 groups with N >= 5") +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      thm + ggplot2::labs(
        title = "D. Pairwise Niche Equivalency (Schoener's D vs. Null Permutations)",
        x = "Schoener's D Niche Overlap", y = "Permutation Count"
      )
  }

  # Panel E: Horizontal Stacked Variance-Partitioning Bar Chart (Replaces Pie Chart)
  vp_df <- srf_res$variance_partition
  vp_df$Label <- ifelse(vp_df$Variance_Pct >= 5, sprintf("%.1f%%", vp_df$Variance_Pct), "")
  vp_df$Component <- factor(vp_df$Component, levels = rev(vp_df$Component))
  vp_cols <- stats::setNames(c("#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#999999")[seq_len(nrow(vp_df))], levels(vp_df$Component))
  p5 <- ggplot2::ggplot(vp_df, ggplot2::aes(x = "Variance\nPartition", y = Variance_Pct, fill = Component)) +
    ggplot2::geom_col(width = 0.52, color = "white", linewidth = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(ggplot2::aes(label = Label), position = ggplot2::position_stack(vjust = 0.5),
                       size = 3.1, fontface = "bold", color = "white") +
    ggplot2::scale_fill_manual(values = vp_cols) +
    thm + ggplot2::theme(legend.position = "bottom", legend.direction = "vertical") +
    ggplot2::labs(
      title = "E. Multi-Scale Variance Partitioning",
      subtitle = "Relative contribution of [a] Edaphic/Outcrop, [b] Climate & [c] Spatial MEMs",
      x = NULL, y = "Percentage of Explained Niche Variance (%)", fill = "Variance Component"
    )

  tmp_pdf <- paste0(out_pdf, ".tmp")
  pdf(tmp_pdf, width = 13, height = 8.5)
  tryCatch({
    gridExtra::grid.arrange(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 2, 5), c(3, 4, 5)), widths = c(1.15, 1.15, 0.95))
  }, finally = {
    dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, out_pdf)

  # Export auxiliary Step 06 figures via generate_step06_pub_figures if available,
  # then atomically render dynamic top-RF-importance fig06_ecological_niche_phenology.{pdf,png}
  if (exists("generate_step06_pub_figures", mode = "function")) {
    tryCatch(
      generate_step06_pub_figures(df, srf_res, niche_res, NULL, out_pdf),
      error = function(e) warning("generate_step06_pub_figures warning: ", e$message)
    )
  }
  tryCatch(
    export_dynamic_fig06_ecological_niche_phenology(df, srf_res, niche_res, out_pdf),
    error = function(e) warning("export_dynamic_fig06_ecological_niche_phenology warning: ", e$message)
  )
}

# ------------------------------------------------------------------------------
# 7. Main Execution Pipeline
# ------------------------------------------------------------------------------
main <- function() {
  set.seed(42)
  opts <- parse_args_robust()
  message("=== Starting Multimodal Spatial RF & Ecological Niche Analysis ===")
  
  if (!file.exists(opts$vouchers)) stop("Missing vouchers file: ", opts$vouchers)
  vouchers_df <- readr::read_csv(opts$vouchers, show_col_types = FALSE)
  vouchers_df$catalogNumber <- as.character(vouchers_df$catalogNumber)
  sp_cand <- rep("Unknown", nrow(vouchers_df))
  for (sc in intersect(c("species_standardized", "scientificName", "species_raw"), names(vouchers_df))) {
    v <- standardize_packera_taxon(vouchers_df[[sc]])
    m <- sp_cand == "Unknown" & !is.na(v) & v != "Unknown"
    sp_cand[m] <- v[m]
  }
  vouchers_df$species_standardized <- sp_cand

  # Merge Morphometrics if present (resolve morphometric_flags.csv <-> morphometrics_misidentification_flags.csv)
  if (!file.exists(opts$morphometrics)) {
    sibling_morph <- file.path(dirname(opts$vouchers), "morphometric_flags.csv")
    if (file.exists(sibling_morph)) {
      opts$morphometrics <- sibling_morph
    } else if (file.exists("data/tables/morphometric_flags.csv")) {
      opts$morphometrics <- "data/tables/morphometric_flags.csv"
    }
  }
  if (file.exists(opts$morphometrics)) {
    morph_df <- readr::read_csv(opts$morphometrics, show_col_types = FALSE)
    morph_df$catalogNumber <- as.character(morph_df$catalogNumber)
    if ("morphometric_outlier_flag" %in% names(morph_df) && !("misidentification_flag" %in% names(morph_df))) {
      morph_df$misidentification_flag <- morph_df$morphometric_outlier_flag
    } else if ("misidentification_flag" %in% names(morph_df) && !("morphometric_outlier_flag" %in% names(morph_df))) {
      morph_df$morphometric_outlier_flag <- morph_df$misidentification_flag
    }
    m_cols <- intersect(names(morph_df), c("catalogNumber", "plant_individual_id", "morph_tier_source", "n_tier1_leaves", "n_tier2_leaves", "cda_predicted_taxon", "cda_posterior_prob", "can1", "can2", "can3", "gmm_cluster", "morphometric_outlier_flag", "misidentification_flag"))
    morph_sub <- morph_df[!duplicated(morph_df$catalogNumber), m_cols, drop = FALSE]
    vouchers_df <- dplyr::left_join(vouchers_df, morph_sub, by = "catalogNumber")
  }

  # Merge Cleanlab Vision Audit if present
  if (file.exists(opts$vision_audit)) {
    vis_df <- readr::read_csv(opts$vision_audit, show_col_types = FALSE)
    vis_df$catalogNumber <- as.character(vis_df$catalogNumber)
    v_cols <- intersect(names(vis_df), c("catalogNumber", "predicted_label", "confidence_predicted_class", "c_error", "is_label_corrupted"))
    names(vis_df)[names(vis_df) == "predicted_label"] <- "vision_predicted_label"
    vis_sub <- vis_df[!duplicated(vis_df$catalogNumber), intersect(names(vis_df), c("catalogNumber", "vision_predicted_label", "c_error", "is_label_corrupted")), drop = FALSE]
    vouchers_df <- dplyr::left_join(vouchers_df, vis_sub, by = "catalogNumber")
  }

  # Ingest Voucher Phenological States (Computer Vision Anthesis Ingestion via left_join)
  if (file.exists(opts$pheno_states)) {
    message(sprintf("Ingesting voucher phenological states from %s...", opts$pheno_states))
    pheno_states_df <- readr::read_csv(opts$pheno_states, show_col_types = FALSE)
    pheno_states_df$catalogNumber <- as.character(pheno_states_df$catalogNumber)
    vouchers_df$catalogNumber <- as.character(vouchers_df$catalogNumber)
    p_cols <- intersect(names(pheno_states_df), c("catalogNumber", "total_capitula", "n_anthesis", "n_fruit", "n_bud", "voucher_phenological_state", "dominant_pappus_ratio"))
    pheno_states_df <- pheno_states_df[!duplicated(pheno_states_df$catalogNumber), p_cols, drop = FALSE]
    vouchers_df <- dplyr::left_join(vouchers_df, pheno_states_df, by = "catalogNumber")
    message(sprintf("Left-joined %d vouchers with phenological states.", nrow(vouchers_df)))
  } else {
    message(sprintf("[INFO] Phenological states CSV not found at %s. Gracefully falling back to calendar + Darwin Core text.", opts$pheno_states))
    vouchers_df$voucher_phenological_state <- NA_character_
    vouchers_df$total_capitula <- NA_real_
    vouchers_df$n_anthesis <- NA_real_
    vouchers_df$n_fruit <- NA_real_
    vouchers_df$n_bud <- NA_real_
    vouchers_df$dominant_pappus_ratio <- NA_real_
  }

  # Filter strictly for primary duplicate vouchers (exsiccatae deduplication)
  if ("is_primary_duplicate" %in% names(vouchers_df)) {
    n_pre_filter <- nrow(vouchers_df)
    vouchers_df <- vouchers_df[toupper(as.character(vouchers_df$is_primary_duplicate)) %in% c("TRUE", "T", "1"), ]
    n_post_filter <- nrow(vouchers_df)
    message(sprintf("Exsiccatae Stratification Filter: Retaining %d primary vouchers (%d duplicate sheets excluded to prevent spatial pseudoreplication).",
                    n_post_filter, n_pre_filter - n_post_filter))
  }

  # Filter valid coordinates
  vouchers_df <- vouchers_df[!is.na(vouchers_df$latitude) & !is.na(vouchers_df$longitude), ]
  message(sprintf("Processing %d georeferenced vouchers across target taxa...", nrow(vouchers_df)))

  demo_mode <- isTRUE(opts$demo)

  # Stage 1: Environmental Layer Extraction (250m SoilGrids & WorldClim v2.1)
  vouchers_df <- extract_environmental_layers(vouchers_df, opts$env_dir, demo = demo_mode)

  # Stage 1b: 10m/50m-Buffered USDA SSURGO Vector Pedology, 30m POLARIS/gSSURGO, 100m SOLUS100, SGMC Lithology & 30m 3DEP Micro-Topography
  if (exists("load_and_update_fine_grain_cache", mode = "function")) {
    fine_cache_csv <- file.path(opts$env_dir, "fine_grain_edaphic_cache.csv")
    vouchers_df <- load_and_update_fine_grain_cache(
      vouchers_df,
      env_dir = opts$env_dir,
      ssurgo_cache_dir = opts$ssurgo_cache,
      batch_size = opts$batch_size,
      demo = demo_mode,
      query_ssurgo_fn = query_ssurgo_sda,
      cache_path = fine_cache_csv
    )
  } else {
    vouchers_df <- query_ssurgo_sda(vouchers_df, cache_dir = opts$ssurgo_cache, batch_size = opts$batch_size, demo = demo_mode)
    if (exists("extract_sgmc_lithology", mode = "function")) {
      vouchers_df <- extract_sgmc_lithology(vouchers_df, env_dir = opts$env_dir, demo = demo_mode)
    }
    if (exists("extract_uncertainty_buffered_predictors", mode = "function")) {
      vouchers_df <- extract_uncertainty_buffered_predictors(vouchers_df, env_dir = opts$env_dir, demo = demo_mode)
    }
  }

  edaphic_synth <- synthesize_multiscale_edaphics(
    vouchers_df,
    ssurgo_out = opts$ssurgo_table,
    contingency_out = opts$contingency_table,
    demo = demo_mode
  )

  # Stage 1c: Latitudinal Spring Baseline & Phenological Anomalies Engine
  pheno_res <- compute_phenological_anomalies(
    vouchers_df,
    output_csv = opts$pheno_anomalies,
    summary_csv = opts$pheno_summary,
    demo = demo_mode
  )
  vouchers_df <- pheno_res$df
  export_phenological_latitudinal_figures(vouchers_df, pheno_res$baseline, opts$pheno_plot)

  # Impute median delta_doy per species if NA for complete Spatial RF matrices
  for (tx in unique(vouchers_df$species_standardized)) {
    tx_mask <- vouchers_df$species_standardized == tx
    med_d <- stats::median(vouchers_df$delta_doy[tx_mask], na.rm = TRUE)
    if (is.na(med_d)) med_d <- 0.0
    na_sub <- tx_mask & is.na(vouchers_df$delta_doy)
    if (any(na_sub)) vouchers_df$delta_doy[na_sub] <- med_d
  }

  env_vars <- c(
    "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm",
    "bdricm", "resdept_r", "cfvo", "awc", "fine_soil_cec", "alr_sand_clay", "alr_silt_clay",
    "is_azonal_lithology", "azonal_contact_dist_m", "azonal_proximity_index",
    "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer",
    "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm", "topo_tpi_30m", "topo_twi_30m",
    "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality",
    "delta_doy"
  )
  env_vars <- intersect(env_vars, names(vouchers_df))
  for (ev in env_vars) {
    if (ev %in% names(vouchers_df)) {
      vouchers_df[[ev]] <- suppressWarnings(as.numeric(vouchers_df[[ev]]))
      for (tx in unique(vouchers_df$species_standardized)) {
        tx_m <- vouchers_df$species_standardized == tx
        med_v <- stats::median(vouchers_df[[ev]][tx_m], na.rm = TRUE)
        if (is.na(med_v)) med_v <- stats::median(vouchers_df[[ev]], na.rm = TRUE)
        if (is.na(med_v)) med_v <- 0.0
        na_m <- tx_m & is.na(vouchers_df[[ev]])
        if (any(na_m)) vouchers_df[[ev]][na_m] <- med_v
      }
    }
  }

  # Stage 2: Resolve single-taxon / intraspecific analysis_group & Cross-Modal Consensus Checks
  merged_df <- resolve_analysis_group(vouchers_df, intraspecific_mode = isTRUE(opts$intraspecific_mode))
  vouchers_df <- execute_crossmodal_consensus(merged_df, intraspecific_mode = isTRUE(opts$intraspecific_mode))

  # Stage 3: Moran's Eigenvector Maps (MEM1..MEM5) & Authority-Anchored Spatial RF
  srf_res <- compute_spatial_rf_mems(vouchers_df, env_vars, n_trees = opts$n_trees, intraspecific_mode = isTRUE(opts$intraspecific_mode))
  vouchers_df$spatial_block_cv_acc <- round(srf_res$Spatial_Block_CV_Accuracy, 4)
  vouchers_df$spatial_autocorr_range_km <- srf_res$spatial_autocorr_range_km
  vouchers_df$rf_predicted_taxon <- srf_res$rf_predicted_taxon
  vouchers_df$rf_predicted_prob <- srf_res$rf_predicted_prob
  vouchers_df$rf_given_taxon_prob <- srf_res$rf_given_taxon_prob
  vouchers_df$rf_prediction_source <- srf_res$rf_prediction_source
  vouchers_df$spatial_block_id <- srf_res$spatial_block_id

  # Stage 4: Partitioned PCA-env Niche Equivalency & Background Similarity Tests
  niche_res <- run_warren_niche_identity_tests(vouchers_df, env_vars, n_perm = opts$permutations, intraspecific_mode = isTRUE(opts$intraspecific_mode))

  # Stage 5: Export Conflict Flags, Summary Report & Figures (Atomic write-then-replace)
  if (write_csv_atomic(vouchers_df, opts$output_flags, demo = demo_mode)) {
    message(sprintf("Saved %d multimodal conflict flags to %s.", nrow(vouchers_df), opts$output_flags))
  }

  if (write_csv_atomic(niche_res$summary, opts$output_summary, demo = demo_mode)) {
    message(sprintf("Saved Warren's Niche Identity summary to %s.", opts$output_summary))
  }

  export_spatial_rf_figures(vouchers_df, srf_res, niche_res, opts$output_plot)
  tryCatch({
    if (file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R")
      format_predictor_labels <<- format_step06_predictor_labels
      generate_step06_pub_figures(vouchers_df, srf_res, niche_res, pheno_res, opts$output_plot)
    }
    export_dynamic_fig06_ecological_niche_phenology(vouchers_df, srf_res, niche_res, opts$output_plot)
  }, error = function(e) warning("Step 06 publication figures export warning: ", e$message))
  message("=== Multimodal Spatial RF & Niche Modeling Analysis Complete ===")
}

if (!interactive() && sys.nframe() == 0L) main()
