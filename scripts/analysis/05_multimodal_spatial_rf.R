#!/usr/bin/env Rscript
# ==============================================================================
# Script: 05_multimodal_spatial_rf.R
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

for (.pf_p in c("scripts/morphometrics/pub_figures_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "morphometrics", "pub_figures_utils.R"))) {
  if (!exists("atomic_save_dual_plot", mode = "function") && !is.na(.pf_p) && nzchar(.pf_p) && file.exists(.pf_p)) { source(.pf_p); break }
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
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
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

WC_BIO_LOOKUP <- c(
  "1"  = "bio1_temp_mean",
  "2"  = "bio2_mean_diurnal_range",
  "3"  = "bio3_isothermality",
  "4"  = "bio4_temp_seasonality",
  "5"  = "bio5_max_temp_warmest_month",
  "6"  = "bio6_min_temp_coldest_month",
  "7"  = "bio7_temp_annual_range",
  "8"  = "bio8_mean_temp_wettest_qtr",
  "9"  = "bio9_mean_temp_driest_qtr",
  "10" = "bio10_mean_temp_warmest_qtr",
  "11" = "bio11_mean_temp_coldest_qtr",
  "12" = "bio12_precip_annual",
  "13" = "bio13_precip_wettest_month",
  "14" = "bio14_precip_driest_month",
  "15" = "bio15_precip_seasonality",
  "16" = "bio16_precip_wettest_qtr",
  "17" = "bio17_precip_driest_qtr",
  "18" = "bio18_precip_warmest_qtr",
  "19" = "bio19_precip_coldest_qtr"
)

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
  "bio2_mean_diurnal_range"          = "Mean Diurnal Temp Range (BIO2, deg C)",
  "bio3_isothermality"               = "Isothermality (BIO3, BIO2/BIO7 x 100)",
  "bio4_temp_seasonality"            = "Temp Seasonality (BIO4, SD x 100)",
  "bio5_max_temp_warmest_month"      = "Max Temp Warmest Month (BIO5, deg C)",
  "bio6_min_temp_coldest_month"      = "Min Temp Coldest Month (BIO6, deg C)",
  "bio7_temp_annual_range"           = "Temp Annual Range (BIO7, deg C)",
  "bio8_mean_temp_wettest_qtr"       = "Mean Temp Wettest Quarter (BIO8, deg C)",
  "bio9_mean_temp_driest_qtr"        = "Mean Temp Driest Quarter (BIO9, deg C)",
  "bio10_mean_temp_warmest_qtr"      = "Mean Temp Warmest Quarter (BIO10, deg C)",
  "bio11_mean_temp_coldest_qtr"      = "Mean Temp Coldest Quarter (BIO11, deg C)",
  "bio12_precip_annual"              = "Annual Precipitation (BIO12, mm)",
  "bio13_precip_wettest_month"       = "Wettest-Month Precip (BIO13, mm)",
  "bio14_precip_driest_month"        = "Driest-Month Precip (BIO14, mm)",
  "bio15_precip_seasonality"         = "Precipitation Seasonality (BIO15, CV)",
  "bio16_precip_wettest_qtr"         = "Wettest-Quarter Precip (BIO16, mm)",
  "bio17_precip_driest_qtr"          = "Driest-Quarter Precip (BIO17, mm)",
  "bio18_precip_warmest_qtr"         = "Warmest-Quarter Precip (BIO18, mm)",
  "bio19_precip_coldest_qtr"         = "Coldest-Quarter Precip (BIO19, mm)",
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
        compute_empirical_bayes_authority_weights(
          df,
          mode = "continuous_prior_posterior",
          target_modality = "ecological"
        ),
        error = function(e1) {
          compute_empirical_bayes_authority_weights(
            df = df,
            taxon_col = "species_standardized",
            tier_col = "determiner_tier",
            mode = "continuous_prior_posterior",
            target_modality = "ecological"
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
  has_eco <- ("rf_predicted_taxon" %in% names(df)) || ("edaphic_best_fit_taxon" %in% names(df))
  obs_succ <- rep(0.0, n)
  obs_trials <- rep(0.0, n)
  if (has_cda) {
    cda_v <- as.character(df$cda_predicted_taxon)
    ok_c <- !is.na(cda_v) & nzchar(cda_v) & sp_vec != "Unknown"
    obs_succ[ok_c] <- obs_succ[ok_c] + as.numeric(cda_v[ok_c] == sp_vec[ok_c])
    obs_trials[ok_c] <- obs_trials[ok_c] + 1.0
  }
  if (has_eco) {
    eco_col <- if ("rf_predicted_taxon" %in% names(df)) "rf_predicted_taxon" else "edaphic_best_fit_taxon"
    eco_v <- as.character(df[[eco_col]])
    ok_e <- !is.na(eco_v) & nzchar(eco_v) & sp_vec != "Unknown"
    obs_succ[ok_e] <- obs_succ[ok_e] + as.numeric(eco_v[ok_e] == sp_vec[ok_e])
    obs_trials[ok_e] <- obs_trials[ok_e] + 1.0
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
  if (exists(".identify_type_specimen_vec", mode = "function")) {
    ti <- .identify_type_specimen_vec(df)
    w_out[ti$type_rank == 2L] <- 3.0
    w_out[ti$type_rank == 1L] <- pmax(w_out[ti$type_rank == 1L], 2.8)
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
    optparse::make_option(c("--ssurgo-batch-size"), dest = "batch_size", type = "integer",
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
    raw_args <- commandArgs(trailingOnly = TRUE)
    raw_args <- sub("^--batch-size(?==|$)", "--ssurgo-batch-size", raw_args, perl = TRUE)
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list), args = raw_args)
    names(opts) <- gsub("-", "_", names(opts))
    return(opts)
  }
  return(list(
    vouchers = "data/tables/curated_vouchers.csv", morphometrics = "data/tables/morphometrics_misidentification_flags.csv",
    pheno_states = "data/tables/voucher_phenological_states.csv",
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
      ep <- sub("^[×-]+", "", ep)
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
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(df, file = tmp_path, na = "", quote = "auto", logical01 = FALSE)
  } else {
    readr::write_csv(df, tmp_path)
  }
  file.rename(tmp_path, target_path)
  return(invisible(TRUE))
}

# ------------------------------------------------------------------------------
# 2. Environmental Extraction (SoilGrids 250m 0-15cm Rhizosphere & Uncertainty-Buffered WorldClim v2.1)
# ------------------------------------------------------------------------------
extract_rast_buffered_mean <- function(rast_obj, pts_vect, radii_m) {
  r_num <- suppressWarnings(as.numeric(radii_m))
  r_num[is.na(r_num) | !is.finite(r_num)] <- 250.0
  r_num <- pmax(100.0, pmin(1000.0, r_num))

  crds <- tryCatch(terra::crds(pts_vect), error = function(e) NULL)
  inv_idx <- NULL
  pts_use <- pts_vect
  r_use <- r_num
  if (!is.null(crds) && is.matrix(crds) && nrow(crds) == length(r_num) && nrow(crds) > 1L) {
    u_key <- paste(round(crds[, 1], 6), round(crds[, 2], 6), round(r_num, 1), sep = "_")
    u_keys <- unique(u_key)
    if (length(u_keys) < nrow(crds)) {
      first_idx <- match(u_keys, u_key)
      inv_idx <- match(u_key, u_keys)
      pts_use <- pts_vect[first_idx]
      r_use <- r_num[first_idx]
    }
  }

  pt_vals <- terra::extract(rast_obj, pts_use)
  buf_vals <- tryCatch({
    pts_buf <- terra::buffer(pts_use, width = r_use)
    terra::extract(rast_obj, pts_buf, fun = mean, na.rm = TRUE)
  }, error = function(e) NULL)
  if (!is.null(buf_vals)) {
    for (col_nm in setdiff(names(pt_vals), "ID")) {
      if (col_nm %in% names(buf_vals)) {
        bv <- suppressWarnings(as.numeric(buf_vals[[col_nm]]))
        pv <- suppressWarnings(as.numeric(pt_vals[[col_nm]]))
        pt_vals[[col_nm]] <- ifelse(is.finite(bv), bv, pv)
      }
    }
  }
  if (!is.null(inv_idx)) {
    pt_vals <- pt_vals[inv_idx, , drop = FALSE]
    pt_vals$ID <- seq_along(inv_idx)
    rownames(pt_vals) <- NULL
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

coalesce_legacy_wc_bio_cols <- function(df, ensure_all_19 = TRUE) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  n_r <- nrow(df)
  if ("bio18_precip_warmest_qt" %in% names(df)) {
    v_qt <- suppressWarnings(as.numeric(df$bio18_precip_warmest_qt))
    v_qtr <- if ("bio18_precip_warmest_qtr" %in% names(df)) suppressWarnings(as.numeric(df$bio18_precip_warmest_qtr)) else rep(NA_real_, n_r)
    df$bio18_precip_warmest_qtr <- ifelse(is.finite(v_qtr), v_qtr, v_qt)
  }
  if ("bio19_precip_coldest_qt" %in% names(df)) {
    v_qt <- suppressWarnings(as.numeric(df$bio19_precip_coldest_qt))
    v_qtr <- if ("bio19_precip_coldest_qtr" %in% names(df)) suppressWarnings(as.numeric(df$bio19_precip_coldest_qtr)) else rep(NA_real_, n_r)
    df$bio19_precip_coldest_qtr <- ifelse(is.finite(v_qtr), v_qtr, v_qt)
  }
  if (isTRUE(ensure_all_19)) {
    for (b_col in unname(WC_BIO_LOOKUP)) {
      if (!(b_col %in% names(df))) {
        df[[b_col]] <- rep(NA_real_, n_r)
      } else {
        df[[b_col]] <- suppressWarnings(as.numeric(df[[b_col]]))
      }
    }
  }
  df
}

assign_wc_bio_from_extracted <- function(df, wc_vals) {
  if (is.null(wc_vals) || !is.data.frame(wc_vals)) return(coalesce_legacy_wc_bio_cols(df))
  val_nms <- setdiff(names(wc_vals), "ID")
  for (b_idx in seq_len(19L)) {
    target_col <- WC_BIO_LOOKUP[[as.character(b_idx)]]
    cand_hit <- character(0)
    if (target_col %in% val_nms) {
      cand_hit <- target_col
    } else if (b_idx == 18L && "bio18_precip_warmest_qt" %in% val_nms) {
      cand_hit <- "bio18_precip_warmest_qt"
    } else if (b_idx == 19L && "bio19_precip_coldest_qt" %in% val_nms) {
      cand_hit <- "bio19_precip_coldest_qt"
    } else {
      pat_exact <- sprintf("(?:^|[._-])bio[_.-]?0?%d$", b_idx)
      m_ex <- grep(pat_exact, val_nms, value = TRUE, ignore.case = TRUE, perl = TRUE)
      if (length(m_ex) > 0L) {
        cand_hit <- m_ex[1]
      } else {
        pat_pref <- sprintf("^bio0?%d(?:_|$)", b_idx)
        m_pr <- grep(pat_pref, val_nms, value = TRUE, ignore.case = TRUE, perl = TRUE)
        if (length(m_pr) > 0L) cand_hit <- m_pr[1]
      }
    }
    if (length(cand_hit) == 0L && length(val_nms) == 19L && all(grepl("bio|lyr|band", val_nms, ignore.case = TRUE))) {
      cand_hit <- val_nms[b_idx]
    }
    if (length(cand_hit) > 0L) {
      v_ext <- suppressWarnings(as.numeric(wc_vals[[cand_hit[1]]]))
      if (target_col %in% names(df)) {
        v_cur <- suppressWarnings(as.numeric(df[[target_col]]))
        df[[target_col]] <- ifelse(is.finite(v_ext), v_ext, v_cur)
      } else {
        df[[target_col]] <- v_ext
      }
    }
  }
  coalesce_legacy_wc_bio_cols(df)
}

discover_candidate_eco_predictors <- function(df, include_pheno = TRUE, min_finite = 3L) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(character(0))
  df <- coalesce_legacy_wc_bio_cols(df, ensure_all_19 = FALSE)
  nms <- names(df)
  wc_cols <- intersect(unname(WC_BIO_LOOKUP), nms)
  pheno_cols <- if (isTRUE(include_pheno)) intersect("delta_doy", nms) else character(0)

  fine_pattern <- paste0(
    "^(fine_|ssurgo_|polaris_|solus_|bdricm$|resdept_r$|cfvo$|awc$|",
    "alr_|clr_|is_azonal_lithology$|azonal_|",
    "ph_q10_buffer$|ph_q90_buffer$|min_bedrock_depth_buffer_cm$|max_cfvo_buffer_pct$|",
    "azonal_lithology_presence_buffer$)"
  )
  fine_cols <- grep(fine_pattern, nms, value = TRUE, ignore.case = TRUE)

  has_finite_col <- function(col_nm) {
    if (!(col_nm %in% nms)) return(FALSE)
    v <- suppressWarnings(as.numeric(df[[col_nm]]))
    sum(is.finite(v)) >= min(min_finite, max(2L, nrow(df)))
  }

  sg_fallback <- character(0)
  if (!has_finite_col("fine_soil_ph") && "soil_ph" %in% nms) sg_fallback <- c(sg_fallback, "soil_ph")
  if (!has_finite_col("fine_soil_sand") && "soil_sand" %in% nms) sg_fallback <- c(sg_fallback, "soil_sand")
  if (!has_finite_col("fine_soil_cec") && "soil_cec" %in% nms) sg_fallback <- c(sg_fallback, "soil_cec")
  if (!has_finite_col("polaris_bd_30m") && "soil_bulk_density" %in% nms && length(fine_cols) == 0L) {
    sg_fallback <- c(sg_fallback, "soil_bulk_density")
  }

  topo_cols <- grep("^topo_", nms, value = TRUE, ignore.case = TRUE)

  exclude_cols <- c(
    "catalogNumber", "plant_individual_id", "latitude", "longitude",
    "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
    "coordinateUncertaintyInMeters", "buffer_radius_m",
    "high_coord_uncertainty_flag", "is_training_eligible_coord",
    "edaphic_best_fit_taxon", "edaphic_z_score", "climatic_z_score",
    "soil_silt_pct", "bio18_precip_warmest_qt", "bio19_precip_coldest_qt"
  )
  cand <- setdiff(unique(c(fine_cols, sg_fallback, topo_cols, wc_cols, pheno_cols)), exclude_cols)

  min_obs <- min(as.integer(min_finite), max(2L, nrow(df)))
  keep <- vapply(cand, function(cl) {
    raw_col <- df[[cl]]
    if (!(is.numeric(raw_col) || is.integer(raw_col) || is.logical(raw_col))) return(FALSE)
    v <- suppressWarnings(as.numeric(raw_col))
    fin <- is.finite(v)
    if (sum(fin) < min_obs) return(FALSE)
    sd_v <- stats::sd(v[fin])
    is.finite(sd_v) && sd_v > 1e-8
  }, logical(1))
  cand[keep]
}

crop_and_compress_na_rasters <- function(rast_obj = NULL, env_dir = NULL, na_bbox = c(-130.0, -50.0, 14.0, 66.0)) {
  if (!requireNamespace("terra", quietly = TRUE)) return(rast_obj)
  na_ext <- terra::ext(na_bbox[1], na_bbox[2], na_bbox[3], na_bbox[4])
  gdal_opts <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

  .crop_single_tif_on_disk <- function(tif_path) {
    if (is.null(tif_path) || !nzchar(tif_path) || !file.exists(tif_path)) return(tif_path)
    if (grepl("\\.crop_tmp\\.tif$", tif_path, ignore.case = TRUE)) return(tif_path)
    r_disk <- tryCatch(terra::rast(tif_path), error = function(e) NULL)
    if (is.null(r_disk)) return(tif_path)
    ex <- tryCatch(terra::ext(r_disk), error = function(e) NULL)
    if (is.null(ex)) return(tif_path)
    is_global_uncropped <- (terra::xmin(ex) < -150.0) || (terra::xmax(ex) > 0.0) ||
                           (terra::ymin(ex) < 0.0) || (terra::ymax(ex) > 75.0)
    if (!isTRUE(is_global_uncropped)) return(tif_path)

    tmp_tif <- paste0(tif_path, ".crop_tmp.tif")
    ok_crop <- tryCatch({
      r_cropped <- terra::crop(r_disk, na_ext)
      tryCatch(
        terra::writeRaster(r_cropped, filename = tmp_tif, overwrite = TRUE, gdal = gdal_opts),
        error = function(e_pred) {
          terra::writeRaster(r_cropped, filename = tmp_tif, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"))
        }
      )
      TRUE
    }, error = function(e) FALSE)
    rm(r_disk)
    if (isTRUE(ok_crop) && file.exists(tmp_tif)) {
      unlink(tif_path)
      file.rename(tmp_tif, tif_path)
      message(sprintf("[Raster Cropping & DEFLATE] Cropped global GeoTIFF to North America [-130, -50, 14, 66] & compressed: %s", basename(tif_path)))
    } else if (file.exists(tmp_tif)) {
      unlink(tmp_tif)
    }
    tif_path
  }

  cand_files <- character(0)
  if (!is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir)) {
    cand_files <- c(
      cand_files,
      list.files(env_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    )
  }
  if (is.character(rast_obj)) {
    cand_files <- c(cand_files, rast_obj[file.exists(rast_obj)])
  } else if (inherits(rast_obj, "SpatRaster")) {
    srcs <- tryCatch(terra::sources(rast_obj), error = function(e) character(0))
    srcs <- srcs[nzchar(srcs) & file.exists(srcs)]
    cand_files <- c(cand_files, srcs)
  } else if (is.list(rast_obj)) {
    for (el in rast_obj) {
      if (inherits(el, "SpatRaster")) {
        srcs <- tryCatch(terra::sources(el), error = function(e) character(0))
        srcs <- srcs[nzchar(srcs) & file.exists(srcs)]
        cand_files <- c(cand_files, srcs)
      }
    }
  }

  cand_files <- unique(cand_files[nzchar(cand_files) & file.exists(cand_files)])
  if (length(cand_files) > 0L) {
    for (cf in cand_files) {
      tryCatch(.crop_single_tif_on_disk(cf), error = function(e) NULL)
    }
  }

  if (inherits(rast_obj, "SpatRaster")) {
    srcs <- tryCatch(terra::sources(rast_obj), error = function(e) character(0))
    srcs <- srcs[nzchar(srcs) & file.exists(srcs)]
    if (length(srcs) > 0L && all(file.exists(srcs))) {
      rast_obj <- tryCatch(terra::rast(srcs), error = function(e) terra::crop(rast_obj, na_ext))
    }
    ex_curr <- tryCatch(terra::ext(rast_obj), error = function(e) NULL)
    if (!is.null(ex_curr) && (terra::xmin(ex_curr) < -135.0 || terra::xmax(ex_curr) > -45.0 ||
                              terra::ymin(ex_curr) < 10.0 || terra::ymax(ex_curr) > 70.0)) {
      rast_obj <- tryCatch(terra::crop(rast_obj, na_ext), error = function(e) rast_obj)
    }
    return(rast_obj)
  } else if (is.list(rast_obj)) {
    return(lapply(rast_obj, function(el) {
      if (inherits(el, "SpatRaster")) crop_and_compress_na_rasters(el, env_dir = NULL, na_bbox = na_bbox) else el
    }))
  }
  rast_obj
}

discover_env_rasters <- function(env_dir, recursive = TRUE, crop_global = TRUE) {
  if (is.null(env_dir) || !nzchar(env_dir) || !dir.exists(env_dir)) return(character(0))
  if (isTRUE(crop_global) && requireNamespace("terra", quietly = TRUE)) {
    tryCatch(crop_and_compress_na_rasters(env_dir = env_dir), error = function(e) NULL)
  }
  f_vec <- list.files(env_dir, pattern = "\\.(tif|tiff|grd)$", full.names = TRUE, recursive = isTRUE(recursive), ignore.case = TRUE)
  f_vec[!grepl("\\.crop_tmp\\.tif$", f_vec, ignore.case = TRUE)]
}

extract_environmental_layers <- function(df, env_dir, demo = FALSE) {
  message("Integrating SoilGrids 250m (0-15cm Rhizosphere) & Uncertainty-Buffered WorldClim v2.1 (19 Bioclimatic) layers...")
  if (!isTRUE(demo) && !is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir) && requireNamespace("terra", quietly = TRUE)) {
    tryCatch(crop_and_compress_na_rasters(env_dir = env_dir), error = function(e) NULL)
  }
  if (exists("flag_coordinate_uncertainty", mode = "function")) {
    df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0)
  } else {
    if (!("coordinateUncertaintyInMeters" %in% names(df))) df$coordinateUncertaintyInMeters <- NA_real_
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    df$high_coord_uncertainty_flag <- !is.na(u_m) & u_m > 10000.0
    df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
    df$buffer_radius_m <- ifelse(is.na(u_m) | !is.finite(u_m) | u_m <= 0, 250.0, pmax(100.0, pmin(1000.0, u_m)))
  }

  # Cache-first SQLite environmental_cache pre-hydration in canonical non-demo runs
  req_wc_sg_cols <- c(unname(WC_BIO_LOOKUP), "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")
  full_df_ref <- NULL
  uncached_env_idx <- seq_len(nrow(df))
  is_canon_env <- !isTRUE(demo) && nrow(df) > 0L && is_canonical_prod_path(env_dir) &&
                  exists("query_environmental_cache_sqlite", mode = "function")
  if (is_canon_env) {
    ll_pre <- ifelse(
      is.finite(suppressWarnings(as.numeric(df$latitude))) & is.finite(suppressWarnings(as.numeric(df$longitude))),
      sprintf("%.5f_%.5f", round(as.numeric(df$latitude), 5), round(as.numeric(df$longitude), 5)),
      NA_character_
    )
    sql_pre <- tryCatch(
      query_environmental_cache_sqlite(
        catalog_numbers = as.character(df$catalogNumber),
        lat_lon_keys = ll_pre,
        cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv"),
        demo = demo
      ),
      error = function(e) NULL
    )
    if (!is.null(sql_pre) && is.data.frame(sql_pre) && nrow(sql_pre) > 0L) {
      m_cat_p <- if ("catalogNumber" %in% names(sql_pre)) match(as.character(df$catalogNumber), as.character(sql_pre$catalogNumber)) else rep(NA_integer_, nrow(df))
      m_ll_p <- if ("lat_lon_key" %in% names(sql_pre)) match(ll_pre, as.character(sql_pre$lat_lon_key)) else rep(NA_integer_, nrow(df))
      m_hit_p <- ifelse(!is.na(m_cat_p), m_cat_p, m_ll_p)
      hyd_p <- intersect(c(req_wc_sg_cols, "bio18_precip_warmest_qt", "bio19_precip_coldest_qt"), names(sql_pre))
      for (hc in hyd_p) {
        v_s <- suppressWarnings(as.numeric(sql_pre[[hc]][m_hit_p]))
        if (hc %in% names(df)) {
          v_c <- suppressWarnings(as.numeric(df[[hc]]))
          df[[hc]] <- ifelse(is.finite(v_c), v_c, v_s)
        } else {
          df[[hc]] <- v_s
        }
      }
      df <- coalesce_legacy_wc_bio_cols(df, ensure_all_19 = TRUE)
      if (all(req_wc_sg_cols %in% names(df))) {
        fin_mat <- vapply(req_wc_sg_cols, function(cl) is.finite(suppressWarnings(as.numeric(df[[cl]]))), logical(nrow(df)))
        if (!is.matrix(fin_mat)) fin_mat <- matrix(fin_mat, nrow = nrow(df))
        valid_ll <- is.finite(suppressWarnings(as.numeric(df$latitude))) & is.finite(suppressWarnings(as.numeric(df$longitude)))
        uncached_mask <- valid_ll & (rowSums(!fin_mat) > 0L)
        if (!any(uncached_mask) && any(valid_ll)) {
          message(sprintf("Hydrated all %d WorldClim v2.1 (bio1..bio19) & SoilGrids 250m records directly from SQLite environmental_cache (0 uncached rasters to extract).", sum(valid_ll)))
          return(df)
        }
        if (sum(uncached_mask) < nrow(df) && sum(uncached_mask) > 0L) {
          message(sprintf("SQLite environmental_cache hit for %d / %d records; extracting rasters for %d uncached coordinates...", nrow(df) - sum(uncached_mask), nrow(df), sum(uncached_mask)))
          full_df_ref <- df
          uncached_env_idx <- which(uncached_mask)
          df <- df[uncached_env_idx, , drop = FALSE]
        }
      }
    }
  }
  .finalize_env_extract <- function(extracted_sub) {
    if (!is.null(full_df_ref)) {
      out_full <- full_df_ref
      for (cl in setdiff(names(extracted_sub), c("catalogNumber", "latitude", "longitude"))) {
        if (!(cl %in% names(out_full))) out_full[[cl]] <- rep(NA, nrow(out_full))
        out_full[[cl]][uncached_env_idx] <- extracted_sub[[cl]]
      }
    } else {
      out_full <- extracted_sub
    }
    if (is_canon_env && exists("upsert_environmental_cache_sqlite", mode = "function")) {
      tryCatch({
        up_cols <- intersect(c("catalogNumber", "latitude", "longitude", req_wc_sg_cols), names(extracted_sub))
        if (all(c("catalogNumber", "latitude", "longitude") %in% up_cols)) {
          up_df <- extracted_sub[, up_cols, drop = FALSE]
          up_df$lat_lon_key <- sprintf("%.5f_%.5f", round(as.numeric(up_df$latitude), 5), round(as.numeric(up_df$longitude), 5))
          upsert_environmental_cache_sqlite(up_df, cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv"), demo = demo)
        }
      }, error = function(e) NULL)
    }
    out_full
  }

  # Check if GeoTIFF rasters exist directly in env_dir
  raster_files <- discover_env_rasters(env_dir, recursive = FALSE, crop_global = TRUE)
  extract_err <- NULL
  if (length(raster_files) > 0 && requireNamespace("terra", quietly = TRUE)) {
    extracted_df <- tryCatch({
      raster_files <- crop_and_compress_na_rasters(raster_files, env_dir = env_dir)
      env_stack <- crop_and_compress_na_rasters(terra::rast(raster_files), env_dir = env_dir)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      extracted <- extract_rast_buffered_mean(env_stack, pts, df$buffer_radius_m)
      for (v in setdiff(names(extracted), "ID")) {
        if (!(v %in% unname(WC_BIO_LOOKUP))) df[[v]] <- extracted[[v]]
      }
      df <- assign_wc_bio_from_extracted(df, extracted)
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
      if (!("soil_ph" %in% names(df)) && any(grepl("phh2o", names(extracted), ignore.case = TRUE))) {
        df$soil_ph <- extract_soilgrids_rhizosphere_col(extracted, "phh2o")
      }
      if (!("soil_cec" %in% names(df)) && any(grepl("cec", names(extracted), ignore.case = TRUE))) {
        df$soil_cec <- extract_soilgrids_rhizosphere_col(extracted, "cec")
      }
      if (!("soil_sand" %in% names(df)) && any(grepl("sand", names(extracted), ignore.case = TRUE))) {
        df$soil_sand <- extract_soilgrids_rhizosphere_col(extracted, "sand")
      }
      if (!("soil_bulk_density" %in% names(df)) && any(grepl("bdod", names(extracted), ignore.case = TRUE))) {
        df$soil_bulk_density <- extract_soilgrids_rhizosphere_col(extracted, "bdod")
      }
      message(sprintf("Extracted %d raster layers (with adaptive buffer radius r_i in [100, 1000]m) from %s.", length(raster_files), env_dir))
      df
    }, error = function(e) {
      extract_err <<- e$message
      NULL
    })
    if (!is.null(extracted_df)) return(.finalize_env_extract(extracted_df))
  }

  sub_files <- discover_env_rasters(env_dir, recursive = TRUE, crop_global = TRUE)
  wc_files <- sub_files[grepl("bio[_.-]?[0-9]+", basename(sub_files), ignore.case = TRUE)]
  sg_files <- sub_files[grepl("phh2o|cec|sand|bdod", basename(sub_files), ignore.case = TRUE)]
  if (length(wc_files) > 0 && length(sg_files) > 0 && requireNamespace("terra", quietly = TRUE)) {
    extracted_cached <- tryCatch({
      wc_files <- crop_and_compress_na_rasters(wc_files, env_dir = env_dir)
      sg_files <- crop_and_compress_na_rasters(sg_files, env_dir = env_dir)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      wc_rast <- crop_and_compress_na_rasters(terra::rast(wc_files), env_dir = env_dir)
      sg_rast <- crop_and_compress_na_rasters(terra::rast(sg_files), env_dir = env_dir)
      wc_vals <- extract_rast_buffered_mean(wc_rast, pts, df$buffer_radius_m)
      sg_vals <- extract_rast_buffered_mean(sg_rast, pts, df$buffer_radius_m)
      df <- assign_wc_bio_from_extracted(df, wc_vals)
      df$soil_ph <- extract_soilgrids_rhizosphere_col(sg_vals, "phh2o")
      df$soil_cec <- extract_soilgrids_rhizosphere_col(sg_vals, "cec")
      df$soil_sand <- extract_soilgrids_rhizosphere_col(sg_vals, "sand")
      df$soil_bulk_density <- extract_soilgrids_rhizosphere_col(sg_vals, "bdod")
      message(sprintf("Extracted %d WorldClim (bio1..bio19) and %d SoilGrids cached rasters (0-15cm rhizosphere & buffered) from %s.", length(wc_files), length(sg_files), env_dir))
      df
    }, error = function(e) {
      extract_err <<- e$message
      NULL
    })
    if (!is.null(extracted_cached)) return(.finalize_env_extract(extracted_cached))
  }

  # Automated GeoTIFF Acquisition via geodata when local rasters are absent and --demo is FALSE
  if (!isTRUE(demo) && requireNamespace("geodata", quietly = TRUE) && requireNamespace("terra", quietly = TRUE) &&
      (is_canonical_prod_path(env_dir) || identical(Sys.getenv("PACKERA_ENABLE_GEODATA_DOWNLOAD"), "TRUE"))) {
    extracted_geo <- tryCatch({
      message(sprintf("Attempting automated WorldClim v2.1 (bio1..bio19) & SoilGrids 250m (0-5cm & 5-15cm) acquisition via geodata into '%s'...", env_dir))
      dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
      wc_rast <- geodata::worldclim_global(var = "bio", res = 2.5, path = env_dir)
      wc_rast <- crop_and_compress_na_rasters(wc_rast, env_dir = env_dir)
      sg_0_5 <- lapply(c("phh2o", "cec", "sand", "bdod"), function(v) tryCatch(geodata::soil_world(var = v, depth = 5, path = env_dir), error = function(e) NULL))
      sg_5_15 <- lapply(c("phh2o", "cec", "sand", "bdod"), function(v) geodata::soil_world(var = v, depth = 15, path = env_dir))
      sg_list <- crop_and_compress_na_rasters(c(Filter(Negate(is.null), sg_0_5), sg_5_15), env_dir = env_dir)
      sg_rast <- crop_and_compress_na_rasters(terra::rast(sg_list), env_dir = env_dir)
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      wc_vals <- extract_rast_buffered_mean(wc_rast, pts, df$buffer_radius_m)
      sg_vals <- extract_rast_buffered_mean(sg_rast, pts, df$buffer_radius_m)
      df <- assign_wc_bio_from_extracted(df, wc_vals)
      df$soil_ph <- extract_soilgrids_rhizosphere_col(sg_vals, "phh2o")
      df$soil_cec <- extract_soilgrids_rhizosphere_col(sg_vals, "cec")
      df$soil_sand <- extract_soilgrids_rhizosphere_col(sg_vals, "sand")
      df$soil_bulk_density <- extract_soilgrids_rhizosphere_col(sg_vals, "bdod")
      df
    }, error = function(e) {
      extract_err <<- e$message
      NULL
    })
    if (!is.null(extracted_geo)) return(.finalize_env_extract(extracted_geo))
  }
  if (!is.null(full_df_ref)) df <- full_df_ref

  # Offline SQLite / CSV environmental_cache hydration when --demo is FALSE
  if (!isTRUE(demo) && nrow(df) > 0L) {
    if (exists("query_environmental_cache_sqlite", mode = "function")) {
      ll_keys <- ifelse(
        is.finite(suppressWarnings(as.numeric(df$latitude))) & is.finite(suppressWarnings(as.numeric(df$longitude))),
        sprintf("%.5f_%.5f", round(as.numeric(df$latitude), 5), round(as.numeric(df$longitude), 5)),
        NA_character_
      )
      sql_env <- tryCatch(
        query_environmental_cache_sqlite(
          catalog_numbers = as.character(df$catalogNumber),
          lat_lon_keys = ll_keys,
          cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv"),
          demo = demo
        ),
        error = function(e) NULL
      )
      if (!is.null(sql_env) && is.data.frame(sql_env) && nrow(sql_env) > 0L) {
        m_cat <- if ("catalogNumber" %in% names(sql_env)) match(as.character(df$catalogNumber), as.character(sql_env$catalogNumber)) else rep(NA_integer_, nrow(df))
        m_ll <- if ("lat_lon_key" %in% names(sql_env)) match(ll_keys, as.character(sql_env$lat_lon_key)) else rep(NA_integer_, nrow(df))
        m_hit <- ifelse(!is.na(m_cat), m_cat, m_ll)
        hyd_cols <- intersect(
          c(unname(WC_BIO_LOOKUP), "bio18_precip_warmest_qt", "bio19_precip_coldest_qt",
            "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density",
            "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_soil_cec", "fine_bedrock_depth_cm"),
          names(sql_env)
        )
        for (hc in hyd_cols) {
          v_sql <- suppressWarnings(as.numeric(sql_env[[hc]][m_hit]))
          if (hc %in% names(df)) {
            v_cur <- suppressWarnings(as.numeric(df[[hc]]))
            df[[hc]] <- ifelse(is.finite(v_cur), v_cur, v_sql)
          } else {
            df[[hc]] <- v_sql
          }
        }
      }
    }

    df <- coalesce_legacy_wc_bio_cols(df, ensure_all_19 = TRUE)
    has_cached_clim <- ("bio1_temp_mean" %in% names(df)) && any(is.finite(suppressWarnings(as.numeric(df$bio1_temp_mean))))
    has_cached_soil <- any(vapply(intersect(c("soil_ph", "fine_soil_ph"), names(df)), function(cl) {
      any(is.finite(suppressWarnings(as.numeric(df[[cl]]))))
    }, logical(1)))

    if (has_cached_clim && has_cached_soil) {
      if (!("soil_ph" %in% names(df))) df$soil_ph <- rep(NA_real_, nrow(df))
      if ("fine_soil_ph" %in% names(df)) {
        cur_ph <- suppressWarnings(as.numeric(df$soil_ph))
        fin_ph <- suppressWarnings(as.numeric(df$fine_soil_ph))
        df$soil_ph <- ifelse(is.finite(cur_ph), cur_ph, fin_ph)
      }
      if (!("soil_sand" %in% names(df))) df$soil_sand <- rep(NA_real_, nrow(df))
      if ("fine_soil_sand" %in% names(df)) {
        cur_s <- suppressWarnings(as.numeric(df$soil_sand))
        fin_s <- suppressWarnings(as.numeric(df$fine_soil_sand))
        df$soil_sand <- ifelse(is.finite(cur_s), cur_s, fin_s)
      }
      if (!("soil_cec" %in% names(df))) df$soil_cec <- rep(NA_real_, nrow(df))
      if ("fine_soil_cec" %in% names(df)) {
        cur_c <- suppressWarnings(as.numeric(df$soil_cec))
        fin_c <- suppressWarnings(as.numeric(df$fine_soil_cec))
        df$soil_cec <- ifelse(is.finite(cur_c), cur_c, fin_c)
      }
      if (!("soil_bulk_density" %in% names(df))) df$soil_bulk_density <- rep(NA_real_, nrow(df))
      message("Hydrated WorldClim v2.1 & SoilGrids 250m covariates directly from SQLite environmental_cache (uncached bio* left as NA_real_).")
      return(df)
    }
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

  message("[DEMO MODE] Generating synthetic SoilGrids (0-15cm rhizosphere) & all 19 WorldClim bioclimatic layers...")
  set.seed(42)
  lat <- as.numeric(df$latitude); lon <- as.numeric(df$longitude); n <- nrow(df)
  reg <- if ("regional_group" %in% names(df)) as.character(df$regional_group) else rep("Unknown", n)

  is_flatrock <- isTRUE(demo) & (grepl("Flatrock", reg, ignore.case = TRUE) | (lon > -84 & lon < -79 & lat > 33 & lat < 36.5))
  is_sandhill <- isTRUE(demo) & (grepl("Sandhill|Coastal", reg, ignore.case = TRUE) | (lon > -82 & lon < -75 & lat < 36.5))
  is_prairie  <- isTRUE(demo) & (grepl("Midwest|Prairie", reg, ignore.case = TRUE) | (lon < -88))

  ph_0_5 <- ifelse(is_sandhill, 4.6 + 0.3 * sin(lat), ifelse(is_flatrock, 5.1 + 0.25 * cos(lon),
            ifelse(is_prairie, 7.0 + 0.3 * sin(lat), 6.3 + 0.3 * cos(lat)))) + stats::rnorm(n, 0, 0.15)
  ph_5_15 <- ph_0_5 + 0.15 + stats::rnorm(n, 0, 0.10)
  df$soil_ph <- round(pmax(3.8, pmin(8.5, (5.0 * ph_0_5 + 10.0 * ph_5_15) / 15.0)), 2)

  cec_0_5 <- ifelse(is_sandhill, 5.6 + 1.2 * abs(sin(lat)), ifelse(is_flatrock, 12.0 + 2.0 * cos(lat),
             ifelse(is_prairie, 27.2 + 3.5 * sin(lon), 18.6 + 2.5 * cos(lat)))) + stats::rnorm(n, 0, 1.0)
  cec_5_15 <- cec_0_5 - 0.8 + stats::rnorm(n, 0, 0.8)
  df$soil_cec <- round(pmax(1.5, pmin(50.0, (5.0 * cec_0_5 + 10.0 * cec_5_15) / 15.0)), 1)

  sand_0_5 <- ifelse(is_sandhill, 83.0 + 5.0 * sin(lon), ifelse(is_flatrock, 53.0 + 6.0 * cos(lat),
              ifelse(is_prairie, 24.5 + 4.5 * sin(lat), 38.5 + 5.0 * cos(lon)))) + stats::rnorm(n, 0, 3.0)
  sand_5_15 <- sand_0_5 - 1.5 + stats::rnorm(n, 0, 2.5)
  df$soil_sand <- round(pmax(5.0, pmin(98.0, (5.0 * sand_0_5 + 10.0 * sand_5_15) / 15.0)), 1)

  bd_0_5 <- ifelse(is_sandhill, 1.44, ifelse(is_flatrock, 1.28, 1.20)) + stats::rnorm(n, 0, 0.04)
  bd_5_15 <- bd_0_5 + 0.06 + stats::rnorm(n, 0, 0.04)
  df$soil_bulk_density <- round((5.0 * bd_0_5 + 10.0 * bd_5_15) / 15.0, 2)

  # All 19 WorldClim v2.1 Bioclimatic Predictors (WC_BIO_LOOKUP)
  df$bio1_temp_mean            <- round(24.5 - 0.72 * (lat - 28.0) - 0.05 * abs(lon + 82) + stats::rnorm(n, 0, 0.4), 2)
  df$bio2_mean_diurnal_range   <- round(pmax(4.5, 11.4 + 0.11 * abs(lon + 80.0) + stats::rnorm(n, 0, 0.45)), 2)
  df$bio4_temp_seasonality     <- round(pmax(80.0, 520.0 + 28.5 * (lat - 30.0) + stats::rnorm(n, 0, 15.0)), 1)
  seas_amp                     <- (df$bio4_temp_seasonality / 100.0) * 1.25
  df$bio5_max_temp_warmest_month <- round(df$bio1_temp_mean + seas_amp + 0.55 * df$bio2_mean_diurnal_range + stats::rnorm(n, 0, 0.35), 2)
  df$bio6_min_temp_coldest_month <- round(df$bio1_temp_mean - seas_amp - 0.55 * df$bio2_mean_diurnal_range + stats::rnorm(n, 0, 0.35), 2)
  df$bio7_temp_annual_range    <- round(pmax(df$bio2_mean_diurnal_range + 1.0, df$bio5_max_temp_warmest_month - df$bio6_min_temp_coldest_month), 2)
  df$bio3_isothermality        <- round(100.0 * df$bio2_mean_diurnal_range / pmax(df$bio7_temp_annual_range, 1.0), 2)
  df$bio10_mean_temp_warmest_qtr <- round(df$bio1_temp_mean + 0.95 * (df$bio4_temp_seasonality / 100.0) + stats::rnorm(n, 0, 0.3), 2)
  df$bio11_mean_temp_coldest_qtr <- round(df$bio1_temp_mean - 0.95 * (df$bio4_temp_seasonality / 100.0) + stats::rnorm(n, 0, 0.3), 2)
  df$bio8_mean_temp_wettest_qtr  <- round(0.65 * df$bio10_mean_temp_warmest_qtr + 0.35 * df$bio1_temp_mean + stats::rnorm(n, 0, 0.6), 2)
  df$bio9_mean_temp_driest_qtr   <- round(0.45 * df$bio11_mean_temp_coldest_qtr + 0.55 * df$bio1_temp_mean + stats::rnorm(n, 0, 0.6), 2)

  df$bio12_precip_annual       <- round(pmax(120.0, 1450.0 - 15.0 * (lat - 30.0) + 12.0 * (lon + 85.0) + stats::rnorm(n, 0, 35.0)), 1)
  df$bio15_precip_seasonality  <- round(pmax(5.0, 22.0 + 0.8 * (lat - 30.0) + stats::rnorm(n, 0, 1.5)), 1)
  df$bio13_precip_wettest_month <- round(pmax(15.0, df$bio12_precip_annual * (0.115 + 0.0012 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 6.0)), 1)
  df$bio14_precip_driest_month  <- round(pmax(1.0, df$bio12_precip_annual * pmax(0.015, 0.068 - 0.0009 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 4.0)), 1)
  df$bio16_precip_wettest_qtr  <- round(pmax(df$bio13_precip_wettest_month, df$bio12_precip_annual * (0.31 + 0.0020 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 12.0)), 1)
  df$bio17_precip_driest_qtr   <- round(pmax(df$bio14_precip_driest_month, df$bio12_precip_annual * pmax(0.06, 0.21 - 0.0018 * df$bio15_precip_seasonality) + stats::rnorm(n, 0, 10.0)), 1)
  df$bio18_precip_warmest_qtr  <- round(pmax(10.0, df$bio12_precip_annual * (0.32 + 0.02 * sin(lat)) + stats::rnorm(n, 0, 15.0)), 1)
  df$bio19_precip_coldest_qtr  <- round(pmax(10.0, df$bio12_precip_annual * (0.22 - 0.01 * cos(lon)) + stats::rnorm(n, 0, 12.0)), 1)
  return(coalesce_legacy_wc_bio_cols(df, ensure_all_19 = TRUE))
}

# ------------------------------------------------------------------------------
# 2b. USDA SSURGO Vector Pedology Query & Caching Engine (soilDB / SDA Integration)
# ------------------------------------------------------------------------------
classify_ssurgo_outcrop <- function(
  muname = "",
  compname = "",
  taxsuborder = "",
  taxgreatgroup = "",
  compkind = "",
  lat = NA,
  lon = NA,
  is_demo = FALSE,
  demo = NULL,
  reskind = "",
  soil_ph = NA_real_,
  soil_sand = NA_real_,
  soil_cec = NA_real_,
  bedrock_depth_cm = NA_real_,
  rock_outcrop_pct = NA_real_
) {
  # Coordinate-free SSURGO outcrop & lithic classification from mapunit/component
  # text descriptors and diagnostic physical soil/bedrock indicators.
  # Regional lat/lon coordinate-box inference is restricted strictly to demo == TRUE.
  eff_demo <- if (!is.null(demo)) isTRUE(demo) else isTRUE(is_demo)
  desc <- tolower(trimws(paste(muname, compname, taxsuborder, taxgreatgroup, compkind, reskind)))
  lithic_series_regex <- "(?i)(rock outcrop|wake|louisburg|ashlar|rion|wateree|wilkes|iredell|mecklenburg|poindexter|enon|goldston|manteo|boulders|granite|diabase)"
  is_lithic_series <- grepl(lithic_series_regex, paste(muname, compname), perl = TRUE)
  is_misc_rock <- grepl("miscellaneous area", tolower(as.character(compkind))) &&
                  grepl("rock|rubble|outcrop|boulder", desc)

  ph_num   <- suppressWarnings(as.numeric(soil_ph)[1])
  sand_num <- suppressWarnings(as.numeric(soil_sand)[1])
  cec_num  <- suppressWarnings(as.numeric(soil_cec)[1])
  bed_num  <- suppressWarnings(as.numeric(bedrock_depth_cm)[1])
  ro_num   <- suppressWarnings(as.numeric(rock_outcrop_pct)[1])
  has_shallow_phys <- (is.finite(ro_num) && ro_num >= 15.0) || (is.finite(bed_num) && bed_num <= 50.0)

  is_serp_txt  <- grepl("chrome|trego|conowingo|aldino|serpentine|ultramafic|peridotite|dunite|soapstone", desc)
  is_mafic_txt <- grepl("iredell|mecklenburg|wilkes|enon|poindexter|diabase|gabbro|basalt|amphibolite|trap rock|greenstone|mafic", desc)
  is_lime_txt  <- grepl("gladeville|barfield|cedars|talbott|colbert|limestone|dolomite|dolostone|calcareous|alvar|karst|marble|chalk|rendolls|calciustolls|calciudolls|glade", desc)
  is_sand_txt  <- grepl("ramsey|lily|sewanee|goldston|manteo|wallen|dekab|sandstone|quartzite|siliceous|chert|novaculite|shale|argillite", desc)
  is_gran_txt  <- grepl("wake|louisburg|flatrock|granit|gneiss|felsic|catao|hempstead|ashlar|rion|wateree|pegmatite|rhyolite", desc) ||
                  (grepl("appling", desc) && grepl("flatrock", desc))

  demo_gran <- FALSE
  demo_lime <- FALSE
  if (isTRUE(eff_demo)) {
    la <- suppressWarnings(as.numeric(lat)[1])
    lo <- suppressWarnings(as.numeric(lon)[1])
    demo_gran <- grepl("dystrudepts|rock outcrop", desc) &&
                 !is.na(lo) && lo > -84.5 && lo < -77.5 && !is.na(la) && la > 32.5 && la < 37.0
    demo_lime <- grepl("hapludolls|rendolls", desc) &&
                 !is.na(lo) && lo > -88.5 && lo < -83.0 && !is.na(la) && la > 34.5 && la < 37.5
  }

  is_outcrop <- grepl("rock outcrop|flatrock|glade|lithic|paralithic|bedrock|rubble|boulder", desc) ||
                is_lithic_series || is_misc_rock || has_shallow_phys ||
                is_serp_txt || is_mafic_txt || is_lime_txt || is_sand_txt || is_gran_txt ||
                demo_gran || demo_lime

  outcrop_class <- if (is_lime_txt || demo_lime) {
    "Limestone glade"
  } else if (is_sand_txt) {
    "Sandstone/Siliceous"
  } else if (is_serp_txt) {
    "Serpentine"
  } else if (is_mafic_txt) {
    "Diabase/Mafic glade"
  } else if (is_gran_txt || demo_gran) {
    "Granitic flatrock"
  } else if (is_outcrop) {
    if (is.finite(ph_num) && ph_num >= 6.8) {
      "Limestone glade"
    } else if (is.finite(cec_num) && cec_num >= 20.0 && (!is.finite(sand_num) || sand_num < 50.0)) {
      "Diabase/Mafic glade"
    } else if (is.finite(sand_num) && sand_num >= 60.0 && (!is.finite(ph_num) || ph_num <= 5.5)) {
      "Granitic flatrock"
    } else {
      "Sandstone/Siliceous"
    }
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
  if (!isTRUE(demo) && exists("query_environmental_cache_sqlite", mode = "function")) {
    b_ll <- sprintf("%.5f_%.5f", batch_df$latitude, batch_df$longitude)
    sql_ss <- tryCatch(query_environmental_cache_sqlite(catalog_numbers = batch_df$catalogNumber, lat_lon_keys = b_ll, cache_path = cache_dir, demo = demo), error = function(e) NULL)
    if (!is.null(sql_ss) && nrow(sql_ss) > 0L && any(c("musym", "outcrop_class") %in% names(sql_ss))) {
      m_cat <- match(as.character(batch_df$catalogNumber), as.character(sql_ss$catalogNumber))
      m_ll <- if ("lat_lon_key" %in% names(sql_ss)) match(b_ll, as.character(sql_ss$lat_lon_key)) else rep(NA_integer_, nrow(batch_df))
      m_hit <- ifelse(!is.na(m_cat), m_cat, m_ll)
      if (all(!is.na(m_hit))) {
        hit_df <- sql_ss[m_hit, , drop = FALSE]; hit_df$catalogNumber <- as.character(batch_df$catalogNumber)
        hit_df$latitude <- batch_df$latitude; hit_df$longitude <- batch_df$longitude
        if (("musym" %in% names(hit_df) && any(!is.na(hit_df$musym))) || ("ssurgo_rock_outcrop_pct" %in% names(hit_df) && any(!is.na(hit_df$ssurgo_rock_outcrop_pct)))) {
          message(sprintf("Loaded SSURGO batch %d directly from SQLite environmental_cache.", batch_idx))
          return(normalize_ssurgo_cols(hit_df))
        }
      }
    }
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
              agg_df <- aggregate_ssurgo_tabular_to_points(sub_pts, buf_tab, classify_fn = classify_ssurgo_outcrop, is_demo = isTRUE(demo))
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
                  agg_df <- aggregate_ssurgo_tabular_to_points(sub_pts, tab_joined, classify_fn = classify_ssurgo_outcrop, is_demo = isTRUE(demo))
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
    
    is_flatrock_zone <- isTRUE(demo) && (grepl("Flatrock", reg, ignore.case = TRUE) ||
                        (lon > -84.5 && lon < -78.5 && lat > 33.2 && lat < 36.8))
    is_sandhill_zone <- isTRUE(demo) && (grepl("Sandhill|Coastal", reg, ignore.case = TRUE) ||
                        (lon > -82.0 && lon < -75.5 && lat < 36.5))
    is_prairie_zone  <- isTRUE(demo) && (grepl("Midwest|Prairie", reg, ignore.case = TRUE) || (lon < -88.0))
    is_appalachian   <- isTRUE(demo) && (grepl("Appalachian", reg, ignore.case = TRUE) ||
                        (lon >= -84.0 && lon <= -79.0 && lat >= 35.0 && lat <= 39.0))
    
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
    
    cl <- classify_ssurgo_outcrop(munames[i], "", suborders[i], greatgroups[i], "", lat, lon, is_demo = TRUE)
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
    batch_results <- list()
    rem_conus <- conus_df

    # Single-batch SQLite environmental_cache pre-hydration across all CONUS points
    if (!isTRUE(demo) && exists("query_environmental_cache_sqlite", mode = "function")) {
      c_ll <- sprintf("%.5f_%.5f", round(as.numeric(conus_df$latitude), 5), round(as.numeric(conus_df$longitude), 5))
      sql_all <- tryCatch(
        query_environmental_cache_sqlite(
          catalog_numbers = as.character(conus_df$catalogNumber),
          lat_lon_keys = c_ll,
          cache_path = cache_dir,
          demo = demo
        ),
        error = function(e) NULL
      )
      if (!is.null(sql_all) && nrow(sql_all) > 0L && any(c("musym", "ssurgo_rock_outcrop_pct") %in% names(sql_all))) {
        m_cat <- if ("catalogNumber" %in% names(sql_all)) match(as.character(conus_df$catalogNumber), as.character(sql_all$catalogNumber)) else rep(NA_integer_, n_conus)
        m_ll <- if ("lat_lon_key" %in% names(sql_all)) match(c_ll, as.character(sql_all$lat_lon_key)) else rep(NA_integer_, n_conus)
        m_hit <- ifelse(!is.na(m_cat), m_cat, m_ll)
        has_ss_hit <- !is.na(m_hit)
        if (any(has_ss_hit)) {
          cand_hit <- sql_all[m_hit[has_ss_hit], , drop = FALSE]
          valid_ss <- (if ("musym" %in% names(cand_hit)) !is.na(cand_hit$musym) & nzchar(as.character(cand_hit$musym)) else rep(FALSE, nrow(cand_hit))) |
                      (if ("ssurgo_rock_outcrop_pct" %in% names(cand_hit)) !is.na(cand_hit$ssurgo_rock_outcrop_pct) else rep(FALSE, nrow(cand_hit)))
          if (any(valid_ss)) {
            hit_idx <- which(has_ss_hit)[valid_ss]
            hit_rows <- cand_hit[valid_ss, , drop = FALSE]
            hit_rows$catalogNumber <- as.character(conus_df$catalogNumber[hit_idx])
            hit_rows$latitude <- conus_df$latitude[hit_idx]
            hit_rows$longitude <- conus_df$longitude[hit_idx]
            batch_results[[1L]] <- normalize_ssurgo_cols(hit_rows)
            rem_conus <- conus_df[-hit_idx, , drop = FALSE]
            message(sprintf("Pre-hydrated %d / %d CONUS SSURGO records directly from SQLite environmental_cache.", length(hit_idx), n_conus))
          }
        }
      }
    }

    if (nrow(rem_conus) > 0L) {
      batch_indices <- split(seq_len(nrow(rem_conus)), ceiling(seq_len(nrow(rem_conus)) / batch_size))
      run_one_ss_batch <- function(b_idx) {
        idx <- batch_indices[[b_idx]]
        sub_df <- rem_conus[idx, , drop = FALSE]
        query_ssurgo_sda_batch(sub_df, cache_dir, b_idx, demo = demo, batch_size = batch_size)
      }
      rem_batches <- if (!isTRUE(demo) && length(batch_indices) >= 2L && exists("packera_mclapply", mode = "function")) {
        packera_mclapply(seq_along(batch_indices), run_one_ss_batch, mc.cores = min(4L, length(batch_indices)), mc.preschedule = FALSE)
      } else {
        lapply(seq_along(batch_indices), run_one_ss_batch)
      }
      batch_results <- c(batch_results, rem_batches)
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

# Source modular phenology (model_latitudinal_spring_baseline, hopkins_adjusted_doy, delta_doy, fill = NA) and ecotype/type-locality helpers
step05_dir <- if (dir.exists("scripts/analysis")) {
  "scripts/analysis"
} else {
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
}
source(file.path(step05_dir, "step05_phenology_ecotypes.R"))

execute_crossmodal_consensus <- function(
  df,
  intraspecific_mode = FALSE,
  rf_res = NULL,
  primary_env_vars = NULL,
  primary_clim_vars = NULL,
  primary_edaph_vars = NULL,
  rf_var_weights = NULL
) {
  message("Executing 3-way cross-modal consensus checks (Morpho + Fine-Grain Edaphic + Phenology)...")
  if (exists("flag_coordinate_uncertainty", mode = "function")) {
    df <- flag_coordinate_uncertainty(df, max_uncertainty_m = 10000.0)
  } else {
    if (!("coordinateUncertaintyInMeters" %in% names(df))) df$coordinateUncertaintyInMeters <- NA_real_
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    df$high_coord_uncertainty_flag <- !is.na(u_m) & u_m > 10000.0
    df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
    df$buffer_radius_m <- ifelse(is.na(u_m) | !is.finite(u_m) | u_m <= 0, 250.0, pmax(100.0, pmin(1000.0, u_m)))
  }

  df <- coalesce_legacy_wc_bio_cols(df, ensure_all_19 = FALSE)
  df <- resolve_analysis_group(
    df,
    intraspecific_mode = intraspecific_mode,
    primary_env_vars = primary_env_vars,
    rf_var_weights = rf_var_weights
  )
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  active_taxa <- if (is_intra) sort(unique(stats::na.omit(df$analysis_group))) else get_active_packera_taxa(df$species_standardized)
  if (length(active_taxa) == 0L) active_taxa <- sort(unique(as.character(stats::na.omit(df$analysis_group))))
  grp_vec <- if (is_intra) df$analysis_group else df$species_standardized

  n <- nrow(df)

  # Coalesce fine-grain soil aliases from SoilGrids/SSURGO if missing without fabricating synthetic constants
  if (!("fine_soil_ph" %in% names(df))) df$fine_soil_ph <- if ("soil_ph" %in% names(df)) df$soil_ph else rep(NA_real_, n)
  if (!("fine_soil_sand" %in% names(df))) df$fine_soil_sand <- if ("soil_sand" %in% names(df)) df$soil_sand else rep(NA_real_, n)
  if (!("fine_soil_clay" %in% names(df))) {
    df$fine_soil_clay <- if (any(is.finite(as.numeric(df$fine_soil_sand)))) {
      ifelse(is.finite(as.numeric(df$fine_soil_sand)), pmax(2.0, pmin(65.0, 100.0 - as.numeric(df$fine_soil_sand) - 22.0)), NA_real_)
    } else {
      rep(NA_real_, n)
    }
  }
  if (!("fine_bedrock_depth_cm" %in% names(df))) {
    df$fine_bedrock_depth_cm <- if ("ssurgo_bedrock_depth_cm" %in% names(df)) df$ssurgo_bedrock_depth_cm else rep(NA_real_, n)
  }
  if (!("ssurgo_rock_outcrop_pct" %in% names(df))) {
    df$ssurgo_rock_outcrop_pct <- if ("is_rock_outcrop" %in% names(df)) ifelse(as.logical(df$is_rock_outcrop), 45.0, 0.0) else rep(NA_real_, n)
  }
  if (!("fine_soil_cec" %in% names(df))) df$fine_soil_cec <- if ("soil_cec" %in% names(df)) df$soil_cec else rep(NA_real_, n)

  # Phenological circular stats (cohort circular mean fallback instead of hardcoded 120) & Empirical Bayes Authority Weights
  raw_doy <- if ("doy" %in% names(df)) suppressWarnings(as.numeric(df$doy)) else rep(NA_real_, n)
  emp_doy_ref <- compute_circular_mean_doy(raw_doy)
  if (!is.finite(emp_doy_ref)) {
    cfg_b <- read_config_phenology_bounds()
    emp_doy_ref <- 0.5 * (cfg_b$min_flowering_doy + cfg_b$max_flowering_doy)
  }
  doy <- ifelse(!is.finite(raw_doy) | raw_doy < 1 | raw_doy > 366, round(emp_doy_ref), raw_doy)
  theta <- 2 * pi * doy / 365.25
  df$pheno_theta <- theta
  df$authority_weight <- resolve_step06_authority_weights(df)
  df$spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df) else rep(1.0, nrow(df))
  df$case_wts <- df$authority_weight * df$spatial_balance_weight

  if ((is.null(primary_env_vars) || length(primary_env_vars) == 0L) && !is.null(rf_res) && is.list(rf_res) && !is.null(rf_res$importance)) {
    sel_auto <- select_primary_predictors_from_rf(
      rf_res = rf_res,
      vif_res = rf_res$collinearity_screening,
      df_env = df,
      intraspecific_mode = is_intra
    )
    if (length(sel_auto$primary_env_vars) > 0L) {
      primary_env_vars <- sel_auto$primary_env_vars
      if (is.null(primary_clim_vars) || length(primary_clim_vars) == 0L) primary_clim_vars <- sel_auto$primary_clim_vars
      if (is.null(primary_edaph_vars) || length(primary_edaph_vars) == 0L) primary_edaph_vars <- sel_auto$primary_edaph_vars
      if (is.null(rf_var_weights) || length(rf_var_weights) == 0L) rf_var_weights <- sel_auto$rf_var_weights
    }
  }

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
  # Strictly on Tier_1_Gold, Tier_2_Silver & promoted Bronze anchors with coordinateUncertaintyInMeters <= 10000
  elig_coord <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, nrow(df))
  elig_coord[is.na(elig_coord)] <- TRUE
  is_prom_anc <- rep(FALSE, n)
  if ("is_promoted_anchor" %in% names(df)) { pa <- as.logical(df$is_promoted_anchor); pa[is.na(pa)] <- FALSE; is_prom_anc <- is_prom_anc | pa }
  if ("is_promoted_bronze" %in% names(df)) { pb <- as.logical(df$is_promoted_bronze); pb[is.na(pb)] <- FALSE; is_prom_anc <- is_prom_anc | pb }
  det_t_vec <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n)
  is_anchor <- ((det_t_vec %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_prom_anc) & elig_coord & grp_vec %in% active_taxa
  anchor_df <- df[is_anchor, ]
  if (nrow(anchor_df) < min(10L, n)) anchor_df <- df[elig_coord & grp_vec %in% active_taxa, ]
  if (nrow(anchor_df) < min(10L, n)) anchor_df <- df[grp_vec %in% active_taxa, ]

  cand_eco <- discover_candidate_eco_predictors(df, include_pheno = TRUE, min_finite = 3L)
  eco_cols <- if (!is.null(primary_env_vars) && length(intersect(primary_env_vars, cand_eco)) > 0L) {
    intersect(as.character(primary_env_vars), cand_eco)
  } else {
    cand_eco
  }
  if (length(eco_cols) == 0L) {
    fallback_c <- intersect(c("fine_soil_ph", "fine_soil_sand", "fine_soil_cec", "soil_ph", "soil_sand", "soil_cec", unname(WC_BIO_LOOKUP)), names(df))
    eco_cols <- fallback_c[vapply(fallback_c, function(cl) any(is.finite(suppressWarnings(as.numeric(df[[cl]])))), logical(1))]
  }
  min_sd_map <- stats::setNames(vapply(eco_cols, function(cl) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    vf <- v[is.finite(v)]
    if (length(vf) < 2L) return(1e-3)
    disp_g <- max(c(stats::sd(vf), stats::IQR(vf) / 1.349), na.rm = TRUE)
    if (!is.finite(disp_g) || disp_g <= 1e-8) 1e-3 else max(0.15 * disp_g, 1e-4)
  }, numeric(1)), eco_cols)

  taxa_profiles <- lapply(active_taxa, function(tx) {
    sub_df <- anchor_df[(if (is_intra) anchor_df$analysis_group else anchor_df$species_standardized) == tx, ]
    if (nrow(sub_df) < 3) sub_df <- df[grp_vec == tx, ]
    w_sub <- sub_df$authority_weight * sub_df$spatial_balance_weight
    sin_m <- w_mean(sin(sub_df$pheno_theta), w_sub)
    cos_m <- w_mean(cos(sub_df$pheno_theta), w_sub)
    eco_m <- sapply(eco_cols, function(cl) w_mean(as.numeric(sub_df[[cl]]), w_sub))
    eco_s <- sapply(eco_cols, function(cl) w_sd(as.numeric(sub_df[[cl]]), w_sub, ifelse(cl %in% names(min_sd_map), min_sd_map[[cl]], 1e-3)))
    list(
      mean_theta = atan2(sin_m, cos_m), R_bar = sqrt(sin_m^2 + cos_m^2),
      eco_mean = eco_m, eco_sd = eco_s
    )
  })
  names(taxa_profiles) <- active_taxa
  for (cl in eco_cols) {
    bar_sd <- mean(sapply(taxa_profiles, function(p) p$eco_sd[[cl]]), na.rm = TRUE)
    for (tx in active_taxa) taxa_profiles[[tx]]$eco_sd[[cl]] <- 0.5 * taxa_profiles[[tx]]$eco_sd[[cl]] + 0.5 * bar_sd
  }

  z_engine <- compute_specimen_edaphic_z(
    df = df,
    primary_edaph_vars = primary_edaph_vars,
    primary_clim_vars = primary_clim_vars,
    primary_env_vars = eco_cols,
    rf_var_weights = rf_var_weights,
    intraspecific_mode = is_intra
  )

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

  # Precompute decorrelated PCA-env + Ledoit-Wolf shrinkage Mahalanobis likelihoods
  # and extract optional Spatial Block-CV out-of-fold Random Forest posteriors
  K_tx <- length(active_taxa)
  maha_prob_mat <- matrix(1.0 / max(1L, K_tx), nrow = n, ncol = max(1L, K_tx),
                          dimnames = list(NULL, if (K_tx > 0L) active_taxa else "Unknown"))
  if (n > 0L && K_tx >= 2L && length(eco_cols) >= 1L) {
    p_eco <- length(eco_cols)
    E_all <- matrix(0.0, nrow = n, ncol = p_eco, dimnames = list(NULL, eco_cols))
    E_anc <- matrix(0.0, nrow = nrow(anchor_df), ncol = p_eco, dimnames = list(NULL, eco_cols))
    mu_scale <- numeric(p_eco); sd_scale <- numeric(p_eco)
    w_anc_all <- anchor_df$authority_weight * anchor_df$spatial_balance_weight
    for (j in seq_len(p_eco)) {
      cl <- eco_cols[j]
      v_all <- suppressWarnings(as.numeric(df[[cl]]))
      v_anc <- suppressWarnings(as.numeric(anchor_df[[cl]]))
      med_j <- stats::median(v_anc[is.finite(v_anc)], na.rm = TRUE)
      if (!is.finite(med_j)) med_j <- stats::median(v_all[is.finite(v_all)], na.rm = TRUE)
      if (!is.finite(med_j)) med_j <- 0.0
      v_all[!is.finite(v_all)] <- med_j
      v_anc[!is.finite(v_anc)] <- med_j
      E_all[, j] <- v_all
      E_anc[, j] <- v_anc
      min_s <- if (cl %in% names(min_sd_map)) min_sd_map[[cl]] else 0.25
      mu_j <- w_mean(v_anc, w_anc_all)
      if (!is.finite(mu_j)) mu_j <- mean(v_all)
      mu_scale[j] <- mu_j
      sd_scale[j] <- max(w_sd(v_anc, w_anc_all, min_s), min_s, 1e-4)
    }
    Z_raw_all <- sweep(sweep(E_all, 2L, mu_scale, "-"), 2L, sd_scale, "/")
    Z_raw_anc <- sweep(sweep(E_anc, 2L, mu_scale, "-"), 2L, sd_scale, "/")
    Z_raw_all[!is.finite(Z_raw_all)] <- 0.0
    Z_raw_anc[!is.finite(Z_raw_anc)] <- 0.0

    if (p_eco >= 2L && nrow(Z_raw_anc) >= 3L) {
      pca_env <- tryCatch(stats::prcomp(Z_raw_anc, center = TRUE, scale. = FALSE), error = function(e) NULL)
    } else {
      pca_env <- NULL
    }
    if (!is.null(pca_env) && ncol(pca_env$rotation) >= 1L) {
      eig_v <- pca_env$sdev^2
      cum_v <- if (sum(eig_v) > 1e-12) cumsum(eig_v) / sum(eig_v) else rep(1.0, length(eig_v))
      k_90 <- which(cum_v >= 0.90)[1]
      if (is.na(k_90)) k_90 <- ncol(pca_env$rotation)
      k_max <- max(1L, min(6L, ncol(pca_env$rotation), max(1L, nrow(Z_raw_anc) - K_tx)))
      k_pca <- max(1L, min(k_90, k_max))
      rot_k <- pca_env$rotation[, seq_len(k_pca), drop = FALSE]
      Z_all <- sweep(Z_raw_all, 2L, pca_env$center, "-") %*% rot_k
      Z_anc <- sweep(Z_raw_anc, 2L, pca_env$center, "-") %*% rot_k
    } else {
      k_pca <- ncol(Z_raw_all)
      Z_all <- Z_raw_all
      Z_anc <- Z_raw_anc
    }

    mu_pca <- matrix(0.0, nrow = K_tx, ncol = k_pca, dimnames = list(active_taxa, NULL))
    res_chunks <- list()
    w_tot_pool <- 0.0
    n_obs_pool <- 0L
    for (tx in active_taxa) {
      idx_a <- which(!is.na(anc_grp_vec) & anc_grp_vec == tx)
      if (length(idx_a) >= 2L) {
        Zk <- Z_anc[idx_a, , drop = FALSE]
        wk <- pmax(1e-3, w_anc_all[idx_a])
      } else {
        idx_f <- which(!is.na(grp_vec) & grp_vec == tx)
        Zk <- Z_all[idx_f, , drop = FALSE]
        wk <- pmax(1e-3, df$case_wts[idx_f])
      }
      if (nrow(Zk) >= 1L) {
        sw <- sum(wk)
        mk <- colSums(Zk * wk) / max(sw, 1e-8)
        mu_pca[tx, ] <- mk
        dk <- sweep(Zk, 2L, mk, "-")
        res_chunks[[length(res_chunks) + 1L]] <- dk * sqrt(wk)
        w_tot_pool <- w_tot_pool + sw
        n_obs_pool <- n_obs_pool + nrow(Zk)
      }
    }
    X_res_pca <- if (length(res_chunks) > 0L) do.call(rbind, res_chunks) else Z_anc
    S_w <- crossprod(X_res_pca) / max(w_tot_pool - K_tx, 1.0)
    S_w[!is.finite(S_w)] <- 0.0
    tr_sw <- sum(diag(S_w))
    target_var <- if (is.finite(tr_sw) && tr_sw > 1e-8) (tr_sw / k_pca) else 1.0
    if (k_pca == 1L) {
      lam_lw <- 0.10
    } else {
      tr_sw2 <- sum(S_w * S_w)
      denom_lw <- (max(n_obs_pool, 3L) + 1.0 - 2.0 / k_pca) * (tr_sw2 - (tr_sw * tr_sw) / k_pca)
      num_lw <- (1.0 - 2.0 / k_pca) * tr_sw2 + (tr_sw * tr_sw)
      lam_raw <- if (is.finite(denom_lw) && denom_lw > 1e-12) (num_lw / denom_lw) else 0.15
      lam_lw <- min(0.50, max(0.05, lam_raw))
    }
    S_reg_pca <- (1.0 - lam_lw) * S_w + lam_lw * diag(target_var, nrow = k_pca, ncol = k_pca)
    S_inv_pca <- tryCatch(
      solve(S_reg_pca),
      error = function(e) diag(1.0 / pmax(diag(S_reg_pca), 1e-4), nrow = k_pca, ncol = k_pca)
    )

    d2_maha <- matrix(0.0, nrow = n, ncol = K_tx, dimnames = list(NULL, active_taxa))
    for (tx in active_taxa) {
      diff_tx <- sweep(Z_all, 2L, mu_pca[tx, ], "-")
      d2_maha[, tx] <- pmax(0.0, rowSums((diff_tx %*% S_inv_pca) * diff_tx))
    }
    nu_df <- max(3.0, k_pca + 2.0)
    ll_maha <- -0.5 * (nu_df + k_pca) * log1p(d2_maha / nu_df)
    ll_max <- apply(ll_maha, 1L, max)
    p_maha_unnorm <- exp(sweep(ll_maha, 1L, ll_max, "-"))
    maha_prob_mat <- sweep(p_maha_unnorm, 1L, pmax(rowSums(p_maha_unnorm), 1e-12), "/")
  }

  rf_prob_mat <- NULL
  rf_acc_est <- NA_real_
  if (!is.null(rf_res) && is.list(rf_res)) {
    for (acc_nm in c("Spatial_Block_CV_Accuracy", "Spatial_Block_CV_Acc", "spatial_block_cv_acc", "block_cv_acc", "OOB_Accuracy")) {
      if (!is.null(rf_res[[acc_nm]]) && is.finite(suppressWarnings(as.numeric(rf_res[[acc_nm]][1])))) {
        rf_acc_est <- suppressWarnings(as.numeric(rf_res[[acc_nm]][1]))
        break
      }
    }
    cand_mat <- if (!is.null(rf_res$oof_prob_matrix)) rf_res$oof_prob_matrix else rf_res$full_prob_mat
    if (is.matrix(cand_mat) && nrow(cand_mat) == n && K_tx >= 2L && length(intersect(colnames(cand_mat), active_taxa)) >= 1L) {
      rf_prob_mat <- matrix(NA_real_, nrow = n, ncol = K_tx, dimnames = list(NULL, active_taxa))
      common_tx <- intersect(colnames(cand_mat), active_taxa)
      for (ctx in common_tx) {
        rf_prob_mat[, ctx] <- suppressWarnings(as.numeric(cand_mat[, ctx]))
      }
    } else {
      rf_p_tx <- if (!is.null(rf_res$rf_predicted_taxon)) rf_res$rf_predicted_taxon else rf_res$predicted_taxon
      rf_p_pr <- if (!is.null(rf_res$rf_predicted_prob)) rf_res$rf_predicted_prob else rf_res$predicted_prob
      rf_g_pr <- rf_res$rf_given_taxon_prob
      if (!is.null(rf_p_tx) && !is.null(rf_p_pr) && length(rf_p_tx) == n && length(rf_p_pr) == n && K_tx >= 2L) {
        rf_prob_mat <- matrix(NA_real_, nrow = n, ncol = K_tx, dimnames = list(NULL, active_taxa))
        for (ri in seq_len(n)) {
          ptx <- as.character(rf_p_tx[ri]); ppr <- suppressWarnings(as.numeric(rf_p_pr[ri]))
          if (!is.na(ptx) && (ptx %in% active_taxa) && is.finite(ppr)) {
            ppr <- min(0.999, max(1e-4, ppr))
            gtx <- as.character(grp_vec[ri])
            gpr <- if (!is.null(rf_g_pr) && length(rf_g_pr) == n) suppressWarnings(as.numeric(rf_g_pr[ri])) else NA_real_
            if (!is.na(gtx) && (gtx %in% active_taxa) && gtx != ptx && is.finite(gpr)) {
              gpr <- min(1.0 - ppr, max(0.0, gpr))
              rem_tx <- setdiff(active_taxa, c(ptx, gtx))
              rem_p <- if (length(rem_tx) > 0L) max(0.0, 1.0 - ppr - gpr) / length(rem_tx) else 0.0
              rf_prob_mat[ri, ] <- rem_p
              rf_prob_mat[ri, ptx] <- ppr
              rf_prob_mat[ri, gtx] <- gpr
            } else {
              rem_p <- (1.0 - ppr) / max(1L, K_tx - 1L)
              rf_prob_mat[ri, ] <- rem_p
              rf_prob_mat[ri, ptx] <- ppr
            }
          }
        }
      }
    }
  } else if (K_tx >= 2L && all(c("rf_predicted_taxon", "rf_predicted_prob") %in% names(df))) {
    if ("spatial_block_cv_acc" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$spatial_block_cv_acc))))) {
      rf_acc_est <- stats::median(suppressWarnings(as.numeric(df$spatial_block_cv_acc)), na.rm = TRUE)
    }
    rf_prob_mat <- matrix(NA_real_, nrow = n, ncol = K_tx, dimnames = list(NULL, active_taxa))
    for (ri in seq_len(n)) {
      ptx <- as.character(df$rf_predicted_taxon[ri]); ppr <- suppressWarnings(as.numeric(df$rf_predicted_prob[ri]))
      if (!is.na(ptx) && (ptx %in% active_taxa) && is.finite(ppr)) {
        ppr <- min(0.999, max(1e-4, ppr))
        gtx <- as.character(grp_vec[ri])
        gpr <- if ("rf_given_taxon_prob" %in% names(df)) suppressWarnings(as.numeric(df$rf_given_taxon_prob[ri])) else NA_real_
        if (!is.na(gtx) && (gtx %in% active_taxa) && gtx != ptx && is.finite(gpr)) {
          gpr <- min(1.0 - ppr, max(0.0, gpr))
          rem_tx <- setdiff(active_taxa, c(ptx, gtx))
          rem_p <- if (length(rem_tx) > 0L) max(0.0, 1.0 - ppr - gpr) / length(rem_tx) else 0.0
          rf_prob_mat[ri, ] <- rem_p
          rf_prob_mat[ri, ptx] <- ppr
          rf_prob_mat[ri, gtx] <- gpr
        } else {
          rf_prob_mat[ri, ] <- (1.0 - ppr) / max(1L, K_tx - 1L)
          rf_prob_mat[ri, ptx] <- ppr
        }
      }
    }
  }

  chance_p <- 1.0 / max(2.0, K_tx)
  if (!is.finite(rf_acc_est)) rf_acc_est <- 0.75
  w_rf_pool <- if (rf_acc_est > chance_p) {
    min(0.85, max(0.35, (rf_acc_est - chance_p) / max(1e-6, 1.0 - chance_p)))
  } else {
    max(0.10, 0.35 * max(0.0, rf_acc_est / max(chance_p, 1e-6)))
  }

  # Evaluate modal concordance per specimen
  conf_flags <- logical(n); triage_cats <- character(n); concordances <- numeric(n); discord_notes <- character(n)
  edaph_best_vec <- character(n); edaph_z_vec <- numeric(n); clim_z_vec <- rep(NA_real_, n); outcrop_spec_flags <- logical(n)

  for (i in seq_len(n)) {
    given_sp <- grp_vec[i]
    morph_sp <- if ("cda_predicted_taxon" %in% names(df) && !is.na(df$cda_predicted_taxon[i])) df$cda_predicted_taxon[i] else given_sp
    if (is_intra) {
      if (!(morph_sp %in% active_taxa)) morph_sp <- given_sp
    }
    cda_prob <- if ("cda_posterior_prob" %in% names(df) && !is.na(df$cda_posterior_prob[i])) as.numeric(df$cda_posterior_prob[i]) else NA_real_
    is_m_out <- if ("morphometric_outlier_flag" %in% names(df) && !is.na(df$morphometric_outlier_flag[i])) isTRUE(as.logical(df$morphometric_outlier_flag[i])) else FALSE

    # Combine Spatial Block-CV OOF RF posterior with decorrelated PCA-env Mahalanobis likelihood
    if (K_tx == 0L) {
      edaph_sp <- given_sp
    } else if (K_tx == 1L) {
      edaph_sp <- active_taxa[1L]
    } else {
      p_maha_i <- pmax(as.numeric(maha_prob_mat[i, ]), 1e-6)
      p_maha_i <- p_maha_i / sum(p_maha_i)
      has_rf_i <- !is.null(rf_prob_mat) && any(is.finite(rf_prob_mat[i, ]))
      if (has_rf_i) {
        p_rf_i <- as.numeric(rf_prob_mat[i, ])
        na_rf <- !is.finite(p_rf_i)
        if (any(na_rf)) {
          rem_mass <- max(0.0, 1.0 - sum(p_rf_i[!na_rf]))
          p_rf_i[na_rf] <- rem_mass / sum(na_rf)
        }
        p_rf_i <- pmax(p_rf_i, 1e-6)
        p_rf_i <- p_rf_i / sum(p_rf_i)
        edaph_scores <- w_rf_pool * log(p_rf_i) + (1.0 - w_rf_pool) * log(p_maha_i)
      } else {
        edaph_scores <- log(p_maha_i)
      }
      edaph_sp <- if (all(!is.finite(edaph_scores))) given_sp else active_taxa[which.max(edaph_scores)][1]
    }
    edaph_best_vec[i] <- edaph_sp
    edaph_z_given <- z_engine$eval_specimen_z(i, given_sp, "edaphic")
    clim_z_given  <- z_engine$eval_specimen_z(i, given_sp, "climatic")
    edaph_z_vec[i] <- round(edaph_z_given, 3)
    clim_z_vec[i]  <- if (is.finite(clim_z_given)) round(clim_z_given, 3) else NA_real_

    alt_parent <- if (morph_sp != given_sp) morph_sp else edaph_sp
    edaph_z_alt <- z_engine$eval_specimen_z(i, alt_parent, "edaphic")

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

    # Agreement matrix across Morphometrics, Fine-Grain Edaphic/Climatic Niche, and Reproductive Phenology
    matches <- c(morph_sp == given_sp, edaph_sp == given_sp, pheno_agree)
    concordance <- sum(matches, na.rm = TRUE) / 3.0
    concordances[i] <- round(concordance, 2)

    # Multi-evidence decision triage
    if (morph_sp != given_sp && edaph_sp != given_sp && (edaph_sp == morph_sp || !pheno_agree || is_m_out || (!is.na(cda_prob) && cda_prob >= 0.70))) {
      conf_flags[i] <- TRUE; triage_cats[i] <- "Severe_Triple_Stream_Conflict"
      discord_notes[i] <- sprintf("Given %s contradicted by Morphology (%s) and Fine-Grain Edaphic (%s)",
                                  given_sp, morph_sp, edaph_sp)
    } else if (edaph_z_given > 3.0 && morph_sp == given_sp) {
      conf_flags[i] <- TRUE; triage_cats[i] <- "Edaphic_Envelope_Violation"
      discord_notes[i] <- sprintf(
        "Specimen outside physiological primary-edaphic tolerance (%s; edaphic_z=%.2f)",
        paste(utils::head(z_engine$primary_edaph_vars, 4L), collapse = ","),
        edaph_z_given
      )
    } else if (!pheno_agree && morph_sp == given_sp && edaph_sp == given_sp) {
      conf_flags[i] <- TRUE; triage_cats[i] <- "Phenological_Anomaly_Outlier"
      discord_notes[i] <- sprintf("Flowering DOY %d diverges from seasonal peak (ang_diff=%.2f rad)", doy[i], ang_diff)
    } else if (morph_sp != given_sp || edaph_sp != given_sp) {
      if ((morph_sp != edaph_sp || (!is.na(cda_prob) && cda_prob >= 0.40 && cda_prob <= 0.65)) &&
          (edaph_z_given <= 2.5 && edaph_z_alt <= 2.5)) {
        conf_flags[i] <- TRUE; triage_cats[i] <- "Putative_Hybrid_Zone_Intergrade"
        discord_notes[i] <- sprintf("Intermediate morphometrics & fine-grain edaphic overlap between %s and %s", given_sp, alt_parent)
      } else if (morph_sp != given_sp && (edaph_sp == morph_sp || is_m_out || (!is.na(cda_prob) && cda_prob >= 0.70))) {
        conf_flags[i] <- TRUE; triage_cats[i] <- "Severe_Triple_Stream_Conflict"
        discord_notes[i] <- sprintf("Cross-modal misidentification signal favoring %s over given %s", alt_parent, given_sp)
      } else {
        conf_flags[i] <- FALSE; triage_cats[i] <- "Clean_MultiModal_Consensus"
        discord_notes[i] <- "Concordant across primary morphological, phenological, and ecological axes"
      }
    } else {
      conf_flags[i] <- FALSE; triage_cats[i] <- "Clean_MultiModal_Consensus"
      discord_notes[i] <- "Concordant across morphological, phenological, and ecological axes"
    }
  }

  df$multimodal_concordance <- concordances
  df$multimodal_conflict_flag <- conf_flags
  df$triage_category <- triage_cats
  df$conflict_discordance_notes <- discord_notes
  df$edaphic_best_fit_taxon <- edaph_best_vec
  df$edaphic_z_score <- edaph_z_vec
  df$climatic_z_score <- clim_z_vec
  df$primary_env_predictors <- rep(paste(z_engine$primary_env_vars, collapse = ";"), n)
  df$primary_clim_predictors <- rep(paste(z_engine$primary_clim_vars, collapse = ";"), n)
  df$primary_edaph_predictors <- rep(paste(z_engine$primary_edaph_vars, collapse = ";"), n)
  df$substrate_specialist_flag <- outcrop_spec_flags
  df$piedmont_outcrop_specialist_flag <- outcrop_spec_flags
  type_eco_ctrl <- compute_type_locality_ecological_control(
    df = df,
    primary_env_vars = z_engine$primary_env_vars,
    primary_edaph_vars = z_engine$primary_edaph_vars,
    rf_var_weights = z_engine$rf_var_weights,
    intraspecific_mode = is_intra
  )
  df$type_edaphic_z_score <- type_eco_ctrl$type_edaphic_z_score
  df$nearest_type_niche_taxon <- type_eco_ctrl$nearest_type_niche_taxon
  df$matches_type_lithology <- type_eco_ctrl$matches_type_lithology
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
  if (n > 2500L) {
    sub_idx <- unique(round(seq(1L, n, length.out = 2500L)))
    coords <- coords[sub_idx, , drop = FALSE]
    X_env <- X_env[sub_idx, , drop = FALSE]
  }
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

compute_spatial_rf_mems <- function(df, env_vars, n_trees = 500, intraspecific_mode = FALSE, vif_res = NULL) {
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
  if (n <= 900L && requireNamespace("adespatial", quietly = TRUE)) {
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
    # Fast O(N^2) double-centering + truncated Lanczos (RSpectra::eigs_sym) for top-5 positive MEM eigenfunctions
    D <- as.matrix(dist(coords))
    k_nn <- min(8L, n)
    thresh <- max(apply(D, 1, function(r) sort(r, partial = k_nn)[k_nn]))
    W <- exp(-D / (thresh + 1e-5)) * (D <= thresh); diag(W) <- 0
    rm_w <- rowMeans(W)
    M <- sweep(sweep(W, 1L, rm_w, "-"), 2L, rm_w, "-") + mean(rm_w)
    eig <- if (n >= 20L && requireNamespace("RSpectra", quietly = TRUE)) {
      tryCatch(RSpectra::eigs_sym(M, k = min(8L, n - 1L), which = "LA"), error = function(e) eigen(M, symmetric = TRUE))
    } else {
      eigen(M, symmetric = TRUE)
    }
    pos_idx <- which(eig$values > 1e-4)
    if (length(pos_idx) == 0L) pos_idx <- 1L
    num_mems <- min(length(pos_idx), 5L)
    MEMs <- eig$vectors[, pos_idx[seq_len(num_mems)], drop = FALSE]
    colnames(MEMs) <- paste0("MEM", seq_len(num_mems))
  }

  # Prepare feature matrices & 3-Block Collinearity Pruning (Dormann et al. 2013: |r| <= 0.70, VIF <= 5.0)
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
  anc_grp_for_vif <- resolve_active_anchor_groups(df, intraspecific_mode = intraspecific_mode)
  col_screen <- if (!is.null(vif_res) && is.list(vif_res) && (!is.null(vif_res$kept_vars) || !is.null(vif_res$selected_variables))) {
    vif_res
  } else if (exists("prune_collinear_predictors", mode = "function")) {
    prune_collinear_predictors(X_env, cor.threshold = 0.70, vif.threshold = 5.0)
  } else if (exists("screen_collinear_predictors_by_block", mode = "function")) {
    screen_collinear_predictors_by_block(X_env, cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = anc_grp_for_vif)
  } else if (exists("screen_collinear_predictors", mode = "function")) {
    screen_collinear_predictors(X_env, cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = anc_grp_for_vif)
  } else {
    ac <- if (requireNamespace("spatialRF", quietly = TRUE)) tryCatch(spatialRF::auto_cor(x = as.data.frame(X_env), cor.threshold = 0.70), error = function(e) NULL) else NULL
    X_s <- if (!is.null(ac) && !is.null(ac$selected.variables)) X_env[, ac$selected.variables, drop = FALSE] else X_env
    av <- if (requireNamespace("spatialRF", quietly = TRUE)) tryCatch(spatialRF::auto_vif(x = as.data.frame(X_s), vif.threshold = 5.0), error = function(e) NULL) else NULL
    s_v <- if (!is.null(av) && !is.null(av$selected.variables)) av$selected.variables else colnames(X_s)
    list(selected_variables = s_v, kept_vars = s_v, pruned_variables = setdiff(colnames(X_env), s_v))
  }
  kept_v <- if (!is.null(col_screen$kept_vars)) col_screen$kept_vars else col_screen$selected_variables
  if (length(intersect(colnames(X_env), kept_v)) >= 1L) {
    env_vars <- intersect(colnames(X_env), unique(kept_v))
    X_env <- X_env[, env_vars, drop = FALSE]
  }
  X_all <- cbind(X_env, MEMs)
  y <- factor(df$analysis_group)

  # Authority-Anchored Training Cohort Selection:
  # Strictly Tier_1_Gold + Tier_2_Silver + promoted Bronze anchors with coordinateUncertaintyInMeters <= 10000
  elig_coord <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, n)
  elig_coord[is.na(elig_coord)] <- TRUE
  if ("coordinateUncertaintyInMeters" %in% names(df)) {
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    elig_coord <- elig_coord & (is.na(u_m) | u_m <= 10000.0)
  }
  det_tier <- if (exists(".standardize_determiner_tier_vec", mode = "function")) {
    .standardize_determiner_tier_vec(df)
  } else if ("determiner_tier" %in% names(df)) {
    as.character(df$determiner_tier)
  } else {
    rep("Tier_2_Silver", nrow(df))
  }
  type_info_rf <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    .identify_type_specimen_vec(df)
  } else {
    list(is_type_specimen = rep(FALSE, nrow(df)))
  }
  is_promoted <- rep(FALSE, n)
  if ("is_promoted_anchor" %in% names(df)) { pa <- as.logical(df$is_promoted_anchor); pa[is.na(pa)] <- FALSE; is_promoted <- is_promoted | pa }
  if ("is_promoted_bronze" %in% names(df)) { pb <- as.logical(df$is_promoted_bronze); pb[is.na(pb)] <- FALSE; is_promoted <- is_promoted | pb }
  is_gold_silver <- (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted | type_info_rf$is_type_specimen
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
  case_wts[anchor_idx] <- pmin(4.0, pmax(0.25, case_wts[anchor_idx]))
  case_wts[anchor_idx] <- case_wts[anchor_idx] / mean(case_wts[anchor_idx], na.rm = TRUE)

  X_env_anc <- X_env[anchor_idx, , drop = FALSE]
  X_all_anc <- X_all[anchor_idx, , drop = FALSE]
  wts_anc   <- case_wts[anchor_idx]

  rf_threads <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(32L) else 1L
  if (requireNamespace("ranger", quietly = TRUE)) {
    df_env_anc <- data.frame(y = y_anchor, as.data.frame(X_env_anc), check.names = FALSE)
    df_all_anc <- data.frame(y = y_anchor, as.data.frame(X_all_anc), check.names = FALSE)
    rf_nonspatial <- ranger::ranger(y ~ ., data = df_env_anc, num.trees = n_trees, importance = "permutation", probability = FALSE, case.weights = wts_anc, num.threads = rf_threads, seed = 42)
    oob_acc_nonspatial <- 1 - rf_nonspatial$prediction.error
    rf_spatial <- ranger::ranger(y ~ ., data = df_all_anc, num.trees = n_trees, importance = "permutation", probability = FALSE, case.weights = wts_anc, num.threads = rf_threads, seed = 42)
    oob_acc_spatial <- 1 - rf_spatial$prediction.error
    rf_prob_anchor <- ranger::ranger(y ~ ., data = df_env_anc, num.trees = n_trees, probability = TRUE, case.weights = wts_anc, num.threads = rf_threads, seed = 42)
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
                                 probability = TRUE, case.weights = wts_anc[tr_sub], num.threads = rf_threads, seed = 42 + as.integer(f_id))
        pr_mat <- stats::predict(fit_cv, data = te_df, num.threads = rf_threads)$predictions
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
      pass_pr <- stats::predict(rf_prob_anchor, data = pass_df, num.threads = rf_threads)$predictions
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
    oof_prob_matrix = full_prob_mat,
    full_prob_mat = full_prob_mat,
    block_cv_acc = spatial_block_cv_acc,
    predicted_taxon = full_preds,
    predicted_prob = round(rf_pred_prob, 4),
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
.prepare_partition_scores_glob <- function(bg_mat, rf_var_weights = NULL) {
  if (is.null(bg_mat) || ncol(bg_mat) == 0 || nrow(bg_mat) < 6) return(NULL)
  for (j in seq_len(ncol(bg_mat))) {
    col_v <- bg_mat[, j]
    if (any(!is.finite(col_v))) {
      med_j <- stats::median(col_v[is.finite(col_v)], na.rm = TRUE)
      if (!is.finite(med_j)) med_j <- 0.0
      bg_mat[!is.finite(col_v), j] <- med_j
    }
  }
  active_cols <- which(apply(bg_mat, 2, stats::sd, na.rm = TRUE) > 1e-8)
  if (length(active_cols) == 0) return(NULL)
  bg_sub <- bg_mat[, active_cols, drop = FALSE]
  use_rf_w <- !is.null(rf_var_weights) && !is.null(colnames(bg_sub)) && any(colnames(bg_sub) %in% names(rf_var_weights))
  if (use_rf_w) {
    w_vec <- rep(1.0 / ncol(bg_sub), ncol(bg_sub))
    names(w_vec) <- colnames(bg_sub)
    ov_nm <- intersect(colnames(bg_sub), names(rf_var_weights))
    w_raw <- as.numeric(rf_var_weights[ov_nm])
    w_raw[!is.finite(w_raw) | w_raw <= 0] <- 1e-6
    w_vec[ov_nm] <- w_raw
    w_vec <- (w_vec / sum(w_vec)) * ncol(bg_sub)
    bg_z <- scale(bg_sub, center = TRUE, scale = TRUE)
    bg_z <- sweep(bg_z, 2L, sqrt(w_vec), `*`)
  } else {
    bg_z <- NULL
  }
  scores_glob <- tryCatch({
    if (ncol(bg_sub) >= 2) {
      if (!is.null(bg_z)) {
        as.data.frame(stats::prcomp(bg_z, center = FALSE, scale. = FALSE)$x[, 1:2, drop = FALSE])
      } else {
        as.data.frame(stats::prcomp(bg_sub, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE])
      }
    } else {
      z1 <- as.numeric(scale(bg_sub[, 1]))
      data.frame(PC1 = z1, PC2 = 0.05 * seq_along(z1) / length(z1))
    }
  }, error = function(e) NULL)
  list(bg_sub = bg_sub, scores_glob = scores_glob)
}

.run_ecospat_test_chunk <- function(z1, z2, rep_n, test_kind) {
  if (identical(test_kind, "equiv")) {
    eq <- tryCatch(
      suppressMessages(suppressWarnings(ecospat::ecospat.niche.equivalency.test(z1, z2, rep = rep_n, overlap.alternative = "lower", ncores = 1))),
      error = function(e) suppressMessages(suppressWarnings(ecospat::ecospat.niche.equivalency.test(z1, z2, rep = rep_n, ncores = 1)))
    )
    list(D = as.numeric(eq$sim$D), I = as.numeric(eq$sim$I))
  } else {
    bg_sim <- tryCatch(
      suppressMessages(suppressWarnings(ecospat::ecospat.niche.similarity.test(z1, z2, rep = rep_n, overlap.alternative = "lower", rand.type = 2, ncores = 1))),
      error = function(e) tryCatch(
        suppressMessages(suppressWarnings(ecospat::ecospat.niche.similarity.test(z1, z2, rep = rep_n, rand.type = 2, ncores = 1))),
        error = function(e2) NULL
      )
    )
    if (!is.null(bg_sim) && !is.null(bg_sim$sim$D)) {
      list(D = as.numeric(bg_sim$sim$D), I = as.numeric(bg_sim$sim$I))
    } else {
      list(D = numeric(0), I = numeric(0))
    }
  }
}

.evaluate_partition_jobs_parallel <- function(job_specs, grp_col, n_perm) {
  n_jobs <- length(job_specs)
  if (n_jobs == 0L) return(list())
  can_par <- (n_perm >= 99L) && exists("packera_mclapply", mode = "function") &&
             exists("resolve_packera_workers", mode = "function") && (resolve_packera_workers(32L) > 1L)
  if (!can_par) {
    return(lapply(job_specs, function(jb) {
      compute_pca_env_partition_overlap(
        jb$bg_mat, grp_col, jb$sp1, jb$sp2,
        n_perm = n_perm, rf_var_weights = jb$rf_var_weights, precomputed_prep = jb$prep
      )
    }))
  }

  # Build z1, z2, ov grids once per (pair, partition) job, then fan out equivalency + similarity chunks across 32 cores
  grid_states <- vector("list", n_jobs)
  flat_tasks <- list()
  max_w <- resolve_packera_workers(32L)
  chunks_per_test <- max(1L, min(6L, floor(max_w / max(2L * n_jobs, 1L))))

  for (j_i in seq_len(n_jobs)) {
    jb <- job_specs[[j_i]]
    prep <- if (!is.null(jb$prep)) jb$prep else .prepare_partition_scores_glob(jb$bg_mat, jb$rf_var_weights)
    idx1 <- which(grp_col == jb$sp1)
    idx2 <- which(grp_col == jb$sp2)
    g_st <- NULL
    if (!is.null(prep) && !is.null(prep$scores_glob) && requireNamespace("ecospat", quietly = TRUE) &&
        length(idx1) >= 5L && length(idx2) >= 5L) {
      g_st <- tryCatch({
        sc <- prep$scores_glob
        z1 <- ecospat::ecospat.grid.clim.dyn(glob = sc, glob1 = sc, sp = sc[idx1, 1:2, drop = FALSE], R = 100)
        z2 <- ecospat::ecospat.grid.clim.dyn(glob = sc, glob1 = sc, sp = sc[idx2, 1:2, drop = FALSE], R = 100)
        ov <- ecospat::ecospat.niche.overlap(z1, z2, cor = TRUE)
        if (is.finite(as.numeric(ov$D)) && is.finite(as.numeric(ov$I))) {
          list(z1 = z1, z2 = z2, D_emp = as.numeric(ov$D), I_emp = as.numeric(ov$I))
        } else NULL
      }, error = function(e) NULL)
    }
    grid_states[[j_i]] <- list(prep = prep, g_st = g_st)
    if (!is.null(g_st)) {
      base_rep <- floor(n_perm / chunks_per_test)
      rem_rep <- n_perm - base_rep * chunks_per_test
      rep_vec <- rep(base_rep, chunks_per_test)
      if (rem_rep > 0L) rep_vec[seq_len(rem_rep)] <- rep_vec[seq_len(rem_rep)] + 1L
      rep_vec <- rep_vec[rep_vec > 0L]
      for (tk in c("equiv", "bg")) {
        for (rv in rep_vec) {
          flat_tasks[[length(flat_tasks) + 1L]] <- list(job_idx = j_i, kind = tk, rep_n = as.integer(rv))
        }
      }
    }
  }

  chunk_outs <- list()
  if (length(flat_tasks) > 0L) {
    base_seed <- sample.int(100000000L, 1L)
    n_workers <- resolve_packera_workers(32L, length(flat_tasks))
    chunk_outs <- packera_mclapply(
      seq_along(flat_tasks),
      function(t_i) {
        set.seed(as.integer((base_seed + t_i * 104729L) %% 2147483647L))
        tk_obj <- flat_tasks[[t_i]]
        gs <- grid_states[[tk_obj$job_idx]]$g_st
        tryCatch(.run_ecospat_test_chunk(gs$z1, gs$z2, tk_obj$rep_n, tk_obj$kind), error = function(e) NULL)
      },
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )
  }

  out_res <- vector("list", n_jobs)
  for (j_i in seq_len(n_jobs)) {
    jb <- job_specs[[j_i]]
    gs <- grid_states[[j_i]]$g_st
    prep <- grid_states[[j_i]]$prep
    if (!is.null(gs) && length(chunk_outs) == length(flat_tasks)) {
      eq_d <- numeric(0); eq_i <- numeric(0); bg_d <- numeric(0); bg_i <- numeric(0)
      ok_chunks <- TRUE
      for (t_i in seq_along(flat_tasks)) {
        if (flat_tasks[[t_i]]$job_idx == j_i) {
          co <- chunk_outs[[t_i]]
          if (is.null(co) || (identical(flat_tasks[[t_i]]$kind, "equiv") && length(co$D) == 0L)) {
            ok_chunks <- FALSE
            break
          }
          if (identical(flat_tasks[[t_i]]$kind, "equiv")) {
            eq_d <- c(eq_d, co$D); eq_i <- c(eq_i, co$I)
          } else {
            bg_d <- c(bg_d, co$D); bg_i <- c(bg_i, co$I)
          }
        }
      }
      if (ok_chunks && length(eq_d) == n_perm) {
        p_bg_d <- if (length(bg_d) > 0L) (sum(bg_d <= gs$D_emp, na.rm = TRUE) + 1) / (length(bg_d) + 1) else NA_real_
        p_bg_i <- if (length(bg_i) > 0L) (sum(bg_i <= gs$I_emp, na.rm = TRUE) + 1) / (length(bg_i) + 1) else NA_real_
        if (is.na(p_bg_d) || !is.finite(p_bg_d)) {
          # Run fast analytical background fallback if ecospat similarity returned empty
          fb <- compute_pca_env_partition_overlap(
            jb$bg_mat, grp_col, jb$sp1, jb$sp2,
            n_perm = n_perm, rf_var_weights = jb$rf_var_weights, precomputed_prep = prep
          )
          out_res[[j_i]] <- fb
        } else {
          out_res[[j_i]] <- list(
            D_emp = gs$D_emp,
            I_emp = gs$I_emp,
            p_equiv_D = (sum(eq_d <= gs$D_emp, na.rm = TRUE) + 1) / (length(eq_d) + 1),
            p_equiv_I = (sum(eq_i <= gs$I_emp, na.rm = TRUE) + 1) / (length(eq_i) + 1),
            p_bg_D = p_bg_d,
            p_bg_I = p_bg_i,
            D_null = eq_d,
            I_null = eq_i
          )
        }
        next
      }
    }
    out_res[[j_i]] <- compute_pca_env_partition_overlap(
      jb$bg_mat, grp_col, jb$sp1, jb$sp2,
      n_perm = n_perm, rf_var_weights = jb$rf_var_weights, precomputed_prep = prep
    )
  }
  out_res
}

compute_pca_env_partition_overlap <- function(bg_mat, grp_col, sp1, sp2, n_perm = 100, rf_var_weights = NULL, precomputed_prep = NULL) {
  prep <- if (!is.null(precomputed_prep)) precomputed_prep else .prepare_partition_scores_glob(bg_mat, rf_var_weights)
  if (is.null(prep)) {
    return(list(D_emp = NA_real_, I_emp = NA_real_, p_equiv_D = NA_real_, p_equiv_I = NA_real_,
                p_bg_D = NA_real_, p_bg_I = NA_real_, D_null = rep(NA_real_, n_perm), I_null = rep(NA_real_, n_perm)))
  }
  bg_sub <- prep$bg_sub
  scores_glob <- prep$scores_glob

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

compute_edaphic_schoener_d <- function(
  df,
  primary_edaph_vars = NULL,
  rf_var_weights = NULL,
  n_perm = 100,
  intraspecific_mode = FALSE
) {
  if (!("analysis_group" %in% names(df)) || isTRUE(intraspecific_mode)) {
    df <- resolve_analysis_group(
      df,
      intraspecific_mode = intraspecific_mode,
      primary_env_vars = primary_edaph_vars,
      rf_var_weights = rf_var_weights
    )
  }
  elig_mask <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, nrow(df))
  elig_mask[is.na(elig_mask)] <- TRUE
  if ("coordinateUncertaintyInMeters" %in% names(df)) {
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    elig_mask <- elig_mask & (is.na(u_m) | u_m <= 10000.0)
  }
  df_calib <- if (sum(elig_mask) >= 10L) df[elig_mask, , drop = FALSE] else df
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df_calib$species_standardized))) < 2L)
  grp_col <- if (is_intra) df_calib$analysis_group else df_calib$species_standardized
  cand_groups <- if (is_intra) sort(unique(stats::na.omit(grp_col))) else get_active_packera_taxa(grp_col)
  if (length(cand_groups) < 2L) cand_groups <- sort(unique(stats::na.omit(df_calib$analysis_group)))
  pairs <- if (length(cand_groups) >= 2L) utils::combn(cand_groups, 2L, simplify = FALSE) else list()

  edaph_cols <- intersect(as.character(primary_edaph_vars), names(df_calib))
  if (length(edaph_cols) == 0L) {
    cand_all <- discover_candidate_eco_predictors(df_calib, include_pheno = FALSE, min_finite = 3L)
    edaph_cols <- grep("^(fine_|ssurgo_|polaris_|solus_|soil_|bdricm|resdept_r|cfvo|awc|alr_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)", cand_all, value = TRUE, ignore.case = TRUE)
  }
  if (length(edaph_cols) == 0L || length(pairs) == 0L) {
    return(data.frame(
      Pair = character(0), Taxon1 = character(0), Taxon2 = character(0),
      D_edaph = numeric(0), I_edaph = numeric(0),
      p_equiv_D_edaph = numeric(0), p_bg_D_edaph = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  X_ed <- as.matrix(df_calib[, edaph_cols, drop = FALSE])
  prep_ed <- .prepare_partition_scores_glob(X_ed, rf_var_weights)
  ed_jobs <- list()
  ed_meta <- list()
  for (pr in pairs) {
    sp1 <- pr[1]; sp2 <- pr[2]
    if (sum(grp_col == sp1, na.rm = TRUE) < 5L || sum(grp_col == sp2, na.rm = TRUE) < 5L) next
    pair_name <- paste0(gsub("Packera ", "P.", sp1), " vs ", gsub("Packera ", "P.", sp2))
    ed_jobs[[length(ed_jobs) + 1L]] <- list(bg_mat = X_ed, sp1 = sp1, sp2 = sp2, rf_var_weights = rf_var_weights, prep = prep_ed)
    ed_meta[[length(ed_meta) + 1L]] <- list(pair_name = pair_name, sp1 = sp1, sp2 = sp2)
  }
  if (length(ed_jobs) == 0L) {
    return(data.frame(
      Pair = character(0), Taxon1 = character(0), Taxon2 = character(0),
      D_edaph = numeric(0), I_edaph = numeric(0),
      p_equiv_D_edaph = numeric(0), p_bg_D_edaph = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  ed_outs <- .evaluate_partition_jobs_parallel(ed_jobs, grp_col, n_perm)
  out_list <- vector("list", length(ed_jobs))
  for (k in seq_along(ed_jobs)) {
    ov <- ed_outs[[k]]
    md <- ed_meta[[k]]
    out_list[[k]] <- data.frame(
      Pair = md$pair_name, Taxon1 = md$sp1, Taxon2 = md$sp2,
      D_edaph = round(ov$D_emp, 4), I_edaph = round(ov$I_emp, 4),
      p_equiv_D_edaph = round(ov$p_equiv_D, 4), p_bg_D_edaph = round(ov$p_bg_D, 4),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out_list)
}

run_warren_niche_identity_tests <- function(
  df,
  env_vars = NULL,
  n_perm = 100,
  intraspecific_mode = FALSE,
  primary_env_vars = NULL,
  primary_clim_vars = NULL,
  primary_edaph_vars = NULL,
  rf_var_weights = NULL
) {
  message(sprintf("Executing Partitioned PCA-env Niche Equivalency & Background Similarity Tests (%d permutations per pair)...", n_perm))
  if (!("analysis_group" %in% names(df)) || isTRUE(intraspecific_mode)) {
    df <- resolve_analysis_group(
      df,
      intraspecific_mode = intraspecific_mode,
      primary_env_vars = primary_env_vars,
      rf_var_weights = rf_var_weights
    )
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

  if (is.null(env_vars) || length(env_vars) == 0L) {
    env_vars <- discover_candidate_eco_predictors(df_calib, include_pheno = TRUE, min_finite = 3L)
  }
  all_req_vars <- unique(c(env_vars, primary_env_vars, primary_clim_vars, primary_edaph_vars))
  env_vars <- intersect(all_req_vars, names(df_calib))
  X_raw <- as.matrix(df_calib[, env_vars, drop = FALSE])

  has_primary_override <- (!is.null(primary_env_vars) && length(intersect(primary_env_vars, colnames(X_raw))) > 0L) ||
    (!is.null(primary_clim_vars) && length(intersect(primary_clim_vars, colnames(X_raw))) > 0L) ||
    (!is.null(primary_edaph_vars) && length(intersect(primary_edaph_vars, colnames(X_raw))) > 0L)

  if (has_primary_override) {
    clim_vars <- intersect(as.character(primary_clim_vars), colnames(X_raw))
    edaph_vars <- intersect(as.character(primary_edaph_vars), colnames(X_raw))
    joint_vars <- intersect(as.character(primary_env_vars), colnames(X_raw))
    if (length(clim_vars) == 0L) clim_vars <- grep("^(bio[0-9]+|delta_doy|doy)", joint_vars, value = TRUE, ignore.case = TRUE)
    if (length(edaph_vars) == 0L) edaph_vars <- setdiff(joint_vars, clim_vars)
    if (length(joint_vars) == 0L) joint_vars <- unique(c(clim_vars, edaph_vars))
    block_screen <- list(
      climatic_variables = clim_vars,
      edaphic_variables = edaph_vars,
      selected_variables = joint_vars,
      kept_vars = joint_vars
    )
  } else {
    anc_grp <- resolve_active_anchor_groups(df_calib, intraspecific_mode = is_intra)
    block_screen <- if (exists("screen_collinear_predictors_by_block", mode = "function")) {
      screen_collinear_predictors_by_block(X_raw, cor.threshold = 0.70, vif.threshold = 5.0, anchor_groups = anc_grp)
    } else {
      clim_c <- grep("^(bio|delta_doy|doy)", env_vars, value = TRUE, ignore.case = TRUE)
      edaph_c <- setdiff(env_vars, clim_c)
      list(climatic_variables = clim_c, edaphic_variables = edaph_c, selected_variables = env_vars, kept_vars = env_vars)
    }
    clim_vars <- intersect(block_screen$climatic_variables, colnames(X_raw))
    if (length(clim_vars) == 0) clim_vars <- intersect(grep("^(bio|delta_doy|doy)", colnames(X_raw), value = TRUE, ignore.case = TRUE), colnames(X_raw))
    edaph_vars <- intersect(unique(c(block_screen$edaphic_variables, block_screen$topo_variables)), colnames(X_raw))
    if (length(edaph_vars) == 0) edaph_vars <- setdiff(colnames(X_raw), clim_vars)
    joint_vars <- unique(c(clim_vars, edaph_vars))
  }

  X_clim_part <- if (length(clim_vars) > 0) X_raw[, clim_vars, drop = FALSE] else X_raw
  X_edaph_part <- if (length(edaph_vars) > 0) X_raw[, edaph_vars, drop = FALSE] else X_raw
  X_joint_part <- X_raw[, joint_vars, drop = FALSE]
  prep_clim <- .prepare_partition_scores_glob(X_clim_part, rf_var_weights)
  prep_edaph <- .prepare_partition_scores_glob(X_edaph_part, rf_var_weights)
  prep_joint <- .prepare_partition_scores_glob(X_joint_part, rf_var_weights)

  valid_pair_meta <- list()
  part_jobs <- list()
  for (pr in pairs) {
    sp1 <- pr[1]; sp2 <- pr[2]; pair_name <- paste0(gsub("Packera ", "P.", sp1), " vs ", gsub("Packera ", "P.", sp2))
    n1 <- sum(grp_col == sp1, na.rm = TRUE)
    n2 <- sum(grp_col == sp2, na.rm = TRUE)
    if (n1 < 5 || n2 < 5) next
    valid_pair_meta[[length(valid_pair_meta) + 1L]] <- list(sp1 = sp1, sp2 = sp2, pair_name = pair_name, n1 = n1, n2 = n2)
    part_jobs[[length(part_jobs) + 1L]] <- list(bg_mat = X_clim_part, sp1 = sp1, sp2 = sp2, rf_var_weights = rf_var_weights, prep = prep_clim)
    part_jobs[[length(part_jobs) + 1L]] <- list(bg_mat = X_edaph_part, sp1 = sp1, sp2 = sp2, rf_var_weights = rf_var_weights, prep = prep_edaph)
    part_jobs[[length(part_jobs) + 1L]] <- list(bg_mat = X_joint_part, sp1 = sp1, sp2 = sp2, rf_var_weights = rf_var_weights, prep = prep_joint)
  }

  part_outs <- .evaluate_partition_jobs_parallel(part_jobs, grp_col, n_perm)

  for (p_i in seq_along(valid_pair_meta)) {
    md <- valid_pair_meta[[p_i]]
    sp1 <- md$sp1; sp2 <- md$sp2; pair_name <- md$pair_name; n1 <- md$n1; n2 <- md$n2
    clim_res  <- part_outs[[(p_i - 1L) * 3L + 1L]]
    edaph_res <- part_outs[[(p_i - 1L) * 3L + 2L]]
    joint_res <- part_outs[[(p_i - 1L) * 3L + 3L]]

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

# Source modular Figure 05, Spatial RF plotting (fill = NA), and MaxEnt/SRF SDM projection helpers
source(file.path(step05_dir, "step05_sdm_maxent.R"))


# ------------------------------------------------------------------------------
# 7. Main Execution Pipeline
# ------------------------------------------------------------------------------
left_join_by_voucher_individual <- function(base_df, mod_df, keep_cols) {
  if (is.null(mod_df) || nrow(mod_df) == 0L || !("catalogNumber" %in% names(mod_df))) return(base_df)
  base_df$catalogNumber <- as.character(base_df$catalogNumber)
  mod_df$catalogNumber <- as.character(mod_df$catalogNumber)
  k_cols <- intersect(names(mod_df), keep_cols)
  val_cols <- setdiff(k_cols, c("catalogNumber", "plant_individual_id"))
  overlap_cols <- intersect(val_cols, names(base_df))
  add_cols <- setdiff(val_cols, overlap_cols)

  has_mod_pid <- "plant_individual_id" %in% k_cols
  if (has_mod_pid) {
    mpid <- ifelse(is.na(mod_df$plant_individual_id), 0L, suppressWarnings(as.integer(mod_df$plant_individual_id)))
    mpid[is.na(mpid)] <- 0L
    mod_df$plant_individual_id <- mpid
    mod_sub <- mod_df[!duplicated(paste0(mod_df$catalogNumber, "__", mod_df$plant_individual_id)), k_cols, drop = FALSE]
  } else {
    mod_sub <- mod_df[!duplicated(mod_df$catalogNumber), k_cols, drop = FALSE]
  }

  has_base_pid <- "plant_individual_id" %in% names(base_df)
  if (has_base_pid && has_mod_pid) {
    bpid <- ifelse(is.na(base_df$plant_individual_id), 0L, suppressWarnings(as.integer(base_df$plant_individual_id)))
    bpid[is.na(bpid)] <- 0L
    base_df$plant_individual_id <- bpid

    b_counts <- table(base_df$catalogNumber)
    single_zero_cats <- unique(base_df$catalogNumber[base_df$plant_individual_id == 0L & as.integer(b_counts[base_df$catalogNumber]) == 1L])
    mod_nz <- mod_sub[mod_sub$catalogNumber %in% single_zero_cats & mod_sub$plant_individual_id != 0L, c("catalogNumber", "plant_individual_id"), drop = FALSE]
    if (nrow(mod_nz) > 0L) {
      keep_b <- !(base_df$catalogNumber %in% unique(mod_nz$catalogNumber))
      exp_b <- dplyr::left_join(base_df[!keep_b, setdiff(names(base_df), "plant_individual_id"), drop = FALSE], mod_nz, by = "catalogNumber")
      base_df <- dplyr::bind_rows(base_df[keep_b, , drop = FALSE], exp_b[, names(base_df), drop = FALSE])
    }

    join_sub <- mod_sub[, unique(c("catalogNumber", "plant_individual_id", add_cols)), drop = FALSE]
    out_df <- dplyr::left_join(base_df, join_sub, by = c("catalogNumber", "plant_individual_id"))

    m_counts <- table(mod_sub$catalogNumber)
    sheet_lvl_sub <- mod_sub[mod_sub$plant_individual_id == 0L & as.integer(m_counts[mod_sub$catalogNumber]) == 1L, , drop = FALSE]
    idx_exact <- match(paste0(out_df$catalogNumber, "__", out_df$plant_individual_id), paste0(mod_sub$catalogNumber, "__", mod_sub$plant_individual_id))
    idx_sheet <- if (nrow(sheet_lvl_sub) > 0L) match(out_df$catalogNumber, sheet_lvl_sub$catalogNumber) else rep(NA_integer_, nrow(out_df))
    b_cnt2 <- table(out_df$catalogNumber)
    s_cats <- intersect(names(b_cnt2[b_cnt2 == 1L]), names(m_counts[m_counts == 1L]))
    idx_single_cat <- ifelse(out_df$catalogNumber %in% s_cats, match(out_df$catalogNumber, mod_sub$catalogNumber), NA_integer_)

    for (cl in add_cols) {
      na_m <- is.na(out_df[[cl]]) & is.na(idx_exact)
      if (any(na_m & !is.na(idx_sheet))) { hit <- na_m & !is.na(idx_sheet); out_df[[cl]][hit] <- sheet_lvl_sub[[cl]][idx_sheet[hit]] }
      na_m2 <- is.na(out_df[[cl]]) & is.na(idx_exact) & !is.na(idx_single_cat)
      if (any(na_m2)) out_df[[cl]][na_m2] <- mod_sub[[cl]][idx_single_cat[na_m2]]
    }
    for (cl in overlap_cols) {
      na_o <- is.na(out_df[[cl]])
      if (any(na_o & !is.na(idx_exact))) { h1 <- na_o & !is.na(idx_exact); out_df[[cl]][h1] <- mod_sub[[cl]][idx_exact[h1]] }
      na_o2 <- is.na(out_df[[cl]])
      if (any(na_o2 & !is.na(idx_sheet))) { h2 <- na_o2 & !is.na(idx_sheet); out_df[[cl]][h2] <- sheet_lvl_sub[[cl]][idx_sheet[h2]] }
      na_o3 <- is.na(out_df[[cl]]) & !is.na(idx_single_cat)
      if (any(na_o3)) out_df[[cl]][na_o3] <- mod_sub[[cl]][idx_single_cat[na_o3]]
    }
    return(out_df)
  } else {
    join_cols <- unique(c("catalogNumber", if (has_mod_pid && !has_base_pid) "plant_individual_id" else character(0), add_cols))
    out_df <- dplyr::left_join(base_df, mod_sub[, join_cols, drop = FALSE], by = "catalogNumber")
    if (length(overlap_cols) > 0L) {
      idx_c <- if (has_mod_pid && "plant_individual_id" %in% names(out_df)) {
        match(paste0(out_df$catalogNumber, "__", out_df$plant_individual_id), paste0(mod_sub$catalogNumber, "__", mod_sub$plant_individual_id))
      } else match(out_df$catalogNumber, mod_sub$catalogNumber)
      for (cl in overlap_cols) {
        na_o <- is.na(out_df[[cl]]) & !is.na(idx_c)
        if (any(na_o)) out_df[[cl]][na_o] <- mod_sub[[cl]][idx_c[na_o]]
      }
    }
    return(out_df)
  }
}

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
  sql_master <- if (!isTRUE(opts$demo) && exists("try_fetch_master_view_for_cohort", mode = "function")) {
    tryCatch(try_fetch_master_view_for_cohort(vouchers_df = vouchers_df, vouchers_path = opts$vouchers, demo = isTRUE(opts$demo), require_step04 = TRUE), error = function(e) NULL)
  } else NULL
  if (!is.null(sql_master) && nrow(sql_master) > 0L) {
    vouchers_df <- sql_master
    vouchers_df$species_standardized <- standardize_packera_taxon(ifelse(!is.na(vouchers_df$species_standardized) & nzchar(vouchers_df$species_standardized), vouchers_df$species_standardized, vouchers_df$scientificName))
    message(sprintf("Hydrated %d specimen rows across Steps 01-04 directly from SQLite v_integrative_triage_master.", nrow(vouchers_df)))
  } else {
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
    if ("cda_predicted" %in% names(morph_df) && !("cda_predicted_taxon" %in% names(morph_df))) {
      morph_df$cda_predicted_taxon <- morph_df$cda_predicted
    }
    if ("is_outlier" %in% names(morph_df) && !("morphometric_outlier_flag" %in% names(morph_df))) {
      morph_df$morphometric_outlier_flag <- as.logical(morph_df$is_outlier)
    }
    if ("morphometric_outlier_flag" %in% names(morph_df) && !("misidentification_flag" %in% names(morph_df))) {
      morph_df$misidentification_flag <- morph_df$morphometric_outlier_flag
    } else if ("misidentification_flag" %in% names(morph_df) && !("morphometric_outlier_flag" %in% names(morph_df))) {
      morph_df$morphometric_outlier_flag <- morph_df$misidentification_flag
    }
    if ("plant_height_mm" %in% names(morph_df) && !("plant_height_cm" %in% names(morph_df))) {
      morph_df$plant_height_cm <- round(suppressWarnings(as.numeric(morph_df$plant_height_mm)) / 10.0, 2)
    }
    if ("is_promoted_anchor" %in% names(morph_df) && !("is_promoted_bronze" %in% names(morph_df))) {
      morph_df$is_promoted_bronze <- as.logical(morph_df$is_promoted_anchor)
    }
    dyn_pc_cols <- grep("^(Sym_|Asym_)?PC[0-9]+$", names(morph_df), value = TRUE)
    m_cols <- intersect(names(morph_df), unique(c(
      "catalogNumber", "plant_individual_id", "morph_tier_source", "n_tier1_leaves", "n_tier2_leaves",
      "gmm_cluster", "gmm_uncertainty", "is_novel_morphotype",
      "cda_predicted_taxon", "cda_predicted", "cda_posterior_prob",
      "cda_second_taxon", "cda_second_posterior_prob", "cda_posterior_margin",
      "loocv_accuracy", "mahalanobis_d2", "mahalanobis_pval",
      "is_outlier", "morphometric_outlier_flag", "misidentification_flag",
      "can1", "can2", "can3", dyn_pc_cols,
      "efd_harmonic_opt", "efd_k95_pcs", "selected_morph_features",
      "solidity", "q25_solidity", "dissection_index", "dissected_leaf_fraction", "basal_leaf_count", "heteroblasty_index", "is_dissected",
      "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
      "plant_height_mm", "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem",
      "capitula_count", "involucre_height_mm", "involucre_width_mm", "AR_inv",
      "ecoregion_cda_predicted", "ecoregion_cda_posterior_prob",
      "type_status", "is_type_specimen", "has_empirical_type_anchor",
      "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance", "is_typical_morphotype",
      "is_promoted_bronze", "is_promoted_anchor", "authority_weight"
    )))
    vouchers_df <- left_join_by_voucher_individual(vouchers_df, morph_df, m_cols)
  }

  # Ingest Voucher Phenological States (Computer Vision Anthesis Ingestion via left_join)
  if (file.exists(opts$pheno_states)) {
    message(sprintf("Ingesting voucher phenological states from %s...", opts$pheno_states))
    pheno_states_df <- readr::read_csv(opts$pheno_states, show_col_types = FALSE)
    pheno_states_df$catalogNumber <- as.character(pheno_states_df$catalogNumber)
    vouchers_df$catalogNumber <- as.character(vouchers_df$catalogNumber)
    p_cols <- intersect(names(pheno_states_df), c("catalogNumber", "total_capitula", "n_anthesis", "n_fruit", "n_fruiting", "n_bud", "n_senescent", "voucher_phenological_state", "reproductive_stage_confidence", "dominant_pappus_ratio"))
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
  }

  # Filter strictly for primary duplicate vouchers (exsiccatae deduplication)
  if ("is_primary_duplicate" %in% names(vouchers_df)) {
    n_pre_filter <- nrow(vouchers_df)
    vouchers_df <- vouchers_df[is.na(vouchers_df$is_primary_duplicate) | toupper(as.character(vouchers_df$is_primary_duplicate)) %in% c("TRUE", "T", "1"), ]
    n_post_filter <- nrow(vouchers_df)
    message(sprintf("Exsiccatae Stratification Filter: Retaining %d primary vouchers (%d duplicate sheets excluded to prevent spatial pseudoreplication).",
                    n_post_filter, n_pre_filter - n_post_filter))
  }

  # Filter valid coordinates
  vouchers_df <- vouchers_df[!is.na(vouchers_df$latitude) & !is.na(vouchers_df$longitude), ]
  message(sprintf("Processing %d georeferenced vouchers across target taxa...", nrow(vouchers_df)))

  demo_mode <- isTRUE(opts$demo)

  # Stage 1: Environmental Layer Extraction (250m SoilGrids & all 19 WorldClim v2.1 bio1..bio19)
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

  # Stage 1c: Cohort-Adaptive Latitudinal Baseline & Circular Anthesis Phenology Engine
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

  # Dynamically discover all 19 bio1..bio19, fine_*, ssurgo_*, soil_*, topo_*, and pheno predictors
  cand_eco <- discover_candidate_eco_predictors(vouchers_df, include_pheno = TRUE, min_finite = 3L)
  for (ev in cand_eco) {
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
  env_vars <- cand_eco

  # Stage 2: Empirical 3-Block Collinearity Pruning (|r| < 0.70, VIF < 5.0) & Spatial RF Importance Gating
  is_intra_run <- isTRUE(opts$intraspecific_mode) || (length(unique(stats::na.omit(vouchers_df$species_standardized))) < 2L)
  if (is_intra_run) {
    vouchers_df <- resolve_analysis_group(vouchers_df, intraspecific_mode = TRUE)
  } else {
    vouchers_df$analysis_group <- vouchers_df$species_standardized
  }

  anc_grp <- resolve_active_anchor_groups(vouchers_df, intraspecific_mode = is_intra_run)
  vif_res <- screen_collinear_predictors_by_block(
    vouchers_df[, cand_eco, drop = FALSE],
    cor.threshold = 0.70,
    vif.threshold = 5.0,
    anchor_groups = anc_grp
  )
  kept_rf_vars <- if (length(vif_res$kept_vars) > 0L) vif_res$kept_vars else cand_eco

  srf_res <- compute_spatial_rf_mems(
    vouchers_df,
    kept_rf_vars,
    n_trees = opts$n_trees,
    intraspecific_mode = is_intra_run,
    vif_res = vif_res
  )
  vouchers_df$spatial_block_cv_acc <- round(srf_res$Spatial_Block_CV_Accuracy, 4)
  vouchers_df$spatial_autocorr_range_km <- srf_res$spatial_autocorr_range_km
  vouchers_df$rf_predicted_taxon <- srf_res$rf_predicted_taxon
  vouchers_df$rf_predicted_prob <- srf_res$rf_predicted_prob
  vouchers_df$rf_given_taxon_prob <- srf_res$rf_given_taxon_prob
  vouchers_df$rf_prediction_source <- srf_res$rf_prediction_source
  vouchers_df$spatial_block_id <- srf_res$spatial_block_id

  prim_sel <- select_primary_predictors_from_rf(
    rf_res = srf_res,
    vif_res = vif_res,
    df_env = vouchers_df,
    intraspecific_mode = is_intra_run
  )
  primary_env_vars <- prim_sel$primary_env_vars
  primary_clim_vars <- prim_sel$primary_clim_vars
  primary_edaph_vars <- prim_sel$primary_edaph_vars
  rf_var_weights <- prim_sel$rf_var_weights
  message(sprintf(
    "[Primary Predictor Selection] Retained %d empirical primary predictors (clim=%d, edaph=%d, topo=%d): %s",
    length(primary_env_vars), length(primary_clim_vars), length(primary_edaph_vars), length(prim_sel$primary_topo_vars),
    paste(primary_env_vars, collapse = ", ")
  ))

  # Stage 3: Pass Empirically Validated Primary Predictors & RF Importance Weights to Ecotype Partitioning, Niche Overlap & Cross-Modal Consensus
  if (is_intra_run) {
    vouchers_df <- resolve_analysis_group(
      vouchers_df,
      intraspecific_mode = TRUE,
      primary_env_vars = primary_env_vars,
      rf_var_weights = rf_var_weights
    )
  }
  vouchers_df <- execute_crossmodal_consensus(
    vouchers_df,
    intraspecific_mode = is_intra_run,
    rf_res = srf_res,
    primary_env_vars = primary_env_vars,
    primary_clim_vars = primary_clim_vars,
    primary_edaph_vars = primary_edaph_vars,
    rf_var_weights = rf_var_weights
  )
  vouchers_df$primary_env_predictors <- paste(primary_env_vars, collapse = ";")
  vouchers_df$primary_clim_predictors <- paste(primary_clim_vars, collapse = ";")
  vouchers_df$primary_edaph_predictors <- paste(primary_edaph_vars, collapse = ";")

  # Stage 4: Partitioned PCA-env Niche Equivalency & Background Similarity Tests on Empirically Selected Primary Predictors
  niche_res <- run_warren_niche_identity_tests(
    vouchers_df,
    env_vars = kept_rf_vars,
    n_perm = opts$permutations,
    intraspecific_mode = is_intra_run,
    primary_env_vars = primary_env_vars,
    primary_clim_vars = primary_clim_vars,
    primary_edaph_vars = primary_edaph_vars,
    rf_var_weights = rf_var_weights
  )
  edaph_sum_cols <- c("Pair", "Taxon1", "Taxon2", "D_edaph", "I_edaph", "p_equiv_D_edaph", "p_bg_D_edaph")
  edaph_schoener_df <- if (length(primary_edaph_vars) > 0L && !is.null(niche_res$summary) &&
                           is.data.frame(niche_res$summary) && all(edaph_sum_cols %in% names(niche_res$summary))) {
    niche_res$summary[, edaph_sum_cols, drop = FALSE]
  } else {
    compute_edaphic_schoener_d(
      vouchers_df,
      primary_edaph_vars = primary_edaph_vars,
      rf_var_weights = rf_var_weights,
      n_perm = opts$permutations,
      intraspecific_mode = is_intra_run
    )
  }
  niche_res$edaphic_schoener_d <- edaph_schoener_df

  # Stage 5: Export Conflict Flags, Summary Report & Figures (Atomic write-then-replace)
  if (write_csv_atomic(vouchers_df, opts$output_flags, demo = demo_mode)) {
    message(sprintf("Saved %d multimodal conflict flags to %s.", nrow(vouchers_df), opts$output_flags))
  }
  if (!isTRUE(demo_mode) && exists("persist_step06_to_sqlite", mode = "function")) {
    tryCatch(persist_step06_to_sqlite(vouchers_df, opts$vouchers, demo = demo_mode), error = function(e) warning("SQLite Step 06 persistence note: ", e$message))
  }

  niche_summary_df <- if (!is.null(niche_res$summary) && is.data.frame(niche_res$summary)) {
    niche_res$summary
  } else {
    data.frame(Pair = character(0), Taxon1 = character(0), Taxon2 = character(0), Schoener_D = numeric(0), p_value_D = numeric(0), stringsAsFactors = FALSE)
  }
  if (write_csv_atomic(niche_summary_df, opts$output_summary, demo = demo_mode)) {
    message(sprintf("Saved Warren's Niche Identity summary to %s.", opts$output_summary))
  }

  export_spatial_rf_figures(vouchers_df, srf_res, niche_res, opts$output_plot, pheno_res = pheno_res)

  sdm_out_dir <- dirname(dirname(opts$output_plot))
  if (is.null(sdm_out_dir) || !nzchar(sdm_out_dir) || identical(sdm_out_dir, ".")) sdm_out_dir <- "outputs"
  focal_tx_use <- if (!is.null(opts$focal_taxon) && nzchar(trimws(opts$focal_taxon))) {
    standardize_packera_taxon(opts$focal_taxon)
  } else {
    "Packera dubia"
  }
  tryCatch(
    fit_and_project_maxent_sdm_maps(
      df = vouchers_df,
      srf_res = srf_res,
      env_dir = opts$env_dir,
      out_dir = sdm_out_dir,
      kept_vars = kept_rf_vars,
      primary_clim_vars = primary_clim_vars,
      primary_edaph_vars = primary_edaph_vars,
      focal_taxon = focal_tx_use,
      intraspecific_mode = is_intra_run,
      demo = demo_mode
    ),
    error = function(e) warning("Step 06 MaxEnt/Ensemble SDM projection warning: ", e$message)
  )

  message("=== Multimodal Spatial RF & Niche Modeling Analysis Complete ===")
}

if (!interactive() && sys.nframe() == 0L) main()
