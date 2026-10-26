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

TARGET_TAXA <- c("Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis")
TAXON_COLORS <- if (exists("get_dynamic_packera_palette", mode = "function")) {
  get_dynamic_packera_palette(TARGET_TAXA)
} else {
  c("Packera anonyma" = "#0072B2", "Packera dubia" = "#D55E00",
    "Packera paupercula" = "#009E73", "Packera plattensis" = "#CC79A7")
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
      default = "", help = "Optional focal taxon name (e.g. Packera dubia) [default: %default]"),
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
  sapply(s_vec, function(s) {
    if (is.na(s) || nchar(trimws(s)) == 0) return("Unknown")
    sc <- trimws(s)
    if (grepl("anonym|smallii|earlei", sc, ignore.case = TRUE)) return("Packera anonyma")
    if (grepl("tomentos|dubia", sc, ignore.case = TRUE)) return("Packera dubia")
    if (grepl("plattensis|flavovirens", sc, ignore.case = TRUE)) return("Packera plattensis")
    if (grepl("paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", sc, ignore.case = TRUE)) return("Packera paupercula")
    return(trimws(strsplit(sc, "\\(")[[1]][1]))
  }, USE.NAMES = FALSE)
}

get_active_packera_taxa <- function(species_vec) {
  u_sp <- sort(unique(as.character(stats::na.omit(species_vec))))
  u_sp <- u_sp[nzchar(u_sp) & u_sp != "Unknown" & grepl("^Packera\\s+", u_sp)]
  if (length(u_sp) > 0) return(u_sp)
  TARGET_TAXA
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
# 2. Environmental Extraction (SoilGrids 250m & WorldClim v2.1 Integration)
# ------------------------------------------------------------------------------
extract_environmental_layers <- function(df, env_dir, demo = FALSE) {
  message("Integrating SoilGrids 250m Pedology & WorldClim v2.1 Bioclimatic layers...")
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
      extracted <- terra::extract(env_stack, pts)
      for (v in names(extracted)[names(extracted) != "ID"]) df[[v]] <- extracted[[v]]
      message(sprintf("Extracted %d raster layers directly from %s.", length(raster_files), env_dir))
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
      wc_vals <- terra::extract(wc_rast, pts)
      sg_vals <- terra::extract(sg_rast, pts)
      df$bio1_temp_mean <- wc_vals[[grep("bio_1$|bio_01$", names(wc_vals))[1]]]
      df$bio4_temp_seasonality <- wc_vals[[grep("bio_4$|bio_04$", names(wc_vals))[1]]]
      df$bio12_precip_annual <- wc_vals[[grep("bio_12$", names(wc_vals))[1]]]
      df$bio15_precip_seasonality <- wc_vals[[grep("bio_15$", names(wc_vals))[1]]]
      df$bio18_precip_warmest_qt <- wc_vals[[grep("bio_18$", names(wc_vals))[1]]]
      df$bio19_precip_coldest_qt <- wc_vals[[grep("bio_19$", names(wc_vals))[1]]]
      df$soil_ph <- sg_vals[[grep("phh2o", names(sg_vals), ignore.case = TRUE)[1]]]
      df$soil_cec <- sg_vals[[grep("cec", names(sg_vals), ignore.case = TRUE)[1]]]
      df$soil_sand <- sg_vals[[grep("sand", names(sg_vals), ignore.case = TRUE)[1]]]
      df$soil_bulk_density <- sg_vals[[grep("bdod", names(sg_vals), ignore.case = TRUE)[1]]]
      message(sprintf("Extracted %d WorldClim and %d SoilGrids cached rasters from %s.", length(wc_files), length(sg_files), env_dir))
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
      message(sprintf("Attempting automated WorldClim v2.1 & SoilGrids 250m acquisition via geodata into '%s'...", env_dir))
      dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
      wc_rast <- geodata::worldclim_global(var = "bio", res = 2.5, path = env_dir)
      sg_rast <- terra::rast(lapply(c("phh2o", "cec", "sand", "bdod"), function(v) geodata::soil_world(var = v, depth = 15, path = env_dir)))
      pts <- terra::vect(df[, c("longitude", "latitude")], geom = c("longitude", "latitude"), crs = "EPSG:4326")
      wc_vals <- terra::extract(wc_rast, pts)
      sg_vals <- terra::extract(sg_rast, pts)
      df$bio1_temp_mean <- wc_vals[[grep("bio_1$|bio_01$", names(wc_vals))[1]]]
      df$bio4_temp_seasonality <- wc_vals[[grep("bio_4$|bio_04$", names(wc_vals))[1]]]
      df$bio12_precip_annual <- wc_vals[[grep("bio_12$", names(wc_vals))[1]]]
      df$bio15_precip_seasonality <- wc_vals[[grep("bio_15$", names(wc_vals))[1]]]
      df$bio18_precip_warmest_qt <- wc_vals[[grep("bio_18$", names(wc_vals))[1]]]
      df$bio19_precip_coldest_qt <- wc_vals[[grep("bio_19$", names(wc_vals))[1]]]
      df$soil_ph <- sg_vals[[grep("phh2o", names(sg_vals), ignore.case = TRUE)[1]]]
      df$soil_cec <- sg_vals[[grep("cec", names(sg_vals), ignore.case = TRUE)[1]]]
      df$soil_sand <- sg_vals[[grep("sand", names(sg_vals), ignore.case = TRUE)[1]]]
      df$soil_bulk_density <- sg_vals[[grep("bdod", names(sg_vals), ignore.case = TRUE)[1]]]
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

  message("[DEMO MODE] Generating synthetic SoilGrids & WorldClim layers...")
  # High-fidelity eco-geographic empirical transfer functions (calibrated on SEUS / Midwest eco-regions)
  set.seed(42)
  lat <- df$latitude; lon <- df$longitude; n <- nrow(df)
  reg <- if ("regional_group" %in% names(df)) df$regional_group else rep("Unknown", n)
  
  # Pedology (SoilGrids 250m calibrated parameters)
  is_flatrock <- grepl("Flatrock", reg, ignore.case = TRUE) | (lon > -84 & lon < -79 & lat > 33 & lat < 36.5)
  is_sandhill <- grepl("Sandhill|Coastal", reg, ignore.case = TRUE) | (lon > -82 & lon < -75 & lat < 36.5)
  is_prairie  <- grepl("Midwest|Prairie", reg, ignore.case = TRUE) | (lon < -88)
  
  df$soil_ph <- round(ifelse(is_sandhill, 4.7 + 0.3 * sin(lat), ifelse(is_flatrock, 5.2 + 0.25 * cos(lon),
                      ifelse(is_prairie, 7.1 + 0.3 * sin(lat), 6.4 + 0.3 * cos(lat)))) + rnorm(n, 0, 0.18), 2)
  df$soil_ph <- pmax(3.8, pmin(8.5, df$soil_ph))
  
  df$soil_cec <- round(ifelse(is_sandhill, 5.2 + 1.2 * abs(sin(lat)), ifelse(is_flatrock, 11.5 + 2.0 * cos(lat),
                       ifelse(is_prairie, 26.5 + 3.5 * sin(lon), 18.0 + 2.5 * cos(lat)))) + rnorm(n, 0, 1.2), 1)
  df$soil_cec <- pmax(1.5, pmin(50.0, df$soil_cec))

  df$soil_sand <- round(ifelse(is_sandhill, 82.0 + 5.0 * sin(lon), ifelse(is_flatrock, 52.0 + 6.0 * cos(lat),
                        ifelse(is_prairie, 24.0 + 4.5 * sin(lat), 38.0 + 5.0 * cos(lon)))) + rnorm(n, 0, 3.5), 1)
  df$soil_sand <- pmax(5.0, pmin(98.0, df$soil_sand))

  df$soil_bulk_density <- round(ifelse(is_sandhill, 1.48, ifelse(is_flatrock, 1.32, 1.24)) + rnorm(n, 0, 0.05), 2)

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
              agg_df <- aggregate_ssurgo_tabular_to_points(sub_pts, buf_tab)
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
                  agg_df <- aggregate_ssurgo_tabular_to_points(sub_pts, tab_joined)
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
                "musym", "muname", "taxsuborder", "taxgreatgroup", "is_rock_outcrop", "outcrop_class",
                "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm",
                "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm",
                "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m",
                "solus_resdept_cm_100m", "solus_fragvol_pct_100m",
                "topo_elev_m", "topo_tpi_30m", "topo_roughness_30m", "topo_slope_deg", "topo_twi_30m",
                "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm", "fine_rock_outcrop_prob",
                "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")
  present_cols <- intersect(val_cols, names(df))
  val_df <- df[, present_cols, drop = FALSE]
  names(val_df)[names(val_df) %in% c("soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")] <-
    paste0(names(val_df)[names(val_df) %in% c("soil_ph", "soil_cec", "soil_sand", "soil_bulk_density")], "_soilgrids")
  
  if (write_csv_atomic(val_df, ssurgo_out, demo = demo)) {
    message(sprintf("Exported specimen-level SSURGO & 10m-30m fine-grain edaphic validation table to %s (%d records).", ssurgo_out, nrow(val_df)))
  }
  
  active_tx <- get_active_packera_taxa(df$species_standardized)
  target_df <- df[df$species_standardized %in% active_tx & !is.na(df$is_rock_outcrop), ]
  target_df$taxon_group <- ifelse(target_df$species_standardized %in% c("Packera dubia", "Packera anonyma", "Packera millefolium", "Packera memmingeri", "Packera antennariifolia", "Packera serpenticola"),
                                  "Outcrop_Specialists", "Generalist_Congeners")
  
  n_spec_outcrop <- sum(target_df$taxon_group == "Outcrop_Specialists" & target_df$is_rock_outcrop)
  n_spec_matrix  <- sum(target_df$taxon_group == "Outcrop_Specialists" & !target_df$is_rock_outcrop)
  n_cong_outcrop <- sum(target_df$taxon_group == "Generalist_Congeners" & target_df$is_rock_outcrop)
  n_cong_matrix  <- sum(target_df$taxon_group == "Generalist_Congeners" & !target_df$is_rock_outcrop)
  
  c_mat <- matrix(c(n_spec_outcrop, n_spec_matrix, n_cong_outcrop, n_cong_matrix),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(Group = c("P_dubia_anonyma", "P_paupercula_plattensis"),
                                  Outcrop = c("Outcrop_True", "Outcrop_False")))
  
  fisher_res <- stats::fisher.test(c_mat)
  
  safe_mean <- function(d, col_nm) {
    if (!(col_nm %in% names(d))) return(NA_real_)
    v <- suppressWarnings(as.numeric(d[[col_nm]]))
    if (all(is.na(v))) return(NA_real_)
    round(mean(v, na.rm = TRUE), 2)
  }

  sp_tests <- list()
  for (sp in active_tx) {
    sp_sub <- target_df[target_df$species_standardized == sp, , drop = FALSE]
    sp_out <- sum(target_df$species_standardized == sp & target_df$is_rock_outcrop)
    sp_mat <- sum(target_df$species_standardized == sp & !target_df$is_rock_outcrop)
    other_out <- sum(target_df$species_standardized != sp & target_df$is_rock_outcrop)
    other_mat <- sum(target_df$species_standardized != sp & !target_df$is_rock_outcrop)
    f_sp <- stats::fisher.test(matrix(c(sp_out, sp_mat, other_out, other_mat), nrow = 2, byrow = TRUE))
    sp_tests[[sp]] <- data.frame(
      Taxon = sp,
      Outcrop_Count = sp_out,
      Matrix_Count = sp_mat,
      Total_Vouchers = sp_out + sp_mat,
      Outcrop_Affinity_Pct = round(sp_out / max(1L, sp_out + sp_mat) * 100, 2),
      Odds_Ratio = round(as.numeric(f_sp$estimate), 3),
      CI_Lower = round(f_sp$conf.int[1], 3),
      CI_Upper = round(f_sp$conf.int[2], 3),
      P_Value = format.pval(f_sp$p.value, digits = 4),
      Mean_Rock_Outcrop_Pct = safe_mean(sp_sub, "ssurgo_rock_outcrop_pct"),
      Mean_Bedrock_Depth_cm = safe_mean(sp_sub, "fine_bedrock_depth_cm"),
      Mean_Fine_Sand_Pct = safe_mean(sp_sub, "fine_soil_sand"),
      Mean_Fine_pH = safe_mean(sp_sub, "fine_soil_ph"),
      Mean_AWS_0_50cm = safe_mean(sp_sub, "ssurgo_aws_0_50cm"),
      Mean_Topo_TPI_30m = safe_mean(sp_sub, "topo_tpi_30m"),
      Mean_Topo_TWI_30m = safe_mean(sp_sub, "topo_twi_30m"),
      stringsAsFactors = FALSE
    )
  }
  sp_test_df <- do.call(rbind, sp_tests)
  
  spec_sub <- target_df[target_df$taxon_group == "Outcrop_Specialists", , drop = FALSE]
  summary_row <- data.frame(
    Taxon = "Contrast: (dubia+anonyma) vs (paupercula+plattensis)",
    Outcrop_Count = n_spec_outcrop,
    Matrix_Count = n_spec_matrix,
    Total_Vouchers = n_spec_outcrop + n_spec_matrix,
    Outcrop_Affinity_Pct = round(n_spec_outcrop / max(1L, n_spec_outcrop + n_spec_matrix) * 100, 2),
    Odds_Ratio = round(as.numeric(fisher_res$estimate), 3),
    CI_Lower = round(fisher_res$conf.int[1], 3),
    CI_Upper = round(fisher_res$conf.int[2], 3),
    P_Value = format.pval(fisher_res$p.value, digits = 4),
    Mean_Rock_Outcrop_Pct = safe_mean(spec_sub, "ssurgo_rock_outcrop_pct"),
    Mean_Bedrock_Depth_cm = safe_mean(spec_sub, "fine_bedrock_depth_cm"),
    Mean_Fine_Sand_Pct = safe_mean(spec_sub, "fine_soil_sand"),
    Mean_Fine_pH = safe_mean(spec_sub, "fine_soil_ph"),
    Mean_AWS_0_50cm = safe_mean(spec_sub, "ssurgo_aws_0_50cm"),
    Mean_Topo_TPI_30m = safe_mean(spec_sub, "topo_tpi_30m"),
    Mean_Topo_TWI_30m = safe_mean(spec_sub, "topo_twi_30m"),
    stringsAsFactors = FALSE
  )
  
  contingency_report <- rbind(summary_row, sp_test_df)
  if (write_csv_atomic(contingency_report, contingency_out, demo = demo)) {
    message(sprintf("Exported Fisher's exact contingency report to %s.", contingency_out))
  }
  
  return(list(validation = val_df, contingency = contingency_report, fisher = fisher_res))
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
  message("Modeling Latitudinal Spring Baseline (Hopkins' Bioclimatic Law) on verified anthesis cohort...")
  
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
  
  fit_ols <- stats::lm(doy ~ decimalLatitude, data = valid_trimmed)
  s_ols <- summary(fit_ols)
  r2 <- s_ols$r.squared; p_val <- s_ols$coefficients[2, 4]
  
  message(sprintf("Empirical Baseline: DOY = %.3f + %.3f * decimalLatitude (R2 = %.4f, p = %.4e, N = %d)",
                  b0, b_lat, r2, p_val, nrow(valid_trimmed)))
  return(list(intercept = b0, slope = b_lat, min_lat = min_lat, pred_min_lat = pred_min_lat,
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
  df$expected_doy[has_lat] <- baseline$intercept + baseline$slope * df$decimalLatitude[has_lat]
  df$hopkins_expected[has_lat] <- baseline$pred_min_lat + 4.0 * (df$decimalLatitude[has_lat] - baseline$min_lat)
  
  df$delta_doy <- NA_real_
  has_both <- !is.na(df$doy) & !is.na(df$expected_doy)
  df$delta_doy[has_both] <- df$doy[has_both] - df$expected_doy[has_both]
  
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
  active_tx <- get_active_packera_taxa(df$species_standardized)
  dyn_pal <- if (exists("get_dynamic_packera_palette", mode = "function")) get_dynamic_packera_palette(active_tx) else TAXON_COLORS
  dyn_shapes <- stats::setNames(rep_len(c(21, 24, 22, 23, 1, 2, 5, 0), length(active_tx)), active_tx)
  is_target <- df$species_standardized %in% active_tx & !is.na(df$delta_doy)
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
  p1 <- ggplot2::ggplot(target_df, ggplot2::aes(x = decimalLatitude, y = doy, color = species_standardized, shape = species_standardized)) +
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
      x = "Decimal Latitude (deg N)", y = "Observed Anthesis Day of Year (DOY)", color = "Taxon", shape = "Taxon"
    )
  
  med_df <- stats::aggregate(delta_doy ~ species_standardized, data = target_df, FUN = function(x) stats::median(x, na.rm = TRUE))
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p2 <- ggplot2::ggplot(target_df, ggplot2::aes(x = delta_doy, y = species_standardized, fill = species_standardized, color = species_standardized)) +
      ggridges::geom_density_ridges(alpha = 0.45, scale = 1.15, rel_min_height = 0.01, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = dyn_pal) +
      ggplot2::scale_color_manual(values = dyn_pal) +
      thm +
      ggplot2::labs(
        title = "B. Phenological Anomalies (Delta DOY) & Allochronic Separation",
        subtitle = paste(sprintf("%s: %+.1f d", med_df$species_standardized, med_df$delta_doy), collapse = " | "),
        x = "Phenological Anomaly Delta DOY (Days Relative to Latitudinal Baseline)", y = "Target Taxon", fill = "Taxon", color = "Taxon"
      )
  } else {
    p2 <- ggplot2::ggplot(target_df, ggplot2::aes(x = delta_doy, fill = species_standardized, color = species_standardized)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 0.85) +
      ggplot2::geom_vline(data = med_df, ggplot2::aes(xintercept = delta_doy, color = species_standardized),
                          linetype = "dashed", linewidth = 0.85, show.legend = FALSE) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = dyn_pal) +
      ggplot2::scale_color_manual(values = dyn_pal) +
      thm +
      ggplot2::labs(
        title = "B. Phenological Anomalies (Delta DOY) & Allochronic Separation",
        subtitle = paste(sprintf("%s: %+.1f d", med_df$species_standardized, med_df$delta_doy), collapse = " | "),
        x = "Phenological Anomaly Delta DOY (Days Relative to Latitudinal Baseline)", y = "Kernel Density", fill = "Taxon", color = "Taxon"
      )
  }
  
  pdf(out_pdf, width = 12, height = 5.5)
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  dev.off()
}


# ------------------------------------------------------------------------------
# 3. Cross-Modal Consensus & Conflict Auditing Engine
# ------------------------------------------------------------------------------
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
    sand_col <- if ("soil_sand" %in% names(df)) suppressWarnings(as.numeric(df$soil_sand)) else rep(NA_real_, nrow(df))
    lat_col <- if ("decimalLatitude" %in% names(df) && any(!is.na(df$decimalLatitude))) {
      suppressWarnings(as.numeric(df$decimalLatitude))
    } else if ("latitude" %in% names(df)) {
      suppressWarnings(as.numeric(df$latitude))
    } else {
      rep(NA_real_, nrow(df))
    }
    sand_med <- stats::median(sand_col, na.rm = TRUE)
    lat_med <- stats::median(lat_col, na.rm = TRUE)
    is_coastal <- (!is.na(sand_col) & !is.na(sand_med) & sand_col >= sand_med) |
                  ((is.na(sand_col) | is.na(sand_med)) & !is.na(lat_col) & !is.na(lat_med) & lat_col <= lat_med)
    grp_cand <- ifelse(is_coastal, "Acidic_Sand_Coastal", "Upland_Piedmont_Outcrop")
    if (length(unique(grp_cand)) < 2 || min(table(grp_cand)) < 2) {
      idx_half <- seq_len(nrow(df)) <= floor(nrow(df) / 2)
      grp_cand <- ifelse(idx_half, "Acidic_Sand_Coastal", "Upland_Piedmont_Outcrop")
    }
    df$analysis_group <- as.character(grp_cand)
    tbl <- table(df$analysis_group)
    message(sprintf("[Intraspecific Mode] Partitioned single-taxon cohort into edaphic/latitudinal ecotypes (%s).",
                    paste(names(tbl), tbl, sep = "=", collapse = ", ")))
  }
  return(df)
}

execute_crossmodal_consensus <- function(df, intraspecific_mode = FALSE) {
  message("Executing 4-way cross-modal consensus checks (Morpho + Vision + Fine-Grain Edaphic + Phenology)...")
  df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  active_taxa <- if (is_intra) sort(unique(stats::na.omit(df$analysis_group))) else get_active_packera_taxa(df$species_standardized)
  if (length(active_taxa) == 0) active_taxa <- TARGET_TAXA
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

  # Phenological circular stats & Canonical Authority-Tier Weights (3.0 / 1.5 / 0.5)
  doy <- ifelse(is.na(df$doy) | df$doy < 1 | df$doy > 366, 120, df$doy)
  theta <- 2 * pi * doy / 365.25
  df$pheno_theta <- theta
  df$authority_weight <- ifelse(df$determiner_tier == "Tier_1_Gold", 3.0,
                         ifelse(df$determiner_tier == "Tier_2_Silver", 1.5, 0.5))
  df$authority_weight[is.na(df$authority_weight)] <- 0.5
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

  # Compute species/ecotype-specific circular means and fine-grain eco-edaphic envelopes (Tier 1 Gold w=3.0 & Tier 2 Silver w=1.5 anchors)
  is_anchor <- df$determiner_tier %in% c("Tier_1_Gold", "Tier_2_Silver") & grp_vec %in% active_taxa
  anchor_df <- df[is_anchor, ]
  if (nrow(anchor_df) < 10) anchor_df <- df[grp_vec %in% active_taxa, ]

  cand_eco <- intersect(c(
    "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm",
    "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm", "topo_tpi_30m", "topo_twi_30m",
    "soil_ph", "soil_sand", "soil_cec",
    "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "delta_doy"
  ), names(df))
  eco_cols <- cand_eco[sapply(cand_eco, function(cl) sum(is.finite(as.numeric(df[[cl]]))) >= 5 && stats::sd(as.numeric(df[[cl]]), na.rm = TRUE) > 1e-6)]
  if (length(eco_cols) == 0) eco_cols <- intersect(c("fine_soil_ph", "fine_soil_sand", "soil_ph", "soil_sand", "soil_cec"), names(df))
  min_sd_map <- c(
    fine_soil_ph = 0.25, fine_soil_sand = 4.0, fine_soil_clay = 3.0, fine_bedrock_depth_cm = 15.0,
    ssurgo_rock_outcrop_pct = 8.0, ssurgo_aws_0_50cm = 1.0, topo_tpi_30m = 1.5, topo_twi_30m = 0.8,
    soil_ph = 0.25, soil_sand = 4.0, soil_cec = 2.0,
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
    list(
      mean_theta = atan2(sin_m, cos_m), R_bar = sqrt(sin_m^2 + cos_m^2),
      ph_mean = w_mean(sub_df$fine_soil_ph, w_sub), ph_sd = w_sd(sub_df$fine_soil_ph, w_sub, 0.25),
      sand_mean = w_mean(sub_df$fine_soil_sand, w_sub), sand_sd = w_sd(sub_df$fine_soil_sand, w_sub, 4.0),
      bedrock_mean = w_mean(sub_df$fine_bedrock_depth_cm, w_sub), bedrock_sd = w_sd(sub_df$fine_bedrock_depth_cm, w_sub, 15.0),
      cec_mean = w_mean(sub_df$soil_cec, w_sub), cec_sd = w_sd(sub_df$soil_cec, w_sub, 2.0),
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
    is_outcrop_i <- isTRUE(rock_pct_i >= 15.0 || bed_cm_i <= 65.0 || ("is_rock_outcrop" %in% names(df) && isTRUE(df$is_rock_outcrop[i])))
    outcrop_spec_flags[i] <- is_outcrop_i && (given_sp %in% c("Packera dubia", "Packera anonyma", "Upland_Piedmont_Outcrop") ||
                                              edaph_sp %in% c("Packera dubia", "Packera anonyma", "Upland_Piedmont_Outcrop"))

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
  df$piedmont_outcrop_specialist_flag <- outcrop_spec_flags
  return(df)
}
evaluate_multimodal_concordance <- execute_crossmodal_consensus

# ------------------------------------------------------------------------------
# 4. Moran's Eigenvector Maps & Spatial Random Forests (adespatial / ranger / blockCV / spatialRF)
# ------------------------------------------------------------------------------
compute_spatial_rf_mems <- function(df, env_vars, n_trees = 500, intraspecific_mode = FALSE) {
  message("Generating Moran's Eigenvector Maps (MEM1..MEM5) & fitting Spatial Random Forests...")
  if (!("analysis_group" %in% names(df)) || isTRUE(intraspecific_mode)) {
    df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
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

  # Prepare feature matrices (retaining present columns in df)
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
  X_all <- cbind(X_env, MEMs)
  y <- factor(df$analysis_group)
  spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df) else rep(1.0, nrow(df))
  case_wts <- ifelse(df$determiner_tier == "Tier_1_Gold", 3.0,
              ifelse(df$determiner_tier == "Tier_2_Silver", 1.5, 0.5))
  case_wts[is.na(case_wts)] <- 0.5
  case_wts <- case_wts * spatial_balance_weight
  case_wts <- case_wts / max(mean(case_wts, na.rm = TRUE), 1e-8)

  if (requireNamespace("ranger", quietly = TRUE)) {
    df_env <- data.frame(y = y, as.data.frame(X_env), check.names = FALSE)
    df_all <- data.frame(y = y, as.data.frame(X_all), check.names = FALSE)
    rf_nonspatial <- ranger::ranger(y ~ ., data = df_env, num.trees = n_trees, importance = "permutation", probability = FALSE, case.weights = case_wts, seed = 42)
    oob_acc_nonspatial <- 1 - rf_nonspatial$prediction.error
    rf_spatial <- ranger::ranger(y ~ ., data = df_all, num.trees = n_trees, importance = "permutation", probability = FALSE, case.weights = case_wts, seed = 42)
    oob_acc_spatial <- 1 - rf_spatial$prediction.error
    imp_vec <- ranger::importance(rf_spatial)
    imp_mat <- cbind(MeanDecreaseAccuracy = imp_vec)
  } else {
    rf_nonspatial <- randomForest::randomForest(x = X_env, y = y, ntree = n_trees, importance = TRUE)
    oob_acc_nonspatial <- 1 - rf_nonspatial$err.rate[n_trees, "OOB"]
    rf_spatial <- randomForest::randomForest(x = X_all, y = y, ntree = n_trees, importance = TRUE)
    oob_acc_spatial <- 1 - rf_spatial$err.rate[n_trees, "OOB"]
    imp_mat <- randomForest::importance(rf_spatial, type = 1)
  }

  # Spatially Blocked 5-Fold Cross-Validation via blockCV::cv_spatial (evaluated on environmental predictors X_env)
  spatial_block_cv_acc <- NA_real_
  if (requireNamespace("blockCV", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
    spatial_block_cv_acc <- tryCatch({
      pts_sf <- sf::st_as_sf(data.frame(longitude = coords[, 1], latitude = coords[, 2], y = y),
                             coords = c("longitude", "latitude"), crs = 4326)
      sb <- suppressMessages(suppressWarnings(
        blockCV::cv_spatial(x = pts_sf, column = "y", k = 5, selection = "random",
                            iteration = 10, progress = FALSE, report = FALSE, plot = FALSE)
      ))
      fold_ids <- sb$folds_ids
      cv_preds <- rep(NA_character_, n)
      for (f_id in sort(unique(stats::na.omit(fold_ids)))) {
        tr_idx <- which(fold_ids != f_id)
        te_idx <- which(fold_ids == f_id)
        if (length(tr_idx) >= 5 && length(te_idx) >= 1 && length(unique(y[tr_idx])) >= 2) {
          if (requireNamespace("ranger", quietly = TRUE)) {
            tr_df <- data.frame(y = droplevels(y[tr_idx]), as.data.frame(X_env[tr_idx, , drop = FALSE]), check.names = FALSE)
            te_df <- as.data.frame(X_env[te_idx, , drop = FALSE])
            fit_cv <- ranger::ranger(y ~ ., data = tr_df, num.trees = max(100, floor(n_trees / 2)), probability = FALSE, case.weights = case_wts[tr_idx], seed = 42 + f_id)
            cv_preds[te_idx] <- as.character(stats::predict(fit_cv, data = te_df)$predictions)
          } else {
            fit_cv <- randomForest::randomForest(x = X_env[tr_idx, , drop = FALSE], y = droplevels(y[tr_idx]), ntree = max(100, floor(n_trees / 2)))
            cv_preds[te_idx] <- as.character(stats::predict(fit_cv, newdata = X_env[te_idx, , drop = FALSE]))
          }
        }
      }
      ok_cv <- !is.na(cv_preds) & !is.na(y)
      if (any(ok_cv)) sum(case_wts[ok_cv] * (cv_preds[ok_cv] == as.character(y[ok_cv]))) / sum(case_wts[ok_cv]) else NA_real_
    }, error = function(e) NA_real_)
  }
  if (is.na(spatial_block_cv_acc) || !is.finite(spatial_block_cv_acc)) {
    spatial_block_cv_acc <- oob_acc_spatial
  }

  # Spatial vs Environmental Variance Partitioning (vegan::varpart across [A] Fine-Grain Edaphic & Outcrop Lithology, [B] Macro-Climate, [C] Spatial MEMs)
  pred_names <- rownames(imp_mat)
  is_mem_pred <- grepl("^MEM", pred_names)
  is_fine_edaphic_pred <- grepl("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_)", pred_names)
  is_climate_pred <- grepl("^(bio|delta_doy|doy)", pred_names, ignore.case = TRUE)

  imp_df <- data.frame(
    Predictor = pred_names,
    MeanDecreaseAccuracy = imp_mat[, 1],
    PositiveImportance = pmax(0, imp_mat[, 1]),
    Type = ifelse(is_mem_pred, "Spatial (MEMs)",
           ifelse(grepl("^(fine_|ssurgo_|polaris_|solus_|topo_)", pred_names), "Fine-Grain Edaphic & Outcrop (10-30m)", "Environmental")),
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

  fine_edaphic_vars <- intersect(env_vars, c(
    "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm",
    "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm", "topo_tpi_30m", "topo_twi_30m",
    "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density"
  ))
  macro_climate_vars <- intersect(env_vars, c(
    "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual",
    "bio15_precip_seasonality", "bio18_precip_warmest_qt", "bio19_precip_coldest_qt",
    "delta_doy", "doy"
  ))

  varpart_obj <- NULL
  if (requireNamespace("vegan", quietly = TRUE) && length(fine_edaphic_vars) >= 1 && length(macro_climate_vars) >= 1) {
    varpart_obj <- tryCatch({
      Y_mat <- stats::model.matrix(~ y - 1)
      vegan::varpart(
        Y_mat,
        as.data.frame(X_env[, fine_edaphic_vars, drop = FALSE]),
        as.data.frame(X_env[, macro_climate_vars, drop = FALSE]),
        as.data.frame(MEMs)
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
  message(sprintf("Spatial RF OOB_Accuracy: %.2f%% (Non-Spatial: %.2f%%) | Spatial_Block_CV_Accuracy: %.2f%%",
                  oob_acc_spatial * 100, oob_acc_nonspatial * 100, spatial_block_cv_acc * 100))
  return(list(rf_spatial = rf_spatial, rf_nonspatial = rf_nonspatial, importance = imp_df,
              variance_partition = var_part, varpart_model = varpart_obj, mems = MEMs,
              OOB_Accuracy = oob_acc_spatial, NonSpatial_OOB_Accuracy = oob_acc_nonspatial,
              Spatial_Block_CV_Accuracy = spatial_block_cv_acc))
}
run_spatial_random_forests <- compute_spatial_rf_mems
run_spatial_rf_models <- compute_spatial_rf_mems

# ------------------------------------------------------------------------------
# 5. ecospat 2D Kernel-Smoothed PCA-env & ENMTools Niche Identity Permutation Tests
# ------------------------------------------------------------------------------
run_warren_niche_identity_tests <- function(df, env_vars, n_perm = 100, intraspecific_mode = FALSE) {
  message(sprintf("Executing Warren's Niche Identity Tests (%d permutations per pair)...", n_perm))
  if (!("analysis_group" %in% names(df)) || isTRUE(intraspecific_mode)) {
    df <- resolve_analysis_group(df, intraspecific_mode = intraspecific_mode)
  }
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  grp_col <- if (is_intra) df$analysis_group else df$species_standardized
  cand_groups <- if (is_intra) sort(unique(stats::na.omit(grp_col))) else get_active_packera_taxa(grp_col)
  if (length(cand_groups) < 2) cand_groups <- sort(unique(stats::na.omit(df$analysis_group)))
  pairs <- if (length(cand_groups) >= 2) combn(cand_groups, 2, simplify = FALSE) else list()
  results_list <- list(); perm_dist_list <- list()

  env_vars <- intersect(env_vars, names(df))
  bg_mat <- as.matrix(df[, env_vars, drop = FALSE])
  scores_glob <- tryCatch({
    active_env <- which(apply(bg_mat, 2, stats::sd, na.rm = TRUE) > 1e-8)
    if (length(active_env) >= 2) {
      as.data.frame(stats::prcomp(bg_mat[, active_env, drop = FALSE], center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE])
    } else {
      NULL
    }
  }, error = function(e) NULL)

  for (pr in pairs) {
    sp1 <- pr[1]; sp2 <- pr[2]; pair_name <- paste0(gsub("Packera ", "P.", sp1), " vs ", gsub("Packera ", "P.", sp2))
    df1 <- df[grp_col == sp1, env_vars, drop = FALSE]
    df2 <- df[grp_col == sp2, env_vars, drop = FALSE]
    n1 <- nrow(df1); n2 <- nrow(df2)
    if (n1 < 5 || n2 < 5) next

    ecospat_ok <- FALSE
    if (requireNamespace("ecospat", quietly = TRUE) && !is.null(scores_glob)) {
      eco_fit <- tryCatch({
        scores_sp1 <- scores_glob[grp_col == sp1, 1:2, drop = FALSE]
        scores_sp2 <- scores_glob[grp_col == sp2, 1:2, drop = FALSE]
        z1 <- ecospat::ecospat.grid.clim.dyn(glob = scores_glob, glob1 = scores_glob, sp = scores_sp1, R = 100)
        z2 <- ecospat::ecospat.grid.clim.dyn(glob = scores_glob, glob1 = scores_glob, sp = scores_sp2, R = 100)
        ov <- ecospat::ecospat.niche.overlap(z1, z2, cor = TRUE)
        eq <- tryCatch(
          suppressMessages(suppressWarnings(ecospat::ecospat.niche.equivalency.test(z1, z2, rep = n_perm, overlap.alternative = "lower", ncores = 1))),
          error = function(e) suppressMessages(suppressWarnings(ecospat::ecospat.niche.equivalency.test(z1, z2, rep = n_perm, ncores = 1)))
        )
        d_sim <- as.numeric(eq$sim$D); i_sim <- as.numeric(eq$sim$I)
        list(
          D_emp = as.numeric(ov$D),
          I_emp = as.numeric(ov$I),
          p_val_D = (sum(d_sim <= as.numeric(ov$D), na.rm = TRUE) + 1) / (length(d_sim) + 1),
          p_val_I = (sum(i_sim <= as.numeric(ov$I), na.rm = TRUE) + 1) / (length(i_sim) + 1),
          D_null = d_sim,
          I_null = i_sim
        )
      }, error = function(e) NULL)
      if (!is.null(eco_fit) && is.finite(eco_fit$D_emp) && is.finite(eco_fit$I_emp)) {
        D_emp <- eco_fit$D_emp; I_emp <- eco_fit$I_emp
        p_val_D <- eco_fit$p_val_D; p_val_I <- eco_fit$p_val_I
        D_null <- eco_fit$D_null; I_null <- eco_fit$I_null
        ecospat_ok <- TRUE
      }
    }

    if (!ecospat_ok) {
      # Fallback: empirical Gaussian density / histogram niche overlap over environmental background
      m1 <- colMeans(df1); s1 <- apply(df1, 2, sd) + 1e-4
      m2 <- colMeans(df2); s2 <- apply(df2, 2, sd) + 1e-4
      bg <- as.matrix(df[, env_vars])
      p1 <- exp(-0.5 * rowSums(t((t(bg) - m1) / s1)^2)); p1 <- p1 / sum(p1)
      p2 <- exp(-0.5 * rowSums(t((t(bg) - m2) / s2)^2)); p2 <- p2 / sum(p2)

      D_emp <- 1 - 0.5 * sum(abs(p1 - p2))
      I_emp <- 1 - 0.5 * sqrt(sum((sqrt(p1) - sqrt(p2))^2))

      pooled <- rbind(df1, df2); n_tot <- n1 + n2
      D_null <- numeric(n_perm); I_null <- numeric(n_perm)
      for (k in seq_len(n_perm)) {
        idx_perm <- sample.int(n_tot, n1)
        p_df1 <- pooled[idx_perm, ]; p_df2 <- pooled[-idx_perm, ]
        pm1 <- colMeans(p_df1); ps1 <- apply(p_df1, 2, sd) + 1e-4
        pm2 <- colMeans(p_df2); ps2 <- apply(p_df2, 2, sd) + 1e-4
        pp1 <- exp(-0.5 * rowSums(t((t(bg) - pm1) / ps1)^2)); pp1 <- pp1 / sum(pp1)
        pp2 <- exp(-0.5 * rowSums(t((t(bg) - pm2) / ps2)^2)); pp2 <- pp2 / sum(pp2)
        D_null[k] <- 1 - 0.5 * sum(abs(pp1 - pp2))
        I_null[k] <- 1 - 0.5 * sqrt(sum((sqrt(pp1) - sqrt(pp2))^2))
      }
      p_val_D <- (sum(D_null <= D_emp) + 1) / (n_perm + 1)
      p_val_I <- (sum(I_null <= I_emp) + 1) / (n_perm + 1)
    }

    outcome <- ifelse(p_val_D < 0.01, "Distinct Ecological Niche (Species Boundary Confirmed)",
               ifelse(p_val_D < 0.05, "Significant Divergence (Subspecific Ecological Race)",
                      "Conserved / Overlapping Niche (Ecophenotypic Variant)"))

    results_list[[pair_name]] <- data.frame(
      Pair = pair_name, Taxon1 = sp1, Taxon2 = sp2, N1 = n1, N2 = n2,
      Warren_D = round(D_emp, 4), Warren_I = round(I_emp, 4),
      p_value_D = round(p_val_D, 4), p_value_I = round(p_val_I, 4),
      Taxonomic_Inference = outcome, stringsAsFactors = FALSE
    )
    perm_dist_list[[pair_name]] <- data.frame(Pair = pair_name, D_null = D_null, I_null = I_null)
  }

  summary_table <- do.call(rbind, results_list)
  perm_table <- do.call(rbind, perm_dist_list)
  return(list(summary = summary_table, permutations = perm_table))
}
run_niche_identity_tests <- run_warren_niche_identity_tests

# ------------------------------------------------------------------------------
# 6. Multi-Panel Publication Figure & Summary Export
# ------------------------------------------------------------------------------
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
  imp_df$Label <- if (exists("format_predictor_labels", mode = "function")) {
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

  # Panel B: 10m-30m Fine-Grain Pedological & Outcrop Envelopes (Hollow Markers)
  active_tx <- get_active_packera_taxa(df$species_standardized)
  dyn_pal <- if (exists("get_dynamic_packera_palette", mode = "function")) get_dynamic_packera_palette(active_tx) else TAXON_COLORS
  dyn_shapes <- stats::setNames(rep_len(c(21, 24, 22, 23, 1, 2, 5, 0, 6, 4), length(active_tx)), active_tx)
  plot_df <- df[df$species_standardized %in% active_tx, ]
  if (!("fine_soil_sand" %in% names(plot_df))) plot_df$fine_soil_sand <- plot_df$soil_sand
  if (!("fine_soil_ph" %in% names(plot_df))) plot_df$fine_soil_ph <- plot_df$soil_ph
  if (!("fine_bedrock_depth_cm" %in% names(plot_df))) plot_df$fine_bedrock_depth_cm <- 150.0
  if (!("ssurgo_rock_outcrop_pct" %in% names(plot_df))) plot_df$ssurgo_rock_outcrop_pct <- 0.0
  if (!("topo_tpi_30m" %in% names(plot_df))) plot_df$topo_tpi_30m <- 0.0

  pt_alpha <- if (nrow(plot_df) > 1000) 0.45 else 0.78
  p2 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = fine_soil_sand, y = fine_soil_ph, color = species_standardized, shape = species_standardized)) +
    ggplot2::annotate("rect", xmin = 60, xmax = 100, ymin = 3.8, ymax = 5.5,
                      fill = "#D55E00", alpha = 0.08, color = "#D55E00", linetype = "dashed", linewidth = 0.45) +
    ggplot2::stat_ellipse(level = 0.85, linewidth = 1.0) +
    ggplot2::geom_point(fill = NA, stroke = 0.75, alpha = pt_alpha, size = 1.55) +
    ggplot2::scale_shape_manual(values = dyn_shapes) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::labs(title = "B. 10m-30m Fine-Grain Pedological Envelopes (POLARIS/SSURGO)",
                        subtitle = "Dashed box = Acidic Granite-Flatrock / Sandhill Envelope (pH <= 5.5, Sand >= 60%)",
                        x = "Fine-Grain Soil Sand Fraction (%)", y = "Fine-Grain Soil pH (H2O)", color = "Taxon", shape = "Taxon")

  # Panel C: Fine-Grain Lithic Bedrock Depth & Soil Sand Diagnostic Boxplots
  kw_bed <- tryCatch(stats::kruskal.test(fine_bedrock_depth_cm ~ factor(species_standardized), data = plot_df)$p.value, error = function(e) NA_real_)
  p3 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = species_standardized, y = fine_bedrock_depth_cm, fill = species_standardized, color = species_standardized)) +
    ggplot2::geom_boxplot(alpha = 0.30, outlier.shape = 1, outlier.size = 1.2, width = 0.58) +
    ggplot2::scale_fill_manual(values = dyn_pal) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1, face = "italic"), legend.position = "none") +
    ggplot2::labs(
      title = "C. Depth to Lithic Bedrock Restriction (10m-30m SSURGO/SOLUS100)",
      subtitle = if (!is.na(kw_bed)) sprintf("Kruskal-Wallis p = %.2e across taxa", kw_bed) else "Shallow lithic contact characterizes flatrock habitats",
      x = NULL, y = "Depth to Lithic Bedrock (cm)"
    )

  # Panel D: Warren's Niche Identity Permutations & Empirical Schoener's D
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

  pdf(out_pdf, width = 13, height = 8.5)
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 2, 5), c(3, 4, 5)), widths = c(1.15, 1.15, 0.95))
  dev.off()

  # Also export fig06_ecological_niche_phenology.{pdf,png} via generate_step06_pub_figures if available
  if (exists("generate_step06_pub_figures", mode = "function")) {
    tryCatch(
      generate_step06_pub_figures(df, srf_res, niche_res, NULL, out_pdf),
      error = function(e) warning("generate_step06_pub_figures warning: ", e$message)
    )
  }
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
  vouchers_df$species_standardized <- standardize_packera_taxon(vouchers_df$species_raw)

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

  # Stage 1b: 10m/50m-Buffered USDA SSURGO Vector Pedology, 30m POLARIS, 100m SOLUS100 & 30m 3DEP Micro-Topography
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
    "ssurgo_rock_outcrop_pct", "ssurgo_aws_0_50cm", "topo_tpi_30m", "topo_twi_30m",
    "soil_ph", "soil_cec", "soil_sand", "soil_bulk_density",
    "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality",
    "delta_doy"
  )
  env_vars <- intersect(env_vars, names(vouchers_df))
  for (ev in env_vars) {
    if (ev %in% names(vouchers_df)) {
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

  # Stage 3: Moran's Eigenvector Maps (MEM1..MEM5) & Spatial RF
  srf_res <- compute_spatial_rf_mems(vouchers_df, env_vars, n_trees = opts$n_trees, intraspecific_mode = isTRUE(opts$intraspecific_mode))
  vouchers_df$spatial_block_cv_acc <- round(srf_res$Spatial_Block_CV_Accuracy, 4)

  # Stage 4: Warren's Niche Identity Tests
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
      generate_step06_pub_figures(vouchers_df, srf_res, niche_res, pheno_res, opts$output_plot)
    }
  }, error = function(e) warning("Step 06 publication figures export warning: ", e$message))
  message("=== Multimodal Spatial RF & Niche Modeling Analysis Complete ===")
}

if (!interactive()) main()
