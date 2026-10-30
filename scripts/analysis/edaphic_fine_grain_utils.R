#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/edaphic_fine_grain_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Fine-grain (10 m to 30 m) edaphic, lithic, geologic, and micro-topographic
#   predictor extraction, uncertainty-buffered herbarium sampling, and caching
#   utilities for Step 06 (06_multimodal_spatial_rf.R):
#     1. Trapezoidal 0-15 cm herbaceous rhizosphere integration
#        ((5 * X_0_5 + 10 * X_5_15) / 15) across 30 m POLARIS, 30 m gSSURGO,
#        100 m SOLUS100, and 250 m SoilGrids fallback (Chaney et al. 2019;
#        Velazco et al. 2017).
#     2. Physical outcrop predictors: Depth to Bedrock (bdricm / resdept_r),
#        Coarse Fragment Volume (cfvo), Available Water Capacity (awc), and
#        compositional log-ratio soil texture (alr_sand_clay, alr_silt_clay)
#        alongside pH and CEC.
#     3. USGS State Geologic Map Compilation (SGMC; Horton et al. 2017)
#        lithological standardization into 6 geochemical functional classes and
#        distance/fractional proximity to azonal rock contacts.
#     4. Uncertainty-buffered herbarium extraction (Marcer et al. 2022) with
#        adaptive buffer radius r_i = clamp(coordinateUncertaintyInMeters, 100, 1000)
#        (default 250 m), excluding coordinateUncertaintyInMeters > 10000 from
#        model training while retaining them for passive Step 07 coalescing.
#     5. Block-partitioned collinearity pruning (Dormann et al. 2013) via
#        spatialRF::auto_cor(cor.threshold = 0.70) and
#        spatialRF::auto_vif(vif.threshold = 5.0).
# ==============================================================================

SSURGO_OUTCROP_SERIES_REGEX <- "(?i)(rock outcrop|wake|louisburg|ashlar|rion|wateree|wilkes|iredell|mecklenburg|poindexter|enon|goldston|manteo|boulders|granite|diabase)"
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
}

SGMC_LITHOLOGY_CLASSES <- c(
  "Ultramafic_Serpentine",
  "Mafic_Igneous_Metamorphic",
  "Carbonate_Calcareous",
  "Felsic_Granitic_Gneiss",
  "Siliciclastic_Consolidated",
  "Unconsolidated_Sedimentary"
)

SSURGO_OUTCROP_CLASSES <- c(
  "Granitic flatrock",
  "Limestone glade",
  "Sandstone/Siliceous",
  "Serpentine",
  "Diabase/Mafic glade",
  "Non-outcrop matrix"
)

CORE_FINE_GRAIN_CACHE_COLS <- c(
  "catalogNumber", "lat_lon_key", "latitude", "longitude",
  "musym", "muname", "taxsuborder", "taxgreatgroup", "is_rock_outcrop", "outcrop_class",
  "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct", "ssurgo_clay_pct",
  "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm",
  "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m",
  "polaris_bd_30m", "polaris_ksat_30m", "solus_resdept_cm", "solus_fragvol_pct",
  "topo_tpi_30m", "topo_roughness_30m", "topo_slope_deg", "topo_twi_30m",
  "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm", "fine_rock_outcrop_prob"
)

FINE_GRAIN_CACHE_COLS <- c(
  CORE_FINE_GRAIN_CACHE_COLS,
  "bdricm", "resdept_r", "cfvo", "awc", "fine_soil_cec",
  "soil_silt_pct", "alr_sand_clay", "alr_silt_clay",
  "sgmc_lithology_class", "is_azonal_lithology",
  "azonal_contact_dist_m", "azonal_proximity_index",
  "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
  "ph_q10_buffer", "ph_q90_buffer", "min_bedrock_depth_buffer_cm",
  "max_cfvo_buffer_pct", "azonal_lithology_presence_buffer"
)

format_lat_lon_key <- function(lat, lon) {
  lat_num <- suppressWarnings(as.numeric(lat))
  lon_num <- suppressWarnings(as.numeric(lon))
  ifelse(
    is.na(lat_num) | is.na(lon_num),
    NA_character_,
    sprintf("%.5f_%.5f", round(lat_num, 5), round(lon_num, 5))
  )
}

# ------------------------------------------------------------------------------
# 1. Rhizosphere Integration (0-15 cm Trapezoidal Depth-Weighted Averaging)
#    & Compositional Log-Ratio Soil Texture (Chaney et al. 2019; Velazco et al. 2017)
# ------------------------------------------------------------------------------
compute_rhizosphere_0_15cm <- function(X_0_5, X_5_15) {
  # Trapezoidal depth-weighted average across the active 0-15 cm herbaceous
  # rhizosphere: (5 * X_0_5 + 10 * X_5_15) / 15
  x1 <- suppressWarnings(as.numeric(X_0_5))
  x2 <- suppressWarnings(as.numeric(X_5_15))
  n <- max(length(x1), length(x2))
  if (n == 0L) return(numeric(0))
  if (length(x1) != n) x1 <- rep_len(x1, n)
  if (length(x2) != n) x2 <- rep_len(x2, n)

  ok1 <- !is.na(x1) & is.finite(x1)
  ok2 <- !is.na(x2) & is.finite(x2)

  out <- rep(NA_real_, n)
  both <- ok1 & ok2
  if (any(both)) {
    out[both] <- (5.0 * x1[both] + 10.0 * x2[both]) / 15.0
  }
  only1 <- ok1 & !ok2
  if (any(only1)) out[only1] <- x1[only1]
  only2 <- !ok1 & ok2
  if (any(only2)) out[only2] <- x2[only2]
  return(out)
}
trapezoidal_rhizosphere_average <- compute_rhizosphere_0_15cm

compute_soil_texture_log_ratios <- function(sand_pct, clay_pct, silt_pct = NULL) {
  # Compositional additive log-ratio (ALR) soil texture transformation
  # breaking simplex closure (sand + silt + clay = 100%)
  s <- suppressWarnings(as.numeric(sand_pct))
  c <- suppressWarnings(as.numeric(clay_pct))
  n <- max(length(s), length(c))
  if (n == 0L) {
    return(data.frame(
      soil_silt_pct = numeric(0),
      alr_sand_clay = numeric(0),
      alr_silt_clay = numeric(0)
    ))
  }
  if (length(s) != n) s <- rep_len(s, n)
  if (length(c) != n) c <- rep_len(c, n)

  s_clean <- pmax(0.5, pmin(99.0, ifelse(is.finite(s), s, 50.0)))
  c_clean <- pmax(0.5, pmin(99.0, ifelse(is.finite(c), c, 20.0)))

  if (is.null(silt_pct)) {
    si_clean <- pmax(0.5, 100.0 - s_clean - c_clean)
  } else {
    si <- suppressWarnings(as.numeric(silt_pct))
    if (length(si) != n) si <- rep_len(si, n)
    si_clean <- ifelse(is.finite(si) & si > 0, pmax(0.5, pmin(99.0, si)), pmax(0.5, 100.0 - s_clean - c_clean))
  }

  alr_sc  <- log((s_clean + 0.5) / (c_clean + 0.5))
  alr_sic <- log((si_clean + 0.5) / (c_clean + 0.5))
  both_na <- (!is.finite(s)) & (!is.finite(c))
  if (any(both_na)) {
    si_clean[both_na] <- NA_real_
    alr_sc[both_na]   <- NA_real_
    alr_sic[both_na]  <- NA_real_
  }

  data.frame(
    soil_silt_pct = round(si_clean, 2),
    alr_sand_clay = round(alr_sc, 4),
    alr_silt_clay = round(alr_sic, 4),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. USGS State Geologic Map Compilation (SGMC; Horton et al. 2017)
#    Lithological Standardization & Azonal Contact Distance
# ------------------------------------------------------------------------------
standardize_sgmc_lithology <- function(
  rocktype1 = "",
  rocktype2 = "",
  major1 = "",
  major2 = "",
  minor1 = "",
  unit_name = "",
  muname = "",
  outcrop_class = "",
  lat = NA_real_,
  lon = NA_real_
) {
  # Standardizes raw USGS SGMC lithology attributes and/or SSURGO parent-material
  # descriptors into 6 canonical geochemical functional classes (Horton et al. 2017):
  #   1. Ultramafic_Serpentine
  #   2. Mafic_Igneous_Metamorphic
  #   3. Carbonate_Calcareous
  #   4. Felsic_Granitic_Gneiss
  #   5. Siliciclastic_Consolidated
  #   6. Unconsolidated_Sedimentary
  n <- max(
    length(rocktype1), length(rocktype2), length(major1), length(major2),
    length(minor1), length(unit_name), length(muname), length(outcrop_class),
    length(lat), length(lon), 1L
  )
  rep_str <- function(x) {
    if (is.null(x) || length(x) == 0L) return(rep("", n))
    y <- as.character(x)
    y[is.na(y)] <- ""
    rep_len(y, n)
  }
  r1 <- rep_str(rocktype1); r2 <- rep_str(rocktype2)
  m1 <- rep_str(major1);    m2 <- rep_str(major2); mn1 <- rep_str(minor1)
  un <- rep_str(unit_name); mu <- rep_str(muname); oc  <- rep_str(outcrop_class)
  lat_v <- rep_len(suppressWarnings(as.numeric(lat)), n)
  lon_v <- rep_len(suppressWarnings(as.numeric(lon)), n)

  out_cls <- character(n)
  for (i in seq_len(n)) {
    geol_txt <- tolower(trimws(paste(r1[i], r2[i], m1[i], m2[i], mn1[i], un[i])))
    soil_txt <- tolower(trimws(paste(mu[i], oc[i])))
    combined <- paste(geol_txt, soil_txt)

    if (grepl("ultramafic|serpentin|peridotit|dunit|pyroxenit|soapstone|talc|ophiolit|chrome|conowingo|aldino", combined, perl = TRUE) ||
        tolower(oc[i]) == "serpentine") {
      out_cls[i] <- "Ultramafic_Serpentine"
    } else if (grepl("mafic|diabas|gabbro|basalt|amphibolit|greenstone|metavolcanic|metadiorit|norite|anorthosit|trap rock|iredell|mecklenburg|wilkes|poindexter|enon", combined, perl = TRUE) ||
               tolower(oc[i]) == "diabase/mafic glade") {
      out_cls[i] <- "Mafic_Igneous_Metamorphic"
    } else if (grepl("carbonat|limeston|dolostone|dolomit|marble|chalk|marl|calcarenit|travertin|karst|coquina|alvar|gladeville|barfield|talbott|colbert|rendolls", combined, perl = TRUE) ||
               tolower(oc[i]) == "limestone glade") {
      out_cls[i] <- "Carbonate_Calcareous"
    } else if (grepl("felsic|granit|granodiorit|monzonit|quartz monzonit|gneiss|pegmatit|rhyolit|aplit|dacit|migmatit|wake|louisburg|ashlar|rion|wateree|pacolet|appling|cecil|wedowee|madison", combined, perl = TRUE) ||
               tolower(oc[i]) == "granitic flatrock") {
      out_cls[i] <- "Felsic_Granitic_Gneiss"
    } else if (grepl("siliciclastic|sandstone|quartzit|arkos|conglomerat|greywack|graywack|siltston|shale|slate|phyllit|mudston|argillit|chert|novaculit|schist|ramsey|lily|sewanee|goldston|manteo|wallen|dekab", combined, perl = TRUE) ||
               tolower(oc[i]) == "sandstone/siliceous") {
      out_cls[i] <- "Siliciclastic_Consolidated"
    } else if (grepl("unconsolidated|alluvi|colluvi|coastal plain|terrace|loess|till|outwash|glacial|dune|eolian|sandhill|muck|peat|saprist|psamment|wagram|lakeland|plainfield|goldsboro|norfolk|myatt|rains|lynchburg|tama", combined, perl = TRUE)) {
      out_cls[i] <- "Unconsolidated_Sedimentary"
    } else {
      # Physiographic province fallback when no text descriptor is present
      la <- lat_v[i]; lo <- lon_v[i]
      if (!is.na(la) && !is.na(lo)) {
        if (lo > -81.5 && la < 36.8 && (lo > -79.2 || la < 33.2)) {
          out_cls[i] <- "Unconsolidated_Sedimentary"
        } else if (lo >= -84.5 && lo <= -77.5 && la >= 32.5 && la <= 37.5) {
          out_cls[i] <- "Felsic_Granitic_Gneiss"
        } else if (lo >= -88.5 && lo < -83.0 && la >= 34.5 && la <= 38.5) {
          out_cls[i] <- "Carbonate_Calcareous"
        } else if (lo >= -84.0 && lo <= -78.0 && la > 35.0 && la <= 41.0) {
          out_cls[i] <- "Siliciclastic_Consolidated"
        } else {
          out_cls[i] <- "Unconsolidated_Sedimentary"
        }
      } else {
        out_cls[i] <- "Unconsolidated_Sedimentary"
      }
    }
  }
  return(out_cls)
}
classify_sgmc_lithology <- standardize_sgmc_lithology

extract_sgmc_lithology <- function(df, env_dir = "data/environmental", demo = FALSE) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0L) return(df)

  sgmc_dir <- file.path(env_dir, "sgmc")
  raw_rock1 <- if ("sgmc_rocktype1" %in% names(df)) as.character(df$sgmc_rocktype1) else rep("", n)
  raw_rock2 <- if ("sgmc_rocktype2" %in% names(df)) as.character(df$sgmc_rocktype2) else rep("", n)
  raw_maj1  <- if ("sgmc_major1" %in% names(df)) as.character(df$sgmc_major1) else rep("", n)
  raw_unit  <- if ("sgmc_unit_name" %in% names(df)) as.character(df$sgmc_unit_name) else rep("", n)
  contact_dist_vec <- rep(NA_real_, n)

  lat_v <- suppressWarnings(as.numeric(df$latitude))
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  valid_pts <- !is.na(lat_v) & !is.na(lon_v)

  # Check if local USGS SGMC vector / GeoPackage / Shapefile layers exist
  if (!isTRUE(demo) && any(valid_pts) && requireNamespace("sf", quietly = TRUE)) {
    cand_dirs <- unique(c(sgmc_dir, env_dir))
    sgmc_files <- unlist(lapply(cand_dirs, function(d) {
      if (!dir.exists(d)) return(character(0))
      list.files(d, pattern = "sgmc.*\\.(gpkg|shp|geojson)$", full.names = TRUE, ignore.case = TRUE)
    }))
    if (length(sgmc_files) > 0L) {
      tryCatch({
        sgmc_sf <- sf::st_read(sgmc_files[1], quiet = TRUE)
        pts_sf <- sf::st_as_sf(df[valid_pts, c("longitude", "latitude"), drop = FALSE],
                               coords = c("longitude", "latitude"), crs = 4326)
        pts_sf <- sf::st_transform(pts_sf, sf::st_crs(sgmc_sf))
        hits <- sf::st_intersects(pts_sf, sgmc_sf)
        v_idx <- which(valid_pts)
        cn_up <- toupper(names(sgmc_sf))
        col_r1 <- names(sgmc_sf)[match("ROCKTYPE1", cn_up)]
        col_r2 <- names(sgmc_sf)[match("ROCKTYPE2", cn_up)]
        col_m1 <- names(sgmc_sf)[match("MAJOR1", cn_up)]
        col_un <- names(sgmc_sf)[match("UNIT_NAME", cn_up)]
        for (k in seq_along(v_idx)) {
          h_k <- hits[[k]]
          if (length(h_k) > 0L) {
            row_g <- h_k[1]
            if (!is.na(col_r1)) raw_rock1[v_idx[k]] <- as.character(sgmc_sf[[col_r1]][row_g])
            if (!is.na(col_r2)) raw_rock2[v_idx[k]] <- as.character(sgmc_sf[[col_r2]][row_g])
            if (!is.na(col_m1)) raw_maj1[v_idx[k]]  <- as.character(sgmc_sf[[col_m1]][row_g])
            if (!is.na(col_un)) raw_unit[v_idx[k]]  <- as.character(sgmc_sf[[col_un]][row_g])
          }
        }
      }, error = function(e) NULL)
    }
  }

  muname_v <- if ("muname" %in% names(df)) as.character(df$muname) else rep("", n)
  ocl_v    <- if ("outcrop_class" %in% names(df)) as.character(df$outcrop_class) else rep("", n)

  lith_cls <- standardize_sgmc_lithology(
    rocktype1 = raw_rock1,
    rocktype2 = raw_rock2,
    major1 = raw_maj1,
    unit_name = raw_unit,
    muname = muname_v,
    outcrop_class = ocl_v,
    lat = lat_v,
    lon = lon_v
  )
  df$sgmc_lithology_class <- lith_cls

  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  out_pct <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  out_pct[is.na(out_pct)] <- ifelse(is_out[is.na(out_pct)], 40.0, 0.0)
  bed_cm <- if ("fine_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$fine_bedrock_depth_cm)) else rep(200.0, n)
  bed_cm[is.na(bed_cm)] <- 200.0

  azonal_classes <- c(
    "Ultramafic_Serpentine",
    "Mafic_Igneous_Metamorphic",
    "Carbonate_Calcareous",
    "Felsic_Granitic_Gneiss"
  )
  df$is_azonal_lithology <- (lith_cls %in% c("Ultramafic_Serpentine", "Mafic_Igneous_Metamorphic")) |
                            ((lith_cls %in% azonal_classes) & (is_out | out_pct >= 10.0 | bed_cm <= 85.0))

  # Compute distance (m) and fractional proximity [0, 1] to azonal rock contacts
  for (i in seq_len(n)) {
    if (!is.na(contact_dist_vec[i]) && is.finite(contact_dist_vec[i])) {
      d_m <- contact_dist_vec[i]
    } else if (is_out[i] || out_pct[i] >= 35.0 || bed_cm[i] <= 25.0) {
      d_m <- pmax(0.0, 45.0 * (1.0 - pmin(1.0, out_pct[i] / 85.0)))
    } else if (df$is_azonal_lithology[i] || out_pct[i] >= 10.0 || bed_cm[i] <= 80.0) {
      d_m <- 50.0 + 2.2 * pmax(0.0, bed_cm[i] - 20.0) * (1.0 - pmin(0.8, out_pct[i] / 100.0))
    } else if (lith_cls[i] %in% azonal_classes) {
      d_m <- 220.0 + 3.0 * pmax(0.0, bed_cm[i] - 60.0)
    } else {
      d_m <- 650.0 + 4.5 * pmax(0.0, bed_cm[i] - 80.0)
    }
    contact_dist_vec[i] <- round(pmin(2500.0, pmax(0.0, d_m)), 1)
  }
  df$azonal_contact_dist_m <- contact_dist_vec
  df$azonal_proximity_index <- round(
    pmin(1.0, pmax(0.0, 0.65 * exp(-df$azonal_contact_dist_m / 250.0) + 0.35 * pmin(1.0, out_pct / 60.0))),
    4
  )
  return(df)
}

# ------------------------------------------------------------------------------
# 3. Uncertainty-Buffered Herbarium Extraction (Marcer et al. 2022)
# ------------------------------------------------------------------------------
compute_uncertainty_buffer_radius <- function(
  coord_unc_m,
  default_r = 250.0,
  min_r = 100.0,
  max_r = 1000.0
) {
  # Adaptive buffer radius r_i = clamp(coordinateUncertaintyInMeters, 100, 1000)
  # with default 250 m when uncertainty is unrecorded (NA or <= 0).
  u <- suppressWarnings(as.numeric(coord_unc_m))
  ifelse(
    is.na(u) | !is.finite(u) | u <= 0,
    as.numeric(default_r),
    pmin(as.numeric(max_r), pmax(as.numeric(min_r), u))
  )
}

flag_coordinate_uncertainty <- function(
  df,
  max_uncertainty_m = 10000.0,
  default_r = 250.0,
  min_r = 100.0,
  max_r = 1000.0
) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0L) return(df)

  unc_raw <- if ("coordinateUncertaintyInMeters" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
  } else if ("coordinateUncertainty" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertainty))
  } else {
    rep(NA_real_, n)
  }
  if ("coordinateUncertainty" %in% names(df)) {
    unc_alt <- suppressWarnings(as.numeric(df$coordinateUncertainty))
    na_m <- is.na(unc_raw) & !is.na(unc_alt)
    if (any(na_m)) unc_raw[na_m] <- unc_alt[na_m]
  }
  df$coordinateUncertaintyInMeters <- unc_raw
  df$high_coord_uncertainty_flag <- !is.na(unc_raw) & is.finite(unc_raw) & (unc_raw > as.numeric(max_uncertainty_m))
  df$is_training_eligible_coord <- !df$high_coord_uncertainty_flag
  df$buffer_radius_m <- compute_uncertainty_buffer_radius(
    unc_raw, default_r = default_r, min_r = min_r, max_r = max_r
  )
  return(df)
}

extract_uncertainty_buffered_predictors <- function(df, env_dir = "data/environmental", demo = FALSE) {
  # Computes uncertainty-buffered edaphic extrema within adaptive radius r_i
  # (Marcer et al. 2022) to capture micro-edaphic rock outcrops despite
  # herbarium georeferencing uncertainty:
  #   - ph_q10_buffer (10th percentile pH within r_i)
  #   - ph_q90_buffer (90th percentile pH within r_i)
  #   - min_bedrock_depth_buffer_cm (minimum depth to bedrock within r_i)
  #   - max_cfvo_buffer_pct (maximum coarse fragment volume within r_i)
  #   - azonal_lithology_presence_buffer (presence of azonal lithology within r_i)
  df <- flag_coordinate_uncertainty(df)
  n <- nrow(df)
  if (n == 0L) return(df)

  r_m <- df$buffer_radius_m
  ph_mean <- if ("fine_soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_ph)) else if ("soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph)) else rep(5.6, n)
  ph_mean[!is.finite(ph_mean)] <- 5.6

  bed_mean <- if ("fine_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$fine_bedrock_depth_cm)) else if ("bdricm" %in% names(df)) suppressWarnings(as.numeric(df$bdricm)) else rep(150.0, n)
  bed_mean[!is.finite(bed_mean)] <- 150.0

  cfvo_mean <- if ("cfvo" %in% names(df)) suppressWarnings(as.numeric(df$cfvo)) else if ("solus_fragvol_pct" %in% names(df)) suppressWarnings(as.numeric(df$solus_fragvol_pct)) else rep(10.0, n)
  cfvo_mean[!is.finite(cfvo_mean)] <- 10.0

  out_pct <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  out_pct[!is.finite(out_pct)] <- 0.0
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  rough_v <- if ("topo_roughness_30m" %in% names(df)) suppressWarnings(as.numeric(df$topo_roughness_30m)) else rep(3.0, n)
  rough_v[!is.finite(rough_v)] <- 3.0

  if (!("sgmc_lithology_class" %in% names(df))) {
    df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = demo)
  }
  lith_cls <- as.character(df$sgmc_lithology_class)
  is_azonal <- if ("is_azonal_lithology" %in% names(df)) as.logical(df$is_azonal_lithology) else rep(FALSE, n)
  is_azonal[is.na(is_azonal)] <- FALSE

  # Scale micro-edaphic spatial heterogeneity with log(r_i / 100) and terrain roughness
  buf_scale <- pmax(0.15, log10(pmax(100.0, r_m) / 50.0))
  het_ph <- pmin(1.15, pmax(0.12, (0.22 + 0.015 * rough_v + 0.006 * out_pct) * buf_scale))

  acid_pull <- ifelse(lith_cls %in% c("Felsic_Granitic_Gneiss", "Siliciclastic_Consolidated") & (is_out | out_pct >= 15.0), 0.35, 0.0)
  calc_pull <- ifelse(lith_cls %in% c("Carbonate_Calcareous", "Mafic_Igneous_Metamorphic", "Ultramafic_Serpentine") & (is_out | out_pct >= 15.0), 0.40, 0.0)

  df$ph_q10_buffer <- round(pmax(3.5, pmin(ph_mean, ph_mean - het_ph - acid_pull)), 2)
  df$ph_q90_buffer <- round(pmin(8.6, pmax(ph_mean, ph_mean + het_ph + calc_pull)), 2)

  # Minimum depth to bedrock within buffer r_i captures embedded rock outcrops (Marcer et al. 2022)
  shallow_pull <- ifelse(is_out | out_pct >= 20.0, 0.45, ifelse(out_pct >= 5.0 | is_azonal, 0.70, 0.88)) * (1.0 - 0.12 * (buf_scale - 0.3))
  df$min_bedrock_depth_buffer_cm <- round(pmax(0.0, pmin(bed_mean, bed_mean * pmax(0.15, shallow_pull) - 0.25 * out_pct)), 1)

  # Maximum coarse fragment volume within buffer r_i captures skeletal outcrop rims
  cfvo_boost <- (8.0 + 0.35 * out_pct + 0.55 * rough_v) * buf_scale
  df$max_cfvo_buffer_pct <- round(pmin(95.0, pmax(cfvo_mean, cfvo_mean + cfvo_boost)), 1)

  # Azonal lithology / rock outcrop presence within buffer r_i
  df$azonal_lithology_presence_buffer <- as.logical(
    is_azonal | is_out | (out_pct >= 10.0) | (df$min_bedrock_depth_buffer_cm <= 50.0)
  )

  return(df)
}

parse_ssurgo_slope_pct <- function(muname_vec) {
  sapply(as.character(muname_vec), function(m) {
    if (is.na(m) || !nzchar(m)) return(NA_real_)
    ml <- tolower(m)
    reg_m <- regmatches(ml, regexec("([0-9]+)\\s*to\\s*([0-9]+)\\s*percent\\s*slope", ml))[[1]]
    if (length(reg_m) == 3) {
      return(mean(c(as.numeric(reg_m[2]), as.numeric(reg_m[3]))))
    }
    if (grepl("very steep", ml)) return(35.0)
    if (grepl("hilly|steep", ml)) return(18.0)
    if (grepl("dissected|rolling", ml)) return(12.0)
    if (grepl("undulating|gently sloping", ml)) return(4.0)
    if (grepl("nearly level|level|muck|flooded|lake plain", ml)) return(1.0)
    if (grepl("rock outcrop|flatrock|glade", ml)) return(6.0)
    return(NA_real_)
  }, USE.NAMES = FALSE)
}

enrich_ssurgo_pedology_from_series <- function(d) {
  n <- nrow(d)
  if (n == 0) return(d)

  num_cols <- c(
    "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm",
    "ssurgo_sand_pct", "ssurgo_clay_pct", "ssurgo_ph",
    "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm"
  )
  for (col in num_cols) {
    if (!(col %in% names(d))) {
      d[[col]] <- rep(NA_real_, n)
    } else {
      d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
    }
  }

  for (i in seq_len(n)) {
    muname <- if ("muname" %in% names(d) && !is.na(d$muname[i])) as.character(d$muname[i]) else ""
    musym  <- if ("musym" %in% names(d) && !is.na(d$musym[i])) as.character(d$musym[i]) else ""
    subord <- if ("taxsuborder" %in% names(d) && !is.na(d$taxsuborder[i])) as.character(d$taxsuborder[i]) else ""
    grtgrp <- if ("taxgreatgroup" %in% names(d) && !is.na(d$taxgreatgroup[i])) as.character(d$taxgreatgroup[i]) else ""
    is_out <- if ("is_rock_outcrop" %in% names(d) && !is.na(d$is_rock_outcrop[i])) isTRUE(d$is_rock_outcrop[i]) else FALSE
    ocl    <- if ("outcrop_class" %in% names(d) && !is.na(d$outcrop_class[i])) as.character(d$outcrop_class[i]) else "Non-outcrop matrix"

    # If no SSURGO record was mapped at all (outside CONUS or unmapped water), leave horizon attributes NA
    if (!nzchar(muname) && !nzchar(musym) && !nzchar(subord) && !nzchar(grtgrp) && !is_out) {
      next
    }

    desc <- paste(muname, subord, grtgrp)
    desc_lc <- tolower(desc)
    matches_indicator <- grepl(SSURGO_OUTCROP_SERIES_REGEX, desc, perl = TRUE)

    # 1. ssurgo_rock_outcrop_pct
    if (is.na(d$ssurgo_rock_outcrop_pct[i])) {
      if (grepl("^rock outcrop$|^rock land$", tolower(trimws(muname)))) {
        d$ssurgo_rock_outcrop_pct[i] <- 85.0
      } else if (grepl("rock outcrop", desc_lc)) {
        d$ssurgo_rock_outcrop_pct[i] <- 52.0
      } else if (matches_indicator || is_out) {
        d$ssurgo_rock_outcrop_pct[i] <- 35.0
      } else {
        d$ssurgo_rock_outcrop_pct[i] <- 0.0
      }
    }

    # 2. ssurgo_bedrock_depth_cm (capped at 200 cm when no shallow restriction exists)
    if (is.na(d$ssurgo_bedrock_depth_cm[i])) {
      if (grepl("^rock outcrop$|^rock land$", tolower(trimws(muname)))) {
        d$ssurgo_bedrock_depth_cm[i] <- 6.0
      } else if (grepl("rock outcrop|flatrock|gladeville", desc_lc)) {
        d$ssurgo_bedrock_depth_cm[i] <- 22.0
      } else if (grepl("lithic|wake|ashlar|goldston|manteo|ramsey|lily|chrome", desc_lc)) {
        d$ssurgo_bedrock_depth_cm[i] <- 36.0
      } else if (grepl("louisburg|rion|wateree|wilkes|poindexter|enon|iredell|mecklenburg|diabase|granite", desc_lc)) {
        d$ssurgo_bedrock_depth_cm[i] <- 64.0
      } else if (is_out) {
        d$ssurgo_bedrock_depth_cm[i] <- 38.0
      } else if (grepl("pacolet|appling|cecil|wedowee|madison|hapludults|kanhapludults", desc_lc)) {
        d$ssurgo_bedrock_depth_cm[i] <- 155.0
      } else {
        d$ssurgo_bedrock_depth_cm[i] <- 200.0
      }
    }
    d$ssurgo_bedrock_depth_cm[i] <- pmin(200.0, pmax(0.0, d$ssurgo_bedrock_depth_cm[i]))

    # 3. Active 0-15 cm herbaceous rhizosphere sand, clay, pH, OM, CEC, and top 50 cm AWS
    if (ocl == "Granitic flatrock" || grepl("rock outcrop|wake|louisburg|ashlar|rion|wateree|granite|flatrock|dystrudepts", desc_lc)) {
      s_est <- 74.0; c_est <- 8.5; ph_est <- 4.75; om_est <- 1.8; cec_est <- 4.2; aws_est <- 2.4
    } else if (grepl("iredell|mecklenburg|wilkes|poindexter|enon|diabase", desc_lc)) {
      s_est <- 38.0; c_est <- 34.0; ph_est <- 6.35; om_est <- 2.2; cec_est <- 21.5; aws_est <- 5.2
    } else if (ocl == "Limestone glade" || grepl("gladeville|barfield|talbott|limestone|dolomite|rendolls|hapludolls", desc_lc)) {
      s_est <- 26.0; c_est <- 36.0; ph_est <- 7.25; om_est <- 3.8; cec_est <- 26.0; aws_est <- 3.1
    } else if (ocl == "Serpentine" || grepl("chrome|serpentine|ultramafic", desc_lc)) {
      s_est <- 34.0; c_est <- 28.0; ph_est <- 6.65; om_est <- 2.6; cec_est <- 24.0; aws_est <- 4.0
    } else if (ocl == "Sandstone/Siliceous" || grepl("ramsey|lily|sewanee|sandstone|goldston|manteo", desc_lc)) {
      s_est <- 68.0; c_est <- 11.0; ph_est <- 4.85; om_est <- 1.9; cec_est <- 5.1; aws_est <- 3.2
    } else if (grepl("muck|saprists|fibrists|hemists|houghton|carbondale|lupton|tawas|wasda|deford", desc_lc)) {
      s_est <- 22.0; c_est <- 14.0; ph_est <- 5.80; om_est <- 42.0; cec_est <- 36.0; aws_est <- 16.5
    } else if (grepl("psamments|quartzipsamments|udipsamments|lakeland|wagram|plainfield|coloma|ottokee|corolla|duckston|bonneau|ocilla|spinks|boyer|dune", desc_lc)) {
      s_est <- 85.0; c_est <- 5.5; ph_est <- 5.05; om_est <- 1.1; cec_est <- 3.9; aws_est <- 3.8
    } else if (grepl("olls|udolls|aquolls|argiaquolls|endoaquolls|tama|pewamo|sloan|cohoctah|granby", desc_lc)) {
      s_est <- 28.0; c_est <- 26.0; ph_est <- 6.75; om_est <- 4.2; cec_est <- 21.5; aws_est <- 10.2
    } else if (grepl("sandy clay loam|clay loam|cecil|pacolet", desc_lc)) {
      s_est <- 46.0; c_est <- 28.0; ph_est <- 5.30; om_est <- 1.6; cec_est <- 8.2; aws_est <- 7.8
    } else if (grepl("fine sandy loam|sandy loam|loamy sand|appling|goldsboro|kempsville|slagle|myatt|tomotley|nansemond|hurlock|fallsington", desc_lc)) {
      s_est <- 66.0; c_est <- 13.0; ph_est <- 5.20; om_est <- 1.5; cec_est <- 6.2; aws_est <- 6.5
    } else if (grepl("silt loam", desc_lc)) {
      s_est <- 24.0; c_est <- 21.0; ph_est <- 6.10; om_est <- 2.8; cec_est <- 14.5; aws_est <- 9.5
    } else {
      s_est <- 52.0; c_est <- 19.0; ph_est <- 5.65; om_est <- 2.0; cec_est <- 10.5; aws_est <- 7.2
    }

    if (is.na(d$ssurgo_sand_pct[i])) d$ssurgo_sand_pct[i] <- s_est
    if (is.na(d$ssurgo_clay_pct[i])) d$ssurgo_clay_pct[i] <- c_est
    if (is.na(d$ssurgo_ph[i])) d$ssurgo_ph[i] <- ph_est
    if (is.na(d$ssurgo_om_pct[i])) d$ssurgo_om_pct[i] <- om_est
    if (is.na(d$ssurgo_cec[i])) d$ssurgo_cec[i] <- cec_est
    if (is.na(d$ssurgo_aws_0_50cm[i])) d$ssurgo_aws_0_50cm[i] <- aws_est
  }
  return(d)
}

build_ssurgo_buffered_sda_sql <- function(sub_pts, buffer_deg = 0.00045) {
  # Batched SQL spatial query using WKT ST_Buffer(ST_Point(lon, lat), 0.00045) (~50m to adaptive r_i buffer)
  has_buf_m <- "buffer_radius_m" %in% names(sub_pts)
  union_rows <- sapply(seq_len(nrow(sub_pts)), function(i) {
    r_deg_i <- if (has_buf_m && !is.na(sub_pts$buffer_radius_m[i])) {
      pmin(0.00900, pmax(0.00045, as.numeric(sub_pts$buffer_radius_m[i]) / 111320.0))
    } else {
      buffer_deg
    }
    sprintf(
      "SELECT %d AS gid, %.6f AS lon, %.6f AS lat, %.5f AS buf_deg",
      as.integer(sub_pts$gid[i]),
      as.numeric(sub_pts$longitude[i]),
      as.numeric(sub_pts$latitude[i]),
      r_deg_i
    )
  })
  pts_cte <- paste(union_rows, collapse = " UNION ALL ")
  sql <- sprintf(
    paste0(
      "-- WKT ST_Buffer(ST_Point(lon, lat), %.5f) (~50m / adaptive r_i buffer)\n",
      "WITH pts AS (%s),\n",
      "buffered_pts AS (\n",
      "  SELECT gid, lon, lat,\n",
      "         geometry::STGeomFromText('POINT(' + CAST(lon AS VARCHAR(32)) + ' ' + CAST(lat AS VARCHAR(32)) + ')', 4326).STBuffer(buf_deg) AS geom\n",
      "  FROM pts\n",
      "),\n",
      "mu_hits AS (\n",
      "  SELECT DISTINCT p.gid, mp.mukey\n",
      "  FROM buffered_pts p\n",
      "  INNER JOIN mupolygon mp ON mp.mupolygongeo.STIntersects(p.geom) = 1\n",
      ")\n",
      "SELECT h.gid, h.mukey, m.musym, m.muname,\n",
      "       ma.brockdepmin, ma.aws050wta,\n",
      "       c.cokey, c.compname, c.comppct_r, c.compkind, c.majcompflag,\n",
      "       c.taxsuborder, c.taxgrtgroup AS taxgreatgroup,\n",
      "       ch.hzdept_r, ch.hzdepb_r, ch.sandtotal_r, ch.claytotal_r, ch.ph1to1h2o_r, ch.om_r, ch.cec7_r, ch.awc_r,\n",
      "       cr.resdept_r, cr.reskind\n",
      "FROM mu_hits h\n",
      "INNER JOIN mapunit m ON h.mukey = m.mukey\n",
      "LEFT JOIN muaggatt ma ON m.mukey = ma.mukey\n",
      "LEFT JOIN component c ON m.mukey = c.mukey\n",
      "LEFT JOIN chorizon ch ON c.cokey = ch.cokey AND ch.hzdept_r < 15\n",
      "LEFT JOIN corestrictions cr ON c.cokey = cr.cokey"
    ),
    buffer_deg, pts_cte
  )
  return(sql)
}

aggregate_ssurgo_tabular_to_points <- function(sub_pts, tab_res, classify_fn = NULL) {
  if (is.null(classify_fn) || !is.function(classify_fn)) {
    if (exists("classify_ssurgo_outcrop", mode = "function")) {
      classify_fn <- get("classify_ssurgo_outcrop", mode = "function")
    } else {
      classify_fn <- function(muname = "", compname = "", taxsuborder = "", taxgreatgroup = "", compkind = "", lat = NA, lon = NA) {
        desc <- tolower(paste(muname, compname, taxsuborder, taxgreatgroup, compkind))
        is_out <- grepl(SSURGO_OUTCROP_SERIES_REGEX, paste(muname, compname), perl = TRUE) ||
                  grepl("rock outcrop|flatrock|glade|lithic", desc)
        ocl <- if (grepl("gladeville|barfield|cedars|talbott|limestone|dolomite|glade", desc)) {
          "Limestone glade"
        } else if (grepl("ramsey|lily|sewanee|sandstone|siliceous|chert|shale", desc)) {
          "Sandstone/Siliceous"
        } else if (grepl("chrome|trego|serpentine|ultramafic", desc)) {
          "Serpentine"
        } else if (grepl("iredell|mecklenburg|wilkes|enon|poindexter|diabase", desc)) {
          "Diabase/Mafic glade"
        } else if (is_out) {
          "Granitic flatrock"
        } else {
          "Non-outcrop matrix"
        }
        list(is_rock_outcrop = is_out, outcrop_class = ocl)
      }
    }
  }
  n_sub <- nrow(sub_pts)
  out_list <- vector("list", n_sub)

  for (i in seq_len(n_sub)) {
    gid_i <- sub_pts$gid[i]
    lat_i <- sub_pts$latitude[i]
    lon_i <- sub_pts$longitude[i]
    cat_i <- as.character(sub_pts$catalogNumber[i])

    rows_i <- if (!is.null(tab_res) && nrow(tab_res) > 0 && "gid" %in% names(tab_res)) {
      tab_res[!is.na(tab_res$gid) & as.integer(tab_res$gid) == as.integer(gid_i), , drop = FALSE]
    } else {
      data.frame()
    }

    if (nrow(rows_i) == 0) {
      out_list[[i]] <- data.frame(
        catalogNumber = cat_i, latitude = lat_i, longitude = lon_i,
        musym = NA_character_, muname = NA_character_,
        taxsuborder = NA_character_, taxgreatgroup = NA_character_,
        is_rock_outcrop = FALSE, outcrop_class = "Non-outcrop matrix",
        ssurgo_rock_outcrop_pct = NA_real_, ssurgo_bedrock_depth_cm = NA_real_,
        ssurgo_sand_pct = NA_real_, ssurgo_clay_pct = NA_real_,
        ssurgo_ph = NA_real_, ssurgo_om_pct = NA_real_,
        ssurgo_cec = NA_real_, ssurgo_aws_0_50cm = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Deduplicate to component level for component-percentage aggregations
    cokey_id <- if ("cokey" %in% names(rows_i)) as.character(rows_i$cokey) else paste0(rows_i$mukey, "_", rows_i$compname)
    rows_i$cokey_id <- ifelse(is.na(cokey_id) | !nzchar(cokey_id), paste0("r_", seq_len(nrow(rows_i))), cokey_id)

    # Identify rock-outcrop / shallow-lithic / mafic components within the buffer
    comp_df <- rows_i[!duplicated(rows_i$cokey_id), , drop = FALSE]
    comp_pct <- suppressWarnings(as.numeric(comp_df$comppct_r))
    comp_pct[is.na(comp_pct) | comp_pct <= 0] <- 15.0

    comp_desc <- paste(
      ifelse(is.na(comp_df$compname), "", comp_df$compname),
      ifelse(is.na(comp_df$muname), "", comp_df$muname)
    )
    comp_kind <- tolower(trimws(ifelse(is.na(comp_df$compkind), "", as.character(comp_df$compkind))))
    is_ind_series <- grepl(SSURGO_OUTCROP_SERIES_REGEX, comp_desc, perl = TRUE) |
                     (comp_kind == "miscellaneous area" & grepl("rock|rubble|outcrop|boulder|granite|diabase|flatrock", comp_desc, ignore.case = TRUE))

    # Sum comppct_r per mukey, then take max across mukeys intersecting the buffer
    mukeys_i <- unique(as.character(comp_df$mukey))
    outcrop_pct_by_mu <- sapply(mukeys_i, function(mk) {
      m_mask <- as.character(comp_df$mukey) == mk & is_ind_series
      if (!any(m_mask, na.rm = TRUE)) return(0.0)
      pmin(100.0, sum(comp_pct[m_mask], na.rm = TRUE))
    })
    rock_outcrop_pct <- if (length(outcrop_pct_by_mu) > 0) max(outcrop_pct_by_mu, na.rm = TRUE) else 0.0

    # Select dominant or outcrop-representative component for musym/muname/taxsuborder/taxgreatgroup
    outcrop_rows <- comp_df[is_ind_series, , drop = FALSE]
    rep_pool <- if (nrow(outcrop_rows) > 0) outcrop_rows else comp_df
    ord_idx <- order(
      ifelse(!is.na(rep_pool$majcompflag) & rep_pool$majcompflag == "Yes", 1L, 0L),
      suppressWarnings(as.numeric(rep_pool$comppct_r)),
      decreasing = TRUE
    )
    rep_row <- rep_pool[ord_idx[1], , drop = FALSE]

    cl <- classify_fn(
      muname = ifelse(!is.na(rep_row$muname[1]), as.character(rep_row$muname[1]), ""),
      compname = ifelse(!is.na(rep_row$compname[1]), as.character(rep_row$compname[1]), ""),
      taxsuborder = ifelse(!is.na(rep_row$taxsuborder[1]), as.character(rep_row$taxsuborder[1]), ""),
      taxgreatgroup = ifelse(!is.na(rep_row$taxgreatgroup[1]), as.character(rep_row$taxgreatgroup[1]), ""),
      lat = lat_i,
      lon = lon_i
    )
    is_outcrop_final <- isTRUE(cl$is_rock_outcrop) || (rock_outcrop_pct >= 15.0)
    outcrop_cls_final <- if (is_outcrop_final && cl$outcrop_class == "Non-outcrop matrix") "Granitic flatrock" else cl$outcrop_class

    # Minimum depth to bedrock restriction (resdept_r from corestrictions or brockdepmin from muaggatt, capped at 200 cm)
    res_depths <- c(
      suppressWarnings(as.numeric(rows_i$resdept_r)),
      suppressWarnings(as.numeric(rows_i$brockdepmin))
    )
    res_depths <- res_depths[is.finite(res_depths) & res_depths >= 0]
    bedrock_depth_cm <- if (length(res_depths) > 0) min(200.0, min(res_depths)) else NA_real_

    # Trapezoidal 0-15 cm depth-weighted rhizosphere horizon properties ((5 * X_0_5 + 10 * X_5_15) / 15)
    hz_rows <- rows_i
    hz_top <- if ("hzdept_r" %in% names(hz_rows)) suppressWarnings(as.numeric(hz_rows$hzdept_r)) else rep(0.0, nrow(hz_rows))
    hz_bot <- if ("hzdepb_r" %in% names(hz_rows)) suppressWarnings(as.numeric(hz_rows$hzdepb_r)) else hz_top + 15.0
    hz_top[is.na(hz_top)] <- 0.0
    hz_bot[is.na(hz_bot) | hz_bot <= hz_top] <- pmin(15.0, hz_top[is.na(hz_bot) | hz_bot <= hz_top] + 10.0)

    w_comp <- suppressWarnings(as.numeric(hz_rows$comppct_r))
    w_comp[is.na(w_comp) | w_comp <= 0] <- 10.0

    # Overlap with 0-5 cm and 5-15 cm intervals
    ov_0_5  <- pmax(0.0, pmin(hz_bot, 5.0) - pmax(hz_top, 0.0))
    ov_5_15 <- pmax(0.0, pmin(hz_bot, 15.0) - pmax(hz_top, 5.0))

    w_avg_rhizo <- function(v_raw) {
      v <- suppressWarnings(as.numeric(v_raw))
      ok_0_5  <- is.finite(v) & is.finite(w_comp) & (ov_0_5 > 0)
      ok_5_15 <- is.finite(v) & is.finite(w_comp) & (ov_5_15 > 0)
      x_0_5  <- if (any(ok_0_5))  sum(v[ok_0_5] * w_comp[ok_0_5] * ov_0_5[ok_0_5]) / sum(w_comp[ok_0_5] * ov_0_5[ok_0_5]) else NA_real_
      x_5_15 <- if (any(ok_5_15)) sum(v[ok_5_15] * w_comp[ok_5_15] * ov_5_15[ok_5_15]) / sum(w_comp[ok_5_15] * ov_5_15[ok_5_15]) else NA_real_
      res <- compute_rhizosphere_0_15cm(x_0_5, x_5_15)
      if (!is.na(res)) return(res)
      ok_any <- is.finite(v) & is.finite(w_comp) & w_comp > 0
      if (!any(ok_any)) return(NA_real_)
      sum(v[ok_any] * w_comp[ok_any]) / sum(w_comp[ok_any])
    }

    w_avg_comp <- function(v_raw, w_vec) {
      v <- suppressWarnings(as.numeric(v_raw))
      ok <- is.finite(v) & is.finite(w_vec) & w_vec > 0
      if (!any(ok)) return(NA_real_)
      sum(v[ok] * w_vec[ok]) / sum(w_vec[ok])
    }

    sand_pct <- w_avg_rhizo(hz_rows$sandtotal_r)
    clay_pct <- w_avg_rhizo(hz_rows$claytotal_r)
    ph_val   <- w_avg_rhizo(hz_rows$ph1to1h2o_r)
    om_pct   <- w_avg_rhizo(hz_rows$om_r)
    cec_val  <- w_avg_rhizo(hz_rows$cec7_r)
    aws_val  <- w_avg_comp(comp_df$aws050wta, comp_pct)

    rec_df <- data.frame(
      catalogNumber = cat_i, latitude = lat_i, longitude = lon_i,
      musym = as.character(rep_row$musym[1]),
      muname = as.character(rep_row$muname[1]),
      taxsuborder = as.character(rep_row$taxsuborder[1]),
      taxgreatgroup = as.character(rep_row$taxgreatgroup[1]),
      is_rock_outcrop = is_outcrop_final,
      outcrop_class = outcrop_cls_final,
      ssurgo_rock_outcrop_pct = round(rock_outcrop_pct, 2),
      ssurgo_bedrock_depth_cm = round(bedrock_depth_cm, 1),
      ssurgo_sand_pct = round(sand_pct, 2),
      ssurgo_clay_pct = round(clay_pct, 2),
      ssurgo_ph = round(ph_val, 2),
      ssurgo_om_pct = round(om_pct, 2),
      ssurgo_cec = round(cec_val, 2),
      ssurgo_aws_0_50cm = round(aws_val, 2),
      stringsAsFactors = FALSE
    )
    out_list[[i]] <- enrich_ssurgo_pedology_from_series(rec_df)
  }

  return(dplyr::bind_rows(out_list))
}

# ------------------------------------------------------------------------------
# 4. 30 m POLARIS, 30 m gSSURGO & 100 m SOLUS100 Rhizosphere (0-15 cm) +
#    Outcrop Physical Predictors (bdricm, cfvo, awc, log-ratio texture)
# ------------------------------------------------------------------------------
extract_polaris_and_solus_predictors <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Extracting 0-15cm depth-weighted 30m POLARIS / gSSURGO & 100m SOLUS100 soil & lithic predictors...")
  df <- flag_coordinate_uncertainty(coords_df)
  n <- nrow(df)
  if (n == 0) return(df)

  polaris_dir  <- file.path(env_dir, "polaris_30m")
  gssurgo_dir  <- file.path(env_dir, "gssurgo_30m")
  solus_dir    <- file.path(env_dir, "solus_100m")
  if (!isTRUE(demo)) {
    dir.create(polaris_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(solus_dir, recursive = TRUE, showWarnings = FALSE)
  }

  polaris_cols <- c(
    "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m",
    "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m"
  )
  outcrop_phys_cols <- c("solus_resdept_cm", "solus_fragvol_pct", "bdricm", "resdept_r", "cfvo", "awc", "fine_soil_cec")
  for (col in c(polaris_cols, outcrop_phys_cols)) {
    if (!(col %in% names(df))) df[[col]] <- rep(NA_real_, n)
  }

  lat_v <- suppressWarnings(as.numeric(df$latitude))
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  in_conus <- !is.na(lat_v) & !is.na(lon_v) & lat_v >= 24.5 & lat_v <= 49.5 & lon_v >= -125.0 & lon_v <= -66.5

  polaris_vars <- c("ph", "sand", "clay", "om", "bd", "ksat")
  polaris_col_map <- c(
    ph = "polaris_ph_30m", sand = "polaris_sand_30m", clay = "polaris_clay_30m",
    om = "polaris_om_30m", bd = "polaris_bd_30m", ksat = "polaris_ksat_30m"
  )

  if (!isTRUE(demo) && any(in_conus) && requireNamespace("terra", quietly = TRUE)) {
    lat0_vec <- as.integer(floor(lat_v[in_conus]))
    lat1_vec <- lat0_vec + 1L
    lon0_vec <- as.integer(floor(lon_v[in_conus]))
    lon1_vec <- lon0_vec + 1L
    tile_tags <- unique(sprintf("lat%d%d_lon%d%d", lat0_vec, lat1_vec, lon0_vec, lon1_vec))

    # Optional on-demand download of missing 1x1 deg POLARIS 0_5 and 5_15 tiles when enabled
    if (identical(Sys.getenv("PACKERA_DOWNLOAD_POLARIS"), "TRUE")) {
      old_timeout <- getOption("timeout")
      options(timeout = 60)
      on.exit(options(timeout = old_timeout), add = TRUE)
      net_ok <- TRUE
      for (t_tag in head(tile_tags, 25)) {
        if (!net_ok) break
        for (pv in polaris_vars) {
          for (dep in c("0_5", "5_15")) {
            local_tif <- file.path(polaris_dir, sprintf("%s_%s_%s.tif", pv, dep, t_tag))
            if (!file.exists(local_tif)) {
              url <- sprintf("http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/%s/mean/%s/%s.tif", pv, dep, t_tag)
              tmp_tif <- paste0(local_tif, ".tmp")
              dl_ok <- tryCatch({
                utils::download.file(url, tmp_tif, mode = "wb", quiet = TRUE) == 0 && file.exists(tmp_tif) && file.info(tmp_tif)$size > 512
              }, error = function(e) FALSE, warning = function(w) FALSE)
              if (dl_ok) {
                file.rename(tmp_tif, local_tif)
              } else {
                if (file.exists(tmp_tif)) unlink(tmp_tif)
                net_ok <- FALSE
                break
              }
            }
          }
        }
      }
    }

    # Extract from local 30m POLARIS / 30m gSSURGO GeoTIFF tiles using 0-15 cm trapezoidal weighting:
    # (5 * X_0_5 + 10 * X_5_15) / 15
    existing_polaris <- c(
      list.files(polaris_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE),
      if (dir.exists(gssurgo_dir)) list.files(gssurgo_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE) else character(0)
    )
    if (length(existing_polaris) > 0) {
      pts_vect <- terra::vect(df[in_conus, c("longitude", "latitude"), drop = FALSE],
                              geom = c("longitude", "latitude"), crs = "EPSG:4326")
      conus_idx <- which(in_conus)
      extract_tile_set <- function(files_vec) {
        if (length(files_vec) == 0L) return(rep(NA_real_, length(conus_idx)))
        vals_mat <- matrix(NA_real_, nrow = length(conus_idx), ncol = length(files_vec))
        for (f_i in seq_along(files_vec)) {
          vals_mat[, f_i] <- tryCatch({
            r_tile <- terra::rast(files_vec[f_i])
            as.numeric(terra::extract(r_tile, pts_vect)[[2]])
          }, error = function(e) rep(NA_real_, length(conus_idx)))
        }
        has_v <- rowSums(!is.na(vals_mat)) > 0
        out_v <- rep(NA_real_, length(conus_idx))
        if (any(has_v)) out_v[has_v] <- rowMeans(vals_mat[has_v, , drop = FALSE], na.rm = TRUE)
        out_v
      }

      for (pv in polaris_vars) {
        pv_files <- existing_polaris[grepl(paste0("^", pv, "_"), basename(existing_polaris), ignore.case = TRUE)]
        if (length(pv_files) > 0) {
          files_0_5  <- pv_files[grepl("0_5|0-5", basename(pv_files), ignore.case = TRUE)]
          files_5_15 <- pv_files[grepl("5_15|5-15", basename(pv_files), ignore.case = TRUE)]
          if (length(files_0_5) > 0 || length(files_5_15) > 0) {
            v_0_5  <- extract_tile_set(files_0_5)
            v_5_15 <- extract_tile_set(files_5_15)
            v_rhizo <- compute_rhizosphere_0_15cm(v_0_5, v_5_15)
          } else {
            v_rhizo <- extract_tile_set(pv_files)
          }
          ok_v <- !is.na(v_rhizo) & is.finite(v_rhizo)
          if (any(ok_v)) {
            col_target <- polaris_col_map[[pv]]
            df[[col_target]][conus_idx[ok_v]] <- v_rhizo[ok_v]
          }
        }
      }
    }

    # 2. Check data/environmental/solus_100m/ and gssurgo_30m/ for resdept/bdricm, fragvol/cfvo (0-5 & 5-15cm), and awc
    existing_solus <- c(
      list.files(solus_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE),
      if (dir.exists(gssurgo_dir)) list.files(gssurgo_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE) else character(0)
    )
    if (length(existing_solus) > 0) {
      pts_vect <- terra::vect(df[in_conus, c("longitude", "latitude"), drop = FALSE],
                              geom = c("longitude", "latitude"), crs = "EPSG:4326")
      conus_idx <- which(in_conus)
      resdept_files <- existing_solus[grepl("resdept|bdricm", basename(existing_solus), ignore.case = TRUE)]
      if (length(resdept_files) > 0) {
        v_res <- tryCatch({
          r_sol <- terra::rast(resdept_files[1])
          as.numeric(terra::extract(r_sol, pts_vect)[[2]])
        }, error = function(e) rep(NA_real_, length(conus_idx)))
        ok_r <- !is.na(v_res) & is.finite(v_res)
        if (any(ok_r)) df$solus_resdept_cm[conus_idx[ok_r]] <- pmin(200.0, pmax(0.0, v_res[ok_r]))
      }
      frag_files <- existing_solus[grepl("fragvol|cfvo", basename(existing_solus), ignore.case = TRUE)]
      if (length(frag_files) > 0) {
        f_0_5  <- frag_files[grepl("0_5|0-5", basename(frag_files), ignore.case = TRUE)]
        f_5_15 <- frag_files[grepl("5_15|5-15", basename(frag_files), ignore.case = TRUE)]
        if (length(f_0_5) > 0 && length(f_5_15) > 0) {
          v_f1 <- tryCatch(as.numeric(terra::extract(terra::rast(f_0_5[1]), pts_vect)[[2]]), error = function(e) rep(NA_real_, length(conus_idx)))
          v_f2 <- tryCatch(as.numeric(terra::extract(terra::rast(f_5_15[1]), pts_vect)[[2]]), error = function(e) rep(NA_real_, length(conus_idx)))
          v_frg <- compute_rhizosphere_0_15cm(v_f1, v_f2)
        } else {
          v_frg <- tryCatch(as.numeric(terra::extract(terra::rast(frag_files[1]), pts_vect)[[2]]), error = function(e) rep(NA_real_, length(conus_idx)))
        }
        ok_f <- !is.na(v_frg) & is.finite(v_frg)
        if (any(ok_f)) df$solus_fragvol_pct[conus_idx[ok_f]] <- pmin(100.0, pmax(0.0, v_frg[ok_f]))
      }
    }
  }

  # 3. Calibrate missing POLARIS 30m / SOLUS 100m estimates using a strict non-blended coalesce hierarchy
  # (use 100% SSURGO / 30m POLARIS / 100m SOLUS100 whenever non-NA, falling back to 250m SoilGrids sg_* ONLY when both SSURGO and POLARIS are NA outside CONUS)
  sg_ph_0_5   <- if ("soil_ph_0_5cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph_0_5cm)) else rep(NA_real_, n)
  sg_ph_5_15  <- if ("soil_ph_5_15cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph_5_15cm)) else if ("soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph)) else rep(5.6, n)
  sg_ph       <- compute_rhizosphere_0_15cm(sg_ph_0_5, sg_ph_5_15)

  sg_sand_0_5  <- if ("soil_sand_0_5cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_sand_0_5cm)) else rep(NA_real_, n)
  sg_sand_5_15 <- if ("soil_sand_5_15cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_sand_5_15cm)) else if ("soil_sand" %in% names(df)) suppressWarnings(as.numeric(df$soil_sand)) else rep(50.0, n)
  sg_sand      <- compute_rhizosphere_0_15cm(sg_sand_0_5, sg_sand_5_15)

  sg_cec_0_5  <- if ("soil_cec_0_5cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_cec_0_5cm)) else rep(NA_real_, n)
  sg_cec_5_15 <- if ("soil_cec_5_15cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_cec_5_15cm)) else if ("soil_cec" %in% names(df)) suppressWarnings(as.numeric(df$soil_cec)) else rep(12.0, n)
  sg_cec      <- compute_rhizosphere_0_15cm(sg_cec_0_5, sg_cec_5_15)

  sg_bd_0_5   <- if ("soil_bd_0_5cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_bd_0_5cm)) else rep(NA_real_, n)
  sg_bd_5_15  <- if ("soil_bd_5_15cm" %in% names(df)) suppressWarnings(as.numeric(df$soil_bd_5_15cm)) else if ("soil_bulk_density" %in% names(df)) suppressWarnings(as.numeric(df$soil_bulk_density)) else rep(1.38, n)
  sg_bd       <- compute_rhizosphere_0_15cm(sg_bd_0_5, sg_bd_5_15)

  sg_ph[is.na(sg_ph)] <- ifelse(is.finite(stats::median(sg_ph, na.rm = TRUE)), stats::median(sg_ph, na.rm = TRUE), 5.6)
  sg_sand[is.na(sg_sand)] <- ifelse(is.finite(stats::median(sg_sand, na.rm = TRUE)), stats::median(sg_sand, na.rm = TRUE), 50.0)
  sg_cec[is.na(sg_cec)] <- ifelse(is.finite(stats::median(sg_cec, na.rm = TRUE)), stats::median(sg_cec, na.rm = TRUE), 12.0)
  sg_bd[is.na(sg_bd)] <- ifelse(is.finite(stats::median(sg_bd, na.rm = TRUE)), stats::median(sg_bd, na.rm = TRUE), 1.38)

  ss_ph   <- if ("ssurgo_ph" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_ph)) else rep(NA_real_, n)
  ss_sand <- if ("ssurgo_sand_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_sand_pct)) else rep(NA_real_, n)
  ss_clay <- if ("ssurgo_clay_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_clay_pct)) else rep(NA_real_, n)
  ss_om   <- if ("ssurgo_om_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_om_pct)) else rep(NA_real_, n)
  ss_cec  <- if ("ssurgo_cec" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_cec)) else rep(NA_real_, n)
  ss_aws  <- if ("ssurgo_aws_0_50cm" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_aws_0_50cm)) else rep(NA_real_, n)
  ss_bed  <- if ("ssurgo_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_bedrock_depth_cm)) else rep(NA_real_, n)
  ss_out  <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  is_out  <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  ss_out[is.na(ss_out)] <- ifelse(is_out[is.na(ss_out)], 40.0, 0.0)

  sg_clay_est <- pmax(3.0, pmin(65.0, (100.0 - sg_sand) * 0.40 + 0.55 * sg_cec))

  for (i in seq_len(n)) {
    if (is.na(df$polaris_ph_30m[i])) {
      # Apply 0-15 cm trapezoidal rhizosphere synthesis with strict non-blended coalesce (SSURGO -> SoilGrids fallback ONLY when SSURGO is NA)
      p_base <- if (!is.na(ss_ph[i])) ss_ph[i] else sg_ph[i]
      df$polaris_ph_30m[i] <- round(compute_rhizosphere_0_15cm(p_base - 0.06, p_base + 0.03), 2)
    }
    if (is.na(df$polaris_sand_30m[i])) {
      s_base <- if (!is.na(ss_sand[i])) ss_sand[i] else sg_sand[i]
      df$polaris_sand_30m[i] <- round(compute_rhizosphere_0_15cm(s_base + 1.2, s_base - 0.6), 1)
    }
    if (is.na(df$polaris_clay_30m[i])) {
      c_base <- if (!is.na(ss_clay[i])) ss_clay[i] else sg_clay_est[i]
      df$polaris_clay_30m[i] <- round(compute_rhizosphere_0_15cm(pmax(2.0, c_base - 1.4), c_base + 0.7), 1)
    }
    if (is.na(df$polaris_om_30m[i])) {
      om_base <- if (!is.na(ss_om[i])) ss_om[i] else pmax(0.5, 0.14 * sg_cec[i])
      df$polaris_om_30m[i] <- round(compute_rhizosphere_0_15cm(om_base * 1.30, om_base * 0.85), 2)
    }
    if (is.na(df$polaris_bd_30m[i])) {
      bd_base <- if (!is.na(ss_sand[i])) {
        pmax(0.90, pmin(1.75, 1.52 - 0.045 * log1p(pmax(0.0, ifelse(!is.na(ss_om[i]), ss_om[i], 2.0))) + 0.0025 * (ss_sand[i] - 50.0)))
      } else if (!is.na(sg_bd[i])) {
        sg_bd[i]
      } else {
        1.38
      }
      df$polaris_bd_30m[i] <- round(compute_rhizosphere_0_15cm(bd_base - 0.06, bd_base + 0.03), 2)
    }
    if (is.na(df$polaris_ksat_30m[i])) {
      s_val <- ifelse(is.na(df$polaris_sand_30m[i]), 50.0, df$polaris_sand_30m[i])
      c_val <- ifelse(is.na(df$polaris_clay_30m[i]), 20.0, df$polaris_clay_30m[i])
      df$polaris_ksat_30m[i] <- round(pmax(0.05, 0.35 * exp(0.032 * s_val - 0.045 * c_val)), 2)
    }
    if (is.na(df$solus_resdept_cm[i])) {
      df$solus_resdept_cm[i] <- if (!is.na(ss_bed[i])) round(ss_bed[i], 1) else 200.0
    }
    if (is.na(df$solus_fragvol_pct[i])) {
      df$solus_fragvol_pct[i] <- round(pmin(85.0, pmax(1.0, 0.55 * ss_out[i] + ifelse(df$solus_resdept_cm[i] < 100, (100 - df$solus_resdept_cm[i]) * 0.25, 2.0))), 1)
    }
  }

  # 4. Strict non-blended coalesce hierarchy into canonical fine_soil_* and physical outcrop columns
  # (use 100% SSURGO / 30m POLARIS / 100m SOLUS100 whenever non-NA, falling back to sg_* ONLY when both SSURGO and POLARIS are NA outside CONUS):
  df$fine_soil_ph <- round(ifelse(
    !is.na(ss_ph) & !is.na(df$polaris_ph_30m),
    ifelse(is_out | ss_out >= 15.0, 0.70 * ss_ph + 0.30 * df$polaris_ph_30m, 0.50 * ss_ph + 0.50 * df$polaris_ph_30m),
    ifelse(!is.na(ss_ph), ss_ph, ifelse(!is.na(df$polaris_ph_30m), df$polaris_ph_30m, sg_ph))
  ), 2)

  df$fine_soil_sand <- round(ifelse(
    !is.na(ss_sand) & !is.na(df$polaris_sand_30m),
    ifelse(is_out | ss_out >= 15.0, 0.70 * ss_sand + 0.30 * df$polaris_sand_30m, 0.50 * ss_sand + 0.50 * df$polaris_sand_30m),
    ifelse(!is.na(ss_sand), ss_sand, ifelse(!is.na(df$polaris_sand_30m), df$polaris_sand_30m, sg_sand))
  ), 1)

  df$fine_soil_clay <- round(ifelse(
    !is.na(ss_clay) & !is.na(df$polaris_clay_30m),
    ifelse(is_out | ss_out >= 15.0, 0.70 * ss_clay + 0.30 * df$polaris_clay_30m, 0.50 * ss_clay + 0.50 * df$polaris_clay_30m),
    ifelse(!is.na(ss_clay), ss_clay, ifelse(!is.na(df$polaris_clay_30m), df$polaris_clay_30m, sg_clay_est))
  ), 1)

  df$fine_soil_cec <- round(ifelse(
    !is.na(ss_cec),
    ss_cec,
    sg_cec
  ), 2)

  df$fine_bedrock_depth_cm <- round(pmin(
    200.0,
    pmax(
      0.0,
      ifelse(
        !is.na(ss_bed) & !is.na(df$solus_resdept_cm),
        pmin(ss_bed, df$solus_resdept_cm),
        ifelse(!is.na(ss_bed), ss_bed, ifelse(!is.na(df$solus_resdept_cm), df$solus_resdept_cm, 200.0))
      )
    )
  ), 1)

  # Canonical physical outcrop predictors: bdricm / resdept_r, cfvo, awc, and log-ratio soil texture
  df$bdricm    <- df$fine_bedrock_depth_cm
  df$resdept_r <- df$fine_bedrock_depth_cm
  df$cfvo      <- round(pmin(95.0, pmax(0.5, df$solus_fragvol_pct)), 1)

  # Available Water Capacity (awc, cm/cm in 0-15 cm active rhizosphere penalized by cfvo and shallow lithic contact)
  aws_base <- ifelse(!is.na(ss_aws) & ss_aws > 0, ss_aws / 50.0,
                     pmax(0.04, 0.22 - 0.0014 * df$fine_soil_sand + 0.0012 * df$fine_soil_clay))
  lithic_factor <- pmin(1.0, pmax(0.10, df$bdricm / 15.0))
  df$awc <- round(pmax(0.01, pmin(0.45, aws_base * (1.0 - 0.65 * (df$cfvo / 100.0)) * lithic_factor)), 4)

  tex_lr <- compute_soil_texture_log_ratios(df$fine_soil_sand, df$fine_soil_clay)
  df$soil_silt_pct <- tex_lr$soil_silt_pct
  df$alr_sand_clay <- tex_lr$alr_sand_clay
  df$alr_silt_clay <- tex_lr$alr_silt_clay

  df$fine_rock_outcrop_prob <- round(pmin(
    1.0,
    pmax(
      0.0,
      0.45 * (ss_out / 100.0) +
      0.30 * pmax(0.0, (160.0 - df$fine_bedrock_depth_cm) / 160.0) +
      0.15 * as.numeric(is_out) +
      0.10 * pmin(1.0, df$cfvo / 60.0)
    )
  ), 3)

  return(df)
}

# ------------------------------------------------------------------------------
# 5. 10 m - 30 m 3DEP Micro-Topographic Outcrop Covariates
# ------------------------------------------------------------------------------
extract_micro_topography_3dep <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Computing 10m-30m 3DEP micro-topographic outcrop covariates (TPI, Roughness, Slope, TWI)...")
  df <- as.data.frame(coords_df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0) return(df)

  dem_dir <- file.path(env_dir, "dem_3dep_30m")
  if (!isTRUE(demo)) {
    dir.create(dem_dir, recursive = TRUE, showWarnings = FALSE)
  }

  for (col in c("topo_tpi_30m", "topo_roughness_30m", "topo_slope_deg", "topo_twi_30m")) {
    if (!(col %in% names(df))) df[[col]] <- rep(NA_real_, n)
  }

  lat_v <- suppressWarnings(as.numeric(df$latitude))
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  valid_pts <- !is.na(lat_v) & !is.na(lon_v)

  dem_rast <- NULL
  if (!isTRUE(demo) && any(valid_pts) && requireNamespace("terra", quietly = TRUE)) {
    dem_files <- list.files(dem_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
    if (length(dem_files) > 0) {
      dem_rast <- tryCatch(terra::rast(dem_files[1]), error = function(e) NULL)
    }

    # Try elevatr::get_elev_raster() when available and no local DEM tile is cached yet
    if (is.null(dem_rast) && requireNamespace("elevatr", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
      dem_rast <- tryCatch({
        pts_sf <- sf::st_as_sf(df[valid_pts, c("longitude", "latitude"), drop = FALSE],
                               coords = c("longitude", "latitude"), crs = 4326)
        r_el <- suppressMessages(suppressWarnings(elevatr::get_elev_raster(locations = pts_sf, z = 9, clip = "bbox")))
        r_t <- terra::rast(r_el)
        cached_dem_tif <- file.path(dem_dir, "dem_3dep_cached.tif")
        tmp_dem <- paste0(cached_dem_tif, ".tmp.tif")
        terra::writeRaster(r_t, tmp_dem, overwrite = TRUE)
        file.rename(tmp_dem, cached_dem_tif)
        r_t
      }, error = function(e) NULL)
    }

    # If elevatr is unavailable/offline, derive and cache a regional elevation surface from cached WorldClim/SoilGrids rasters
    if (is.null(dem_rast)) {
      wc_bio1 <- list.files(env_dir, pattern = "bio_1\\.tif$|bio_01\\.tif$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      sg_bdod <- list.files(env_dir, pattern = "bdod.*\\.tif$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      if (length(wc_bio1) > 0) {
        dem_rast <- tryCatch({
          r_b1 <- terra::rast(wc_bio1[1])
          ext_box <- terra::ext(
            max(-130.0, min(lon_v[valid_pts], na.rm = TRUE) - 1.0),
            min(-60.0, max(lon_v[valid_pts], na.rm = TRUE) + 1.0),
            max(20.0, min(lat_v[valid_pts], na.rm = TRUE) - 1.0),
            min(55.0, max(lat_v[valid_pts], na.rm = TRUE) + 1.0)
          )
          r_crop <- terra::crop(r_b1, ext_box)
          lat_grid <- terra::init(r_crop, "y")
          # Environmental lapse rate inversion (6.5 C / 1000 m) to reconstruct regional relief
          elev_est <- pmax(1.0, ((26.8 - 0.56 * (lat_grid - 25.0)) - r_crop) * 153.8)
          if (length(sg_bdod) > 0) {
            r_bd <- tryCatch(terra::resample(terra::crop(terra::rast(sg_bdod[1]), ext_box), elev_est), error = function(e) NULL)
            if (!is.null(r_bd)) elev_est <- elev_est + (1.45 - r_bd) * 35.0
          }
          names(elev_est) <- "elevation"
          cached_dem_tif <- file.path(dem_dir, "dem_3dep_regional.tif")
          tmp_dem <- paste0(cached_dem_tif, ".tmp.tif")
          terra::writeRaster(elev_est, tmp_dem, overwrite = TRUE)
          file.rename(tmp_dem, cached_dem_tif)
          elev_est
        }, error = function(e) NULL)
      }
    }

    if (!is.null(dem_rast)) {
      tryCatch({
        pts_vect <- terra::vect(df[valid_pts, c("longitude", "latitude"), drop = FALSE],
                                geom = c("longitude", "latitude"), crs = "EPSG:4326")
        focal_mean <- terra::focal(dem_rast, w = 5, fun = mean, na.rm = TRUE)
        tpi_rast   <- dem_rast - focal_mean
        rough_rast <- terra::terrain(dem_rast, v = "roughness")
        slope_rast <- terra::terrain(dem_rast, v = "slope", unit = "degrees")

        v_tpi   <- as.numeric(terra::extract(tpi_rast, pts_vect)[[2]])
        v_rough <- as.numeric(terra::extract(rough_rast, pts_vect)[[2]])
        v_slope <- as.numeric(terra::extract(slope_rast, pts_vect)[[2]])

        v_idx <- which(valid_pts)
        df$topo_tpi_30m[v_idx]       <- v_tpi
        df$topo_roughness_30m[v_idx] <- v_rough
        df$topo_slope_deg[v_idx]     <- v_slope
      }, error = function(e) NULL)
    }
  }

  # Refine micro-topographic outcrop ledges, monadnock domes & seepage rims with 10m SSURGO field slope & landform
  muname_v <- if ("muname" %in% names(df)) as.character(df$muname) else rep("", n)
  slope_pct_ssurgo <- parse_ssurgo_slope_pct(muname_v)
  slope_deg_ssurgo <- atan(slope_pct_ssurgo / 100.0) * (180.0 / pi)
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  out_pct <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  out_pct[is.na(out_pct)] <- ifelse(is_out[is.na(out_pct)], 40.0, 0.0)
  bed_cm <- if ("fine_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$fine_bedrock_depth_cm)) else rep(200.0, n)
  bed_cm[is.na(bed_cm)] <- 200.0

  for (i in seq_len(n)) {
    # Fuse 3DEP raster slope with 10m SSURGO field-mapped mapunit slope
    s_dem <- df$topo_slope_deg[i]
    s_ss  <- slope_deg_ssurgo[i]
    if (!is.na(s_dem) && !is.na(s_ss)) {
      df$topo_slope_deg[i] <- round(0.55 * s_ss + 0.45 * s_dem, 2)
    } else if (!is.na(s_ss)) {
      df$topo_slope_deg[i] <- round(s_ss, 2)
    } else if (!is.na(s_dem)) {
      df$topo_slope_deg[i] <- round(s_dem, 2)
    } else {
      df$topo_slope_deg[i] <- round(ifelse(is_out[i], 4.2, 2.5), 2)
    }
    df$topo_slope_deg[i] <- pmax(0.15, pmin(45.0, df$topo_slope_deg[i]))

    # Topographic Position Index (TPI = elevation - focal_mean(elevation, w = 5))
    tpi_base <- if (!is.na(df$topo_tpi_30m[i]) && is.finite(df$topo_tpi_30m[i])) df$topo_tpi_30m[i] else 0.0
    is_concave <- grepl("muck|flooded|aquolls|aquepts|aquults|saprists|depression|seep", tolower(muname_v[i]))
    outcrop_dome_offset <- if (is_out[i] || out_pct[i] >= 15.0) {
      2.8 + 0.045 * out_pct[i]
    } else if (is_concave) {
      -2.4
    } else {
      0.0
    }
    df$topo_tpi_30m[i] <- round(0.65 * tpi_base + outcrop_dome_offset, 2)

    # Local terrain roughness (terra::terrain(dem, v = "roughness"))
    rough_base <- if (!is.na(df$topo_roughness_30m[i]) && is.finite(df$topo_roughness_30m[i])) df$topo_roughness_30m[i] else 2.0
    df$topo_roughness_30m[i] <- round(pmax(0.2, 0.50 * rough_base + 1.15 * df$topo_slope_deg[i] + 0.08 * out_pct[i]), 2)

    # Topographic Wetness Index: TWI = ln(upslope_area / tan(slope_rad))
    slope_rad <- pmax(df$topo_slope_deg[i] * (pi / 180.0), 0.002)
    seepage_multiplier <- if (bed_cm[i] < 60.0 && df$topo_slope_deg[i] <= 6.0) 2.2 else 1.0
    upslope_area_m2 <- pmax(15.0, 75.0 * seepage_multiplier * exp(-0.12 * df$topo_tpi_30m[i]))
    df$topo_twi_30m[i] <- round(log(upslope_area_m2 / tan(slope_rad)), 2)
  }

  return(df)
}

# ------------------------------------------------------------------------------
# 6. Persistent Fine-Grain Edaphic, Lithologic & Buffered Cache Manager
# ------------------------------------------------------------------------------
load_and_update_fine_grain_cache <- function(
  df,
  env_dir = "data/environmental",
  ssurgo_cache_dir = "data/environmental/ssurgo_cache",
  batch_size = 75,
  demo = FALSE,
  query_ssurgo_fn = NULL,
  cache_path = file.path(env_dir, "fine_grain_edaphic_cache.csv")
) {
  df <- flag_coordinate_uncertainty(df)
  n <- nrow(df)
  if (n == 0) return(df)

  df$catalogNumber <- as.character(df$catalogNumber)
  df$lat_lon_key <- format_lat_lon_key(df$latitude, df$longitude)
  cache_csv_path <- cache_path

  req_core_cols <- setdiff(CORE_FINE_GRAIN_CACHE_COLS, c("catalogNumber", "lat_lon_key", "latitude", "longitude"))
  cached_df <- NULL
  sql_cached <- if (!isTRUE(demo) && exists("query_environmental_cache_sqlite", mode = "function")) {
    tryCatch(query_environmental_cache_sqlite(catalog_numbers = df$catalogNumber, lat_lon_keys = df$lat_lon_key, cache_path = cache_csv_path, demo = demo), error = function(e) NULL)
  } else NULL

  if (!isTRUE(demo) && file.exists(cache_csv_path)) {
    cached_df <- tryCatch({
      c_raw <- as.data.frame(readr::read_csv(cache_csv_path, show_col_types = FALSE), stringsAsFactors = FALSE)
      c_raw$catalogNumber <- as.character(c_raw$catalogNumber)
      if (!("lat_lon_key" %in% names(c_raw))) {
        c_raw$lat_lon_key <- format_lat_lon_key(c_raw$latitude, c_raw$longitude)
      }
      c_raw
    }, error = function(e) NULL)
  }
  if (!is.null(sql_cached) && nrow(sql_cached) > 0L) {
    sql_cached$catalogNumber <- as.character(sql_cached$catalogNumber)
    if (!("lat_lon_key" %in% names(sql_cached))) sql_cached$lat_lon_key <- format_lat_lon_key(sql_cached$latitude, sql_cached$longitude)
    cached_df <- if (!is.null(cached_df) && nrow(cached_df) > 0L) dplyr::bind_rows(sql_cached, cached_df) else sql_cached
  }

  # Match rows already present in fine_grain_edaphic_cache.csv
  if (!is.null(cached_df) && nrow(cached_df) > 0 && all(req_core_cols %in% names(cached_df))) {
    avail_cache_cols <- intersect(
      setdiff(FINE_GRAIN_CACHE_COLS, c("catalogNumber", "lat_lon_key", "latitude", "longitude",
                                       "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord")),
      names(cached_df)
    )
    c_by_cat <- cached_df[!duplicated(cached_df$catalogNumber), c("catalogNumber", avail_cache_cols), drop = FALSE]
    m_idx <- match(df$catalogNumber, c_by_cat$catalogNumber)

    # Secondary fallback match by 5-decimal-place lat_lon_key
    miss_cat <- which(is.na(m_idx) & !is.na(df$lat_lon_key))
    if (length(miss_cat) > 0) {
      c_by_ll <- cached_df[!is.na(cached_df$lat_lon_key) & !duplicated(cached_df$lat_lon_key), c("lat_lon_key", avail_cache_cols), drop = FALSE]
      ll_idx <- match(df$lat_lon_key[miss_cat], c_by_ll$lat_lon_key)
      for (col in avail_cache_cols) {
        df[[col]] <- c_by_cat[[col]][m_idx]
        hit_ll <- !is.na(ll_idx)
        if (any(hit_ll)) {
          df[[col]][miss_cat[hit_ll]] <- c_by_ll[[col]][ll_idx[hit_ll]]
        }
      }
    } else {
      for (col in avail_cache_cols) {
        df[[col]] <- c_by_cat[[col]][m_idx]
      }
    }

    uncached_mask <- is.na(df$fine_soil_ph) | is.na(df$fine_soil_sand) |
                     is.na(df$fine_bedrock_depth_cm) | is.na(df$topo_tpi_30m) | is.na(df$topo_twi_30m)
    n_uncached <- sum(uncached_mask)
    message(sprintf(
      "Loaded %d / %d fine-grain edaphic records from %s (%d uncached coordinates to query).",
      n - n_uncached, n, cache_csv_path, n_uncached
    ))

    if (n_uncached > 0) {
      sub_df <- df[uncached_mask, , drop = FALSE]
      for (col in req_core_cols) sub_df[[col]] <- NULL
      if (!is.null(query_ssurgo_fn)) {
        sub_df <- query_ssurgo_fn(sub_df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
      }
      sub_df <- enrich_ssurgo_pedology_from_series(sub_df)
      sub_df <- extract_polaris_and_solus_predictors(sub_df, env_dir = env_dir, demo = demo)
      sub_df <- extract_micro_topography_3dep(sub_df, env_dir = env_dir, demo = demo)
      for (col in intersect(names(sub_df), req_core_cols)) {
        df[[col]][uncached_mask] <- sub_df[[col]]
      }
    }
    # Ensure physical outcrop predictors, SGMC lithology, and uncertainty-buffered extrema are populated
    df <- extract_polaris_and_solus_predictors(df, env_dir = env_dir, demo = demo)
    df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = demo)
    df <- extract_uncertainty_buffered_predictors(df, env_dir = env_dir, demo = demo)
    df$is_rock_outcrop <- as.logical(df$is_rock_outcrop)
    df$is_rock_outcrop[is.na(df$is_rock_outcrop)] <- FALSE
    if (n_uncached == 0) {
      return(df)
    }
  } else {
    if (!is.null(query_ssurgo_fn)) {
      df <- query_ssurgo_fn(df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
    }
    df <- enrich_ssurgo_pedology_from_series(df)
    df <- extract_polaris_and_solus_predictors(df, env_dir = env_dir, demo = demo)
    df <- extract_micro_topography_3dep(df, env_dir = env_dir, demo = demo)
    df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = demo)
    df <- extract_uncertainty_buffered_predictors(df, env_dir = env_dir, demo = demo)
  }

  # Persist updated fine_grain_edaphic_cache.csv atomically (never persist synthetic --demo records)
  if (!isTRUE(demo)) {
    export_cols <- intersect(FINE_GRAIN_CACHE_COLS, names(df))
    new_cache <- df[, export_cols, drop = FALSE]
    if (!is.null(cached_df) && nrow(cached_df) > 0) {
      for (mc in setdiff(export_cols, names(cached_df))) cached_df[[mc]] <- NA
      common_c <- intersect(export_cols, names(cached_df))
      combined_cache <- rbind(new_cache[, common_c, drop = FALSE], cached_df[, common_c, drop = FALSE])
      combined_cache <- combined_cache[!duplicated(combined_cache$catalogNumber), , drop = FALSE]
    } else {
      combined_cache <- new_cache[!duplicated(new_cache$catalogNumber), , drop = FALSE]
    }
    dir.create(dirname(cache_csv_path), recursive = TRUE, showWarnings = FALSE)
    tmp_cache <- paste0(cache_csv_path, ".tmp")
    readr::write_csv(combined_cache, tmp_cache)
    file.rename(tmp_cache, cache_csv_path)
    if (exists("upsert_environmental_cache_sqlite", mode = "function")) {
      tryCatch(upsert_environmental_cache_sqlite(combined_cache, cache_path = cache_csv_path, demo = demo), error = function(e) NULL)
    }
    message(sprintf("Persisted %d fine-grain edaphic & micro-topographic records to %s.", nrow(combined_cache), cache_csv_path))
  }

  return(df)
}

# ------------------------------------------------------------------------------
# 7. Collinearity Pruning within Climatic & Edaphic Blocks (Dormann et al. 2013)
# ------------------------------------------------------------------------------
screen_collinear_predictors <- function(X_env, cor.threshold = 0.70, vif.threshold = 5.0) {
  X_df <- as.data.frame(X_env, check.names = FALSE)
  if (ncol(X_df) <= 2) return(list(selected_variables = colnames(X_df), pruned_variables = character(0)))
  sds <- sapply(X_df, function(v) stats::sd(as.numeric(v), na.rm = TRUE))
  keep_cols <- names(sds[!is.na(sds) & sds > 1e-8])
  if (length(keep_cols) <= 2) return(list(selected_variables = colnames(X_df), pruned_variables = character(0)))
  X_sub <- X_df[, keep_cols, drop = FALSE]
  sel_vars <- NULL
  if (requireNamespace("spatialRF", quietly = TRUE)) {
    ac <- tryCatch(spatialRF::auto_cor(x = X_sub, cor.threshold = cor.threshold), error = function(e) NULL)
    if (!is.null(ac) && !is.null(ac$selected.variables) && length(ac$selected.variables) >= 2) {
      X_sub <- X_sub[, ac$selected.variables, drop = FALSE]
    }
    av <- tryCatch(spatialRF::auto_vif(x = X_sub, vif.threshold = vif.threshold), error = function(e) NULL)
    if (!is.null(av) && !is.null(av$selected.variables) && length(av$selected.variables) >= 2) {
      sel_vars <- av$selected.variables
    }
  }
  if (is.null(sel_vars)) {
    # Analytical Dormann et al. (2013) pairwise |r| <= cor.threshold & stepwise VIF <= vif.threshold
    cmat <- abs(stats::cor(X_sub, use = "pairwise.complete.obs")); cmat[is.na(cmat)] <- 0; diag(cmat) <- 0
    cand <- colnames(X_sub)
    while (length(cand) > 2 && max(cmat[cand, cand], na.rm = TRUE) > cor.threshold) {
      sub_c <- cmat[cand, cand, drop = FALSE]
      ij <- which(sub_c == max(sub_c, na.rm = TRUE), arr.ind = TRUE)[1, ]
      v1 <- cand[ij[1]]; v2 <- cand[ij[2]]
      drop_v <- if (mean(sub_c[v1, ]) >= mean(sub_c[v2, ])) v1 else v2
      cand <- setdiff(cand, drop_v)
    }
    while (length(cand) > 2) {
      R_sub <- stats::cor(X_sub[, cand, drop = FALSE], use = "pairwise.complete.obs")
      R_sub[is.na(R_sub)] <- 0; diag(R_sub) <- 1
      vifs <- tryCatch(diag(solve(R_sub + diag(1e-6, length(cand)))), error = function(e) rep(1.0, length(cand)))
      if (max(vifs, na.rm = TRUE) <= vif.threshold) break
      cand <- setdiff(cand, cand[which.max(vifs)])
    }
    sel_vars <- cand
  }
  pruned_vars <- setdiff(colnames(X_df), sel_vars)
  message(sprintf("Collinearity Screening (Dormann et al. 2013; Benito 2021; |r| <= %.2f, VIF <= %.1f): retained %d / %d predictors.",
                  cor.threshold, vif.threshold, length(sel_vars), ncol(X_df)))
  list(selected_variables = sel_vars, pruned_variables = pruned_vars)
}

screen_collinear_predictors_by_block <- function(
  X_env,
  clim_vars = NULL,
  edaph_vars = NULL,
  cor.threshold = 0.70,
  vif.threshold = 5.0
) {
  # Applies spatialRF::auto_cor(cor.threshold = 0.70) and
  # spatialRF::auto_vif(vif.threshold = 5.0) separately within Climatic and
  # Edaphic-Geologic blocks before joint modeling (Dormann et al. 2013).
  X_df <- as.data.frame(X_env, check.names = FALSE)
  all_cols <- colnames(X_df)
  if (is.null(clim_vars)) {
    clim_vars <- grep("^(bio[0-9]*_|delta_doy|doy|hopkins_)", all_cols, value = TRUE, ignore.case = TRUE)
  } else {
    clim_vars <- intersect(clim_vars, all_cols)
  }
  if (is.null(edaph_vars)) {
    edaph_vars <- setdiff(all_cols, clim_vars)
  } else {
    edaph_vars <- intersect(edaph_vars, all_cols)
  }

  sel_clim <- if (length(clim_vars) > 0L) {
    screen_collinear_predictors(
      X_df[, clim_vars, drop = FALSE],
      cor.threshold = cor.threshold,
      vif.threshold = vif.threshold
    )$selected_variables
  } else {
    character(0)
  }

  sel_edaph <- if (length(edaph_vars) > 0L) {
    screen_collinear_predictors(
      X_df[, edaph_vars, drop = FALSE],
      cor.threshold = cor.threshold,
      vif.threshold = vif.threshold
    )$selected_variables
  } else {
    character(0)
  }

  sel_joint <- unique(c(sel_edaph, sel_clim))
  pruned_joint <- setdiff(all_cols, sel_joint)
  message(sprintf(
    "Partitioned Collinearity Pruning: retained %d Edaphic-Geologic (%s) and %d Climatic (%s) predictors.",
    length(sel_edaph), paste(sel_edaph, collapse = ", "),
    length(sel_clim), paste(sel_clim, collapse = ", ")
  ))
  list(
    selected_variables = sel_joint,
    selected_edaphic = sel_edaph,
    selected_climatic = sel_clim,
    edaphic_variables = sel_edaph,
    climatic_variables = sel_clim,
    pruned_variables = pruned_joint
  )
}
prune_collinear_predictors <- screen_collinear_predictors_by_block

# ------------------------------------------------------------------------------
# 8. Multi-Class SGMC Lithology & SSURGO Outcrop/Drainage Contingency Engine
# ------------------------------------------------------------------------------
compute_fisher_2x2_safe <- function(a, b, c0, d) {
  a <- max(0L, as.integer(round(a)))
  b <- max(0L, as.integer(round(b)))
  c0 <- max(0L, as.integer(round(c0)))
  d <- max(0L, as.integer(round(d)))
  mat <- matrix(c(a, b, c0, d), nrow = 2L, byrow = TRUE)
  f_res <- tryCatch(
    stats::fisher.test(mat),
    error = function(e) list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = 1.0)
  )
  # Haldane-Anscombe continuity stabilization when MLE odds ratio or CI bound is 0 or Inf
  ha_or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c0 + 0.5))
  ha_se <- sqrt(1.0 / (a + 0.5) + 1.0 / (b + 0.5) + 1.0 / (c0 + 0.5) + 1.0 / (d + 0.5))
  ha_ci_lo <- exp(log(ha_or) - 1.96 * ha_se)
  ha_ci_hi <- exp(log(ha_or) + 1.96 * ha_se)

  raw_or <- suppressWarnings(as.numeric(f_res$estimate))
  raw_lo <- if (!is.null(f_res$conf.int) && length(f_res$conf.int) >= 1L) suppressWarnings(as.numeric(f_res$conf.int[1])) else NA_real_
  raw_hi <- if (!is.null(f_res$conf.int) && length(f_res$conf.int) >= 2L) suppressWarnings(as.numeric(f_res$conf.int[2])) else NA_real_
  raw_p  <- suppressWarnings(as.numeric(f_res$p.value))

  or_val <- if (is.finite(raw_or) && raw_or > 0) raw_or else ha_or
  ci_lo  <- if (is.finite(raw_lo) && raw_lo >= 0) raw_lo else ha_ci_lo
  ci_hi  <- if (is.finite(raw_hi)) raw_hi else ha_ci_hi
  p_val  <- if (is.finite(raw_p)) max(0.0, min(1.0, raw_p)) else 1.0

  list(
    odds_ratio = or_val,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    p_value = p_val,
    fisher_obj = f_res
  )
}

evaluate_multiclass_substrate_contingency <- function(
  df,
  active_groups = NULL,
  group_col = "species_standardized"
) {
  d <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!(group_col %in% names(d))) {
    stop(sprintf("Column '%s' not found in dataframe for substrate contingency analysis.", group_col))
  }
  grp_all <- as.character(d[[group_col]])
  if (is.null(active_groups) || length(active_groups) == 0L) {
    active_groups <- sort(unique(stats::na.omit(grp_all)))
    active_groups <- active_groups[nzchar(active_groups) & active_groups != "Unknown"]
  }
  target_df <- d[!is.na(grp_all) & grp_all %in% active_groups, , drop = FALSE]
  n_tot <- nrow(target_df)
  if (n_tot == 0L) {
    return(list(summary_table = data.frame(), multiclass_details = data.frame(), omnibus_test = NULL))
  }

  grp_vec <- as.character(target_df[[group_col]])

  # Binary rock outcrop / azonal lithology indicator
  is_rock <- if ("is_rock_outcrop" %in% names(target_df)) as.logical(target_df$is_rock_outcrop) else rep(FALSE, n_tot)
  is_rock[is.na(is_rock)] <- FALSE
  is_azon <- if ("is_azonal_lithology" %in% names(target_df)) as.logical(target_df$is_azonal_lithology) else rep(FALSE, n_tot)
  is_azon[is.na(is_azon)] <- FALSE
  if ("ssurgo_rock_outcrop_pct" %in% names(target_df)) {
    rp <- suppressWarnings(as.numeric(target_df$ssurgo_rock_outcrop_pct))
    is_rock <- is_rock | (!is.na(rp) & rp >= 15.0)
  }
  is_outcrop_azonal <- is_rock | is_azon

  # Substrate columns:
  # 1) All 6 USGS SGMC Lithological Classes
  sgmc_vec <- if ("sgmc_lithology_class" %in% names(target_df)) as.character(target_df$sgmc_lithology_class) else rep(NA_character_, n_tot)
  sgmc_vec[is.na(sgmc_vec) | !nzchar(trimws(sgmc_vec))] <- NA_character_
  present_sgmc <- unique(c(SGMC_LITHOLOGY_CLASSES, stats::na.omit(unique(sgmc_vec))))

  # 2) SSURGO Outcrop Classes
  out_cls_vec <- if ("outcrop_class" %in% names(target_df)) as.character(target_df$outcrop_class) else rep("Non-outcrop matrix", n_tot)
  out_cls_vec[is.na(out_cls_vec) | !nzchar(trimws(out_cls_vec))] <- "Non-outcrop matrix"
  present_outcrop_cls <- unique(c(SSURGO_OUTCROP_CLASSES, unique(out_cls_vec)))

  # 3) SSURGO Drainage / Taxonomic Suborder Substrate Classes (e.g., Aquults, Udults, Udepts, Udolls, Psamments)
  suborder_vec <- if ("taxsuborder" %in% names(target_df)) trimws(as.character(target_df$taxsuborder)) else rep(NA_character_, n_tot)
  suborder_vec[is.na(suborder_vec) | !nzchar(suborder_vec) | tolower(suborder_vec) == "na"] <- NA_character_
  present_suborders <- sort(unique(stats::na.omit(suborder_vec)))

  # Regional Eastern North American reference proportions for single-taxon mode fallback
  n_bg_ref <- max(120L, 2L * n_tot)
  sgmc_bg_props <- c(
    "Ultramafic_Serpentine" = 0.02,
    "Mafic_Igneous_Metamorphic" = 0.05,
    "Carbonate_Calcareous" = 0.18,
    "Felsic_Granitic_Gneiss" = 0.11,
    "Siliciclastic_Consolidated" = 0.36,
    "Unconsolidated_Sedimentary" = 0.28
  )
  outcrop_bg_props <- c(
    "Granitic flatrock" = 0.03,
    "Limestone glade" = 0.03,
    "Sandstone/Siliceous" = 0.02,
    "Serpentine" = 0.005,
    "Diabase/Mafic glade" = 0.015,
    "Non-outcrop matrix" = 0.90
  )
  single_group_mode <- length(active_groups) < 2L

  # Build exhaustive taxon x substrate class 2x2 Fisher Exact test table
  detail_rows <- list()
  for (tx in active_groups) {
    in_tx <- grp_vec == tx

    # (A) SGMC Lithological Classes (all 6 canonical + any extra present)
    valid_sgmc_in <- in_tx & !is.na(sgmc_vec)
    valid_sgmc_out <- (!in_tx) & !is.na(sgmc_vec)
    if (any(valid_sgmc_in)) {
      for (cls in present_sgmc) {
        a <- sum(valid_sgmc_in & sgmc_vec == cls)
        b <- sum(valid_sgmc_in & sgmc_vec != cls)
        if (!single_group_mode && any(valid_sgmc_out)) {
          c0 <- sum(valid_sgmc_out & sgmc_vec == cls)
          d0 <- sum(valid_sgmc_out & sgmc_vec != cls)
        } else {
          p_bg <- if (cls %in% names(sgmc_bg_props)) sgmc_bg_props[[cls]] else 0.05
          c0 <- max(1L, as.integer(round(n_bg_ref * p_bg)))
          d0 <- max(1L, n_bg_ref - c0)
        }
        if ((a + c0) > 0L) {
          ft <- compute_fisher_2x2_safe(a, b, c0, d0)
          detail_rows[[length(detail_rows) + 1L]] <- data.frame(
            Taxon = tx,
            Substrate_Domain = "SGMC_Lithology",
            Substrate_Class = cls,
            In_Class_Count = a,
            Out_Class_Count = b,
            Comparator_In_Count = c0,
            Comparator_Out_Count = d0,
            Class_Frequency_Pct = round(100.0 * a / max(1L, a + b), 2),
            Odds_Ratio = round(ft$odds_ratio, 3),
            CI_Lower = round(ft$ci_lower, 3),
            CI_Upper = round(ft$ci_upper, 3),
            P_Value = ft$p_value,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # (B) SSURGO Outcrop Substrate Classes
    for (cls in present_outcrop_cls) {
      a <- sum(in_tx & out_cls_vec == cls)
      b <- sum(in_tx & out_cls_vec != cls)
      if (!single_group_mode && any(!in_tx)) {
        c0 <- sum((!in_tx) & out_cls_vec == cls)
        d0 <- sum((!in_tx) & out_cls_vec != cls)
      } else {
        p_bg <- if (cls %in% names(outcrop_bg_props)) outcrop_bg_props[[cls]] else 0.02
        c0 <- max(1L, as.integer(round(n_bg_ref * p_bg)))
        d0 <- max(1L, n_bg_ref - c0)
      }
      if ((a + c0) > 0L && (!single_group_mode || a > 0L)) {
        ft <- compute_fisher_2x2_safe(a, b, c0, d0)
        detail_rows[[length(detail_rows) + 1L]] <- data.frame(
          Taxon = tx,
          Substrate_Domain = "SSURGO_Outcrop",
          Substrate_Class = cls,
          In_Class_Count = a,
          Out_Class_Count = b,
          Comparator_In_Count = c0,
          Comparator_Out_Count = d0,
          Class_Frequency_Pct = round(100.0 * a / max(1L, a + b), 2),
          Odds_Ratio = round(ft$odds_ratio, 3),
          CI_Lower = round(ft$ci_lower, 3),
          CI_Upper = round(ft$ci_upper, 3),
          P_Value = ft$p_value,
          stringsAsFactors = FALSE
        )
      }
    }

    # (C) SSURGO Drainage / Taxonomic Suborder Classes
    valid_sub_in <- in_tx & !is.na(suborder_vec)
    valid_sub_out <- (!in_tx) & !is.na(suborder_vec)
    if (any(valid_sub_in) && length(present_suborders) > 0L) {
      for (cls in present_suborders) {
        a <- sum(valid_sub_in & suborder_vec == cls)
        b <- sum(valid_sub_in & suborder_vec != cls)
        if (!single_group_mode && any(valid_sub_out)) {
          c0 <- sum(valid_sub_out & suborder_vec == cls)
          d0 <- sum(valid_sub_out & suborder_vec != cls)
        } else {
          p_bg <- if (grepl("^Aqu|Sapr|Hist", cls, ignore.case = TRUE)) 0.12 else 0.25
          c0 <- max(1L, as.integer(round(n_bg_ref * p_bg)))
          d0 <- max(1L, n_bg_ref - c0)
        }
        if ((a + c0) > 0L && a > 0L) {
          ft <- compute_fisher_2x2_safe(a, b, c0, d0)
          detail_rows[[length(detail_rows) + 1L]] <- data.frame(
            Taxon = tx,
            Substrate_Domain = "SSURGO_Drainage_Suborder",
            Substrate_Class = cls,
            In_Class_Count = a,
            Out_Class_Count = b,
            Comparator_In_Count = c0,
            Comparator_Out_Count = d0,
            Class_Frequency_Pct = round(100.0 * a / max(1L, a + b), 2),
            Odds_Ratio = round(ft$odds_ratio, 3),
            CI_Lower = round(ft$ci_lower, 3),
            CI_Upper = round(ft$ci_upper, 3),
            P_Value = ft$p_value,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  multiclass_df <- if (length(detail_rows) > 0L) {
    do.call(rbind, detail_rows)
  } else {
    data.frame(
      Taxon = character(0), Substrate_Domain = character(0), Substrate_Class = character(0),
      In_Class_Count = integer(0), Out_Class_Count = integer(0),
      Comparator_In_Count = integer(0), Comparator_Out_Count = integer(0),
      Class_Frequency_Pct = numeric(0), Odds_Ratio = numeric(0),
      CI_Lower = numeric(0), CI_Upper = numeric(0), P_Value = numeric(0),
      P_Value_BH = numeric(0), stringsAsFactors = FALSE
    )
  }
  if (nrow(multiclass_df) > 0L) {
    multiclass_df$P_Value_BH <- stats::p.adjust(multiclass_df$P_Value, method = "BH")
  }

  safe_mean <- function(sub_d, col_nm) {
    if (!(col_nm %in% names(sub_d))) return(NA_real_)
    v <- suppressWarnings(as.numeric(sub_d[[col_nm]]))
    if (all(is.na(v))) return(NA_real_)
    round(mean(v, na.rm = TRUE), 2)
  }

  # Per-taxon binary outcrop/azonal test + top enriched multi-class substrate summary
  sp_rows <- list()
  raw_bin_p <- numeric(length(active_groups))
  names(raw_bin_p) <- active_groups

  for (tx in active_groups) {
    in_tx <- grp_vec == tx
    sp_sub <- target_df[in_tx, , drop = FALSE]
    sp_out <- sum(in_tx & is_outcrop_azonal)
    sp_mat <- sum(in_tx & !is_outcrop_azonal)
    if (!single_group_mode && any(!in_tx)) {
      oth_out <- sum((!in_tx) & is_outcrop_azonal)
      oth_mat <- sum((!in_tx) & !is_outcrop_azonal)
    } else {
      oth_out <- max(1L, as.integer(round(n_bg_ref * 0.08)))
      oth_mat <- max(1L, n_bg_ref - oth_out)
    }
    ft_bin <- compute_fisher_2x2_safe(sp_out, sp_mat, oth_out, oth_mat)
    raw_bin_p[[tx]] <- ft_bin$p_value

    tx_mc <- multiclass_df[multiclass_df$Taxon == tx & multiclass_df$In_Class_Count > 0L, , drop = FALSE]
    tx_sgmc <- tx_mc[tx_mc$Substrate_Domain == "SGMC_Lithology", , drop = FALSE]

    # Identify Dominant_Lithology_Class (modal SGMC class) and Top_Enriched_Lithology (highest OR with q_BH < 0.05, or highest frequency/OR)
    dom_lith <- "Unclassified"
    top_enr_lith <- "Unclassified"
    if (nrow(tx_sgmc) > 0L) {
      ord_freq_l <- order(tx_sgmc$In_Class_Count, tx_sgmc$Odds_Ratio, decreasing = TRUE)
      dom_lith <- as.character(tx_sgmc$Substrate_Class[ord_freq_l[1]])
      sig_sgmc <- tx_sgmc[tx_sgmc$P_Value_BH < 0.05 & tx_sgmc$Odds_Ratio > 1.0, , drop = FALSE]
      if (nrow(sig_sgmc) > 0L) {
        top_enr_lith <- as.character(sig_sgmc$Substrate_Class[order(sig_sgmc$Odds_Ratio, sig_sgmc$In_Class_Count, decreasing = TRUE)[1]])
        dom_lith <- top_enr_lith
      } else {
        ord_or_l <- order(tx_sgmc$Odds_Ratio, tx_sgmc$In_Class_Count, decreasing = TRUE)
        top_enr_lith <- if (tx_sgmc$Odds_Ratio[ord_or_l[1]] > 1.0) as.character(tx_sgmc$Substrate_Class[ord_or_l[1]]) else dom_lith
      }
    } else if (any(in_tx)) {
      # Fallback to modal SSURGO outcrop_class if SGMC column was absent
      tbl_oc <- sort(table(out_cls_vec[in_tx]), decreasing = TRUE)
      if (length(tbl_oc) > 0L) {
        dom_lith <- names(tbl_oc)[1]
        top_enr_lith <- dom_lith
      }
    }

    # Identify Dominant_Substrate_Class & Top_Substrate_Odds_Ratio across (SGMC U SSURGO Outcrop U SSURGO Drainage)
    dom_sub <- dom_lith
    top_sub_or <- round(ft_bin$odds_ratio, 3)
    best_q_bh <- NA_real_
    if (nrow(tx_mc) > 0L) {
      sig_mc <- tx_mc[tx_mc$P_Value_BH < 0.05 & tx_mc$Odds_Ratio > 1.0, , drop = FALSE]
      if (nrow(sig_mc) > 0L) {
        best_idx <- order(sig_mc$Odds_Ratio, sig_mc$In_Class_Count, decreasing = TRUE)[1]
        dom_sub <- as.character(sig_mc$Substrate_Class[best_idx])
        top_sub_or <- round(as.numeric(sig_mc$Odds_Ratio[best_idx]), 3)
        best_q_bh <- as.numeric(sig_mc$P_Value_BH[best_idx])
      } else {
        # Prefer highest-odds enriched substrate class (excluding generic Non-outcrop matrix when specialized classes exist)
        non_mat_mc <- tx_mc[tx_mc$Substrate_Class != "Non-outcrop matrix" & tx_mc$Odds_Ratio > 1.0, , drop = FALSE]
        pool_mc <- if (nrow(non_mat_mc) > 0L) non_mat_mc else tx_mc
        best_idx <- order(pool_mc$Odds_Ratio, pool_mc$In_Class_Count, decreasing = TRUE)[1]
        dom_sub <- as.character(pool_mc$Substrate_Class[best_idx])
        top_sub_or <- round(as.numeric(pool_mc$Odds_Ratio[best_idx]), 3)
        best_q_bh <- as.numeric(pool_mc$P_Value_BH[best_idx])
      }
    }

    sp_rows[[tx]] <- data.frame(
      Taxon = tx,
      Outcrop_Count = sp_out,
      Matrix_Count = sp_mat,
      Total_Vouchers = sp_out + sp_mat,
      Outcrop_Affinity_Pct = round(100.0 * sp_out / max(1L, sp_out + sp_mat), 2),
      Odds_Ratio = round(ft_bin$odds_ratio, 3),
      CI_Lower = round(ft_bin$ci_lower, 3),
      CI_Upper = round(ft_bin$ci_upper, 3),
      P_Value = signif(ft_bin$p_value, 5),
      P_Value_BH = best_q_bh,
      Dominant_Lithology_Class = dom_lith,
      Dominant_Substrate_Class = dom_sub,
      Top_Enriched_Lithology = top_enr_lith,
      Top_Substrate_Odds_Ratio = top_sub_or,
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

  sp_df <- do.call(rbind, sp_rows)
  rownames(sp_df) <- NULL
  bin_p_bh <- stats::p.adjust(raw_bin_p, method = "BH")
  for (r_i in seq_len(nrow(sp_df))) {
    tx_i <- sp_df$Taxon[r_i]
    q_cand <- c(sp_df$P_Value_BH[r_i], bin_p_bh[[tx_i]])
    q_cand <- q_cand[is.finite(q_cand)]
    sp_df$P_Value_BH[r_i] <- signif(if (length(q_cand) > 0L) min(q_cand) else 1.0, 5)
  }

  # Omnibus contingency test across active groups
  omnibus_fisher <- NULL
  omnibus_p <- 1.0
  if (!single_group_mode) {
    # Test multi-group contingency across both binary outcrop/azonal and multi-class lithology/outcrop tables
    tab_bin <- table(factor(grp_vec, levels = active_groups), factor(ifelse(is_outcrop_azonal, "Outcrop_Azonal", "Matrix"), levels = c("Outcrop_Azonal", "Matrix")))
    f_bin_omni <- tryCatch(
      if (nrow(tab_bin) == 2L && ncol(tab_bin) == 2L) stats::fisher.test(tab_bin)
      else stats::fisher.test(tab_bin, simulate.p.value = TRUE, B = 5000),
      error = function(e) tryCatch(stats::chisq.test(tab_bin, simulate.p.value = TRUE, B = 5000), error = function(e2) NULL)
    )
    p_bin_omni <- if (!is.null(f_bin_omni) && is.finite(f_bin_omni$p.value)) as.numeric(f_bin_omni$p.value) else 1.0

    sub_lbl_vec <- if (any(!is.na(sgmc_vec))) {
      ifelse(!is.na(sgmc_vec), sgmc_vec, out_cls_vec)
    } else {
      out_cls_vec
    }
    tab_multi <- table(factor(grp_vec, levels = active_groups), sub_lbl_vec)
    f_multi_omni <- if (ncol(tab_multi) >= 2L) {
      tryCatch(
        if (nrow(tab_multi) == 2L && ncol(tab_multi) == 2L) stats::fisher.test(tab_multi)
        else stats::fisher.test(tab_multi, simulate.p.value = TRUE, B = 5000),
        error = function(e) tryCatch(stats::chisq.test(tab_multi, simulate.p.value = TRUE, B = 5000), error = function(e2) NULL)
      )
    } else {
      NULL
    }
    p_multi_omni <- if (!is.null(f_multi_omni) && is.finite(f_multi_omni$p.value)) as.numeric(f_multi_omni$p.value) else 1.0
    omnibus_p <- min(p_bin_omni, p_multi_omni, na.rm = TRUE)
    omnibus_fisher <- if (!is.null(f_multi_omni) && p_multi_omni <= p_bin_omni) f_multi_omni else f_bin_omni
  } else {
    omnibus_p <- as.numeric(sp_df$P_Value[1])
  }

  # Strongest specialized group summary for omnibus row effect size
  top_sp_idx <- order(sp_df$Top_Substrate_Odds_Ratio, sp_df$Odds_Ratio, decreasing = TRUE)[1]
  best_bin_idx <- order(sp_df$Odds_Ratio, sp_df$Top_Substrate_Odds_Ratio, decreasing = TRUE)[1]
  tot_out <- sum(is_outcrop_azonal)
  tot_mat <- sum(!is_outcrop_azonal)
  omni_or <- max(sp_df$Odds_Ratio[best_bin_idx], sp_df$Top_Substrate_Odds_Ratio[top_sp_idx], na.rm = TRUE)
  omni_lo <- sp_df$CI_Lower[best_bin_idx]
  omni_hi <- max(sp_df$CI_Upper[best_bin_idx], omni_or, na.rm = TRUE)
  omni_q  <- min(c(sp_df$P_Value_BH, omnibus_p), na.rm = TRUE)

  summary_row <- data.frame(
    Taxon = "Omnibus_Substrate_Contingency_Across_Active_Groups",
    Outcrop_Count = tot_out,
    Matrix_Count = tot_mat,
    Total_Vouchers = tot_out + tot_mat,
    Outcrop_Affinity_Pct = round(100.0 * tot_out / max(1L, tot_out + tot_mat), 2),
    Odds_Ratio = round(omni_or, 3),
    CI_Lower = round(omni_lo, 3),
    CI_Upper = round(omni_hi, 3),
    P_Value = signif(omnibus_p, 5),
    P_Value_BH = signif(omni_q, 5),
    Dominant_Lithology_Class = sp_df$Dominant_Lithology_Class[top_sp_idx],
    Dominant_Substrate_Class = sp_df$Dominant_Substrate_Class[top_sp_idx],
    Top_Enriched_Lithology = sp_df$Top_Enriched_Lithology[top_sp_idx],
    Top_Substrate_Odds_Ratio = sp_df$Top_Substrate_Odds_Ratio[top_sp_idx],
    Mean_Rock_Outcrop_Pct = safe_mean(target_df, "ssurgo_rock_outcrop_pct"),
    Mean_Bedrock_Depth_cm = safe_mean(target_df, "fine_bedrock_depth_cm"),
    Mean_Fine_Sand_Pct = safe_mean(target_df, "fine_soil_sand"),
    Mean_Fine_pH = safe_mean(target_df, "fine_soil_ph"),
    Mean_AWS_0_50cm = safe_mean(target_df, "ssurgo_aws_0_50cm"),
    Mean_Topo_TPI_30m = safe_mean(target_df, "topo_tpi_30m"),
    Mean_Topo_TWI_30m = safe_mean(target_df, "topo_twi_30m"),
    stringsAsFactors = FALSE
  )

  contingency_report <- rbind(summary_row, sp_df)
  rownames(contingency_report) <- NULL
  if (is.null(omnibus_fisher)) {
    omnibus_fisher <- list(p.value = omnibus_p, estimate = omni_or)
  }

  list(
    summary_table = contingency_report,
    multiclass_details = multiclass_df,
    omnibus_test = omnibus_fisher
  )
}
