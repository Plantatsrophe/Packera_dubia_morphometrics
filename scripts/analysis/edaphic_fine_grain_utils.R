#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/edaphic_fine_grain_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Fine-grain (10 m to 30 m) edaphic, lithic, geologic, and micro-topographic
#   predictor extraction, uncertainty-buffered herbarium sampling, and caching
#   utilities for Step 05 (05_multimodal_spatial_rf.R):
#     1. Trapezoidal 0-15 cm herbaceous rhizosphere integration
#        ((5 * X_0_5 + 10 * X_5_15) / 15) across 30 m POLARIS, 30 m gSSURGO,
#        100 m SOLUS100, and 250 m SoilGrids fallback (Chaney et al. 2019;
#        Velazco et al. 2017).
#     2. Harmonized cross-border 0-15 cm rhizosphere predictors (Mod et al. 2016;
#        Nauman et al. 2024): fine_soil_ph, orthonormal isometric log-ratio soil
#        texture (ilr_texture_1, ilr_texture_2; Egozcue et al. 2003),
#        fine_soil_cec, fine_cfvo_pct, and fine_bedrock_depth_cm.
#     3. 30 m 3DEP micro-topographic predictors: DEM slope (topo_slope_deg),
#        Topographic Position Index (topo_tpi_30m), Heat Load Index
#        (topo_hli_30m; McCune & Keon 2002), and D8 hydrological Topographic
#        Wetness Index (topo_twi_30m; Beven & Kirkby 1979).
#     4. USGS State Geologic Map Compilation (SGMC; Horton et al. 2017)
#        lithological standardization into 6 geochemical functional classes for
#        active-cohort Fisher's exact substrate contingency testing.
#     5. Block-partitioned collinearity pruning (Dormann et al. 2013) via
#        spatialRF::auto_cor(cor.threshold = 0.70) and
#        spatialRF::auto_vif(vif.threshold = 5.0).
# ==============================================================================

SSURGO_OUTCROP_SERIES_REGEX <- "(?i)(\\brock outcrop\\b|\\brock land\\b|\\bflatrock\\b|\\bglade\\b|\\bboulders\\b|\\bwake\\b|\\bashlar\\b|\\bgladeville\\b|\\bbarfield\\b|\\bgoldston\\b|\\bmanteo\\b|\\bramsey\\b|\\blily\\b|\\bchrome\\b|\\bconowingo\\b)"
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
  "polaris_bd_30m", "polaris_ksat_30m",
  "solus_resdept_cm", "solus_fragvol_pct", "solus_sand_pct", "solus_clay_pct", "solus_silt_pct", "solus_cec",
  "topo_tpi_30m", "topo_slope_deg", "topo_hli_30m", "topo_twi_30m", "elevation_m",
  "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "soil_silt_pct",
  "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm"
)

FINE_GRAIN_CACHE_COLS <- c(
  CORE_FINE_GRAIN_CACHE_COLS,
  "sgmc_lithology_class", "is_azonal_lithology",
  "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
  "is_clim_eligible_coord", "is_edaphic_eligible_coord", "edaphic_coord_weight",
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

compute_soil_texture_log_ratios <- function(sand_pct, clay_pct, silt_pct = NULL) {
  # Orthonormal isometric log-ratio (ILR; Egozcue et al. 2003) soil texture
  # transformation breaking 3-part simplex closure (sand + silt + clay = 100%).
  # Strictly returns NA_real_ whenever sand or clay is missing (never imputing
  # synthetic 50% sand / 20% clay constants).
  s <- suppressWarnings(as.numeric(sand_pct))
  c <- suppressWarnings(as.numeric(clay_pct))
  n <- max(length(s), length(c))
  if (n == 0L) {
    return(data.frame(
      soil_silt_pct = numeric(0),
      ilr_texture_1 = numeric(0),
      ilr_texture_2 = numeric(0)
    ))
  }
  if (length(s) != n) s <- rep_len(s, n)
  if (length(c) != n) c <- rep_len(c, n)
  si_raw <- if (is.null(silt_pct)) rep(NA_real_, n) else rep_len(suppressWarnings(as.numeric(silt_pct)), n)

  valid <- is.finite(s) & is.finite(c) & (s >= 0.0) & (c >= 0.0)
  si_out <- rep(NA_real_, n)
  ilr1 <- rep(NA_real_, n)
  ilr2 <- rep(NA_real_, n)

  if (any(valid)) {
    sv <- pmax(0.5, s[valid])
    cv <- pmax(0.5, c[valid])
    siv_in <- si_raw[valid]
    siv <- ifelse(is.finite(siv_in) & siv_in >= 0.0, pmax(0.5, siv_in), pmax(0.5, 100.0 - sv - cv))
    tot <- sv + cv + siv
    sv_c <- 100.0 * (sv / tot)
    cv_c <- 100.0 * (cv / tot)
    siv_c <- 100.0 * (siv / tot)
    si_out[valid] <- round(siv_c, 2)
    ilr1[valid] <- round(sqrt(2.0 / 3.0) * log(sv_c / sqrt(siv_c * cv_c)), 4)
    ilr2[valid] <- round(sqrt(0.5) * log(siv_c / cv_c), 4)
  }

  data.frame(
    soil_silt_pct = si_out,
    ilr_texture_1 = ilr1,
    ilr_texture_2 = ilr2,
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
  lon = NA_real_,
  is_demo = FALSE,
  demo = NULL,
  generalize = "",
  rock_type = "",
  compname = "",
  taxsuborder = "",
  taxgreatgroup = "",
  reskind = "",
  soil_ph = NA_real_,
  soil_sand = NA_real_,
  soil_cec = NA_real_,
  bedrock_depth_cm = NA_real_,
  rock_outcrop_pct = NA_real_
) {
  # Standardizes raw USGS SGMC lithology attributes and/or SSURGO parent-material
  # descriptors into 6 canonical geochemical functional classes (Horton et al. 2017):
  #   1. Ultramafic_Serpentine
  #   2. Mafic_Igneous_Metamorphic
  #   3. Carbonate_Calcareous
  #   4. Felsic_Granitic_Gneiss
  #   5. Siliciclastic_Consolidated
  #   6. Unconsolidated_Sedimentary
  # In non-demo runs (is_demo = FALSE), never uses coordinate bounding boxes;
  # returns "Unknown_Lithology" when both text descriptors and diagnostic physical
  # soil indicators are absent.
  eff_demo <- if (!is.null(demo)) isTRUE(demo) else isTRUE(is_demo)
  n <- max(
    length(rocktype1), length(rocktype2), length(major1), length(major2),
    length(minor1), length(unit_name), length(generalize), length(rock_type),
    length(muname), length(compname), length(taxsuborder), length(taxgreatgroup),
    length(reskind), length(outcrop_class), length(soil_ph), length(soil_sand),
    length(soil_cec), length(bedrock_depth_cm), length(rock_outcrop_pct),
    length(lat), length(lon), 1L
  )
  rep_str <- function(x) {
    if (is.null(x) || length(x) == 0L) return(rep("", n))
    y <- as.character(x)
    y[is.na(y)] <- ""
    rep_len(y, n)
  }
  rep_num <- function(x) {
    if (is.null(x) || length(x) == 0L) return(rep(NA_real_, n))
    rep_len(suppressWarnings(as.numeric(x)), n)
  }
  r1 <- rep_str(rocktype1); r2 <- rep_str(rocktype2)
  m1 <- rep_str(major1);    m2 <- rep_str(major2); mn1 <- rep_str(minor1)
  un <- rep_str(unit_name); gn <- rep_str(generalize); rt <- rep_str(rock_type)
  mu <- rep_str(muname);    cn <- rep_str(compname)
  so <- rep_str(taxsuborder); gg <- rep_str(taxgreatgroup); rk <- rep_str(reskind)
  oc <- rep_str(outcrop_class)

  ph_v   <- rep_num(soil_ph)
  sand_v <- rep_num(soil_sand)
  cec_v  <- rep_num(soil_cec)
  bed_v  <- rep_num(bedrock_depth_cm)
  ro_v   <- rep_num(rock_outcrop_pct)
  lat_v  <- rep_num(lat)
  lon_v  <- rep_num(lon)

  geol_txt <- tolower(trimws(paste(r1, r2, m1, m2, mn1, un, gn, rt)))
  oc_clean <- tolower(trimws(oc))
  oc_diag <- ifelse(oc_clean == "non-outcrop matrix", "", oc_clean)
  soil_txt <- tolower(trimws(paste(mu, cn, so, gg, rk, oc_diag)))
  combined <- trimws(paste(geol_txt, soil_txt))

  m_ultra <- grepl("ultramafic|serpentin|peridotit|dunit|pyroxenit|soapstone|talc|ophiolit|chrome|conowingo|aldino", combined, perl = TRUE) | (oc_clean == "serpentine")
  m_mafic <- !m_ultra & (grepl("mafic|diabas|gabbro|basalt|amphibolit|greenstone|metavolcanic|metadiorit|norite|anorthosit|trap rock|iredell|mecklenburg|wilkes|poindexter|enon", combined, perl = TRUE) | (oc_clean == "diabase/mafic glade"))
  m_carb  <- !(m_ultra | m_mafic) & (grepl("carbonat|limeston|dolostone|dolomit|marble|chalk|marl|calcarenit|travertin|karst|coquina|alvar|gladeville|barfield|talbott|colbert|rendolls|calciustolls|calciudolls", combined, perl = TRUE) | (oc_clean == "limestone glade"))
  m_fels  <- !(m_ultra | m_mafic | m_carb) & (grepl("felsic|crystalline|granit|granodiorit|monzonit|quartz monzonit|gneiss|pegmatit|rhyolit|aplit|dacit|migmatit|wake|louisburg|ashlar|rion|wateree|pacolet|appling|cecil|wedowee|madison", combined, perl = TRUE) | (oc_clean == "granitic flatrock"))
  m_silic <- !(m_ultra | m_mafic | m_carb | m_fels) & (grepl("siliciclastic|\\bconsolidated\\b|sandstone|quartzit|arkos|conglomerat|greywack|graywack|siltston|shale|slate|phyllit|mudston|argillit|chert|novaculit|schist|ramsey|lily|sewanee|goldston|manteo|wallen|dekab", combined, perl = TRUE) | (oc_clean == "sandstone/siliceous"))
  m_uncon <- !(m_ultra | m_mafic | m_carb | m_fels | m_silic) & grepl("unconsolidated|alluvi|colluvi|coastal plain|terrace|loess|till|outwash|glacial|dune|eolian|sandhill|muck|peat|saprist|psamment|aquult|aquoll|aquept|fluvent|wagram|lakeland|plainfield|goldsboro|norfolk|myatt|rains|lynchburg|tama", combined, perl = TRUE)

  out_cls <- character(n)
  out_cls[m_ultra] <- "Ultramafic_Serpentine"
  out_cls[m_mafic] <- "Mafic_Igneous_Metamorphic"
  out_cls[m_carb]  <- "Carbonate_Calcareous"
  out_cls[m_fels]  <- "Felsic_Granitic_Gneiss"
  out_cls[m_silic] <- "Siliciclastic_Consolidated"
  out_cls[m_uncon] <- "Unconsolidated_Sedimentary"

  rem_idx <- which(!nzchar(out_cls))
  if (length(rem_idx) > 0L) {
    rem_lithic_txt <- grepl("lithic|paralithic|bedrock|rock outcrop", combined[rem_idx], perl = TRUE)
    for (k in seq_along(rem_idx)) {
      i <- rem_idx[k]
      # Diagnostic physical soil/bedrock inference (strictly coordinate-free)
      ph_i <- ph_v[i]; sand_i <- sand_v[i]; cec_i <- cec_v[i]
      bed_i <- bed_v[i]; ro_i <- ro_v[i]
      is_shallow_lithic <- (is.finite(bed_i) && bed_i <= 50.0) ||
                           (is.finite(ro_i) && ro_i >= 15.0) ||
                           rem_lithic_txt[k]
      if (is_shallow_lithic && is.finite(ph_i) && ph_i >= 7.0) {
        out_cls[i] <- "Carbonate_Calcareous"
      } else if (is_shallow_lithic && is.finite(ph_i) && ph_i >= 6.2 &&
                 is.finite(cec_i) && cec_i >= 20.0 && (!is.finite(sand_i) || sand_i < 45.0)) {
        out_cls[i] <- "Mafic_Igneous_Metamorphic"
      } else if (is_shallow_lithic && is.finite(ph_i) && ph_i <= 5.3 &&
                 is.finite(sand_i) && sand_i >= 60.0) {
        out_cls[i] <- "Felsic_Granitic_Gneiss"
      } else if (is_shallow_lithic) {
        out_cls[i] <- "Siliciclastic_Consolidated"
      } else if (is.finite(bed_i) && bed_i >= 150.0 &&
                 (!is.finite(ro_i) || ro_i < 5.0) &&
                 ((is.finite(sand_i) && sand_i >= 70.0) || nzchar(combined[i]))) {
        out_cls[i] <- "Unconsolidated_Sedimentary"
      } else if (isTRUE(eff_demo)) {
        # Synthetic demonstration fallback gated strictly behind demo == TRUE
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
      } else {
        out_cls[i] <- "Unknown_Lithology"
      }
    }
  }
  return(out_cls)
}

extract_sgmc_lithology <- function(df, env_dir = "data/environmental", demo = FALSE) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0L) return(df)

  sgmc_dir <- file.path(env_dir, "sgmc")
  raw_rock1 <- if ("sgmc_rocktype1" %in% names(df)) as.character(df$sgmc_rocktype1) else if ("rocktype1" %in% names(df)) as.character(df$rocktype1) else rep("", n)
  raw_rock2 <- if ("sgmc_rocktype2" %in% names(df)) as.character(df$sgmc_rocktype2) else if ("rocktype2" %in% names(df)) as.character(df$rocktype2) else rep("", n)
  raw_maj1  <- if ("sgmc_major1" %in% names(df)) as.character(df$sgmc_major1) else if ("major1" %in% names(df)) as.character(df$major1) else rep("", n)
  raw_unit  <- if ("sgmc_unit_name" %in% names(df)) as.character(df$sgmc_unit_name) else if ("unit_name" %in% names(df)) as.character(df$unit_name) else rep("", n)
  raw_gen   <- if ("sgmc_generalize" %in% names(df)) as.character(df$sgmc_generalize) else if ("GENERALIZE" %in% names(df)) as.character(df$GENERALIZE) else rep("", n)
  raw_rtype <- if ("sgmc_rock_type" %in% names(df)) as.character(df$sgmc_rock_type) else if ("ROCK_TYPE" %in% names(df)) as.character(df$ROCK_TYPE) else rep("", n)
  prev_sgmc <- if ("sgmc_lithology_class" %in% names(df)) as.character(df$sgmc_lithology_class) else rep("", n)
  prev_sgmc[is.na(prev_sgmc) | tolower(trimws(prev_sgmc)) %in% c("na", "unknown", "unknown_lithology")] <- ""
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
        pick_sgmc_col <- function(cands) {
          for (cc in cands) {
            idx_c <- match(cc, cn_up)
            if (!is.na(idx_c)) return(names(sgmc_sf)[idx_c])
          }
          NA_character_
        }
        col_r1 <- pick_sgmc_col(c("ROCKTYPE1", "MAJOR1"))
        col_r2 <- pick_sgmc_col(c("ROCKTYPE2", "MAJOR2", "MINOR1"))
        col_m1 <- pick_sgmc_col(c("MAJOR1", "ROCKTYPE1"))
        col_un <- pick_sgmc_col(c("UNIT_NAME", "UNITNAME"))
        col_gn <- pick_sgmc_col(c("GENERALIZE", "GENERALIZED_LITH", "GEN_LITH"))
        col_rt <- pick_sgmc_col(c("ROCK_TYPE", "ORIG_LABEL"))
        for (k in seq_along(v_idx)) {
          h_k <- hits[[k]]
          if (length(h_k) > 0L) {
            row_g <- h_k[1]
            if (!is.na(col_r1)) raw_rock1[v_idx[k]] <- as.character(sgmc_sf[[col_r1]][row_g])
            if (!is.na(col_r2)) raw_rock2[v_idx[k]] <- as.character(sgmc_sf[[col_r2]][row_g])
            if (!is.na(col_m1)) raw_maj1[v_idx[k]]  <- as.character(sgmc_sf[[col_m1]][row_g])
            if (!is.na(col_un)) raw_unit[v_idx[k]]  <- as.character(sgmc_sf[[col_un]][row_g])
            if (!is.na(col_gn)) raw_gen[v_idx[k]]   <- as.character(sgmc_sf[[col_gn]][row_g])
            if (!is.na(col_rt)) raw_rtype[v_idx[k]] <- as.character(sgmc_sf[[col_rt]][row_g])
          }
        }
      }, error = function(e) NULL)
    }
  }

  muname_v <- if ("muname" %in% names(df)) as.character(df$muname) else rep("", n)
  comp_v   <- if ("compname" %in% names(df)) as.character(df$compname) else rep("", n)
  subord_v <- if ("taxsuborder" %in% names(df)) as.character(df$taxsuborder) else rep("", n)
  grtgrp_v <- if ("taxgreatgroup" %in% names(df)) as.character(df$taxgreatgroup) else rep("", n)
  resk_v   <- if ("reskind" %in% names(df)) as.character(df$reskind) else rep("", n)
  ocl_v    <- if ("outcrop_class" %in% names(df)) as.character(df$outcrop_class) else rep("", n)
  ph_col   <- if ("fine_soil_ph" %in% names(df)) df$fine_soil_ph else if ("soil_ph" %in% names(df)) df$soil_ph else rep(NA_real_, n)
  sand_col <- if ("fine_soil_sand" %in% names(df)) df$fine_soil_sand else if ("soil_sand" %in% names(df)) df$soil_sand else rep(NA_real_, n)
  cec_col  <- if ("fine_soil_cec" %in% names(df)) df$fine_soil_cec else if ("soil_cec" %in% names(df)) df$soil_cec else rep(NA_real_, n)
  bed_col  <- if ("fine_bedrock_depth_cm" %in% names(df)) df$fine_bedrock_depth_cm else if ("ssurgo_bedrock_depth_cm" %in% names(df)) df$ssurgo_bedrock_depth_cm else rep(NA_real_, n)
  ro_col   <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) df$ssurgo_rock_outcrop_pct else rep(NA_real_, n)

  gen_combined <- ifelse(nzchar(trimws(raw_gen)), raw_gen, prev_sgmc)
  lith_cls <- standardize_sgmc_lithology(
    rocktype1 = raw_rock1,
    rocktype2 = raw_rock2,
    major1 = raw_maj1,
    unit_name = raw_unit,
    muname = muname_v,
    outcrop_class = ocl_v,
    lat = lat_v,
    lon = lon_v,
    is_demo = isTRUE(demo),
    generalize = gen_combined,
    rock_type = raw_rtype,
    compname = comp_v,
    taxsuborder = subord_v,
    taxgreatgroup = grtgrp_v,
    reskind = resk_v,
    soil_ph = ph_col,
    soil_sand = sand_col,
    soil_cec = cec_col,
    bedrock_depth_cm = bed_col,
    rock_outcrop_pct = ro_col
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
  df$is_azonal_lithology <- (lith_cls == "Ultramafic_Serpentine") |
                            ((lith_cls %in% azonal_classes) & (is_out | out_pct >= 15.0 | bed_cm <= 40.0))
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
  max_edaphic_uncertainty_m = 2000.0,
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
  df$is_clim_eligible_coord <- !df$high_coord_uncertainty_flag
  yr_raw <- if ("year" %in% names(df)) suppressWarnings(as.numeric(df$year)) else rep(NA_real_, n)
  has_year_col <- ("year" %in% names(df)) && any(is.finite(yr_raw))
  modern_gps_fallback <- if (has_year_col) (is.finite(yr_raw) & yr_raw >= 2000) else rep(TRUE, n)
  df$is_edaphic_eligible_coord <- (is.finite(unc_raw) & unc_raw > 0 & (unc_raw <= as.numeric(max_edaphic_uncertainty_m))) |
                                  ((is.na(unc_raw) | !is.finite(unc_raw)) & modern_gps_fallback)
  df$edaphic_coord_weight <- ifelse(
    is.finite(unc_raw) & unc_raw > 0,
    ifelse(unc_raw <= 500.0, 1.0, ifelse(unc_raw <= as.numeric(max_edaphic_uncertainty_m), 0.70, 0.0)),
    ifelse(modern_gps_fallback, 0.60, 0.0)
  )
  df$buffer_radius_m <- compute_uncertainty_buffer_radius(
    unc_raw, default_r = default_r, min_r = min_r, max_r = max_r
  )
  return(df)
}

.list_real_fine_rasters <- function(dir_path, exclude_pattern = "^dem_3dep_regional\\.tif$") {
  if (is.null(dir_path) || !nzchar(dir_path) || !dir.exists(dir_path)) return(character(0))
  fls <- list.files(dir_path, pattern = "\\.(tif|tiff|hgt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fls) == 0L) return(character(0))
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
    fls <- fls[!grepl(exclude_pattern, basename(fls), ignore.case = TRUE)]
  }
  fls[file.exists(fls) & file.info(fls)$size > 256]
}

.project_vect_to_rast <- function(v_geom, r_rast) {
  if (is.null(v_geom) || is.null(r_rast)) return(v_geom)
  r_crs <- tryCatch(terra::crs(r_rast), error = function(e) "")
  v_crs <- tryCatch(terra::crs(v_geom), error = function(e) "")
  if (nzchar(r_crs) && nzchar(v_crs) && !isTRUE(tryCatch(terra::same.crs(v_geom, r_rast), error = function(e) FALSE))) {
    v_geom <- tryCatch(terra::project(v_geom, r_crs), error = function(e) v_geom)
  }
  v_geom
}

.rescale_solus_ph_vec <- function(v_raw) {
  v <- suppressWarnings(as.numeric(v_raw))
  fin <- v[is.finite(v)]
  if (length(fin) == 0L) return(v)
  med_v <- stats::median(fin, na.rm = TRUE)
  if (is.finite(med_v) && med_v > 100.0) {
    v <- v / 100.0
  } else if (is.finite(med_v) && med_v > 14.5) {
    v <- v / 10.0
  }
  v[!is.finite(v) | v < 2.0 | v > 11.5] <- NA_real_
  v
}

extract_uncertainty_buffered_predictors <- function(df, env_dir = "data/environmental", demo = FALSE) {
  # Computes uncertainty-buffered edaphic extrema within adaptive radius r_i
  # strictly from real local 30 m / 100 m GeoTIFF rasters when demo == FALSE.
  # When fine-grain rasters are absent on disk in non-demo runs, buffer columns
  # remain NA rather than fabricating deterministic algebraic clones.
  df <- flag_coordinate_uncertainty(df)
  n <- nrow(df)
  if (n == 0L) return(df)

  if (!("sgmc_lithology_class" %in% names(df))) {
    df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = demo)
  }

  df$ph_q10_buffer <- rep(NA_real_, n)
  df$ph_q90_buffer <- rep(NA_real_, n)
  df$min_bedrock_depth_buffer_cm <- rep(NA_real_, n)
  df$max_cfvo_buffer_pct <- rep(NA_real_, n)
  df$azonal_lithology_presence_buffer <- rep(NA, n)

  if (!isTRUE(demo)) {
    pol_files <- c(
      .list_real_fine_rasters(file.path(env_dir, "polaris_30m")),
      .list_real_fine_rasters(file.path(env_dir, "gssurgo_30m"))
    )
    sol_files <- c(
      .list_real_fine_rasters(file.path(env_dir, "solus_100m")),
      .list_real_fine_rasters(file.path(env_dir, "gssurgo_30m"))
    )
    lat_v <- suppressWarnings(as.numeric(df$latitude))
    lon_v <- suppressWarnings(as.numeric(df$longitude))
    valid_pts <- is.finite(lat_v) & is.finite(lon_v) & df$is_edaphic_eligible_coord
    if (any(valid_pts) && requireNamespace("terra", quietly = TRUE) && (length(pol_files) > 0L || length(sol_files) > 0L)) {
      v_idx <- which(valid_pts)
      pts_v <- tryCatch(
        terra::vect(df[v_idx, c("longitude", "latitude"), drop = FALSE], geom = c("longitude", "latitude"), crs = "EPSG:4326"),
        error = function(e) NULL
      )
      buf_v <- if (!is.null(pts_v)) {
        tryCatch(terra::buffer(pts_v, width = pmax(100.0, pmin(1000.0, as.numeric(df$buffer_radius_m[v_idx])))), error = function(e) NULL)
      } else NULL

      if (!is.null(buf_v)) {
        ph_tifs <- c(
          pol_files[grepl("^ph_", basename(pol_files), ignore.case = TRUE)],
          sol_files[grepl("^ph1to1h2o.*_5_cm_p|^ph1to1h2o|^ph_", basename(sol_files), ignore.case = TRUE)]
        )
        if (length(ph_tifs) > 0L) {
          tryCatch({
            r_ph <- terra::rast(ph_tifs[1L])
            buf_ph <- .project_vect_to_rast(buf_v, r_ph)
            ex_q10 <- .rescale_solus_ph_vec(terra::extract(r_ph, buf_ph, fun = function(x) stats::quantile(x, probs = 0.10, na.rm = TRUE))[[2]])
            ex_q90 <- .rescale_solus_ph_vec(terra::extract(r_ph, buf_ph, fun = function(x) stats::quantile(x, probs = 0.90, na.rm = TRUE))[[2]])
            df$ph_q10_buffer[v_idx] <- round(as.numeric(ex_q10), 2)
            df$ph_q90_buffer[v_idx] <- round(as.numeric(ex_q90), 2)
          }, error = function(e) NULL)
        }
        res_tifs <- sol_files[grepl("anylithicdpt|resdept|bdricm", basename(sol_files), ignore.case = TRUE)]
        if (length(res_tifs) > 0L) {
          min_accum <- rep(Inf, length(v_idx))
          for (rt_p in res_tifs) {
            tryCatch({
              r_bd <- terra::rast(rt_p)
              buf_bd <- .project_vect_to_rast(buf_v, r_bd)
              ex_m <- suppressWarnings(as.numeric(terra::extract(r_bd, buf_bd, fun = function(x) mean(x, na.rm = TRUE))[[2]]))
              ok_m <- is.finite(ex_m) & ex_m >= 0.0
              min_accum[ok_m] <- pmin(min_accum[ok_m], ex_m[ok_m])
            }, error = function(e) NULL)
          }
          min_accum[!is.finite(min_accum)] <- NA_real_
          df$min_bedrock_depth_buffer_cm[v_idx] <- round(pmin(200.0, pmax(0.0, min_accum)), 1)
        }
        cf_tifs <- sol_files[grepl("fragvol_0_cm_p|fragvol_5_cm_p|fragvol_15_cm_p|fragvol|cfvo", basename(sol_files), ignore.case = TRUE)]
        if (length(cf_tifs) > 0L) {
          max_cf <- rep(-Inf, length(v_idx))
          for (cf_p in utils::head(cf_tifs, 3L)) {
            tryCatch({
              r_cf <- terra::rast(cf_p)
              buf_cf <- .project_vect_to_rast(buf_v, r_cf)
              ex_c <- suppressWarnings(as.numeric(terra::extract(r_cf, buf_cf, fun = function(x) max(x, na.rm = TRUE))[[2]]))
              ok_c <- is.finite(ex_c) & ex_c >= 0.0
              max_cf[ok_c] <- pmax(max_cf[ok_c], ex_c[ok_c])
            }, error = function(e) NULL)
          }
          max_cf[!is.finite(max_cf)] <- NA_real_
          df$max_cfvo_buffer_pct[v_idx] <- round(pmin(100.0, pmax(0.0, max_cf)), 1)
        }
      }
    }
    return(df)
  }

  # Explicit demonstration simulation mode ONLY (demo == TRUE)
  r_m <- df$buffer_radius_m
  ph_mean <- if ("fine_soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_ph)) else if ("soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph)) else rep(5.6, n)
  ph_mean[!is.finite(ph_mean)] <- 5.6
  bed_mean <- if ("fine_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$fine_bedrock_depth_cm)) else rep(150.0, n)
  bed_mean[!is.finite(bed_mean)] <- 150.0
  cfvo_mean <- if ("fine_cfvo_pct" %in% names(df)) suppressWarnings(as.numeric(df$fine_cfvo_pct)) else if ("solus_fragvol_pct" %in% names(df)) suppressWarnings(as.numeric(df$solus_fragvol_pct)) else rep(10.0, n)
  cfvo_mean[!is.finite(cfvo_mean)] <- 10.0
  out_pct <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  out_pct[!is.finite(out_pct)] <- 0.0
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  rough_v <- if ("topo_slope_deg" %in% names(df)) suppressWarnings(as.numeric(df$topo_slope_deg)) else rep(3.0, n)
  rough_v[!is.finite(rough_v)] <- 3.0
  lith_cls <- as.character(df$sgmc_lithology_class)
  is_azonal <- if ("is_azonal_lithology" %in% names(df)) as.logical(df$is_azonal_lithology) else rep(FALSE, n)
  is_azonal[is.na(is_azonal)] <- FALSE

  buf_scale <- pmax(0.15, log10(pmax(100.0, r_m) / 50.0))
  het_ph <- pmin(1.15, pmax(0.12, (0.22 + 0.015 * rough_v + 0.006 * out_pct) * buf_scale))
  acid_pull <- ifelse(lith_cls %in% c("Felsic_Granitic_Gneiss", "Siliciclastic_Consolidated") & (is_out | out_pct >= 15.0), 0.35, 0.0)
  calc_pull <- ifelse(lith_cls %in% c("Carbonate_Calcareous", "Mafic_Igneous_Metamorphic", "Ultramafic_Serpentine") & (is_out | out_pct >= 15.0), 0.40, 0.0)

  df$ph_q10_buffer <- round(pmax(3.5, pmin(ph_mean, ph_mean - het_ph - acid_pull)), 2)
  df$ph_q90_buffer <- round(pmin(8.6, pmax(ph_mean, ph_mean + het_ph + calc_pull)), 2)
  shallow_pull <- ifelse(is_out | out_pct >= 20.0, 0.45, ifelse(out_pct >= 5.0 | is_azonal, 0.70, 0.88)) * (1.0 - 0.12 * (buf_scale - 0.3))
  df$min_bedrock_depth_buffer_cm <- round(pmax(0.0, pmin(bed_mean, bed_mean * pmax(0.15, shallow_pull) - 0.25 * out_pct)), 1)
  cfvo_boost <- (8.0 + 0.35 * out_pct + 0.55 * rough_v) * buf_scale
  df$max_cfvo_buffer_pct <- round(pmin(95.0, pmax(cfvo_mean, cfvo_mean + cfvo_boost)), 1)
  df$azonal_lithology_presence_buffer <- as.logical(
    is_azonal | is_out | (out_pct >= 15.0) | (df$min_bedrock_depth_buffer_cm <= 40.0)
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

enrich_ssurgo_pedology_from_series <- function(d, demo = FALSE) {
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

  # In production (demo == FALSE), never invent laboratory soil chemistry or
  # bedrock depth point-masses from regex matches on soil series names; preserve
  # measured SSURGO tabular horizon values or leave NA_real_ for continuous rasters.
  if (!isTRUE(demo)) {
    for (i in seq_len(n)) {
      if (!is.na(d$ssurgo_bedrock_depth_cm[i]) && is.finite(d$ssurgo_bedrock_depth_cm[i])) {
        d$ssurgo_bedrock_depth_cm[i] <- pmin(200.0, pmax(0.0, d$ssurgo_bedrock_depth_cm[i]))
      }
      muname <- if ("muname" %in% names(d) && !is.na(d$muname[i])) as.character(d$muname[i]) else ""
      is_out <- if ("is_rock_outcrop" %in% names(d) && !is.na(d$is_rock_outcrop[i])) isTRUE(d$is_rock_outcrop[i]) else FALSE
      if (is.na(d$ssurgo_rock_outcrop_pct[i]) && (nzchar(muname) || is_out)) {
        if (grepl("^rock outcrop$|^rock land$", tolower(trimws(muname)))) {
          d$ssurgo_rock_outcrop_pct[i] <- 85.0
        } else if (grepl("rock outcrop", tolower(muname))) {
          d$ssurgo_rock_outcrop_pct[i] <- 52.0
        } else if (!is_out) {
          d$ssurgo_rock_outcrop_pct[i] <- 0.0
        }
      }
    }
    return(d)
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

aggregate_ssurgo_tabular_to_points <- function(sub_pts, tab_res, classify_fn = NULL, is_demo = FALSE) {
  if (is.null(classify_fn) || !is.function(classify_fn)) {
    if (exists("classify_ssurgo_outcrop", mode = "function")) {
      classify_fn <- get("classify_ssurgo_outcrop", mode = "function")
    } else {
      classify_fn <- function(
        muname = "", compname = "", taxsuborder = "", taxgreatgroup = "", compkind = "",
        lat = NA, lon = NA, is_demo = FALSE, demo = NULL, reskind = "",
        soil_ph = NA_real_, soil_sand = NA_real_, soil_cec = NA_real_,
        bedrock_depth_cm = NA_real_, rock_outcrop_pct = NA_real_
      ) {
        desc <- tolower(paste(muname, compname, taxsuborder, taxgreatgroup, compkind, reskind))
        is_out <- grepl(SSURGO_OUTCROP_SERIES_REGEX, paste(muname, compname), perl = TRUE) ||
                  grepl("rock outcrop|flatrock|glade|lithic", desc) ||
                  (is.finite(rock_outcrop_pct) && rock_outcrop_pct >= 15.0)
        ocl <- if (grepl("gladeville|barfield|cedars|talbott|limestone|dolomite|dolostone|calcareous|alvar|karst|rendolls|glade", desc)) {
          "Limestone glade"
        } else if (grepl("ramsey|lily|sewanee|goldston|manteo|wallen|dekab|sandstone|quartzite|siliceous|chert|shale", desc)) {
          "Sandstone/Siliceous"
        } else if (grepl("chrome|trego|conowingo|aldino|serpentine|ultramafic|peridotite", desc)) {
          "Serpentine"
        } else if (grepl("iredell|mecklenburg|wilkes|enon|poindexter|diabase|gabbro|basalt|amphibolite|mafic", desc)) {
          "Diabase/Mafic glade"
        } else if (grepl("wake|louisburg|ashlar|rion|wateree|granit|gneiss|felsic|flatrock", desc)) {
          "Granitic flatrock"
        } else if (is_out) {
          if (is.finite(soil_ph) && soil_ph >= 6.8) {
            "Limestone glade"
          } else if (is.finite(soil_cec) && soil_cec >= 20.0 && (!is.finite(soil_sand) || soil_sand < 50.0)) {
            "Diabase/Mafic glade"
          } else if (is.finite(soil_sand) && soil_sand >= 60.0 && (!is.finite(soil_ph) || soil_ph <= 5.5)) {
            "Granitic flatrock"
          } else {
            "Sandstone/Siliceous"
          }
        } else {
          "Non-outcrop matrix"
        }
        list(is_rock_outcrop = is_out, outcrop_class = ocl)
      }
    }
  }
  fn_formals <- names(formals(classify_fn))
  has_dots <- "..." %in% fn_formals

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

    reskind_val <- if ("reskind" %in% names(rows_i)) {
      rk_v <- stats::na.omit(unique(trimws(as.character(rows_i$reskind))))
      rk_v <- rk_v[nzchar(rk_v)]
      if (length(rk_v) > 0L) paste(rk_v, collapse = "; ") else ""
    } else ""

    call_args <- list(
      muname = ifelse(!is.na(rep_row$muname[1]), as.character(rep_row$muname[1]), ""),
      compname = ifelse(!is.na(rep_row$compname[1]), as.character(rep_row$compname[1]), ""),
      taxsuborder = ifelse(!is.na(rep_row$taxsuborder[1]), as.character(rep_row$taxsuborder[1]), ""),
      taxgreatgroup = ifelse(!is.na(rep_row$taxgreatgroup[1]), as.character(rep_row$taxgreatgroup[1]), ""),
      lat = lat_i,
      lon = lon_i
    )
    extra_args <- list(
      compkind = if ("compkind" %in% names(rep_row) && !is.na(rep_row$compkind[1])) as.character(rep_row$compkind[1]) else "",
      is_demo = isTRUE(is_demo),
      reskind = reskind_val,
      soil_ph = ph_val,
      soil_sand = sand_pct,
      soil_cec = cec_val,
      bedrock_depth_cm = bedrock_depth_cm,
      rock_outcrop_pct = rock_outcrop_pct
    )
    for (arg_nm in names(extra_args)) {
      if (has_dots || (arg_nm %in% fn_formals)) {
        call_args[[arg_nm]] <- extra_args[[arg_nm]]
      }
    }
    cl <- do.call(classify_fn, call_args)
    is_outcrop_final <- isTRUE(cl$is_rock_outcrop) || (rock_outcrop_pct >= 15.0)
    outcrop_cls_final <- if (is_outcrop_final && cl$outcrop_class == "Non-outcrop matrix") {
      if (is.finite(ph_val) && ph_val >= 6.8) {
        "Limestone glade"
      } else if (is.finite(cec_val) && cec_val >= 20.0 && (!is.finite(sand_pct) || sand_pct < 50.0)) {
        "Diabase/Mafic glade"
      } else if (is.finite(sand_pct) && sand_pct >= 60.0 && (!is.finite(ph_val) || ph_val <= 5.5)) {
        "Granitic flatrock"
      } else {
        "Sandstone/Siliceous"
      }
    } else {
      cl$outcrop_class
    }

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
# 4. Harmonized Cross-Border 0-15 cm Rhizosphere Edaphics
#    (30 m POLARIS / 100 m SOLUS100 / Measured SSURGO / 250 m SoilGrids)
# ------------------------------------------------------------------------------
.normalize_soilgrids_prop <- function(v_raw, prop_type = "pct") {
  v <- suppressWarnings(as.numeric(v_raw))
  fin <- v[is.finite(v)]
  if (length(fin) == 0L) return(rep(NA_real_, length(v)))
  med_v <- stats::median(fin, na.rm = TRUE)
  if (identical(prop_type, "ph")) {
    if (is.finite(med_v) && med_v > 100.0) v <- v / 100.0
    else if (is.finite(med_v) && med_v > 14.5) v <- v / 10.0
    v[!is.finite(v) | v < 2.0 | v > 11.5] <- NA_real_
  } else if (identical(prop_type, "pct")) {
    # ISRIC SoilGrids v2.0 stores sand/clay/silt in g/kg (0..1000) and cfvo in cm3/dm3*10 (0..1000)
    if (is.finite(med_v) && med_v > 105.0) v <- v / 10.0
    v[!is.finite(v) | v < 0.0 | v > 100.0] <- NA_real_
  } else if (identical(prop_type, "cec")) {
    # ISRIC SoilGrids v2.0 stores cec in mmol(c)/kg (= cmol(+)/kg * 10); SOLUS100 cec7 is also x10
    if (is.finite(med_v) && med_v > 90.0) v <- v / 10.0
    v[!is.finite(v) | v < 0.0 | v > 250.0] <- NA_real_
  } else if (identical(prop_type, "bd")) {
    # ISRIC SoilGrids v2.0 stores bdod in cg/cm3 (= g/cm3 * 100)
    if (is.finite(med_v) && med_v > 10.0) v <- v / 100.0
    v[!is.finite(v) | v < 0.1 | v > 2.65] <- NA_real_
  }
  v
}

extract_polaris_and_solus_predictors <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Extracting harmonized 0-15cm rhizosphere soil & lithic predictors (POLARIS 30m / SOLUS100 / SSURGO / SoilGrids 250m)...")
  df <- flag_coordinate_uncertainty(coords_df)
  n <- nrow(df)
  if (n == 0L) return(df)

  polaris_dir <- file.path(env_dir, "polaris_30m")
  gssurgo_dir <- file.path(env_dir, "gssurgo_30m")
  solus_dir   <- file.path(env_dir, "solus_100m")
  if (!isTRUE(demo)) {
    dir.create(polaris_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(solus_dir, recursive = TRUE, showWarnings = FALSE)
  }

  polaris_cols <- c(
    "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m",
    "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m"
  )
  solus_cols <- c(
    "solus_resdept_cm", "solus_fragvol_pct",
    "solus_sand_pct", "solus_clay_pct", "solus_silt_pct", "solus_cec"
  )
  canon_soil_cols <- c(
    "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "soil_silt_pct",
    "ilr_texture_1", "ilr_texture_2", "fine_soil_cec", "fine_cfvo_pct", "fine_bedrock_depth_cm"
  )
  for (col in c(polaris_cols, solus_cols, canon_soil_cols)) {
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

  solus_ph_vec <- rep(NA_real_, n)
  if (!isTRUE(demo) && any(in_conus) && requireNamespace("terra", quietly = TRUE)) {
    conus_idx <- which(in_conus)
    pts_vect <- terra::vect(df[conus_idx, c("longitude", "latitude"), drop = FALSE],
                            geom = c("longitude", "latitude"), crs = "EPSG:4326")

    existing_polaris <- c(
      list.files(polaris_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE),
      if (dir.exists(gssurgo_dir)) list.files(gssurgo_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE) else character(0)
    )
    if (length(existing_polaris) > 0L) {
      extract_tile_set <- function(files_vec) {
        if (length(files_vec) == 0L) return(rep(NA_real_, length(conus_idx)))
        vals_mat <- matrix(NA_real_, nrow = length(conus_idx), ncol = length(files_vec))
        for (f_i in seq_along(files_vec)) {
          vals_mat[, f_i] <- tryCatch({
            r_tile <- terra::rast(files_vec[f_i])
            pts_tile <- .project_vect_to_rast(pts_vect, r_tile)
            as.numeric(terra::extract(r_tile, pts_tile)[[2]])
          }, error = function(e) rep(NA_real_, length(conus_idx)))
        }
        has_v <- rowSums(is.finite(vals_mat)) > 0L
        out_v <- rep(NA_real_, length(conus_idx))
        if (any(has_v)) out_v[has_v] <- rowMeans(vals_mat[has_v, , drop = FALSE], na.rm = TRUE)
        out_v
      }

      for (pv in polaris_vars) {
        pv_files <- existing_polaris[grepl(paste0("^", pv, "_"), basename(existing_polaris), ignore.case = TRUE)]
        if (length(pv_files) > 0L) {
          files_0_5  <- pv_files[grepl("0_5|0-5", basename(pv_files), ignore.case = TRUE)]
          files_5_15 <- pv_files[grepl("5_15|5-15", basename(pv_files), ignore.case = TRUE)]
          v_rhizo <- if (length(files_0_5) > 0L || length(files_5_15) > 0L) {
            compute_rhizosphere_0_15cm(extract_tile_set(files_0_5), extract_tile_set(files_5_15))
          } else {
            extract_tile_set(pv_files)
          }
          ok_v <- is.finite(v_rhizo)
          if (any(ok_v)) df[[polaris_col_map[[pv]]]][conus_idx[ok_v]] <- v_rhizo[ok_v]
        }
      }
    }

    existing_solus <- c(
      list.files(solus_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE),
      if (dir.exists(gssurgo_dir)) list.files(gssurgo_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE) else character(0)
    )
    if (length(existing_solus) > 0L) {
      ext_solus_1 <- function(r_path) {
        tryCatch({
          r_sol <- terra::rast(r_path)
          pts_proj <- .project_vect_to_rast(pts_vect, r_sol)
          as.numeric(terra::extract(r_sol, pts_proj)[[2]])
        }, error = function(e) rep(NA_real_, length(conus_idx)))
      }
      ext_solus_0_15 <- function(cand_files, transform_fn = identity) {
        if (length(cand_files) == 0L) return(rep(NA_real_, length(conus_idx)))
        f_0_5  <- cand_files[grepl("0_5|0-5", basename(cand_files), ignore.case = TRUE)]
        f_5_15 <- cand_files[grepl("5_15|5-15", basename(cand_files), ignore.case = TRUE)]
        f_0cm  <- cand_files[grepl("(^|_)0_cm(_|\\.)", basename(cand_files), ignore.case = TRUE)]
        f_5cm  <- cand_files[grepl("(^|_)5_cm(_|\\.)", basename(cand_files), ignore.case = TRUE)]
        f_15cm <- cand_files[grepl("(^|_)15_cm(_|\\.)", basename(cand_files), ignore.case = TRUE)]
        if (length(f_0_5) > 0L && length(f_5_15) > 0L) {
          compute_rhizosphere_0_15cm(transform_fn(ext_solus_1(f_0_5[1L])), transform_fn(ext_solus_1(f_5_15[1L])))
        } else if (length(f_0cm) > 0L && length(f_5cm) > 0L && length(f_15cm) > 0L) {
          x0  <- transform_fn(ext_solus_1(f_0cm[1L]))
          x5  <- transform_fn(ext_solus_1(f_5cm[1L]))
          x15 <- transform_fn(ext_solus_1(f_15cm[1L]))
          compute_rhizosphere_0_15cm(0.5 * (x0 + x5), 0.5 * (x5 + x15))
        } else {
          transform_fn(ext_solus_1(cand_files[1L]))
        }
      }

      resdept_files <- existing_solus[grepl("anylithicdpt|resdept", basename(existing_solus), ignore.case = TRUE)]
      if (length(resdept_files) > 0L) {
        v_res <- rep(Inf, length(conus_idx))
        for (rf_p in resdept_files) {
          vr_i <- ext_solus_1(rf_p)
          ok_i <- is.finite(vr_i) & vr_i >= 0.0
          v_res[ok_i] <- pmin(v_res[ok_i], vr_i[ok_i])
        }
        v_res[!is.finite(v_res)] <- NA_real_
        ok_r <- is.finite(v_res)
        if (any(ok_r)) df$solus_resdept_cm[conus_idx[ok_r]] <- round(pmin(200.0, pmax(0.0, v_res[ok_r])), 1)
      }

      v_frg <- ext_solus_0_15(existing_solus[grepl("fragvol|cfvo", basename(existing_solus), ignore.case = TRUE)])
      ok_f <- is.finite(v_frg) & v_frg >= 0.0
      if (any(ok_f)) df$solus_fragvol_pct[conus_idx[ok_f]] <- round(pmin(100.0, pmax(0.0, v_frg[ok_f])), 1)

      v_phs <- ext_solus_0_15(existing_solus[grepl("^ph1to1h2o|^ph_solus", basename(existing_solus), ignore.case = TRUE)], .rescale_solus_ph_vec)
      ok_p <- is.finite(v_phs)
      if (any(ok_p)) solus_ph_vec[conus_idx[ok_p]] <- round(v_phs[ok_p], 2)

      v_snd <- ext_solus_0_15(existing_solus[grepl("^sandtotal", basename(existing_solus), ignore.case = TRUE)], function(z) .normalize_soilgrids_prop(z, "pct"))
      ok_sn <- is.finite(v_snd)
      if (any(ok_sn)) df$solus_sand_pct[conus_idx[ok_sn]] <- round(v_snd[ok_sn], 1)

      v_cly <- ext_solus_0_15(existing_solus[grepl("^claytotal", basename(existing_solus), ignore.case = TRUE)], function(z) .normalize_soilgrids_prop(z, "pct"))
      ok_cl <- is.finite(v_cly)
      if (any(ok_cl)) df$solus_clay_pct[conus_idx[ok_cl]] <- round(v_cly[ok_cl], 1)

      v_slt <- ext_solus_0_15(existing_solus[grepl("^silttotal", basename(existing_solus), ignore.case = TRUE)], function(z) .normalize_soilgrids_prop(z, "pct"))
      ok_si <- is.finite(v_slt)
      if (any(ok_si)) df$solus_silt_pct[conus_idx[ok_si]] <- round(v_slt[ok_si], 1)

      v_cec <- ext_solus_0_15(existing_solus[grepl("^cec7", basename(existing_solus), ignore.case = TRUE)], function(z) .normalize_soilgrids_prop(z / 10.0, "cec"))
      ok_ce <- is.finite(v_cec)
      if (any(ok_ce)) df$solus_cec[conus_idx[ok_ce]] <- round(v_cec[ok_ce], 2)
    }
  }

  # SoilGrids 250 m cross-border 0-15 cm rhizosphere columns (strictly preserving NA_real_ when missing in production)
  get_sg_0_15 <- function(col_0_5, col_5_15, col_base, prop_type) {
    x05 <- if (col_0_5 %in% names(df)) suppressWarnings(as.numeric(df[[col_0_5]])) else rep(NA_real_, n)
    x515 <- if (col_5_15 %in% names(df)) suppressWarnings(as.numeric(df[[col_5_15]])) else if (col_base %in% names(df)) suppressWarnings(as.numeric(df[[col_base]])) else rep(NA_real_, n)
    .normalize_soilgrids_prop(compute_rhizosphere_0_15cm(x05, x515), prop_type = prop_type)
  }
  sg_ph   <- get_sg_0_15("soil_ph_0_5cm",   "soil_ph_5_15cm",   "soil_ph",   "ph")
  sg_sand <- get_sg_0_15("soil_sand_0_5cm", "soil_sand_5_15cm", "soil_sand", "pct")
  sg_clay <- get_sg_0_15("soil_clay_0_5cm", "soil_clay_5_15cm", "soil_clay", "pct")
  sg_silt <- get_sg_0_15("soil_silt_0_5cm", "soil_silt_5_15cm", "soil_silt", "pct")
  sg_cec  <- get_sg_0_15("soil_cec_0_5cm",  "soil_cec_5_15cm",  "soil_cec",  "cec")
  sg_cfvo <- get_sg_0_15("soil_cfvo_0_5cm", "soil_cfvo_5_15cm", "soil_cfvo", "pct")
  sg_bd   <- get_sg_0_15("soil_bd_0_5cm",   "soil_bd_5_15cm",   "soil_bulk_density", "bd")

  ss_ph   <- if ("ssurgo_ph" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_ph)) else rep(NA_real_, n)
  ss_sand <- if ("ssurgo_sand_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_sand_pct)) else rep(NA_real_, n)
  ss_clay <- if ("ssurgo_clay_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_clay_pct)) else rep(NA_real_, n)
  ss_cec  <- if ("ssurgo_cec" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_cec)) else rep(NA_real_, n)
  ss_bed  <- if ("ssurgo_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_bedrock_depth_cm)) else rep(NA_real_, n)
  ss_out  <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  is_out  <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  ss_out[is.na(ss_out)] <- ifelse(is_out[is.na(ss_out)], 40.0, 0.0)

  if (isTRUE(demo)) {
    for (i in seq_len(n)) {
      p_b <- if (is.finite(ss_ph[i])) ss_ph[i] else if (is.finite(sg_ph[i])) sg_ph[i] else 5.6
      s_b <- if (is.finite(ss_sand[i])) ss_sand[i] else if (is.finite(sg_sand[i])) sg_sand[i] else 50.0
      c_b <- if (is.finite(ss_clay[i])) ss_clay[i] else if (is.finite(sg_clay[i])) sg_clay[i] else pmax(5.0, pmin(55.0, (100.0 - s_b) * 0.42))
      if (is.na(df$polaris_ph_30m[i])) df$polaris_ph_30m[i] <- round(compute_rhizosphere_0_15cm(p_b - 0.06, p_b + 0.03), 2)
      if (is.na(df$polaris_sand_30m[i])) df$polaris_sand_30m[i] <- round(compute_rhizosphere_0_15cm(s_b + 1.2, s_b - 0.6), 1)
      if (is.na(df$polaris_clay_30m[i])) df$polaris_clay_30m[i] <- round(compute_rhizosphere_0_15cm(pmax(2.0, c_b - 1.4), c_b + 0.7), 1)
      if (is.na(df$solus_resdept_cm[i])) df$solus_resdept_cm[i] <- if (is.finite(ss_bed[i])) round(ss_bed[i], 1) else 160.0
      if (is.na(df$solus_fragvol_pct[i])) {
        df$solus_fragvol_pct[i] <- round(pmin(85.0, pmax(1.0, 0.55 * ss_out[i] + ifelse(df$solus_resdept_cm[i] < 100, (100 - df$solus_resdept_cm[i]) * 0.25, 4.0))), 1)
      }
      if (is.na(df$solus_cec[i])) df$solus_cec[i] <- if (is.finite(ss_cec[i])) round(ss_cec[i], 2) else if (is.finite(sg_cec[i])) round(sg_cec[i], 2) else 11.5
    }
  }

  coalesce_num <- function(...) {
    vecs <- list(...)
    out <- rep(NA_real_, n)
    for (v in vecs) {
      vn <- suppressWarnings(as.numeric(v))
      fill <- !is.finite(out) & is.finite(vn)
      if (any(fill)) out[fill] <- vn[fill]
    }
    out
  }

  # Strict non-blended cross-border coalesce hierarchy:
  # CONUS fine-grain rasters (POLARIS 30m -> SOLUS100 100m) -> measured SSURGO -> global SoilGrids 250m -> NA_real_
  df$fine_soil_ph   <- round(coalesce_num(df$polaris_ph_30m, solus_ph_vec, ss_ph, sg_ph), 2)
  df$fine_soil_sand <- round(coalesce_num(df$polaris_sand_30m, df$solus_sand_pct, ss_sand, sg_sand), 1)
  df$fine_soil_clay <- round(coalesce_num(df$polaris_clay_30m, df$solus_clay_pct, ss_clay, sg_clay), 1)
  silt_obs          <- coalesce_num(df$solus_silt_pct, sg_silt)
  df$fine_soil_cec  <- round(coalesce_num(df$solus_cec, ss_cec, sg_cec), 2)
  df$fine_cfvo_pct  <- round(pmin(100.0, pmax(0.0, coalesce_num(df$solus_fragvol_pct, sg_cfvo))), 1)

  bed_comb <- ifelse(
    is.finite(ss_bed) & is.finite(df$solus_resdept_cm),
    pmin(ss_bed, df$solus_resdept_cm),
    coalesce_num(df$solus_resdept_cm, ss_bed)
  )
  df$fine_bedrock_depth_cm <- ifelse(is.finite(bed_comb), round(pmin(200.0, pmax(0.0, bed_comb)), 1), NA_real_)

  tex_lr <- compute_soil_texture_log_ratios(df$fine_soil_sand, df$fine_soil_clay, silt_pct = silt_obs)
  df$soil_silt_pct <- tex_lr$soil_silt_pct
  df$ilr_texture_1 <- tex_lr$ilr_texture_1
  df$ilr_texture_2 <- tex_lr$ilr_texture_2

  return(df)
}

# ------------------------------------------------------------------------------
# 5. 30 m 3DEP Micro-Topography: DEM Slope, TPI, McCune & Keon (2002) HLI & D8 TWI
# ------------------------------------------------------------------------------
.compute_mccune_keon_hli <- function(lat_deg, slope_deg, aspect_deg) {
  # McCune & Keon (2002) J. Veg. Sci. 13:603-606 Eq. 3 (30-60 deg N potential direct incident radiation / heat load)
  rad <- pi / 180.0
  L <- pmin(60.0, pmax(15.0, abs(as.numeric(lat_deg)))) * rad
  S <- pmin(85.0, pmax(0.0, as.numeric(slope_deg))) * rad
  A_deg <- as.numeric(aspect_deg) %% 360.0
  # Folded aspect around SW (225 deg) in radians on [0, pi]
  F_rad <- abs(pi - abs(A_deg * rad - (225.0 * rad)))
  ln_hli <- -1.467 +
            1.582 * cos(L) * cos(S) -
            1.500 * cos(F_rad) * sin(S) * sin(L) -
            0.262 * sin(L) * sin(S) +
            0.607 * sin(F_rad) * sin(S)
  hli <- exp(ln_hli)
  ok <- is.finite(lat_deg) & is.finite(slope_deg)
  # Flat cells (slope == 0) have undefined aspect; set F_rad = pi/2 neutral aspect
  flat <- ok & (!is.finite(aspect_deg) | slope_deg < 1e-4)
  if (any(flat)) {
    Lf <- L[flat]
    hli[flat] <- exp(-1.467 + 1.582 * cos(Lf))
  }
  hli[!ok] <- NA_real_
  pmin(1.50, pmax(0.05, hli))
}

extract_micro_topography_3dep <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Computing 30m 3DEP micro-topographic covariates (Slope, TPI, HLI, D8 Hydrological TWI)...")
  df <- as.data.frame(coords_df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0L) return(df)

  dem_dir <- file.path(env_dir, "dem_3dep_30m")
  if (!isTRUE(demo)) {
    dir.create(dem_dir, recursive = TRUE, showWarnings = FALSE)
  }

  for (col in c("topo_tpi_30m", "topo_slope_deg", "topo_hli_30m", "topo_twi_30m", "elevation_m")) {
    if (!(col %in% names(df))) df[[col]] <- rep(NA_real_, n)
  }

  lat_v <- suppressWarnings(as.numeric(df$latitude))
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  valid_pts <- is.finite(lat_v) & is.finite(lon_v)

  stale_regional_dem <- file.path(dem_dir, "dem_3dep_regional.tif")
  if (!isTRUE(demo) && file.exists(stale_regional_dem)) {
    tryCatch(unlink(stale_regional_dem), error = function(e) NULL)
  }

  if (!isTRUE(demo) && any(valid_pts) && requireNamespace("terra", quietly = TRUE)) {
    dem_files <- .list_real_fine_rasters(dem_dir, exclude_pattern = "^dem_3dep_regional\\.tif$")
    if (length(dem_files) > 0L) {
      v_idx <- which(valid_pts)
      stems_up <- toupper(sub("\\.(tif|tiff|hgt)$", "", basename(dem_files), ignore.case = TRUE))
      is_grid_tile <- grepl("^[NS][0-9]{2}[EW][0-9]{3}$", stems_up)

      extract_from_single_dem <- function(r_dem, sub_rows_idx) {
        if (is.null(r_dem) || length(sub_rows_idx) == 0L) return(invisible(NULL))
        tryCatch({
          pts_sub <- terra::vect(df[sub_rows_idx, c("longitude", "latitude"), drop = FALSE],
                                 geom = c("longitude", "latitude"), crs = "EPSG:4326")
          pts_proj <- .project_vect_to_rast(pts_sub, r_dem)
          ext_pts <- terra::ext(pts_proj)
          res_xy <- terra::res(r_dem)
          is_lonlat <- isTRUE(tryCatch(terra::is.lonlat(r_dem), error = function(e) TRUE))
          # Pad >= 65 cells (~2 km on 30 m DEM) so D8 upslope catchment accumulation is physically grounded
          pad_x <- max(res_xy[1] * 65.0, if (is_lonlat) 0.018 else 1800.0)
          pad_y <- max(res_xy[2] * 65.0, if (is_lonlat) 0.018 else 1800.0)
          crop_win <- terra::ext(ext_pts[1] - pad_x, ext_pts[2] + pad_x, ext_pts[3] - pad_y, ext_pts[4] + pad_y)
          r_sub <- tryCatch(terra::crop(r_dem, crop_win), error = function(e) NULL)
          if (is.null(r_sub) || terra::ncell(r_sub) < 9L) return(invisible(NULL))

          focal_mean <- terra::focal(r_sub, w = 5, fun = mean, na.rm = TRUE)
          tpi_rast   <- r_sub - focal_mean
          slope_rast <- terra::terrain(r_sub, v = "slope", unit = "degrees")
          asp_rast   <- terra::terrain(r_sub, v = "aspect", unit = "degrees")
          fdir_rast  <- tryCatch(terra::terrain(r_sub, v = "flowdir"), error = function(e) NULL)
          facc_rast  <- if (!is.null(fdir_rast)) tryCatch(terra::flowAccumulation(fdir_rast), error = function(e) NULL) else NULL

          v_elev  <- as.numeric(terra::extract(r_sub, pts_proj)[[2]])
          v_tpi   <- as.numeric(terra::extract(tpi_rast, pts_proj)[[2]])
          v_slope <- as.numeric(terra::extract(slope_rast, pts_proj)[[2]])
          v_asp   <- as.numeric(terra::extract(asp_rast, pts_proj)[[2]])
          v_facc  <- if (!is.null(facc_rast)) as.numeric(terra::extract(facc_rast, pts_proj)[[2]]) else rep(NA_real_, length(sub_rows_idx))

          lat_sub <- lat_v[sub_rows_idx]
          cell_w_m <- if (is_lonlat) {
            dx_m <- res_xy[1] * 111320.0 * cos(pmin(80.0, abs(lat_sub)) * pi / 180.0)
            dy_m <- res_xy[2] * 110574.0
            sqrt(pmax(1.0, dx_m * dy_m))
          } else {
            sqrt(pmax(1.0, res_xy[1] * res_xy[2]))
          }

          ok_e <- is.finite(v_elev) & is.na(df$elevation_m[sub_rows_idx])
          if (any(ok_e)) df$elevation_m[sub_rows_idx[ok_e]] <<- round(v_elev[ok_e], 1)

          ok_t <- is.finite(v_tpi)
          if (any(ok_t)) df$topo_tpi_30m[sub_rows_idx[ok_t]] <<- round(v_tpi[ok_t], 2)

          ok_s <- is.finite(v_slope)
          if (any(ok_s)) {
            sl_clamped <- pmax(0.05, pmin(75.0, v_slope[ok_s]))
            df$topo_slope_deg[sub_rows_idx[ok_s]] <<- round(sl_clamped, 2)
            hli_sub <- .compute_mccune_keon_hli(lat_sub[ok_s], sl_clamped, v_asp[ok_s])
            ok_h <- is.finite(hli_sub)
            if (any(ok_h)) df$topo_hli_30m[sub_rows_idx[ok_s][ok_h]] <<- round(hli_sub[ok_h], 4)
          }

          ok_w <- is.finite(v_facc) & is.finite(v_slope)
          if (any(ok_w)) {
            a_spec <- (pmax(0.0, v_facc[ok_w]) + 1.0) * cell_w_m[ok_w]
            tan_b  <- tan(pmax(v_slope[ok_w], 0.10) * (pi / 180.0))
            twi_v  <- log(pmax(1e-3, a_spec / pmax(tan_b, 1e-4)))
            df$topo_twi_30m[sub_rows_idx[ok_w]] <<- round(pmin(25.0, pmax(0.5, twi_v)), 2)
          }
        }, error = function(e) NULL)
      }

      if (any(is_grid_tile)) {
        tile_map <- stats::setNames(dem_files[is_grid_tile], stems_up[is_grid_tile])
        lat_fl <- as.integer(floor(lat_v[v_idx]))
        lon_fl <- as.integer(floor(lon_v[v_idx]))
        pt_tiles <- sprintf("%s%02d%s%03d", ifelse(lat_fl >= 0L, "N", "S"), abs(lat_fl), ifelse(lon_fl >= 0L, "E", "W"), abs(lon_fl))
        for (tid in unique(pt_tiles)) {
          if (tid %in% names(tile_map)) {
            rows_tid <- v_idx[pt_tiles == tid]
            r_tile <- tryCatch(terra::rast(tile_map[[tid]]), error = function(e) NULL)
            extract_from_single_dem(r_tile, rows_tid)
          }
        }
      }

      rem_idx <- v_idx[is.na(df$topo_tpi_30m[v_idx])]
      if (length(rem_idx) > 0L) {
        cand_fb_dems <- c(dem_files[!is_grid_tile], dem_files[is_grid_tile])
        for (od_p in cand_fb_dems) {
          r_od <- tryCatch(terra::rast(od_p), error = function(e) NULL)
          if (is.null(r_od)) next
          ex_od <- tryCatch(terra::ext(r_od), error = function(e) NULL)
          if (!is.null(ex_od)) {
            in_ex <- rem_idx[lon_v[rem_idx] >= ex_od[1] & lon_v[rem_idx] <= ex_od[2] &
                             lat_v[rem_idx] >= ex_od[3] & lat_v[rem_idx] <= ex_od[4]]
            if (length(in_ex) > 0L) extract_from_single_dem(r_od, in_ex)
          } else {
            extract_from_single_dem(r_od, rem_idx)
          }
          rem_idx <- v_idx[is.na(df$topo_tpi_30m[v_idx])]
          if (length(rem_idx) == 0L) break
        }
      }
    }
  }

  muname_v <- if ("muname" %in% names(df)) as.character(df$muname) else rep("", n)
  slope_pct_ssurgo <- parse_ssurgo_slope_pct(muname_v)
  slope_deg_ssurgo <- atan(slope_pct_ssurgo / 100.0) * (180.0 / pi)

  if (!isTRUE(demo)) {
    # Use SSURGO mapunit slope ONLY as a fallback when DEM slope is unavailable (never 50/50 blending)
    miss_s <- !is.finite(df$topo_slope_deg) & is.finite(slope_deg_ssurgo)
    if (any(miss_s)) {
      df$topo_slope_deg[miss_s] <- round(pmax(0.15, pmin(60.0, slope_deg_ssurgo[miss_s])), 2)
    }
    return(df)
  }

  # Explicit demonstration simulation mode ONLY (demo == TRUE)
  is_out <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  out_pct <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  out_pct[is.na(out_pct)] <- ifelse(is_out[is.na(out_pct)], 40.0, 0.0)
  for (i in seq_len(n)) {
    s_ss <- slope_deg_ssurgo[i]
    df$topo_slope_deg[i] <- pmax(0.15, pmin(45.0, round(if (!is.na(s_ss)) s_ss else ifelse(is_out[i], 4.2, 2.5), 2)))
    is_concave <- grepl("muck|flooded|aquolls|aquepts|aquults|saprists|depression|seep", tolower(muname_v[i]))
    df$topo_tpi_30m[i] <- round(if (is_out[i] || out_pct[i] >= 15.0) 2.8 + 0.045 * out_pct[i] else if (is_concave) -2.4 else 0.0, 2)
    lat_i <- if (is.finite(lat_v[i])) lat_v[i] else 35.5
    asp_sim <- if (is_out[i]) 215.0 else 135.0
    df$topo_hli_30m[i] <- round(.compute_mccune_keon_hli(lat_i, df$topo_slope_deg[i], asp_sim), 4)
    tan_b <- tan(pmax(df$topo_slope_deg[i], 0.15) * (pi / 180.0))
    a_sim <- if (is_concave) 450.0 else if (is_out[i]) 45.0 else 120.0
    df$topo_twi_30m[i] <- round(log(a_sim / tan_b), 2)
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
  if (n == 0L) return(df)

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
      c_raw <- if (requireNamespace("data.table", quietly = TRUE)) {
        as.data.frame(data.table::fread(cache_csv_path, showProgress = FALSE, na.strings = c("", "NA")), stringsAsFactors = FALSE)
      } else {
        as.data.frame(readr::read_csv(cache_csv_path, show_col_types = FALSE), stringsAsFactors = FALSE)
      }
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
    cached_df <- if (!is.null(cached_df) && nrow(cached_df) > 0L) {
      if (exists(".bind_rows_schema_safe", mode = "function")) {
        .bind_rows_schema_safe(sql_cached, cached_df)
      } else {
        for (.cn in intersect(names(sql_cached), names(cached_df))) {
          if (!identical(typeof(sql_cached[[.cn]]), typeof(cached_df[[.cn]]))) {
            if (is.logical(cached_df[[.cn]]) && (is.numeric(sql_cached[[.cn]]) || is.integer(sql_cached[[.cn]]))) {
              sql_cached[[.cn]] <- ifelse(is.na(sql_cached[[.cn]]), NA, sql_cached[[.cn]] != 0)
            } else if ((is.numeric(sql_cached[[.cn]]) || is.integer(sql_cached[[.cn]])) && (is.numeric(cached_df[[.cn]]) || is.integer(cached_df[[.cn]]))) {
              sql_cached[[.cn]] <- suppressWarnings(as.numeric(sql_cached[[.cn]]))
              cached_df[[.cn]] <- suppressWarnings(as.numeric(cached_df[[.cn]]))
            } else {
              sql_cached[[.cn]] <- as.character(sql_cached[[.cn]])
              cached_df[[.cn]] <- as.character(cached_df[[.cn]])
            }
          }
        }
        dplyr::bind_rows(sql_cached, cached_df)
      }
    } else sql_cached
  }

  pol_disk_fls <- c(.list_real_fine_rasters(file.path(env_dir, "polaris_30m")), .list_real_fine_rasters(file.path(env_dir, "gssurgo_30m")))
  sol_disk_fls <- c(.list_real_fine_rasters(file.path(env_dir, "solus_100m")), .list_real_fine_rasters(file.path(env_dir, "gssurgo_30m")))
  has_pol_disk <- length(pol_disk_fls) > 0L
  has_sol_disk <- length(sol_disk_fls) > 0L
  has_sol_ph_disk <- any(grepl("^ph1to1h2o|^ph_", basename(sol_disk_fls), ignore.case = TRUE))
  stale_dem_p <- file.path(env_dir, "dem_3dep_30m", "dem_3dep_regional.tif")
  if (!isTRUE(demo) && file.exists(stale_dem_p)) tryCatch(unlink(stale_dem_p), error = function(e) NULL)
  has_dem_disk <- length(.list_real_fine_rasters(file.path(env_dir, "dem_3dep_30m"), exclude_pattern = "^dem_3dep_regional\\.tif$")) > 0L

  if (!is.null(cached_df) && nrow(cached_df) > 0L && any(c("fine_soil_ph", "polaris_ph_30m", "solus_resdept_cm") %in% names(cached_df))) {
    avail_cache_cols <- intersect(
      setdiff(FINE_GRAIN_CACHE_COLS, c("catalogNumber", "lat_lon_key", "latitude", "longitude",
                                       "buffer_radius_m", "high_coord_uncertainty_flag", "is_training_eligible_coord",
                                       "is_clim_eligible_coord", "is_edaphic_eligible_coord", "edaphic_coord_weight")),
      names(cached_df)
    )
    c_by_cat <- cached_df[!duplicated(cached_df$catalogNumber), c("catalogNumber", avail_cache_cols), drop = FALSE]
    m_idx <- match(df$catalogNumber, c_by_cat$catalogNumber)

    miss_cat <- which(is.na(m_idx) & !is.na(df$lat_lon_key))
    if (length(miss_cat) > 0L) {
      c_by_ll <- cached_df[!is.na(cached_df$lat_lon_key) & !duplicated(cached_df$lat_lon_key), c("lat_lon_key", avail_cache_cols), drop = FALSE]
      ll_idx <- match(df$lat_lon_key[miss_cat], c_by_ll$lat_lon_key)
      for (col in avail_cache_cols) {
        df[[col]] <- c_by_cat[[col]][m_idx]
        hit_ll <- !is.na(ll_idx)
        if (any(hit_ll)) df[[col]][miss_cat[hit_ll]] <- c_by_ll[[col]][ll_idx[hit_ll]]
      }
    } else {
      for (col in avail_cache_cols) df[[col]] <- c_by_cat[[col]][m_idx]
    }

    uncached_mask <- is.na(df$fine_soil_ph) |
                     !(all(c("ilr_texture_1", "ilr_texture_2") %in% names(df))) |
                     (has_dem_disk & (!("topo_hli_30m" %in% names(df)) | is.na(df$topo_hli_30m) | is.na(df$topo_twi_30m)))
    n_uncached <- sum(uncached_mask)
    message(sprintf(
      "Loaded %d / %d fine-grain edaphic records from %s (%d coordinates requiring refresh).",
      n - n_uncached, n, cache_csv_path, n_uncached
    ))

    if (n_uncached > 0L) {
      sub_df <- df[uncached_mask, , drop = FALSE]
      for (col in req_core_cols) sub_df[[col]] <- NULL
      if (!is.null(query_ssurgo_fn)) {
        sub_df <- query_ssurgo_fn(sub_df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
      }
      sub_df <- enrich_ssurgo_pedology_from_series(sub_df, demo = demo)
      sub_df <- extract_polaris_and_solus_predictors(sub_df, env_dir = env_dir, demo = demo)
      sub_df <- extract_micro_topography_3dep(sub_df, env_dir = env_dir, demo = demo)
      for (col in setdiff(names(sub_df), c("catalogNumber", "lat_lon_key", "latitude", "longitude"))) {
        if (!(col %in% names(df))) df[[col]] <- rep(NA, n)
        df[[col]][uncached_mask] <- sub_df[[col]]
      }
    }

    need_sgmc <- !("sgmc_lithology_class" %in% names(df)) | is.na(df$sgmc_lithology_class) |
                 !nzchar(trimws(as.character(df$sgmc_lithology_class)))
    if (any(need_sgmc)) {
      g_sub <- extract_sgmc_lithology(df[need_sgmc, , drop = FALSE], env_dir = env_dir, demo = demo)
      for (col in c("sgmc_lithology_class", "is_azonal_lithology")) {
        if (col %in% names(g_sub)) {
          if (!(col %in% names(df))) df[[col]] <- rep(NA, n)
          df[[col]][need_sgmc] <- g_sub[[col]]
        }
      }
    }

    sanitized_stale <- FALSE
    if (!isTRUE(demo)) {
      if (!has_pol_disk) {
        pol_only_cols <- c("polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m")
        if (!has_sol_ph_disk) pol_only_cols <- c(pol_only_cols, "ph_q10_buffer", "ph_q90_buffer")
        pol_c <- intersect(pol_only_cols, names(df))
        if (any(vapply(pol_c, function(cc) any(!is.na(df[[cc]])), logical(1)))) {
          for (cc in pol_c) df[[cc]] <- rep(NA_real_, n)
          sanitized_stale <- TRUE
        }
      }
      if (!has_sol_disk) {
        sol_c <- intersect(c("solus_resdept_cm", "solus_fragvol_pct", "solus_sand_pct", "solus_clay_pct", "solus_silt_pct", "solus_cec", "min_bedrock_depth_buffer_cm", "max_cfvo_buffer_pct"), names(df))
        if (any(vapply(sol_c, function(cc) any(!is.na(df[[cc]])), logical(1)))) {
          for (cc in sol_c) df[[cc]] <- rep(NA_real_, n)
          sanitized_stale <- TRUE
        }
      }
      if (!has_dem_disk) {
        dem_c <- intersect(c("topo_tpi_30m", "topo_hli_30m", "topo_twi_30m"), names(df))
        if (any(vapply(dem_c, function(cc) any(!is.na(df[[cc]])), logical(1)))) {
          for (cc in dem_c) df[[cc]] <- rep(NA_real_, n)
          df <- extract_micro_topography_3dep(df, env_dir = env_dir, demo = FALSE)
          sanitized_stale <- TRUE
        }
      }
      if (isTRUE(sanitized_stale)) {
        df <- extract_polaris_and_solus_predictors(df, env_dir = env_dir, demo = FALSE)
        df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = FALSE)
      }
    }

    if (all(c("fine_soil_sand", "fine_soil_clay") %in% names(df))) {
      miss_ilr <- !("ilr_texture_1" %in% names(df)) | !("ilr_texture_2" %in% names(df)) |
                  (is.na(df$ilr_texture_1) & is.finite(suppressWarnings(as.numeric(df$fine_soil_sand))) & is.finite(suppressWarnings(as.numeric(df$fine_soil_clay))))
      if (any(miss_ilr)) {
        silt_in <- if ("soil_silt_pct" %in% names(df)) df$soil_silt_pct else NULL
        tex_hyd <- compute_soil_texture_log_ratios(df$fine_soil_sand, df$fine_soil_clay, silt_pct = silt_in)
        if (!("soil_silt_pct" %in% names(df))) df$soil_silt_pct <- tex_hyd$soil_silt_pct
        else df$soil_silt_pct[miss_ilr] <- tex_hyd$soil_silt_pct[miss_ilr]
        if (!("ilr_texture_1" %in% names(df))) df$ilr_texture_1 <- tex_hyd$ilr_texture_1
        else df$ilr_texture_1[miss_ilr] <- tex_hyd$ilr_texture_1[miss_ilr]
        if (!("ilr_texture_2" %in% names(df))) df$ilr_texture_2 <- tex_hyd$ilr_texture_2
        else df$ilr_texture_2[miss_ilr] <- tex_hyd$ilr_texture_2[miss_ilr]
      }
    }

    df$is_rock_outcrop <- as.logical(df$is_rock_outcrop)
    df$is_rock_outcrop[is.na(df$is_rock_outcrop)] <- FALSE
    if (!sanitized_stale && n_uncached == 0L && !any(need_sgmc)) {
      return(df)
    }
  } else {
    if (!is.null(query_ssurgo_fn)) {
      df <- query_ssurgo_fn(df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
    }
    df <- enrich_ssurgo_pedology_from_series(df, demo = demo)
    df <- extract_polaris_and_solus_predictors(df, env_dir = env_dir, demo = demo)
    df <- extract_micro_topography_3dep(df, env_dir = env_dir, demo = demo)
    df <- extract_sgmc_lithology(df, env_dir = env_dir, demo = demo)
    df <- extract_uncertainty_buffered_predictors(df, env_dir = env_dir, demo = demo)
  }

  if (!isTRUE(demo)) {
    export_cols <- intersect(FINE_GRAIN_CACHE_COLS, names(df))
    new_cache <- df[, export_cols, drop = FALSE]
    if (!is.null(cached_df) && nrow(cached_df) > 0L) {
      for (mc in setdiff(export_cols, names(cached_df))) cached_df[[mc]] <- NA
      common_c <- intersect(export_cols, names(cached_df))
      combined_cache <- rbind(new_cache[, common_c, drop = FALSE], cached_df[, common_c, drop = FALSE])
      combined_cache <- combined_cache[!duplicated(combined_cache$catalogNumber), , drop = FALSE]
    } else {
      combined_cache <- new_cache[!duplicated(new_cache$catalogNumber), , drop = FALSE]
    }
    dir.create(dirname(cache_csv_path), recursive = TRUE, showWarnings = FALSE)
    tmp_cache <- paste0(cache_csv_path, ".tmp")
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fwrite(combined_cache, file = tmp_cache, na = "", quote = "auto", logical01 = FALSE)
    } else {
      readr::write_csv(combined_cache, tmp_cache)
    }
    file.rename(tmp_cache, cache_csv_path)
    if (exists("upsert_environmental_cache_sqlite", mode = "function")) {
      tryCatch(upsert_environmental_cache_sqlite(combined_cache, cache_path = cache_csv_path, demo = demo), error = function(e) NULL)
    }
    message(sprintf("Persisted %d fine-grain edaphic & micro-topographic records to %s.", nrow(combined_cache), cache_csv_path))
  }

  return(df)
}

# ------------------------------------------------------------------------------
# 7. Collinearity Pruning Across Climatic, Edaphic & Topographic Blocks
#    (Dormann et al. 2013; Empirical Anchor Separation Tie-Breaking)
# ------------------------------------------------------------------------------
compute_anchor_separation_scores <- function(X_df, anchor_groups = NULL) {
  X_d <- as.data.frame(X_df, check.names = FALSE)
  cols <- colnames(X_d)
  if (length(cols) == 0L) return(numeric(0))
  n_rows <- nrow(X_d)
  has_grp <- !is.null(anchor_groups) && (length(anchor_groups) == n_rows)
  g_all <- if (has_grp) trimws(as.character(anchor_groups)) else character(0)

  scores <- vapply(cols, function(cl) {
    v <- suppressWarnings(as.numeric(X_d[[cl]]))
    fin <- is.finite(v)
    fin_frac <- if (n_rows > 0L) sum(fin) / n_rows else 0.0
    if (sum(fin) < 3L) return(0.0)
    sd_v <- stats::sd(v[fin])
    if (!is.finite(sd_v) || sd_v <= 1e-8) return(0.0)

    if (has_grp) {
      ok <- fin & !is.na(g_all) & nzchar(g_all) & !(tolower(g_all) %in% c("unknown", "na", "none"))
      if (sum(ok) >= 4L) {
        tbl_g <- table(g_all[ok])
        valid_g <- names(tbl_g[tbl_g >= 2L])
        if (length(valid_g) >= 2L) {
          m_g <- ok & (g_all %in% valid_g)
          v_g <- v[m_g]
          f_g <- factor(g_all[m_g])
          kw_stat <- tryCatch(
            suppressWarnings(as.numeric(stats::kruskal.test(v_g ~ f_g)$statistic[1])),
            error = function(e) NA_real_
          )
          f_stat <- tryCatch({
            af <- summary(stats::aov(v_g ~ f_g))[[1]]
            suppressWarnings(as.numeric(af[["F value"]][1]))
          }, error = function(e) NA_real_)
          if (is.finite(kw_stat) && kw_stat >= 0.0) {
            f_tie <- if (is.finite(f_stat) && f_stat >= 0.0) 1e-4 * f_stat else 0.0
            return(kw_stat + f_tie + 1e-6 * fin_frac)
          }
          if (is.finite(f_stat) && f_stat >= 0.0) {
            return(f_stat + 1e-6 * fin_frac)
          }
        }
      }
    }
    fin_frac
  }, numeric(1))
  names(scores) <- cols
  scores
}

screen_collinear_predictors <- function(
  X_env,
  cor.threshold = 0.70,
  vif.threshold = 5.0,
  anchor_groups = NULL,
  var_scores = NULL
) {
  X_df <- as.data.frame(X_env, check.names = FALSE)
  all_cols <- colnames(X_df)
  if (length(all_cols) == 0L) {
    return(list(
      selected_variables = character(0),
      kept_vars = character(0),
      pruned_variables = character(0),
      anchor_scores = numeric(0)
    ))
  }
  if (length(all_cols) == 1L) {
    return(list(
      selected_variables = all_cols,
      kept_vars = all_cols,
      pruned_variables = character(0),
      anchor_scores = stats::setNames(1.0, all_cols)
    ))
  }

  sds <- vapply(X_df, function(v) stats::sd(suppressWarnings(as.numeric(v)), na.rm = TRUE), numeric(1))
  keep_cols <- names(sds[is.finite(sds) & sds > 1e-8])
  if (length(keep_cols) == 0L) {
    return(list(
      selected_variables = all_cols[1],
      kept_vars = all_cols[1],
      pruned_variables = setdiff(all_cols, all_cols[1]),
      anchor_scores = stats::setNames(rep(0.0, length(all_cols)), all_cols)
    ))
  }
  if (length(keep_cols) == 1L) {
    return(list(
      selected_variables = keep_cols,
      kept_vars = keep_cols,
      pruned_variables = setdiff(all_cols, keep_cols),
      anchor_scores = stats::setNames(1.0, keep_cols)
    ))
  }

  scores_sub <- if (!is.null(var_scores) && all(keep_cols %in% names(var_scores))) {
    var_scores[keep_cols]
  } else {
    compute_anchor_separation_scores(X_df[, keep_cols, drop = FALSE], anchor_groups = anchor_groups)
  }
  scores_sub[!is.finite(scores_sub)] <- 0.0
  ranked_cols <- keep_cols[order(scores_sub[keep_cols], decreasing = TRUE)]

  X_imp <- as.data.frame(lapply(ranked_cols, function(cl) {
    v <- suppressWarnings(as.numeric(X_df[[cl]]))
    if (any(!is.finite(v))) {
      med_v <- stats::median(v[is.finite(v)], na.rm = TRUE)
      v[!is.finite(v)] <- if (is.finite(med_v)) med_v else 0.0
    }
    v
  }), stringsAsFactors = FALSE)
  colnames(X_imp) <- ranked_cols

  sel_vars <- NULL
  if (requireNamespace("spatialRF", quietly = TRUE) && ncol(X_imp) >= 3L) {
    ac <- tryCatch(
      spatialRF::auto_cor(x = X_imp, cor.threshold = cor.threshold, preference.order = ranked_cols),
      error = function(e) tryCatch(spatialRF::auto_cor(x = X_imp, cor.threshold = cor.threshold), error = function(e2) NULL)
    )
    X_ac <- if (!is.null(ac) && !is.null(ac$selected.variables) && length(ac$selected.variables) >= 1L) {
      X_imp[, ac$selected.variables, drop = FALSE]
    } else {
      X_imp
    }
    if (ncol(X_ac) >= 2L) {
      pref_ac <- intersect(ranked_cols, colnames(X_ac))
      av <- tryCatch(
        spatialRF::auto_vif(x = X_ac, vif.threshold = vif.threshold, preference.order = pref_ac),
        error = function(e) tryCatch(spatialRF::auto_vif(x = X_ac, vif.threshold = vif.threshold), error = function(e2) NULL)
      )
      if (!is.null(av) && !is.null(av$selected.variables) && length(av$selected.variables) >= 1L) {
        sel_vars <- av$selected.variables
      }
    } else if (ncol(X_ac) == 1L) {
      sel_vars <- colnames(X_ac)
    }
  }

  # Verify |r| <= cor.threshold and VIF <= vif.threshold via analytical Dormann et al. (2013)
  # tie-broken strictly by empirical anchor separation score (Kruskal-Wallis chi^2 / ANOVA F)
  cmat <- abs(stats::cor(X_imp, use = "pairwise.complete.obs"))
  cmat[!is.finite(cmat)] <- 0.0
  diag(cmat) <- 0.0
  needs_analytical <- is.null(sel_vars) || length(sel_vars) == 0L
  if (!needs_analytical && length(sel_vars) >= 2L) {
    max_r_sel <- max(cmat[sel_vars, sel_vars, drop = FALSE], na.rm = TRUE)
    if (is.finite(max_r_sel) && max_r_sel > (cor.threshold + 1e-6)) {
      needs_analytical <- TRUE
    }
  }

  if (needs_analytical) {
    cand <- ranked_cols
    while (length(cand) > 1L && max(cmat[cand, cand, drop = FALSE], na.rm = TRUE) > cor.threshold) {
      sub_c <- cmat[cand, cand, drop = FALSE]
      ij <- which(sub_c == max(sub_c, na.rm = TRUE), arr.ind = TRUE)[1, ]
      v1 <- cand[ij[1]]
      v2 <- cand[ij[2]]
      s1 <- scores_sub[[v1]]
      s2 <- scores_sub[[v2]]
      if (is.finite(s1) && is.finite(s2) && abs(s1 - s2) > 1e-9) {
        drop_v <- if (s1 < s2) v1 else v2
      } else {
        drop_v <- if (mean(sub_c[v1, ]) >= mean(sub_c[v2, ])) v1 else v2
      }
      cand <- setdiff(cand, drop_v)
    }
    while (length(cand) > 2L) {
      R_sub <- stats::cor(X_imp[, cand, drop = FALSE], use = "pairwise.complete.obs")
      R_sub[!is.finite(R_sub)] <- 0.0
      diag(R_sub) <- 1.0
      vifs <- tryCatch(diag(solve(R_sub + diag(1e-6, length(cand)))), error = function(e) rep(1.0, length(cand)))
      names(vifs) <- cand
      max_vif <- max(vifs, na.rm = TRUE)
      if (!is.finite(max_vif) || max_vif <= vif.threshold) break
      over_v <- names(vifs[vifs >= (max_vif - 1e-9)])
      drop_v <- if (length(over_v) > 1L) {
        over_v[order(scores_sub[over_v], decreasing = FALSE)][1]
      } else {
        over_v[1]
      }
      cand <- setdiff(cand, drop_v)
    }
    sel_vars <- cand
  }

  pruned_vars <- setdiff(all_cols, sel_vars)
  message(sprintf(
    "Collinearity Screening (Dormann et al. 2013; Benito 2021; |r| <= %.2f, VIF <= %.1f): retained %d / %d predictors.",
    cor.threshold, vif.threshold, length(sel_vars), ncol(X_df)
  ))
  list(
    selected_variables = sel_vars,
    kept_vars = sel_vars,
    pruned_variables = pruned_vars,
    anchor_scores = scores_sub
  )
}

screen_collinear_predictors_by_block <- function(
  X_env,
  clim_vars = NULL,
  edaph_vars = NULL,
  cor.threshold = 0.70,
  vif.threshold = 5.0,
  topo_vars = NULL,
  anchor_groups = NULL
) {
  # Applies spatialRF::auto_cor(cor.threshold = 0.70) and
  # spatialRF::auto_vif(vif.threshold = 5.0) separately within 3 ecological blocks:
  #   1. Climatic block (bio1..bio19 + phenology if passed)
  #   2. Edaphic-Lithic block (fine_* + ssurgo_* + polaris_* + solus_* + soil_* + lithic)
  #   3. Topographic block (topo_*)
  # Within-block tie-breaking ranks variables empirically by Kruskal-Wallis chi^2 /
  # ANOVA F separating active reference anchors (anchor_groups) in the current run.
  X_df <- as.data.frame(X_env, check.names = FALSE)
  all_cols <- colnames(X_df)
  if (length(all_cols) == 0L) {
    return(list(
      selected_variables = character(0),
      kept_vars = character(0),
      selected_edaphic = character(0),
      selected_climatic = character(0),
      selected_topo = character(0),
      edaphic_variables = character(0),
      climatic_variables = character(0),
      topo_variables = character(0),
      pruned_variables = character(0),
      anchor_scores = numeric(0)
    ))
  }

  if (is.null(clim_vars)) {
    clim_vars <- grep("^(bio[0-9]*_|delta_doy|doy|hopkins_|sin_doy|cos_doy)", all_cols, value = TRUE, ignore.case = TRUE)
  } else {
    clim_vars <- intersect(clim_vars, all_cols)
  }

  if (is.null(topo_vars)) {
    if (is.null(edaph_vars)) {
      topo_vars <- setdiff(
        grep("^(topo_|elevation$|elev_)", all_cols, value = TRUE, ignore.case = TRUE),
        clim_vars
      )
    } else {
      topo_vars <- setdiff(
        grep("^(topo_|elevation$|elev_)", all_cols, value = TRUE, ignore.case = TRUE),
        unique(c(clim_vars, edaph_vars))
      )
    }
  } else {
    topo_vars <- setdiff(intersect(topo_vars, all_cols), clim_vars)
  }

  if (is.null(edaph_vars)) {
    edaph_vars <- setdiff(all_cols, unique(c(clim_vars, topo_vars)))
  } else {
    edaph_vars <- setdiff(intersect(edaph_vars, all_cols), unique(c(clim_vars, topo_vars)))
  }

  scores_all <- compute_anchor_separation_scores(X_df, anchor_groups = anchor_groups)

  sel_clim <- if (length(clim_vars) > 0L) {
    screen_collinear_predictors(
      X_df[, clim_vars, drop = FALSE],
      cor.threshold = cor.threshold,
      vif.threshold = vif.threshold,
      anchor_groups = anchor_groups,
      var_scores = scores_all
    )$selected_variables
  } else {
    character(0)
  }

  sel_edaph <- if (length(edaph_vars) > 0L) {
    screen_collinear_predictors(
      X_df[, edaph_vars, drop = FALSE],
      cor.threshold = cor.threshold,
      vif.threshold = vif.threshold,
      anchor_groups = anchor_groups,
      var_scores = scores_all
    )$selected_variables
  } else {
    character(0)
  }

  sel_topo <- if (length(topo_vars) > 0L) {
    screen_collinear_predictors(
      X_df[, topo_vars, drop = FALSE],
      cor.threshold = cor.threshold,
      vif.threshold = vif.threshold,
      anchor_groups = anchor_groups,
      var_scores = scores_all
    )$selected_variables
  } else {
    character(0)
  }

  sel_joint <- unique(c(sel_edaph, sel_topo, sel_clim))
  pruned_joint <- setdiff(all_cols, sel_joint)
  message(sprintf(
    "Partitioned 3-Block Collinearity Pruning: retained %d Edaphic-Lithic (%s), %d Topographic (%s), and %d Climatic (%s) predictors.",
    length(sel_edaph), paste(sel_edaph, collapse = ", "),
    length(sel_topo), paste(sel_topo, collapse = ", "),
    length(sel_clim), paste(sel_clim, collapse = ", ")
  ))
  list(
    selected_variables = sel_joint,
    kept_vars = sel_joint,
    selected_edaphic = sel_edaph,
    selected_climatic = sel_clim,
    selected_topo = sel_topo,
    edaphic_variables = sel_edaph,
    climatic_variables = sel_clim,
    topo_variables = sel_topo,
    pruned_variables = pruned_joint,
    anchor_scores = scores_all
  )
}

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
  in_grp_mask <- !is.na(grp_all) & (grp_all %in% active_groups)
  if ("is_confident_label_noise" %in% names(d)) {
    cl_n <- as.logical(d$is_confident_label_noise); cl_n[is.na(cl_n)] <- FALSE
    is_t <- if ("is_type_specimen" %in% names(d)) as.logical(d$is_type_specimen) else rep(FALSE, nrow(d))
    is_t[is.na(is_t)] <- FALSE
    in_grp_mask <- in_grp_mask & (!cl_n | is_t)
  }
  if ("determiner_tier" %in% names(d)) {
    dt_v <- as.character(d$determiner_tier); dt_v[is.na(dt_v)] <- "Tier_3_Bronze"
    pr_v <- rep(FALSE, nrow(d))
    if ("is_promoted_bronze" %in% names(d)) { p1 <- as.logical(d$is_promoted_bronze); p1[is.na(p1)] <- FALSE; pr_v <- pr_v | p1 }
    is_t2 <- if ("is_type_specimen" %in% names(d)) as.logical(d$is_type_specimen) else rep(FALSE, nrow(d))
    is_t2[is.na(is_t2)] <- FALSE
    act_mask <- in_grp_mask & (grepl("Tier_1|Gold|Tier_2|Silver", dt_v, ignore.case = TRUE) | pr_v | is_t2)
    if (sum(act_mask) >= 6L) in_grp_mask <- act_mask
  }
  target_df <- d[in_grp_mask, , drop = FALSE]
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
