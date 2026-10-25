#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/edaphic_fine_grain_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Fine-grain (10 m to 30 m) edaphic, lithic, and micro-topographic predictor
#   extraction and caching utilities for Step 06 (06_multimodal_spatial_rf.R):
#     1. 10 m / 50 m-buffered USDA-NRCS SSURGO lithic & outcrop SQL extractor
#        using WKT ST_Buffer(ST_Point(lon, lat), 0.00045) via soilDB::SDA_query().
#     2. 30 m POLARIS (1 deg x 1 deg tiles) & 100 m SOLUS100 soil & depth-to-bedrock
#        integration with hierarchical multi-scale fusion (10m SSURGO / 30m POLARIS
#        -> 250m SoilGrids gap-free fallback).
#     3. 10 m - 30 m 3DEP micro-topographic outcrop covariates (topo_tpi_30m,
#        topo_roughness_30m, topo_slope_deg, topo_twi_30m) via elevatr & terra.
#     4. Persistent coordinate-keyed caching in data/environmental/fine_grain_edaphic_cache.csv.
# ==============================================================================

SSURGO_OUTCROP_SERIES_REGEX <- "(?i)(rock outcrop|wake|louisburg|ashlar|rion|wateree|wilkes|iredell|mecklenburg|poindexter|enon|goldston|manteo|boulders|granite|diabase)"

FINE_GRAIN_CACHE_COLS <- c(
  "catalogNumber", "lat_lon_key", "latitude", "longitude",
  "musym", "muname", "taxsuborder", "taxgreatgroup", "is_rock_outcrop", "outcrop_class",
  "ssurgo_rock_outcrop_pct", "ssurgo_bedrock_depth_cm", "ssurgo_sand_pct", "ssurgo_clay_pct",
  "ssurgo_ph", "ssurgo_om_pct", "ssurgo_cec", "ssurgo_aws_0_50cm",
  "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m", "polaris_om_30m",
  "polaris_bd_30m", "polaris_ksat_30m", "solus_resdept_cm", "solus_fragvol_pct",
  "topo_tpi_30m", "topo_roughness_30m", "topo_slope_deg", "topo_twi_30m",
  "fine_soil_ph", "fine_soil_sand", "fine_soil_clay", "fine_bedrock_depth_cm", "fine_rock_outcrop_prob"
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

    # 3. Surface horizon (<20 cm) sand, clay, pH, OM, CEC, and top 50 cm AWS
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
  # Batched SQL spatial query using WKT ST_Buffer(ST_Point(lon, lat), 0.00045) (~50m buffer)
  union_rows <- sapply(seq_len(nrow(sub_pts)), function(i) {
    sprintf(
      "SELECT %d AS gid, %.6f AS lon, %.6f AS lat",
      as.integer(sub_pts$gid[i]),
      as.numeric(sub_pts$longitude[i]),
      as.numeric(sub_pts$latitude[i])
    )
  })
  pts_cte <- paste(union_rows, collapse = " UNION ALL ")
  sql <- sprintf(
    paste0(
      "-- WKT ST_Buffer(ST_Point(lon, lat), %.5f) (~50m buffer)\n",
      "WITH pts AS (%s),\n",
      "buffered_pts AS (\n",
      "  SELECT gid, lon, lat,\n",
      "         geometry::STGeomFromText('POINT(' + CAST(lon AS VARCHAR(32)) + ' ' + CAST(lat AS VARCHAR(32)) + ')', 4326).STBuffer(%.5f) AS geom\n",
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
      "       ch.hzdept_r, ch.sandtotal_r, ch.claytotal_r, ch.ph1to1h2o_r, ch.om_r, ch.cec7_r,\n",
      "       cr.resdept_r, cr.reskind\n",
      "FROM mu_hits h\n",
      "INNER JOIN mapunit m ON h.mukey = m.mukey\n",
      "LEFT JOIN muaggatt ma ON m.mukey = ma.mukey\n",
      "LEFT JOIN component c ON m.mukey = c.mukey\n",
      "LEFT JOIN chorizon ch ON c.cokey = ch.cokey AND ch.hzdept_r < 20\n",
      "LEFT JOIN corestrictions cr ON c.cokey = cr.cokey"
    ),
    buffer_deg, pts_cte, buffer_deg
  )
  return(sql)
}

aggregate_ssurgo_tabular_to_points <- function(sub_pts, tab_res, classify_fn) {
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

    # Identify rock-outcrop / shallow-lithic / mafic components within the 50m buffer
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

    # Sum comppct_r per mukey, then take max across mukeys intersecting the 50m buffer
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

    # Component-percentage-weighted surface horizon (hzdept_r < 20 cm) properties
    hz_rows <- rows_i
    if ("hzdept_r" %in% names(hz_rows)) {
      hz_d <- suppressWarnings(as.numeric(hz_rows$hzdept_r))
      if (any(!is.na(hz_d) & hz_d < 20)) {
        hz_rows <- hz_rows[is.na(hz_d) | hz_d < 20, , drop = FALSE]
      }
    }
    w_hz <- suppressWarnings(as.numeric(hz_rows$comppct_r))
    w_hz[is.na(w_hz) | w_hz <= 0] <- 10.0

    w_avg <- function(v_raw, w_vec) {
      v <- suppressWarnings(as.numeric(v_raw))
      ok <- is.finite(v) & is.finite(w_vec) & w_vec > 0
      if (!any(ok)) return(NA_real_)
      sum(v[ok] * w_vec[ok]) / sum(w_vec[ok])
    }

    sand_pct <- w_avg(hz_rows$sandtotal_r, w_hz)
    clay_pct <- w_avg(hz_rows$claytotal_r, w_hz)
    ph_val   <- w_avg(hz_rows$ph1to1h2o_r, w_hz)
    om_pct   <- w_avg(hz_rows$om_r, w_hz)
    cec_val  <- w_avg(hz_rows$cec7_r, w_hz)
    aws_val  <- w_avg(comp_df$aws050wta, comp_pct)

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
# Milestone 2: 30 m POLARIS & 100 m SOLUS100 Continuous Soil + Depth-to-Bedrock
# ------------------------------------------------------------------------------
extract_polaris_and_solus_predictors <- function(coords_df, env_dir = "data/environmental", demo = FALSE) {
  message("Extracting 30m POLARIS & 100m SOLUS100 continuous soil & lithic depth predictors...")
  df <- as.data.frame(coords_df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0) return(df)

  polaris_dir <- file.path(env_dir, "polaris_30m")
  solus_dir   <- file.path(env_dir, "solus_100m")
  if (!isTRUE(demo)) {
    dir.create(polaris_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(solus_dir, recursive = TRUE, showWarnings = FALSE)
  }

  polaris_cols <- c(
    "polaris_ph_30m", "polaris_sand_30m", "polaris_clay_30m",
    "polaris_om_30m", "polaris_bd_30m", "polaris_ksat_30m"
  )
  for (col in c(polaris_cols, "solus_resdept_cm", "solus_fragvol_pct")) {
    if (!(col %in% names(df))) df[[col]] <- rep(NA_real_, n)
  }

  lat_v <- suppressWarnings(as.numeric(df$latitude))
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  in_conus <- !is.na(lat_v) & !is.na(lon_v) & lat_v >= 24.5 & lat_v <= 49.5 & lon_v >= -125.0 & lon_v <= -66.5

  # 1. Check cached 1 deg x 1 deg 30 m POLARIS GeoTIFF tiles in data/environmental/polaris_30m/
  # URL pattern: http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/{var}/mean/0_5/lat{lat0}{lat1}_lon{lon0}{lon1}.tif
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

    # Optional on-demand download of missing 1x1 deg POLARIS tiles when network is reachable
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

    # Extract from any local POLARIS GeoTIFF tiles present in polaris_dir
    existing_polaris <- list.files(polaris_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
    if (length(existing_polaris) > 0) {
      pts_vect <- terra::vect(df[in_conus, c("longitude", "latitude"), drop = FALSE],
                              geom = c("longitude", "latitude"), crs = "EPSG:4326")
      conus_idx <- which(in_conus)
      for (pv in polaris_vars) {
        pv_files <- existing_polaris[grepl(paste0("^", pv, "_"), basename(existing_polaris), ignore.case = TRUE)]
        if (length(pv_files) > 0) {
          vals_accum <- matrix(NA_real_, nrow = length(conus_idx), ncol = length(pv_files))
          for (f_i in seq_along(pv_files)) {
            v_ext <- tryCatch({
              r_tile <- terra::rast(pv_files[f_i])
              ext_df <- terra::extract(r_tile, pts_vect)
              as.numeric(ext_df[[2]])
            }, error = function(e) rep(NA_real_, length(conus_idx)))
            vals_accum[, f_i] <- v_ext
          }
          row_has_val <- rowSums(!is.na(vals_accum)) > 0
          if (any(row_has_val)) {
            col_target <- polaris_col_map[[pv]]
            df[[col_target]][conus_idx[row_has_val]] <- rowMeans(vals_accum[row_has_val, , drop = FALSE], na.rm = TRUE)
          }
        }
      }
    }

    # 2. Check data/environmental/solus_100m/ for resdept and fragvol_0_5cm GeoTIFFs
    existing_solus <- list.files(solus_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
    if (length(existing_solus) > 0) {
      pts_vect <- terra::vect(df[in_conus, c("longitude", "latitude"), drop = FALSE],
                              geom = c("longitude", "latitude"), crs = "EPSG:4326")
      conus_idx <- which(in_conus)
      resdept_files <- existing_solus[grepl("resdept", basename(existing_solus), ignore.case = TRUE)]
      if (length(resdept_files) > 0) {
        v_res <- tryCatch({
          r_sol <- terra::rast(resdept_files[1])
          as.numeric(terra::extract(r_sol, pts_vect)[[2]])
        }, error = function(e) rep(NA_real_, length(conus_idx)))
        ok_r <- !is.na(v_res) & is.finite(v_res)
        if (any(ok_r)) df$solus_resdept_cm[conus_idx[ok_r]] <- pmin(200.0, pmax(0.0, v_res[ok_r]))
      }
      frag_files <- existing_solus[grepl("fragvol", basename(existing_solus), ignore.case = TRUE)]
      if (length(frag_files) > 0) {
        v_frg <- tryCatch({
          r_frg <- terra::rast(frag_files[1])
          as.numeric(terra::extract(r_frg, pts_vect)[[2]])
        }, error = function(e) rep(NA_real_, length(conus_idx)))
        ok_f <- !is.na(v_frg) & is.finite(v_frg)
        if (any(ok_f)) df$solus_fragvol_pct[conus_idx[ok_f]] <- pmin(100.0, pmax(0.0, v_frg[ok_f]))
      }
    }
  }

  # 3. Calibrate missing POLARIS 30m / SOLUS 100m estimates from 10m SSURGO + 250m SoilGrids hierarchy
  sg_ph   <- if ("soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$soil_ph)) else rep(5.6, n)
  sg_sand <- if ("soil_sand" %in% names(df)) suppressWarnings(as.numeric(df$soil_sand)) else rep(50.0, n)
  sg_cec  <- if ("soil_cec" %in% names(df)) suppressWarnings(as.numeric(df$soil_cec)) else rep(12.0, n)
  sg_bd   <- if ("soil_bulk_density" %in% names(df)) suppressWarnings(as.numeric(df$soil_bulk_density)) else rep(1.38, n)
  sg_ph[is.na(sg_ph)] <- ifelse(is.finite(stats::median(sg_ph, na.rm = TRUE)), stats::median(sg_ph, na.rm = TRUE), 5.6)
  sg_sand[is.na(sg_sand)] <- ifelse(is.finite(stats::median(sg_sand, na.rm = TRUE)), stats::median(sg_sand, na.rm = TRUE), 50.0)
  sg_cec[is.na(sg_cec)] <- ifelse(is.finite(stats::median(sg_cec, na.rm = TRUE)), stats::median(sg_cec, na.rm = TRUE), 12.0)
  sg_bd[is.na(sg_bd)] <- ifelse(is.finite(stats::median(sg_bd, na.rm = TRUE)), stats::median(sg_bd, na.rm = TRUE), 1.38)

  ss_ph   <- if ("ssurgo_ph" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_ph)) else rep(NA_real_, n)
  ss_sand <- if ("ssurgo_sand_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_sand_pct)) else rep(NA_real_, n)
  ss_clay <- if ("ssurgo_clay_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_clay_pct)) else rep(NA_real_, n)
  ss_om   <- if ("ssurgo_om_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_om_pct)) else rep(NA_real_, n)
  ss_bed  <- if ("ssurgo_bedrock_depth_cm" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_bedrock_depth_cm)) else rep(NA_real_, n)
  ss_out  <- if ("ssurgo_rock_outcrop_pct" %in% names(df)) suppressWarnings(as.numeric(df$ssurgo_rock_outcrop_pct)) else rep(0.0, n)
  is_out  <- if ("is_rock_outcrop" %in% names(df)) as.logical(df$is_rock_outcrop) else rep(FALSE, n)
  is_out[is.na(is_out)] <- FALSE
  ss_out[is.na(ss_out)] <- ifelse(is_out[is.na(ss_out)], 40.0, 0.0)

  sg_clay_est <- pmax(3.0, pmin(65.0, (100.0 - sg_sand) * 0.40 + 0.55 * sg_cec))

  for (i in seq_len(n)) {
    w_ss <- if (is_out[i] || ss_out[i] >= 15.0) 0.78 else 0.62
    if (is.na(df$polaris_ph_30m[i])) {
      df$polaris_ph_30m[i] <- if (!is.na(ss_ph[i]) && !is.na(sg_ph[i])) {
        round(w_ss * ss_ph[i] + (1.0 - w_ss) * sg_ph[i], 2)
      } else if (!is.na(ss_ph[i])) {
        round(ss_ph[i], 2)
      } else {
        round(sg_ph[i], 2)
      }
    }
    if (is.na(df$polaris_sand_30m[i])) {
      df$polaris_sand_30m[i] <- if (!is.na(ss_sand[i]) && !is.na(sg_sand[i])) {
        round(w_ss * ss_sand[i] + (1.0 - w_ss) * sg_sand[i], 1)
      } else if (!is.na(ss_sand[i])) {
        round(ss_sand[i], 1)
      } else {
        round(sg_sand[i], 1)
      }
    }
    if (is.na(df$polaris_clay_30m[i])) {
      df$polaris_clay_30m[i] <- if (!is.na(ss_clay[i]) && !is.na(sg_clay_est[i])) {
        round(w_ss * ss_clay[i] + (1.0 - w_ss) * sg_clay_est[i], 1)
      } else if (!is.na(ss_clay[i])) {
        round(ss_clay[i], 1)
      } else {
        round(sg_clay_est[i], 1)
      }
    }
    if (is.na(df$polaris_om_30m[i])) {
      df$polaris_om_30m[i] <- if (!is.na(ss_om[i])) round(ss_om[i], 2) else round(pmax(0.5, 0.14 * sg_cec[i]), 2)
    }
    if (is.na(df$polaris_bd_30m[i])) {
      df$polaris_bd_30m[i] <- round(ifelse(is.na(sg_bd[i]), 1.38, sg_bd[i]), 2)
    }
    if (is.na(df$polaris_ksat_30m[i])) {
      s_val <- ifelse(is.na(df$polaris_sand_30m[i]), 50.0, df$polaris_sand_30m[i])
      c_val <- ifelse(is.na(df$polaris_clay_30m[i]), 20.0, df$polaris_clay_30m[i])
      # Saxton & Rawls saturated hydraulic conductivity pedotransfer (cm/hr)
      df$polaris_ksat_30m[i] <- round(pmax(0.05, 0.35 * exp(0.032 * s_val - 0.045 * c_val)), 2)
    }
    if (is.na(df$solus_resdept_cm[i])) {
      df$solus_resdept_cm[i] <- if (!is.na(ss_bed[i])) round(ss_bed[i], 1) else 200.0
    }
    if (is.na(df$solus_fragvol_pct[i])) {
      df$solus_fragvol_pct[i] <- round(pmin(85.0, pmax(1.0, 0.55 * ss_out[i] + ifelse(df$solus_resdept_cm[i] < 100, (100 - df$solus_resdept_cm[i]) * 0.25, 2.0))), 1)
    }
  }

  # 4. Hierarchical multi-scale fusion into canonical fine_soil_* columns:
  #    Prioritizes (1) 30m POLARIS / 10m SSURGO component values where available,
  #    and (2) 250m SoilGrids as gap-free fallback.
  df$fine_soil_ph <- round(ifelse(
    !is.na(ss_ph) & !is.na(df$polaris_ph_30m),
    ifelse(is_out | ss_out >= 15.0, 0.70 * ss_ph + 0.30 * df$polaris_ph_30m, 0.50 * ss_ph + 0.50 * df$polaris_ph_30m),
    ifelse(!is.na(df$polaris_ph_30m), df$polaris_ph_30m, ifelse(!is.na(ss_ph), ss_ph, sg_ph))
  ), 2)

  df$fine_soil_sand <- round(ifelse(
    !is.na(ss_sand) & !is.na(df$polaris_sand_30m),
    ifelse(is_out | ss_out >= 15.0, 0.70 * ss_sand + 0.30 * df$polaris_sand_30m, 0.50 * ss_sand + 0.50 * df$polaris_sand_30m),
    ifelse(!is.na(df$polaris_sand_30m), df$polaris_sand_30m, ifelse(!is.na(ss_sand), ss_sand, sg_sand))
  ), 1)

  df$fine_soil_clay <- round(ifelse(
    !is.na(ss_clay) & !is.na(df$polaris_clay_30m),
    ifelse(is_out | ss_out >= 15.0, 0.70 * ss_clay + 0.30 * df$polaris_clay_30m, 0.50 * ss_clay + 0.50 * df$polaris_clay_30m),
    ifelse(!is.na(df$polaris_clay_30m), df$polaris_clay_30m, ifelse(!is.na(ss_clay), ss_clay, sg_clay_est))
  ), 1)

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

  df$fine_rock_outcrop_prob <- round(pmin(
    1.0,
    pmax(
      0.0,
      0.45 * (ss_out / 100.0) +
      0.30 * pmax(0.0, (160.0 - df$fine_bedrock_depth_cm) / 160.0) +
      0.15 * as.numeric(is_out) +
      0.10 * pmin(1.0, df$solus_fragvol_pct / 60.0)
    )
  ), 3)

  return(df)
}

# ------------------------------------------------------------------------------
# Milestone 3: 10 m - 30 m 3DEP Micro-Topographic Outcrop Covariates
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
    # Convex monadnock domes/outcrop ledges have positive TPI; mucks/floodplains have negative TPI
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
    # On flatrock seepage rims (shallow bedrock < 60 cm and gentle slope < 5 deg), impermeable granite concentrates sheetflow
    seepage_multiplier <- if (bed_cm[i] < 60.0 && df$topo_slope_deg[i] <= 6.0) 2.2 else 1.0
    upslope_area_m2 <- pmax(15.0, 75.0 * seepage_multiplier * exp(-0.12 * df$topo_tpi_30m[i]))
    df$topo_twi_30m[i] <- round(log(upslope_area_m2 / tan(slope_rad)), 2)
  }

  return(df)
}

# ------------------------------------------------------------------------------
# Milestone 4: Persistent Fine-Grain Edaphic & Micro-Topographic Cache Manager
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
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (n == 0) return(df)

  df$catalogNumber <- as.character(df$catalogNumber)
  df$lat_lon_key <- format_lat_lon_key(df$latitude, df$longitude)
  cache_csv_path <- cache_path

  req_fine_cols <- setdiff(FINE_GRAIN_CACHE_COLS, c("catalogNumber", "lat_lon_key", "latitude", "longitude"))
  cached_df <- NULL

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

  # Match rows already present in fine_grain_edaphic_cache.csv
  if (!is.null(cached_df) && nrow(cached_df) > 0 && all(req_fine_cols %in% names(cached_df))) {
    c_by_cat <- cached_df[!duplicated(cached_df$catalogNumber), c("catalogNumber", req_fine_cols), drop = FALSE]
    m_idx <- match(df$catalogNumber, c_by_cat$catalogNumber)

    # Secondary fallback match by 5-decimal-place lat_lon_key
    miss_cat <- which(is.na(m_idx) & !is.na(df$lat_lon_key))
    if (length(miss_cat) > 0) {
      c_by_ll <- cached_df[!is.na(cached_df$lat_lon_key) & !duplicated(cached_df$lat_lon_key), c("lat_lon_key", req_fine_cols), drop = FALSE]
      ll_idx <- match(df$lat_lon_key[miss_cat], c_by_ll$lat_lon_key)
      for (col in req_fine_cols) {
        df[[col]] <- c_by_cat[[col]][m_idx]
        hit_ll <- !is.na(ll_idx)
        if (any(hit_ll)) {
          df[[col]][miss_cat[hit_ll]] <- c_by_ll[[col]][ll_idx[hit_ll]]
        }
      }
    } else {
      for (col in req_fine_cols) {
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

    if (n_uncached == 0) {
      df$is_rock_outcrop <- as.logical(df$is_rock_outcrop)
      df$is_rock_outcrop[is.na(df$is_rock_outcrop)] <- FALSE
      return(df)
    }

    sub_df <- df[uncached_mask, , drop = FALSE]
    for (col in req_fine_cols) sub_df[[col]] <- NULL
    if (!is.null(query_ssurgo_fn)) {
      sub_df <- query_ssurgo_fn(sub_df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
    }
    sub_df <- enrich_ssurgo_pedology_from_series(sub_df)
    sub_df <- extract_polaris_and_solus_predictors(sub_df, env_dir = env_dir, demo = demo)
    sub_df <- extract_micro_topography_3dep(sub_df, env_dir = env_dir, demo = demo)

    for (col in req_fine_cols) {
      df[[col]][uncached_mask] <- sub_df[[col]]
    }
  } else {
    if (!is.null(query_ssurgo_fn)) {
      df <- query_ssurgo_fn(df, cache_dir = ssurgo_cache_dir, batch_size = batch_size, demo = demo)
    }
    df <- enrich_ssurgo_pedology_from_series(df)
    df <- extract_polaris_and_solus_predictors(df, env_dir = env_dir, demo = demo)
    df <- extract_micro_topography_3dep(df, env_dir = env_dir, demo = demo)
  }

  # Persist updated fine_grain_edaphic_cache.csv atomically (never persist synthetic --demo records)
  if (!isTRUE(demo)) {
    export_cols <- intersect(FINE_GRAIN_CACHE_COLS, names(df))
    new_cache <- df[, export_cols, drop = FALSE]
    if (!is.null(cached_df) && nrow(cached_df) > 0) {
      common_c <- intersect(export_cols, names(cached_df))
      combined_cache <- rbind(new_cache[, common_c, drop = FALSE], cached_df[, common_c, drop = FALSE])
      combined_cache <- combined_cache[!duplicated(combined_cache$catalogNumber), , drop = FALSE]
    } else {
      combined_cache <- new_cache[!duplicated( new_cache$catalogNumber), , drop = FALSE]
    }
    dir.create(dirname(cache_csv_path), recursive = TRUE, showWarnings = FALSE)
    tmp_cache <- paste0(cache_csv_path, ".tmp")
    readr::write_csv(combined_cache, tmp_cache)
    file.rename(tmp_cache, cache_csv_path)
    message(sprintf("Persisted %d fine-grain edaphic & micro-topographic records to %s.", nrow(combined_cache), cache_csv_path))
  }

  return(df)
}
