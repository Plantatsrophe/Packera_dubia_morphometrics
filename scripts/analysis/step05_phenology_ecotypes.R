# ==============================================================================
# Step 05 Module: Phenology, Circular DOY, & Type-Locality Ecotype Controls
# Extracted from 05_multimodal_spatial_rf.R
# ==============================================================================

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

read_config_phenology_bounds <- function(config_path = "config/config.yaml") {
  min_doy <- 60.0
  max_doy <- 220.0
  if (!is.null(config_path) && nzchar(config_path) && file.exists(config_path)) {
    lines <- tryCatch(readLines(config_path, warn = FALSE), error = function(e) character(0))
    if (length(lines) > 0L) {
      m_min <- grep("^\\s*min_flowering_doy\\s*:", lines, value = TRUE)
      if (length(m_min) > 0L) {
        val <- suppressWarnings(as.numeric(sub("#.*$", "", sub("^[^:]+:\\s*", "", m_min[1]))))
        if (is.finite(val) && val >= 1.0 && val <= 366.0) min_doy <- val
      }
      m_max <- grep("^\\s*max_flowering_doy\\s*:", lines, value = TRUE)
      if (length(m_max) > 0L) {
        val <- suppressWarnings(as.numeric(sub("#.*$", "", sub("^[^:]+:\\s*", "", m_max[1]))))
        if (is.finite(val) && val >= 1.0 && val <= 366.0) max_doy <- val
      }
    }
  }
  list(min_flowering_doy = min_doy, max_flowering_doy = max_doy)
}

compute_circular_mean_doy <- function(doy_vec) {
  v <- suppressWarnings(as.numeric(doy_vec))
  v <- v[is.finite(v) & v >= 1.0 & v <= 366.0]
  if (length(v) == 0L) return(NA_real_)
  theta <- 2.0 * pi * (v / 365.25)
  sin_m <- mean(sin(theta))
  cos_m <- mean(cos(theta))
  if (!is.finite(sin_m) || !is.finite(cos_m) || sqrt(sin_m^2 + cos_m^2) < 1e-8) {
    return(as.numeric(stats::median(v)))
  }
  ang <- atan2(sin_m, cos_m)
  if (ang <= 0.0) ang <- ang + 2.0 * pi
  c_doy <- ang * 365.25 / (2.0 * pi)
  if (c_doy <= 0.0) c_doy <- c_doy + 365.25
  as.numeric(c_doy)
}

circular_signed_doy_offset <- function(doy_vec, center_doy) {
  v <- suppressWarnings(as.numeric(doy_vec))
  if (is.null(center_doy) || length(center_doy) == 0L || !is.finite(center_doy[1])) {
    return(rep(NA_real_, length(v)))
  }
  d_theta <- 2.0 * pi * ((v - as.numeric(center_doy[1])) / 365.25)
  atan2(sin(d_theta), cos(d_theta)) * 365.25 / (2.0 * pi)
}

wrap_doy_365 <- function(doy_vec) {
  v <- suppressWarnings(as.numeric(doy_vec))
  w <- v %% 365.25
  idx_le <- which(is.finite(w) & w <= 0.0)
  if (length(idx_le) > 0L) w[idx_le] <- w[idx_le] + 365.25
  w
}

model_latitudinal_spring_baseline <- function(flowering_df, circ_mean_doy = NULL, lat_bounds = NULL) {
  message("Modeling Cohort-Adaptive Latitudinal Anthesis Baseline (Hopkins' Bioclimatic Law, Hopkins 1918) on verified anthesis cohort...")

  lat_col <- if ("decimalLatitude" %in% names(flowering_df) && any(!is.na(flowering_df$decimalLatitude))) {
    suppressWarnings(as.numeric(flowering_df$decimalLatitude))
  } else if ("latitude" %in% names(flowering_df)) {
    suppressWarnings(as.numeric(flowering_df$latitude))
  } else {
    rep(NA_real_, nrow(flowering_df))
  }
  flowering_df$decimalLatitude <- lat_col
  doy_num <- suppressWarnings(as.numeric(flowering_df$doy))
  flowering_df$doy <- doy_num

  keep_valid <- !is.na(flowering_df$decimalLatitude) & is.finite(flowering_df$decimalLatitude) &
                flowering_df$decimalLatitude >= -90.0 & flowering_df$decimalLatitude <= 90.0 &
                !is.na(flowering_df$doy) & is.finite(flowering_df$doy) &
                flowering_df$doy >= 1.0 & flowering_df$doy <= 366.0
  if (!is.null(lat_bounds) && length(lat_bounds) == 2L && all(is.finite(lat_bounds))) {
    keep_valid <- keep_valid &
                  flowering_df$decimalLatitude >= min(lat_bounds) &
                  flowering_df$decimalLatitude <= max(lat_bounds)
  }
  valid <- flowering_df[keep_valid, , drop = FALSE]

  ref_circ_doy <- if (!is.null(circ_mean_doy) && length(circ_mean_doy) >= 1L && is.finite(circ_mean_doy[1])) {
    as.numeric(circ_mean_doy[1])
  } else {
    compute_circular_mean_doy(valid$doy)
  }

  if (nrow(valid) == 0L) {
    message("[WARN] Zero georeferenced anthesis records with valid DOY available for latitudinal baseline.")
    return(list(
      intercept = ref_circ_doy, slope = 0.0,
      b_elev = 0.0, ref_elev = 0.0, b_mat = 0.0, ref_mat = 0.0,
      min_lat = NA_real_, pred_min_lat = ref_circ_doy,
      circ_mean_doy = ref_circ_doy,
      hopkins_expected = numeric(0), r_squared = 0.0, p_value = NA_real_,
      n = 0L, model = NULL
    ))
  }

  # Unwrap circular DOY continuously around cohort circular mean DOY so winter/wrap-around
  # and late-summer alpine cohorts regress smoothly without calendar boundary discontinuities
  valid$doy_offset <- circular_signed_doy_offset(valid$doy, ref_circ_doy)
  valid$doy_unwrapped <- ref_circ_doy + valid$doy_offset

  if (nrow(valid) >= 10L) {
    p1 <- as.numeric(stats::quantile(valid$doy_offset, 0.01, na.rm = TRUE))
    p99 <- as.numeric(stats::quantile(valid$doy_offset, 0.99, na.rm = TRUE))
    valid_trimmed <- valid[valid$doy_offset >= p1 & valid$doy_offset <= p99, , drop = FALSE]
    if (nrow(valid_trimmed) < 3L) valid_trimmed <- valid
  } else {
    valid_trimmed <- valid
  }

  min_lat <- min(valid_trimmed$decimalLatitude, na.rm = TRUE)
  lat_sd <- if (nrow(valid_trimmed) >= 2L) stats::sd(valid_trimmed$decimalLatitude, na.rm = TRUE) else 0.0
  if (nrow(valid_trimmed) < 3L || !is.finite(lat_sd) || lat_sd <= 1e-6) {
    b0 <- ref_circ_doy
    b_lat <- 0.0
    pred_min_lat <- ref_circ_doy
    hopkins_expected <- pred_min_lat + 4.0 * (valid_trimmed$decimalLatitude - min_lat)
    return(list(
      intercept = b0, slope = b_lat,
      b_elev = 0.0, ref_elev = 0.0, b_mat = 0.0, ref_mat = 0.0,
      min_lat = min_lat, pred_min_lat = pred_min_lat,
      circ_mean_doy = ref_circ_doy,
      hopkins_expected = hopkins_expected, r_squared = 0.0, p_value = NA_real_,
      n = nrow(valid_trimmed), model = NULL
    ))
  }

  fit <- if (requireNamespace("MASS", quietly = TRUE)) {
    tryCatch(
      MASS::rlm(doy_unwrapped ~ decimalLatitude, data = valid_trimmed),
      error = function(e) stats::lm(doy_unwrapped ~ decimalLatitude, data = valid_trimmed)
    )
  } else {
    stats::lm(doy_unwrapped ~ decimalLatitude, data = valid_trimmed)
  }
  coefs <- stats::coef(fit)
  b0 <- as.numeric(coefs[1])
  b_lat <- if (length(coefs) >= 2L && is.finite(coefs[2])) as.numeric(coefs[2]) else 0.0
  pred_min_lat <- b0 + b_lat * min_lat
  decimalLatitude <- valid_trimmed$decimalLatitude
  hopkins_expected <- pred_min_lat + 4.0 * (decimalLatitude - min_lat)
  valid_trimmed$hopkins_expected <- hopkins_expected

  # Multivariate bioclimatic adjustment for latitude, elevation, and MAT (Hopkins 1918) when covariates are available
  elev_col <- if ("topo_elev_m" %in% names(valid_trimmed)) {
    suppressWarnings(as.numeric(valid_trimmed$topo_elev_m))
  } else if ("elevation" %in% names(valid_trimmed)) {
    suppressWarnings(as.numeric(valid_trimmed$elevation))
  } else {
    rep(NA_real_, nrow(valid_trimmed))
  }
  mat_col <- if ("bio1_temp_mean" %in% names(valid_trimmed)) {
    suppressWarnings(as.numeric(valid_trimmed$bio1_temp_mean))
  } else {
    rep(NA_real_, nrow(valid_trimmed))
  }
  b_elev <- 0.0; ref_elev <- 0.0; b_mat <- 0.0; ref_mat <- 0.0
  if (sum(is.finite(elev_col)) >= 10L && stats::sd(elev_col, na.rm = TRUE) > 5.0) {
    ref_elev <- stats::median(elev_col, na.rm = TRUE)
    b_elev <- 0.0125 # Hopkins' Bioclimatic Law: +1.25 days per +100 m elevation (+0.0125 d/m)
  }
  if (sum(is.finite(mat_col)) >= 10L && stats::sd(mat_col, na.rm = TRUE) > 0.5) {
    ref_mat <- stats::median(mat_col, na.rm = TRUE)
    fit_mat <- tryCatch(
      stats::lm(I(doy_unwrapped - (b0 + b_lat * decimalLatitude)) ~ mat_col, data = valid_trimmed),
      error = function(e) NULL
    )
    if (!is.null(fit_mat) && length(stats::coef(fit_mat)) >= 2L && is.finite(stats::coef(fit_mat)[2])) {
      b_mat <- max(-2.5, min(2.5, as.numeric(stats::coef(fit_mat)[2])))
    }
  }

  fit_ols <- stats::lm(doy_unwrapped ~ decimalLatitude, data = valid_trimmed)
  s_ols <- summary(fit_ols)
  r2 <- if (!is.null(s_ols$r.squared) && is.finite(s_ols$r.squared)) s_ols$r.squared else 0.0
  p_val <- if (nrow(s_ols$coefficients) >= 2L && ncol(s_ols$coefficients) >= 4L) s_ols$coefficients[2, 4] else NA_real_

  message(sprintf("Empirical Baseline (Hopkins 1918): DOY = %.3f + %.3f * decimalLatitude (circ_mean = %.1f, b_elev = %.4f d/m, b_mat = %.3f d/C, R2 = %.4f, p = %.4e, N = %d)",
                  b0, b_lat, ref_circ_doy, b_elev, b_mat, r2, p_val, nrow(valid_trimmed)))
  return(list(
    intercept = b0, slope = b_lat,
    b_elev = b_elev, ref_elev = ref_elev, b_mat = b_mat, ref_mat = ref_mat,
    min_lat = min_lat, pred_min_lat = pred_min_lat,
    circ_mean_doy = ref_circ_doy,
    hopkins_expected = hopkins_expected, r_squared = r2, p_value = p_val,
    n = nrow(valid_trimmed), model = fit_ols
  ))
}

compute_phenological_anomalies <- function(df, output_csv = NULL, summary_csv = NULL, demo = FALSE, config_path = "config/config.yaml") {
  message("Computing cohort-adaptive circular phenological anomalies (Delta DOY)...")
  df$doy <- extract_doy_robust(df)
  lat_val <- if ("decimalLatitude" %in% names(df) && any(!is.na(df$decimalLatitude))) {
    suppressWarnings(as.numeric(df$decimalLatitude))
  } else if ("latitude" %in% names(df)) {
    suppressWarnings(as.numeric(df$latitude))
  } else {
    rep(NA_real_, nrow(df))
  }
  df$decimalLatitude <- lat_val
  df$latitude <- lat_val

  valid_doy <- !is.na(df$doy) & is.finite(df$doy) & df$doy >= 1.0 & df$doy <= 366.0
  valid_lat <- !is.na(df$decimalLatitude) & is.finite(df$decimalLatitude) &
               df$decimalLatitude >= -90.0 & df$decimalLatitude <= 90.0

  # ----------------------------------------------------------------------------
  # 3-Tier Cohort-Adaptive Anthesis Verification Filter
  # ----------------------------------------------------------------------------
  # Darwin Core reproductiveCondition check: exclude fruiting-only / sterile-only sheets
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

  # Computer Vision check: voucher_phenological_state == "anthesis"
  if ("voucher_phenological_state" %in% names(df)) {
    v_state <- tolower(trimws(as.character(df$voucher_phenological_state)))
    v_state[is.na(v_state)] <- ""
  } else {
    v_state <- rep("", nrow(df))
  }

  if ("total_capitula" %in% names(df)) {
    tot_caps <- suppressWarnings(as.numeric(df$total_capitula))
    tot_caps[is.na(tot_caps)] <- 0
  } else if ("capitula_count" %in% names(df)) {
    tot_caps <- suppressWarnings(as.numeric(df$capitula_count))
    tot_caps[is.na(tot_caps)] <- 0
  } else {
    tot_caps <- rep(0, nrow(df))
  }

  n_anth_vec <- if ("n_anthesis" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis)) else rep(0, nrow(df))
  n_anth_vec[is.na(n_anth_vec)] <- 0
  has_valid_cv <- ((v_state != "" & v_state != "unknown" & v_state != "na") & tot_caps > 0) | (n_anth_vec > 0)
  is_cv_anth_label <- grepl("anthesis", v_state) | (n_anth_vec > 0)
  raw_anthesis_cand <- valid_doy & valid_lat & pass_dwc &
                       ((has_valid_cv & is_cv_anth_label) | (!has_valid_cv))

  # Empirical cohort latitude bounds: [min(lat) - 3.0, max(lat) + 3.0]
  lat_pool <- if (any(raw_anthesis_cand)) df$decimalLatitude[raw_anthesis_cand] else df$decimalLatitude[valid_lat]
  if (length(lat_pool) > 0L) {
    lat_min_gate <- max(-90.0, min(lat_pool, na.rm = TRUE) - 3.0)
    lat_max_gate <- min(90.0, max(lat_pool, na.rm = TRUE) + 3.0)
  } else {
    lat_min_gate <- -90.0
    lat_max_gate <- 90.0
  }

  # Cohort-adaptive circular DOY gate (99% angular coverage arc +/- 21 days);
  # fall back to config/config.yaml phenology bounds ONLY when N_verified_anthesis < 5
  n_raw_anthesis <- sum(raw_anthesis_cand)
  if (n_raw_anthesis >= 5L) {
    circ_mean_doy <- compute_circular_mean_doy(df$doy[raw_anthesis_cand])
    cand_offsets <- circular_signed_doy_offset(df$doy[raw_anthesis_cand], circ_mean_doy)
    arc_lo <- max(-182.625, as.numeric(stats::quantile(cand_offsets, 0.005, na.rm = TRUE)) - 21.0)
    arc_hi <- min(182.625, as.numeric(stats::quantile(cand_offsets, 0.995, na.rm = TRUE)) + 21.0)
    all_offsets <- circular_signed_doy_offset(df$doy, circ_mean_doy)
    pass_doy_gate <- valid_doy & is.finite(all_offsets) & (all_offsets >= arc_lo) & (all_offsets <= arc_hi)
    message(sprintf("[INFO] Cohort-adaptive circular anthesis gate: circ_mean_doy=%.1f, arc=[%+.1f, %+.1f] d, lat=[%.2f, %.2f] (N_cand=%d).",
                    circ_mean_doy, arc_lo, arc_hi, lat_min_gate, lat_max_gate, n_raw_anthesis))
  } else {
    cfg_bounds <- read_config_phenology_bounds(config_path)
    min_d <- cfg_bounds$min_flowering_doy
    max_d <- cfg_bounds$max_flowering_doy
    if (min_d <= max_d) {
      pass_doy_gate <- valid_doy & (df$doy >= min_d) & (df$doy <= max_d)
    } else {
      pass_doy_gate <- valid_doy & (df$doy >= min_d | df$doy <= max_d)
    }
    circ_mean_doy <- if (any(raw_anthesis_cand)) {
      compute_circular_mean_doy(df$doy[raw_anthesis_cand])
    } else if (any(valid_doy)) {
      compute_circular_mean_doy(df$doy[valid_doy])
    } else {
      compute_circular_mean_doy(c(min_d, max_d))
    }
    message(sprintf("[WARN] Verified anthesis candidates N=%d (< 5); falling back to config phenology bounds DOY [%.0f, %.0f], lat=[%.2f, %.2f].",
                    n_raw_anthesis, min_d, max_d, lat_min_gate, lat_max_gate))
  }

  pass_calendar <- pass_doy_gate & valid_lat &
                   (df$decimalLatitude >= lat_min_gate) & (df$decimalLatitude <= lat_max_gate)

  cv_anthesis <- rep(FALSE, nrow(df))
  cv_fruit <- rep(FALSE, nrow(df))
  cv_sterile <- rep(FALSE, nrow(df))

  cv_anthesis[has_valid_cv & is_cv_anth_label] <- TRUE
  cv_fruit[has_valid_cv & !is_cv_anth_label & grepl("fruit", v_state)] <- TRUE
  cv_sterile[has_valid_cv & !is_cv_anth_label & (v_state %in% c("sterile", "bud", "senescent") | grepl("sterile|bud", v_state))] <- TRUE

  missing_cv <- !has_valid_cv
  cv_anthesis[missing_cv & pass_calendar & pass_dwc] <- TRUE
  cv_fruit[missing_cv & dwc_fruit_only] <- TRUE
  cv_sterile[missing_cv & dwc_sterile_only] <- TRUE

  is_verified_anthesis <- pass_calendar & pass_dwc & cv_anthesis
  is_fruit_excluded <- pass_calendar & (cv_fruit | dwc_fruit_only) & !is_verified_anthesis
  is_sterile_excluded <- pass_calendar & (cv_sterile | dwc_sterile_only) & !is_verified_anthesis & !is_fruit_excluded

  df$is_verified_anthesis <- is_verified_anthesis
  df$anthesis_filter_status <- "Excluded_Calendar_Gate"
  df$anthesis_filter_status[pass_calendar] <- "Excluded_Other"
  df$anthesis_filter_status[is_sterile_excluded] <- "Excluded_Sterile"
  df$anthesis_filter_status[is_fruit_excluded] <- "Excluded_Fruit_Only"
  df$anthesis_filter_status[is_verified_anthesis] <- "Verified_Anthesis"

  n_verified <- sum(is_verified_anthesis)
  n_fruit <- sum(is_fruit_excluded)
  n_sterile <- sum(is_sterile_excluded)
  message(sprintf("[INFO] Filtered phenology vouchers: %d verified anthesis, %d fruiting-only (excluded), %d sterile (excluded).",
                  n_verified, n_fruit, n_sterile))

  anthesis_cohort <- df[df$is_verified_anthesis, , drop = FALSE]
  if (nrow(anthesis_cohort) < 5L) {
    message("[WARN] Verified anthesis cohort has < 5 vouchers; using calendar-gated cohort as baseline fallback.")
    anthesis_cohort <- df[pass_calendar, , drop = FALSE]
    if (nrow(anthesis_cohort) == 0L) anthesis_cohort <- df[valid_doy & valid_lat, , drop = FALSE]
  }
  baseline <- model_latitudinal_spring_baseline(
    anthesis_cohort,
    circ_mean_doy = circ_mean_doy,
    lat_bounds = c(lat_min_gate, lat_max_gate)
  )

  df$expected_doy <- NA_real_
  df$hopkins_expected <- NA_real_
  has_lat <- valid_lat
  elev_v <- if ("topo_elev_m" %in% names(df)) {
    suppressWarnings(as.numeric(df$topo_elev_m))
  } else if ("elevation" %in% names(df)) {
    suppressWarnings(as.numeric(df$elevation))
  } else {
    rep(NA_real_, nrow(df))
  }
  mat_v <- if ("bio1_temp_mean" %in% names(df)) {
    suppressWarnings(as.numeric(df$bio1_temp_mean))
  } else {
    rep(NA_real_, nrow(df))
  }
  elev_adj <- ifelse(is.finite(elev_v) & !is.null(baseline$b_elev) & baseline$b_elev != 0,
                     baseline$b_elev * (elev_v - baseline$ref_elev), 0.0)
  mat_adj <- ifelse(is.finite(mat_v) & !is.null(baseline$b_mat) & baseline$b_mat != 0,
                    baseline$b_mat * (mat_v - baseline$ref_mat), 0.0)

  if (any(has_lat) && is.finite(baseline$intercept) && is.finite(baseline$slope)) {
    exp_unwrapped <- baseline$intercept + baseline$slope * df$decimalLatitude[has_lat] +
                     elev_adj[has_lat] + mat_adj[has_lat]
    df$expected_doy[has_lat] <- wrap_doy_365(exp_unwrapped)
    if (is.finite(baseline$min_lat) && is.finite(baseline$pred_min_lat)) {
      hop_unwrapped <- baseline$pred_min_lat + 4.0 * (df$decimalLatitude[has_lat] - baseline$min_lat) +
                       elev_adj[has_lat]
      df$hopkins_expected[has_lat] <- wrap_doy_365(hop_unwrapped)
    }
  }

  df$delta_doy <- NA_real_
  has_both <- valid_doy & !is.na(df$expected_doy) & is.finite(df$expected_doy)
  if (any(has_both)) {
    df$delta_doy[has_both] <- circular_signed_doy_offset(df$doy[has_both], df$expected_doy[has_both])
  }

  # Hopkins-adjusted flowering DOY & circular coordinates (sin_doy, cos_doy) standardized to empirical cohort median latitude
  med_cohort_lat <- if (any(has_lat)) stats::median(df$decimalLatitude[has_lat], na.rm = TRUE) else NA_real_
  ref_doy_base <- if (is.finite(med_cohort_lat) && is.finite(baseline$intercept) && is.finite(baseline$slope)) {
    wrap_doy_365(baseline$intercept + baseline$slope * med_cohort_lat)
  } else {
    circ_mean_doy
  }
  df$hopkins_adjusted_doy <- ifelse(
    has_both & is.finite(ref_doy_base),
    wrap_doy_365(ref_doy_base + df$delta_doy),
    ifelse(valid_doy, df$doy, NA_real_)
  )
  valid_hdoy <- !is.na(df$hopkins_adjusted_doy) & is.finite(df$hopkins_adjusted_doy)
  df$sin_doy <- ifelse(valid_hdoy, round(sin(2.0 * pi * df$hopkins_adjusted_doy / 365.25), 6), NA_real_)
  df$cos_doy <- ifelse(valid_hdoy, round(cos(2.0 * pi * df$hopkins_adjusted_doy / 365.25), 6), NA_real_)

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
    export_df <- df[, intersect(keep_cols, names(df)), drop = FALSE]
    if (write_csv_atomic(export_df, output_csv, demo = demo)) {
      message(sprintf("Exported %d voucher phenological anomalies to %s.", nrow(export_df), output_csv))
    }
  }

  # Statistical divergence testing across target taxa strictly on verified anthesis cohort
  active_tx <- get_active_packera_taxa(df$species_standardized)
  target_flowering <- df[df$species_standardized %in% active_tx &
                         !is.na(df$delta_doy) & df$is_verified_anthesis, , drop = FALSE]

  if (nrow(target_flowering) >= 5L && length(unique(target_flowering$species_standardized)) > 1L) {
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
    P_Value = if (!is.na(baseline$p_value)) sprintf("%.4e", baseline$p_value) else "NA",
    Interpretation = sprintf("Empirical anthesis advance slope = +%.2f days/degree N (Hopkins' Law)", baseline$slope),
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
    if (length(sub_vals) > 0L) {
      m <- mean(sub_vals); s <- stats::sd(sub_vals); med <- stats::median(sub_vals); iq <- stats::IQR(sub_vals)
      n_tx <- length(sub_vals); se <- if (is.finite(s) && n_tx > 0L) s / sqrt(n_tx) else NA_real_
      summary_records[[length(summary_records) + 1L]] <- data.frame(
        Analysis_Section = "Taxon_Distribution", Taxon_or_Comparison = tx, Sample_Size = n_tx,
        Mean_Delta_DOY = round(m, 3), SD_Delta_DOY = round(s, 3), Median_Delta_DOY = round(med, 3),
        IQR_Delta_DOY = round(iq, 3),
        CI_Lower_95 = if (is.finite(se)) round(m - 1.96 * se, 3) else NA_real_,
        CI_Upper_95 = if (is.finite(se)) round(m + 1.96 * se, 3) else NA_real_,
        P_Value = "",
        Interpretation = if (m < -3) "Early blooming relative to latitude" else if (m > 3) "Late blooming relative to latitude" else "Conforms to latitudinal cline",
        stringsAsFactors = FALSE
      )
    }
  }
  if (nrow(tukey_mat) > 0L) {
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
      summary_records[[length(summary_records) + 1L]] <- data.frame(
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

evaluate_circular_phenology <- compute_phenological_anomalies

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
  if (length(active_tx) == 0L) active_tx <- get_active_packera_taxa(df$species_standardized)
  is_target <- df$plot_group %in% active_tx & !is.na(df$delta_doy)
  target_df <- if ("is_verified_anthesis" %in% names(df)) {
    df[is_target & df$is_verified_anthesis, , drop = FALSE]
  } else {
    df[is_target & !is.na(df$doy) & df$doy >= 1.0 & df$doy <= 366.0, , drop = FALSE]
  }
  if (nrow(target_df) == 0L) target_df <- df[is_target, , drop = FALSE]
  if (nrow(target_df) == 0L) {
    message("[WARN] Zero records with valid delta_doy available for phenology figure; skipping plot export.")
    return(invisible(NULL))
  }

  thm <- ggplot2::theme_bw(base_size = 9.5) +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                        plot.subtitle = ggplot2::element_text(size = 7.5, color = "grey25"),
                        panel.grid.minor = ggplot2::element_blank())

  lat_col <- if ("decimalLatitude" %in% names(target_df) && any(!is.na(target_df$decimalLatitude))) target_df$decimalLatitude else target_df$latitude
  target_df$decimalLatitude <- suppressWarnings(as.numeric(lat_col))
  target_df <- target_df[is.finite(target_df$decimalLatitude) & is.finite(target_df$doy), , drop = FALSE]
  if (nrow(target_df) == 0L) return(invisible(NULL))
  keep_na_lat <- is.finite(target_df$decimalLatitude) & target_df$decimalLatitude >= 22.0 & target_df$decimalLatitude <= 62.0
  if (sum(keep_na_lat) >= 3L) target_df <- target_df[keep_na_lat, , drop = FALSE]

  active_tx <- sort(unique(as.character(target_df$plot_group)))
  fg_tx <- active_tx
  bg_tx <- character(0)
  if (length(active_tx) > 6L) {
    tier_info <- if (exists("resolve_taxon_display_tier", mode = "function")) {
      tryCatch(resolve_taxon_display_tier(target_df, target_df$plot_group), error = function(e) NULL)
    } else {
      NULL
    }
    if (is.list(tier_info) && !is.null(tier_info$foreground) && length(intersect(tier_info$foreground, active_tx)) > 0L) {
      fg_tx <- intersect(as.character(tier_info$foreground), active_tx)
    } else if (is.list(tier_info) && !is.null(tier_info$fg_taxa) && length(intersect(tier_info$fg_taxa, active_tx)) > 0L) {
      fg_tx <- intersect(as.character(tier_info$fg_taxa), active_tx)
    } else if (is.character(tier_info) && length(intersect(tier_info, active_tx)) > 0L && length(intersect(tier_info, active_tx)) <= 6L) {
      fg_tx <- intersect(tier_info, active_tx)
    } else {
      dubia_tx <- active_tx[grepl("dubia", active_tx, ignore.case = TRUE)]
      other_tx <- setdiff(active_tx, dubia_tx)
      cnt_other <- vapply(other_tx, function(tx) sum(target_df$plot_group == tx, na.rm = TRUE), integer(1))
      other_ord <- other_tx[order(-cnt_other, other_tx)]
      n_keep_oth <- max(0L, 6L - length(dubia_tx))
      fg_tx <- c(dubia_tx, utils::head(other_ord, n_keep_oth))
    }
    bg_tx <- setdiff(active_tx, fg_tx)
  }

  dyn_pal <- resolve_dynamic_group_palette(active_tx)
  if (length(bg_tx) > 0L) {
    fg_pal <- resolve_dynamic_group_palette(fg_tx)
    for (tx in fg_tx) dyn_pal[[tx]] <- fg_pal[[tx]]
    for (tx in bg_tx) dyn_pal[[tx]] <- "#B0B0B0"
  }
  dyn_shapes <- stats::setNames(rep_len(c(21, 24, 22, 23, 1, 2, 5, 0), length(active_tx)), active_tx)

  min_lat <- min(target_df$decimalLatitude, na.rm = TRUE)
  max_lat <- max(target_df$decimalLatitude, na.rm = TRUE)
  b0_use <- if (!is.null(baseline$intercept) && is.finite(baseline$intercept)) baseline$intercept else compute_circular_mean_doy(target_df$doy)
  b1_use <- if (!is.null(baseline$slope) && is.finite(baseline$slope)) baseline$slope else 0.0
  pred_min_lat <- b0_use + b1_use * min_lat
  lat_seq <- seq(min_lat, max_lat, length.out = 100)
  pred_line <- data.frame(
    latitude = lat_seq,
    decimalLatitude = lat_seq,
    doy = pmin(365.0, pmax(1.0, b0_use + b1_use * lat_seq)),
    hopkins_expected = pmin(365.0, pmax(1.0, pred_min_lat + 4.0 * (lat_seq - min_lat)))
  )

  pt_alpha <- if (nrow(target_df) > 1000L) 0.45 else 0.75
  if (length(active_tx) > 6L) {
    bg_df <- target_df[target_df$plot_group %in% bg_tx, , drop = FALSE]
    fg_df <- target_df[target_df$plot_group %in% fg_tx, , drop = FALSE]
    fg_df$plot_group <- factor(fg_df$plot_group, levels = fg_tx)
    fg_cnts <- table(fg_df$plot_group)
    reg_tx <- names(fg_cnts)[fg_cnts >= 3L]
    reg_df <- fg_df[fg_df$plot_group %in% reg_tx, , drop = FALSE]

    p1 <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = bg_df,
        ggplot2::aes(x = decimalLatitude, y = doy),
        shape = 1, color = "#B0B0B0", fill = NA, stroke = 0.65, alpha = 0.35, size = 1.45
      ) +
      ggplot2::geom_point(
        data = fg_df,
        ggplot2::aes(x = decimalLatitude, y = doy, color = plot_group, shape = plot_group),
        fill = NA, stroke = 0.80, alpha = pt_alpha, size = 1.70
      )
    if (nrow(reg_df) >= 3L) {
      p1 <- p1 + ggplot2::geom_smooth(
        data = reg_df,
        ggplot2::aes(x = decimalLatitude, y = doy, color = plot_group),
        method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.75, linetype = "solid", show.legend = FALSE
      )
    }
    p1 <- p1 +
      ggplot2::geom_line(data = pred_line, ggplot2::aes(x = decimalLatitude, y = doy),
                         color = "black", linewidth = 1.1, linetype = "solid", inherit.aes = FALSE) +
      ggplot2::geom_line(data = pred_line, ggplot2::aes(x = decimalLatitude, y = hopkins_expected),
                         color = "#D55E00", linewidth = 0.95, linetype = "dashed", inherit.aes = FALSE) +
      ggplot2::scale_color_manual(values = dyn_pal[fg_tx], breaks = fg_tx) +
      ggplot2::scale_shape_manual(values = dyn_shapes[fg_tx], breaks = fg_tx) +
      thm +
      ggplot2::labs(
        title = sprintf("A. Verified Latitudinal Anthesis Cline: DOY ~ Latitude (Slope = +%.2f d/deg, R2 = %.2f)",
                        b1_use, ifelse(is.finite(baseline$r_squared), baseline$r_squared, 0.0)),
        subtitle = "Solid black = Empirical baseline; Dashed vermilion = Hopkins' Law (+4.0 d/deg N); Colored lines = Foreground species clines",
        x = "Decimal Latitude (deg N)", y = "Observed Anthesis Day of Year (DOY)", color = "Foreground Taxa", shape = "Foreground Taxa"
      )
  } else {
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
        title = sprintf("A. Verified Latitudinal Anthesis Cline: DOY ~ Latitude (Slope = +%.2f d/deg, R2 = %.2f)",
                        b1_use, ifelse(is.finite(baseline$r_squared), baseline$r_squared, 0.0)),
        subtitle = "Solid black = Empirical robust regression; Dashed vermilion = Hopkins' Law (+4.0 d / deg N)",
        x = "Decimal Latitude (deg N)", y = "Observed Anthesis Day of Year (DOY)", color = "Group", shape = "Group"
      )
  }
  if ("is_type_specimen" %in% names(target_df) && any(as.logical(target_df$is_type_specimen) == TRUE, na.rm = TRUE)) {
    p1 <- p1 + ggplot2::geom_point(
      data = target_df[which(as.logical(target_df$is_type_specimen) == TRUE), , drop = FALSE],
      ggplot2::aes(x = decimalLatitude, y = doy),
      shape = 5, fill = NA, stroke = 1.3, size = 2.8, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
    )
  }

  df_ridge <- target_df[is.finite(target_df$delta_doy), , drop = FALSE]
  med_df <- stats::aggregate(delta_doy ~ plot_group, data = df_ridge, FUN = function(x) stats::median(x, na.rm = TRUE))
  ord_lvls <- as.character(med_df$plot_group[order(med_df$delta_doy)])
  df_ridge$plot_group <- factor(as.character(df_ridge$plot_group), levels = ord_lvls)
  med_df <- med_df[order(med_df$delta_doy), , drop = FALSE]
  med_df$plot_group <- factor(as.character(med_df$plot_group), levels = ord_lvls)
  med_sub_df <- if (nrow(med_df) > 6L) med_df[as.character(med_df$plot_group) %in% fg_tx, , drop = FALSE] else med_df
  sub_p2_txt <- paste(sprintf("%s: %+.1f d", gsub("^Packera\\s+", "P. ", as.character(med_sub_df$plot_group)), med_sub_df$delta_doy), collapse = " | ")

  if (requireNamespace("ggridges", quietly = TRUE)) {
    p2 <- ggplot2::ggplot(df_ridge, ggplot2::aes(x = delta_doy, y = plot_group, fill = plot_group, color = plot_group)) +
      ggridges::geom_density_ridges(alpha = 0.45, scale = 1.15, rel_min_height = 0.01, linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = c(-14, 14), linetype = "dotted", color = "#D55E00", linewidth = 0.78) +
      ggplot2::scale_fill_manual(values = dyn_pal) +
      ggplot2::scale_color_manual(values = dyn_pal) +
      thm +
      ggplot2::theme(legend.position = if (length(active_tx) > 6L) "none" else "right") +
      ggplot2::labs(
        title = "B. Phenological Anomalies (Delta DOY) & Allochronic Separation",
        subtitle = sub_p2_txt,
        x = "Phenological Anomaly Delta DOY (Days Relative to Latitudinal Baseline; Dotted = +/-14 d)", y = "Target Group (Ordered by Median Delta DOY)", fill = "Group", color = "Group"
      )
  } else {
    p2 <- ggplot2::ggplot(df_ridge, ggplot2::aes(x = delta_doy, fill = plot_group, color = plot_group)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 0.85) +
      ggplot2::geom_vline(data = med_df, ggplot2::aes(xintercept = delta_doy, color = plot_group),
                          linetype = "dashed", linewidth = 0.85, show.legend = FALSE) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = c(-14, 14), linetype = "dotted", color = "#D55E00", linewidth = 0.78) +
      ggplot2::scale_fill_manual(values = dyn_pal) +
      ggplot2::scale_color_manual(values = dyn_pal) +
      thm +
      ggplot2::theme(legend.position = if (length(active_tx) > 6L) "none" else "right") +
      ggplot2::labs(
        title = "B. Phenological Anomalies (Delta DOY) & Allochronic Separation",
        subtitle = sub_p2_txt,
        x = "Phenological Anomaly Delta DOY (Days Relative to Latitudinal Baseline; Dotted = +/-14 d)", y = "Kernel Density", fill = "Group", color = "Group"
      )
  }

  fig_w <- 12.5
  fig_h <- max(5.5, 0.32 * length(active_tx) + 2.2)
  grob_pheno <- gridExtra::arrangeGrob(p1, p2, ncol = 2)

  tmp_pdf <- paste0(out_pdf, ".tmp")
  grDevices::pdf(tmp_pdf, width = fig_w, height = fig_h)
  tryCatch({
    grid::grid.draw(grob_pheno)
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, out_pdf)

  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  if (identical(out_png, out_pdf)) out_png <- paste0(out_pdf, ".png")
  tmp_png <- paste0(out_png, ".tmp")
  grDevices::png(tmp_png, width = fig_w, height = fig_h, units = "in", res = 300)
  tryCatch({
    grid::grid.draw(grob_pheno)
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(tmp_png)) file.rename(tmp_png, out_png)
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

compute_geo_knn_kappa <- function(coords_lon_lat, cluster_vec, k_nn = 5L, max_dist_km = 180.0) {
  lon <- suppressWarnings(as.numeric(coords_lon_lat[, 1]))
  lat <- suppressWarnings(as.numeric(coords_lon_lat[, 2]))
  cl <- as.integer(factor(cluster_vec))
  ok <- is.finite(lon) & is.finite(lat) & !is.na(cl)
  if (sum(ok) < 8L || length(unique(cl[ok])) < 2L) return(NA_real_)
  lon <- lon[ok]; lat <- lat[ok]; cl <- cl[ok]
  n <- length(cl)
  k_use <- min(as.integer(k_nn), n - 1L)
  if (k_use < 1L) return(NA_real_)

  rad <- pi / 180.0
  lat_r <- lat * rad; lon_r <- lon * rad
  dlat <- outer(lat_r, lat_r, "-")
  dlon <- outer(lon_r, lon_r, "-")
  hav <- sin(dlat * 0.5)^2 + outer(cos(lat_r), cos(lat_r), "*") * sin(dlon * 0.5)^2
  dist_km <- matrix(2.0 * 6371.0 * asin(pmin(1.0, sqrt(pmax(0.0, hav)))), nrow = n, ncol = n)
  diag(dist_km) <- Inf

  nn_maj <- integer(n)
  sympatric_conflict <- logical(n)
  for (i in seq_len(n)) {
    d_i <- dist_km[i, ]
    ord_i <- order(d_i)[seq_len(k_use)]
    near_i <- ord_i[is.finite(d_i[ord_i]) & d_i[ord_i] <= max_dist_km]
    if (length(near_i) == 0L) near_i <- ord_i[1L]
    tbl_i <- table(cl[near_i])
    nn_maj[i] <- as.integer(names(tbl_i)[which.max(tbl_i)][1L])
    close_10km <- which(d_i <= 10.0)
    sympatric_conflict[i] <- length(close_10km) > 0L && any(cl[close_10km] != cl[i])
  }

  p_obs <- mean(nn_maj == cl)
  lvls <- sort(unique(cl))
  p_exp <- sum(vapply(lvls, function(lv) mean(cl == lv) * mean(nn_maj == lv), numeric(1)))
  kappa_geo <- if ((1.0 - p_exp) > 1e-6) (p_obs - p_exp) / (1.0 - p_exp) else 0.0
  attr(kappa_geo, "sympatric_10km_rate") <- mean(sympatric_conflict)
  kappa_geo
}

evaluate_phenotypic_ecotype_coupling <- function(df, cluster_vec) {
  n <- nrow(df)
  cl_fac <- factor(cluster_vec)
  if (n < 10L || length(levels(cl_fac)) < 2L) {
    return(list(has_independent_modality = FALSE, pass = FALSE, morph_p = NA_real_, morph_r2_adj = NA_real_, pheno_p = NA_real_, pheno_diff_days = 0.0))
  }
  morph_cands <- unique(c(
    grep("^(Sym_PC|ECT_PC)[1-4]$", names(df), value = TRUE),
    intersect(c("can1", "can2", "blade_length_to_width_ratio", "petiole_to_blade_ratio", "solidity", "dissection_index"), names(df))
  ))
  morph_cols <- morph_cands[vapply(morph_cands, function(cn) {
    v <- suppressWarnings(as.numeric(df[[cn]]))
    sum(is.finite(v)) >= 8L && stats::sd(v[is.finite(v)]) > 1e-6
  }, logical(1))]

  morph_p <- NA_real_
  morph_r2_adj <- NA_real_
  if (length(morph_cols) >= 1L) {
    Y_m <- scale(as.matrix(df[, morph_cols, drop = FALSE]))
    for (j in seq_len(ncol(Y_m))) {
      bad <- !is.finite(Y_m[, j])
      if (any(bad)) Y_m[bad, j] <- 0.0
    }
    lat_v <- if ("latitude" %in% names(df)) suppressWarnings(as.numeric(df$latitude)) else rep(0.0, n)
    lon_v <- if ("longitude" %in% names(df)) suppressWarnings(as.numeric(df$longitude)) else rep(0.0, n)
    lat_v[!is.finite(lat_v)] <- stats::median(lat_v, na.rm = TRUE)
    lon_v[!is.finite(lon_v)] <- stats::median(lon_v, na.rm = TRUE)

    # Partial multivariate F-test of morphology ~ cluster | (latitude + longitude)
    X0 <- stats::model.matrix(~ lat_v + lon_v)
    X1 <- stats::model.matrix(~ lat_v + lon_v + cl_fac)
    rss0 <- sum( qr.resid(qr(X0), Y_m)^2 )
    rss1 <- sum( qr.resid(qr(X1), Y_m)^2 )
    df1 <- ncol(X1) - ncol(X0)
    df2 <- max(1L, n - ncol(X1))
    if (rss0 > 1e-8 && df1 >= 1L && df2 >= 4L) {
      r2_part <- max(0.0, (rss0 - rss1) / rss0)
      morph_r2_adj <- 1.0 - (1.0 - r2_part) * ((n - ncol(X0)) / df2)
      f_val <- ((rss0 - rss1) / df1) / max(rss1 / df2, 1e-8)
      morph_p <- stats::pf(f_val, df1 = df1 * ncol(Y_m), df2 = df2 * ncol(Y_m), lower.tail = FALSE)
    }
  }

  pheno_p <- NA_real_
  pheno_diff <- 0.0
  has_pheno <- "delta_doy" %in% names(df) && sum(is.finite(suppressWarnings(as.numeric(df$delta_doy)))) >= 10L
  if (has_pheno) {
    dd <- suppressWarnings(as.numeric(df$delta_doy))
    ok_d <- is.finite(dd) & !is.na(cl_fac)
    if (sum(ok_d) >= 10L && length(unique(cl_fac[ok_d])) >= 2L) {
      meds <- tapply(dd[ok_d], cl_fac[ok_d], stats::median, na.rm = TRUE)
      pheno_diff <- max(meds, na.rm = TRUE) - min(meds, na.rm = TRUE)
      pheno_p <- tryCatch(stats::kruskal.test(dd[ok_d] ~ droplevels(cl_fac[ok_d]))$p.value, error = function(e) NA_real_)
    }
  }

  has_indep <- (length(morph_cols) >= 1L) || has_pheno
  pass_morph <- is.finite(morph_p) && is.finite(morph_r2_adj) && (morph_p < 0.01) && (morph_r2_adj >= 0.05)
  pass_pheno <- is.finite(pheno_p) && (pheno_p < 0.01) && (pheno_diff >= 10.0)
  list(
    has_independent_modality = has_indep,
    pass = isTRUE(pass_morph || pass_pheno),
    morph_p = morph_p,
    morph_r2_adj = morph_r2_adj,
    pheno_p = pheno_p,
    pheno_diff_days = pheno_diff
  )
}

compute_isolation_by_environment <- function(df, edaph_vars = character(0), clim_vars = character(0), mem_mat = NULL, n_perm = 99L) {
  # Non-circular Isolation-by-Environment (IBE; Wang & Bradburd 2014) variance decomposition
  # testing continuous Morphological + Phenological divergence as a function of
  # [A] Edaphics and [B] Macro-Climate conditioned on [C] Spatial MEMs.
  n <- nrow(df)
  empty_ibe <- data.frame(
    Modality = character(0), N_Specimens = integer(0),
    Pure_Edaphic_AdjR2 = numeric(0), Pure_Climate_AdjR2 = numeric(0),
    Pure_Space_MEM_AdjR2 = numeric(0), Total_Explained_AdjR2 = numeric(0),
    Partial_RDA_F = numeric(0), Partial_RDA_p = numeric(0),
    IBE_Inference = character(0), stringsAsFactors = FALSE
  )
  if (is.null(df) || n < 8L) return(empty_ibe)

  ed_cols <- intersect(as.character(edaph_vars), names(df))
  cl_cols <- intersect(setdiff(as.character(clim_vars), "delta_doy"), names(df))
  clean_block <- function(cols) {
    if (length(cols) == 0L) return(NULL)
    M <- matrix(0.0, nrow = n, ncol = length(cols), dimnames = list(NULL, cols))
    keep <- logical(length(cols))
    for (j in seq_along(cols)) {
      v <- suppressWarnings(as.numeric(df[[cols[j]]]))
      fin <- is.finite(v)
      if (sum(fin) >= 5L && stats::sd(v[fin]) > 1e-8) {
        v[!fin] <- stats::median(v[fin], na.rm = TRUE)
        M[, j] <- as.numeric(scale(v))
        keep[j] <- TRUE
      }
    }
    if (!any(keep)) return(NULL)
    M[, keep, drop = FALSE]
  }

  X_ed <- clean_block(ed_cols)
  X_cl <- clean_block(cl_cols)
  if (is.null(X_ed) && is.null(X_cl)) return(empty_ibe)

  if (is.null(mem_mat) || nrow(as.matrix(mem_mat)) != n) {
    if (all(c("longitude", "latitude") %in% names(df))) {
      lon <- suppressWarnings(as.numeric(df$longitude)); lat <- suppressWarnings(as.numeric(df$latitude))
      lon[!is.finite(lon)] <- stats::median(lon, na.rm = TRUE)
      lat[!is.finite(lat)] <- stats::median(lat, na.rm = TRUE)
      mem_mat <- scale(cbind(MEM_Lon = lon, MEM_Lat = lat))
    } else {
      mem_mat <- matrix(0.0, nrow = n, ncol = 1L, dimnames = list(NULL, "MEM1"))
    }
  } else {
    mem_mat <- scale(as.matrix(mem_mat))
    mem_mat[!is.finite(mem_mat)] <- 0.0
  }

  eval_response_block <- function(Y_raw, mod_name) {
    if (is.null(Y_raw) || ncol(Y_raw) == 0L) return(NULL)
    Y <- scale(as.matrix(Y_raw))
    Y[!is.finite(Y)] <- 0.0
    X_env <- cbind(X_ed, X_cl)

    adj_r2_ols <- function(X_mat) {
      if (is.null(X_mat) || ncol(X_mat) == 0L) return(0.0)
      X_d <- cbind(1.0, X_mat)
      qr_x <- qr(X_d)
      p_k <- qr_x$rank - 1L
      if (p_k < 1L || (n - qr_x$rank) < 2L) return(0.0)
      rss <- sum(qr.resid(qr_x, Y)^2)
      tss <- sum(sweep(Y, 2L, colMeans(Y), "-")^2)
      if (tss <= 1e-10) return(0.0)
      r2 <- max(0.0, 1.0 - rss / tss)
      1.0 - (1.0 - r2) * ((n - 1.0) / max(1.0, n - p_k - 1.0))
    }

    r2_all <- adj_r2_ols(cbind(X_ed, X_cl, mem_mat))
    r2_no_ed <- adj_r2_ols(cbind(X_cl, mem_mat))
    r2_no_cl <- adj_r2_ols(cbind(X_ed, mem_mat))
    r2_no_sp <- adj_r2_ols(cbind(X_ed, X_cl))

    pure_ed <- max(0.0, r2_all - r2_no_ed)
    pure_cl <- max(0.0, r2_all - r2_no_cl)
    pure_sp <- max(0.0, r2_all - r2_no_sp)

    f_obs <- NA_real_; p_perm <- NA_real_
    if (requireNamespace("vegan", quietly = TRUE) && ncol(X_env) >= 1L && (n - ncol(X_env) - ncol(mem_mat)) >= 4L) {
      rda_res <- tryCatch({
        fit_r <- vegan::rda(Y ~ X_env + Condition(mem_mat))
        an_r <- suppressMessages(suppressWarnings(vegan::anova.cca(fit_r, permutations = as.integer(max(19L, n_perm)))))
        list(F_val = as.numeric(an_r$F[1L]), p_val = as.numeric(an_r$`Pr(>F)`[1L]))
      }, error = function(e) NULL)
      if (!is.null(rda_res)) {
        f_obs <- rda_res$F_val; p_perm <- rda_res$p_val
      }
    }
    if (!is.finite(p_perm)) {
      X0 <- cbind(1.0, mem_mat)
      X1 <- cbind(1.0, mem_mat, X_env)
      q0 <- qr(X0); q1 <- qr(X1)
      rss0 <- sum(qr.resid(q0, Y)^2); rss1 <- sum(qr.resid(q1, Y)^2)
      df1 <- max(1L, q1$rank - q0$rank); df2 <- max(1L, n - q1$rank)
      f_obs <- ((rss0 - rss1) / df1) / max(rss1 / df2, 1e-8)
      Y_res0 <- qr.resid(q0, Y)
      f_null <- numeric(n_perm)
      for (b in seq_len(n_perm)) {
        Yp <- Y_res0[sample.int(n), , drop = FALSE]
        rss0_b <- sum(qr.resid(q0, Yp)^2)
        rss1_b <- sum(qr.resid(q1, Yp)^2)
        f_null[b] <- ((rss0_b - rss1_b) / df1) / max(rss1_b / df2, 1e-8)
      }
      p_perm <- (sum(f_null >= f_obs) + 1.0) / (n_perm + 1.0)
    }

    pure_env <- pure_ed + pure_cl
    inf_lbl <- if (is.finite(p_perm) && p_perm < 0.01 && pure_env >= 0.05) {
      "Significant Isolation by Environment (Adaptive Divergence Candidate)"
    } else if (is.finite(p_perm) && p_perm < 0.05) {
      "Weak Environmental Cline (Plasticity / Gradual Ecotone)"
    } else {
      "Null IBE (Panmictic / Facultative Niche Breadth)"
    }

    data.frame(
      Modality = mod_name,
      N_Specimens = as.integer(n),
      Pure_Edaphic_AdjR2 = round(pure_ed, 4),
      Pure_Climate_AdjR2 = round(pure_cl, 4),
      Pure_Space_MEM_AdjR2 = round(pure_sp, 4),
      Total_Explained_AdjR2 = round(max(0.0, r2_all), 4),
      Partial_RDA_F = round(f_obs, 4),
      Partial_RDA_p = round(p_perm, 4),
      IBE_Inference = inf_lbl,
      stringsAsFactors = FALSE
    )
  }

  rows_out <- list()
  morph_cols <- unique(c(
    grep("^(Sym_PC|ECT_PC)[1-4]$", names(df), value = TRUE),
    intersect(c("blade_length_to_width_ratio", "petiole_to_blade_ratio", "solidity", "dissection_index"), names(df))
  ))
  Y_morph <- clean_block(morph_cols)
  if (!is.null(Y_morph)) rows_out[[length(rows_out) + 1L]] <- eval_response_block(Y_morph, "Foliar_Morphology_IBE")

  Y_pheno <- clean_block(intersect("delta_doy", names(df)))
  if (!is.null(Y_pheno)) rows_out[[length(rows_out) + 1L]] <- eval_response_block(Y_pheno, "Anthesis_Phenology_IBE")

  if (length(rows_out) == 0L) return(empty_ibe)
  do.call(rbind, rows_out)
}

resolve_analysis_group <- function(
  df,
  intraspecific_mode = FALSE,
  primary_env_vars = NULL,
  rf_var_weights = NULL,
  allow_unimodal_null = NULL
) {
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  if (!is_intra) {
    df$analysis_group <- as.character(df$species_standardized)
    df$is_missing_gmm_cluster <- rep(FALSE, nrow(df))
    return(df)
  }

  has_step04_gmm <- "gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster) & nzchar(trimws(as.character(df$gmm_cluster))))
  used_gmm <- FALSE
  if (has_step04_gmm) {
    miss_gmm <- is.na(df$gmm_cluster) | !nzchar(trimws(as.character(df$gmm_cluster)))
    g_raw <- ifelse(miss_gmm, 1L, suppressWarnings(as.integer(df$gmm_cluster)))
    g_raw[is.na(g_raw)] <- 1L
    grp_cand <- paste0("Morphotype_K", g_raw)
    tbl <- table(grp_cand)
    if (length(tbl) >= 2 && sum(tbl >= 5) >= 2) {
      df$analysis_group <- as.character(grp_cand)
      df$is_missing_gmm_cluster <- as.logical(miss_gmm)
      used_gmm <- TRUE
      message(sprintf("[Intraspecific Mode] Using independent Step 04 GMM morphotypes as analysis_group (%s).",
                      paste(names(tbl), tbl, sep = "=", collapse = ", ")))
    }
  }

  if (!used_gmm) {
    df$is_missing_gmm_cluster <- rep(FALSE, nrow(df))
    n_rows <- nrow(df)
    if (n_rows == 0L) {
      df$analysis_group <- character(0)
      return(df)
    }

    sp_std <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Packera dubia", n_rows)
    valid_sp <- sp_std[!is.na(sp_std) & nzchar(trimws(sp_std)) & sp_std != "Unknown"]
    canon_sp <- if (length(valid_sp) > 0L) valid_sp[1L] else "Packera dubia"

    has_coords <- all(c("latitude", "longitude") %in% names(df)) &&
                  sum(is.finite(suppressWarnings(as.numeric(df$latitude))) & is.finite(suppressWarnings(as.numeric(df$longitude)))) >= 10L
    if (is.null(allow_unimodal_null)) {
      env_flag <- Sys.getenv("PACKERA_ALLOW_UNIMODAL_ECOTYPE", unset = "")
      allow_unimodal_null <- if (nzchar(env_flag)) {
        !(tolower(env_flag) %in% c("0", "false", "no"))
      } else {
        isTRUE(has_step04_gmm || has_coords)
      }
    }

    # Exclude coordinate columns AND binary/synthetic indicator columns from continuous environmental PCA
    exclude_cols <- c(
      "catalogNumber", "plant_individual_id", "latitude", "longitude",
      "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
      "coordinateUncertaintyInMeters", "buffer_radius_m",
      "high_coord_uncertainty_flag", "is_training_eligible_coord",
      "is_clim_eligible_coord", "is_edaphic_eligible_coord", "edaphic_coord_weight",
      "edaphic_best_fit_taxon", "edaphic_z_score", "climatic_z_score", "soil_silt_pct",
      "is_rock_outcrop", "is_azonal_lithology", "azonal_lithology_presence_buffer",
      "azonal_contact_dist_m", "resdept_r", "bdricm"
    )
    cand_cols <- if (!is.null(primary_env_vars) && length(primary_env_vars) > 0L) {
      setdiff(intersect(as.character(primary_env_vars), names(df)), exclude_cols)
    } else {
      character(0)
    }
    if (length(cand_cols) == 0L) {
      env_pattern <- paste0(
        "^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bio[0-9]*_|cfvo$|awc$|",
        "alr_|clr_|azonal_proximity_index$|ph_q10_buffer$|ph_q90_buffer$|",
        "min_bedrock_depth_buffer_cm$|max_cfvo_buffer_pct$|elevation$|elev_|clim_|edaph_|env_)"
      )
      cand_cols <- setdiff(grep(env_pattern, names(df), value = TRUE, ignore.case = TRUE), exclude_cols)
    }

    num_env <- list()
    min_fin <- min(5L, max(2L, n_rows))
    for (cl in cand_cols) {
      if (is.numeric(df[[cl]]) || is.integer(df[[cl]])) {
        v <- suppressWarnings(as.numeric(df[[cl]]))
        u_vals <- unique(round(v[is.finite(v)], 6))
        # Require strictly continuous/multi-valued environmental predictors (>= 4 distinct finite values when n >= 10)
        min_distinct <- if (n_rows >= 10L) 4L else 2L
        if (sum(is.finite(v)) >= min_fin && length(u_vals) >= min_distinct) {
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
    if ("fine_soil_ph" %in% names(num_env) && "ssurgo_ph" %in% names(num_env)) num_env[["ssurgo_ph"]] <- NULL
    if ("fine_soil_sand" %in% names(num_env) && "ssurgo_sand_pct" %in% names(num_env)) num_env[["ssurgo_sand_pct"]] <- NULL
    if ("fine_soil_clay" %in% names(num_env) && "ssurgo_clay_pct" %in% names(num_env)) num_env[["ssurgo_clay_pct"]] <- NULL
    if ("fine_soil_cec" %in% names(num_env) && "ssurgo_cec" %in% names(num_env)) num_env[["ssurgo_cec"]] <- NULL
    if ("awc" %in% names(num_env) && "ssurgo_aws_0_50cm" %in% names(num_env)) num_env[["ssurgo_aws_0_50cm"]] <- NULL

    if (length(num_env) == 0L || n_rows < 2L) {
      if (isTRUE(allow_unimodal_null)) {
        df$analysis_group <- rep(canon_sp, n_rows)
        message(sprintf("[Intraspecific Mode] Unimodal environmental variance; retaining single taxon '%s' (K=1).", canon_sp))
      } else {
        idx_half <- seq_len(n_rows) <= max(1L, floor(n_rows / 2))
        df$analysis_group <- ifelse(idx_half, "Ecotype_Cluster_1", "Ecotype_Cluster_2")
      }
      return(df)
    }

    E_mat <- do.call(cbind, num_env)
    # Prune collinear variables before PCA-env so redundant blocks cannot dominate PC1 (Dormann et al. 2013)
    if (ncol(E_mat) >= 3L && exists("screen_collinear_predictors", mode = "function")) {
      pruned_e <- tryCatch(screen_collinear_predictors(E_mat, cor.threshold = 0.70, vif.threshold = 5.0)$kept_vars, error = function(e) colnames(E_mat))
      if (length(intersect(colnames(E_mat), pruned_e)) >= 2L) {
        E_mat <- E_mat[, intersect(colnames(E_mat), pruned_e), drop = FALSE]
      }
    }
    # Never weight E_scaled by rf_var_weights inside resolve_analysis_group (prevents post-clustering circularity)
    E_scaled <- scale(E_mat, center = TRUE, scale = TRUE)
    E_scaled[!is.finite(E_scaled)] <- 0.0

    if (ncol(E_scaled) >= 2L) {
      pca_fit <- stats::prcomp(E_scaled, center = FALSE, scale. = FALSE)
      n_pcs <- min(2L, ncol(pca_fit$x))
      pca_scores <- as.matrix(pca_fit$x[, seq_len(n_pcs), drop = FALSE])
      if (ncol(pca_scores) == 1L) pca_scores <- cbind(PC1 = pca_scores[, 1], PC2 = rep(0.0, n_rows))
      pc1_load <- abs(pca_fit$rotation[, 1])
      top_var <- names(which.max(pc1_load))[1]
    } else {
      pca_scores <- cbind(PC1 = as.numeric(E_scaled[, 1]), PC2 = rep(0.0, n_rows))
      top_var <- colnames(E_scaled)[1]
    }

    n_distinct <- nrow(unique(round(pca_scores, 8)))
    if (n_distinct < 2L) {
      if (isTRUE(allow_unimodal_null)) {
        df$analysis_group <- rep(canon_sp, n_rows)
      } else {
        idx_half <- seq_len(n_rows) <= max(1L, floor(n_rows / 2))
        df$analysis_group <- ifelse(idx_half, "Ecotype_Cluster_1", "Ecotype_Cluster_2")
      }
      return(df)
    }

    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    set.seed(42)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)

    D_pca <- as.matrix(stats::dist(pca_scores))
    km2 <- stats::kmeans(pca_scores, centers = 2, nstart = 25)
    best_cl <- km2$cluster
    best_k <- 2L
    best_tot_within <- km2$tot.withinss
    tot_ss <- km2$totss
    sil2 <- compute_mean_silhouette_base(km2$cluster, D_pca)
    best_sil <- sil2

    if (n_rows >= 30L && n_distinct >= 3L) {
      km3 <- tryCatch(stats::kmeans(pca_scores, centers = 3, nstart = 25), error = function(e) NULL)
      if (!is.null(km3) && min(table(km3$cluster)) >= 5L) {
        sil3 <- compute_mean_silhouette_base(km3$cluster, D_pca)
        if (is.finite(sil3) && is.finite(sil2) && sil3 > (sil2 + 0.05)) {
          best_cl <- km3$cluster
          best_k <- 3L
          best_tot_within <- km3$tot.withinss
          best_sil <- sil3
        }
      }
    }

    # Non-Circular 3-Gate Validation of Candidate K >= 2 vs Unimodal K = 1 Null:
    #   Gate 1 (Continuous PCA-env Multimodality): silhouette >= 0.42 & variance reduction (1 - W_K/W_1) >= 0.50
    #   Gate 2 (Geographic Cohesion / Anti-Sympatric-Mixing): 5-NN spatial cluster kappa >= 0.55 & <= 15% 10km sympatric mixing
    #   Gate 3 (Independent Phenotypic Coupling when Step 04 morphology returned G = 1): significant morphological or allochronic divergence
    if (isTRUE(allow_unimodal_null)) {
      var_red <- if (is.finite(tot_ss) && tot_ss > 1e-8) (1.0 - best_tot_within / tot_ss) else 0.0
      pass_env_gap <- is.finite(best_sil) && (best_sil >= 0.42) && (var_red >= 0.50) && (min(table(best_cl)) >= 5L)

      geo_kappa <- if (has_coords) compute_geo_knn_kappa(df[, c("longitude", "latitude"), drop = FALSE], best_cl, k_nn = 5L) else NA_real_
      symp_10km <- if (!is.null(attr(geo_kappa, "sympatric_10km_rate"))) as.numeric(attr(geo_kappa, "sympatric_10km_rate")) else 0.0
      pass_geo <- !has_coords || (is.finite(geo_kappa) && geo_kappa >= 0.55 && symp_10km <= 0.15)

      pheno_cpl <- evaluate_phenotypic_ecotype_coupling(df, best_cl)
      # If Step 04 evaluated leaf morphology and found G = 1 unimodal morphology, require independent phenotypic coupling
      pass_pheno_gate <- if (isTRUE(has_step04_gmm) || isTRUE(pheno_cpl$has_independent_modality && has_step04_gmm)) {
        isTRUE(pheno_cpl$pass)
      } else {
        !isTRUE(pheno_cpl$has_independent_modality) || isTRUE(pheno_cpl$pass) || (is.finite(geo_kappa) && geo_kappa >= 0.75 && best_sil >= 0.55)
      }

      if (!(pass_env_gap && pass_geo && pass_pheno_gate)) {
        df$analysis_group <- rep(canon_sp, n_rows)
        attr(df, "ecotype_diagnostics") <- list(
          optimal_k = 1L,
          status = "Unimodal_Single_Lineage_Facultative_Breadth",
          candidate_k = best_k,
          pca_silhouette = round(best_sil, 4),
          variance_reduction = round(var_red, 4),
          geo_knn_kappa = round(geo_kappa, 4),
          sympatric_10km_rate = round(symp_10km, 4),
          morph_coupling_p = pheno_cpl$morph_p,
          morph_coupling_r2_adj = pheno_cpl$morph_r2_adj,
          pheno_coupling_p = pheno_cpl$pheno_p
        )
        message(sprintf(
          paste0(
            "[Intraspecific Mode: Non-Circular Ecotype Gate] Retained K=1 single taxon '%s' ",
            "(rejected sympatric/unimodal k-means split: sil=%.3f, var_red=%.3f, geo_kappa=%s, symp_10km=%.2f, morph_p=%s, pheno_p=%s)."
          ),
          canon_sp,
          ifelse(is.finite(best_sil), best_sil, NA_real_),
          var_red,
          ifelse(is.finite(geo_kappa), sprintf("%.3f", geo_kappa), "NA"),
          symp_10km,
          ifelse(is.finite(pheno_cpl$morph_p), sprintf("%.3g", pheno_cpl$morph_p), "NA"),
          ifelse(is.finite(pheno_cpl$pheno_p), sprintf("%.3g", pheno_cpl$pheno_p), "NA")
        ))
        return(df)
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
    message(sprintf("[Intraspecific Mode] Validated allopatric/divergent PCA-env ecotypes (K=%d, Axis-1 driver=%s: %s).",
                    best_k, clean_top_var, paste(names(tbl), tbl, sep = "=", collapse = ", ")))
  }
  return(df)
}

resolve_active_anchor_groups <- function(df, intraspecific_mode = FALSE) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(character(0))
  n_r <- nrow(df)
  is_intra <- isTRUE(intraspecific_mode) ||
    (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  g_raw <- if (is_intra && "analysis_group" %in% names(df)) {
    as.character(df$analysis_group)
  } else if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("analysis_group" %in% names(df)) {
    as.character(df$analysis_group)
  } else {
    rep(NA_character_, n_r)
  }
  elig <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, n_r)
  elig[is.na(elig)] <- TRUE
  if ("coordinateUncertaintyInMeters" %in% names(df)) {
    u_m <- suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
    elig <- elig & (is.na(u_m) | u_m <= 10000.0)
  }
  det_t <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n_r)
  is_prom <- rep(FALSE, n_r)
  if ("is_promoted_anchor" %in% names(df)) {
    pa <- as.logical(df$is_promoted_anchor); pa[is.na(pa)] <- FALSE; is_prom <- is_prom | pa
  }
  if ("is_promoted_bronze" %in% names(df)) {
    pb <- as.logical(df$is_promoted_bronze); pb[is.na(pb)] <- FALSE; is_prom <- is_prom | pb
  }
  miss_gmm <- if (is_intra && "is_missing_gmm_cluster" %in% names(df)) {
    as.logical(df$is_missing_gmm_cluster)
  } else if (is_intra && "analysis_group" %in% names(df) && all(grepl("^Morphotype_K", stats::na.omit(df$analysis_group))) && "gmm_cluster" %in% names(df)) {
    is.na(df$gmm_cluster) | !nzchar(trimws(as.character(df$gmm_cluster)))
  } else {
    rep(FALSE, n_r)
  }
  miss_gmm[is.na(miss_gmm)] <- FALSE
  base_anc <- ((det_t %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_prom) & elig & !is.na(g_raw) & nzchar(g_raw) & (g_raw != "Unknown")
  is_anc <- base_anc & !miss_gmm
  if (sum(is_anc) < 6L || length(unique(g_raw[is_anc])) < 2L) {
    is_anc <- (elig & !is.na(g_raw) & nzchar(g_raw) & (g_raw != "Unknown")) & !miss_gmm
  }
  if (sum(is_anc) < 6L || length(unique(g_raw[is_anc])) < 2L) {
    is_anc <- base_anc
  }
  if (sum(is_anc) < 6L || length(unique(g_raw[is_anc])) < 2L) {
    is_anc <- elig & !is.na(g_raw) & nzchar(g_raw) & (g_raw != "Unknown")
  }
  if (sum(is_anc) < 4L || length(unique(g_raw[is_anc])) < 2L) {
    is_anc <- !is.na(g_raw) & nzchar(g_raw) & (g_raw != "Unknown")
  }
  ifelse(is_anc, g_raw, NA_character_)
}

select_primary_predictors_from_rf <- function(
  rf_res = NULL,
  vif_res = NULL,
  df_env = NULL,
  cum_imp_target = 0.85,
  intraspecific_mode = FALSE
) {
  clim_pat <- "^(bio[0-9]*_|delta_doy$|doy$|hopkins_|sin_doy$|cos_doy$)"
  topo_pat <- "^(topo_|elevation$|elev_)"

  cand_vars <- character(0)
  if (!is.null(vif_res) && is.list(vif_res)) {
    v_kept <- if (!is.null(vif_res$kept_vars)) vif_res$kept_vars else vif_res$selected_variables
    if (length(v_kept) > 0L) cand_vars <- as.character(v_kept)
  }
  if (length(cand_vars) == 0L && !is.null(rf_res) && is.list(rf_res) && !is.null(rf_res$importance)) {
    cand_vars <- as.character(rf_res$importance$Predictor)
  }
  if (length(cand_vars) == 0L && !is.null(df_env) && is.data.frame(df_env)) {
    cand_vars <- discover_candidate_eco_predictors(df_env, include_pheno = TRUE, min_finite = 3L)
  }
  cand_vars <- unique(cand_vars[!is.na(cand_vars) & nzchar(cand_vars) & !grepl("^MEM[0-9]*$", cand_vars)])
  if (!is.null(df_env) && is.data.frame(df_env)) {
    cand_vars <- intersect(cand_vars, names(df_env))
    cand_vars <- cand_vars[vapply(cand_vars, function(cl) {
      v <- suppressWarnings(as.numeric(df_env[[cl]]))
      fin <- is.finite(v)
      sum(fin) >= min(3L, max(2L, nrow(df_env))) && is.finite(stats::sd(v[fin])) && stats::sd(v[fin]) > 1e-8
    }, logical(1))]
  }
  if (length(cand_vars) == 0L) {
    return(list(
      primary_env_vars = character(0),
      primary_clim_vars = character(0),
      primary_edaph_vars = character(0),
      primary_topo_vars = character(0),
      rf_var_weights = numeric(0),
      selection_method = "empty"
    ))
  }

  mda_map <- stats::setNames(rep(NA_real_, length(cand_vars)), cand_vars)
  if (!is.null(rf_res) && is.list(rf_res) && !is.null(rf_res$importance) &&
      all(c("Predictor", "MeanDecreaseAccuracy") %in% names(rf_res$importance))) {
    imp_df <- rf_res$importance
    hit_m <- match(cand_vars, as.character(imp_df$Predictor))
    ok_m <- !is.na(hit_m)
    mda_map[ok_m] <- suppressWarnings(as.numeric(imp_df$MeanDecreaseAccuracy[hit_m[ok_m]]))
  }

  pos_vars <- names(mda_map[is.finite(mda_map) & mda_map > 0.0])
  sel_method <- "rf_permutation_85pct"
  selected_vars <- character(0)
  raw_weights <- numeric(0)

  if (length(pos_vars) >= 1L) {
    ord_pos <- pos_vars[order(mda_map[pos_vars], decreasing = TRUE)]
    pos_vals <- mda_map[ord_pos]
    cum_frac <- cumsum(pos_vals) / max(sum(pos_vals), 1e-12)
    k_cut <- which(cum_frac >= cum_imp_target)[1]
    if (is.na(k_cut)) k_cut <- length(ord_pos)
    k_keep <- min(length(ord_pos), max(min(2L, length(ord_pos)), as.integer(k_cut)))
    selected_vars <- ord_pos[seq_len(k_keep)]
    raw_weights <- pos_vals[seq_len(k_keep)]
  } else {
    # Fallback when RF permutation importance has zero positive features or rf_res is unavailable:
    # univariate Kruskal-Wallis / ANOVA BH-adjusted p < 0.05 across active reference anchors
    sel_method <- "univariate_bh_fdr_fallback"
    anc_grp <- if (!is.null(df_env) && is.data.frame(df_env)) {
      resolve_active_anchor_groups(df_env, intraspecific_mode = intraspecific_mode)
    } else {
      character(0)
    }
    p_vals <- stats::setNames(rep(1.0, length(cand_vars)), cand_vars)
    sep_scores <- if (!is.null(df_env) && is.data.frame(df_env) && exists("compute_anchor_separation_scores", mode = "function")) {
      compute_anchor_separation_scores(df_env[, cand_vars, drop = FALSE], anchor_groups = anc_grp)
    } else {
      stats::setNames(seq(length(cand_vars), 1L), cand_vars)
    }
    if (!is.null(df_env) && is.data.frame(df_env) && length(anc_grp) == nrow(df_env)) {
      for (cl in cand_vars) {
        v <- suppressWarnings(as.numeric(df_env[[cl]]))
        ok <- is.finite(v) & !is.na(anc_grp) & nzchar(anc_grp)
        if (sum(ok) >= 4L && length(unique(anc_grp[ok])) >= 2L) {
          f_g <- factor(anc_grp[ok])
          pv_kw <- tryCatch(suppressWarnings(as.numeric(stats::kruskal.test(v[ok] ~ f_g)$p.value)), error = function(e) NA_real_)
          pv_av <- tryCatch({
            af <- summary(stats::aov(v[ok] ~ f_g))[[1]]
            suppressWarnings(as.numeric(af[["Pr(>F)"]][1]))
          }, error = function(e) NA_real_)
          pv <- min(c(pv_kw, pv_av), na.rm = TRUE)
          if (is.finite(pv)) p_vals[[cl]] <- pv
        }
      }
    }
    p_adj <- stats::p.adjust(p_vals, method = "BH")
    sig_vars <- names(p_adj[is.finite(p_adj) & p_adj < 0.05])
    if (length(sig_vars) >= 1L) {
      selected_vars <- sig_vars[order(sep_scores[sig_vars], decreasing = TRUE)]
    } else {
      sel_method <- "anchor_separation_top_fallback"
      ord_sep <- names(sort(sep_scores, decreasing = TRUE))
      selected_vars <- ord_sep[seq_len(min(length(ord_sep), max(2L, min(6L, length(ord_sep)))))]
    }
    raw_weights <- pmax(as.numeric(sep_scores[selected_vars]), 1e-4)
    names(raw_weights) <- selected_vars
  }

  primary_clim <- grep(clim_pat, selected_vars, value = TRUE, ignore.case = TRUE)
  primary_topo <- setdiff(grep(topo_pat, selected_vars, value = TRUE, ignore.case = TRUE), primary_clim)
  primary_edaph <- setdiff(selected_vars, unique(c(primary_clim, primary_topo)))

  # If positive-importance predictors existed in a block (climatic or edaphic) just outside the 85% cutoff,
  # retain the top positive predictor of that block so partitioned climatic vs edaphic tests remain empirically grounded
  cand_clim_pos <- intersect(grep(clim_pat, cand_vars, value = TRUE, ignore.case = TRUE), pos_vars)
  if (length(primary_clim) == 0L && length(cand_clim_pos) > 0L) {
    best_c <- cand_clim_pos[which.max(mda_map[cand_clim_pos])][1]
    primary_clim <- best_c
    selected_vars <- unique(c(selected_vars, best_c))
    raw_weights[[best_c]] <- max(mda_map[[best_c]], 1e-4)
  }
  cand_edaph_all <- setdiff(cand_vars, unique(c(
    grep(clim_pat, cand_vars, value = TRUE, ignore.case = TRUE),
    grep(topo_pat, cand_vars, value = TRUE, ignore.case = TRUE)
  )))
  cand_edaph_pos <- intersect(cand_edaph_all, pos_vars)
  if (length(primary_edaph) == 0L && length(cand_edaph_pos) > 0L) {
    best_e <- cand_edaph_pos[which.max(mda_map[cand_edaph_pos])][1]
    primary_edaph <- best_e
    selected_vars <- unique(c(selected_vars, best_e))
    raw_weights[[best_e]] <- max(mda_map[[best_e]], 1e-4)
  }

  raw_weights <- raw_weights[selected_vars]
  raw_weights[!is.finite(raw_weights) | raw_weights <= 0.0] <- 1e-4
  rf_var_weights <- raw_weights / sum(raw_weights)
  names(rf_var_weights) <- selected_vars

  list(
    primary_env_vars = selected_vars,
    primary_clim_vars = primary_clim,
    primary_edaph_vars = primary_edaph,
    primary_topo_vars = primary_topo,
    rf_var_weights = rf_var_weights,
    selection_method = sel_method
  )
}

compute_specimen_edaphic_z <- function(
  df,
  primary_edaph_vars = NULL,
  primary_clim_vars = NULL,
  primary_env_vars = NULL,
  rf_var_weights = NULL,
  intraspecific_mode = FALSE
) {
  n <- nrow(df)
  if (is.null(df) || !is.data.frame(df) || n == 0L) {
    return(list(
      edaphic_z_score = numeric(0),
      climatic_z_score = numeric(0),
      env_z_score = numeric(0),
      primary_edaph_vars = character(0),
      primary_clim_vars = character(0),
      primary_env_vars = character(0),
      eval_specimen_z = function(idx, target_tx, mode = "edaphic") 0.0
    ))
  }

  is_intra <- isTRUE(intraspecific_mode) ||
    (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  if (is_intra && !("analysis_group" %in% names(df))) {
    df <- resolve_analysis_group(df, intraspecific_mode = TRUE, primary_env_vars = primary_env_vars, rf_var_weights = rf_var_weights)
  }
  grp_vec <- if (is_intra && "analysis_group" %in% names(df)) as.character(df$analysis_group) else as.character(df$species_standardized)
  active_taxa <- sort(unique(grp_vec[!is.na(grp_vec) & nzchar(grp_vec) & grp_vec != "Unknown"]))

  clim_pat <- "^(bio[0-9]*_|delta_doy$|doy$|hopkins_|sin_doy$|cos_doy$)"
  cand_all <- discover_candidate_eco_predictors(df, include_pheno = TRUE, min_finite = 3L)

  if (is.null(primary_env_vars) || length(primary_env_vars) == 0L) {
    primary_env_vars <- unique(c(as.character(primary_edaph_vars), as.character(primary_clim_vars)))
    if (length(primary_env_vars) == 0L) primary_env_vars <- cand_all
  }
  primary_env_vars <- intersect(unique(as.character(primary_env_vars)), names(df))
  primary_env_vars <- primary_env_vars[vapply(primary_env_vars, function(cl) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    sum(is.finite(v)) >= min(3L, max(2L, n))
  }, logical(1))]

  if (is.null(primary_clim_vars) || length(primary_clim_vars) == 0L) {
    primary_clim_vars <- grep(clim_pat, primary_env_vars, value = TRUE, ignore.case = TRUE)
    if (length(primary_clim_vars) == 0L) {
      primary_clim_vars <- grep(clim_pat, cand_all, value = TRUE, ignore.case = TRUE)
    }
  } else {
    primary_clim_vars <- intersect(unique(as.character(primary_clim_vars)), names(df))
  }

  if (is.null(primary_edaph_vars) || length(primary_edaph_vars) == 0L) {
    primary_edaph_vars <- setdiff(primary_env_vars, grep(clim_pat, primary_env_vars, value = TRUE, ignore.case = TRUE))
    if (length(primary_edaph_vars) == 0L) {
      primary_edaph_vars <- setdiff(cand_all, grep(clim_pat, cand_all, value = TRUE, ignore.case = TRUE))
    }
  } else {
    primary_edaph_vars <- intersect(unique(as.character(primary_edaph_vars)), names(df))
  }

  all_z_vars <- unique(c(primary_env_vars, primary_edaph_vars, primary_clim_vars))
  all_z_vars <- all_z_vars[vapply(all_z_vars, function(cl) {
    (cl %in% names(df)) && sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= min(3L, max(2L, n))
  }, logical(1))]

  auth_w <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else resolve_step06_authority_weights(df)
  auth_w[!is.finite(auth_w) | auth_w <= 0] <- 1.0
  spat_w <- if ("spatial_balance_weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$spatial_balance_weight))
  } else if (exists("compute_spatial_balance_weights", mode = "function")) {
    compute_spatial_balance_weights(df)
  } else {
    rep(1.0, n)
  }
  spat_w[!is.finite(spat_w) | spat_w <= 0] <- 1.0
  row_wts <- auth_w * spat_w

  elig_coord <- if ("is_training_eligible_coord" %in% names(df)) as.logical(df$is_training_eligible_coord) else rep(TRUE, n)
  elig_coord[is.na(elig_coord)] <- TRUE
  det_tier <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n)
  is_prom <- rep(FALSE, n)
  if ("is_promoted_anchor" %in% names(df)) {
    pa <- as.logical(df$is_promoted_anchor); pa[is.na(pa)] <- FALSE; is_prom <- is_prom | pa
  }
  if ("is_promoted_bronze" %in% names(df)) {
    pb <- as.logical(df$is_promoted_bronze); pb[is.na(pb)] <- FALSE; is_prom <- is_prom | pb
  }
  is_anchor <- ((det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_prom) & elig_coord & (grp_vec %in% active_taxa)
  if (sum(is_anchor) < min(6L, n)) is_anchor <- elig_coord & (grp_vec %in% active_taxa)
  if (sum(is_anchor) < min(4L, n)) is_anchor <- grp_vec %in% active_taxa

  w_mean_z <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) return(mean(x[is.finite(x)], na.rm = TRUE))
    sum(x[ok] * w[ok]) / sum(w[ok])
  }
  w_sd_z <- function(x, w, floor_sd) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (sum(ok) < 2L) return(max(stats::sd(x[is.finite(x)], na.rm = TRUE), floor_sd, na.rm = TRUE))
    wm <- sum(x[ok] * w[ok]) / sum(w[ok])
    denom <- sum(w[ok]) - (sum(w[ok]^2) / sum(w[ok]))
    wv <- if (denom > 1e-6) sum(w[ok] * (x[ok] - wm)^2) / denom else stats::var(x[ok])
    max(sqrt(max(wv, 0.0, na.rm = TRUE)), floor_sd, na.rm = TRUE)
  }

  # Data-driven cohort dispersion floor per variable (15% of pooled SD/IQR across cohort, avoiding unit hardcoding)
  scale_floors <- stats::setNames(vapply(all_z_vars, function(cl) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    vf <- v[is.finite(v)]
    if (length(vf) < 2L) return(1e-3)
    sd_g <- stats::sd(vf)
    iqr_g <- stats::IQR(vf) / 1.349
    disp_g <- max(c(sd_g, iqr_g), na.rm = TRUE)
    if (!is.finite(disp_g) || disp_g <= 1e-8) 1e-3 else max(0.15 * disp_g, 1e-4)
  }, numeric(1)), all_z_vars)

  z_profiles <- list()
  if (length(active_taxa) > 0L && length(all_z_vars) > 0L) {
    for (tx in active_taxa) {
      idx_tx <- which(is_anchor & grp_vec == tx)
      if (length(idx_tx) < 3L) idx_tx <- which(elig_coord & grp_vec == tx)
      if (length(idx_tx) < 2L) idx_tx <- which(grp_vec == tx)
      w_tx <- row_wts[idx_tx]
      mu_v <- stats::setNames(numeric(length(all_z_vars)), all_z_vars)
      sd_v <- stats::setNames(numeric(length(all_z_vars)), all_z_vars)
      for (cl in all_z_vars) {
        x_tx <- suppressWarnings(as.numeric(df[[cl]][idx_tx]))
        mu_v[[cl]] <- w_mean_z(x_tx, w_tx)
        sd_v[[cl]] <- w_sd_z(x_tx, w_tx, scale_floors[[cl]])
      }
      z_profiles[[tx]] <- list(mu = mu_v, sd = sd_v)
    }
    for (cl in all_z_vars) {
      bar_sd <- mean(vapply(z_profiles, function(pr) pr$sd[[cl]], numeric(1)), na.rm = TRUE)
      if (is.finite(bar_sd) && bar_sd > 0) {
        for (tx in active_taxa) {
          z_profiles[[tx]]$sd[[cl]] <- max(0.5 * z_profiles[[tx]]$sd[[cl]] + 0.5 * bar_sd, scale_floors[[cl]])
        }
      }
    }
  }

  resolve_subset_weights <- function(vars_sub) {
    if (length(vars_sub) == 0L) return(numeric(0))
    w_out <- stats::setNames(rep(1.0 / length(vars_sub), length(vars_sub)), vars_sub)
    if (!is.null(rf_var_weights) && length(rf_var_weights) > 0L && !is.null(names(rf_var_weights))) {
      hit <- intersect(vars_sub, names(rf_var_weights))
      if (length(hit) > 0L) {
        med_pos <- stats::median(as.numeric(rf_var_weights[rf_var_weights > 0]), na.rm = TRUE)
        if (!is.finite(med_pos) || med_pos <= 0) med_pos <- 1.0 / length(vars_sub)
        w_raw <- stats::setNames(rep(0.25 * med_pos, length(vars_sub)), vars_sub)
        w_hit <- suppressWarnings(as.numeric(rf_var_weights[hit]))
        w_hit[!is.finite(w_hit) | w_hit <= 0] <- 0.25 * med_pos
        w_raw[hit] <- w_hit
        w_out <- w_raw / sum(w_raw)
      }
    }
    w_out
  }

  edaph_use <- if (length(primary_edaph_vars) > 0L) intersect(primary_edaph_vars, all_z_vars) else intersect(primary_env_vars, all_z_vars)
  clim_use  <- intersect(primary_clim_vars, all_z_vars)
  env_use   <- if (length(primary_env_vars) > 0L) intersect(primary_env_vars, all_z_vars) else all_z_vars

  w_edaph <- resolve_subset_weights(edaph_use)
  w_clim  <- resolve_subset_weights(clim_use)
  w_env   <- resolve_subset_weights(env_use)

  eval_specimen_z <- function(idx, target_tx, mode = "edaphic") {
    if (is.na(target_tx) || !(as.character(target_tx) %in% names(z_profiles))) {
      return(if (identical(mode, "climatic")) NA_real_ else 0.0)
    }
    vars_m <- if (identical(mode, "climatic")) clim_use else if (identical(mode, "env")) env_use else edaph_use
    wts_m  <- if (identical(mode, "climatic")) w_clim else if (identical(mode, "env")) w_env else w_edaph
    if (length(vars_m) == 0L) return(if (identical(mode, "climatic")) NA_real_ else 0.0)
    pr <- z_profiles[[as.character(target_tx)]]
    x_row <- vapply(vars_m, function(cl) suppressWarnings(as.numeric(df[[cl]][idx])), numeric(1))
    mu_r <- pr$mu[vars_m]
    sd_r <- pmax(pr$sd[vars_m], scale_floors[vars_m], 1e-4)
    z_sq <- ((x_row - mu_r) / sd_r)^2
    ok <- is.finite(z_sq) & is.finite(wts_m) & (wts_m > 0)
    if (!any(ok)) return(if (identical(mode, "climatic")) NA_real_ else 0.0)
    sqrt(sum(wts_m[ok] * z_sq[ok]) / sum(wts_m[ok]))
  }

  edaph_z_vec <- vapply(seq_len(n), function(i) eval_specimen_z(i, grp_vec[i], "edaphic"), numeric(1))
  clim_z_vec  <- vapply(seq_len(n), function(i) eval_specimen_z(i, grp_vec[i], "climatic"), numeric(1))
  env_z_vec   <- vapply(seq_len(n), function(i) eval_specimen_z(i, grp_vec[i], "env"), numeric(1))

  list(
    edaphic_z_score = round(edaph_z_vec, 4),
    climatic_z_score = round(clim_z_vec, 4),
    env_z_score = round(env_z_vec, 4),
    primary_edaph_vars = edaph_use,
    primary_clim_vars = clim_use,
    primary_env_vars = env_use,
    rf_var_weights = w_env,
    eval_specimen_z = eval_specimen_z
  )
}

compute_type_locality_ecological_control <- function(
  df,
  primary_env_vars = NULL,
  primary_edaph_vars = NULL,
  rf_var_weights = NULL,
  intraspecific_mode = FALSE,
  max_radius_km = 150.0
) {
  n <- if (is.null(df) || !is.data.frame(df)) 0L else nrow(df)
  if (n == 0L) {
    return(list(
      type_edaphic_z_score = numeric(0),
      nearest_type_niche_taxon = character(0),
      matches_type_lithology = logical(0)
    ))
  }

  t_info <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    .identify_type_specimen_vec(df)
  } else {
    ts_raw <- if ("type_status" %in% names(df)) as.character(df$type_status) else rep("", n)
    ts_raw[is.na(ts_raw)] <- ""
    is_t <- grepl("\\b(holotype|lectotype|neotype|epitype|isotype|isolectotype|isoneotype|syntype|isosyntype|paratype|isoparatype|type|cotype)s?\\b", ts_raw, ignore.case = TRUE) &
      !grepl("not a type|non-type|false type", ts_raw, ignore.case = TRUE)
    if ("is_type_specimen" %in% names(df)) {
      b <- as.logical(df$is_type_specimen); b[is.na(b)] <- FALSE
      is_t <- is_t | b
    }
    list(type_rank = ifelse(is_t, 2L, 0L), is_type_specimen = is_t)
  }
  is_type_vec <- as.logical(t_info$is_type_specimen)
  type_rank_vec <- as.integer(t_info$type_rank)

  is_intra <- isTRUE(intraspecific_mode) ||
    (("species_standardized" %in% names(df)) && length(unique(stats::na.omit(df$species_standardized))) < 2L)
  grp_vec <- if (is_intra && "analysis_group" %in% names(df)) {
    as.character(df$analysis_group)
  } else if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("analysis_group" %in% names(df)) {
    as.character(df$analysis_group)
  } else {
    rep("Unknown", n)
  }
  active_taxa <- sort(unique(grp_vec[!is.na(grp_vec) & nzchar(grp_vec) & grp_vec != "Unknown"]))
  if (length(active_taxa) == 0L) active_taxa <- "Unknown"

  cand_all <- discover_candidate_eco_predictors(df, include_pheno = TRUE, min_finite = 2L)
  env_use <- if (!is.null(primary_env_vars) && length(primary_env_vars) > 0L) {
    intersect(as.character(primary_env_vars), names(df))
  } else {
    cand_all
  }
  if (length(env_use) == 0L) {
    env_use <- intersect(c("fine_soil_ph", "fine_soil_sand", "fine_soil_cec", "soil_ph", "soil_sand", "soil_cec", "bio1_temp_mean", "bio12_precip_annual"), names(df))
  }
  clim_pat <- "^(bio[0-9]*_|delta_doy$|doy$|hopkins_|sin_doy$|cos_doy$)"
  edaph_use <- if (!is.null(primary_edaph_vars) && length(primary_edaph_vars) > 0L) {
    intersect(as.character(primary_edaph_vars), names(df))
  } else {
    setdiff(env_use, grep(clim_pat, env_use, value = TRUE, ignore.case = TRUE))
  }
  if (length(edaph_use) == 0L) edaph_use <- env_use
  all_vars <- unique(c(env_use, edaph_use))
  all_vars <- all_vars[vapply(all_vars, function(cl) {
    (cl %in% names(df)) && sum(is.finite(suppressWarnings(as.numeric(df[[cl]])))) >= 1L
  }, logical(1))]

  lat_v <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else if ("latitude" %in% names(df)) suppressWarnings(as.numeric(df$latitude)) else rep(NA_real_, n)
  lon_v <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else if ("longitude" %in% names(df)) suppressWarnings(as.numeric(df$longitude)) else rep(NA_real_, n)
  has_xy <- is.finite(lat_v) & is.finite(lon_v)
  unc_v <- if ("coordinateUncertaintyInMeters" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
  } else if ("coordinateUncertainty" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertainty))
  } else {
    rep(0.0, n)
  }
  unc_v[!is.finite(unc_v) | unc_v < 0.0] <- 0.0
  is_precise_edaph <- unc_v <= 2000.0
  is_precise_clim <- unc_v <= 10000.0
  reg_v <- if ("regional_group" %in% names(df)) trimws(as.character(df$regional_group)) else rep("", n)
  reg_v[is.na(reg_v) | tolower(reg_v) %in% c("na", "unknown", "unknown_ecoregion")] <- ""
  det_t <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n)
  is_gold_v <- grepl("tier_?[01]|gold|type", det_t, ignore.case = TRUE) | is_type_vec
  is_anc_v <- is_gold_v | grepl("tier_?2|silver", det_t, ignore.case = TRUE)
  if ("is_promoted_bronze" %in% names(df)) {
    pb <- as.logical(df$is_promoted_bronze); pb[is.na(pb)] <- FALSE; is_anc_v <- is_anc_v | pb
  }

  haversine_km <- function(lat1, lon1, lat2_vec, lon2_vec) {
    rad <- pi / 180.0
    dlat <- (lat2_vec - lat1) * rad
    dlon <- (lon2_vec - lon1) * rad
    a <- sin(dlat / 2.0)^2 + cos(lat1 * rad) * cos(lat2_vec * rad) * sin(dlon / 2.0)^2
    6371.0 * (2.0 * atan2(sqrt(pmax(0.0, a)), sqrt(pmax(0.0, 1.0 - a))))
  }

  scale_floors <- stats::setNames(vapply(all_vars, function(cl) {
    v <- suppressWarnings(as.numeric(df[[cl]]))
    vf <- v[is.finite(v)]
    if (length(vf) < 2L) return(1.0)
    disp_g <- max(c(stats::mad(vf, na.rm = TRUE), stats::sd(vf, na.rm = TRUE), stats::IQR(vf, na.rm = TRUE) / 1.349), na.rm = TRUE)
    if (!is.finite(disp_g) || disp_g <= 1e-8) 1.0 else max(0.20 * disp_g, 1e-4)
  }, numeric(1)), all_vars)

  lith_col <- if ("sgmc_lithology_class" %in% names(df)) {
    "sgmc_lithology_class"
  } else if ("Dominant_Lithology_Class" %in% names(df)) {
    "Dominant_Lithology_Class"
  } else {
    NULL
  }
  lith_vec <- if (!is.null(lith_col)) trimws(as.character(df[[lith_col]])) else rep("", n)
  lith_vec[is.na(lith_vec) | tolower(lith_vec) %in% c("na", "unknown", "none", "unmapped")] <- ""

  type_profiles <- list()
  for (tx in active_taxa) {
    tx_rows <- which(grp_vec == tx)
    type_rows <- tx_rows[is_type_vec[tx_rows]]
    loc_rows <- integer(0)
    loc_wts <- numeric(0)
    if (length(type_rows) > 0L) {
      near_gold <- integer(0)
      gold_cand <- setdiff(tx_rows[is_gold_v[tx_rows]], type_rows)
      if (length(gold_cand) > 0L) {
        type_regs <- unique(reg_v[type_rows[nzchar(reg_v[type_rows])]])
        gold_cand_eco <- gold_cand[is_precise_edaph[gold_cand]]
        if (length(gold_cand_eco) == 0L) gold_cand_eco <- gold_cand
        near_eco <- if (length(type_regs) > 0L) gold_cand_eco[reg_v[gold_cand_eco] %in% type_regs] else integer(0)
        if (length(near_eco) == 0L && length(type_regs) > 0L) {
          near_eco <- gold_cand[reg_v[gold_cand] %in% type_regs]
        }
        near_xy <- integer(0)
        type_xy <- type_rows[has_xy[type_rows]]
        cand_xy <- gold_cand[has_xy[gold_cand] & is_precise_edaph[gold_cand]]
        if (length(cand_xy) == 0L) cand_xy <- gold_cand[has_xy[gold_cand]]
        if (length(type_xy) > 0L && length(cand_xy) > 0L) {
          min_d <- rep(Inf, length(cand_xy))
          for (ti in type_xy) {
            min_d <- pmin(min_d, haversine_km(lat_v[ti], lon_v[ti], lat_v[cand_xy], lon_v[cand_xy]))
          }
          near_xy <- cand_xy[is.finite(min_d) & (min_d <= max_radius_km)]
        }
        near_gold <- unique(c(near_eco, near_xy))
      }
      loc_rows <- c(type_rows, near_gold)
      loc_wts <- c(ifelse(type_rank_vec[type_rows] == 2L, 3.5, 2.5), rep(1.0, length(near_gold)))
    } else {
      anc_tx <- tx_rows[is_anc_v[tx_rows]]
      loc_rows <- if (length(anc_tx) >= 1L) anc_tx else tx_rows
      loc_wts <- rep(1.0, length(loc_rows))
    }

    mu_tx <- stats::setNames(numeric(length(all_vars)), all_vars)
    sd_tx <- stats::setNames(numeric(length(all_vars)), all_vars)
    for (cl in all_vars) {
      col_v <- suppressWarnings(as.numeric(df[[cl]]))
      is_edaph_col <- cl %in% edaph_use
      type_rows_precise <- if (is_edaph_col) type_rows[is_precise_edaph[type_rows]] else type_rows[is_precise_clim[type_rows]]
      loc_rows_precise <- if (is_edaph_col) loc_rows[is_precise_edaph[loc_rows]] else loc_rows[is_precise_clim[loc_rows]]
      if (length(loc_rows_precise) == 0L || !any(is.finite(col_v[loc_rows_precise]))) {
        loc_rows_precise <- loc_rows
      }
      t_vals <- if (length(type_rows_precise) > 0L) col_v[type_rows_precise] else numeric(0)
      t_vals <- t_vals[is.finite(t_vals)]
      l_vals <- if (length(loc_rows_precise) > 0L) col_v[loc_rows_precise] else numeric(0)
      l_vals <- l_vals[is.finite(l_vals)]
      x_all_tx <- col_v[tx_rows]

      if (length(t_vals) > 0L && length(l_vals) > 0L) {
        mu_tx[[cl]] <- 0.65 * stats::median(t_vals) + 0.35 * stats::median(l_vals)
      } else if (length(t_vals) > 0L) {
        mu_tx[[cl]] <- stats::median(t_vals)
      } else if (length(l_vals) > 0L) {
        mu_tx[[cl]] <- stats::median(l_vals)
      } else {
        t_fallback <- col_v[type_rows]
        t_fallback <- t_fallback[is.finite(t_fallback)]
        if (length(t_fallback) > 0L) {
          mu_tx[[cl]] <- stats::median(t_fallback)
        } else {
          mu_tx[[cl]] <- stats::median(x_all_tx[is.finite(x_all_tx)], na.rm = TRUE)
        }
      }
      if (!is.finite(mu_tx[[cl]])) mu_tx[[cl]] <- 0.0

      anc_tx_rows <- tx_rows[is_anc_v[tx_rows]]
      x_anc_tx <- col_v[anc_tx_rows]
      pool_x <- if (length(l_vals) >= 2L) {
        l_vals
      } else if (sum(is.finite(x_anc_tx)) >= 2L) {
        x_anc_tx[is.finite(x_anc_tx)]
      } else {
        x_all_tx[is.finite(x_all_tx)]
      }
      mad_v <- if (length(pool_x) >= 2L) stats::mad(pool_x, center = mu_tx[[cl]], na.rm = TRUE) else NA_real_
      sd_v <- if (length(pool_x) >= 2L) stats::sd(pool_x, na.rm = TRUE) else NA_real_
      disp_tx <- if (is.finite(mad_v) && mad_v > 1e-6) mad_v else sd_v
      sd_tx[[cl]] <- max(c(disp_tx, scale_floors[[cl]]), na.rm = TRUE)
      if (!is.finite(sd_tx[[cl]]) || sd_tx[[cl]] <= 1e-6) sd_tx[[cl]] <- scale_floors[[cl]]
    }

    type_lith_rows <- type_rows[is_precise_edaph[type_rows]]
    lith_ref <- unique(lith_vec[type_lith_rows[nzchar(lith_vec[type_lith_rows])]])
    if (length(lith_ref) == 0L) {
      loc_lith_precise <- loc_rows[is_precise_edaph[loc_rows]]
      l_cand <- lith_vec[loc_lith_precise[nzchar(lith_vec[loc_lith_precise])]]
      if (length(l_cand) == 0L) {
        l_cand <- lith_vec[loc_rows[nzchar(lith_vec[loc_rows])]]
      }
      if (length(l_cand) == 0L) {
        l_cand <- lith_vec[type_rows[nzchar(lith_vec[type_rows])]]
      }
      if (length(l_cand) > 0L) {
        tb_l <- table(l_cand)
        lith_ref <- names(tb_l[tb_l == max(tb_l)])
      }
    }
    type_profiles[[tx]] <- list(mu = mu_tx, sd = sd_tx, lith_ref = lith_ref)
  }

  get_wts <- function(vars_s) {
    if (length(vars_s) == 0L) return(numeric(0))
    w <- stats::setNames(rep(1.0 / length(vars_s), length(vars_s)), vars_s)
    if (!is.null(rf_var_weights) && length(rf_var_weights) > 0L && !is.null(names(rf_var_weights))) {
      hit <- intersect(vars_s, names(rf_var_weights))
      if (length(hit) > 0L) {
        w_pos <- suppressWarnings(as.numeric(rf_var_weights[hit]))
        w_pos[!is.finite(w_pos) | w_pos <= 0] <- 1e-4
        w[hit] <- w_pos
        w <- w / sum(w)
      }
    }
    w
  }
  edaph_v <- intersect(edaph_use, all_vars)
  env_v <- intersect(env_use, all_vars)
  w_ed <- get_wts(edaph_v)
  w_en <- get_wts(env_v)

  calc_z <- function(idx, tx_name, vars_s, wts_s) {
    if (length(vars_s) == 0L || is.na(tx_name) || !(tx_name %in% names(type_profiles))) return(0.0)
    pr <- type_profiles[[tx_name]]
    xv <- vapply(vars_s, function(cl) suppressWarnings(as.numeric(df[[cl]][idx])), numeric(1))
    z2 <- ((xv - pr$mu[vars_s]) / pmax(pr$sd[vars_s], scale_floors[vars_s], 1e-4))^2
    ok <- is.finite(z2) & is.finite(wts_s) & (wts_s > 0)
    if (!any(ok)) return(0.0)
    sqrt(sum(wts_s[ok] * z2[ok]) / sum(wts_s[ok]))
  }

  type_edaphic_z <- numeric(n)
  nearest_type_niche <- character(n)
  matches_lith <- rep(TRUE, n)

  for (i in seq_len(n)) {
    g_sp <- grp_vec[i]
    if (is.na(g_sp) || !(g_sp %in% active_taxa)) g_sp <- active_taxa[1L]
    z_ed_vec <- vapply(active_taxa, function(tx) calc_z(i, tx, edaph_v, w_ed), numeric(1))
    z_en_vec <- vapply(active_taxa, function(tx) calc_z(i, tx, env_v, w_en), numeric(1))
    comb_z <- 0.60 * z_ed_vec + 0.40 * z_en_vec
    g_pos <- match(g_sp, active_taxa)
    type_edaphic_z[i] <- round(z_ed_vec[if (is.na(g_pos)) 1L else g_pos], 4)
    nearest_type_niche[i] <- active_taxa[which.min(comb_z)][1L]

    pr_g <- type_profiles[[g_sp]]
    if (!is.null(pr_g) && length(pr_g$lith_ref) > 0L && nzchar(lith_vec[i])) {
      matches_lith[i] <- lith_vec[i] %in% pr_g$lith_ref
    } else {
      matches_lith[i] <- type_edaphic_z[i] <= 2.5
    }
  }

  list(
    type_edaphic_z_score = type_edaphic_z,
    nearest_type_niche_taxon = nearest_type_niche,
    matches_type_lithology = as.logical(matches_lith)
  )
}

