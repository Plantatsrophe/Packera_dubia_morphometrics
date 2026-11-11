# ==============================================================================
# Module: dissected_leaf_utils.R - Packera Morphometrics Pipeline (NCU)
# Literature-backed utilities for positional basal-organ homology filtering,
# phase-invariant Dissection Index (DI),
# solidity-stratified EFA harmonic power calibration, multimodal regime-aware
# MCD outlier screening, heteroblastic OTU summary statistics, MorphoTools2
# collinearity pruning (|r| < 0.75), and O(N*P) trace Euclidean PERMANOVA.
#
# Key References:
#   - Kincaid & Schneider (1983) Can. J. Bot. 61:2333-2342; McLellan & Endler
#     (1998) Syst. Biol. 47:264-281; Hearn (2009) Taxon 58:934-954 (DI & EFA)
#   - Claude (2008); Bonhomme et al. (2014) J. Stat. Softw. 56(13) (99% power)
#   - Hubert & Debruyne (2010) WIREs Comput. Stat. 2:36-43 (multimodal MCD)
#   - Chitwood & Sinha (2016); Trock (2006) FNA Vol. 20 (Packera heteroblasty)
#   - Slenker et al. (2022) Appl. Plant Sci. 10:e11495 (MorphoTools2 |r| < 0.75)
# ==============================================================================

compute_polygon_area_perim_pure <- function(coo) {
  if (is.null(coo) || !is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 3L) {
    return(list(area = NA_real_, perim = NA_real_, circularity = NA_real_, solidity = NA_real_))
  }
  x <- as.numeric(coo[, 1L]); y <- as.numeric(coo[, 2L])
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    return(list(area = NA_real_, perim = NA_real_, circularity = NA_real_, solidity = NA_real_))
  }
  xn <- c(x[-1L], x[1L]); yn <- c(y[-1L], y[1L])
  area <- 0.5 * abs(sum(x * yn - xn * y))
  perim <- sum(sqrt((xn - x)^2 + (yn - y)^2))
  circ <- if (is.finite(perim) && perim > 1e-8 && is.finite(area) && area > 0) {
    min(1.0, (4.0 * pi * area) / (perim^2))
  } else {
    NA_real_
  }
  sol <- NA_real_
  if (is.finite(area) && area > 1e-8) {
    ch <- tryCatch(grDevices::chull(x, y), error = function(e) integer(0))
    if (length(ch) >= 3L) {
      hx <- x[ch]; hy <- y[ch]
      ha <- 0.5 * abs(sum(hx * c(hy[-1L], hy[1L]) - c(hx[-1L], hx[1L]) * hy))
      if (is.finite(ha) && ha > 1e-8) sol <- min(1.0, area / ha)
    }
  }
  list(area = area, perim = perim, circularity = circ, solidity = sol)
}

compute_dissection_index_from_coo <- function(coo = NULL, area_val = NULL, perim_val = NULL) {
  a <- if (!is.null(area_val) && length(area_val) == 1L) suppressWarnings(as.numeric(area_val)) else NA_real_
  p <- if (!is.null(perim_val) && length(perim_val) == 1L) suppressWarnings(as.numeric(perim_val)) else NA_real_
  if ((!is.finite(a) || a <= 1e-8 || !is.finite(p) || p <= 1e-8) && !is.null(coo)) {
    gp <- compute_polygon_area_perim_pure(coo)
    a <- gp$area; p <- gp$perim
  }
  if (!is.finite(a) || a <= 1e-8 || !is.finite(p) || p <= 1e-8) return(NA_real_)
  round(max(1.0, p / (2.0 * sqrt(pi * a))), 4)
}

classify_and_partition_leaf_zones <- function(
  fac_df,
  min_cy_plant = 0.40,
  min_cy_sheet = 0.42,
  fallback_band = 0.15
) {
  if (is.null(fac_df) || nrow(fac_df) == 0L) return(fac_df)
  n <- nrow(fac_df)
  if (!("plant_individual_id" %in% names(fac_df))) fac_df$plant_individual_id <- 0L
  fac_df$plant_individual_id <- ifelse(
    is.na(fac_df$plant_individual_id),
    0L,
    suppressWarnings(as.integer(fac_df$plant_individual_id))
  )
  fac_df$plant_individual_id[is.na(fac_df$plant_individual_id)] <- 0L

  raw_zone <- if ("leaf_zone" %in% names(fac_df)) tolower(trimws(as.character(fac_df$leaf_zone))) else rep(NA_character_, n)
  raw_zone[!(raw_zone %in% c("basal", "cauline"))] <- NA_character_
  zone_conf <- if ("leaf_zone_confidence" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$leaf_zone_confidence)) else rep(NA_real_, n)
  d_crown <- if ("dist_from_crown_norm" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$dist_from_crown_norm)) else rep(NA_real_, n)
  cy_p <- if ("cy_rel_plant_bbox" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$cy_rel_plant_bbox)) else rep(NA_real_, n)
  cy_s <- if ("cy_rel_sheet" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$cy_rel_sheet)) else rep(NA_real_, n)
  pbr <- if ("petiole_to_blade_ratio" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$petiole_to_blade_ratio)) else rep(NA_real_, n)
  area_px <- if ("area_px" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$area_px)) else rep(NA_real_, n)
  area_mm2 <- if ("area_mm2" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$area_mm2)) else rep(NA_real_, n)

  leaf_zone <- raw_zone
  is_fallback_vec <- rep(FALSE, n)
  is_caul_only_vec <- rep(FALSE, n)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)

  for (g_idx in idx_by_grp) {
    lam_ok <- g_idx
    if (length(lam_ok) == 0L) next
    all_explicit_cauline <- all(!is.na(raw_zone[lam_ok]) & raw_zone[lam_ok] == "cauline")

    use_px_pool <- any(is.finite(area_px[lam_ok]) & area_px[lam_ok] > 0.0)
    a_pool <- if (use_px_pool) {
      ifelse(is.finite(area_px[lam_ok]) & area_px[lam_ok] > 0.0, area_px[lam_ok], NA_real_)
    } else {
      ifelse(is.finite(area_mm2[lam_ok]) & area_mm2[lam_ok] > 0.0, area_mm2[lam_ok], NA_real_)
    }
    max_pid_a <- if (any(is.finite(a_pool))) max(a_pool, na.rm = TRUE) else 1.0

    for (ii in g_idx) {
      if (!is.na(leaf_zone[ii])) {
        if (!is.finite(zone_conf[ii])) zone_conf[ii] <- 0.90
        next
      }
      pos_ii <- match(ii, lam_ok)
      rel_a <- if (!is.na(pos_ii) && is.finite(a_pool[pos_ii]) && max_pid_a > 0.0) a_pool[pos_ii] / max_pid_a else 1.0
      if (is.finite(d_crown[ii])) {
        s_b <- if (d_crown[ii] <= 0.42) 0.76 else if (d_crown[ii] > 0.45) 0.22 else 0.50
        if (is.finite(cy_p[ii])) s_b <- s_b + (if (cy_p[ii] >= min_cy_plant) 0.12 else -0.12)
        if (is.finite(pbr[ii])) s_b <- s_b + (if (pbr[ii] >= 0.12) 0.14 else if (pbr[ii] < 0.08 && d_crown[ii] > 0.38) -0.14 else 0.0)
        if (is.finite(rel_a)) s_b <- s_b + (if (rel_a < 0.45 && d_crown[ii] > 0.38) -0.12 else if (rel_a >= 0.65) 0.08 else 0.0)
        p_b <- min(1.0, max(0.0, s_b))
        leaf_zone[ii] <- if (p_b >= 0.50) "basal" else "cauline"
        zone_conf[ii] <- round(max(0.55, abs(p_b - 0.50) * 2.0), 4)
      } else if (is.finite(cy_p[ii]) || is.finite(cy_s[ii])) {
        in_prim <- if (is.finite(cy_p[ii])) (cy_p[ii] >= min_cy_plant) else (cy_s[ii] >= min_cy_sheet)
        leaf_zone[ii] <- if (isTRUE(in_prim)) "basal" else "cauline"
        zone_conf[ii] <- if (is.finite(zone_conf[ii])) zone_conf[ii] else 0.80
      } else {
        is_caul_heur <- is.finite(pbr[ii]) && (pbr[ii] < 0.06) && is.finite(rel_a) && (rel_a < 0.40) && (length(lam_ok) >= 2L)
        leaf_zone[ii] <- if (is_caul_heur) "cauline" else "basal"
        zone_conf[ii] <- if (is.finite(zone_conf[ii])) zone_conf[ii] else 0.75
      }
    }

    has_basal <- any(leaf_zone[lam_ok] == "basal", na.rm = TRUE)
    if (!has_basal) {
      if (all_explicit_cauline) {
        is_caul_only_vec[lam_ok] <- TRUE
        is_fallback_vec[lam_ok] <- TRUE
      } else {
        has_cy_p <- is.finite(cy_p[lam_ok])
        has_cy_s <- !has_cy_p & is.finite(cy_s[lam_ok])
        has_dc <- is.finite(d_crown[lam_ok])
        if (any(has_cy_p)) {
          c_max <- max(cy_p[lam_ok[has_cy_p]], na.rm = TRUE)
          fb_hit <- has_cy_p & (cy_p[lam_ok] >= (c_max - fallback_band))
        } else if (any(has_cy_s)) {
          c_max <- max(cy_s[lam_ok[has_cy_s]], na.rm = TRUE)
          fb_hit <- has_cy_s & (cy_s[lam_ok] >= (c_max - fallback_band))
        } else if (any(has_dc)) {
          d_min <- min(d_crown[lam_ok[has_dc]], na.rm = TRUE)
          fb_hit <- has_dc & (d_crown[lam_ok] <= (d_min + fallback_band))
        } else {
          fb_hit <- seq_along(lam_ok) == which.max(ifelse(is.finite(a_pool), a_pool, 0.0))
        }
        if (any(fb_hit)) {
          leaf_zone[lam_ok[fb_hit]] <- "basal"
          is_fallback_vec[lam_ok[fb_hit]] <- TRUE
        }
      }
    }
  }

  fac_df$leaf_zone <- leaf_zone
  fac_df$leaf_zone_confidence <- zone_conf
  fac_df$dist_from_crown_norm <- d_crown
  fac_df$is_basal_fallback <- is_fallback_vec
  fac_df$is_cauline_only_individual <- is_caul_only_vec
  fac_df
}

filter_basal_lamina_contours <- function(
  fac_df,
  coo_list,
  min_cy_plant = 0.40,
  min_cy_sheet = 0.42,
  fallback_band = 0.15,
  retain_cauline_tagged = FALSE
) {
  if (is.null(fac_df) || nrow(fac_df) == 0L || length(coo_list) == 0L) {
    return(list(fac_df = fac_df, coo_list = coo_list, stats = list()))
  }
  n <- nrow(fac_df)
  di_vec <- if ("dissection_index" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$dissection_index)) else rep(NA_real_, n)
  sol_vec <- if ("solidity" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$solidity)) else rep(NA_real_, n)
  circ_vec <- if ("circularity" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$circularity)) else rep(NA_real_, n)

  miss_geom <- which(!is.finite(di_vec) | !is.finite(sol_vec) | !is.finite(circ_vec))
  if (length(miss_geom) > 0L) {
    for (i in miss_geom) {
      gp <- compute_polygon_area_perim_pure(coo_list[[i]])
      if (!is.finite(sol_vec[i])) sol_vec[i] <- round(gp$solidity, 4)
      if (!is.finite(circ_vec[i])) circ_vec[i] <- round(gp$circularity, 4)
      if (!is.finite(di_vec[i])) di_vec[i] <- compute_dissection_index_from_coo(area_val = gp$area, perim_val = gp$perim)
    }
  }

  fac_df$solidity <- sol_vec
  fac_df$circularity <- circ_vec
  fac_df$dissection_index <- di_vec
  raw_dis <- if ("is_dissected" %in% names(fac_df)) as.logical(fac_df$is_dissected) else rep(NA, n)
  fac_df$is_dissected <- ifelse(
    !is.na(raw_dis),
    raw_dis | (is.finite(sol_vec) & sol_vec < 0.72),
    (is.finite(sol_vec) & sol_vec < 0.72)
  )

  fac_df <- classify_and_partition_leaf_zones(
    fac_df,
    min_cy_plant = min_cy_plant,
    min_cy_sheet = min_cy_sheet,
    fallback_band = fallback_band
  )
  is_fallback_vec <- as.logical(fac_df$is_basal_fallback)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)
  n_fb_individuals <- sum(vapply(idx_by_grp, function(g_idx) any(is_fallback_vec[g_idx] & !fac_df$is_cauline_only_individual[g_idx]), logical(1)))

  basal_keep <- (fac_df$leaf_zone == "basal") | is_fallback_vec
  for (g_idx in idx_by_grp) {
    if (!any(basal_keep[g_idx]) && length(g_idx) > 0L) {
      basal_keep[g_idx[1L]] <- TRUE
    }
  }
  if (!any(basal_keep)) basal_keep <- rep(TRUE, n)

  n_area_drop <- 0L
  n_caul_drop <- sum(!basal_keep)
  keep_mask <- if (isTRUE(retain_cauline_tagged)) rep(TRUE, n) else basal_keep

  if (n_caul_drop > 0L || n_fb_individuals > 0L) {
    message(sprintf(
      paste0(
        "Positional Homology & Basal Lamina Gate: classified %d / %d primary basal blades ",
        "(cauline-partitioned = %d; lowest-lamina fallback individuals = %d)."
      ),
      sum(basal_keep), n, n_caul_drop, n_fb_individuals
    ))
  }
  fac_sub <- fac_df[keep_mask, , drop = FALSE]
  coo_sub <- coo_list[keep_mask]
  rownames(fac_sub) <- names(coo_sub)
  list(
    fac_df = fac_sub,
    coo_list = coo_sub,
    stats = list(
      n_input = n,
      n_retained = sum(keep_mask),
      n_basal_retained = sum(basal_keep),
      n_cauline_excluded = n_caul_drop,
      n_subfloor_excluded = n_area_drop,
      n_fallback_individuals = n_fb_individuals
    )
  )
}

.extract_harmonic_opt_from_cal <- function(cal_res, min_h = 6L, max_h = 24L) {
  if (is.null(cal_res) || is.null(cal_res$q) || length(cal_res$q) == 0L) return(NULL)
  q_mat <- as.matrix(cal_res$q)
  q_vec <- if ("50%" %in% rownames(q_mat)) as.numeric(q_mat["50%", ]) else apply(q_mat, 2, stats::median, na.rm = TRUE)
  if (!any(is.finite(q_vec))) return(NULL)
  if (max(q_vec, na.rm = TRUE) <= 1.05) q_vec <- q_vec * 100.0
  h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(q_mat))))
  if (length(h_nums) != length(q_vec) || any(is.na(h_nums))) h_nums <- seq_along(q_vec)
  hit_idx <- which(is.finite(q_vec) & q_vec >= 99.0)
  pick_i <- if (length(hit_idx) > 0L) hit_idx[1L] else which.max(h_nums)
  h_val <- max(as.integer(min_h), min(as.integer(max_h), as.integer(h_nums[pick_i])))
  col_m <- match(h_val, h_nums)
  cum_p <- as.numeric(q_vec[if (!is.na(col_m)) col_m else pick_i])
  list(h_opt = as.integer(h_val), cum_pct = cum_p)
}

calibrate_dissection_aware_harmonics <- function(
  out_obj = NULL,
  harm_mat = NULL,
  solidity_vec = NULL,
  nb_h = 12L,
  min_h = 6L,
  max_h = 24L
) {
  fb_h <- max(as.integer(min_h), min(as.integer(max_h), as.integer(if (is.null(nb_h) || is.na(nb_h)) 12L else nb_h)))
  h_opt <- fb_h
  cum_pct <- NA_real_

  if (!is.null(out_obj) && requireNamespace("Momocs", quietly = TRUE) && !is.null(out_obj$coo)) {
    n_tot <- length(out_obj$coo)
    sol <- if (!is.null(solidity_vec) && length(solidity_vec) == n_tot) {
      suppressWarnings(as.numeric(solidity_vec))
    } else if (!is.null(out_obj$fac) && "solidity" %in% names(out_obj$fac)) {
      suppressWarnings(as.numeric(out_obj$fac$solidity))
    } else {
      rep(NA_real_, n_tot)
    }

    strat_idx <- if (n_tot > 240L && sum(is.finite(sol)) >= 20L) {
      qs <- stats::quantile(sol[is.finite(sol)], probs = c(0, 0.25, 0.50, 0.75, 1.0), na.rm = TRUE, names = FALSE)
      bins <- findInterval(ifelse(is.finite(sol), sol, stats::median(sol, na.rm = TRUE)), unique(qs), rightmost.closed = TRUE)
      unlist(lapply(sort(unique(bins)), function(b) {
        ib <- which(bins == b)
        if (length(ib) <= 60L) ib else ib[round(seq(1, length(ib), length.out = 60L))]
      }), use.names = FALSE)
    } else if (n_tot > 240L) {
      round(seq(1, n_tot, length.out = 240L))
    } else {
      seq_len(n_tot)
    }
    strat_idx <- sort(unique(pmax(1L, pmin(n_tot, as.integer(strat_idx)))))

    cal_all <- tryCatch({
      sub_all <- if (length(strat_idx) < n_tot) Momocs::slice(out_obj, strat_idx) else out_obj
      Momocs::calibrate_harmonicpower_efourier(sub_all, nb.h = as.integer(max_h), plot = FALSE)
    }, error = function(e) NULL)
    parsed_all <- .extract_harmonic_opt_from_cal(cal_all, min_h = min_h, max_h = max_h)
    if (!is.null(parsed_all)) {
      h_opt <- parsed_all$h_opt
      cum_pct <- parsed_all$cum_pct
    }

    # Ensure dissected / lowest-solidity quartile also achieves >= 99% cumulative harmonic power
    dis_idx <- which(is.finite(sol) & sol < 0.72)
    if (length(dis_idx) >= 5L) {
      if (length(dis_idx) > 80L) {
        ord_d <- dis_idx[order(sol[dis_idx])]
        dis_idx <- sort(ord_d[round(seq(1, length(ord_d), length.out = 80L))])
      }
      cal_dis <- tryCatch({
        sub_dis <- Momocs::slice(out_obj, dis_idx)
        Momocs::calibrate_harmonicpower_efourier(sub_dis, nb.h = as.integer(max_h), plot = FALSE)
      }, error = function(e) NULL)
      parsed_dis <- .extract_harmonic_opt_from_cal(cal_dis, min_h = min_h, max_h = max_h)
      if (!is.null(parsed_dis) && parsed_dis$h_opt > h_opt) {
        h_opt <- parsed_dis$h_opt
        cum_pct <- min(c(cum_pct, parsed_dis$cum_pct), na.rm = TRUE)
      }
    }
  }

  list(h_opt = as.integer(h_opt), cum_pct = as.numeric(cum_pct))
}

compute_dissection_aware_roughness <- function(coe, solidity_vec = NULL, pct_asymmetry_vec = NULL) {
  coe <- as.matrix(coe)
  n_rows <- nrow(coe)
  if (n_rows == 0L || ncol(coe) == 0L) {
    return(list(harmonic_roughness = numeric(0), is_contour_outlier = logical(0)))
  }
  all_h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(coe))))
  max_h <- if (any(is.finite(all_h_nums))) max(all_h_nums, na.rm = TRUE) else 12L
  # Measure high-frequency segmentation jaggedness above the biological pinna-lobe band (n = 4..8)
  min_high_n <- max(9L, min(14L, max_h - 3L))
  if (max_h < 9L) min_high_n <- max(6L, max_h - 1L)

  low_h_cols <- colnames(coe)[!is.na(all_h_nums) & all_h_nums >= 1L & all_h_nums <= 3L]
  high_h_cols <- colnames(coe)[!is.na(all_h_nums) & all_h_nums >= min_high_n]
  high_n <- as.numeric(gsub("[^0-9]", "", high_h_cols))

  low_power <- if (length(low_h_cols) > 0L) rowSums(coe[, low_h_cols, drop = FALSE]^2, na.rm = TRUE) else rep(1.0, n_rows)
  high_power <- if (length(high_h_cols) > 0L) {
    rowSums(sweep(coe[, high_h_cols, drop = FALSE]^2, 2L, (high_n / 2.2)^2, "*"), na.rm = TRUE)
  } else {
    rep(0.0, n_rows)
  }
  harmonic_roughness <- high_power / pmax(low_power, 1e-12)

  pct_asym <- if (!is.null(pct_asymmetry_vec) && length(pct_asymmetry_vec) == n_rows) {
    as.numeric(pct_asymmetry_vec)
  } else {
    rep(0.0, n_rows)
  }
  sol <- if (!is.null(solidity_vec) && length(solidity_vec) == n_rows) {
    suppressWarnings(as.numeric(solidity_vec))
  } else {
    rep(1.0, n_rows)
  }
  sol[!is.finite(sol)] <- 1.0

  # Bilaterally symmetric deeply dissected leaves naturally carry finite energy in n >= 9;
  # scale roughness ceiling smoothly up to 0.28 only when bilateral symmetry passes (pct_asymmetry <= 0.20)
  rough_ceil <- ifelse(
    sol < 0.72 & is.finite(pct_asym) & pct_asym <= 0.20,
    pmin(0.28, 0.15 + 0.24 * (0.72 - sol)),
    0.15
  )
  is_out <- if (n_rows >= 10L) (pct_asym > 0.22) | (harmonic_roughness > rough_ceil) else rep(FALSE, n_rows)
  is_out[is.na(is_out)] <- FALSE
  if (all(is_out)) is_out <- rep(FALSE, n_rows)
  list(harmonic_roughness = harmonic_roughness, is_contour_outlier = is_out)
}

screen_regime_stratified_mcd_outliers <- function(
  pc_scores,
  fit_pos = NULL,
  solidity_vec = NULL,
  is_dissected_vec = NULL,
  k_mcd = 5L,
  p_cutoff = 0.001,
  seed = 42L
) {
  pc_scores <- as.matrix(pc_scores)
  n <- nrow(pc_scores)
  mcd_out <- rep(FALSE, n)
  if (n < 15L || !requireNamespace("MASS", quietly = TRUE)) return(mcd_out)
  k_eff <- max(2L, min(as.integer(k_mcd), ncol(pc_scores)))
  pc_sub <- pc_scores[, seq_len(k_eff), drop = FALSE]

  fit_mask <- rep(FALSE, n)
  if (!is.null(fit_pos) && length(fit_pos) > 0L) {
    valid_fp <- fit_pos[!is.na(fit_pos) & fit_pos >= 1L & fit_pos <= n]
    fit_mask[valid_fp] <- TRUE
  } else {
    fit_mask[] <- TRUE
  }
  if (sum(fit_mask) < 20L) return(mcd_out)

  sol <- if (!is.null(solidity_vec) && length(solidity_vec) == n) suppressWarnings(as.numeric(solidity_vec)) else rep(0.90, n)
  sol[!is.finite(sol)] <- 0.90
  is_dis <- if (!is.null(is_dissected_vec) && length(is_dissected_vec) == n) {
    (!is.na(is_dissected_vec) & as.logical(is_dissected_vec)) | (sol < 0.72)
  } else {
    sol < 0.72
  }

  # Partition into morphological regimes so entire-leaf dominance (>50%) never collapses dissected modes (Hubert & Debruyne 2010)
  regimes <- ifelse(!is_dis, "Entire", "Dissected")
  if (sum(is_dis & (sol < 0.50)) >= 60L && sum(is_dis & (sol >= 0.50)) >= 60L) {
    regimes[is_dis & (sol < 0.50)] <- "Deep_Pinnatifid"
    regimes[is_dis & (sol >= 0.50)] <- "Moderate_Lobed"
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)

  for (reg_nm in unique(regimes)) {
    reg_idx <- which(regimes == reg_nm)
    reg_fit <- which((regimes == reg_nm) & fit_mask)
    if (length(reg_fit) < max(12L, k_eff + 5L)) reg_fit <- reg_idx
    k_reg <- min(k_eff, max(2L, floor((length(reg_fit) - 4L) / 2L)))
    if (length(reg_fit) < (k_reg + 5L) || length(reg_idx) < (k_reg + 5L)) next
    X_reg_all <- pc_sub[reg_idx, seq_len(k_reg), drop = FALSE]
    X_reg_fit <- pc_sub[reg_fit, seq_len(k_reg), drop = FALSE]
    v_ok <- apply(X_reg_fit, 2L, stats::var, na.rm = TRUE) > 1e-10
    if (sum(v_ok) < 2L) next
    X_reg_all <- X_reg_all[, v_ok, drop = FALSE]
    X_reg_fit <- X_reg_fit[, v_ok, drop = FALSE]
    df_reg <- ncol(X_reg_fit)
    var_floor <- pmax(1e-4, 0.05 * apply(X_reg_fit, 2L, stats::var, na.rm = TRUE))
    var_floor[!is.finite(var_floor)] <- 1e-4
    if (nrow(X_reg_fit) >= 30L) {
      q_used <- min(nrow(X_reg_fit) - 1L, max(floor(0.88 * nrow(X_reg_fit)), floor((nrow(X_reg_fit) + df_reg + 1L) / 2L)))
      if (nrow(X_reg_fit) > 1000L && exists(".fit_fast_mcd_cov", mode = "function")) {
        mcd_fit <- tryCatch(.fit_fast_mcd_cov(X_reg_fit, q_used = q_used, seed = as.integer(seed)), error = function(e) NULL)
      } else if (nrow(X_reg_fit) > 1000L) {
        set.seed(as.integer(seed))
        ord_r <- order(X_reg_fit[, 1L])
        sub_r <- sort(ord_r[unique(round(seq(1, nrow(X_reg_fit), length.out = 800L)))])
        q_sub <- min(length(sub_r) - 1L, max(floor((q_used / nrow(X_reg_fit)) * length(sub_r)), floor((length(sub_r) + df_reg + 1L) / 2L)))
        sub_fit <- tryCatch(MASS::cov.mcd(X_reg_fit[sub_r, , drop = FALSE], quantile.used = q_sub), error = function(e) NULL)
        if (!is.null(sub_fit) && !is.null(sub_fit$cov)) {
          ctr_c <- sub_fit$center
          cov_c <- sub_fit$cov + diag(1e-7, df_reg)
          for (c_it in 1:2) {
            d2_c <- tryCatch(stats::mahalanobis(X_reg_fit, center = ctr_c, cov = cov_c), error = function(e) NULL)
            if (is.null(d2_c) || any(!is.finite(d2_c))) break
            h_idx <- order(d2_c)[seq_len(q_used)]
            X_h <- X_reg_fit[h_idx, , drop = FALSE]
            ctr_c <- colMeans(X_h)
            cov_c <- stats::cov(X_h) + diag(1e-7, df_reg)
          }
          mcd_fit <- list(center = ctr_c, cov = cov_c)
        } else {
          mcd_fit <- NULL
        }
      } else {
        set.seed(as.integer(seed))
        mcd_fit <- tryCatch(MASS::cov.mcd(X_reg_fit, quantile.used = q_used), error = function(e) NULL)
      }
    } else {
      mcd_fit <- NULL
      q_used <- nrow(X_reg_fit)
    }
    if (!is.null(mcd_fit) && !is.null(mcd_fit$cov)) {
      alpha_q <- min(0.999, max(0.50, q_used / nrow(X_reg_fit)))
      denom_c <- stats::pgamma(stats::qchisq(alpha_q, df = df_reg) / 2.0, shape = df_reg / 2.0 + 1.0)
      c_fac <- if (is.finite(denom_c) && denom_c > 1e-8) (alpha_q / denom_c) else 1.0
      m_center <- mcd_fit$center
      cov_reg <- (c_fac * mcd_fit$cov) + diag(var_floor, nrow = df_reg, ncol = df_reg)
    } else {
      m_center <- colMeans(X_reg_fit, na.rm = TRUE)
      cov_reg <- stats::cov(X_reg_fit) + diag(var_floor, nrow = df_reg, ncol = df_reg)
    }
    d2 <- tryCatch(
      stats::mahalanobis(X_reg_all, center = m_center, cov = cov_reg),
      error = function(e) rep(0.0, nrow(X_reg_all))
    )
    thr <- stats::qchisq(1.0 - p_cutoff, df = df_reg)
    out_r <- is.finite(d2) & (d2 > thr)
    if (any(out_r) && sum(!out_r) >= 5L) {
      mcd_out[reg_idx[out_r]] <- TRUE
    }
  }
  if (sum(!mcd_out) < 5L) mcd_out[] <- FALSE
  mcd_out
}

summarize_otu_dissection_metrics <- function(grp_df) {
  if ("leaf_zone" %in% names(grp_df)) {
    lz <- tolower(trimws(as.character(grp_df$leaf_zone)))
    if (any(!is.na(lz) & lz == "basal")) {
      grp_df <- grp_df[!is.na(lz) & lz == "basal", , drop = FALSE]
    }
  }
  n_l <- nrow(grp_df)
  med_num <- function(col_nm, digits = 4L) {
    if (!(col_nm %in% names(grp_df))) return(NA_real_)
    v <- suppressWarnings(as.numeric(grp_df[[col_nm]]))
    v <- v[is.finite(v)]
    if (length(v) == 0L) NA_real_ else round(stats::median(v), digits)
  }
  sol_v <- if ("solidity" %in% names(grp_df)) suppressWarnings(as.numeric(grp_df$solidity)) else numeric(0)
  sol_fin <- sol_v[is.finite(sol_v)]
  q25_sol <- if (length(sol_fin) > 0L) {
    round(as.numeric(stats::quantile(sol_fin, probs = 0.25, na.rm = TRUE, names = FALSE)), 4)
  } else {
    NA_real_
  }
  dis_raw <- if ("is_dissected" %in% names(grp_df)) as.logical(grp_df$is_dissected) else rep(NA, n_l)
  dis_leaf <- ifelse(!is.na(dis_raw), dis_raw, FALSE) | (if (length(sol_v) == n_l) (is.finite(sol_v) & sol_v < 0.72) else FALSE)
  dis_frac <- if (n_l > 0L && (any(!is.na(dis_raw)) || length(sol_fin) > 0L)) round(mean(dis_leaf, na.rm = TRUE), 4) else NA_real_

  list(
    basal_leaf_count = as.integer(n_l),
    solidity = med_num("solidity", 4L),
    q25_solidity = q25_sol,
    dissection_index = med_num("dissection_index", 4L),
    dissected_leaf_fraction = dis_frac,
    aspect_ratio = med_num("aspect_ratio", 4L),
    circularity = med_num("circularity", 4L),
    area_px = med_num("area_px", 1L),
    sym_var = med_num("sym_var", 6L),
    asym_var = med_num("asym_var", 6L),
    pct_asymmetry = med_num("pct_asymmetry", 6L),
    harmonic_roughness = med_num("harmonic_roughness", 6L),
    leaf_quality_weight = med_num("leaf_quality_weight", 6L)
  )
}

read_numpy_ect_r <- function(npy_path, num_thresh = NULL, num_dirs = NULL) {
  if (is.null(npy_path) || is.na(npy_path) || !nzchar(trimws(as.character(npy_path))) || !file.exists(npy_path)) {
    return(NULL)
  }
  fsize <- file.info(npy_path)$size
  if (is.na(fsize) || fsize < 80L) return(NULL)
  con <- file(npy_path, "rb")
  on.exit(close(con), add = TRUE)
  # Support full 32x360 float32/float64 LM3 ECT matrices (up to 512 KB)
  raw_b <- readBin(con, what = "raw", n = min(as.integer(fsize), 524288L))
  if (length(raw_b) < 16L || as.integer(raw_b[1L]) != 147L) return(NULL)
  if (rawToChar(raw_b[2:6]) != "NUMPY") return(NULL)
  ver_maj <- as.integer(raw_b[7L])
  if (ver_maj == 1L) {
    hlen <- as.integer(raw_b[9L]) + 256L * as.integer(raw_b[10L])
    data_start <- 10L + hlen
  } else {
    hlen <- as.integer(raw_b[9L]) + 256L * as.integer(raw_b[10L]) + 65536L * as.integer(raw_b[11L]) + 16777216L * as.integer(raw_b[12L])
    data_start <- 12L + hlen
  }
  if (data_start >= length(raw_b)) return(NULL)
  hdr_str <- rawToChar(raw_b[seq.int(if (ver_maj == 1L) 11L else 13L, data_start)])
  is_f8 <- grepl("<f8|float64", hdr_str)
  is_fortran <- grepl("'fortran_order':\\s*True", hdr_str)
  elem_sz <- if (is_f8) 8L else 4L
  nr_use <- if (!is.null(num_thresh) && is.finite(num_thresh)) as.integer(num_thresh) else 32L
  nc_use <- if (!is.null(num_dirs) && is.finite(num_dirs)) as.integer(num_dirs) else 360L
  sh_m <- regmatches(hdr_str, regexec("'shape'\\s*:\\s*\\(\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\)", hdr_str))[[1L]]
  if (length(sh_m) >= 3L) {
    nr_hdr <- suppressWarnings(as.integer(sh_m[2L]))
    nc_hdr <- suppressWarnings(as.integer(sh_m[3L]))
    if (is.finite(nr_hdr) && nr_hdr >= 4L && is.finite(nc_hdr) && nc_hdr >= 4L) {
      nr_use <- nr_hdr
      nc_use <- nc_hdr
    }
  }
  n_need <- as.integer(nr_use) * as.integer(nc_use)
  payload <- raw_b[seq.int(data_start + 1L, length(raw_b))]
  if (length(payload) < n_need * elem_sz) return(NULL)
  vals <- readBin(payload, what = "numeric", n = n_need, size = elem_sz, endian = "little")
  if (length(vals) != n_need || any(!is.finite(vals))) return(NULL)
  matrix(vals, nrow = as.integer(nr_use), ncol = as.integer(nc_use), byrow = !is_fortran)
}

.resample_ect_directions_r <- function(ect_mat, target_nr = 32L, target_nc = 360L) {
  if (is.null(ect_mat) || !is.matrix(ect_mat)) return(NULL)
  if (nrow(ect_mat) == target_nr && ncol(ect_mat) == target_nc) return(ect_mat)
  # Periodic linear interpolation along phi in [0, 2*pi) preserving Column 1 = phi = 0 (lamina_tip)
  if (nrow(ect_mat) != target_nr) {
    x_old <- seq(0.0, 1.0, length.out = nrow(ect_mat))
    x_new <- seq(0.0, 1.0, length.out = target_nr)
    ect_r <- matrix(0.0, nrow = target_nr, ncol = ncol(ect_mat))
    for (j in seq_len(ncol(ect_mat))) {
      ect_r[, j] <- stats::approx(x_old, ect_mat[, j], xout = x_new, rule = 2)$y
    }
    ect_mat <- ect_r
  }
  if (ncol(ect_mat) != target_nc) {
    nc_old <- ncol(ect_mat)
    phi_ext <- seq(0.0, 1.0, length.out = nc_old + 1L)
    phi_new <- seq(0.0, 1.0 - 1.0 / target_nc, length.out = target_nc)
    ect_c <- matrix(0.0, nrow = target_nr, ncol = target_nc)
    for (i in seq_len(target_nr)) {
      y_ext <- c(ect_mat[i, ], ect_mat[i, 1L])
      ect_c[i, ] <- stats::approx(phi_ext, y_ext, xout = phi_new, rule = 2)$y
    }
    ect_mat <- ect_c
  }
  ect_mat
}

compute_cycle_ect_r <- function(coo, num_thresh = 32L, num_dirs = 360L) {
  if (is.null(coo) || !is.matrix(coo) || nrow(coo) < 6L || ncol(coo) < 2L) return(NULL)
  xy <- matrix(as.numeric(coo[, 1:2]), ncol = 2L)
  if (any(!is.finite(xy))) return(NULL)
  # Orient tip-up: when coo is already horizontally aligned by Step 03 (+X = lamina_tip or -X = lamina_base),
  # use the exact horizontal midrib axis; otherwise project from vertex 1 (lamina_base) to apical midline cap
  if (isTRUE(attr(coo, "lm3_xax_aligned"))) {
    cx0 <- mean(xy[, 1L]); p_base <- xy[1L, ]
    u_tip <- if (p_base[1L] <= cx0) c(1.0, 0.0) else c(-1.0, 0.0)
    v_len <- 1.0
    v_axis <- u_tip
  } else {
    p_base <- xy[1L, ]
    d_base <- sqrt((xy[, 1L] - p_base[1L])^2 + (xy[, 2L] - p_base[2L])^2)
    cap_i <- which(d_base >= 0.985 * max(d_base, na.rm = TRUE))
    v_axis <- c(stats::median(xy[cap_i, 1L]) - p_base[1L], stats::median(xy[cap_i, 2L]) - p_base[2L])
  }
  v_len <- sqrt(sum(v_axis^2))
  if (!is.finite(v_len) || v_len <= 1e-8) {
    u_tip <- c(1.0, 0.0)
  } else {
    u_tip <- v_axis / v_len
  }
  # In Momocs Cartesian space (+X = lamina_tip, +Y = anatomical left, -Y = anatomical right),
  # 90-deg clockwise from u_tip = (ux, uy) is u_right = (uy, -ux), matching Step 02 tip_up_directions().
  u_right <- c(u_tip[2L], -u_tip[1L])
  cx <- mean(xy[, 1L]); cy <- mean(xy[, 2L])
  xy_c <- cbind(xy[, 1L] - cx, xy[, 2L] - cy)
  radii <- sqrt(xy_c[, 1L]^2 + xy_c[, 2L]^2)
  r_max <- max(radii, na.rm = TRUE)
  if (!is.finite(r_max) || r_max <= 1e-8) return(NULL)
  xy_n <- xy_c / r_max
  # Sweep phi in [0, 2*pi) starting at phi = 0 (lamina_tip) and ascending clockwise toward u_right
  phis <- seq(0.0, 2.0 * pi * (1.0 - 1.0 / num_dirs), length.out = num_dirs)
  thresh <- seq(-1.0, 1.0, length.out = num_thresh)
  n_v <- nrow(xy_n)
  next_idx <- c(seq.int(2L, n_v), 1L)
  ect_mat <- matrix(0.0, nrow = num_thresh, ncol = num_dirs)
  for (j in seq_len(num_dirs)) {
    phi_j <- phis[j]
    dir_vec <- cos(phi_j) * u_tip + sin(phi_j) * u_right
    h_v <- as.numeric(xy_n[, 1L] * dir_vec[1L] + xy_n[, 2L] * dir_vec[2L])
    h_e <- pmax(h_v, h_v[next_idx])
    # Exact Euler characteristic of EmbeddedGraph cycle: chi(phi_j, t_i) = #V(h <= t_i) - #E(max(h_u, h_v) <= t_i)
    n_v_le <- findInterval(thresh, sort(h_v))
    n_e_le <- findInterval(thresh, sort(h_e))
    ect_mat[, j] <- as.numeric(n_v_le - n_e_le)
  }
  ect_mat
}

.summarize_and_symmetrize_ect_matrix <- function(ect_mat) {
  if (is.null(ect_mat) || !is.matrix(ect_mat) || nrow(ect_mat) < 4L || ncol(ect_mat) < 4L) {
    return(NULL)
  }
  nr <- nrow(ect_mat); nc <- ncol(ect_mat)
  # Step 02 tip_up_directions() and compute_cycle_ect_r() place phi = 0 (lamina_tip)
  # strictly at Column 1; never re-roll onto wider lateral axes on reniform (W > L) blades!
  # Column 1 is phi = 0 (North/Tip); Columns 2..nc reflect via nc + 2 - j
  refl_cols <- c(1L, seq.int(nc, 2L, by = -1L))
  ect_refl <- ect_mat[, refl_cols, drop = FALSE]
  ect_sym <- 0.5 * (ect_mat + ect_refl)
  tot_norm2 <- sum(ect_mat^2, na.rm = TRUE)
  asym_val <- if (is.finite(tot_norm2) && tot_norm2 > 1e-8) {
    min(1.0, sum((ect_mat - ect_refl)^2, na.rm = TRUE) / (2.0 * tot_norm2))
  } else {
    0.0
  }
  # Total variation along each direction column j: convex shape has TV_j = 2; extra lobes/sinuses increase TV_j and max chi > 1
  tv_cols <- colSums(abs(diff(rbind(0.0, ect_mat, 0.0))), na.rm = TRUE)
  excess_tv <- mean(pmax(0.0, (tv_cols - 2.0) / 2.0), na.rm = TRUE)
  multi_comp <- mean(pmax(0.0, ect_mat - 1.0), na.rm = TRUE)
  lobe_comp <- round(as.numeric(excess_tv + multi_comp), 4)
  half_nc <- floor(nc / 2L) + 1L
  sym_half <- ect_sym[, seq_len(half_nc), drop = FALSE]
  sym_vec <- as.numeric(sign(sym_half) * log1p(abs(sym_half)))
  list(
    ect_mat = ect_mat,
    ect_sym = ect_sym,
    sym_vec = sym_vec,
    lobe_complexity = lobe_comp,
    bilateral_asym = round(as.numeric(asym_val), 5)
  )
}

compute_leaf_ect_features <- function(efa_df, coo_list = NULL, contours_dir = NULL) {
  if (is.null(efa_df) || nrow(efa_df) == 0L) return(efa_df)
  n_leaves <- nrow(efa_df)
  ect_dirs <- unique(c(
    if (!is.null(contours_dir) && nzchar(contours_dir)) file.path(dirname(contours_dir), "ect") else character(0),
    "data/ect"
  ))
  ect_dirs <- ect_dirs[dir.exists(ect_dirs)]

  raw_mats <- vector("list", n_leaves)
  ect_mats <- vector("list", n_leaves)
  sym_vecs <- vector("list", n_leaves)
  lobe_comp <- rep(NA_real_, n_leaves)
  bilat_asym <- rep(NA_real_, n_leaves)

  for (i in seq_len(n_leaves)) {
    cat_i <- trimws(as.character(efa_df$catalogNumber[i]))
    pid_i <- if ("plant_individual_id" %in% names(efa_df) && !is.na(efa_df$plant_individual_id[i])) as.integer(efa_df$plant_individual_id[i]) else 0L
    lid_i <- if ("leaf_id" %in% names(efa_df) && !is.na(efa_df$leaf_id[i])) as.integer(efa_df$leaf_id[i]) else i
    sid_i <- if ("shape_id" %in% names(efa_df) && !is.na(efa_df$shape_id[i])) trimws(as.character(efa_df$shape_id[i])) else ""
    ep_i <- if ("ect_path" %in% names(efa_df)) as.character(efa_df$ect_path[i]) else NA_character_
    cand_paths <- character(0)
    if (!is.na(ep_i) && nzchar(trimws(ep_i))) cand_paths <- c(cand_paths, trimws(ep_i))
    for (ed in ect_dirs) {
      if (nzchar(sid_i)) cand_paths <- c(cand_paths, file.path(ed, sprintf("%s_ect.npy", sid_i)))
      if (pid_i > 0L) cand_paths <- c(cand_paths, file.path(ed, sprintf("%s_p%d_leaf%d_ect.npy", cat_i, pid_i, lid_i)))
      cand_paths <- c(cand_paths, file.path(ed, sprintf("%s_leaf%d_ect.npy", cat_i, lid_i)))
    }
    for (cp in unique(cand_paths)) {
      if (file.exists(cp)) {
        m_try <- read_numpy_ect_r(cp)
        if (!is.null(m_try)) {
          raw_mats[[i]] <- m_try
          break
        }
      }
    }
  }

  # Resolve canonical (num_thresh, num_dirs) from loaded .npy matrices (default 32 x 360 per LM3_settings.yaml)
  loaded_dims <- Filter(Negate(is.null), lapply(raw_mats, dim))
  target_nr <- if (length(loaded_dims) > 0L) as.integer(stats::median(vapply(loaded_dims, `[`, integer(1L), 1L))) else 32L
  target_nc <- if (length(loaded_dims) > 0L) as.integer(stats::median(vapply(loaded_dims, `[`, integer(1L), 2L))) else 360L

  for (i in seq_len(n_leaves)) {
    mat_i <- raw_mats[[i]]
    if (is.null(mat_i) && !is.null(coo_list)) {
      sid_i <- if ("shape_id" %in% names(efa_df) && !is.na(efa_df$shape_id[i])) as.character(efa_df$shape_id[i]) else rownames(efa_df)[i]
      coo_i <- if (!is.null(sid_i) && sid_i %in% names(coo_list)) coo_list[[sid_i]] else if (i <= length(coo_list)) coo_list[[i]] else NULL
      mat_i <- compute_cycle_ect_r(coo_i, num_thresh = target_nr, num_dirs = target_nc)
    }
    if (!is.null(mat_i)) {
      mat_i <- .resample_ect_directions_r(mat_i, target_nr = target_nr, target_nc = target_nc)
      s_res <- .summarize_and_symmetrize_ect_matrix(mat_i)
      if (!is.null(s_res)) {
        ect_mats[[i]] <- s_res$ect_sym
        sym_vecs[[i]] <- s_res$sym_vec
        lobe_comp[i] <- s_res$lobe_complexity
        bilat_asym[i] <- s_res$bilateral_asym
      }
    }
  }

  efa_df$ect_lobe_complexity <- lobe_comp
  efa_df$ect_bilateral_asym <- bilat_asym
  efa_df$ECT_PC1 <- NA_real_
  efa_df$ECT_PC2 <- NA_real_
  efa_df$ECT_PC3 <- NA_real_

  valid_idx <- which(!vapply(sym_vecs, is.null, logical(1)))
  if (length(valid_idx) >= 2L) {
    M_ect <- do.call(rbind, sym_vecs[valid_idx])
    storage.mode(M_ect) <- "double"
    M_ect[!is.finite(M_ect)] <- 0.0
    lz <- if ("leaf_zone" %in% names(efa_df)) tolower(trimws(as.character(efa_df$leaf_zone[valid_idx]))) else rep("basal", length(valid_idx))
    out_v <- if ("is_contour_outlier" %in% names(efa_df)) as.logical(efa_df$is_contour_outlier[valid_idx]) else rep(FALSE, length(valid_idx))
    typ_v <- if ("is_type_specimen" %in% names(efa_df)) (!is.na(efa_df$is_type_specimen[valid_idx]) & as.logical(efa_df$is_type_specimen[valid_idx])) else rep(FALSE, length(valid_idx))
    out_v[is.na(out_v) | typ_v] <- FALSE
    wt_v <- if ("authority_weight" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$authority_weight[valid_idx])) else rep(1.0, length(valid_idx))
    wt_v[!is.finite(wt_v) | wt_v < 0.0] <- 0.0
    sb_v <- if ("spatial_balance_weight" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$spatial_balance_weight[valid_idx])) else rep(1.0, length(valid_idx))
    sb_v[!is.finite(sb_v) | sb_v <= 0.0] <- 1.0
    eff_wt_v <- wt_v * sb_v

    bas_pos <- which((is.na(lz) | lz == "basal") & !out_v)
    if (length(bas_pos) < 2L) bas_pos <- which(is.na(lz) | lz == "basal")
    if (length(bas_pos) < 2L) bas_pos <- seq_along(valid_idx)
    ref_pos <- bas_pos[wt_v[bas_pos] > 0.0]
    if (length(ref_pos) < 2L) ref_pos <- bas_pos

    M_fit <- M_ect[ref_pos, , drop = FALSE]
    w_fit <- if (any(eff_wt_v[ref_pos] > 0.0)) eff_wt_v[ref_pos] else wt_v[ref_pos]
    if (sum(w_fit, na.rm = TRUE) <= 1e-12) w_fit <- rep(1.0, nrow(M_fit))
    w_norm <- w_fit / sum(w_fit)
    col_v <- apply(M_fit, 2L, stats::var, na.rm = TRUE)
    act_c <- which(is.finite(col_v) & col_v > 1e-10)
    if (length(act_c) >= 1L) {
      mu_fit <- as.numeric(crossprod(w_norm, M_fit[, act_c, drop = FALSE]))
      X_fit <- sweep(M_fit[, act_c, drop = FALSE], 2L, mu_fit, "-")
      Z_fit <- sweep(X_fit, 1L, sqrt(w_norm), "*")
      svd_e <- tryCatch(svd(Z_fit, nu = 0L, nv = min(3L, nrow(Z_fit), ncol(Z_fit))), error = function(e) NULL)
      if (!is.null(svd_e) && !is.null(svd_e$v) && ncol(svd_e$v) >= 1L) {
        X_all <- sweep(M_ect[, act_c, drop = FALSE], 2L, mu_fit, "-")
        sc_all <- X_all %*% svd_e$v
        sc_ref <- sc_all[ref_pos, 1L]
        lc_ref <- lobe_comp[valid_idx[ref_pos]]
        # Align ECT_PC1 sign positively with ect_lobe_complexity strictly on active basal reference anchors
        if (sum(is.finite(lc_ref)) >= 3L && stats::sd(sc_ref, na.rm = TRUE) > 1e-8 && stats::sd(lc_ref, na.rm = TRUE) > 1e-8) {
          r_sign <- suppressWarnings(stats::cor(sc_ref, lc_ref, use = "complete.obs"))
          if (is.finite(r_sign) && r_sign < 0.0) {
            sc_all[, 1L] <- -sc_all[, 1L]
            svd_e$v[, 1L] <- -svd_e$v[, 1L]
          }
        }
        efa_df$ECT_PC1[valid_idx] <- round(sc_all[, 1L], 6)
        if (ncol(sc_all) >= 2L) efa_df$ECT_PC2[valid_idx] <- round(sc_all[, 2L], 6) else efa_df$ECT_PC2[valid_idx] <- 0.0
        if (ncol(sc_all) >= 3L) efa_df$ECT_PC3[valid_idx] <- round(sc_all[, 3L], 6) else efa_df$ECT_PC3[valid_idx] <- 0.0
        denom_w <- max(1e-8, 1.0 - sum(w_norm^2))
        tot_wvar <- sum(colSums(Z_fit^2, na.rm = TRUE)) / denom_w
        ev_vals <- (svd_e$d^2) / denom_w
        attr(efa_df, "ect_basal_var_exp") <- round(ev_vals / max(tot_wvar, 1e-12) * 100.0, 2)
      }
    }
    attr(efa_df, "ect_sym_dim") <- ncol(M_ect)
  }
  attr(efa_df, "ect_sym_mats") <- ect_mats
  attr(efa_df, "ect_sym_vecs") <- sym_vecs
  attr(efa_df, "ect_coo_list") <- coo_list
  efa_df
}

compute_separate_cauline_otu_track <- function(efa_df, specimen_df, sym_names = NULL) {
  if (is.null(specimen_df) || nrow(specimen_df) == 0L) return(specimen_df)
  n_spec <- nrow(specimen_df)
  cauline_leaf_count <- integer(n_spec)
  cauline_solidity_med <- rep(NA_real_, n_spec)
  cauline_dissection_index_med <- rep(NA_real_, n_spec)
  cauline_aspect_ratio_med <- rep(NA_real_, n_spec)
  cauline_blade_length_mm <- rep(NA_real_, n_spec)
  basal_vs_cauline_length_ratio <- rep(NA_real_, n_spec)
  basal_vs_cauline_solidity_delta <- rep(NA_real_, n_spec)
  heteroblasty_efa_distance <- rep(NA_real_, n_spec)
  Cauline_PC1 <- rep(NA_real_, n_spec)
  Cauline_PC2 <- rep(NA_real_, n_spec)
  Cauline_PC3 <- rep(NA_real_, n_spec)
  ECT_PC1 <- rep(NA_real_, n_spec)
  ECT_PC2 <- rep(NA_real_, n_spec)
  ECT_PC3 <- rep(NA_real_, n_spec)
  Cauline_ECT_PC1 <- rep(NA_real_, n_spec)
  Cauline_ECT_PC2 <- rep(NA_real_, n_spec)
  ect_lobe_complexity <- rep(NA_real_, n_spec)
  cauline_ect_lobe_complexity <- rep(NA_real_, n_spec)
  heteroblasty_ect_distance <- rep(NA_real_, n_spec)

  if (!is.null(efa_df) && nrow(efa_df) > 0L) {
    if (is.null(sym_names) || length(sym_names) == 0L) {
      sym_names <- grep("^[AD][0-9]+$", names(efa_df), value = TRUE)
    }
    sym_names <- intersect(as.character(sym_names), names(efa_df))
    ect_sym_vecs <- attr(efa_df, "ect_sym_vecs")
    pids_leaf <- if ("plant_individual_id" %in% names(efa_df)) ifelse(is.na(efa_df$plant_individual_id), 0L, as.integer(efa_df$plant_individual_id)) else rep(0L, nrow(efa_df))
    pids_spec <- if ("plant_individual_id" %in% names(specimen_df)) ifelse(is.na(specimen_df$plant_individual_id), 0L, as.integer(specimen_df$plant_individual_id)) else rep(0L, n_spec)
    leaf_keys <- paste0(trimws(as.character(efa_df$catalogNumber)), "__", pids_leaf)
    spec_keys <- paste0(trimws(as.character(specimen_df$catalogNumber)), "__", pids_spec)
    idx_by_spec <- split(seq_len(nrow(efa_df)), leaf_keys)

    med_col <- function(df_sub, cn, digits = 4L, pos_only = FALSE) {
      if (!(cn %in% names(df_sub))) return(NA_real_)
      v <- suppressWarnings(as.numeric(df_sub[[cn]]))
      v <- if (pos_only) v[is.finite(v) & v > 0.0] else v[is.finite(v)]
      if (length(v) == 0L) NA_real_ else round(stats::median(v), digits)
    }

    caul_harm_list <- vector("list", n_spec)
    caul_ect_list <- vector("list", n_spec)
    for (i in seq_len(n_spec)) {
      r_idx <- idx_by_spec[[spec_keys[i]]]
      if (length(r_idx) == 0L) next
      grp <- efa_df[r_idx, , drop = FALSE]
      keep_sub <- rep(TRUE, nrow(grp))

      is_type_grp <- if ("is_type_specimen" %in% names(grp)) (!is.na(grp$is_type_specimen) & as.logical(grp$is_type_specimen)) else rep(FALSE, nrow(grp))
      if ("is_contour_outlier" %in% names(grp) && any((!as.logical(grp$is_contour_outlier[keep_sub])) | is_type_grp[keep_sub], na.rm = TRUE)) {
        keep_sub <- keep_sub & ((!is.na(grp$is_contour_outlier) & !as.logical(grp$is_contour_outlier)) | is_type_grp)
      }
      grp <- grp[keep_sub, , drop = FALSE]
      r_idx_kept <- r_idx[keep_sub]
      lz <- if ("leaf_zone" %in% names(grp)) tolower(trimws(as.character(grp$leaf_zone))) else rep("basal", nrow(grp))
      lz[is.na(lz) | !(lz %in% c("basal", "cauline"))] <- "basal"
      basal_rows <- grp[lz == "basal", , drop = FALSE]
      cauline_rows <- grp[lz == "cauline", , drop = FALSE]
      bas_orig_idx <- r_idx_kept[lz == "basal"]
      caul_orig_idx <- r_idx_kept[lz == "cauline"]
      n_bas <- nrow(basal_rows)
      n_caul <- nrow(cauline_rows)
      cauline_leaf_count[i] <- as.integer(n_caul)
      if ("basal_leaf_count" %in% names(specimen_df)) specimen_df$basal_leaf_count[i] <- as.integer(n_bas)
      if ("leaf_count" %in% names(specimen_df)) specimen_df$leaf_count[i] <- as.integer(n_bas + n_caul)
      if (n_bas == 0L && n_caul >= 1L && ("morph_tier_source" %in% names(specimen_df))) {
        specimen_df$morph_tier_source[i] <- "Cauline_Only_Fallback"
      }

      prim_rows <- if (n_bas >= 1L) basal_rows else cauline_rows
      if (nrow(prim_rows) >= 1L) {
        ECT_PC1[i] <- med_col(prim_rows, "ECT_PC1", 6L)
        ECT_PC2[i] <- med_col(prim_rows, "ECT_PC2", 6L)
        ECT_PC3[i] <- med_col(prim_rows, "ECT_PC3", 6L)
        ect_lobe_complexity[i] <- med_col(prim_rows, "ect_lobe_complexity", 4L)
        if ("n_lamina_components" %in% names(prim_rows) && any(!is.na(prim_rows$n_lamina_components) & as.integer(prim_rows$n_lamina_components) > 1L)) {
          if ("solidity" %in% names(specimen_df)) specimen_df$solidity[i] <- med_col(prim_rows, "solidity", 4L)
          if ("q25_solidity" %in% names(specimen_df) && "solidity" %in% names(prim_rows)) {
            sv_p <- suppressWarnings(as.numeric(prim_rows$solidity)); sv_p <- sv_p[is.finite(sv_p)]
            if (length(sv_p) >= 1L) specimen_df$q25_solidity[i] <- round(as.numeric(stats::quantile(sv_p, probs = 0.25, names = FALSE)), 4)
          }
          if ("dissection_index" %in% names(specimen_df)) specimen_df$dissection_index[i] <- med_col(prim_rows, "dissection_index", 4L)
          if ("dissected_leaf_fraction" %in% names(specimen_df) && "is_dissected" %in% names(prim_rows)) {
            dv_p <- as.logical(prim_rows$is_dissected); dv_p <- dv_p[!is.na(dv_p)]
            if (length(dv_p) >= 1L) specimen_df$dissected_leaf_fraction[i] <- round(mean(dv_p), 4)
          }
        }
      }

      if (n_caul >= 1L) {
        cauline_solidity_med[i] <- med_col(cauline_rows, "solidity", 4L)
        cauline_dissection_index_med[i] <- med_col(cauline_rows, "dissection_index", 4L)
        cauline_aspect_ratio_med[i] <- med_col(cauline_rows, "aspect_ratio", 4L)
        cauline_ect_lobe_complexity[i] <- med_col(cauline_rows, "ect_lobe_complexity", 4L)
        c_bl <- med_col(cauline_rows, "blade_length_mm", 3L, pos_only = TRUE)
        c_cs <- med_col(cauline_rows, "centroid_size", 3L, pos_only = TRUE)
        c_cs_px <- med_col(cauline_rows, "centroid_size_px", 3L, pos_only = TRUE)
        c_apx <- med_col(cauline_rows, "area_px", 1L, pos_only = TRUE)
        cauline_blade_length_mm[i] <- if (is.finite(c_bl)) c_bl else NA_real_

        b_bl <- if (n_bas >= 1L) med_col(basal_rows, "blade_length_mm", 4L, pos_only = TRUE) else NA_real_
        if (!is.finite(b_bl) && ("blade_length_mm" %in% names(specimen_df))) {
          v_b <- suppressWarnings(as.numeric(specimen_df$blade_length_mm[i]))
          if (is.finite(v_b) && v_b > 0.0) b_bl <- v_b
        }
        b_cs <- if (n_bas >= 1L) med_col(basal_rows, "centroid_size", 4L, pos_only = TRUE) else NA_real_
        if (!is.finite(b_cs) && ("centroid_size" %in% names(specimen_df))) {
          v_c <- suppressWarnings(as.numeric(specimen_df$centroid_size[i]))
          if (is.finite(v_c) && v_c > 0.0) b_cs <- v_c
        }
        b_cs_px <- if (n_bas >= 1L) med_col(basal_rows, "centroid_size_px", 3L, pos_only = TRUE) else NA_real_
        b_apx <- if (n_bas >= 1L) med_col(basal_rows, "area_px", 1L, pos_only = TRUE) else NA_real_
        if (is.finite(b_bl) && is.finite(c_bl) && c_bl > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_bl / c_bl, 4)
        } else if (is.finite(b_cs) && is.finite(c_cs) && c_cs > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_cs / c_cs, 4)
        } else if (is.finite(b_cs_px) && is.finite(c_cs_px) && c_cs_px > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_cs_px / c_cs_px, 4)
        } else if (is.finite(b_apx) && is.finite(c_apx) && c_apx > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(sqrt(b_apx / c_apx), 4)
        }

        b_sol <- if ("solidity" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$solidity[i])) else med_col(basal_rows, "solidity", 4L)
        if (is.finite(b_sol) && is.finite(cauline_solidity_med[i])) {
          basal_vs_cauline_solidity_delta[i] <- round(b_sol - cauline_solidity_med[i], 4)
        }

        if (length(sym_names) >= 1L) {
          caul_mat_i <- as.matrix(cauline_rows[, sym_names, drop = FALSE])
          storage.mode(caul_mat_i) <- "double"
          if (any(complete.cases(caul_mat_i))) {
            caul_med_h <- apply(caul_mat_i, 2L, stats::median, na.rm = TRUE)
            caul_med_h[!is.finite(caul_med_h)] <- 0.0
            caul_harm_list[[i]] <- caul_med_h
            if (n_bas >= 1L) {
              bas_mat_i <- as.matrix(basal_rows[, sym_names, drop = FALSE])
              storage.mode(bas_mat_i) <- "double"
              if (any(complete.cases(bas_mat_i))) {
                bas_med_h <- apply(bas_mat_i, 2L, stats::median, na.rm = TRUE)
                bas_med_h[!is.finite(bas_med_h)] <- 0.0
                d_efa <- sqrt(sum((bas_med_h - caul_med_h)^2, na.rm = TRUE))
                if (is.finite(d_efa)) {
                  heteroblasty_efa_distance[i] <- round(d_efa, 6)
                }
              }
            }
          }
        }

        if (!is.null(ect_sym_vecs) && length(caul_orig_idx) >= 1L) {
          c_vecs <- Filter(Negate(is.null), ect_sym_vecs[caul_orig_idx])
          if (length(c_vecs) >= 1L) {
            c_em <- do.call(rbind, c_vecs)
            c_med_e <- apply(c_em, 2L, stats::median, na.rm = TRUE)
            c_med_e[!is.finite(c_med_e)] <- 0.0
            caul_ect_list[[i]] <- c_med_e
            if (n_bas >= 1L && length(bas_orig_idx) >= 1L) {
              b_vecs <- Filter(Negate(is.null), ect_sym_vecs[bas_orig_idx])
              if (length(b_vecs) >= 1L) {
                b_em <- do.call(rbind, b_vecs)
                b_med_e <- apply(b_em, 2L, stats::median, na.rm = TRUE)
                b_med_e[!is.finite(b_med_e)] <- 0.0
                d_ect <- sqrt(mean((b_med_e - c_med_e)^2, na.rm = TRUE))
                if (is.finite(d_ect)) heteroblasty_ect_distance[i] <- round(d_ect, 6)
              }
            }
          }
        } else if (all(c("ECT_PC1", "ECT_PC2", "ECT_PC3") %in% names(cauline_rows))) {
          c_epc <- c(
            med_col(cauline_rows, "ECT_PC1", 6L),
            med_col(cauline_rows, "ECT_PC2", 6L),
            med_col(cauline_rows, "ECT_PC3", 6L)
          )
          c_epc[!is.finite(c_epc)] <- 0.0
          caul_ect_list[[i]] <- c_epc
          if (n_bas >= 1L) {
            b_epc <- c(
              med_col(basal_rows, "ECT_PC1", 6L),
              med_col(basal_rows, "ECT_PC2", 6L),
              med_col(basal_rows, "ECT_PC3", 6L)
            )
            b_epc[!is.finite(b_epc)] <- 0.0
            ect_dim_ref <- attr(efa_df, "ect_sym_dim")
            ect_dim_use <- if (!is.null(ect_dim_ref) && is.finite(ect_dim_ref) && ect_dim_ref >= 3L) as.numeric(ect_dim_ref) else 5792.0
            d_ect <- sqrt(sum((b_epc - c_epc)^2, na.rm = TRUE) / ect_dim_use)
            if (is.finite(d_ect)) heteroblasty_ect_distance[i] <- round(d_ect, 6)
          }
        }
      }
    }

    auth_spec_v <- if ("authority_weight" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$authority_weight)) else rep(1.0, n_spec)
    auth_spec_v[!is.finite(auth_spec_v) | auth_spec_v < 0.0] <- 0.0
    sb_spec_v <- if ("spatial_balance_weight" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$spatial_balance_weight)) else rep(1.0, n_spec)
    sb_spec_v[!is.finite(sb_spec_v) | sb_spec_v <= 0.0] <- 1.0
    eff_spec_v <- auth_spec_v * sb_spec_v

    caul_spec_idx <- which(cauline_leaf_count >= 1L & !vapply(caul_harm_list, is.null, logical(1)))
    if (length(caul_spec_idx) >= 1L && length(sym_names) >= 1L) {
      M_caul <- do.call(rbind, caul_harm_list[caul_spec_idx])
      storage.mode(M_caul) <- "double"
      ref_c_pos <- which(auth_spec_v[caul_spec_idx] > 0.0)
      if (length(ref_c_pos) < 3L) ref_c_pos <- seq_len(nrow(M_caul))
      for (jc in seq_len(ncol(M_caul))) {
        bad_c <- !is.finite(M_caul[, jc])
        if (any(bad_c)) {
          ok_r <- intersect(ref_c_pos, which(!bad_c))
          M_caul[bad_c, jc] <- if (length(ok_r) > 0L) stats::median(M_caul[ok_r, jc], na.rm = TRUE) else if (any(!bad_c)) stats::median(M_caul[!bad_c, jc], na.rm = TRUE) else 0.0
        }
      }
      c_scores <- NULL
      if (length(ref_c_pos) >= 3L) {
        M_ref_c <- M_caul[ref_c_pos, , drop = FALSE]
        c_vars <- apply(M_ref_c, 2L, stats::var, na.rm = TRUE)
        act_c <- which(is.finite(c_vars) & c_vars > 1e-12)
        if (length(act_c) >= 1L) {
          w_c <- eff_spec_v[caul_spec_idx[ref_c_pos]]
          if (sum(w_c, na.rm = TRUE) <= 1e-12) w_c <- rep(1.0, length(ref_c_pos))
          wt_c <- stats::cov.wt(M_ref_c[, act_c, drop = FALSE], wt = w_c, center = TRUE)
          eig_c <- eigen(wt_c$cov, symmetric = TRUE)
          c_scores <- sweep(M_caul[, act_c, drop = FALSE], 2L, wt_c$center, "-") %*% eig_c$vectors
        }
      }
      if (is.null(c_scores)) {
        all_sym_m <- as.matrix(efa_df[, sym_names, drop = FALSE])
        storage.mode(all_sym_m) <- "double"
        all_sym_m[!is.finite(all_sym_m)] <- 0.0
        auth_leaf_v <- if ("authority_weight" %in% names(efa_df)) suppressWarnings(as.numeric(efa_df$authority_weight)) else rep(1.0, nrow(efa_df))
        ref_l_pos <- which(is.finite(auth_leaf_v) & auth_leaf_v > 0.0)
        if (length(ref_l_pos) < 2L) ref_l_pos <- seq_len(nrow(all_sym_m))
        if (length(ref_l_pos) >= 2L) {
          c_vars <- apply(all_sym_m[ref_l_pos, , drop = FALSE], 2L, stats::var, na.rm = TRUE)
          act_c <- which(is.finite(c_vars) & c_vars > 1e-12)
          if (length(act_c) >= 1L) {
            mu_all <- colMeans(all_sym_m[ref_l_pos, act_c, drop = FALSE])
            eig_all <- eigen(stats::cov(all_sym_m[ref_l_pos, act_c, drop = FALSE]), symmetric = TRUE)
            c_scores <- sweep(M_caul[, act_c, drop = FALSE], 2L, mu_all, "-") %*% eig_all$vectors
          }
        }
      }
      if (!is.null(c_scores) && nrow(c_scores) == length(caul_spec_idx)) {
        Cauline_PC1[caul_spec_idx] <- if (ncol(c_scores) >= 1L) round(c_scores[, 1L], 6) else 0.0
        Cauline_PC2[caul_spec_idx] <- if (ncol(c_scores) >= 2L) round(c_scores[, 2L], 6) else 0.0
        Cauline_PC3[caul_spec_idx] <- if (ncol(c_scores) >= 3L) round(c_scores[, 3L], 6) else 0.0
      }
    }

    caul_ect_idx <- which(cauline_leaf_count >= 1L & !vapply(caul_ect_list, is.null, logical(1)))
    if (length(caul_ect_idx) >= 1L) {
      M_cect <- do.call(rbind, caul_ect_list[caul_ect_idx])
      storage.mode(M_cect) <- "double"
      M_cect[!is.finite(M_cect)] <- 0.0
      ce_scores <- NULL
      ref_ce_pos <- which(auth_spec_v[caul_ect_idx] > 0.0)
      if (length(ref_ce_pos) < 2L) ref_ce_pos <- seq_len(nrow(M_cect))
      if (length(ref_ce_pos) >= 2L) {
        M_ref_ce <- M_cect[ref_ce_pos, , drop = FALSE]
        ce_vars <- apply(M_ref_ce, 2L, stats::var, na.rm = TRUE)
        act_ce <- which(is.finite(ce_vars) & ce_vars > 1e-12)
        if (length(act_ce) >= 1L) {
          w_ce <- eff_spec_v[caul_ect_idx[ref_ce_pos]]
          if (sum(w_ce, na.rm = TRUE) <= 1e-12) w_ce <- rep(1.0, length(ref_ce_pos))
          w_ce_n <- w_ce / sum(w_ce)
          mu_ce <- as.numeric(crossprod(w_ce_n, M_ref_ce[, act_ce, drop = FALSE]))
          X_ce_ref <- sweep(M_ref_ce[, act_ce, drop = FALSE], 2L, mu_ce, "-")
          Z_ce_ref <- sweep(X_ce_ref, 1L, sqrt(w_ce_n), "*")
          svd_ce <- tryCatch(svd(Z_ce_ref, nu = 0L, nv = min(2L, nrow(Z_ce_ref), ncol(Z_ce_ref))), error = function(e) NULL)
          if (!is.null(svd_ce) && !is.null(svd_ce$v) && ncol(svd_ce$v) >= 1L) {
            ce_scores <- sweep(M_cect[, act_ce, drop = FALSE], 2L, mu_ce, "-") %*% svd_ce$v
            denom_ce <- max(1e-8, 1.0 - sum(w_ce_n^2))
            ev_ce <- (svd_ce$d^2) / denom_ce
            tot_ce <- sum(colSums(Z_ce_ref^2, na.rm = TRUE)) / denom_ce
            attr(specimen_df, "ect_cauline_var_exp") <- round(ev_ce / max(tot_ce, 1e-12) * 100.0, 2)
          }
        }
      }
      if (!is.null(ce_scores) && nrow(ce_scores) == length(caul_ect_idx)) {
        Cauline_ECT_PC1[caul_ect_idx] <- if (ncol(ce_scores) >= 1L) round(ce_scores[, 1L], 6) else 0.0
        Cauline_ECT_PC2[caul_ect_idx] <- if (ncol(ce_scores) >= 2L) round(ce_scores[, 2L], 6) else 0.0
      } else {
        Cauline_ECT_PC1[caul_ect_idx] <- 0.0
        Cauline_ECT_PC2[caul_ect_idx] <- 0.0
      }
    }
  }

  specimen_df$cauline_leaf_count <- as.integer(cauline_leaf_count)
  specimen_df$cauline_solidity_med <- as.numeric(cauline_solidity_med)
  specimen_df$cauline_dissection_index_med <- as.numeric(cauline_dissection_index_med)
  specimen_df$cauline_aspect_ratio_med <- as.numeric(cauline_aspect_ratio_med)
  specimen_df$cauline_blade_length_mm <- as.numeric(cauline_blade_length_mm)
  specimen_df$basal_vs_cauline_length_ratio <- as.numeric(basal_vs_cauline_length_ratio)
  specimen_df$basal_vs_cauline_solidity_delta <- as.numeric(basal_vs_cauline_solidity_delta)
  specimen_df$heteroblasty_efa_distance <- as.numeric(heteroblasty_efa_distance)
  specimen_df$Cauline_PC1 <- as.numeric(Cauline_PC1)
  specimen_df$Cauline_PC2 <- as.numeric(Cauline_PC2)
  specimen_df$Cauline_PC3 <- as.numeric(Cauline_PC3)
  specimen_df$ECT_PC1 <- as.numeric(ECT_PC1)
  specimen_df$ECT_PC2 <- as.numeric(ECT_PC2)
  specimen_df$ECT_PC3 <- as.numeric(ECT_PC3)
  specimen_df$Cauline_ECT_PC1 <- as.numeric(Cauline_ECT_PC1)
  specimen_df$Cauline_ECT_PC2 <- as.numeric(Cauline_ECT_PC2)
  specimen_df$ect_lobe_complexity <- as.numeric(ect_lobe_complexity)
  specimen_df$cauline_ect_lobe_complexity <- as.numeric(cauline_ect_lobe_complexity)
  specimen_df$heteroblasty_ect_distance <- as.numeric(heteroblasty_ect_distance)
  attr(specimen_df, "ect_basal_var_exp") <- attr(efa_df, "ect_basal_var_exp")
  specimen_df
}

filter_biologically_variable_dissection_scalars <- function(df, cand_scalars, active_mask = NULL) {
  cand_scalars <- intersect(as.character(cand_scalars), names(df))
  if (length(cand_scalars) == 0L || nrow(df) < 5L) return(character(0))
  act <- if (!is.null(active_mask) && length(active_mask) == nrow(df) && sum(active_mask, na.rm = TRUE) >= 5L) {
    which(as.logical(active_mask))
  } else {
    seq_len(nrow(df))
  }
  cauline_het_scalars <- c(
    "heteroblasty_efa_distance", "basal_vs_cauline_length_ratio",
    "basal_vs_cauline_solidity_delta", "cauline_dissection_index_med",
    "cauline_solidity_med", "Cauline_PC1",
    "Cauline_ECT_PC1", "cauline_ect_lobe_complexity", "heteroblasty_ect_distance"
  )
  keep <- character(0)
  for (sc in cand_scalars) {
    raw_nm <- sub("_std$", "", sc)
    raw_v <- if (raw_nm %in% names(df)) suppressWarnings(as.numeric(df[[raw_nm]][act])) else suppressWarnings(as.numeric(df[[sc]][act]))
    raw_v <- raw_v[is.finite(raw_v)]
    if (length(raw_v) < 5L) next
    if ((raw_nm %in% cauline_het_scalars) && (length(raw_v) < ceiling(0.50 * length(act)))) next
    sd_v <- stats::sd(raw_v)
    q_span <- as.numeric(diff(stats::quantile(raw_v, probs = c(0.10, 0.90), names = FALSE)))
    ok <- if (raw_nm %in% c("solidity", "q25_solidity", "cauline_solidity_med")) {
      (is.finite(sd_v) && sd_v >= 0.045) || any(raw_v < 0.75) || (is.finite(q_span) && q_span >= 0.10)
    } else if (raw_nm %in% c("dissection_index", "cauline_dissection_index_med")) {
      (is.finite(sd_v) && sd_v >= 0.12) || any(raw_v > 1.35)
    } else if (raw_nm %in% c("within_rosette_efa_dispersion", "heteroblasty_efa_distance", "heteroblasty_ect_distance")) {
      (is.finite(sd_v) && sd_v >= 0.015) && (sum(raw_v > 0) >= 5L)
    } else if (raw_nm %in% c("basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "Cauline_PC1",
                             "ECT_PC1", "ECT_PC2", "Cauline_ECT_PC1", "ect_lobe_complexity", "cauline_ect_lobe_complexity")) {
      is.finite(sd_v) && sd_v >= 0.015
    } else {
      is.finite(sd_v) && sd_v > 1e-6
    }
    if (isTRUE(ok)) keep <- c(keep, sc)
  }
  keep
}

prune_collinear_dissection_scalars <- function(
  df,
  base_cols,
  cand_scalars,
  active_mask = NULL,
  group_vec = NULL,
  max_r = 0.75
) {
  base_cols <- intersect(as.character(base_cols), names(df))
  cand_scalars <- setdiff(intersect(as.character(cand_scalars), names(df)), base_cols)
  if (length(cand_scalars) == 0L || nrow(df) < 5L) return(character(0))
  act <- if (!is.null(active_mask) && length(active_mask) == nrow(df) && sum(active_mask, na.rm = TRUE) >= 5L) {
    which(as.logical(active_mask))
  } else {
    seq_len(nrow(df))
  }

  # Rank candidate scalars empirically by univariate ANOVA F on active reference anchors (MorphoTools2 protocol)
  f_scores <- vapply(cand_scalars, function(sc) {
    v <- suppressWarnings(as.numeric(df[[sc]][act]))
    fin <- is.finite(v)
    if (sum(fin) < 5L || stats::sd(v[fin]) <= 1e-8) return(-1.0)
    if (!is.null(group_vec) && length(group_vec) == nrow(df)) {
      g <- droplevels(factor(as.character(group_vec[act][fin])))
      if (nlevels(g) >= 2L && min(table(g)) >= 2L) {
        fv <- tryCatch(summary(stats::aov(v[fin] ~ g))[[1]][["F value"]][1L], error = function(e) NA_real_)
        if (is.finite(fv)) return(as.numeric(fv))
      }
    }
    stats::var(v[fin])
  }, numeric(1))

  valid_cands <- cand_scalars[f_scores >= 0.0]
  if (length(valid_cands) == 0L) return(character(0))
  ord_cands <- valid_cands[order(f_scores[valid_cands], decreasing = TRUE)]

  accepted <- base_cols
  kept_scalars <- character(0)
  for (sc in ord_cands) {
    v_sc <- suppressWarnings(as.numeric(df[[sc]][act]))
    if (sum(is.finite(v_sc)) < 5L || stats::sd(v_sc, na.rm = TRUE) <= 1e-8) next
    too_corr <- FALSE
    for (ac in accepted) {
      v_ac <- suppressWarnings(as.numeric(df[[ac]][act]))
      ok <- is.finite(v_sc) & is.finite(v_ac)
      if (sum(ok) >= 5L && stats::sd(v_ac[ok]) > 1e-8 && stats::sd(v_sc[ok]) > 1e-8) {
        r_val <- suppressWarnings(abs(stats::cor(v_sc[ok], v_ac[ok])))
        if (is.finite(r_val) && r_val >= max_r) {
          too_corr <- TRUE
          break
        }
      }
    }
    if (!too_corr) {
      accepted <- c(accepted, sc)
      kept_scalars <- c(kept_scalars, sc)
    }
  }
  kept_scalars
}

fast_euclidean_permanova_trace <- function(X_mat, grp_vec, permutations = 999L, seed = 42L) {
  X <- as.matrix(X_mat)
  storage.mode(X) <- "double"
  g_fac <- droplevels(factor(as.character(grp_vec)))
  ok <- stats::complete.cases(X) & !is.na(g_fac)
  X <- X[ok, , drop = FALSE]
  g_fac <- droplevels(g_fac[ok])
  N <- nrow(X)
  g <- nlevels(g_fac)
  if (N <= g || g < 2L || ncol(X) < 1L) {
    return(list(F = NA_real_, R2 = NA_real_, p_value = NA_real_, df_group = max(0L, g - 1L), df_res = max(0L, N - g)))
  }

  X_c <- scale(X, center = TRUE, scale = FALSE)
  ss_tot <- sum(X_c^2)
  if (!is.finite(ss_tot) || ss_tot <= 1e-12) {
    return(list(F = 0.0, R2 = 0.0, p_value = 1.0, df_group = g - 1L, df_res = N - g))
  }

  g_int <- as.integer(g_fac)
  n_vec <- as.numeric(tabulate(g_int, nbins = g))
  grp_sums <- rowsum(X_c, g_int, reorder = TRUE)
  ss_between <- sum(rowSums(grp_sums^2) / n_vec)
  ss_within <- max(0.0, ss_tot - ss_between)
  df1 <- as.numeric(g - 1L)
  df2 <- as.numeric(N - g)
  r2_obs <- ss_between / ss_tot
  f_obs <- if (ss_within > 1e-12) (ss_between / df1) / (ss_within / df2) else Inf

  B <- max(0L, as.integer(permutations))
  p_val <- NA_real_
  if (B >= 1L) {
    set.seed(as.integer(seed))
    tol <- ss_between - 1e-12 * max(1.0, ss_tot)
    if (N * B >= 250000L && exists("packera_mclapply", mode = "function")) {
      perm_list <- replicate(B, sample.int(N), simplify = FALSE)
      n_w <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(16L, B) else 4L
      chunks <- split(seq_len(B), cut(seq_len(B), breaks = min(n_w, B), labels = FALSE))
      cnts <- packera_mclapply(chunks, function(idx_b) {
        c_loc <- 0L
        for (bi in idx_b) {
          sums_b <- rowsum(X_c[perm_list[[bi]], , drop = FALSE], g_int, reorder = TRUE)
          if (sum(rowSums(sums_b^2) / n_vec) >= tol) c_loc <- c_loc + 1L
        }
        c_loc
      }, mc.cores = length(chunks))
      ge_cnt <- sum(unlist(cnts))
    } else {
      ge_cnt <- 0L
      for (b in seq_len(B)) {
        g_perm <- sample.int(N)
        sums_b <- rowsum(X_c[g_perm, , drop = FALSE], g_int, reorder = TRUE)
        ss_b <- sum(rowSums(sums_b^2) / n_vec)
        if (ss_b >= tol) ge_cnt <- ge_cnt + 1L
      }
    }
    p_val <- (1.0 + ge_cnt) / (1.0 + B)
  }

  list(
    F = as.numeric(f_obs),
    R2 = as.numeric(r2_obs),
    p_value = as.numeric(p_val),
    SumOfSqs = as.numeric(ss_between),
    Residual_SS = as.numeric(ss_within),
    Total_SS = as.numeric(ss_tot),
    df_group = as.integer(df1),
    df_res = as.integer(df2)
  )
}
