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
  ba_cm2 <- if ("blade_area_cm2" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$blade_area_cm2)) else rep(NA_real_, n)
  min_ba <- if ("effective_min_blade_area_cm2" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$effective_min_blade_area_cm2)) else rep(NA_real_, n)
  area_ok <- ifelse(is.finite(ba_cm2) & is.finite(min_ba) & min_ba > 0, ba_cm2 >= min_ba, TRUE)

  leaf_zone <- raw_zone
  is_fallback_vec <- rep(FALSE, n)
  is_caul_only_vec <- rep(FALSE, n)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)

  for (g_idx in idx_by_grp) {
    lam_ok <- g_idx[area_ok[g_idx]]
    if (length(lam_ok) == 0L) lam_ok <- g_idx
    if (length(lam_ok) == 0L) next
    all_explicit_cauline <- all(!is.na(raw_zone[lam_ok]) & raw_zone[lam_ok] == "cauline")

    a_pool <- ifelse(is.finite(area_px[lam_ok]) & area_px[lam_ok] > 0, area_px[lam_ok],
                     ifelse(is.finite(ba_cm2[lam_ok]) & ba_cm2[lam_ok] > 0, ba_cm2[lam_ok] * 100.0, 1.0))
    max_pid_a <- max(a_pool, 1.0, na.rm = TRUE)

    for (ii in g_idx) {
      if (!is.na(leaf_zone[ii])) {
        if (!is.finite(zone_conf[ii])) zone_conf[ii] <- 0.90
        next
      }
      rel_a <- if (is.finite(area_px[ii]) && area_px[ii] > 0) area_px[ii] / max_pid_a else if (is.finite(ba_cm2[ii]) && ba_cm2[ii] > 0) (ba_cm2[ii] * 100.0) / max_pid_a else 1.0
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
          fb_hit <- seq_along(lam_ok) == which.max(a_pool)
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
  fac_df$area_ok_flag <- area_ok
  fac_df
}

partition_basal_and_cauline_contours <- classify_and_partition_leaf_zones

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
    raw_dis | (is.finite(sol_vec) & sol_vec < 0.68),
    (is.finite(sol_vec) & sol_vec < 0.68)
  )

  fac_df <- classify_and_partition_leaf_zones(
    fac_df,
    min_cy_plant = min_cy_plant,
    min_cy_sheet = min_cy_sheet,
    fallback_band = fallback_band
  )
  area_ok <- as.logical(fac_df$area_ok_flag)
  fac_df$area_ok_flag <- NULL
  is_fallback_vec <- as.logical(fac_df$is_basal_fallback)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)
  n_fb_individuals <- sum(vapply(idx_by_grp, function(g_idx) any(is_fallback_vec[g_idx] & !fac_df$is_cauline_only_individual[g_idx]), logical(1)))

  basal_keep <- area_ok & ((fac_df$leaf_zone == "basal") | is_fallback_vec)
  for (g_idx in idx_by_grp) {
    if (!any(basal_keep[g_idx])) {
      lam_ok <- g_idx[area_ok[g_idx]]
      if (length(lam_ok) == 0L) lam_ok <- g_idx
      if (length(lam_ok) > 0L) basal_keep[lam_ok[1L]] <- TRUE
    }
  }
  if (!any(basal_keep)) basal_keep <- rep(TRUE, n)

  n_area_drop <- sum(!area_ok & !basal_keep)
  n_caul_drop <- sum(area_ok & !basal_keep)
  keep_mask <- if (isTRUE(retain_cauline_tagged)) {
    km <- area_ok | basal_keep
    if (!any(km)) rep(TRUE, n) else km
  } else {
    basal_keep
  }

  if (n_area_drop > 0L || n_caul_drop > 0L || n_fb_individuals > 0L) {
    message(sprintf(
      paste0(
        "Positional Homology & Basal Lamina Gate: retained %d / %d Tier 1 basal blades ",
        "(excluded %d upper cauline, %d sub-floor blade fragments; ",
        "lowest-lamina fallback individuals = %d)."
      ),
      sum(basal_keep), n, n_caul_drop, n_area_drop, n_fb_individuals
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
    ucs_score = med_num("ucs_score", 4L),
    area_px = med_num("area_px", 1L),
    sym_var = med_num("sym_var", 6L),
    asym_var = med_num("asym_var", 6L),
    pct_asymmetry = med_num("pct_asymmetry", 6L),
    harmonic_roughness = med_num("harmonic_roughness", 6L),
    leaf_quality_weight = med_num("leaf_quality_weight", 6L)
  )
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

  if (!is.null(efa_df) && nrow(efa_df) > 0L) {
    if (is.null(sym_names) || length(sym_names) == 0L) {
      sym_names <- grep("^[AD][0-9]+$", names(efa_df), value = TRUE)
    }
    sym_names <- intersect(as.character(sym_names), names(efa_df))
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
    for (i in seq_len(n_spec)) {
      r_idx <- idx_by_spec[[spec_keys[i]]]
      if (length(r_idx) == 0L) next
      grp <- efa_df[r_idx, , drop = FALSE]
      if ("tier_assignment" %in% names(grp) && any(grp$tier_assignment == "Tier_1_Direct", na.rm = TRUE)) {
        grp <- grp[grp$tier_assignment == "Tier_1_Direct", , drop = FALSE]
      }
      if ("is_contour_outlier" %in% names(grp) && any(!as.logical(grp$is_contour_outlier), na.rm = TRUE)) {
        grp <- grp[!is.na(grp$is_contour_outlier) & !as.logical(grp$is_contour_outlier), , drop = FALSE]
      }
      lz <- if ("leaf_zone" %in% names(grp)) tolower(trimws(as.character(grp$leaf_zone))) else rep("basal", nrow(grp))
      lz[is.na(lz) | !(lz %in% c("basal", "cauline"))] <- "basal"
      basal_rows <- grp[lz == "basal", , drop = FALSE]
      cauline_rows <- grp[lz == "cauline", , drop = FALSE]
      n_bas <- nrow(basal_rows)
      n_caul <- nrow(cauline_rows)
      cauline_leaf_count[i] <- as.integer(n_caul)
      if ("basal_leaf_count" %in% names(specimen_df)) specimen_df$basal_leaf_count[i] <- as.integer(n_bas)
      if ("leaf_count" %in% names(specimen_df)) specimen_df$leaf_count[i] <- as.integer(n_bas + n_caul)
      if (n_bas == 0L && n_caul >= 1L && ("morph_tier_source" %in% names(specimen_df))) {
        specimen_df$morph_tier_source[i] <- "Cauline_Only_Fallback"
      }

      if (n_caul >= 1L) {
        cauline_solidity_med[i] <- med_col(cauline_rows, "solidity", 4L)
        cauline_dissection_index_med[i] <- med_col(cauline_rows, "dissection_index", 4L)
        cauline_aspect_ratio_med[i] <- med_col(cauline_rows, "aspect_ratio", 4L)
        c_bl <- med_col(cauline_rows, "blade_length_mm", 3L, pos_only = TRUE)
        c_cs <- med_col(cauline_rows, "centroid_size", 3L, pos_only = TRUE)
        cauline_blade_length_mm[i] <- if (is.finite(c_bl)) c_bl else c_cs

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
        if (is.finite(b_bl) && is.finite(c_bl) && c_bl > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_bl / c_bl, 4)
        } else if (is.finite(b_cs) && is.finite(c_cs) && c_cs > 1e-6) {
          basal_vs_cauline_length_ratio[i] <- round(b_cs / c_cs, 4)
        }

        b_sol <- if ("solidity" %in% names(specimen_df)) suppressWarnings(as.numeric(specimen_df$solidity[i])) else med_col(basal_rows, "solidity", 4L)
        if (is.finite(b_sol) && is.finite(cauline_solidity_med[i])) {
          basal_vs_cauline_solidity_delta[i] <- round(b_sol - cauline_solidity_med[i], 4)
        }

        if (length(sym_names) >= 1L) {
          caul_mat_i <- as.matrix(cauline_rows[, sym_names, drop = FALSE])
          storage.mode(caul_mat_i) <- "double"
          caul_med_h <- apply(caul_mat_i, 2L, stats::median, na.rm = TRUE)
          caul_med_h[!is.finite(caul_med_h)] <- 0.0
          caul_harm_list[[i]] <- caul_med_h
          if (n_bas >= 1L) {
            bas_mat_i <- as.matrix(basal_rows[, sym_names, drop = FALSE])
            storage.mode(bas_mat_i) <- "double"
            bas_med_h <- apply(bas_mat_i, 2L, stats::median, na.rm = TRUE)
            bas_med_h[!is.finite(bas_med_h)] <- 0.0
            d_efa <- sqrt(sum((bas_med_h - caul_med_h)^2, na.rm = TRUE))
            if (is.finite(d_efa)) {
              heteroblasty_efa_distance[i] <- round(d_efa, 6)
              if ("heteroblasty_index" %in% names(specimen_df)) {
                specimen_df$heteroblasty_index[i] <- round(d_efa, 6)
              }
            }
          }
        }
      }
    }

    caul_spec_idx <- which(cauline_leaf_count >= 1L & !vapply(caul_harm_list, is.null, logical(1)))
    if (length(caul_spec_idx) >= 1L && length(sym_names) >= 1L) {
      M_caul <- do.call(rbind, caul_harm_list[caul_spec_idx])
      storage.mode(M_caul) <- "double"
      for (jc in seq_len(ncol(M_caul))) {
        bad_c <- !is.finite(M_caul[, jc])
        if (any(bad_c)) M_caul[bad_c, jc] <- if (any(!bad_c)) stats::median(M_caul[!bad_c, jc], na.rm = TRUE) else 0.0
      }
      c_scores <- NULL
      if (nrow(M_caul) >= 3L) {
        c_vars <- apply(M_caul, 2L, stats::var, na.rm = TRUE)
        act_c <- which(is.finite(c_vars) & c_vars > 1e-12)
        if (length(act_c) >= 1L) {
          M_sub <- scale(M_caul[, act_c, drop = FALSE], center = TRUE, scale = FALSE)
          eig_c <- eigen(stats::cov(M_sub), symmetric = TRUE)
          c_scores <- M_sub %*% eig_c$vectors
        }
      }
      if (is.null(c_scores)) {
        all_sym_m <- as.matrix(efa_df[, sym_names, drop = FALSE])
        storage.mode(all_sym_m) <- "double"
        all_sym_m[!is.finite(all_sym_m)] <- 0.0
        if (nrow(all_sym_m) >= 2L) {
          c_vars <- apply(all_sym_m, 2L, stats::var, na.rm = TRUE)
          act_c <- which(is.finite(c_vars) & c_vars > 1e-12)
          if (length(act_c) >= 1L) {
            mu_all <- colMeans(all_sym_m[, act_c, drop = FALSE])
            eig_all <- eigen(stats::cov(all_sym_m[, act_c, drop = FALSE]), symmetric = TRUE)
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
    "cauline_solidity_med", "Cauline_PC1"
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
    } else if (raw_nm %in% c("heteroblasty_index", "heteroblasty_efa_distance")) {
      (is.finite(sd_v) && sd_v >= 0.015) && (sum(raw_v > 0) >= 5L)
    } else if (raw_nm %in% c("basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "Cauline_PC1")) {
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
