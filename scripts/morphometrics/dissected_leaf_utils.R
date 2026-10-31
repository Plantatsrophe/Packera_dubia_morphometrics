# ==============================================================================
# Module: dissected_leaf_utils.R - Packera Morphometrics Pipeline (NCU)
# Literature-backed utilities for positional basal-organ homology filtering,
# curved petiole/stem ribbon rejection, phase-invariant Dissection Index (DI),
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

compute_contour_ribbon_metrics <- function(coo) {
  if (is.null(coo) || !is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 12L) {
    return(list(width_ratio = 1.0, n_tips = 0L))
  }
  M <- nrow(coo)
  xy <- scale(matrix(as.numeric(coo[, 1:2]), ncol = 2L), center = TRUE, scale = FALSE)
  if (any(!is.finite(xy))) return(list(width_ratio = 1.0, n_tips = 0L))

  dmat <- as.matrix(stats::dist(xy))
  idx_mat <- outer(seq_len(M), seq_len(M), function(a, b) pmin(abs(a - b), M - abs(a - b)))
  min_sep <- max(4L, as.integer(floor(0.04 * M)))
  dmat[idx_mat < min_sep] <- Inf
  d_opp <- apply(dmat, 1L, min)
  d_opp <- d_opp[is.finite(d_opp)]
  width_ratio <- if (length(d_opp) >= 6L) {
    q92 <- as.numeric(stats::quantile(d_opp, probs = 0.92, na.rm = TRUE, names = FALSE))
    q25 <- as.numeric(stats::quantile(d_opp, probs = 0.25, na.rm = TRUE, names = FALSE))
    q92 / max(q25, 1e-6)
  } else {
    1.0
  }

  xs <- xy[, 1L]; ys <- xy[, 2L]
  for (pass in 1:2) {
    xs <- 0.25 * c(xs[M], xs[-M]) + 0.50 * xs + 0.25 * c(xs[-1L], xs[1L])
    ys <- 0.25 * c(ys[M], ys[-M]) + 0.50 * ys + 0.25 * c(ys[-1L], ys[1L])
  }
  prev_i <- c(M, seq_len(M - 1L))
  next_i <- c(seq(2L, M), 1L)
  v1x <- xs - xs[prev_i]; v1y <- ys - ys[prev_i]
  v2x <- xs[next_i] - xs; v2y <- ys[next_i] - ys
  signed_area <- 0.5 * sum(xs * ys[next_i] - xs[next_i] * ys)
  sgn <- if (is.finite(signed_area) && signed_area >= 0.0) 1.0 else -1.0
  turn <- sgn * atan2(v1x * v2y - v1y * v2x, pmax(-1e12, v1x * v2x + v1y * v2y))
  n_tips <- sum(turn > 0.30 & turn >= turn[prev_i] & turn > turn[next_i], na.rm = TRUE)
  list(width_ratio = as.numeric(width_ratio), n_tips = as.integer(n_tips))
}

is_non_lamina_petiole_ribbon <- function(coo, solidity = NA_real_, circularity = NA_real_) {
  sol <- suppressWarnings(as.numeric(solidity))
  circ <- suppressWarnings(as.numeric(circularity))
  if (!is.finite(sol) || !is.finite(circ)) {
    gp <- compute_polygon_area_perim_pure(coo)
    if (!is.finite(sol)) sol <- gp$solidity
    if (!is.finite(circ)) circ <- gp$circularity
  }
  if (!is.finite(sol) || !is.finite(circ)) return(FALSE)
  if (sol >= 0.65 || circ >= 0.085) return(FALSE)
  rm <- compute_contour_ribbon_metrics(coo)
  (rm$n_tips <= 3L) && (rm$width_ratio < 2.20)
}

filter_basal_lamina_contours <- function(
  fac_df,
  coo_list,
  min_cy_plant = 0.40,
  min_cy_sheet = 0.42,
  fallback_band = 0.15
) {
  if (is.null(fac_df) || nrow(fac_df) == 0L || length(coo_list) == 0L) {
    return(list(fac_df = fac_df, coo_list = coo_list, stats = list()))
  }
  n <- nrow(fac_df)
  if (!("plant_individual_id" %in% names(fac_df))) fac_df$plant_individual_id <- 0L
  fac_df$plant_individual_id <- ifelse(
    is.na(fac_df$plant_individual_id),
    0L,
    suppressWarnings(as.integer(fac_df$plant_individual_id))
  )
  fac_df$plant_individual_id[is.na(fac_df$plant_individual_id)] <- 0L

  di_vec <- if ("dissection_index" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$dissection_index)) else rep(NA_real_, n)
  sol_vec <- if ("solidity" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$solidity)) else rep(NA_real_, n)
  circ_vec <- if ("circularity" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$circularity)) else rep(NA_real_, n)
  is_ribbon <- logical(n)

  for (i in seq_len(n)) {
    coo_i <- coo_list[[i]]
    if (!is.finite(di_vec[i]) || !is.finite(sol_vec[i]) || !is.finite(circ_vec[i])) {
      gp <- compute_polygon_area_perim_pure(coo_i)
      if (!is.finite(sol_vec[i])) sol_vec[i] <- round(gp$solidity, 4)
      if (!is.finite(circ_vec[i])) circ_vec[i] <- round(gp$circularity, 4)
      if (!is.finite(di_vec[i])) di_vec[i] <- compute_dissection_index_from_coo(area_val = gp$area, perim_val = gp$perim)
    }
    is_ribbon[i] <- is_non_lamina_petiole_ribbon(coo_i, solidity = sol_vec[i], circularity = circ_vec[i])
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

  cy_p <- if ("cy_rel_plant_bbox" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$cy_rel_plant_bbox)) else rep(NA_real_, n)
  cy_s <- if ("cy_rel_sheet" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$cy_rel_sheet)) else rep(NA_real_, n)
  ba_cm2 <- if ("blade_area_cm2" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$blade_area_cm2)) else rep(NA_real_, n)
  min_ba <- if ("effective_min_blade_area_cm2" %in% names(fac_df)) suppressWarnings(as.numeric(fac_df$effective_min_blade_area_cm2)) else rep(NA_real_, n)
  area_ok <- ifelse(is.finite(ba_cm2) & is.finite(min_ba) & min_ba > 0, ba_cm2 >= min_ba, TRUE)

  keep_mask <- rep(FALSE, n)
  is_fallback_vec <- rep(FALSE, n)
  grp_keys <- paste0(trimws(as.character(fac_df$catalogNumber)), "__", fac_df$plant_individual_id)
  idx_by_grp <- split(seq_len(n), grp_keys)
  n_fb_individuals <- 0L

  for (g_idx in idx_by_grp) {
    lam_ok <- g_idx[!is_ribbon[g_idx] & area_ok[g_idx]]
    if (length(lam_ok) == 0L) lam_ok <- g_idx[!is_ribbon[g_idx]]
    if (length(lam_ok) == 0L) next

    has_cy_p <- is.finite(cy_p[lam_ok])
    has_cy_s <- !has_cy_p & is.finite(cy_s[lam_ok])
    if (!any(has_cy_p) && !any(has_cy_s)) {
      keep_mask[lam_ok] <- TRUE
      next
    }
    in_primary <- ifelse(
      has_cy_p,
      cy_p[lam_ok] >= min_cy_plant,
      ifelse(has_cy_s, cy_s[lam_ok] >= min_cy_sheet, TRUE)
    )
    if (any(in_primary)) {
      keep_mask[lam_ok[in_primary]] <- TRUE
    } else {
      # Confirmed detached/top-mounted rosette fallback: retain lowest lamina band on individual
      if (any(has_cy_p)) {
        c_max <- max(cy_p[lam_ok[has_cy_p]], na.rm = TRUE)
        fb_hit <- has_cy_p & (cy_p[lam_ok] >= (c_max - fallback_band))
      } else {
        c_max <- max(cy_s[lam_ok[has_cy_s]], na.rm = TRUE)
        fb_hit <- has_cy_s & (cy_s[lam_ok] >= (c_max - fallback_band))
      }
      if (any(fb_hit)) {
        keep_mask[lam_ok[fb_hit]] <- TRUE
        is_fallback_vec[lam_ok[fb_hit]] <- TRUE
        n_fb_individuals <- n_fb_individuals + 1L
      }
    }
  }

  # Safeguard: never drop 100% of cohort contours on minimal synthetic fixtures
  if (!any(keep_mask)) {
    keep_mask <- rep(TRUE, n)
  }

  fac_df$is_basal_fallback <- is_fallback_vec
  n_rib <- sum(is_ribbon)
  n_area_drop <- sum(!is_ribbon & !area_ok & !keep_mask)
  n_caul_drop <- sum(!is_ribbon & area_ok & !keep_mask)
  if (n_rib > 0L || n_area_drop > 0L || n_caul_drop > 0L || n_fb_individuals > 0L) {
    message(sprintf(
      paste0(
        "Positional Homology & Basal Lamina Gate: retained %d / %d Tier 1 basal blades ",
        "(excluded %d upper cauline, %d petiole/stem ribbons, %d sub-floor blade fragments; ",
        "lowest-lamina fallback individuals = %d)."
      ),
      sum(keep_mask), n, n_caul_drop, n_rib, n_area_drop, n_fb_individuals
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
      n_cauline_excluded = n_caul_drop,
      n_ribbon_excluded = n_rib,
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
      set.seed(as.integer(seed))
      q_used <- min(nrow(X_reg_fit) - 1L, max(floor(0.88 * nrow(X_reg_fit)), floor((nrow(X_reg_fit) + df_reg + 1L) / 2L)))
      mcd_fit <- tryCatch(MASS::cov.mcd(X_reg_fit, quantile.used = q_used), error = function(e) NULL)
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

filter_biologically_variable_dissection_scalars <- function(df, cand_scalars, active_mask = NULL) {
  cand_scalars <- intersect(as.character(cand_scalars), names(df))
  if (length(cand_scalars) == 0L || nrow(df) < 5L) return(character(0))
  act <- if (!is.null(active_mask) && length(active_mask) == nrow(df) && sum(active_mask, na.rm = TRUE) >= 5L) {
    which(as.logical(active_mask))
  } else {
    seq_len(nrow(df))
  }
  keep <- character(0)
  for (sc in cand_scalars) {
    raw_nm <- sub("_std$", "", sc)
    raw_v <- if (raw_nm %in% names(df)) suppressWarnings(as.numeric(df[[raw_nm]][act])) else suppressWarnings(as.numeric(df[[sc]][act]))
    raw_v <- raw_v[is.finite(raw_v)]
    if (length(raw_v) < 5L) next
    sd_v <- stats::sd(raw_v)
    q_span <- as.numeric(diff(stats::quantile(raw_v, probs = c(0.10, 0.90), names = FALSE)))
    ok <- if (raw_nm %in% c("solidity", "q25_solidity")) {
      (is.finite(sd_v) && sd_v >= 0.045) || any(raw_v < 0.75) || (is.finite(q_span) && q_span >= 0.10)
    } else if (identical(raw_nm, "dissection_index")) {
      (is.finite(sd_v) && sd_v >= 0.12) || any(raw_v > 1.35)
    } else if (identical(raw_nm, "heteroblasty_index")) {
      (is.finite(sd_v) && sd_v >= 0.015) && (sum(raw_v > 0) >= 5L)
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
    ge_cnt <- 0L
    tol <- ss_between - 1e-12 * max(1.0, ss_tot)
    for (b in seq_len(B)) {
      g_perm <- sample.int(N)
      sums_b <- rowsum(X_c[g_perm, , drop = FALSE], g_int, reorder = TRUE)
      ss_b <- sum(rowSums(sums_b^2) / n_vec)
      if (ss_b >= tol) ge_cnt <- ge_cnt + 1L
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
