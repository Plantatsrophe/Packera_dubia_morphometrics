#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/plasticity_mitigation_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Strictly modality-independent phenotypic plasticity mitigation utilities for
#   Step 04 unsupervised GMM clustering and Canonical Discriminant Analysis:
#     1. Specimen-intrinsic capitulum phenophase covariate ingestion from Step 02
#        (`voucher_phenological_states.csv`) and within-taxon ontogenetic/fruiting
#        residualization free of geographic coordinates, collection dates, or
#        ecoregion labels (preserving Padial et al. 2010 modality independence).
#     2. Ray & Lindsay (2005) / Hennig (2010) multivariate Gaussian ridgeline
#        density unimodality diagnosis & skewed-tail component merging in pure
#        morphological space, allowing genuine G = 1 unimodal species nulls while
#        preserving true multimodal morphological discontinuities.
# ==============================================================================

merge_voucher_phenophase_covariates <- function(df, harmonics_path = NULL, vouchers_path = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  n <- nrow(df)

  cand_dirs <- unique(stats::na.omit(c(
    if (!is.null(harmonics_path) && nzchar(harmonics_path)) dirname(harmonics_path) else NA_character_,
    if (!is.null(vouchers_path) && nzchar(vouchers_path)) dirname(vouchers_path) else NA_character_
  )))
  for (d_cand in cand_dirs) {
    ph_csv <- file.path(d_cand, "voucher_phenological_states.csv")
    if (file.exists(ph_csv) && isTRUE(file.info(ph_csv)$size > 0)) {
      ph_df <- tryCatch(utils::read.csv(ph_csv, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(ph_df) && nrow(ph_df) > 0L && "catalogNumber" %in% names(ph_df)) {
        ph_pid <- if ("plant_individual_id" %in% names(ph_df)) {
          ifelse(is.na(ph_df$plant_individual_id), 0L, as.integer(ph_df$plant_individual_id))
        } else {
          rep(0L, nrow(ph_df))
        }
        df_pid <- if ("plant_individual_id" %in% names(df)) {
          ifelse(is.na(df$plant_individual_id), 0L, as.integer(df$plant_individual_id))
        } else {
          rep(0L, n)
        }
        k_ph <- paste(trimws(as.character(ph_df$catalogNumber)), ph_pid, sep = "__")
        k_df <- paste(trimws(as.character(df$catalogNumber)), df_pid, sep = "__")
        n_pids_sheet <- tapply(ph_pid, trimws(as.character(ph_df$catalogNumber)), function(z) length(unique(z)))
        single_pid_sheet <- !(trimws(as.character(df$catalogNumber)) %in% names(n_pids_sheet[n_pids_sheet > 1L]))
        m_exact <- match(k_df, k_ph)
        m_sheet <- match(paste(trimws(as.character(df$catalogNumber)), 0L, sep = "__"), k_ph)
        m_use <- ifelse(!is.na(m_exact), m_exact, ifelse(single_pid_sheet, m_sheet, NA_integer_))
        ph_cols <- intersect(
          c("n_anthesis", "n_bud", "n_fruiting", "n_senescent", "total_capitula", "phenological_state", "is_verified_anthesis"),
          names(ph_df)
        )
        for (pc_nm in ph_cols) {
          if (!(pc_nm %in% names(df))) {
            df[[pc_nm]] <- ph_df[[pc_nm]][m_use]
          } else {
            miss_m <- is.na(df[[pc_nm]]) | !nzchar(trimws(as.character(df[[pc_nm]])))
            if (any(miss_m)) df[[pc_nm]][miss_m] <- ph_df[[pc_nm]][m_use[miss_m]]
          }
        }
      }
      break
    }
  }
  df
}

residualize_ontogenetic_phenophase_plasticity <- function(
    df,
    target_cols,
    ref_mask = NULL,
    intraspecific_mode = FALSE,
    min_partial_r2 = 0.04
) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 12L) return(df)
  n <- nrow(df)
  target_cols <- intersect(as.character(target_cols), names(df))
  if (length(target_cols) == 0L) return(df)

  auth_w <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n)
  auth_w[!is.finite(auth_w) | auth_w < 0.0] <- 0.0
  sb_w <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, n)
  sb_w[!is.finite(sb_w) | sb_w <= 0.0] <- 1.0
  eff_w <- auth_w * sb_w

  if (is.null(ref_mask) || length(ref_mask) != n) {
    is_t1 <- if ("morph_tier_source" %in% names(df)) (!is.na(df$morph_tier_source) & df$morph_tier_source == "Basal_Rosette") else rep(TRUE, n)
    ref_mask <- (auth_w > 0.0) & is_t1
  }
  ref_idx <- which(as.logical(ref_mask))
  if (length(ref_idx) < 12L) return(df)

  # Strictly specimen-intrinsic developmental covariates (zero geographic coordinates, DOY, or ecoregion labels)
  cov_list <- list()

  # 1. Post-anthesis fruiting/senescent capitulum fraction from Step 02 single-individual reproductive detections
  n_fr <- rep(0.0, n)
  has_fr <- FALSE
  for (fc in intersect(c("n_fruiting", "n_senescent"), names(df))) {
    vv <- suppressWarnings(as.numeric(df[[fc]]))
    if (any(is.finite(vv) & vv > 0.0)) {
      has_fr <- TRUE
      n_fr <- n_fr + ifelse(is.finite(vv) & vv > 0.0, vv, 0.0)
    }
  }
  tot_cap <- if ("total_capitula" %in% names(df)) suppressWarnings(as.numeric(df$total_capitula)) else if ("capitula_count" %in% names(df)) suppressWarnings(as.numeric(df$capitula_count)) else rep(NA_real_, n)
  if (has_fr) {
    has_cap_obs <- is.finite(tot_cap) & (tot_cap >= 1.0)
    fr_frac <- ifelse(has_cap_obs, pmin(1.0, pmax(0.0, n_fr / pmax(tot_cap, n_fr))), NA_real_)
    ok_fr_ref <- ref_idx[is.finite(fr_frac[ref_idx])]
    if (length(ok_fr_ref) >= 12L && stats::sd(fr_frac[ok_fr_ref], na.rm = TRUE) > 0.05) {
      cov_list[["fruiting_frac"]] <- fr_frac
    }
  }

  # 2. Ruler-calibrated foliar size covariate for non-allometry-corrected architectural traits
  has_log_cs <- FALSE
  log_cs_v <- rep(NA_real_, n)
  if ("log_centsize" %in% names(df)) {
    cal_ok <- if ("is_calibrated_mm" %in% names(df)) (!is.na(df$is_calibrated_mm) & as.logical(df$is_calibrated_mm)) else rep(TRUE, n)
    vv <- ifelse(cal_ok, suppressWarnings(as.numeric(df$log_centsize)), NA_real_)
    ok_cs <- ref_idx[is.finite(vv[ref_idx])]
    if (length(ok_cs) >= 12L && stats::sd(vv[ok_cs]) > 0.05) {
      med_cs <- stats::median(vv[ok_cs])
      log_cs_v <- ifelse(is.finite(vv), vv - med_cs, NA_real_)
      has_log_cs <- TRUE
    }
  }

  if (length(cov_list) == 0L && !has_log_cs) return(df)

  sp_raw <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Packera_sp", n)
  sp_raw[is.na(sp_raw) | !nzchar(trimws(sp_raw))] <- "Packera_sp"

  for (col_nm in target_cols) {
    yv <- suppressWarnings(as.numeric(df[[col_nm]]))
    if (sum(is.finite(yv[ref_idx])) < 12L) next
    is_pc_col <- grepl("^PC[0-9]+$", col_nm)
    if (is_pc_col && !isTRUE(intraspecific_mode)) next

    use_covs <- cov_list
    if (!is_pc_col && has_log_cs && !isTRUE(attr(df, "allometry_corrected_scalars"))) {
      use_covs[["log_centsize_dev"]] <- log_cs_v
    }
    if (length(use_covs) == 0L) next

    C_full <- do.call(cbind, use_covs)
    fin_sub <- is.finite(yv[ref_idx]) & stats::complete.cases(C_full[ref_idx, , drop = FALSE])
    if (sum(fin_sub) < 12L) next

    sub_rows <- ref_idx[fin_sub]
    sub_sp <- droplevels(factor(sp_raw[sub_rows]))
    sp_cnts <- table(sub_sp)
    keep_sp <- names(sp_cnts[sp_cnts >= 3L])
    ok_sp <- as.character(sub_sp) %in% keep_sp
    if (sum(ok_sp) < 12L) next

    sub_rows <- sub_rows[ok_sp]
    sub_sp <- droplevels(sub_sp[ok_sp])
    w_sub <- eff_w[sub_rows]
    if (sum(w_sub) <= 1e-12) w_sub <- rep(1.0, length(sub_rows))

    w_sp_mean <- function(v_in) {
      unname(vapply(split(seq_along(v_in), sub_sp), function(ii) stats::weighted.mean(v_in[ii], w = w_sub[ii], na.rm = TRUE), numeric(1))[as.character(sub_sp)])
    }
    y_dm <- as.numeric(yv[sub_rows] - w_sp_mean(yv[sub_rows]))
    ss_y <- sum(w_sub * (y_dm^2))
    if (!is.finite(ss_y) || ss_y <= 1e-8) next

    C_dm <- apply(C_full[sub_rows, , drop = FALSE], 2L, function(cv) as.numeric(cv - w_sp_mean(cv)))
    if (!is.matrix(C_dm)) C_dm <- matrix(C_dm, ncol = ncol(C_full))
    keep_c <- which(apply(C_dm, 2L, stats::sd) > 1e-5)
    if (length(keep_c) == 0L) next
    C_dm <- C_dm[, keep_c, drop = FALSE]
    C_eval <- C_full[, keep_c, drop = FALSE]

    fit_w <- tryCatch(stats::lm.wfit(x = C_dm, y = y_dm, w = w_sub), error = function(e) NULL)
    if (is.null(fit_w)) next
    beta_w <- fit_w$coefficients
    beta_w[!is.finite(beta_w)] <- 0.0
    rss_w <- sum(w_sub * (fit_w$residuals^2))
    r2_w <- max(0.0, 1.0 - (rss_w / ss_y))
    df1 <- max(1L, fit_w$rank)
    df2 <- length(sub_rows) - nlevels(sub_sp) - df1
    if (df2 < 5L) next
    f_stat <- ((ss_y - rss_w) / df1) / max(rss_w / df2, 1e-12)
    p_val <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)

    if (is.finite(r2_w) && r2_w >= min_partial_r2 && is.finite(p_val) && p_val < 0.05) {
      c_ref_means <- vapply(seq_len(ncol(C_eval)), function(j) {
        ok_j <- is.finite(C_eval[sub_rows, j])
        if (any(ok_j)) stats::weighted.mean(C_eval[sub_rows[ok_j], j], w = w_sub[ok_j]) else 0.0
      }, numeric(1))
      fin_all <- which(is.finite(yv) & stats::complete.cases(C_eval))
      if (length(fin_all) > 0L) {
        C_cent <- sweep(C_eval[fin_all, , drop = FALSE], 2L, c_ref_means, "-")
        adj_vals <- yv[fin_all] - as.numeric(C_cent %*% beta_w)
        df[[col_nm]][fin_all] <- round(adj_vals, 6)
      }
    }
  }
  df
}

.compute_ridgeline_valley_ratio <- function(pi_j, mu_j, Sig_j, pi_l, mu_l, Sig_l, n_grid = 101L) {
  p <- length(mu_j)
  reg_j <- Sig_j + diag(1e-5, nrow = p, ncol = p)
  reg_l <- Sig_l + diag(1e-5, nrow = p, ncol = p)
  inv_j <- tryCatch(solve(reg_j), error = function(e) diag(x = 1.0 / pmax(diag(as.matrix(reg_j)), 1e-5), nrow = p, ncol = p))
  inv_l <- tryCatch(solve(reg_l), error = function(e) diag(x = 1.0 / pmax(diag(as.matrix(reg_l)), 1e-5), nrow = p, ncol = p))
  ldet_j <- as.numeric(determinant(reg_j, logarithm = TRUE)$modulus)
  ldet_l <- as.numeric(determinant(reg_l, logarithm = TRUE)$modulus)

  w_j <- pi_j / max(pi_j + pi_l, 1e-12)
  w_l <- 1.0 - w_j
  alphas <- seq(0.0, 1.0, length.out = n_grid)
  dens_v <- numeric(n_grid)

  for (idx_a in seq_along(alphas)) {
    al <- alphas[idx_a]
    A_al <- (1.0 - al) * inv_j + al * inv_l
    rhs <- (1.0 - al) * (inv_j %*% mu_j) + al * (inv_l %*% mu_l)
    x_al <- tryCatch(as.numeric(solve(A_al, rhs)), error = function(e) (1.0 - al) * mu_j + al * mu_l)
    dj <- x_al - mu_j
    dl <- x_al - mu_l
    log_fj <- -0.5 * (p * log(2.0 * pi) + ldet_j + sum(dj * as.numeric(inv_j %*% dj)))
    log_fl <- -0.5 * (p * log(2.0 * pi) + ldet_l + sum(dl * as.numeric(inv_l %*% dl)))
    m_log <- max(log_fj, log_fl)
    dens_v[idx_a] <- exp(m_log) * (w_j * exp(log_fj - m_log) + w_l * exp(log_fl - m_log))
  }

  is_peak <- logical(n_grid)
  is_peak[1L] <- dens_v[1L] >= dens_v[2L]
  is_peak[n_grid] <- dens_v[n_grid] >= dens_v[n_grid - 1L]
  if (n_grid >= 3L) {
    for (ii in 2L:(n_grid - 1L)) {
      is_peak[ii] <- (dens_v[ii] >= dens_v[ii - 1L]) && (dens_v[ii] >= dens_v[ii + 1L]) &&
        ((dens_v[ii] > dens_v[ii - 1L]) || (dens_v[ii] > dens_v[ii + 1L]))
    }
  }
  peak_idx <- which(is_peak)
  if (length(peak_idx) <= 1L) {
    return(list(n_modes = 1L, ridge_ratio = 1.0))
  }
  p_left <- peak_idx[1L]
  p_right <- peak_idx[length(peak_idx)]
  if (p_right - p_left <= 1L) {
    return(list(n_modes = 1L, ridge_ratio = 1.0))
  }
  f_saddle <- min(dens_v[p_left:p_right])
  f_small_peak <- min(dens_v[p_left], dens_v[p_right])
  ratio <- if (f_small_peak > 1e-300) min(1.0, f_saddle / f_small_peak) else 1.0
  list(n_modes = length(peak_idx), ridge_ratio = as.numeric(ratio))
}

merge_plastic_gmm_components_ridgeline <- function(
    X_eval,
    X_core,
    probs_eval,
    probs_core,
    intraspecific_mode = FALSE,
    ridge_intraspecific_max = 0.45,
    ridge_interspecific_max = 0.90
) {
  if (is.null(probs_eval) || !is.matrix(probs_eval) || ncol(probs_eval) < 2L) {
    return(list(
      best_k = 1L,
      classification = rep(1L, nrow(X_eval)),
      uncertainty = rep(0.0, nrow(X_eval)),
      probs_eval = matrix(1.0, nrow = nrow(X_eval), ncol = 1L),
      merged_any = FALSE
    ))
  }

  Z_core <- as.matrix(probs_core)
  Z_eval <- as.matrix(probs_eval)
  p <- ncol(X_core)
  merged_any <- FALSE
  ridge_thr <- if (isTRUE(intraspecific_mode)) as.numeric(ridge_intraspecific_max) else as.numeric(ridge_interspecific_max)

  repeat {
    G_cur <- ncol(Z_core)
    if (G_cur < 2L) break

    nk <- colSums(Z_core, na.rm = TRUE)
    valid_k <- which(nk > 1e-6)
    if (length(valid_k) < G_cur) {
      if (length(valid_k) == 0L) valid_k <- 1L
      Z_core <- Z_core[, valid_k, drop = FALSE]
      Z_eval <- Z_eval[, valid_k, drop = FALSE]
      Z_core <- sweep(Z_core, 1L, pmax(rowSums(Z_core), 1e-12), "/")
      Z_eval <- sweep(Z_eval, 1L, pmax(rowSums(Z_eval), 1e-12), "/")
      merged_any <- TRUE
      next
    }

    pis <- nk / sum(nk)
    mus <- lapply(seq_len(G_cur), function(k) {
      wk <- Z_core[, k]
      colSums(X_core * wk) / max(sum(wk), 1e-12)
    })
    Sigs <- lapply(seq_len(G_cur), function(k) {
      wk <- Z_core[, k]
      Xc_k <- sweep(X_core, 2L, mus[[k]], "-")
      Sk <- crossprod(Xc_k * sqrt(wk)) / max(sum(wk), 1.0)
      Sk + diag(0.02 * mean(diag(as.matrix(Sk))) + 1e-4, nrow = p, ncol = p)
    })

    best_pair <- NULL
    best_ratio <- -1.0

    for (j in seq_len(G_cur - 1L)) {
      for (l in (j + 1L):G_cur) {
        rl <- .compute_ridgeline_valley_ratio(pis[j], mus[[j]], Sigs[[j]], pis[l], mus[[l]], Sigs[[l]])
        should_merge <- (rl$n_modes <= 1L) || (rl$ridge_ratio >= ridge_thr)
        if (should_merge && rl$ridge_ratio > best_ratio) {
          best_ratio <- rl$ridge_ratio
          best_pair <- c(j, l)
        }
      }
    }

    if (is.null(best_pair)) break

    j_keep <- best_pair[1L]
    l_drop <- best_pair[2L]
    Z_core[, j_keep] <- Z_core[, j_keep] + Z_core[, l_drop]
    Z_core <- Z_core[, -l_drop, drop = FALSE]
    Z_eval[, j_keep] <- Z_eval[, j_keep] + Z_eval[, l_drop]
    Z_eval <- Z_eval[, -l_drop, drop = FALSE]
    merged_any <- TRUE
  }

  best_k <- ncol(Z_eval)
  cls_out <- if (best_k == 1L) rep(1L, nrow(Z_eval)) else as.integer(max.col(Z_eval, ties.method = "first"))
  unc_out <- if (best_k == 1L) rep(0.0, nrow(Z_eval)) else round(pmax(0.0, 1.0 - apply(Z_eval, 1L, max)), 6)
  list(
    best_k = as.integer(best_k),
    classification = cls_out,
    uncertainty = unc_out,
    probs_eval = Z_eval,
    merged_any = merged_any
  )
}
