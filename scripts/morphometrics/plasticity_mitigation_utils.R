#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/plasticity_mitigation_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Non-circular phenotypic plasticity mitigation utilities for Step 04
#   unsupervised GMM clustering and Canonical Discriminant Analysis:
#     1. Phenophase & seasonal heteroblasty covariate ingestion and within-deme
#        ontogenetic/seasonal residualization (Collyer & Adams 2018; Klingenberg 2016).
#     2. Within-population (conspecific, same-ecoregion <= 15 km deme) plasticity
#        covariance whitening (Q_ST-proxy eigen-shrinkage) on active anchors.
#     3. Spatial Principal Component Analysis (sPCA; Jombart et al. 2008) global
#        Moran's I permutation test & cluster spatial segregation test.
#     4. Ray & Lindsay (2005) / Hennig (2010) multivariate Gaussian ridgeline
#        unimodality diagnosis & sympatric-plasticity component merging, allowing
#        genuine G = 1 unimodal species nulls while preserving geographically
#        segregated ecotypes and deeply bimodal sympatric cryptic lineages.
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
        m_exact <- match(k_df, k_ph)
        # Fallback to sheet-level pid = 0 when individual-level entry is absent
        m_sheet <- match(paste(trimws(as.character(df$catalogNumber)), 0L, sep = "__"), k_ph)
        m_use <- ifelse(!is.na(m_exact), m_exact, m_sheet)
        ph_cols <- intersect(
          c("n_anthesis", "n_bud", "n_fruiting", "n_fruit", "n_senescent", "total_capitula", "phenological_state", "is_verified_anthesis"),
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

.build_ecological_deme_ids <- function(df, idx_vec, deme_radius_km = 15.0) {
  m <- length(idx_vec)
  if (m == 0L) return(character(0))
  sub_df <- df[idx_vec, , drop = FALSE]

  sp_v <- if ("species_standardized" %in% names(sub_df)) as.character(sub_df$species_standardized) else rep("Packera_sp", m)
  sp_v[is.na(sp_v) | !nzchar(trimws(sp_v))] <- "Packera_sp"
  rg_v <- if ("regional_group" %in% names(sub_df)) trimws(as.character(sub_df$regional_group)) else rep("", m)
  rg_v[is.na(rg_v) | tolower(rg_v) %in% c("na", "unknown", "none")] <- ""

  lat_v <- rep(NA_real_, m)
  lon_v <- rep(NA_real_, m)
  for (lc in c("decimalLatitude", "latitude")) {
    if (lc %in% names(sub_df)) {
      vv <- suppressWarnings(as.numeric(sub_df[[lc]]))
      ok <- !is.finite(lat_v) & is.finite(vv)
      lat_v[ok] <- vv[ok]
    }
  }
  for (lc in c("decimalLongitude", "longitude")) {
    if (lc %in% names(sub_df)) {
      vv <- suppressWarnings(as.numeric(sub_df[[lc]]))
      ok <- !is.finite(lon_v) & is.finite(vv)
      lon_v[ok] <- vv[ok]
    }
  }
  has_c <- is.finite(lat_v) & is.finite(lon_v) & (abs(lat_v) <= 90.0) & (abs(lon_v) <= 180.0)

  out_ids <- sprintf("ISOLATE_%d", seq_len(m))
  grp_keys <- paste(sp_v, rg_v, sep = "||")
  for (gk in unique(grp_keys[has_c])) {
    loc <- which(has_c & grp_keys == gk)
    ng <- length(loc)
    if (ng == 1L) {
      out_ids[loc] <- paste0(gk, "__D1")
      next
    }
    lat_r <- lat_v[loc] * (pi / 180.0)
    lon_r <- lon_v[loc] * (pi / 180.0)
    comp <- seq_len(ng)
    for (a in seq_len(ng - 1L)) {
      b_seq <- (a + 1L):ng
      dlat <- lat_r[b_seq] - lat_r[a]
      dlon <- lon_r[b_seq] - lon_r[a]
      hav <- sin(dlat * 0.5)^2 + cos(lat_r[a]) * cos(lat_r[b_seq]) * (sin(dlon * 0.5)^2)
      dkm <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hav))))
      hit <- b_seq[is.finite(dkm) & (dkm <= deme_radius_km)]
      if (length(hit) > 0L) {
        old_c <- unique(c(comp[a], comp[hit]))
        comp[comp %in% old_c] <- min(old_c)
      }
    }
    out_ids[loc] <- paste0(gk, "__D", match(comp, sort(unique(comp))))
  }
  out_ids
}

residualize_within_deme_phenophase_plasticity <- function(
    df,
    target_cols,
    ref_mask = NULL,
    intraspecific_mode = FALSE,
    deme_radius_km = 25.0,
    min_partial_r2 = 0.04
) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 12L) return(df)
  n <- nrow(df)
  target_cols <- intersect(as.character(target_cols), names(df))
  if (length(target_cols) == 0L) return(df)

  if (is.null(ref_mask) || length(ref_mask) != n) {
    auth_w <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n)
    is_t1 <- if ("morph_tier_source" %in% names(df)) (!is.na(df$morph_tier_source) & df$morph_tier_source == "Tier_1_Pure") else rep(TRUE, n)
    ref_mask <- is.finite(auth_w) & (auth_w > 0.0) & is_t1
  }
  ref_idx <- which(as.logical(ref_mask))
  if (length(ref_idx) < 12L) return(df)

  # Build specimen-level developmental/seasonal plasticity covariates
  cov_list <- list()

  # 1. Fruiting / post-anthesis capitulum fraction (captures pappus flaring in AR_inv & leaf senescence)
  n_fr <- rep(0.0, n)
  has_fr <- FALSE
  for (fc in intersect(c("n_fruiting", "n_fruit", "n_senescent"), names(df))) {
    vv <- suppressWarnings(as.numeric(df[[fc]]))
    if (any(is.finite(vv) & vv > 0.0)) {
      has_fr <- TRUE
      n_fr <- n_fr + ifelse(is.finite(vv) & vv > 0.0, vv, 0.0)
    }
  }
  tot_cap <- if ("total_capitula" %in% names(df)) suppressWarnings(as.numeric(df$total_capitula)) else if ("capitula_count" %in% names(df)) suppressWarnings(as.numeric(df$capitula_count)) else rep(NA_real_, n)
  if (has_fr) {
    denom_c <- ifelse(is.finite(tot_cap) & tot_cap >= 1.0, pmax(tot_cap, n_fr), pmax(1.0, n_fr))
    fr_frac <- pmin(1.0, pmax(0.0, n_fr / denom_c))
    if (stats::sd(fr_frac[ref_idx], na.rm = TRUE) > 0.05) {
      cov_list[["fruiting_frac"]] <- fr_frac
    }
  }

  # 2. Hopkins-detrended seasonal collection offset & late-summer bolting/senescence hinge (> +30 days)
  if ("doy" %in% names(df)) {
    doy_v <- suppressWarnings(as.numeric(df$doy))
    lat_c <- intersect(c("decimalLatitude", "latitude"), names(df))
    lat_v <- if (length(lat_c) > 0L) suppressWarnings(as.numeric(df[[lat_c[1L]]])) else rep(NA_real_, n)
    ok_doy_ref <- ref_idx[is.finite(doy_v[ref_idx]) & doy_v[ref_idx] >= 1.0 & doy_v[ref_idx] <= 366.0]
    if (length(ok_doy_ref) >= 12L) {
      med_doy <- stats::median(doy_v[ok_doy_ref], na.rm = TRUE)
      ok_lat_ref <- ok_doy_ref[is.finite(lat_v[ok_doy_ref])]
      med_lat <- if (length(ok_lat_ref) >= 5L) stats::median(lat_v[ok_lat_ref], na.rm = TRUE) else 35.0
      exp_doy <- med_doy + ifelse(is.finite(lat_v), 4.0 * (lat_v - med_lat), 0.0)
      hop_off <- ifelse(is.finite(doy_v) & doy_v >= 1.0 & doy_v <= 366.0, doy_v - exp_doy, 0.0)
      # Only activate seasonal plasticity covariate when genuine wide phenological spread exists (SD > 12 days),
      # preventing overfitting on narrow uniform-spring or synthetic benchmark cohorts
      if (stats::sd(hop_off[ok_doy_ref], na.rm = TRUE) > 12.0) {
        cov_list[["hopkins_offset"]] <- hop_off / 30.0
        late_hinge <- pmax(0.0, hop_off - 30.0) / 30.0
        if (sum(late_hinge[ok_doy_ref] > 0.0) >= 3L && stats::sd(late_hinge[ok_doy_ref]) > 0.05) {
          cov_list[["late_season_hinge"]] <- late_hinge
        }
      }
    }
  }

  # 3. Foliar size covariate for scalar architectural traits (which test_and_correct_allometry skips)
  has_log_cs <- FALSE
  log_cs_v <- rep(0.0, n)
  if ("log_centsize" %in% names(df)) {
    vv <- suppressWarnings(as.numeric(df$log_centsize))
    ok_cs <- ref_idx[is.finite(vv[ref_idx])]
    if (length(ok_cs) >= 12L && stats::sd(vv[ok_cs]) > 0.05) {
      med_cs <- stats::median(vv[ok_cs])
      log_cs_v <- ifelse(is.finite(vv), vv - med_cs, 0.0)
      has_log_cs <- TRUE
    }
  }

  if (length(cov_list) == 0L && !has_log_cs) return(df)

  deme_ref <- .build_ecological_deme_ids(df, ref_idx, deme_radius_km = deme_radius_km)
  deme_tab <- table(deme_ref)
  multi_demes <- names(deme_tab[deme_tab >= 2L])
  in_multi <- deme_ref %in% multi_demes
  df_within <- sum(deme_tab[multi_demes] - 1L)
  if (df_within < 10L) return(df)

  ref_multi_idx <- ref_idx[in_multi]
  deme_f <- factor(deme_ref[in_multi])

  center_by_deme <- function(vec_sub) {
    d_means <- tapply(vec_sub, deme_f, mean, na.rm = TRUE)
    as.numeric(vec_sub - unname(d_means[as.character(deme_f)]))
  }

  for (col_nm in target_cols) {
    yv <- suppressWarnings(as.numeric(df[[col_nm]]))
    if (sum(is.finite(yv[ref_multi_idx])) < 12L) next
    is_pc_col <- grepl("^(Sym_)?PC[0-9]+$", col_nm)
    if (is_pc_col && !isTRUE(intraspecific_mode)) next

    use_covs <- cov_list
    if (!is_pc_col && has_log_cs) {
      use_covs[["log_centsize_dev"]] <- log_cs_v
    }
    if (length(use_covs) == 0L) next

    C_full <- do.call(cbind, use_covs)
    fin_sub <- is.finite(yv[ref_multi_idx]) & stats::complete.cases(C_full[ref_multi_idx, , drop = FALSE])
    if (sum(fin_sub) < 12L) next

    sub_rows <- ref_multi_idx[fin_sub]
    sub_deme <- droplevels(factor(deme_ref[in_multi][fin_sub]))
    d_cnts <- table(sub_deme)
    keep_d <- names(d_cnts[d_cnts >= 2L])
    ok_d <- as.character(sub_deme) %in% keep_d
    if (sum(ok_d) < 10L || sum(d_cnts[keep_d] - 1L) < 8L) next

    sub_rows <- sub_rows[ok_d]
    sub_deme <- droplevels(sub_deme[ok_d])
    y_dm <- as.numeric(yv[sub_rows] - unname(tapply(yv[sub_rows], sub_deme, mean)[as.character(sub_deme)]))
    ss_y <- sum(y_dm^2)
    if (!is.finite(ss_y) || ss_y <= 1e-8) next

    C_dm <- apply(C_full[sub_rows, , drop = FALSE], 2L, function(cv) {
      as.numeric(cv - unname(tapply(cv, sub_deme, mean)[as.character(sub_deme)]))
    })
    if (!is.matrix(C_dm)) C_dm <- matrix(C_dm, ncol = ncol(C_full))
    keep_c <- which(apply(C_dm, 2L, stats::sd) > 1e-5)
    if (length(keep_c) == 0L) next
    C_dm <- C_dm[, keep_c, drop = FALSE]
    C_eval <- C_full[, keep_c, drop = FALSE]

    fit_w <- tryCatch(stats::lm.fit(x = C_dm, y = y_dm), error = function(e) NULL)
    if (is.null(fit_w)) next
    beta_w <- fit_w$coefficients
    beta_w[!is.finite(beta_w)] <- 0.0
    rss_w <- sum(fit_w$residuals^2)
    r2_w <- max(0.0, 1.0 - (rss_w / ss_y))
    df1 <- max(1L, fit_w$rank)
    df2 <- sum(d_cnts[keep_d] - 1L) - df1
    if (df2 < 5L) next
    f_stat <- ((ss_y - rss_w) / df1) / max(rss_w / df2, 1e-12)
    p_val <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)

    if (is.finite(r2_w) && r2_w >= min_partial_r2 && is.finite(p_val) && p_val < 0.05) {
      c_ref_means <- colMeans(C_eval[ref_idx, , drop = FALSE], na.rm = TRUE)
      c_ref_means[!is.finite(c_ref_means)] <- 0.0
      fin_all <- which(is.finite(yv))
      C_cent <- sweep(C_eval[fin_all, , drop = FALSE], 2L, c_ref_means, "-")
      C_cent[!is.finite(C_cent)] <- 0.0
      adj_vals <- yv[fin_all] - as.numeric(C_cent %*% beta_w)
      df[[col_nm]][fin_all] <- round(adj_vals, 6)
    }
  }
  df
}

whiten_intraspecific_deme_plasticity <- function(
    X,
    df,
    ref_idx,
    deme_radius_km = 15.0,
    max_compression = 3.0
) {
  p <- ncol(X)
  if (is.null(X) || p < 1L || length(ref_idx) < max(10L, p + 4L)) {
    return(list(X = X, applied = FALSE, qst_proxy = rep(NA_real_, max(1L, p))))
  }

  deme_ids <- .build_ecological_deme_ids(df, ref_idx, deme_radius_km = deme_radius_km)
  is_geo_deme <- !grepl("^ISOLATE_", deme_ids)
  if (sum(is_geo_deme) < max(10L, p + 4L)) {
    return(list(X = X, applied = FALSE, qst_proxy = rep(NA_real_, p)))
  }

  geo_pos <- which(is_geo_deme)
  geo_ref_idx <- ref_idx[geo_pos]
  geo_demes <- deme_ids[geo_pos]
  d_tab <- table(geo_demes)
  multi_d <- names(d_tab[d_tab >= 2L])
  df_w <- sum(d_tab[multi_d] - 1L)
  n_unique_demes <- length(d_tab)

  # Require >= 3 distinct geographic demes AND >= max(6, p) within-deme degrees of freedom
  if (n_unique_demes < 3L || length(multi_d) < 2L || df_w < max(6L, p)) {
    return(list(X = X, applied = FALSE, qst_proxy = rep(NA_real_, p)))
  }

  X_ref <- X[geo_ref_idx, , drop = FALSE]
  sp_wt <- if ("spatial_balance_weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$spatial_balance_weight[geo_ref_idx]))
  } else {
    rep(1.0, length(geo_ref_idx))
  }
  sp_wt[!is.finite(sp_wt) | sp_wt <= 0.0] <- 1.0

  # Compute deme centroids and pooled within-deme residual matrix R_w
  u_demes <- names(d_tab)
  deme_means <- matrix(0.0, nrow = length(u_demes), ncol = p, dimnames = list(u_demes, colnames(X)))
  deme_wts <- numeric(length(u_demes))
  for (d_i in seq_along(u_demes)) {
    dk <- u_demes[d_i]
    m_k <- which(geo_demes == dk)
    wk <- sp_wt[m_k]
    deme_means[d_i, ] <- colSums(X_ref[m_k, , drop = FALSE] * wk) / max(sum(wk), 1e-12)
    # Weight each deme centroid gently by sqrt(n_k) * mean(spatial_balance_weight) to resist uneven herbarium density
    deme_wts[d_i] <- mean(wk) * sqrt(length(m_k))
  }

  in_m <- which(geo_demes %in% multi_d)
  R_w <- X_ref[in_m, , drop = FALSE] - deme_means[geo_demes[in_m], , drop = FALSE]
  S_w_raw <- crossprod(R_w) / max(df_w, 1L)
  S_w <- if (exists("compute_ledoit_wolf_cov", mode = "function")) {
    compute_ledoit_wolf_cov(S_w_raw, X_centered = R_w, n_obs = nrow(R_w))
  } else {
    S_w_raw + diag(0.05 * mean(diag(S_w_raw)) + 1e-4, p)
  }

  grand_mu <- colSums(deme_means * deme_wts) / max(sum(deme_wts), 1e-12)
  D_cent <- sweep(deme_means, 2L, grand_mu, "-")
  S_b <- crossprod(D_cent * sqrt(deme_wts)) / max(sum(deme_wts) - (sum(deme_wts^2) / max(sum(deme_wts), 1e-12)), 1e-6)

  eig_w <- eigen(S_w, symmetric = TRUE)
  U_w <- eig_w$vectors
  lam_w <- pmax(eig_w$values, 1e-6)
  var_b_proj <- pmax(0.0, diag(t(U_w) %*% S_b %*% U_w))
  # Subtract expected sampling noise (lam_w / harmonic_mean_deme_n) from between-deme variance
  harm_n <- length(u_demes) / sum(1.0 / as.numeric(d_tab))
  var_b_net <- pmax(0.0, var_b_proj - (lam_w / max(harm_n, 1.5)))
  eta <- var_b_net / lam_w
  qst_proxy <- eta / (eta + 2.0)

  floor_g <- 1.0 / max(1.0, as.numeric(max_compression))
  gain_vec <- pmin(1.0, pmax(floor_g, sqrt(eta / (eta + 1.0))))

  T_mat <- U_w %*% diag(gain_vec, nrow = p, ncol = p) %*% t(U_w)
  X_out <- sweep(X, 2L, grand_mu, "-") %*% T_mat
  X_out <- sweep(X_out, 2L, grand_mu, "+")
  colnames(X_out) <- colnames(X)
  list(X = X_out, applied = TRUE, qst_proxy = round(qst_proxy, 4), gain = round(gain_vec, 4))
}

compute_cluster_spatial_segregation_p <- function(
    cls_vec,
    df_sub,
    k_nn = 6L,
    max_dist_km = 300.0,
    permutations = 199L,
    seed = 42L
) {
  valid <- which(!is.na(cls_vec))
  if (length(valid) < 8L || length(unique(cls_vec[valid])) < 2L) return(1.0)
  sub_df <- df_sub[valid, , drop = FALSE]
  cl_v <- as.integer(cls_vec[valid])

  # 1. Check EPA Level III Ecoregion contingency if >= 2 ecoregions exist
  p_eco <- 1.0
  if ("regional_group" %in% names(sub_df)) {
    rg <- trimws(as.character(sub_df$regional_group))
    rg[is.na(rg) | !nzchar(rg) | tolower(rg) %in% c("na", "unknown", "none")] <- NA_character_
    ok_rg <- !is.na(rg)
    if (sum(ok_rg) >= 8L && length(unique(rg[ok_rg])) >= 2L && length(unique(cl_v[ok_rg])) >= 2L) {
      tab_rg <- table(cl_v[ok_rg], rg[ok_rg])
      p_eco <- tryCatch({
        set.seed(seed)
        stats::fisher.test(tab_rg, simulate.p.value = TRUE, B = max(499L, as.integer(permutations)))$p.value
      }, error = function(e) 1.0)
    }
  }

  # 2. Spatial k-NN / inverse-distance join-count concordance (Moran-type test on cluster assignments)
  lat_v <- rep(NA_real_, nrow(sub_df))
  lon_v <- rep(NA_real_, nrow(sub_df))
  for (lc in c("decimalLatitude", "latitude")) {
    if (lc %in% names(sub_df)) {
      vv <- suppressWarnings(as.numeric(sub_df[[lc]]))
      m <- !is.finite(lat_v) & is.finite(vv)
      lat_v[m] <- vv[m]
    }
  }
  for (lc in c("decimalLongitude", "longitude")) {
    if (lc %in% names(sub_df)) {
      vv <- suppressWarnings(as.numeric(sub_df[[lc]]))
      m <- !is.finite(lon_v) & is.finite(vv)
      lon_v[m] <- vv[m]
    }
  }
  ok_c <- which(is.finite(lat_v) & is.finite(lon_v) & abs(lat_v) <= 90.0 & abs(lon_v) <= 180.0)
  if (length(ok_c) < 8L || length(unique(cl_v[ok_c])) < 2L) return(p_eco)

  lat_r <- lat_v[ok_c] * (pi / 180.0)
  lon_r <- lon_v[ok_c] * (pi / 180.0)
  cl_c <- cl_v[ok_c]
  nc <- length(ok_c)

  dlat <- outer(lat_r, lat_r, "-")
  dlon <- outer(lon_r, lon_r, "-")
  hav <- sin(dlat * 0.5)^2 + (cos(lat_r) %o% cos(lat_r)) * (sin(dlon * 0.5)^2)
  dist_km <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hav))))
  diag(dist_km) <- Inf

  # Require >= 2 distinct spatial locations (> 5 km apart)
  if (max(ifelse(is.finite(dist_km), dist_km, 0.0)) <= 5.0) return(p_eco)

  W <- matrix(0.0, nrow = nc, ncol = nc)
  k_use <- min(as.integer(k_nn), nc - 1L)
  bw_km <- max(25.0, stats::median(apply(dist_km, 1L, min, na.rm = TRUE), na.rm = TRUE) * 2.0)
  for (i in seq_len(nc)) {
    di <- dist_km[i, ]
    # Exclude exact 0-km same-site co-collections from dominating if > 1 km neighbors exist,
    # or include all within max_dist_km so local neighborhood homogeneity is tested across sites
    ord_i <- order(di)[seq_len(k_use)]
    ord_i <- ord_i[is.finite(di[ord_i]) & di[ord_i] <= max_dist_km]
    if (length(ord_i) > 0L) {
      wi <- exp(-di[ord_i] / bw_km)
      W[i, ord_i] <- wi / max(sum(wi), 1e-12)
    }
  }
  W_sym <- 0.5 * (W + t(W))
  if (sum(W_sym) <= 1e-8) return(p_eco)

  same_mat <- outer(cl_c, cl_c, "==")
  obs_stat <- sum(W_sym * same_mat) / sum(W_sym)

  set.seed(seed)
  n_perm <- max(99L, as.integer(permutations))
  ge_cnt <- 0L
  for (b in seq_len(n_perm)) {
    p_cl <- sample(cl_c)
    perm_stat <- sum(W_sym * outer(p_cl, p_cl, "==")) / sum(W_sym)
    if (perm_stat >= obs_stat - 1e-12) ge_cnt <- ge_cnt + 1L
  }
  p_geo <- (ge_cnt + 1.0) / (n_perm + 1.0)
  min(p_eco, p_geo)
}

compute_spca_global_test <- function(
    X_core,
    df_core,
    k_nn = 8L,
    max_dist_km = 300.0,
    permutations = 199L,
    seed = 42L
) {
  n <- nrow(X_core)
  p <- ncol(X_core)
  if (is.null(X_core) || n < 10L || p < 1L || is.null(df_core)) {
    return(list(lambda_max = NA_real_, p_value = 1.0, scores = X_core))
  }
  lat_v <- rep(NA_real_, n)
  lon_v <- rep(NA_real_, n)
  for (lc in c("decimalLatitude", "latitude")) {
    if (lc %in% names(df_core)) {
      vv <- suppressWarnings(as.numeric(df_core[[lc]]))
      m <- !is.finite(lat_v) & is.finite(vv)
      lat_v[m] <- vv[m]
    }
  }
  for (lc in c("decimalLongitude", "longitude")) {
    if (lc %in% names(df_core)) {
      vv <- suppressWarnings(as.numeric(df_core[[lc]]))
      m <- !is.finite(lon_v) & is.finite(vv)
      lon_v[m] <- vv[m]
    }
  }
  ok <- which(is.finite(lat_v) & is.finite(lon_v) & abs(lat_v) <= 90.0 & abs(lon_v) <= 180.0)
  if (length(ok) < 10L) return(list(lambda_max = NA_real_, p_value = 1.0, scores = X_core))

  lat_r <- lat_v[ok] * (pi / 180.0)
  lon_r <- lon_v[ok] * (pi / 180.0)
  ng <- length(ok)
  dlat <- outer(lat_r, lat_r, "-")
  dlon <- outer(lon_r, lon_r, "-")
  hav <- sin(dlat * 0.5)^2 + (cos(lat_r) %o% cos(lat_r)) * (sin(dlon * 0.5)^2)
  dist_km <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hav))))
  diag(dist_km) <- Inf
  if (max(ifelse(is.finite(dist_km), dist_km, 0.0)) <= 5.0) {
    return(list(lambda_max = NA_real_, p_value = 1.0, scores = X_core))
  }

  L <- matrix(0.0, nrow = ng, ncol = ng)
  k_use <- min(as.integer(k_nn), ng - 1L)
  bw_km <- max(30.0, stats::median(apply(dist_km, 1L, min, na.rm = TRUE), na.rm = TRUE) * 2.0)
  for (i in seq_len(ng)) {
    di <- dist_km[i, ]
    ord_i <- order(di)[seq_len(k_use)]
    ord_i <- ord_i[is.finite(di[ord_i]) & di[ord_i] <= max_dist_km]
    if (length(ord_i) > 0L) {
      wi <- exp(-di[ord_i] / bw_km)
      L[i, ord_i] <- wi / max(sum(wi), 1e-12)
    }
  }
  L_sym <- 0.5 * (L + t(L))
  if (sum(L_sym) <= 1e-8) return(list(lambda_max = NA_real_, p_value = 1.0, scores = X_core))

  Xc <- scale(X_core[ok, , drop = FALSE], center = TRUE, scale = TRUE)
  Xc[!is.finite(Xc)] <- 0.0
  C_spca <- (t(Xc) %*% L_sym %*% Xc) / ng
  C_spca <- 0.5 * (C_spca + t(C_spca))
  eig_s <- eigen(C_spca, symmetric = TRUE)
  obs_lmax <- max(0.0, eig_s$values[1L])

  set.seed(seed)
  n_perm <- max(99L, as.integer(permutations))
  ge_cnt <- 0L
  for (b in seq_len(n_perm)) {
    idx_p <- sample.int(ng)
    Xp <- Xc[idx_p, , drop = FALSE]
    Cp <- (t(Xp) %*% L_sym %*% Xp) / ng
    Cp <- 0.5 * (Cp + t(Cp))
    l_p <- max(eigen(Cp, symmetric = TRUE, only.values = TRUE)$values)
    if (l_p >= obs_lmax - 1e-12) ge_cnt <- ge_cnt + 1L
  }
  p_global <- (ge_cnt + 1.0) / (n_perm + 1.0)
  list(lambda_max = round(obs_lmax, 5), p_value = round(p_global, 5), eigenvectors = eig_s$vectors, eigenvalues = eig_s$values)
}

.compute_ridgeline_valley_ratio <- function(pi_j, mu_j, Sig_j, pi_l, mu_l, Sig_l, n_grid = 101L) {
  p <- length(mu_j)
  reg_j <- Sig_j + diag(1e-5, p)
  reg_l <- Sig_l + diag(1e-5, p)
  inv_j <- tryCatch(solve(reg_j), error = function(e) diag(1.0 / pmax(diag(reg_j), 1e-5), p))
  inv_l <- tryCatch(solve(reg_l), error = function(e) diag(1.0 / pmax(diag(reg_l), 1e-5), p))
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

  # Identify local maxima along ridgeline [0, 1]
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
  # Take the leftmost and rightmost prominent peaks (near component centers alpha=0 and alpha=1)
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
    df_core = NULL,
    intraspecific_mode = FALSE,
    ridge_deep_gap = 0.20,
    ridge_sympatric_max = 0.65,
    permutations = 199L,
    seed = 42L
) {
  if (is.null(probs_eval) || !is.matrix(probs_eval) || ncol(probs_eval) < 2L) {
    k_one <- 1L
    return(list(
      best_k = k_one,
      classification = rep(1L, nrow(X_eval)),
      uncertainty = rep(0.0, nrow(X_eval)),
      probs_eval = matrix(1.0, nrow = nrow(X_eval), ncol = 1L),
      merged_any = FALSE
    ))
  }

  Z_core <- as.matrix(probs_core)
  Z_eval <- as.matrix(probs_eval)
  p <- ncol(X_core)
  n_c <- nrow(X_core)
  merged_any <- FALSE

  repeat {
    G_cur <- ncol(Z_core)
    if (G_cur < 2L) break

    nk <- colSums(Z_core, na.rm = TRUE)
    # Drop empty / degenerate zero-weight components immediately
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
      Sk + diag(0.02 * mean(diag(Sk)) + 1e-4, p)
    })

    best_pair <- NULL
    best_ratio <- -1.0
    cls_core <- max.col(Z_core, ties.method = "first")

    for (j in seq_len(G_cur - 1L)) {
      for (l in (j + 1L):G_cur) {
        rl <- .compute_ridgeline_valley_ratio(pis[j], mus[[j]], Sigs[[j]], pis[l], mus[[l]], Sigs[[l]])
        should_merge <- FALSE

        # Rule 1 (Hennig 2010; Cadena et al. 2018): Strictly unimodal ridgeline or negligible saddle
        if (rl$n_modes <= 1L || rl$ridge_ratio >= 0.92) {
          should_merge <- TRUE
        } else if (isTRUE(intraspecific_mode) && rl$ridge_ratio > ridge_deep_gap) {
          # Rule 2 (Intraspecific Plasticity Gate): Shallow saddle (> 0.20) is merged when the two
          # candidate clusters are geographically/ecologically sympatric (lack spatial segregation),
          # or when saddle ratio exceeds ridge_sympatric_max (0.65)
          if (rl$ridge_ratio > ridge_sympatric_max) {
            should_merge <- TRUE
          } else if (!is.null(df_core) && nrow(df_core) == n_c) {
            pair_idx <- which(cls_core %in% c(j, l))
            if (length(pair_idx) >= 8L && sum(cls_core == j) >= 2L && sum(cls_core == l) >= 2L) {
              p_seg <- compute_cluster_spatial_segregation_p(
                cls_core[pair_idx],
                df_core[pair_idx, , drop = FALSE],
                permutations = permutations,
                seed = seed + j * 31L + l
              )
              if (!is.finite(p_seg) || p_seg > 0.05) {
                should_merge <- TRUE
              }
            } else {
              should_merge <- TRUE
            }
          }
        }

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
