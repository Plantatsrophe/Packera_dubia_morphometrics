# ==============================================================================
# Module: scripts/analysis/step06_spatial_landscape_genomics.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Modular Phase 2 Landscape Genomics Helpers for Step 06:
#     1. Small-sample (`median_n < 5`) statistical guardrail explicitly disabling
#        Fst outlier scans (OutFLANK, BayeScan), SFS models (dadi, fastsimcoal2),
#        and 1D equilibrium clines (HZAR).
#     2. Wang (2013) MMRR (`D_geno ~ D_geo_log + D_eco_step05 + D_morph_step04`,
#        9,999 matrix permutations), read-only post-stream lookup.
#     3. `conStruct`-style spatial vs. non-spatial allelic covariance decay.
#     4. `FEEMS`-style Delaunay/KNN graph-Laplacian effective migration corridors.
# ==============================================================================

.haversine_km_mat <- function(lat, lon) {
  n <- length(lat)
  dmat <- matrix(0.0, nrow = n, ncol = n)
  r_earth <- 6371.0088
  rlat <- lat * pi / 180.0
  rlon <- lon * pi / 180.0
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      dlat <- rlat[j] - rlat[i]
      dlon <- rlon[j] - rlon[i]
      a <- sin(dlat / 2)^2 + cos(rlat[i]) * cos(rlat[j]) * sin(dlon / 2)^2
      dkm <- 2.0 * r_earth * asin(min(1.0, sqrt(max(0.0, a))))
      dmat[i, j] <- dkm
      dmat[j, i] <- dkm
    }
  }
  dmat
}

.zscore_vec <- function(v) {
  sd_v <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(sd_v) || sd_v <= 1e-12) return(rep(0.0, length(v)))
  (v - mean(v, na.rm = TRUE)) / sd_v
}

.fetch_posthoc_step04_step05_covariates <- function(pass_df, db_path, run_id, demo) {
  n <- nrow(pass_df)
  out_morph <- matrix(NA_real_, nrow = n, ncol = 5L)
  out_eco <- matrix(NA_real_, nrow = n, ncol = 6L)

  if (!exists("connect_packera_db", mode = "function") || isTRUE(demo) || !nzchar(db_path) || !file.exists(db_path)) {
    return(list(morph = out_morph, eco = out_eco))
  }
  conn <- tryCatch(connect_packera_db(db_path = db_path, read_only = TRUE, demo = FALSE), error = function(e) NULL)
  if (is.null(conn)) return(list(morph = out_morph, eco = out_eco))
  on.exit(disconnect_packera_db(conn), add = TRUE)

  mv <- tryCatch(fetch_run_master_view(conn, run_id = run_id), error = function(e) NULL)
  if (is.null(mv) || nrow(mv) == 0L) return(list(morph = out_morph, eco = out_eco))

  mcols <- intersect(paste0("PC", 1:5), names(mv))
  ecols <- intersect(c("bio1_temp_mean", "bio12_precip_annual", "bio4_temp_seasonality", "soil_ph", "soil_clay_pct", "elevation_m"), names(mv))

  cats <- trimws(as.character(pass_df$catalogNumber))
  mv_cats <- trimws(as.character(mv$catalogNumber))
  for (i in seq_len(n)) {
    if (is.na(cats[i]) || !nzchar(cats[i])) next
    hits <- which(mv_cats == cats[i])
    if (length(hits) == 0L) next
    if (length(mcols) > 0L) {
      m_sub <- apply(as.matrix(mv[hits, mcols, drop = FALSE]), 2L, function(x) stats::median(as.numeric(x), na.rm = TRUE))
      out_morph[i, seq_along(mcols)] <- m_sub
    }
    if (length(ecols) > 0L) {
      e_sub <- apply(as.matrix(mv[hits, ecols, drop = FALSE]), 2L, function(x) stats::median(as.numeric(x), na.rm = TRUE))
      out_eco[i, seq_along(ecols)] <- e_sub
    }
  }
  list(morph = out_morph, eco = out_eco)
}

run_step06_landscape_genomics <- function(
  pass_df,
  d_geno_mat,
  db_path = "",
  run_id = "canonical",
  demo = FALSE,
  n_perm = 9999L
) {
  n <- nrow(pass_df)
  pop_tab <- table(as.character(pass_df$population_id[!is.na(pass_df$population_id) & nzchar(pass_df$population_id)]))
  median_n_pop <- if (length(pop_tab) > 0L) stats::median(as.numeric(pop_tab)) else 1.0
  small_n_guard <- isTRUE(median_n_pop < 5.0)
  if (small_n_guard) {
    message(sprintf(
      "[Step 06 Landscape] small_n_guardrail_active=TRUE (median n/pop=%.1f < 5): disabling Fst outliers (OutFLANK/BayeScan), SFS demography (dadi/fsc2), and HZAR clines.",
      median_n_pop
    ))
  }

  if (requireNamespace("conStruct", quietly = TRUE)) {
    message("[Step 06 Landscape] Detected installed R package `conStruct`.")
  }

  lat <- suppressWarnings(as.numeric(pass_df$decimalLatitude))
  lon <- suppressWarnings(as.numeric(pass_df$decimalLongitude))
  geo_idx <- which(is.finite(lat) & is.finite(lon) & abs(lat) <= 90 & abs(lon) <= 180)

  if (length(geo_idx) < 4L) {
    summary_df <- data.frame(
      n_georef_samples = length(geo_idx),
      median_n_per_pop = median_n_pop,
      small_n_guardrail_active = small_n_guard,
      excluded_large_n_methods = if (small_n_guard) "WeirCockerham_FST_OutFLANK_BayeScan;SFS_dadi_fsc2_TajimaD;1D_cline_HZAR" else "none",
      mmrr_r2 = NA_real_, mmrr_p_overall = NA_real_,
      mmrr_beta_ibd = NA_real_, mmrr_p_ibd = NA_real_,
      mmrr_beta_ibe = NA_real_, mmrr_p_ibe = NA_real_,
      mmrr_beta_morph = NA_real_, mmrr_p_morph = NA_real_,
      construct_r2_spatial = NA_real_, construct_r2_nonspatial = NA_real_,
      feems_mean_log_weight = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(summary_df = summary_df, feems_edges_df = data.frame()))
  }

  d_geo <- .haversine_km_mat(lat[geo_idx], lon[geo_idx])
  d_geo_log <- log1p(d_geo)
  dg_sub <- d_geno_mat[geo_idx, geo_idx, drop = FALSE]

  covs <- .fetch_posthoc_step04_step05_covariates(pass_df[geo_idx, , drop = FALSE], db_path, run_id, demo)
  # Geographic latitude/longitude climate proxy fallback when sample lacks linked BioClim/SSURGO row
  eco_mat <- covs$eco
  miss_eco <- rowSums(is.finite(eco_mat)) == 0L
  if (any(miss_eco)) {
    eco_mat[miss_eco, 1L] <- abs(lat[geo_idx][miss_eco])
    eco_mat[miss_eco, 2L] <- lon[geo_idx][miss_eco]
  }
  for (j in seq_len(ncol(eco_mat))) {
    fin <- is.finite(eco_mat[, j])
    if (sum(fin) >= 2L && stats::sd(eco_mat[fin, j]) > 1e-9) {
      eco_mat[fin, j] <- (eco_mat[fin, j] - mean(eco_mat[fin, j])) / stats::sd(eco_mat[fin, j])
      eco_mat[!fin, j] <- 0.0
    } else {
      eco_mat[, j] <- 0.0
    }
  }
  d_eco <- as.matrix(stats::dist(eco_mat))

  morph_mat <- covs$morph
  has_morph <- sum(rowSums(is.finite(morph_mat)) > 0L) >= 4L
  if (has_morph) {
    for (j in seq_len(ncol(morph_mat))) {
      fin <- is.finite(morph_mat[, j])
      if (any(fin)) morph_mat[!fin, j] <- mean(morph_mat[fin, j]) else morph_mat[, j] <- 0.0
    }
    d_morph <- as.matrix(stats::dist(morph_mat))
  } else {
    d_morph <- matrix(0.0, nrow = length(geo_idx), ncol = length(geo_idx))
  }

  # Wang (2013) MMRR across lower-triangle vectors with 9,999 Mantel matrix permutations
  lt <- lower.tri(dg_sub)
  y_v <- .zscore_vec(dg_sub[lt])
  x1_v <- .zscore_vec(d_geo_log[lt])
  x2_v <- .zscore_vec(d_eco[lt])
  x3_v <- .zscore_vec(d_morph[lt])

  X_mmrr <- cbind(IBD = x1_v, IBE = x2_v)
  if (has_morph && stats::sd(x3_v) > 1e-8) X_mmrr <- cbind(X_mmrr, Morph = x3_v)
  fit0 <- stats::lm.fit(X_mmrr, y_v)
  betas <- fit0$coefficients
  b_ibd <- unname(betas["IBD"])
  b_ibe <- unname(betas["IBE"])
  b_morph <- if ("Morph" %in% names(betas)) unname(betas["Morph"]) else NA_real_
  ss_res <- sum(fit0$residuals^2)
  r2_obs <- max(0.0, 1.0 - ss_res / max(1e-10, sum((y_v - mean(y_v))^2)))

  # Fast matrix permutation using precomputed QR decomposition of predictors
  q_obj <- qr(X_mmrr)
  eff_perm <- min(as.integer(n_perm), 9999L)
  if (nrow(dg_sub) <= 6L) eff_perm <- min(eff_perm, 719L)
  set.seed(202611L)
  m <- length(geo_idx)
  ge_r2 <- 0L
  ge_ibd <- 0L
  ge_ibe <- 0L
  ge_morph <- 0L
  for (pm in seq_len(eff_perm)) {
    perm <- sample.int(m)
    yp <- .zscore_vec(dg_sub[perm, perm][lt])
    bp <- qr.coef(q_obj, yp)
    rp <- yp - X_mmrr %*% bp
    r2_p <- max(0.0, 1.0 - sum(rp^2) / max(1e-10, sum((yp - mean(yp))^2)))
    if (r2_p >= r2_obs) ge_r2 <- ge_r2 + 1L
    if (abs(bp["IBD"]) >= abs(b_ibd)) ge_ibd <- ge_ibd + 1L
    if (abs(bp["IBE"]) >= abs(b_ibe)) ge_ibe <- ge_ibe + 1L
    if ("Morph" %in% names(betas) && abs(bp["Morph"]) >= abs(b_morph)) ge_morph <- ge_morph + 1L
  }
  p_r2 <- (ge_r2 + 1.0) / (eff_perm + 1.0)
  p_ibd <- (ge_ibd + 1.0) / (eff_perm + 1.0)
  p_ibe <- (ge_ibe + 1.0) / (eff_perm + 1.0)
  p_morph <- if (is.finite(b_morph)) (ge_morph + 1.0) / (eff_perm + 1.0) else NA_real_

  # conStruct-inspired spatial vs non-spatial allelic covariance decay
  sim_geno <- pmax(0.0, 1.0 - dg_sub[lt])
  d_km_v <- d_geo[lt]
  decay_scale <- max(25.0, stats::median(d_km_v, na.rm = TRUE))
  pred_sp <- exp(-d_km_v / decay_scale)
  r2_sp <- max(0.0, stats::cor(sim_geno, pred_sp, use = "complete.obs")^2)
  r2_nonsp <- max(0.0, stats::cor(sim_geno, x2_v, use = "complete.obs")^2)

  # FEEMS-style KNN/Delaunay graph edge residual migration weights
  k_nn <- min(3L, m - 1L)
  fit_ibd <- stats::lm(dg_sub[lt] ~ d_geo_log[lt])
  resid_mat <- matrix(0.0, nrow = m, ncol = m)
  resid_mat[lt] <- fit_ibd$residuals
  resid_mat <- resid_mat + t(resid_mat)

  edge_rows <- list()
  seen_e <- character(0)
  for (a in seq_len(m)) {
    ord_n <- order(d_geo[a, ])
    nbs <- setdiff(ord_n, a)[seq_len(k_nn)]
    for (b in nbs) {
      ekey <- paste(min(a, b), max(a, b), sep = "_")
      if (ekey %in% seen_e) next
      seen_e <- c(seen_e, ekey)
      # Negative genetic residual vs geographic expectation -> effective migration corridor (+weight)
      log_w <- max(-2.5, min(2.5, -resid_mat[a, b] / max(0.02, stats::sd(resid_mat[lt], na.rm = TRUE))))
      edge_rows[[length(edge_rows) + 1L]] <- data.frame(
        sample_a = pass_df$genetic_sample_id[geo_idx[a]],
        sample_b = pass_df$genetic_sample_id[geo_idx[b]],
        lon_a = lon[geo_idx[a]], lat_a = lat[geo_idx[a]],
        lon_b = lon[geo_idx[b]], lat_b = lat[geo_idx[b]],
        geo_km = round(d_geo[a, b], 2),
        d_geno = round(dg_sub[a, b], 4),
        feems_log_migration = round(log_w, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  feems_df <- do.call(rbind, edge_rows)

  summary_df <- data.frame(
    n_georef_samples = length(geo_idx),
    median_n_per_pop = median_n_pop,
    small_n_guardrail_active = small_n_guard,
    excluded_large_n_methods = if (small_n_guard) "WeirCockerham_FST_OutFLANK_BayeScan;SFS_dadi_fsc2_TajimaD;1D_cline_HZAR" else "none",
    mmrr_r2 = round(r2_obs, 4), mmrr_p_overall = round(p_r2, 4),
    mmrr_beta_ibd = round(b_ibd, 4), mmrr_p_ibd = round(p_ibd, 4),
    mmrr_beta_ibe = round(b_ibe, 4), mmrr_p_ibe = round(p_ibe, 4),
    mmrr_beta_morph = round(b_morph, 4), mmrr_p_morph = round(p_morph, 4),
    construct_r2_spatial = round(r2_sp, 4),
    construct_r2_nonspatial = round(r2_nonsp, 4),
    feems_mean_log_weight = round(mean(feems_df$feems_log_migration, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  list(summary_df = summary_df, feems_edges_df = feems_df)
}
