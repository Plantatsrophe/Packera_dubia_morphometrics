# ==============================================================================
# Module: scripts/analysis/step06_spatial_landscape_genomics.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Modular Phase 2 Landscape Genomics Helpers for Step 06:
#     1. Small-sample (`median_n < 5`) statistical guardrail disabling Fst outliers.
#     2. Nested 5-model Wang (2013) MMRR with Gaussian AICc, Akaike weights,
#        partial IBE|IBD permutations, and centroid/rooftop coordinate quarantine.
#     3. `conStruct`-style spatial vs. non-spatial allelic covariance decay.
#     4. `FEEMS`-style kNN/MST graph-Laplacian migration corridors with LOOCV.
# ==============================================================================

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
}

.haversine_km_mat <- function(lat, lon) {
  n <- length(lat)
  dmat <- matrix(0.0, nrow = n, ncol = n)
  rlat <- lat * pi / 180.0; rlon <- lon * pi / 180.0
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      dlat <- rlat[j] - rlat[i]; dlon <- rlon[j] - rlon[i]
      a <- sin(dlat / 2)^2 + cos(rlat[i]) * cos(rlat[j]) * sin(dlon / 2)^2
      dkm <- 2.0 * 6371.0088 * asin(min(1.0, sqrt(max(0.0, a))))
      dmat[i, j] <- dkm; dmat[j, i] <- dkm
    }
  }
  dmat
}

.filter_step06_valid_geo_idx <- function(pass_df, lat, lon, db_path = "") {
  n <- length(lat)
  if (n == 0L) return(integer(0))
  valid <- is.finite(lat) & is.finite(lon) & abs(lat) <= 90 & abs(lon) <= 180 &
    !(abs(lat) < 1e-4 & abs(lon) < 1e-4) & !(lat > 15.0 & lon > 40.0)
  if ("is_centroid_or_inst_coord" %in% names(pass_df)) {
    cf <- pass_df$is_centroid_or_inst_coord
    c_bad <- if (is.logical(cf)) !is.na(cf) & cf else tolower(trimws(as.character(cf))) %in% c("1", "true", "t", "yes")
    valid <- valid & !c_bad
  }
  cats <- if ("catalogNumber" %in% names(pass_df)) trimws(as.character(pass_df$catalogNumber)) else rep("", n)
  cats[is.na(cats)] <- ""
  u_cats <- unique(cats[nzchar(cats)])
  if (length(u_cats) > 0L && !is.null(db_path) && nzchar(db_path) && file.exists(db_path) &&
      requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    conn <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO), error = function(e) NULL)
    if (!is.null(conn)) {
      v_cols <- tryCatch(DBI::dbGetQuery(conn, "PRAGMA table_info(vouchers);")$name, error = function(e) character(0))
      if ("is_centroid_or_inst_coord" %in% v_cols) {
        q_c <- paste(vapply(u_cats, function(s) as.character(DBI::dbQuoteString(conn, s)), character(1)), collapse = ", ")
        bad_v <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT catalogNumber FROM vouchers WHERE is_centroid_or_inst_coord = 1 AND catalogNumber IN (%s);", q_c))$catalogNumber, error = function(e) character(0))
        if (length(bad_v) > 0L) valid <- valid & !(cats %in% trimws(as.character(bad_v)))
      }
      tryCatch(DBI::dbDisconnect(conn), error = function(e) NULL)
    }
  }
  inst_lat <- c(35.9108, 36.0014, 35.7877, 40.8623, 42.3785, 38.8889, 38.6128, 41.8663, 39.9571, 32.7397, 30.2850, 43.0742, 44.9845, 35.9544, 33.9480, 37.2707, 30.4133, 40.4433, 39.5097, 42.8186, 45.3920, 45.5569, 43.6677)
  inst_lon <- c(-79.0505, -78.9382, -78.6703, -73.8803, -71.1160, -77.0260, -90.2594, -87.6170, -75.1714, -97.3631, -97.7340, -89.4012, -93.1832, -83.9295, -83.3773, -76.7153, -91.1800, -79.9496, -84.7345, -103.0005, -75.7160, -73.5575, -79.3948)
  for (i in which(valid & !nzchar(cats))) {
    rl1 <- lat[i] * pi / 180.0; rn1 <- lon[i] * pi / 180.0
    rl2 <- inst_lat * pi / 180.0; rn2 <- inst_lon * pi / 180.0
    a <- sin((rl2 - rl1) / 2)^2 + cos(rl1) * cos(rl2) * sin((rn2 - rn1) / 2)^2
    if (min(2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, a))))) <= 0.25) valid[i] <- FALSE
  }
  which(valid)
}

.zscore_vec <- function(v) {
  sd_v <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(sd_v) || sd_v <= 1e-12) return(rep(0.0, length(v)))
  (v - mean(v, na.rm = TRUE)) / sd_v
}

.scale_feature_dist <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    fin <- is.finite(mat[, j])
    if (sum(fin) >= 2L && stats::sd(mat[fin, j]) > 1e-9) {
      mat[fin, j] <- (mat[fin, j] - mean(mat[fin, j])) / stats::sd(mat[fin, j])
      mat[!fin, j] <- 0.0
    } else mat[, j] <- 0.0
  }
  as.matrix(stats::dist(mat))
}

.agg_by_catalog <- function(df, cols, cats) {
  out <- matrix(NA_real_, nrow = length(cats), ncol = length(cols))
  if (length(cols) == 0L || is.null(df) || nrow(df) == 0L) return(out)
  dcats <- trimws(as.character(df$catalogNumber))
  for (i in seq_along(cats)) {
    if (is.na(cats[i]) || !nzchar(cats[i])) next
    hits <- which(dcats == cats[i])
    if (length(hits) > 0L) out[i, ] <- apply(as.matrix(df[hits, cols, drop = FALSE]), 2L, function(x) stats::median(as.numeric(x), na.rm = TRUE))
  }
  out
}

.build_landscape_summary_df <- function(
  n_geo, med_n, guard, excl, r2 = NA_real_, p_all = NA_real_, b_ibd = NA_real_, p_ibd = NA_real_,
  b_ibe = NA_real_, p_ibe = NA_real_, b_m = NA_real_, p_m = NA_real_,
  r2_sp = NA_real_, r2_nonsp = NA_real_, feems_w = NA_real_, feems_lam = NA_real_
) {
  data.frame(
    n_georef_samples = n_geo, median_n_per_pop = med_n, small_n_guardrail_active = guard, excluded_large_n_methods = excl,
    mmrr_r2 = round(r2, 4), mmrr_p_overall = round(p_all, 4), mmrr_beta_ibd = round(b_ibd, 4), mmrr_p_ibd = round(p_ibd, 4),
    mmrr_beta_ibe = round(b_ibe, 4), mmrr_p_ibe = round(p_ibe, 4), mmrr_beta_morph = round(b_m, 4), mmrr_p_morph = round(p_m, 4),
    construct_r2_spatial = round(r2_sp, 4), construct_r2_nonspatial = round(r2_nonsp, 4),
    feems_mean_log_weight = round(feems_w, 4), feems_cv_lambda = feems_lam, stringsAsFactors = FALSE
  )
}

.fetch_posthoc_step04_step05_covariates <- function(pass_df, db_path, run_id, demo) {
  n <- nrow(pass_df); pc_cands <- paste0("PC", 1:15)
  eco_cands <- c("bio1_temp_mean", "bio12_precip_annual", "bio4_temp_seasonality", "soil_ph", "soil_clay_pct", "elevation_m")
  out_morph <- matrix(NA_real_, nrow = n, ncol = 5L); out_eco <- matrix(NA_real_, nrow = n, ncol = length(eco_cands))
  if (!exists("connect_packera_db", mode = "function") || isTRUE(demo) || !nzchar(db_path) || !file.exists(db_path)) {
    return(list(morph = out_morph, eco = out_eco))
  }
  conn <- tryCatch(connect_packera_db(db_path = db_path, read_only = TRUE, demo = FALSE), error = function(e) NULL)
  if (is.null(conn)) return(list(morph = out_morph, eco = out_eco))
  on.exit(disconnect_packera_db(conn), add = TRUE)

  mv <- if (nzchar(run_id) && exists("fetch_run_master_view", mode = "function")) tryCatch(fetch_run_master_view(conn, run_id = run_id), error = function(e) NULL) else NULL
  cats <- trimws(as.character(pass_df$catalogNumber))
  if (!is.null(mv) && nrow(mv) > 0L) {
    mcols <- intersect(pc_cands, names(mv)); ecols <- intersect(eco_cands, names(mv))
    if (length(mcols) > 0L) out_morph <- .agg_by_catalog(mv, mcols, cats)
    if (length(ecols) > 0L) out_eco[, seq_along(ecols)] <- .agg_by_catalog(mv, ecols, cats)
    return(list(morph = out_morph, eco = out_eco))
  }
  if (DBI::dbExistsTable(conn, "otu_morphometrics")) {
    otu_df <- tryCatch({
      if (nzchar(run_id)) DBI::dbGetQuery(conn, "SELECT * FROM otu_morphometrics WHERE run_id = ?", params = list(run_id))
      else DBI::dbGetQuery(conn, "SELECT * FROM otu_morphometrics")
    }, error = function(e) NULL)
    mcols <- if (!is.null(otu_df) && nrow(otu_df) > 0L) intersect(pc_cands, names(otu_df)) else character(0)
    if (length(mcols) > 0L) out_morph <- .agg_by_catalog(otu_df, mcols, cats)
  }
  if (DBI::dbExistsTable(conn, "environmental_cache") && all(c("decimalLatitude", "decimalLongitude") %in% names(pass_df))) {
    ec_df <- tryCatch(DBI::dbGetQuery(conn, "SELECT * FROM environmental_cache"), error = function(e) NULL)
    if (!is.null(ec_df) && nrow(ec_df) > 0L && "lat_lon_key" %in% names(ec_df)) {
      ecols <- intersect(eco_cands, names(ec_df))
      if (length(ecols) > 0L) {
        lat_p <- suppressWarnings(as.numeric(pass_df$decimalLatitude)); lon_p <- suppressWarnings(as.numeric(pass_df$decimalLongitude))
        keys <- ifelse(is.finite(lat_p) & is.finite(lon_p), sprintf("%.4f_%.4f", round(lat_p, 4), round(lon_p, 4)), NA_character_)
        idx_e <- match(keys, as.character(ec_df$lat_lon_key)); ok_e <- which(!is.na(idx_e))
        if (length(ok_e) > 0L) for (j in seq_along(ecols)) out_eco[ok_e, j] <- suppressWarnings(as.numeric(ec_df[[ecols[j]]][idx_e[ok_e]]))
      }
    }
  }
  list(morph = out_morph, eco = out_eco)
}

.fit_mmrr_model_table <- function(dg_sub, d_geo_log, d_eco, d_morph, use_morph, eff_perm) {
  ut <- upper.tri(dg_sub); n_pairs <- sum(ut); m <- nrow(dg_sub)
  y_v <- .zscore_vec(dg_sub[ut]); x1_v <- .zscore_vec(d_geo_log[ut]); x2_v <- .zscore_vec(d_eco[ut])
  x3_v <- if (use_morph) .zscore_vec(d_morph[ut]) else numeric(0)
  sst <- max(1e-12, sum((y_v - mean(y_v))^2))
  specs <- list(Null = matrix(1.0, nrow = n_pairs, ncol = 1L, dimnames = list(NULL, "Intercept")),
                IBD_only = cbind(Intercept = 1.0, IBD = x1_v), IBE_only = cbind(Intercept = 1.0, IBE = x2_v),
                IBD_plus_IBE = cbind(Intercept = 1.0, IBD = x1_v, IBE = x2_v))
  if (use_morph) specs[["IBD_plus_IBE_plus_Morph"]] <- cbind(Intercept = 1.0, IBD = x1_v, IBE = x2_v, Morph = x3_v)

  m_names <- names(specs); k_vec <- vapply(specs, ncol, integer(1))
  fits <- lapply(specs, function(X) stats::lm.fit(X, y_v))
  rss_vec <- vapply(fits, function(f) sum(f$residuals^2), numeric(1))
  r2_vec <- pmax(0.0, 1.0 - rss_vec / sst)
  adj_r2_vec <- vapply(seq_along(m_names), function(i) if (k_vec[i] <= 1L) 0.0 else 1.0 - (1.0 - r2_vec[i]) * (n_pairs - 1.0) / max(1.0, n_pairs - k_vec[i]), numeric(1))
  aicc_vec <- vapply(seq_along(m_names), function(i) {
    kp <- k_vec[i]
    n_pairs * log(max(1e-12, rss_vec[i] / n_pairs)) + 2.0 * (kp + 1.0) + (2.0 * (kp + 1.0) * (kp + 2.0)) / max(1.0, n_pairs - kp - 2.0)
  }, numeric(1))
  delta_aicc <- aicc_vec - min(aicc_vec); w_aicc <- exp(-0.5 * delta_aicc); w_aicc <- w_aicc / sum(w_aicc)
  full_name <- m_names[length(m_names)]; fb <- fits[[full_name]]$coefficients
  b_ibd <- unname(fb["IBD"]); b_ibe <- unname(fb["IBE"]); b_morph <- if ("Morph" %in% names(fb)) unname(fb["Morph"]) else NA_real_

  q_list <- lapply(specs, qr); q_ibd <- q_list[["IBD_only"]]
  e_y <- fits[["IBD_only"]]$residuals; ssr_y_ibd <- max(1e-12, sum(e_y^2))
  r2_ibd <- r2_vec[match("IBD_only", m_names)]; delta_r2_obs <- r2_vec[match("IBD_plus_IBE", m_names)] - r2_ibd
  perms_mm <- replicate(eff_perm, sample.int(m), simplify = FALSE)
  eval_mm_chunk <- function(idx_v) {
    gm <- integer(length(m_names)); gibd <- 0L; gibe <- 0L; gmp <- 0L; gpi <- 0L
    for (pm in idx_v) {
      perm <- perms_mm[[pm]]; yp <- .zscore_vec(dg_sub[perm, perm][ut]); sst_p <- max(1e-12, sum((yp - mean(yp))^2))
      for (j in seq_along(m_names)[-1L]) {
        bp <- qr.coef(q_list[[j]], yp); bp[is.na(bp)] <- 0.0
        if (max(0.0, 1.0 - sum((yp - specs[[j]] %*% bp)^2) / sst_p) >= r2_vec[j]) gm[j] <- gm[j] + 1L
        if (m_names[j] == full_name) {
          if (is.finite(b_ibd) && is.finite(bp["IBD"]) && abs(bp["IBD"]) >= abs(b_ibd)) gibd <- gibd + 1L
          if (is.finite(b_ibe) && is.finite(bp["IBE"]) && abs(bp["IBE"]) >= abs(b_ibe)) gibe <- gibe + 1L
          if ("Morph" %in% names(bp) && is.finite(b_morph) && is.finite(bp["Morph"]) && abs(bp["Morph"]) >= abs(b_morph)) gmp <- gmp + 1L
        }
      }
      x2_p <- .zscore_vec(d_eco[perm, perm][ut]); c_x2 <- qr.coef(q_ibd, x2_p); c_x2[is.na(c_x2)] <- 0.0
      e_x2 <- x2_p - specs[["IBD_only"]] %*% c_x2; denom_x2 <- sum(e_x2^2)
      dr2_p <- if (denom_x2 > 1e-12) ((sum(e_y * e_x2))^2 / (denom_x2 * ssr_y_ibd)) * (1.0 - r2_ibd) else 0.0
      if (dr2_p >= max(0.0, delta_r2_obs) - 1e-12) gpi <- gpi + 1L
    }
    c(gm, gibd, gibe, gmp, gpi)
  }
  nw_mm <- if (eff_perm >= 99L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, ceiling(eff_perm / 16L)) else 1L
  mm_tot <- if (nw_mm > 1L) Reduce("+", packera_mclapply(parallel::splitIndices(eff_perm, nw_mm), eval_mm_chunk, mc.cores = nw_mm)) else eval_mm_chunk(seq_len(eff_perm))
  nk <- length(m_names); ge_mod <- mm_tot[seq_len(nk)]; ge_ibd <- mm_tot[nk + 1L]; ge_ibe <- mm_tot[nk + 2L]; ge_morph <- mm_tot[nk + 3L]; ge_part_ibe <- mm_tot[nk + 4L]

  perm_p <- vapply(seq_along(m_names), function(j) if (j == 1L) 1.0 else (ge_mod[j] + 1.0) / (eff_perm + 1.0), numeric(1))
  p_part_ibe <- (ge_part_ibe + 1.0) / (eff_perm + 1.0)
  models_df <- data.frame(
    model_name = m_names, n_samples = rep(as.integer(m), length(m_names)), n_pairs = rep(as.integer(n_pairs), length(m_names)),
    k_params = as.integer(k_vec), r2 = round(r2_vec, 4), adj_r2 = round(adj_r2_vec, 4), aicc = round(aicc_vec, 4),
    delta_aicc = round(delta_aicc, 4), akaike_weight = round(w_aicc, 4),
    delta_r2_ibe_given_ibd = ifelse(m_names == "IBD_plus_IBE", round(delta_r2_obs, 4), NA_real_),
    perm_p_value = round(perm_p, 4), partial_ibe_p_value = ifelse(m_names == "IBD_plus_IBE", round(p_part_ibe, 4), NA_real_),
    is_selected_model = seq_along(m_names) == which.min(aicc_vec), stringsAsFactors = FALSE
  )
  list(
    models_df = models_df, r2_full = r2_vec[length(m_names)], p_full = perm_p[length(m_names)],
    b_ibd = b_ibd, p_ibd = if (is.finite(b_ibd)) (ge_ibd + 1.0) / (eff_perm + 1.0) else NA_real_,
    b_ibe = b_ibe, p_ibe = if (is.finite(b_ibe)) (ge_ibe + 1.0) / (eff_perm + 1.0) else NA_real_,
    b_morph = b_morph, p_morph = if (is.finite(b_morph)) (ge_morph + 1.0) / (eff_perm + 1.0) else NA_real_
  )
}

.select_feems_lambda_and_edges <- function(pass_df, geo_idx, lat, lon, dg_sub, d_geo, d_geo_log) {
  m <- length(geo_idx); ut <- upper.tri(dg_sub); k_nn <- min(3L, m - 1L); A_sp <- matrix(0.0, nrow = m, ncol = m)
  add_edge <- function(u, v) { w <- 1.0 / (1.0 + d_geo_log[u, v]); A_sp[u, v] <<- w; A_sp[v, u] <<- w }
  for (a in seq_len(m)) for (b in setdiff(order(d_geo[a, ]), a)[seq_len(k_nn)]) add_edge(a, b)
  in_tree <- rep(FALSE, m); in_tree[1L] <- TRUE
  for (step_mst in seq_len(m - 1L)) {
    pos <- which.min(d_geo[in_tree, !in_tree, drop = FALSE])
    u_idx <- which(in_tree)[((pos - 1L) %% sum(in_tree)) + 1L]
    v_idx <- which(!in_tree)[((pos - 1L) %/% sum(in_tree)) + 1L]
    add_edge(u_idx, v_idx); in_tree[v_idx] <- TRUE
  }
  L_sp <- diag(rowSums(A_sp)) - A_sp
  lambdas <- c(1e-3, 1e-2, 1e-1, 1.0, 10.0)
  cv_mse <- vapply(lambdas, function(lam) {
    sq_errs <- numeric(0)
    for (v in seq_len(m)) {
      tr <- setdiff(seq_len(m), v); ut_tr <- upper.tri(d_geo_log[tr, tr, drop = FALSE])
      fit_tr <- stats::lm.fit(cbind(1.0, d_geo_log[tr, tr, drop = FALSE][ut_tr]), dg_sub[tr, tr, drop = FALSE][ut_tr])
      b_tr <- ifelse(is.na(fit_tr$coefficients), 0.0, fit_tr$coefficients)
      res_tr <- matrix(0.0, nrow = m - 1L, ncol = m - 1L); res_tr[ut_tr] <- fit_tr$residuals; res_tr <- res_tr + t(res_tr)
      u_tr <- rowSums(res_tr) / max(1.0, m - 2.0); L_tr <- diag(rowSums(A_sp[tr, tr, drop = FALSE])) - A_sp[tr, tr, drop = FALSE]
      u_sm <- tryCatch(as.numeric(solve(diag(m - 1L) + lam * L_tr, u_tr)), error = function(e) u_tr)
      w_v <- A_sp[v, tr]; w_v <- if (sum(w_v) > 0.0) w_v / sum(w_v) else rep(1.0 / (m - 1L), m - 1L)
      sq_errs <- c(sq_errs, (dg_sub[v, tr] - (b_tr[1L] + b_tr[2L] * d_geo_log[v, tr] + (sum(w_v * u_sm) + u_sm)))^2)
    }
    mean(sq_errs, na.rm = TRUE)
  }, numeric(1))
  feems_cv_lambda <- lambdas[which.min(cv_mse)]

  fit_ibd <- stats::lm(dg_sub[ut] ~ d_geo_log[ut])
  resid_mat <- matrix(0.0, nrow = m, ncol = m); resid_mat[ut] <- fit_ibd$residuals; resid_mat <- resid_mat + t(resid_mat)
  S_sp <- tryCatch(solve(diag(m) + feems_cv_lambda * L_sp), error = function(e) diag(m))
  resid_smooth <- S_sp %*% resid_mat %*% t(S_sp); denom <- max(0.02, stats::sd(resid_smooth[ut], na.rm = TRUE))
  idx_e <- which(A_sp > 0.0 & upper.tri(A_sp), arr.ind = TRUE); ia <- idx_e[, 1L]; ib <- idx_e[, 2L]
  edges_df <- data.frame(
    sample_a = pass_df$genetic_sample_id[geo_idx[ia]], sample_b = pass_df$genetic_sample_id[geo_idx[ib]],
    lon_a = lon[geo_idx[ia]], lat_a = lat[geo_idx[ia]], lon_b = lon[geo_idx[ib]], lat_b = lat[geo_idx[ib]],
    geo_km = round(d_geo[idx_e], 2), d_geno = round(dg_sub[idx_e], 4), feems_log_migration = round(pmax(-2.5, pmin(2.5, -resid_smooth[idx_e] / denom)), 4),
    stringsAsFactors = FALSE
  )
  list(edges_df = edges_df, cv_lambda = feems_cv_lambda)
}

run_step06_landscape_genomics <- function(pass_df, d_geno_mat, db_path = "", run_id = "", demo = FALSE, seed = 202611L, n_perm = 999L, out_dir = "") {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (is.null(old_seed)) { if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }
    else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))

  pop_tab <- table(as.character(pass_df$population_id[!is.na(pass_df$population_id) & nzchar(pass_df$population_id)]))
  median_n_pop <- if (length(pop_tab) > 0L) stats::median(as.numeric(pop_tab)) else 1.0
  small_n_guard <- isTRUE(median_n_pop < 5.0)
  excl_str <- if (small_n_guard) "WeirCockerham_FST_OutFLANK_BayeScan;SFS_dadi_fsc2_TajimaD;1D_cline_HZAR" else "none"
  if (small_n_guard) message(sprintf("[Step 06 Landscape] small_n_guardrail_active=TRUE (median n/pop=%.1f < 5): disabling Fst outliers (OutFLANK/BayeScan), SFS demography (dadi/fsc2), and HZAR clines.", median_n_pop))
  if (requireNamespace("conStruct", quietly = TRUE)) message("[Step 06 Landscape] Detected installed R package `conStruct`.")

  lat <- suppressWarnings(as.numeric(pass_df$decimalLatitude)); lon <- suppressWarnings(as.numeric(pass_df$decimalLongitude))
  geo_idx <- .filter_step06_valid_geo_idx(pass_df, lat, lon, db_path = db_path)
  if (length(geo_idx) < 4L) {
    empty_mmrr <- data.frame(
      model_name = character(0), n_samples = integer(0), n_pairs = integer(0), k_params = integer(0),
      r2 = numeric(0), adj_r2 = numeric(0), aicc = numeric(0), delta_aicc = numeric(0), akaike_weight = numeric(0),
      delta_r2_ibe_given_ibd = numeric(0), perm_p_value = numeric(0), partial_ibe_p_value = numeric(0),
      is_selected_model = logical(0), stringsAsFactors = FALSE
    )
    if (nzchar(out_dir)) { dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); utils::write.csv(empty_mmrr, file.path(out_dir, "mmrr_model_selection.csv"), row.names = FALSE) }
    return(list(summary_df = .build_landscape_summary_df(length(geo_idx), median_n_pop, small_n_guard, excl_str), feems_edges_df = data.frame(), mmrr_models_df = empty_mmrr))
  }

  d_geo <- .haversine_km_mat(lat[geo_idx], lon[geo_idx]); d_geo_log <- log1p(d_geo)
  dg_sub <- d_geno_mat[geo_idx, geo_idx, drop = FALSE]; ut <- upper.tri(dg_sub)
  covs <- .fetch_posthoc_step04_step05_covariates(pass_df[geo_idx, , drop = FALSE], db_path, run_id, demo)
  eco_mat <- covs$eco; miss_eco <- rowSums(is.finite(eco_mat)) == 0L
  if (any(miss_eco)) { eco_mat[miss_eco, 1L] <- abs(lat[geo_idx][miss_eco]); eco_mat[miss_eco, 2L] <- lon[geo_idx][miss_eco] }
  d_eco <- .scale_feature_dist(eco_mat)

  morph_mat <- covs$morph; has_morph <- sum(rowSums(is.finite(morph_mat)) > 0L) >= 4L
  d_morph <- if (has_morph) .scale_feature_dist(morph_mat) else matrix(0.0, nrow = length(geo_idx), ncol = length(geo_idx))
  use_morph <- has_morph && stats::sd(d_morph[ut], na.rm = TRUE) > 1e-8
  eff_perm <- min(max(1L, as.integer(n_perm)), 9999L); if (nrow(dg_sub) <= 6L) eff_perm <- min(eff_perm, 719L)
  mmrr_res <- .fit_mmrr_model_table(dg_sub, d_geo_log, d_eco, d_morph, use_morph, eff_perm)
  if (nzchar(out_dir)) { dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); utils::write.csv(mmrr_res$models_df, file.path(out_dir, "mmrr_model_selection.csv"), row.names = FALSE) }

  sim_geno <- pmax(0.0, 1.0 - dg_sub[ut]); d_km_v <- d_geo[ut]; decay_scale <- max(25.0, stats::median(d_km_v, na.rm = TRUE))
  r2_sp <- suppressWarnings(stats::cor(sim_geno, exp(-d_km_v / decay_scale), use = "complete.obs")^2)
  r2_nonsp <- suppressWarnings(stats::cor(sim_geno, .zscore_vec(d_eco[ut]), use = "complete.obs")^2)
  r2_sp <- if (is.finite(r2_sp)) max(0.0, r2_sp) else 0.0; r2_nonsp <- if (is.finite(r2_nonsp)) max(0.0, r2_nonsp) else 0.0
  feems_res <- .select_feems_lambda_and_edges(pass_df, geo_idx, lat, lon, dg_sub, d_geo, d_geo_log)
  summary_df <- .build_landscape_summary_df(
    length(geo_idx), median_n_pop, small_n_guard, excl_str, mmrr_res$r2_full, mmrr_res$p_full,
    mmrr_res$b_ibd, mmrr_res$p_ibd, mmrr_res$b_ibe, mmrr_res$p_ibe, mmrr_res$b_morph, mmrr_res$p_morph,
    r2_sp, r2_nonsp, mean(feems_res$edges_df$feems_log_migration, na.rm = TRUE), feems_res$cv_lambda
  )
  list(summary_df = summary_df, feems_edges_df = feems_res$edges_df, mmrr_models_df = mmrr_res$models_df)
}
