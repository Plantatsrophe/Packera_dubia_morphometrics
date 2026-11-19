# ==============================================================================
# Step 05 Module: Circular von Mises Hopkins AICc Models & Phase Divergence Tests
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

doy_to_rad <- function(doy) 2.0 * pi * (as.numeric(doy) / 365.25)
wrap_rad <- function(d) atan2(sin(d), cos(d))

estimate_von_mises_kappa <- function(R_bar, n = Inf) {
  rb <- suppressWarnings(as.numeric(R_bar)[1L])
  n_val <- suppressWarnings(as.numeric(n)[1L])
  if (!is.finite(rb) || rb <= 0.0) return(0.0)
  rb <- min(max(rb, 0.0), 1.0 - 1e-8)

  # Best & Fisher (1981) initializer + 6 Newton-Raphson iterations on A1(k) = R_bar
  k <- if (rb < 0.53) {
    2.0 * rb + rb^3 + (5.0 / 6.0) * rb^5
  } else if (rb < 0.85) {
    -0.4 + 1.39 * rb + 0.43 / (1.0 - rb)
  } else {
    1.0 / max(1e-8, rb^3 - 4.0 * rb^2 + 3.0 * rb)
  }
  for (iter in seq_len(6L)) {
    if (!is.finite(k) || k <= 1e-8) { k <- 1e-6; break }
    i0 <- besselI(k, nu = 0.0, expon.scaled = TRUE)
    i1 <- besselI(k, nu = 1.0, expon.scaled = TRUE)
    if (!is.finite(i0) || i0 <= 1e-300) break
    a1 <- i1 / i0
    k_next <- k - (a1 - rb) / max(1e-12, 1.0 - a1 * a1 - (a1 / k))
    if (!is.finite(k_next)) break
    k <- min(500.0, max(1e-8, k_next))
  }
  # Mardia & Jupp (2000) bias correction for n < 15
  if (is.finite(n_val) && n_val >= 2.0 && n_val < 15.0) {
    k <- if (k < 2.0) max(k - (2.0 / max(1e-6, n_val * k)), 0.0) else max(((n_val - 1.0)^3 * k) / (n_val^3 + n_val), 0.0)
  } else if (is.finite(n_val) && n_val < 2.0) {
    k <- 0.0
  }
  min(500.0, max(0.0, k))
}

loglik_von_mises <- function(theta, mu_hat, kappa) {
  th <- as.numeric(theta); th <- th[is.finite(th)]
  n_obs <- length(th)
  if (n_obs == 0L) return(0.0)
  mu_val <- if (is.finite(mu_hat)) as.numeric(mu_hat)[1L] else 0.0
  kap <- if (is.finite(kappa)) as.numeric(kappa)[1L] else 0.0
  if (kap <= 1e-8) return(-n_obs * log(2.0 * pi))
  i0_sc <- max(1e-300, besselI(kap, nu = 0.0, expon.scaled = TRUE))
  sum(kap * cos(th - mu_val) - log(2.0 * pi) - (log(i0_sc) + kap))
}

run_rayleigh_uniformity_test <- function(theta) {
  th <- as.numeric(theta); th <- th[is.finite(th)]
  n <- length(th)
  if (n == 0L) {
    return(list(n = 0L, C_bar = NA_real_, S_bar = NA_real_, R_bar = NA_real_, Z = NA_real_, mu_hat = NA_real_, p_value = 1.0))
  }
  c_bar <- mean(cos(th)); s_bar <- mean(sin(th))
  r_bar <- sqrt(c_bar^2 + s_bar^2)
  r_n <- n * r_bar
  # Greenwood & Durand (1955) Rayleigh p-value approximation
  radicand <- max(0.0, 1.0 + 4.0 * n + 4.0 * (n - r_n) * (n + r_n))
  p_val <- min(1.0, max(0.0, exp(sqrt(radicand) - (1.0 + 2.0 * n))))
  list(n = as.integer(n), C_bar = c_bar, S_bar = s_bar, R_bar = r_bar, Z = n * (r_bar^2), mu_hat = atan2(s_bar, c_bar), p_value = p_val)
}

.extract_circ_cov <- function(df, candidates, lo = -Inf, hi = Inf) {
  hits <- intersect(candidates, names(df))
  if (length(hits) == 0L) return(rep(NA_real_, nrow(df)))
  v <- suppressWarnings(as.numeric(df[[hits[1L]]]))
  v[!is.finite(v) | v < lo | v > hi] <- NA_real_
  v
}

.eval_vm_profile <- function(theta_obs, off_days) {
  r <- theta_obs - (2.0 * pi / 365.25) * off_days
  c_r <- mean(cos(r)); s_r <- mean(sin(r))
  mu_0 <- atan2(s_r, c_r)
  r_bar <- sqrt(c_r^2 + s_r^2)
  kap <- estimate_von_mises_kappa(r_bar, length(r))
  list(mu_0 = mu_0, kappa = kap, R_bar = r_bar, loglik = loglik_von_mises(r, mu_0, kap))
}

select_circular_hopkins_baseline_model <- function(flowering_df) {
  df_in <- as.data.frame(flowering_df, check.names = FALSE, stringsAsFactors = FALSE)
  lat_raw <- .extract_circ_cov(df_in, c("decimalLatitude", "latitude"), -90.0, 90.0)
  doy_raw <- .extract_circ_cov(df_in, c("doy"), 1.0, 366.0)
  ok <- is.finite(lat_raw) & is.finite(doy_raw)
  df_v <- df_in[ok, , drop = FALSE]
  lat_v <- lat_raw[ok]
  doy_v <- doy_raw[ok]
  n <- length(lat_v)

  if (n < 2L) {
    ray0 <- run_rayleigh_uniformity_test(if (n == 1L) doy_to_rad(doy_v) else numeric(0))
    mu_fb <- if (n == 1L) wrap_rad(doy_to_rad(doy_v[1L])) else 0.0
    tbl0 <- data.frame(
      Model = "M0_Constant_Phase", k_params = 2L, N = as.integer(n), LogLik = NA_real_, Kappa = 0.0,
      Slope_Lat_d_per_deg = 0.0, AICc = NA_real_, Delta_AICc = NA_real_, Akaike_Weight = NA_real_,
      Delta_AICc_vs_Hopkins1918 = NA_real_, stringsAsFactors = FALSE
    )
    return(list(
      selected_model = "M0_Constant_Phase", n = as.integer(n), mu_0 = mu_fb, kappa = 0.0,
      b_lat = 0.0, b_lon = 0.0, b_elev = 0.0, b_mat = 0.0, b_year = 0.0,
      lat_ref = if (n == 1L) lat_v[1L] else 0.0, lon_ref = 0.0, elev_ref = 0.0, mat_ref = 0.0, year_ref = 1990.0,
      rayleigh_z = ray0$Z, rayleigh_p = ray0$p_value, r_bar = ray0$R_bar, is_uniform_null = TRUE,
      aicc_weight = NA_real_, delta_aicc_vs_hopkins1918 = NA_real_, model_selection_table = tbl0,
      predict_doy_rad = function(newdata_df) {
        nd_lat <- .extract_circ_cov(as.data.frame(newdata_df), c("decimalLatitude", "latitude"), -90.0, 90.0)
        ifelse(is.finite(nd_lat), wrap_rad(mu_fb), NA_real_)
      }
    ))
  }

  lon_v <- .extract_circ_cov(df_v, c("decimalLongitude", "longitude"), -180.0, 180.0)
  elev_v <- .extract_circ_cov(df_v, c("elevation_m", "elevation", "minimumElevationInMeters"), -500.0, 9000.0)
  mat_v <- .extract_circ_cov(df_v, c("bio1_temp_mean", "mat"), -60.0, 60.0)
  yr_v <- .extract_circ_cov(df_v, c("year"), 1750.0, 2100.0)
  theta_obs <- doy_to_rad(doy_v)
  ray <- run_rayleigh_uniformity_test(theta_obs)
  is_uniform_null <- isTRUE(is.finite(ray$p_value) && is.finite(ray$R_bar) && ray$p_value >= 0.05 && ray$R_bar < 0.15)

  med_or <- function(x, fb = 0.0) { f <- x[is.finite(x)]; if (length(f) > 0L) stats::median(f) else fb }
  lat_ref <- med_or(lat_v, 0.0); lon_ref <- med_or(lon_v, 0.0); elev_ref <- med_or(elev_v, 0.0)
  mat_ref <- med_or(mat_v, 0.0); yr_ref <- med_or(yr_v, 1990.0)

  lat_c <- lat_v - lat_ref
  lon_c <- ifelse(is.finite(lon_v), lon_v - lon_ref, 0.0)
  elev_c <- ifelse(is.finite(elev_v), elev_v - elev_ref, 0.0)
  mat_c <- ifelse(is.finite(mat_v), mat_v - mat_ref, 0.0)
  yr_c <- ifelse(is.finite(yr_v), yr_v - yr_ref, 0.0)

  m_list <- list()
  add_cand <- function(nm, k_p, prf, bl = 0.0, blo = 0.0, be = 0.0, bm = 0.0, by = 0.0) {
    m_list[[length(m_list) + 1L]] <<- list(
      name = nm, k = as.integer(k_p), mu_0 = prf$mu_0, kappa = prf$kappa, loglik = prf$loglik,
      b_lat = bl, b_lon = blo, b_elev = be, b_mat = bm, b_year = by
    )
  }

  has_elev <- sum(is.finite(elev_v)) >= max(5L, as.integer(ceiling(0.50 * n))) &&
              is.finite(stats::sd(elev_v, na.rm = TRUE)) && stats::sd(elev_v, na.rm = TRUE) > 5.0
  has_lon <- sum(is.finite(lon_v)) >= 10L &&
             is.finite(stats::sd(lon_v, na.rm = TRUE)) && stats::sd(lon_v, na.rm = TRUE) > 0.5
  has_mat <- !has_elev && sum(is.finite(mat_v)) >= max(10L, as.integer(ceiling(0.50 * n))) &&
             is.finite(stats::sd(mat_v, na.rm = TRUE)) && stats::sd(mat_v, na.rm = TRUE) > 0.5
  has_yr <- sum(is.finite(yr_v)) >= max(15L, as.integer(ceiling(0.80 * n))) &&
            is.finite(stats::sd(yr_v, na.rm = TRUE)) && stats::sd(yr_v, na.rm = TRUE) >= 5.0
  be_m1 <- if (has_elev) 4.00 / 121.92 else 0.0
  blo_m1 <- if (has_lon) 0.80 else 0.0

  add_cand("M0_Constant_Phase", 2L, .eval_vm_profile(theta_obs, rep(0.0, n)))
  add_cand("M1_Hopkins_1918_Fixed", 2L, .eval_vm_profile(theta_obs, 4.00 * lat_c + be_m1 * elev_c + blo_m1 * lon_c),
           bl = 4.00, blo = blo_m1, be = be_m1)

  b_lat_m2 <- 0.0
  lat_sd <- stats::sd(lat_v, na.rm = TRUE)
  if (n >= 10L && is.finite(lat_sd) && lat_sd > 0.05) {
    base_off <- blo_m1 * lon_c + be_m1 * elev_c
    g_lat <- seq(-15.0, 25.0, by = 0.25)
    b_init <- g_lat[which.max(vapply(g_lat, function(bl) .eval_vm_profile(theta_obs, bl * lat_c + base_off)$loglik, numeric(1)))]
    opt2 <- stats::optimize(
      function(bl) .eval_vm_profile(theta_obs, bl * lat_c + base_off)$loglik,
      interval = c(max(-15.0, b_init - 1.0), min(25.0, b_init + 1.0)), maximum = TRUE
    )
    b_lat_m2 <- as.numeric(opt2$maximum)
    add_cand("M2_Empirical_Lat", 3L, .eval_vm_profile(theta_obs, b_lat_m2 * lat_c + base_off), bl = b_lat_m2, blo = blo_m1, be = be_m1)
  }

  calc_aicc <- function(ll, kp) if (n > kp + 1L && is.finite(ll)) -2.0 * ll + 2.0 * kp + (2.0 * kp * (kp + 1.0)) / (n - kp - 1.0) else Inf
  fit_vm_sub <- function(fixed_off, cov_specs, k_val, starts_fn) {
    lbs <- c(-15.0, vapply(cov_specs, function(s) s$lb, numeric(1)))
    ubs <- c(25.0,  vapply(cov_specs, function(s) s$ub, numeric(1)))
    calc_off <- function(pc) {
      off <- fixed_off + pc[1L] * lat_c
      for (idx in seq_along(cov_specs)) off <- off + pc[idx + 1L] * cov_specs[[idx]]$x
      off
    }
    eval_fn <- function(par) {
      pc <- pmin(ubs, pmax(lbs, par))
      .eval_vm_profile(theta_obs, calc_off(pc))$loglik - 1e-4 * sum((par - pc)^2)
    }
    starts <- starts_fn(cov_specs)
    best_p <- pmin(ubs, pmax(lbs, starts[[1L]]))
    best_ll <- .eval_vm_profile(theta_obs, calc_off(best_p))$loglik
    for (st in starts) {
      res_o <- tryCatch(stats::optim(st, eval_fn, method = "BFGS", control = list(fnscale = -1, maxit = 200L)), error = function(e) NULL)
      if (!is.null(res_o) && is.finite(res_o$value)) {
        pc <- pmin(ubs, pmax(lbs, res_o$par))
        ll_ex <- .eval_vm_profile(theta_obs, calc_off(pc))$loglik
        if (is.finite(ll_ex) && ll_ex > best_ll) { best_ll <- ll_ex; best_p <- pc }
      }
    }
    b_ex <- c(elev = be_m1, mat = 0.0, year = 0.0)
    for (idx in seq_along(cov_specs)) b_ex[cov_specs[[idx]]$nm] <- best_p[idx + 1L]
    list(k = as.integer(k_val), prf = .eval_vm_profile(theta_obs, calc_off(best_p)), aicc = calc_aicc(best_ll, k_val),
         bl = best_p[1L], blo = blo_m1, be = b_ex["elev"], bm = b_ex["mat"], by = b_ex["year"])
  }

  if (n >= 18L && is.finite(lat_sd) && lat_sd > 0.05 && (has_elev || has_mat || has_yr)) {
    yr_sp <- list(nm = "year", x = yr_c, lb = -1.5, ub = 1.5, p_init = -0.08)
    tc_sp <- if (has_elev) list(nm = "elev", x = elev_c, lb = -0.05, ub = 0.10, p_init = 4.0 / 121.92) else if (has_mat) list(nm = "mat", x = mat_c, lb = -15.0, ub = 5.0, p_init = -2.0) else NULL
    sub_vars <- list()
    if (has_yr && n > 5L) sub_vars[[1L]] <- fit_vm_sub(blo_m1 * lon_c + be_m1 * elev_c, list(yr_sp), 4L, function(cs) list(c(b_lat_m2, 0.0), c(4.0, -0.08)))
    if (n >= 25L && !is.null(tc_sp)) {
      cs_3b <- if (has_yr) list(tc_sp, yr_sp) else list(tc_sp)
      sub_vars[[length(sub_vars) + 1L]] <- fit_vm_sub(blo_m1 * lon_c, cs_3b, 3L + length(cs_3b), function(cs) list(c(b_lat_m2, rep(0.0, length(cs))), c(4.0, vapply(cs, function(s) s$p_init, numeric(1)))))
    }
    if (length(sub_vars) > 0L) {
      w3 <- sub_vars[[which.min(vapply(sub_vars, function(v) v$aicc, numeric(1)))]]
      add_cand("M3_Empirical_TopoClimate_Year", w3$k, w3$prf, bl = w3$bl, blo = w3$blo, be = w3$be, bm = w3$bm, by = w3$by)
    }
  }

  aiccs <- vapply(m_list, function(m) {
    if (n > m$k + 1L && is.finite(m$loglik)) -2.0 * m$loglik + 2.0 * m$k + (2.0 * m$k * (m$k + 1.0)) / (n - m$k - 1.0) else Inf
  }, numeric(1))
  fin_aicc <- is.finite(aiccs)
  best_idx <- if (any(fin_aicc)) which.min(aiccs) else 1L
  min_aicc <- if (any(fin_aicc)) min(aiccs[fin_aicc]) else NA_real_
  delta_aicc <- if (is.finite(min_aicc)) aiccs - min_aicc else rep(NA_real_, length(m_list))
  raw_w <- ifelse(is.finite(delta_aicc), exp(-0.5 * pmin(delta_aicc, 500.0)), 0.0)
  aw <- if (sum(raw_w) > 0.0) raw_w / sum(raw_w) else rep(NA_real_, length(m_list))
  m1_idx <- which(vapply(m_list, function(m) m$name == "M1_Hopkins_1918_Fixed", logical(1)))[1L]
  aicc_m1 <- if (!is.na(m1_idx)) aiccs[m1_idx] else NA_real_
  d_m1 <- if (is.finite(aicc_m1)) aiccs - aicc_m1 else rep(NA_real_, length(m_list))

  sel <- m_list[[best_idx]]
  sel_tbl <- data.frame(
    Model = vapply(m_list, function(m) m$name, character(1)),
    k_params = vapply(m_list, function(m) m$k, integer(1)),
    N = rep(as.integer(n), length(m_list)),
    LogLik = round(vapply(m_list, function(m) m$loglik, numeric(1)), 4),
    Kappa = round(vapply(m_list, function(m) m$kappa, numeric(1)), 4),
    Slope_Lat_d_per_deg = round(vapply(m_list, function(m) m$b_lat, numeric(1)), 4),
    AICc = round(aiccs, 4), Delta_AICc = round(delta_aicc, 4),
    Akaike_Weight = round(aw, 6), Delta_AICc_vs_Hopkins1918 = round(d_m1, 4),
    stringsAsFactors = FALSE
  )

  rad_per_day <- 2.0 * pi / 365.25
  predict_doy_rad <- function(newdata_df) {
    nd <- as.data.frame(newdata_df, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(nd) == 0L) return(numeric(0))
    qlat <- .extract_circ_cov(nd, c("decimalLatitude", "latitude"), -90.0, 90.0)
    qlon <- .extract_circ_cov(nd, c("decimalLongitude", "longitude"), -180.0, 180.0)
    qelev <- .extract_circ_cov(nd, c("elevation_m", "elevation", "minimumElevationInMeters"), -500.0, 9000.0)
    qmat <- .extract_circ_cov(nd, c("bio1_temp_mean", "mat"), -60.0, 60.0)
    qyr <- .extract_circ_cov(nd, c("year"), 1750.0, 2100.0)
    off_q <- sel$b_lat * (qlat - lat_ref) +
             ifelse(is.finite(qlon), sel$b_lon * (qlon - lon_ref), 0.0) +
             ifelse(is.finite(qelev), sel$b_elev * (qelev - elev_ref), 0.0) +
             ifelse(is.finite(qmat), sel$b_mat * (qmat - mat_ref), 0.0) +
             ifelse(is.finite(qyr), sel$b_year * (qyr - yr_ref), 0.0)
    ifelse(is.finite(qlat), wrap_rad(sel$mu_0 + rad_per_day * off_q), NA_real_)
  }

  list(
    selected_model = sel$name, n = as.integer(n), mu_0 = sel$mu_0, kappa = sel$kappa,
    b_lat = sel$b_lat, b_lon = sel$b_lon, b_elev = sel$b_elev, b_mat = sel$b_mat, b_year = sel$b_year,
    lat_ref = lat_ref, lon_ref = lon_ref, elev_ref = elev_ref, mat_ref = mat_ref, year_ref = yr_ref,
    rayleigh_z = ray$Z, rayleigh_p = ray$p_value, r_bar = ray$R_bar, is_uniform_null = is_uniform_null,
    aicc_weight = as.numeric(aw[best_idx]), delta_aicc_vs_hopkins1918 = as.numeric(d_m1[best_idx]),
    model_selection_table = sel_tbl, predict_doy_rad = predict_doy_rad
  )
}

run_circular_multigroup_phase_test <- function(theta_vec, group_vec, B = 999L, seed = 42L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))

  th_raw <- suppressWarnings(as.numeric(theta_vec))
  gr_raw <- as.character(group_vec)
  keep <- is.finite(th_raw) & !is.na(gr_raw) & nzchar(gr_raw)
  th_f <- th_raw[keep]; gr_f <- gr_raw[keep]
  tab <- table(gr_f); valid_grps <- names(tab)[tab >= 3L]
  keep_g <- gr_f %in% valid_grps
  th <- th_f[keep_g]; grp <- factor(gr_f[keep_g], levels = valid_grps)
  s <- nlevels(grp); n_tot <- length(th)
  if (s < 2L || n_tot < (s + 2L)) {
    return(list(watson_f = NA_real_, watson_p = NA_real_, perm_p = NA_real_, n_groups = as.integer(s), n_obs = as.integer(n_tot), groups = valid_grps))
  }

  mat_cs <- cbind(cos(th), sin(th))
  g_idx <- as.integer(grp)
  rs_obs <- rowsum(mat_cs, g_idx, reorder = FALSE)
  sum_rg <- sum(sqrt(rs_obs[, 1L]^2 + rs_obs[, 2L]^2))
  r_tot <- sqrt(sum(mat_cs[, 1L])^2 + sum(mat_cs[, 2L])^2)
  kap_w <- estimate_von_mises_kappa(min(1.0 - 1e-8, max(0.0, sum_rg / n_tot)), n_tot)
  if (kap_w > 1e-4 && (n_tot - sum_rg) > 1e-10) {
    watson_f <- max(0.0, (1.0 + 3.0 / (8.0 * kap_w)) * ((n_tot - s) * (sum_rg - r_tot)) / ((s - 1.0) * (n_tot - sum_rg)))
    watson_p <- stats::pf(watson_f, df1 = s - 1.0, df2 = n_tot - s, lower.tail = FALSE)
  } else {
    watson_f <- 0.0; watson_p <- 1.0
  }

  obs_stat <- sum_rg - r_tot
  B_int <- max(1L, as.integer(B))
  perms <- replicate(B_int, sample.int(n_tot), simplify = FALSE)
  eval_idx <- function(idx_v) sum(vapply(idx_v, function(b) {
    rs_b <- rowsum(mat_cs, g_idx[perms[[b]]], reorder = FALSE)
    st <- sum(sqrt(rs_b[, 1L]^2 + rs_b[, 2L]^2)) - r_tot
    is.finite(st) && st >= (obs_stat - 1e-12)
  }, logical(1)))
  nw_c <- if (B_int >= 99L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, ceiling(B_int / 16L)) else 1L
  n_ge <- if (nw_c > 1L) sum(unlist(packera_mclapply(parallel::splitIndices(B_int, nw_c), eval_idx, mc.cores = nw_c), use.names = FALSE)) else eval_idx(seq_len(B_int))
  perm_p <- (n_ge + 1.0) / (B_int + 1.0)
  list(watson_f = watson_f, watson_p = watson_p, perm_p = perm_p, n_groups = as.integer(s), n_obs = as.integer(n_tot), groups = valid_grps)
}
