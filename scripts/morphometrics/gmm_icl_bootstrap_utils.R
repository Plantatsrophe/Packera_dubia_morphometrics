# Non-Gaussian Tail & Over-Splitting Guards: ICL Profiling and Parametric Bootstrap LRT
# Protocol Reference: UNC-BOT-SOP-2026-04-REV5 (Biernacki, Celeux & Govaert 2000; McLachlan 1987)

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.save_random_seed <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(had_seed = had_seed, old_seed = if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL)
}

.restore_random_seed <- function(seed_state) {
  if (isTRUE(seed_state$had_seed)) assign(".Random.seed", seed_state$old_seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
}

.ensure_mclust_attached <- function() {
  has_mclust <- requireNamespace("mclust", quietly = TRUE)
  if (has_mclust && !("package:mclust" %in% search())) tryCatch(suppressPackageStartupMessages(attachNamespace("mclust")), error = function(e) NULL)
  has_mclust
}

.kmeans_soft_responsibilities <- function(X_mat, k_val, seed_i) {
  n_r <- nrow(X_mat); set.seed(seed_i)
  km <- tryCatch(stats::kmeans(X_mat, centers = k_val, nstart = 5L), error = function(e) NULL)
  if (is.null(km) || is.null(km$centers)) return(matrix(1.0 / k_val, nrow = n_r, ncol = k_val))
  d2 <- matrix(0.0, nrow = n_r, ncol = k_val)
  for (j in seq_len(k_val)) d2[, j] <- rowSums(sweep(X_mat, 2L, km$centers[j, ], "-")^2)
  w <- exp(-0.5 * sweep(d2, 1L, apply(d2, 1L, min), "-")); rs <- rowSums(w); rs[!is.finite(rs) | rs <= 0.0] <- 1.0
  w / rs
}

compute_gmm_icl_profile <- function(X_core, bic_table, modelNames = NULL, pure_res = NULL, seed = 42L) {
  seed_state <- .save_random_seed(); on.exit(.restore_random_seed(seed_state), add = TRUE)
  X_core <- as.matrix(X_core); storage.mode(X_core) <- "double"; n_core <- nrow(X_core)
  seed_i <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  k_vec <- as.integer(bic_table$K); bic_vec <- as.numeric(bic_table$BIC)
  icl_vec <- rep(NA_real_, length(k_vec)); ent_vec <- rep(NA_real_, length(k_vec))
  bic_model_cols <- grep("^BIC_", names(bic_table), value = TRUE); has_mclust <- .ensure_mclust_attached()
  for (i in seq_along(k_vec)) {
    k_val <- k_vec[i]; bic_val <- bic_vec[i]; if (!is.finite(k_val) || !is.finite(bic_val)) next
    if (k_val <= 1L) { ent_vec[i] <- 0.0; icl_vec[i] <- bic_val; next }
    win_model_k <- modelNames
    if (length(bic_model_cols) > 0L) {
      row_bics <- as.numeric(unlist(bic_table[i, bic_model_cols, drop = TRUE])); row_bics[!is.finite(row_bics)] <- NA_real_
      if (any(!is.na(row_bics))) win_model_k <- sub("^BIC_", "", bic_model_cols[which.max(row_bics)])
    }
    z_raw <- NULL
    if (has_mclust) {
      set.seed(seed_i); fit_k <- tryCatch(mclust::Mclust(X_core, G = k_val, modelNames = win_model_k, verbose = FALSE), error = function(e) NULL)
      if (is.null(fit_k) && !identical(win_model_k, modelNames)) { set.seed(seed_i); fit_k <- tryCatch(mclust::Mclust(X_core, G = k_val, modelNames = modelNames, verbose = FALSE), error = function(e) NULL) }
      if (!is.null(fit_k) && !is.null(fit_k$z)) z_raw <- fit_k$z
    }
    if (is.null(z_raw) && !is.null(pure_res) && !is.null(pure_res$best) && identical(as.integer(pure_res$best$k), k_val) && !is.null(pure_res$best$z) && nrow(pure_res$best$z) == n_core) z_raw <- pure_res$best$z
    if (is.null(z_raw) && exists("fit_gmm_em_pure", mode = "function")) {
      pure_k <- tryCatch(fit_gmm_em_pure(X_core, max_k = k_val, seed = seed_i), error = function(e) NULL)
      if (!is.null(pure_k) && !is.null(pure_k$best) && identical(as.integer(pure_k$best$k), k_val) && nrow(pure_k$best$z) == n_core) z_raw <- pure_k$best$z
    }
    if (is.null(z_raw)) z_raw <- .kmeans_soft_responsibilities(X_core, k_val, seed_i)
    z_mat <- matrix(as.numeric(z_raw), nrow = n_core, ncol = k_val); z_mat[!is.finite(z_mat)] <- 1.0 / k_val
    ent_k <- -sum(z_mat * log(pmax(z_mat, 1e-12))); ent_vec[i] <- ent_k; icl_vec[i] <- bic_val - 2.0 * ent_k
  }
  idx_1 <- match(1L, k_vec); icl_1 <- if (!is.na(idx_1) && is.finite(icl_vec[idx_1])) icl_vec[idx_1] else icl_vec[1L]
  data.frame(K = bic_table$K, ICL = round(icl_vec, 2), Delta_ICL_vs_Null = round(icl_vec - icl_1, 2), Classification_Entropy = round(ent_vec, 3), stringsAsFactors = FALSE)
}

.fit_g_pair <- function(X_mat, g_null, g_alt, mn_null, mn_alt, has_mclust, fit_seed) {
  g_null <- as.integer(g_null); g_alt <- as.integer(g_alt)
  if (has_mclust) {
    fit_one <- function(G_v, mn_v) {
      set.seed(fit_seed); f <- tryCatch(mclust::Mclust(X_mat, G = G_v, modelNames = mn_v, verbose = FALSE), error = function(e) NULL)
      if (is.null(f) && !is.null(mn_v)) { set.seed(fit_seed); f <- tryCatch(mclust::Mclust(X_mat, G = G_v, modelNames = NULL, verbose = FALSE), error = function(e) NULL) }
      f
    }
    f_null <- fit_one(g_null, mn_null); f_alt <- fit_one(g_alt, mn_alt)
    if (!is.null(f_null) && !is.null(f_alt) && is.finite(as.numeric(f_null$loglik[1L])) && is.finite(as.numeric(f_alt$loglik[1L]))) {
      return(list(fit_null = f_null, fit_alt = f_alt, ok = TRUE))
    }
  }
  if (exists("fit_gmm_em_pure", mode = "function")) {
    pe <- tryCatch(fit_gmm_em_pure(X_mat, max_k = max(g_null, g_alt), seed = fit_seed), error = function(e) NULL)
    if (!is.null(pe) && !is.null(pe$bic_table)) {
      ll1 <- as.numeric(pe$bic_table$LogLik[match(g_null, pe$bic_table$K)]); ll2 <- as.numeric(pe$bic_table$LogLik[match(g_alt, pe$bic_table$K)])
      if (is.finite(ll1) && is.finite(ll2)) {
        p_d <- ncol(X_mat); cov_x <- if (p_d == 1L) matrix(stats::var(as.numeric(X_mat)), 1L, 1L) else stats::cov(X_mat)
        f_null <- if (g_null == 1L) list(loglik = ll1, modelName = "pure_EM", parameters = list(pro = 1.0, mean = colMeans(X_mat), variance = list(sigma = cov_x))) else {
          set.seed(fit_seed); km <- tryCatch(stats::kmeans(X_mat, centers = g_null, nstart = 5L), error = function(e) NULL)
          mu_m <- if (!is.null(km)) t(km$centers) else matrix(rep(colMeans(X_mat), g_null), nrow = p_d, ncol = g_null)
          pro_v <- if (!is.null(km)) as.numeric(table(factor(km$cluster, levels = seq_len(g_null)))) / nrow(X_mat) else rep(1.0 / g_null, g_null)
          sig_arr <- array(0.0, dim = c(p_d, p_d, g_null)); for (kc in seq_len(g_null)) sig_arr[, , kc] <- cov_x
          list(loglik = ll1, modelName = "pure_EM", parameters = list(pro = pro_v, mean = mu_m, variance = list(sigma = sig_arr)))
        }
        return(list(fit_null = f_null, fit_alt = list(loglik = ll2, modelName = "pure_EM"), ok = TRUE))
      }
    }
  }
  list(fit_null = NULL, fit_alt = NULL, ok = FALSE)
}

run_gmm_pair_bootstrap_lrt <- function(X_core, g_null = 1L, g_alt = 2L, modelNames = NULL, B = 199L, seed = 42L) {
  seed_state <- .save_random_seed(); on.exit(.restore_random_seed(seed_state), add = TRUE)
  X_core <- as.matrix(X_core); storage.mode(X_core) <- "double"; n_c <- nrow(X_core); p_c <- ncol(X_core); B <- max(1L, as.integer(B))
  g_null <- max(1L, as.integer(g_null)); g_alt <- max(g_null + 1L, as.integer(g_alt))
  seed_base <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  if (n_c < max(3L, g_alt + 1L) || p_c < 1L) return(list(lrt_stat = 0.0, p_value = 1.0, n_boot_valid = 0L, null_model = "None", alt_model = "None"))
  has_mclust <- .ensure_mclust_attached(); mn_use <- if (p_c == 1L) NULL else modelNames
  obs_pair <- .fit_g_pair(X_core, g_null, g_alt, mn_use, mn_use, has_mclust, seed_base)
  if (!isTRUE(obs_pair$ok)) return(list(lrt_stat = 0.0, p_value = 1.0, n_boot_valid = 0L, null_model = "None", alt_model = "None"))
  fit_null <- obs_pair$fit_null; fit_alt <- obs_pair$fit_alt
  obs_lrt <- max(0.0, 2.0 * (as.numeric(fit_alt$loglik) - as.numeric(fit_null$loglik)))
  pro_raw <- if (!is.null(fit_null$parameters$pro)) as.numeric(fit_null$parameters$pro) else rep(1.0 / g_null, g_null)
  if (length(pro_raw) != g_null || any(!is.finite(pro_raw)) || sum(pmax(pro_raw, 0.0)) <= 0.0) pro_raw <- rep(1.0 / g_null, g_null)
  pro_0 <- pmax(pro_raw, 0.0); pro_0 <- pro_0 / sum(pro_0)
  mean_raw <- fit_null$parameters$mean
  mu_mat <- if (!is.null(mean_raw) && length(as.numeric(mean_raw)) == p_c * g_null && all(is.finite(as.numeric(mean_raw)))) matrix(as.numeric(mean_raw), nrow = p_c, ncol = g_null) else matrix(rep(colMeans(X_core), g_null), nrow = p_c, ncol = g_null)
  sig_raw <- if (!is.null(fit_null$parameters$variance$sigma)) fit_null$parameters$variance$sigma else NULL
  sq_raw <- if (!is.null(fit_null$parameters$variance$sigmasq)) as.numeric(fit_null$parameters$variance$sigmasq) else NULL
  cov_fb <- if (p_c == 1L) matrix(max(stats::var(as.numeric(X_core)), 1e-6), 1L, 1L) else stats::cov(X_core); cov_fb <- 0.5 * (cov_fb + t(cov_fb))
  chol_list <- lapply(seq_len(g_null), function(k) {
    Sig_k <- if (!is.null(sig_raw) && length(dim(sig_raw)) == 3L && dim(sig_raw)[3L] >= k && identical(dim(sig_raw)[1:2], c(p_c, p_c))) matrix(as.numeric(sig_raw[, , k]), nrow = p_c, ncol = p_c)
    else if (!is.null(sig_raw) && is.matrix(sig_raw) && identical(dim(sig_raw), c(p_c, p_c))) matrix(as.numeric(sig_raw), nrow = p_c, ncol = p_c)
    else if (p_c == 1L && !is.null(sq_raw) && length(sq_raw) >= 1L) matrix(as.numeric(sq_raw[min(k, length(sq_raw))]), 1L, 1L) else cov_fb
    if (any(!is.finite(Sig_k))) Sig_k <- cov_fb; Sig_k <- 0.5 * (Sig_k + t(Sig_k))
    tryCatch(chol(Sig_k + diag(1e-8, p_c)), error = function(e) tryCatch(chol(cov_fb + diag(1e-4, p_c)), error = function(e2) diag(1.0, p_c)))
  })
  mn_boot_null <- if (!is.null(fit_null$modelName) && nzchar(fit_null$modelName) && !identical(fit_null$modelName, "pure_EM")) fit_null$modelName else mn_use
  mn_boot_alt <- if (!is.null(fit_alt$modelName) && nzchar(fit_alt$modelName) && !identical(fit_alt$modelName, "pure_EM")) fit_alt$modelName else mn_use
  boot_worker <- function(b) {
    b_seed <- seed_base + as.integer(b); set.seed(b_seed)
    if (g_null == 1L) X_b <- sweep(matrix(stats::rnorm(n_c * p_c), nrow = n_c, ncol = p_c) %*% chol_list[[1L]], 2L, mu_mat[, 1L], "+") else {
      comp_b <- sample.int(g_null, size = n_c, replace = TRUE, prob = pro_0); X_b <- matrix(0.0, nrow = n_c, ncol = p_c)
      for (k in seq_len(g_null)) { ik <- which(comp_b == k); if (length(ik) > 0L) X_b[ik, ] <- sweep(matrix(stats::rnorm(length(ik) * p_c), nrow = length(ik), ncol = p_c) %*% chol_list[[k]], 2L, mu_mat[, k], "+") }
    }
    bp <- tryCatch(.fit_g_pair(X_b, g_null, g_alt, mn_boot_null, mn_boot_alt, has_mclust, b_seed), error = function(e) list(ok = FALSE))
    if (!isTRUE(bp$ok)) return(NA_real_)
    ll1 <- as.numeric(bp$fit_null$loglik); ll2 <- as.numeric(bp$fit_alt$loglik)
    if (!is.finite(ll1) || !is.finite(ll2)) NA_real_ else max(0.0, 2.0 * (ll2 - ll1))
  }
  lrt_list <- if (exists("packera_mclapply", mode = "function")) packera_mclapply(seq_len(B), boot_worker, mc.cores = if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(16L, B) else 1L) else lapply(seq_len(B), boot_worker)
  lrt_boot <- as.numeric(unlist(lrt_list)); n_valid <- sum(is.finite(lrt_boot))
  boot_p <- if (n_valid > 0L) (1.0 + sum(lrt_boot >= obs_lrt, na.rm = TRUE)) / (n_valid + 1.0) else 1.0
  list(lrt_stat = round(obs_lrt, 3), p_value = round(boot_p, 4), n_boot_valid = as.integer(n_valid), null_model = as.character(fit_null$modelName), alt_model = as.character(fit_alt$modelName))
}

run_gmm_k1_vs_k2_bootstrap_lrt <- function(X_core, modelNames = NULL, B = 199L, seed = 42L) {
  run_gmm_pair_bootstrap_lrt(X_core, g_null = 1L, g_alt = 2L, modelNames = modelNames, B = B, seed = seed)
}

.refit_gmm_stepped <- function(X_core, X_pcs, k_sup, modelNames, seed) {
  seed_state <- .save_random_seed(); on.exit(.restore_random_seed(seed_state), add = TRUE)
  seed_i <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  if (.ensure_mclust_attached()) {
    set.seed(seed_i); fit_s <- tryCatch(mclust::Mclust(X_core, G = k_sup, modelNames = modelNames, verbose = FALSE), error = function(e) NULL)
    if (is.null(fit_s) && !is.null(modelNames)) { set.seed(seed_i); fit_s <- tryCatch(mclust::Mclust(X_core, G = k_sup, modelNames = NULL, verbose = FALSE), error = function(e) NULL) }
    if (!is.null(fit_s)) { pr_s <- stats::predict(fit_s, newdata = X_pcs); return(list(cls = as.integer(pr_s$classification), unc = if (!is.null(pr_s$uncertainty)) as.numeric(pr_s$uncertainty) else 1.0 - apply(as.matrix(pr_s$z), 1L, max))) }
  }
  if (exists("fit_gmm_em_pure", mode = "function")) {
    pe_s <- tryCatch(fit_gmm_em_pure(X_core, max_k = k_sup, X_eval = X_pcs, seed = seed_i), error = function(e) NULL)
    if (!is.null(pe_s) && !is.null(pe_s$best) && length(pe_s$best$classification) == nrow(X_pcs)) return(list(cls = as.integer(pe_s$best$classification), unc = as.numeric(pe_s$best$uncertainty)))
  }
  set.seed(seed_i); km <- tryCatch(stats::kmeans(X_core, centers = k_sup, nstart = 5L), error = function(e) NULL)
  if (!is.null(km) && !is.null(km$centers)) {
    d2 <- vapply(seq_len(k_sup), function(j) rowSums(sweep(X_pcs, 2L, km$centers[j, ], "-")^2), numeric(nrow(X_pcs)))
    return(list(cls = apply(d2, 1L, which.min), unc = rep(0.0, nrow(X_pcs))))
  }
  list(cls = rep(1L, nrow(X_pcs)), unc = rep(0.0, nrow(X_pcs)))
}

apply_gmm_icl_bootstrap_gate <- function(X_core, X_pcs, bic_table, best_k, best_name, cls_sub, unc_sub,
                                         modelNames = NULL, pure_res = NULL, is_intra_gmm = FALSE, seed = 42L) {
  icl_prof <- compute_gmm_icl_profile(X_core, bic_table, modelNames = modelNames, pure_res = pure_res, seed = seed)
  bic_table$ICL <- icl_prof$ICL; bic_table$Delta_ICL_vs_Null <- icl_prof$Delta_ICL_vs_Null
  bic_table$Boot_LRT_Stat_K1_vs_K2 <- NA_real_; bic_table$Boot_LRT_P_K1_vs_K2 <- NA_real_
  bic_table$Boot_LRT_Stat_prev <- NA_real_; bic_table$Boot_LRT_P_prev <- NA_real_
  boot_p <- NA_real_; did_icl_collapse <- FALSE
  idx_1 <- match(1L, bic_table$K); idx_2 <- match(2L, bic_table$K)
  if (!is.na(idx_1) && !is.na(idx_2) && is.finite(bic_table$BIC[idx_1]) && is.finite(bic_table$BIC[idx_2])) {
    b12 <- run_gmm_pair_bootstrap_lrt(X_core, g_null = 1L, g_alt = 2L, modelNames = modelNames, B = 199L, seed = seed); boot_p <- b12$p_value
    bic_table$Boot_LRT_Stat_K1_vs_K2[idx_1] <- 0.0; bic_table$Boot_LRT_P_K1_vs_K2[idx_1] <- 1.0
    bic_table$Boot_LRT_Stat_K1_vs_K2[idx_2] <- b12$lrt_stat; bic_table$Boot_LRT_P_K1_vs_K2[idx_2] <- boot_p
    bic_table$Boot_LRT_Stat_prev[idx_2] <- b12$lrt_stat; bic_table$Boot_LRT_P_prev[idx_2] <- boot_p
  }
  fails_gate <- function(delta_icl, bp) if (isTRUE(is_intra_gmm)) !(is.finite(delta_icl) && delta_icl > 0.0 && is.finite(bp) && bp < 0.05) else is.finite(delta_icl) && delta_icl <= 0.0 && is.finite(bp) && bp >= 0.05
  if (best_k >= 2L) {
    delta_icl_best <- as.numeric(bic_table$Delta_ICL_vs_Null[match(best_k, bic_table$K)])
    if (fails_gate(delta_icl_best, boot_p)) {
      message(sprintf("[Model Selection Gate] Collapsed GMM K=%d -> K=1 (Delta_ICL=%.2f, Boot_LRT_p=%.4f).", as.integer(best_k), delta_icl_best, as.numeric(boot_p)))
      best_k <- 1L; cls_sub <- rep(1L, nrow(X_pcs)); unc_sub <- rep(0.0, nrow(X_pcs)); best_name <- paste0(best_name, "_ICLBootCollapsed"); did_icl_collapse <- TRUE
    }
  }
  if (!did_icl_collapse && best_k > 2L) {
    k_orig <- as.integer(best_k); k_sup <- 2L; chain_ok <- TRUE
    for (k_seq in seq.int(3L, k_orig)) {
      rk <- match(k_seq, bic_table$K); rp <- match(k_seq - 1L, bic_table$K)
      if (is.na(rk) || is.na(rp) || !is.finite(bic_table$BIC[rk]) || !is.finite(bic_table$BIC[rp])) { chain_ok <- FALSE; next }
      b_seq <- run_gmm_pair_bootstrap_lrt(X_core, g_null = k_seq - 1L, g_alt = k_seq, modelNames = modelNames, B = 199L, seed = seed)
      bic_table$Boot_LRT_Stat_prev[rk] <- b_seq$lrt_stat; bic_table$Boot_LRT_P_prev[rk] <- b_seq$p_value
      delta_icl_prev <- as.numeric(bic_table$ICL[rk]) - as.numeric(bic_table$ICL[rp])
      if (chain_ok && !fails_gate(delta_icl_prev, b_seq$p_value)) k_sup <- k_seq else chain_ok <- FALSE
    }
    if (k_sup < k_orig) {
      message(sprintf("[Model Selection Gate] Sequential bootstrap LRT stepped GMM K=%d -> K=%d.", k_orig, k_sup))
      best_k <- as.integer(k_sup); best_name <- paste0(best_name, "_SeqBootStepped")
      rf <- .refit_gmm_stepped(X_core, X_pcs, k_sup, modelNames, seed); cls_sub <- rf$cls; unc_sub <- rf$unc
    }
  }
  list(bic_table = bic_table, best_k = best_k, best_name = best_name, cls_sub = cls_sub, unc_sub = unc_sub, did_icl_collapse = did_icl_collapse)
}
