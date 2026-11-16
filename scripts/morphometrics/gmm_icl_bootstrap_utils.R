# Non-Gaussian Tail & Over-Splitting Guards: ICL Profiling and Parametric Bootstrap LRT
# Protocol Reference: UNC-BOT-SOP-2026-04-REV5 (Biernacki, Celeux & Govaert 2000; McLachlan 1987)

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.save_random_seed <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(had_seed = had_seed, old_seed = if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL)
}

.restore_random_seed <- function(seed_state) {
  if (isTRUE(seed_state$had_seed)) {
    assign(".Random.seed", seed_state$old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

.ensure_mclust_attached <- function() {
  has_mclust <- requireNamespace("mclust", quietly = TRUE)
  if (has_mclust && !("package:mclust" %in% search())) {
    tryCatch(suppressPackageStartupMessages(attachNamespace("mclust")), error = function(e) NULL)
  }
  has_mclust
}

.kmeans_soft_responsibilities <- function(X_mat, k_val, seed_i) {
  n_r <- nrow(X_mat); set.seed(seed_i)
  km <- tryCatch(stats::kmeans(X_mat, centers = k_val, nstart = 5L), error = function(e) NULL)
  if (is.null(km) || is.null(km$centers)) return(matrix(1.0 / k_val, nrow = n_r, ncol = k_val))
  d2 <- matrix(0.0, nrow = n_r, ncol = k_val)
  for (j in seq_len(k_val)) d2[, j] <- rowSums(sweep(X_mat, 2L, km$centers[j, ], "-")^2)
  w <- exp(-0.5 * sweep(d2, 1L, apply(d2, 1L, min), "-"))
  rs <- rowSums(w); rs[!is.finite(rs) | rs <= 0.0] <- 1.0
  w / rs
}

compute_gmm_icl_profile <- function(X_core, bic_table, modelNames = NULL, pure_res = NULL, seed = 42L) {
  seed_state <- .save_random_seed()
  on.exit(.restore_random_seed(seed_state), add = TRUE)

  X_core <- as.matrix(X_core); storage.mode(X_core) <- "double"
  n_core <- nrow(X_core)
  seed_i <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  k_vec <- as.integer(bic_table$K); bic_vec <- as.numeric(bic_table$BIC)
  icl_vec <- rep(NA_real_, length(k_vec)); ent_vec <- rep(NA_real_, length(k_vec))
  bic_model_cols <- grep("^BIC_", names(bic_table), value = TRUE)
  has_mclust <- .ensure_mclust_attached()

  for (i in seq_along(k_vec)) {
    k_val <- k_vec[i]; bic_val <- bic_vec[i]
    if (!is.finite(k_val) || !is.finite(bic_val)) next
    if (k_val <= 1L) { ent_vec[i] <- 0.0; icl_vec[i] <- bic_val; next }

    win_model_k <- modelNames
    if (length(bic_model_cols) > 0L) {
      row_bics <- as.numeric(unlist(bic_table[i, bic_model_cols, drop = TRUE]))
      row_bics[!is.finite(row_bics)] <- NA_real_
      if (any(!is.na(row_bics))) win_model_k <- sub("^BIC_", "", bic_model_cols[which.max(row_bics)])
    }

    z_raw <- NULL
    if (has_mclust) {
      set.seed(seed_i)
      fit_k <- tryCatch(mclust::Mclust(X_core, G = k_val, modelNames = win_model_k, verbose = FALSE), error = function(e) NULL)
      if (is.null(fit_k) && !identical(win_model_k, modelNames)) {
        set.seed(seed_i)
        fit_k <- tryCatch(mclust::Mclust(X_core, G = k_val, modelNames = modelNames, verbose = FALSE), error = function(e) NULL)
      }
      if (!is.null(fit_k) && !is.null(fit_k$z)) z_raw <- fit_k$z
    }
    if (is.null(z_raw) && !is.null(pure_res) && !is.null(pure_res$best) &&
        identical(as.integer(pure_res$best$k), k_val) && !is.null(pure_res$best$z) && nrow(pure_res$best$z) == n_core) {
      z_raw <- pure_res$best$z
    }
    if (is.null(z_raw) && exists("fit_gmm_em_pure", mode = "function")) {
      pure_k <- tryCatch(fit_gmm_em_pure(X_core, max_k = k_val, seed = seed_i), error = function(e) NULL)
      if (!is.null(pure_k) && !is.null(pure_k$best) && identical(as.integer(pure_k$best$k), k_val) && nrow(pure_k$best$z) == n_core) {
        z_raw <- pure_k$best$z
      }
    }
    if (is.null(z_raw)) z_raw <- .kmeans_soft_responsibilities(X_core, k_val, seed_i)

    z_mat <- matrix(as.numeric(z_raw), nrow = n_core, ncol = k_val)
    z_mat[!is.finite(z_mat)] <- 1.0 / k_val
    ent_k <- -sum(z_mat * log(pmax(z_mat, 1e-12)))
    ent_vec[i] <- ent_k; icl_vec[i] <- bic_val - 2.0 * ent_k
  }

  idx_1 <- match(1L, k_vec)
  icl_1 <- if (!is.na(idx_1) && is.finite(icl_vec[idx_1])) icl_vec[idx_1] else icl_vec[1L]
  data.frame(
    K = bic_table$K, ICL = round(icl_vec, 2),
    Delta_ICL_vs_Null = round(icl_vec - icl_1, 2),
    Classification_Entropy = round(ent_vec, 3), stringsAsFactors = FALSE
  )
}

.fit_g1_g2_pair <- function(X_mat, mn_1, mn_2, has_mclust, fit_seed) {
  set.seed(fit_seed)
  if (has_mclust) {
    f1 <- tryCatch(mclust::Mclust(X_mat, G = 1L, modelNames = mn_1, verbose = FALSE), error = function(e) NULL)
    if (is.null(f1) && !is.null(mn_1)) { set.seed(fit_seed); f1 <- tryCatch(mclust::Mclust(X_mat, G = 1L, modelNames = NULL, verbose = FALSE), error = function(e) NULL) }
    set.seed(fit_seed)
    f2 <- tryCatch(mclust::Mclust(X_mat, G = 2L, modelNames = mn_2, verbose = FALSE), error = function(e) NULL)
    if (is.null(f2) && !is.null(mn_2)) { set.seed(fit_seed); f2 <- tryCatch(mclust::Mclust(X_mat, G = 2L, modelNames = NULL, verbose = FALSE), error = function(e) NULL) }
    if (!is.null(f1) && !is.null(f2) && !is.null(f1$loglik) && !is.null(f2$loglik) &&
        is.finite(as.numeric(f1$loglik[1L])) && is.finite(as.numeric(f2$loglik[1L]))) {
      return(list(fit1 = f1, fit2 = f2, ok = TRUE))
    }
  }
  if (exists("fit_gmm_em_pure", mode = "function")) {
    pe <- tryCatch(fit_gmm_em_pure(X_mat, max_k = 2L, seed = fit_seed), error = function(e) NULL)
    if (!is.null(pe) && !is.null(pe$bic_table) && nrow(pe$bic_table) >= 2L) {
      ll1 <- as.numeric(pe$bic_table$LogLik[match(1L, pe$bic_table$K)])
      ll2 <- as.numeric(pe$bic_table$LogLik[match(2L, pe$bic_table$K)])
      if (is.finite(ll1) && is.finite(ll2)) {
        f1 <- list(loglik = ll1, modelName = "pure_EM", parameters = list(mean = colMeans(X_mat), variance = list(sigma = stats::cov(X_mat))))
        return(list(fit1 = f1, fit2 = list(loglik = ll2, modelName = "pure_EM"), ok = TRUE))
      }
    }
  }
  list(fit1 = NULL, fit2 = NULL, ok = FALSE)
}

run_gmm_k1_vs_k2_bootstrap_lrt <- function(X_core, modelNames = NULL, B = 199L, seed = 42L) {
  seed_state <- .save_random_seed()
  on.exit(.restore_random_seed(seed_state), add = TRUE)

  X_core <- as.matrix(X_core); storage.mode(X_core) <- "double"
  n_c <- nrow(X_core); p_c <- ncol(X_core); B <- max(1L, as.integer(B))
  seed_base <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  if (n_c < 3L || p_c < 1L) return(list(lrt_stat = 0.0, p_value = 1.0, n_boot_valid = 0L, null_model = "None", alt_model = "None"))

  has_mclust <- .ensure_mclust_attached()
  mn_use <- if (p_c == 1L) NULL else modelNames
  obs_pair <- .fit_g1_g2_pair(X_core, mn_use, mn_use, has_mclust, seed_base)
  if (!isTRUE(obs_pair$ok)) return(list(lrt_stat = 0.0, p_value = 1.0, n_boot_valid = 0L, null_model = "None", alt_model = "None"))

  fit_g1 <- obs_pair$fit1; fit_g2 <- obs_pair$fit2
  obs_lrt <- max(0.0, 2.0 * (as.numeric(fit_g2$loglik) - as.numeric(fit_g1$loglik)))
  mu_raw <- if (!is.null(fit_g1$parameters$mean)) as.numeric(fit_g1$parameters$mean) else colMeans(X_core)
  mu_0 <- if (length(mu_raw) == p_c && all(is.finite(mu_raw))) mu_raw else colMeans(X_core)

  sig_raw <- if (!is.null(fit_g1$parameters$variance$sigma)) fit_g1$parameters$variance$sigma else NULL
  Sigma_0 <- if (length(dim(sig_raw)) == 3L) {
    matrix(as.numeric(sig_raw[, , 1L]), nrow = p_c, ncol = p_c)
  } else if (is.matrix(sig_raw)) {
    matrix(as.numeric(sig_raw), nrow = nrow(sig_raw), ncol = ncol(sig_raw))
  } else if (!is.null(fit_g1$parameters$variance$sigmasq) && p_c == 1L) {
    matrix(as.numeric(fit_g1$parameters$variance$sigmasq[1L]), 1L, 1L)
  } else NULL
  if (is.null(Sigma_0) || !identical(dim(Sigma_0), c(p_c, p_c)) || any(!is.finite(Sigma_0))) {
    Sigma_0 <- if (p_c == 1L) matrix(stats::var(as.numeric(X_core)), 1L, 1L) else stats::cov(X_core)
  }
  Sigma_0 <- 0.5 * (Sigma_0 + t(Sigma_0))
  chol_R <- tryCatch(chol(Sigma_0 + diag(1e-8, p_c)), error = function(e) {
    cov_fb <- if (p_c == 1L) matrix(stats::var(as.numeric(X_core)), 1L, 1L) else stats::cov(X_core)
    tryCatch(chol(cov_fb + diag(1e-4, p_c)), error = function(e2) diag(1.0, p_c))
  })

  mn_boot_1 <- if (n_c > 250L && !is.null(fit_g1$modelName) && nzchar(fit_g1$modelName) && !identical(fit_g1$modelName, "pure_EM")) fit_g1$modelName else mn_use
  mn_boot_2 <- if (n_c > 250L && !is.null(fit_g2$modelName) && nzchar(fit_g2$modelName) && !identical(fit_g2$modelName, "pure_EM")) fit_g2$modelName else mn_use

  boot_worker <- function(b) {
    b_seed <- seed_base + as.integer(b)
    set.seed(b_seed)
    Z_b <- matrix(stats::rnorm(n_c * p_c), nrow = n_c, ncol = p_c)
    X_b <- sweep(Z_b %*% chol_R, 2L, mu_0, "+")
    b_pair <- tryCatch(.fit_g1_g2_pair(X_b, mn_boot_1, mn_boot_2, has_mclust, b_seed), error = function(e) list(ok = FALSE))
    if (!isTRUE(b_pair$ok)) return(NA_real_)
    ll1_b <- as.numeric(b_pair$fit1$loglik); ll2_b <- as.numeric(b_pair$fit2$loglik)
    if (!is.finite(ll1_b) || !is.finite(ll2_b)) return(NA_real_)
    max(0.0, 2.0 * (ll2_b - ll1_b))
  }

  lrt_list <- if (exists("packera_mclapply", mode = "function")) {
    n_workers <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(16L, B) else 1L
    packera_mclapply(seq_len(B), boot_worker, mc.cores = n_workers)
  } else {
    lapply(seq_len(B), boot_worker)
  }
  lrt_boot <- as.numeric(unlist(lrt_list))
  n_valid <- sum(is.finite(lrt_boot))
  boot_p <- if (n_valid > 0L) (1.0 + sum(lrt_boot >= obs_lrt, na.rm = TRUE)) / (n_valid + 1.0) else 1.0

  list(lrt_stat = round(obs_lrt, 3), p_value = round(boot_p, 4), n_boot_valid = as.integer(n_valid),
       null_model = as.character(fit_g1$modelName), alt_model = as.character(fit_g2$modelName))
}

apply_gmm_icl_bootstrap_gate <- function(X_core, X_pcs, bic_table, best_k, best_name, cls_sub, unc_sub,
                                         modelNames = NULL, pure_res = NULL, is_intra_gmm = FALSE, seed = 42L) {
  icl_prof <- compute_gmm_icl_profile(X_core, bic_table, modelNames = modelNames, pure_res = pure_res, seed = seed)
  bic_table$ICL <- icl_prof$ICL
  bic_table$Delta_ICL_vs_Null <- icl_prof$Delta_ICL_vs_Null

  boot_stat <- NA_real_; boot_p <- NA_real_; did_icl_collapse <- FALSE
  bic_table$Boot_LRT_Stat_K1_vs_K2 <- NA_real_
  bic_table$Boot_LRT_P_K1_vs_K2 <- NA_real_

  if (nrow(bic_table) >= 2L && is.finite(bic_table$BIC[1L]) && is.finite(bic_table$BIC[2L])) {
    boot_res <- run_gmm_k1_vs_k2_bootstrap_lrt(X_core, modelNames = modelNames, B = 199L, seed = seed)
    boot_stat <- boot_res$lrt_stat; boot_p <- boot_res$p_value
    bic_table$Boot_LRT_Stat_K1_vs_K2[bic_table$K == 1L] <- 0.0
    bic_table$Boot_LRT_P_K1_vs_K2[bic_table$K == 1L] <- 1.0
    bic_table$Boot_LRT_Stat_K1_vs_K2[bic_table$K == 2L] <- boot_stat
    bic_table$Boot_LRT_P_K1_vs_K2[bic_table$K == 2L] <- boot_p
  }

  if (best_k >= 2L) {
    delta_icl_best <- as.numeric(bic_table$Delta_ICL_vs_Null[match(best_k, bic_table$K)])
    collapse_to_1 <- if (isTRUE(is_intra_gmm)) {
      !(is.finite(delta_icl_best) && delta_icl_best > 0.0 && is.finite(boot_p) && boot_p < 0.05)
    } else {
      is.finite(delta_icl_best) && delta_icl_best <= 0.0 && is.finite(boot_p) && boot_p >= 0.05
    }
    if (isTRUE(collapse_to_1)) {
      orig_k <- as.integer(best_k)
      message(sprintf("[Model Selection Gate] Collapsed GMM K=%d -> K=1 (Delta_ICL=%.2f, Boot_LRT_p=%.4f).", orig_k, delta_icl_best, as.numeric(boot_p)))
      best_k <- 1L; cls_sub <- rep(1L, nrow(X_pcs)); unc_sub <- rep(0.0, nrow(X_pcs))
      best_name <- paste0(best_name, "_ICLBootCollapsed"); did_icl_collapse <- TRUE
    }
  }

  list(bic_table = bic_table, best_k = best_k, best_name = best_name,
       cls_sub = cls_sub, unc_sub = unc_sub, did_icl_collapse = did_icl_collapse)
}
