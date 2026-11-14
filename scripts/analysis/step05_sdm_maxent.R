# ==============================================================================
# Module: scripts/analysis/step05_sdm_maxent.R
# Step 05 Spatial RF Niche Figure Delegation + MaxEnt/Ensemble SDM Engine
# ==============================================================================

.ensure_step05_fig_modules <- function() {
  if (!exists("render_step05a_ecological_niche_png", mode = "function")) {
    for (.f in c("scripts/figures/theme_packera_pub.R", "scripts/figures/step05_figures.R")) {
      if (file.exists(.f)) source(.f)
    }
  }
}

export_dynamic_fig05_ecological_niche_phenology <- function(df, srf_res = NULL, niche_res = NULL, out_plot_path = NULL, pheno_res = NULL) {
  .ensure_step05_fig_modules()
  out_png <- if (!is.null(out_plot_path) && nzchar(out_plot_path)) normalize_png_path(out_plot_path) else "outputs/figures/fig05a_ecological_niche_and_rf_importance.png"
  render_step05a_ecological_niche_png(df = df, srf_res = srf_res, niche_res = niche_res, out_plot_path = out_png)
}

export_spatial_rf_figures <- function(df, srf_res, niche_res, out_pdf, pheno_res = NULL) {
  .ensure_step05_fig_modules()
  out_png <- if (!is.null(out_pdf) && nzchar(out_pdf)) normalize_png_path(out_pdf) else "outputs/figures/fig05a_ecological_niche_and_rf_importance.png"
  message(sprintf("Generating 300-DPI Spatial RF & Ecological Niche plate -> %s", out_png))
  render_step05a_ecological_niche_png(df = df, srf_res = srf_res, niche_res = niche_res, out_plot_path = out_png)
}

# ------------------------------------------------------------------------------
# 6b. Default MaxEnt (maxnet / Authority-Weighted IPP-Cloglog GLM) & Spatial RF Ensemble SDM
#     Projection, Continuous Boyce Index (Hirzel et al. 2006), & Habitat Maps
# ------------------------------------------------------------------------------
compute_continuous_boyce_index <- function(suitability_all, suitability_pres, n_bins = 20L) {
  s_all <- suppressWarnings(as.numeric(suitability_all))
  s_all <- s_all[is.finite(s_all)]
  s_pres <- suppressWarnings(as.numeric(suitability_pres))
  s_pres <- s_pres[is.finite(s_pres)]
  empty_curve <- data.frame(suitability = numeric(0), pe_ratio = numeric(0), stringsAsFactors = FALSE)
  if (length(s_pres) < 2L || length(s_all) < 4L) {
    res <- NA_real_
    attr(res, "boyce_curve") <- empty_curve
    return(res)
  }
  s_min <- min(c(s_all, s_pres), na.rm = TRUE)
  s_max <- max(c(s_all, s_pres), na.rm = TRUE)
  if (!is.finite(s_min) || !is.finite(s_max) || (s_max - s_min) < 1e-8) {
    res <- 0.0
    attr(res, "boyce_curve") <- empty_curve
    return(res)
  }
  n_b <- max(5L, as.integer(n_bins))
  win_w <- max((s_max - s_min) / 10.0, 1e-4)
  centers <- seq(s_min + 0.5 * win_w, s_max - 0.5 * win_w, length.out = n_b)
  pe_ratio <- rep(NA_real_, n_b)
  n_p <- length(s_pres)
  n_a <- length(s_all)
  for (i in seq_len(n_b)) {
    lo <- centers[i] - 0.5 * win_w
    hi <- centers[i] + 0.5 * win_w
    p_i <- sum(s_pres >= lo & s_pres <= hi) / n_p
    e_i <- sum(s_all >= lo & s_all <= hi) / n_a
    if (e_i > 0) pe_ratio[i] <- p_i / e_i
  }
  ok <- is.finite(pe_ratio) & is.finite(centers)
  curve_df <- data.frame(suitability = centers[ok], pe_ratio = pe_ratio[ok], stringsAsFactors = FALSE)
  if (sum(ok) < 3L || stats::sd(pe_ratio[ok]) < 1e-10) {
    res <- 0.0
    attr(res, "boyce_curve") <- curve_df
    return(res)
  }
  cbi <- suppressWarnings(stats::cor(centers[ok], pe_ratio[ok], method = "spearman", use = "complete.obs"))
  if (!is.finite(cbi)) cbi <- 0.0
  res <- round(as.numeric(cbi), 4)
  attr(res, "boyce_curve") <- curve_df
  res
}

.compute_roc_auc_curve <- function(y01, probs, wts = NULL) {
  y_v <- suppressWarnings(as.numeric(y01))
  p_v <- suppressWarnings(as.numeric(probs))
  w_v <- if (is.null(wts) || length(wts) != length(y_v)) rep(1.0, length(y_v)) else suppressWarnings(as.numeric(wts))
  w_v[!is.finite(w_v) | w_v <= 0] <- 1.0
  ok <- is.finite(y_v) & is.finite(p_v) & (y_v %in% c(0, 1))
  if (sum(y_v[ok] == 1) < 1L || sum(y_v[ok] == 0) < 1L) {
    return(list(auc = NA_real_, fpr = c(0, 1), tpr = c(0, 1)))
  }
  y_v <- y_v[ok]; p_v <- p_v[ok]; w_v <- w_v[ok]
  ord <- order(p_v, decreasing = TRUE)
  y_s <- y_v[ord]; p_s <- p_v[ord]; w_s <- w_v[ord]
  pos_w <- sum(w_s[y_s == 1])
  neg_w <- sum(w_s[y_s == 0])
  if (pos_w <= 0 || neg_w <= 0) return(list(auc = 0.5, fpr = c(0, 1), tpr = c(0, 1)))
  tp_cum <- c(0, cumsum(ifelse(y_s == 1, w_s, 0)) / pos_w)
  fp_cum <- c(0, cumsum(ifelse(y_s == 0, w_s, 0)) / neg_w)
  dups <- !duplicated(c(Inf, p_s), fromLast = TRUE)
  tpr <- tp_cum[dups]
  fpr <- fp_cum[dups]
  if ( utils::tail(fpr, 1) < 1.0 || utils::tail(tpr, 1) < 1.0 ) {
    fpr <- c(fpr, 1.0); tpr <- c(tpr, 1.0)
  }
  auc_val <- sum(0.5 * (tpr[-1] + tpr[-length(tpr)]) * diff(fpr))
  list(auc = round(min(1.0, max(0.0, auc_val)), 4), fpr = fpr, tpr = tpr)
}

.fit_single_maxent_model <- function(p_vec, X_df, w_vec = NULL) {
  if (requireNamespace("glmnet", quietly = TRUE)) {
    tryCatch(glmnet::glmnet.control(factory = TRUE), error = function(e) NULL)
    on.exit(tryCatch(glmnet::glmnet.control(factory = TRUE), error = function(e) NULL), add = TRUE)
  }

  p_num <- as.integer(as.numeric(p_vec) > 0.5)
  X_mat <- as.data.frame(X_df)
  pred_cols <- colnames(X_mat)
  n_rows <- nrow(X_mat)
  for (cn in pred_cols) {
    v <- suppressWarnings(as.numeric(X_mat[[cn]]))
    med_v <- stats::median(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(med_v)) med_v <- 0.0
    v[!is.finite(v)] <- med_v
    X_mat[[cn]] <- v
  }
  w_use <- if (is.null(w_vec) || length(w_vec) != n_rows) rep(1.0, n_rows) else as.numeric(w_vec)
  w_use[!is.finite(w_use) | w_use <= 0] <- 1.0
  n_pres <- sum(p_num == 1L)
  n_bg_raw <- sum(p_num == 0L)

  has_nonuniform_weights <- (length(w_use[p_num == 1L]) >= 2L) && (stats::sd(w_use[p_num == 1L], na.rm = TRUE) > 1e-4)

  # Inhomogeneous Poisson Process (IPP) quadrature augmentation (Phillips et al. 2017,
  # addsamplestobackground = TRUE): include unique presence environments in the background
  # quadrature set with weight 100 so density-ratio estimation is strictly concave and
  # never diverges on small or separable target-group background sets.
  pres_idx <- which(p_num == 1L)
  bg_idx <- which(p_num == 0L)
  X_pres <- X_mat[pres_idx, , drop = FALSE]
  X_bg <- X_mat[bg_idx, , drop = FALSE]
  w_pres_raw <- w_use[pres_idx]
  w_pres_norm <- w_pres_raw / max(mean(w_pres_raw, na.rm = TRUE), 1e-6)

  if (nrow(X_pres) > 0L) {
    bg_keys <- if (nrow(X_bg) > 0L) do.call(paste, c(lapply(X_bg, function(z) round(as.numeric(z), 6)), sep = "\r")) else character(0)
    pr_keys <- do.call(paste, c(lapply(X_pres, function(z) round(as.numeric(z), 6)), sep = "\r"))
    miss_in_bg <- !(pr_keys %in% bg_keys) & !duplicated(pr_keys)
    X_bg_aug <- if (any(miss_in_bg)) rbind(X_bg, X_pres[miss_in_bg, , drop = FALSE]) else X_bg
  } else {
    X_bg_aug <- X_bg
  }
  X_ipp <- rbind(X_pres, X_bg_aug)
  p_ipp <- c(rep(1L, nrow(X_pres)), rep(0L, nrow(X_bg_aug)))
  w_ipp <- c(w_pres_norm, rep(100.0, nrow(X_bg_aug)))
  n_bg_ipp <- nrow(X_bg_aug)

  # Unified [0, 1]-Rescaled, Collinearity-Pruned Inhomogeneous Poisson Process (IPP) MaxEnt Engine
  # (Phillips & Dudik 2008; Renner & Warton 2013; Phillips et al. 2017):
  #   - Avoids calling maxnet::maxnet() directly because maxnet mutates global glmnet.control(fdev = 0)
  #     and passes unscaled quadratic (~10^6) + 98-knot collinear hinge columns with standardize = FALSE,
  #     which causes coordinate-descent non-convergence on collinear raster surfaces.
  #   - Rescales all linear, quadratic, and 10-knot piecewise-linear hinge features strictly onto [0, 1],
  #     prunes redundant collinear basis functions (|r| > 0.985), and solves the exact Phillips & Dudik (2008)
  #     class-penalized IPP objective via standardized single-pass glmnet in milliseconds.
  vmin_ipp <- vapply(pred_cols, function(cn) min(as.numeric(X_ipp[[cn]]), na.rm = TRUE), numeric(1))
  vmax_ipp <- vapply(pred_cols, function(cn) max(as.numeric(X_ipp[[cn]]), na.rm = TRUE), numeric(1))
  vrng_ipp <- pmax(vmax_ipp - vmin_ipp, 1e-6)
  med_ipp  <- vapply(pred_cols, function(cn) stats::median(as.numeric(X_mat[[cn]]), na.rm = TRUE), numeric(1))
  hinge_knots <- seq(1.0 / 11.0, 10.0 / 11.0, length.out = 10L)

  .build_ipp_maxent_basis <- function(df_in, cl_mode = "lq") {
    nr <- nrow(df_in)
    np <- length(pred_cols)
    U <- matrix(0.0, nrow = nr, ncol = np, dimnames = list(NULL, paste0("l__", pred_cols)))
    for (j in seq_len(np)) {
      cn <- pred_cols[j]
      v <- if (cn %in% names(df_in)) suppressWarnings(as.numeric(df_in[[cn]])) else rep(med_ipp[j], nr)
      v[!is.finite(v)] <- med_ipp[j]
      U[, j] <- pmin(1.0, pmax(0.0, (v - vmin_ipp[j]) / vrng_ipp[j]))
    }
    mats <- list(U)
    ftypes <- rep("l", np)
    if (grepl("q", cl_mode, fixed = TRUE)) {
      Uq <- U * U
      colnames(Uq) <- paste0("q__", pred_cols)
      mats[[length(mats) + 1L]] <- Uq
      ftypes <- c(ftypes, rep("q", np))
    }
    if (grepl("h", cl_mode, fixed = TRUE)) {
      nk <- length(hinge_knots)
      H_mat <- matrix(0.0, nrow = nr, ncol = 2L * nk * np)
      h_nms <- character(2L * nk * np)
      col_ptr <- 1L
      for (j in seq_len(np)) {
        uj <- U[, j]
        cn <- pred_cols[j]
        for (ki in seq_len(nk)) {
          tk <- hinge_knots[ki]
          H_mat[, col_ptr] <- pmax(0.0, (uj - tk) / (1.0 - tk))
          h_nms[col_ptr] <- sprintf("hf_%s_%02d", cn, ki)
          col_ptr <- col_ptr + 1L
          H_mat[, col_ptr] <- pmax(0.0, (tk - uj) / tk)
          h_nms[col_ptr] <- sprintf("hr_%s_%02d", cn, ki)
          col_ptr <- col_ptr + 1L
        }
      }
      colnames(H_mat) <- h_nms
      mats[[length(mats) + 1L]] <- H_mat
      ftypes <- c(ftypes, rep("h", 2L * nk * np))
    }
    B_out <- do.call(cbind, mats)
    attr(B_out, "ftypes") <- ftypes
    B_out
  }

  cls_cascade <- if (n_pres >= 15L && n_bg_raw >= 60L) c("lqh", "lq", "l") else if (n_pres >= 8L) c("lq", "l") else c("l")

  if (requireNamespace("glmnet", quietly = TRUE) && n_pres >= 3L && n_bg_ipp >= 4L && length(pred_cols) >= 1L) {
    for (cl_opt in cls_cascade) {
      w_fit <- tryCatch({
        glmnet::glmnet.control(factory = TRUE)
        mm_raw <- .build_ipp_maxent_basis(X_ipp, cl_mode = cl_opt)
        ft_raw <- attr(mm_raw, "ftypes")
        col_sds <- apply(mm_raw, 2L, stats::sd)
        keep_v <- which(is.finite(col_sds) & col_sds > 1e-5)
        if (length(keep_v) < 1L) stop("zero-variance basis")
        mm <- mm_raw[, keep_v, drop = FALSE]
        ft <- ft_raw[keep_v]

        # Deduplicate near-collinear basis columns (|r| > 0.985), prioritizing strongest presence-background contrast
        if (ncol(mm) > 2L) {
          mu_pr <- colMeans(mm[p_ipp == 1L, , drop = FALSE])
          mu_bg <- colMeans(mm[p_ipp == 0L, , drop = FALSE])
          sep_score <- abs(mu_pr - mu_bg)
          ord_c <- order(sep_score, decreasing = TRUE)
          C_mat <- suppressWarnings(abs(stats::cor(mm[, ord_c, drop = FALSE])))
          C_mat[!is.finite(C_mat)] <- 0.0
          keep_ord <- logical(ncol(mm))
          for (ci in seq_len(ncol(mm))) {
            if (ci == 1L || all(C_mat[ci, which(keep_ord)] <= 0.985)) {
              keep_ord[ci] <- TRUE
            }
          }
          sel_idx <- sort(ord_c[keep_ord])
          mm <- mm[, sel_idx, drop = FALSE]
          ft <- ft[sel_idx]
        }

        # Phillips & Dudik (2008, Table 1 & Eq. 4) empirical feature-dispersion L1 regularization schedule
        c_lq <- stats::approx(x = c(0, 10, 15, 30, 100), y = c(1.00, 0.80, 0.50, 0.25, 0.05), xout = n_pres, rule = 2)$y
        c_mult <- ifelse(ft == "h", 0.50, c_lq)
        mm_pr <- mm[p_ipp == 1L, , drop = FALSE]
        wp_norm <- w_ipp[p_ipp == 1L] / max(sum(w_ipp[p_ipp == 1L]), 1e-8)
        mu_p <- as.numeric(crossprod(wp_norm, mm_pr))
        var_p <- as.numeric(crossprod(wp_norm, sweep(mm_pr, 2L, mu_p, "-")^2))
        reg <- c_mult * sqrt(pmax(var_p, 0.0) + (1.0 / max(n_pres, 1L))) / sqrt(max(n_pres, 1L))
        reg[!is.finite(reg) | reg <= 1e-6] <- 1e-3

        # Rescale penalty.factor by weighted column SD so standardize = TRUE optimizes exact IPP objective with O(1) conditioning
        w_tot <- sum(w_ipp)
        w_rel <- w_ipp / max(w_tot, 1e-8)
        mu_all <- as.numeric(crossprod(w_rel, mm))
        sd_all <- sqrt(pmax(as.numeric(crossprod(w_rel, sweep(mm, 2L, mu_all, "-")^2)), 1e-8))
        pf <- pmax(1e-4, reg / sd_all)
        pf <- pf / max(mean(pf), 1e-8)

        sum_p_eff <- sum(w_ipp[p_ipp == 1L])
        lam_seq <- (10^seq(3.2, 0.0, length.out = 50L)) * mean(reg) * (sum_p_eff / w_tot)
        pad_added <- FALSE
        if (ncol(mm) == 1L) {
          mm <- cbind(mm, `__zero_pad__` = 0.0)
          pf <- c(pf, 1.0)
          pad_added <- TRUE
        }
        gn <- suppressWarnings(glmnet::glmnet(
          x = mm,
          y = as.factor(p_ipp),
          family = "binomial",
          standardize = TRUE,
          penalty.factor = pf,
          lambda = lam_seq,
          weights = w_ipp,
          thresh = 1e-5,
          maxit = 1500L
        ))
        if (is.null(gn$beta) || ncol(gn$beta) < 1L) stop("glmnet returned empty path")
        bb <- as.numeric(gn$beta[, ncol(gn$beta)])
        names(bb) <- rownames(gn$beta)
        if (pad_added) bb <- bb[names(bb) != "__zero_pad__"]
        betas <- bb[is.finite(bb) & abs(bb) > 0.0]
        mm_bg <- mm[p_ipp == 0L, , drop = FALSE]
        lp_bg <- if (length(betas) > 0L) as.numeric(mm_bg[, names(betas), drop = FALSE] %*% betas) else rep(0.0, nrow(mm_bg))
        lp_shift <- max(lp_bg, na.rm = TRUE)
        if (!is.finite(lp_shift)) lp_shift <- 0.0
        rr <- exp(pmin(50.0, pmax(-50.0, lp_bg - lp_shift)))
        sum_rr <- max(sum(rr), 1e-12)
        raw <- rr / sum_rr
        ent <- -sum(raw * ifelse(raw > 1e-15, log(raw), 0.0))
        alp <- -log(sum_rr) - lp_shift
        list(cl_opt = cl_opt, betas = betas, entropy = ent, alpha = alp)
      }, error = function(e) NULL)

      if (!is.null(w_fit)) {
        cl_saved <- w_fit$cl_opt
        betas_saved <- w_fit$betas
        ent_saved <- w_fit$entropy
        alp_saved <- w_fit$alpha
        eng_lbl <- if (has_nonuniform_weights) sprintf("weighted_maxnet_glmnet (%s/cloglog)", cl_saved) else sprintf("maxnet_ipp_glmnet (%s/cloglog)", cl_saved)
        return(list(
          engine = eng_lbl,
          model = w_fit,
          predict_cloglog = function(new_df) {
            nd <- as.data.frame(new_df)
            if (nrow(nd) == 0L) return(numeric(0))
            if (length(betas_saved) == 0L) {
              pr0 <- rep(1.0 - exp(-exp(ent_saved + alp_saved)), nrow(nd))
              pr0[!is.finite(pr0)] <- 0.5
              return(pmin(1.0, pmax(0.0, pr0)))
            }
            B_new <- .build_ipp_maxent_basis(nd, cl_mode = cl_saved)
            use_b <- intersect(names(betas_saved), colnames(B_new))
            if (length(use_b) == 0L) return(rep(0.5, nrow(nd)))
            lp <- as.numeric(B_new[, use_b, drop = FALSE] %*% betas_saved[use_b]) + alp_saved
            pr <- 1.0 - exp(-exp(pmin(10.0, pmax(-20.0, ent_saved + lp))))
            pr[!is.finite(pr)] <- 0.0
            pmin(1.0, pmax(0.0, pr))
          }
        ))
      }
    }
  }

  # 2. Secondary Fallback Engine: Standardized IPP Cloglog Cross-Validated GLM (glmnet::cv.glmnet)
  Phi <- .build_ipp_maxent_basis(X_ipp, cl_mode = if (n_pres >= 8L) "lq" else "l")
  if (requireNamespace("glmnet", quietly = TRUE) && ncol(Phi) >= 2L && n_pres >= 4L && n_bg_ipp >= 4L) {
    glmnet::glmnet.control(factory = TRUE)
    sum_p_eff <- sum(w_ipp[p_ipp == 1L])
    lam_fb <- (10^seq(3.0, 0.0, length.out = 40L)) * (sum_p_eff / sum(w_ipp))
    nf_fb <- min(5L, max(3L, as.integer(floor(min(n_pres, n_bg_ipp) / 2L))))
    fid_ipp <- integer(length(p_ipp))
    fid_ipp[p_ipp == 1L] <- rep_len(seq_len(nf_fb), sum(p_ipp == 1L))
    fid_ipp[p_ipp == 0L] <- rep_len(seq_len(nf_fb), sum(p_ipp == 0L))
    cv_fit <- tryCatch(
      suppressWarnings(glmnet::cv.glmnet(
        x = Phi, y = as.factor(p_ipp), family = "binomial", standardize = TRUE,
        weights = w_ipp, alpha = 0.85, lambda = lam_fb, foldid = fid_ipp, thresh = 1e-5, maxit = 1000L
      )),
      error = function(e) NULL
    )
    if (!is.null(cv_fit)) {
      s_use <- if (is.finite(cv_fit$lambda.1se)) cv_fit$lambda.1se else cv_fit$lambda.min
      eta_tr <- as.numeric(stats::predict(cv_fit, newx = Phi, s = s_use, type = "link"))
      med_eta_p <- stats::median(eta_tr[p_ipp == 1L], na.rm = TRUE)
      if (!is.finite(med_eta_p)) med_eta_p <- 0.0
      c_ent <- -med_eta_p
      cl_fb <- if (n_pres >= 8L) "lq" else "l"
      return(list(
        engine = sprintf("glmnet_maxent (%s/cloglog)", cl_fb),
        model = cv_fit,
        predict_cloglog = function(new_df) {
          B_new <- .build_ipp_maxent_basis(as.data.frame(new_df), cl_mode = cl_fb)
          eta <- tryCatch(as.numeric(stats::predict(cv_fit, newx = B_new, s = s_use, type = "link")) + c_ent, error = function(e) rep(0.0, nrow(B_new)))
          eta <- pmin(10.0, pmax(-15.0, eta))
          pr <- 1.0 - exp(-exp(eta))
          pr[!is.finite(pr)] <- 0.0
          pmin(1.0, pmax(0.0, pr))
        }
      ))
    }
  }

  # 3. Ridge-stabilized weighted cloglog GLM fallback (zero external package dependency)
  X_des <- cbind(Intercept = 1.0, Phi)
  w_glm <- c(w_pres_norm, rep(2.0 * max(sum(w_pres_norm), 1e-6) / max(n_bg_ipp, 1L), n_bg_ipp))
  W_diag <- w_glm / max(mean(w_glm), 1e-6)
  P_pen <- diag(c(0.0, rep(0.25, ncol(Phi))), nrow = ncol(X_des), ncol = ncol(X_des))
  y_target <- ifelse(p_ipp == 1L, 0.85, -1.25)
  beta_hat <- tryCatch(
    solve(crossprod(X_des * sqrt(W_diag)) + P_pen, crossprod(X_des, W_diag * y_target)),
    error = function(e) rep(0.0, ncol(X_des))
  )
  beta_hat <- as.numeric(beta_hat)
  eta_tr <- as.numeric(X_des %*% beta_hat)
  med_eta_p <- stats::median(eta_tr[p_ipp == 1L], na.rm = TRUE)
  if (!is.finite(med_eta_p)) med_eta_p <- 0.0
  c_ent <- -med_eta_p
  cl_fb2 <- if (n_pres >= 8L) "lq" else "l"

  list(
    engine = sprintf("weighted_glm_maxent (%s/cloglog)", cl_fb2),
    model = list(beta = beta_hat, c_ent = c_ent),
    predict_cloglog = function(new_df) {
      B_new <- cbind(Intercept = 1.0, .build_ipp_maxent_basis(as.data.frame(new_df), cl_mode = cl_fb2))
      eta <- pmin(10.0, pmax(-15.0, as.numeric(B_new %*% beta_hat) + c_ent))
      pr <- 1.0 - exp(-exp(eta))
      pr[!is.finite(pr)] <- 0.0
      pmin(1.0, pmax(0.0, pr))
    }
  )
}

.fit_single_srf_binary_model <- function(p_vec, X_df, w_vec = NULL, n_trees = 250L) {
  p_num <- as.integer(as.numeric(p_vec) > 0.5)
  X_mat <- as.data.frame(X_df)
  pred_cols <- colnames(X_mat)
  n_rows <- nrow(X_mat)
  w_use <- if (is.null(w_vec) || length(w_vec) != n_rows) rep(1.0, n_rows) else as.numeric(w_vec)
  w_use[!is.finite(w_use) | w_use <= 0] <- 1.0
  # Inverse class-frequency weighting (Valavi et al. 2021)
  n_p <- max(sum(p_num == 1L), 1L)
  n_b <- max(sum(p_num == 0L), 1L)
  w_rf <- w_use * ifelse(p_num == 1L, n_rows / (2.0 * n_p), n_rows / (2.0 * n_b))
  y_fac <- factor(ifelse(p_num == 1L, "Presence", "Background"), levels = c("Background", "Presence"))

  if (requireNamespace("ranger", quietly = TRUE) && sum(p_num == 1L) >= 2L && sum(p_num == 0L) >= 2L) {
    max_w <- if (n_rows <= 2500L) 2L else 16L
    rf_threads <- if (exists("resolve_packera_workers", mode = "function")) min(max_w, resolve_packera_workers(max_w)) else 1L
    tr_df <- data.frame(y = y_fac, X_mat, check.names = FALSE)
    fit_rg <- tryCatch(
      ranger::ranger(y ~ ., data = tr_df, num.trees = n_trees, probability = TRUE, case.weights = w_rf, num.threads = rf_threads, seed = 42),
      error = function(e) NULL
    )
    if (!is.null(fit_rg)) {
      return(list(
        model = fit_rg,
        predict_prob = function(new_df) {
          nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]
          if (nrow(nd) == 0L) return(numeric(0))
          for (cn in pred_cols) {
            v <- suppressWarnings(as.numeric(nd[[cn]]))
            v[!is.finite(v)] <- stats::median(X_mat[[cn]], na.rm = TRUE)
            nd[[cn]] <- v
          }
          pr_m <- tryCatch(stats::predict(fit_rg, data = nd, num.threads = rf_threads)$predictions, error = function(e) NULL)
          if (is.matrix(pr_m) && "Presence" %in% colnames(pr_m)) {
            pr <- as.numeric(pr_m[, "Presence"])
          } else {
            pr <- rep(0.5, nrow(nd))
          }
          pr[!is.finite(pr)] <- 0.0
          pmin(1.0, pmax(0.0, pr))
        }
      ))
    }
  }
  if (requireNamespace("randomForest", quietly = TRUE) && sum(p_num == 1L) >= 2L && sum(p_num == 0L) >= 2L) {
    fit_rf <- tryCatch(randomForest::randomForest(x = X_mat, y = y_fac, ntree = n_trees), error = function(e) NULL)
    if (!is.null(fit_rf)) {
      return(list(
        model = fit_rf,
        predict_prob = function(new_df) {
          nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]
          pr_m <- tryCatch(stats::predict(fit_rf, newdata = nd, type = "prob"), error = function(e) NULL)
          pr <- if (is.matrix(pr_m) && "Presence" %in% colnames(pr_m)) as.numeric(pr_m[, "Presence"]) else rep(0.5, nrow(nd))
          pr[!is.finite(pr)] <- 0.0
          pmin(1.0, pmax(0.0, pr))
        }
      ))
    }
  }
  # Fallback to MaxEnt cloglog if neither ranger nor randomForest is available
  mx_fb <- .fit_single_maxent_model(p_vec, X_df, w_vec)
  list(model = mx_fb$model, predict_prob = mx_fb$predict_cloglog)
}

.build_regional_sdm_raster_stack <- function(df, env_dir, kept_vars, ncols = 110L, nrows = 90L, demo = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE)) return(NULL)
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  lat_v <- suppressWarnings(as.numeric(df$latitude))
  ok_ll <- is.finite(lon_v) & is.finite(lat_v)
  if (!any(ok_ll)) return(NULL)
  if (sum(ok_ll) >= 20L) {
    lon_rng <- as.numeric(stats::quantile(lon_v[ok_ll], probs = c(0.005, 0.995), na.rm = TRUE))
    lat_rng <- as.numeric(stats::quantile(lat_v[ok_ll], probs = c(0.005, 0.995), na.rm = TRUE))
  } else {
    lon_rng <- range(lon_v[ok_ll], na.rm = TRUE)
    lat_rng <- range(lat_v[ok_ll], na.rm = TRUE)
  }
  pad_x <- max(0.75, 0.10 * diff(lon_rng))
  pad_y <- max(0.75, 0.10 * diff(lat_rng))
  med_lat <- stats::median(lat_v[ok_ll], na.rm = TRUE)
  min_lat_clamp <- if (is.finite(med_lat) && med_lat >= 28.0) 24.5 else 14.0
  xmin <- max(-130.0, lon_rng[1] - pad_x)
  xmax <- min(-50.0,  lon_rng[2] + pad_x)
  ymin <- max(min_lat_clamp, lat_rng[1] - pad_y)
  ymax <- min(66.0,   lat_rng[2] + pad_y)
  if (xmax - xmin < 1.0) { xmin <- xmin - 0.75; xmax <- xmax + 0.75 }
  if (ymax - ymin < 1.0) { ymin <- ymin - 0.75; ymax <- ymax + 0.75 }

  study_ext <- terra::ext(xmin, xmax, ymin, ymax)
  tmpl_rast <- terra::rast(ext = study_ext, ncols = ncols, nrows = nrows, crs = "EPSG:4326")

  .crop_resample_to_tmpl <- function(r_src, nm_out, scale_div = 1.0) {
    tryCatch({
      r_p <- r_src
      if (!isTRUE(tryCatch(terra::same.crs(r_p, tmpl_rast), error = function(e) FALSE))) {
        ext_src <- terra::project(terra::as.polygons(study_ext, crs = "EPSG:4326"), terra::crs(r_p))
        r_c <- terra::crop(r_p, ext_src)
        rs <- terra::project(r_c, tmpl_rast, method = "bilinear")
      } else {
        rc <- terra::crop(r_p, study_ext)
        rs <- terra::resample(rc, tmpl_rast, method = "bilinear")
      }
      if (scale_div != 1.0) rs <- rs / scale_div
      names(rs) <- nm_out
      rs
    }, error = function(e) NULL)
  }

  loaded_layers <- list()
  raw_aux_rasters <- list()
  if (!is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir)) {
    tif_files <- discover_env_rasters(env_dir, recursive = TRUE, crop_global = TRUE)
    if (length(tif_files) > 0L) {
      for (tf in tif_files) {
        bn <- tolower(basename(tf))
        r_t <- tryCatch(terra::rast(tf), error = function(e) NULL)
        if (is.null(r_t)) next
        nms_t <- names(r_t)
        for (li in seq_along(nms_t)) {
          lnm <- nms_t[li]
          mapped_nm <- lnm
          scale_div <- 1.0
          if (!(mapped_nm %in% kept_vars)) {
            m_bio <- regmatches(lnm, regexec("bio[_.-]?([0-9]{1,2})$", tolower(lnm)))[[1]]
            if (length(m_bio) == 2L && m_bio[2] %in% names(WC_BIO_LOOKUP)) {
              mapped_nm <- WC_BIO_LOOKUP[[m_bio[2]]]
            } else if (grepl("^ph1to1h2o_5_cm", bn)) {
              mapped_nm <- "fine_soil_ph"; scale_div <- 100.0
            } else if (grepl("phh2o.*5-15|phh2o.*5_15", paste(bn, tolower(lnm)))) {
              mapped_nm <- if ("fine_soil_ph" %in% kept_vars && !("fine_soil_ph" %in% names(loaded_layers))) "fine_soil_ph" else "soil_ph"
              v_med <- stats::median(terra::values(r_t[[li]]), na.rm = TRUE)
              if (is.finite(v_med) && v_med > 14.5) scale_div <- 10.0
            } else if (grepl("^cec7_5_cm", bn)) {
              mapped_nm <- "fine_soil_cec"; scale_div <- 10.0
            } else if (grepl("cec.*5-15|cec.*5_15", paste(bn, tolower(lnm)))) {
              mapped_nm <- if ("fine_soil_cec" %in% kept_vars && !("fine_soil_cec" %in% names(loaded_layers))) "fine_soil_cec" else "soil_cec"
              v_med <- stats::median(terra::values(r_t[[li]]), na.rm = TRUE)
              if (is.finite(v_med) && v_med > 90.0) scale_div <- 10.0
            } else if (grepl("^fragvol_5_cm", bn)) {
              mapped_nm <- "fine_cfvo_pct"
              if (is.null(raw_aux_rasters[["cfvo"]])) raw_aux_rasters[["cfvo"]] <- .crop_resample_to_tmpl(r_t[[li]], "cfvo", 1.0)
            } else if (grepl("cfvo.*5-15|cfvo.*5_15", paste(bn, tolower(lnm)))) {
              mapped_nm <- if ("fine_cfvo_pct" %in% kept_vars && !("fine_cfvo_pct" %in% names(loaded_layers))) "fine_cfvo_pct" else "soil_cfvo"
              v_med <- stats::median(terra::values(r_t[[li]]), na.rm = TRUE)
              if (is.finite(v_med) && v_med > 105.0) scale_div <- 10.0
              if (is.null(raw_aux_rasters[["cfvo"]])) raw_aux_rasters[["cfvo"]] <- .crop_resample_to_tmpl(r_t[[li]], "cfvo", scale_div)
            } else if (grepl("^anylithicdpt|^resdept", bn)) {
              mapped_nm <- "fine_bedrock_depth_cm"
              if (is.null(raw_aux_rasters[["bedrock"]])) raw_aux_rasters[["bedrock"]] <- .crop_resample_to_tmpl(r_t[[li]], "bedrock", 1.0)
            } else if (grepl("^sandtotal_5_cm|sand.*5-15|sand.*5_15", paste(bn, tolower(lnm)))) {
              raw_aux_rasters[["sand"]] <- .crop_resample_to_tmpl(r_t[[li]], "sand", if (grepl("sandtotal", bn)) 1.0 else { vm <- stats::median(terra::values(r_t[[li]]), na.rm = TRUE); if (is.finite(vm) && vm > 105) 10.0 else 1.0 })
            } else if (grepl("^claytotal_5_cm|clay.*5-15|clay.*5_15", paste(bn, tolower(lnm)))) {
              raw_aux_rasters[["clay"]] <- .crop_resample_to_tmpl(r_t[[li]], "clay", if (grepl("claytotal", bn)) 1.0 else { vm <- stats::median(terra::values(r_t[[li]]), na.rm = TRUE); if (is.finite(vm) && vm > 105) 10.0 else 1.0 })
            } else if (grepl("^silttotal_5_cm|silt.*5-15|silt.*5_15", paste(bn, tolower(lnm)))) {
              raw_aux_rasters[["silt"]] <- .crop_resample_to_tmpl(r_t[[li]], "silt", if (grepl("silttotal", bn)) 1.0 else { vm <- stats::median(terra::values(r_t[[li]]), na.rm = TRUE); if (is.finite(vm) && vm > 105) 10.0 else 1.0 })
            }
          }
          if (mapped_nm %in% kept_vars && !(mapped_nm %in% names(loaded_layers))) {
            r_sub <- .crop_resample_to_tmpl(r_t[[li]], mapped_nm, scale_div = scale_div)
            if (!is.null(r_sub)) loaded_layers[[mapped_nm]] <- r_sub
          }
        }
      }
    }
  }

  # Derive continuous ilr_texture_1 / ilr_texture_2 rasters when sand + clay rasters exist
  if ((("ilr_texture_1" %in% kept_vars && !("ilr_texture_1" %in% names(loaded_layers))) ||
       ("ilr_texture_2" %in% kept_vars && !("ilr_texture_2" %in% names(loaded_layers)))) &&
      !is.null(raw_aux_rasters[["sand"]]) && !is.null(raw_aux_rasters[["clay"]])) {
    sv <- terra::values(raw_aux_rasters[["sand"]])[, 1]
    cv <- terra::values(raw_aux_rasters[["clay"]])[, 1]
    siv <- if (!is.null(raw_aux_rasters[["silt"]])) terra::values(raw_aux_rasters[["silt"]])[, 1] else NULL
    tex_g <- compute_soil_texture_log_ratios(sv, cv, silt_pct = siv)
    if ("ilr_texture_1" %in% kept_vars && !("ilr_texture_1" %in% names(loaded_layers))) {
      r_i1 <- terra::rast(tmpl_rast); terra::values(r_i1) <- tex_g$ilr_texture_1; names(r_i1) <- "ilr_texture_1"
      loaded_layers[["ilr_texture_1"]] <- r_i1
    }
    if ("ilr_texture_2" %in% kept_vars && !("ilr_texture_2" %in% names(loaded_layers))) {
      r_i2 <- terra::rast(tmpl_rast); terra::values(r_i2) <- tex_g$ilr_texture_2; names(r_i2) <- "ilr_texture_2"
      loaded_layers[["ilr_texture_2"]] <- r_i2
    }
  }

  # Derive continuous lithic_outcrop_index raster from cached bedrock (+ optional cfvo) rasters
  if ("lithic_outcrop_index" %in% kept_vars && !("lithic_outcrop_index" %in% names(loaded_layers)) && !is.null(raw_aux_rasters[["bedrock"]])) {
    bd_v <- pmin(pmax(as.numeric(terra::values(raw_aux_rasters[["bedrock"]])[, 1]), 0.0), 150.0)
    cfvo_mult <- if (!is.null(raw_aux_rasters[["cfvo"]])) {
      cv_v <- pmin(pmax(as.numeric(terra::values(raw_aux_rasters[["cfvo"]])[, 1]), 0.0), 100.0)
      ifelse(is.finite(cv_v), 0.7 + 0.3 * (cv_v / 100.0), 1.0)
    } else rep(1.0, length(bd_v))
    loi_v <- 100.0 * pmax(0.0, 1.0 - bd_v / 150.0) * cfvo_mult
    r_loi <- terra::rast(tmpl_rast); terra::values(r_loi) <- loi_v; names(r_loi) <- "lithic_outcrop_index"
    loaded_layers[["lithic_outcrop_index"]] <- r_loi
  }

  # Enforce domain-coverage integrity when demo == FALSE: never fabricate continental rasters via point IDW!
  if (!isTRUE(demo)) {
    if (length(loaded_layers) == 0L) {
      warning("Skipping Step 05 continuous SDM raster map projection: zero real GeoTIFF layers matched kept_vars in env_dir (IDW raster fabrication is disabled when demo == FALSE).")
      return(NULL)
    }
    bio_nms <- names(loaded_layers)[grepl("^bio[0-9]+_", names(loaded_layers))]
    if (length(bio_nms) > 0L) {
      land_ref <- loaded_layers[[bio_nms[1L]]]
      for (nm_l in setdiff(names(loaded_layers), bio_nms)) {
        rl_f <- tryCatch({
          rf <- terra::focal(loaded_layers[[nm_l]], w = 5, fun = mean, na.rm = TRUE, na.policy = "only")
          names(rf) <- nm_l
          terra::mask(rf, land_ref)
        }, error = function(e) NULL)
        if (!is.null(rl_f)) loaded_layers[[nm_l]] <- rl_f
      }
    }
    cov_counts <- vapply(loaded_layers, function(rl) sum(is.finite(terra::values(rl)[, 1])), integer(1))
    ref_land_cells <- max(cov_counts, 1L)
    valid_cov_nms <- names(cov_counts)[cov_counts >= floor(0.70 * ref_land_cells)]
    loaded_layers <- loaded_layers[intersect(kept_vars, valid_cov_nms)]
    if (length(loaded_layers) < 2L) {
      warning(sprintf("Skipping Step 05 continuous SDM raster map projection: only %d raster layer(s) cover >=70%% of study domain in '%s' (demo == FALSE).", length(loaded_layers), env_dir))
      return(NULL)
    }
    return(terra::rast(loaded_layers))
  }

  land_mask_rast <- if (length(loaded_layers) > 0L) loaded_layers[[1L]] else NULL

  # Demonstration mode ONLY (demo == TRUE): interpolate missing demo rasters across template_rast
  miss_vars <- setdiff(kept_vars, names(loaded_layers))
  if (length(miss_vars) > 0L) {
    u_key <- paste(round(lon_v[ok_ll], 4), round(lat_v[ok_ll], 4), sep = "_")
    u_idx <- which(ok_ll)[!duplicated(u_key)]
    if (length(u_idx) > 450L) {
      u_idx <- u_idx[unique(round(seq(1L, length(u_idx), length.out = 450L)))]
    }
    vx <- lon_v[u_idx]
    vy <- lat_v[u_idx]
    V_mat <- matrix(0.0, nrow = length(u_idx), ncol = length(miss_vars), dimnames = list(NULL, miss_vars))
    for (j in seq_along(miss_vars)) {
      vn <- miss_vars[j]
      vv <- if (vn %in% names(df)) suppressWarnings(as.numeric(df[[vn]][u_idx])) else rep(0.0, length(u_idx))
      med_v <- stats::median(vv[is.finite(vv)], na.rm = TRUE)
      if (!is.finite(med_v)) med_v <- 0.0
      vv[!is.finite(vv)] <- med_v
      V_mat[, j] <- vv
    }
    g_xy <- terra::xyFromCell(tmpl_rast, seq_len(terra::ncell(tmpl_rast)))
    D2 <- outer(g_xy[, 1], vx, "-")^2 + outer(g_xy[, 2], vy, "-")^2 + 0.0025
    W <- 1.0 / (D2^1.15)
    W <- W / pmax(rowSums(W), 1e-12)
    G_mat <- W %*% V_mat
    for (j in seq_along(miss_vars)) {
      vn <- miss_vars[j]
      r_idw <- terra::rast(tmpl_rast)
      terra::values(r_idw) <- as.numeric(G_mat[, j])
      if (!is.null(land_mask_rast)) {
        r_idw <- tryCatch(terra::mask(r_idw, land_mask_rast), error = function(e) r_idw)
      }
      names(r_idw) <- vn
      loaded_layers[[vn]] <- r_idw
    }
  }

  ord_layers <- loaded_layers[kept_vars]
  terra::rast(ord_layers)
}

fit_and_project_maxent_sdm_maps <- function(
  df,
  srf_res,
  env_dir,
  out_dir = "outputs",
  kept_vars = NULL,
  primary_clim_vars = NULL,
  primary_edaph_vars = NULL,
  focal_taxon = "Packera dubia",
  intraspecific_mode = FALSE,
  demo = FALSE
) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  if (is.null(out_dir) || !nzchar(out_dir)) out_dir <- "outputs"
  message("=== Fitting MaxEnt (cloglog) & Spatial RF Ensemble SDMs & Projecting Suitability Rasters ===")

  # Resolve non-MEM environmental predictors
  if (is.null(kept_vars) || length(kept_vars) == 0L) {
    if (!is.null(srf_res) && !is.null(srf_res$collinearity_screening$kept_vars)) {
      kept_vars <- srf_res$collinearity_screening$kept_vars
    } else {
      imp_csv <- file.path(out_dir, "reports", "srf_feature_importance.csv")
      if (file.exists(imp_csv)) {
        imp_tbl <- tryCatch(utils::read.csv(imp_csv, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(imp_tbl) && "Predictor" %in% names(imp_tbl)) {
          kept_vars <- as.character(imp_tbl$Predictor)
        }
      }
      if (is.null(kept_vars) || length(kept_vars) == 0L) {
        cand_eco <- discover_candidate_eco_predictors(df, include_pheno = FALSE, min_finite = 2L)
        if (length(cand_eco) > 12L && exists("screen_collinear_predictors_by_block", mode = "function")) {
          kept_vars <- tryCatch(
            screen_collinear_predictors_by_block(df[, cand_eco, drop = FALSE], cor.threshold = 0.70, vif.threshold = 5.0)$kept_vars,
            error = function(e) NULL
          )
        }
        if (is.null(kept_vars) || length(kept_vars) == 0L) {
          clim_c <- cand_eco[grepl("^bio[0-9]+", cand_eco)]
          soil_c <- setdiff(cand_eco, clim_c)
          kept_vars <- c(utils::head(clim_c, 6L), utils::head(soil_c, 6L))
        }
      }
    }
  }
  kept_vars <- intersect(kept_vars[!grepl("^MEM", kept_vars)], names(df))
  if (length(kept_vars) > 14L) {
    kv_clim <- kept_vars[grepl("^bio[0-9]+", kept_vars)]
    kv_edaph <- setdiff(kept_vars, kv_clim)
    kept_vars <- c(utils::head(kv_clim, 7L), utils::head(kv_edaph, 7L))
  }
  if (length(kept_vars) == 0L) {
    kept_vars <- intersect(c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_ph", "soil_ph", "fine_soil_sand", "soil_sand"), names(df))
  }
  if (length(kept_vars) == 0L) return(invisible(NULL))

  # Ensure numeric finite predictor values in df
  for (kv in kept_vars) {
    v <- suppressWarnings(as.numeric(df[[kv]]))
    med_v <- stats::median(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(med_v)) med_v <- 0.0
    v[!is.finite(v)] <- med_v
    df[[kv]] <- v
  }

  # Identify top climatic and edaphic predictors
  imp_df <- if (!is.null(srf_res) && is.data.frame(srf_res$importance)) srf_res$importance else NULL
  edaph_rx <- "^(fine_|ilr_texture_|ssurgo_|polaris_|solus_|topo_|soil_|sgmc_|is_azonal|lithic_outcrop|min_bedrock|ph_q|max_cfvo)"
  clim_cand <- if (!is.null(primary_clim_vars) && length(intersect(primary_clim_vars, kept_vars)) > 0L) {
    intersect(primary_clim_vars, kept_vars)
  } else {
    kept_vars[!grepl(edaph_rx, kept_vars)]
  }
  edaph_cand <- if (!is.null(primary_edaph_vars) && length(intersect(primary_edaph_vars, kept_vars)) > 0L) {
    intersect(primary_edaph_vars, kept_vars)
  } else {
    kept_vars[grepl(edaph_rx, kept_vars)]
  }
  if (!is.null(imp_df) && all(c("Predictor", "MeanDecreaseAccuracy") %in% names(imp_df))) {
    imp_map <- stats::setNames(as.numeric(imp_df$MeanDecreaseAccuracy), as.character(imp_df$Predictor))
    if (length(clim_cand) > 1L) clim_cand <- clim_cand[order(ifelse(clim_cand %in% names(imp_map), imp_map[clim_cand], -Inf), decreasing = TRUE)]
    if (length(edaph_cand) > 1L) edaph_cand <- edaph_cand[order(ifelse(edaph_cand %in% names(imp_map), imp_map[edaph_cand], -Inf), decreasing = TRUE)]
  }
  top_clim_pred <- if (length(clim_cand) > 0L) clim_cand[1L] else kept_vars[1L]
  top_edaph_pred <- if (length(edaph_cand) > 0L) edaph_cand[1L] else "None_Retained"

  # Resolve authority prior weights (preserving 0.0 on unpromoted Tier_3_Bronze & Cleanlab noise!)
  raw_auth_w <- if ("Combined_Prior_Weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$Combined_Prior_Weight))
  } else if ("authority_weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$authority_weight))
  } else {
    resolve_ecological_authority_weights(df)
  }
  raw_auth_w[!is.finite(raw_auth_w) | raw_auth_w < 0.0] <- 0.0
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    .identify_type_specimen_vec(df)
  } else {
    is_t <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, nrow(df))
    is_t[is.na(is_t)] <- FALSE
    list(is_type_specimen = is_t, type_rank = ifelse(is_t, 2L, 0L))
  }
  cl_sdm <- if ("is_confident_label_noise" %in% names(df)) as.logical(df$is_confident_label_noise) else rep(FALSE, nrow(df))
  cl_sdm[is.na(cl_sdm)] <- FALSE
  raw_auth_w[cl_sdm & !type_info$is_type_specimen] <- 0.0
  pos_m <- raw_auth_w > 0.0
  w_max_curr <- if (any(pos_m)) max(raw_auth_w[pos_m], 1e-6) else 1.0
  auth_w_01 <- ifelse(pos_m, pmin(1.0, pmax(0.20, raw_auth_w / w_max_curr)), 0.0)
  auth_w_01[type_info$is_type_specimen] <- 1.00
  df$sdm_prior_weight_01 <- auth_w_01

  # Determine training-eligible coordinates (coordinateUncertaintyInMeters <= 10000 & finite lat/lon)
  lon_all <- suppressWarnings(as.numeric(df$longitude))
  lat_all <- suppressWarnings(as.numeric(df$latitude))
  unc_col <- if ("coordinateUncertaintyInMeters" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
  } else if ("coordinateUncertainty" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertainty))
  } else {
    rep(NA_real_, nrow(df))
  }
  elig_coord_only <- is.finite(lon_all) & is.finite(lat_all) & (is.na(unc_col) | unc_col <= 10000.0 | type_info$is_type_specimen)
  if (sum(elig_coord_only) < 4L) elig_coord_only <- is.finite(lon_all) & is.finite(lat_all)
  elig_mask <- elig_coord_only & (df$sdm_prior_weight_01 > 0.0)
  if (sum(elig_mask) < 4L) {
    elig_mask <- elig_coord_only & (!cl_sdm | type_info$is_type_specimen)
    df$sdm_prior_weight_01[elig_mask & df$sdm_prior_weight_01 <= 0.0] <- 0.50
  }

  # Build regional environmental SpatRaster stack across the study bounding box (no IDW when demo == FALSE)
  env_stack <- .build_regional_sdm_raster_stack(df[elig_coord_only, , drop = FALSE], env_dir = env_dir, kept_vars = kept_vars, demo = demo)
  if (!is.null(env_stack) && length(intersect(kept_vars, names(env_stack))) >= 2L) {
    kept_vars <- intersect(kept_vars, names(env_stack))
    if (!(top_clim_pred %in% kept_vars)) top_clim_pred <- kept_vars[1L]
    edaph_in_stack <- intersect(edaph_cand, kept_vars)
    top_edaph_pred <- if (length(edaph_in_stack) > 0L) edaph_in_stack[1L] else if (length(edaph_cand) > 0L) edaph_cand[1L] else "None_Retained"
  }

  # Generate supplemental random background cells from env_stack so target-group + random background >= 1,000 points
  set.seed(42)
  bg_pool_df <- NULL
  if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE)) {
    n_cells_tot <- terra::ncell(env_stack)
    s_idx <- sample.int(n_cells_tot, size = min(1200L, n_cells_tot), replace = FALSE)
    s_xy <- terra::xyFromCell(env_stack, s_idx)
    s_vals <- as.data.frame(terra::values(env_stack)[s_idx, , drop = FALSE])
    colnames(s_vals) <- names(env_stack)
    bg_pool_df <- data.frame(
      longitude = as.numeric(s_xy[, 1]),
      latitude = as.numeric(s_xy[, 2]),
      s_vals[, kept_vars, drop = FALSE],
      sdm_prior_weight_01 = 1.0,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    ok_bg <- stats::complete.cases(bg_pool_df[, c("longitude", "latitude", kept_vars), drop = FALSE])
    bg_pool_df <- bg_pool_df[ok_bg, , drop = FALSE]
  }

  # Determine target taxa / morphotypes to model
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2L)
  grp_vec <- if (is_intra && "analysis_group" %in% names(df)) {
    as.character(df$analysis_group)
  } else {
    as.character(df$species_standardized)
  }
  grp_vec[is.na(grp_vec) | !nzchar(grp_vec)] <- "Unknown"
  tbl_elig <- table(grp_vec[elig_mask & grp_vec != "Unknown"])
  target_groups <- names(tbl_elig[tbl_elig >= 5L])
  if (length(target_groups) == 0L) target_groups <- names(tbl_elig[tbl_elig >= 2L])
  if (length(target_groups) == 0L && any(grp_vec != "Unknown")) target_groups <- unique(grp_vec[grp_vec != "Unknown"])

  # Also include focal_taxon explicitly if present in species_standardized with >= 3 records and not already in target_groups
  if (!is.null(focal_taxon) && nzchar(focal_taxon) && !(focal_taxon %in% target_groups)) {
    if (sum(df$species_standardized == focal_taxon & elig_mask, na.rm = TRUE) >= 3L) {
      target_groups <- c(focal_taxon, target_groups)
    }
  }
  # Order focal_taxon first
  if (!is.null(focal_taxon) && focal_taxon %in% target_groups) {
    target_groups <- c(focal_taxon, setdiff(target_groups, focal_taxon))
  }

  block_deg_use <- if (!is.null(srf_res) && is.finite(srf_res$spatial_autocorr_range_km)) {
    max(0.25, min(2.5, as.numeric(srf_res$spatial_autocorr_range_km) / 111.0))
  } else {
    0.50
  }

  raster_dir <- file.path(out_dir, "rasters")
  write_prod_rasters <- !(isTRUE(demo) && is_canonical_prod_path(raster_dir))
  if (write_prod_rasters) dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
  gdal_deflate <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

  perf_rows <- list()
  focal_artifacts <- NULL
  all_sdm_artifacts <- list()

  for (tg in target_groups) {
    pres_mask <- if (tg %in% df$species_standardized && !(tg %in% grp_vec)) {
      elig_mask & (df$species_standardized == tg)
    } else {
      elig_mask & (grp_vec == tg)
    }
    pres_df <- df[pres_mask, c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
    n_pres <- nrow(pres_df)
    if (n_pres < 2L) next

    # Target-group background points (Phillips et al. 2009): all other georeferenced Packera localities + supplemental raster cells up to >= 1,000
    tg_bg_mask <- elig_coord_only & !(if (tg %in% df$species_standardized && !(tg %in% grp_vec)) (df$species_standardized == tg) else (grp_vec == tg))
    tg_bg_df <- df[tg_bg_mask, c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
    if (nrow(tg_bg_df) > 0L) tg_bg_df$sdm_prior_weight_01 <- 1.0
    n_need_supp <- max(0L, 1000L - nrow(tg_bg_df))
    if (n_need_supp > 0L && !is.null(bg_pool_df) && nrow(bg_pool_df) > 0L) {
      take_n <- min(n_need_supp, nrow(bg_pool_df))
      supp_df <- bg_pool_df[seq_len(take_n), c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
      bg_df <- rbind(tg_bg_df, supp_df)
    } else if (nrow(tg_bg_df) > 0L) {
      bg_df <- tg_bg_df
    } else if (!is.null(bg_pool_df) && nrow(bg_pool_df) > 0L) {
      bg_df <- bg_pool_df
    } else {
      next
    }
    n_bg <- nrow(bg_df)

    comb_df <- rbind(pres_df, bg_df)
    p_vec <- c(rep(1L, n_pres), rep(0L, n_bg))
    w_vec <- as.numeric(comb_df$sdm_prior_weight_01)
    X_comb <- comb_df[, kept_vars, drop = FALSE]
    coords_comb <- as.matrix(comb_df[, c("longitude", "latitude")])

    # Spatial-Block Cross-Validation (k = 5 spatial blocks)
    k_folds <- min(5L, max(2L, n_pres))
    block_ids <- assign_spatial_grid_blocks(coords_comb, block_deg = block_deg_use, k = k_folds)
    # Ensure presence points span at least 2 folds when possible
    if (length(unique(block_ids[p_vec == 1L])) < 2L && n_pres >= 2L) {
      block_ids[p_vec == 1L] <- ((seq_len(n_pres) - 1L) %% min(k_folds, n_pres)) + 1L
    }

    oof_mx <- rep(NA_real_, length(p_vec))
    oof_rf <- rep(NA_real_, length(p_vec))
    u_folds <- sort(unique(block_ids))
    for (fid in u_folds) {
      tr_idx <- which(block_ids != fid)
      te_idx <- which(block_ids == fid)
      if (sum(p_vec[tr_idx] == 1L) >= 2L && sum(p_vec[tr_idx] == 0L) >= 2L && length(te_idx) >= 1L) {
        mx_cv <- .fit_single_maxent_model(p_vec[tr_idx], X_comb[tr_idx, , drop = FALSE], w_vec[tr_idx])
        rf_cv <- .fit_single_srf_binary_model(p_vec[tr_idx], X_comb[tr_idx, , drop = FALSE], w_vec[tr_idx], n_trees = 200L)
        oof_mx[te_idx] <- mx_cv$predict_cloglog(X_comb[te_idx, , drop = FALSE])
        oof_rf[te_idx] <- rf_cv$predict_prob(X_comb[te_idx, , drop = FALSE])
      }
    }

    # Fit full models on all presence + target-group background data
    mx_full <- .fit_single_maxent_model(p_vec, X_comb, w_vec)
    rf_full <- .fit_single_srf_binary_model(p_vec, X_comb, w_vec, n_trees = 300L)

    # Fill any fold that lacked training presences with full model predictions
    na_mx <- !is.finite(oof_mx)
    if (any(na_mx)) oof_mx[na_mx] <- mx_full$predict_cloglog(X_comb[na_mx, , drop = FALSE])
    na_rf <- !is.finite(oof_rf)
    if (any(na_rf)) oof_rf[na_rf] <- rf_full$predict_prob(X_comb[na_rf, , drop = FALSE])

    roc_mx <- .compute_roc_auc_curve(p_vec, oof_mx, w_vec)
    roc_rf <- .compute_roc_auc_curve(p_vec, oof_rf, w_vec)
    auc_mx <- if (is.finite(roc_mx$auc)) roc_mx$auc else 0.5
    auc_rf <- if (is.finite(roc_rf$auc)) roc_rf$auc else 0.5

    # AUC_block-weighted ensemble (default 0.5 * MaxEnt + 0.5 * SpatialRF weighted by AUC_block)
    w_ens_mx <- if ((auc_mx + auc_rf) > 1e-6) auc_mx / (auc_mx + auc_rf) else 0.5
    w_ens_rf <- 1.0 - w_ens_mx
    oof_ens <- pmin(1.0, pmax(0.0, w_ens_mx * oof_mx + w_ens_rf * oof_rf))
    roc_ens <- .compute_roc_auc_curve(p_vec, oof_ens, w_vec)
    auc_ens <- if (is.finite(roc_ens$auc)) roc_ens$auc else round(0.5 * (auc_mx + auc_rf), 4)

    cbi_mx <- compute_continuous_boyce_index(oof_mx, oof_mx[p_vec == 1L], n_bins = 20L)
    cbi_rf <- compute_continuous_boyce_index(oof_rf, oof_rf[p_vec == 1L], n_bins = 20L)
    cbi_ens <- compute_continuous_boyce_index(oof_ens, oof_ens[p_vec == 1L], n_bins = 20L)

    perf_rows[[length(perf_rows) + 1L]] <- data.frame(
      Taxon_or_Morphotype = tg,
      N_Presences = as.integer(n_pres),
      N_Background = as.integer(n_bg),
      MaxEnt_Engine = as.character(mx_full$engine),
      AUC_Block_MaxEnt = round(auc_mx, 4),
      CBI_MaxEnt = round(as.numeric(cbi_mx), 4),
      AUC_Block_SpatialRF = round(auc_rf, 4),
      CBI_SpatialRF = round(as.numeric(cbi_rf), 4),
      AUC_Block_Ensemble = round(auc_ens, 4),
      CBI_Ensemble = round(as.numeric(cbi_ens), 4),
      Top_Climatic_Predictor = as.character(top_clim_pred),
      Top_Edaphic_Predictor = as.character(top_edaph_pred),
      stringsAsFactors = FALSE
    )

    # Project continuous habitat suitability across env_stack via terra::predict()
    rast_mx <- NULL
    rast_rf <- NULL
    rast_ens <- NULL
    if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE)) {
      rast_mx <- tryCatch(
        terra::predict(env_stack, mx_full, fun = function(model, data, ...) model$predict_cloglog(data), na.rm = TRUE),
        error = function(e) NULL
      )
      rast_rf <- tryCatch(
        terra::predict(env_stack, rf_full, fun = function(model, data, ...) model$predict_prob(data), na.rm = TRUE),
        error = function(e) NULL
      )
      if (!is.null(rast_mx) && !is.null(rast_rf)) {
        names(rast_mx) <- "maxent_cloglog"
        names(rast_rf) <- "spatial_rf_suitability"
        rast_ens <- w_ens_mx * rast_mx + w_ens_rf * rast_rf
        names(rast_ens) <- "ensemble_suitability"

        safe_name <- gsub("^_+|_+$", "", tolower(gsub("[^A-Za-z0-9]+", "_", tg)))
        mx_tif_out <- file.path(raster_dir, sprintf("sdm_%s_maxent_cloglog.tif", safe_name))
        ens_tif_out <- file.path(raster_dir, sprintf("sdm_%s_ensemble.tif", safe_name))

        .write_deflate_tif <- function(r_obj, path_target) {
          tmp_t <- paste0(path_target, ".tmp.tif")
          dir.create(dirname(path_target), recursive = TRUE, showWarnings = FALSE)
          tryCatch(
            terra::writeRaster(r_obj, filename = tmp_t, overwrite = TRUE, gdal = gdal_deflate),
            error = function(e) terra::writeRaster(r_obj, filename = tmp_t, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"))
          )
          if (file.exists(tmp_t)) file.rename(tmp_t, path_target)
        }

        if (write_prod_rasters) {
          tryCatch(.write_deflate_tif(rast_mx, mx_tif_out), error = function(e) NULL)
          tryCatch(.write_deflate_tif(rast_ens, ens_tif_out), error = function(e) NULL)
        }
      }
    }

    tg_art <- list(
      taxon = tg,
      pres_mask = pres_mask,
      mx_full = mx_full,
      rf_full = rf_full,
      rast_mx = rast_mx,
      rast_rf = rast_rf,
      rast_ens = rast_ens,
      roc_mx = roc_mx,
      roc_rf = roc_rf,
      roc_ens = roc_ens,
      cbi_mx = cbi_mx,
      cbi_rf = cbi_rf,
      cbi_ens = cbi_ens,
      X_comb = X_comb,
      p_vec = p_vec
    )
    all_sdm_artifacts[[tg]] <- tg_art
    if (is.null(focal_artifacts)) {
      focal_artifacts <- tg_art
    }
  }

  perf_df <- if (length(perf_rows) > 0L) {
    do.call(rbind, perf_rows)
  } else {
    data.frame(
      Taxon_or_Morphotype = character(0), N_Presences = integer(0), N_Background = integer(0),
      MaxEnt_Engine = character(0), AUC_Block_MaxEnt = numeric(0), CBI_MaxEnt = numeric(0),
      AUC_Block_SpatialRF = numeric(0), CBI_SpatialRF = numeric(0), AUC_Block_Ensemble = numeric(0),
      CBI_Ensemble = numeric(0), Top_Climatic_Predictor = character(0), Top_Edaphic_Predictor = character(0),
      stringsAsFactors = FALSE
    )
  }

  perf_csv_out <- file.path(out_dir, "reports", "sdm_model_performance.csv")
  if (write_csv_atomic(perf_df, perf_csv_out, demo = demo)) {
    message(sprintf("Saved SDM performance metrics (%d modeled taxa/morphotypes) to %s.", nrow(perf_df), perf_csv_out))
  }


  # Render canonical 300-DPI PNG figure fig05c_sdm_suitability_and_distribution_map.png
  fig_dir <- file.path(out_dir, "figures")
  if (!(isTRUE(demo) && is_canonical_prod_path(fig_dir))) {
    .ensure_step05_fig_modules()
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(focal_artifacts) && !is.null(focal_artifacts$rast_mx) && !is.null(focal_artifacts$rast_ens)) {
      resp_vars <- intersect(unique(c(top_clim_pred, top_edaph_pred, clim_cand, edaph_cand, kept_vars)), colnames(focal_artifacts$X_comb))
      resp_vars <- resp_vars[seq_len(min(3L, length(resp_vars)))]
      resp_lbls <- format_step05_predictor_labels(resp_vars)

      pres_mat <- focal_artifacts$X_comb[focal_artifacts$p_vec == 1L, , drop = FALSE]
      if (nrow(pres_mat) == 0L) pres_mat <- focal_artifacts$X_comb
      med_row <- apply(pres_mat, 2, stats::median, na.rm = TRUE)
      med_row_df <- as.data.frame(as.list(med_row), check.names = FALSE)
      u_grid <- seq(0, 1, length.out = 100L)
      resp_rows <- lapply(seq_along(resp_vars), function(vi) {
        vn <- resp_vars[vi]
        v_raw <- as.numeric(focal_artifacts$X_comb[[vn]])
        q_lo <- as.numeric(stats::quantile(v_raw, probs = 0.05, na.rm = TRUE))
        q_hi <- as.numeric(stats::quantile(v_raw, probs = 0.95, na.rm = TRUE))
        if (!is.finite(q_lo) || !is.finite(q_hi) || abs(q_hi - q_lo) < 1e-8) {
          q_lo <- min(v_raw, na.rm = TRUE); q_hi <- max(v_raw, na.rm = TRUE)
        }
        eval_df <- med_row_df[rep(1L, 100L), , drop = FALSE]
        x_seq <- q_lo + u_grid * (q_hi - q_lo)
        eval_df[[vn]] <- x_seq
        y_resp <- focal_artifacts$mx_full$predict_cloglog(eval_df)
        data.frame(
          predictor = rep(resp_lbls[vi], 100L),
          x_val = x_seq,
          suitability = as.numeric(y_resp),
          stringsAsFactors = FALSE
        )
      })
      resp_df_plot <- if (length(resp_rows) > 0L) do.call(rbind, resp_rows) else NULL

      roc_parts <- list()
      if (!is.null(focal_artifacts$roc_mx)) roc_parts[[1L]] <- data.frame(fpr = focal_artifacts$roc_mx$fpr, tpr = focal_artifacts$roc_mx$tpr, model = "MaxEnt", stringsAsFactors = FALSE)
      if (!is.null(focal_artifacts$roc_rf)) roc_parts[[2L]] <- data.frame(fpr = focal_artifacts$roc_rf$fpr, tpr = focal_artifacts$roc_rf$tpr, model = "Spatial RF", stringsAsFactors = FALSE)
      if (!is.null(focal_artifacts$roc_ens)) roc_parts[[3L]] <- data.frame(fpr = focal_artifacts$roc_ens$fpr, tpr = focal_artifacts$roc_ens$tpr, model = "Ensemble", stringsAsFactors = FALSE)
      roc_df_plot <- if (length(roc_parts) > 0L) do.call(rbind, Filter(Negate(is.null), roc_parts)) else NULL
      boyce_curve_df <- attr(focal_artifacts$cbi_ens, "boyce_curve")

      valid_multi_arts <- Filter(function(a) !is.null(a$rast_ens), all_sdm_artifacts)
      multi_rast_list <- if (length(valid_multi_arts) > 1L) {
        stats::setNames(lapply(valid_multi_arts, function(a) a$rast_ens), vapply(valid_multi_arts, function(a) as.character(a$taxon), character(1)))
      } else NULL

      render_step05c_sdm_suitability_png(
        df = df,
        focal_rasters = list(maxent = focal_artifacts$rast_mx, ensemble = focal_artifacts$rast_ens),
        multi_rasters = multi_rast_list,
        perf_df = perf_df,
        resp_df = resp_df_plot,
        roc_df = roc_df_plot,
        boyce_df = boyce_curve_df,
        figures_dir = fig_dir,
        focal_taxon = focal_artifacts$taxon
      )
      message(sprintf("Exported 300-DPI Ensemble SDM & Distribution plate -> %s", file.path(fig_dir, "fig05c_sdm_suitability_and_distribution_map.png")))
    } else {
      render_step05c_sdm_suitability_png(
        df = df,
        figures_dir = fig_dir,
        focal_taxon = focal_taxon
      )
    }
  }

  invisible(list(performance = perf_df, focal = focal_artifacts))
}
