if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/morphometrics/step04_allometry_merger_utils.R
# Step 04 Config, Calibration Gating, EM Fallback, Trait Mergers & Batch Audit
# ==============================================================================

load_step04_morphometrics_config <- function(cfg_path = "config/config.yaml") {
  out <- list(correct_allometry = "auto", allometry_r2_threshold = 0.10, num_pcs = 5L, max_k = 8L, random_seed = 42L)
  eff_cfg <- if (!is.null(cfg_path) && length(cfg_path) >= 1L && !is.na(cfg_path[1L]) && nzchar(as.character(cfg_path[1L]))) {
    as.character(cfg_path[1L])
  } else {
    e1 <- Sys.getenv("PLANTMATICS_CONFIG_PATH", unset = ""); e2 <- Sys.getenv("PACKERA_CONFIG_PATH", unset = "")
    if (nzchar(e1)) e1 else if (nzchar(e2)) e2 else "config/config.yaml"
  }
  if (!nzchar(eff_cfg) || !file.exists(eff_cfg)) return(out)
  m_map <- NULL
  if (requireNamespace("yaml", quietly = TRUE)) {
    raw_y <- tryCatch(yaml::read_yaml(eff_cfg), error = function(e) NULL)
    if (is.list(raw_y) && is.list(raw_y$morphometrics)) m_map <- raw_y$morphometrics
  }
  if (is.null(m_map)) {
    lines <- tryCatch(readLines(eff_cfg, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
    m_map <- list(); in_m <- FALSE; m_ind <- NA_integer_; sub_ind <- NA_integer_
    for (ln in lines) {
      if (grepl("^[[:space:]]*(#|$)", ln)) next
      cln <- sub("[[:space:]]+#.*$", "", ln); if (!grepl("^[[:space:]]*[A-Za-z0-9_]+:", cln)) next
      ind <- nchar(sub("^([[:space:]]*).*$", "\\1", cln))
      if (ind == 0L) { if (trimws(sub(":.*$", "", cln)) == "morphometrics") { in_m <- TRUE; next }; if (in_m) break else next }
      if (!in_m) next
      if (is.na(m_ind)) m_ind <- ind
      if (!is.na(sub_ind) && ind > sub_ind) next else sub_ind <- NA_integer_
      key <- trimws(sub(":.*$", "", cln)); rest <- trimws(sub("^[^:]+:", "", cln))
      if (!nzchar(rest)) { sub_ind <- ind; next }
      if (ind == m_ind) { if (grepl("^[\"'].*[\"']$", rest)) rest <- substr(rest, 2L, nchar(rest) - 1L); m_map[[key]] <- rest }
    }
  }
  if (!is.null(m_map$correct_allometry) && length(m_map$correct_allometry) >= 1L) {
    ca <- m_map$correct_allometry[1L]
    out$correct_allometry <- if (is.logical(ca)) (if (isTRUE(ca)) "true" else "false") else { s <- tolower(trimws(as.character(ca))); if (nzchar(s)) s else out$correct_allometry }
  }
  r2_v <- suppressWarnings(as.numeric(m_map$allometry_r2_threshold[1L]))
  if (length(r2_v) == 1L && is.finite(r2_v) && r2_v >= 0.0 && r2_v <= 1.0) out$allometry_r2_threshold <- as.numeric(r2_v)
  np_v <- suppressWarnings(as.integer(round(as.numeric(m_map$num_pcs[1L]))))
  if (length(np_v) == 1L && is.finite(np_v) && np_v >= 2L) out$num_pcs <- np_v
  mk_v <- suppressWarnings(as.integer(round(as.numeric(m_map$max_k[1L]))))
  if (length(mk_v) == 1L && is.finite(mk_v) && mk_v >= 2L) out$max_k <- mk_v
  rs_raw <- if (!is.null(m_map$random_seed)) m_map$random_seed else m_map$seed
  rs_v <- suppressWarnings(as.integer(round(as.numeric(rs_raw[1L]))))
  if (length(rs_v) == 1L && is.finite(rs_v)) out$random_seed <- rs_v
  out
}

is_step04_calibrated_mm <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(logical(0))
  n <- nrow(df); has_cal <- "is_calibrated_mm" %in% names(df); has_ss <- "scale_source" %in% names(df)
  if (!has_cal && !has_ss) return(rep(TRUE, n))
  cal_ok <- if (has_cal) {
    raw_c <- df$is_calibrated_mm
    if (is.character(raw_c) || is.factor(raw_c)) {
      lc <- tolower(trimws(as.character(raw_c))); out_c <- rep(FALSE, n); out_c[lc %in% c("true", "t", "1", "1.0", "yes")] <- TRUE; out_c
    } else !is.na(raw_c) & as.logical(raw_c)
  } else rep(TRUE, n)
  if (has_ss) {
    ss <- tolower(trimws(as.character(df$scale_source))); ss[is.na(ss)] <- ""
    cal_ok <- cal_ok & !(ss %in% c("mp_estimated", "predicted_from_megapixels", "uncalibrated", "tick_lattice_unboxed"))
  }
  as.logical(cal_ok)
}

.parse_step04_cli_fallback <- function(raw_args, defaults) {
  opts <- defaults; if (length(raw_args) == 0L) return(opts)
  expanded <- unlist(lapply(as.character(raw_args), function(tok) {
    if (grepl("^--[A-Za-z0-9_-]+=", tok)) { eq <- regexpr("=", tok, fixed = TRUE)[1L]; c(substr(tok, 1L, eq - 1L), substr(tok, eq + 1L, nchar(tok))) } else tok
  }), use.names = FALSE)
  str_map <- c("-v"="vouchers", "--vouchers"="vouchers", "-e"="harmonics", "--harmonics"="harmonics", "-m"="reproductive", "--reproductive"="reproductive",
               "--correct-allometry"="correct_allometry", "--allometry-report"="allometry_report", "--allometry-pairwise"="allometry_pairwise",
               "--allometry-plot"="allometry_plot", "-f"="output_flags", "--output-flags"="output_flags", "-p"="output_plot", "--output-plot"="output_plot",
               "-r"="output_report", "--output-report"="output_report", "-b"="output_batch_audit", "--output-batch-audit"="output_batch_audit", "--focal-taxon"="focal_taxon")
  int_map <- c("-k"="max_k", "--max-k"="max_k", "-d"="num_pcs", "--num-pcs"="num_pcs", "--permutations"="permutations", "--seed"="seed")
  bool_map <- c("--outline-only"="outline_only", "--use-clustvarsel"="use_clustvarsel", "--intraspecific-mode"="intraspecific_mode")
  i <- 1L; n <- length(expanded)
  while (i <= n) {
    a <- expanded[i]; nxt <- if (i < n) expanded[i + 1L] else ""
    if (a %in% names(bool_map)) { opts[[bool_map[[a]]]] <- TRUE; i <- i + 1L }
    else if (a %in% names(str_map) && i < n) { opts[[str_map[[a]]]] <- nxt; i <- i + 2L }
    else if (a %in% names(int_map) && i < n) { opts[[int_map[[a]]]] <- as.integer(nxt); i <- i + 2L }
    else i <- i + 1L
  }
  opts
}

.compute_exact_mclust_params <- function(K_vec, p_dim, bic_mat = NULL, default_model = "VVI", explicit_n_params = NULL) {
  K_int <- pmax(1L, as.integer(K_vec)); p <- max(1L, as.integer(p_dim))
  if (!is.null(explicit_n_params) && length(explicit_n_params) == length(K_int) && all(is.finite(explicit_n_params))) return(as.numeric(explicit_n_params))
  row_models <- rep(sub("(_RidgelineMerged|_ICLBootCollapsed|_SeqBootStepped)+$", "", as.character(default_model)[1L]), length(K_int))
  if (!is.null(bic_mat) && is.matrix(bic_mat) && !is.null(colnames(bic_mat))) {
    for (idx in seq_along(K_int)) {
      r_nm <- as.character(K_int[idx])
      if (r_nm %in% rownames(bic_mat)) { rv <- suppressWarnings(as.numeric(bic_mat[r_nm, ])); if (any(is.finite(rv))) row_models[idx] <- colnames(bic_mat)[which.max(ifelse(is.finite(rv), rv, -Inf))] }
    }
  }
  has_mclust <- requireNamespace("mclust", quietly = TRUE)
  vapply(seq_along(K_int), function(idx) {
    G <- K_int[idx]; mn <- toupper(trimws(row_models[idx]))
    if (has_mclust) { np_m <- tryCatch(as.numeric(mclust::nMclustParams(modelName = mn, d = p, G = G)), error = function(e) NA_real_); if (is.finite(np_m) && np_m >= 1.0) return(np_m) }
    cov_p <- switch(mn, "E"=1.0, "V"=as.numeric(G), "EII"=1.0, "VII"=as.numeric(G), "EEI"=as.numeric(p), "VEI"=as.numeric(G+p-1L),
                    "EVI"=as.numeric(1L+G*(p-1L)), "VVI"=as.numeric(G*p), "EEE"=as.numeric(p*(p+1L)/2.0), "VVV"=as.numeric(G*p*(p+1L)/2.0),
                    "FULL_COVARIANCE_EM"=as.numeric(G*p*(p+1L)/2.0), as.numeric(G*p))
    (G - 1.0) + G * p + cov_p
  }, numeric(1))
}

fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150, X_eval = NULL, seed = 42L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({ if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
  if (is.null(X_eval)) X_eval <- X
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL; base_seed <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  bic_table <- data.frame(K = seq_len(max_k), LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)
  for (k in seq_len(max_k)) {
    set.seed(base_seed + k)
    if (k == 1L) { mu <- matrix(colMeans(X), nrow = 1L); sigma <- list(stats::cov(X) + diag(1e-6, p)); pi_k <- c(1.0) } else {
      km <- stats::kmeans(X, centers = k, nstart = 10L, iter.max = 30L); mu <- km$centers; pi_k <- as.vector(table(factor(km$cluster, levels = seq_len(k)))) / N
      sigma <- lapply(seq_len(k), function(c) { sx <- X[km$cluster == c, , drop = FALSE]; if (nrow(sx) > p) stats::cov(sx) + diag(1e-5, p) else diag(1e-2, p) })
    }
    loglik_old <- -Inf; resp <- matrix(0, nrow = N, ncol = k)
    for (iter in seq_len(max_iter)) {
      log_dens <- matrix(0, nrow = N, ncol = k)
      for (c in seq_len(k)) {
        dx <- t(t(X) - mu[c, ]); si <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        log_dens[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + rowSums((dx %*% si) * dx))
      }
      max_l <- apply(log_dens, 1L, max); lse <- max_l + log(rowSums(exp(log_dens - max_l))); resp <- exp(log_dens - lse); loglik <- sum(lse)
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik; N_k <- colSums(resp); pi_k <- N_k / N
      for (c in seq_len(k)) if (N_k[c] > 1e-6) { mu[c, ] <- colSums(resp[, c] * X) / N_k[c]; dx <- t(t(X) - mu[c, ]); sigma[[c]] <- (t(dx) %*% (resp[, c] * dx)) / N_k[c] + diag(1e-5, p) }
    }
    n_params <- (k - 1L) + k * p + k * (p * (p + 1L) / 2); bic_val <- 2 * loglik - n_params * log(N)
    bic_table$LogLik[k] <- round(loglik, 2); bic_table$n_params[k] <- n_params; bic_table$BIC[k] <- round(bic_val, 2)
    if (bic_val > best_bic) {
      best_bic <- bic_val; ld_ev <- matrix(0, nrow = nrow(X_eval), ncol = k)
      for (c in seq_len(k)) { dx_e <- t(t(X_eval) - mu[c, ]); si_e <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]]))); ld_ev[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + rowSums((dx_e %*% si_e) * dx_e)) }
      ml_e <- apply(ld_ev, 1L, max); resp_ev <- exp(ld_ev - (ml_e + log(rowSums(exp(ld_ev - ml_e)))))
      best_model <- list(k = k, classification = apply(resp_ev, 1L, which.max), uncertainty = 1 - apply(resp_ev, 1L, max), z = resp_ev, bic = bic_val)
    }
  }
  list(best = best_model, bic_table = bic_table)
}

.backfill_cols_by_individual <- function(agg_df, src_pid, src_cat, cols, single_cat, idx_pid, idx_cat, as_num = FALSE) {
  for (cl in cols) {
    cur <- if (cl %in% names(agg_df)) (if (as_num) suppressWarnings(as.numeric(agg_df[[cl]])) else agg_df[[cl]]) else rep(if (as_num) NA_real_ else NA, nrow(agg_df))
    if (!is.null(src_pid)) { hit_p <- is.na(cur) & !is.na(idx_pid); if (any(hit_p)) cur[hit_p] <- src_pid[[cl]][idx_pid[hit_p]] }
    hit_c <- is.na(cur) & !is.na(idx_cat) & single_cat; if (any(hit_c)) cur[hit_c] <- src_cat[[cl]][idx_cat[hit_c]]
    agg_df[[cl]] <- cur
  }
  agg_df
}

.resolve_ind_join_keys <- function(agg_df, src_df) {
  src_df$catalogNumber <- as.character(src_df$catalogNumber); agg_df$catalogNumber <- as.character(agg_df$catalogNumber)
  has_pid <- ("plant_individual_id" %in% names(src_df)) && ("plant_individual_id" %in% names(agg_df))
  pid_clean <- function(v) { p <- ifelse(is.na(v), 0L, suppressWarnings(as.integer(v))); p[is.na(p)] <- 0L; p }
  if (!has_pid) return(list(has_pid = FALSE, spid = NULL, apid = NULL, single_cat = rep(TRUE, nrow(agg_df))))
  spid <- pid_clean(src_df$plant_individual_id); apid <- pid_clean(agg_df$plant_individual_id)
  n_pids <- tapply(spid, src_df$catalogNumber, function(z) length(unique(z)))
  list(has_pid = TRUE, spid = spid, apid = apid, single_cat = !(agg_df$catalogNumber %in% names(n_pids[n_pids > 1L])))
}

merge_manifest_petiole_traits_by_individual <- function(agg_df, man_cand) {
  if (is.null(agg_df) || nrow(agg_df) == 0L || is.null(man_cand) || !nzchar(man_cand) || !file.exists(man_cand)) return(agg_df)
  mm_cols <- c("petiole_length_mm", "blade_length_mm")
  phot_cols <- c("blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta", "blade_cie_a_delta", "blade_cie_b_delta", "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta")
  target_pcols <- c(mm_cols, "petiole_to_blade_ratio", phot_cols)
  need_cols <- target_pcols[vapply(target_pcols, function(cl) !(cl %in% names(agg_df)) || any(is.na(agg_df[[cl]])), logical(1))]
  if (length(need_cols) == 0L) return(agg_df)
  man_df <- tryCatch(utils::read.csv(man_cand, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(man_df) || nrow(man_df) == 0L || !("catalogNumber" %in% names(man_df))) return(agg_df)
  p_cols <- intersect(target_pcols, names(man_df)); if (length(p_cols) == 0L) return(agg_df)
  cal_m <- is_step04_calibrated_mm(man_df)
  for (mc in intersect(mm_cols, p_cols)) { mv <- suppressWarnings(as.numeric(man_df[[mc]])); mv[!cal_m] <- NA_real_; man_df[[mc]] <- mv }
  if ("leaf_zone" %in% names(man_df)) {
    lz_low <- tolower(trimws(as.character(man_df$leaf_zone))); lz_low[is.na(lz_low) | !nzchar(lz_low)] <- "basal"
    is_basal_primary <- lz_low != "cauline"
    if (any(is_basal_primary)) {
      mp_tmp <- if ("plant_individual_id" %in% names(man_df)) ifelse(is.na(man_df$plant_individual_id), 0L, suppressWarnings(as.integer(man_df$plant_individual_id))) else rep(0L, nrow(man_df))
      mp_tmp[is.na(mp_tmp)] <- 0L; grp_k <- if ("plant_individual_id" %in% names(man_df)) paste0(man_df$catalogNumber, "__", mp_tmp) else as.character(man_df$catalogNumber)
      has_basal_grp <- tapply(is_basal_primary, grp_k, any); man_df <- man_df[which(is_basal_primary | !has_basal_grp[grp_k]), , drop = FALSE]
    }
  }
  q75_fn <- function(x) { v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]; if (length(v) > 0L) round(stats::quantile(v, probs = 0.75, na.rm = TRUE, names = FALSE), 4) else NA_real_ }
  med_fn <- function(x) { v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]; if (length(v) > 0L) round(stats::median(v), 4) else NA_real_ }
  jk <- .resolve_ind_join_keys(agg_df, man_df); q75_cols <- setdiff(p_cols, phot_cols); med_cols <- intersect(p_cols, phot_cols)
  agg_by_fn <- function(by_list) {
    parts <- list()
    if (length(q75_cols) > 0L) parts[[length(parts) + 1L]] <- stats::aggregate(man_df[, q75_cols, drop = FALSE], by = by_list, FUN = q75_fn)
    if (length(med_cols) > 0L) parts[[length(parts) + 1L]] <- stats::aggregate(man_df[, med_cols, drop = FALSE], by = by_list, FUN = med_fn)
    if (length(parts) == 1L) parts[[1L]] else merge(parts[[1L]], parts[[2L]], by = names(by_list), all = TRUE, sort = FALSE)
  }
  agg_cat <- agg_by_fn(list(catalogNumber = as.character(man_df$catalogNumber)))
  agg_pid <- if (jk$has_pid) agg_by_fn(list(catalogNumber = as.character(man_df$catalogNumber), plant_individual_id = jk$spid)) else NULL
  idx_pid <- if (!is.null(agg_pid)) match(paste0(agg_df$catalogNumber, "__", jk$apid), paste0(agg_pid$catalogNumber, "__", agg_pid$plant_individual_id)) else rep(NA_integer_, nrow(agg_df))
  .backfill_cols_by_individual(agg_df, agg_pid, agg_cat, intersect(need_cols, p_cols), jk$single_cat, idx_pid, match(as.character(agg_df$catalogNumber), agg_cat$catalogNumber), as_num = TRUE)
}

merge_reproductive_traits_by_individual <- function(agg_df, repro_path) {
  if (is.null(agg_df) || nrow(agg_df) == 0L || is.null(repro_path) || !nzchar(repro_path) || !file.exists(repro_path)) return(agg_df)
  repro_df <- tryCatch(utils::read.csv(repro_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(repro_df) || nrow(repro_df) == 0L || !("catalogNumber" %in% names(repro_df))) return(agg_df)
  if ("plant_height_mm" %in% names(repro_df) && !("plant_height_cm" %in% names(repro_df))) repro_df$plant_height_cm <- round(suppressWarnings(as.numeric(repro_df$plant_height_mm)) / 10.0, 2)
  if (!("capitula_count" %in% names(repro_df)) && ("capitulum_count" %in% names(repro_df))) repro_df$capitula_count <- repro_df$capitulum_count
  if (any(c("capitula_count", "corymb_head_count_est") %in% names(repro_df))) {
    cc_v <- if ("capitula_count" %in% names(repro_df)) suppressWarnings(as.numeric(repro_df$capitula_count)) else rep(NA_real_, nrow(repro_df))
    ce_v <- if ("corymb_head_count_est" %in% names(repro_df)) suppressWarnings(as.numeric(repro_df$corymb_head_count_est)) else rep(NA_real_, nrow(repro_df))
    cc_v[!is.finite(cc_v)] <- NA_real_; ce_v[!is.finite(ce_v)] <- NA_real_; cc_eff <- ifelse(is.finite(cc_v) | is.finite(ce_v), pmax(cc_v, ce_v, na.rm = TRUE), NA_real_)
    l1_eff <- ifelse(is.finite(cc_eff) & cc_eff >= 1, round(log1p(cc_eff), 5), NA_real_)
    if (!("log1p_capitulum_count" %in% names(repro_df))) repro_df$log1p_capitulum_count <- l1_eff else { cur_l1 <- suppressWarnings(as.numeric(repro_df$log1p_capitulum_count)); cur_l1[is.finite(cc_eff)] <- l1_eff[is.finite(cc_eff)]; repro_df$log1p_capitulum_count <- cur_l1 }
  }
  cand_r <- c("AR_inv", "capitula_count", "corymb_head_count_est", "anthesis_head_count", "fruiting_head_count", "Flower_MANY", "Seed_Fruit_MANY",
              "log1p_capitulum_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm",
              "plant_height_mm", "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem", "petiole_to_blade_ratio")
  r_val_cols <- intersect(names(repro_df), cand_r)
  r_val_cols <- r_val_cols[vapply(r_val_cols, function(cl) !(cl %in% names(agg_df)) || any(is.na(agg_df[[cl]])), logical(1))]
  if (length(r_val_cols) == 0L) return(agg_df)
  jk <- .resolve_ind_join_keys(agg_df, repro_df)
  r_sub_pid <- if (jk$has_pid) repro_df[!duplicated(paste0(as.character(repro_df$catalogNumber), "__", jk$spid)), , drop = FALSE] else NULL
  idx_pid <- if (!is.null(r_sub_pid)) {
    rp_s <- ifelse(is.na(r_sub_pid$plant_individual_id), 0L, suppressWarnings(as.integer(r_sub_pid$plant_individual_id))); rp_s[is.na(rp_s)] <- 0L
    match(paste0(agg_df$catalogNumber, "__", jk$apid), paste0(as.character(r_sub_pid$catalogNumber), "__", rp_s))
  } else rep(NA_integer_, nrow(agg_df))
  r_sub_cat <- repro_df[!duplicated(repro_df$catalogNumber), , drop = FALSE]
  .backfill_cols_by_individual(agg_df, r_sub_pid, r_sub_cat, r_val_cols, jk$single_cat, idx_pid, match(as.character(agg_df$catalogNumber), as.character(r_sub_cat$catalogNumber)), as_num = FALSE)
}

.update_step04_allometry_pcs <- function(morpho_df, sub_df, gold_df, valid_idx, sym_names, active_cols, pca_resids, ref_mean_cs, eff_weight) {
  all_sym_idx <- which(stats::complete.cases(morpho_df[, sym_names, drop = FALSE]))
  if (length(all_sym_idx) == 0L) return(morpho_df)
  proj_pcs <- scale(as.matrix(morpho_df[all_sym_idx, sym_names[active_cols], drop = FALSE]), center = pca_resids$center, scale = FALSE) %*% pca_resids$rotation
  pc_cols <- grep("^PC[0-9]+$", names(morpho_df), value = TRUE)
  if (length(pc_cols) > 1L) pc_cols <- pc_cols[order(as.integer(gsub("[^0-9]", "", pc_cols)))]
  for (p in seq_along(pc_cols)) morpho_df[all_sym_idx, pc_cols[p]] <- if (p <= ncol(proj_pcs)) round(proj_pcs[, p], 6) else 0.0
  morpho_df
}

.update_step04_allometry_scalars <- function(morpho_df, sub_df, gold_df, valid_idx, ref_mean_cs) {
  sc_cands <- intersect(c("ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "solidity", "q25_solidity", "dissection_index", "aspect_ratio", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio"), names(morpho_df))
  for (sc_nm in sc_cands) {
    yv_g <- suppressWarnings(as.numeric(gold_df[[sc_nm]])); ok_g <- is.finite(yv_g) & is.finite(gold_df$log_centsize)
    if (sum(ok_g) < 10L || stats::sd(yv_g[ok_g]) <= 1e-6) next
    g_sub <- gold_df[ok_g, , drop = FALSE]; g_sub$.y_sc <- yv_g[ok_g]; g_sub$scientificName <- droplevels(factor(g_sub$scientificName)); p_int_sc <- 1.0
    if (nlevels(g_sub$scientificName) >= 2L) {
      fi <- tryCatch(stats::lm(.y_sc ~ log_centsize * scientificName, data = g_sub, weights = eff_weight), error = function(e) NULL)
      if (!is.null(fi)) { ai <- stats::anova(fi); ii <- grep(":", rownames(ai))[1L]; if (!is.na(ii) && is.finite(ai[["Pr(>F)"]][ii])) p_int_sc <- ai[["Pr(>F)"]][ii] }
    }
    if (p_int_sc > 0.05) {
      fs <- tryCatch(if (nlevels(g_sub$scientificName) >= 2L) stats::lm(.y_sc ~ log_centsize + scientificName, data = g_sub, weights = eff_weight) else stats::lm(.y_sc ~ log_centsize, data = g_sub, weights = eff_weight), error = function(e) NULL)
      if (!is.null(fs)) {
        b_sc <- unname(stats::coef(fs)["log_centsize"]); p_s_sc <- tryCatch(summary(fs)$coefficients["log_centsize", "Pr(>|t|)"], error = function(e) 1.0)
        if (is.finite(b_sc) && is.finite(p_s_sc) && p_s_sc < 0.05) { yv_a <- suppressWarnings(as.numeric(morpho_df[[sc_nm]][valid_idx])); fv <- is.finite(yv_a); if (any(fv)) morpho_df[[sc_nm]][valid_idx[fv]] <- round(yv_a[fv] - b_sc * (sub_df$log_centsize[fv] - ref_mean_cs), 6) }
      }
    }
  }
  attr(morpho_df, "allometry_corrected_scalars") <- TRUE; morpho_df
}

audit_institutional_batch_effects <- function(morpho_df, output_csv = "outputs/reports/institutional_batch_effect_audit.csv", min_obs = 15, selected_features = NULL) {
  if (!("institutionCode" %in% names(morpho_df))) return(NULL)
  dyn_sel <- if (!is.null(selected_features) && length(selected_features) > 0L) intersect(as.character(selected_features), names(morpho_df)) else discover_sorted_morph_cols(morpho_df, "^PC[0-9]+$")
  test_feats <- unique(c(intersect(c("PC1", "PC2"), names(morpho_df)), dyn_sel))
  test_feats <- test_feats[vapply(test_feats, function(fc) is.numeric(morpho_df[[fc]]) && sum(is.finite(morpho_df[[fc]])) >= 10L, logical(1))]
  valid_df <- morpho_df[!is.na(morpho_df$institutionCode) & nzchar(trimws(morpho_df$institutionCode)), , drop = FALSE]
  clean_b <- (!("authority_weight" %in% names(valid_df)) | (!is.na(valid_df$authority_weight) & valid_df$authority_weight > 0)) &
    (!("is_passive" %in% names(valid_df)) | is.na(valid_df$is_passive) | !as.logical(valid_df$is_passive)) &
    (!("morph_tier_source" %in% names(valid_df)) | is.na(valid_df$morph_tier_source) | !(valid_df$morph_tier_source %in% c("Cauline_Only_Fallback", "Transitional_Passive")))
  if (sum(clean_b, na.rm = TRUE) >= 10L && length(unique(valid_df$institutionCode[clean_b])) >= 2L) valid_df <- valid_df[clean_b, , drop = FALSE]
  if (nrow(valid_df) < 10L || length(test_feats) == 0L) return(NULL)
  inst_counts <- table(valid_df$institutionCode); rep_insts <- names(inst_counts[inst_counts >= min_obs])
  audit_subset <- if (length(rep_insts) >= 2L) valid_df[valid_df$institutionCode %in% rep_insts, , drop = FALSE] else valid_df
  audit_subset$institutionCode <- droplevels(factor(audit_subset$institutionCode))
  if (!("species_standardized" %in% names(audit_subset))) audit_subset$species_standardized <- if ("scientificName" %in% names(audit_subset)) audit_subset$scientificName else "Unassigned_Packera"
  audit_subset$species_standardized <- droplevels(factor(audit_subset$species_standardized))
  if (nlevels(audit_subset$institutionCode) < 2L) return(NULL)
  fits <- list(); rows_out <- list(); raw_ps <- numeric(0); pre_eb_map <- attr(morpho_df, "pre_eb_batch_partial_eta_sq")
  for (feat in test_feats) {
    sub_f <- audit_subset[is.finite(audit_subset[[feat]]), , drop = FALSE]
    sub_f$institutionCode <- droplevels(sub_f$institutionCode); sub_f$species_standardized <- droplevels(sub_f$species_standardized)
    if (nrow(sub_f) < 6L || nlevels(sub_f$institutionCode) < 2L || stats::var(sub_f[[feat]], na.rm = TRUE) <= 1e-12) next
    w_bf <- (if ("authority_weight" %in% names(sub_f)) suppressWarnings(as.numeric(sub_f$authority_weight)) else rep(1.0, nrow(sub_f))) * (if ("spatial_balance_weight" %in% names(sub_f)) suppressWarnings(as.numeric(sub_f$spatial_balance_weight)) else rep(1.0, nrow(sub_f)))
    w_bf[!is.finite(w_bf) | w_bf <= 0] <- 1.0
    rg_raw <- if ("regional_group" %in% names(sub_f)) trimws(as.character(sub_f$regional_group)) else rep("", nrow(sub_f))
    rg_ok <- !is.na(rg_raw) & nzchar(rg_raw) & !(rg_raw %in% c("Unassigned", "Unknown", "NA")); rg_clean <- unique(rg_raw[rg_ok])
    sub_f$.rg_clean <- factor(ifelse(rg_ok, rg_raw, "Unassigned")); has_sp <- nlevels(sub_f$species_standardized) >= 2L
    cand_pairs <- list()
    if (length(rg_clean) >= 2L) {
      rhs_r <- if (has_sp) "species_standardized + .rg_clean" else ".rg_clean"
      cand_pairs[[1L]] <- list(red = stats::as.formula(paste(feat, "~", rhs_r)), full = stats::as.formula(paste(feat, "~", rhs_r, "+ institutionCode")), min_df = 3L)
    }
    if (has_sp) cand_pairs[[length(cand_pairs) + 1L]] <- list(red = stats::as.formula(paste(feat, "~ species_standardized")), full = stats::as.formula(paste(feat, "~ species_standardized + institutionCode")), min_df = 2L)
    cand_pairs[[length(cand_pairs) + 1L]] <- list(red = stats::as.formula(paste(feat, "~ 1")), full = stats::as.formula(paste(feat, "~ institutionCode")), min_df = 1L)
    sel_pair <- NULL
    for (cp in cand_pairs) if (qr(stats::model.matrix(cp$full, data = sub_f))$rank > qr(stats::model.matrix(cp$red, data = sub_f))$rank && (nrow(sub_f) - qr(stats::model.matrix(cp$full, data = sub_f))$rank) >= cp$min_df) { sel_pair <- cp; break }
    if (is.null(sel_pair)) next
    fit_red <- stats::lm(sel_pair$red, data = sub_f, weights = w_bf); fit_full <- stats::lm(sel_pair$full, data = sub_f, weights = w_bf)
    aov_cmp <- stats::anova(fit_red, fit_full); s_f <- summary(fit_full)
    df_i <- as.integer(aov_cmp$Df[2L]); df_r <- as.integer(aov_cmp$Res.Df[2L]); ss_i <- max(0.0, as.numeric(aov_cmp[["Sum of Sq"]][2L])); ss_r <- max(0.0, as.numeric(aov_cmp$RSS[2L]))
    if (is.na(df_i) || df_i <= 0L || is.na(df_r) || df_r <= 0L) next
    f_val <- as.numeric(aov_cmp$F[2L]); pv_f <- as.numeric(aov_cmp[["Pr(>F)"]][2L]); if (!is.finite(f_val)) f_val <- 0.0; if (!is.finite(pv_f)) pv_f <- 1.0
    r2_f <- ss_i / max(ss_i + ss_r, 1e-12); pre_eta_f <- if (!is.null(pre_eb_map) && feat %in% names(pre_eb_map) && is.finite(pre_eb_map[[feat]])) as.numeric(pre_eb_map[[feat]]) else r2_f
    fits[[feat]] <- fit_full; raw_ps <- c(raw_ps, pv_f)
    rows_out[[length(rows_out) + 1L]] <- data.frame(
      Response = c(feat, feat), Term = c("institutionCode", "Residuals"), Df = c(df_i, df_r), Sum_Sq = round(c(ss_i, ss_r), 4),
      Mean_Sq = round(c(ss_i / df_i, ss_r / df_r), 4), F_value = c(round(f_val, 4), NA), p_value = c(signif(pv_f, 4), NA), P_value = c(signif(pv_f, 4), NA),
      p_value_bh = c(NA_real_, NA_real_), R_squared = round(c(r2_f, NA), 4), Adj_R_squared = round(c(s_f$adj.r.squared, NA), 4),
      partial_eta_sq = round(c(r2_f, NA), 4), pre_eb_partial_eta_sq = round(c(pre_eta_f, NA), 4),
      Batch_Effect_Status = c(ifelse(r2_f < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA), stringsAsFactors = FALSE)
  }
  if (length(rows_out) == 0L) return(NULL)
  anova_df <- do.call(rbind, rows_out); inst_rows <- which(anova_df$Term == "institutionCode"); q_bh <- stats::p.adjust(raw_ps, method = "BH")
  anova_df$p_value_bh[inst_rows] <- signif(q_bh, 4)
  anova_df$Batch_Effect_Status[inst_rows] <- ifelse(anova_df$partial_eta_sq[inst_rows] >= 0.05 & q_bh < 0.05, "Exceeds 0.05 (q_BH < 0.05)", ifelse(anova_df$partial_eta_sq[inst_rows] < 0.05, "Pass (<0.05)", "Pass (q_BH >= 0.05)"))
  if (!is.null(output_csv) && nzchar(output_csv)) { dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE); tmp_ba <- paste0(output_csv, ".", Sys.getpid(), ".tmp"); write.csv(anova_df, file = tmp_ba, row.names = FALSE); file.rename(tmp_ba, output_csv) }
  invisible(list(lm_pc1 = fits[["PC1"]], lm_pc2 = fits[["PC2"]], models = fits, anova_table = anova_df))
}
