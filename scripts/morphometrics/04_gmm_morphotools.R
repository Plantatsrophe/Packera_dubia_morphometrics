#!/usr/bin/env Rscript
# Script: 04_gmm_morphotools.R - Packera dubia Morphometrics Pipeline (NCU)
suppressPackageStartupMessages({
  for (pkg in c("MASS", "mclust", "clustvarsel", "MorphoTools2", "RRPP", "vegan", "dplyr", "readr", "tibble", "ggplot2", "optparse")) {
    if (requireNamespace(pkg, quietly = TRUE)) library(pkg, character.only = TRUE)
  }
})

# 1. CLI Argument Parsing
parse_args_robust <- function() {
  option_list <- list(
    optparse::make_option(c("-v", "--vouchers"), type = "character", default = "data/tables/curated_vouchers.csv"),
    optparse::make_option(c("-e", "--harmonics"), type = "character", default = "data/tables/leaf_efa_harmonics.csv"),
    optparse::make_option(c("-m", "--reproductive"), type = "character", default = "data/tables/voucher_reproductive_metrics.csv"),
    optparse::make_option(c("--outline-only"), action = "store_true", default = FALSE),
    optparse::make_option(c("--correct-allometry"), type = "character", default = "auto"),
    optparse::make_option(c("--allometry-report"), type = "character", default = "outputs/reports/allometry_significance_audit.csv"),
    optparse::make_option(c("--allometry-pairwise"), type = "character", default = "outputs/reports/allometry_pairwise_trajectories.csv"),
    optparse::make_option(c("--allometry-plot"), type = "character", default = "outputs/figures/shape_allometry_regression.pdf"),
    optparse::make_option(c("--use-clustvarsel"), action = "store_true", default = FALSE),
    optparse::make_option(c("-f", "--output-flags"), type = "character", default = "data/tables/morphometric_flags.csv"),
    optparse::make_option(c("-p", "--output-plot"), type = "character", default = "outputs/figures/cda_passive_projection.pdf"),
    optparse::make_option(c("-r", "--output-report"), type = "character", default = "outputs/reports/gmm_bayes_factors_summary.csv"),
    optparse::make_option(c("-b", "--output-batch-audit"), type = "character", default = "outputs/reports/institutional_batch_effect_audit.csv"),
    optparse::make_option(c("-k", "--max-k"), type = "integer", default = 8), optparse::make_option(c("-d", "--num-pcs"), type = "integer", default = 5),
    optparse::make_option(c("--permutations"), type = "integer", default = 999), optparse::make_option(c("--seed"), type = "integer", default = 42),
    optparse::make_option(c("--intraspecific-mode"), action = "store_true", default = FALSE), optparse::make_option(c("--focal-taxon"), type = "character", default = "")
  )
  if (requireNamespace("optparse", quietly = TRUE)) { opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list)); names(opts) <- gsub("-", "_", names(opts)); return(opts) }
  list(vouchers = "data/tables/curated_vouchers.csv", harmonics = "data/tables/leaf_efa_harmonics.csv", reproductive = "data/tables/voucher_reproductive_metrics.csv",
       outline_only = FALSE, correct_allometry = "auto", allometry_report = "outputs/reports/allometry_significance_audit.csv",
       allometry_pairwise = "outputs/reports/allometry_pairwise_trajectories.csv", allometry_plot = "outputs/figures/shape_allometry_regression.pdf",
       use_clustvarsel = FALSE, output_flags = "data/tables/morphometric_flags.csv", output_plot = "outputs/figures/cda_passive_projection.pdf",
       output_report = "outputs/reports/gmm_bayes_factors_summary.csv", output_batch_audit = "outputs/reports/institutional_batch_effect_audit.csv",
       max_k = 8, num_pcs = 5, permutations = 999, seed = 42, intraspecific_mode = FALSE, focal_taxon = "")
}

# 2. Circular Phenology & Taxonomic Concept Standardization
calculate_circular_phenology <- function(df) {
  if ("doy" %in% names(df)) {
    doy_val <- suppressWarnings(as.numeric(df$doy)); valid <- !is.na(doy_val) & doy_val >= 1 & doy_val <= 366
    df$pheno_sin <- NA_real_; df$pheno_cos <- NA_real_
    df$pheno_sin[valid] <- round(sin(2 * pi * doy_val[valid] / 365.25), 6); df$pheno_cos[valid] <- round(cos(2 * pi * doy_val[valid] / 365.25), 6)
  }
  return(df)
}

if (file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
for (.aw_p in c("scripts/morphometrics/authority_weighting_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "authority_weighting_utils.R"))) {
  if (!is.na(.aw_p) && nzchar(.aw_p) && file.exists(.aw_p)) { source(.aw_p); break }
}
.base_std_taxon <- if (exists("standardize_taxon_name", mode = "function")) standardize_taxon_name else NULL
standardize_packera_taxon <- function(species_vec) {
  if (!is.null(.base_std_taxon)) { res <- .base_std_taxon(species_vec); return(ifelse(is.na(res), "Unknown", res)) }
  as.character(species_vec)
}

# Empirical Bayes / Dawid-Skene hierarchical active/passive reference authority weights (Slenker et al. 2022; Dawid & Skene 1979)
compute_authority_precision_weights <- function(df, active_mask = NULL) compute_empirical_bayes_authority_weights(df, mode = "active_passive", active_mask = active_mask)

# 2b. Foliar Allometry Multivariate Testing & Residual Correction
generate_allometry_diagnostic_plot <- function(sub_df, r2_size, p_size, r2_int, output_pdf, p_int = 1.0) {
  if (!exists("render_allometry_diagnostic_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
  if (exists("render_allometry_diagnostic_pdf", mode = "function")) render_allometry_diagnostic_pdf(sub_df, r2_size, p_size, r2_int, output_pdf, p_int = p_int)
  message("Shape allometry regression plot saved -> ", output_pdf)
}

test_and_correct_allometry <- function(morpho_df, correct_mode = "auto", permutations = 999, seed = 42, r2_threshold = 0.10,
                                       output_report = "outputs/reports/allometry_significance_audit.csv",
                                       output_plot = "outputs/figures/shape_allometry_regression.pdf",
                                       output_pairwise = "outputs/reports/allometry_pairwise_trajectories.csv") {
  message("=== Foliar Allometry Assessment (Multivariate Size vs. Shape Covariation) ===")
  sym_names <- intersect(c(paste0("A", seq_len(12)), paste0("D", seq_len(12))), names(morpho_df))
  if (!"log_centsize" %in% names(morpho_df) || length(sym_names) < 4) {
    message("Warning: log_centsize or symmetric harmonics not available. Skipping allometry test.")
    return(list(morpho_matrix = NULL, corrected = FALSE, anova_table = NULL, df = morpho_df))
  }
  cal_mask <- if ("is_calibrated_mm" %in% names(morpho_df)) (!is.na(morpho_df$is_calibrated_mm) & as.logical(morpho_df$is_calibrated_mm)) else rep(TRUE, nrow(morpho_df))
  valid_idx <- which(cal_mask & !is.na(morpho_df$log_centsize) & is.finite(morpho_df$log_centsize) & complete.cases(morpho_df[, sym_names]))
  if (length(valid_idx) < 10) return(list(morpho_matrix = NULL, corrected = FALSE, anova_table = NULL, df = morpho_df))
  sub_df <- morpho_df[valid_idx, , drop = FALSE]; sub_df$scientificName <- factor(standardize_packera_taxon(as.character(if ("species_standardized" %in% names(sub_df)) sub_df$species_standardized else sub_df$scientificName)))
  is_tier1_pure <- if ("morph_tier_source" %in% names(sub_df)) (!is.na(sub_df$morph_tier_source) & sub_df$morph_tier_source == "Tier_1_Pure") else rep(TRUE, nrow(sub_df))
  anc_allom <- select_adaptive_reference_anchors(sub_df, feature_matrix = as.matrix(sub_df[, sym_names[seq_len(min(4L, length(sym_names)))], drop = FALSE]), taxon_col = "scientificName")
  sub_df$is_promoted_bronze <- anc_allom$is_promoted_bronze; sub_df$authority_weight <- compute_authority_precision_weights(sub_df)
  sub_df$spatial_balance_weight <- compute_spatial_balance_weights(sub_df); sub_df$eff_weight <- sub_df$authority_weight
  gold_mask <- anc_allom$active_mask & !is.na(sub_df$scientificName) & is_tier1_pure
  min_sp <- min(2L, length(unique(stats::na.omit(sub_df$scientificName))))
  gold_df <- if (sum(gold_mask) >= 6 && length(unique(stats::na.omit(sub_df$scientificName[gold_mask]))) >= min_sp) sub_df[gold_mask, , drop = FALSE] else if (sum(is_tier1_pure) >= 6) sub_df[is_tier1_pure, , drop = FALSE] else sub_df
  if (all(gold_df$eff_weight <= 0)) gold_df$eff_weight <- 1.0
  gold_df$scientificName <- droplevels(factor(gold_df$scientificName)); Y <- as.matrix(gold_df[, sym_names, drop = FALSE]); set.seed(seed); anova_df <- NULL; pw_df <- NULL
  if (requireNamespace("RRPP", quietly = TRUE) && nlevels(gold_df$scientificName) >= 2) {
    rdf <- RRPP::rrpp.data.frame(Y = Y, log_centsize = gold_df$log_centsize, scientificName = gold_df$scientificName)
    fit_rrpp <- tryCatch(RRPP::lm.rrpp(Y ~ log_centsize * scientificName, data = rdf, iter = permutations, print.progress = FALSE), error = function(e) NULL)
    if (!is.null(fit_rrpp)) {
      anova_df <- as.data.frame(anova(fit_rrpp)$table); anova_df$Term <- rownames(anova_df); rownames(anova_df) <- NULL
      if ("SS" %in% names(anova_df)) anova_df$SumOfSqs <- anova_df$SS; if ("Rsq" %in% names(anova_df)) anova_df$R2 <- anova_df$Rsq
      if ("Pr(>F)" %in% names(anova_df)) anova_df$p_value <- anova_df[["Pr(>F)"]]
      pw <- tryCatch(RRPP::pairwise(fit_rrpp, groups = gold_df$scientificName, covariate = gold_df$log_centsize), error = function(e) NULL)
      if (!is.null(pw)) {
        vc_s <- tryCatch(as.data.frame(summary(pw, test.type = "VC", angle.type = "deg")$summary.table), error = function(e) NULL)
        dt_s <- tryCatch(as.data.frame(summary(pw, test.type = "dist")$summary.table), error = function(e) NULL)
        if (!is.null(vc_s) && nrow(vc_s) > 0) {
          pw_df <- data.frame(Pair = rownames(vc_s), r = round(vc_s[[1]], 4), angle_deg = round(vc_s[[2]], 3), Z_angle = round(vc_s[["Z"]], 4), p_angle = round(vc_s[[ncol(vc_s)]], 5), stringsAsFactors = FALSE)
          if (!is.null(dt_s) && nrow(dt_s) == nrow(vc_s)) { pw_df$slope_dist <- round(dt_s[[1]], 4); pw_df$p_dist <- round(dt_s[[ncol(dt_s)]], 5) }
        }
      }
    }
  }
  if (is.null(anova_df) && requireNamespace("vegan", quietly = TRUE) && nlevels(gold_df$scientificName) >= 2) {
    allometry_fit <- tryCatch(vegan::adonis2(Y ~ log_centsize * scientificName, data = gold_df, permutations = permutations, method = "euclidean", by = "terms"), error = function(e) NULL)
    if (!is.null(allometry_fit)) {
      anova_df <- as.data.frame(allometry_fit); anova_df$Term <- rownames(anova_df); rownames(anova_df) <- NULL
      if ("F.Model" %in% names(anova_df)) anova_df$F <- anova_df$F.Model; if ("Pr(>F)" %in% names(anova_df)) anova_df$p_value <- anova_df[["Pr(>F)"]]
    }
  }
  if (is.null(anova_df)) {
    Y_cent <- scale(Y, center = TRUE, scale = FALSE); SS_tot <- sum(Y_cent^2); N_obs <- nrow(gold_df)
    fit1 <- stats::lm(Y_cent ~ log_centsize, data = gold_df, weights = eff_weight); SS_size <- sum(stats::fitted(fit1)^2); df_size <- fit1$rank - 1
    if (nlevels(gold_df$scientificName) >= 2) {
      fit2 <- stats::lm(Y_cent ~ log_centsize + scientificName, data = gold_df, weights = eff_weight); SS_12 <- sum(stats::fitted(fit2)^2); SS_taxa <- SS_12 - SS_size; df_taxa <- fit2$rank - fit1$rank
      fit3 <- stats::lm(Y_cent ~ log_centsize * scientificName, data = gold_df, weights = eff_weight); SS_123 <- sum(stats::fitted(fit3)^2); SS_int <- SS_123 - SS_12; df_int <- fit3$rank - fit2$rank
      SS_res <- SS_tot - SS_123; df_res <- N_obs - fit3$rank; MS_res <- SS_res / max(df_res, 1)
      F_size <- (SS_size / max(df_size, 1)) / MS_res; F_taxa <- (SS_taxa / max(df_taxa, 1)) / MS_res; F_int <- (SS_int / max(df_int, 1)) / MS_res
      p_taxa <- round(stats::pf(F_taxa, df_taxa, df_res, lower.tail = FALSE), 5); p_int <- round(stats::pf(F_int, df_int, df_res, lower.tail = FALSE), 5)
    } else {
      SS_taxa <- 0.0; df_taxa <- 0L; SS_int <- 0.0; df_int <- 0L; SS_res <- SS_tot - SS_size; df_res <- N_obs - fit1$rank; MS_res <- SS_res / max(df_res, 1)
      F_size <- (SS_size / max(df_size, 1)) / MS_res; F_taxa <- NA_real_; F_int <- NA_real_; p_taxa <- NA_real_; p_int <- NA_real_
    }
    anova_df <- data.frame(
      Term = c("log_centsize", "scientificName", "log_centsize:scientificName", "Residual", "Total"),
      Df = c(df_size, df_taxa, df_int, df_res, N_obs - 1), SumOfSqs = round(c(SS_size, SS_taxa, SS_int, SS_res, SS_tot), 4),
      R2 = round(c(SS_size / SS_tot, SS_taxa / SS_tot, SS_int / SS_tot, SS_res / SS_tot, 1.0), 6),
      F = c(round(F_size, 4), round(F_taxa, 4), round(F_int, 4), NA_real_, NA_real_),
      p_value = c(round(stats::pf(F_size, df_size, df_res, lower.tail = FALSE), 5), p_taxa, p_int, NA_real_, NA_real_), stringsAsFactors = FALSE)
  }
  if (!is.null(output_pairwise) && nzchar(output_pairwise) && !is.null(pw_df)) {
    dir.create(dirname(output_pairwise), recursive = TRUE, showWarnings = FALSE); tmp_pw <- paste0(output_pairwise, ".tmp"); write.csv(pw_df, file = tmp_pw, row.names = FALSE); file.rename(tmp_pw, output_pairwise)
  }
  anova_df$Term <- gsub("^Residuals$", "Residual", anova_df$Term)
  if ("Total" %in% anova_df$Term) anova_df$R2[anova_df$Term == "Total"] <- 1.0
  anova_df <- anova_df[, intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(anova_df)), drop = FALSE]
  if (!is.null(output_report) && nzchar(output_report)) {
    dir.create(dirname(output_report), recursive = TRUE, showWarnings = FALSE); tmp_rep <- paste0(output_report, ".tmp"); write.csv(anova_df, file = tmp_rep, row.names = FALSE); file.rename(tmp_rep, output_report)
    message("Allometry significance audit report exported -> ", output_report)
  }
  size_row <- anova_df[grepl("^log_centsize$", anova_df$Term, ignore.case = TRUE), ]
  r2_size <- if (nrow(size_row) > 0) size_row$R2[1] else 0.0; p_size <- if (nrow(size_row) > 0 && !is.na(size_row$p_value[1])) size_row$p_value[1] else 0.0
  int_row <- anova_df[grepl("log_centsize:scientificName|scientificName:log_centsize", anova_df$Term, ignore.case = TRUE), ]
  r2_int <- if (nrow(int_row) > 0) r2_int <- int_row$R2[1] else 0.0; p_int <- if (nrow(int_row) > 0 && !is.na(int_row$p_value[1])) int_row$p_value[1] else 1.0
  authority_weight <- gold_df$authority_weight; spatial_balance_weight <- gold_df$spatial_balance_weight; eff_weight <- gold_df$eff_weight
  # Test homogeneity of allometric slopes (PC ~ log_centsize * taxon, Klingenberg 2016) on active reference anchors
  if (nlevels(gold_df$scientificName) >= 2 && "PC1" %in% names(gold_df)) {
    pc_slope_fit <- tryCatch(stats::lm(PC1 ~ log_centsize * scientificName, data = gold_df, weights = eff_weight), error = function(e) NULL)
    if (!is.null(pc_slope_fit) && is.na(p_int)) { aov_pc <- stats::anova(pc_slope_fit); p_int <- aov_pc[["Pr(>F)"]][grep(":", rownames(aov_pc))[1]] }
  }
  correct_flag <- tolower(trimws(as.character(correct_mode))); slopes_homogeneous <- is.na(p_int) || (p_int > 0.05)
  message(sprintf("Klingenberg (2016) Homogeneity of Allometric Slopes (PC ~ log_centsize * taxon): R2_int = %.4f, p_int = %.4f (homogeneous = %s)", r2_int, p_int, slopes_homogeneous))
  should_correct <- (correct_flag %in% c("true", "1", "yes")) || (correct_flag == "auto" && !is.na(r2_size) && r2_size >= r2_threshold && slopes_homogeneous)
  fit_within <- if (nlevels(gold_df$scientificName) >= 2) stats::lm(Y ~ log_centsize + scientificName, data = gold_df, weights = eff_weight) else stats::lm(Y ~ log_centsize, data = gold_df, weights = eff_weight)
  beta_within <- stats::coef(fit_within)["log_centsize", ]; beta_within[is.na(beta_within)] <- 0.0
  log_centsize <- gold_df$log_centsize
  Y_adj <- Y - outer(log_centsize - mean(log_centsize, na.rm = TRUE), beta_within)
  ref_mean_cs <- stats::weighted.mean(gold_df$log_centsize, w = eff_weight, na.rm = TRUE)
  Y_ref_adj <- Y - outer(gold_df$log_centsize - ref_mean_cs, beta_within)
  active_cols <- which(apply(Y_ref_adj, 2, stats::var, na.rm = TRUE) > 1e-8); if (length(active_cols) == 0L) active_cols <- seq_len(ncol(Y_ref_adj))
  wt_res_cov <- stats::cov.wt(Y_ref_adj[, active_cols, drop = FALSE], wt = eff_weight, center = TRUE); eig_res <- eigen(wt_res_cov$cov, symmetric = TRUE)
  colnames(eig_res$vectors) <- paste0("PC", seq_len(ncol(eig_res$vectors))); rownames(eig_res$vectors) <- sym_names[active_cols]
  Y_all <- as.matrix(sub_df[, sym_names, drop = FALSE]); log_centsize <- sub_df$log_centsize; Y_adj <- Y_all - outer(log_centsize - ref_mean_cs, beta_within)
  res_scores <- sweep(Y_adj[, active_cols, drop = FALSE], 2, wt_res_cov$center, "-") %*% eig_res$vectors
  pca_resids <- structure(list(sdev = sqrt(pmax(eig_res$values, 0)), rotation = eig_res$vectors, center = wt_res_cov$center, scale = FALSE, x = res_scores), class = "prcomp")
  sub_df$PC1_raw <- sub_df$PC1; sub_df$PC2_raw <- sub_df$PC2; sub_df$PC1_corrected <- round(pca_resids$x[, 1], 6); sub_df$PC2_corrected <- round(pca_resids$x[, 2], 6)
  if (should_correct) {
    message(sprintf("[INFO] Foliar allometry detected (R2 = %.2f, p = %.2f). Analyzing allometry-free shape residuals.", r2_size, p_size))
    morpho_matrix_for_gmm <- Y_adj; morpho_df[valid_idx, sym_names] <- round(Y_adj, 6)
    all_sym_idx <- which(complete.cases(morpho_df[, sym_names]))
    proj_pcs <- scale(as.matrix(morpho_df[all_sym_idx, sym_names[active_cols], drop = FALSE]), center = pca_resids$center, scale = FALSE) %*% pca_resids$rotation
    for (p in seq_len(min(5, ncol(proj_pcs)))) morpho_df[all_sym_idx, paste0("PC", p)] <- round(proj_pcs[, p], 6)
  } else {
    message(sprintf("[INFO] Foliar allometry negligible (R2 = %.2f). Analyzing raw symmetric harmonics.", r2_size))
    morpho_matrix_for_gmm <- Y_all
  }
  if (!is.null(output_plot) && nzchar(output_plot)) generate_allometry_diagnostic_plot(sub_df, r2_size, p_size, r2_int, output_plot, p_int = p_int)
  return(list(morpho_matrix = morpho_matrix_for_gmm, corrected = should_correct, beta_within = beta_within,
              r2_size = r2_size, p_size = p_size, r2_interaction = r2_int, p_interaction = p_int, anova_table = anova_df, pairwise_table = pw_df, df = morpho_df))
}
evaluate_ontogenetic_allometry <- test_and_correct_allometry

# Subspace Selection: Horn's Parallel Analysis / Broken-Stick & Williams & Titus (1988) Rule (P <= floor(N_min / 3))
select_morphometric_subspace <- function(X, group_vec = NULL, min_p = 2L) {
  M <- ncol(X); N <- nrow(X); if (M <= min_p) return(colnames(X))
  pc_idx <- grep("^PC[0-9]+$", colnames(X))
  if (length(pc_idx) >= 3L) {
    v_obs <- apply(X[, pc_idx, drop = FALSE], 2, stats::var, na.rm = TRUE); v_prop <- v_obs / max(sum(v_obs), 1e-12)
    k_m <- length(pc_idx); bs_prop <- sapply(seq_len(k_m), function(j) sum(1 / (j:k_m)) / k_m)
    horn_thr <- (1 + sqrt(k_m / max(N, k_m + 1)))^2 / k_m
    keep_pc <- which(v_prop >= pmin(bs_prop * 0.55, horn_thr * 0.55))
    if (length(keep_pc) < min_p) keep_pc <- seq_len(min(min_p, k_m))
    X <- X[, c(pc_idx[keep_pc], setdiff(seq_len(M), pc_idx)), drop = FALSE]
  }
  N_min <- if (!is.null(group_vec) && length(unique(stats::na.omit(group_vec))) >= 2) min(as.integer(table(droplevels(factor(group_vec))))) else N
  P_max <- if (N_min >= 6L) max(min_p, floor(N_min / 3)) else ncol(X) # Williams & Titus (1988): P <= floor(N_min / 3)
  sel <- colnames(X)[seq_len(min(ncol(X), P_max))]
  message(sprintf("Subspace Selection (Horn Parallel Analysis / Broken-Stick & Williams & Titus 1988 P <= floor(N_min / 3)): retained P = %d / %d features (N_min = %d).", length(sel), M, N_min))
  sel
}

# 3. Gaussian Mixture Modeling & Bayes Factor Testing (Seed = 42)
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150, X_eval = NULL) {
  if (is.null(X_eval)) X_eval <- X
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL
  bic_table <- data.frame(K = seq_len(max_k), LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)
  for (k in seq_len(max_k)) {
    set.seed(42 + k)
    if (k == 1) {
      mu <- matrix(colMeans(X), nrow = 1); sigma <- list(cov(X) + diag(1e-6, p)); pi_k <- c(1.0)
    } else {
      km <- stats::kmeans(X, centers = k, nstart = 10, iter.max = 30); mu <- km$centers
      pi_k <- as.vector(table(factor(km$cluster, levels = seq_len(k)))) / N
      sigma <- lapply(seq_len(k), function(c) { sub_x <- X[km$cluster == c, , drop = FALSE]; if (nrow(sub_x) > p) cov(sub_x) + diag(1e-5, p) else diag(1e-2, p) })
    }
    loglik_old <- -Inf; resp <- matrix(0, nrow = N, ncol = k)
    for (iter in seq_len(max_iter)) {
      log_dens <- matrix(0, nrow = N, ncol = k)
      for (c in seq_len(k)) {
        diff_x <- t(t(X) - mu[c, ]); sig_inv <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        quad <- rowSums((diff_x %*% sig_inv) * diff_x)
        log_dens[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + quad)
      }
      max_l <- apply(log_dens, 1, max); log_sum_exp <- max_l + log(rowSums(exp(log_dens - max_l)))
      resp <- exp(log_dens - log_sum_exp); loglik <- sum(log_sum_exp)
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik; N_k <- colSums(resp); pi_k <- N_k / N
      for (c in seq_len(k)) if (N_k[c] > 1e-6) {
        mu[c, ] <- colSums(resp[, c] * X) / N_k[c]; diff_x <- t(t(X) - mu[c, ])
        sigma[[c]] <- (t(diff_x) %*% (resp[, c] * diff_x)) / N_k[c] + diag(1e-5, p)
      }
    }
    n_params <- (k - 1) + k * p + k * (p * (p + 1) / 2); bic_val <- 2 * loglik - n_params * log(N)
    bic_table$LogLik[k] <- round(loglik, 2); bic_table$n_params[k] <- n_params; bic_table$BIC[k] <- round(bic_val, 2)
    if (bic_val > best_bic) {
      best_bic <- bic_val
      ld_ev <- matrix(0, nrow = nrow(X_eval), ncol = k)
      for (c in seq_len(k)) {
        dx_e <- t(t(X_eval) - mu[c, ]); si_e <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        ld_ev[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + rowSums((dx_e %*% si_e) * dx_e))
      }
      ml_e <- apply(ld_ev, 1, max); resp_ev <- exp(ld_ev - (ml_e + log(rowSums(exp(ld_ev - ml_e)))))
      best_model <- list(k = k, classification = apply(resp_ev, 1, which.max), uncertainty = 1 - apply(resp_ev, 1, max), bic = bic_val)
    }
  }
  return(list(best = best_model, bic_table = bic_table))
}

run_gmm_cluster_analysis <- function(df, feature_cols, max_k = 8, use_clustvarsel = FALSE, intraspecific_mode = FALSE) {
  pc_guard_cols <- intersect(paste0("PC", seq_len(5)), feature_cols)
  req_cols <- if (length(pc_guard_cols) > 0) pc_guard_cols else feature_cols
  valid_rows <- which(complete.cases(df[, req_cols, drop = FALSE])); N_valid <- length(valid_rows)
  if (N_valid < 3) {
    return(list(classification = rep(NA_integer_, nrow(df)), uncertainty = rep(NA_real_, nrow(df)), best_k = 1, model_name = "None",
                selected_features = feature_cols, bic_table = data.frame(K = seq_len(max_k), BIC = NA_real_, Two_Delta_BIC = NA_real_, Evidence = "Inconclusive")))
  }
  sub_v <- df[valid_rows, , drop = FALSE]
  is_tier1_pure <- if ("morph_tier_source" %in% names(sub_v)) (!is.na(sub_v$morph_tier_source) & sub_v$morph_tier_source == "Tier_1_Pure") else rep(TRUE, N_valid)
  anc_gmm <- select_adaptive_reference_anchors(sub_v, feature_matrix = as.matrix(sub_v[, req_cols, drop = FALSE]), taxon_col = "species_standardized")
  ref_cand <- anc_gmm$active_mask & is_tier1_pure
  if (sum(ref_cand, na.rm = TRUE) < 6L) ref_cand <- if (sum(is_tier1_pure) >= 6L) is_tier1_pure else rep(TRUE, N_valid)
  X <- as.matrix(sub_v[, feature_cols, drop = FALSE])
  for (j in seq_len(ncol(X))) if (any(is.na(X[, j]))) X[is.na(X[, j]), j] <- if (any(!is.na(X[ref_cand, j]))) mean(X[ref_cand, j], na.rm = TRUE) else 0.0
  selected_cols <- select_morphometric_subspace(X[ref_cand, , drop = FALSE]); X <- X[, selected_cols, drop = FALSE]
  set.seed(42); best_k <- 1; best_name <- "Full_Covariance_EM"; bic_table <- NULL; cls_sub <- NULL; unc_sub <- NULL
  X_mah <- X[, if (length(intersect(pc_guard_cols, colnames(X))) > 0) intersect(pc_guard_cols, colnames(X)) else colnames(X), drop = FALSE]
  X_mah_ref <- X_mah[ref_cand, , drop = FALSE]
  q_mcd <- max(floor(0.85 * nrow(X_mah_ref)), floor((nrow(X_mah_ref) + ncol(X_mah) + 1L) / 2L))
  rob_cov <- if (nrow(X_mah_ref) >= 30 && requireNamespace("MASS", quietly = TRUE)) tryCatch(MASS::cov.rob(X_mah_ref, quantile.used = min(nrow(X_mah_ref) - 1L, q_mcd)), error = function(e) NULL) else NULL
  m_center <- if (!is.null(rob_cov)) rob_cov$center else colMeans(X_mah_ref); m_cov <- if (!is.null(rob_cov)) rob_cov$cov else stats::cov(X_mah_ref)
  mahal_d2 <- tryCatch(stats::mahalanobis(X_mah, center = m_center, cov = m_cov + diag(1e-6, ncol(X_mah))), error = function(e) rep(0.0, nrow(X)))
  mahal_outliers <- (nrow(X_mah_ref) >= 10) & !is.na(mahal_d2) & (mahal_d2 > stats::qchisq(0.995, df = ncol(X_mah)))
  anchor_mask <- ref_cand & !mahal_outliers
  core_idx <- if (sum(anchor_mask, na.rm = TRUE) >= 6) which(anchor_mask) else which(ref_cand)
  if (requireNamespace("clustvarsel", quietly = TRUE) && (isTRUE(use_clustvarsel) || ncol(X) >= 4)) {
    cvs <- tryCatch({ set.seed(42); clustvarsel::clustvarsel(X[if (length(core_idx) > 350L) core_idx[seq_len(350L)] else core_idx, , drop = FALSE], G = seq_len(min(max_k, 5L)), verbose = FALSE) }, error = function(e) NULL)
    if (!is.null(cvs) && length(cvs$subset) >= 2) { selected_cols <- colnames(X)[cvs$subset]; X <- X[, cvs$subset, drop = FALSE] }
  }
  X_pcs <- X; X_core <- X_pcs[core_idx, , drop = FALSE]; eff_k <- max(2L, min(if (nrow(X_core) > 100 && max_k >= 7L) 6L else max_k, nrow(X_core) - 1L))
  modelNames <- if ((nrow(X_core) / ncol(X_pcs)) < 5.0) c("EII", "VII", "EEI", "VEI", "EVI", "VVI") else NULL
  is_intra_gmm <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2L)
  if (requireNamespace("mclust", quietly = TRUE)) {
    gmm_fit <- tryCatch(mclust::Mclust(X_core, G = seq_len(eff_k), modelNames = modelNames, verbose = FALSE), error = function(e) NULL)
    if (!is.null(gmm_fit)) {
      bic_vals <- apply(gmm_fit$BIC, 1, function(row) if (all(is.na(row))) NA_real_ else max(row, na.rm = TRUE))
      bic_table <- data.frame(K = seq_len(eff_k), BIC = round(bic_vals[seq_len(eff_k)], 2))
      use_fit <- if (is_intra_gmm && gmm_fit$G < 2L) tryCatch(mclust::Mclust(X_core, G = 2:eff_k, modelNames = modelNames, verbose = FALSE), error = function(e) NULL) else gmm_fit
      if (!is.null(use_fit)) {
        best_k <- use_fit$G; best_name <- use_fit$modelName; gmm_pred <- predict(use_fit, newdata = X_pcs)
        cls_sub <- as.integer(gmm_pred$classification)
        unc_sub <- if (!is.null(gmm_pred$uncertainty)) gmm_pred$uncertainty else 1.0 - apply(as.matrix(gmm_pred$z), 1, max)
      }
    }
  }
  if (is.null(bic_table)) {
    pure_res <- fit_gmm_em_pure(X_core, max_k = eff_k, X_eval = X_pcs); bic_table <- pure_res$bic_table
    if (!is_intra_gmm || pure_res$best$k >= 2L) { best_k <- pure_res$best$k; cls_sub <- pure_res$best$classification; unc_sub <- pure_res$best$uncertainty }
  }
  if (is.null(cls_sub) || (is_intra_gmm && length(unique(stats::na.omit(cls_sub))) < 2L)) {
    set.seed(42); km <- stats::kmeans(X_core, centers = 2L, nstart = 20L); best_k <- 2L; best_name <- "kmeans_K2_fallback"
    d1 <- sqrt(rowSums(sweep(X_pcs, 2, km$centers[1, ], "-")^2)); d2 <- sqrt(rowSums(sweep(X_pcs, 2, km$centers[2, ], "-")^2))
    cls_sub <- ifelse(d1 <= d2, 1L, 2L); unc_sub <- pmin(d1, d2) / pmax(d1 + d2, 1e-9)
  }
  cls <- rep(NA_integer_, nrow(df)); unc <- rep(NA_real_, nrow(df)); cls[valid_rows] <- cls_sub; unc[valid_rows] <- unc_sub
  bic_table$Delta_BIC_prev <- c(0, diff(bic_table$BIC)); bic_table$Two_Delta_BIC <- round(bic_table$BIC - bic_table$BIC[1], 2)
  bic_table$BayesFactor_vs_Null <- bic_table$Two_Delta_BIC
  bic_table$Evidence <- sapply(seq_len(nrow(bic_table)), function(idx_b) { # Kass & Raftery (1995) 2 * Delta_BIC vs G = 1
    if (idx_b == 1L) return("Null Model (G = 1)"); val <- bic_table$Two_Delta_BIC[idx_b]
    if (is.na(val) || val < 0) return("Negative (favors G = 1)")
    if (val < 2) return("Weak / Inconclusive"); if (val < 6) return("Positive Evidence"); if (val < 10) return("Strong Evidence")
    return("Decisive Evidence (Species Boundary)")
  })
  return(list(classification = cls, uncertainty = unc, best_k = best_k, model_name = best_name, selected_features = selected_cols, bic_table = bic_table))
}

# 4. Canonical Discriminant Analysis with Passive Sample Projection (Track 1: Unsupervised Morphotype_K1..Kn or Multi-Species)
run_canonical_discriminant_analysis <- function(df, feature_cols, target_taxa, gmm_res = NULL, intraspecific_mode = FALSE) {
  message("Configuring Canonical Discriminant Analysis in MorphoTools2 Architecture (Tier_3_Bronze weight = 0.0)...")
  df$authority_weight <- compute_authority_precision_weights(df); df$spatial_balance_weight <- compute_spatial_balance_weights(df)
  authority_weight <- df$authority_weight; spatial_balance_weight <- df$spatial_balance_weight; eff_weight <- authority_weight
  is_tier1_pure <- if ("morph_tier_source" %in% names(df)) (!is.na(df$morph_tier_source) & df$morph_tier_source == "Tier_1_Pure") else rep(TRUE, nrow(df))
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  pc_req <- intersect(paste0("PC", seq_len(5)), feature_cols)
  has_pcs <- complete.cases(df[, if (length(pc_req) > 0) pc_req else feature_cols, drop = FALSE])
  if (is_intra) {
    cl_u <- if (!is.null(gmm_res) && length(unique(stats::na.omit(gmm_res$classification))) >= 2L) as.integer(gmm_res$classification) else NULL
    if (is.null(cl_u)) {
      set.seed(42); km_fb <- stats::kmeans(as.matrix(df[has_pcs, if (length(pc_req) > 0) pc_req else feature_cols, drop = FALSE]), centers = 2L, nstart = 20L)
      cl_u <- rep(NA_integer_, nrow(df)); cl_u[has_pcs] <- as.integer(km_fb$cluster)
    }
    df$grouping_col <- ifelse(is.na(cl_u), NA_character_, paste0("Morphotype_K", cl_u))
    unc_v <- if (!is.null(gmm_res) && !is.null(gmm_res$uncertainty)) gmm_res$uncertainty else rep(0.0, nrow(df))
    u_cut <- stats::quantile(unc_v, 0.70, na.rm = TRUE)
    is_active <- (df$authority_weight > 0) & is_tier1_pure & !is.na(df$grouping_col) & !is.na(unc_v) & (unc_v <= u_cut) & has_pcs
    if (sum(table(df$grouping_col[is_active]) >= 2L) < 2L) is_active <- (df$authority_weight > 0) & is_tier1_pure & !is.na(df$grouping_col) & has_pcs
    if (sum(table(df$grouping_col[is_active]) >= 2L) < 2L) is_active <- is_tier1_pure & !is.na(df$grouping_col) & has_pcs
    if (sum(table(df$grouping_col[is_active]) >= 2L) < 2L) is_active <- !is.na(df$grouping_col) & has_pcs
    act_tab <- table(df$grouping_col[is_active]); target_taxa <- sort(names(act_tab[act_tab >= 2L]))
    if (length(target_taxa) < 2L) target_taxa <- sort(unique(stats::na.omit(df$grouping_col[is_active])))
    is_active <- is_active & (df$grouping_col %in% target_taxa); is_train_cand <- is_active
  } else {
    df$grouping_col <- df$species_standardized
    anc_cda <- select_adaptive_reference_anchors(df, feature_matrix = as.matrix(df[, if (length(pc_req) > 0) pc_req else feature_cols, drop = FALSE]), taxon_col = "grouping_col")
    df$is_promoted_bronze <- anc_cda$is_promoted_bronze; df$anchor_tier_source <- anc_cda$anchor_tier_source
    df$authority_weight <- compute_authority_precision_weights(df); authority_weight <- df$authority_weight; eff_weight <- authority_weight
    is_train_cand <- anc_cda$active_mask & is_tier1_pure & (df$grouping_col %in% target_taxa) & has_pcs
    tr_tab <- table(factor(df$grouping_col[is_train_cand], levels = target_taxa))
    if (sum(tr_tab >= 2L) >= 2L) { target_taxa <- names(tr_tab[tr_tab >= 2L]); is_train_cand <- is_train_cand & (df$grouping_col %in% target_taxa)
    } else if (any(tr_tab < 2L)) {
      if (length(unique(df$grouping_col[is_train_cand])) >= 2L) target_taxa <- intersect(target_taxa, unique(df$grouping_col[is_train_cand]))
      fallback_cand <- anc_cda$active_mask & (df$species_standardized %in% target_taxa) & has_pcs
      is_train_cand <- if (all(table(factor(df$grouping_col[fallback_cand], levels = target_taxa)) >= 2L)) fallback_cand else (df$grouping_col %in% target_taxa & has_pcs)
    }
    is_active <- is_train_cand
  }
  X_full <- as.matrix(df[, feature_cols, drop = FALSE])
  for (j in seq_len(ncol(X_full))) if (any(is.na(X_full[, j]))) X_full[is.na(X_full[, j]), j] <- if (any(!is.na(X_full[is_train_cand, j]))) mean(X_full[is_train_cand, j], na.rm = TRUE) else 0.0
  feature_cols <- select_morphometric_subspace(X_full[is_train_cand, , drop = FALSE], group_vec = df$grouping_col[is_train_cand])
  X <- X_full[, feature_cols, drop = FALSE]; p <- ncol(X)
  for (grp in target_taxa) {
    cand_k <- which(is_train_cand & df$grouping_col == grp); n_k <- length(cand_k)
    if (n_k >= 5) {
      sub_x <- X[cand_k, , drop = FALSE]
      d2_k <- tryCatch(stats::mahalanobis(sub_x, colMeans(sub_x), stats::cov(sub_x) + diag(1e-4, p)), error = function(e) rep(0, n_k))
      drop_k <- !is.na(d2_k) & (d2_k > stats::qchisq(0.95, df = p))
      if (sum(!drop_k) >= 2) is_train_cand[cand_k[drop_k]] <- FALSE
    }
  }
  if (!is_intra) is_active <- is_train_cand
  active_idx <- which(is_train_cand); passive_idx <- which(!is_active)
  message(sprintf("Active Anchors (Complete): %d | Passive Projected: %d", length(active_idx), length(passive_idx)))
  if (length(active_idx) < length(target_taxa)) stop("Insufficient complete active Tier 1 Gold anchors for CDA.")

  X_a <- X[active_idx, , drop = FALSE]; y_a <- df$grouping_col[active_idx]; eff_weight <- eff_weight[active_idx]; if (all(eff_weight <= 0)) eff_weight <- rep(1.0, length(active_idx)); w_a <- eff_weight; N_a <- nrow(X_a)
  X_mat <- X_a; grp_factor <- droplevels(factor(y_a))
  lda_cv <- if (requireNamespace("MASS", quietly = TRUE) && nlevels(grp_factor) >= 2) {
    tryCatch(MASS::lda(x = X_mat, grouping = grp_factor, prior = rep(1.0 / nlevels(grp_factor), nlevels(grp_factor)), CV = TRUE), error = function(e) NULL)
  } else NULL
  loocv_accuracy <- if (!is.null(lda_cv) && !is.null(lda_cv$class)) round(mean(as.character(lda_cv$class) == as.character(grp_factor), na.rm = TRUE), 4) else NA_real_
  message(sprintf("CDA Leave-One-Out Cross-Validation Accuracy (loocv_accuracy): %.4f", loocv_accuracy))

  g_levels <- target_taxa; g <- length(g_levels); grand_mean <- colSums(X_a * w_a) / sum(w_a)
  B <- matrix(0, p, p); W <- matrix(0, p, p)
  group_means <- matrix(0, g, p, dimnames = list(g_levels, feature_cols)); group_counts <- integer(g)
  for (k in seq_along(g_levels)) {
    grp <- g_levels[k]; k_mask <- (y_a == grp); sub_x <- X_a[k_mask, , drop = FALSE]; w_k <- eff_weight[k_mask]; group_counts[k] <- nrow(sub_x)
    if (nrow(sub_x) > 0) {
      m_k <- colSums(sub_x * w_k) / sum(w_k); group_means[k, ] <- m_k; diff_m <- m_k - grand_mean
      B <- B + sum(w_k) * (diff_m %*% t(diff_m)); diff_x <- t(t(sub_x) - m_k); W <- W + (t(diff_x) %*% (diff_x * w_k))
    }
  }

  S_reg <- (W / max(sum(w_a) - g, 1)) + diag(1e-7, p); eig <- eigen(solve(S_reg, B))
  real_idx <- which(abs(Im(eig$values)) < 1e-6); real_vals <- Re(eig$values[real_idx]); real_vecs <- Re(eig$vectors[, real_idx, drop = FALSE])
  order_idx <- order(real_vals, decreasing = TRUE); num_axes <- max(1L, min(3L, ncol(X_mat), g - 1L))
  eig_vals <- real_vals[order_idx][seq_len(num_axes)]; eig_vecs <- real_vecs[, order_idx[seq_len(num_axes)], drop = FALSE]
  for (j in seq_len(num_axes)) {
    v <- eig_vecs[, j]; s <- as.numeric(sqrt(t(v) %*% S_reg %*% v))
    if (s > 1e-8) eig_vecs[, j] <- v / s
  }
  var_pct <- round((eig_vals / max(sum(eig_vals), 1e-12)) * 100, 2)
  message("=== CDA Canonical Variates Eigenvalues & Variance ===")
  for (j in seq_len(num_axes)) message(sprintf("  Can%d: Eigenvalue = %7.4f (%5.2f%% variation)", j, eig_vals[j], var_pct[j]))

  Z_active <- t(t(X_a) - grand_mean) %*% eig_vecs; struct_corr <- stats::cor(X_a, Z_active); X_proj <- X
  for (col_idx in seq_len(p)) {
    na_mask <- is.na(X_proj[, col_idx]); if (any(na_mask)) X_proj[na_mask, col_idx] <- grand_mean[col_idx]
  }
  Z_all <- t(t(X_proj) - grand_mean) %*% eig_vecs; canonical_centroids <- t(t(group_means) - grand_mean) %*% eig_vecs
  priors <- rep(1/g, g); post_probs <- matrix(0, nrow(df), g, dimnames = list(NULL, g_levels))
  for (k in seq_along(g_levels)) {
    diff_z <- t(t(Z_all) - canonical_centroids[k, ]); post_probs[, k] <- log(priors[k]) - 0.5 * rowSums(diff_z^2)
  }
  post_probs <- exp(post_probs - apply(post_probs, 1, max, na.rm = TRUE))
  post_probs <- post_probs / pmax(rowSums(post_probs, na.rm = TRUE), 1e-12); pred_idx <- max.col(post_probs, ties.method = "first")
  posterior_prob <- apply(post_probs, 1, max)
  second_posterior_prob <- if (g >= 2) apply(post_probs, 1, function(r) sort(r, decreasing = TRUE)[2]) else rep(0.0, nrow(df))
  posterior_margin <- posterior_prob - second_posterior_prob
  orth_against <- function(vec, basis_mat) {
    v <- as.numeric(vec) - mean(vec, na.rm = TRUE)
    for (b_col in seq_len(ncol(basis_mat))) { b <- as.numeric(basis_mat[, b_col]); v <- v - (sum(v * b, na.rm = TRUE) / max(sum(b^2, na.rm = TRUE), 1e-12)) * b }
    v
  }
  if (ncol(Z_all) == 1L) {
    src1 <- if ("PC1" %in% colnames(X_proj)) X_proj[, "PC1"] else X_proj[, 1]; can2_v <- orth_against(src1, Z_all)
    if (stats::sd(can2_v, na.rm = TRUE) < 1e-8 && ncol(X_proj) >= 2) can2_v <- orth_against(X_proj[, 2], Z_all)
    Z_all <- cbind(Z_all, can2_v)
  }
  if (ncol(Z_all) == 2L) {
    src3 <- if ("PC3" %in% colnames(X_proj)) X_proj[, "PC3"] else if ("PC2" %in% colnames(X_proj)) X_proj[, "PC2"] else X_proj[, min(3L, ncol(X_proj))]; can3_v <- orth_against(src3, Z_all)
    if (stats::sd(can3_v, na.rm = TRUE) < 1e-8 && ncol(X_proj) >= 2) can3_v <- orth_against(X_proj[, 2], Z_all)
    Z_all <- cbind(Z_all, can3_v)
  }
  colnames(Z_all) <- paste0("Can", seq_len(ncol(Z_all)))
  if (ncol(canonical_centroids) < ncol(Z_all)) {
    for (extra_c in (ncol(canonical_centroids) + 1L):ncol(Z_all)) {
      canonical_centroids <- cbind(canonical_centroids, sapply(g_levels, function(grp) mean(Z_all[active_idx[y_a == grp], extra_c], na.rm = TRUE)))
    }
  }
  colnames(canonical_centroids) <- colnames(Z_all); rownames(canonical_centroids) <- g_levels
  return(list(eigenvalues = eig_vals, variance_pct = var_pct, canonical_scores = Z_all, centroids = canonical_centroids,
              predicted_taxon = g_levels[pred_idx], posterior_prob = posterior_prob, second_posterior_prob = second_posterior_prob, post_probs = post_probs,
              posterior_margin = posterior_margin, is_active = is_active, is_promoted_bronze = df$is_promoted_bronze, anchor_tier_source = df$anchor_tier_source,
              authority_weight = df$authority_weight, structure_loadings = struct_corr, loocv_accuracy = loocv_accuracy, grouping_col = df$grouping_col))
}
run_cda_with_passive_projection <- run_canonical_discriminant_analysis

# 5. Outlier Flagging & Herbarium Misidentification Auditing
audit_misidentifications <- function(df, cda_res, gmm_res) {
  message("Auditing specimens for morphospace outliers & label discordances...")
  n <- nrow(df); flags <- logical(n); triage <- character(n); reasons <- character(n)
  grp_vec <- if (!is.null(cda_res$grouping_col)) cda_res$grouping_col else df$species_standardized
  is_intra_audit <- length(unique(stats::na.omit(df$species_standardized))) < 2 || any(grepl("^Morphotype_K|^LatTertile_", grp_vec))
  pc_c <- intersect(paste0("PC", seq_len(5)), names(df))
  mahal_out <- if (is_intra_audit && length(pc_c) >= 2 && nrow(df) >= 10) {
    Xm <- as.matrix(df[, pc_c, drop = FALSE]); d2 <- tryCatch(stats::mahalanobis(Xm, colMeans(Xm, na.rm = TRUE), stats::cov(Xm, use = "complete.obs") + diag(1e-6, ncol(Xm))), error = function(e) rep(0, n))
    !is.na(d2) & (d2 > stats::qchisq(0.975, df = ncol(Xm)))
  } else rep(FALSE, n)
  for (i in seq_len(nrow(df))) {
    raw_sp <- as.character(grp_vec[i]); if (is.na(raw_sp) || !nzchar(raw_sp)) raw_sp <- as.character(df$species_standardized[i])
    pred_sp <- as.character(cda_res$predicted_taxon[i]); if (is.na(pred_sp) || !nzchar(pred_sp)) pred_sp <- raw_sp
    post_p <- cda_res$posterior_prob[i]; is_act <- isTRUE(cda_res$is_active[i])
    if (is_act) {
      if (raw_sp == pred_sp) { flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Verified_Tier1_Gold_Anchor" }
      else if (post_p >= 0.85) { flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- sprintf("Tier_1_Gold_Morphological_Discordance_to_%s", gsub(" ", "_", pred_sp)) }
      else { flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Tier_1_Gold_Borderline_Variant" }
    } else {
      if (raw_sp == pred_sp && post_p >= 0.75) { flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Congruent_Bronze_Determination" }
      else if (raw_sp != pred_sp && post_p >= 0.70) {
        flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- sprintf("%s_Discordant_to_%s", gsub(" ", "_", raw_sp), gsub(" ", "_", pred_sp))
      } else if (post_p < 0.60) { flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Morphological_Intermediate_Ambiguous" }
      else { flags[i] <- (raw_sp != pred_sp); triage[i] <- if (raw_sp != pred_sp) "MEDIUM" else "LOW"
             reasons[i] <- if (raw_sp != pred_sp) "Candidate_Discordance_Moderate_Posterior" else "Minor_Uncertainty" }
    }
    if (is_intra_audit && !flags[i] && isTRUE(mahal_out[i])) { flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- "Intraspecific_Mahalanobis_Morphospace_Outlier" }
  }
  if (!is.null(cda_res$is_promoted_bronze)) df$is_promoted_bronze <- cda_res$is_promoted_bronze
  if (!is.null(cda_res$anchor_tier_source)) df$anchor_tier_source <- cda_res$anchor_tier_source
  df$authority_weight <- if (!is.null(cda_res$authority_weight)) cda_res$authority_weight else compute_authority_precision_weights(df)
  df$gmm_cluster <- gmm_res$classification; df$gmm_uncertainty <- round(gmm_res$uncertainty, 4); df$loocv_accuracy <- cda_res$loocv_accuracy
  df$cda_predicted_taxon <- cda_res$predicted_taxon; df$cda_posterior_prob <- round(cda_res$posterior_prob, 4)
  df$cda_second_posterior_prob <- round(cda_res$second_posterior_prob, 4); df$cda_posterior_margin <- round(cda_res$posterior_margin, 4)
  df$can1 <- round(cda_res$canonical_scores[, 1], 5)
  df$can2 <- if (ncol(cda_res$canonical_scores) >= 2) round(cda_res$canonical_scores[, 2], 5) else 0.0
  df$can3 <- if (ncol(cda_res$canonical_scores) >= 3) round(cda_res$canonical_scores[, 3], 5) else 0.0
  df$is_passive <- !cda_res$is_active; df$morphometric_outlier_flag <- flags; df$misidentification_flag <- flags
  df$triage_priority <- triage; df$discordance_reason <- reasons
  if ("AR_inv" %in% names(df)) df$reproductive_trait_status <- ifelse(is.na(df$AR_inv), "Unexpanded_or_Absent", "Measured")
  message(sprintf("Audit Complete: %d morphospace outliers/discordances flagged (%d HIGH, %d MEDIUM)", sum(flags), sum(triage == "HIGH"), sum(triage == "MEDIUM")))
  return(df)
}

# 5b. Institutional Batch-Effect Audit (Taxon-Adjusted Marginal ANCOVA & Partial R-Squared)
audit_institutional_batch_effects <- function(morpho_df, output_csv = "outputs/reports/institutional_batch_effect_audit.csv", min_obs = 15) {
  if (!("institutionCode" %in% names(morpho_df))) return(NULL)
  valid_df <- morpho_df[!is.na(morpho_df$institutionCode) & trimws(morpho_df$institutionCode) != "" & !is.na(morpho_df$PC1) & !is.na(morpho_df$PC2), ]
  if (nrow(valid_df) < 10) return(NULL)
  inst_counts <- table(valid_df$institutionCode); rep_insts <- names(inst_counts[inst_counts >= min_obs])
  audit_subset <- if (length(rep_insts) >= 2) valid_df[valid_df$institutionCode %in% rep_insts, ] else valid_df
  audit_subset$institutionCode <- droplevels(factor(audit_subset$institutionCode))
  if (!("species_standardized" %in% names(audit_subset))) audit_subset$species_standardized <- if ("scientificName" %in% names(audit_subset)) audit_subset$scientificName else "Unassigned_Packera"
  audit_subset$species_standardized <- droplevels(factor(audit_subset$species_standardized))
  if (nlevels(audit_subset$institutionCode) < 2) return(NULL)
  has_multi_taxa <- nlevels(audit_subset$species_standardized) >= 2
  lm_pc1 <- if (has_multi_taxa) lm(PC1 ~ species_standardized + institutionCode, data = audit_subset) else lm(PC1 ~ institutionCode, data = audit_subset)
  lm_pc2 <- if (has_multi_taxa) lm(PC2 ~ species_standardized + institutionCode, data = audit_subset) else lm(PC2 ~ institutionCode, data = audit_subset)
  s1 <- summary(lm_pc1); s2 <- summary(lm_pc2); aov1 <- anova(lm_pc1); aov2 <- anova(lm_pc2); i1 <- which(rownames(aov1) == "institutionCode")[1]; r1 <- which(rownames(aov1) == "Residuals")[1]
  i2 <- which(rownames(aov2) == "institutionCode")[1]; r2 <- which(rownames(aov2) == "Residuals")[1]; r2_pc1 <- aov1$`Sum Sq`[i1] / max(aov1$`Sum Sq`[i1] + aov1$`Sum Sq`[r1], 1e-12); r2_pc2 <- aov2$`Sum Sq`[i2] / max(aov2$`Sum Sq`[i2] + aov2$`Sum Sq`[r2], 1e-12)
  anova_df <- data.frame(
    Response = c("PC1", "PC1", "PC2", "PC2"), Term = c("institutionCode", "Residuals", "institutionCode", "Residuals"),
    Df = c(aov1$Df[i1], aov1$Df[r1], aov2$Df[i2], aov2$Df[r2]), Sum_Sq = round(c(aov1$`Sum Sq`[i1], aov1$`Sum Sq`[r1], aov2$`Sum Sq`[i2], aov2$`Sum Sq`[r2]), 4),
    Mean_Sq = round(c(aov1$`Mean Sq`[i1], aov1$`Mean Sq`[r1], aov2$`Mean Sq`[i2], aov2$`Mean Sq`[r2]), 4),
    F_value = c(round(aov1$`F value`[i1], 4), NA, round(aov2$`F value`[i2], 4), NA), P_value = c(signif(aov1$`Pr(>F)`[i1], 4), NA, signif(aov2$`Pr(>F)`[i2], 4), NA),
    R_squared = round(c(r2_pc1, NA, r2_pc2, NA), 4), Adj_R_squared = round(c(s1$adj.r.squared, NA, s2$adj.r.squared, NA), 4),
    Batch_Effect_Status = c(ifelse(r2_pc1 < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA, ifelse(r2_pc2 < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA), stringsAsFactors = FALSE)
  if (!is.null(output_csv) && nzchar(output_csv)) {
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE); tmp_ba <- paste0(output_csv, ".tmp"); write.csv(anova_df, file = tmp_ba, row.names = FALSE); file.rename(tmp_ba, output_csv)
  }
  return(invisible(list(lm_pc1 = lm_pc1, lm_pc2 = lm_pc2, anova_table = anova_df)))
}

# 6. Publication-Quality 2D CDA & 2D GMM Biplot Vector Graphic (PDF)
generate_cda_biplot_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL, manova_res = NULL) {
  if (!exists("render_cda_biplot_4panel_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
  if (exists("render_cda_biplot_4panel_pdf", mode = "function")) render_cda_biplot_4panel_pdf(audit_df, cda_res, gmm_res, output_pdf, feature_cols, manova_res = manova_res)
}

# 7. Main Workflow Orchestrator
run_gmm_morphotools_pipeline <- function(opts) {
  if (!file.exists(opts$harmonics)) stop("Harmonics file not found: ", opts$harmonics)
  efa_df <- read.csv(opts$harmonics, stringsAsFactors = FALSE)
  vouchers_df <- if (file.exists(opts$vouchers)) calculate_circular_phenology(read.csv(opts$vouchers, stringsAsFactors = FALSE)) else NULL
  pca_cols <- paste0("PC", seq_len(opts$num_pcs)); missing_pca <- setdiff(pca_cols, names(efa_df))
  if (length(missing_pca) > 0) stop("Missing PCA columns in EFA table: ", paste(missing_pca, collapse = ", "))
  closed_mask <- grepl("^tier_?[12]", efa_df$assigned_tier, ignore.case = TRUE)
  efa_df$assigned_tier <- ifelse(grepl("^tier_?1", efa_df$assigned_tier, ignore.case = TRUE), "Tier_1_Direct", ifelse(grepl("^tier_?2", efa_df$assigned_tier, ignore.case = TRUE), "Tier_2_Reflected", efa_df$assigned_tier))
  closed_df <- efa_df[closed_mask, ]; valid_idx <- which(complete.cases(closed_df[, pca_cols])); closed_df <- closed_df[valid_idx, ]
  # Ensure primary GMM and CDA run on specimen-level table (1 row = 1 voucher collection event; group_by(catalogNumber) & plant_individual_id)
  if (!"plant_individual_id" %in% names(closed_df)) closed_df$plant_individual_id <- 0L
  spec_keys <- paste0(closed_df$catalogNumber, "_", closed_df$plant_individual_id)
  if (any(duplicated(spec_keys))) {
    pc_sym_cols <- intersect(c(pca_cols, paste0("PC", seq_len(5)), paste0("Sym_PC", seq_len(5))), names(closed_df))
    pet_cols <- intersect(c("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio"), names(closed_df))
    num_cols <- intersect(c(pc_sym_cols, pet_cols, "plant_height_cm", "scape_to_basal_leaf_ratio", "foliar_variance", "leaf_count", "aspect_ratio", "solidity", "ucs_score", "area_px", "centroid_size", "log_centsize", "harmonic_roughness", "leaf_quality_weight"), names(closed_df))
    other_num <- setdiff(num_cols, c(pc_sym_cols, pet_cols)); has_w <- "leaf_quality_weight" %in% names(closed_df)
    meta_cols <- setdiff(names(closed_df), c(num_cols, "catalogNumber", "plant_individual_id", "shape_id", "leaf_id", "mask_source", "n_tier1_leaves", "n_tier2_leaves", "morph_tier_source"))
    grp_list <- split(closed_df, list(closed_df$catalogNumber, closed_df$plant_individual_id), drop = TRUE)
    closed_df <- do.call(rbind, lapply(grp_list, function(grp) {
      n_t1 <- sum(grp$assigned_tier == "Tier_1_Direct", na.rm = TRUE); n_t2 <- sum(grp$assigned_tier == "Tier_2_Reflected", na.rm = TRUE)
      if (n_t1 >= 1L) { grp <- grp[grp$assigned_tier == "Tier_1_Direct", , drop = FALSE]; m_src <- "Tier_1_Pure" } else { m_src <- "Tier_2_Fallback" }
      leaf_quality_weight <- if (has_w) suppressWarnings(as.numeric(grp$leaf_quality_weight)) else rep(1.0, nrow(grp))
      leaf_quality_weight[is.na(leaf_quality_weight) | leaf_quality_weight <= 0] <- 1.0
      is_contour_outlier <- if ("is_contour_outlier" %in% names(grp)) as.logical(grp$is_contour_outlier) else rep(FALSE, nrow(grp)); is_contour_outlier[is.na(is_contour_outlier)] <- FALSE
      if (nrow(grp) >= 3) {
        pc_mat <- as.matrix(grp[, intersect(paste0("PC", seq_len(5)), names(grp)), drop = FALSE]); med_pc <- apply(pc_mat, 2, stats::median, na.rm = TRUE)
        dist <- sqrt(rowSums(sweep(pc_mat, 2, med_pc, "-")^2, na.rm = TRUE)); iqr_d <- stats::IQR(dist, na.rm = TRUE); med_d <- stats::median(dist, na.rm = TRUE)
        drop_leaf <- ((dist > 2.0 * iqr_d & iqr_d > 0) | (dist > 2.5 * med_d)) & (is_contour_outlier | nrow(grp) >= 3)
        if (sum(!drop_leaf) < 2) drop_leaf[utils::head(order(leaf_quality_weight, -dist, decreasing = TRUE), 2)] <- FALSE
        grp <- grp[!drop_leaf, , drop = FALSE]
      }
      out <- grp[1, , drop = FALSE]
      for (col in c(pc_sym_cols, other_num)) out[[col]] <- stats::median(grp[[col]], na.rm = TRUE)
      for (col in pet_cols) out[[col]] <- stats::quantile(suppressWarnings(as.numeric(grp[[col]])), probs = 0.75, na.rm = TRUE, names = FALSE)
      if ("leaf_count" %in% names(out)) out$leaf_count <- nrow(grp)
      out$n_tier1_leaves <- n_t1; out$n_tier2_leaves <- n_t2; out$morph_tier_source <- m_src
      for (col in meta_cols) { v <- stats::na.omit(grp[[col]]); out[[col]] <- if (length(v) > 0) v[1] else grp[[col]][1] }
      out
    })); rownames(closed_df) <- NULL
  } else {
    lc <- if ("leaf_count" %in% names(closed_df) && !any(is.na(closed_df$leaf_count))) as.integer(closed_df$leaf_count) else rep(1L, nrow(closed_df))
    if (!"n_tier1_leaves" %in% names(closed_df)) closed_df$n_tier1_leaves <- ifelse(closed_df$assigned_tier == "Tier_1_Direct", lc, 0L)
    if (!"n_tier2_leaves" %in% names(closed_df)) closed_df$n_tier2_leaves <- ifelse(closed_df$assigned_tier == "Tier_2_Reflected", lc, 0L)
    if (!"morph_tier_source" %in% names(closed_df)) closed_df$morph_tier_source <- ifelse(closed_df$n_tier1_leaves >= 1L | closed_df$assigned_tier == "Tier_1_Direct", "Tier_1_Pure", "Tier_2_Fallback")
  }
  raw_sp_col <- if ("species_raw" %in% names(closed_df) && any(!is.na(closed_df$species_raw))) closed_df$species_raw else if ("scientificName" %in% names(closed_df)) closed_df$scientificName else closed_df$species_standardized; closed_df$species_standardized <- standardize_packera_taxon(raw_sp_col); closed_df$authority_weight <- compute_authority_precision_weights(closed_df)
  if (!is.null(vouchers_df)) {
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), , drop = FALSE]; v_idx <- match(as.character(closed_df$catalogNumber), as.character(v_sub$catalogNumber))
    for (col in intersect(names(v_sub), c("institutionCode", "county", "stateProvince", "decimalLatitude", "decimalLongitude", "latitude", "longitude", "doy", "pheno_sin", "pheno_cos", "regional_group", "spatial_balance_weight"))) if (!col %in% names(closed_df)) closed_df[[col]] <- v_sub[[col]][v_idx]
  }
  if ("regional_group" %in% names(closed_df) && exists("standardize_epa_l3_label_r", mode = "function")) closed_df$regional_group <- standardize_epa_l3_label_r(closed_df$regional_group)
  closed_df$spatial_balance_weight <- compute_spatial_balance_weights(closed_df)
  allom_res <- test_and_correct_allometry(closed_df, correct_mode = opts$correct_allometry, permutations = opts$permutations, seed = opts$seed, output_report = opts$allometry_report, output_plot = opts$allometry_plot, output_pairwise = opts$allometry_pairwise)
  if (!is.null(allom_res$df)) closed_df <- allom_res$df
  include_repro <- FALSE; feature_cols <- pca_cols
  if (!isTRUE(opts$outline_only)) {
    if ("capitulum_aspect_ratio" %in% names(closed_df)) closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$capitulum_aspect_ratio)) else if ("AR_inv" %in% names(closed_df)) closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$AR_inv))
    if (!is.null(opts$reproductive) && file.exists(opts$reproductive)) {
      repro_df <- read.csv(opts$reproductive, stringsAsFactors = FALSE)
      if ("capitulum_aspect_ratio" %in% names(repro_df) && !"AR_inv" %in% names(repro_df)) repro_df$AR_inv <- suppressWarnings(as.numeric(repro_df$capitulum_aspect_ratio))
      r_cols <- c("catalogNumber", setdiff(intersect(names(repro_df), c("AR_inv", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm", "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem", "petiole_to_blade_ratio")), names(closed_df)))
      if (length(r_cols) > 1L) closed_df <- merge(closed_df, repro_df[!duplicated(repro_df$catalogNumber), r_cols, drop = FALSE], by = "catalogNumber", all.x = TRUE)
    }
    man_cand <- file.path(dirname(opts$harmonics), "extracted_leaf_manifest.csv")
    if (!"petiole_to_blade_ratio" %in% names(closed_df) && file.exists(man_cand)) {
      man_df <- tryCatch(read.csv(man_cand, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(man_df) && all(c("catalogNumber", "petiole_to_blade_ratio") %in% names(man_df))) {
        p_cols <- intersect(c("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio"), names(man_df))
        agg_m <- stats::aggregate(man_df[, p_cols, drop = FALSE], by = list(catalogNumber = man_df$catalogNumber), FUN = function(x) stats::quantile(suppressWarnings(as.numeric(x)), probs = 0.75, na.rm = TRUE, names = FALSE))
        closed_df <- merge(closed_df, agg_m[, c("catalogNumber", setdiff(p_cols, names(closed_df))), drop = FALSE], by = "catalogNumber", all.x = TRUE)
      }
    }
    for (tr_col in c("AR_inv", "petiole_to_blade_ratio", "plant_height_cm")) if (tr_col %in% names(closed_df)) {
      closed_df[[tr_col]] <- suppressWarnings(as.numeric(closed_df[[tr_col]])); valid_tr <- which(!is.na(closed_df[[tr_col]]) & is.finite(closed_df[[tr_col]]))
      ref_tr <- intersect(valid_tr, which(closed_df$authority_weight > 0)); if (length(ref_tr) < 3L) ref_tr <- valid_tr
      if (length(valid_tr) >= 5) {
        std_col <- paste0(tr_col, "_std"); closed_df[[std_col]] <- NA_real_; m_tr <- mean(closed_df[[tr_col]][ref_tr]); s_tr <- stats::sd(closed_df[[tr_col]][ref_tr])
        closed_df[[std_col]][valid_tr] <- if (is.finite(s_tr) && s_tr > 1e-6) round((closed_df[[tr_col]][valid_tr] - m_tr) / s_tr, 5) else 0.0
        include_repro <- TRUE; feature_cols <- unique(c(feature_cols, std_col))
      }
    }
  }
  cand_taxa <- sort(unique(closed_df$species_standardized[!is.na(closed_df$species_standardized) & grepl("^Packera\\s+", closed_df$species_standardized)]))
  tx_n <- table(factor(closed_df$species_standardized, levels = cand_taxa)); target_taxa <- if (sum(tx_n >= 2L) >= 2L) names(tx_n[tx_n >= 2L]) else if (length(cand_taxa) > 0L) cand_taxa else sort(unique(stats::na.omit(closed_df$species_standardized)))
  is_intra_mode <- isTRUE(opts$intraspecific_mode) || (length(target_taxa) < 2L)
  gmm_res <- run_gmm_cluster_analysis(closed_df, feature_cols, max_k = opts$max_k, use_clustvarsel = isTRUE(opts$use_clustvarsel), intraspecific_mode = is_intra_mode)
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE); tmp_bic <- paste0(opts$output_report, ".tmp"); write.csv(gmm_res$bic_table, file = tmp_bic, row.names = FALSE); file.rename(tmp_bic, opts$output_report)
  veg_extra <- intersect(c("log_centsize", "aspect_ratio", "solidity"), names(closed_df)); ref_std_idx <- which(closed_df$authority_weight > 0); if (length(ref_std_idx) < 3L) ref_std_idx <- seq_len(nrow(closed_df))
  for (vc in veg_extra) if (all(is.finite(closed_df[[vc]])) && stats::sd(closed_df[[vc]][ref_std_idx]) > 1e-6) closed_df[[paste0(vc, "_std")]] <- (closed_df[[vc]] - mean(closed_df[[vc]][ref_std_idx])) / stats::sd(closed_df[[vc]][ref_std_idx])
  cda_cols <- c(feature_cols, if (!isTRUE(opts$outline_only)) intersect(paste0(veg_extra, "_std"), names(closed_df)) else character(0))
  cda_unsupervised <- run_cda_with_passive_projection(closed_df, cda_cols, target_taxa, gmm_res = gmm_res, intraspecific_mode = is_intra_mode); cda_res <- cda_unsupervised
  audited_df <- audit_misidentifications(closed_df, cda_res, gmm_res); audit_institutional_batch_effects(audited_df, output_csv = opts$output_batch_audit)
  rep_dir <- dirname(opts$output_report); eco_sum_csv <- file.path(rep_dir, "ecoregion_supervised_cda_summary.csv")
  eco_cda_res <- if (exists("run_supervised_ecoregion_cda", mode = "function")) run_supervised_ecoregion_cda(audited_df, cda_cols, gmm_res = gmm_res, output_summary_csv = eco_sum_csv, permutations = opts$permutations, seed = opts$seed) else NULL
  if (!is.null(eco_cda_res)) {
    audited_df$regional_group <- eco_cda_res$ecoregion_vec; audited_df$ecoregion_cda_predicted <- eco_cda_res$ecoregion_cda_predicted
    audited_df$ecoregion_cda_posterior_prob <- eco_cda_res$ecoregion_cda_posterior_prob; audited_df$ecoregion_can1 <- eco_cda_res$ecoregion_can1; audited_df$ecoregion_can2 <- eco_cda_res$ecoregion_can2
  }
  manova_res <- if (exists("run_multivariate_manova_audit", mode = "function")) run_multivariate_manova_audit(audited_df, cda_cols, gmm_res, output_group_csv = file.path(rep_dir, "manova_group_significance.csv"), output_pairwise_csv = file.path(rep_dir, "manova_pairwise_posthoc.csv")) else NULL
  lead_cols <- intersect(c("catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier", "n_tier1_leaves", "n_tier2_leaves", "morph_tier_source", "pixels_per_mm", "is_calibrated_mm",
                           "species_raw", "species_standardized", "determiner_tier", "authority_weight", "spatial_balance_weight", "cda_predicted_taxon", "cda_posterior_prob", "cda_second_posterior_prob",
                           "cda_posterior_margin", "loocv_accuracy", "is_passive", "morphometric_outlier_flag", "misidentification_flag", "triage_priority", "discordance_reason", "gmm_cluster", "gmm_uncertainty",
                           "can1", "can2", "can3", "ecoregion_cda_predicted", "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2", "pheno_sin", "pheno_cos", "PC1", "PC2", "PC3", "PC4", "PC5",
                           "centroid_size", "log_centsize", "solidity", "ucs_score", "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "petiole_to_blade_ratio_std",
                           "AR_inv", "AR_inv_std", "plant_height_cm", "plant_height_cm_std", "scape_to_basal_leaf_ratio", "is_bent_stem", "reproductive_trait_status", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm"), names(audited_df))
  out_df <- audited_df[, c(lead_cols, setdiff(names(audited_df), lead_cols))]; dir.create(dirname(opts$output_flags), recursive = TRUE, showWarnings = FALSE); tmp_flags <- paste0(opts$output_flags, ".tmp"); write.csv(out_df, file = tmp_flags, row.names = FALSE, na = ""); file.rename(tmp_flags, opts$output_flags)
  generate_cda_biplot_pdf(audited_df, cda_res, gmm_res, opts$output_plot, feature_cols, manova_res = manova_res)
  tryCatch({ if (exists("generate_step04_pub_figures", mode = "function")) generate_step04_pub_figures(audited_df, cda_res, gmm_res, NULL, opts$output_plot, eco_cda_res = eco_cda_res, manova_res = manova_res) }, error = function(e) warning("Step 04 publication figures export warning: ", e$message))
  return(invisible(out_df))
}
if (sys.nframe() == 0) run_gmm_morphotools_pipeline(parse_args_robust())
