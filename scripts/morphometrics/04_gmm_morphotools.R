#!/usr/bin/env Rscript
# Script: 04_gmm_morphotools.R - Packera Morphometrics Pipeline (NCU)
suppressPackageStartupMessages({
  for (pkg in c("MASS", "mclust", "clustvarsel", "MorphoTools2", "RRPP", "vegan", "dplyr", "readr", "tibble", "ggplot2", "optparse")) {
    if (requireNamespace(pkg, quietly = TRUE)) library(pkg, character.only = TRUE)
  }
})

# 1. CLI Argument Parsing
parse_args_robust <- function() {
  option_list <- list(
    optparse::make_option(c("-v", "--vouchers"), type = "character", default = "data/tables/curated_vouchers.csv"), optparse::make_option(c("-e", "--harmonics"), type = "character", default = "data/tables/leaf_efa_harmonics.csv"),
    optparse::make_option(c("-m", "--reproductive"), type = "character", default = "data/tables/voucher_reproductive_metrics.csv"), optparse::make_option(c("--outline-only"), action = "store_true", default = FALSE),
    optparse::make_option(c("--correct-allometry"), type = "character", default = "auto"), optparse::make_option(c("--allometry-report"), type = "character", default = "outputs/reports/allometry_significance_audit.csv"),
    optparse::make_option(c("--allometry-pairwise"), type = "character", default = "outputs/reports/allometry_pairwise_trajectories.csv"), optparse::make_option(c("--allometry-plot"), type = "character", default = "outputs/figures/fig04a_foliar_allometry_audit.png"),
    optparse::make_option(c("--use-clustvarsel"), action = "store_true", default = FALSE), optparse::make_option(c("-f", "--output-flags"), type = "character", default = "data/tables/morphometric_flags.csv"),
    optparse::make_option(c("-p", "--output-plot"), type = "character", default = "outputs/figures/fig04b_gmm_cda_morphometrics.png"), optparse::make_option(c("-r", "--output-report"), type = "character", default = "outputs/reports/gmm_bayes_factors_summary.csv"),
    optparse::make_option(c("-b", "--output-batch-audit"), type = "character", default = "outputs/reports/institutional_batch_effect_audit.csv"), optparse::make_option(c("-k", "--max-k"), type = "integer", default = 8),
    optparse::make_option(c("-d", "--num-pcs"), type = "integer", default = 5), optparse::make_option(c("--permutations"), type = "integer", default = 999),
    optparse::make_option(c("--seed"), type = "integer", default = 42), optparse::make_option(c("--intraspecific-mode"), action = "store_true", default = FALSE), optparse::make_option(c("--focal-taxon"), type = "character", default = "")
  )
  if (requireNamespace("optparse", quietly = TRUE)) { opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list)); names(opts) <- gsub("-", "_", names(opts)); return(opts) }
  defs <- list(vouchers = "data/tables/curated_vouchers.csv", harmonics = "data/tables/leaf_efa_harmonics.csv", reproductive = "data/tables/voucher_reproductive_metrics.csv", outline_only = FALSE, correct_allometry = "auto",
       allometry_report = "outputs/reports/allometry_significance_audit.csv", allometry_pairwise = "outputs/reports/allometry_pairwise_trajectories.csv", allometry_plot = "outputs/figures/fig04a_foliar_allometry_audit.png",
       use_clustvarsel = FALSE, output_flags = "data/tables/morphometric_flags.csv", output_plot = "outputs/figures/fig04b_gmm_cda_morphometrics.png", output_report = "outputs/reports/gmm_bayes_factors_summary.csv",
       output_batch_audit = "outputs/reports/institutional_batch_effect_audit.csv", max_k = 8L, num_pcs = 5L, permutations = 999L, seed = 42L, intraspecific_mode = FALSE, focal_taxon = "")
  if (exists(".parse_step04_cli_fallback", mode = "function")) .parse_step04_cli_fallback(commandArgs(trailingOnly = TRUE), defs) else defs
}

# 2. Taxonomic Concept Standardization & Dynamic Column Discovery

for (.pf_p in c("scripts/morphometrics/pub_figures_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "pub_figures_utils.R"))) {
  if (!exists("compute_spatial_balance_weights", mode = "function") && !is.na(.pf_p) && nzchar(.pf_p) && file.exists(.pf_p)) { source(.pf_p); break }
}
for (.aw_p in c("scripts/morphometrics/authority_weighting_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "authority_weighting_utils.R"))) if (!exists("compute_empirical_bayes_authority_weights", mode = "function") && !is.na(.aw_p) && nzchar(.aw_p) && file.exists(.aw_p)) { source(.aw_p); break }
for (.pm_p in c("scripts/morphometrics/plasticity_mitigation_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "plasticity_mitigation_utils.R"))) if (!exists("merge_plastic_gmm_components_ridgeline", mode = "function") && !is.na(.pm_p) && nzchar(.pm_p) && file.exists(.pm_p)) { source(.pm_p); break }
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!exists("fetch_step04_morphometrics_from_db", mode = "function") && !is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
}
.base_std_taxon <- if (exists("standardize_taxon_name", mode = "function")) standardize_taxon_name else NULL
standardize_packera_taxon <- function(species_vec) {
  if (!is.null(.base_std_taxon)) { res <- .base_std_taxon(species_vec); return(ifelse(is.na(res), "Unknown", res)) }
  as.character(species_vec)
}

discover_sorted_morph_cols <- function(x, pattern, require_finite = TRUE) {
  nms <- if (is.data.frame(x) || is.matrix(x)) colnames(x) else as.character(x); hits <- grep(pattern, nms, value = TRUE)
  if (length(hits) > 1L) hits <- hits[order(suppressWarnings(as.integer(gsub("[^0-9]", "", hits))), hits)]
  if (isTRUE(require_finite) && (is.data.frame(x) || is.matrix(x)) && length(hits) > 0L) hits <- hits[vapply(hits, function(cl) any(is.finite(suppressWarnings(as.numeric(if (is.data.frame(x)) x[[cl]] else x[, cl])))), logical(1))]
  hits
}

# 2b. Foliar Allometry Multivariate Testing & Residual Correction
generate_allometry_diagnostic_plot <- function(sub_df, r2_size, p_size, r2_int, output_pdf, p_int = 1.0) {
  if (!exists("render_allometry_diagnostic_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
  if (exists("render_allometry_diagnostic_pdf", mode = "function")) render_allometry_diagnostic_pdf(sub_df, r2_size, p_size, r2_int, output_pdf, p_int = p_int)
  message("Shape allometry regression plot saved -> ", output_pdf)
}

test_and_correct_allometry <- function(morpho_df, correct_mode = "auto", permutations = 999, seed = 42, r2_threshold = 0.10,
                                       output_report = "outputs/reports/allometry_significance_audit.csv",
                                       output_plot = "outputs/figures/fig04a_foliar_allometry_audit.png",
                                       output_pairwise = "outputs/reports/allometry_pairwise_trajectories.csv") {
  message("=== Foliar Allometry Assessment (Multivariate Size vs. Shape Covariation) ===")
  sym_names <- c(discover_sorted_morph_cols(morpho_df, "^A[0-9]+$"), discover_sorted_morph_cols(morpho_df, "^D[0-9]+$"))
  if (!"log_centsize" %in% names(morpho_df) || length(sym_names) < 4) {
    message("Warning: log_centsize or symmetric harmonics not available. Skipping allometry test.")
    return(list(morpho_matrix = NULL, corrected = FALSE, anova_table = NULL, df = morpho_df))
  }
  cal_mask <- if ("is_calibrated_mm" %in% names(morpho_df)) (!is.na(morpho_df$is_calibrated_mm) & as.logical(morpho_df$is_calibrated_mm)) else rep(TRUE, nrow(morpho_df))
  if ("scale_source" %in% names(morpho_df)) cal_mask <- cal_mask & !(tolower(trimws(as.character(morpho_df$scale_source))) %in% c("mp_estimated", "predicted_from_megapixels", "uncalibrated"))
  valid_idx <- which(cal_mask & !is.na(morpho_df$log_centsize) & is.finite(morpho_df$log_centsize) & complete.cases(morpho_df[, sym_names]))
  if (length(valid_idx) < 10) return(list(morpho_matrix = NULL, corrected = FALSE, anova_table = NULL, df = morpho_df))
  sub_df <- morpho_df[valid_idx, , drop = FALSE]; sub_df$scientificName <- factor(standardize_packera_taxon(as.character(if ("species_standardized" %in% names(sub_df)) sub_df$species_standardized else sub_df$scientificName)))
  is_tier1_pure <- if ("morph_tier_source" %in% names(sub_df)) (is.na(sub_df$morph_tier_source) | sub_df$morph_tier_source != "Cauline_Only_Fallback") else rep(TRUE, nrow(sub_df))
  anc_allom <- select_adaptive_reference_anchors(sub_df, feature_matrix = as.matrix(sub_df[, sym_names[seq_len(min(4L, length(sym_names)))], drop = FALSE]), taxon_col = "scientificName")
  sub_df$is_promoted_bronze <- anc_allom$is_promoted_bronze; sub_df$authority_weight <- compute_authority_precision_weights(sub_df)
  gold_mask <- anc_allom$active_mask & !is.na(sub_df$scientificName) & is_tier1_pure
  sub_df$spatial_balance_weight <- compute_spatial_balance_weights(sub_df, active_mask = gold_mask, recompute = TRUE)
  sub_df$eff_weight <- sub_df$authority_weight * sub_df$spatial_balance_weight
  min_sp <- min(2L, length(unique(stats::na.omit(sub_df$scientificName))))
  gold_df <- if (sum(gold_mask) >= 6 && length(unique(stats::na.omit(sub_df$scientificName[gold_mask]))) >= min_sp) sub_df[gold_mask, , drop = FALSE] else if (sum(is_tier1_pure) >= 6) sub_df[is_tier1_pure, , drop = FALSE] else sub_df
  if (all(gold_df$eff_weight <= 0)) gold_df$eff_weight <- 1.0
  gold_df$scientificName <- droplevels(factor(gold_df$scientificName)); Y <- as.matrix(gold_df[, sym_names, drop = FALSE]); set.seed(seed); anova_df <- NULL; pw_df <- NULL
  if (requireNamespace("RRPP", quietly = TRUE) && nlevels(gold_df$scientificName) >= 2) {
    rr_par <- if (nrow(gold_df) >= 120L && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(16L) else FALSE
    rdf <- RRPP::rrpp.data.frame(Y = Y, log_centsize = gold_df$log_centsize, scientificName = gold_df$scientificName)
    fit_rrpp <- tryCatch(RRPP::lm.rrpp(Y ~ log_centsize * scientificName, data = rdf, weights = gold_df$eff_weight, iter = permutations, turbo = TRUE, Parallel = rr_par, print.progress = FALSE), error = function(e) tryCatch(RRPP::lm.rrpp(Y ~ log_centsize * scientificName, data = rdf, iter = permutations, turbo = TRUE, Parallel = rr_par, print.progress = FALSE), error = function(e2) NULL))
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
  if (!is.null(output_pairwise) && nzchar(output_pairwise) && !is.null(pw_df)) { dir.create(dirname(output_pairwise), recursive = TRUE, showWarnings = FALSE); tmp_pw <- paste0(output_pairwise, ".", Sys.getpid(), ".tmp"); write.csv(pw_df, file = tmp_pw, row.names = FALSE); file.rename(tmp_pw, output_pairwise) }
  anova_df$Term <- gsub("^Residuals$", "Residual", anova_df$Term)
  if ("Total" %in% anova_df$Term) anova_df$R2[anova_df$Term == "Total"] <- 1.0
  anova_df <- anova_df[, intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(anova_df)), drop = FALSE]
  if (!is.null(output_report) && nzchar(output_report)) {
    dir.create(dirname(output_report), recursive = TRUE, showWarnings = FALSE); tmp_rep <- paste0(output_report, ".", Sys.getpid(), ".tmp"); write.csv(anova_df, file = tmp_rep, row.names = FALSE); file.rename(tmp_rep, output_report)
    message("Allometry significance audit report exported -> ", output_report)
  }
  size_row <- anova_df[grepl("^log_centsize$", anova_df$Term, ignore.case = TRUE), ]
  r2_size <- if (nrow(size_row) > 0 && is.finite(size_row$R2[1])) size_row$R2[1] else 0.0; p_size <- if (nrow(size_row) > 0 && !is.na(size_row$p_value[1]) && is.finite(size_row$p_value[1])) size_row$p_value[1] else 1.0
  int_row <- anova_df[grepl("log_centsize:scientificName|scientificName:log_centsize", anova_df$Term, ignore.case = TRUE), ]
  r2_int <- if (nrow(int_row) > 0 && is.finite(int_row$R2[1])) int_row$R2[1] else 0.0; p_int <- if (nrow(int_row) > 0) int_row$p_value[1] else NA_real_
  authority_weight <- gold_df$authority_weight; spatial_balance_weight <- gold_df$spatial_balance_weight; eff_weight <- gold_df$eff_weight
  # Fallback to univariate PC1 ~ log_centsize * scientificName only when multivariate interaction p_int is unavailable
  used_pc1_fb <- FALSE
  if (nlevels(gold_df$scientificName) >= 2 && (is.na(p_int) || !is.finite(p_int)) && "PC1" %in% names(gold_df)) {
    pc_slope_fit <- tryCatch(stats::lm(PC1 ~ log_centsize * scientificName, data = gold_df, weights = eff_weight), error = function(e) NULL)
    if (!is.null(pc_slope_fit)) { aov_pc <- stats::anova(pc_slope_fit); p_int <- aov_pc[["Pr(>F)"]][grep(":", rownames(aov_pc))[1]]; used_pc1_fb <- TRUE }
  }
  if (is.na(p_int) || !is.finite(p_int)) p_int <- 1.0
  correct_flag <- tolower(trimws(as.character(correct_mode))); slopes_homogeneous <- p_int > 0.05
  message(sprintf("Klingenberg (2016) Homogeneity of Allometric Slopes (%s): R2_int = %.4f, p_int = %.4f (homogeneous = %s)", ifelse(used_pc1_fb, "PC1 ~ log_centsize * taxon fallback", "Y ~ log_centsize * taxon"), r2_int, p_int, slopes_homogeneous))
  should_correct <- (correct_flag %in% c("true", "1", "yes")) || (correct_flag == "auto" && !is.na(r2_size) && r2_size >= r2_threshold && !is.na(p_size) && p_size < 0.05 && slopes_homogeneous)
  fit_within <- if (nlevels(gold_df$scientificName) >= 2) stats::lm(Y ~ log_centsize + scientificName, data = gold_df, weights = eff_weight) else stats::lm(Y ~ log_centsize, data = gold_df, weights = eff_weight)
  beta_within <- stats::coef(fit_within)["log_centsize", ]; beta_within[is.na(beta_within)] <- 0.0
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
    message(sprintf("[INFO] Foliar allometry detected (R2 = %.2f, p = %.4f). Analyzing allometry-free shape residuals.", r2_size, p_size))
    morpho_matrix_for_gmm <- Y_adj; morpho_df[valid_idx, sym_names] <- round(Y_adj, 6)
    morpho_df <- .update_step04_allometry_pcs(morpho_df, sub_df, gold_df, valid_idx, sym_names, active_cols, pca_resids, ref_mean_cs, eff_weight)
    morpho_df <- .update_step04_allometry_scalars(morpho_df, sub_df, gold_df, valid_idx, ref_mean_cs)
  } else {
    message(sprintf("[INFO] Foliar allometry negligible (R2 = %.2f). Analyzing raw symmetric harmonics.", r2_size))
    morpho_matrix_for_gmm <- Y_all
  }
  if (!is.null(output_plot) && nzchar(output_plot)) generate_allometry_diagnostic_plot(sub_df, r2_size, p_size, r2_int, output_plot, p_int = p_int)
  return(list(morpho_matrix = morpho_matrix_for_gmm, corrected = should_correct, beta_within = beta_within,
              r2_size = r2_size, p_size = p_size, r2_interaction = r2_int, p_interaction = p_int, anova_table = anova_df, pairwise_table = pw_df, df = morpho_df))
}
# Subspace Selection: Horn's Parallel Analysis / Broken-Stick & Williams & Titus (1988) Rule (P <= floor(N_min / 3))
# 3. Gaussian Mixture Modeling & deviance-scale 2*ln(BF) Testing (fit_gmm_em_pure in pub_synthesis_plates.R)
run_gmm_cluster_analysis <- function(df, feature_cols, max_k = 8, use_clustvarsel = FALSE, intraspecific_mode = FALSE, seed = 42L) {
  feature_cols <- intersect(as.character(feature_cols), names(df))
  pc_guard_cols <- grep("^PC[0-9]+$", feature_cols, value = TRUE)
  req_cols <- if (length(pc_guard_cols) > 0) pc_guard_cols else feature_cols
  valid_rows <- which(complete.cases(df[, req_cols, drop = FALSE])); N_valid <- length(valid_rows)
  if (N_valid < 3) {
    return(list(classification = rep(NA_integer_, nrow(df)), uncertainty = rep(NA_real_, nrow(df)), best_k = 1, model_name = "None",
                selected_features = feature_cols, bic_table = data.frame(K = seq_len(max_k), BIC = NA_real_, Two_Delta_BIC = NA_real_, Evidence = "Inconclusive")))
  }
  sub_v <- df[valid_rows, , drop = FALSE]
  is_tier1_pure <- if ("morph_tier_source" %in% names(sub_v)) (is.na(sub_v$morph_tier_source) | sub_v$morph_tier_source != "Cauline_Only_Fallback") else rep(TRUE, N_valid)
  anc_gmm <- select_adaptive_reference_anchors(sub_v, feature_matrix = as.matrix(sub_v[, req_cols, drop = FALSE]), taxon_col = "species_standardized")
  ref_cand <- anc_gmm$active_mask & !anc_gmm$is_confident_label_noise & is_tier1_pure
  if (sum(ref_cand, na.rm = TRUE) < 6L) ref_cand <- if (sum(is_tier1_pure) >= 6L) is_tier1_pure else rep(TRUE, N_valid)
  keep_comp <- vapply(feature_cols, function(fc) (fc %in% pc_guard_cols) || (mean(is.finite(suppressWarnings(as.numeric(sub_v[[fc]][ref_cand])))) >= 0.70), logical(1))
  if (sum(keep_comp) >= 2L) feature_cols <- feature_cols[keep_comp]
  X <- as.matrix(sub_v[, feature_cols, drop = FALSE]); storage.mode(X) <- "double"
  for (j in seq_len(ncol(X))) if (any(is.na(X[, j]))) X[is.na(X[, j]), j] <- if (any(!is.na(X[ref_cand, j]))) stats::median(X[ref_cand, j], na.rm = TRUE) else 0.0
  if (exists("screen_gmm_point_mass_columns", mode = "function")) X <- X[, screen_gmm_point_mass_columns(X, ref_cand), drop = FALSE]
  selected_cols <- select_morphometric_subspace(X[ref_cand, , drop = FALSE], aux_df = sub_v[ref_cand, , drop = FALSE]); X <- X[, selected_cols, drop = FALSE]
  seed_i <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L); set.seed(seed_i); best_k <- 1; best_name <- "Full_Covariance_EM"; bic_table <- NULL; bic_raw_mat <- NULL; cls_sub <- NULL; unc_sub <- NULL; did_merge <- FALSE
  X_mah <- X[, if (length(intersect(pc_guard_cols, colnames(X))) > 0) intersect(pc_guard_cols, colnames(X)) else colnames(X), drop = FALSE]
  mahal_outliers <- screen_regime_stratified_mcd_outliers(X_mah, fit_pos = which(ref_cand), solidity_vec = sub_v$solidity, is_dissected_vec = sub_v$is_dissected, k_mcd = ncol(X_mah), p_cutoff = 0.005, seed = seed_i)
  anchor_mask <- ref_cand & !mahal_outliers
  core_idx <- if (sum(anchor_mask, na.rm = TRUE) >= 6) which(anchor_mask) else which(ref_cand)
  if (requireNamespace("clustvarsel", quietly = TRUE) && (isTRUE(use_clustvarsel) || ncol(X) >= 4)) {
    cvs <- if (exists("run_bounded_clustvarsel", mode = "function")) run_bounded_clustvarsel(X, core_idx, max_k = max_k, seed = seed_i) else tryCatch({ set.seed(seed_i); clustvarsel::clustvarsel(X[if (length(core_idx) > 350L) core_idx[seq_len(350L)] else core_idx, , drop = FALSE], G = seq_len(min(max_k, 5L)), itermax = 100, parallel = FALSE, verbose = FALSE) }, error = function(e) NULL)
    if (!is.null(cvs) && length(cvs$subset) >= 2) { selected_cols <- colnames(X)[cvs$subset]; X <- X[, cvs$subset, drop = FALSE] }
  }
  is_intra_gmm <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2L)
  X_pcs <- X; X_core <- X_pcs[core_idx, , drop = FALSE]; eff_k <- max(2L, min(if (nrow(X_core) > 100 && max_k >= 7L) 6L else max_k, nrow(X_core) - 1L))
  modelNames <- if ((nrow(X_core) / ncol(X_pcs)) < 5.0) c("EII", "VII", "EEI", "VEI", "EVI", "VVI") else NULL; z_eval_mat <- NULL
  if (requireNamespace("mclust", quietly = TRUE)) {
    set.seed(seed_i); fit_mclust_fn <- function(G_arg) if (exists("fit_parallel_mclust_gmm", mode = "function")) fit_parallel_mclust_gmm(X_core, G_arg, modelNames) else mclust::Mclust(X_core, G = G_arg, modelNames = modelNames, verbose = FALSE)
    gmm_fit <- tryCatch(fit_mclust_fn(seq_len(eff_k)), error = function(e) NULL)
    if (!is.null(gmm_fit)) {
      bic_raw_mat <- as.matrix(gmm_fit$BIC); bic_vals <- apply(bic_raw_mat, 1, function(row) if (all(is.na(row))) NA_real_ else max(row, na.rm = TRUE))
      bic_table <- data.frame(K = seq_len(eff_k), BIC = round(bic_vals[seq_len(eff_k)], 2))
      if (!is.null(colnames(bic_raw_mat))) for (mn in colnames(bic_raw_mat)) bic_table[[paste0("BIC_", mn)]] <- round(as.numeric(bic_raw_mat[seq_len(eff_k), mn]), 2)
      best_k <- gmm_fit$G; best_name <- gmm_fit$modelName; gmm_pred <- predict(gmm_fit, newdata = X_pcs); z_eval_mat <- as.matrix(gmm_pred$z)
      cls_sub <- as.integer(gmm_pred$classification); unc_sub <- if (!is.null(gmm_pred$uncertainty)) gmm_pred$uncertainty else 1.0 - apply(z_eval_mat, 1, max)
    }
  }
  if (is.null(bic_table)) {
    pure_res <- fit_gmm_em_pure(X_core, max_k = eff_k, X_eval = X_pcs, seed = seed_i); bic_table <- pure_res$bic_table
    best_k <- pure_res$best$k; cls_sub <- pure_res$best$classification; unc_sub <- pure_res$best$uncertainty; z_eval_mat <- pure_res$best$z
  }
  if (!is.null(z_eval_mat) && ncol(z_eval_mat) >= 2L && exists("merge_plastic_gmm_components_ridgeline", mode = "function")) {
    mrg <- merge_plastic_gmm_components_ridgeline(X_pcs, X_core, z_eval_mat, z_eval_mat[core_idx, , drop = FALSE], intraspecific_mode = is_intra_gmm)
    if (isTRUE(mrg$merged_any)) { did_merge <- TRUE; message(sprintf("[Plasticity Gate] Morphological ridgeline unimodality test merged GMM K=%d -> K=%d.", best_k, mrg$best_k)); best_k <- mrg$best_k; cls_sub <- mrg$classification; unc_sub <- mrg$uncertainty; best_name <- paste0(best_name, "_RidgelineMerged") }
  }
  if (is.null(cls_sub)) { best_k <- 1L; best_name <- "Unimodal_G1_Null"; cls_sub <- rep(1L, nrow(X_pcs)); unc_sub <- rep(0.0, nrow(X_pcs)) }
  relab <- .canonical_relabel_gmm_clusters(sub_v, cls_sub, best_k, ref_cand, core_idx); cls_sub <- relab$cls_sub; cls <- rep(NA_integer_, nrow(df)); unc <- rep(NA_real_, nrow(df)); cls[valid_rows] <- cls_sub; unc[valid_rows] <- unc_sub
  n_core_g <- max(2.0, nrow(X_core)); n_eff_g <- max(2.0, n_core_g * min(1.0, max(0.10, relab$n_eff_ratio))); p_g <- ncol(X_pcs)
  k_par_vec <- if (exists(".compute_exact_mclust_params", mode = "function")) .compute_exact_mclust_params(bic_table$K, p_g, bic_mat = bic_raw_mat, default_model = best_name, explicit_n_params = if ("n_params" %in% names(bic_table)) bic_table$n_params else NULL) else (bic_table$K * p_g + (bic_table$K - 1L) + bic_table$K * p_g)
  bic_neff <- bic_table$BIC + k_par_vec * (log(n_core_g) - log(n_eff_g))
  # Note: mclust BIC = 2*logLik - k*log(n) is already on the deviance 2*ln(BF) scale (Kass & Raftery 1995); BayesFactor_vs_Null stores Kish-N_eff-adjusted 2*ln(BF_{K,1}).
  bic_table$Delta_BIC_prev <- c(0, diff(bic_table$BIC)); bic_table$Two_Delta_BIC <- round(bic_table$BIC - bic_table$BIC[1], 2)
  bic_table$BayesFactor_vs_Null <- round(bic_neff - bic_neff[1], 2)
  bic_table$Evidence <- sapply(seq_len(nrow(bic_table)), function(idx_b) { # Kass & Raftery (1995) deviance-scale 2*ln(BF) = Delta_BIC_mclust vs G = 1
    if (idx_b == 1L) return("Null Model (G = 1)"); val <- bic_table$Two_Delta_BIC[idx_b]
    ev_str <- if (is.na(val) || val < 0) "Negative (favors G = 1)" else if (val < 2) "Weak / Inconclusive" else if (val < 6) "Positive Evidence" else if (val < 10) "Strong Evidence" else "Decisive Evidence (Species Boundary)"
    if (did_merge && bic_table$K[idx_b] > best_k) paste0(ev_str, sprintf(" [Ridgeline Merged -> K=%d]", best_k)) else ev_str
  })
  return(list(classification = cls, uncertainty = unc, best_k = best_k, model_name = best_name, selected_features = selected_cols, bic_table = bic_table))
}

# 4. Canonical Discriminant Analysis with Passive Sample Projection (Track 1: Unsupervised Morphotype_K1..Kn or Multi-Species)
run_canonical_discriminant_analysis <- function(df, feature_cols, target_taxa, gmm_res = NULL, intraspecific_mode = FALSE) {
  message("Configuring Weighted Ledoit-Wolf Shrinkage CDA (Tier_0_Type_Anchor / Tier_1_Gold / Tier_2_Silver / Promoted_Bronze vs. passive projection)...")
  feature_cols <- intersect(as.character(feature_cols), names(df))
  df$authority_weight <- compute_authority_precision_weights(df)
  authority_weight <- df$authority_weight
  is_tier1_pure <- if ("morph_tier_source" %in% names(df)) (is.na(df$morph_tier_source) | df$morph_tier_source != "Cauline_Only_Fallback") else rep(TRUE, nrow(df))
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  pc_req <- grep("^PC[0-9]+$", feature_cols, value = TRUE)
  has_pcs <- complete.cases(df[, if (length(pc_req) > 0) pc_req else feature_cols, drop = FALSE])
  screened_by_adaptive_anchors <- FALSE; df$is_confident_label_noise <- rep(FALSE, nrow(df))
  if (is_intra) {
    cl_u <- if (!is.null(gmm_res) && !is.null(gmm_res$classification)) as.integer(gmm_res$classification) else rep(1L, nrow(df))
    cl_u[has_pcs & is.na(cl_u)] <- 1L; n_intra_k <- length(unique(stats::na.omit(cl_u[has_pcs])))
    df$grouping_col <- ifelse(is.na(cl_u), NA_character_, paste0("Morphotype_K", cl_u)); min_grp_k <- min(2L, max(1L, n_intra_k))
    is_active <- (df$authority_weight > 0) & is_tier1_pure & !is.na(df$grouping_col) & has_pcs
    if (sum(table(df$grouping_col[is_active]) >= 2L) < min_grp_k) is_active <- is_tier1_pure & !is.na(df$grouping_col) & has_pcs
    if (sum(table(df$grouping_col[is_active]) >= 2L) < min_grp_k) is_active <- !is.na(df$grouping_col) & has_pcs
    act_tab <- table(df$grouping_col[is_active]); target_taxa <- sort(names(act_tab[act_tab >= 2L]))
    if (length(target_taxa) < min_grp_k) target_taxa <- sort(unique(stats::na.omit(df$grouping_col[is_active])))
    is_active <- is_active & (df$grouping_col %in% target_taxa); is_train_cand <- is_active
  } else {
    df$grouping_col <- df$species_standardized; used_unscreened_fb <- FALSE
    anc_cda <- select_adaptive_reference_anchors(df, feature_matrix = as.matrix(df[, if (length(pc_req) > 0) pc_req else feature_cols, drop = FALSE]), taxon_col = "grouping_col")
    df$is_promoted_bronze <- anc_cda$is_promoted_bronze; df$anchor_tier_source <- anc_cda$anchor_tier_source; df$is_confident_label_noise <- as.logical(anc_cda$is_confident_label_noise)
    df$authority_weight <- compute_authority_precision_weights(df); authority_weight <- df$authority_weight
    is_train_cand <- anc_cda$active_mask & is_tier1_pure & (df$grouping_col %in% target_taxa) & has_pcs
    tr_tab <- table(factor(df$grouping_col[is_train_cand], levels = target_taxa))
    if (sum(tr_tab >= 2L) >= 2L) { target_taxa <- names(tr_tab[tr_tab >= 2L]); is_train_cand <- is_train_cand & (df$grouping_col %in% target_taxa)
    } else if (any(tr_tab < 2L)) {
      if (length(unique(df$grouping_col[is_train_cand])) >= 2L) target_taxa <- intersect(target_taxa, unique(df$grouping_col[is_train_cand]))
      fallback_cand <- anc_cda$active_mask & (df$species_standardized %in% target_taxa) & has_pcs
      if (all(table(factor(df$grouping_col[fallback_cand], levels = target_taxa)) >= 2L)) {
        is_train_cand <- fallback_cand
      } else {
        is_train_cand <- (df$grouping_col %in% target_taxa) & has_pcs; used_unscreened_fb <- TRUE
      }
    }
    screened_by_adaptive_anchors <- (!used_unscreened_fb) && (sum(is_train_cand) >= 5L)
    is_active <- is_train_cand
  }
  keep_cda_comp <- vapply(feature_cols, function(fc) (fc %in% pc_req) || (mean(is.finite(suppressWarnings(as.numeric(df[[fc]][is_train_cand])))) >= 0.70), logical(1))
  if (sum(keep_cda_comp) >= 2L) feature_cols <- feature_cols[keep_cda_comp]
  X_full <- as.matrix(df[, feature_cols, drop = FALSE]); storage.mode(X_full) <- "double"
  for (j in seq_len(ncol(X_full))) if (any(is.na(X_full[, j]))) X_full[is.na(X_full[, j]), j] <- if (any(!is.na(X_full[is_train_cand, j]))) stats::median(X_full[is_train_cand, j], na.rm = TRUE) else 0.0
  feature_cols <- select_morphometric_subspace(X_full[is_train_cand, , drop = FALSE], group_vec = df$grouping_col[is_train_cand], aux_df = df[is_train_cand, , drop = FALSE], allow_supervised_f = FALSE)
  X <- X_full[, feature_cols, drop = FALSE]; p <- ncol(X)
  if (!screened_by_adaptive_anchors) {
    for (grp in target_taxa) {
      cand_k <- which(is_train_cand & df$grouping_col == grp); n_k <- length(cand_k)
      if (n_k >= 5) {
        sub_x <- X[cand_k, , drop = FALSE]
        d2_k <- tryCatch(stats::mahalanobis(sub_x, colMeans(sub_x), stats::cov(sub_x) + diag(1e-4, p)), error = function(e) rep(0, n_k))
        drop_k <- !is.na(d2_k) & (d2_k > stats::qchisq(0.95, df = p)); drop_k[.identify_type_specimen_vec(df[cand_k, , drop = FALSE])$is_type_specimen] <- FALSE
        if (sum(!drop_k) >= 2) is_train_cand[cand_k[drop_k]] <- FALSE
      }
    }
  }
  if (!is_intra) is_active <- is_train_cand
  df$spatial_balance_weight <- compute_spatial_balance_weights(df, active_mask = is_train_cand, recompute = TRUE)
  spatial_balance_weight <- df$spatial_balance_weight; eff_weight <- authority_weight * spatial_balance_weight
  active_idx <- which(is_train_cand); passive_idx <- which(!is_active)
  message(sprintf("Active Anchors (Complete): %d | Passive Projected: %d", length(active_idx), length(passive_idx)))
  if (length(active_idx) < length(target_taxa)) stop("Insufficient complete active reference anchors for CDA.")

  X_a <- X[active_idx, , drop = FALSE]; y_a <- df$grouping_col[active_idx]; eff_weight <- eff_weight[active_idx]; if (all(eff_weight <= 0)) eff_weight <- rep(1.0, length(active_idx)); w_a <- eff_weight; N_a <- nrow(X_a)
  X_mat <- X_a; grp_factor <- droplevels(factor(y_a))
  loocv_accuracy <- if (is_intra && all(grepl("^Morphotype_K", levels(grp_factor)))) NA_real_ else tryCatch(compute_weighted_shrinkage_cda_loocv(X_a, grp_factor, w_a, coords_df = df[active_idx, , drop = FALSE], X_cand = X_full[active_idx, , drop = FALSE]), error = function(e) NA_real_)
  if (is.na(loocv_accuracy) && !is_intra && requireNamespace("MASS", quietly = TRUE) && nlevels(grp_factor) >= 2) {
    lda_cv <- tryCatch(MASS::lda(x = X_mat, grouping = grp_factor, prior = rep(1.0 / nlevels(grp_factor), nlevels(grp_factor)), CV = TRUE), error = function(e) NULL)
    if (!is.null(lda_cv) && !is.null(lda_cv$class)) loocv_accuracy <- round(mean(as.character(lda_cv$class) == as.character(grp_factor), na.rm = TRUE), 4)
  }
  message(sprintf("CDA Gathering/Deme-Blocked Morphological CV Accuracy (loocv_accuracy): %s", ifelse(is.na(loocv_accuracy), "NA (post-clustering descriptive ordination)", sprintf("%.4f", loocv_accuracy))))

  g_levels <- target_taxa; g <- length(g_levels); grand_mean <- colSums(X_a * w_a) / sum(w_a)
  if (g == 1L) {
    X_cent1 <- sweep(X_a, 2L, grand_mean, "-"); S_w1 <- crossprod(X_cent1 * sqrt(w_a)) / max(sum(w_a) - 1.0, 1.0)
    S_reg1 <- compute_ledoit_wolf_cov(S_w1, X_centered = X_cent1, n_obs = N_a, weights = w_a); eig1 <- eigen(S_reg1, symmetric = TRUE)
    X_p1 <- X; for (ci in seq_len(p)) if (any(is.na(X_p1[, ci]))) X_p1[is.na(X_p1[, ci]), ci] <- grand_mean[ci]
    Z1 <- sweep(X_p1, 2L, grand_mean, "-") %*% eig1$vectors[, seq_len(min(3L, p)), drop = FALSE]
    while (ncol(Z1) < 3L) Z1 <- cbind(Z1, 0.0); colnames(Z1) <- c("Can1", "Can2", "Can3")
    vp1 <- round(100.0 * pmax(eig1$values, 0) / max(sum(pmax(eig1$values, 0)), 1e-12), 2); sc1 <- suppressWarnings(stats::cor(X_a, Z1[active_idx, 1:2, drop = FALSE])); sc1[!is.finite(sc1)] <- 0.0
    c_mat1 <- matrix(0.0, nrow = 1L, ncol = 3L, dimnames = list(g_levels, colnames(Z1))); pp1 <- matrix(1.0, nrow = nrow(df), ncol = 1L, dimnames = list(NULL, g_levels))
    return(list(eigenvalues = pmax(eig1$values, 0), variance_pct = vp1, canonical_scores = Z1, centroids = c_mat1, canonical_scores_full = Z1[, 1, drop = FALSE], centroids_full = c_mat1[, 1, drop = FALSE],
                num_canonical_axes = 1L, shrinkage_lambda = attr(S_reg1, "shrinkage_lambda"), selected_features = feature_cols, predicted_taxon = rep(g_levels[1L], nrow(df)), posterior_prob = rep(1.0, nrow(df)),
                second_taxon = rep(NA_character_, nrow(df)), second_posterior_prob = rep(0.0, nrow(df)), post_probs = pp1, posterior_margin = rep(1.0, nrow(df)), is_active = is_active,
                is_promoted_bronze = df$is_promoted_bronze, is_confident_label_noise = df$is_confident_label_noise, anchor_tier_source = df$anchor_tier_source, authority_weight = df$authority_weight,
                structure_loadings = sc1, loocv_accuracy = loocv_accuracy, grouping_col = df$grouping_col, X_proj = X_p1, S_reg = S_reg1))
  }
  B <- matrix(0, p, p); W <- matrix(0, p, p); X_centered <- matrix(0, N_a, p, dimnames = dimnames(X_a))
  group_means <- matrix(0, g, p, dimnames = list(g_levels, feature_cols)); group_counts <- integer(g)
  for (k in seq_along(g_levels)) {
    grp <- g_levels[k]; k_mask <- (y_a == grp); sub_x <- X_a[k_mask, , drop = FALSE]; w_k <- eff_weight[k_mask]; group_counts[k] <- nrow(sub_x)
    if (nrow(sub_x) > 0) {
      m_k <- colSums(sub_x * w_k) / sum(w_k); group_means[k, ] <- m_k; diff_m <- m_k - grand_mean
      B <- B + sum(w_k) * (diff_m %*% t(diff_m)); diff_x <- t(t(sub_x) - m_k); W <- W + (t(diff_x) %*% (diff_x * w_k))
      X_centered[k_mask, ] <- diff_x
    }
  }

  denom_w <- max(sum(w_a) - sum(vapply(g_levels, function(g_k) { wk <- w_a[y_a == g_k]; if (sum(wk) > 0) sum(wk^2) / sum(wk) else 1.0 }, numeric(1))), 1.0)
  S_w <- W / denom_w
  S_reg <- compute_ledoit_wolf_cov(S_w, X_centered = X_centered, n_obs = N_a, weights = w_a)
  k_full <- max(1L, min(ncol(X_mat), g - 1L))
  eig_cda <- solve_shrinkage_cda_eigen(S_reg, B, max_axes = k_full)
  eig_vals <- eig_cda$values; eig_vecs_full <- eig_cda$vectors
  var_pct <- round((eig_vals / max(sum(eig_vals), 1e-12)) * 100, 2)
  message(sprintf("=== CDA Canonical Variates (K_full = %d, Ledoit-Wolf lambda = %.5f) ===", k_full, attr(S_reg, "shrinkage_lambda")))
  for (j in seq_len(k_full)) message(sprintf("  Can%d: Eigenvalue = %7.4f (%5.2f%% variation)", j, eig_vals[j], var_pct[j]))

  Z_anchors_full <- sweep(X_a, 2L, grand_mean, "-") %*% eig_vecs_full
  struct_corr <- suppressWarnings(stats::cor(X_a, Z_anchors_full)); struct_corr[!is.finite(struct_corr)] <- 0.0
  X_proj <- X
  for (col_idx in seq_len(p)) { na_mask <- is.na(X_proj[, col_idx]); if (any(na_mask)) X_proj[na_mask, col_idx] <- grand_mean[col_idx] }
  Z_all_full <- sweep(X_proj, 2L, grand_mean, "-") %*% eig_vecs_full
  centroids_full <- sweep(group_means, 2L, grand_mean, "-") %*% eig_vecs_full
  colnames(Z_all_full) <- colnames(Z_anchors_full) <- colnames(centroids_full) <- paste0("Can", seq_len(k_full))
  rownames(centroids_full) <- g_levels

  priors <- rep(1/g, g); post_probs <- matrix(0.0, nrow = nrow(df), ncol = g, dimnames = list(NULL, g_levels))
  for (k in seq_along(g_levels)) { diff_z <- sweep(Z_all_full, 2L, centroids_full[k, ], "-"); post_probs[, k] <- log(priors[k]) - 0.5 * rowSums(diff_z^2) }
  post_probs <- exp(post_probs - apply(post_probs, 1, max, na.rm = TRUE))
  post_probs <- post_probs / pmax(rowSums(post_probs, na.rm = TRUE), 1e-12); pred_idx <- max.col(post_probs, ties.method = "first")
  posterior_prob <- apply(post_probs, 1, max)
  second_idx <- if (g >= 2L) vapply(seq_len(nrow(post_probs)), function(i) { oth <- setdiff(seq_len(g), pred_idx[i]); oth[which.max(post_probs[i, oth])] }, integer(1)) else rep(NA_integer_, nrow(df))
  second_taxon <- if (g >= 2L) colnames(post_probs)[second_idx] else rep(NA_character_, nrow(df))
  second_posterior_prob <- if (g >= 2L) post_probs[cbind(seq_len(nrow(post_probs)), second_idx)] else rep(0.0, nrow(df))
  posterior_margin <- posterior_prob - second_posterior_prob

  num_ord <- min(3L, k_full); Z_all <- Z_all_full[, seq_len(num_ord), drop = FALSE]; canonical_centroids <- centroids_full[, seq_len(num_ord), drop = FALSE]
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
              canonical_scores_full = Z_all_full, centroids_full = centroids_full, num_canonical_axes = k_full,
              shrinkage_lambda = attr(S_reg, "shrinkage_lambda"), selected_features = feature_cols,
              predicted_taxon = g_levels[pred_idx], posterior_prob = posterior_prob, second_taxon = second_taxon,
              second_posterior_prob = second_posterior_prob, post_probs = post_probs,
              posterior_margin = posterior_margin, is_active = is_active, is_promoted_bronze = df$is_promoted_bronze, is_confident_label_noise = df$is_confident_label_noise, anchor_tier_source = df$anchor_tier_source,
              authority_weight = df$authority_weight, structure_loadings = struct_corr, loocv_accuracy = loocv_accuracy, grouping_col = df$grouping_col, X_proj = X_proj, S_reg = S_reg))
}

# 5. Outlier Flagging & Herbarium Misidentification Auditing
audit_misidentifications <- function(df, cda_res, gmm_res = NULL, selected_features = NULL) {
  message("Auditing specimens for morphospace outliers & label discordances...")
  n <- nrow(df); flags <- logical(n); triage <- character(n); reasons <- character(n)
  grp_vec <- if (!is.null(cda_res$grouping_col)) cda_res$grouping_col else df$species_standardized
  is_intra_audit <- length(unique(stats::na.omit(df$species_standardized))) < 2 || any(grepl("^Morphotype_K", grp_vec))
  sel_feats <- unique(c(selected_features, cda_res$selected_features, gmm_res$selected_features))
  if (length(sel_feats) == 0L) sel_feats <- discover_sorted_morph_cols(df, "^PC[0-9]+$")
  sel_feats <- intersect(as.character(sel_feats), names(df))
  pc_c <- sel_feats[vapply(sel_feats, function(fc) is.numeric(df[[fc]]) && sum(is.finite(df[[fc]])) >= 10L, logical(1))]
  if (length(pc_c) < 2L) pc_c <- discover_sorted_morph_cols(df, "^PC[0-9]+$")
  mahal_out <- if (is_intra_audit && length(pc_c) >= 2 && nrow(df) >= 10) {
    Xm <- as.matrix(df[, pc_c, drop = FALSE]); act_m <- which(as.logical(cda_res$is_active)); if (length(act_m) < 5L) act_m <- seq_len(n)
    for (jc in seq_len(ncol(Xm))) Xm[!is.finite(Xm[, jc]), jc] <- stats::median(Xm[act_m, jc], na.rm = TRUE)
    rob_f <- if (exists(".fit_fast_mcd_cov", mode = "function")) .fit_fast_mcd_cov(Xm[act_m, , drop = FALSE]) else NULL
    mu_m <- if (!is.null(rob_f)) rob_f$center else apply(Xm[act_m, , drop = FALSE], 2, stats::median, na.rm = TRUE)
    S_m <- if (!is.null(rob_f)) (rob_f$cov + diag(1e-5, ncol(Xm))) else (stats::cov(Xm[act_m, , drop = FALSE]) + diag(1e-5, ncol(Xm)))
    d2 <- tryCatch(stats::mahalanobis(Xm, mu_m, S_m), error = function(e) rep(0, n))
    out_m <- !is.na(d2) & (d2 > stats::qchisq(0.975, df = ncol(Xm))) & stats::complete.cases(df[, pc_c, drop = FALSE])
    out_m[.identify_type_specimen_vec(df)$is_type_specimen] <- FALSE; out_m
  } else rep(FALSE, n)
  for (i in seq_len(nrow(df))) {
    raw_sp <- as.character(grp_vec[i]); if (is.na(raw_sp) || !nzchar(raw_sp)) raw_sp <- as.character(df$species_standardized[i])
    pred_sp <- as.character(cda_res$predicted_taxon[i]); post_p <- cda_res$posterior_prob[i]; is_act <- isTRUE(cda_res$is_active[i])
    ar_i <- .format_step04_anchor_audit_reason(i, df, cda_res, raw_sp, pred_sp, post_p, is_act)
    flags[i] <- ar_i$flag; triage[i] <- ar_i$triage; reasons[i] <- ar_i$reason
    if (is_intra_audit && !flags[i] && isTRUE(mahal_out[i])) { flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- "Intraspecific_Mahalanobis_Morphospace_Outlier" }
  }
  tc <- compute_type_specimen_morph_control(df, cda_res$X_proj, cda_res$S_reg, cda_res, gmm_res, flags = flags, triage = triage, reasons = reasons)
  flags <- tc$flags; triage <- tc$triage; reasons <- tc$reasons
  for (cn in c("type_status", "is_type_specimen", "has_empirical_type_anchor", "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance", "is_typical_morphotype", "min_congener_mahal_d2", "min_congener_taxon", "open_set_novelty_p")) df[[cn]] <- tc[[cn]]
  p_audit <- max(2L, if (!is.null(cda_res$X_proj)) ncol(cda_res$X_proj) else length(pc_c))
  is_type_audit <- as.logical(df$is_type_specimen); is_type_audit[is.na(is_type_audit)] <- FALSE
  is_osr_novel <- !is.na(df$open_set_novelty_p) & (df$open_set_novelty_p < 0.01)
  misid_flags <- !is_type_audit & ((!is_intra_audit & flags & !is_osr_novel) | (is_intra_audit & is.finite(df$min_congener_mahal_d2) & is.finite(df$type_morph_mahal_d2) & (df$min_congener_mahal_d2 < stats::qchisq(0.95, df = p_audit)) & (df$type_morph_mahal_d2 > stats::qchisq(0.99, df = p_audit))))
  outlier_flags <- !is_type_audit & (is_osr_novel | mahal_out | (!is.na(df$type_morph_mahal_d2) & (df$type_morph_mahal_d2 > stats::qchisq(0.995, df = p_audit)) & !misid_flags))
  novel_only_idx <- which(outlier_flags & !misid_flags & (triage %in% c("CLEAN", "LOW")))
  if (length(novel_only_idx) > 0L) { triage[novel_only_idx] <- "HIGH"; reasons[novel_only_idx] <- "Open_Set_Morphospace_Outlier_Novelty" }
  if (!is.null(cda_res$is_promoted_bronze)) df$is_promoted_bronze <- as.logical(cda_res$is_promoted_bronze)
  df$is_confident_label_noise <- if (!is.null(cda_res$is_confident_label_noise)) as.logical(cda_res$is_confident_label_noise) else rep(FALSE, n)
  if (!is.null(cda_res$anchor_tier_source)) df$anchor_tier_source <- cda_res$anchor_tier_source
  df$authority_weight <- if (!is.null(cda_res$authority_weight)) cda_res$authority_weight else compute_authority_precision_weights(df)
  if (!is.null(gmm_res$classification)) df$gmm_cluster <- gmm_res$classification
  if (!is.null(gmm_res$uncertainty)) df$gmm_uncertainty <- round(gmm_res$uncertainty, 4)
  df$loocv_accuracy <- cda_res$loocv_accuracy
  df$selected_morph_features <- paste(sel_feats, collapse = ";")
  df$cda_predicted_taxon <- cda_res$predicted_taxon; df$cda_posterior_prob <- round(cda_res$posterior_prob, 4)
  df$cda_second_taxon <- if (!is.null(cda_res$second_taxon)) as.character(cda_res$second_taxon) else NA_character_
  df$cda_second_posterior_prob <- round(cda_res$second_posterior_prob, 4); df$cda_posterior_margin <- round(cda_res$posterior_margin, 4)
  df$can1 <- round(cda_res$canonical_scores[, 1], 5)
  df$can2 <- if (ncol(cda_res$canonical_scores) >= 2) round(cda_res$canonical_scores[, 2], 5) else 0.0
  df$can3 <- if (ncol(cda_res$canonical_scores) >= 3) round(cda_res$canonical_scores[, 3], 5) else 0.0
  df$is_passive <- !cda_res$is_active; df$morphometric_outlier_flag <- as.logical(outlier_flags); df$misidentification_flag <- as.logical(misid_flags)
  df$triage_priority <- triage; df$discordance_reason <- reasons
  if ("AR_inv" %in% names(df)) df$reproductive_trait_status <- ifelse(is.na(df$AR_inv), "Unexpanded_or_Absent", "Measured")
  message(sprintf("Audit Complete: %d morphospace outliers/discordances flagged (%d HIGH, %d MEDIUM)", sum(flags), sum(triage == "HIGH"), sum(triage == "MEDIUM")))
  return(df)
}

# 5b. Institutional Batch-Effect Audit (Taxon-Adjusted Marginal ANCOVA, Partial R-Squared & BH FDR)
audit_institutional_batch_effects <- function(morpho_df, output_csv = "outputs/reports/institutional_batch_effect_audit.csv", min_obs = 15, selected_features = NULL) {
  if (!("institutionCode" %in% names(morpho_df))) return(NULL)
  dyn_sel <- if (!is.null(selected_features) && length(selected_features) > 0L) intersect(as.character(selected_features), names(morpho_df)) else discover_sorted_morph_cols(morpho_df, "^PC[0-9]+$")
  test_feats <- unique(c(intersect(c("PC1", "PC2"), names(morpho_df)), dyn_sel))
  test_feats <- test_feats[vapply(test_feats, function(fc) is.numeric(morpho_df[[fc]]) && sum(is.finite(morpho_df[[fc]])) >= 10L, logical(1))]
  valid_df <- morpho_df[!is.na(morpho_df$institutionCode) & trimws(morpho_df$institutionCode) != "", , drop = FALSE]
  clean_b <- (!("authority_weight" %in% names(valid_df)) | (!is.na(valid_df$authority_weight) & valid_df$authority_weight > 0)) & (!("morph_tier_source" %in% names(valid_df)) | is.na(valid_df$morph_tier_source) | valid_df$morph_tier_source == "Basal_Rosette")
  if (sum(clean_b, na.rm = TRUE) >= 10L && length(unique(valid_df$institutionCode[clean_b])) >= 2L) valid_df <- valid_df[clean_b, , drop = FALSE]
  if (nrow(valid_df) < 10 || length(test_feats) == 0L) return(NULL)
  inst_counts <- table(valid_df$institutionCode); rep_insts <- names(inst_counts[inst_counts >= min_obs])
  audit_subset <- if (length(rep_insts) >= 2) valid_df[valid_df$institutionCode %in% rep_insts, , drop = FALSE] else valid_df
  audit_subset$institutionCode <- droplevels(factor(audit_subset$institutionCode))
  if (!("species_standardized" %in% names(audit_subset))) audit_subset$species_standardized <- if ("scientificName" %in% names(audit_subset)) audit_subset$scientificName else "Unassigned_Packera"
  audit_subset$species_standardized <- droplevels(factor(audit_subset$species_standardized))
  if (nlevels(audit_subset$institutionCode) < 2) return(NULL)
  fits <- list(); rows_out <- list(); raw_ps <- numeric(0)
  for (feat in test_feats) {
    sub_f <- audit_subset[is.finite(audit_subset[[feat]]), , drop = FALSE]; sub_f$institutionCode <- droplevels(sub_f$institutionCode); sub_f$species_standardized <- droplevels(sub_f$species_standardized)
    if (nrow(sub_f) < 6L || nlevels(sub_f$institutionCode) < 2L || stats::var(sub_f[[feat]], na.rm = TRUE) <= 1e-12) next
    fml <- if (nlevels(sub_f$species_standardized) >= 2L) stats::as.formula(paste(feat, "~ species_standardized + institutionCode")) else stats::as.formula(paste(feat, "~ institutionCode"))
    w_bf <- (if ("authority_weight" %in% names(sub_f)) suppressWarnings(as.numeric(sub_f$authority_weight)) else rep(1.0, nrow(sub_f))) * (if ("spatial_balance_weight" %in% names(sub_f)) suppressWarnings(as.numeric(sub_f$spatial_balance_weight)) else rep(1.0, nrow(sub_f))); w_bf[!is.finite(w_bf) | w_bf <= 0] <- 1.0
    fit_f <- stats::lm(fml, data = sub_f, weights = w_bf); s_f <- summary(fit_f); aov_f <- stats::anova(fit_f)
    ii <- which(rownames(aov_f) == "institutionCode")[1]; ri <- which(rownames(aov_f) == "Residuals")[1]; if (is.na(ii) || is.na(ri)) next
    r2_f <- aov_f$`Sum Sq`[ii] / max(aov_f$`Sum Sq`[ii] + aov_f$`Sum Sq`[ri], 1e-12); pv_f <- as.numeric(aov_f$`Pr(>F)`[ii])
    fits[[feat]] <- fit_f; raw_ps <- c(raw_ps, pv_f)
    rows_out[[length(rows_out) + 1L]] <- data.frame(
      Response = c(feat, feat), Term = c("institutionCode", "Residuals"), Df = c(aov_f$Df[ii], aov_f$Df[ri]),
      Sum_Sq = round(c(aov_f$`Sum Sq`[ii], aov_f$`Sum Sq`[ri]), 4), Mean_Sq = round(c(aov_f$`Mean Sq`[ii], aov_f$`Mean Sq`[ri]), 4),
      F_value = c(round(aov_f$`F value`[ii], 4), NA), p_value = c(signif(pv_f, 4), NA), P_value = c(signif(pv_f, 4), NA), p_value_bh = c(NA_real_, NA_real_),
      R_squared = round(c(r2_f, NA), 4), Adj_R_squared = round(c(s_f$adj.r.squared, NA), 4),
      Batch_Effect_Status = c(ifelse(r2_f < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA), stringsAsFactors = FALSE)
  }
  if (length(rows_out) == 0L) return(NULL)
  anova_df <- do.call(rbind, rows_out); anova_df$p_value_bh[anova_df$Term == "institutionCode"] <- signif(stats::p.adjust(raw_ps, method = "BH"), 4)
  if (!is.null(output_csv) && nzchar(output_csv)) { dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE); tmp_ba <- paste0(output_csv, ".", Sys.getpid(), ".tmp"); write.csv(anova_df, file = tmp_ba, row.names = FALSE); file.rename(tmp_ba, output_csv) }
  return(invisible(list(lm_pc1 = fits[["PC1"]], lm_pc2 = fits[["PC2"]], models = fits, anova_table = anova_df)))
}

# 6. Publication-Quality 2D CDA & 2D GMM Biplot Vector Graphic (PDF)
generate_cda_biplot_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL, manova_res = NULL) {
  if (!exists("render_cda_biplot_4panel_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
  if (exists("render_cda_biplot_4panel_pdf", mode = "function")) render_cda_biplot_4panel_pdf(audit_df, cda_res, gmm_res, output_pdf, feature_cols, manova_res = manova_res)
}

# 7. Main Workflow Orchestrator
.read_csv_s04 <- function(p) if (requireNamespace("data.table", quietly = TRUE) && isTRUE(file.info(p)$size > 1048576)) as.data.frame(data.table::fread(p, data.table = FALSE, showProgress = FALSE)) else read.csv(p, stringsAsFactors = FALSE)
.write_csv_s04 <- function(d, p, na = "NA") if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(d, file = p, na = na, quote = "auto", logical01 = FALSE) else write.csv(d, file = p, row.names = FALSE, na = na)
run_gmm_morphotools_pipeline <- function(opts) {
  vouchers_df <- if (file.exists(opts$vouchers)) .read_csv_s04(opts$vouchers) else NULL
  db_efa <- if (exists("fetch_step04_morphometrics_from_db", mode = "function")) fetch_step04_morphometrics_from_db(opts$vouchers, opts$harmonics, active_catalogs = if (!is.null(vouchers_df) && "catalogNumber" %in% names(vouchers_df)) vouchers_df$catalogNumber else NULL) else NULL
  if (!is.null(db_efa) && nrow(db_efa) > 0L) {
    na_dyn <- grep("^PC[0-9]+$|^[ABCD][0-9]+$", names(db_efa), value = TRUE)
    drop_na_dyn <- na_dyn[vapply(na_dyn, function(cl) !any(is.finite(suppressWarnings(as.numeric(db_efa[[cl]])))), logical(1))]
    if (length(drop_na_dyn) > 0L) db_efa <- db_efa[, setdiff(names(db_efa), drop_na_dyn), drop = FALSE]
  }
  db_pcs <- if (!is.null(db_efa) && nrow(db_efa) > 0L) discover_sorted_morph_cols(db_efa, "^PC[0-9]+$") else character(0)
  if (!is.null(db_efa) && nrow(db_efa) > 0L && length(db_pcs) >= 2L) { efa_df <- db_efa } else {
    if (!file.exists(opts$harmonics)) stop("Harmonics file not found: ", opts$harmonics)
    efa_df <- .read_csv_s04(opts$harmonics)
  }
  pc_cols <- discover_sorted_morph_cols(efa_df, "^PC[0-9]+$"); ect_pc_cols <- discover_sorted_morph_cols(efa_df, "^ECT_PC[0-9]+$")
  pca_cols <- if (length(pc_cols) >= 2L) pc_cols else ect_pc_cols
  if (length(pca_cols) < 2L) stop("Missing PCA columns in EFA/ECT table (need >= 2 PC/ECT_PC axes).")
  has_primary_pc <- complete.cases(efa_df[, pca_cols[seq_len(min(2L, length(pca_cols)))], drop = FALSE])
  has_ect_topo <- (if (length(ect_pc_cols) >= 2L) complete.cases(efa_df[, ect_pc_cols[1:2], drop = FALSE]) else rep(FALSE, nrow(efa_df))) | (if (all(c("solidity", "dissection_index") %in% names(efa_df))) is.finite(suppressWarnings(as.numeric(efa_df$solidity))) & is.finite(suppressWarnings(as.numeric(efa_df$dissection_index))) else rep(FALSE, nrow(efa_df)))
  valid_idx <- which(has_primary_pc | has_ect_topo); closed_df <- efa_df[valid_idx, , drop = FALSE]
  # Ensure primary GMM and CDA run on specimen-level table (1 row = 1 voucher collection event; group_by(catalogNumber) & plant_individual_id)
  if (!"plant_individual_id" %in% names(closed_df)) closed_df$plant_individual_id <- 0L
  spec_keys <- paste0(closed_df$catalogNumber, "_", closed_df$plant_individual_id)
  if (any(duplicated(spec_keys))) {
    disagg_leaf_df <- closed_df; pc_sym_cols <- unique(c(pc_cols, ect_pc_cols)); pet_cols <- intersect(c("petiole_length_mm", "petiole_width_mm", "blade_length_mm", "blade_width_mm", "petiole_to_blade_ratio"), names(closed_df))
    lin_med_cols <- intersect(c("median_area_mm2", "median_perimeter_mm", "median_length_mm", "median_width_mm", "median_circularity", "median_aspect_ratio", "median_solidity", "median_dissection_index", "area_mm2", "perimeter_mm", "length_mm", "width_mm", "circularity"), names(closed_df))
    num_cols <- intersect(c(pc_sym_cols, pet_cols, lin_med_cols, "plant_height_cm", "scape_to_basal_leaf_ratio", "foliar_variance", "leaf_count", "aspect_ratio", "solidity", "dissection_index", "ect_lobe_complexity", "area_px", "centroid_size", "log_centsize", "harmonic_roughness", "leaf_quality_weight"), names(closed_df))
    other_num <- setdiff(num_cols, c(pc_sym_cols, pet_cols)); has_w <- "leaf_quality_weight" %in% names(closed_df); has_lz <- "leaf_zone" %in% names(closed_df)
    meta_cols <- setdiff(names(closed_df), c(num_cols, "catalogNumber", "plant_individual_id", "shape_id", "leaf_id", "mask_source", "morph_tier_source"))
    grp_list <- split(closed_df, list(closed_df$catalogNumber, closed_df$plant_individual_id), drop = TRUE)
    agg_rows <- (if (length(grp_list) >= 200L && exists("packera_mclapply", mode = "function")) packera_mclapply else lapply)(grp_list, function(grp) {
      m_src <- "Basal_Rosette"
      if (has_lz) { lz_g <- tolower(trimws(as.character(grp$leaf_zone))); if (any(!is.na(lz_g) & lz_g == "basal")) grp <- grp[!is.na(lz_g) & lz_g == "basal", , drop = FALSE] else if (all(!is.na(lz_g) & lz_g == "cauline")) m_src <- "Cauline_Only_Fallback" }
      leaf_quality_weight <- if (has_w) suppressWarnings(as.numeric(grp$leaf_quality_weight)) else rep(1.0, nrow(grp))
      leaf_quality_weight[is.na(leaf_quality_weight) | leaf_quality_weight <= 0] <- 1.0
      is_type_g <- .identify_type_specimen_vec(grp)$is_type_specimen
      is_contour_outlier <- if ("is_contour_outlier" %in% names(grp)) as.logical(grp$is_contour_outlier) else rep(FALSE, nrow(grp)); is_contour_outlier[is.na(is_contour_outlier) | is_type_g] <- FALSE
      if (nrow(grp) >= 3 && !any(is_type_g)) {
        pc_mat <- as.matrix(grp[, intersect(pca_cols, names(grp)), drop = FALSE]); med_pc <- apply(pc_mat, 2, stats::median, na.rm = TRUE)
        dist <- sqrt(rowSums(sweep(pc_mat, 2, med_pc, "-")^2, na.rm = TRUE)); iqr_d <- stats::IQR(dist, na.rm = TRUE); med_d <- stats::median(dist, na.rm = TRUE)
        drop_leaf <- ((dist > 2.0 * iqr_d & iqr_d > 0) | (dist > 2.5 * med_d)) & (is_contour_outlier | nrow(grp) >= 3)
        if (sum(!drop_leaf) < 2) drop_leaf[utils::head(order(leaf_quality_weight, -dist, decreasing = TRUE), 2)] <- FALSE
        grp <- grp[!drop_leaf, , drop = FALSE]
      } else if (any(!is_contour_outlier)) { grp <- grp[!is_contour_outlier, , drop = FALSE] }
      out <- grp[1, , drop = FALSE]; for (col in c(pc_sym_cols, other_num)) { vv <- suppressWarnings(as.numeric(grp[[col]])); vv <- vv[is.finite(vv)]; out[[col]] <- if (length(vv) > 0L) stats::median(vv) else NA_real_ }
      dm_sub <- summarize_otu_dissection_metrics(grp)
      for (d_nm in c("basal_leaf_count", "q25_solidity", "dissection_index", "dissected_leaf_fraction")) if (!is.na(dm_sub[[d_nm]])) out[[d_nm]] <- dm_sub[[d_nm]]
      for (col in pet_cols) out[[col]] <- stats::quantile(suppressWarnings(as.numeric(grp[[col]])), probs = 0.75, na.rm = TRUE, names = FALSE)
      if ("leaf_count" %in% names(out)) out$leaf_count <- nrow(grp)
      out$morph_tier_source <- m_src
      for (col in meta_cols) { v <- stats::na.omit(grp[[col]]); out[[col]] <- if (length(v) > 0) v[1] else grp[[col]][1] }
      out
    })
    closed_df <- if (requireNamespace("data.table", quietly = TRUE) && length(agg_rows) >= 200L) as.data.frame(data.table::rbindlist(agg_rows, use.names = TRUE, fill = TRUE)) else do.call(rbind, agg_rows); rownames(closed_df) <- NULL
    if (!("cauline_leaf_count" %in% names(closed_df)) && exists("compute_separate_cauline_otu_track", mode = "function")) closed_df <- compute_separate_cauline_otu_track(disagg_leaf_df, closed_df)
  } else {
    if (!"morph_tier_source" %in% names(closed_df)) closed_df$morph_tier_source <- "Basal_Rosette"
  }
  if (!is.null(vouchers_df)) {
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), , drop = FALSE]; v_idx <- match(as.character(closed_df$catalogNumber), as.character(v_sub$catalogNumber))
    for (col in intersect(names(v_sub), c("scientificName", "species_standardized", "species_raw", "identifiedBy", "determiner_tier", "type_status", "is_type_specimen", "institutionCode", "county", "stateProvince", "decimalLatitude", "decimalLongitude", "doy", "regional_group", "spatial_balance_weight"))) { vv <- v_sub[[col]][v_idx]; if (!col %in% names(closed_df)) closed_df[[col]] <- vv else { m_v <- !is.na(vv) & nzchar(trimws(as.character(vv))) & (col %in% c("scientificName", "species_standardized", "species_raw", "determiner_tier", "type_status", "is_type_specimen", "regional_group") | is.na(closed_df[[col]])); closed_df[[col]][m_v] <- vv[m_v] } }
  }
  sp_cand <- rep("Unknown", nrow(closed_df)); for (sc in intersect(c("species_standardized", "scientificName", "species_raw"), names(closed_df))) { v <- standardize_packera_taxon(closed_df[[sc]]); m <- sp_cand == "Unknown" & !is.na(v) & v != "Unknown"; sp_cand[m] <- v[m] }
  closed_df$species_standardized <- sp_cand; closed_df$authority_weight <- compute_authority_precision_weights(closed_df)
  if ("regional_group" %in% names(closed_df) && exists("standardize_epa_l3_label_r", mode = "function")) closed_df$regional_group <- standardize_epa_l3_label_r(closed_df$regional_group)
  closed_df$spatial_balance_weight <- compute_spatial_balance_weights(closed_df, active_mask = (closed_df$authority_weight > 0), recompute = TRUE)
  if (!isTRUE(opts$outline_only)) {
    if ("AR_inv" %in% names(closed_df)) closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$AR_inv))
    closed_df <- merge_reproductive_traits_by_individual(closed_df, opts$reproductive)
    closed_df <- merge_manifest_petiole_traits_by_individual(closed_df, file.path(dirname(opts$harmonics), "extracted_leaf_manifest.csv"))
  }
  allom_res <- test_and_correct_allometry(closed_df, correct_mode = opts$correct_allometry, permutations = opts$permutations, seed = opts$seed, output_report = opts$allometry_report, output_plot = opts$allometry_plot, output_pairwise = opts$allometry_pairwise)
  if (!is.null(allom_res$df)) closed_df <- allom_res$df
  pc_cols <- discover_sorted_morph_cols(closed_df, "^PC[0-9]+$"); ect_pc_cols <- discover_sorted_morph_cols(closed_df, "^ECT_PC[0-9]+$"); pca_cols <- if (length(pc_cols) >= 2L) pc_cols else ect_pc_cols
  if (exists("assign_spatial_gathering_demes", mode = "function")) closed_df <- assign_spatial_gathering_demes(closed_df, deme_radius_km = 10.0)
  feature_cols <- pca_cols
  if (!isTRUE(opts$outline_only)) {
    for (tr_col in c("AR_inv", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio", "petiole_length_mm", "plant_height_cm")) if (tr_col %in% names(closed_df)) {
      closed_df[[tr_col]] <- suppressWarnings(as.numeric(closed_df[[tr_col]])); valid_tr <- which(!is.na(closed_df[[tr_col]]) & is.finite(closed_df[[tr_col]]))
      ref_tr <- intersect(valid_tr, which(closed_df$authority_weight > 0)); if (length(ref_tr) < 3L) ref_tr <- valid_tr
      if (length(valid_tr) >= 5) {
        std_col <- paste0(tr_col, "_std"); closed_df[[std_col]] <- NA_real_; m_tr <- mean(closed_df[[tr_col]][ref_tr]); s_tr <- stats::sd(closed_df[[tr_col]][ref_tr])
        closed_df[[std_col]][valid_tr] <- if (is.finite(s_tr) && s_tr > 1e-6) round((closed_df[[tr_col]][valid_tr] - m_tr) / s_tr, 5) else 0.0
        # Restrict core multivariate feature matrices strictly to scale-invariant ratios (excluding ruler-calibrated mm/cm traits that spike on uncalibrated sheets)
        if (tr_col %in% c("AR_inv", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio")) {
          act_all <- which(closed_df$authority_weight > 0); if (length(act_all) < 3L) act_all <- seq_len(nrow(closed_df))
          if ((length(intersect(valid_tr, act_all)) / max(1L, length(act_all))) >= 0.70) feature_cols <- unique(c(feature_cols, std_col))
        }
      }
    }
  }
  caul_het_raw <- c("basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "heteroblasty_efa_distance", "heteroblasty_ect_distance", "cauline_dissection_index_med", "cauline_solidity_med", "cauline_ect_lobe_complexity", "Cauline_PC1", "Cauline_ECT_PC1")
  veg_all <- intersect(c("log_centsize", "aspect_ratio", "solidity", "q25_solidity", "dissection_index", "within_rosette_efa_dispersion", "ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", caul_het_raw), names(closed_df))
  ref_std_idx <- which(closed_df$authority_weight > 0); if (length(ref_std_idx) < 3L) ref_std_idx <- seq_len(nrow(closed_df))
  cal_ok_s04 <- if ("is_calibrated_mm" %in% names(closed_df)) (!is.na(closed_df$is_calibrated_mm) & as.logical(closed_df$is_calibrated_mm)) else rep(TRUE, nrow(closed_df))
  if ("scale_source" %in% names(closed_df)) cal_ok_s04 <- cal_ok_s04 & !(tolower(trimws(as.character(closed_df$scale_source))) %in% c("mp_estimated", "predicted_from_megapixels", "uncalibrated"))
  for (vc in veg_all) {
    vv <- suppressWarnings(as.numeric(closed_df[[vc]])); if (vc == "log_centsize") vv[!cal_ok_s04] <- NA_real_; closed_df[[vc]] <- vv
    v_fin <- which(is.finite(vv)); v_ref <- intersect(v_fin, ref_std_idx); if (length(v_ref) < 3L) v_ref <- v_fin
    min_req <- if (vc %in% caul_het_raw) max(5L, ceiling(0.50 * length(ref_std_idx))) else 5L
    if (length(v_ref) >= min_req && is.finite(stats::sd(vv[v_ref])) && stats::sd(vv[v_ref]) > 1e-6) {
      closed_df[[paste0(vc, "_std")]] <- NA_real_; closed_df[[paste0(vc, "_std")]][v_fin] <- round((vv[v_fin] - mean(vv[v_ref])) / stats::sd(vv[v_ref]), 5)
    }
  }
  if (!isTRUE(opts$outline_only)) {
    dis_cands <- filter_biologically_variable_dissection_scalars(closed_df, intersect(c("solidity_std", "q25_solidity_std", "dissection_index_std", "within_rosette_efa_dispersion_std", "ECT_PC1_std", "ECT_PC2_std", "ECT_PC3_std", "ect_lobe_complexity_std", paste0(caul_het_raw, "_std")), names(closed_df)), active_mask = (closed_df$authority_weight > 0))
    kept_dis <- prune_collinear_dissection_scalars(closed_df, base_cols = feature_cols, cand_scalars = dis_cands, active_mask = (closed_df$authority_weight > 0), group_vec = NULL, max_r = 0.75)
    if (length(kept_dis) > 0L) feature_cols <- unique(c(feature_cols, kept_dis))
  }
  cand_taxa <- sort(unique(closed_df$species_standardized[!is.na(closed_df$species_standardized) & grepl("^Packera\\s+", closed_df$species_standardized)]))
  tx_n <- table(factor(closed_df$species_standardized, levels = cand_taxa)); target_taxa <- if (sum(tx_n >= 2L) >= 2L) names(tx_n[tx_n >= 2L]) else if (length(cand_taxa) > 0L) cand_taxa else sort(unique(stats::na.omit(closed_df$species_standardized)))
  is_intra_mode <- isTRUE(opts$intraspecific_mode) || (length(target_taxa) < 2L)
  if (exists("merge_voucher_phenophase_covariates", mode = "function")) closed_df <- merge_voucher_phenophase_covariates(closed_df, opts$harmonics, opts$vouchers)
  std_plast_cols <- setdiff(grep("_std$", names(closed_df), value = TRUE), c("log_centsize_std", "petiole_length_mm_std", "plant_height_cm_std"))
  if (exists("residualize_ontogenetic_phenophase_plasticity", mode = "function")) closed_df <- residualize_ontogenetic_phenophase_plasticity(closed_df, target_cols = unique(c(if (is_intra_mode) pca_cols else character(0), std_plast_cols)), ref_mask = (closed_df$authority_weight > 0) & (is.na(closed_df$morph_tier_source) | closed_df$morph_tier_source != "Cauline_Only_Fallback"), intraspecific_mode = is_intra_mode)
  gmm_res <- run_gmm_cluster_analysis(closed_df, feature_cols, max_k = opts$max_k, use_clustvarsel = isTRUE(opts$use_clustvarsel), intraspecific_mode = is_intra_mode, seed = opts$seed)
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE); tmp_bic <- paste0(opts$output_report, ".", Sys.getpid(), ".tmp"); write.csv(gmm_res$bic_table, file = tmp_bic, row.names = FALSE); file.rename(tmp_bic, opts$output_report)
  veg_cda <- if (!isTRUE(opts$outline_only)) prune_collinear_dissection_scalars(closed_df, base_cols = feature_cols, cand_scalars = intersect(c("aspect_ratio_std"), names(closed_df)), active_mask = (closed_df$authority_weight > 0), group_vec = NULL, max_r = 0.75) else character(0)
  cda_cols <- unique(c(feature_cols, veg_cda))
  cda_res <- run_canonical_discriminant_analysis(closed_df, cda_cols, target_taxa, gmm_res = gmm_res, intraspecific_mode = is_intra_mode)
  if (exists(".project_shared_topology_rows_step04", mode = "function")) { tp_proj <- .project_shared_topology_rows_step04(closed_df, cda_res, gmm_res, pca_cols); cda_res <- tp_proj$cda_res; gmm_res <- tp_proj$gmm_res }
  selected_features <- unique(c(cda_res$selected_features, gmm_res$selected_features)); if (length(selected_features) == 0L) selected_features <- feature_cols
  audited_df <- audit_misidentifications(closed_df, cda_res, gmm_res, selected_features = selected_features)
  audit_institutional_batch_effects(audited_df, output_csv = opts$output_batch_audit, selected_features = selected_features)
  rep_dir <- dirname(opts$output_report); eco_sum_csv <- file.path(rep_dir, "ecoregion_supervised_cda_summary.csv")
  eco_cda_res <- if (exists("run_supervised_ecoregion_cda", mode = "function")) run_supervised_ecoregion_cda(audited_df, cda_cols, gmm_res = gmm_res, output_summary_csv = eco_sum_csv, permutations = opts$permutations, seed = opts$seed) else NULL
  if (!is.null(eco_cda_res)) { audited_df$regional_group <- eco_cda_res$ecoregion_vec; audited_df$ecoregion_cda_predicted <- eco_cda_res$ecoregion_cda_predicted; audited_df$ecoregion_cda_posterior_prob <- eco_cda_res$ecoregion_cda_posterior_prob; audited_df$ecoregion_can1 <- eco_cda_res$ecoregion_can1; audited_df$ecoregion_can2 <- eco_cda_res$ecoregion_can2; if (is_intra_mode && all(is.na(audited_df$loocv_accuracy)) && !is.null(eco_cda_res$loocv_accuracy)) audited_df$loocv_accuracy <- eco_cda_res$loocv_accuracy }
  manova_res <- if (exists("run_multivariate_manova_audit", mode = "function")) run_multivariate_manova_audit(audited_df, cda_cols, gmm_res, output_group_csv = file.path(rep_dir, "manova_group_significance.csv"), output_pairwise_csv = file.path(rep_dir, "manova_pairwise_posthoc.csv"), include_gmm_clusters = !is_intra_mode) else NULL
  lead_cols <- intersect(c("catalogNumber", "plant_individual_id", "leaf_id", "basal_leaf_count", "cauline_leaf_count", "morph_tier_source", "pixels_per_mm", "is_calibrated_mm", "species_raw", "species_standardized", "identifiedBy", "determiner_tier", "type_status", "is_type_specimen", "is_promoted_bronze", "has_empirical_type_anchor", "authority_weight", "spatial_balance_weight", "population_deme_id", "deme_n_specimens", "cda_predicted_taxon", "cda_posterior_prob", "cda_second_taxon", "cda_second_posterior_prob", "cda_posterior_margin", "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance", "is_typical_morphotype", "min_congener_mahal_d2", "min_congener_taxon", "open_set_novelty_p", "is_confident_label_noise", "loocv_accuracy", "selected_morph_features", "is_passive", "morphometric_outlier_flag", "misidentification_flag", "triage_priority", "discordance_reason", "gmm_cluster", "gmm_uncertainty", "can1", "can2", "can3", "ecoregion_cda_predicted", "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2", pc_cols, "Cauline_PC1", "Cauline_PC2", "Cauline_PC3", "ECT_PC1", "ECT_PC2", "ECT_PC3", "Cauline_ECT_PC1", "Cauline_ECT_PC2", "ect_lobe_complexity", "cauline_ect_lobe_complexity", "heteroblasty_ect_distance", "PC1_raw_individual", "PC2_raw_individual", "centroid_size", "log_centsize", "solidity", "q25_solidity", "dissection_index", "dissected_leaf_fraction", "within_rosette_efa_dispersion", "heteroblasty_efa_distance", "cauline_solidity_med", "cauline_dissection_index_med", "cauline_aspect_ratio_med", "cauline_blade_length_mm", "basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "petiole_length_mm", "petiole_length_mm_std", "blade_length_mm", "petiole_to_blade_ratio", "petiole_to_blade_ratio_std", "AR_inv", "AR_inv_std", "plant_height_cm", "plant_height_cm_std", "scape_to_basal_leaf_ratio", "scape_to_basal_leaf_ratio_std", "is_bent_stem", "reproductive_trait_status", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm"), names(audited_df))
  out_df <- audited_df[, c(lead_cols, setdiff(names(audited_df), lead_cols))]; dir.create(dirname(opts$output_flags), recursive = TRUE, showWarnings = FALSE); tmp_flags <- paste0(opts$output_flags, ".", Sys.getpid(), ".tmp"); .write_csv_s04(out_df, tmp_flags, na = ""); file.rename(tmp_flags, opts$output_flags)
  if (exists("persist_step04_to_sqlite", mode = "function")) tryCatch(persist_step04_to_sqlite(out_df, opts$vouchers, focal_taxon = opts$focal_taxon, intraspecific_mode = is_intra_mode), error = function(e) warning("SQLite Step 04 persistence note: ", e$message))
  generate_cda_biplot_pdf(audited_df, cda_res, gmm_res, opts$output_plot, selected_features, manova_res = manova_res)
  tryCatch({ if (exists("generate_step04_pub_figures", mode = "function")) generate_step04_pub_figures(audited_df, cda_res, gmm_res, NULL, opts$output_plot, eco_cda_res = eco_cda_res, manova_res = manova_res, batch_csv = opts$output_batch_audit) }, error = function(e) warning("Step 04 publication figures export warning: ", e$message))
  return(invisible(out_df))
}
if (sys.nframe() == 0) run_gmm_morphotools_pipeline(parse_args_robust())
