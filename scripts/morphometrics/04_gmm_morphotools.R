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
    optparse::make_option(c("-v", "--vouchers"), type = "character", default = "data/tables/curated_vouchers.csv", help = "Vouchers metadata CSV"),
    optparse::make_option(c("-e", "--harmonics"), type = "character", default = "data/tables/leaf_efa_harmonics.csv", help = "Leaf EFA harmonics CSV"),
    optparse::make_option(c("-m", "--reproductive"), type = "character", default = "data/tables/voucher_reproductive_metrics.csv", help = "Voucher reproductive metrics CSV"),
    optparse::make_option(c("--outline-only"), action = "store_true", default = FALSE, help = "Evaluate vegetative leaf outlines only"),
    optparse::make_option(c("--correct-allometry"), type = "character", default = "auto", help = "Foliar allometry correction [auto, true, false]"),
    optparse::make_option(c("--allometry-report"), type = "character", default = "outputs/reports/allometry_significance_audit.csv", help = "Allometry test audit CSV"),
    optparse::make_option(c("--allometry-pairwise"), type = "character", default = "outputs/reports/allometry_pairwise_trajectories.csv", help = "Pairwise allometric trajectories CSV"),
    optparse::make_option(c("--allometry-plot"), type = "character", default = "outputs/figures/shape_allometry_regression.pdf", help = "Allometry regression plot PDF"),
    optparse::make_option(c("--use-clustvarsel"), action = "store_true", default = FALSE, help = "Use clustvarsel feature selection before Mclust"),
    optparse::make_option(c("-f", "--output-flags"), type = "character", default = "data/tables/morphometric_flags.csv", help = "Flags output CSV"),
    optparse::make_option(c("-p", "--output-plot"), type = "character", default = "outputs/figures/cda_passive_projection.pdf", help = "Output CDA PDF"),
    optparse::make_option(c("-r", "--output-report"), type = "character", default = "outputs/reports/gmm_bayes_factors_summary.csv", help = "BIC report CSV"),
    optparse::make_option(c("-b", "--output-batch-audit"), type = "character", default = "outputs/reports/institutional_batch_effect_audit.csv", help = "Batch audit CSV"),
    optparse::make_option(c("-k", "--max-k"), type = "integer", default = 8, help = "Max mixture components for GMM"),
    optparse::make_option(c("-d", "--num-pcs"), type = "integer", default = 5, help = "Number of PCA dimensions"),
    optparse::make_option(c("--permutations"), type = "integer", default = 999, help = "PERMANOVA permutations"),
    optparse::make_option(c("--seed"), type = "integer", default = 42, help = "Random seed"),
    optparse::make_option(c("--intraspecific-mode"), action = "store_true", default = FALSE, help = "Enable single-taxon intraspecific mode"),
    optparse::make_option(c("--focal-taxon"), type = "character", default = "", help = "Focal taxon name")
  )
  if (requireNamespace("optparse", quietly = TRUE)) {
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list))
    names(opts) <- gsub("-", "_", names(opts)); return(opts)
  }
  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(vouchers = "data/tables/curated_vouchers.csv", harmonics = "data/tables/leaf_efa_harmonics.csv", reproductive = "data/tables/voucher_reproductive_metrics.csv",
               outline_only = FALSE, correct_allometry = "auto", allometry_report = "outputs/reports/allometry_significance_audit.csv",
               allometry_pairwise = "outputs/reports/allometry_pairwise_trajectories.csv", allometry_plot = "outputs/figures/shape_allometry_regression.pdf",
               use_clustvarsel = FALSE, output_flags = "data/tables/morphometric_flags.csv", output_plot = "outputs/figures/cda_passive_projection.pdf",
               output_report = "outputs/reports/gmm_bayes_factors_summary.csv", output_batch_audit = "outputs/reports/institutional_batch_effect_audit.csv",
               max_k = 8, num_pcs = 5, permutations = 999, seed = 42, intraspecific_mode = FALSE, focal_taxon = "")
  i <- 1
  while (i <= length(raw_args)) {
    a <- raw_args[i]
    if (a %in% c("-v", "--vouchers") && i < length(raw_args)) { opts$vouchers <- raw_args[i+1]; i <- i+2 } else if (a %in% c("-e", "--harmonics") && i < length(raw_args)) { opts$harmonics <- raw_args[i+1]; i <- i+2
    } else if (a %in% c("-m", "--reproductive") && i < length(raw_args)) { opts$reproductive <- raw_args[i+1]; i <- i+2 } else if (a == "--outline-only") { opts$outline_only <- TRUE; i <- i+1
    } else if (a == "--use-clustvarsel") { opts$use_clustvarsel <- TRUE; i <- i+1 } else if (a == "--intraspecific-mode") { opts$intraspecific_mode <- TRUE; i <- i+1
    } else if (a == "--focal-taxon" && i < length(raw_args)) { opts$focal_taxon <- raw_args[i+1]; i <- i+2 } else if (a == "--correct-allometry" && i < length(raw_args)) { opts$correct_allometry <- raw_args[i+1]; i <- i+2
    } else if (a == "--allometry-report" && i < length(raw_args)) { opts$allometry_report <- raw_args[i+1]; i <- i+2 } else if (a == "--allometry-pairwise" && i < length(raw_args)) { opts$allometry_pairwise <- raw_args[i+1]; i <- i+2
    } else if (a == "--allometry-plot" && i < length(raw_args)) { opts$allometry_plot <- raw_args[i+1]; i <- i+2 } else if (a %in% c("-f", "--output-flags") && i < length(raw_args)) { opts$output_flags <- raw_args[i+1]; i <- i+2
    } else if (a %in% c("-p", "--output-plot") && i < length(raw_args)) { opts$output_plot <- raw_args[i+1]; i <- i+2 } else if (a %in% c("-r", "--output-report") && i < length(raw_args)) { opts$output_report <- raw_args[i+1]; i <- i+2
    } else if (a %in% c("-b", "--output-batch-audit") && i < length(raw_args)) { opts$output_batch_audit <- raw_args[i+1]; i <- i+2 } else if (a %in% c("-k", "--max-k") && i < length(raw_args)) { opts$max_k <- as.integer(raw_args[i+1]); i <- i+2
    } else if (a %in% c("-d", "--num-pcs") && i < length(raw_args)) { opts$num_pcs <- as.integer(raw_args[i+1]); i <- i+2 } else if (a == "--permutations" && i < length(raw_args)) { opts$permutations <- as.integer(raw_args[i+1]); i <- i+2
    } else if (a == "--seed" && i < length(raw_args)) { opts$seed <- as.integer(raw_args[i+1]); i <- i+2 } else { i <- i+1 }
  }
  return(opts)
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

standardize_packera_taxon <- function(species_vec) {
  sapply(species_vec, function(s) {
    if (is.na(s) || nchar(trimws(s)) == 0) return("Unknown")
    s_clean <- trimws(s)
    if (grepl("anonym|smallii|earlei", s_clean, ignore.case = TRUE)) return("Packera anonyma")
    if (grepl("tomentos|dubia", s_clean, ignore.case = TRUE)) return("Packera dubia")
    if (grepl("plattensis|flavovirens", s_clean, ignore.case = TRUE)) return("Packera plattensis")
    if (grepl("paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", s_clean, ignore.case = TRUE)) return("Packera paupercula")
    return(trimws(strsplit(s_clean, "\\(")[[1]][1]))
  }, USE.NAMES = FALSE)
}

# 2b. Foliar Allometry Multivariate Testing & Residual Correction
generate_allometry_diagnostic_plot <- function(sub_df, r2_size, p_size, r2_int, output_pdf) {
  if (!exists("render_allometry_diagnostic_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) {
    source("scripts/morphometrics/pub_figures_utils.R")
  }
  if (exists("render_allometry_diagnostic_pdf", mode = "function")) {
    render_allometry_diagnostic_pdf(sub_df, r2_size, p_size, r2_int, output_pdf)
  }
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
  sub_df$authority_weight <- if ("authority_weight" %in% names(sub_df)) suppressWarnings(as.numeric(sub_df$authority_weight)) else ifelse(sub_df$determiner_tier == "Tier_1_Gold", 2.0, ifelse(sub_df$determiner_tier == "Tier_2_Silver", 1.25, 0.75))
  sub_df$authority_weight[is.na(sub_df$authority_weight) | sub_df$authority_weight <= 0] <- 0.75
  gold_mask <- if ("determiner_tier" %in% names(sub_df)) (sub_df$determiner_tier %in% c("Tier_1_Gold", "Tier_2_Silver") & !is.na(sub_df$scientificName)) else rep(TRUE, nrow(sub_df))
  gold_df <- if (sum(gold_mask) >= 10 && length(unique(sub_df$scientificName[gold_mask])) >= 2) sub_df[gold_mask, , drop = FALSE] else sub_df
  gold_df$scientificName <- droplevels(factor(gold_df$scientificName))
  Y <- as.matrix(gold_df[, sym_names, drop = FALSE]); set.seed(seed); anova_df <- NULL; pw_df <- NULL

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
          pw_df <- data.frame(Pair = rownames(vc_s), r = round(vc_s[[1]], 4), angle_deg = round(vc_s[[2]], 3),
                              Z_angle = round(vc_s[["Z"]], 4), p_angle = round(vc_s[[ncol(vc_s)]], 5), stringsAsFactors = FALSE)
          if (!is.null(dt_s) && nrow(dt_s) == nrow(vc_s)) { pw_df$slope_dist <- round(dt_s[[1]], 4); pw_df$p_dist <- round(dt_s[[ncol(dt_s)]], 5) }
        }
      }
    }
  }
  if (is.null(anova_df) && requireNamespace("vegan", quietly = TRUE)) {
    allometry_fit <- tryCatch(vegan::adonis2(Y ~ log_centsize * scientificName, data = gold_df, permutations = permutations, method = "euclidean", by = "terms"), error = function(e) NULL)
    if (!is.null(allometry_fit)) {
      anova_df <- as.data.frame(allometry_fit); anova_df$Term <- rownames(anova_df); rownames(anova_df) <- NULL
      if ("F.Model" %in% names(anova_df)) anova_df$F <- anova_df$F.Model; if ("Pr(>F)" %in% names(anova_df)) anova_df$p_value <- anova_df[["Pr(>F)"]]
    }
  }
  if (is.null(anova_df)) {
    Y_cent <- scale(Y, center = TRUE, scale = FALSE); SS_tot <- sum(Y_cent^2); N_obs <- nrow(gold_df)
    fit1 <- stats::lm(Y_cent ~ log_centsize, data = gold_df, weights = authority_weight); SS_size <- sum(stats::fitted(fit1)^2); df_size <- fit1$rank - 1
    fit2 <- stats::lm(Y_cent ~ log_centsize + scientificName, data = gold_df, weights = authority_weight); SS_12 <- sum(stats::fitted(fit2)^2); SS_taxa <- SS_12 - SS_size; df_taxa <- fit2$rank - fit1$rank
    fit3 <- stats::lm(Y_cent ~ log_centsize * scientificName, data = gold_df, weights = authority_weight); SS_123 <- sum(stats::fitted(fit3)^2); SS_int <- SS_123 - SS_12; df_int <- fit3$rank - fit2$rank
    SS_res <- SS_tot - SS_123; df_res <- N_obs - fit3$rank; MS_res <- SS_res / max(df_res, 1)
    F_size <- (SS_size / max(df_size, 1)) / MS_res; F_taxa <- (SS_taxa / max(df_taxa, 1)) / MS_res; F_int <- (SS_int / max(df_int, 1)) / MS_res
    anova_df <- data.frame(
      Term = c("log_centsize", "scientificName", "log_centsize:scientificName", "Residual", "Total"),
      Df = c(df_size, df_taxa, df_int, df_res, N_obs - 1), SumOfSqs = round(c(SS_size, SS_taxa, SS_int, SS_res, SS_tot), 4),
      R2 = round(c(SS_size / SS_tot, SS_taxa / SS_tot, SS_int / SS_tot, SS_res / SS_tot, 1.0), 6),
      F = c(round(F_size, 4), round(F_taxa, 4), round(F_int, 4), NA_real_, NA_real_),
      p_value = c(round(stats::pf(F_size, df_size, df_res, lower.tail = FALSE), 5), round(stats::pf(F_taxa, df_taxa, df_res, lower.tail = FALSE), 5),
                  round(stats::pf(F_int, df_int, df_res, lower.tail = FALSE), 5), NA_real_, NA_real_), stringsAsFactors = FALSE)
  }
  if (!is.null(output_pairwise) && nzchar(output_pairwise) && !is.null(pw_df)) {
    dir.create(dirname(output_pairwise), recursive = TRUE, showWarnings = FALSE)
    tmp_pw <- paste0(output_pairwise, ".tmp"); write.csv(pw_df, file = tmp_pw, row.names = FALSE); file.rename(tmp_pw, output_pairwise)
  }
  anova_df$Term <- gsub("^Residuals$", "Residual", anova_df$Term)
  if ("Total" %in% anova_df$Term) anova_df$R2[anova_df$Term == "Total"] <- 1.0
  anova_df <- anova_df[, intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(anova_df)), drop = FALSE]
  if (!is.null(output_report) && nzchar(output_report)) {
    dir.create(dirname(output_report), recursive = TRUE, showWarnings = FALSE)
    tmp_rep <- paste0(output_report, ".tmp"); write.csv(anova_df, file = tmp_rep, row.names = FALSE); file.rename(tmp_rep, output_report)
    message("Allometry significance audit report exported -> ", output_report)
  }
  size_row <- anova_df[grepl("^log_centsize$", anova_df$Term, ignore.case = TRUE), ]
  r2_size <- if (nrow(size_row) > 0) size_row$R2[1] else 0.0
  p_size <- if (nrow(size_row) > 0 && !is.na(size_row$p_value[1])) size_row$p_value[1] else 0.0
  int_row <- anova_df[grepl("log_centsize:scientificName|scientificName:log_centsize", anova_df$Term, ignore.case = TRUE), ]
  r2_int <- if (nrow(int_row) > 0) int_row$R2[1] else 0.0
  p_int <- if (nrow(int_row) > 0 && !is.na(int_row$p_value[1])) int_row$p_value[1] else 1.0
  correct_flag <- tolower(trimws(as.character(correct_mode))); slopes_homogeneous <- is.na(p_int) || (p_int > 0.05)
  should_correct <- (correct_flag %in% c("true", "1", "yes")) || (correct_flag == "auto" && !is.na(r2_size) && r2_size >= r2_threshold && slopes_homogeneous)
  authority_weight <- gold_df$authority_weight
  fit_within <- if (nlevels(gold_df$scientificName) >= 2) stats::lm(Y ~ log_centsize + scientificName, data = gold_df, weights = authority_weight) else stats::lm(Y ~ log_centsize, data = gold_df, weights = authority_weight)
  beta_within <- stats::coef(fit_within)["log_centsize", ]; beta_within[is.na(beta_within)] <- 0.0
  Y <- as.matrix(sub_df[, sym_names, drop = FALSE]); log_centsize <- sub_df$log_centsize
  Y_adj <- Y - outer(log_centsize - mean(log_centsize, na.rm = TRUE), beta_within)
  active_cols <- which(apply(Y_adj, 2, stats::var, na.rm = TRUE) > 1e-8)
  pca_resids <- stats::prcomp(Y_adj[, active_cols, drop = FALSE], center = TRUE, scale. = FALSE)
  sub_df$PC1_raw <- sub_df$PC1; sub_df$PC2_raw <- sub_df$PC2
  sub_df$PC1_corrected <- round(pca_resids$x[, 1], 6); sub_df$PC2_corrected <- round(pca_resids$x[, 2], 6)
  if (should_correct) {
    message(sprintf("[INFO] Foliar allometry detected (R2 = %.2f, p = %.2f). Analyzing allometry-free shape residuals.", r2_size, p_size))
    morpho_matrix_for_gmm <- Y_adj; morpho_df[valid_idx, sym_names] <- round(Y_adj, 6)
    all_sym_idx <- which(complete.cases(morpho_df[, sym_names]))
    proj_pcs <- scale(as.matrix(morpho_df[all_sym_idx, sym_names[active_cols], drop = FALSE]), center = pca_resids$center, scale = FALSE) %*% pca_resids$rotation
    for (p in seq_len(min(5, ncol(proj_pcs)))) morpho_df[all_sym_idx, paste0("PC", p)] <- round(proj_pcs[, p], 6)
  } else {
    message(sprintf("[INFO] Foliar allometry negligible (R2 = %.2f). Analyzing raw symmetric harmonics.", r2_size))
    morpho_matrix_for_gmm <- Y
  }
  if (!is.null(output_plot) && nzchar(output_plot)) generate_allometry_diagnostic_plot(sub_df, r2_size, p_size, r2_int, output_plot)
  return(list(morpho_matrix = morpho_matrix_for_gmm, corrected = should_correct, beta_within = beta_within,
              r2_size = r2_size, p_size = p_size, r2_interaction = r2_int, p_interaction = p_int, anova_table = anova_df, pairwise_table = pw_df, df = morpho_df))
}
evaluate_ontogenetic_allometry <- test_and_correct_allometry

# ------------------------------------------------------------------------------
# 3. Gaussian Mixture Modeling & Bayes Factor Testing (Seed = 42)
# ------------------------------------------------------------------------------
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150) {
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
    if (bic_val > best_bic) best_bic <- bic_val; if (bic_val == best_bic) best_model <- list(k = k, classification = apply(resp, 1, which.max), uncertainty = 1 - apply(resp, 1, max), bic = bic_val)
  }
  return(list(best = best_model, bic_table = bic_table))
}

run_gmm_cluster_analysis <- function(df, feature_cols, max_k = 8, use_clustvarsel = FALSE) {
  valid_rows <- which(complete.cases(df[, feature_cols])); N_valid <- length(valid_rows)
  if (N_valid < 3) {
    message("Warning: Too few complete observations for GMM clustering.")
    return(list(classification = rep(NA_integer_, nrow(df)), uncertainty = rep(NA_real_, nrow(df)), best_k = 1, model_name = "None",
                selected_features = feature_cols, bic_table = data.frame(K = seq_len(max_k), BIC = NA_real_, Two_Delta_BIC = NA_real_, Evidence = "Inconclusive")))
  }
  X <- as.matrix(df[valid_rows, feature_cols, drop = FALSE]); selected_cols <- colnames(X)
  message(sprintf("Fitting Label-Blind Gaussian Mixture Models (K = 1 to %d, N = %d, seed = 42)...", max_k, N_valid))
  set.seed(42); best_k <- 1; best_name <- "Full_Covariance_EM"; bic_table <- NULL; cls_sub <- NULL; unc_sub <- NULL
  pc_guard_cols <- intersect(paste0("PC", seq_len(5)), colnames(X))
  X_mah <- X[, if (length(pc_guard_cols) > 0) pc_guard_cols else colnames(X), drop = FALSE]
  rob_cov <- if (nrow(X) >= 30 && requireNamespace("MASS", quietly = TRUE)) tryCatch(MASS::cov.rob(X_mah), error = function(e) NULL) else NULL
  m_center <- if (!is.null(rob_cov)) rob_cov$center else colMeans(X_mah)
  m_cov <- if (!is.null(rob_cov)) rob_cov$cov else stats::cov(X_mah)
  mahal_d2 <- tryCatch(stats::mahalanobis(X_mah, center = m_center, cov = m_cov + diag(1e-6, ncol(X_mah))), error = function(e) rep(0.0, nrow(X)))
  mahal_outliers <- (nrow(X) >= 30) & !is.na(mahal_d2) & (mahal_d2 > stats::qchisq(0.995, df = ncol(X_mah)))
  anchor_mask <- if ("determiner_tier" %in% names(df)) (df$determiner_tier[valid_rows] %in% c("Tier_1_Gold", "Tier_2_Silver")) & !mahal_outliers else !mahal_outliers
  core_idx <- if (sum(anchor_mask, na.rm = TRUE) >= 20) which(anchor_mask) else if (any(mahal_outliers) && sum(!mahal_outliers) >= max(10, ncol(X) + 2)) which(!mahal_outliers) else seq_len(nrow(X))
  if (requireNamespace("clustvarsel", quietly = TRUE) && (isTRUE(use_clustvarsel) || ncol(X) >= 4)) {
    cvs <- tryCatch(clustvarsel::clustvarsel(X[core_idx, , drop = FALSE], G = seq_len(max_k), verbose = FALSE), error = function(e) NULL)
    if (!is.null(cvs) && length(cvs$subset) >= 2) {
      selected_cols <- colnames(X)[cvs$subset]; X <- X[, cvs$subset, drop = FALSE]
      message("clustvarsel selected informative subset: ", paste(selected_cols, collapse = ", "))
    }
  }
  X_pcs <- X; X_core <- X_pcs[core_idx, , drop = FALSE]; n_samples <- nrow(X_pcs); eff_k <- if (nrow(X_core) > 100 && max_k >= 7L) 6L else max_k
  modelNames <- if ((nrow(X_core) / ncol(X_pcs)) < 5.0) c("EII", "VII", "EEI", "VEI", "EVI", "VVI") else NULL
  if (requireNamespace("mclust", quietly = TRUE)) {
    gmm_fit <- tryCatch(mclust::Mclust(X_core, G = seq_len(eff_k), modelNames = modelNames, verbose = FALSE), error = function(e) NULL)
    if (!is.null(gmm_fit)) {
      best_k <- gmm_fit$G; best_name <- gmm_fit$modelName
      gmm_pred <- predict(gmm_fit, newdata = X_pcs)
      cls_sub <- as.integer(gmm_pred$classification)
      unc_sub <- if (!is.null(gmm_pred$uncertainty)) gmm_pred$uncertainty else 1.0 - apply(as.matrix(gmm_pred$z), 1, max)
      bic_vals <- apply(gmm_fit$BIC, 1, function(row) if (all(is.na(row))) NA_real_ else max(row, na.rm = TRUE))
      bic_table <- data.frame(K = seq_len(eff_k), BIC = round(bic_vals[seq_len(eff_k)], 2))
    }
  }
  if (is.null(bic_table)) {
    pure_res <- fit_gmm_em_pure(X, max_k = max_k); best_k <- pure_res$best$k
    cls_sub <- pure_res$best$classification; unc_sub <- pure_res$best$uncertainty; bic_table <- pure_res$bic_table
  }
  cls <- rep(NA_integer_, nrow(df)); unc <- rep(NA_real_, nrow(df)); cls[valid_rows] <- cls_sub; unc[valid_rows] <- unc_sub
  bic_table$Delta_BIC_prev <- c(0, diff(bic_table$BIC)); bic_table$Two_Delta_BIC <- round(bic_table$Delta_BIC_prev, 2)
  bic_table$BayesFactor_vs_Null <- round(bic_table$BIC - bic_table$BIC[1], 2)
  bic_table$Evidence <- sapply(bic_table$Two_Delta_BIC, function(val) {
    if (is.na(val) || val < 0) return("Negative (favors K-1)")
    if (val < 2) return("Weak / Inconclusive"); if (val < 6) return("Positive Evidence"); if (val < 10) return("Strong Evidence")
    return("Decisive Evidence (Species Boundary)")
  })
  return(list(classification = cls, uncertainty = unc, best_k = best_k, model_name = best_name, selected_features = selected_cols, bic_table = bic_table))
}

# 4. Canonical Discriminant Analysis with Passive Sample Projection
run_canonical_discriminant_analysis <- function(df, feature_cols, target_taxa, gmm_res = NULL, intraspecific_mode = FALSE) {
  message("Configuring Canonical Discriminant Analysis in MorphoTools2 Architecture...")
  df$authority_weight <- ifelse(df$determiner_tier == "Tier_1_Gold", 2.0, ifelse(df$determiner_tier == "Tier_2_Silver", 1.25, 0.75))
  df$authority_weight[is.na(df$authority_weight)] <- 0.75
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2)
  if (is_intra) {
    if (!is.null(gmm_res) && !is.null(gmm_res$best_k) && gmm_res$best_k >= 2 && length(unique(stats::na.omit(gmm_res$classification))) >= 2) {
      df$grouping_col <- ifelse(is.na(gmm_res$classification), NA_character_, paste0("Morphotype_K", gmm_res$classification))
      u_cut <- stats::quantile(gmm_res$uncertainty, 0.60, na.rm = TRUE)
      is_active <- !is.na(df$grouping_col) & !is.na(gmm_res$uncertainty) & (gmm_res$uncertainty <= u_cut) & complete.cases(df[, feature_cols])
      if (length(unique(df$grouping_col[is_active])) < 2) is_active <- !is.na(df$grouping_col) & complete.cases(df[, feature_cols])
    } else {
      rg <- if ("regional_group" %in% names(df)) trimws(as.character(df$regional_group)) else rep("", nrow(df))
      st <- if ("stateProvince" %in% names(df)) trimws(as.character(df$stateProvince)) else rep("", nrow(df))
      if (length(unique(rg[!is.na(rg) & nzchar(rg)])) >= 2) df$grouping_col <- ifelse(is.na(rg) | !nzchar(rg), "Other_Region", rg)
      else if (length(unique(st[!is.na(st) & nzchar(st)])) >= 2) df$grouping_col <- ifelse(is.na(st) | !nzchar(st), "Other_State", st)
      else if ("latitude" %in% names(df) && length(unique(stats::na.omit(df$latitude))) >= 3) {
        q_l <- unique(stats::quantile(df$latitude, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
        df$grouping_col <- paste0("LatTertile_", as.integer(cut(df$latitude, breaks = q_l, include.lowest = TRUE)))
      } else {
        p1 <- if ("PC1" %in% names(df)) df$PC1 else df[[feature_cols[1]]]
        df$grouping_col <- ifelse(p1 >= stats::median(p1, na.rm = TRUE), "Morphotype_K2", "Morphotype_K1")
      }
      g_act <- (df$determiner_tier == "Tier_1_Gold") & !is.na(df$grouping_col) & complete.cases(df[, feature_cols])
      is_active <- if (length(unique(df$grouping_col[g_act])) >= 2) g_act else (!is.na(df$grouping_col) & complete.cases(df[, feature_cols]))
    }
    target_taxa <- sort(unique(stats::na.omit(df$grouping_col[is_active])))
    is_train_cand <- is_active
  } else {
    df$grouping_col <- df$species_standardized
    is_active <- (df$determiner_tier == "Tier_1_Gold") & (df$grouping_col %in% target_taxa) & complete.cases(df[, feature_cols])
    is_train_cand <- df$determiner_tier %in% c("Tier_1_Gold", "Tier_2_Silver") & df$species_standardized %in% target_taxa & complete.cases(df[, feature_cols])
    if (length(unique(df$grouping_col[is_train_cand])) >= 2) target_taxa <- intersect(target_taxa, unique(df$grouping_col[is_train_cand]))
    if (any(table(factor(df$grouping_col[is_train_cand], levels = target_taxa)) < 2)) {
      is_train_cand <- if (all(table(factor(df$grouping_col[is_active], levels = target_taxa)) >= 2)) is_active else (df$grouping_col %in% target_taxa & complete.cases(df[, feature_cols]))
    }
  }
  X <- as.matrix(df[, feature_cols, drop = FALSE]); p <- ncol(X)
  for (grp in target_taxa) {
    cand_k <- which(is_train_cand & df$grouping_col == grp); n_k <- length(cand_k)
    if (n_k >= 6) {
      sub_x <- X[cand_k, , drop = FALSE]
      d2_k <- tryCatch(stats::mahalanobis(sub_x, colMeans(sub_x), stats::cov(sub_x) + diag(1e-4, p)), error = function(e) rep(0, n_k))
      drop_k <- !is.na(d2_k) & (d2_k > stats::qchisq(0.95, df = p))
      if (sum(!drop_k) >= 2) is_train_cand[cand_k[drop_k]] <- FALSE
    }
  }
  active_idx <- which(is_train_cand); passive_idx <- which(!is_active)
  message(sprintf("Active Anchors (Complete): %d | Passive Projected: %d", length(active_idx), length(passive_idx)))
  if (length(active_idx) < length(target_taxa)) stop("Insufficient complete active Tier 1 Gold anchors for CDA.")

  X_a <- X[active_idx, , drop = FALSE]; y_a <- df$grouping_col[active_idx]; w_a <- df$authority_weight[active_idx]; N_a <- nrow(X_a)
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
    grp <- g_levels[k]; k_mask <- (y_a == grp); sub_x <- X_a[k_mask, , drop = FALSE]; w_k <- w_a[k_mask]; group_counts[k] <- nrow(sub_x)
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
  priors <- rep(1.0 / g, g); post_probs <- matrix(0, nrow(df), g, dimnames = list(NULL, g_levels))
  for (k in seq_along(g_levels)) {
    diff_z <- t(t(Z_all) - canonical_centroids[k, ]); post_probs[, k] <- priors[k] * exp(-0.5 * rowSums(diff_z^2))
  }
  post_probs <- post_probs / rowSums(post_probs); pred_idx <- apply(post_probs, 1, which.max)
  posterior_prob <- apply(post_probs, 1, max)
  second_posterior_prob <- if (g >= 2) apply(post_probs, 1, function(r) sort(r, decreasing = TRUE)[2]) else rep(0.0, nrow(df))
  posterior_margin <- posterior_prob - second_posterior_prob
  orth_against <- function(vec, basis_mat) {
    v <- as.numeric(vec) - mean(vec, na.rm = TRUE)
    for (b_col in seq_len(ncol(basis_mat))) {
      b <- as.numeric(basis_mat[, b_col]); v <- v - (sum(v * b, na.rm = TRUE) / max(sum(b^2, na.rm = TRUE), 1e-12)) * b
    }
    v
  }
  if (ncol(Z_all) == 1L) {
    src1 <- if ("PC1" %in% colnames(X_proj)) X_proj[, "PC1"] else X_proj[, 1]
    can2_v <- orth_against(src1, Z_all)
    if (stats::sd(can2_v, na.rm = TRUE) < 1e-8 && ncol(X_proj) >= 2) can2_v <- orth_against(X_proj[, 2], Z_all)
    Z_all <- cbind(Z_all, can2_v)
  }
  if (ncol(Z_all) == 2L) {
    src3 <- if ("PC3" %in% colnames(X_proj)) X_proj[, "PC3"] else if ("PC2" %in% colnames(X_proj)) X_proj[, "PC2"] else X_proj[, min(3L, ncol(X_proj))]
    can3_v <- orth_against(src3, Z_all)
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
              predicted_taxon = g_levels[pred_idx], posterior_prob = posterior_prob, second_posterior_prob = second_posterior_prob,
              posterior_margin = posterior_margin, is_active = is_active, structure_loadings = struct_corr, loocv_accuracy = loocv_accuracy, grouping_col = df$grouping_col))
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
    raw_sp <- grp_vec[i]; pred_sp <- cda_res$predicted_taxon[i]; post_p <- cda_res$posterior_prob[i]; is_act <- cda_res$is_active[i]
    if (is_act) {
      if (raw_sp == pred_sp) { flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Verified_Tier1_Gold_Anchor" }
      else if (post_p >= 0.85) { flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- sprintf("Tier_1_Gold_Morphological_Discordance_to_%s", gsub(" ", "_", pred_sp)) }
      else { flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Tier_1_Gold_Borderline_Variant" }
    } else {
      if (raw_sp == pred_sp && post_p >= 0.75) { flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Congruent_Bronze_Determination" }
      else if (raw_sp != pred_sp && post_p >= 0.70) {
        flags[i] <- TRUE; triage[i] <- "HIGH"
        reasons[i] <- if (grepl("dubia|tomentos", raw_sp, ignore.case = TRUE) && pred_sp == "Packera anonyma") "Glabrescent_P_dubia_Misidentified_as_P_anonyma" else sprintf("Tier_3_Bronze_Misidentified_as_%s", gsub(" ", "_", pred_sp))
      } else if (post_p < 0.60) { flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Morphological_Intermediate_Ambiguous" }
      else { flags[i] <- (raw_sp != pred_sp); triage[i] <- if (raw_sp != pred_sp) "MEDIUM" else "LOW"
             reasons[i] <- if (raw_sp != pred_sp) "Candidate_Discordance_Moderate_Posterior" else "Minor_Uncertainty" }
    }
    if (is_intra_audit && !flags[i] && isTRUE(mahal_out[i])) { flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- "Intraspecific_Mahalanobis_Morphospace_Outlier" }
  }
  df$authority_weight <- ifelse(df$determiner_tier == "Tier_1_Gold", 2.0, ifelse(df$determiner_tier == "Tier_2_Silver", 1.25, 0.75))
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

# 5b. Institutional Batch-Effect Audit (ANOVA & Variance Explained)
audit_institutional_batch_effects <- function(morpho_df, output_csv = "outputs/reports/institutional_batch_effect_audit.csv", min_obs = 15) {
  message("=== Institutional Batch-Effect Audit (Digitization Rig Diagnostic) ===")
  if (!("institutionCode" %in% names(morpho_df))) return(NULL)
  valid_df <- morpho_df[!is.na(morpho_df$institutionCode) & trimws(morpho_df$institutionCode) != "" & !is.na(morpho_df$PC1) & !is.na(morpho_df$PC2), ]
  if (nrow(valid_df) < 10) return(NULL)
  inst_counts <- table(valid_df$institutionCode); rep_insts <- names(inst_counts[inst_counts >= min_obs])
  audit_subset <- if (length(rep_insts) >= 2) valid_df[valid_df$institutionCode %in% rep_insts, ] else valid_df
  audit_subset$institutionCode <- droplevels(factor(audit_subset$institutionCode))
  if (nlevels(audit_subset$institutionCode) < 2) { message("  [INFO] Fewer than 2 distinct institutions represented; skipping batch-effect ANOVA."); return(NULL) }
  lm_pc1 <- lm(PC1 ~ institutionCode, data = audit_subset); lm_pc2 <- lm(PC2 ~ institutionCode, data = audit_subset)
  s1 <- summary(lm_pc1); s2 <- summary(lm_pc2); aov1 <- anova(lm_pc1); aov2 <- anova(lm_pc2); r2_pc1 <- s1$r.squared; r2_pc2 <- s2$r.squared
  if (r2_pc1 >= 0.05 || r2_pc2 >= 0.05) warning(sprintf("Institutional batch effect R^2 exceeds 0.05 threshold (PC1 R^2=%.4f, PC2 R^2=%.4f)", r2_pc1, r2_pc2))
  else message("  [PASS] Assertion satisfied: Institutional identity explains < 5% variance (R^2 < 0.05).")
  anova_df <- data.frame(
    Response = c("PC1", "PC1", "PC2", "PC2"), Term = c("institutionCode", "Residuals", "institutionCode", "Residuals"),
    Df = c(aov1$Df[1], aov1$Df[2], aov2$Df[1], aov2$Df[2]), Sum_Sq = round(c(aov1$`Sum Sq`[1], aov1$`Sum Sq`[2], aov2$`Sum Sq`[1], aov2$`Sum Sq`[2]), 4),
    Mean_Sq = round(c(aov1$`Mean Sq`[1], aov1$`Mean Sq`[2], aov2$`Mean Sq`[1], aov2$`Mean Sq`[2]), 4),
    F_value = c(round(aov1$`F value`[1], 4), NA, round(aov2$`F value`[1], 4), NA), P_value = c(signif(aov1$`Pr(>F)`[1], 4), NA, signif(aov2$`Pr(>F)`[1], 4), NA),
    R_squared = round(c(r2_pc1, NA, r2_pc2, NA), 4), Adj_R_squared = round(c(s1$adj.r.squared, NA, s2$adj.r.squared, NA), 4),
    Batch_Effect_Status = c(ifelse(r2_pc1 < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA, ifelse(r2_pc2 < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA), stringsAsFactors = FALSE)
  if (!is.null(output_csv) && nzchar(output_csv)) {
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE); tmp_ba <- paste0(output_csv, ".tmp"); write.csv(anova_df, file = tmp_ba, row.names = FALSE); file.rename(tmp_ba, output_csv)
  }
  return(invisible(list(lm_pc1 = lm_pc1, lm_pc2 = lm_pc2, anova_table = anova_df)))
}

# 6. Publication-Quality 3D CDA & 3D GMM Biplot Vector Graphic (PDF)
generate_cda_biplot_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL) {
  if (!exists("render_cda_biplot_4panel_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) {
    source("scripts/morphometrics/pub_figures_utils.R")
  }
  if (exists("render_cda_biplot_4panel_pdf", mode = "function")) {
    render_cda_biplot_4panel_pdf(audit_df, cda_res, gmm_res, output_pdf, feature_cols)
  }
  message("CDA biplot graphic successfully saved to ", output_pdf)
}

# 7. Main Workflow Orchestrator
run_gmm_morphotools_pipeline <- function(opts) {
  message("==================================================================")
  message("Starting Packera Morphometrics GMM & Passive MorphoTools2 CDA")
  message("Vouchers: ", opts$vouchers, " | Harmonics: ", opts$harmonics, " | Allometry Mode: ", opts$correct_allometry)
  message("==================================================================")
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
    message("Aggregating to specimen-level profile (1 row = 1 voucher collection event)...")
    pc_sym_cols <- intersect(c(pca_cols, paste0("PC", seq_len(5)), paste0("Sym_PC", seq_len(5))), names(closed_df))
    num_cols <- intersect(c(pc_sym_cols, "foliar_variance", "leaf_count", "aspect_ratio", "solidity", "ucs_score", "area_px", "centroid_size", "log_centsize", "harmonic_roughness", "leaf_quality_weight"), names(closed_df))
    other_num <- setdiff(num_cols, pc_sym_cols); has_w <- "leaf_quality_weight" %in% names(closed_df)
    meta_cols <- setdiff(names(closed_df), c(num_cols, "catalogNumber", "plant_individual_id", "shape_id", "leaf_id", "mask_source"))
    grp_list <- split(closed_df, list(closed_df$catalogNumber, closed_df$plant_individual_id), drop = TRUE)
    closed_df <- do.call(rbind, lapply(grp_list, function(grp) {
      leaf_quality_weight <- if (has_w) suppressWarnings(as.numeric(grp$leaf_quality_weight)) else rep(1.0, nrow(grp))
      leaf_quality_weight[is.na(leaf_quality_weight) | leaf_quality_weight <= 0] <- 1.0
      is_contour_outlier <- if ("is_contour_outlier" %in% names(grp)) as.logical(grp$is_contour_outlier) else rep(FALSE, nrow(grp)); is_contour_outlier[is.na(is_contour_outlier)] <- FALSE
      if (nrow(grp) >= 3) {
        pc_mat <- as.matrix(grp[, intersect(paste0("PC", seq_len(5)), names(grp)), drop = FALSE]); med_pc <- apply(pc_mat, 2, stats::median, na.rm = TRUE)
        dist <- sqrt(rowSums(sweep(pc_mat, 2, med_pc, "-")^2, na.rm = TRUE)); iqr_d <- stats::IQR(dist, na.rm = TRUE); med_d <- stats::median(dist, na.rm = TRUE)
        drop_leaf <- ((dist > 2.0 * iqr_d & iqr_d > 0) | (dist > 2.5 * med_d)) & (is_contour_outlier | nrow(grp) >= 3)
        if (sum(!drop_leaf) < 2) drop_leaf[utils::head(order(leaf_quality_weight, -dist, decreasing = TRUE), 2)] <- FALSE
        grp <- grp[!drop_leaf, , drop = FALSE]; leaf_quality_weight <- leaf_quality_weight[!drop_leaf]
      }
      out <- grp[1, , drop = FALSE]
      for (col in pc_sym_cols) out[[col]] <- if (has_w) stats::weighted.mean(grp[[col]], w = leaf_quality_weight, na.rm = TRUE) else mean(grp[[col]], na.rm = TRUE)
      for (col in other_num) out[[col]] <- stats::median(grp[[col]], na.rm = TRUE)
      if ("leaf_count" %in% names(out)) out$leaf_count <- nrow(grp)
      for (col in meta_cols) { v <- stats::na.omit(grp[[col]]); out[[col]] <- if (length(v) > 0) v[1] else grp[[col]][1] }
      out
    })); rownames(closed_df) <- NULL
  }
  message(sprintf("Specimen-level collection units for morphometric modeling (1 row = 1 voucher): %d", nrow(closed_df)))
  closed_df$species_standardized <- standardize_packera_taxon(closed_df$species_raw)
  closed_df$authority_weight <- ifelse(closed_df$determiner_tier == "Tier_1_Gold", 2.0, ifelse(closed_df$determiner_tier == "Tier_2_Silver", 1.25, 0.75))
  if (!is.null(vouchers_df)) {
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), , drop = FALSE]; v_idx <- match(as.character(closed_df$catalogNumber), as.character(v_sub$catalogNumber))
    for (col in intersect(names(v_sub), c("institutionCode", "county", "stateProvince", "latitude", "longitude", "doy", "pheno_sin", "pheno_cos", "regional_group"))) {
      if (!col %in% names(closed_df)) closed_df[[col]] <- v_sub[[col]][v_idx]
    }
  }
  allom_res <- test_and_correct_allometry(closed_df, correct_mode = opts$correct_allometry, permutations = opts$permutations,
                                          seed = opts$seed, output_report = opts$allometry_report, output_plot = opts$allometry_plot, output_pairwise = opts$allometry_pairwise)
  if (!is.null(allom_res$df)) closed_df <- allom_res$df
  include_repro <- FALSE; feature_cols <- pca_cols
  if (!isTRUE(opts$outline_only)) {
    if ("capitulum_aspect_ratio" %in% names(closed_df)) closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$capitulum_aspect_ratio))
    else if ("AR_inv" %in% names(closed_df)) closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$AR_inv))
    else if (!is.null(opts$reproductive) && file.exists(opts$reproductive)) {
      repro_df <- read.csv(opts$reproductive, stringsAsFactors = FALSE)
      if ("capitulum_aspect_ratio" %in% names(repro_df)) repro_df$AR_inv <- suppressWarnings(as.numeric(repro_df$capitulum_aspect_ratio))
      r_cols <- intersect(names(repro_df), c("catalogNumber", "AR_inv", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm"))
      closed_df <- merge(closed_df, repro_df[!duplicated(repro_df$catalogNumber), r_cols, drop = FALSE], by = "catalogNumber", all.x = TRUE)
    }
    if ("AR_inv" %in% names(closed_df)) {
      valid_repro <- which(!is.na(closed_df$AR_inv) & is.finite(closed_df$AR_inv))
      if (length(valid_repro) >= 5) {
        closed_df$AR_inv_std <- NA_real_; m_ar <- mean(closed_df$AR_inv[valid_repro]); s_ar <- stats::sd(closed_df$AR_inv[valid_repro])
        closed_df$AR_inv_std[valid_repro] <- if (is.finite(s_ar) && s_ar > 1e-6) round((closed_df$AR_inv[valid_repro] - m_ar) / s_ar, 5) else 0.0
        include_repro <- TRUE; feature_cols <- c(pca_cols, "AR_inv_std")
        message(sprintf("Integrated macro-reproductive trait (AR_inv): %d / %d vouchers with valid capitula (mean=%.3f, sd=%.3f)", length(valid_repro), nrow(closed_df), m_ar, s_ar))
      }
    }
  }
  target_taxa <- c("Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula")
  gmm_res <- run_gmm_cluster_analysis(closed_df, feature_cols, max_k = opts$max_k, use_clustvarsel = isTRUE(opts$use_clustvarsel))
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE)
  tmp_bic <- paste0(opts$output_report, ".tmp"); write.csv(gmm_res$bic_table, file = tmp_bic, row.names = FALSE); file.rename(tmp_bic, opts$output_report)
  message("Bayes Factor (2ΔBIC) summary report saved: ", opts$output_report)
  veg_extra <- intersect(c("log_centsize", "aspect_ratio", "solidity"), names(closed_df))
  for (vc in veg_extra) if (all(is.finite(closed_df[[vc]])) && stats::sd(closed_df[[vc]]) > 1e-6) closed_df[[paste0(vc, "_std")]] <- as.numeric(scale(closed_df[[vc]]))
  cda_cols <- c(feature_cols, if (!isTRUE(opts$outline_only)) intersect(paste0(veg_extra, "_std"), names(closed_df)) else character(0))
  is_intra_mode <- isTRUE(opts$intraspecific_mode) || (length(unique(closed_df$species_standardized)) < 2)
  cda_res <- run_cda_with_passive_projection(closed_df, cda_cols, target_taxa, gmm_res = gmm_res, intraspecific_mode = is_intra_mode)
  audited_df <- audit_misidentifications(closed_df, cda_res, gmm_res)
  audit_institutional_batch_effects(audited_df, output_csv = opts$output_batch_audit)
  lead_cols <- intersect(c("catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier", "pixels_per_mm", "is_calibrated_mm",
                           "species_raw", "species_standardized", "determiner_tier", "authority_weight", "cda_predicted_taxon", "cda_posterior_prob",
                           "cda_second_posterior_prob", "cda_posterior_margin", "loocv_accuracy", "is_passive", "morphometric_outlier_flag", "misidentification_flag",
                           "triage_priority", "discordance_reason", "gmm_cluster", "gmm_uncertainty", "can1", "can2", "can3", "pheno_sin", "pheno_cos",
                           "PC1", "PC2", "PC3", "PC4", "PC5", "centroid_size", "log_centsize", "solidity", "ucs_score", "AR_inv", "AR_inv_std",
                           "reproductive_trait_status", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm"), names(audited_df))
  out_df <- audited_df[, c(lead_cols, setdiff(names(audited_df), lead_cols))]
  dir.create(dirname(opts$output_flags), recursive = TRUE, showWarnings = FALSE)
  tmp_flags <- paste0(opts$output_flags, ".tmp"); write.csv(out_df, file = tmp_flags, row.names = FALSE, na = ""); file.rename(tmp_flags, opts$output_flags)
  message("Master morphometric flags table exported: ", opts$output_flags, " (Rows: ", nrow(out_df), ")")
  generate_cda_biplot_pdf(audited_df, cda_res, gmm_res, opts$output_plot, feature_cols)
  tryCatch({
    if (file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R"); generate_step04_pub_figures(audited_df, cda_res, gmm_res, NULL, opts$output_plot)
    }
  }, error = function(e) warning("Step 04 publication figures export warning: ", e$message))
  message("Morphometric Pipeline Complete Successfully."); return(invisible(out_df))
}

if (sys.nframe() == 0) run_gmm_morphotools_pipeline(parse_args_robust())
