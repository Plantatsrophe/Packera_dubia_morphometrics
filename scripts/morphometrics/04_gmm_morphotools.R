#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_gmm_morphotools.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("mclust", quietly = TRUE)) library(mclust)
  if (requireNamespace("MorphoTools2", quietly = TRUE)) library(MorphoTools2)
  if (requireNamespace("vegan", quietly = TRUE)) library(vegan)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("tibble", quietly = TRUE)) library(tibble)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})

# ------------------------------------------------------------------------------
# 1. CLI Argument Parsing
# ------------------------------------------------------------------------------
parse_args_robust <- function() {
  option_list <- list(
    optparse::make_option(c("-v", "--vouchers"), type = "character", default = "data/tables/curated_vouchers.csv", help = "Vouchers metadata CSV [default: %default]"),
    optparse::make_option(c("-e", "--harmonics"), type = "character", default = "data/tables/leaf_efa_harmonics.csv", help = "Leaf EFA harmonics CSV [default: %default]"),
    optparse::make_option(c("-m", "--reproductive"), type = "character", default = "data/tables/voucher_reproductive_metrics.csv", help = "Voucher reproductive metrics CSV [default: %default]"),
    optparse::make_option(c("--outline-only"), action = "store_true", default = FALSE, help = "Evaluate vegetative leaf outlines only [default: %default]"),
    optparse::make_option(c("--correct-allometry"), type = "character", default = "auto", help = "Foliar allometry correction [auto, true, false] [default: %default]"),
    optparse::make_option(c("--allometry-report"), type = "character", default = "outputs/reports/allometry_significance_audit.csv", help = "Allometry test audit CSV [default: %default]"),
    optparse::make_option(c("--allometry-plot"), type = "character", default = "outputs/figures/shape_allometry_regression.pdf", help = "Allometry regression plot PDF [default: %default]"),
    optparse::make_option(c("-f", "--output-flags"), type = "character", default = "data/tables/morphometric_flags.csv", help = "Flags output CSV [default: %default]"),
    optparse::make_option(c("-p", "--output-plot"), type = "character", default = "outputs/figures/cda_passive_projection.pdf", help = "Output CDA PDF [default: %default]"),
    optparse::make_option(c("-r", "--output-report"), type = "character", default = "outputs/reports/gmm_bayes_factors_summary.csv", help = "BIC report CSV [default: %default]"),
    optparse::make_option(c("-b", "--output-batch-audit"), type = "character", default = "outputs/reports/institutional_batch_effect_audit.csv", help = "Batch audit CSV [default: %default]"),
    optparse::make_option(c("-k", "--max-k"), type = "integer", default = 8, help = "Max mixture components for GMM [default: %default]"),
    optparse::make_option(c("-d", "--num-pcs"), type = "integer", default = 5, help = "Number of PCA dimensions [default: %default]"),
    optparse::make_option(c("--permutations"), type = "integer", default = 999, help = "PERMANOVA permutations [default: %default]"),
    optparse::make_option(c("--seed"), type = "integer", default = 42, help = "Random seed [default: %default]")
  )
  if (requireNamespace("optparse", quietly = TRUE)) {
    return(optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list)))
  }
  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(vouchers = "data/tables/curated_vouchers.csv", harmonics = "data/tables/leaf_efa_harmonics.csv",
               reproductive = "data/tables/voucher_reproductive_metrics.csv", outline_only = FALSE,
               correct_allometry = "auto", allometry_report = "outputs/reports/allometry_significance_audit.csv",
               allometry_plot = "outputs/figures/shape_allometry_regression.pdf", output_flags = "data/tables/morphometric_flags.csv",
               output_plot = "outputs/figures/cda_passive_projection.pdf", output_report = "outputs/reports/gmm_bayes_factors_summary.csv",
               output_batch_audit = "outputs/reports/institutional_batch_effect_audit.csv", max_k = 8, num_pcs = 5, permutations = 999, seed = 42)
  i <- 1
  while (i <= length(raw_args)) {
    a <- raw_args[i]
    if (a %in% c("-v", "--vouchers") && i < length(raw_args)) { opts$vouchers <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-e", "--harmonics") && i < length(raw_args)) { opts$harmonics <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-m", "--reproductive") && i < length(raw_args)) { opts$reproductive <- raw_args[i+1]; i <- i+2 }
    else if (a == "--outline-only") { opts$outline_only <- TRUE; i <- i+1 }
    else if (a == "--correct-allometry" && i < length(raw_args)) { opts$correct_allometry <- raw_args[i+1]; i <- i+2 }
    else if (a == "--allometry-report" && i < length(raw_args)) { opts$allometry_report <- raw_args[i+1]; i <- i+2 }
    else if (a == "--allometry-plot" && i < length(raw_args)) { opts$allometry_plot <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-f", "--output-flags") && i < length(raw_args)) { opts$output_flags <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-p", "--output-plot") && i < length(raw_args)) { opts$output_plot <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-r", "--output-report") && i < length(raw_args)) { opts$output_report <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-b", "--output-batch-audit") && i < length(raw_args)) { opts$output_batch_audit <- raw_args[i+1]; i <- i+2 }
    else if (a %in% c("-k", "--max-k") && i < length(raw_args)) { opts$max_k <- as.integer(raw_args[i+1]); i <- i+2 }
    else if (a %in% c("-d", "--num-pcs") && i < length(raw_args)) { opts$num_pcs <- as.integer(raw_args[i+1]); i <- i+2 }
    else if (a == "--permutations" && i < length(raw_args)) { opts$permutations <- as.integer(raw_args[i+1]); i <- i+2 }
    else if (a == "--seed" && i < length(raw_args)) { opts$seed <- as.integer(raw_args[i+1]); i <- i+2 }
    else { i <- i+1 }
  }
  return(opts)
}

# ------------------------------------------------------------------------------
# 2. Circular Phenology & Taxonomic Concept Standardization
# ------------------------------------------------------------------------------
calculate_circular_phenology <- function(df) {
  if ("doy" %in% names(df)) {
    doy_val <- suppressWarnings(as.numeric(df$doy)); valid <- !is.na(doy_val) & doy_val >= 1 & doy_val <= 366
    df$pheno_sin <- NA_real_; df$pheno_cos <- NA_real_
    df$pheno_sin[valid] <- round(sin(2 * pi * doy_val[valid] / 365.25), 6)
    df$pheno_cos[valid] <- round(cos(2 * pi * doy_val[valid] / 365.25), 6)
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

# ------------------------------------------------------------------------------
# 2b. Foliar Allometry Multivariate Testing & Residual Correction
# ------------------------------------------------------------------------------
generate_allometry_diagnostic_plot <- function(sub_df, r2_size, p_size, r2_int, output_pdf) {
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  taxon_cols <- c("Packera anonyma" = "#E69F00", "Packera dubia" = "#009E73", "Packera plattensis" = "#56B4E9", "Packera paupercula" = "#CC79A7")
  taxa_in_data <- intersect(names(taxon_cols), unique(as.character(sub_df$scientificName)))
  if (length(taxa_in_data) == 0) taxa_in_data <- unique(as.character(sub_df$scientificName))
  pdf(output_pdf, width = 12, height = 5.5, pointsize = 11)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), oma = c(0, 0, 2, 0), family = "sans")

  # Panel A: Scatterplot of PC1 vs log(CentroidSize) with linear trendlines
  plot(sub_df$log_centsize, sub_df$PC1_raw, type = "n", xlab = expression(log(Centroid~Size)~"[pixels]"), ylab = "Raw Symmetric Morphospace PC1",
       main = sprintf("A. Foliar Allometry (R2 = %.3f, p = %.3f)", r2_size, p_size), font.main = 2)
  grid(col = "gray90", lty = "solid")
  for (sp in taxa_in_data) {
    sp_idx <- which(sub_df$scientificName == sp)
    if (length(sp_idx) > 0) {
      pt_col <- if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50"
      points(sub_df$log_centsize[sp_idx], sub_df$PC1_raw[sp_idx], pch = 21, bg = grDevices::adjustcolor(pt_col, alpha.f = 0.65), col = pt_col, cex = 1.1)
      if (length(sp_idx) >= 3) {
        sp_lm <- tryCatch(stats::lm(PC1_raw ~ log_centsize, data = sub_df[sp_idx, ]), error = function(e) NULL)
        if (!is.null(sp_lm)) abline(sp_lm, col = pt_col, lwd = 1.8)
      }
    }
  }
  glob_lm <- tryCatch(stats::lm(PC1_raw ~ log_centsize, data = sub_df), error = function(e) NULL)
  if (!is.null(glob_lm)) abline(glob_lm, col = "black", lwd = 2.2, lty = 2)
  legend("topleft", legend = c(taxa_in_data, "Global Trendline"),
         col = c(sapply(taxa_in_data, function(sp) if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50"), "black"),
         pch = c(rep(21, length(taxa_in_data)), NA), lty = c(rep(1, length(taxa_in_data)), 2),
         pt.bg = c(sapply(taxa_in_data, function(sp) if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50"), NA), bty = "n", cex = 0.75)

  # Panel B: Morphospace Biplot Before vs. After Allometric Correction
  x_lim <- range(c(sub_df$PC1_raw, sub_df$PC1_corrected), na.rm = TRUE); y_lim <- range(c(sub_df$PC2_raw, sub_df$PC2_corrected), na.rm = TRUE)
  plot(sub_df$PC1_raw, sub_df$PC2_raw, type = "n", xlim = x_lim, ylim = y_lim, xlab = "Morphospace PC1", ylab = "Morphospace PC2",
       main = "B. Morphospace Before vs. After Allometric Correction", font.main = 2)
  grid(col = "gray90", lty = "solid")
  arrows(sub_df$PC1_raw, sub_df$PC2_raw, sub_df$PC1_corrected, sub_df$PC2_corrected, length = 0.05, col = grDevices::adjustcolor("gray65", alpha.f = 0.6), lwd = 0.8)
  points(sub_df$PC1_raw, sub_df$PC2_raw, pch = 1, col = grDevices::adjustcolor("gray40", alpha.f = 0.6), cex = 0.9)
  for (sp in taxa_in_data) {
    sp_idx <- which(sub_df$scientificName == sp)
    if (length(sp_idx) > 0) {
      pt_col <- if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50"
      points(sub_df$PC1_corrected[sp_idx], sub_df$PC2_corrected[sp_idx], pch = 19, col = grDevices::adjustcolor(pt_col, alpha.f = 0.85), cex = 1.1)
    }
  }
  legend("bottomright", legend = c("Before (Raw Harmonic EFA)", "After (Allometric Residuals)", "Shift Vector"),
         pch = c(1, 19, NA), lty = c(NA, NA, 1), col = c("gray40", "#009E73", "gray65"), bty = "n", cex = 0.8)
  title("Packera dubia Complex: Foliar Size-Shape Allometry & Morphospace Correction", outer = TRUE, cex.main = 1.3)
  dev.off()
  message("Shape allometry regression plot saved -> ", output_pdf)
}

test_and_correct_allometry <- function(morpho_df, correct_mode = "auto", permutations = 999, seed = 42,
                                       output_report = "outputs/reports/allometry_significance_audit.csv",
                                       output_plot = "outputs/figures/shape_allometry_regression.pdf") {
  message("=== Foliar Allometry Assessment (Multivariate Size vs. Shape Covariation) ===")
  sym_names <- intersect(c(paste0("A", 1:12), paste0("D", 1:12)), names(morpho_df))
  if (!"log_centsize" %in% names(morpho_df) || length(sym_names) < 4) {
    message("Warning: log_centsize or symmetric harmonics not available. Skipping allometry test.")
    return(list(morpho_matrix = NULL, corrected = FALSE, anova_table = NULL, df = morpho_df))
  }
  valid_idx <- which(!is.na(morpho_df$log_centsize) & is.finite(morpho_df$log_centsize) & complete.cases(morpho_df[, sym_names]))
  if (length(valid_idx) < 10) return(list(morpho_matrix = NULL, corrected = FALSE, anova_table = NULL, df = morpho_df))

  sub_df <- morpho_df[valid_idx, , drop = FALSE]; sub_df$scientificName <- factor(sub_df$scientificName)
  Y <- as.matrix(sub_df[, sym_names, drop = FALSE]); set.seed(seed); allometry_fit <- NULL

  if (requireNamespace("vegan", quietly = TRUE)) {
    allometry_fit <- tryCatch(vegan::adonis2(Y ~ log_centsize * scientificName, data = sub_df,
                                             permutations = permutations, method = "euclidean", by = "terms"),
                              error = function(e) NULL)
  }
  if (!is.null(allometry_fit)) {
    anova_df <- as.data.frame(allometry_fit); anova_df$Term <- rownames(anova_df); rownames(anova_df) <- NULL
    if ("F.Model" %in% names(anova_df)) anova_df$F <- anova_df$F.Model
    if ("Pr(>F)" %in% names(anova_df)) anova_df$p_value <- anova_df[["Pr(>F)"]]
  } else {
    Y_cent <- scale(Y, center = TRUE, scale = FALSE); SS_tot <- sum(Y_cent^2); N_obs <- nrow(sub_df)
    fit1 <- stats::lm(Y_cent ~ log_centsize, data = sub_df); SS_size <- sum(stats::fitted(fit1)^2); df_size <- fit1$rank - 1
    fit2 <- stats::lm(Y_cent ~ log_centsize + scientificName, data = sub_df); SS_12 <- sum(stats::fitted(fit2)^2); SS_taxa <- SS_12 - SS_size; df_taxa <- fit2$rank - fit1$rank
    fit3 <- stats::lm(Y_cent ~ log_centsize * scientificName, data = sub_df); SS_123 <- sum(stats::fitted(fit3)^2); SS_int <- SS_123 - SS_12; df_int <- fit3$rank - fit2$rank
    SS_res <- SS_tot - SS_123; df_res <- N_obs - fit3$rank; MS_res <- SS_res / max(df_res, 1)
    F_size <- (SS_size / max(df_size, 1)) / MS_res; F_taxa <- (SS_taxa / max(df_taxa, 1)) / MS_res; F_int <- (SS_int / max(df_int, 1)) / MS_res
    anova_df <- data.frame(
      Term = c("log_centsize", "scientificName", "log_centsize:scientificName", "Residual", "Total"),
      Df = c(df_size, df_taxa, df_int, df_res, N_obs - 1),
      SumOfSqs = round(c(SS_size, SS_taxa, SS_int, SS_res, SS_tot), 4),
      R2 = round(c(SS_size / SS_tot, SS_taxa / SS_tot, SS_int / SS_tot, SS_res / SS_tot, 1.0), 6),
      F = c(round(F_size, 4), round(F_taxa, 4), round(F_int, 4), NA_real_, NA_real_),
      p_value = c(round(stats::pf(F_size, df_size, df_res, lower.tail = FALSE), 5),
                  round(stats::pf(F_taxa, df_taxa, df_res, lower.tail = FALSE), 5),
                  round(stats::pf(F_int, df_int, df_res, lower.tail = FALSE), 5), NA_real_, NA_real_),
      stringsAsFactors = FALSE
    )
  }
  report_cols <- intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(anova_df))
  anova_df <- anova_df[, report_cols, drop = FALSE]
  if (!is.null(output_report) && nzchar(output_report)) {
    dir.create(dirname(output_report), recursive = TRUE, showWarnings = FALSE)
    write.csv(anova_df, file = output_report, row.names = FALSE)
    message("Allometry significance audit report exported -> ", output_report)
  }

  size_row <- anova_df[grepl("^log_centsize$", anova_df$Term, ignore.case = TRUE), ]
  r2_size <- if (nrow(size_row) > 0) size_row$R2[1] else 0.0
  p_size <- if (nrow(size_row) > 0 && !is.na(size_row$p_value[1])) size_row$p_value[1] else 0.0
  int_row <- anova_df[grepl("log_centsize:scientificName|scientificName:log_centsize", anova_df$Term, ignore.case = TRUE), ]
  r2_int <- if (nrow(int_row) > 0) int_row$R2[1] else 0.0

  correct_flag <- tolower(trimws(as.character(correct_mode)))
  should_correct <- (correct_flag %in% c("true", "1", "yes")) || (correct_flag == "auto" && !is.na(r2_size) && r2_size >= 0.10)

  sub_df$PC1_raw <- sub_df$PC1; sub_df$PC2_raw <- sub_df$PC2
  fit_lm <- stats::lm(Y ~ log_centsize, data = sub_df)
  shape_residuals <- stats::residuals(fit_lm)
  pca_resids <- stats::prcomp(shape_residuals, center = TRUE, scale. = TRUE)
  sub_df$PC1_corrected <- round(pca_resids$x[, 1], 6); sub_df$PC2_corrected <- round(pca_resids$x[, 2], 6)

  if (should_correct) {
    message(sprintf("[INFO] Foliar allometry detected (R2 = %.2f, p = %.2f). Analyzing allometry-free shape residuals.", r2_size, p_size))
    morpho_matrix_for_gmm <- shape_residuals
    for (p in 1:min(5, ncol(pca_resids$x))) {
      morpho_df[[paste0("PC", p)]] <- NA_real_; morpho_df[valid_idx, paste0("PC", p)] <- round(pca_resids$x[, p], 6)
    }
  } else {
    message(sprintf("[INFO] Foliar allometry negligible (R2 = %.2f). Analyzing raw symmetric harmonics.", r2_size))
    morpho_matrix_for_gmm <- Y
  }

  if (!is.null(output_plot) && nzchar(output_plot)) generate_allometry_diagnostic_plot(sub_df, r2_size, p_size, r2_int, output_plot)
  return(list(morpho_matrix = morpho_matrix_for_gmm, corrected = should_correct,
              r2_size = r2_size, p_size = p_size, r2_interaction = r2_int, anova_table = anova_df, df = morpho_df))
}

# ------------------------------------------------------------------------------
# 3. Gaussian Mixture Modeling & Bayes Factor Testing (Seed = 42)
# ------------------------------------------------------------------------------
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150) {
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL
  bic_table <- data.frame(K = 1:max_k, LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)
  for (k in 1:max_k) {
    set.seed(42 + k)
    if (k == 1) {
      mu <- matrix(colMeans(X), nrow = 1); sigma <- list(cov(X) + diag(1e-6, p)); pi_k <- c(1.0)
    } else {
      km <- stats::kmeans(X, centers = k, nstart = 10, iter.max = 30); mu <- km$centers
      pi_k <- as.vector(table(factor(km$cluster, levels = 1:k))) / N
      sigma <- lapply(1:k, function(c) { sub_x <- X[km$cluster == c, , drop = FALSE]; if (nrow(sub_x) > p) cov(sub_x) + diag(1e-5, p) else diag(1e-2, p) })
    }
    loglik_old <- -Inf; resp <- matrix(0, nrow = N, ncol = k)
    for (iter in 1:max_iter) {
      log_dens <- matrix(0, nrow = N, ncol = k)
      for (c in 1:k) {
        diff_x <- t(t(X) - mu[c, ]); sig_inv <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        quad <- rowSums((diff_x %*% sig_inv) * diff_x)
        log_dens[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + quad)
      }
      max_l <- apply(log_dens, 1, max); log_sum_exp <- max_l + log(rowSums(exp(log_dens - max_l)))
      resp <- exp(log_dens - log_sum_exp); loglik <- sum(log_sum_exp)
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik; N_k <- colSums(resp); pi_k <- N_k / N
      for (c in 1:k) {
        if (N_k[c] > 1e-6) {
          mu[c, ] <- colSums(resp[, c] * X) / N_k[c]; diff_x <- t(t(X) - mu[c, ])
          sigma[[c]] <- (t(diff_x) %*% (resp[, c] * diff_x)) / N_k[c] + diag(1e-5, p)
        }
      }
    }
    n_params <- (k - 1) + k * p + k * (p * (p + 1) / 2); bic_val <- 2 * loglik - n_params * log(N)
    bic_table$LogLik[k] <- round(loglik, 2); bic_table$n_params[k] <- n_params; bic_table$BIC[k] <- round(bic_val, 2)
    if (bic_val > best_bic) {
      best_bic <- bic_val; best_model <- list(k = k, classification = apply(resp, 1, which.max), uncertainty = 1 - apply(resp, 1, max), bic = bic_val)
    }
  }
  return(list(best = best_model, bic_table = bic_table))
}

run_gmm_cluster_analysis <- function(df, feature_cols, max_k = 8) {
  valid_rows <- which(complete.cases(df[, feature_cols])); N_valid <- length(valid_rows)
  if (N_valid < 3) {
    message("Warning: Too few complete observations for GMM clustering.")
    return(list(classification = rep(NA_integer_, nrow(df)), uncertainty = rep(NA_real_, nrow(df)),
                best_k = 1, model_name = "None",
                bic_table = data.frame(K = 1:max_k, BIC = NA_real_, Two_Delta_BIC = NA_real_, Evidence = "Inconclusive")))
  }
  X <- as.matrix(df[valid_rows, feature_cols, drop = FALSE])
  message(sprintf("Fitting Label-Blind Gaussian Mixture Models (K = 1 to %d, N = %d, seed = 42)...", max_k, N_valid))
  set.seed(42); best_k <- 1; best_name <- "Full_Covariance_EM"; bic_table <- NULL; cls_sub <- NULL; unc_sub <- NULL

  if (requireNamespace("mclust", quietly = TRUE)) {
    mc <- tryCatch(mclust::Mclust(X, G = 1:max_k, verbose = FALSE), error = function(e) NULL)
    if (!is.null(mc)) {
      best_k <- mc$G; best_name <- mc$modelName; cls_sub <- mc$classification; unc_sub <- mc$uncertainty
      bic_vals <- apply(mc$BIC, 1, function(row) if (all(is.na(row))) NA_real_ else max(row, na.rm = TRUE))
      bic_table <- data.frame(K = 1:max_k, BIC = round(bic_vals[1:max_k], 2))
    }
  }
  if (is.null(bic_table)) {
    pure_res <- fit_gmm_em_pure(X, max_k = max_k); best_k <- pure_res$best$k
    cls_sub <- pure_res$best$classification; unc_sub <- pure_res$best$uncertainty; bic_table <- pure_res$bic_table
  }

  cls <- rep(NA_integer_, nrow(df)); unc <- rep(NA_real_, nrow(df))
  cls[valid_rows] <- cls_sub; unc[valid_rows] <- unc_sub
  bic_table$Delta_BIC_prev <- c(0, diff(bic_table$BIC)); bic_table$Two_Delta_BIC <- round(bic_table$Delta_BIC_prev, 2)
  bic_table$BayesFactor_vs_Null <- round(bic_table$BIC - bic_table$BIC[1], 2)
  bic_table$Evidence <- sapply(bic_table$Two_Delta_BIC, function(val) {
    if (is.na(val) || val < 0) return("Negative (favors K-1)")
    if (val < 2) return("Weak / Inconclusive")
    if (val < 6) return("Positive Evidence")
    if (val < 10) return("Strong Evidence")
    return("Decisive Evidence (Species Boundary)")
  })
  return(list(classification = cls, uncertainty = unc, best_k = best_k, model_name = best_name, bic_table = bic_table))
}

# ------------------------------------------------------------------------------
# 4. Canonical Discriminant Analysis with Passive Sample Projection
# ------------------------------------------------------------------------------
run_cda_with_passive_projection <- function(df, feature_cols, target_taxa) {
  message("Configuring Canonical Discriminant Analysis in MorphoTools2 Architecture...")
  is_active <- (df$determiner_tier == "Tier_1_Gold") & (df$species_standardized %in% target_taxa) & complete.cases(df[, feature_cols])
  active_idx <- which(is_active); passive_idx <- which(!is_active)
  message(sprintf("Active Anchors (Tier 1 Gold, Complete): %d | Passive Projected: %d", length(active_idx), length(passive_idx)))
  if (length(active_idx) < length(target_taxa)) stop("Insufficient complete active Tier 1 Gold anchors for CDA.")

  X <- as.matrix(df[, feature_cols, drop = FALSE]); X_a <- X[active_idx, , drop = FALSE]
  y_a <- df$species_standardized[active_idx]; N_a <- nrow(X_a); p <- ncol(X_a)
  g_levels <- target_taxa; g <- length(g_levels); grand_mean <- colMeans(X_a)
  B <- matrix(0, p, p); W <- matrix(0, p, p)
  group_means <- matrix(0, g, p, dimnames = list(g_levels, feature_cols)); group_counts <- integer(g)

  for (k in seq_along(g_levels)) {
    grp <- g_levels[k]; sub_x <- X_a[y_a == grp, , drop = FALSE]; group_counts[k] <- nrow(sub_x)
    if (nrow(sub_x) > 0) {
      m_k <- colMeans(sub_x); group_means[k, ] <- m_k; diff_m <- m_k - grand_mean
      B <- B + group_counts[k] * (diff_m %*% t(diff_m)); diff_x <- t(t(sub_x) - m_k); W <- W + (t(diff_x) %*% diff_x)
    }
  }

  S_reg <- (W / max(N_a - g, 1)) + diag(1e-7, p); eig <- eigen(solve(S_reg, B))
  real_idx <- which(abs(Im(eig$values)) < 1e-6); real_vals <- Re(eig$values[real_idx]); real_vecs <- Re(eig$vectors[, real_idx, drop = FALSE])
  order_idx <- order(real_vals, decreasing = TRUE); num_axes <- min(g - 1, p)
  eig_vals <- real_vals[order_idx][1:num_axes]; eig_vecs <- real_vecs[, order_idx[1:num_axes], drop = FALSE]

  for (j in seq_len(num_axes)) {
    v <- eig_vecs[, j]; s <- as.numeric(sqrt(t(v) %*% S_reg %*% v))
    if (s > 1e-8) eig_vecs[, j] <- v / s
  }
  var_pct <- round((eig_vals / sum(eig_vals)) * 100, 2)
  message("=== CDA Canonical Variates Eigenvalues & Variance ===")
  for (j in seq_len(num_axes)) message(sprintf("  Can%d: Eigenvalue = %7.4f (%5.2f%% variation)", j, eig_vals[j], var_pct[j]))

  Z_active <- t(t(X_a) - grand_mean) %*% eig_vecs; struct_corr <- stats::cor(X_a, Z_active)
  X_proj <- X
  for (col_idx in seq_len(p)) {
    na_mask <- is.na(X_proj[, col_idx])
    if (any(na_mask)) X_proj[na_mask, col_idx] <- grand_mean[col_idx]
  }

  Z_all <- t(t(X_proj) - grand_mean) %*% eig_vecs; canonical_centroids <- t(t(group_means) - grand_mean) %*% eig_vecs
  priors <- group_counts / sum(group_counts); post_probs <- matrix(0, nrow(df), g, dimnames = list(NULL, g_levels))
  for (k in seq_along(g_levels)) {
    diff_z <- t(t(Z_all) - canonical_centroids[k, ])
    post_probs[, k] <- priors[k] * exp(-0.5 * rowSums(diff_z^2))
  }
  post_probs <- post_probs / rowSums(post_probs); pred_idx <- apply(post_probs, 1, which.max)

  return(list(eigenvalues = eig_vals, variance_pct = var_pct, canonical_scores = Z_all,
              centroids = canonical_centroids, predicted_taxon = g_levels[pred_idx],
              posterior_prob = apply(post_probs, 1, max), is_active = is_active, structure_loadings = struct_corr))
}

# ------------------------------------------------------------------------------
# 5. Outlier Flagging & Herbarium Misidentification Auditing
# ------------------------------------------------------------------------------
audit_misidentifications <- function(df, cda_res, gmm_res) {
  message("Auditing specimens for morphospace outliers & label discordances...")
  n <- nrow(df); flags <- logical(n); triage <- character(n); reasons <- character(n)

  for (i in seq_len(n)) {
    raw_sp <- df$species_standardized[i]; pred_sp <- cda_res$predicted_taxon[i]
    post_p <- cda_res$posterior_prob[i]; is_act <- cda_res$is_active[i]
    if (is_act) {
      if (raw_sp == pred_sp) { flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Verified_Tier1_Gold_Anchor" }
      else if (post_p >= 0.85) { flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- sprintf("Tier_1_Gold_Morphological_Discordance_to_%s", gsub(" ", "_", pred_sp)) }
      else { flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Tier_1_Gold_Borderline_Variant" }
    } else {
      if (raw_sp == pred_sp && post_p >= 0.75) { flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Congruent_Bronze_Determination" }
      else if (raw_sp != pred_sp && post_p >= 0.70) {
        flags[i] <- TRUE; triage[i] <- "HIGH"
        if (grepl("dubia|tomentos", raw_sp, ignore.case = TRUE) && pred_sp == "Packera anonyma") {
          reasons[i] <- "Glabrescent_P_dubia_Misidentified_as_P_anonyma"
        } else { reasons[i] <- sprintf("Tier_3_Bronze_Misidentified_as_%s", gsub(" ", "_", pred_sp)) }
      } else if (post_p < 0.60) { flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Morphological_Intermediate_Ambiguous" }
      else { flags[i] <- (raw_sp != pred_sp); triage[i] <- if (raw_sp != pred_sp) "MEDIUM" else "LOW"
             reasons[i] <- if (raw_sp != pred_sp) "Candidate_Discordance_Moderate_Posterior" else "Minor_Uncertainty" }
    }
  }

  df$gmm_cluster <- gmm_res$classification; df$gmm_uncertainty <- round(gmm_res$uncertainty, 4)
  df$cda_predicted_taxon <- cda_res$predicted_taxon; df$cda_posterior_prob <- round(cda_res$posterior_prob, 4)
  df$can1 <- round(cda_res$canonical_scores[, 1], 5)
  df$can2 <- if (ncol(cda_res$canonical_scores) >= 2) round(cda_res$canonical_scores[, 2], 5) else 0.0
  df$is_passive <- !cda_res$is_active; df$morphometric_outlier_flag <- flags
  df$triage_priority <- triage; df$discordance_reason <- reasons
  if ("AR_inv" %in% names(df)) df$reproductive_trait_status <- ifelse(is.na(df$AR_inv), "Unexpanded_or_Absent", "Measured")
  message(sprintf("Audit Complete: %d morphospace outliers/discordances flagged (%d HIGH, %d MEDIUM)",
                  sum(flags), sum(triage == "HIGH"), sum(triage == "MEDIUM")))
  return(df)
}

# ------------------------------------------------------------------------------
# 5b. Institutional Batch-Effect Audit (ANOVA & Variance Explained)
# ------------------------------------------------------------------------------
audit_institutional_batch_effects <- function(morpho_df, output_csv = "outputs/reports/institutional_batch_effect_audit.csv", min_obs = 15) {
  message("=== Institutional Batch-Effect Audit (Digitization Rig Diagnostic) ===")
  if (!("institutionCode" %in% names(morpho_df))) return(NULL)
  valid_df <- morpho_df[!is.na(morpho_df$institutionCode) & trimws(morpho_df$institutionCode) != "", ]
  valid_df <- valid_df[!is.na(valid_df$PC1) & !is.na(valid_df$PC2), ]
  if (nrow(valid_df) < 10) return(NULL)

  inst_counts <- table(valid_df$institutionCode); rep_insts <- names(inst_counts[inst_counts >= min_obs])
  audit_subset <- if (length(rep_insts) >= 2) valid_df[valid_df$institutionCode %in% rep_insts, ] else valid_df
  audit_subset$institutionCode <- factor(audit_subset$institutionCode)
  lm_pc1 <- lm(PC1 ~ institutionCode, data = audit_subset); lm_pc2 <- lm(PC2 ~ institutionCode, data = audit_subset)
  s1 <- summary(lm_pc1); s2 <- summary(lm_pc2); aov1 <- anova(lm_pc1); aov2 <- anova(lm_pc2)
  r2_pc1 <- s1$r.squared; r2_pc2 <- s2$r.squared

  if (r2_pc1 >= 0.05 || r2_pc2 >= 0.05) {
    warning(sprintf("Institutional batch effect R^2 exceeds 0.05 threshold (PC1 R^2=%.4f, PC2 R^2=%.4f)", r2_pc1, r2_pc2))
  } else {
    message("  [PASS] Assertion satisfied: Institutional identity explains < 5% variance (R^2 < 0.05).")
  }

  anova_df <- data.frame(
    Response = c("PC1", "PC1", "PC2", "PC2"), Term = c("institutionCode", "Residuals", "institutionCode", "Residuals"),
    Df = c(aov1$Df[1], aov1$Df[2], aov2$Df[1], aov2$Df[2]),
    Sum_Sq = round(c(aov1$`Sum Sq`[1], aov1$`Sum Sq`[2], aov2$`Sum Sq`[1], aov2$`Sum Sq`[2]), 4),
    Mean_Sq = round(c(aov1$`Mean Sq`[1], aov1$`Mean Sq`[2], aov2$`Mean Sq`[1], aov2$`Mean Sq`[2]), 4),
    F_value = c(round(aov1$`F value`[1], 4), NA, round(aov2$`F value`[1], 4), NA),
    P_value = c(signif(aov1$`Pr(>F)`[1], 4), NA, signif(aov2$`Pr(>F)`[1], 4), NA),
    R_squared = round(c(r2_pc1, NA, r2_pc2, NA), 4), Adj_R_squared = round(c(s1$adj.r.squared, NA, s2$adj.r.squared, NA), 4),
    Batch_Effect_Status = c(ifelse(r2_pc1 < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA, ifelse(r2_pc2 < 0.05, "Pass (<0.05)", "Exceeds 0.05"), NA),
    stringsAsFactors = FALSE
  )
  if (!is.null(output_csv) && nzchar(output_csv)) {
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
    write.csv(anova_df, file = output_csv, row.names = FALSE)
  }
  return(invisible(list(lm_pc1 = lm_pc1, lm_pc2 = lm_pc2, anova_table = anova_df)))
}

# ------------------------------------------------------------------------------
# 6. Publication-Quality CDA Biplot Vector Graphic (PDF)
# ------------------------------------------------------------------------------
generate_cda_biplot_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL) {
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(output_pdf, width = 12, height = 10, pointsize = 11)
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.2, 1.5), oma = c(1, 1, 2, 1), family = "sans")

  taxon_cols <- c("Packera anonyma" = "#E69F00", "Packera dubia" = "#009E73", "Packera plattensis" = "#56B4E9", "Packera paupercula" = "#CC79A7")
  cluster_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  is_joint <- !is.null(feature_cols) && ("AR_inv_std" %in% feature_cols)
  cda_title <- if (is_joint) "A. Joint Vegetative + Reproductive CDA (MorphoTools2)" else "A. CDA with Passive Sample Projection (MorphoTools2)"

  plot(audit_df$can1, audit_df$can2, type = "n", xlab = sprintf("Canonical Axis 1 (%.1f%% Between-Group Var)", cda_res$variance_pct[1]),
       ylab = sprintf("Canonical Axis 2 (%.1f%% Between-Group Var)", cda_res$variance_pct[2]), main = cda_title, font.main = 2)
  grid(col = "gray90", lty = "solid")
  active_sub <- audit_df[!audit_df$is_passive, ]
  for (sp in names(taxon_cols)) {
    sp_idx <- which(active_sub$species_standardized == sp)
    if (length(sp_idx) > 0) points(active_sub$can1[sp_idx], active_sub$can2[sp_idx], pch = 21, bg = grDevices::adjustcolor(taxon_cols[sp], alpha.f = 0.65), col = taxon_cols[sp], cex = 0.9)
  }
  passive_sub <- audit_df[audit_df$is_passive, ]
  if (nrow(passive_sub) > 0) {
    clean_p <- passive_sub[!passive_sub$morphometric_outlier_flag, ]; flagged_p <- passive_sub[passive_sub$morphometric_outlier_flag, ]
    if (nrow(clean_p) > 0) points(clean_p$can1, clean_p$can2, pch = 5, col = grDevices::adjustcolor("gray40", alpha.f = 0.5), cex = 0.8)
    if (nrow(flagged_p) > 0) points(flagged_p$can1, flagged_p$can2, pch = 23, bg = "#D55E00", col = "black", cex = 1.1, lwd = 1.2)
  }
  for (k in seq_len(nrow(cda_res$centroids))) {
    points(cda_res$centroids[k, 1], cda_res$centroids[k, 2], pch = 3, col = "black", cex = 2.2, lwd = 3)
    text(cda_res$centroids[k, 1], cda_res$centroids[k, 2], labels = rownames(cda_res$centroids)[k], pos = 3, cex = 0.8, font = 4)
  }
  legend("topleft", legend = c("Tier 1 Gold Anchors", "Tier 3 Passive Congruent", "Flagged Morphospace Outliers"),
         pch = c(21, 5, 23), pt.bg = c("#009E73", NA, "#D55E00"), col = c("#009E73", "gray40", "black"), bty = "n", cex = 0.8)

  col_vec <- ifelse(!is.na(audit_df$gmm_cluster), grDevices::adjustcolor(cluster_cols[(audit_df$gmm_cluster - 1) %% length(cluster_cols) + 1], alpha.f = 0.6), grDevices::adjustcolor("gray70", alpha.f = 0.3))
  plot(audit_df$PC1, audit_df$PC2, pch = 20, col = col_vec, xlab = "Morphospace PC1", ylab = "Morphospace PC2",
       main = sprintf("B. Label-Blind GMM Clusters (Optimal K = %d)", gmm_res$best_k), font.main = 2)
  grid(col = "gray90", lty = "solid")
  legend("bottomright", legend = paste("Cluster", 1:gmm_res$best_k), col = cluster_cols[1:gmm_res$best_k], pch = 20, bty = "n", cex = 0.8, ncol = 2)

  bic_df <- gmm_res$bic_table; delta_bic <- bic_df$Two_Delta_BIC; delta_bic[1] <- 0
  bar_cols <- ifelse(delta_bic >= 10, "#0072B2", ifelse(delta_bic >= 6, "#56B4E9", ifelse(delta_bic >= 2, "#F0E442", "#999999")))
  bp <- barplot(delta_bic, names.arg = paste0("K=", bic_df$K), col = bar_cols, border = "white",
                xlab = "Number of Mixture Components (K)", ylab = "Bayes Factor (2ΔBIC vs K-1)",
                main = "C. Species Boundary Evidence (Kass & Raftery 1995)", font.main = 2)
  abline(h = c(2, 6, 10), lty = 2, col = c("gray60", "gray40", "red3"))

  hist(audit_df$cda_posterior_prob[audit_df$morphometric_outlier_flag], breaks = 15, col = grDevices::adjustcolor("#D55E00", alpha.f = 0.6), border = "white",
       xlab = "CDA Posterior Classification Confidence", ylab = "Specimen Count", main = "D. Posterior Confidence of Flagged Outliers", font.main = 2)
  hist(audit_df$cda_posterior_prob[!audit_df$morphometric_outlier_flag], breaks = 20, col = grDevices::adjustcolor("#009E73", alpha.f = 0.35), border = "white", add = TRUE)
  legend("topleft", legend = c("Flagged Outliers", "Congruent Vouchers"), fill = c(grDevices::adjustcolor("#D55E00", alpha.f = 0.6), grDevices::adjustcolor("#009E73", alpha.f = 0.35)), bty = "n", cex = 0.8)

  title(if (is_joint) "Packera dubia Complex: Joint Morphometrics & Passive CDA" else "Packera dubia Complex: Morphometrics & Passive CDA", outer = TRUE, cex.main = 1.4)
  dev.off()
  message("CDA biplot graphic successfully saved to ", output_pdf)
}

# ------------------------------------------------------------------------------
# 7. Main Workflow Orchestrator
# ------------------------------------------------------------------------------
run_gmm_morphotools_pipeline <- function(opts) {
  message("==================================================================")
  message("Starting Packera Morphometrics GMM & Passive MorphoTools2 CDA")
  message("Vouchers: ", opts$vouchers, " | Harmonics: ", opts$harmonics)
  message("Allometry Mode: ", opts$correct_allometry, " | Report: ", opts$allometry_report)
  message("==================================================================")
  if (!file.exists(opts$harmonics)) stop("Harmonics file not found: ", opts$harmonics)
  efa_df <- read.csv(opts$harmonics, stringsAsFactors = FALSE)
  vouchers_df <- if (file.exists(opts$vouchers)) read.csv(opts$vouchers, stringsAsFactors = FALSE) else NULL
  if (!is.null(vouchers_df)) vouchers_df <- calculate_circular_phenology(vouchers_df)

  pca_cols <- paste0("PC", 1:opts$num_pcs); missing_pca <- setdiff(pca_cols, names(efa_df))
  if (length(missing_pca) > 0) stop("Missing PCA columns in EFA table: ", paste(missing_pca, collapse = ", "))

  closed_df <- efa_df[efa_df$assigned_tier %in% c("Tier_1_Direct", "Tier_2_Reflected"), ]
  valid_idx <- which(complete.cases(closed_df[, pca_cols])); closed_df <- closed_df[valid_idx, ]

  # Ensure primary GMM and CDA run on specimen-level median table (1 row = 1 voucher collection event)
  if (any(duplicated(closed_df$catalogNumber))) {
    message("Aggregating to specimen-level median profile (1 row = 1 voucher collection event)...")
    num_cols <- intersect(c(pca_cols, "foliar_variance", "leaf_count", "aspect_ratio", "area_px", "centroid_size", "log_centsize"), names(closed_df))
    meta_cols <- setdiff(names(closed_df), c(num_cols, "shape_id", "leaf_id", "mask_source"))
    closed_df <- closed_df %>%
      dplyr::group_by(catalogNumber) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(num_cols), ~ stats::median(.x, na.rm = TRUE)),
                       dplyr::across(dplyr::all_of(meta_cols[meta_cols != "catalogNumber"]), ~ stats::na.omit(.x)[1]),
                       .groups = "drop")
  }
  message(sprintf("Specimen-level collection units for morphometric modeling (1 row = 1 voucher): %d", nrow(closed_df)))
  closed_df$species_standardized <- standardize_packera_taxon(closed_df$species_raw)

  # 0. Foliar Allometry Multivariate Testing & Shape Residual Correction
  allom_res <- test_and_correct_allometry(closed_df, correct_mode = opts$correct_allometry, permutations = opts$permutations,
                                          seed = opts$seed, output_report = opts$allometry_report, output_plot = opts$allometry_plot)
  if (!is.null(allom_res$df)) closed_df <- allom_res$df

  include_repro <- FALSE; feature_cols <- pca_cols
  if (!isTRUE(opts$outline_only)) {
    if ("capitulum_aspect_ratio" %in% names(closed_df)) {
      closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$capitulum_aspect_ratio))
    } else if ("AR_inv" %in% names(closed_df)) {
      closed_df$AR_inv <- suppressWarnings(as.numeric(closed_df$AR_inv))
    } else if (!is.null(opts$reproductive) && file.exists(opts$reproductive)) {
      repro_df <- read.csv(opts$reproductive, stringsAsFactors = FALSE)
      if ("capitulum_aspect_ratio" %in% names(repro_df)) repro_df$AR_inv <- suppressWarnings(as.numeric(repro_df$capitulum_aspect_ratio))
      r_cols <- intersect(names(repro_df), c("catalogNumber", "AR_inv", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm"))
      r_sub <- repro_df[!duplicated(repro_df$catalogNumber), r_cols, drop = FALSE]
      closed_df <- merge(closed_df, r_sub, by = "catalogNumber", all.x = TRUE)
    }

    if ("AR_inv" %in% names(closed_df)) {
      valid_repro <- which(!is.na(closed_df$AR_inv) & is.finite(closed_df$AR_inv))
      if (length(valid_repro) >= 5) {
        closed_df$AR_inv_std <- NA_real_; m_ar <- mean(closed_df$AR_inv[valid_repro]); s_ar <- stats::sd(closed_df$AR_inv[valid_repro])
        closed_df$AR_inv_std[valid_repro] <- if (is.finite(s_ar) && s_ar > 1e-6) round((closed_df$AR_inv[valid_repro] - m_ar) / s_ar, 5) else 0.0
        include_repro <- TRUE; feature_cols <- c(pca_cols, "AR_inv_std")
        message(sprintf("Integrated macro-reproductive trait (AR_inv): %d / %d vouchers with valid capitula (mean=%.3f, sd=%.3f)",
                        length(valid_repro), nrow(closed_df), m_ar, s_ar))
      }
    }
  }

  target_taxa <- c("Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula")
  gmm_res <- run_gmm_cluster_analysis(closed_df, feature_cols, max_k = opts$max_k)
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE)
  write.csv(gmm_res$bic_table, file = opts$output_report, row.names = FALSE)
  message("Bayes Factor (2ΔBIC) summary report saved: ", opts$output_report)

  cda_res <- run_cda_with_passive_projection(closed_df, feature_cols, target_taxa)
  audited_df <- audit_misidentifications(closed_df, cda_res, gmm_res)

  if (!is.null(vouchers_df)) {
    v_meta <- intersect(names(vouchers_df), c("catalogNumber", "institutionCode", "county", "stateProvince", "latitude", "longitude", "doy", "pheno_sin", "pheno_cos", "regional_group"))
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), v_meta, drop = FALSE]
    audited_df <- merge(audited_df, v_sub, by = "catalogNumber", all.x = TRUE, suffixes = c("", "_meta"))
  }
  audit_institutional_batch_effects(audited_df, output_csv = opts$output_batch_audit)

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "species_raw", "species_standardized",
                 "determiner_tier", "cda_predicted_taxon", "cda_posterior_prob", "is_passive",
                 "morphometric_outlier_flag", "triage_priority", "discordance_reason", "gmm_cluster",
                 "gmm_uncertainty", "can1", "can2", "pheno_sin", "pheno_cos", "PC1", "PC2", "PC3", "PC4", "PC5",
                 "centroid_size", "log_centsize",
                 "AR_inv", "AR_inv_std", "reproductive_trait_status", "capitula_count",
                 "involucre_height_px", "involucre_width_px", "involucre_height_mm", "involucre_width_mm")
  lead_cols <- intersect(lead_cols, names(audited_df))
  out_df <- audited_df[, c(lead_cols, setdiff(names(audited_df), lead_cols))]

  dir.create(dirname(opts$output_flags), recursive = TRUE, showWarnings = FALSE)
  write.csv(out_df, file = opts$output_flags, row.names = FALSE, na = "")
  message("Master morphometric flags table exported: ", opts$output_flags, " (Rows: ", nrow(out_df), ")")
  generate_cda_biplot_pdf(audited_df, cda_res, gmm_res, opts$output_plot, feature_cols)
  message("==================================================================")
  message("Morphometric Pipeline Complete Successfully.")
  message("==================================================================")
  return(invisible(out_df))
}

if (sys.nframe() == 0) {
  opts <- parse_args_robust()
  run_gmm_morphotools_pipeline(opts)
}
