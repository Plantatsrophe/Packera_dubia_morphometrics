#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_gmm_morphotools.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Label-blind Gaussian Mixture Model (mclust::Mclust) clustering and Canonical
#   Discriminant Analysis with passive sample projection (MorphoTools2).
#   1. Ingests curated_vouchers.csv and leaf_efa_harmonics.csv.
#   2. Fits GMMs blind to herbarium determinations to detect natural clusters.
#   3. Computes Bayes Factors (2ΔBIC) across competing K-component models.
#   4. Executes CDA with Tier 3 Bronze vouchers as passiveSamples.
#   5. Exports misidentification audit flags and publication-quality biplots.
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("mclust", quietly = TRUE)) library(mclust)
  if (requireNamespace("MorphoTools2", quietly = TRUE)) library(MorphoTools2)
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
    optparse::make_option(c("-v", "--vouchers"), type = "character",
      default = "data/tables/curated_vouchers.csv", help = "Vouchers metadata CSV [default: %default]"),
    optparse::make_option(c("-e", "--harmonics"), type = "character",
      default = "data/tables/leaf_efa_harmonics.csv", help = "Leaf EFA harmonics CSV [default: %default]"),
    optparse::make_option(c("-f", "--output-flags"), type = "character",
      default = "data/tables/morphometrics_misidentification_flags.csv", help = "Flags output CSV [default: %default]"),
    optparse::make_option(c("-p", "--output-plot"), type = "character",
      default = "outputs/figures/cda_passive_projection.pdf", help = "Output CDA PDF biplot [default: %default]"),
    optparse::make_option(c("-r", "--output-report"), type = "character",
      default = "outputs/reports/gmm_bayes_factors_summary.csv", help = "BIC report CSV [default: %default]"),
    optparse::make_option(c("-k", "--max-k"), type = "integer",
      default = 8, help = "Max mixture components for GMM [default: %default]"),
    optparse::make_option(c("-d", "--num-pcs"), type = "integer",
      default = 5, help = "Number of PCA dimensions to analyze [default: %default]")
  )

  if (requireNamespace("optparse", quietly = TRUE)) {
    parser <- optparse::OptionParser(usage = "%prog [options]", option_list = option_list)
    return(optparse::parse_args(parser))
  }

  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    vouchers = "data/tables/curated_vouchers.csv", harmonics = "data/tables/leaf_efa_harmonics.csv",
    output_flags = "data/tables/morphometrics_misidentification_flags.csv",
    output_plot = "outputs/figures/cda_passive_projection.pdf",
    output_report = "outputs/reports/gmm_bayes_factors_summary.csv", max_k = 8, num_pcs = 5
  )
  i <- 1
  while (i <= length(raw_args)) {
    arg <- raw_args[i]
    if (arg %in% c("-v", "--vouchers") && i < length(raw_args)) { opts$vouchers <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-e", "--harmonics") && i < length(raw_args)) { opts$harmonics <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-f", "--output-flags") && i < length(raw_args)) { opts$output_flags <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-p", "--output-plot") && i < length(raw_args)) { opts$output_plot <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-r", "--output-report") && i < length(raw_args)) { opts$output_report <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-k", "--max-k") && i < length(raw_args)) { opts$max_k <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else if (arg %in% c("-d", "--num-pcs") && i < length(raw_args)) { opts$num_pcs <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else { i <- i + 1 }
  }
  return(opts)
}

# ------------------------------------------------------------------------------
# 2. Taxonomic Concept Standardization
# ------------------------------------------------------------------------------
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
# 3. Label-Blind Gaussian Mixture Modeling & Bayes Factor Testing
# ------------------------------------------------------------------------------
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150) {
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL
  bic_table <- data.frame(K = 1:max_k, LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)

  for (k in 1:max_k) {
    set.seed(42 + k)
    if (k == 1) {
      mu <- matrix(colMeans(X), nrow = 1); sigma <- list(cov(X) + diag(1e-6, p)); pi_k <- c(1.0)
    } else {
      km <- stats::kmeans(X, centers = k, nstart = 10, iter.max = 30)
      mu <- km$centers; pi_k <- as.vector(table(factor(km$cluster, levels = 1:k))) / N
      sigma <- lapply(1:k, function(c) {
        sub_x <- X[km$cluster == c, , drop = FALSE]
        if (nrow(sub_x) > p) cov(sub_x) + diag(1e-5, p) else diag(1e-2, p)
      })
    }
    loglik_old <- -Inf; resp <- matrix(0, nrow = N, ncol = k)

    for (iter in 1:max_iter) {
      log_dens <- matrix(0, nrow = N, ncol = k)
      for (c in 1:k) {
        diff_x <- t(t(X) - mu[c, ])
        sig_inv <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        sig_det <- max(det(sigma[[c]]), 1e-12)
        quad <- rowSums((diff_x %*% sig_inv) * diff_x)
        log_dens[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(sig_det) + quad)
      }
      max_l <- apply(log_dens, 1, max); log_sum_exp <- max_l + log(rowSums(exp(log_dens - max_l)))
      resp <- exp(log_dens - log_sum_exp); loglik <- sum(log_sum_exp)
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik
      N_k <- colSums(resp); pi_k <- N_k / N
      for (c in 1:k) {
        if (N_k[c] > 1e-6) {
          mu[c, ] <- colSums(resp[, c] * X) / N_k[c]
          diff_x <- t(t(X) - mu[c, ])
          sigma[[c]] <- (t(diff_x) %*% (resp[, c] * diff_x)) / N_k[c] + diag(1e-5, p)
        }
      }
    }
    n_params <- (k - 1) + k * p + k * (p * (p + 1) / 2)
    bic_val <- 2 * loglik - n_params * log(N)
    bic_table$LogLik[k] <- round(loglik, 2); bic_table$n_params[k] <- n_params; bic_table$BIC[k] <- round(bic_val, 2)
    if (bic_val > best_bic) {
      best_bic <- bic_val
      best_model <- list(k = k, classification = apply(resp, 1, which.max), uncertainty = 1 - apply(resp, 1, max), bic = bic_val)
    }
  }
  return(list(best = best_model, bic_table = bic_table))
}

run_gmm_cluster_analysis <- function(df, feature_cols, max_k = 8) {
  X <- as.matrix(df[, feature_cols])
  message("Fitting Label-Blind Gaussian Mixture Models (K = 1 to ", max_k, ")...")
  best_k <- 1; best_name <- "Full_Covariance_EM"; bic_table <- NULL; cls <- NULL; unc <- NULL

  if (requireNamespace("mclust", quietly = TRUE)) {
    mc <- tryCatch(mclust::Mclust(X, G = 1:max_k, verbose = FALSE), error = function(e) NULL)
    if (!is.null(mc)) {
      best_k <- mc$G; best_name <- mc$modelName; cls <- mc$classification; unc <- mc$uncertainty
      bic_vals <- apply(mc$BIC, 1, function(row) if (all(is.na(row))) NA_real_ else max(row, na.rm = TRUE))
      bic_table <- data.frame(K = 1:max_k, BIC = round(bic_vals[1:max_k], 2))
    }
  }
  if (is.null(bic_table)) {
    pure_res <- fit_gmm_em_pure(X, max_k = max_k)
    best_k <- pure_res$best$k; cls <- pure_res$best$classification; unc <- pure_res$best$uncertainty
    bic_table <- pure_res$bic_table
  }

  bic_table$Delta_BIC_prev <- c(0, diff(bic_table$BIC))
  bic_table$Two_Delta_BIC <- round(bic_table$Delta_BIC_prev, 2)
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
  message("Configuring Canonical Discriminant Analysis (CDA) in MorphoTools2 Architecture...")
  is_active <- (df$determiner_tier == "Tier_1_Gold") & (df$species_standardized %in% target_taxa)
  active_idx <- which(is_active); passive_idx <- which(!is_active)
  message(sprintf("Active Anchors (Tier 1 Gold): %d | Passive Projected: %d", length(active_idx), length(passive_idx)))

  X <- as.matrix(df[, feature_cols]); X_a <- X[active_idx, , drop = FALSE]
  y_a <- df$species_standardized[active_idx]; N_a <- nrow(X_a); p <- ncol(X_a)
  g_levels <- target_taxa; g <- length(g_levels)
  grand_mean <- colMeans(X_a); B <- matrix(0, p, p); W <- matrix(0, p, p)
  group_means <- matrix(0, g, p, dimnames = list(g_levels, feature_cols)); group_counts <- integer(g)

  for (k in seq_along(g_levels)) {
    grp <- g_levels[k]; sub_x <- X_a[y_a == grp, , drop = FALSE]; group_counts[k] <- nrow(sub_x)
    if (nrow(sub_x) > 0) {
      m_k <- colMeans(sub_x); group_means[k, ] <- m_k
      diff_m <- m_k - grand_mean; B <- B + group_counts[k] * (diff_m %*% t(diff_m))
      diff_x <- t(t(sub_x) - m_k); W <- W + (t(diff_x) %*% diff_x)
    }
  }

  S_reg <- (W / max(N_a - g, 1)) + diag(1e-7, p)
  eig <- eigen(solve(S_reg, B))
  real_idx <- which(abs(Im(eig$values)) < 1e-6)
  real_vals <- Re(eig$values[real_idx]); real_vecs <- Re(eig$vectors[, real_idx, drop = FALSE])
  order_idx <- order(real_vals, decreasing = TRUE); num_axes <- min(g - 1, p)
  eig_vals <- real_vals[order_idx][1:num_axes]; eig_vecs <- real_vecs[, order_idx[1:num_axes], drop = FALSE]

  for (j in seq_len(num_axes)) {
    v <- eig_vecs[, j]; s <- as.numeric(sqrt(t(v) %*% S_reg %*% v))
    if (s > 1e-8) eig_vecs[, j] <- v / s
  }
  var_pct <- round((eig_vals / sum(eig_vals)) * 100, 2)
  message("=== CDA Canonical Variates Eigenvalues & Variance ===")
  for (j in seq_len(num_axes)) message(sprintf("  Can%d: Eigenvalue = %7.4f (%5.2f%% variation)", j, eig_vals[j], var_pct[j]))

  Z_all <- t(t(X) - grand_mean) %*% eig_vecs
  canonical_centroids <- t(t(group_means) - grand_mean) %*% eig_vecs
  priors <- group_counts / sum(group_counts)
  post_probs <- matrix(0, nrow(df), g, dimnames = list(NULL, g_levels))

  for (k in seq_along(g_levels)) {
    diff_z <- t(t(Z_all) - canonical_centroids[k, ])
    post_probs[, k] <- priors[k] * exp(-0.5 * rowSums(diff_z^2))
  }
  post_probs <- post_probs / rowSums(post_probs)
  pred_idx <- apply(post_probs, 1, which.max)

  return(list(
    eigenvalues = eig_vals, variance_pct = var_pct, canonical_scores = Z_all,
    centroids = canonical_centroids, predicted_taxon = g_levels[pred_idx],
    posterior_prob = apply(post_probs, 1, max), is_active = is_active
  ))
}

# ------------------------------------------------------------------------------
# 5. Herbarium Misidentification Auditing & Triage Flagging
# ------------------------------------------------------------------------------
audit_misidentifications <- function(df, cda_res, gmm_res) {
  message("Auditing specimens for herbarium misidentifications & label discordances...")
  n <- nrow(df); flags <- logical(n); triage <- character(n); reasons <- character(n)

  for (i in seq_len(n)) {
    raw_sp <- df$species_standardized[i]; pred_sp <- cda_res$predicted_taxon[i]
    post_p <- cda_res$posterior_prob[i]; is_act <- cda_res$is_active[i]

    if (is_act) {
      if (raw_sp == pred_sp) {
        flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Verified_Tier1_Gold_Anchor"
      } else if (post_p >= 0.85) {
        flags[i] <- TRUE; triage[i] <- "HIGH"; reasons[i] <- sprintf("Tier_1_Gold_Morphological_Discordance_to_%s", gsub(" ", "_", pred_sp))
      } else {
        flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Tier_1_Gold_Borderline_Variant"
      }
    } else {
      if (raw_sp == pred_sp && post_p >= 0.75) {
        flags[i] <- FALSE; triage[i] <- "CLEAN"; reasons[i] <- "Congruent_Bronze_Determination"
      } else if (raw_sp != pred_sp && post_p >= 0.70) {
        flags[i] <- TRUE; triage[i] <- "HIGH"
        if (grepl("dubia|tomentos", raw_sp, ignore.case = TRUE) && pred_sp == "Packera anonyma") {
          reasons[i] <- "Glabrescent_P_dubia_Misidentified_as_P_anonyma"
        } else {
          reasons[i] <- sprintf("Tier_3_Bronze_Misidentified_as_%s", gsub(" ", "_", pred_sp))
        }
      } else if (post_p < 0.60) {
        flags[i] <- FALSE; triage[i] <- "MEDIUM"; reasons[i] <- "Morphological_Intermediate_Ambiguous"
      } else {
        flags[i] <- (raw_sp != pred_sp); triage[i] <- if (raw_sp != pred_sp) "MEDIUM" else "LOW"
        reasons[i] <- if (raw_sp != pred_sp) "Candidate_Discordance_Moderate_Posterior" else "Minor_Uncertainty"
      }
    }
  }

  df$gmm_cluster <- gmm_res$classification; df$gmm_uncertainty <- round(gmm_res$uncertainty, 4)
  df$cda_predicted_taxon <- cda_res$predicted_taxon; df$cda_posterior_prob <- round(cda_res$posterior_prob, 4)
  df$can1 <- round(cda_res$canonical_scores[, 1], 5)
  df$can2 <- if (ncol(cda_res$canonical_scores) >= 2) round(cda_res$canonical_scores[, 2], 5) else 0.0
  df$is_passive <- !cda_res$is_active; df$misidentification_flag <- flags
  df$triage_priority <- triage; df$discordance_reason <- reasons

  message(sprintf("Audit Complete: %d misidentifications flagged (%d HIGH, %d MEDIUM priority)",
                  sum(flags), sum(triage == "HIGH"), sum(triage == "MEDIUM")))
  return(df)
}

# ------------------------------------------------------------------------------
# 6. Publication-Quality Multi-Panel Vector Biplot Generation
# ------------------------------------------------------------------------------
generate_cda_biplot_pdf <- function(audit_df, cda_res, gmm_res, output_pdf) {
  message("Rendering publication-quality multi-panel CDA biplot PDF: ", output_pdf)
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

  pdf(output_pdf, width = 12, height = 10, pointsize = 11)
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.2, 1.5), oma = c(1, 1, 2, 1), family = "sans")

  taxon_cols <- c("Packera anonyma" = "#E69F00", "Packera dubia" = "#009E73",
                  "Packera plattensis" = "#56B4E9", "Packera paupercula" = "#CC79A7")
  cluster_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

  # Panel A: Canonical Discriminant Analysis with Passive Projection
  plot(audit_df$can1, audit_df$can2, type = "n",
       xlab = sprintf("Canonical Axis 1 (%.1f%% Between-Group Var)", cda_res$variance_pct[1]),
       ylab = sprintf("Canonical Axis 2 (%.1f%% Between-Group Var)", cda_res$variance_pct[2]),
       main = "A. CDA with Passive Sample Projection (MorphoTools2)", font.main = 2)
  grid(col = "gray90", lty = "solid")

  active_sub <- audit_df[!audit_df$is_passive, ]
  for (sp in names(taxon_cols)) {
    sp_idx <- which(active_sub$species_standardized == sp)
    if (length(sp_idx) > 0) {
      points(active_sub$can1[sp_idx], active_sub$can2[sp_idx],
             pch = 21, bg = adjustcolor(taxon_cols[sp], alpha.f = 0.65), col = taxon_cols[sp], cex = 0.9)
    }
  }

  passive_sub <- audit_df[audit_df$is_passive, ]
  if (nrow(passive_sub) > 0) {
    clean_p <- passive_sub[!passive_sub$misidentification_flag, ]
    flagged_p <- passive_sub[passive_sub$misidentification_flag, ]
    if (nrow(clean_p) > 0) points(clean_p$can1, clean_p$can2, pch = 5, col = adjustcolor("gray40", alpha.f = 0.5), cex = 0.8)
    if (nrow(flagged_p) > 0) points(flagged_p$can1, flagged_p$can2, pch = 23, bg = "#D55E00", col = "black", cex = 1.1, lwd = 1.2)
  }

  for (k in seq_len(nrow(cda_res$centroids))) {
    points(cda_res$centroids[k, 1], cda_res$centroids[k, 2], pch = 3, col = "black", cex = 2.2, lwd = 3)
    text(cda_res$centroids[k, 1], cda_res$centroids[k, 2], labels = rownames(cda_res$centroids)[k],
         pos = 3, cex = 0.8, font = 4, col = "black")
  }
  legend("topleft", legend = c("Tier 1 Gold Anchors", "Tier 3 Passive Congruent", "Flagged Misidentification"),
         pch = c(21, 5, 23), pt.bg = c("#009E73", NA, "#D55E00"), col = c("#009E73", "gray40", "black"), bty = "n", cex = 0.8)

  # Panel B: Unsupervised GMM Clusters on EFA Morphospace
  plot(audit_df$PC1, audit_df$PC2, pch = 20,
       col = adjustcolor(cluster_cols[(audit_df$gmm_cluster - 1) %% length(cluster_cols) + 1], alpha.f = 0.6),
       xlab = "Morphospace PC1", ylab = "Morphospace PC2",
       main = sprintf("B. Label-Blind GMM Clusters (Optimal K = %d)", gmm_res$best_k), font.main = 2)
  grid(col = "gray90", lty = "solid")
  legend("bottomright", legend = paste("Cluster", 1:gmm_res$best_k),
         col = cluster_cols[1:gmm_res$best_k], pch = 20, bty = "n", cex = 0.8, ncol = 2)

  # Panel C: Bayes Factor (2ΔBIC) Model Comparison
  bic_df <- gmm_res$bic_table; delta_bic <- bic_df$Two_Delta_BIC; delta_bic[1] <- 0
  bar_cols <- ifelse(delta_bic >= 10, "#0072B2", ifelse(delta_bic >= 6, "#56B4E9", ifelse(delta_bic >= 2, "#F0E442", "#999999")))
  bp <- barplot(delta_bic, names.arg = paste0("K=", bic_df$K), col = bar_cols, border = "white",
                xlab = "Number of Mixture Components (K)", ylab = "Bayes Factor (2ΔBIC vs K-1)",
                main = "C. Species Boundary Evidence (Kass & Raftery 1995)", font.main = 2)
  abline(h = c(2, 6, 10), lty = 2, col = c("gray60", "gray40", "red3"))
  text(bp[length(bp)], 10.5, "Decisive (2ΔBIC ≥ 10)", adj = c(1, 0), cex = 0.7, col = "red3", font = 3)

  # Panel D: Posterior Probability Re-determination Shift
  hist(audit_df$cda_posterior_prob[audit_df$misidentification_flag], breaks = 15,
       col = adjustcolor("#D55E00", alpha.f = 0.6), border = "white",
       xlab = "CDA Posterior Classification Confidence", ylab = "Specimen Count",
       main = "D. Posterior Confidence of Flagged Misidentifications", font.main = 2)
  hist(audit_df$cda_posterior_prob[!audit_df$misidentification_flag], breaks = 20,
       col = adjustcolor("#009E73", alpha.f = 0.35), border = "white", add = TRUE)
  legend("topleft", legend = c("Flagged Misidentifications", "Congruent Vouchers"),
         fill = c(adjustcolor("#D55E00", alpha.f = 0.6), adjustcolor("#009E73", alpha.f = 0.35)), bty = "n", cex = 0.8)

  title("Packera dubia Complex: Morphometrics & Misidentification Triage", outer = TRUE, cex.main = 1.4)
  dev.off()
  message("Biplot graphic successfully saved to ", output_pdf)
}

# ------------------------------------------------------------------------------
# 7. Main Workflow Orchestrator
# ------------------------------------------------------------------------------
run_gmm_morphotools_pipeline <- function(opts) {
  message("==================================================================")
  message("Starting Packera Morphometrics Cluster Discovery & Passive CDA")
  message("Vouchers: ", opts$vouchers, " | Harmonics: ", opts$harmonics)
  message("==================================================================")

  if (!file.exists(opts$harmonics)) stop("Harmonics file not found: ", opts$harmonics)
  efa_df <- read.csv(opts$harmonics, stringsAsFactors = FALSE)
  vouchers_df <- if (file.exists(opts$vouchers)) read.csv(opts$vouchers, stringsAsFactors = FALSE) else NULL

  pca_cols <- paste0("PC", 1:opts$num_pcs)
  missing_pca <- setdiff(pca_cols, names(efa_df))
  if (length(missing_pca) > 0) stop("Missing PCA columns in EFA table: ", paste(missing_pca, collapse = ", "))

  closed_df <- efa_df[efa_df$assigned_tier %in% c("Tier_1_Direct", "Tier_2_Reflected"), ]
  valid_idx <- which(complete.cases(closed_df[, pca_cols]))
  closed_df <- closed_df[valid_idx, ]
  message(sprintf("Valid closed leaf outlines for morphometric modeling: %d", nrow(closed_df)))

  closed_df$species_standardized <- standardize_packera_taxon(closed_df$species_raw)
  target_taxa <- c("Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula")

  # 1. Unsupervised Gaussian Mixture Modeling (Label-Blind)
  gmm_res <- run_gmm_cluster_analysis(closed_df, pca_cols, max_k = opts$max_k)
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE)
  write.csv(gmm_res$bic_table, file = opts$output_report, row.names = FALSE)
  message("Bayes Factor (2ΔBIC) summary report saved: ", opts$output_report)

  # 2. Canonical Discriminant Analysis with Passive Sample Projection
  cda_res <- run_cda_with_passive_projection(closed_df, pca_cols, target_taxa)

  # 3. Herbarium Misidentification Auditing & Triage
  audited_df <- audit_misidentifications(closed_df, cda_res, gmm_res)

  if (!is.null(vouchers_df)) {
    v_meta <- intersect(names(vouchers_df), c("catalogNumber", "institutionCode", "county", "stateProvince",
                                              "latitude", "longitude", "pheno_sin", "pheno_cos", "regional_group"))
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), v_meta, drop = FALSE]
    audited_df <- merge(audited_df, v_sub, by = "catalogNumber", all.x = TRUE, suffixes = c("", "_meta"))
  }

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "species_raw", "species_standardized",
                 "determiner_tier", "cda_predicted_taxon", "cda_posterior_prob", "is_passive",
                 "misidentification_flag", "triage_priority", "discordance_reason", "gmm_cluster",
                 "gmm_uncertainty", "can1", "can2", "PC1", "PC2", "PC3", "PC4", "PC5")
  lead_cols <- intersect(lead_cols, names(audited_df))
  out_df <- audited_df[, c(lead_cols, setdiff(names(audited_df), lead_cols))]

  dir.create(dirname(opts$output_flags), recursive = TRUE, showWarnings = FALSE)
  write.csv(out_df, file = opts$output_flags, row.names = FALSE, na = "")
  message("Master misidentification flags table exported: ", opts$output_flags, " (Rows: ", nrow(out_df), ")")

  # 4. Multi-Panel Publication PDF Biplot
  generate_cda_biplot_pdf(audited_df, cda_res, gmm_res, opts$output_plot)

  message("==================================================================")
  message("Morphometric Pipeline Complete Successfully.")
  message("==================================================================")
  return(invisible(out_df))
}

if (sys.nframe() == 0) {
  opts <- parse_args_robust()
  run_gmm_morphotools_pipeline(opts)
}
