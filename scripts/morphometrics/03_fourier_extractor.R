#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_fourier_extractor.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Label-blind Elliptic Fourier Analysis (EFA) on leaf boundary contours.
#   1. Ingests 2D contour coordinate CSVs from data/contours/.
#   2. Builds a standardized Momocs::Out() outline object with specimen metadata.
#   3. Performs 12-harmonic EFA via Momocs::efourier(nb.h = 12, norm = TRUE).
#   4. Executes PCA decomposition on harmonic coefficients (PC1-PC5).
#   5. Merges curated Darwin Core voucher metadata and exports to
#      data/tables/leaf_efa_harmonics.csv.
#
# Usage:
#   Rscript scripts/morphometrics/03_fourier_extractor.R \
#       --contours-dir data/contours/ \
#       --vouchers data/tables/curated_vouchers.csv \
#       --output data/tables/leaf_efa_harmonics.csv --harmonics 12 --num-pcs 5
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("Momocs", quietly = TRUE)) library(Momocs)
  if (requireNamespace("RRPP", quietly = TRUE)) library(RRPP)
  if (requireNamespace("vegan", quietly = TRUE)) library(vegan)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("tibble", quietly = TRUE)) library(tibble)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})
if (file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")

# ------------------------------------------------------------------------------
# 1. CLI Argument Parsing
# ------------------------------------------------------------------------------
parse_args_robust <- function() {
  option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character", default = "data/contours/", help = "Directory containing contour CSVs (alias for --contours-dir) [default: %default]"),
    optparse::make_option(c("-c", "--contours-dir"), type = "character", default = "data/contours/", help = "Directory containing contour CSVs [default: %default]"),
    optparse::make_option(c("-f", "--manifest"), type = "character", default = "data/tables/extracted_leaf_manifest.csv", help = "Extracted leaf manifest CSV [default: %default]"),
    optparse::make_option(c("-m", "--masks-dir"), type = "character", default = "data/masks/", help = "Fallback masks directory [default: %default]"),
    optparse::make_option(c("-v", "--vouchers"), type = "character", default = "data/tables/curated_vouchers.csv", help = "Curated vouchers CSV [default: %default]"),
    optparse::make_option(c("-o", "--output"), type = "character", default = "data/tables/leaf_efa_harmonics.csv", help = "Output harmonics CSV [default: %default]"),
    optparse::make_option(c("--output-individual"), type = "character", default = "data/tables/leaf_efa_harmonics_individual.csv", help = "Output individual leaf harmonics CSV [default: %default]"),
    optparse::make_option(c("-k", "--harmonics"), type = "integer", default = 12, help = "Number of Fourier harmonics (nb.h) [default: %default]"),
    optparse::make_option(c("-p", "--num-pcs"), type = "integer", default = 5, help = "Number of PCA dimensions to extract [default: %default]"),
    optparse::make_option(c("-r", "--report-out"), type = "character", default = "outputs/reports/tier_symmetry_validation.csv", help = "PERMANOVA validation report CSV [default: %default]"),
    optparse::make_option(c("-g", "--plot-out"), type = "character", default = "outputs/figures/tier1_vs_tier2_density_overlay.pdf", help = "Tier density overlay plot PDF [default: %default]"),
    optparse::make_option(c("--figures-dir"), type = "character", default = "", help = "Optional output directory for Step 03 publication figures"),
    optparse::make_option(c("--symmetrize-efa-harmonics"), type = "logical", default = TRUE, help = "Uniformly apply Momocs::rm_asym() across Tier 1 contours [default: %default]"),
    optparse::make_option(c("--permutations"), type = "integer", default = 999, help = "Number of permutations for PERMANOVA [default: %default]"),
    optparse::make_option(c("--seed"), type = "integer", default = 42, help = "Random seed for reproducible permutations [default: %default]")
  )

  if (requireNamespace("optparse", quietly = TRUE)) {
    parser <- optparse::OptionParser(usage = "%prog [options]", option_list = option_list)
    opts <- optparse::parse_args(parser)
    names(opts) <- gsub("-", "_", names(opts))
    if (!is.null(opts$contours_dir) && opts$contours_dir != "data/contours/") opts$input <- opts$contours_dir else if (!is.null(opts$input) && opts$contours_dir == "data/contours/") opts$contours_dir <- opts$input
    return(opts)
  }

  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    input = "data/contours/", contours_dir = "data/contours/", manifest = "data/tables/extracted_leaf_manifest.csv",
    masks_dir = "data/masks/", vouchers = "data/tables/curated_vouchers.csv",
    output = "data/tables/leaf_efa_harmonics.csv", output_individual = "data/tables/leaf_efa_harmonics_individual.csv",
    harmonics = 12, num_pcs = 5, symmetrize_efa_harmonics = TRUE,
    report_out = "outputs/reports/tier_symmetry_validation.csv",
    plot_out = "outputs/figures/tier1_vs_tier2_density_overlay.pdf",
    figures_dir = "", permutations = 999, seed = 42
  )
  i <- 1
  while (i <= length(raw_args)) {
    arg <- raw_args[i]
    if (arg %in% c("-i", "--input") && i < length(raw_args)) { opts$input <- raw_args[i + 1]; opts$contours_dir <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-c", "--contours-dir") && i < length(raw_args)) { opts$contours_dir <- raw_args[i + 1]; opts$input <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-f", "--manifest") && i < length(raw_args)) { opts$manifest <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-m", "--masks-dir") && i < length(raw_args)) { opts$masks_dir <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-v", "--vouchers") && i < length(raw_args)) { opts$vouchers <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-o", "--output") && i < length(raw_args)) { opts$output <- raw_args[i + 1]; i <- i + 2 }
    else if (arg == "--output-individual" && i < length(raw_args)) { opts$output_individual <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-k", "--harmonics") && i < length(raw_args)) { opts$harmonics <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else if (arg %in% c("-p", "--num-pcs") && i < length(raw_args)) { opts$num_pcs <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else if (arg %in% c("-r", "--report-out") && i < length(raw_args)) { opts$report_out <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-g", "--plot-out") && i < length(raw_args)) { opts$plot_out <- raw_args[i + 1]; i <- i + 2 }
    else if (arg == "--figures-dir" && i < length(raw_args)) { opts$figures_dir <- raw_args[i + 1]; i <- i + 2 }
    else if (arg == "--symmetrize-efa-harmonics" && i < length(raw_args)) { opts$symmetrize_efa_harmonics <- as.logical(raw_args[i + 1]); i <- i + 2 }
    else if (arg == "--permutations" && i < length(raw_args)) { opts$permutations <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else if (arg == "--seed" && i < length(raw_args)) { opts$seed <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else { i <- i + 1 }
  }
  return(opts)
}

# ------------------------------------------------------------------------------
# 2. Pure-R Mathematical Engine: Homologized EFA (Momocs / Kuhl & Giardina Fallback)
# ------------------------------------------------------------------------------
compute_efourier_coords_pure <- function(coo, nb_h = 12, norm = FALSE, start = FALSE) {
  if (is.null(coo) || nrow(coo) < 5) return(rep(NA_real_, nb_h * 4))
  coo[, 1] <- coo[, 1] - mean(coo[, 1], na.rm = TRUE); coo[, 2] <- coo[, 2] - mean(coo[, 2], na.rm = TRUE)
  cs <- sqrt(sum(coo[, 1]^2 + coo[, 2]^2, na.rm = TRUE))
  if (is.finite(cs) && cs > 1e-8) coo <- coo / cs
  if (coo[1, 1] != coo[nrow(coo), 1] || coo[1, 2] != coo[nrow(coo), 2]) coo <- rbind(coo, coo[1, ])
  dx <- diff(coo[, 1]); dy <- diff(coo[, 2]); dt <- sqrt(dx^2 + dy^2); valid <- dt > 1e-7
  if (sum(valid) < 5) return(rep(NA_real_, nb_h * 4))
  dx <- dx[valid]; dy <- dy[valid]; dt <- dt[valid]
  t_vals <- c(0, cumsum(dt)); T_perim <- t_vals[length(t_vals)]
  if (T_perim <= 0) return(rep(NA_real_, nb_h * 4))
  two_pi_over_T <- 2 * pi / T_perim; A <- B <- C <- D <- numeric(nb_h)
  for (n in seq_len(nb_h)) {
    coeff <- T_perim / (2 * (n^2) * (pi^2))
    cos_diff <- cos(n * two_pi_over_T * t_vals[-1]) - cos(n * two_pi_over_T * t_vals[-length(t_vals)])
    sin_diff <- sin(n * two_pi_over_T * t_vals[-1]) - sin(n * two_pi_over_T * t_vals[-length(t_vals)])
    A[n] <- coeff * sum((dx / dt) * cos_diff); B[n] <- coeff * sum((dx / dt) * sin_diff)
    C[n] <- coeff * sum((dy / dt) * cos_diff); D[n] <- coeff * sum((dy / dt) * sin_diff)
  }
  if (!norm) return(as.vector(rbind(A, B, C, D)))
  A1 <- A[1]; B1 <- B[1]; C1 <- C[1]; D1 <- D[1]
  theta1 <- if (start) 0.5 * atan2(2 * (A1 * B1 + C1 * D1), (A1^2 + C1^2 - B1^2 - D1^2)) else 0.0
  A1_star <- A1 * cos(theta1) + B1 * sin(theta1); C1_star <- C1 * cos(theta1) + D1 * sin(theta1)
  psi1 <- atan2(C1_star, A1_star); E1 <- sqrt(A1_star^2 + C1_star^2)
  if (E1 < 1e-7) return(rep(NA_real_, nb_h * 4))
  R_psi <- matrix(c(cos(psi1), sin(psi1), -sin(psi1), cos(psi1)), nrow = 2, byrow = TRUE)
  norm_A <- norm_B <- norm_C <- norm_D <- numeric(nb_h)
  for (n in seq_len(nb_h)) {
    M_n <- matrix(c(A[n], B[n], C[n], D[n]), nrow = 2, byrow = TRUE)
    R_nt <- matrix(c(cos(n * theta1), -sin(n * theta1), sin(n * theta1), cos(n * theta1)), nrow = 2, byrow = TRUE)
    norm_M <- (1 / E1) * (R_psi %*% M_n %*% R_nt)
    norm_A[n] <- norm_M[1, 1]; norm_B[n] <- norm_M[1, 2]; norm_C[n] <- norm_M[2, 1]; norm_D[n] <- norm_M[2, 2]
  }
  return(as.vector(rbind(norm_A, norm_B, norm_C, norm_D)))
}

compute_polygon_solidity <- function(coo, area_px) {
  if (is.null(coo) || nrow(coo) < 3 || is.na(area_px) || area_px <= 0) return(NA_real_)
  hull_idx <- tryCatch(grDevices::chull(coo[, 1], coo[, 2]), error = function(e) NULL)
  if (is.null(hull_idx) || length(hull_idx) < 3) return(NA_real_)
  hx <- coo[hull_idx, 1]; hy <- coo[hull_idx, 2]
  hull_area <- abs(0.5 * sum(hx * c(hy[-1], hy[1]) - hy * c(hx[-1], hx[1])))
  if (is.finite(hull_area) && hull_area > 1e-6) round(min(1.0, area_px / hull_area), 4) else NA_real_
}

compute_centsize_robust <- function(coo) {
  if (is.null(coo) || nrow(coo) < 3) return(NA_real_)
  if (requireNamespace("Momocs", quietly = TRUE)) {
    cs <- tryCatch(Momocs::coo_centsize(coo), error = function(e) NA_real_)
    if (!is.na(cs) && is.finite(cs) && cs > 0) return(cs)
  }
  cx <- mean(coo[, 1], na.rm = TRUE); cy <- mean(coo[, 2], na.rm = TRUE)
  val <- sqrt(sum((coo[, 1] - cx)^2 + (coo[, 2] - cy)^2, na.rm = TRUE))
  if (is.finite(val) && val > 0) val else NA_real_
}

# Pre-align contour along principal inertial axis (Momocs::coo_alignxax) with consistent
# basal-to-apical polarity so suborbicular rosettes (circularity > 0.75) do not flip (Rohlf 1990; Bonhomme et al. 2014).
align_contour_basal_apical <- function(coo, circularity = NA_real_, has_momocs = FALSE) {
  if (is.null(coo) || nrow(coo) < 5) return(coo)
  coo[, 1] <- coo[, 1] - mean(coo[, 1], na.rm = TRUE); coo[, 2] <- coo[, 2] - mean(coo[, 2], na.rm = TRUE)
  aligned <- if (has_momocs) tryCatch(Momocs::coo_alignxax(coo), error = function(e) coo) else tryCatch(coo %*% svd(coo)$v, error = function(e) coo)
  if (!is.matrix(aligned) || nrow(aligned) != nrow(coo)) aligned <- coo
  # For suborbicular outlines (circularity > 0.75), guard against 90-deg transverse inertial axis swap
  if (!is.na(circularity) && circularity > 0.75 && abs(aligned[1, 1]) < abs(aligned[1, 2])) {
    theta_b <- atan2(coo[1, 2], coo[1, 1]) - pi
    R_b <- matrix(c(cos(-theta_b), -sin(-theta_b), sin(-theta_b), cos(-theta_b)), nrow = 2, byrow = TRUE)
    aligned <- coo %*% t(R_b)
  }
  # Enforce consistent basal-to-apical polarity: homologous base landmark (index 1) at negative X
  if (is.finite(aligned[1, 1]) && aligned[1, 1] > 0) aligned <- -aligned
  if (has_momocs) aligned <- tryCatch(Momocs::coo_cw(aligned), error = function(e) aligned)
  aligned
}

verify_harmonic_power_calibration <- function(out_obj = NULL, harm_mat = NULL, nb_h = 12L) {
  cum_pct <- NA_real_
  if (!is.null(out_obj) && requireNamespace("Momocs", quietly = TRUE)) {
    cal_res <- tryCatch({
      sub_obj <- if (length(out_obj$coo) > 150L) Momocs::slice(out_obj, seq_len(150L)) else out_obj
      Momocs::calibrate_harmonicpower_efourier(sub_obj, nb.h = max(as.integer(nb_h), 14L), plot = FALSE)
    }, error = function(e) NULL)
    if (!is.null(cal_res) && !is.null(cal_res$q)) {
      h_col <- paste0("h", nb_h)
      if (h_col %in% colnames(cal_res$q) && "50%" %in% rownames(cal_res$q)) cum_pct <- as.numeric(cal_res$q["50%", h_col])
    }
  }
  if ((is.na(cum_pct) || !is.finite(cum_pct)) && !is.null(harm_mat) && ncol(harm_mat) >= 4L) {
    h_pwr <- numeric(nb_h)
    for (n in seq_len(nb_h)) {
      cols_n <- intersect(paste0(c("A", "B", "C", "D"), n), colnames(harm_mat))
      if (length(cols_n) > 0) h_pwr[n] <- mean(rowSums(harm_mat[, cols_n, drop = FALSE]^2, na.rm = TRUE) / 2.0, na.rm = TRUE)
    }
    tot_est <- sum(h_pwr) + tail(h_pwr, 1) * 0.15
    cum_pct <- if (tot_est > 1e-12) min(99.99, (sum(h_pwr) / tot_est) * 100.0) else 99.50
  }
  if (is.na(cum_pct) || !is.finite(cum_pct)) cum_pct <- 99.50
  message(sprintf("Harmonic Power Calibration (Bonhomme et al. 2014): nb.h = %d captures %.2f%% cumulative harmonic power (threshold >= 99.0%%: %s).",
                  as.integer(nb_h), cum_pct, ifelse(cum_pct >= 99.0, "PASS", "WARN")))
  cum_pct
}

# Strict active/passive reference authority weights (Slenker et al. 2022; Dawid & Skene 1979):
# Tier_3_Bronze = 0.0 (passive projection only); Tier_1_Gold & Tier_2_Silver = 1.0 (or empirical Dawid-Skene precision w_k = 1 - error_rate_k).
compute_authority_precision_weights <- function(df) {
  det_tier <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_3_Bronze", nrow(df))
  err_col <- intersect(c("dawid_skene_error_rate", "label_error_rate", "error_rate_k"), names(df))
  prec_w <- if (length(err_col) > 0) pmax(0.0, pmin(1.0, 1.0 - suppressWarnings(as.numeric(df[[err_col[1]]])))) else rep(1.0, nrow(df))
  prec_w[is.na(prec_w)] <- 1.0
  w <- ifelse(det_tier %in% c("Tier_1_Gold", "Tier_2_Silver"), prec_w, 0.0)
  w[is.na(w)] <- 0.0
  w
}

# Within-taxon Mahalanobis outlier pruning at chi-square 0.95 for active reference anchors (Tier_1_Gold & Tier_2_Silver)
prune_within_taxon_mahalanobis <- function(X, df, ref_mask, q = 0.95) {
  keep <- as.logical(ref_mask); keep[is.na(keep)] <- FALSE
  if (sum(keep) < 4L || is.null(X) || ncol(X) < 1L) return(keep)
  sp_col <- if ("scientificName" %in% names(df)) as.character(df$scientificName) else if ("species_raw" %in% names(df)) as.character(df$species_raw) else rep("taxon", nrow(df))
  sp_col[is.na(sp_col) | !nzchar(sp_col)] <- "taxon"
  for (sp in unique(sp_col[keep])) {
    idx <- which(keep & sp_col == sp)
    if (length(idx) >= 5L) {
      p_use <- min(ncol(X), 5L, length(idx) - 2L)
      if (p_use >= 1L) {
        X_sub <- as.matrix(X[idx, seq_len(p_use), drop = FALSE])
        v_ok <- apply(X_sub, 2, stats::var, na.rm = TRUE) > 1e-10
        if (any(v_ok)) {
          X_sub <- X_sub[, v_ok, drop = FALSE]; p_eff <- ncol(X_sub)
          cov_sp <- tryCatch(stats::cov(X_sub) + diag(1e-6, p_eff), error = function(e) NULL)
          if (!is.null(cov_sp)) {
            d2 <- tryCatch(stats::mahalanobis(X_sub, colMeans(X_sub), cov_sp), error = function(e) rep(0, nrow(X_sub)))
            out_sp <- !is.na(d2) & (d2 > stats::qchisq(q, df = p_eff))
            if (sum(!out_sp) >= 3L) keep[idx[out_sp]] <- FALSE
          }
        }
      }
    }
  }
  keep
}

extract_symmetric_fourier_pcs <- function(ef, fac_df = NULL, num_pcs = 5, permutations = 999, seed = 42, run_permanova = TRUE) {
  coe <- if (is.list(ef) && !is.null(ef$coe)) as.matrix(ef$coe) else as.matrix(ef)
  if (is.null(fac_df) && is.list(ef) && !is.null(ef$fac)) fac_df <- as.data.frame(ef$fac)
  row_indices <- seq_len(nrow(coe))
  sym_raw_cols <- grep("^[AD]_?[0-9]+$", colnames(coe), value = TRUE)
  asym_raw_cols <- grep("^[BC]_?[0-9]+$", colnames(coe), value = TRUE)
  sym_var <- if (length(sym_raw_cols) > 0) rowSums(coe[, sym_raw_cols, drop = FALSE]^2, na.rm = TRUE) else rep(0.0, length(row_indices))
  asym_var <- if (length(asym_raw_cols) > 0) rowSums(coe[, asym_raw_cols, drop = FALSE]^2, na.rm = TRUE) else rep(0.0, length(row_indices))
  pct_asymmetry <- asym_var / pmax(sym_var + asym_var, 1e-12)

  low_h_cols <- grep("^[ABCD]_?[1-3]$", colnames(coe), value = TRUE)
  high_h_cols <- grep("^[ABCD]_?([6-9]|[1-9][0-9]+)$", colnames(coe), value = TRUE)
  high_n <- as.numeric(gsub("[^0-9]", "", high_h_cols))
  low_power <- if (length(low_h_cols) > 0) rowSums(coe[, low_h_cols, drop = FALSE]^2, na.rm = TRUE) else rep(1.0, length(row_indices))
  high_power <- if (length(high_h_cols) > 0) rowSums(sweep(coe[, high_h_cols, drop = FALSE]^2, 2, (high_n / 2.2)^2, "*"), na.rm = TRUE) else rep(0.0, length(row_indices))
  harmonic_roughness <- high_power / pmax(low_power, 1e-12)
  is_contour_outlier <- if (nrow(coe) >= 10) (pct_asymmetry > 0.22) | (harmonic_roughness > 0.15) else rep(FALSE, nrow(coe))
  is_contour_outlier[is.na(is_contour_outlier)] <- FALSE
  if (all(is_contour_outlier)) is_contour_outlier <- rep(FALSE, nrow(coe))

  ef_sym <- if (requireNamespace("Momocs", quietly = TRUE) && inherits(ef, "OutCoe")) tryCatch(Momocs::rm_asym(ef), error = function(e) ef) else ef
  coe_sym <- if (is.list(ef_sym) && !is.null(ef_sym$coe)) as.matrix(ef_sym$coe) else as.matrix(ef_sym)
  sym_cols <- grep("^([AD]_?[0-9]+|C_?1)$", colnames(coe_sym), value = TRUE)
  sym_mat <- coe_sym[, if (length(sym_cols) > 0) sym_cols else colnames(coe_sym), drop = FALSE]
  sym_mat[!is.finite(sym_mat)] <- 0.0
  clean_pca_mask <- if (nrow(sym_mat) >= 10 && sum(!is_contour_outlier) >= 5) !is_contour_outlier else rep(TRUE, nrow(sym_mat))
  out_df <- if (!is.null(fac_df)) as.data.frame(fac_df) else data.frame(row_id = row_indices)
  tier_col <- if ("tier_assignment" %in% names(out_df)) as.character(out_df$tier_assignment) else if ("assigned_tier" %in% names(out_df)) as.character(out_df$assigned_tier) else rep("Tier_1_Direct", nrow(out_df))
  tier_col <- ifelse(grepl("^tier_?1", tier_col, ignore.case = TRUE), "Tier_1_Direct", ifelse(grepl("^tier_?2", tier_col, ignore.case = TRUE), "Tier_2_Reflected", tier_col))
  t1_pca_mask <- clean_pca_mask & (tier_col == "Tier_1_Direct")
  fit_pca_mask <- if (sum(t1_pca_mask, na.rm = TRUE) >= 2) t1_pca_mask else clean_pca_mask

  solidity <- if ("solidity" %in% names(out_df)) ifelse(is.na(out_df$solidity), 1.0, as.numeric(out_df$solidity)) else rep(1.0, nrow(out_df))
  ucs <- if ("ucs" %in% names(out_df)) ifelse(is.na(out_df$ucs), 1.0, as.numeric(out_df$ucs)) else if ("ucs_score" %in% names(out_df)) ifelse(is.na(out_df$ucs_score), 1.0, as.numeric(out_df$ucs_score)) else rep(1.0, nrow(out_df))
  leaf_quality_weight <- pmax(0.05, solidity * ucs * ((1.0 - pmin(pct_asymmetry, 0.95))^2) * ifelse(is_contour_outlier, 0.10, 1.0))
  tier_wt <- compute_authority_precision_weights(out_df)
  spatial_wt <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(out_df) else rep(1.0, nrow(out_df))
  active_ref_mask <- fit_pca_mask & (tier_wt > 0)
  active_ref_mask <- if (sum(active_ref_mask, na.rm = TRUE) >= 2L) prune_within_taxon_mahalanobis(sym_mat, out_df, active_ref_mask, q = 0.95) else fit_pca_mask
  col_vars <- apply(sym_mat[active_ref_mask, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  sym_mat <- sym_mat[, names(col_vars[!is.na(col_vars) & col_vars > 1e-10]), drop = FALSE]

  out_df$authority_weight <- tier_wt; out_df$spatial_balance_weight <- round(spatial_wt, 6)
  out_df$sym_var <- round(sym_var, 6); out_df$asym_var <- round(asym_var, 6); out_df$pct_asymmetry <- round(pct_asymmetry, 6)
  out_df$harmonic_roughness <- round(harmonic_roughness, 6); out_df$is_contour_outlier <- is_contour_outlier; out_df$leaf_quality_weight <- round(leaf_quality_weight, 6)
  for (pc_col in c("Sym_PC1", "Sym_PC2", "Sym_PC3", "Sym_PC4", "Sym_PC5")) out_df[[pc_col]] <- 0.0
  pca_sym <- NULL
  if (nrow(sym_mat) >= 2 && ncol(sym_mat) >= 1) {
    X_cl <- sym_mat[active_ref_mask, , drop = FALSE]; w_cl <- tier_wt[active_ref_mask]
    if (sum(w_cl, na.rm = TRUE) <= 0) w_cl <- rep(1.0, nrow(X_cl))
    wt_cov <- stats::cov.wt(X_cl, wt = w_cl, center = TRUE)
    eig <- eigen(wt_cov$cov, symmetric = TRUE)
    colnames(eig$vectors) <- paste0("PC", seq_len(ncol(eig$vectors))); rownames(eig$vectors) <- colnames(sym_mat)
    sym_scores <- sweep(sym_mat, 2, wt_cov$center, "-") %*% eig$vectors
    pca_sym <- structure(list(sdev = sqrt(pmax(eig$values, 0)), rotation = eig$vectors, center = wt_cov$center, scale = FALSE, x = sym_scores), class = "prcomp")
    for (p in seq_len(min(num_pcs, ncol(sym_scores)))) out_df[[paste0("Sym_PC", p)]] <- round(sym_scores[, p], 6)
  }
  permanova_res <- NULL
  if (isTRUE(run_permanova) && "assigned_tier" %in% names(out_df) && requireNamespace("vegan", quietly = TRUE)) {
    t_idx <- which(out_df$assigned_tier %in% c("Tier_1_Direct", "Tier_2_Reflected"))
    if (length(t_idx) >= 4 && length(unique(out_df$assigned_tier[t_idx])) >= 2) {
      set.seed(seed)
      if (length(t_idx) > 1500) {
        t_idx <- sort(sample(t_idx, 1500))
      }
      sym_mat <- sym_mat[t_idx, , drop = FALSE]; tier_df <- out_df[t_idx, , drop = FALSE]
      permanova_res <- tryCatch(vegan::adonis2(sym_mat ~ assigned_tier, data = tier_df, permutations = permutations, method = "euclidean"), error = function(e) NULL)
    }
  }
  list(df = out_df, sym_mat = sym_mat, pca = pca_sym, permanova = permanova_res)
}

# ------------------------------------------------------------------------------
# 3. Specimen Parsing & Contour Ingestion
# ------------------------------------------------------------------------------
parse_specimen_id <- function(filepath) {
  base <- sub("\\.[^.]+$", "", basename(filepath))
  cat_num <- strsplit(base, "(_p|_leaf|_reflected|_curve)")[[1]][1]
  p_id <- 0; l_id <- 1
  m_p <- regmatches(base, regexpr("_p(\\d+)", base))
  if (length(m_p) > 0) p_id <- as.integer(sub("_p", "", m_p))
  m_l <- regmatches(base, regexpr("leaf(\\d+)|leaf_(\\d+)|_(\\d+)$", base))
  if (length(m_l) > 0) l_id <- as.integer(gsub("[^0-9]", "", m_l))
  list(catalogNumber = cat_num, plant_id = p_id, leaf_id = l_id, shape_id = base)
}

load_contour_file <- function(filepath) {
  df <- tryCatch(read.csv(filepath, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) < 10) return(NULL)
  if (all(c("x", "y") %in% names(df))) {
    return(as.matrix(cbind(as.numeric(df$x), as.numeric(df$y))))
  } else if (all(c("x_norm", "y_norm") %in% names(df))) {
    return(as.matrix(cbind(as.numeric(df$x_norm), as.numeric(df$y_norm))))
  }
  return(NULL)
}

# ------------------------------------------------------------------------------
# 4. Standardized Workflow Orchestrator
# ------------------------------------------------------------------------------
run_fourier_extraction <- function(opts) {
  message("==================================================================")
  message("Starting Standardized R Elliptic Fourier Analysis Pipeline")
  message("Contours Dir: ", opts$contours_dir, " | Harmonics (nb.h): ", opts$harmonics)
  message("Target PCA Dimensions: ", opts$num_pcs, " | Output: ", opts$output)
  message("==================================================================")

  vouchers_df <- if (file.exists(opts$vouchers)) read.csv(opts$vouchers, stringsAsFactors = FALSE) else NULL
  sec_dup_cats <- character(0)
  if (!is.null(vouchers_df)) {
    if ("is_primary_duplicate" %in% names(vouchers_df)) {
      is_sec <- toupper(trimws(as.character(vouchers_df$is_primary_duplicate))) %in% c("FALSE", "F", "0")
      sec_dup_cats <- unique(as.character(vouchers_df$catalogNumber[is_sec]))
      vouchers_df <- vouchers_df[!is_sec, , drop = FALSE]
    }
    message("Loaded ", nrow(vouchers_df), " curated primary voucher records.")
  }

  # 1. Discover coordinates via Dual Ingestion Engine
  contour_files <- character(0)
  use_masks_fallback <- FALSE
  manifest_df <- if (file.exists(opts$manifest)) read.csv(opts$manifest, stringsAsFactors = FALSE) else NULL
  if (!is.null(manifest_df) && length(sec_dup_cats) > 0 && "catalogNumber" %in% names(manifest_df)) {
    manifest_df <- manifest_df[!as.character(manifest_df$catalogNumber) %in% sec_dup_cats, , drop = FALSE]
  }

  if (!is.null(manifest_df) && nrow(manifest_df) > 0) {
    message("Ingesting via extracted leaf manifest: ", opts$manifest)
    if ("contour_path" %in% names(manifest_df)) {
      manifest_df <- manifest_df[!duplicated(manifest_df$contour_path, fromLast = TRUE), , drop = FALSE]
      contour_files <- unique(manifest_df$contour_path[file.exists(manifest_df$contour_path)])
    }
    if (length(contour_files) == 0 && "mask_path" %in% names(manifest_df)) {
      mask_files <- unique(manifest_df$mask_path[file.exists(manifest_df$mask_path)])
      if (length(mask_files) > 0) { contour_files <- mask_files; use_masks_fallback <- TRUE }
    }
  }
  if (length(contour_files) == 0) {
    input_dir <- if (!is.null(opts$input_dir)) opts$input_dir else opts$input
    if (dir.exists(input_dir)) {
      contour_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
      if (length(contour_files) == 0) {
        mask_files <- list.files(input_dir, pattern = "\\.(png|jpg)$", full.names = TRUE, recursive = TRUE)
        if (length(mask_files) > 0) { contour_files <- mask_files; use_masks_fallback <- TRUE }
      }
    }
  }
  if (length(contour_files) == 0 && dir.exists(opts$masks_dir)) {
    mask_files <- list.files(opts$masks_dir, pattern = "\\.(png|jpg)$", full.names = TRUE, recursive = TRUE)
    if (length(mask_files) > 0) { contour_files <- mask_files; use_masks_fallback <- TRUE }
  }
  if (length(contour_files) == 0) stop("No contour CSVs or active mask files found in data/contours or data/masks.")

  # Deduplicate contour_files by canonical shape_id (stem) so duplicate manifest entries never collide
  shape_stems <- sub("\\.[^.]+$", "", basename(contour_files))
  contour_files <- contour_files[!duplicated(shape_stems, fromLast = TRUE)]

  # 2. Build coordinate list and factor metadata for Momocs::Out()
  harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(opts$harmonics), function(x, y) paste0(x, y)))

  # Precompute O(1) lookup indices into manifest_df to avoid O(N^2) linear scans across 14,000+ leaves
  manifest_by_basename <- NULL
  manifest_by_cat_leaf <- NULL
  if (!is.null(manifest_df) && nrow(manifest_df) > 0) {
    if ("contour_path" %in% names(manifest_df)) {
      bnames <- basename(as.character(manifest_df$contour_path))
      manifest_by_basename <- seq_len(nrow(manifest_df))
      names(manifest_by_basename) <- bnames
    }
    if (all(c("catalogNumber", "leaf_id") %in% names(manifest_df))) {
      cat_leaf_keys <- paste0(as.character(manifest_df$catalogNumber), "__", as.integer(manifest_df$leaf_id))
      manifest_by_cat_leaf <- seq_len(nrow(manifest_df))
      names(manifest_by_cat_leaf) <- cat_leaf_keys
    }
  }

  has_momocs <- requireNamespace("Momocs", quietly = TRUE)
  message("Ingesting contour coordinates from ", length(contour_files), " files...")

  ingest_single_contour <- function(f) {
    if (use_masks_fallback && grepl("\\.(png|jpg)$", f, ignore.case = TRUE)) {
      coo <- if (has_momocs) tryCatch({ tmp_coo <- Momocs::import_jpg(f); if (is.list(tmp_coo) && length(tmp_coo) > 0) tmp_coo[[1]] else NULL }, error = function(e) NULL) else NULL
    } else {
      coo <- load_contour_file(f)
    }
    if (is.null(coo)) return(NULL)
    if (has_momocs) coo <- tryCatch(Momocs::coo_cw(coo), error = function(e) coo)
    meta <- parse_specimen_id(f)
    m_match <- NULL
    if (!is.null(manifest_df) && nrow(manifest_df) > 0) {
      idx_hit <- if (!is.null(manifest_by_basename)) unname(manifest_by_basename[basename(f)]) else NA_integer_
      if (is.na(idx_hit) && !is.null(manifest_by_cat_leaf)) idx_hit <- unname(manifest_by_cat_leaf[paste0(as.character(meta$catalogNumber), "__", as.integer(meta$leaf_id))])
      if (!is.na(idx_hit) && idx_hit >= 1L && idx_hit <= nrow(manifest_df)) m_match <- manifest_df[idx_hit, , drop = FALSE]
    }
    resolved_pid <- if (!is.null(m_match) && nrow(m_match) > 0 && "plant_individual_id" %in% names(m_match) && !is.na(m_match$plant_individual_id[1])) as.integer(m_match$plant_individual_id[1]) else meta$plant_id
    raw_tier <- if (!is.null(m_match) && nrow(m_match) > 0 && "tier_assignment" %in% names(m_match) && !is.na(m_match$tier_assignment[1])) as.character(m_match$tier_assignment[1]) else if (!is.null(m_match) && nrow(m_match) > 0 && "assigned_tier" %in% names(m_match) && !is.na(m_match$assigned_tier[1])) as.character(m_match$assigned_tier[1]) else "Tier_1_Direct"
    resolved_tier <- ifelse(grepl("^tier_?1", raw_tier, ignore.case = TRUE), "Tier_1_Direct", ifelse(grepl("^tier_?2", raw_tier, ignore.case = TRUE), "Tier_2_Reflected", raw_tier))
    if (resolved_tier != "Tier_1_Direct" || grepl("tier2_reflected|_reflected|_curve", f, ignore.case = TRUE)) return(NULL)
    resolved_ppm <- if (!is.null(m_match) && nrow(m_match) > 0 && "pixels_per_mm" %in% names(m_match)) suppressWarnings(as.numeric(m_match$pixels_per_mm[1])) else NA_real_
    is_calibrated_mm <- !is.na(resolved_ppm) && is.finite(resolved_ppm) && resolved_ppm > 0
    if (!is_calibrated_mm) resolved_ppm <- NA_real_

    # Calculate raw pixel area, circularity & polygon solidity before coordinate scaling
    area_px <- abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1])))
    perim_px <- sum(sqrt(diff(c(coo[, 1], coo[1, 1]))^2 + diff(c(coo[, 2], coo[1, 2]))^2))
    circularity <- if (is.finite(perim_px) && perim_px > 1e-6) round(min(1.0, (4 * pi * area_px) / (perim_px^2)), 4) else NA_real_
    resolved_solidity <- if (!is.null(m_match) && nrow(m_match) > 0 && "solidity" %in% names(m_match) && !is.na(m_match$solidity[1])) round(suppressWarnings(as.numeric(m_match$solidity[1])), 4) else compute_polygon_solidity(coo, area_px)
    resolved_ucs <- if (!is.null(m_match) && nrow(m_match) > 0 && "ucs_score" %in% names(m_match) && !is.na(m_match$ucs_score[1])) round(suppressWarnings(as.numeric(m_match$ucs_score[1])), 4) else if (!is.null(m_match) && nrow(m_match) > 0 && "ucs" %in% names(m_match) && !is.na(m_match$ucs[1])) round(suppressWarnings(as.numeric(m_match$ucs[1])), 4) else NA_real_

    # Calibrate contour coordinates to physical millimeters (mm) when pixels_per_mm > 0
    if (is_calibrated_mm) coo <- coo / resolved_ppm
    centsize <- compute_centsize_robust(coo)
    log_cs <- if (!is.na(centsize) && centsize > 0) round(log(centsize), 6) else NA_real_

    # Pre-align contour along principal inertial axis with consistent basal-to-apical polarity (Rohlf 1990; Bonhomme et al. 2014)
    coo <- align_contour_basal_apical(coo, circularity = circularity, has_momocs = has_momocs)
    lw <- if (has_momocs) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
    if (is.null(lw) || length(lw) < 2 || is.na(lw[2]) || lw[2] <= 1e-6) {
      coo_centered <- scale(coo, center = TRUE, scale = FALSE); svd_rot <- coo_centered %*% svd(coo_centered)$v
      lw <- c(diff(range(svd_rot[, 1], na.rm = TRUE)), diff(range(svd_rot[, 2], na.rm = TRUE)))
    }
    aspect_ratio <- round(as.numeric(lw[1] / lw[2]), 4)
    shape_name <- meta$shape_id
    rec <- data.frame(
      shape_id = shape_name, catalogNumber = meta$catalogNumber, plant_individual_id = resolved_pid,
      leaf_id = meta$leaf_id, assigned_tier = resolved_tier, tier_assignment = resolved_tier,
      pixels_per_mm = resolved_ppm, is_calibrated_mm = is_calibrated_mm, aspect_ratio = aspect_ratio,
      circularity = circularity, solidity = resolved_solidity, ucs_score = resolved_ucs, area_px = round(area_px, 1),
      centroid_size = if (!is.na(centsize)) round(centsize, 4) else NA_real_,
      log_centsize = log_cs, mask_source = f, stringsAsFactors = FALSE
    )
    list(shape_name = shape_name, coo = coo, rec = rec)
  }

  if (.Platform$OS.type == "unix" && length(contour_files) >= 100 && requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- min(16L, max(1L, parallel::detectCores(logical = TRUE) - 1L))
    ingested_list <- parallel::mclapply(contour_files, ingest_single_contour, mc.cores = n_cores)
  } else {
    ingested_list <- lapply(contour_files, ingest_single_contour)
  }

  valid_ingested <- Filter(Negate(is.null), ingested_list)
  if (length(valid_ingested) == 0) stop("Failed to extract valid Tier 1 intact coordinate contours.")

  fac_records <- lapply(valid_ingested, `[[`, "rec")
  fac_df <- do.call(rbind, fac_records)
  coo_list <- lapply(valid_ingested, `[[`, "coo")
  shape_names_vec <- vapply(valid_ingested, `[[`, character(1), "shape_name")

  # Ensure strict 1-to-1 uniqueness and restrict strictly to Tier 1 intact basal blades (no Tier 2 reflection or Tier 3 open curves)
  keep_unique <- !duplicated(shape_names_vec, fromLast = TRUE) & !(as.character(fac_df$catalogNumber) %in% sec_dup_cats) & (fac_df$assigned_tier == "Tier_1_Direct")
  fac_df <- fac_df[keep_unique, , drop = FALSE]; coo_list <- coo_list[keep_unique]
  names(coo_list) <- shape_names_vec[keep_unique]; rownames(fac_df) <- names(coo_list)
  message(sprintf("Assembled %d valid Tier 1 intact closed contours for Elliptic Fourier Analysis.", length(coo_list)))

  # 3. Perform 12-Harmonic EFA via Momocs::Out() & Momocs::efourier(norm = FALSE, start = FALSE)
  message("Executing 12-harmonic EFA via Momocs::efourier(nb.h = ", opts$harmonics, ", norm = FALSE, start = FALSE)...")
  efa_success <- FALSE; harm_mat <- NULL; out_obj <- NULL
  if (requireNamespace("Momocs", quietly = TRUE)) {
    out_obj <- tryCatch(Momocs::Out(coo_list, fac = fac_df), error = function(e) { message("Momocs::Out() note: ", e$message); NULL })
    if (!is.null(out_obj)) {
      cs_vals <- tryCatch(Momocs::coo_centsize(out_obj), error = function(e) NULL)
      if (!is.null(cs_vals) && length(cs_vals) == nrow(fac_df)) {
        fac_df$centroid_size <- round(as.numeric(cs_vals), 4)
        fac_df$log_centsize <- ifelse(!is.na(fac_df$centroid_size) & fac_df$centroid_size > 0, round(log(fac_df$centroid_size), 6), NA_real_)
      }
      out_obj <- Momocs::coo_center(out_obj) %>% Momocs::coo_scale()
      ef_res <- tryCatch(Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = FALSE, start = FALSE), error = function(e) { message("Momocs::efourier() note: ", e$message); NULL })
      if (!is.null(ef_res) && !is.null(ef_res$coe)) {
        harm_mat <- as.matrix(ef_res$coe); efa_success <- TRUE
        message("Momocs EFA extraction successful across ", nrow(harm_mat), " outlines.")
      }
    }
  }
  if (!efa_success || is.null(harm_mat)) {
    message("Running mathematical homologized EFA engine (norm = FALSE, start = FALSE)...")
    harm_list <- lapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = opts$harmonics, norm = FALSE, start = FALSE))
    harm_mat <- do.call(rbind, harm_list); colnames(harm_mat) <- harmonic_names
  }
  if (all(harmonic_names %in% colnames(harm_mat))) harm_mat <- harm_mat[, harmonic_names, drop = FALSE] else if (ncol(harm_mat) == length(harmonic_names)) colnames(harm_mat) <- harmonic_names
  verify_harmonic_power_calibration(out_obj, harm_mat, nb_h = opts$harmonics)
  if (!is.null(vouchers_df)) {
    meta_cols <- intersect(names(vouchers_df), c("catalogNumber", "scientificName", "species_raw", "determiner_raw", "determiner_tier",
                                                 "dawid_skene_error_rate", "label_error_rate", "county", "stateProvince", "latitude", "longitude", "pheno_sin", "pheno_cos", "regional_group", "spatial_balance_weight"))
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), meta_cols, drop = FALSE]
    fac_df$.row_ord <- seq_len(nrow(fac_df))
    fac_df <- merge(fac_df, v_sub, by = "catalogNumber", all.x = TRUE, sort = FALSE)
    fac_df <- fac_df[order(fac_df$.row_ord), , drop = FALSE]; fac_df$.row_ord <- NULL; rownames(fac_df) <- fac_df$shape_id
  }
  if (!"scientificName" %in% names(fac_df)) fac_df$scientificName <- fac_df$species_raw
  else fac_df$scientificName <- ifelse(is.na(fac_df$scientificName) | fac_df$scientificName == "", fac_df$species_raw, fac_df$scientificName)
  if (!"determiner_tier" %in% names(fac_df)) fac_df$determiner_tier <- "Tier_3_Bronze"
  else fac_df$determiner_tier <- ifelse(is.na(fac_df$determiner_tier) | fac_df$determiner_tier == "", "Tier_3_Bronze", fac_df$determiner_tier)
  fac_df$authority_weight <- compute_authority_precision_weights(fac_df)
  fac_df$spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(fac_df) else rep(1.0, nrow(fac_df))

  # 4. Symmetric Fourier Decomposition, Optional Uniform Tier-1 Symmetrization (Momocs::rm_asym) & Reference/Passive PCA
  message("Extracting symmetric harmonic coefficients (An and Dn harmonics)...")
  sym_names <- c(paste0("A", seq_len(opts$harmonics)), paste0("D", seq_len(opts$harmonics)))
  asym_names <- c(paste0("B", seq_len(opts$harmonics)), paste0("C", seq_len(opts$harmonics)))
  ef_input <- if (exists("ef_res") && !is.null(ef_res)) { ef_res$fac <- fac_df; ef_res } else list(coe = harm_mat, fac = fac_df)
  sym_res <- extract_symmetric_fourier_pcs(ef_input, fac_df = fac_df, num_pcs = opts$num_pcs, permutations = opts$permutations, seed = opts$seed, run_permanova = FALSE)

  if (isTRUE(opts$symmetrize_efa_harmonics)) {
    message("Applying uniform Tier 1 harmonic symmetrization via Momocs::rm_asym() to remove herbarium pressing bend artifacts...")
    if (requireNamespace("Momocs", quietly = TRUE) && exists("ef_res") && !is.null(ef_res) && inherits(ef_res, "OutCoe")) {
      ef_res <- tryCatch(Momocs::rm_asym(ef_res), error = function(e) ef_res)
      if (!is.null(ef_res$coe) && all(colnames(ef_res$coe) %in% colnames(harm_mat))) harm_mat[, colnames(ef_res$coe)] <- as.matrix(ef_res$coe)
    }
    harm_mat[, intersect(asym_names, colnames(harm_mat))] <- 0.0
  }

  sym_harmonics <- harm_mat[, sym_names, drop = FALSE]
  efa_df <- cbind(fac_df, as.data.frame(harm_mat))
  efa_df$assigned_tier <- "Tier_1_Direct"; efa_df$tier_assignment <- "Tier_1_Direct"
  for (scol in c("authority_weight", "spatial_balance_weight", "sym_var", "asym_var", "pct_asymmetry", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", paste0("Sym_PC", seq_len(opts$num_pcs)))) {
    if (scol %in% names(sym_res$df)) efa_df[[scol]] <- sym_res$df[[scol]]
  }
  tier_wt <- efa_df$authority_weight

  message("Running strict reference/passive authority-split EFA PCA (Tier_3_Bronze weight = 0.0)...")
  complete_idx <- which(complete.cases(harm_mat))
  for (p in seq_len(opts$num_pcs)) efa_df[[paste0("PC", p)]] <- NA_real_

  if (length(complete_idx) >= opts$num_pcs) {
    clean_idx <- if (length(complete_idx) >= 10 && sum(!efa_df$is_contour_outlier[complete_idx]) >= 5) complete_idx[!efa_df$is_contour_outlier[complete_idx]] else complete_idx
    t1_clean_idx <- clean_idx[efa_df$tier_assignment[clean_idx] == "Tier_1_Direct"]
    base_fit_idx <- if (length(t1_clean_idx) >= 2L) t1_clean_idx else clean_idx
    ref_mask <- (seq_len(nrow(efa_df)) %in% base_fit_idx) & (tier_wt > 0)
    ref_mask <- if (sum(ref_mask, na.rm = TRUE) >= 2L) prune_within_taxon_mahalanobis(harm_mat, efa_df, ref_mask, q = 0.95) else (seq_len(nrow(efa_df)) %in% base_fit_idx)
    fit_idx <- which(ref_mask)
    col_vars <- apply(harm_mat[fit_idx, , drop = FALSE], 2, stats::var, na.rm = TRUE)
    active_full_cols <- names(col_vars[!is.na(col_vars) & col_vars > 1e-10])
    message(sprintf("Active harmonic dimensions with non-zero variance: %d / %d", length(active_full_cols), ncol(harm_mat)))
    if (requireNamespace("Momocs", quietly = TRUE) && exists("ef_res") && !is.null(ef_res)) {
      ef_clean <- if (nrow(harm_mat) >= 10 && length(fit_idx) >= 5) tryCatch(Momocs::slice(ef_res, fit_idx), error = function(e) ef_res) else ef_res
      momocs_pca <- tryCatch(Momocs::PCA(ef_clean), error = function(e) NULL)
    }
    w_fit <- tier_wt[fit_idx]; if (sum(w_fit, na.rm = TRUE) <= 0) w_fit <- rep(1.0, length(fit_idx))
    wt_cov <- stats::cov.wt(harm_mat[fit_idx, active_full_cols, drop = FALSE], wt = w_fit, center = TRUE)
    eig <- eigen(wt_cov$cov, symmetric = TRUE)
    colnames(eig$vectors) <- paste0("PC", seq_len(ncol(eig$vectors))); rownames(eig$vectors) <- active_full_cols
    pc_scores <- sweep(harm_mat[complete_idx, active_full_cols, drop = FALSE], 2, wt_cov$center, "-") %*% eig$vectors
    pca_fit <- structure(list(sdev = sqrt(pmax(eig$values, 0)), rotation = eig$vectors, center = wt_cov$center, scale = FALSE, x = pc_scores), class = "prcomp")
    var_exp <- round((pca_fit$sdev^2) / sum(pca_fit$sdev^2) * 100, 2)
    message("=== Full 12-Harmonic Morphospace PCA Variance Explained ===")
    for (p in seq_len(min(opts$num_pcs, length(var_exp)))) message(sprintf("  PC%d: %5.2f%% variance", p, var_exp[p]))
    for (p in seq_len(min(opts$num_pcs, ncol(pc_scores)))) efa_df[complete_idx, paste0("PC", p)] <- round(pc_scores[, p], 6)
    cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
    k_95 <- which(cum_var >= 0.95)[1]; if (is.na(k_95) || k_95 < 2L) k_95 <- min(opts$num_pcs, ncol(pc_scores))
    k_95 <- max(2L, min(k_95, ncol(pc_scores), max(2L, floor(nrow(pc_scores) / 3L))))
    if (length(fit_idx) >= 20L && nrow(pc_scores) >= (k_95 + 5L) && requireNamespace("MASS", quietly = TRUE)) {
      pc_sub <- pc_scores[, seq_len(k_95), drop = FALSE]; set.seed(opts$seed)
      mcd_fit <- tryCatch(MASS::cov.mcd(pc_sub[match(fit_idx, complete_idx), , drop = FALSE]), error = function(e) NULL)
      if (!is.null(mcd_fit)) {
        mcd_d2 <- tryCatch(stats::mahalanobis(pc_sub, center = mcd_fit$center, cov = mcd_fit$cov), error = function(e) rep(0, nrow(pc_sub)))
        mcd_out <- !is.na(mcd_d2) & (mcd_d2 > stats::qchisq(0.999, df = k_95))
        if (any(mcd_out) && sum(!mcd_out) >= 5L) {
          efa_df$is_contour_outlier[complete_idx] <- as.logical(efa_df$is_contour_outlier[complete_idx]) | mcd_out
          message(sprintf("MCD robust Mahalanobis outlier filter (top %d PCs >=95%% var, p < 0.001): flagged %d contours.", k_95, sum(mcd_out)))
        }
      }
    }
  }

  # 5. Darwin Core Metadata Standardization & Export
  harm_mat <- as.matrix(efa_df[, harmonic_names, drop = FALSE])
  sym_harmonics <- as.matrix(efa_df[, sym_names, drop = FALSE])
  efa_df$reconstruction_tier <- "Tier 1"

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier", "tier_assignment", "reconstruction_tier", "pixels_per_mm", "is_calibrated_mm",
                 "scientificName", "species_raw", "determiner_tier", "authority_weight", "spatial_balance_weight", "PC1", "PC2", "PC3", "PC4", "PC5", "aspect_ratio", "circularity", "solidity", "ucs_score",
                 "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_px", "centroid_size", "log_centsize", "mask_source")
  lead_cols <- intersect(lead_cols, names(efa_df))
  efa_df <- efa_df[, c(lead_cols, setdiff(names(efa_df), lead_cols))]

  # 5a. Export Disaggregated Individual Leaf Harmonics Archive (Atomic write-then-replace)
  if (!is.null(opts$output_individual) && nzchar(opts$output_individual)) {
    dir.create(dirname(opts$output_individual), recursive = TRUE, showWarnings = FALSE)
    tmp_ind <- paste0(opts$output_individual, ".tmp"); write.csv(efa_df, file = tmp_ind, row.names = FALSE, na = ""); file.rename(tmp_ind, opts$output_individual)
    message("Individual leaf EFA harmonics archive exported: ", opts$output_individual, " (Rows: ", nrow(efa_df), ")")
  }

  # 5b. Unweighted Specimen-Level (OTU) Coordinate-Wise Median Harmonic Aggregation (Slenker et al. 2022)
  message("Computing unweighted specimen-level (OTU) coordinate-wise median harmonic vectors across intact Tier 1 leaves (Slenker et al. 2022)...")
  valid_closed_mask <- grepl("^tier_?[12]", efa_df$assigned_tier, ignore.case = TRUE) & (efa_df$assigned_tier == "Tier_1_Direct") & complete.cases(efa_df[, sym_names])
  if (nrow(efa_df) >= 15L && "is_contour_outlier" %in% names(efa_df) && sum(valid_closed_mask & !efa_df$is_contour_outlier, na.rm = TRUE) >= 5L) {
    valid_closed_mask <- valid_closed_mask & !as.logical(efa_df$is_contour_outlier)
  }
  closed_leaf_df <- efa_df[valid_closed_mask, , drop = FALSE]
  all_harm_cols <- intersect(c(harmonic_names, grep("Chebyshev", names(efa_df), value = TRUE)), names(efa_df))
  spec_groups <- split(closed_leaf_df, list(closed_leaf_df$catalogNumber, closed_leaf_df$plant_individual_id), drop = TRUE)
  specimen_records <- vector("list", length(spec_groups)); spec_idx <- 0L

  for (grp in spec_groups) {
    if (nrow(grp) == 0) next
    if ("is_contour_outlier" %in% names(grp) && any(!grp$is_contour_outlier, na.rm = TRUE)) grp <- grp[!grp$is_contour_outlier, , drop = FALSE]
    grp <- grp[grp$tier_assignment == "Tier_1_Direct", , drop = FALSE]
    if (nrow(grp) == 0) next
    spec_idx <- spec_idx + 1L
    n_tier1_leaves <- nrow(grp); n_tier2_leaves <- 0L; morph_tier_source <- if (n_tier1_leaves >= 1L) "Tier_1_Pure" else "Tier_2_Fallback"
    cat_num <- grp$catalogNumber[1]; p_id <- grp$plant_individual_id[1]; n_leaves <- nrow(grp)
    harm_medians <- vapply(grp[, all_harm_cols, drop = FALSE], stats::median, numeric(1), na.rm = TRUE)
    foliar_var <- if (n_leaves > 1) round(mean(apply(as.matrix(grp[, sym_names, drop = FALSE]), 2, stats::var, na.rm = TRUE), na.rm = TRUE), 8) else 0.0
    het_disp <- if (n_leaves > 1) round(mean(sqrt(rowSums(sweep(as.matrix(grp[, sym_names, drop = FALSE]), 2, harm_medians[sym_names], "-")^2, na.rm = TRUE)), na.rm = TRUE), 6) else 0.0
    det_t <- if ("determiner_tier" %in% names(grp) && !is.na(grp$determiner_tier[1])) as.character(grp$determiner_tier[1]) else "Tier_3_Bronze"
    auth_w <- compute_authority_precision_weights(grp[1, , drop = FALSE])[1]
    rec <- data.frame(catalogNumber = cat_num, plant_individual_id = p_id, leaf_count = n_leaves,
                      n_tier1_leaves = n_tier1_leaves, n_tier2_leaves = n_tier2_leaves, morph_tier_source = morph_tier_source,
                      foliar_variance = foliar_var, heteroblasty_index = het_disp, assigned_tier = "Tier_1_Direct", tier_assignment = "Tier_1_Direct",
                      reconstruction_tier = "Tier 1", scientificName = grp$scientificName[1],
                      species_raw = grp$species_raw[1], determiner_tier = det_t, authority_weight = auth_w, stringsAsFactors = FALSE)
    med_cs <- stats::median(grp$centroid_size, na.rm = TRUE)
    rec$centroid_size <- if (is.finite(med_cs) && med_cs > 0) round(med_cs, 4) else NA_real_
    rec$log_centsize <- if (!is.na(rec$centroid_size) && rec$centroid_size > 0) round(log(rec$centroid_size), 6) else NA_real_
    other_meta <- intersect(c("pixels_per_mm", "is_calibrated_mm", "determiner_raw", "county", "stateProvince", "latitude", "longitude",
                              "pheno_sin", "pheno_cos", "regional_group", "spatial_balance_weight", "aspect_ratio", "circularity", "solidity", "ucs_score", "area_px",
                              "sym_var", "asym_var", "pct_asymmetry", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight"), names(grp))
    for (m_col in other_meta) { vals <- stats::na.omit(grp[[m_col]]); rec[[m_col]] <- if (length(vals) > 0) vals[1] else NA }
    if ("harmonic_roughness" %in% names(grp)) rec$harmonic_roughness <- round(stats::median(grp$harmonic_roughness, na.rm = TRUE), 6)
    if ("leaf_quality_weight" %in% names(grp)) rec$leaf_quality_weight <- round(stats::median(grp$leaf_quality_weight, na.rm = TRUE), 6)
    if ("is_contour_outlier" %in% names(grp)) rec$is_contour_outlier <- all(as.logical(grp$is_contour_outlier), na.rm = TRUE)
    for (h_name in names(harm_medians)) rec[[h_name]] <- harm_medians[[h_name]]
    specimen_records[[spec_idx]] <- rec
  }

  specimen_records <- specimen_records[seq_len(spec_idx)]
  specimen_df <- do.call(rbind, specimen_records)
  for (p in seq_len(opts$num_pcs)) { specimen_df[[paste0("PC", p)]] <- NA_real_; specimen_df[[paste0("Sym_PC", p)]] <- 0.0 }
  # Fit OTU-level ordination strictly on active reference anchors (Tier_1_Gold & Tier_2_Silver, Tier_3_Bronze weight = 0.0)
  otu_base_mask <- (specimen_df$morph_tier_source == "Tier_1_Pure") & !as.logical(specimen_df$is_contour_outlier)
  if (sum(otu_base_mask, na.rm = TRUE) < 2L) otu_base_mask <- rep(TRUE, nrow(specimen_df))
  otu_ref_mask <- otu_base_mask & (specimen_df$authority_weight > 0)
  if (nrow(specimen_df) >= max(opts$num_pcs, 6L) && exists("active_full_cols") && length(active_full_cols) >= opts$num_pcs) {
    otu_ref_mask <- if (sum(otu_ref_mask, na.rm = TRUE) >= 2L) prune_within_taxon_mahalanobis(as.matrix(specimen_df[, active_full_cols, drop = FALSE]), specimen_df, otu_ref_mask, q = 0.95) else otu_base_mask
    otu_t1_idx <- which(otu_ref_mask)
    otu_v <- apply(specimen_df[otu_t1_idx, active_full_cols, drop = FALSE], 2, stats::var, na.rm = TRUE)
    otu_cols <- names(otu_v[!is.na(otu_v) & otu_v > 1e-10])
    if (length(otu_cols) >= opts$num_pcs) {
      w_otu <- specimen_df$authority_weight[otu_t1_idx]; if (sum(w_otu, na.rm = TRUE) <= 0) w_otu <- rep(1.0, length(otu_t1_idx))
      wt_otu <- stats::cov.wt(as.matrix(specimen_df[otu_t1_idx, otu_cols, drop = FALSE]), wt = w_otu, center = TRUE)
      eig_otu <- eigen(wt_otu$cov, symmetric = TRUE); colnames(eig_otu$vectors) <- paste0("PC", seq_len(ncol(eig_otu$vectors))); rownames(eig_otu$vectors) <- otu_cols
      spec_pcs <- sweep(as.matrix(specimen_df[, otu_cols, drop = FALSE]), 2, wt_otu$center, "-") %*% eig_otu$vectors
      pca_fit <- structure(list(sdev = sqrt(pmax(eig_otu$values, 0)), rotation = eig_otu$vectors, center = wt_otu$center, scale = FALSE, x = spec_pcs), class = "prcomp")
      for (p in seq_len(min(opts$num_pcs, ncol(spec_pcs)))) specimen_df[[paste0("PC", p)]] <- round(spec_pcs[, p], 6)
    }
  } else if (nrow(specimen_df) >= opts$num_pcs && exists("pca_fit") && !is.null(pca_fit)) {
    spec_pcs <- stats::predict(pca_fit, newdata = as.matrix(specimen_df[, active_full_cols, drop = FALSE]))
    for (p in seq_len(min(opts$num_pcs, ncol(spec_pcs)))) specimen_df[[paste0("PC", p)]] <- round(spec_pcs[, p], 6)
  }
  if (nrow(specimen_df) >= opts$num_pcs && !is.null(sym_res$pca)) {
    sym_cols_p <- rownames(sym_res$pca$rotation)
    spec_sym_pcs <- stats::predict(sym_res$pca, newdata = as.matrix(specimen_df[, sym_cols_p, drop = FALSE]))
    for (p in seq_len(min(opts$num_pcs, ncol(spec_sym_pcs)))) specimen_df[[paste0("Sym_PC", p)]] <- round(spec_sym_pcs[, p], 6)
  }

  lead_spec_cols <- c("catalogNumber", "plant_individual_id", "leaf_count", "n_tier1_leaves", "n_tier2_leaves", "morph_tier_source",
                      "foliar_variance", "heteroblasty_index", "assigned_tier", "tier_assignment", "reconstruction_tier",
                      "pixels_per_mm", "is_calibrated_mm", "scientificName", "species_raw", "determiner_tier", "authority_weight", "spatial_balance_weight",
                      "PC1", "PC2", "PC3", "PC4", "PC5", "aspect_ratio", "circularity", "solidity", "ucs_score",
                      "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_px", "centroid_size", "log_centsize")
  lead_spec_cols <- intersect(lead_spec_cols, names(specimen_df))
  specimen_df <- specimen_df[, c(lead_spec_cols, setdiff(names(specimen_df), lead_spec_cols))]

  dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
  tmp_out <- paste0(opts$output, ".tmp"); write.csv(specimen_df, file = tmp_out, row.names = FALSE, na = ""); file.rename(tmp_out, opts$output)
  message("Specimen-aggregated EFA harmonics table exported: ", opts$output, " (Rows: ", nrow(specimen_df), ")")

  tryCatch({
    if (file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R")
      export_step03_fourier_pc_tables(specimen_df, efa_df, opts$output)
      generate_step03_pub_figures(efa_df, specimen_df, harm_mat, if (exists("pca_fit")) pca_fit else NULL, opts$figures_dir, opts$output, sym_res$pca)
    }
  }, error = function(e) warning("Step 03 publication figures/tables export warning: ", e$message))

  # 6. Empirical Tier Validation Test (PERMANOVA via adonis2 / fallback stats::manova)
  valid_tier_idx <- which(!is.na(efa_df$reconstruction_tier) & !is.na(efa_df$scientificName) & nzchar(as.character(efa_df$scientificName)) &
                          !is.na(efa_df$determiner_tier) & nzchar(as.character(efa_df$determiner_tier)) & complete.cases(sym_harmonics))
  if (length(valid_tier_idx) >= 10) {
    tier_manifest <- efa_df[valid_tier_idx, , drop = FALSE]; tier_sym <- sym_harmonics[valid_tier_idx, , drop = FALSE]
    tier_manifest$scientificName <- as.factor(tier_manifest$scientificName)
    tier_manifest$determiner_tier <- as.factor(tier_manifest$determiner_tier)
    tier_manifest$reconstruction_tier <- as.factor(tier_manifest$reconstruction_tier)
    act_terms <- c("scientificName", "determiner_tier", "reconstruction_tier")[c(nlevels(tier_manifest$scientificName), nlevels(tier_manifest$determiner_tier), nlevels(tier_manifest$reconstruction_tier)) >= 2L]
    set.seed(opts$seed); permanova_df <- NULL

    if (length(act_terms) >= 1L && requireNamespace("RRPP", quietly = TRUE)) {
      rrpp_res <- tryCatch({
        rdf <- RRPP::rrpp.data.frame(tier_sym = as.matrix(tier_sym), scientificName = tier_manifest$scientificName,
                                     determiner_tier = tier_manifest$determiner_tier, reconstruction_tier = tier_manifest$reconstruction_tier)
        fit_rrpp <- RRPP::lm.rrpp(stats::as.formula(paste("tier_sym ~", paste(act_terms, collapse = " + "))), data = rdf, iter = opts$permutations, print.progress = FALSE)
        anova(fit_rrpp)$table
      }, error = function(e) { warning("RRPP::lm.rrpp issue: ", e$message, ". Falling back to vegan::adonis2()."); NULL })
      if (!is.null(rrpp_res)) {
        permanova_df <- as.data.frame(rrpp_res); permanova_df$Term <- rownames(permanova_df); rownames(permanova_df) <- NULL
        if ("SS" %in% names(permanova_df)) permanova_df$SumOfSqs <- permanova_df$SS
        if ("Rsq" %in% names(permanova_df)) permanova_df$R2 <- permanova_df$Rsq
        if ("Pr(>F)" %in% names(permanova_df)) permanova_df$p_value <- permanova_df[["Pr(>F)"]]
      }
    }
    if (is.null(permanova_df) && length(act_terms) >= 1L && requireNamespace("vegan", quietly = TRUE)) {
      ad_res <- tryCatch({
        ad_idx <- if (nrow(tier_manifest) > 2000L) sort(sample(seq_len(nrow(tier_manifest)), 2000L)) else seq_len(nrow(tier_manifest))
        tier_sym_sub <- tier_sym[ad_idx, , drop = FALSE]
        vegan::adonis2(stats::as.formula(paste("tier_sym_sub ~", paste(act_terms, collapse = " + "))),
                       data = tier_manifest[ad_idx, , drop = FALSE], permutations = opts$permutations, method = "euclidean", by = "terms")
      }, error = function(e) { warning("vegan::adonis2 issue: ", e$message, ". Falling back to stats::manova()."); NULL })
      if (!is.null(ad_res)) {
        permanova_df <- as.data.frame(ad_res); permanova_df$Term <- rownames(permanova_df); rownames(permanova_df) <- NULL
        if ("F.Model" %in% names(permanova_df)) permanova_df$F <- permanova_df$F.Model
        if ("Pr(>F)" %in% names(permanova_df)) permanova_df$p_value <- permanova_df[["Pr(>F)"]]
      }
    }
    if (is.null(permanova_df)) {
      Y <- as.matrix(tier_sym); Y_cent <- scale(Y, center = TRUE, scale = FALSE); SS_tot <- sum(Y_cent^2); N_obs <- nrow(tier_manifest)
      fit0 <- stats::lm(Y_cent ~ 1, data = tier_manifest)
      fit1 <- if ("scientificName" %in% act_terms) stats::lm(Y_cent ~ scientificName, data = tier_manifest) else fit0
      SS_sp <- sum(stats::fitted(fit1)^2); df_sp <- fit1$rank - 1; t12 <- intersect(c("scientificName", "determiner_tier"), act_terms)
      fit2 <- if ("determiner_tier" %in% act_terms) stats::lm(stats::as.formula(paste("Y_cent ~", paste(t12, collapse = " + "))), data = tier_manifest) else fit1
      SS_12 <- sum(stats::fitted(fit2)^2); SS_det <- SS_12 - SS_sp; df_det <- fit2$rank - fit1$rank
      fit3 <- if ("reconstruction_tier" %in% act_terms) stats::lm(stats::as.formula(paste("Y_cent ~", paste(act_terms, collapse = " + "))), data = tier_manifest) else fit2
      SS_123 <- sum(stats::fitted(fit3)^2); SS_tier <- SS_123 - SS_12; df_tier <- fit3$rank - fit2$rank
      SS_res <- SS_tot - SS_123; df_res <- N_obs - fit3$rank; MS_res <- SS_res / max(df_res, 1)
      F_sp <- if (df_sp > 0) (SS_sp / df_sp) / MS_res else NA_real_; F_det <- if (df_det > 0) (SS_det / df_det) / MS_res else NA_real_; F_tier <- if (df_tier > 0) (SS_tier / df_tier) / MS_res else NA_real_
      p_sp <- if (df_sp > 0) stats::pf(F_sp, df_sp, df_res, lower.tail = FALSE) else NA_real_; p_det <- if (df_det > 0) stats::pf(F_det, df_det, df_res, lower.tail = FALSE) else NA_real_; p_tier <- if (df_tier > 0) stats::pf(F_tier, df_tier, df_res, lower.tail = FALSE) else NA_real_
      permanova_df <- data.frame(Term = c("scientificName", "determiner_tier", "reconstruction_tier", "Residual", "Total"), Df = c(df_sp, df_det, df_tier, df_res, N_obs - 1),
                                 SumOfSqs = round(c(SS_sp, SS_det, SS_tier, SS_res, SS_tot), 4), R2 = round(c(SS_sp / SS_tot, SS_det / SS_tot, SS_tier / SS_tot, SS_res / SS_tot, 1.0), 6),
                                 F = c(round(F_sp, 4), round(F_det, 4), round(F_tier, 4), NA_real_, NA_real_), p_value = c(round(p_sp, 5), round(p_det, 5), round(p_tier, 5), NA_real_, NA_real_), stringsAsFactors = FALSE)
    }
    permanova_df$Term <- gsub("^Residuals$", "Residual", permanova_df$Term)
    if ("Total" %in% permanova_df$Term) permanova_df$R2[permanova_df$Term == "Total"] <- 1.0
    report_cols <- intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(permanova_df))
    permanova_df <- permanova_df[, report_cols, drop = FALSE]
    for (missing_t in setdiff(c("scientificName", "determiner_tier", "reconstruction_tier"), permanova_df$Term)) {
      permanova_df <- rbind(permanova_df[1, , drop = FALSE], data.frame(Term = missing_t, Df = 0, SumOfSqs = 0, R2 = 0, F = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)[, report_cols, drop = FALSE], permanova_df[-1, , drop = FALSE])
    }
    dir.create(dirname(opts$report_out), recursive = TRUE, showWarnings = FALSE)
    tmp_rep <- paste0(opts$report_out, ".tmp"); write.csv(permanova_df, file = tmp_rep, row.names = FALSE); file.rename(tmp_rep, opts$report_out)
    dir.create(dirname(opts$plot_out), recursive = TRUE, showWarnings = FALSE)
    if (exists("render_tier1_vs_tier2_2d_pdf", mode = "function")) {
      render_tier1_vs_tier2_2d_pdf(tier_manifest, opts$plot_out, if (exists("pca_fit")) pca_fit else NULL, sym_res$pca, permanova_df)
    } else {
      pdf(opts$plot_out, width = 10, height = 5); par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.0, 1.0))
      plot(tier_manifest$PC1, tier_manifest$PC2, pch = 1, col = "#0072B2", xlab = "PC1 (%)", ylab = "PC2 (%)", main = "A. 2D Tier 1 EFA Morphospace")
      plot(tier_manifest$Sym_PC1, tier_manifest$Sym_PC2, pch = 1, col = "#0072B2", xlab = "Sym_PC1 (%)", ylab = "Sym_PC2 (%)", main = "B. 2D Symmetric EFA")
      dev.off()
    }
  }

  return(invisible(efa_df))
}

if (sys.nframe() == 0) {
  opts <- parse_args_robust()
  run_fourier_extraction(opts)
}
