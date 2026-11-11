#!/usr/bin/env Rscript
# Script: 03_fourier_extractor.R - Packera Morphometrics Pipeline (NCU)
# Label-blind dynamic-harmonic Elliptic Fourier Analysis (Momocs::efourier) & OTU median aggregation.

suppressPackageStartupMessages({
  if (requireNamespace("Momocs", quietly = TRUE)) library(Momocs); if (requireNamespace("RRPP", quietly = TRUE)) library(RRPP)
  if (requireNamespace("vegan", quietly = TRUE)) library(vegan); if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr); if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("tibble", quietly = TRUE)) library(tibble); if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})
for (.pf_p in c("scripts/morphometrics/pub_figures_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "pub_figures_utils.R"))) {
  if (!is.na(.pf_p) && nzchar(.pf_p) && file.exists(.pf_p)) { source(.pf_p); break }
}
for (.aw_p in c("scripts/morphometrics/authority_weighting_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "authority_weighting_utils.R"))) {
  if (!is.na(.aw_p) && nzchar(.aw_p) && file.exists(.aw_p)) { source(.aw_p); break }
}
for (.dl_p in c("scripts/morphometrics/dissected_leaf_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "dissected_leaf_utils.R"))) {
  if (!is.na(.dl_p) && nzchar(.dl_p) && file.exists(.dl_p)) { source(.dl_p); break }
}
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
}
for (.fp_p in c("scripts/morphometrics/fourier_parallel_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "fourier_parallel_utils.R"))) {
  if (!is.na(.fp_p) && nzchar(.fp_p) && file.exists(.fp_p)) { source(.fp_p); break }
}

# 1. CLI Argument Parsing
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

    optparse::make_option(c("--figures-dir"), type = "character", default = "", help = "Optional output directory for Step 03 publication figures"),
    optparse::make_option(c("--symmetrize-efa-harmonics"), type = "logical", default = TRUE, help = "Uniformly apply Momocs::rm_asym() across Tier 1 contours [default: %default]"),
    optparse::make_option(c("--permutations"), type = "integer", default = 999, help = "Number of permutations for PERMANOVA [default: %default]"),
    optparse::make_option(c("--seed"), type = "integer", default = 42, help = "Random seed for reproducible permutations [default: %default]")
  )
  if (requireNamespace("optparse", quietly = TRUE)) {
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list))
    names(opts) <- gsub("-", "_", names(opts))
    if (!is.null(opts$contours_dir) && opts$contours_dir != "data/contours/") opts$input <- opts$contours_dir else if (!is.null(opts$input) && opts$contours_dir == "data/contours/") opts$contours_dir <- opts$input
    return(opts)
  }
  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(input = "data/contours/", contours_dir = "data/contours/", manifest = "data/tables/extracted_leaf_manifest.csv",
               masks_dir = "data/masks/", vouchers = "data/tables/curated_vouchers.csv", output = "data/tables/leaf_efa_harmonics.csv",
               output_individual = "data/tables/leaf_efa_harmonics_individual.csv", harmonics = 12, num_pcs = 5, symmetrize_efa_harmonics = TRUE,
               figures_dir = "", permutations = 999, seed = 42)
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

# 2. Pure-R Mathematical Engine: Homologized EFA (Momocs / Kuhl & Giardina Fallback)
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

verify_harmonic_power_calibration <- function(out_obj = NULL, harm_mat = NULL, nb_h = 12L) {
  fb_h <- max(6L, min(24L, as.integer(if (is.null(nb_h) || is.na(nb_h)) 12L else nb_h)))
  h_opt <- fb_h; cum_pct <- NA_real_
  if (exists("calibrate_dissection_aware_harmonics", mode = "function")) {
    cal_d <- calibrate_dissection_aware_harmonics(out_obj = out_obj, harm_mat = harm_mat, nb_h = nb_h, min_h = 6L, max_h = 24L)
    h_opt <- cal_d$h_opt; cum_pct <- cal_d$cum_pct
  } else if (!is.null(out_obj) && requireNamespace("Momocs", quietly = TRUE)) {
    cal_res <- tryCatch(Momocs::calibrate_harmonicpower_efourier(if (length(out_obj$coo) > 150L) Momocs::slice(out_obj, seq_len(150L)) else out_obj, nb.h = 24L, plot = FALSE), error = function(e) NULL)
    if (!is.null(cal_res) && !is.null(cal_res$q) && length(cal_res$q) > 0L) h_opt <- fb_h
  }
  if ((is.na(cum_pct) || !is.finite(cum_pct)) && !is.null(harm_mat) && ncol(harm_mat) >= 4L) {
    max_hm <- max(suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(harm_mat)))), fb_h, na.rm = TRUE)
    h_pwr <- vapply(seq_len(max_hm), function(n) {
      cols_n <- intersect(paste0(c("A", "B", "C", "D"), n), colnames(harm_mat))
      if (length(cols_n) > 0L) mean(rowSums(harm_mat[, cols_n, drop = FALSE]^2, na.rm = TRUE) / 2.0, na.rm = TRUE) else 0.0
    }, numeric(1))
    tot_est <- sum(h_pwr) + tail(h_pwr, 1) * 0.15
    cum_pct <- if (tot_est > 1e-12) min(99.99, (sum(h_pwr[seq_len(min(h_opt, length(h_pwr)))]) / tot_est) * 100.0) else 99.50
  }
  if (is.na(cum_pct) || !is.finite(cum_pct)) cum_pct <- 99.50
  message(sprintf("Harmonic Power Calibration (Bonhomme et al. 2014): H_opt = %d captures %.2f%% cumulative harmonic power (threshold >= 99.0%%: %s).",
                  as.integer(h_opt), cum_pct, ifelse(cum_pct >= 99.0, "PASS", "WARN")))
  as.integer(h_opt)
}

# Within-taxon Mahalanobis outlier pruning + Adaptive Anchor Promotion for Gold/Silver-poor taxa
prune_within_taxon_mahalanobis <- function(X, df, ref_mask = NULL, q = 0.95, min_anchors_per_taxon = 2L) {
  elig <- if (!is.null(ref_mask) && length(ref_mask) == nrow(df)) { m <- as.logical(ref_mask); m[is.na(m)] <- FALSE; m } else rep(TRUE, nrow(df))
  t_info_all <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else list(is_type_specimen = rep(FALSE, nrow(df)))
  if (sum(elig) < 2L || is.null(X) || ncol(X) < 1L) {
    attr(elig, "is_promoted_bronze") <- rep(FALSE, nrow(df))
    attr(elig, "is_type_specimen") <- as.logical(t_info_all$is_type_specimen)
    attr(elig, "anchor_tier_source") <- ifelse(t_info_all$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted")
    return(elig)
  }
  anc <- select_adaptive_reference_anchors(df[elig, , drop = FALSE], feature_matrix = as.matrix(X[elig, , drop = FALSE]),
                                           taxon_col = "scientificName", min_anchors_per_taxon = min_anchors_per_taxon, mahal_q_gold_silver = q, mahal_q_bronze_promo = 0.80)
  out <- rep(FALSE, nrow(df)); out[elig] <- anc$active_mask
  attr(out, "is_promoted_bronze") <- { p <- rep(FALSE, nrow(df)); p[elig] <- anc$is_promoted_bronze; p }
  attr(out, "is_type_specimen") <- as.logical(t_info_all$is_type_specimen)
  attr(out, "anchor_tier_source") <- { a <- rep("Passive_Unpromoted", nrow(df)); a[elig] <- as.character(anc$anchor_tier_source); a }
  out
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

  out_df <- if (!is.null(fac_df)) as.data.frame(fac_df) else data.frame(row_id = row_indices)
  solidity <- if ("solidity" %in% names(out_df)) ifelse(is.na(out_df$solidity), 1.0, as.numeric(out_df$solidity)) else rep(1.0, nrow(out_df))
  r_chk <- compute_dissection_aware_roughness(coe, solidity_vec = solidity, pct_asymmetry_vec = pct_asymmetry)
  harmonic_roughness <- r_chk$harmonic_roughness
  is_contour_outlier <- r_chk$is_contour_outlier

  ef_sym <- if (requireNamespace("Momocs", quietly = TRUE) && inherits(ef, "OutCoe")) tryCatch(Momocs::rm_asym(ef), error = function(e) ef) else ef
  coe_sym <- if (is.list(ef_sym) && !is.null(ef_sym$coe)) as.matrix(ef_sym$coe) else as.matrix(ef_sym)
  sym_cols <- grep("^([AD]_?[0-9]+|C_?1)$", colnames(coe_sym), value = TRUE)
  sym_mat <- coe_sym[, if (length(sym_cols) > 0) sym_cols else colnames(coe_sym), drop = FALSE]
  sym_mat[!is.finite(sym_mat)] <- 0.0
  clean_pca_mask <- if (nrow(sym_mat) >= 10 && sum(!is_contour_outlier) >= 5) !is_contour_outlier else rep(TRUE, nrow(sym_mat))
  lz_sym <- if ("leaf_zone" %in% names(out_df)) tolower(trimws(as.character(out_df$leaf_zone))) else rep("basal", nrow(out_df))
  bas_pca_mask <- clean_pca_mask & (is.na(lz_sym) | lz_sym == "basal")
  fit_pca_mask <- if (sum(bas_pca_mask, na.rm = TRUE) >= 2) bas_pca_mask else clean_pca_mask

  arch_w <- if ("archetype_score" %in% names(out_df)) {
    av <- suppressWarnings(as.numeric(out_df$archetype_score))
    ifelse(is.na(av) | !is.finite(av) | av <= 0, 1.0, pmin(pmax(av, 0.20), 1.0))
  } else rep(1.0, nrow(out_df))
  leaf_quality_weight <- pmax(0.05, arch_w * ((1.0 - pmin(pct_asymmetry, 0.95))^2) * ifelse(is_contour_outlier, 0.10, 1.0))
  spatial_wt <- if ("spatial_balance_weight" %in% names(out_df) && any(is.finite(out_df$spatial_balance_weight))) as.numeric(out_df$spatial_balance_weight) else if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(out_df) else rep(1.0, nrow(out_df))
  active_ref_mask <- if (sum(fit_pca_mask, na.rm = TRUE) >= 2L) prune_within_taxon_mahalanobis(sym_mat, out_df, fit_pca_mask, q = 0.95) else fit_pca_mask
  out_df$is_promoted_bronze <- if (!is.null(attr(active_ref_mask, "is_promoted_bronze"))) attr(active_ref_mask, "is_promoted_bronze") else rep(FALSE, nrow(out_df))
  out_df$is_type_specimen <- if (!is.null(attr(active_ref_mask, "is_type_specimen"))) attr(active_ref_mask, "is_type_specimen") else if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(out_df)$is_type_specimen else rep(FALSE, nrow(out_df))
  out_df$anchor_tier_source <- if (!is.null(attr(active_ref_mask, "anchor_tier_source"))) attr(active_ref_mask, "anchor_tier_source") else ifelse(out_df$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted")
  tier_wt <- compute_authority_precision_weights(out_df, active_mask = active_ref_mask)
  if (sum(active_ref_mask, na.rm = TRUE) < 2L) active_ref_mask <- fit_pca_mask
  col_vars <- apply(sym_mat[active_ref_mask, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  sym_mat <- sym_mat[, names(col_vars[!is.na(col_vars) & col_vars > 1e-10]), drop = FALSE]

  out_df$authority_weight <- tier_wt; out_df$spatial_balance_weight <- round(spatial_wt, 6)
  out_df$sym_var <- round(sym_var, 6); out_df$asym_var <- round(asym_var, 6); out_df$pct_asymmetry <- round(pct_asymmetry, 6)
  out_df$harmonic_roughness <- round(harmonic_roughness, 6); out_df$is_contour_outlier <- is_contour_outlier; out_df$leaf_quality_weight <- round(leaf_quality_weight, 6)
  for (pc_col in paste0("Sym_PC", seq_len(5L))) out_df[[pc_col]] <- 0.0
  pca_sym <- NULL; k_95_sym <- 5L
  if (nrow(sym_mat) >= 2 && ncol(sym_mat) >= 1) {
    X_cl <- sym_mat[active_ref_mask, , drop = FALSE]; w_cl <- tier_wt[active_ref_mask]; if (sum(w_cl, na.rm = TRUE) <= 0) w_cl <- rep(1.0, nrow(X_cl))
    wt_cov <- stats::cov.wt(X_cl, wt = w_cl, center = TRUE); eig <- eigen(wt_cov$cov, symmetric = TRUE)
    colnames(eig$vectors) <- paste0("PC", seq_len(ncol(eig$vectors))); rownames(eig$vectors) <- colnames(sym_mat)
    sym_scores <- sweep(sym_mat, 2, wt_cov$center, "-") %*% eig$vectors
    pca_sym <- structure(list(sdev = sqrt(pmax(eig$values, 0)), rotation = eig$vectors, center = wt_cov$center, scale = FALSE, x = sym_scores), class = "prcomp")
    sv2 <- pca_sym$sdev^2; k_95_sym <- which((cumsum(sv2) / max(sum(sv2), 1e-12)) >= 0.95)[1]; if (is.na(k_95_sym)) k_95_sym <- ncol(sym_scores)
    k_sym_exp <- max(5L, min(max(as.integer(num_pcs), as.integer(k_95_sym)), 15L, ncol(sym_scores)))
    for (p in seq_len(k_sym_exp)) out_df[[paste0("Sym_PC", p)]] <- if (p <= ncol(sym_scores)) round(sym_scores[, p], 6) else 0.0
  }
  list(df = out_df, sym_mat = sym_mat, pca = pca_sym, k_95 = as.integer(k_95_sym), permanova = NULL, active_ref_mask = active_ref_mask)
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
  sec_dup_cats <- character(0); active_voucher_cats <- character(0)
  if (!is.null(vouchers_df)) {
    if ("is_primary_duplicate" %in% names(vouchers_df)) {
      is_sec <- toupper(trimws(as.character(vouchers_df$is_primary_duplicate))) %in% c("FALSE", "F", "0")
      sec_dup_cats <- unique(as.character(vouchers_df$catalogNumber[is_sec]))
      vouchers_df <- vouchers_df[!is_sec, , drop = FALSE]
    }
    if ("catalogNumber" %in% names(vouchers_df)) active_voucher_cats <- unique(trimws(as.character(vouchers_df$catalogNumber)))
    message("Loaded ", nrow(vouchers_df), " curated primary voucher records.")
  }

  # 1. Discover coordinates via SQLite extracted_leaves or Dual Filesystem Ingestion Engine
  harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(opts$harmonics), function(x, y) paste0(x, y)))
  db_ing <- if (exists("fetch_step03_contours_from_db", mode = "function")) fetch_step03_contours_from_db(opts$vouchers, opts$manifest, opts$contours_dir, active_catalogs = if (length(active_voucher_cats) > 0L) setdiff(active_voucher_cats, sec_dup_cats) else NULL, align_fn = function(coo, ...) sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-"), centsize_fn = compute_centsize_robust) else NULL
  if (!is.null(db_ing)) {
    coo_list <- db_ing$coords_list; fac_df <- db_ing$meta_df; shape_names_vec <- names(coo_list)
    if (exists("connect_packera_db", mode = "function") && requireNamespace("DBI", quietly = TRUE) && !all(c("leaf_zone", "leaf_zone_confidence", "dist_from_crown_norm", "area_mm2") %in% names(fac_df))) {
      tryCatch({ con_lz <- connect_packera_db(read_only = TRUE); on.exit(tryCatch(DBI::dbDisconnect(con_lz), error = function(e) NULL), add = TRUE)
        if (!is.null(con_lz) && DBI::dbExistsTable(con_lz, "extracted_leaves")) {
          f_el <- DBI::dbListFields(con_lz, "extracted_leaves"); sel_el <- intersect(c("catalogNumber", "plant_individual_id", "leaf_id", "contour_xy_json", "leaf_zone", "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox", "area_mm2", "area_px"), f_el)
          if (all(c("catalogNumber", "leaf_id") %in% sel_el)) {
            el_lz <- DBI::dbGetQuery(con_lz, sprintf("SELECT %s FROM extracted_leaves", paste(setdiff(sel_el, "contour_xy_json"), collapse = ", ")))
            k_el <- paste(el_lz$catalogNumber, ifelse(is.na(el_lz$plant_individual_id), 0L, as.integer(el_lz$plant_individual_id)), as.integer(el_lz$leaf_id), sep = "__")
            m_el <- match(paste(fac_df$catalogNumber, ifelse(is.na(fac_df$plant_individual_id), 0L, as.integer(fac_df$plant_individual_id)), as.integer(fac_df$leaf_id), sep = "__"), k_el)
            for (cn_el in intersect(c("leaf_zone", "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox", "area_mm2", "area_px"), names(el_lz))) if (!(cn_el %in% names(fac_df))) fac_df[[cn_el]] <- el_lz[[cn_el]][m_el]
          }
        }
      }, error = function(e) NULL)
    }
    message("Ingested ", length(coo_list), " Tier 1 intact 120-vertex contours directly from SQLite extracted_leaves.")
  } else {
  contour_files <- character(0); use_masks_fallback <- FALSE
  manifest_df <- if (file.exists(opts$manifest)) read.csv(opts$manifest, stringsAsFactors = FALSE) else NULL
  if (!is.null(manifest_df) && "catalogNumber" %in% names(manifest_df)) {
    if (length(sec_dup_cats) > 0) manifest_df <- manifest_df[!as.character(manifest_df$catalogNumber) %in% sec_dup_cats, , drop = FALSE]
    if (length(active_voucher_cats) > 0) manifest_df <- manifest_df[trimws(as.character(manifest_df$catalogNumber)) %in% active_voucher_cats, , drop = FALSE]
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
  shape_stems <- sub("\\.[^.]+$", "", basename(contour_files))
  contour_files <- contour_files[!duplicated(shape_stems, fromLast = TRUE)]
  if (length(active_voucher_cats) > 0) contour_files <- contour_files[vapply(contour_files, function(f) trimws(as.character(parse_specimen_id(f)$catalogNumber)), character(1)) %in% active_voucher_cats]
  manifest_by_basename <- NULL; manifest_by_cat_pid_leaf <- NULL; manifest_by_cat_leaf <- NULL
  if (!is.null(manifest_df) && nrow(manifest_df) > 0) {
    if ("contour_path" %in% names(manifest_df)) {
      bnames <- basename(as.character(manifest_df$contour_path))
      manifest_by_basename <- seq_len(nrow(manifest_df)); names(manifest_by_basename) <- bnames
    }
    if (all(c("catalogNumber", "leaf_id") %in% names(manifest_df))) {
      pids_m <- if ("plant_individual_id" %in% names(manifest_df)) ifelse(is.na(manifest_df$plant_individual_id), 0L, as.integer(manifest_df$plant_individual_id)) else rep(0L, nrow(manifest_df))
      manifest_by_cat_pid_leaf <- stats::setNames(seq_len(nrow(manifest_df)), paste0(as.character(manifest_df$catalogNumber), "__", pids_m, "__", as.integer(manifest_df$leaf_id)))
      cat_leaf_keys <- paste0(as.character(manifest_df$catalogNumber), "__", as.integer(manifest_df$leaf_id))
      manifest_by_cat_leaf <- seq_len(nrow(manifest_df)); names(manifest_by_cat_leaf) <- cat_leaf_keys
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
      if (is.na(idx_hit) && !is.null(manifest_by_cat_pid_leaf)) idx_hit <- unname(manifest_by_cat_pid_leaf[paste0(as.character(meta$catalogNumber), "__", as.integer(meta$plant_id), "__", as.integer(meta$leaf_id))])
      if (is.na(idx_hit) && !is.null(manifest_by_cat_leaf)) idx_hit <- unname(manifest_by_cat_leaf[paste0(as.character(meta$catalogNumber), "__", as.integer(meta$leaf_id))])
      if (!is.na(idx_hit) && idx_hit >= 1L && idx_hit <= nrow(manifest_df)) m_match <- manifest_df[idx_hit, , drop = FALSE]
    }
    resolved_pid <- if (!is.null(m_match) && nrow(m_match) > 0 && "plant_individual_id" %in% names(m_match) && !is.na(m_match$plant_individual_id[1])) as.integer(m_match$plant_individual_id[1]) else meta$plant_id
    if (!is.null(m_match) && nrow(m_match) > 0L && exists("evaluate_step03_lm3_efa_qc", mode = "function") && !isTRUE(evaluate_step03_lm3_efa_qc(m_match)[1L])) return(NULL)
    raw_scale_src <- if (!is.null(m_match) && nrow(m_match) > 0 && "scale_source" %in% names(m_match) && !is.na(m_match$scale_source[1])) tolower(trimws(as.character(m_match$scale_source[1]))) else ""
    resolved_ppm <- if (!is.null(m_match) && nrow(m_match) > 0 && "pixels_per_mm" %in% names(m_match)) suppressWarnings(as.numeric(m_match$pixels_per_mm[1])) else NA_real_
    is_calibrated_mm <- !is.na(resolved_ppm) && is.finite(resolved_ppm) && resolved_ppm > 0 && !(raw_scale_src %in% c("mp_estimated", "predicted_from_megapixels", "uncalibrated"))
    scale_source <- if (is_calibrated_mm) "ruler_calibrated" else if (!is.na(resolved_ppm) && is.finite(resolved_ppm) && resolved_ppm > 0) "mp_estimated" else "uncalibrated" 

    # Calculate raw pixel area, circularity, Dissection Index & polygon solidity before coordinate scaling
    area_px <- abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1])))
    perim_px <- sum(sqrt(diff(c(coo[, 1], coo[1, 1]))^2 + diff(c(coo[, 2], coo[1, 2]))^2))
    circularity <- if (is.finite(perim_px) && perim_px > 1e-6) round(min(1.0, (4 * pi * area_px) / (perim_px^2)), 4) else NA_real_
    dissection_index <- compute_dissection_index_from_coo(area_val = area_px, perim_val = perim_px)
    resolved_solidity <- if (!is.null(m_match) && nrow(m_match) > 0 && "solidity" %in% names(m_match) && !is.na(m_match$solidity[1])) round(suppressWarnings(as.numeric(m_match$solidity[1])), 4) else compute_polygon_solidity(coo, area_px)
    m_num <- function(cn) if (!is.null(m_match) && nrow(m_match) > 0 && cn %in% names(m_match) && !is.na(m_match[[cn]][1])) suppressWarnings(as.numeric(m_match[[cn]][1])) else NA_real_
    m_dis <- !is.na(resolved_solidity) && is.finite(resolved_solidity) && (resolved_solidity < 0.72)

    # Calibrate contour coordinates to physical millimeters (mm) strictly when ruler_calibrated
    centsize_px <- compute_centsize_robust(coo); if (is_calibrated_mm) coo <- coo / resolved_ppm
    centsize <- if (is_calibrated_mm) compute_centsize_robust(coo) else NA_real_
    log_cs <- if (is_calibrated_mm && !is.na(centsize) && centsize > 0) round(log10(centsize), 6) else NA_real_

    # Center contour coordinates while preserving LeafMachine3 tip-up, lamina_base-anchored clockwise outline
    coo <- if (has_momocs) tryCatch(Momocs::coo_center(coo), error = function(e) sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-")) else sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-")
    lw <- if (has_momocs) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
    if (is.null(lw) || length(lw) < 2 || is.na(lw[2]) || lw[2] <= 1e-6) {
      coo_centered <- scale(coo, center = TRUE, scale = FALSE); svd_rot <- coo_centered %*% svd(coo_centered)$v
      lw <- c(diff(range(svd_rot[, 1], na.rm = TRUE)), diff(range(svd_rot[, 2], na.rm = TRUE)))
    }
    aspect_ratio <- round(as.numeric(lw[1] / lw[2]), 4)
    shape_name <- meta$shape_id
    m_zone <- if (!is.null(m_match) && nrow(m_match) > 0 && "leaf_zone" %in% names(m_match) && !is.na(m_match$leaf_zone[1])) as.character(m_match$leaf_zone[1]) else NA_character_
    m_ect <- if (!is.null(m_match) && nrow(m_match) > 0 && "ect_path" %in% names(m_match) && !is.na(m_match$ect_path[1])) as.character(m_match$ect_path[1]) else NA_character_
    m_amm2 <- if (!is.na(m_num("area_mm2"))) m_num("area_mm2") else if (is_calibrated_mm) round(area_px / (resolved_ppm^2), 2) else NA_real_
    rec <- data.frame(
      shape_id = shape_name, catalogNumber = meta$catalogNumber, plant_individual_id = resolved_pid,
      leaf_id = meta$leaf_id, ect_path = m_ect,
      leaf_zone = m_zone, leaf_zone_confidence = m_num("leaf_zone_confidence"), dist_from_crown_norm = m_num("dist_from_crown_norm"),
      area_mm2 = m_amm2, pixels_per_mm = resolved_ppm, is_calibrated_mm = is_calibrated_mm, scale_source = scale_source, aspect_ratio = aspect_ratio,
      circularity = circularity, solidity = resolved_solidity, dissection_index = dissection_index, archetype_score = m_num("archetype_score"),
      petiole_length_mm = m_num("petiole_length_mm"), blade_length_mm = m_num("blade_length_mm"),
      petiole_to_blade_ratio = m_num("petiole_to_blade_ratio"), is_dissected = m_dis,
      cy_rel_plant_bbox = m_num("cy_rel_plant_bbox"), cy_rel_sheet = m_num("cy_rel_sheet"),
      blade_area_cm2 = m_num("blade_area_cm2"),
      area_px = round(area_px, 1), centroid_size_px = if (!is.na(centsize_px)) round(centsize_px, 2) else NA_real_, centroid_size = if (!is.na(centsize)) round(centsize, 4) else NA_real_,
      log_centsize = log_cs, mask_source = f, stringsAsFactors = FALSE
    )
    list(shape_name = shape_name, coo = coo, rec = rec)
  }

  if (.Platform$OS.type == "unix" && length(contour_files) >= 100 && exists("packera_mclapply", mode = "function")) {
    n_cores <- resolve_packera_workers(32L, length(contour_files))
    ingested_list <- packera_mclapply(contour_files, ingest_single_contour, mc.cores = n_cores)
  } else {
    ingested_list <- lapply(contour_files, ingest_single_contour)
  }

  valid_ingested <- Filter(Negate(is.null), ingested_list)
  if (length(valid_ingested) == 0) stop("Failed to extract valid Tier 1 intact coordinate contours.")

  fac_records <- lapply(valid_ingested, `[[`, "rec")
  fac_df <- do.call(rbind, fac_records)
  coo_list <- lapply(valid_ingested, `[[`, "coo")
  shape_names_vec <- vapply(valid_ingested, `[[`, character(1), "shape_name")
  }

  # Ensure strict 1-to-1 uniqueness across active vouchers
  keep_unique <- !duplicated(shape_names_vec, fromLast = TRUE) & !(as.character(fac_df$catalogNumber) %in% sec_dup_cats) & (length(active_voucher_cats) == 0L | trimws(as.character(fac_df$catalogNumber)) %in% active_voucher_cats)
  fac_df <- fac_df[keep_unique, , drop = FALSE]; coo_list <- coo_list[keep_unique]
  names(coo_list) <- shape_names_vec[keep_unique]; rownames(fac_df) <- names(coo_list)
  had_explicit_lz <- ("leaf_zone" %in% names(fac_df) && any(!is.na(fac_df$leaf_zone) & nzchar(trimws(as.character(fac_df$leaf_zone))))) || ("dist_from_crown_norm" %in% names(fac_df) && any(is.finite(suppressWarnings(as.numeric(fac_df$dist_from_crown_norm)))))
  if (exists("filter_basal_lamina_contours", mode = "function")) {
    bg_res <- filter_basal_lamina_contours(fac_df, coo_list, retain_cauline_tagged = TRUE); fac_df <- bg_res$fac_df; coo_list <- bg_res$coo_list
  }
  if (exists("align_contour_horizontal_xax", mode = "function")) coo_list <- lapply(coo_list, align_contour_horizontal_xax)
  message(sprintf("Assembled %d valid Tier 1 intact closed contours for Elliptic Fourier Analysis.", length(coo_list)))

  # 3. Dynamic Harmonic Power Calibration (H_opt in [6, 24], >=99% power) & Homologized EFA
  efa_success <- FALSE; harm_mat <- NULL; out_obj <- NULL; ef_res <- NULL; h_opt <- max(6L, min(24L, as.integer(opts$harmonics)))
  if (requireNamespace("Momocs", quietly = TRUE)) {
    if (length(coo_list) >= 240L && exists("run_parallel_momocs_efa", mode = "function") &&
        exists("resolve_packera_workers", mode = "function") && resolve_packera_workers(32L) > 1L) {
      p_efa <- tryCatch(run_parallel_momocs_efa(coo_list, fac_df, opts$harmonics, verify_harmonic_power_calibration), error = function(e) NULL)
      if (!is.null(p_efa) && !is.null(p_efa$harm_mat)) {
        out_obj <- p_efa$out_obj; fac_df <- p_efa$fac_df; h_opt <- p_efa$h_opt; ef_res <- p_efa$ef_res; harm_mat <- p_efa$harm_mat; efa_success <- TRUE
        harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(h_opt), function(x, y) paste0(x, y)))
        message("Multi-threaded Momocs EFA extraction successful across ", nrow(harm_mat), " outlines (H_opt = ", h_opt, ").")
      }
    }
    if (!efa_success) {
      out_obj <- tryCatch(Momocs::Out(coo_list, fac = fac_df), error = function(e) { message("Momocs::Out() note: ", e$message); NULL })
      if (!is.null(out_obj)) {
        cs_vals <- tryCatch(Momocs::coo_centsize(out_obj), error = function(e) NULL)
        if (!is.null(cs_vals) && length(cs_vals) == nrow(fac_df)) {
          cal_ok <- if ("is_calibrated_mm" %in% names(fac_df)) (!is.na(fac_df$is_calibrated_mm) & as.logical(fac_df$is_calibrated_mm)) else rep(TRUE, nrow(fac_df))
          fac_df$centroid_size <- ifelse(cal_ok, round(as.numeric(cs_vals), 4), NA_real_)
          fac_df$log_centsize <- ifelse(cal_ok & !is.na(fac_df$centroid_size) & fac_df$centroid_size > 0, round(log10(fac_df$centroid_size), 6), NA_real_)
        }
        out_obj <- Momocs::coo_center(out_obj) %>% Momocs::coo_scale()
        h_opt <- verify_harmonic_power_calibration(out_obj, nb_h = opts$harmonics)
        harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(h_opt), function(x, y) paste0(x, y)))
        message("Executing ", h_opt, "-harmonic EFA via Momocs::efourier(nb.h = ", h_opt, ", norm = FALSE, start = FALSE)...")
        # Baseline homologous signature: Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = FALSE, start = FALSE)
        ef_res <- tryCatch(Momocs::efourier(out_obj, nb.h = h_opt, norm = FALSE, start = FALSE), error = function(e) { message("Momocs::efourier() note: ", e$message); NULL })
        if (!is.null(ef_res) && !is.null(ef_res$coe)) {
          harm_mat <- as.matrix(ef_res$coe); efa_success <- TRUE
          message("Momocs EFA extraction successful across ", nrow(harm_mat), " outlines (H_opt = ", h_opt, ").")
        }
      }
    }
  }
  if (!efa_success || is.null(harm_mat)) {
    h_opt <- verify_harmonic_power_calibration(out_obj, nb_h = opts$harmonics)
    harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(h_opt), function(x, y) paste0(x, y)))
    message("Running mathematical homologized EFA engine (nb_h = ", h_opt, ", norm = FALSE, start = FALSE)...")
    harm_list <- if (length(coo_list) >= 240L && exists("packera_mclapply", mode = "function")) packera_mclapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = h_opt, norm = FALSE, start = FALSE), mc.cores = resolve_packera_workers(32L, length(coo_list))) else lapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = h_opt, norm = FALSE, start = FALSE))
    harm_mat <- do.call(rbind, harm_list); colnames(harm_mat) <- harmonic_names
  }
  if (all(harmonic_names %in% colnames(harm_mat))) harm_mat <- harm_mat[, harmonic_names, drop = FALSE] else if (ncol(harm_mat) == length(harmonic_names)) colnames(harm_mat) <- harmonic_names
  if (!is.null(vouchers_df)) {
    meta_cols <- intersect(colnames(vouchers_df), c("catalogNumber", "taxon", "species_standardized", "scientificName", "species_raw", "institutionCode",
      "determiner_raw", "determiner_tier", "authority_tier", "type_status", "dawid_skene_error_rate", "label_error_rate",
      "stateProvince", "county", "locality", "decimalLatitude", "decimalLongitude", "latitude", "longitude",
      "coordinateUncertainty", "year", "month", "day", "eventDate", "pheno_sin", "pheno_cos", "regional_group", "spatial_balance_weight"))
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), meta_cols, drop = FALSE]
    fac_df$.row_ord <- seq_len(nrow(fac_df))
    fac_df <- merge(fac_df, v_sub, by = "catalogNumber", all.x = TRUE, sort = FALSE)
    fac_df <- fac_df[order(fac_df$.row_ord), , drop = FALSE]; fac_df$.row_ord <- NULL; rownames(fac_df) <- fac_df$shape_id
  }
  if (!"species_raw" %in% names(fac_df)) fac_df$species_raw <- if ("taxon" %in% names(fac_df)) fac_df$taxon else NA_character_
  if (!"scientificName" %in% names(fac_df)) fac_df$scientificName <- fac_df$species_raw
  else fac_df$scientificName <- ifelse(is.na(fac_df$scientificName) | fac_df$scientificName == "", fac_df$species_raw, fac_df$scientificName)
  if (exists("standardize_taxon_name", mode = "function")) {
    sp_std <- rep(NA_character_, nrow(fac_df))
    for (sc in intersect(c("species_standardized", "scientificName", "species_raw", "taxon"), names(fac_df))) {
      sv <- standardize_taxon_name(fac_df[[sc]]); hit_s <- is.na(sp_std) & !is.na(sv) & nzchar(sv); sp_std[hit_s] <- sv[hit_s]
    }
    fac_df$species_standardized <- sp_std
    fac_df$scientificName <- ifelse(!is.na(sp_std) & nzchar(sp_std), sp_std, fac_df$scientificName)
  }
  if (!"determiner_tier" %in% names(fac_df)) fac_df$determiner_tier <- if ("authority_tier" %in% names(fac_df)) as.character(fac_df$authority_tier) else "Tier_3_Bronze"
  else fac_df$determiner_tier <- ifelse(is.na(fac_df$determiner_tier) | fac_df$determiner_tier == "", if ("authority_tier" %in% names(fac_df)) ifelse(is.na(fac_df$authority_tier) | fac_df$authority_tier == "", "Tier_3_Bronze", as.character(fac_df$authority_tier)) else "Tier_3_Bronze", fac_df$determiner_tier)
  fac_df$is_type_specimen <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(fac_df)$is_type_specimen else rep(FALSE, nrow(fac_df))
  fac_df$authority_weight <- compute_authority_precision_weights(fac_df)
  fac_df$spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(fac_df) else rep(1.0, nrow(fac_df))

  # 4. Symmetric Fourier Decomposition, Optional Uniform Tier-1 Symmetrization (Momocs::rm_asym) & Reference/Passive PCA
  message("Extracting symmetric harmonic coefficients (A1..A", h_opt, " and D1..D", h_opt, " harmonics)...")
  sym_names <- c(paste0("A", seq_len(h_opt)), paste0("D", seq_len(h_opt)))
  asym_names <- c(paste0("B", seq_len(h_opt)), paste0("C", seq_len(h_opt)))
  ef_input <- if (exists("ef_res") && !is.null(ef_res)) { ef_res$fac <- fac_df; ef_res } else list(coe = harm_mat, fac = fac_df)
  sym_res <- extract_symmetric_fourier_pcs(ef_input, fac_df = fac_df, num_pcs = max(5L, min(15L, as.integer(opts$num_pcs))), permutations = opts$permutations, seed = opts$seed, run_permanova = FALSE)

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
  for (scol in c("is_type_specimen", "anchor_tier_source", "authority_weight", "spatial_balance_weight", "sym_var", "asym_var", "pct_asymmetry", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight")) {
    if (scol %in% names(sym_res$df)) efa_df[[scol]] <- sym_res$df[[scol]]
  }
  tier_wt <- efa_df$authority_weight

  message("Running strict reference/passive authority-split EFA PCA (Tier_3_Bronze weight = 0.0)...")
  complete_idx <- which(complete.cases(harm_mat))
  k_95 <- if (!is.null(sym_res$k_95) && is.finite(sym_res$k_95)) as.integer(sym_res$k_95) else 5L
  k_export <- max(5L, min(k_95, 15L))
  for (p in seq_len(k_export)) { efa_df[[paste0("PC", p)]] <- NA_real_; efa_df[[paste0("Sym_PC", p)]] <- if (paste0("Sym_PC", p) %in% names(sym_res$df)) sym_res$df[[paste0("Sym_PC", p)]] else 0.0 }

  if (length(complete_idx) >= 2L) {
    clean_idx <- if (length(complete_idx) >= 10 && sum(!efa_df$is_contour_outlier[complete_idx]) >= 5) complete_idx[!efa_df$is_contour_outlier[complete_idx]] else complete_idx
    lz_f <- if ("leaf_zone" %in% names(efa_df)) tolower(trimws(as.character(efa_df$leaf_zone[clean_idx]))) else rep("basal", length(clean_idx))
    bas_idx <- clean_idx[is.na(lz_f) | lz_f == "basal"]
    base_fit_idx <- if (length(bas_idx) >= 2L) bas_idx else clean_idx
    base_mask <- seq_len(nrow(efa_df)) %in% base_fit_idx
    if (isTRUE(opts$symmetrize_efa_harmonics) && !is.null(sym_res$active_ref_mask) &&
        length(sym_res$active_ref_mask) == nrow(efa_df) && length(complete_idx) == nrow(efa_df)) {
      ref_mask <- as.logical(sym_res$active_ref_mask) & base_mask
      efa_df$is_promoted_bronze <- if ("is_promoted_bronze" %in% names(sym_res$df)) as.logical(sym_res$df$is_promoted_bronze) else rep(FALSE, nrow(efa_df))
      if (!identical(ref_mask, as.logical(sym_res$active_ref_mask))) {
        efa_df$authority_weight <- compute_authority_precision_weights(efa_df, active_mask = ref_mask); tier_wt <- efa_df$authority_weight
      }
    } else {
      ref_mask <- if (sum(base_mask, na.rm = TRUE) >= 2L) prune_within_taxon_mahalanobis(harm_mat, efa_df, base_mask, q = 0.95) else base_mask
      efa_df$is_promoted_bronze <- if (!is.null(attr(ref_mask, "is_promoted_bronze"))) attr(ref_mask, "is_promoted_bronze") else rep(FALSE, nrow(efa_df))
      efa_df$anchor_tier_source <- if (!is.null(attr(ref_mask, "anchor_tier_source"))) attr(ref_mask, "anchor_tier_source") else ifelse(efa_df$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted")
      efa_df$authority_weight <- compute_authority_precision_weights(efa_df, active_mask = ref_mask); tier_wt <- efa_df$authority_weight
    }
    if (sum(ref_mask, na.rm = TRUE) < 2L) ref_mask <- base_mask
    fit_idx <- which(ref_mask)
    col_vars <- apply(harm_mat[fit_idx, , drop = FALSE], 2, stats::var, na.rm = TRUE)
    active_full_cols <- names(col_vars[!is.na(col_vars) & col_vars > 1e-10])
    if (length(active_full_cols) == 0L) active_full_cols <- sym_names[seq_len(min(length(sym_names), ncol(harm_mat)))]
    message(sprintf("Active harmonic dimensions with non-zero variance: %d / %d", length(active_full_cols), ncol(harm_mat)))
    w_fit <- tier_wt[fit_idx]; if (sum(w_fit, na.rm = TRUE) <= 0) w_fit <- rep(1.0, length(fit_idx))
    wt_cov <- stats::cov.wt(harm_mat[fit_idx, active_full_cols, drop = FALSE], wt = w_fit, center = TRUE)
    eig <- eigen(wt_cov$cov, symmetric = TRUE)
    colnames(eig$vectors) <- paste0("PC", seq_len(ncol(eig$vectors))); rownames(eig$vectors) <- active_full_cols
    pc_scores <- sweep(harm_mat[complete_idx, active_full_cols, drop = FALSE], 2, wt_cov$center, "-") %*% eig$vectors
    pca_fit <- structure(list(sdev = sqrt(pmax(eig$values, 0)), rotation = eig$vectors, center = wt_cov$center, scale = FALSE, x = pc_scores), class = "prcomp")
    tot_ev <- max(sum(pca_fit$sdev^2), 1e-12); var_exp <- round((pca_fit$sdev^2) / tot_ev * 100, 2); cum_var <- cumsum(pca_fit$sdev^2) / tot_ev
    k_95_full <- which(cum_var >= 0.95)[1]; if (is.na(k_95_full) || k_95_full < 1L) k_95_full <- ncol(pc_scores)
    k_95_sym <- if (!is.null(sym_res$k_95) && is.finite(sym_res$k_95)) as.integer(sym_res$k_95) else k_95_full
    k_95 <- as.integer(max(k_95_full, k_95_sym, na.rm = TRUE))
    k_export <- max(5L, min(as.integer(k_95), 15L, ncol(pc_scores)))
    message(sprintf("=== Dynamic %d-Harmonic Morphospace PCA (k_95 = %d, k_export = %d) ===", h_opt, k_95, k_export))
    for (p in seq_len(min(k_export, length(var_exp)))) message(sprintf("  PC%d: %5.2f%% variance", p, var_exp[p]))
    sym_sc_m <- if (!is.null(sym_res$pca) && !is.null(sym_res$pca$x)) sym_res$pca$x else NULL
    for (p in seq_len(k_export)) {
      efa_df[[paste0("PC", p)]] <- NA_real_
      efa_df[complete_idx, paste0("PC", p)] <- if (p <= ncol(pc_scores)) round(pc_scores[, p], 6) else 0.0
      efa_df[[paste0("Sym_PC", p)]] <- if (!is.null(sym_sc_m) && p <= ncol(sym_sc_m)) round(sym_sc_m[, p], 6) else if (paste0("Sym_PC", p) %in% names(sym_res$df)) sym_res$df[[paste0("Sym_PC", p)]] else 0.0
    }
    k_mcd <- max(2L, min(k_95_full, ncol(pc_scores), max(2L, floor(nrow(pc_scores) / 3L))))
    mcd_out <- screen_regime_stratified_mcd_outliers(pc_scores, fit_pos = match(fit_idx, complete_idx), solidity_vec = efa_df$solidity[complete_idx], is_dissected_vec = efa_df$is_dissected[complete_idx], k_mcd = k_mcd, p_cutoff = 0.001, seed = opts$seed)
    if ("is_type_specimen" %in% names(efa_df)) {
      mcd_out[as.logical(efa_df$is_type_specimen[complete_idx])] <- FALSE
    }
    if (any(mcd_out) && sum(!mcd_out) >= 5L) {
      efa_df$is_contour_outlier[complete_idx] <- as.logical(efa_df$is_contour_outlier[complete_idx]) | mcd_out
      message(sprintf("Regime-stratified MCD robust Mahalanobis outlier filter (top %d PCs, k_95=%d, p < 0.001): flagged %d contours.", k_mcd, k_95, sum(mcd_out)))
    }
  }

  # 5. Darwin Core Metadata Standardization & Export
  harm_mat <- as.matrix(efa_df[, harmonic_names, drop = FALSE])
  sym_harmonics <- as.matrix(efa_df[, sym_names, drop = FALSE])
  efa_df$efd_harmonic_opt <- as.integer(h_opt); efa_df$efd_k95_pcs <- as.integer(k_95)
  if (exists("compute_leaf_ect_features", mode = "function")) efa_df <- compute_leaf_ect_features(efa_df, coo_list, opts$contours)
  ect_attrs <- attributes(efa_df)[c("ect_sym_mats", "ect_sym_vecs", "ect_coo_list", "ect_basal_var_exp")]
  exp_pc_cols <- paste0("PC", seq_len(k_export))

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "leaf_zone", "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox", "pixels_per_mm", "is_calibrated_mm",
                 "species_standardized", "scientificName", "species_raw", "determiner_tier", "type_status", "is_type_specimen", "anchor_tier_source", "authority_weight", "spatial_balance_weight", "efd_harmonic_opt", "efd_k95_pcs", exp_pc_cols, "ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "ect_bilateral_asym", "aspect_ratio", "circularity", "solidity", "dissection_index", "archetype_score",
                 "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "is_dissected",
                 "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_px", "centroid_size", "log_centsize", "mask_source")
  lead_cols <- intersect(lead_cols, names(efa_df))
  efa_df <- efa_df[, c(lead_cols, setdiff(names(efa_df), lead_cols))]
  for (anm in intersect(c("ect_sym_mats", "ect_sym_vecs", "ect_coo_list", "ect_basal_var_exp"), names(ect_attrs))) attr(efa_df, anm) <- ect_attrs[[anm]]

  # 5a. Export Disaggregated Individual Leaf Harmonics Archive (Atomic write-then-replace)
  if (!is.null(opts$output_individual) && nzchar(opts$output_individual)) {
    ind_export_df <- if (isTRUE(had_explicit_lz) || !("leaf_zone" %in% names(efa_df)) || !any(efa_df$leaf_zone == "basal", na.rm = TRUE)) efa_df else efa_df[is.na(efa_df$leaf_zone) | efa_df$leaf_zone == "basal", , drop = FALSE]
    write_step03_csv_fast(ind_export_df, opts$output_individual)
    message("Individual leaf EFA harmonics archive exported: ", opts$output_individual, " (Rows: ", nrow(ind_export_df), ")")
  }

  # 5b. Primary Basal Track & Separate Cauline Track Specimen-Level (OTU) Aggregation (Slenker et al. 2022)
  message("Computing primary basal OTU median harmonic vectors and separate cauline morphometric track across intact Tier 1 leaves...")
  valid_closed_mask <- complete.cases(efa_df[, sym_names])
  closed_leaf_df <- efa_df[valid_closed_mask, , drop = FALSE]
  if (!is.null(ect_attrs$ect_sym_vecs)) attr(closed_leaf_df, "ect_sym_vecs") <- ect_attrs$ect_sym_vecs[valid_closed_mask]
  attr(closed_leaf_df, "ect_basal_var_exp") <- ect_attrs$ect_basal_var_exp
  all_harm_cols <- intersect(c(harmonic_names, grep("Chebyshev", names(efa_df), value = TRUE)), names(efa_df))
  specimen_df <- if (exists("aggregate_step03_otu_medians_fast", mode = "function")) {
    agg_tmp <- aggregate_step03_otu_medians_fast(closed_leaf_df, all_harm_cols, sym_names)
    if (exists("apply_step03_heteroblasty_weighted_otu_aggregation", mode = "function")) apply_step03_heteroblasty_weighted_otu_aggregation(agg_tmp, closed_leaf_df, all_harm_cols, sym_names) else agg_tmp
  } else {
    do.call(rbind, lapply(split(closed_leaf_df, paste(closed_leaf_df$catalogNumber, closed_leaf_df$plant_individual_id, sep = "__")), function(grp) {
      lz_g <- tolower(trimws(as.character(grp$leaf_zone))); has_b <- any(!is.na(lz_g) & lz_g == "basal")
      grp_b <- if (has_b) grp[!is.na(lz_g) & lz_g == "basal", , drop = FALSE] else grp; out <- grp_b[1, , drop = FALSE]; w_l <- compute_heteroblasty_leaf_weights(grp_b)
      for (hc in unique(c(all_harm_cols, sym_names))) out[[hc]] <- weighted_otu_stat(grp_b[[hc]], w_l, 6L)
      out$morph_tier_source <- if (!has_b && nrow(grp) >= 1L) "Cauline_Only_Fallback" else "Basal_Rosette"; out
    }))
  }
  if (exists("compute_separate_cauline_otu_track", mode = "function")) specimen_df <- compute_separate_cauline_otu_track(closed_leaf_df, specimen_df, sym_names)
  spec_ect_attrs <- attributes(specimen_df)[c("ect_basal_var_exp", "ect_cauline_var_exp")]
  specimen_df$efd_harmonic_opt <- as.integer(h_opt); specimen_df$efd_k95_pcs <- as.integer(k_95)
  for (p in seq_len(k_export)) { specimen_df[[paste0("PC", p)]] <- NA_real_; specimen_df[[paste0("Sym_PC", p)]] <- 0.0 }
  # Fit OTU-level ordination on Empirical Bayes active reference anchors (with Adaptive Anchor Promotion for 0-Gold/Silver taxa)
  specimen_df$is_type_specimen <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(specimen_df)$is_type_specimen else rep(FALSE, nrow(specimen_df))
  otu_base_mask <- (specimen_df$morph_tier_source != "Cauline_Only_Fallback") & (!as.logical(specimen_df$is_contour_outlier) | as.logical(specimen_df$is_type_specimen))
  if (sum(otu_base_mask, na.rm = TRUE) < 2L) otu_base_mask <- rep(TRUE, nrow(specimen_df))
  otu_h_cols <- intersect(if (exists("active_full_cols") && length(active_full_cols) >= 1L) active_full_cols else sym_names, names(specimen_df))
  otu_ref_mask <- if (sum(otu_base_mask, na.rm = TRUE) >= 2L && length(otu_h_cols) >= 1L) prune_within_taxon_mahalanobis(as.matrix(specimen_df[, otu_h_cols, drop = FALSE]), specimen_df, otu_base_mask, q = 0.95) else otu_base_mask
  specimen_df$is_promoted_bronze <- if (!is.null(attr(otu_ref_mask, "is_promoted_bronze"))) attr(otu_ref_mask, "is_promoted_bronze") else rep(FALSE, nrow(specimen_df))
  specimen_df$anchor_tier_source <- if (!is.null(attr(otu_ref_mask, "anchor_tier_source"))) attr(otu_ref_mask, "anchor_tier_source") else ifelse(specimen_df$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted")
  specimen_df$authority_weight <- compute_authority_precision_weights(specimen_df, active_mask = otu_ref_mask)
  if (sum(otu_ref_mask, na.rm = TRUE) < 2L) otu_ref_mask <- otu_base_mask
  otu_fitted <- FALSE
  if (nrow(specimen_df) >= max(opts$num_pcs, 6L) && exists("active_full_cols") && length(active_full_cols) >= 2L) {
    otu_t1_idx <- which(otu_ref_mask)
    otu_v <- apply(specimen_df[otu_t1_idx, active_full_cols, drop = FALSE], 2, stats::var, na.rm = TRUE)
    otu_cols <- names(otu_v[!is.na(otu_v) & otu_v > 1e-10])
    if (length(otu_cols) >= 2L) {
      w_otu <- specimen_df$authority_weight[otu_t1_idx]; if (sum(w_otu, na.rm = TRUE) <= 0) w_otu <- rep(1.0, length(otu_t1_idx))
      wt_otu <- stats::cov.wt(as.matrix(specimen_df[otu_t1_idx, otu_cols, drop = FALSE]), wt = w_otu, center = TRUE)
      eig_otu <- eigen(wt_otu$cov, symmetric = TRUE); colnames(eig_otu$vectors) <- paste0("PC", seq_len(ncol(eig_otu$vectors))); rownames(eig_otu$vectors) <- otu_cols
      spec_pcs <- sweep(as.matrix(specimen_df[, otu_cols, drop = FALSE]), 2, wt_otu$center, "-") %*% eig_otu$vectors
      pca_fit <- structure(list(sdev = sqrt(pmax(eig_otu$values, 0)), rotation = eig_otu$vectors, center = wt_otu$center, scale = FALSE, x = spec_pcs), class = "prcomp")
      for (p in seq_len(k_export)) specimen_df[[paste0("PC", p)]] <- if (p <= ncol(spec_pcs)) round(spec_pcs[, p], 6) else 0.0
      otu_fitted <- TRUE
    }
  }
  if (!otu_fitted && nrow(specimen_df) >= 1L && exists("pca_fit") && !is.null(pca_fit)) {
    spec_pcs <- stats::predict(pca_fit, newdata = as.matrix(specimen_df[, active_full_cols, drop = FALSE]))
    for (p in seq_len(k_export)) specimen_df[[paste0("PC", p)]] <- if (p <= ncol(spec_pcs)) round(spec_pcs[, p], 6) else 0.0
  }
  if (nrow(specimen_df) >= 1L && !is.null(sym_res$pca)) {
    sym_cols_p <- intersect(rownames(sym_res$pca$rotation), names(specimen_df))
    if (otu_fitted && isTRUE(opts$symmetrize_efa_harmonics)) {
      for (p in seq_len(k_export)) specimen_df[[paste0("Sym_PC", p)]] <- specimen_df[[paste0("PC", p)]]
    } else if (length(sym_cols_p) >= 1L) {
      spec_sym_pcs <- stats::predict(sym_res$pca, newdata = as.matrix(specimen_df[, sym_cols_p, drop = FALSE]))
      for (p in seq_len(k_export)) specimen_df[[paste0("Sym_PC", p)]] <- if (p <= ncol(spec_sym_pcs)) round(spec_sym_pcs[, p], 6) else 0.0
    }
  }

  lead_spec_cols <- intersect(c("catalogNumber", "plant_individual_id", "leaf_count", "basal_leaf_count", "cauline_leaf_count", "morph_tier_source", "foliar_variance", "heteroblasty_index", "heteroblasty_efa_distance", "heteroblasty_ect_distance", "pixels_per_mm", "is_calibrated_mm", "species_standardized", "scientificName", "species_raw", "determiner_tier", "type_status", "is_type_specimen", "anchor_tier_source", "authority_weight", "spatial_balance_weight", "efd_harmonic_opt", "efd_k95_pcs", exp_pc_cols, "Cauline_PC1", "Cauline_PC2", "Cauline_PC3", "ECT_PC1", "ECT_PC2", "ECT_PC3", "Cauline_ECT_PC1", "Cauline_ECT_PC2", "ect_lobe_complexity", "cauline_ect_lobe_complexity", "aspect_ratio", "circularity", "solidity", "q25_solidity", "dissection_index", "dissected_leaf_fraction", "cauline_solidity_med", "cauline_dissection_index_med", "cauline_aspect_ratio_med", "cauline_blade_length_mm", "basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "archetype_score", "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "is_dissected", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_px", "centroid_size", "log_centsize"), names(specimen_df))
  specimen_df <- specimen_df[, c(lead_spec_cols, setdiff(names(specimen_df), lead_spec_cols))]
  for (anm in intersect(c("ect_basal_var_exp", "ect_cauline_var_exp"), names(spec_ect_attrs))) attr(specimen_df, anm) <- spec_ect_attrs[[anm]]
  write_step03_csv_fast(specimen_df, opts$output); message("Specimen-aggregated EFA harmonics table exported: ", opts$output, " (Rows: ", nrow(specimen_df), ")")
  if (exists("persist_step03_to_sqlite", mode = "function")) tryCatch(persist_step03_to_sqlite(efa_df, specimen_df, opts$vouchers), error = function(e) warning("SQLite Step 03 persistence note: ", e$message))
  tryCatch({
    if (!exists("generate_step03_pub_figures", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
    if (exists("generate_step03_pub_figures", mode = "function")) { export_step03_fourier_pc_tables(specimen_df, efa_df, opts$output); generate_step03_pub_figures(efa_df, specimen_df, harm_mat, if (exists("pca_fit")) pca_fit else NULL, opts$figures_dir, opts$output, sym_res$pca) }
  }, error = function(e) warning("Step 03 publication figures/tables export warning: ", e$message))

  return(invisible(efa_df))
}

if (sys.nframe() == 0) run_fourier_extraction(parse_args_robust())
