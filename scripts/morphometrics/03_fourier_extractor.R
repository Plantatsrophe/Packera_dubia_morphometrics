#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_fourier_extractor.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Label-blind Elliptic Fourier Analysis (EFA) on LeafMachine2 leaf silhouettes.
#   Implements 4-tiered extraction hierarchy:
#     - Tier 1: Pristine direct closed outlines (12-harmonic EFA via Momocs)
#     - Tier 2: Hemi-blade bilateral symmetry reflection along midrib
#     - Tier 3: Open Chebyshev orthogonal polynomials (Momocs::opoly)
#     - Tier 4: Whole-rosette dense clump routing to DINOv2 vision embeddings
#   Extracts normalized harmonics (A1-D12), executes PCA (PC1-PC5), merges
#   Darwin Core metadata, and exports data/tables/leaf_efa_harmonics.csv.
#
# Usage:
#   Rscript scripts/morphometrics/03_fourier_extractor.R \
#       --input LM2_Project/Data/output/Packera_dubia_LM2/ \
#       --masks-dir data/masks/ \
#       --vouchers data/tables/curated_vouchers.csv \
#       --output data/tables/leaf_efa_harmonics.csv --harmonics 12 --num-pcs 5
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("Momocs", quietly = TRUE)) library(Momocs)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("tibble", quietly = TRUE)) library(tibble)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})

# ------------------------------------------------------------------------------
# 1. CLI Argument Parsing
# ------------------------------------------------------------------------------
parse_args_robust <- function() {
  option_list <- list(
    optparse::make_option(c("-i", "--input"), type = "character",
      default = "LM2_Project/Data/output/Packera_dubia_LM2/",
      help = "Path to LM2 output directory [default: %default]"),
    optparse::make_option(c("-m", "--masks-dir"), type = "character",
      default = "data/masks/", help = "Path to masks directory [default: %default]"),
    optparse::make_option(c("-v", "--vouchers"), type = "character",
      default = "data/tables/curated_vouchers.csv", help = "Vouchers metadata CSV [default: %default]"),
    optparse::make_option(c("-q", "--qc-table"), type = "character",
      default = "data/tables/leaf_extraction_qc.csv", help = "QC table path [default: %default]"),
    optparse::make_option(c("-o", "--output"), type = "character",
      default = "data/tables/leaf_efa_harmonics.csv", help = "Output CSV [default: %default]"),
    optparse::make_option(c("-k", "--harmonics"), type = "integer",
      default = 12, help = "Number of Fourier harmonics (nb.h) [default: %default]"),
    optparse::make_option(c("-p", "--num-pcs"), type = "integer",
      default = 5, help = "Number of PCA axes to extract [default: %default]")
  )

  if (requireNamespace("optparse", quietly = TRUE)) {
    parser <- optparse::OptionParser(usage = "%prog [options]", option_list = option_list)
    return(optparse::parse_args(parser))
  }

  # Fallback parser if optparse is unavailable
  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    input = "LM2_Project/Data/output/Packera_dubia_LM2/", masks_dir = "data/masks/",
    vouchers = "data/tables/curated_vouchers.csv", qc_table = "data/tables/leaf_extraction_qc.csv",
    output = "data/tables/leaf_efa_harmonics.csv", harmonics = 12, num_pcs = 5
  )
  i <- 1
  while (i <= length(raw_args)) {
    arg <- raw_args[i]
    if (arg %in% c("-i", "--input") && i < length(raw_args)) { opts$input <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-m", "--masks-dir") && i < length(raw_args)) { opts$masks_dir <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-v", "--vouchers") && i < length(raw_args)) { opts$vouchers <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-q", "--qc-table") && i < length(raw_args)) { opts$qc_table <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-o", "--output") && i < length(raw_args)) { opts$output <- raw_args[i + 1]; i <- i + 2 }
    else if (arg %in% c("-k", "--harmonics") && i < length(raw_args)) { opts$harmonics <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else if (arg %in% c("-p", "--num-pcs") && i < length(raw_args)) { opts$num_pcs <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else { i <- i + 1 }
  }
  return(opts)
}

# ------------------------------------------------------------------------------
# 2. Mathematical Engine: Elliptic Fourier Analysis & Orthogonal Polynomials
# ------------------------------------------------------------------------------

#' Computes Normalized Elliptic Fourier Coefficients (Kuhl & Giardina 1982 / Momocs)
compute_efourier_coords <- function(coo, nb_h = 12, norm = TRUE) {
  if (is.null(coo) || nrow(coo) < 5) return(rep(NA_real_, nb_h * 4))
  if (coo[1, 1] != coo[nrow(coo), 1] || coo[1, 2] != coo[nrow(coo), 2]) coo <- rbind(coo, coo[1, ])

  dx <- diff(coo[, 1]); dy <- diff(coo[, 2]); dt <- sqrt(dx^2 + dy^2)
  valid <- dt > 1e-7
  if (sum(valid) < 5) return(rep(NA_real_, nb_h * 4))

  dx <- dx[valid]; dy <- dy[valid]; dt <- dt[valid]
  t_vals <- c(0, cumsum(dt)); T_perim <- t_vals[length(t_vals)]
  if (T_perim <= 0) return(rep(NA_real_, nb_h * 4))

  two_pi_over_T <- 2 * pi / T_perim
  A <- numeric(nb_h); B <- numeric(nb_h); C <- numeric(nb_h); D <- numeric(nb_h)

  for (n in seq_len(nb_h)) {
    coeff <- T_perim / (2 * (n^2) * (pi^2))
    cos_diff <- cos(n * two_pi_over_T * t_vals[-1]) - cos(n * two_pi_over_T * t_vals[-length(t_vals)])
    sin_diff <- sin(n * two_pi_over_T * t_vals[-1]) - sin(n * two_pi_over_T * t_vals[-length(t_vals)])
    A[n] <- coeff * sum((dx / dt) * cos_diff)
    B[n] <- coeff * sum((dx / dt) * sin_diff)
    C[n] <- coeff * sum((dy / dt) * cos_diff)
    D[n] <- coeff * sum((dy / dt) * sin_diff)
  }

  if (!norm) return(as.vector(rbind(A, B, C, D)))

  # Size, rotation, and starting-point normalization
  A1 <- A[1]; B1 <- B[1]; C1 <- C[1]; D1 <- D[1]
  theta1 <- 0.5 * atan2(2 * (A1 * B1 + C1 * D1), (A1^2 + C1^2 - B1^2 - D1^2))
  A1_star <- A1 * cos(theta1) + B1 * sin(theta1)
  C1_star <- C1 * cos(theta1) + D1 * sin(theta1)
  psi1 <- atan2(C1_star, A1_star)
  E1 <- sqrt(A1_star^2 + C1_star^2)
  if (E1 < 1e-7) return(rep(NA_real_, nb_h * 4))

  cos_psi1 <- cos(psi1); sin_psi1 <- sin(psi1)
  R_psi <- matrix(c(cos_psi1, sin_psi1, -sin_psi1, cos_psi1), nrow = 2, byrow = TRUE)

  norm_A <- numeric(nb_h); norm_B <- numeric(nb_h)
  norm_C <- numeric(nb_h); norm_D <- numeric(nb_h)

  for (n in seq_len(nb_h)) {
    M_n <- matrix(c(A[n], B[n], C[n], D[n]), nrow = 2, byrow = TRUE)
    cos_nt <- cos(n * theta1); sin_nt <- sin(n * theta1)
    R_nt <- matrix(c(cos_nt, -sin_nt, sin_nt, cos_nt), nrow = 2, byrow = TRUE)
    norm_M <- (1 / E1) * (R_psi %*% M_n %*% R_nt)
    norm_A[n] <- norm_M[1, 1]; norm_B[n] <- norm_M[1, 2]
    norm_C[n] <- norm_M[2, 1]; norm_D[n] <- norm_M[2, 2]
  }
  return(as.vector(rbind(norm_A, norm_B, norm_C, norm_D)))
}

#' Computes Open Chebyshev Orthogonal Polynomials (Momocs::opoly equivalent)
compute_chebyshev_opoly <- function(x_norm, y_norm, degree = 5) {
  if (length(x_norm) < degree + 1) return(rep(NA_real_, degree + 1))
  u <- 2 * (x_norm - min(x_norm)) / max(max(x_norm) - min(x_norm), 1e-6) - 1
  T_mat <- matrix(0, nrow = length(u), ncol = degree + 1)
  T_mat[, 1] <- 1
  if (degree >= 1) T_mat[, 2] <- u
  if (degree >= 2) {
    for (d in 3:(degree + 1)) T_mat[, d] <- 2 * u * T_mat[, d - 1] - T_mat[, d - 2]
  }
  fit <- tryCatch(lm.fit(T_mat, y_norm), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, degree + 1))
  return(fit$coefficients)
}

#' Extract boundary contour coordinates from binary mask file
extract_contour_coords_from_file <- function(filepath, num_samples = 120) {
  if (!file.exists(filepath)) return(NULL)
  bin <- NULL
  if (requireNamespace("png", quietly = TRUE)) {
    img <- png::readPNG(filepath)
    if (length(dim(img)) == 3) img <- img[, , 1]
    bin <- (img > 0.05)
  } else if (requireNamespace("magick", quietly = TRUE)) {
    im <- magick::image_read(filepath)
    im_gray <- magick::image_convert(im, type = "grayscale")
    raw_dat <- as.numeric(magick::image_data(im_gray, channels = "gray"))
    dim(raw_dat) <- c(magick::image_info(im)$height, magick::image_info(im)$width)
    bin <- (raw_dat > 12)
  }
  if (is.null(bin)) return(NULL)
  h <- nrow(bin); w <- ncol(bin)
  if (sum(bin) < 25) return(NULL)

  padded <- matrix(FALSE, nrow = h + 2, ncol = w + 2)
  padded[2:(h + 1), 2:(w + 1)] <- bin
  edge_x <- c(); edge_y <- c()
  for (c in 2:(w + 1)) {
    ys <- which(padded[, c])
    if (length(ys) > 0) { edge_x <- c(edge_x, c - 1, c - 1); edge_y <- c(edge_y, min(ys) - 1, max(ys) - 1) }
  }
  if (length(edge_x) < 10) return(NULL)
  cx <- mean(edge_x); cy <- mean(edge_y)
  angles <- atan2(edge_y - cy, edge_x - cx)
  coo <- cbind(edge_x[order(angles)], edge_y[order(angles)])
  if (nrow(coo) > num_samples) {
    coo <- coo[as.integer(seq(1, nrow(coo), length.out = num_samples)), ]
  }
  return(coo)
}

# ------------------------------------------------------------------------------
# 3. 4-Tier Hierarchy Orchestrator & Metadata Integration
# ------------------------------------------------------------------------------
run_fourier_extraction <- function(opts) {
  message("==================================================================")
  message("Starting Label-Blind Elliptic Fourier Analysis Pipeline")
  message("Harmonics (nb.h): ", opts$harmonics, " | Target PCA Dims: ", opts$num_pcs)
  message("Masks Dir: ", opts$masks_dir, " | Output: ", opts$output)
  message("==================================================================")

  vouchers_df <- if (file.exists(opts$vouchers)) read.csv(opts$vouchers, stringsAsFactors = FALSE) else NULL
  if (!is.null(vouchers_df)) message("Loaded ", nrow(vouchers_df), " curated voucher records.")

  masks_dir <- opts$masks_dir
  t1_dir <- file.path(masks_dir, "tier1_intact")
  t2_dir <- file.path(masks_dir, "tier2_reflected")
  t3_dir <- file.path(masks_dir, "tier3_open_curves")
  t4_dir <- file.path(masks_dir, "rosettes_dense")

  if (!dir.exists(t1_dir)) {
    t1_dir <- file.path(opts$input, "Plant_Components", "Segmentation_Whole_Leaves")
    if (!dir.exists(t1_dir)) t1_dir <- file.path(opts$input, "crops", "leaves")
  }

  t1_files <- if (dir.exists(t1_dir)) list.files(t1_dir, pattern = "\\.(png|jpg)$", full.names = TRUE) else character(0)
  t2_files <- if (dir.exists(t2_dir)) list.files(t2_dir, pattern = "\\.(png|jpg)$", full.names = TRUE) else character(0)
  t3_files <- if (dir.exists(t3_dir)) list.files(t3_dir, pattern = "\\.csv$", full.names = TRUE) else character(0)
  t4_files <- if (dir.exists(t4_dir)) list.files(t4_dir, pattern = "\\.(png|jpg)$", full.names = TRUE) else character(0)

  message(sprintf("Discovered: %d T1 (Direct), %d T2 (Reflected), %d T3 (Open), %d T4 (Rosette)",
                  length(t1_files), length(t2_files), length(t3_files), length(t4_files)))

  records <- list()
  harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(opts$harmonics), function(x, y) paste0(x, y)))

  parse_specimen_info <- function(filename) {
    base <- sub("\\.[^.]+$", "", basename(filename))
    cat_num <- strsplit(base, "(_p|_leaf|_reflected|_curve)")[[1]][1]
    p_id <- 0; l_id <- 1
    m_p <- regmatches(base, regexpr("_p(\\d+)", base))
    if (length(m_p) > 0) p_id <- as.integer(sub("_p", "", m_p))
    m_l <- regmatches(base, regexpr("leaf(\\d+)|leaf_(\\d+)|_(\\d+)$", base))
    if (length(m_l) > 0) l_id <- as.integer(gsub("[^0-9]", "", m_l))
    list(catalogNumber = cat_num, plant_id = p_id, leaf_id = l_id, base = base)
  }

  # --- TIER 1 DIRECT CLOSED OUTLINES ---
  message("Processing Tier 1 Direct Closed Outlines...")
  for (mf in t1_files) {
    info <- parse_specimen_info(mf)
    coo <- extract_contour_coords_from_file(mf, num_samples = 120)
    h_vals <- if (!is.null(coo)) compute_efourier_coords(coo, nb_h = opts$harmonics, norm = TRUE) else rep(NA_real_, opts$harmonics * 4)
    area_px <- if (!is.null(coo)) abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1]))) else 0
    span_x <- if (!is.null(coo)) diff(range(coo[, 1])) else 0; span_y <- if (!is.null(coo)) diff(range(coo[, 2])) else 0
    aspect_ratio <- if (span_y > 0) round(span_x / span_y, 4) else 1.0

    rec <- list(catalogNumber = info$catalogNumber, plant_individual_id = info$plant_id,
                leaf_id = info$leaf_id, assigned_tier = "Tier_1_Direct",
                aspect_ratio = aspect_ratio, area_px = round(area_px, 1), mask_source = mf)
    for (k in seq_along(harmonic_names)) rec[[harmonic_names[k]]] <- h_vals[k]
    records[[length(records) + 1]] <- rec
  }

  # --- TIER 2 BILATERAL REFLECTED OUTLINES ---
  message("Processing Tier 2 Bilateral Reflected Outlines...")
  for (mf in t2_files) {
    info <- parse_specimen_info(mf)
    coo <- extract_contour_coords_from_file(mf, num_samples = 120)
    h_vals <- if (!is.null(coo)) compute_efourier_coords(coo, nb_h = opts$harmonics, norm = TRUE) else rep(NA_real_, opts$harmonics * 4)
    area_px <- if (!is.null(coo)) abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1]))) else 0
    span_x <- if (!is.null(coo)) diff(range(coo[, 1])) else 0; span_y <- if (!is.null(coo)) diff(range(coo[, 2])) else 0
    aspect_ratio <- if (span_y > 0) round(span_x / span_y, 4) else 1.0

    rec <- list(catalogNumber = info$catalogNumber, plant_individual_id = info$plant_id,
                leaf_id = info$leaf_id, assigned_tier = "Tier_2_Reflected",
                aspect_ratio = aspect_ratio, area_px = round(area_px, 1), mask_source = mf)
    for (k in seq_along(harmonic_names)) rec[[harmonic_names[k]]] <- h_vals[k]
    records[[length(records) + 1]] <- rec
  }

  # --- TIER 3 OPEN MARGIN CURVES ---
  message("Processing Tier 3 Open Margin Curves...")
  for (cf in t3_files) {
    info <- parse_specimen_info(cf)
    curve_data <- tryCatch(read.csv(cf, stringsAsFactors = FALSE), error = function(e) NULL)
    opoly_coeffs <- if (!is.null(curve_data) && nrow(curve_data) >= 10 && "x_norm" %in% names(curve_data)) {
      compute_chebyshev_opoly(curve_data$x_norm, curve_data$y_norm, degree = 5)
    } else rep(NA_real_, 6)

    rec <- list(catalogNumber = info$catalogNumber, plant_individual_id = info$plant_id,
                leaf_id = info$leaf_id, assigned_tier = "Tier_3_OpenCurve",
                aspect_ratio = NA_real_, area_px = NA_real_, mask_source = cf)
    for (k in seq_along(harmonic_names)) rec[[harmonic_names[k]]] <- NA_real_
    for (d in seq_along(opoly_coeffs)) rec[[paste0("Chebyshev_T", d - 1)]] <- round(opoly_coeffs[d], 6)
    records[[length(records) + 1]] <- rec
  }

  # --- TIER 4 WHOLE ROSETTE CLUMP ROUTING ---
  message("Processing Tier 4 Dense Rosette Clump Routing...")
  for (rf in t4_files) {
    info <- parse_specimen_info(rf)
    rec <- list(catalogNumber = info$catalogNumber, plant_individual_id = info$plant_id,
                leaf_id = 0, assigned_tier = "Tier_4_Rosette",
                aspect_ratio = NA_real_, area_px = NA_real_, mask_source = rf)
    for (k in seq_along(harmonic_names)) rec[[harmonic_names[k]]] <- NA_real_
    records[[length(records) + 1]] <- rec
  }

  efa_df <- do.call(rbind, lapply(records, as.data.frame, stringsAsFactors = FALSE))
  message("Total leaves & rosettes harmonized: ", nrow(efa_df))

  # ----------------------------------------------------------------------------
  # 4. Morphospace PCA Decomposition (PC1-PC5)
  # ----------------------------------------------------------------------------
  message("Executing Morphospace PCA Decomposition (PC1-PC5)...")
  closed_idx <- which(efa_df$assigned_tier %in% c("Tier_1_Direct", "Tier_2_Reflected"))
  harm_mat <- as.matrix(efa_df[closed_idx, harmonic_names])

  var_cols <- apply(harm_mat, 2, function(col) stats::var(col, na.rm = TRUE))
  keep_cols <- which(!is.na(var_cols) & var_cols > 1e-8)

  if (length(keep_cols) >= opts$num_pcs && length(closed_idx) > opts$num_pcs) {
    complete_cases <- complete.cases(harm_mat[, keep_cols])
    pca_fit <- stats::prcomp(harm_mat[complete_cases, keep_cols], center = TRUE, scale. = TRUE)
    var_explained <- round((pca_fit$sdev^2) / sum(pca_fit$sdev^2) * 100, 2)
    message("=== Morphospace PCA Variance Explained ===")
    for (p in seq_len(min(opts$num_pcs, length(var_explained)))) {
      message(sprintf("  PC%d: %5.2f%% variance", p, var_explained[p]))
    }
    for (p in seq_len(opts$num_pcs)) efa_df[[paste0("PC", p)]] <- NA_real_
    valid_sub_idx <- closed_idx[complete_cases]
    for (p in seq_len(opts$num_pcs)) {
      if (p <= ncol(pca_fit$x)) efa_df[valid_sub_idx, paste0("PC", p)] <- round(pca_fit$x[, p], 6)
    }
  } else {
    for (p in seq_len(opts$num_pcs)) efa_df[[paste0("PC", p)]] <- NA_real_
  }

  # ----------------------------------------------------------------------------
  # 5. Darwin Core Metadata Join & Export
  # ----------------------------------------------------------------------------
  if (!is.null(vouchers_df)) {
    meta_cols <- intersect(names(vouchers_df), c("catalogNumber", "species_raw", "determiner_raw",
      "determiner_tier", "county", "stateProvince", "latitude", "longitude",
      "pheno_sin", "pheno_cos", "regional_group"))
    vouchers_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), meta_cols, drop = FALSE]
    efa_df <- merge(efa_df, vouchers_sub, by = "catalogNumber", all.x = TRUE)
  }

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier",
                 "species_raw", "determiner_tier", "PC1", "PC2", "PC3", "PC4", "PC5")
  lead_cols <- intersect(lead_cols, names(efa_df))
  efa_df <- efa_df[, c(lead_cols, setdiff(names(efa_df), lead_cols))]

  dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
  write.csv(efa_df, file = opts$output, row.names = FALSE, na = "")
  message("Master EFA harmonics table successfully exported: ", opts$output, " (Rows: ", nrow(efa_df), ")")
  return(invisible(efa_df))
}

if (sys.nframe() == 0) {
  opts <- parse_args_robust()
  run_fourier_extraction(opts)
}
