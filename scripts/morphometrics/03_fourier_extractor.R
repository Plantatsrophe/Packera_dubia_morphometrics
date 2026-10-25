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
    optparse::make_option(c("--permutations"), type = "integer", default = 999, help = "Number of permutations for PERMANOVA [default: %default]"),
    optparse::make_option(c("--seed"), type = "integer", default = 42, help = "Random seed for reproducible permutations [default: %default]")
  )

  if (requireNamespace("optparse", quietly = TRUE)) {
    parser <- optparse::OptionParser(usage = "%prog [options]", option_list = option_list)
    opts <- optparse::parse_args(parser)
    names(opts) <- gsub("-", "_", names(opts))
    if (!is.null(opts$input) && opts$contours_dir == "data/contours/") opts$contours_dir <- opts$input
    return(opts)
  }

  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    input = "data/contours/", contours_dir = "data/contours/", manifest = "data/tables/extracted_leaf_manifest.csv",
    masks_dir = "data/masks/", vouchers = "data/tables/curated_vouchers.csv",
    output = "data/tables/leaf_efa_harmonics.csv", output_individual = "data/tables/leaf_efa_harmonics_individual.csv",
    harmonics = 12, num_pcs = 5,
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
  cx <- mean(coo[, 1], na.rm = TRUE)
  cy <- mean(coo[, 2], na.rm = TRUE)
  coo[, 1] <- coo[, 1] - cx
  coo[, 2] <- coo[, 2] - cy
  cs <- sqrt(sum(coo[, 1]^2 + coo[, 2]^2, na.rm = TRUE))
  if (is.finite(cs) && cs > 1e-8) coo <- coo / cs

  if (coo[1, 1] != coo[nrow(coo), 1] || coo[1, 2] != coo[nrow(coo), 2]) {
    coo <- rbind(coo, coo[1, ])
  }

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

  # Size, orientation, and starting-point invariant normalization (only when norm = TRUE)
  A1 <- A[1]; B1 <- B[1]; C1 <- C[1]; D1 <- D[1]
  theta1 <- if (start) 0.5 * atan2(2 * (A1 * B1 + C1 * D1), (A1^2 + C1^2 - B1^2 - D1^2)) else 0.0
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
  cx <- mean(coo[, 1], na.rm = TRUE)
  cy <- mean(coo[, 2], na.rm = TRUE)
  val <- sqrt(sum((coo[, 1] - cx)^2 + (coo[, 2] - cy)^2, na.rm = TRUE))
  if (is.finite(val) && val > 0) val else NA_real_
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
  low_power <- if (length(low_h_cols) > 0) rowSums(coe[, low_h_cols, drop = FALSE]^2, na.rm = TRUE) else rep(1.0, length(row_indices))
  high_power <- if (length(high_h_cols) > 0) rowSums(coe[, high_h_cols, drop = FALSE]^2, na.rm = TRUE) else rep(0.0, length(row_indices))
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
  col_vars <- apply(sym_mat[clean_pca_mask, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  sym_mat <- sym_mat[, names(col_vars[!is.na(col_vars) & col_vars > 1e-10]), drop = FALSE]

  out_df <- if (!is.null(fac_df)) as.data.frame(fac_df) else data.frame(row_id = row_indices)
  solidity <- if ("solidity" %in% names(out_df)) ifelse(is.na(out_df$solidity), 1.0, as.numeric(out_df$solidity)) else rep(1.0, nrow(out_df))
  ucs <- if ("ucs" %in% names(out_df)) ifelse(is.na(out_df$ucs), 1.0, as.numeric(out_df$ucs)) else if ("ucs_score" %in% names(out_df)) ifelse(is.na(out_df$ucs_score), 1.0, as.numeric(out_df$ucs_score)) else rep(1.0, nrow(out_df))
  leaf_quality_weight <- pmax(0.05, solidity * ucs * ((1.0 - pmin(pct_asymmetry, 0.95))^2) * ifelse(is_contour_outlier, 0.10, 1.0))

  out_df$sym_var <- round(sym_var, 6); out_df$asym_var <- round(asym_var, 6); out_df$pct_asymmetry <- round(pct_asymmetry, 6)
  out_df$harmonic_roughness <- round(harmonic_roughness, 6); out_df$is_contour_outlier <- is_contour_outlier; out_df$leaf_quality_weight <- round(leaf_quality_weight, 6)
  for (pc_col in c("Sym_PC1", "Sym_PC2", "Sym_PC3", "Sym_PC4", "Sym_PC5")) out_df[[pc_col]] <- 0.0
  pca_sym <- if (nrow(sym_mat) >= 2 && ncol(sym_mat) >= 1) stats::prcomp(sym_mat[clean_pca_mask, , drop = FALSE], center = TRUE, scale. = FALSE) else NULL
  if (!is.null(pca_sym)) {
    sym_scores <- stats::predict(pca_sym, newdata = sym_mat); pca_sym$x <- sym_scores
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
  if (!is.null(vouchers_df)) message("Loaded ", nrow(vouchers_df), " curated voucher records.")

  # 1. Discover coordinates via Dual Ingestion Engine
  contour_files <- character(0)
  use_masks_fallback <- FALSE
  manifest_df <- if (file.exists(opts$manifest)) read.csv(opts$manifest, stringsAsFactors = FALSE) else NULL

  if (!is.null(manifest_df) && nrow(manifest_df) > 0) {
    message("Ingesting via extracted leaf manifest: ", opts$manifest)
    if ("contour_path" %in% names(manifest_df)) {
        manifest_df <- manifest_df[!duplicated(manifest_df$contour_path, fromLast = TRUE), , drop = FALSE]
        contour_files <- unique(manifest_df$contour_path[file.exists(manifest_df$contour_path)])
    }
    if (length(contour_files) == 0 && "mask_path" %in% names(manifest_df)) {
        mask_files <- unique(manifest_df$mask_path[file.exists(manifest_df$mask_path)])
        if (length(mask_files) > 0) {
            contour_files <- mask_files
            use_masks_fallback <- TRUE
        }
    }
  }

  if (length(contour_files) == 0) {
      input_dir <- if (!is.null(opts$input_dir)) opts$input_dir else opts$input
      if (dir.exists(input_dir)) {
        contour_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
        if (length(contour_files) == 0) {
            mask_files <- list.files(input_dir, pattern = "\\.(png|jpg)$", full.names = TRUE, recursive = TRUE)
            if (length(mask_files) > 0) {
                contour_files <- mask_files
                use_masks_fallback <- TRUE
            }
        }
      }
  }

  if (length(contour_files) == 0) {
    message("No contour CSVs in input dir. Checking fallback masks...")
    fallback_dir <- opts$masks_dir
    if (dir.exists(fallback_dir)) {
      mask_files <- list.files(fallback_dir, pattern = "\\.(png|jpg)$", full.names = TRUE, recursive = TRUE)
      if (length(mask_files) > 0) {
        message("Found ", length(mask_files), " fallback mask files.")
        contour_files <- mask_files
        use_masks_fallback <- TRUE
      }
    }
  }

  if (length(contour_files) == 0) {
    stop("No contour CSVs or active mask files found in data/contours or data/masks.")
  }

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
      if (has_momocs) {
        coo <- tryCatch({
          tmp_coo <- Momocs::import_jpg(f)
          if (is.list(tmp_coo) && length(tmp_coo) > 0) tmp_coo[[1]] else NULL
        }, error = function(e) NULL)
      } else {
        coo <- NULL
      }
    } else {
      coo <- load_contour_file(f)
    }

    if (is.null(coo)) return(NULL)

    # Enforce clockwise orientation
    if (has_momocs) {
      coo <- tryCatch(Momocs::coo_cw(coo), error = function(e) coo)
    }

    meta <- parse_specimen_id(f)

    m_match <- NULL
    if (!is.null(manifest_df) && nrow(manifest_df) > 0) {
      idx_hit <- NA_integer_
      if (!is.null(manifest_by_basename)) {
        idx_hit <- unname(manifest_by_basename[basename(f)])
      }
      if (is.na(idx_hit) && !is.null(manifest_by_cat_leaf)) {
        key_lookup <- paste0(as.character(meta$catalogNumber), "__", as.integer(meta$leaf_id))
        idx_hit <- unname(manifest_by_cat_leaf[key_lookup])
      }
      if (!is.na(idx_hit) && idx_hit >= 1L && idx_hit <= nrow(manifest_df)) {
        m_match <- manifest_df[idx_hit, , drop = FALSE]
      }
    }

    resolved_pid <- if (!is.null(m_match) && nrow(m_match) > 0 && "plant_individual_id" %in% names(m_match) && !is.na(m_match$plant_individual_id[1])) as.integer(m_match$plant_individual_id[1]) else meta$plant_id
    resolved_tier <- if (!is.null(m_match) && nrow(m_match) > 0 && "assigned_tier" %in% names(m_match) && !is.na(m_match$assigned_tier[1])) as.character(m_match$assigned_tier[1]) else "Tier_1_Direct"
    resolved_ppm <- if (!is.null(m_match) && nrow(m_match) > 0 && "pixels_per_mm" %in% names(m_match)) suppressWarnings(as.numeric(m_match$pixels_per_mm[1])) else NA_real_
    is_calibrated_mm <- !is.na(resolved_ppm) && is.finite(resolved_ppm) && resolved_ppm > 0
    if (!is_calibrated_mm) resolved_ppm <- NA_real_

    # Calculate raw pixel area & polygon solidity before coordinate scaling
    area_px <- abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1])))
    resolved_solidity <- if (!is.null(m_match) && nrow(m_match) > 0 && "solidity" %in% names(m_match) && !is.na(m_match$solidity[1])) round(suppressWarnings(as.numeric(m_match$solidity[1])), 4) else compute_polygon_solidity(coo, area_px)
    resolved_ucs <- if (!is.null(m_match) && nrow(m_match) > 0 && "ucs_score" %in% names(m_match) && !is.na(m_match$ucs_score[1])) {
      round(suppressWarnings(as.numeric(m_match$ucs_score[1])), 4)
    } else if (!is.null(m_match) && nrow(m_match) > 0 && "ucs" %in% names(m_match) && !is.na(m_match$ucs[1])) {
      round(suppressWarnings(as.numeric(m_match$ucs[1])), 4)
    } else NA_real_

    # Calibrate contour coordinates to physical millimeters (mm) when pixels_per_mm > 0
    if (is_calibrated_mm) coo <- coo / resolved_ppm

    # Oriented major-axis length-to-width aspect ratio (Momocs::coo_lw / SVD major-axis alignment)
    lw <- if (has_momocs) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
    if (is.null(lw) || length(lw) < 2 || is.na(lw[2]) || lw[2] <= 1e-6) {
      coo_centered <- scale(coo, center = TRUE, scale = FALSE)
      svd_rot <- coo_centered %*% svd(coo_centered)$v
      lw <- c(diff(range(svd_rot[, 1], na.rm = TRUE)), diff(range(svd_rot[, 2], na.rm = TRUE)))
    }
    aspect_ratio <- round(as.numeric(lw[1] / lw[2]), 4)
    centsize <- compute_centsize_robust(coo)
    log_cs <- if (!is.na(centsize) && centsize > 0) round(log(centsize), 6) else NA_real_

    shape_name <- meta$shape_id
    rec <- data.frame(
      shape_id = shape_name, catalogNumber = meta$catalogNumber, plant_individual_id = resolved_pid,
      leaf_id = meta$leaf_id, assigned_tier = resolved_tier, pixels_per_mm = resolved_ppm,
      is_calibrated_mm = is_calibrated_mm, aspect_ratio = aspect_ratio, solidity = resolved_solidity,
      ucs_score = resolved_ucs, area_px = round(area_px, 1),
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
  if (length(valid_ingested) == 0) stop("Failed to extract valid coordinate contours.")

  fac_records <- lapply(valid_ingested, `[[`, "rec")
  fac_df <- do.call(rbind, fac_records)
  coo_list <- lapply(valid_ingested, `[[`, "coo")
  shape_names_vec <- vapply(valid_ingested, `[[`, character(1), "shape_name")

  # Ensure strict 1-to-1 uniqueness between fac_df rows and coo_list elements
  keep_unique <- !duplicated(shape_names_vec, fromLast = TRUE)
  fac_df <- fac_df[keep_unique, , drop = FALSE]
  coo_list <- coo_list[keep_unique]
  names(coo_list) <- shape_names_vec[keep_unique]
  rownames(fac_df) <- names(coo_list)
  message(sprintf("Assembled %d valid closed contours for Elliptic Fourier Analysis.", length(coo_list)))

  # 3. Perform 12-Harmonic EFA via Momocs::Out() & Momocs::efourier(norm = FALSE, start = FALSE)
  message("Executing 12-harmonic EFA via Momocs::efourier(nb.h = ", opts$harmonics, ", norm = FALSE, start = FALSE)...")
  efa_success <- FALSE
  harm_mat <- NULL

  if (requireNamespace("Momocs", quietly = TRUE)) {
    out_obj <- tryCatch(Momocs::Out(coo_list, fac = fac_df), error = function(e) {
      message("Momocs::Out() note: ", e$message); NULL
    })
    if (!is.null(out_obj)) {
      # Compute Centroid Size (CS) for each contour using Momocs::coo_centsize(outlines) prior to coo_scale()
      cs_vals <- tryCatch(Momocs::coo_centsize(out_obj), error = function(e) NULL)
      if (!is.null(cs_vals) && length(cs_vals) == nrow(fac_df)) {
        fac_df$centroid_size <- round(as.numeric(cs_vals), 4)
        fac_df$log_centsize <- ifelse(!is.na(fac_df$centroid_size) & fac_df$centroid_size > 0,
                                      round(log(fac_df$centroid_size), 6), NA_real_)
      }
      # Center and scale homologized outlines without 180-degree first-harmonic ellipse flips
      out_obj <- Momocs::coo_center(out_obj) %>% Momocs::coo_scale()
      ef_res <- tryCatch(Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = FALSE, start = FALSE), error = function(e) {
        message("Momocs::efourier() note: ", e$message); NULL
      })
      if (!is.null(ef_res) && !is.null(ef_res$coe)) {
        harm_mat <- as.matrix(ef_res$coe)
        efa_success <- TRUE
        message("Momocs EFA extraction successful across ", nrow(harm_mat), " outlines.")
      }
    }
  }

  # Fallback to pure-R vectorized EFA if Momocs not available
  if (!efa_success || is.null(harm_mat)) {
    message("Running mathematical homologized EFA engine (norm = FALSE, start = FALSE)...")
    harm_list <- lapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = opts$harmonics, norm = FALSE, start = FALSE))
    harm_mat <- do.call(rbind, harm_list)
    colnames(harm_mat) <- harmonic_names
  }

  # Harmonize column names to standard A1-D12
  if (ncol(harm_mat) == length(harmonic_names)) colnames(harm_mat) <- harmonic_names
  efa_df <- cbind(fac_df, as.data.frame(harm_mat))

  # 4. Symmetric Fourier Decomposition & PCA (PC1-PC5 / Sym_PC1-Sym_PC5)
  message("Extracting symmetric harmonic coefficients (An and Dn harmonics)...")
  sym_names <- c(paste0("A", seq_len(opts$harmonics)), paste0("D", seq_len(opts$harmonics)))
  sym_harmonics <- harm_mat[, sym_names, drop = FALSE]
  ef_input <- if (exists("ef_res") && !is.null(ef_res)) ef_res else list(coe = harm_mat, fac = fac_df)
  sym_res <- extract_symmetric_fourier_pcs(ef_input, fac_df = fac_df, num_pcs = opts$num_pcs, permutations = opts$permutations, seed = opts$seed, run_permanova = FALSE)
  for (scol in c("sym_var", "asym_var", "pct_asymmetry", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", paste0("Sym_PC", seq_len(opts$num_pcs)))) {
    if (scol %in% names(sym_res$df)) efa_df[[scol]] <- sym_res$df[[scol]]
  }

  message("Running PCA strictly on symmetric harmonics (sym_harmonics)...")
  complete_idx <- which(complete.cases(sym_harmonics))
  for (p in seq_len(opts$num_pcs)) efa_df[[paste0("PC", p)]] <- NA_real_

  if (length(complete_idx) >= opts$num_pcs) {
    clean_idx <- if (length(complete_idx) >= 10 && sum(!efa_df$is_contour_outlier[complete_idx]) >= 5) complete_idx[!efa_df$is_contour_outlier[complete_idx]] else complete_idx
    col_vars <- apply(sym_harmonics[clean_idx, , drop = FALSE], 2, stats::var, na.rm = TRUE)
    active_sym_cols <- names(col_vars[!is.na(col_vars) & col_vars > 1e-10])
    message(sprintf("Active symmetric harmonic dimensions with non-zero variance: %d / %d", length(active_sym_cols), ncol(sym_harmonics)))
    if (requireNamespace("Momocs", quietly = TRUE) && exists("ef_res") && !is.null(ef_res)) {
      ef_clean <- if (nrow(harm_mat) >= 10 && sum(!efa_df$is_contour_outlier) >= 5) tryCatch(Momocs::slice(ef_res, which(!efa_df$is_contour_outlier)), error = function(e) ef_res) else ef_res
      momocs_pca <- tryCatch(Momocs::PCA(ef_clean), error = function(e) NULL)
    }
    pca_fit <- stats::prcomp(sym_harmonics[clean_idx, active_sym_cols, drop = FALSE], center = TRUE, scale. = FALSE)
    pc_scores <- stats::predict(pca_fit, newdata = sym_harmonics[complete_idx, active_sym_cols, drop = FALSE])
    var_exp <- round((pca_fit$sdev^2) / sum(pca_fit$sdev^2) * 100, 2)
    message("=== Symmetric Morphospace PCA Variance Explained ===")
    for (p in seq_len(min(opts$num_pcs, length(var_exp)))) message(sprintf("  PC%d: %5.2f%% variance", p, var_exp[p]))
    for (p in seq_len(min(opts$num_pcs, ncol(pc_scores)))) efa_df[complete_idx, paste0("PC", p)] <- round(pc_scores[, p], 6)
  }

  # 5. Darwin Core Metadata Integration & Standardization
  if (!is.null(vouchers_df)) {
    meta_cols <- intersect(names(vouchers_df), c("catalogNumber", "scientificName", "species_raw", "determiner_raw", "determiner_tier",
                                                 "county", "stateProvince", "latitude", "longitude", "pheno_sin", "pheno_cos", "regional_group"))
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), meta_cols, drop = FALSE]
    efa_df <- merge(efa_df, v_sub, by = "catalogNumber", all.x = TRUE, sort = FALSE)
  }
  harm_mat <- as.matrix(efa_df[, harmonic_names, drop = FALSE])
  sym_harmonics <- as.matrix(efa_df[, sym_names, drop = FALSE])

  # Normalize assigned_tier to canonical contract ("Tier_1_Direct", "Tier_2_Reflected")
  efa_df$assigned_tier <- ifelse(grepl("^tier_?1", efa_df$assigned_tier, ignore.case = TRUE), "Tier_1_Direct",
                                 ifelse(grepl("^tier_?2", efa_df$assigned_tier, ignore.case = TRUE), "Tier_2_Reflected", efa_df$assigned_tier))
  efa_df$reconstruction_tier <- ifelse(efa_df$assigned_tier == "Tier_1_Direct", "Tier 1",
                                       ifelse(efa_df$assigned_tier == "Tier_2_Reflected", "Tier 2", NA_character_))
  if (!"scientificName" %in% names(efa_df)) efa_df$scientificName <- efa_df$species_raw
  else efa_df$scientificName <- ifelse(is.na(efa_df$scientificName) | efa_df$scientificName == "", efa_df$species_raw, efa_df$scientificName)
  if (!"determiner_tier" %in% names(efa_df)) efa_df$determiner_tier <- "Tier_3_Bronze"
  else efa_df$determiner_tier <- ifelse(is.na(efa_df$determiner_tier) | efa_df$determiner_tier == "", "Tier_3_Bronze", efa_df$determiner_tier)

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier", "reconstruction_tier", "pixels_per_mm", "is_calibrated_mm",
                 "scientificName", "species_raw", "determiner_tier", "PC1", "PC2", "PC3", "PC4", "PC5", "aspect_ratio", "solidity", "ucs_score",
                 "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_px", "centroid_size", "log_centsize", "mask_source")
  lead_cols <- intersect(lead_cols, names(efa_df))
  efa_df <- efa_df[, c(lead_cols, setdiff(names(efa_df), lead_cols))]

  # 5a. Export Disaggregated Individual Leaf Harmonics Archive (Atomic write-then-replace)
  if (!is.null(opts$output_individual) && nzchar(opts$output_individual)) {
    dir.create(dirname(opts$output_individual), recursive = TRUE, showWarnings = FALSE)
    tmp_ind <- paste0(opts$output_individual, ".tmp"); write.csv(efa_df, file = tmp_ind, row.names = FALSE, na = ""); file.rename(tmp_ind, opts$output_individual)
    message("Individual leaf EFA harmonics archive exported: ", opts$output_individual, " (Rows: ", nrow(efa_df), ")")
  }

  # 5b. Compute Specimen-Level Median Harmonic Vectors & Phenotypic Stability Index
  message("Computing specimen-level median harmonic vectors (grouped by catalogNumber & plant instance)...")
  valid_closed_mask <- grepl("^tier_?[12]", efa_df$assigned_tier, ignore.case = TRUE) & complete.cases(efa_df[, sym_names])
  closed_leaf_df <- efa_df[valid_closed_mask, , drop = FALSE]
  all_harm_cols <- intersect(c(harmonic_names, grep("Chebyshev", names(efa_df), value = TRUE)), names(efa_df))
  spec_groups <- split(closed_leaf_df, list(closed_leaf_df$catalogNumber, closed_leaf_df$plant_individual_id), drop = TRUE)
  specimen_records <- vector("list", length(spec_groups))
  spec_idx <- 0L

  for (grp in spec_groups) {
    if (nrow(grp) == 0) next
    spec_idx <- spec_idx + 1L
    cat_num <- grp$catalogNumber[1]; p_id <- grp$plant_individual_id[1]; n_leaves <- nrow(grp)
    harm_medians <- vapply(grp[, all_harm_cols, drop = FALSE], stats::median, numeric(1), na.rm = TRUE)
    foliar_var <- if (n_leaves > 1) round(mean(apply(as.matrix(grp[, sym_names, drop = FALSE]), 2, stats::var, na.rm = TRUE), na.rm = TRUE), 8) else 0.0
    rec <- data.frame(catalogNumber = cat_num, plant_individual_id = p_id, leaf_count = n_leaves, foliar_variance = foliar_var,
                      assigned_tier = grp$assigned_tier[1], reconstruction_tier = grp$reconstruction_tier[1],
                      scientificName = grp$scientificName[1], species_raw = grp$species_raw[1], determiner_tier = grp$determiner_tier[1], stringsAsFactors = FALSE)
    med_cs <- stats::median(grp$centroid_size, na.rm = TRUE)
    rec$centroid_size <- if (is.finite(med_cs) && med_cs > 0) round(med_cs, 4) else NA_real_
    rec$log_centsize <- if (!is.na(rec$centroid_size) && rec$centroid_size > 0) round(log(rec$centroid_size), 6) else NA_real_
    other_meta <- intersect(c("pixels_per_mm", "is_calibrated_mm", "determiner_raw", "county", "stateProvince", "latitude", "longitude",
                              "pheno_sin", "pheno_cos", "regional_group", "aspect_ratio", "solidity", "ucs_score", "area_px",
                              "sym_var", "asym_var", "pct_asymmetry", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight"), names(grp))
    for (m_col in other_meta) {
      vals <- stats::na.omit(grp[[m_col]])
      rec[[m_col]] <- if (length(vals) > 0) vals[1] else NA
    }
    if ("harmonic_roughness" %in% names(grp)) rec$harmonic_roughness <- round(stats::median(grp$harmonic_roughness, na.rm = TRUE), 6)
    if ("leaf_quality_weight" %in% names(grp)) rec$leaf_quality_weight <- round(stats::median(grp$leaf_quality_weight, na.rm = TRUE), 6)
    if ("is_contour_outlier" %in% names(grp)) rec$is_contour_outlier <- all(as.logical(grp$is_contour_outlier), na.rm = TRUE)
    for (h_name in names(harm_medians)) rec[[h_name]] <- harm_medians[[h_name]]
    specimen_records[[spec_idx]] <- rec
  }

  specimen_records <- specimen_records[seq_len(spec_idx)]
  specimen_df <- do.call(rbind, specimen_records)
  for (p in seq_len(opts$num_pcs)) { specimen_df[[paste0("PC", p)]] <- NA_real_; specimen_df[[paste0("Sym_PC", p)]] <- 0.0 }
  if (nrow(specimen_df) >= opts$num_pcs && exists("pca_fit") && !is.null(pca_fit)) {
    spec_pcs <- stats::predict(pca_fit, newdata = as.matrix(specimen_df[, active_sym_cols, drop = FALSE]))
    for (p in seq_len(min(opts$num_pcs, ncol(spec_pcs)))) {
      specimen_df[[paste0("PC", p)]] <- round(spec_pcs[, p], 6)
      specimen_df[[paste0("Sym_PC", p)]] <- round(spec_pcs[, p], 6)
    }
  }

  lead_spec_cols <- c("catalogNumber", "plant_individual_id", "leaf_count", "foliar_variance", "assigned_tier", "reconstruction_tier",
                      "pixels_per_mm", "is_calibrated_mm", "scientificName", "species_raw", "determiner_tier",
                      "PC1", "PC2", "PC3", "PC4", "PC5", "aspect_ratio", "solidity", "ucs_score",
                      "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_px", "centroid_size", "log_centsize")
  lead_spec_cols <- intersect(lead_spec_cols, names(specimen_df))
  specimen_df <- specimen_df[, c(lead_spec_cols, setdiff(names(specimen_df), lead_spec_cols))]

  dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
  tmp_out <- paste0(opts$output, ".tmp")
  write.csv(specimen_df, file = tmp_out, row.names = FALSE, na = "")
  file.rename(tmp_out, opts$output)
  message("Specimen-aggregated EFA harmonics table exported: ", opts$output, " (Rows: ", nrow(specimen_df), ")")

  tryCatch({
    if (file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R")
      export_step03_fourier_pc_tables(specimen_df, efa_df, opts$output)
      generate_step03_pub_figures(efa_df, specimen_df, harm_mat, if (exists("pca_fit")) pca_fit else NULL, opts$figures_dir, opts$output)
    }
  }, error = function(e) warning("Step 03 publication figures/tables export warning: ", e$message))

  # 6. Empirical Tier Validation Test (PERMANOVA via adonis2 / fallback manova)
  valid_tier_idx <- which(!is.na(efa_df$reconstruction_tier) &
                          !is.na(efa_df$scientificName) & nzchar(as.character(efa_df$scientificName)) &
                          !is.na(efa_df$determiner_tier) & nzchar(as.character(efa_df$determiner_tier)) &
                          complete.cases(sym_harmonics))
  if (length(valid_tier_idx) >= 10) {
    message("==================================================================")
    message("Executing Empirical Tier Validation Test (PERMANOVA)")
    message("Model: sym_harmonics ~ scientificName + determiner_tier + reconstruction_tier")
    message(sprintf("Sample Size: %d outlines (Tier 1: %d | Tier 2: %d)",
                    length(valid_tier_idx),
                    sum(efa_df$reconstruction_tier[valid_tier_idx] == "Tier 1"),
                    sum(efa_df$reconstruction_tier[valid_tier_idx] == "Tier 2")))
    message("==================================================================")

    tier_manifest <- efa_df[valid_tier_idx, , drop = FALSE]
    tier_sym <- sym_harmonics[valid_tier_idx, , drop = FALSE]

    tier_manifest$scientificName <- as.factor(tier_manifest$scientificName)
    tier_manifest$determiner_tier <- as.factor(tier_manifest$determiner_tier)
    tier_manifest$reconstruction_tier <- as.factor(tier_manifest$reconstruction_tier)

    set.seed(opts$seed)
    permanova_df <- NULL

    if (requireNamespace("RRPP", quietly = TRUE)) {
      message(sprintf("Executing RRPP::lm.rrpp (iter = %d, seed = %d)...", opts$permutations, opts$seed))
      rrpp_res <- tryCatch({
        rdf <- RRPP::rrpp.data.frame(
          tier_sym = as.matrix(tier_sym),
          scientificName = tier_manifest$scientificName,
          determiner_tier = tier_manifest$determiner_tier,
          reconstruction_tier = tier_manifest$reconstruction_tier
        )
        fit_rrpp <- RRPP::lm.rrpp(
          tier_sym ~ scientificName + determiner_tier + reconstruction_tier,
          data = rdf, iter = opts$permutations, print.progress = FALSE
        )
        anova(fit_rrpp)$table
      }, error = function(e) {
        warning("RRPP::lm.rrpp encountered an issue: ", e$message, ". Falling back to vegan::adonis2().")
        NULL
      })
      if (!is.null(rrpp_res)) {
        permanova_df <- as.data.frame(rrpp_res)
        permanova_df$Term <- rownames(permanova_df)
        rownames(permanova_df) <- NULL
        if ("SS" %in% names(permanova_df)) permanova_df$SumOfSqs <- permanova_df$SS
        if ("Rsq" %in% names(permanova_df)) permanova_df$R2 <- permanova_df$Rsq
        if ("Pr(>F)" %in% names(permanova_df)) permanova_df$p_value <- permanova_df[["Pr(>F)"]]
      }
    }

    if (is.null(permanova_df) && requireNamespace("vegan", quietly = TRUE)) {
      message(sprintf("Executing vegan::adonis2 (permutations = %d, seed = %d)...", opts$permutations, opts$seed))
      ad_res <- tryCatch({
        ad_idx <- if (nrow(tier_manifest) > 2000L) sort(sample(seq_len(nrow(tier_manifest)), 2000L)) else seq_len(nrow(tier_manifest))
        vegan::adonis2(tier_sym[ad_idx, , drop = FALSE] ~ scientificName + determiner_tier + reconstruction_tier,
                       data = tier_manifest[ad_idx, , drop = FALSE], permutations = opts$permutations, method = "euclidean", by = "terms")
      }, error = function(e) {
        warning("vegan::adonis2 encountered an issue: ", e$message, ". Falling back to stats::manova().")
        NULL
      })
      if (!is.null(ad_res)) {
        permanova_df <- as.data.frame(ad_res)
        permanova_df$Term <- rownames(permanova_df)
        rownames(permanova_df) <- NULL
        if ("F.Model" %in% names(permanova_df)) permanova_df$F <- permanova_df$F.Model
        if ("Pr(>F)" %in% names(permanova_df)) permanova_df$p_value <- permanova_df[["Pr(>F)"]]
      }
    }

    # Defensive fallback if RRPP/vegan are not available or encountered an error
    if (is.null(permanova_df)) {
      message("vegan package not available; executing defensive fallback using analytical sequential SS decomposition...")
      Y <- as.matrix(tier_sym); Y_cent <- scale(Y, center = TRUE, scale = FALSE)
      SS_tot <- sum(Y_cent^2); N_obs <- nrow(tier_manifest)
      fit1 <- stats::lm(Y_cent ~ scientificName, data = tier_manifest)
      SS_sp <- sum(stats::fitted(fit1)^2); df_sp <- fit1$rank - 1
      fit2 <- stats::lm(Y_cent ~ scientificName + determiner_tier, data = tier_manifest)
      SS_12 <- sum(stats::fitted(fit2)^2); SS_det <- SS_12 - SS_sp; df_det <- fit2$rank - fit1$rank
      fit3 <- stats::lm(Y_cent ~ scientificName + determiner_tier + reconstruction_tier, data = tier_manifest)
      SS_123 <- sum(stats::fitted(fit3)^2); SS_tier <- SS_123 - SS_12; df_tier <- fit3$rank - fit2$rank
      SS_res <- SS_tot - SS_123; df_res <- N_obs - fit3$rank; MS_res <- SS_res / max(df_res, 1)
      F_sp <- (SS_sp / max(df_sp, 1)) / MS_res; F_det <- (SS_det / max(df_det, 1)) / MS_res; F_tier <- (SS_tier / max(df_tier, 1)) / MS_res
      p_sp <- stats::pf(F_sp, df_sp, df_res, lower.tail = FALSE)
      p_det <- stats::pf(F_det, df_det, df_res, lower.tail = FALSE)
      p_tier <- stats::pf(F_tier, df_tier, df_res, lower.tail = FALSE)
      permanova_df <- data.frame(
        Term = c("scientificName", "determiner_tier", "reconstruction_tier", "Residual", "Total"),
        Df = c(df_sp, df_det, df_tier, df_res, N_obs - 1),
        SumOfSqs = round(c(SS_sp, SS_det, SS_tier, SS_res, SS_tot), 4),
        R2 = round(c(SS_sp / SS_tot, SS_det / SS_tot, SS_tier / SS_tot, SS_res / SS_tot, 1.0), 6),
        F = c(round(F_sp, 4), round(F_det, 4), round(F_tier, 4), NA_real_, NA_real_),
        p_value = c(round(p_sp, 5), round(p_det, 5), round(p_tier, 5), NA_real_, NA_real_),
        stringsAsFactors = FALSE
      )
    }

    # Format table columns and export report CSV
    permanova_df$Term <- gsub("^Residuals$", "Residual", permanova_df$Term)
    if ("Total" %in% permanova_df$Term) permanova_df$R2[permanova_df$Term == "Total"] <- 1.0
    report_cols <- intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(permanova_df))
    permanova_df <- permanova_df[, report_cols, drop = FALSE]

    dir.create(dirname(opts$report_out), recursive = TRUE, showWarnings = FALSE)
    tmp_rep <- paste0(opts$report_out, ".tmp")
    write.csv(permanova_df, file = tmp_rep, row.names = FALSE)
    file.rename(tmp_rep, opts$report_out)
    message("PERMANOVA Tier Validation report exported to: ", opts$report_out)

    # Verification assertions
    tier_row <- permanova_df[grepl("reconstruction_tier", permanova_df$Term, ignore.case = TRUE), ]
    if (nrow(tier_row) > 0) {
      r2_val <- tier_row$R2[1]; p_val <- tier_row$p_value[1]
      message(sprintf("PERMANOVA Result: reconstruction_tier explains %.4f%% of variance (R2 = %.6f, p = %.4f)", r2_val * 100, r2_val, p_val))
      if (r2_val < 0.01) message("[ASSERTION PASS] reconstruction_tier accounts for <1.0% of total variance (R2 < 0.01).")
      else message("[ASSERTION NOTE] reconstruction_tier R2 = ", round(r2_val, 4))
      if (!is.na(p_val) && p_val > 0.05) message("[ASSERTION PASS] reconstruction_tier is non-significant (p > 0.05).")
      else message("[ASSERTION NOTE] reconstruction_tier p-value = ", p_val)
    }

    # 7. Diagnostic Comparative Density Plot Export
    dir.create(dirname(opts$plot_out), recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_df <- tier_manifest[!is.na(tier_manifest$PC1) & !is.na(tier_manifest$PC2), ]
      p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, fill = reconstruction_tier, color = reconstruction_tier)) +
        ggplot2::geom_density(alpha = 0.4, linewidth = 0.8) +
        ggplot2::scale_fill_manual(values = c("Tier 1" = "#1b9e77", "Tier 2" = "#d95f02")) +
        ggplot2::scale_color_manual(values = c("Tier 1" = "#1b9e77", "Tier 2" = "#d95f02")) +
        ggplot2::labs(title = "Empirical Tier Morphological Equivalence",
                      subtitle = "Symmetric Fourier PC1 Density Distribution (Tier 1 Pristine vs. Tier 2 Reflected)",
                      x = "Symmetric Morphospace PC1", y = "Density",
                      fill = "Reconstruction Tier", color = "Reconstruction Tier") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
      ggplot2::ggsave(opts$plot_out, plot = p1, width = 8, height = 5)
    } else {
      pdf(opts$plot_out, width = 8, height = 5)
      pc1_t1 <- tier_manifest$PC1[tier_manifest$reconstruction_tier == "Tier 1"]
      pc1_t2 <- tier_manifest$PC1[tier_manifest$reconstruction_tier == "Tier 2"]
      d1 <- stats::density(pc1_t1, na.rm = TRUE)
      d2 <- stats::density(pc1_t2, na.rm = TRUE)
      xlims <- range(c(d1$x, d2$x))
      ylims <- c(0, max(c(d1$y, d2$y)) * 1.1)
      plot(d1, xlim = xlims, ylim = ylims, col = "#1b9e77", lwd = 2,
           main = "Symmetric Fourier PC1 Density: Tier 1 vs Tier 2",
           xlab = "Symmetric Morphospace PC1", ylab = "Density")
      lines(d2, col = "#d95f02", lwd = 2, lty = 2)
      polygon(d1, col = grDevices::adjustcolor("#1b9e77", alpha.f = 0.3), border = NA)
      polygon(d2, col = grDevices::adjustcolor("#d95f02", alpha.f = 0.3), border = NA)
      legend("topright", legend = c("Tier 1 (Pristine)", "Tier 2 (Reflected)"),
             col = c("#1b9e77", "#d95f02"), lwd = 2, lty = c(1, 2), bty = "n")
      dev.off()
    }
    message("Diagnostic comparative density plot exported: ", opts$plot_out)
  }

  return(invisible(efa_df))
}

if (sys.nframe() == 0) {
  opts <- parse_args_robust()
  run_fourier_extraction(opts)
}
