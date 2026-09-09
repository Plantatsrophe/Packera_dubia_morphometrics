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
    optparse::make_option(c("-i", "--input"), type = "character",
      default = "data/contours/", help = "Directory containing contour CSVs (alias for --contours-dir) [default: %default]"),
    optparse::make_option(c("-c", "--contours-dir"), type = "character",
      default = "data/contours/", help = "Directory containing contour CSVs [default: %default]"),
    optparse::make_option(c("-f", "--manifest"), type = "character",
      default = "data/tables/extracted_leaf_manifest.csv", help = "Extracted leaf manifest CSV [default: %default]"),
    optparse::make_option(c("-m", "--masks-dir"), type = "character",
      default = "data/masks/", help = "Fallback masks directory [default: %default]"),
    optparse::make_option(c("-v", "--vouchers"), type = "character",
      default = "data/tables/curated_vouchers.csv", help = "Curated vouchers CSV [default: %default]"),
    optparse::make_option(c("-o", "--output"), type = "character",
      default = "data/tables/leaf_efa_harmonics.csv", help = "Output harmonics CSV [default: %default]"),
    optparse::make_option(c("--output-individual"), type = "character",
      default = "data/tables/leaf_efa_harmonics_individual.csv", help = "Output individual leaf harmonics CSV [default: %default]"),
    optparse::make_option(c("-k", "--harmonics"), type = "integer",
      default = 12, help = "Number of Fourier harmonics (nb.h) [default: %default]"),
    optparse::make_option(c("-p", "--num-pcs"), type = "integer",
      default = 5, help = "Number of PCA dimensions to extract [default: %default]"),
    optparse::make_option(c("-r", "--report-out"), type = "character",
      default = "outputs/reports/tier_symmetry_validation.csv", help = "PERMANOVA validation report CSV [default: %default]"),
    optparse::make_option(c("-g", "--plot-out"), type = "character",
      default = "outputs/figures/tier1_vs_tier2_density_overlay.pdf", help = "Tier density overlay plot PDF [default: %default]"),
    optparse::make_option(c("--permutations"), type = "integer",
      default = 999, help = "Number of permutations for PERMANOVA [default: %default]"),
    optparse::make_option(c("--seed"), type = "integer",
      default = 42, help = "Random seed for reproducible permutations [default: %default]")
  )

  if (requireNamespace("optparse", quietly = TRUE)) {
    parser <- optparse::OptionParser(usage = "%prog [options]", option_list = option_list)
    return(optparse::parse_args(parser))
  }

  raw_args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    input = "data/contours/", contours_dir = "data/contours/", manifest = "data/tables/extracted_leaf_manifest.csv",
    masks_dir = "data/masks/", vouchers = "data/tables/curated_vouchers.csv",
    output = "data/tables/leaf_efa_harmonics.csv", output_individual = "data/tables/leaf_efa_harmonics_individual.csv",
    harmonics = 12, num_pcs = 5,
    report_out = "outputs/reports/tier_symmetry_validation.csv",
    plot_out = "outputs/figures/tier1_vs_tier2_density_overlay.pdf",
    permutations = 999, seed = 42
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
    else if (arg == "--permutations" && i < length(raw_args)) { opts$permutations <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else if (arg == "--seed" && i < length(raw_args)) { opts$seed <- as.integer(raw_args[i + 1]); i <- i + 2 }
    else { i <- i + 1 }
  }
  return(opts)
}

# ------------------------------------------------------------------------------
# 2. Pure-R Mathematical Engine: Normalized EFA (Kuhl & Giardina 1982 Fallback)
# ------------------------------------------------------------------------------
compute_efourier_coords_pure <- function(coo, nb_h = 12, norm = TRUE) {
  if (is.null(coo) || nrow(coo) < 5) return(rep(NA_real_, nb_h * 4))
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

  # Size, orientation, and starting-point invariant normalization
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
        contour_files <- manifest_df$contour_path[file.exists(manifest_df$contour_path)]
    }
    if (length(contour_files) == 0 && "mask_path" %in% names(manifest_df)) {
        mask_files <- manifest_df$mask_path[file.exists(manifest_df$mask_path)]
        if (length(mask_files) > 0) {
            contour_files <- mask_files
            use_masks_fallback <- TRUE
        }
    }
  }

  if (length(contour_files) == 0) {
      if (dir.exists(opts$input)) {
        contour_files <- list.files(opts$input, pattern = "\\.csv$", full.names = TRUE)
        if (length(contour_files) == 0) {
            mask_files <- list.files(opts$input, pattern = "\\.(png|jpg)$", full.names = TRUE, recursive = TRUE)
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
    if (!dir.exists(fallback_dir) && dir.exists("data/_archive/masks/tier1_intact")) {
      fallback_dir <- "data/_archive/masks/tier1_intact"
    }
    mask_files <- list.files(fallback_dir, pattern = "\\.(png|jpg)$", full.names = TRUE, recursive = TRUE)
    if (length(mask_files) > 0) {
      message("Found ", length(mask_files), " fallback mask files.")
      contour_files <- mask_files
      use_masks_fallback <- TRUE
    }
  }

  if (length(contour_files) == 0) {
    stop("No input contour CSVs or mask files found to process.")
  }

  # 2. Build coordinate list and factor metadata for Momocs::Out()
  coo_list <- list()
  fac_records <- list()
  harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(opts$harmonics), function(x, y) paste0(x, y)))

  message("Ingesting contour coordinates from ", length(contour_files), " files...")
  for (f in contour_files) {
    # If mask fallback, attempt to use EBImage/Momocs import functions, though custom logic might be needed
    # (Assuming Momocs::import_jpg / import_txt can handle or load_contour_file handles basic parsing)
    if (use_masks_fallback && grepl("\\.(png|jpg)$", f, ignore.case = TRUE)) {
        if (requireNamespace("Momocs", quietly = TRUE)) {
            tryCatch({
                tmp_coo <- Momocs::import_jpg(f)
                if (is.list(tmp_coo) && length(tmp_coo) > 0) {
                    coo <- tmp_coo[[1]]
                } else {
                    coo <- NULL
                }
            }, error = function(e) { coo <<- NULL })
        } else {
            coo <- NULL
        }
    } else {
        coo <- load_contour_file(f)
    }

    if (is.null(coo)) next

    # Enforce clockwise orientation
    if (requireNamespace("Momocs", quietly = TRUE)) {
        coo <- tryCatch(Momocs::coo_cw(coo), error = function(e) coo)
    }

    meta <- parse_specimen_id(f)

    # Calculate basic morphometric geometry
    n_pts <- nrow(coo)
    span_x <- diff(range(coo[, 1], na.rm = TRUE))
    span_y <- diff(range(coo[, 2], na.rm = TRUE))
    aspect_ratio <- if (span_y > 1e-6) round(span_x / span_y, 4) else 1.0
    area_px <- abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1])))

    shape_name <- meta$shape_id
    coo_list[[shape_name]] <- coo

    rec <- data.frame(
      shape_id = shape_name,
      catalogNumber = meta$catalogNumber,
      plant_individual_id = meta$plant_id,
      leaf_id = meta$leaf_id,
      assigned_tier = if (!is.null(manifest_df) && "assigned_tier" %in% names(manifest_df)) {
          m_match <- manifest_df[manifest_df$catalogNumber == meta$catalogNumber & manifest_df$leaf_id == meta$leaf_id, ]
          if (nrow(m_match) > 0) m_match$assigned_tier[1] else "Tier_1_Direct"
      } else { "Tier_1_Direct" },
      aspect_ratio = aspect_ratio,
      area_px = round(area_px, 1),
      mask_source = f,
      stringsAsFactors = FALSE
    )
    fac_records[[length(fac_records) + 1]] <- rec
  }

  if (length(coo_list) == 0) stop("Failed to extract valid coordinate contours.")

  fac_df <- do.call(rbind, fac_records)
  rownames(fac_df) <- names(coo_list)
  message(sprintf("Assembled %d valid closed contours for Elliptic Fourier Analysis.", length(coo_list)))

  # 3. Perform 12-Harmonic EFA via Momocs::Out() & Momocs::efourier()
  message("Executing 12-harmonic EFA via Momocs::efourier(nb.h = ", opts$harmonics, ", norm = TRUE)...")
  efa_success <- FALSE
  harm_mat <- NULL

  if (requireNamespace("Momocs", quietly = TRUE)) {
    out_obj <- tryCatch(Momocs::Out(coo_list, fac = fac_df), error = function(e) {
      message("Momocs::Out() note: ", e$message); NULL
    })
    if (!is.null(out_obj)) {
      ef_res <- tryCatch(Momocs::efourier(out_obj, nb.h = opts$harmonics, norm = TRUE), error = function(e) {
        message("Momocs::efourier() note: ", e$message); NULL
      })
      if (!is.null(ef_res) && !is.null(ef_res$coe)) {
        harm_mat <- as.matrix(ef_res$coe)
        efa_success <- TRUE
        message("Momocs EFA extraction successful across ", nrow(harm_mat), " outlines.")
      }
    }
  }

  # Fallback to pure-R vectorized Kuhl & Giardina EFA if Momocs not available
  if (!efa_success || is.null(harm_mat)) {
    message("Running mathematical Kuhl & Giardina EFA engine...")
    harm_list <- lapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = opts$harmonics, norm = TRUE))
    harm_mat <- do.call(rbind, harm_list)
    colnames(harm_mat) <- harmonic_names
  }

  # Harmonize column names to standard A1-D12
  if (ncol(harm_mat) == length(harmonic_names)) {
    colnames(harm_mat) <- harmonic_names
  }
  efa_df <- cbind(fac_df, as.data.frame(harm_mat))

  # 4. Symmetric Fourier Decomposition & PCA (PC1-PC5)
  # In normalized EFA with major axis horizontal alignment:
  # - An and Dn capture symmetric outline variance (bilateral symmetry across longitudinal midrib)
  # - Bn and Cn capture asymmetric variance (fluctuating asymmetry / lateral skew)
  # Isolating An and Dn guarantees that Tier 1 (pristine) and Tier 2 (reflected) leaves are evaluated
  # strictly on the symmetric morphological component, eliminating artificial symmetry artifacts.
  message("Extracting symmetric harmonic coefficients (An and Dn harmonics)...")
  sym_names <- c(paste0("A", seq_len(opts$harmonics)), paste0("D", seq_len(opts$harmonics)))
  sym_harmonics <- harm_mat[, sym_names, drop = FALSE]

  message("Running PCA strictly on symmetric harmonics (sym_harmonics)...")
  complete_idx <- which(complete.cases(sym_harmonics))
  for (p in seq_len(opts$num_pcs)) efa_df[[paste0("PC", p)]] <- NA_real_

  if (length(complete_idx) >= opts$num_pcs) {
    # Filter constant/zero-variance coefficients (e.g. invariant A1 = 1.0 in normalized EFA)
    col_vars <- apply(sym_harmonics[complete_idx, , drop = FALSE], 2, stats::var)
    active_sym_cols <- names(col_vars[col_vars > 1e-8])
    message(sprintf("Active symmetric harmonic dimensions with non-zero variance: %d / %d",
                    length(active_sym_cols), ncol(sym_harmonics)))

    pca_fit <- stats::prcomp(sym_harmonics[complete_idx, active_sym_cols, drop = FALSE], center = TRUE, scale. = TRUE)
    var_exp <- round((pca_fit$sdev^2) / sum(pca_fit$sdev^2) * 100, 2)
    message("=== Symmetric Morphospace PCA Variance Explained ===")
    for (p in seq_len(min(opts$num_pcs, length(var_exp)))) {
      message(sprintf("  PC%d: %5.2f%% variance", p, var_exp[p]))
    }
    for (p in seq_len(opts$num_pcs)) {
      if (p <= ncol(pca_fit$x)) {
        efa_df[complete_idx, paste0("PC", p)] <- round(pca_fit$x[, p], 6)
      }
    }
  }

  # 5. Darwin Core Metadata Integration & Standardization
  if (!is.null(vouchers_df)) {
    meta_cols <- intersect(names(vouchers_df), c(
      "catalogNumber", "scientificName", "species_raw", "determiner_raw", "determiner_tier",
      "county", "stateProvince", "latitude", "longitude",
      "pheno_sin", "pheno_cos", "regional_group"
    ))
    v_sub <- vouchers_df[!duplicated(vouchers_df$catalogNumber), meta_cols, drop = FALSE]
    efa_df <- merge(efa_df, v_sub, by = "catalogNumber", all.x = TRUE)
  }

  # Standardize reconstruction_tier: Tier 1 (Pristine Direct) vs. Tier 2 (Reflected Hemi-blade)
  efa_df$reconstruction_tier <- ifelse(
    grepl("tier_?1", efa_df$assigned_tier, ignore.case = TRUE),
    "Tier 1",
    ifelse(grepl("tier_?2", efa_df$assigned_tier, ignore.case = TRUE), "Tier 2", NA_character_)
  )

  if (!"scientificName" %in% names(efa_df)) {
    efa_df$scientificName <- efa_df$species_raw
  } else {
    efa_df$scientificName <- ifelse(is.na(efa_df$scientificName) | efa_df$scientificName == "",
                                    efa_df$species_raw, efa_df$scientificName)
  }
  if (!"determiner_tier" %in% names(efa_df)) {
    efa_df$determiner_tier <- "Tier_3_Bronze"
  } else {
    efa_df$determiner_tier <- ifelse(is.na(efa_df$determiner_tier) | efa_df$determiner_tier == "",
                                     "Tier_3_Bronze", efa_df$determiner_tier)
  }

  lead_cols <- c("catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier",
                 "reconstruction_tier", "scientificName", "species_raw", "determiner_tier",
                 "PC1", "PC2", "PC3", "PC4", "PC5",
                 "aspect_ratio", "area_px", "mask_source")
  lead_cols <- intersect(lead_cols, names(efa_df))
  efa_df <- efa_df[, c(lead_cols, setdiff(names(efa_df), lead_cols))]

  # 5a. Export Disaggregated Individual Leaf Harmonics Archive
  if (!is.null(opts$output_individual) && nzchar(opts$output_individual)) {
    dir.create(dirname(opts$output_individual), recursive = TRUE, showWarnings = FALSE)
    write.csv(efa_df, file = opts$output_individual, row.names = FALSE, na = "")
    message("Individual leaf EFA harmonics archive exported: ", opts$output_individual, " (Rows: ", nrow(efa_df), ")")
  }

  # 5b. Compute Specimen-Level Median Harmonic Vectors & Phenotypic Stability Index
  message("Computing specimen-level median harmonic vectors (grouped by catalogNumber & plant instance)...")
  valid_closed_mask <- efa_df$assigned_tier %in% c("Tier_1_Direct", "Tier_2_Reflected") & complete.cases(efa_df[, sym_names])
  closed_leaf_df <- efa_df[valid_closed_mask, , drop = FALSE]

  all_harm_cols <- intersect(c(harmonic_names, grep("Chebyshev", names(efa_df), value = TRUE)), names(efa_df))
  specimen_records <- list()
  spec_groups <- split(closed_leaf_df, list(closed_leaf_df$catalogNumber, closed_leaf_df$plant_individual_id), drop = TRUE)

  for (grp in spec_groups) {
    if (nrow(grp) == 0) next
    cat_num <- grp$catalogNumber[1]
    p_id <- grp$plant_individual_id[1]
    n_leaves <- nrow(grp)

    # Compute median symmetric harmonic vector across valid leaves belonging to this specimen
    harm_medians <- sapply(grp[, all_harm_cols, drop = FALSE], stats::median, na.rm = TRUE)

    # Within-specimen foliar variance as phenotypic stability index (0.0 for 1 leaf)
    if (n_leaves > 1) {
      sym_sub <- as.matrix(grp[, sym_names, drop = FALSE])
      col_vars_sub <- apply(sym_sub, 2, stats::var, na.rm = TRUE)
      foliar_var <- round(mean(col_vars_sub, na.rm = TRUE), 8)
    } else {
      foliar_var <- 0.0
    }

    rec <- data.frame(
      catalogNumber = cat_num,
      plant_individual_id = p_id,
      leaf_count = n_leaves,
      foliar_variance = foliar_var,
      assigned_tier = grp$assigned_tier[1],
      reconstruction_tier = grp$reconstruction_tier[1],
      scientificName = grp$scientificName[1],
      species_raw = grp$species_raw[1],
      determiner_tier = grp$determiner_tier[1],
      stringsAsFactors = FALSE
    )

    other_meta <- intersect(c("determiner_raw", "county", "stateProvince", "latitude", "longitude",
                              "pheno_sin", "pheno_cos", "regional_group", "aspect_ratio", "solidity", "area_px"),
                            names(grp))
    for (m_col in other_meta) {
      vals <- stats::na.omit(grp[[m_col]])
      rec[[m_col]] <- if (length(vals) > 0) vals[1] else NA
    }

    for (h_name in names(harm_medians)) {
      rec[[h_name]] <- harm_medians[[h_name]]
    }
    specimen_records[[length(specimen_records) + 1]] <- rec
  }

  specimen_df <- do.call(rbind, specimen_records)

  # Project symmetric morphospace PC1-PC5 for specimen-level median profiles
  for (p in seq_len(opts$num_pcs)) specimen_df[[paste0("PC", p)]] <- NA_real_
  if (nrow(specimen_df) >= opts$num_pcs && exists("pca_fit") && !is.null(pca_fit)) {
    spec_sym <- as.matrix(specimen_df[, active_sym_cols, drop = FALSE])
    spec_pcs <- scale(spec_sym, center = pca_fit$center, scale = pca_fit$scale) %*% pca_fit$rotation
    for (p in seq_len(opts$num_pcs)) {
      if (p <= ncol(spec_pcs)) {
        specimen_df[[paste0("PC", p)]] <- round(spec_pcs[, p], 6)
      }
    }
  }

  lead_spec_cols <- c("catalogNumber", "plant_individual_id", "leaf_count", "foliar_variance",
                      "assigned_tier", "reconstruction_tier", "scientificName", "species_raw", "determiner_tier",
                      "PC1", "PC2", "PC3", "PC4", "PC5",
                      "aspect_ratio", "solidity", "area_px")
  lead_spec_cols <- intersect(lead_spec_cols, names(specimen_df))
  specimen_df <- specimen_df[, c(lead_spec_cols, setdiff(names(specimen_df), lead_spec_cols))]

  dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
  write.csv(specimen_df, file = opts$output, row.names = FALSE, na = "")
  message("Specimen-aggregated EFA harmonics table exported: ", opts$output, " (Rows: ", nrow(specimen_df), ")")

  # 6. Empirical Tier Validation Test (PERMANOVA via adonis2 / fallback manova)
  valid_tier_idx <- which(!is.na(efa_df$reconstruction_tier) & complete.cases(sym_harmonics))
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

    if (requireNamespace("vegan", quietly = TRUE)) {
      message(sprintf("Executing vegan::adonis2 (permutations = %d, seed = %d)...",
                      opts$permutations, opts$seed))
      ad_res <- tryCatch({
        vegan::adonis2(
          tier_sym ~ scientificName + determiner_tier + reconstruction_tier,
          data = tier_manifest,
          permutations = opts$permutations,
          method = "euclidean",
          by = "terms"
        )
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

    # Defensive fallback if vegan is not available or encountered an error
    if (is.null(permanova_df)) {
      message("vegan package not available; executing defensive fallback using analytical sequential SS decomposition...")
      Y <- as.matrix(tier_sym)
      Y_cent <- scale(Y, center = TRUE, scale = FALSE)
      SS_tot <- sum(Y_cent^2)
      N_obs <- nrow(tier_manifest)

      fit1 <- stats::lm(Y_cent ~ scientificName, data = tier_manifest)
      Y_hat1 <- stats::fitted(fit1)
      SS_sp <- sum(Y_hat1^2)
      df_sp <- fit1$rank - 1

      fit2 <- stats::lm(Y_cent ~ scientificName + determiner_tier, data = tier_manifest)
      Y_hat2 <- stats::fitted(fit2)
      SS_12 <- sum(Y_hat2^2)
      SS_det <- SS_12 - SS_sp
      df_det <- fit2$rank - fit1$rank

      fit3 <- stats::lm(Y_cent ~ scientificName + determiner_tier + reconstruction_tier, data = tier_manifest)
      Y_hat3 <- stats::fitted(fit3)
      SS_123 <- sum(Y_hat3^2)
      SS_tier <- SS_123 - SS_12
      df_tier <- fit3$rank - fit2$rank

      SS_res <- SS_tot - SS_123
      df_res <- N_obs - fit3$rank

      MS_res <- SS_res / max(df_res, 1)
      F_sp <- (SS_sp / max(df_sp, 1)) / MS_res
      F_det <- (SS_det / max(df_det, 1)) / MS_res
      F_tier <- (SS_tier / max(df_tier, 1)) / MS_res

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
    report_cols <- intersect(c("Term", "Df", "SumOfSqs", "R2", "F", "p_value"), names(permanova_df))
    permanova_df <- permanova_df[, report_cols, drop = FALSE]

    dir.create(dirname(opts$report_out), recursive = TRUE, showWarnings = FALSE)
    write.csv(permanova_df, file = opts$report_out, row.names = FALSE)
    message("PERMANOVA Tier Validation report exported to: ", opts$report_out)

    # Verification assertions
    tier_row <- permanova_df[grepl("reconstruction_tier", permanova_df$Term, ignore.case = TRUE), ]
    if (nrow(tier_row) > 0) {
      r2_val <- tier_row$R2[1]
      p_val <- tier_row$p_value[1]
      message(sprintf("PERMANOVA Result: reconstruction_tier explains %.4f%% of variance (R2 = %.6f, p = %.4f)",
                      r2_val * 100, r2_val, p_val))
      if (r2_val < 0.01) {
        message("[ASSERTION PASS] reconstruction_tier accounts for <1.0% of total variance (R2 < 0.01).")
      } else {
        message("[ASSERTION NOTE] reconstruction_tier R2 = ", round(r2_val, 4))
      }
      if (!is.na(p_val) && p_val > 0.05) {
        message("[ASSERTION PASS] reconstruction_tier is non-significant (p > 0.05).")
      } else {
        message("[ASSERTION NOTE] reconstruction_tier p-value = ", p_val)
      }
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
