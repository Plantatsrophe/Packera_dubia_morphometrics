#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/08_discover_taxon_complexes.R
# Project: Packera Species Delimitation & Morphometrics Pipeline (NCU)
# Description:
#   Unsupervised Quantitative Taxon-Complex Discovery & Integrative Delimitation
#   Engine (Phase 08): (1) Reliability-Weighted Taxon Affinity Matrix (W_ij) with
#   out-of-fold Cohen's kappa weights (w_m = max(0.05, kappa_m)) and missing-
#   modality renormalization; (2) Hopkins-Corrected Circular Anthesis Phenology
#   with BH FDR; (3) 5-Stage Hierarchical Integrative Decision Matrix (Stage 1
#   Misidentification_Reassign -> Stage 2 Hybrid_Introgressant -> Stage 3
#   Putative_Cryptic_Species -> Stage 4 Subspecies_Ecotype/Plasticity -> Stage 5
#   Confirmed_Species); and (4) Leiden/Louvain & Ward.D2 clustering across taxa
#   or Taxon x Morphotype subclusters ('P. <epithet> [K1]', 'P. <epithet> [K2]').
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
  if (requireNamespace("igraph", quietly = TRUE)) library(igraph)
  if (requireNamespace("cluster", quietly = TRUE)) library(cluster)
  if (requireNamespace("MASS", quietly = TRUE)) library(MASS)
})

if (file.exists("scripts/morphometrics/pub_figures_utils.R")) {
  source("scripts/morphometrics/pub_figures_utils.R")
}
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
}

# ------------------------------------------------------------------------------
# 1. CLI Argument Parser & Atomic Export Helpers
# ------------------------------------------------------------------------------
parse_discovery_args <- function() {
  run_dir <- Sys.getenv("PACKERA_RUN_DIR", unset = "")
  rp <- function(rel) if (nzchar(run_dir)) file.path(run_dir, rel) else rel

  defaults <- list(
    morphometrics    = rp("data/tables/morphometric_flags.csv"),
    harmonics        = rp("data/tables/leaf_efa_harmonics.csv"),
    vision_audit     = rp("data/tables/label_noise_audit.csv"),
    vouchers         = rp("data/tables/curated_vouchers.csv"),
    multimodal       = rp("data/tables/multimodal_conflict_flags.csv"),
    focal_taxon      = "",
    min_affinity     = 0.03,
    output_complexes = rp("outputs/reports/empirical_taxon_complexes.csv"),
    output_affinity  = rp("outputs/reports/pairwise_taxon_affinity_matrix.csv"),
    output_neighbors = rp("outputs/reports/focal_auto_neighbors.txt"),
    output_figure    = rp("outputs/figures/fig_empirical_taxon_complex_network.pdf")
  )

  raw_args <- commandArgs(trailingOnly = TRUE)
  if (length(raw_args) == 0) return(defaults)

  opts <- defaults
  i <- 1L
  flag_map <- c(
    "-m" = "morphometrics", "--morphometrics" = "morphometrics",
    "-e" = "harmonics", "--harmonics" = "harmonics",
    "-n" = "vision_audit", "--vision-audit" = "vision_audit",
    "-v" = "vouchers", "--vouchers" = "vouchers",
    "-c" = "multimodal", "--multimodal" = "multimodal",
    "--focal-taxon" = "focal_taxon",
    "--output-complexes" = "output_complexes",
    "--output-affinity" = "output_affinity",
    "--output-neighbors" = "output_neighbors",
    "--output-figure" = "output_figure"
  )
  while (i <= length(raw_args)) {
    a <- raw_args[i]
    if (a %in% names(flag_map) && i < length(raw_args)) {
      opts[[flag_map[[a]]]] <- raw_args[i + 1L]; i <- i + 2L
    } else if (a == "--min-affinity" && i < length(raw_args)) {
      opts$min_affinity <- as.numeric(raw_args[i + 1L]); i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  opts
}

atomic_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  utils::write.csv(df, file = tmp, row.names = FALSE)
  file.rename(tmp, path)
  invisible(path)
}

atomic_write_lines <- function(lines_vec, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  writeLines(as.character(lines_vec), con = tmp)
  file.rename(tmp, path)
  invisible(path)
}

std_taxon <- function(x) {
  if (exists("standardize_taxon_name", mode = "function")) {
    res <- standardize_taxon_name(x)
    return(ifelse(is.na(res), as.character(x), res))
  }
  as.character(x)
}

# ------------------------------------------------------------------------------
# 2. Step A: Load & Filter Tier_1_Pure Gold/Silver Anchor Specimens
# ------------------------------------------------------------------------------
load_and_merge_inputs <- function(opts) {
  load_step08_merged_inputs(opts, std_taxon)
}

select_taxon_anchors <- function(morph_df) {
  is_t1_pure <- rep(TRUE, nrow(morph_df))
  if ("n_tier1_leaves" %in% names(morph_df)) {
    n_t1 <- suppressWarnings(as.numeric(morph_df$n_tier1_leaves))
    if (any(!is.na(n_t1))) is_t1_pure <- !is.na(n_t1) & (n_t1 >= 1)
  } else if ("morph_tier_source" %in% names(morph_df)) {
    mts <- as.character(morph_df$morph_tier_source)
    if (any(!is.na(mts) & nzchar(mts))) is_t1_pure <- !is.na(mts) & (mts == "Tier_1_Pure")
  }
  det_tier <- if ("determiner_tier" %in% names(morph_df)) as.character(morph_df$determiner_tier) else rep("Tier_1_Gold", nrow(morph_df))
  is_gs <- det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")

  taxa <- sort(unique(morph_df$species_standardized))
  anchor_mask <- rep(FALSE, nrow(morph_df))
  n_gs_vec <- stats::setNames(integer(length(taxa)), taxa)
  for (tx in taxa) {
    tx_idx <- which(morph_df$species_standardized == tx)
    gs_pure_idx <- tx_idx[is_t1_pure[tx_idx] & is_gs[tx_idx]]
    n_gs_vec[tx] <- length(gs_pure_idx)
    if (length(gs_pure_idx) >= 5L) anchor_mask[gs_pure_idx] <- TRUE
    else if (length(tx_idx[is_t1_pure[tx_idx]]) >= 3L) anchor_mask[tx_idx[is_t1_pure[tx_idx]]] <- TRUE
    else anchor_mask[tx_idx] <- TRUE
  }
  list(anchor_df = morph_df[anchor_mask, , drop = FALSE], anchor_mask = anchor_mask, n_gold_silver = n_gs_vec)
}

# ------------------------------------------------------------------------------
# 3. Step B: Empirical Cohen's Kappa Reliability Weights (Dawid & Skene 1979)
# ------------------------------------------------------------------------------
compute_cohens_kappa_stat <- function(y_true, y_pred) {
  valid <- !is.na(y_true) & nzchar(y_true) & !is.na(y_pred) & nzchar(y_pred) &
           !(y_pred %in% c("Unknown", "Unresolved", "NA"))
  if (sum(valid) < 2L) return(NA_real_)
  yt <- as.character(y_true[valid])
  yp <- as.character(y_pred[valid])
  lvls <- sort(unique(c(yt, yp)))
  if (length(lvls) < 2L || length(unique(yt)) < 2L) return(max(0.05, min(0.95, mean(yt == yp))))
  tab <- table(factor(yt, levels = lvls), factor(yp, levels = lvls))
  n_tot <- sum(tab)
  p_o <- sum(diag(tab)) / n_tot
  p_e <- sum((rowSums(tab) / n_tot) * (colSums(tab) / n_tot))
  if (p_e >= 1.0 - 1e-12) return(1.0)
  (p_o - p_e) / (1.0 - p_e)
}

compute_modality_kappa_weights <- function(morph_df, vis_df, anchor_mask) {
  # Evaluates out-of-fold Cohen's kappa (w_m = max(0.05, kappa_m)) on
  # Tier_1_Gold + Tier_2_Silver anchors (Dawid & Skene 1979; Kuncheva 2014).
  anch_df <- morph_df[anchor_mask, , drop = FALSE]
  kappa_morph <- if ("cda_predicted_taxon" %in% names(anch_df)) {
    compute_cohens_kappa_stat(anch_df$species_standardized, anch_df$cda_predicted_taxon)
  } else {
    NA_real_
  }
  w_morph <- if (is.finite(kappa_morph)) max(0.05, min(1.0, kappa_morph)) else 0.70

  kappa_vis <- NA_real_
  has_vis <- !is.null(vis_df) && nrow(vis_df) > 0L && "vision_predicted_label" %in% names(vis_df)
  if (has_vis) {
    det_v <- if ("determiner_tier" %in% names(vis_df)) as.character(vis_df$determiner_tier) else rep("Tier_1_Gold", nrow(vis_df))
    v_anch <- vis_df[det_v %in% c("Tier_1_Gold", "Tier_2_Silver"), , drop = FALSE]
    if (nrow(v_anch) < 4L) v_anch <- vis_df
    kappa_vis <- compute_cohens_kappa_stat(v_anch$species_standardized, v_anch$vision_predicted_label)
  }
  w_vis <- if (has_vis) {
    if (is.finite(kappa_vis)) max(0.05, min(1.0, kappa_vis)) else 0.65
  } else {
    NA_real_
  }

  # Dynamic renormalization over present (non-NA) modalities
  present_w <- c(morph = w_morph, vis = w_vis)
  avail <- !is.na(present_w)
  norm_w <- present_w
  norm_w[avail] <- present_w[avail] / sum(present_w[avail])
  norm_w[!avail] <- 0.0
  list(w_morph_raw = w_morph, w_vis_raw = w_vis, w_morph = unname(norm_w["morph"]), w_vis = unname(norm_w["vis"]))
}

# ------------------------------------------------------------------------------
# 4. Step C: Morphometric & Visual Affinity Matrices (S^morph_ij, S^vis_ij)
# ------------------------------------------------------------------------------
compute_morphometric_affinity <- function(morph_df, anchor_df, taxa) {
  n_taxa <- length(taxa)
  sym_cols <- intersect(paste0("Sym_PC", 1:5), names(anchor_df))
  pc_cols <- intersect(paste0("PC", 1:5), names(anchor_df))
  can_cols <- intersect(c("can1", "can2", "can3"), names(anchor_df))

  feat_cols <- if (length(sym_cols) >= 2L) sym_cols else if (length(pc_cols) >= 2L) pc_cols else can_cols
  repro_cols <- intersect(c("involucre_height_mm", "involucre_width_mm", "capitulum_aspect_ratio", "AR_inv_std", "W_inv_std"), names(anchor_df))
  for (rc in repro_cols) {
    vals <- suppressWarnings(as.numeric(anchor_df[[rc]]))
    if (sum(is.finite(vals) & vals > 0) >= 0.5 * nrow(anchor_df) && stats::sd(vals, na.rm = TRUE) > 1e-6) {
      col_nm <- paste0("log_", rc)
      anchor_df[[col_nm]] <- log10(pmax(suppressWarnings(as.numeric(anchor_df[[rc]])), 1e-3))
      morph_df[[col_nm]] <- log10(pmax(suppressWarnings(as.numeric(morph_df[[rc]])), 1e-3))
      feat_cols <- c(feat_cols, col_nm)
    }
  }

  X_anchor <- as.matrix(anchor_df[, feat_cols, drop = FALSE]); storage.mode(X_anchor) <- "double"
  X_all <- as.matrix(morph_df[, feat_cols, drop = FALSE]); storage.mode(X_all) <- "double"
  for (j in seq_len(ncol(X_anchor))) {
    med_j <- stats::median(X_anchor[, j], na.rm = TRUE)
    if (!is.finite(med_j)) med_j <- 0.0
    X_anchor[!is.finite(X_anchor[, j]), j] <- med_j
    X_all[!is.finite(X_all[, j]), j] <- med_j
  }

  p <- ncol(X_anchor)
  centroids <- matrix(0.0, nrow = n_taxa, ncol = p, dimnames = list(taxa, feat_cols))
  S_pooled <- matrix(0.0, nrow = p, ncol = p)
  df_total <- 0L
  for (tx in taxa) {
    sub_x <- X_anchor[anchor_df$species_standardized == tx, , drop = FALSE]
    if (nrow(sub_x) == 0L) sub_x <- X_all[morph_df$species_standardized == tx, , drop = FALSE]
    centroids[tx, ] <- colMeans(sub_x, na.rm = TRUE)
    if (nrow(sub_x) >= 2L) {
      cov_tx <- stats::cov(sub_x); cov_tx[!is.finite(cov_tx)] <- 0.0
      S_pooled <- S_pooled + (nrow(sub_x) - 1L) * cov_tx
      df_total <- df_total + (nrow(sub_x) - 1L)
    }
  }
  S_pooled <- if (df_total > 0L) S_pooled / df_total else stats::cov(X_anchor)
  S_pooled[!is.finite(S_pooled)] <- 0.0
  diag_scale <- max(mean(diag(S_pooled), na.rm = TRUE), 1e-4)
  S_reg <- S_pooled + diag(max(1e-4, 0.02 * diag_scale), nrow = p, ncol = p)
  S_inv <- tryCatch(solve(S_reg), error = function(e) MASS::ginv(S_reg))

  D2_mahal <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) {
    for (j in seq_len(n_taxa)) {
      if (i != j) {
        diff_ij <- centroids[i, ] - centroids[j, ]
        D2_mahal[i, j] <- max(0.0, as.numeric(t(diff_ij) %*% S_inv %*% diff_ij))
      }
    }
  }
  off_diag_d2 <- D2_mahal[upper.tri(D2_mahal)]
  med_d2 <- stats::median(off_diag_d2[off_diag_d2 > 0], na.rm = TRUE)
  if (!is.finite(med_d2) || med_d2 <= 1e-8) med_d2 <- 1.0
  K_mahal <- exp(-D2_mahal / (2.0 * med_d2)); diag(K_mahal) <- 0.0

  C_morph <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]
    idx_i <- which(morph_df$species_standardized == tx_i)
    if (length(idx_i) == 0L) next
    X_i <- X_all[idx_i, , drop = FALSE]
    cda_pred_i <- if ("cda_predicted_taxon" %in% names(morph_df)) as.character(morph_df$cda_predicted_taxon[idx_i]) else rep(tx_i, length(idx_i))

    log_likes <- matrix(0.0, nrow = length(idx_i), ncol = n_taxa)
    for (j in seq_len(n_taxa)) {
      diff_mat <- sweep(X_i, 2, centroids[j, ], "-")
      log_likes[, j] <- -0.5 * pmax(0.0, rowSums((diff_mat %*% S_inv) * diff_mat))
    }
    probs_i <- exp(sweep(log_likes, 1, apply(log_likes, 1, max), "-"))
    probs_i <- sweep(probs_i, 1, pmax(rowSums(probs_i), 1e-12), "/")

    for (j in seq_len(n_taxa)) {
      if (i == j) next
      tx_j <- taxa[j]
      C_morph[i, j] <- mean((!is.na(cda_pred_i) & cda_pred_i == tx_j) | (probs_i[, j] >= 0.15), na.rm = TRUE)
    }
  }
  C_morph_sym <- 0.5 * (C_morph + t(C_morph)); diag(C_morph_sym) <- 0.0
  S_morph <- 0.60 * K_mahal + 0.40 * C_morph_sym; diag(S_morph) <- 0.0
  list(S_morph = S_morph, K_mahal = K_mahal, C_morph_sym = C_morph_sym, D2_mahal = D2_mahal)
}

compute_visual_affinity <- function(vis_df, taxa) {
  n_taxa <- length(taxa)
  C_vis <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  if (is.null(vis_df) || nrow(vis_df) == 0L) return(list(S_vis = C_vis, C_vis = C_vis))

  c_err <- if ("c_error" %in% names(vis_df)) suppressWarnings(as.numeric(vis_df$c_error)) else rep(0.5, nrow(vis_df))
  c_err[!is.finite(c_err)] <- 0.5
  det_tier <- if ("determiner_tier" %in% names(vis_df)) as.character(vis_df$determiner_tier) else rep("Tier_1_Gold", nrow(vis_df))
  is_anchor <- det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")

  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]
    idx_i <- which(vis_df$species_standardized == tx_i)
    if (length(idx_i) == 0L) next
    pred_i <- as.character(vis_df$vision_predicted_label[idx_i])
    cerr_i <- c_err[idx_i]; anch_i <- is_anchor[idx_i]
    w_i <- ifelse(!is.na(pred_i) & pred_i != tx_i, ifelse(anch_i, pmax(1.0 - cerr_i, cerr_i), pmax(cerr_i, 0.25)), 1.0 - cerr_i)
    w_i[!is.finite(w_i) | w_i <= 0] <- 0.5
    for (j in seq_len(n_taxa)) {
      if (i == j) next
      match_j <- !is.na(pred_i) & (pred_i == taxa[j])
      if (any(match_j)) C_vis[i, j] <- sum(w_i[match_j], na.rm = TRUE) / length(idx_i)
    }
  }
  S_vis <- 0.5 * (C_vis + t(C_vis)); diag(S_vis) <- 0.0
  list(S_vis = S_vis, C_vis = C_vis)
}

# ------------------------------------------------------------------------------
# 5. Step D: Hopkins-Corrected Circular Anthesis Phenology (Hopkins 1918; Pewsey et al. 2013)
# ------------------------------------------------------------------------------
compute_hopkins_phenological_affinity <- function(morph_df, taxa) {
  # Gates phenological divergence strictly on active anthesis
  # (capitulum_phenology_classifier.py) and computes shortest-arc circular
  # residuals (atan2) from a Hopkins bioclimatic model:
  #   sin_doy, cos_doy ~ decimalLatitude + elevation + bio1
  n_taxa <- length(taxa)
  S_pheno <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  if (!("doy" %in% names(morph_df))) {
    return(list(S_pheno = S_pheno, has_pheno = FALSE, delta_doy_resid = rep(NA_real_, nrow(morph_df))))
  }
  doy_vec <- suppressWarnings(as.numeric(morph_df$doy))
  anthesis_mask <- is.finite(doy_vec) & doy_vec >= 1 & doy_vec <= 366
  for (pcol in c("voucher_phenological_state", "phenological_Stage", "phenological_state")) {
    if (pcol %in% names(morph_df)) {
      st <- tolower(trimws(as.character(morph_df[[pcol]])))
      has_st <- !is.na(st) & nzchar(st) & st != "na"
      if (any(has_st)) {
        anthesis_mask <- anthesis_mask & (!has_st | st %in% c("anthesis", "flowering", "open_florets", "active_anthesis"))
      }
      break
    }
  }
  if (sum(anthesis_mask) < 4L) {
    return(list(S_pheno = S_pheno, has_pheno = FALSE, delta_doy_resid = rep(NA_real_, nrow(morph_df))))
  }

  lat_col <- intersect(c("decimalLatitude", "latitude"), names(morph_df))
  lat_vec <- if (length(lat_col) > 0L) suppressWarnings(as.numeric(morph_df[[lat_col[1L]]])) else rep(NA_real_, nrow(morph_df))
  elev_col <- intersect(c("topo_elev_m", "elevation_m", "elevation"), names(morph_df))
  elev_vec <- if (length(elev_col) > 0L) suppressWarnings(as.numeric(morph_df[[elev_col[1L]]])) else rep(0.0, nrow(morph_df))
  elev_vec[!is.finite(elev_vec)] <- 0.0
  bio1_col <- intersect(c("bio1_temp_mean", "bio1", "mat_c"), names(morph_df))
  bio1_vec <- if (length(bio1_col) > 0L) suppressWarnings(as.numeric(morph_df[[bio1_col[1L]]])) else rep(NA_real_, nrow(morph_df))

  theta_obs <- 2.0 * pi * (doy_vec / 365.25)
  sin_doy <- sin(theta_obs)
  cos_doy <- cos(theta_obs)

  ref_lat <- if (any(is.finite(lat_vec[anthesis_mask]))) stats::median(lat_vec[anthesis_mask], na.rm = TRUE) else 35.0
  ref_elev <- if (any(is.finite(elev_vec[anthesis_mask]))) stats::median(elev_vec[anthesis_mask], na.rm = TRUE) else 150.0
  ref_bio1 <- if (any(is.finite(bio1_vec[anthesis_mask]))) stats::median(bio1_vec[anthesis_mask], na.rm = TRUE) else 15.0

  lat_shift <- ifelse(is.finite(lat_vec), 4.0 * (lat_vec - ref_lat), 0.0)
  elev_shift <- 1.25 * ((elev_vec - ref_elev) / 100.0)
  bio1_shift <- ifelse(is.finite(bio1_vec), -1.5 * (bio1_vec - ref_bio1), 0.0)
  hopkins_shift <- 0.70 * (lat_shift + elev_shift) + 0.30 * bio1_shift
  adj_doy <- doy_vec - hopkins_shift
  theta_adj <- 2.0 * pi * (adj_doy / 365.25)

  # Circular mean per taxon & shortest-arc atan2 residual (Pewsey et al. 2013)
  mu_theta <- numeric(n_taxa); names(mu_theta) <- taxa
  delta_doy_resid <- rep(NA_real_, nrow(morph_df))
  for (tx in taxa) {
    idx_tx <- which(morph_df$species_standardized == tx & anthesis_mask)
    if (length(idx_tx) > 0L) {
      mu_tx <- atan2(mean(sin(theta_adj[idx_tx]), na.rm = TRUE), mean(cos(theta_adj[idx_tx]), na.rm = TRUE))
      mu_theta[tx] <- mu_tx
      diff_rad <- atan2(sin(theta_adj[idx_tx] - mu_tx), cos(theta_adj[idx_tx] - mu_tx))
      delta_doy_resid[idx_tx] <- abs(diff_rad) * (365.25 / (2.0 * pi))
    } else {
      mu_theta[tx] <- 2.0 * pi * (120.0 / 365.25)
    }
  }

  raw_p <- numeric(0)
  for (i in seq_len(n_taxa)) {
    for (j in seq_len(n_taxa)) {
      if (i != j) {
        arc_rad <- abs(atan2(sin(mu_theta[i] - mu_theta[j]), cos(mu_theta[i] - mu_theta[j])))
        delta_circ <- arc_rad * (365.25 / (2.0 * pi))
        S_pheno[i, j] <- exp(-(delta_circ^2) / (2.0 * (21.0^2)))
        if (i < j) raw_p <- c(raw_p, exp(-(delta_circ / 18.0)))
      }
    }
  }
  pheno_fdr_p <- if (length(raw_p) > 0L) stats::p.adjust(raw_p, method = "BH") else numeric(0)
  diag(S_pheno) <- 0.0
  list(S_pheno = S_pheno, has_pheno = TRUE, sin_doy = sin_doy, cos_doy = cos_doy,
       delta_doy_resid = delta_doy_resid, pheno_fdr_p = pheno_fdr_p)
}

# ------------------------------------------------------------------------------
# 6. Step E: Five-Stage Hierarchical Integrative Decision Matrix
#    (Padial et al. 2010; Schlick-Steiner et al. 2010; Rieseberg 1995; Slenker et al. 2022)
# ------------------------------------------------------------------------------
evaluate_hierarchical_complex_delimitation <- function(morph_df, vis_df, pheno_res) {
  n <- nrow(morph_df)
  sp_given <- as.character(morph_df$species_standardized)
  cda_pred <- if ("cda_predicted_taxon" %in% names(morph_df)) as.character(morph_df$cda_predicted_taxon) else sp_given
  cda_prob <- if ("cda_posterior_prob" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_posterior_prob)) else rep(0.90, n)
  cda_prob[!is.finite(cda_prob)] <- 0.90
  cda_sec <- if ("cda_second_posterior_prob" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_second_posterior_prob)) else pmax(0.0, 1.0 - cda_prob)
  cda_sec[!is.finite(cda_sec)] <- pmax(0.0, 1.0 - cda_prob[!is.finite(cda_sec)])
  cda_margin <- if ("cda_posterior_margin" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_posterior_margin)) else abs(cda_prob - cda_sec)
  cda_margin[!is.finite(cda_margin)] <- abs(cda_prob[!is.finite(cda_margin)] - cda_sec[!is.finite(cda_margin)])

  vis_pred <- sp_given
  c_err <- rep(0.10, n)
  if (!is.null(vis_df) && nrow(vis_df) > 0L && "catalogNumber" %in% names(vis_df)) {
    idx_v <- match(morph_df$catalogNumber, as.character(vis_df$catalogNumber))
    if ("vision_predicted_label" %in% names(vis_df)) {
      vp <- as.character(vis_df$vision_predicted_label)[idx_v]
      vis_pred[!is.na(vp) & nzchar(vp)] <- vp[!is.na(vp) & nzchar(vp)]
    }
    if ("c_error" %in% names(vis_df)) {
      ce <- suppressWarnings(as.numeric(vis_df$c_error))[idx_v]
      c_err[is.finite(ce)] <- ce[is.finite(ce)]
    }
  }

  m_flag <- if ("morphometric_outlier_flag" %in% names(morph_df)) as.logical(morph_df$morphometric_outlier_flag) else rep(FALSE, n)
  m_flag[is.na(m_flag)] <- FALSE
  gmm_cl <- if ("gmm_cluster" %in% names(morph_df)) suppressWarnings(as.integer(morph_df$gmm_cluster)) else rep(1L, n)
  gmm_cl[is.na(gmm_cl)] <- 1L
  n_known_sp <- max(2L, length(unique(sp_given)))
  gmm_bf <- if ("gmm_2_delta_bic" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$gmm_2_delta_bic)) else ifelse(gmm_cl > n_known_sp, 8.0, 0.0)
  gmm_bf[!is.finite(gmm_bf)] <- ifelse(gmm_cl[!is.finite(gmm_bf)] > n_known_sp, 8.0, 0.0)
  d_edaph <- if ("edaphic_schoener_d" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$edaphic_schoener_d)) else rep(0.65, n)
  d_edaph[!is.finite(d_edaph)] <- 0.65
  d_doy <- pheno_res$delta_doy_resid

  is_t1_pure <- rep(TRUE, n)
  if ("n_tier1_leaves" %in% names(morph_df)) {
    nt1 <- suppressWarnings(as.numeric(morph_df$n_tier1_leaves))
    if (any(!is.na(nt1))) is_t1_pure <- !is.na(nt1) & (nt1 >= 1)
  }

  morph_disagrees <- !is.na(cda_pred) & (cda_pred != sp_given) & (cda_prob >= 0.50)
  vis_disagrees <- !is.na(vis_pred) & (vis_pred != sp_given) & (c_err >= 0.55)
  both_alt_agree <- morph_disagrees & vis_disagrees & (cda_pred == vis_pred)
  inter_modal_mosaic <- morph_disagrees & vis_disagrees & (cda_pred != vis_pred)
  single_strong_reassign <- !inter_modal_mosaic & (
    (morph_disagrees & cda_prob >= 0.80 & (!vis_disagrees | cda_pred == vis_pred)) |
    (vis_disagrees & c_err >= 0.85 & (!morph_disagrees | cda_pred == vis_pred)))
  # Stage 1: Misidentification_Reassign (evaluated BEFORE hybridity)
  is_stage1_misid <- both_alt_agree | single_strong_reassign
  # Stage 2: Hybrid_Introgressant (Rieseberg 1995; Slenker et al. 2022)
  is_stage2_hybrid <- !is_stage1_misid & (inter_modal_mosaic | (cda_margin < 0.25 & (morph_disagrees | vis_disagrees) & d_edaph >= 0.30))
  # Stage 3: Putative_Cryptic_Species via Integration by Congruence (Padial et al. 2010)
  has_evol_barrier <- (!is.na(d_doy) & d_doy >= 18.0) | (d_edaph < 0.30)
  is_stage3_cryptic <- !is_stage1_misid & !is_stage2_hybrid & is_t1_pure & (cda_prob < 0.70) & (c_err < 0.75) &
    (m_flag | (gmm_cl > length(unique(sp_given)) & gmm_bf >= 6.0)) & has_evol_barrier
  # Stage 4: Subspecies_Ecotype vs. Ecophenotypic_Plasticity
  is_stage4_ecotype <- !is_stage1_misid & !is_stage2_hybrid & !is_stage3_cryptic &
    ((d_edaph < 0.30) | (!is.na(d_doy) & d_doy >= 14.0)) & (m_flag | gmm_cl > 1L)
  is_stage4_plastic <- !is_stage1_misid & !is_stage2_hybrid & !is_stage3_cryptic & !is_stage4_ecotype & m_flag

  ifelse(is_stage1_misid, "Misidentification_Reassign",
    ifelse(is_stage2_hybrid, "Hybrid_Introgressant",
      ifelse(is_stage3_cryptic, "Putative_Cryptic_Species",
        ifelse(is_stage4_ecotype, "Subspecies_Ecotype",
          ifelse(is_stage4_plastic, "Ecophenotypic_Plasticity", "Confirmed_Species")))))
}

# ------------------------------------------------------------------------------
# 7. Step F: Consensus Affinity Matrix & Unsupervised Graph / Ward Partitioning
# ------------------------------------------------------------------------------
relabel_clusters_sequential <- function(cl_vec, taxa) match(cl_vec, unique(cl_vec[order(taxa)]))

partition_taxon_complexes <- function(W, taxa, n_vouchers_vec, n_gs_vec) {
  n_taxa <- length(taxa); diag(W) <- 0.0
  D_mat <- 1.0 - pmin(pmax(W, 0.0), 1.0); diag(D_mat) <- 0.0
  dist_obj <- stats::as.dist(D_mat)

  off_vals <- W[upper.tri(W)]
  tau <- min(0.03, stats::quantile(off_vals[off_vals > 0], probs = 0.50, na.rm = TRUE))
  if (!is.finite(tau) || tau <= 0) tau <- 0.01
  W_adj <- W; W_adj[W_adj < tau] <- 0.0
  for (i in seq_len(n_taxa)) {
    if (sum(W_adj[i, ] > 0) == 0L && max(W[i, ]) > 0) {
      best_j <- which.max(W[i, ]); W_adj[i, best_j] <- W[i, best_j]; W_adj[best_j, i] <- W[best_j, i]
    }
  }

  leiden_comm <- rep(1L, n_taxa); best_mod <- -Inf
  if (requireNamespace("igraph", quietly = TRUE) && sum(W_adj > 0) > 0) {
    g <- igraph::graph_from_adjacency_matrix(W_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    for (gm in c(0.6, 0.8, 1.0, 1.2)) {
      cl_try <- tryCatch(suppressWarnings(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution = gm, n_iterations = 10)),
        error = function(e) tryCatch(suppressWarnings(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution_parameter = gm, n_iterations = 10)), error = function(e2) NULL))
      if (!is.null(cl_try)) {
        memb <- as.integer(igraph::membership(cl_try)); k_m <- length(unique(memb))
        q_val <- tryCatch(igraph::modularity(g, memb, weights = igraph::E(g)$weight), error = function(e) -1.0)
        if (k_m >= 2L && k_m <= min(8L, n_taxa - 1L) && is.finite(q_val) && q_val > best_mod) { best_mod <- q_val; leiden_comm <- memb }
      }
    }
    if (length(unique(leiden_comm)) < 2L) {
      lv_try <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$weight), error = function(e) NULL)
      if (!is.null(lv_try)) leiden_comm <- as.integer(igraph::membership(lv_try))
    }
  }
  leiden_comm <- relabel_clusters_sequential(leiden_comm, taxa)

  hc <- stats::hclust(dist_obj, method = "ward.D2"); k_max <- min(8L, n_taxa - 1L)
  ward_cluster <- rep(1L, n_taxa); best_sil_mean <- -Inf; best_sil_vec <- rep(0.0, n_taxa)
  if (n_taxa >= 3L && k_max >= 2L) {
    k_candidates <- 2L:k_max
    sil_scores <- vapply(k_candidates, function(kv) mean(cluster::silhouette(stats::cutree(hc, k = kv), dist_obj)[, "sil_width"], na.rm = TRUE), numeric(1))
    best_idx <- which.max(sil_scores); best_k <- k_candidates[best_idx]; best_sil_mean <- sil_scores[best_idx]
    ward_cluster <- relabel_clusters_sequential(as.integer(stats::cutree(hc, k = best_k)), taxa)
    best_sil_vec <- as.numeric(cluster::silhouette(ward_cluster, dist_obj)[, "sil_width"])
  }

  leiden_sil_mean <- -Inf; leiden_sil_vec <- best_sil_vec
  if (length(unique(leiden_comm)) >= 2L && length(unique(leiden_comm)) < n_taxa) {
    l_sil <- cluster::silhouette(leiden_comm, dist_obj)
    leiden_sil_mean <- mean(l_sil[, "sil_width"], na.rm = TRUE); leiden_sil_vec <- as.numeric(l_sil[, "sil_width"])
  }
  if (leiden_sil_mean > best_sil_mean + 0.02) { complex_id <- leiden_comm; sil_width_vec <- leiden_sil_vec }
  else { complex_id <- ward_cluster; sil_width_vec <- best_sil_vec }
  complex_id <- relabel_clusters_sequential(complex_id, taxa)

  complex_labels <- character(n_taxa); deg_vec <- rowSums(W)
  for (cid in sort(unique(complex_id))) {
    m_idx <- which(complex_id == cid)
    ord <- m_idx[order(-deg_vec[m_idx], -n_vouchers_vec[m_idx], taxa[m_idx])]
    top_ep <- gsub("\\s*\\[K([0-9]+)\\]", "_K\\1", sub("^(?:Packera|P\\.)\\s+", "", taxa[ord[seq_len(min(2L, length(ord)))]]))
    complex_labels[m_idx] <- paste0("Complex_", cid, "_", paste(top_ep, collapse = "_"))
  }

  top_n1 <- character(n_taxa); top_a1 <- numeric(n_taxa); top_n2 <- character(n_taxa); top_a2 <- numeric(n_taxa); within_mean_aff <- numeric(n_taxa)
  for (i in seq_len(n_taxa)) {
    w_row <- W[i, ]; w_row[i] <- -1.0; ord_n <- order(w_row, decreasing = TRUE)
    top_n1[i] <- if (length(ord_n) >= 1L) taxa[ord_n[1L]] else ""; top_a1[i] <- if (length(ord_n) >= 1L) round(max(0.0, W[i, ord_n[1L]]), 4) else 0.0
    top_n2[i] <- if (length(ord_n) >= 2L) taxa[ord_n[2L]] else ""; top_a2[i] <- if (length(ord_n) >= 2L) round(max(0.0, W[i, ord_n[2L]]), 4) else 0.0
    peers_i <- which(complex_id == complex_id[i] & seq_len(n_taxa) != i)
    within_mean_aff[i] <- if (length(peers_i) > 0L) round(mean(W[i, peers_i], na.rm = TRUE), 4) else 0.0
  }

  complexes_df <- data.frame(
    species_standardized = taxa, taxon = taxa,
    n_vouchers = as.integer(n_vouchers_vec[taxa]), n_gold_silver_anchors = as.integer(n_gs_vec[taxa]), n_Pure_Gold_Silver = as.integer(n_gs_vec[taxa]),
    complex_id = as.integer(complex_id), consensus_complex_id = as.integer(complex_id), complex_label = complex_labels,
    leiden_community = as.integer(leiden_comm), leiden_complex_id = as.integer(leiden_comm),
    ward_silhouette_cluster = as.integer(ward_cluster), ward_complex_id = as.integer(ward_cluster),
    silhouette_width = round(sil_width_vec, 4), within_complex_mean_affinity = within_mean_aff,
    top_neighbor_1 = top_n1, nearest_neighbor_taxon = top_n1, top_neighbor_1_affinity = top_a1, nearest_neighbor_affinity = top_a1,
    top_neighbor_2 = top_n2, top_neighbor_2_affinity = top_a2, stringsAsFactors = FALSE
  )
  list(complexes_df = complexes_df, hc = hc, dist_obj = dist_obj)
}

# ------------------------------------------------------------------------------
# 8. Step G: Taxon x Morphotype Subcluster Builder & Publication Figure Wrapper
# ------------------------------------------------------------------------------
build_taxon_morphotype_subclusters <- function(morph_df, vis_df = NULL, min_vouchers = 5L) {
  raw_sp <- as.character(morph_df$species_standardized)
  u_sp <- sort(unique(raw_sp[!is.na(raw_sp) & nzchar(raw_sp)]))
  n_sp <- length(u_sp)
  has_gmm <- "gmm_cluster" %in% names(morph_df) && length(unique(stats::na.omit(morph_df$gmm_cluster))) >= 2L
  within_sp_gmm <- if (has_gmm) any(vapply(u_sp, function(sp) {
    ks <- stats::na.omit(as.integer(morph_df$gmm_cluster[raw_sp == sp]))
    length(unique(ks)) >= 2L && sum(table(ks) >= min_vouchers) >= 2L
  }, logical(1))) else FALSE
  if (n_sp > 2L && !within_sp_gmm) return(list(morph_df = morph_df, vis_df = vis_df, subcluster_mode = FALSE))

  k_vec <- if (has_gmm) ifelse(is.na(morph_df$gmm_cluster), 1L, as.integer(morph_df$gmm_cluster)) else rep(1L, nrow(morph_df))
  short_sp <- sub("^Packera\\s+", "P. ", raw_sp)
  cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
  if (sum(table(cand_lbl) >= min_vouchers) < 3L && nrow(morph_df) >= 3L * min_vouchers) {
    fcols <- intersect(c(paste0("Sym_PC", 1:5), paste0("PC", 1:5), "can1", "can2", "can3"), names(morph_df))
    X_m <- scale(as.matrix(morph_df[, fcols, drop = FALSE])); X_m[!is.finite(X_m)] <- 0.0
    for (sp in u_sp) {
      idx_s <- which(raw_sp == sp)
      k_s <- if (n_sp == 1L) min(3L, max(2L, floor(length(idx_s) / min_vouchers))) else min(2L, max(1L, floor(length(idx_s) / min_vouchers)))
      if (k_s >= 2L && length(unique(k_vec[idx_s])) < k_s) {
        km <- tryCatch(stats::kmeans(X_m[idx_s, , drop = FALSE], centers = k_s, nstart = 10L), error = function(e) NULL)
        if (!is.null(km)) k_vec[idx_s] <- as.integer(km$cluster)
      }
    }
    cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
  }
  tab_c <- table(cand_lbl)
  valid_c <- names(tab_c)[tab_c >= min_vouchers]
  if (length(valid_c) < 2L) valid_c <- names(tab_c)[tab_c >= max(1L, min(tab_c))]
  keep_m <- cand_lbl %in% valid_c
  morph_df <- morph_df[keep_m, , drop = FALSE]
  morph_df$base_species <- raw_sp[keep_m]
  morph_df$species_standardized <- cand_lbl[keep_m]
  if ("cda_predicted_taxon" %in% names(morph_df)) {
    cp <- std_taxon(morph_df$cda_predicted_taxon); cp_short <- sub("^Packera\\s+", "P. ", cp)
    morph_df$cda_predicted_taxon <- vapply(seq_len(nrow(morph_df)), function(i) {
      if (!is.na(cp[i]) && cp[i] %in% valid_c) return(cp[i])
      if (!is.na(cp[i]) && cp[i] == morph_df$base_species[i]) return(morph_df$species_standardized[i])
      m_alt <- valid_c[startsWith(valid_c, paste0(cp_short[i], " ["))]
      if (length(m_alt) > 0L) m_alt[1L] else morph_df$species_standardized[i]
    }, character(1))
  }
  if (!is.null(vis_df) && nrow(vis_df) > 0L && "key" %in% names(vis_df)) {
    vis_df <- vis_df[vis_df$key %in% morph_df$key, , drop = FALSE]
    if (nrow(vis_df) > 0L) {
      key_map <- stats::setNames(morph_df$species_standardized, morph_df$key)
      base_map <- stats::setNames(morph_df$base_species, morph_df$key)
      vis_df$species_standardized <- unname(key_map[vis_df$key])
      if ("vision_predicted_label" %in% names(vis_df)) {
        vp <- std_taxon(vis_df$vision_predicted_label); vp_short <- sub("^Packera\\s+", "P. ", vp)
        vis_df$vision_predicted_label <- vapply(seq_len(nrow(vis_df)), function(i) {
          if (!is.na(vp[i]) && vp[i] %in% valid_c) return(vp[i])
          if (!is.na(vp[i]) && vp[i] == unname(base_map[vis_df$key[i]])) return(vis_df$species_standardized[i])
          m_alt <- valid_c[startsWith(valid_c, paste0(vp_short[i], " ["))]
          if (length(m_alt) > 0L) m_alt[1L] else vis_df$species_standardized[i]
        }, character(1))
      }
    }
  }
  message(sprintf("[Subcluster Mode] Built %d Taxon x Morphotype clusters (>= %d vouchers): %s",
                  length(valid_c), min_vouchers, paste(valid_c, collapse = ", ")))
  list(morph_df = morph_df, vis_df = vis_df, subcluster_mode = TRUE)
}

render_complex_discovery_figure <- function(W, part_res, morph_res, vis_res, out_pdf) {
  if (exists("render_step08_complex_discovery_pdf", mode = "function")) {
    return(render_step08_complex_discovery_pdf(W = W, part_res = part_res, out_pdf = out_pdf))
  }
  invisible(out_pdf)
}

# ------------------------------------------------------------------------------
# 9. Main Execution Entrypoint
# ------------------------------------------------------------------------------
main <- function() {
  set.seed(42)
  opts <- parse_discovery_args()
  message("=== Phase 08: Unsupervised Quantitative Taxon-Complex Discovery Engine ===")

  inputs <- load_and_merge_inputs(opts)
  sub_res <- build_taxon_morphotype_subclusters(inputs$morph_df, inputs$vis_df, min_vouchers = 5L)
  morph_df <- sub_res$morph_df
  vis_df <- sub_res$vis_df

  taxa <- sort(unique(morph_df$species_standardized))
  n_taxa <- length(taxa)
  if (n_taxa < 2L) {
    stop("At least 2 distinct Packera taxa or morphotype clusters are required for discovery (found ", n_taxa, ").")
  }

  n_vouchers_vec <- table(morph_df$species_standardized)[taxa]
  anc <- select_taxon_anchors(morph_df)
  kappa_wts <- compute_modality_kappa_weights(morph_df, vis_df, anc$anchor_mask)
  message(sprintf(
    "Evaluated %d vouchers across %d Packera groups (%d Tier_1_Pure Gold/Silver anchors; kappa weights: morph=%.3f, vis=%.3f).",
    nrow(morph_df), n_taxa, nrow(anc$anchor_df), kappa_wts$w_morph, kappa_wts$w_vis
  ))

  morph_res <- compute_morphometric_affinity(morph_df, anc$anchor_df, taxa)
  vis_res <- compute_visual_affinity(vis_df, taxa)
  pheno_res <- compute_hopkins_phenological_affinity(morph_df, taxa)
  specimen_subtypes <- evaluate_hierarchical_complex_delimitation(morph_df, vis_df, pheno_res)

  if (!is.null(vis_df) && nrow(vis_df) > 0L && any(vis_res$S_vis > 0)) {
    W <- kappa_wts$w_morph * morph_res$S_morph + kappa_wts$w_vis * vis_res$S_vis
  } else {
    W <- morph_res$S_morph
  }
  diag(W) <- 0.0
  W <- round(W, 6)

  part_res <- partition_taxon_complexes(W = W, taxa = taxa, n_vouchers_vec = n_vouchers_vec, n_gs_vec = anc$n_gold_silver)
  cdf <- part_res$complexes_df
  hop_ov <- numeric(n_taxa); int_cong <- logical(n_taxa)
  n_misid <- integer(n_taxa); n_hyb <- integer(n_taxa); n_crypt <- integer(n_taxa)
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]; idx_tx <- which(morph_df$species_standardized == tx_i)
    n_misid[i] <- sum(specimen_subtypes[idx_tx] == "Misidentification_Reassign", na.rm = TRUE)
    n_hyb[i] <- sum(specimen_subtypes[idx_tx] == "Hybrid_Introgressant", na.rm = TRUE)
    n_crypt[i] <- sum(specimen_subtypes[idx_tx] == "Putative_Cryptic_Species", na.rm = TRUE)
    nn_tx <- cdf$nearest_neighbor_taxon[i]
    if (nzchar(nn_tx) && nn_tx %in% taxa) {
      hop_ov[i] <- if (isTRUE(pheno_res$has_pheno)) round(pheno_res$S_pheno[tx_i, nn_tx], 4) else NA_real_
      n_mod <- as.integer(morph_res$S_morph[tx_i, nn_tx] >= opts$min_affinity) +
               as.integer(vis_res$S_vis[tx_i, nn_tx] >= opts$min_affinity) +
               as.integer(isTRUE(pheno_res$has_pheno) && pheno_res$S_pheno[tx_i, nn_tx] >= 0.25)
      int_cong[i] <- n_mod >= 2L
    }
  }
  part_res$complexes_df$hopkins_pheno_overlap <- hop_ov
  part_res$complexes_df$integrative_congruence <- int_cong
  part_res$complexes_df$kappa_morph_weight <- round(kappa_wts$w_morph, 4)
  part_res$complexes_df$kappa_vis_weight <- round(kappa_wts$w_vis, 4)
  part_res$complexes_df$n_misidentification_reassign <- n_misid
  part_res$complexes_df$n_hybrid_introgressant <- n_hyb
  part_res$complexes_df$n_putative_cryptic <- n_crypt

  atomic_write_csv(part_res$complexes_df, opts$output_complexes)
  message("Saved empirical taxon complexes -> ", opts$output_complexes)

  aff_df <- data.frame(taxon = taxa, species_standardized = taxa, as.data.frame(W, check.names = FALSE), stringsAsFactors = FALSE, check.names = FALSE)
  atomic_write_csv(aff_df, opts$output_affinity)
  message("Saved pairwise taxon affinity matrix -> ", opts$output_affinity)

  focal_raw <- trimws(as.character(opts$focal_taxon))
  if (nzchar(focal_raw)) {
    focal_std <- std_taxon(focal_raw)
    focal_short <- sub("^Packera\\s+", "P. ", focal_std)
    focal_nodes <- taxa[taxa == focal_std | startsWith(taxa, paste0(focal_short, " ["))]
    if (length(focal_nodes) > 0L) {
      cdf <- part_res$complexes_df
      focal_cids <- unique(cdf$complex_id[cdf$species_standardized %in% focal_nodes])
      same_complex <- cdf$species_standardized[cdf$complex_id %in% focal_cids]
      aff_mask <- apply(morph_res$C_morph_sym[focal_nodes, , drop = FALSE] >= opts$min_affinity |
                        vis_res$S_vis[focal_nodes, , drop = FALSE] >= opts$min_affinity, 2, any)
      selected_neighbors <- sort(setdiff(union(same_complex, taxa[aff_mask]), focal_nodes))
    } else selected_neighbors <- character(0)
    atomic_write_lines(selected_neighbors, opts$output_neighbors)
    message("Saved focal auto-neighbors (", length(selected_neighbors), " taxa) -> ", opts$output_neighbors)
    cat(paste(selected_neighbors, collapse = "\n"), "\n", sep = "")
  } else atomic_write_lines(part_res$complexes_df$species_standardized, opts$output_neighbors)

  render_complex_discovery_figure(W, part_res, morph_res, vis_res, opts$output_figure)
  message("Saved taxon-complex discovery figure -> ", opts$output_figure)
}

if (sys.nframe() == 0L) main()
