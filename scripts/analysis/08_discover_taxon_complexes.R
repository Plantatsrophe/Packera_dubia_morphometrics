#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/08_discover_taxon_complexes.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Unsupervised Quantitative Taxon-Complex Discovery Engine (Phase 08).
#   Computes an empirical N x N Taxon Affinity Matrix (W_ij) fusing:
#     1. Morphospace Overlap (S^morph_ij): Regularized pairwise Mahalanobis
#        affinity kernel (K^mahal_ij) on symmetric Fourier PCs + reproductive
#        traits and symmetrized CDA cross-assignment (C^morph_ij) on Tier_1_Pure
#        Gold/Silver anchors.
#     2. Visual Confusion (S^vis_ij): Symmetrized DINOv2 + Cleanlab out-of-fold
#        visual confusion rates from label_noise_audit.csv.
#   Objectively partitions taxa into empirical species complexes using:
#     - Graph Modularity Maximization (Leiden / Louvain community detection)
#     - Ward.D2 Hierarchical Clustering with Silhouette & Gap Statistic
#   Supports ego-centric --auto-neighbors selection for focal taxon runs and
#   exports a 2-panel publication figure with unfilled/hollow markers.
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
  while (i <= length(raw_args)) {
    a <- raw_args[i]
    if (a %in% c("-m", "--morphometrics") && i < length(raw_args)) {
      opts$morphometrics <- raw_args[i + 1L]; i <- i + 2L
    } else if (a %in% c("-e", "--harmonics") && i < length(raw_args)) {
      opts$harmonics <- raw_args[i + 1L]; i <- i + 2L
    } else if (a %in% c("-n", "--vision-audit") && i < length(raw_args)) {
      opts$vision_audit <- raw_args[i + 1L]; i <- i + 2L
    } else if (a %in% c("-v", "--vouchers") && i < length(raw_args)) {
      opts$vouchers <- raw_args[i + 1L]; i <- i + 2L
    } else if (a == "--focal-taxon" && i < length(raw_args)) {
      opts$focal_taxon <- raw_args[i + 1L]; i <- i + 2L
    } else if (a == "--min-affinity" && i < length(raw_args)) {
      opts$min_affinity <- as.numeric(raw_args[i + 1L]); i <- i + 2L
    } else if (a == "--output-complexes" && i < length(raw_args)) {
      opts$output_complexes <- raw_args[i + 1L]; i <- i + 2L
    } else if (a == "--output-affinity" && i < length(raw_args)) {
      opts$output_affinity <- raw_args[i + 1L]; i <- i + 2L
    } else if (a == "--output-neighbors" && i < length(raw_args)) {
      opts$output_neighbors <- raw_args[i + 1L]; i <- i + 2L
    } else if (a == "--output-figure" && i < length(raw_args)) {
      opts$output_figure <- raw_args[i + 1L]; i <- i + 2L
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
  morph_path <- opts$morphometrics
  if (!file.exists(morph_path)) {
    alt_m <- file.path(dirname(morph_path), "morphometrics_misidentification_flags.csv")
    if (file.exists(alt_m)) morph_path <- alt_m
  }
  if (!file.exists(morph_path)) {
    stop("Morphometrics table not found at: ", opts$morphometrics)
  }
  morph_df <- utils::read.csv(morph_path, stringsAsFactors = FALSE)
  morph_df$catalogNumber <- as.character(morph_df$catalogNumber)

  if (file.exists(opts$harmonics)) {
    harm_df <- tryCatch(utils::read.csv(opts$harmonics, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(harm_df) && nrow(harm_df) > 0 && "catalogNumber" %in% names(harm_df)) {
      harm_df$catalogNumber <- as.character(harm_df$catalogNumber)
      harm_sub <- harm_df[!duplicated(harm_df$catalogNumber), , drop = FALSE]
      add_cols <- setdiff(names(harm_sub), names(morph_df))
      if (length(add_cols) > 0) {
        idx_h <- match(morph_df$catalogNumber, harm_sub$catalogNumber)
        for (cl in add_cols) morph_df[[cl]] <- harm_sub[[cl]][idx_h]
      }
    }
  }

  if (file.exists(opts$vouchers)) {
    v_df <- tryCatch(utils::read.csv(opts$vouchers, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(v_df) && nrow(v_df) > 0 && "catalogNumber" %in% names(v_df)) {
      v_df$catalogNumber <- as.character(v_df$catalogNumber)
      v_sub <- v_df[!duplicated(v_df$catalogNumber), , drop = FALSE]
      add_v <- setdiff(names(v_sub), names(morph_df))
      if (length(add_v) > 0) {
        idx_v <- match(morph_df$catalogNumber, v_sub$catalogNumber)
        for (cl in add_v) morph_df[[cl]] <- v_sub[[cl]][idx_v]
      }
    }
  }

  raw_sp <- if ("species_standardized" %in% names(morph_df)) {
    morph_df$species_standardized
  } else if ("scientificName" %in% names(morph_df)) {
    morph_df$scientificName
  } else if ("species_raw" %in% names(morph_df)) {
    morph_df$species_raw
  } else if ("taxon" %in% names(morph_df)) {
    morph_df$taxon
  } else {
    rep("Unknown", nrow(morph_df))
  }
  morph_df$species_standardized <- std_taxon(raw_sp)
  morph_df <- morph_df[!is.na(morph_df$species_standardized) &
                       nzchar(morph_df$species_standardized) &
                       morph_df$species_standardized != "Unknown", , drop = FALSE]

  if ("cda_predicted_taxon" %in% names(morph_df)) {
    morph_df$cda_predicted_taxon <- std_taxon(morph_df$cda_predicted_taxon)
  }

  vis_df <- NULL
  if (file.exists(opts$vision_audit)) {
    vis_df <- tryCatch(utils::read.csv(opts$vision_audit, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(vis_df) && nrow(vis_df) > 0) {
      vis_sp <- if ("species_standardized" %in% names(vis_df)) {
        vis_df$species_standardized
      } else if ("given_label" %in% names(vis_df)) {
        vis_df$given_label
      } else if ("taxon" %in% names(vis_df)) {
        vis_df$taxon
      } else {
        vis_df$species_raw
      }
      vis_df$species_standardized <- std_taxon(vis_sp)
      pred_col <- if ("vision_predicted_label" %in% names(vis_df)) {
        vis_df$vision_predicted_label
      } else if ("predicted_label" %in% names(vis_df)) {
        vis_df$predicted_label
      } else {
        vis_df$species_standardized
      }
      vis_df$vision_predicted_label <- std_taxon(pred_col)
    }
  }

  list(morph_df = morph_df, vis_df = vis_df)
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

  det_tier <- if ("determiner_tier" %in% names(morph_df)) {
    as.character(morph_df$determiner_tier)
  } else {
    rep("Tier_1_Gold", nrow(morph_df))
  }
  is_gs <- det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")

  taxa <- sort(unique(morph_df$species_standardized))
  anchor_mask <- rep(FALSE, nrow(morph_df))
  n_gs_vec <- integer(length(taxa))
  names(n_gs_vec) <- taxa

  for (tx in taxa) {
    tx_idx <- which(morph_df$species_standardized == tx)
    gs_pure_idx <- tx_idx[is_t1_pure[tx_idx] & is_gs[tx_idx]]
    n_gs_vec[tx] <- length(gs_pure_idx)
    if (length(gs_pure_idx) >= 5L) {
      anchor_mask[gs_pure_idx] <- TRUE
    } else {
      pure_idx <- tx_idx[is_t1_pure[tx_idx]]
      if (length(pure_idx) >= 3L) {
        anchor_mask[pure_idx] <- TRUE
      } else {
        anchor_mask[tx_idx] <- TRUE
      }
    }
  }

  list(
    anchor_df = morph_df[anchor_mask, , drop = FALSE],
    anchor_mask = anchor_mask,
    n_gold_silver = n_gs_vec
  )
}

# ------------------------------------------------------------------------------
# 3. Step B: Morphometric Overlap Matrix (S^morph_ij)
# ------------------------------------------------------------------------------
compute_morphometric_affinity <- function(morph_df, anchor_df, taxa) {
  n_taxa <- length(taxa)
  cand_pcs <- c(
    paste0("Sym_PC", 1:5),
    paste0("PC", 1:5),
    c("can1", "can2", "can3")
  )
  sym_cols <- intersect(paste0("Sym_PC", 1:5), names(anchor_df))
  pc_cols <- intersect(paste0("PC", 1:5), names(anchor_df))
  can_cols <- intersect(c("can1", "can2", "can3"), names(anchor_df))

  feat_cols <- if (length(sym_cols) >= 2L) {
    sym_cols
  } else if (length(pc_cols) >= 2L) {
    pc_cols
  } else if (length(can_cols) >= 2L) {
    can_cols
  } else {
    intersect(cand_pcs, names(anchor_df))
  }

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

  X_anchor <- as.matrix(anchor_df[, feat_cols, drop = FALSE])
  storage.mode(X_anchor) <- "double"
  X_all <- as.matrix(morph_df[, feat_cols, drop = FALSE])
  storage.mode(X_all) <- "double"

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
      cov_tx <- stats::cov(sub_x)
      cov_tx[!is.finite(cov_tx)] <- 0.0
      S_pooled <- S_pooled + (nrow(sub_x) - 1L) * cov_tx
      df_total <- df_total + (nrow(sub_x) - 1L)
    }
  }

  if (df_total > 0L) {
    S_pooled <- S_pooled / df_total
  } else {
    S_pooled <- stats::cov(X_anchor)
  }
  S_pooled[!is.finite(S_pooled)] <- 0.0
  diag_scale <- max(mean(diag(S_pooled), na.rm = TRUE), 1e-4)
  S_reg <- S_pooled + diag(max(1e-4, 0.02 * diag_scale), nrow = p, ncol = p)
  S_inv <- tryCatch(solve(S_reg), error = function(e) MASS::ginv(S_reg))

  D2_mahal <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) {
    for (j in seq_len(n_taxa)) {
      if (i != j) {
        diff_ij <- centroids[i, ] - centroids[j, ]
        d2 <- as.numeric(t(diff_ij) %*% S_inv %*% diff_ij)
        D2_mahal[i, j] <- max(0.0, d2)
      }
    }
  }

  off_diag_d2 <- D2_mahal[upper.tri(D2_mahal)]
  med_d2 <- stats::median(off_diag_d2[off_diag_d2 > 0], na.rm = TRUE)
  if (!is.finite(med_d2) || med_d2 <= 1e-8) med_d2 <- 1.0

  K_mahal <- exp(-D2_mahal / (2.0 * med_d2))
  diag(K_mahal) <- 0.0

  # Empirical CDA cross-assignment matrix C^morph_ij
  C_morph <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]
    idx_i <- which(morph_df$species_standardized == tx_i)
    if (length(idx_i) == 0L) next
    X_i <- X_all[idx_i, , drop = FALSE]
    cda_pred_i <- if ("cda_predicted_taxon" %in% names(morph_df)) {
      as.character(morph_df$cda_predicted_taxon[idx_i])
    } else {
      rep(tx_i, length(idx_i))
    }

    # Compute Gaussian Mahalanobis posteriors P(taxon j | x_m) across taxa
    log_likes <- matrix(0.0, nrow = length(idx_i), ncol = n_taxa)
    for (j in seq_len(n_taxa)) {
      diff_mat <- sweep(X_i, 2, centroids[j, ], "-")
      d2_vec <- rowSums((diff_mat %*% S_inv) * diff_mat)
      log_likes[, j] <- -0.5 * pmax(0.0, d2_vec)
    }
    max_ll <- apply(log_likes, 1, max)
    probs_i <- exp(sweep(log_likes, 1, max_ll, "-"))
    probs_i <- sweep(probs_i, 1, pmax(rowSums(probs_i), 1e-12), "/")

    for (j in seq_len(n_taxa)) {
      if (i == j) next
      tx_j <- taxa[j]
      is_cross <- (!is.na(cda_pred_i) & cda_pred_i == tx_j) | (probs_i[, j] >= 0.15)
      C_morph[i, j] <- mean(is_cross, na.rm = TRUE)
    }
  }

  C_morph_sym <- 0.5 * (C_morph + t(C_morph))
  diag(C_morph_sym) <- 0.0

  S_morph <- 0.60 * K_mahal + 0.40 * C_morph_sym
  diag(S_morph) <- 0.0

  list(
    S_morph = S_morph,
    K_mahal = K_mahal,
    C_morph_sym = C_morph_sym,
    D2_mahal = D2_mahal
  )
}

# ------------------------------------------------------------------------------
# 4. Step C: Cleanlab / DINOv2 Visual Confusion Matrix (S^vis_ij)
# ------------------------------------------------------------------------------
compute_visual_affinity <- function(vis_df, taxa) {
  n_taxa <- length(taxa)
  C_vis <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))

  if (is.null(vis_df) || nrow(vis_df) == 0L) {
    return(list(S_vis = C_vis, C_vis = C_vis))
  }

  c_err <- if ("c_error" %in% names(vis_df)) {
    suppressWarnings(as.numeric(vis_df$c_error))
  } else {
    rep(0.5, nrow(vis_df))
  }
  c_err[!is.finite(c_err)] <- 0.5

  det_tier <- if ("determiner_tier" %in% names(vis_df)) {
    as.character(vis_df$determiner_tier)
  } else {
    rep("Tier_1_Gold", nrow(vis_df))
  }
  is_anchor <- det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")

  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]
    idx_i <- which(vis_df$species_standardized == tx_i)
    if (length(idx_i) == 0L) next
    pred_i <- as.character(vis_df$vision_predicted_label[idx_i])
    cerr_i <- c_err[idx_i]
    anch_i <- is_anchor[idx_i]

    # Weight by (1 - c_error) on anchors or c_error on cross-assigned sheets
    w_i <- ifelse(
      !is.na(pred_i) & pred_i != tx_i,
      ifelse(anch_i, pmax(1.0 - cerr_i, cerr_i), pmax(cerr_i, 0.25)),
      1.0 - cerr_i
    )
    w_i[!is.finite(w_i) | w_i <= 0] <- 0.5

    for (j in seq_len(n_taxa)) {
      if (i == j) next
      tx_j <- taxa[j]
      match_j <- !is.na(pred_i) & (pred_i == tx_j)
      if (any(match_j)) {
        C_vis[i, j] <- sum(w_i[match_j], na.rm = TRUE) / length(idx_i)
      }
    }
  }

  S_vis <- 0.5 * (C_vis + t(C_vis))
  diag(S_vis) <- 0.0
  list(S_vis = S_vis, C_vis = C_vis)
}

# ------------------------------------------------------------------------------
# 5. Step D: Consensus Affinity Matrix & Unsupervised Graph / Ward Partitioning
# ------------------------------------------------------------------------------
relabel_clusters_sequential <- function(cl_vec, taxa) {
  u_ord <- unique(cl_vec[order(taxa)])
  match(cl_vec, u_ord)
}

partition_taxon_complexes <- function(W, taxa, n_vouchers_vec, n_gs_vec) {
  n_taxa <- length(taxa)
  diag(W) <- 0.0
  D_mat <- 1.0 - pmin(pmax(W, 0.0), 1.0)
  diag(D_mat) <- 0.0
  dist_obj <- stats::as.dist(D_mat)

  # Method 1: Graph Modularity Community Detection (Leiden / Louvain)
  off_vals <- W[upper.tri(W)]
  tau <- min(0.03, stats::quantile(off_vals[off_vals > 0], probs = 0.50, na.rm = TRUE))
  if (!is.finite(tau) || tau <= 0) tau <- 0.01
  W_adj <- W
  W_adj[W_adj < tau] <- 0.0
  # Ensure every node keeps at least its strongest edge so isolated singletons don't fragment needlessly
  for (i in seq_len(n_taxa)) {
    if (sum(W_adj[i, ] > 0) == 0L && max(W[i, ]) > 0) {
      best_j <- which.max(W[i, ])
      W_adj[i, best_j] <- W[i, best_j]
      W_adj[best_j, i] <- W[best_j, i]
    }
  }

  leiden_comm <- rep(1L, n_taxa)
  best_mod <- -Inf
  if (requireNamespace("igraph", quietly = TRUE) && sum(W_adj > 0) > 0) {
    g <- igraph::graph_from_adjacency_matrix(W_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    gammas <- c(0.6, 0.8, 1.0, 1.2)
    for (gm in gammas) {
      cl_try <- tryCatch(
        suppressWarnings(igraph::cluster_leiden(
          g,
          objective_function = "modularity",
          weights = igraph::E(g)$weight,
          resolution = gm,
          n_iterations = 10
        )),
        error = function(e) {
          tryCatch(
            suppressWarnings(igraph::cluster_leiden(
              g,
              objective_function = "modularity",
              weights = igraph::E(g)$weight,
              resolution_parameter = gm,
              n_iterations = 10
            )),
            error = function(e2) NULL
          )
        }
      )
      if (!is.null(cl_try)) {
        memb <- as.integer(igraph::membership(cl_try))
        k_m <- length(unique(memb))
        q_val <- tryCatch(igraph::modularity(g, memb, weights = igraph::E(g)$weight), error = function(e) -1.0)
        if (k_m >= 2L && k_m <= min(8L, n_taxa - 1L) && is.finite(q_val) && q_val > best_mod) {
          best_mod <- q_val
          leiden_comm <- memb
        }
      }
    }
    if (length(unique(leiden_comm)) < 2L) {
      lv_try <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$weight), error = function(e) NULL)
      if (!is.null(lv_try)) leiden_comm <- as.integer(igraph::membership(lv_try))
    }
  }
  leiden_comm <- relabel_clusters_sequential(leiden_comm, taxa)

  # Method 2: Ward.D2 Hierarchical Clustering with Silhouette & Gap Statistic
  hc <- stats::hclust(dist_obj, method = "ward.D2")
  k_max <- min(8L, n_taxa - 1L)
  ward_cluster <- rep(1L, n_taxa)
  best_sil_mean <- -Inf
  best_sil_vec <- rep(0.0, n_taxa)

  if (n_taxa >= 3L && k_max >= 2L) {
    k_candidates <- 2L:k_max
    sil_scores <- numeric(length(k_candidates))
    for (idx_k in seq_along(k_candidates)) {
      k_val <- k_candidates[idx_k]
      ct <- stats::cutree(hc, k = k_val)
      sil_obj <- cluster::silhouette(ct, dist_obj)
      sil_scores[idx_k] <- mean(sil_obj[, "sil_width"], na.rm = TRUE)
    }

    # Also compute Tibshirani's Gap Statistic via cluster::clusGap on PCoA coords
    pcoa_coords <- tryCatch(stats::cmdscale(dist_obj, k = min(3L, n_taxa - 1L)), error = function(e) D_mat)
    gap_k <- NA_integer_
    if (requireNamespace("cluster", quietly = TRUE) && nrow(pcoa_coords) >= 4L) {
      gap_fit <- tryCatch({
        fun_ward <- function(x, k) list(cluster = stats::cutree(stats::hclust(stats::dist(x), method = "ward.D2"), k = k))
        cluster::clusGap(pcoa_coords, FUNcluster = fun_ward, K.max = k_max, B = 25, verbose = FALSE)
      }, error = function(e) NULL)
      if (!is.null(gap_fit)) {
        gap_k <- tryCatch(cluster::maxSE(gap_fit$Tab[, "gap"], gap_fit$Tab[, "SE.sim"], method = "firstSEmax"), error = function(e) NA_integer_)
      }
    }

    best_idx <- which.max(sil_scores)
    best_k <- k_candidates[best_idx]
    best_sil_mean <- sil_scores[best_idx]
    ward_cluster <- as.integer(stats::cutree(hc, k = best_k))
    ward_cluster <- relabel_clusters_sequential(ward_cluster, taxa)
    sil_final <- cluster::silhouette(ward_cluster, dist_obj)
    best_sil_vec <- as.numeric(sil_final[, "sil_width"])
  }

  # Evaluate Leiden silhouette width to select consensus complex_id
  leiden_sil_mean <- -Inf
  leiden_sil_vec <- best_sil_vec
  if (length(unique(leiden_comm)) >= 2L && length(unique(leiden_comm)) < n_taxa) {
    l_sil <- cluster::silhouette(leiden_comm, dist_obj)
    leiden_sil_mean <- mean(l_sil[, "sil_width"], na.rm = TRUE)
    leiden_sil_vec <- as.numeric(l_sil[, "sil_width"])
  }

  if (leiden_sil_mean > best_sil_mean + 0.02) {
    complex_id <- leiden_comm
    sil_width_vec <- leiden_sil_vec
  } else {
    complex_id <- ward_cluster
    sil_width_vec <- best_sil_vec
  }
  complex_id <- relabel_clusters_sequential(complex_id, taxa)

  # Auto-name each complex after its two highest-degree / most abundant member species
  complex_labels <- character(n_taxa)
  deg_vec <- rowSums(W)
  for (cid in sort(unique(complex_id))) {
    m_idx <- which(complex_id == cid)
    ord <- m_idx[order(-deg_vec[m_idx], -n_vouchers_vec[m_idx], taxa[m_idx])]
    top_sp <- taxa[ord[seq_len(min(2L, length(ord)))]]
    top_ep <- sub("^Packera\\s+", "", top_sp)
    lbl <- paste0("Complex_", cid, "_", paste(top_ep, collapse = "_"))
    complex_labels[m_idx] <- lbl
  }

  # Top 2 empirical neighbors per taxon
  top_n1 <- character(n_taxa)
  top_a1 <- numeric(n_taxa)
  top_n2 <- character(n_taxa)
  top_a2 <- numeric(n_taxa)
  for (i in seq_len(n_taxa)) {
    w_row <- W[i, ]
    w_row[i] <- -1.0
    ord_n <- order(w_row, decreasing = TRUE)
    top_n1[i] <- if (length(ord_n) >= 1L) taxa[ord_n[1L]] else ""
    top_a1[i] <- if (length(ord_n) >= 1L) round(max(0.0, W[i, ord_n[1L]]), 4) else 0.0
    top_n2[i] <- if (length(ord_n) >= 2L) taxa[ord_n[2L]] else ""
    top_a2[i] <- if (length(ord_n) >= 2L) round(max(0.0, W[i, ord_n[2L]]), 4) else 0.0
  }

  within_mean_aff <- numeric(n_taxa)
  for (i in seq_len(n_taxa)) {
    peers_i <- which(complex_id == complex_id[i] & seq_len(n_taxa) != i)
    within_mean_aff[i] <- if (length(peers_i) > 0L) round(mean(W[i, peers_i], na.rm = TRUE), 4) else 0.0
  }

  complexes_df <- data.frame(
    species_standardized         = taxa,
    taxon                        = taxa,
    n_vouchers                   = as.integer(n_vouchers_vec[taxa]),
    n_gold_silver_anchors        = as.integer(n_gs_vec[taxa]),
    n_Pure_Gold_Silver           = as.integer(n_gs_vec[taxa]),
    complex_id                   = as.integer(complex_id),
    consensus_complex_id         = as.integer(complex_id),
    complex_label                = complex_labels,
    leiden_community             = as.integer(leiden_comm),
    leiden_complex_id            = as.integer(leiden_comm),
    ward_silhouette_cluster      = as.integer(ward_cluster),
    ward_complex_id              = as.integer(ward_cluster),
    silhouette_width             = round(sil_width_vec, 4),
    within_complex_mean_affinity = within_mean_aff,
    top_neighbor_1               = top_n1,
    nearest_neighbor_taxon       = top_n1,
    top_neighbor_1_affinity      = top_a1,
    nearest_neighbor_affinity    = top_a1,
    top_neighbor_2               = top_n2,
    top_neighbor_2_affinity      = top_a2,
    stringsAsFactors             = FALSE
  )

  list(
    complexes_df = complexes_df,
    hc = hc,
    dist_obj = dist_obj
  )
}

# ------------------------------------------------------------------------------
# 6. Step F: 2-Panel Publication Figure (Ward.D2 + Heatmap & 2D PCoA Network)
# ------------------------------------------------------------------------------
render_complex_discovery_figure <- function(W, part_res, morph_res, vis_res, out_pdf) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  cdf <- part_res$complexes_df
  taxa <- cdf$species_standardized
  n_taxa <- length(taxa)
  short_names <- sub("^Packera\\s+", "P. ", taxa)

  dyn_pal <- if (exists("get_dynamic_packera_palette", mode = "function")) {
    get_dynamic_packera_palette(taxa)
  } else {
    stats::setNames(grDevices::hcl.colors(n_taxa, "Dark 3"), taxa)
  }

  complex_pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#882255", "#44AA99")
  hollow_pchs <- c(1L, 2L, 5L, 0L, 6L, 4L, 3L, 8L)

  draw_discovery_panels <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(6.2, 6.2, 3.4, 1.5), family = "sans")

    # Panel A: N x N Empirical Taxon Affinity Matrix Heatmap ordered by Ward.D2 Dendrogram
    ord <- part_res$hc$order
    W_ord <- W[ord, ord, drop = FALSE]
    diag(W_ord) <- 1.0
    lbl_ord <- short_names[ord]
    cid_ord <- cdf$complex_id[ord]

    heat_cols <- grDevices::hcl.colors(64, "YlOrRd", rev = TRUE)
    graphics::image(
      x = seq_len(n_taxa),
      y = seq_len(n_taxa),
      z = t(W_ord[rev(seq_len(n_taxa)), , drop = FALSE]),
      col = heat_cols,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = "A. Empirical Taxon Affinity Matrix (W_ij) & Complexes"
    )
    graphics::axis(1, at = seq_len(n_taxa), labels = lbl_ord, las = 2, cex.axis = 0.78, font = 3)
    graphics::axis(2, at = seq_len(n_taxa), labels = rev(lbl_ord), las = 2, cex.axis = 0.78, font = 3)
    graphics::box()

    # Draw dashed bounding boxes around discovered species complexes along diagonal
    u_cids <- unique(cid_ord)
    for (cid in u_cids) {
      idx_c <- which(cid_ord == cid)
      if (length(idx_c) >= 1L) {
        x0 <- min(idx_c) - 0.5
        x1 <- max(idx_c) + 0.5
        y0 <- (n_taxa - max(idx_c) + 1L) - 0.5
        y1 <- (n_taxa - min(idx_c) + 1L) + 0.5
        b_col <- complex_pal[(cid - 1L) %% length(complex_pal) + 1L]
        graphics::rect(x0, y0, x1, y1, border = b_col, lwd = 2.4, lty = 2)
      }
    }

    # Panel B: 2D PCoA Taxon Affinity Network with Unfilled/Hollow Markers & Convex Hulls
    pcoa_2d <- tryCatch(stats::cmdscale(part_res$dist_obj, k = 2), error = function(e) NULL)
    if (is.null(pcoa_2d) || ncol(pcoa_2d) < 2L) {
      ang <- seq(0, 2 * pi, length.out = n_taxa + 1L)[seq_len(n_taxa)]
      pcoa_2d <- cbind(cos(ang), sin(ang))
    }
    px <- as.numeric(pcoa_2d[, 1])
    py <- as.numeric(pcoa_2d[, 2])
    xr <- range(px, na.rm = TRUE); xpad <- max(0.18 * diff(xr), 0.15)
    yr <- range(py, na.rm = TRUE); ypad <- max(0.18 * diff(yr), 0.15)

    graphics::plot(
      px, py,
      type = "n",
      xlim = c(xr[1] - xpad, xr[2] + xpad),
      ylim = c(yr[1] - ypad, yr[2] + ypad),
      xlab = "PCoA Axis 1 (1 - W_ij Dissimilarity)",
      ylab = "PCoA Axis 2 (1 - W_ij Dissimilarity)",
      main = "B. Unsupervised Taxon-Complex Affinity Network",
      font.main = 2
    )
    graphics::grid(col = "grey88", lty = "solid")
    graphics::abline(h = 0, v = 0, lty = 3, col = "grey72")

    # Draw dashed convex hulls around each discovered species complex
    for (cid in sort(unique(cdf$complex_id))) {
      m <- which(cdf$complex_id == cid)
      h_col <- complex_pal[(cid - 1L) %% length(complex_pal) + 1L]
      if (length(m) >= 3L) {
        hull_idx <- grDevices::chull(px[m], py[m])
        hull_pts <- m[c(hull_idx, hull_idx[1L])]
        graphics::polygon(px[hull_pts], py[hull_pts], border = h_col, col = NA, lty = 2, lwd = 1.8)
      } else if (length(m) == 2L) {
        graphics::segments(px[m[1]], py[m[1]], px[m[2]], py[m[2]], col = h_col, lty = 2, lwd = 2.0)
      }
    }

    # Draw weighted affinity edges W_ij >= 0.05
    for (i in seq_len(n_taxa - 1L)) {
      for (j in (i + 1L):n_taxa) {
        if (W[i, j] >= 0.05) {
          same_c <- cdf$complex_id[i] == cdf$complex_id[j]
          e_col <- if (same_c) {
            grDevices::adjustcolor(complex_pal[(cdf$complex_id[i] - 1L) %% length(complex_pal) + 1L], alpha.f = 0.70)
          } else {
            grDevices::adjustcolor("grey55", alpha.f = 0.45)
          }
          graphics::segments(px[i], py[i], px[j], py[j], col = e_col, lwd = max(0.8, 5.5 * W[i, j]))
        }
      }
    }

    # Plot strictly unfilled/hollow nodes (pch in 1, 2, 5, 0, 6, 4; bg = NA)
    node_cex <- pmax(1.1, pmin(2.6, 0.9 + 0.65 * log10(pmax(cdf$n_vouchers, 1))))
    node_pch <- hollow_pchs[(cdf$complex_id - 1L) %% length(hollow_pchs) + 1L]
    node_col <- dyn_pal[taxa]
    node_col[is.na(node_col)] <- complex_pal[(cdf$complex_id[is.na(node_col)] - 1L) %% length(complex_pal) + 1L]

    graphics::points(
      px, py,
      pch = node_pch,
      bg = NA,
      col = node_col,
      cex = node_cex,
      lwd = 1.8
    )
    graphics::text(px, py, labels = short_names, pos = 3, cex = 0.76, font = 3, col = "black")

    u_cl <- sort(unique(cdf$complex_id))
    cl_names <- vapply(u_cl, function(k) cdf$complex_label[ which(cdf$complex_id == k)[1] ], character(1))
    graphics::legend(
      "topright",
      legend = cl_names,
      col = complex_pal[(u_cl - 1L) %% length(complex_pal) + 1L],
      pch = hollow_pchs[(u_cl - 1L) %% length(hollow_pchs) + 1L],
      pt.bg = NA,
      lty = 2,
      lwd = 1.6,
      cex = 0.72,
      bty = "n"
    )
  }

  pdf_tmp <- paste0(out_pdf, ".", Sys.getpid(), ".tmp")
  grDevices::pdf(pdf_tmp, width = 13, height = 6.2, useDingbats = FALSE)
  tryCatch(draw_discovery_panels(), finally = grDevices::dev.off())
  if (file.exists(pdf_tmp)) file.rename(pdf_tmp, out_pdf)

  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  if (out_png != out_pdf) {
    png_tmp <- paste0(out_png, ".", Sys.getpid(), ".tmp")
    grDevices::png(png_tmp, width = 13, height = 6.2, units = "in", res = 300)
    tryCatch(draw_discovery_panels(), finally = grDevices::dev.off())
    if (file.exists(png_tmp)) file.rename(png_tmp, out_png)
  }
  invisible(out_pdf)
}

# ------------------------------------------------------------------------------
# 7. Main Execution Entrypoint
# ------------------------------------------------------------------------------
main <- function() {
  set.seed(42)
  opts <- parse_discovery_args()
  message("=== Phase 08: Unsupervised Quantitative Taxon-Complex Discovery Engine ===")

  inputs <- load_and_merge_inputs(opts)
  morph_df <- inputs$morph_df
  vis_df <- inputs$vis_df

  taxa <- sort(unique(morph_df$species_standardized))
  n_taxa <- length(taxa)
  if (n_taxa < 2L) {
    stop("At least 2 distinct Packera taxa are required for taxon-complex discovery (found ", n_taxa, ").")
  }

  n_vouchers_vec <- table(morph_df$species_standardized)[taxa]
  anc <- select_taxon_anchors(morph_df)
  message(sprintf("Evaluated %d vouchers across %d Packera taxa (%d Tier_1_Pure Gold/Silver anchors).",
                  nrow(morph_df), n_taxa, nrow(anc$anchor_df)))

  morph_res <- compute_morphometric_affinity(morph_df, anc$anchor_df, taxa)
  vis_res <- compute_visual_affinity(vis_df, taxa)

  if (!is.null(vis_df) && nrow(vis_df) > 0L && any(vis_res$S_vis > 0)) {
    W <- 0.55 * morph_res$S_morph + 0.45 * vis_res$S_vis
  } else {
    W <- morph_res$S_morph
  }
  diag(W) <- 0.0
  W <- round(W, 6)

  part_res <- partition_taxon_complexes(
    W = W,
    taxa = taxa,
    n_vouchers_vec = n_vouchers_vec,
    n_gs_vec = anc$n_gold_silver
  )

  # Export empirical_taxon_complexes.csv and pairwise_taxon_affinity_matrix.csv atomically
  atomic_write_csv(part_res$complexes_df, opts$output_complexes)
  message("Saved empirical taxon complexes -> ", opts$output_complexes)

  aff_df <- data.frame(
    taxon = taxa,
    species_standardized = taxa,
    as.data.frame(W, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  atomic_write_csv(aff_df, opts$output_affinity)
  message("Saved pairwise taxon affinity matrix -> ", opts$output_affinity)

  # Step E: Ego-centric --auto-neighbors selection for focal taxon runs
  focal_raw <- trimws(as.character(opts$focal_taxon))
  if (nzchar(focal_raw)) {
    focal_std <- std_taxon(focal_raw)
    if (focal_std %in% taxa) {
      cdf <- part_res$complexes_df
      focal_cid <- cdf$complex_id[cdf$species_standardized == focal_std][1L]
      same_complex <- cdf$species_standardized[cdf$complex_id == focal_cid]

      cm_focal <- morph_res$C_morph_sym[focal_std, ]
      cv_focal <- vis_res$S_vis[focal_std, ]
      aff_neighbors <- taxa[
        (cm_focal >= opts$min_affinity) |
        (cv_focal >= opts$min_affinity)
      ]
      selected_neighbors <- sort(setdiff(union(same_complex, aff_neighbors), focal_std))
    } else {
      selected_neighbors <- character(0)
    }
    atomic_write_lines(selected_neighbors, opts$output_neighbors)
    message("Saved focal auto-neighbors (", length(selected_neighbors), " taxa) -> ", opts$output_neighbors)
    cat(paste(selected_neighbors, collapse = "\n"), "\n", sep = "")
  } else {
    # Always write default neighbors file listing discovered complex members of top taxon
    cdf <- part_res$complexes_df
    atomic_write_lines(cdf$species_standardized, opts$output_neighbors)
  }

  # Step F: Render 2-panel publication figure
  render_complex_discovery_figure(W, part_res, morph_res, vis_res, opts$output_figure)
  message("Saved taxon-complex discovery figure -> ", opts$output_figure)
}

if (sys.nframe() == 0L) {
  main()
}
