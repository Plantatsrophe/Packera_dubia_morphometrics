#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/authority_weighting_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Unified Empirical Bayes / Dawid-Skene Hierarchical Authority-Tier Weighting
#   and Adaptive Reference Anchor Promotion utilities (Slenker et al. 2022;
#   Dawid & Skene 1979; Gelman et al. 2013; Kaufman & Rousseeuw 1990).
#
#   1. select_adaptive_reference_anchors():
#      - Identifies clean Tier_1_Pure specimens across all active Packera taxa.
#      - Prunes within-taxon Gold + Silver mislabeled outliers via robust MCD /
#        coordinate-wise medoid Mahalanobis distance (chi-sq_0.95,p).
#      - Promotes clean Mahalanobis-screened Tier_3_Bronze vouchers (within
#        chi-sq_0.80,p of the taxon medoid) whenever a species has fewer than
#        min_anchors_per_taxon clean Gold + Silver anchors, preventing rare or
#        under-annotated species from failing PCA, allometry, GMM, or CDA.
#
#   2. compute_empirical_bayes_authority_weights():
#      - Harmonizes authority weights via Beta-Binomial Empirical Bayes
#        shrinkage over hierarchical tier reliability priors (Gold = 0.95,
#        Silver = 0.85, promoted Bronze = 0.65, unpromoted Bronze = 0.35).
#      - Supports both strict "active_passive" projection weighting (active
#        mean = 1.0, unpromoted Bronze = 0.0) and "continuous_prior_posterior"
#        case weighting normalized to [0.5, 3.0].
# ==============================================================================

.standardize_determiner_tier_vec <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  raw <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_3_Bronze", n)
  raw[is.na(raw) | !nzchar(trimws(raw))] <- "Tier_3_Bronze"
  out <- rep("Tier_3_Bronze", n)
  out[grepl("tier_?1|gold", raw, ignore.case = TRUE)] <- "Tier_1_Gold"
  out[grepl("tier_?2|silver", raw, ignore.case = TRUE)] <- "Tier_2_Silver"
  out
}

.resolve_authority_taxon_vec <- function(df, taxon_col = "species_standardized") {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  cand_cols <- unique(c(taxon_col, "species_standardized", "grouping_col", "scientificName", "species_raw"))
  chosen_col <- NULL
  for (cc in cand_cols) {
    if (!is.null(cc) && nzchar(cc) && (cc %in% names(df))) {
      vals <- trimws(as.character(df[[cc]]))
      if (any(!is.na(vals) & nzchar(vals))) {
        chosen_col <- cc
        break
      }
    }
  }
  if (is.null(chosen_col)) return(rep("Unassigned_Packera", n))
  tx <- trimws(as.character(df[[chosen_col]]))
  if (chosen_col %in% c("scientificName", "species_raw") && !any(grepl("^Morphotype_K", tx[!is.na(tx)]))) {
    if (exists("standardize_packera_taxon", mode = "function")) {
      tx <- as.character(standardize_packera_taxon(tx))
    } else if (exists("standardize_taxon_name", mode = "function")) {
      std_tx <- as.character(standardize_taxon_name(tx))
      tx <- ifelse(!is.na(std_tx) & nzchar(std_tx), std_tx, tx)
    }
  }
  tx[is.na(tx) | !nzchar(tx) | tolower(tx) == "unknown"] <- "Unassigned_Packera"
  tx
}

.identify_tier1_pure_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  has_mts <- "morph_tier_source" %in% names(df)
  has_nt1 <- "n_tier1_leaves" %in% names(df)
  if (has_mts) {
    mts <- trimws(as.character(df$morph_tier_source))
    is_t1 <- !is.na(mts) & (mts == "Tier_1_Pure")
    if (has_nt1) {
      nt1 <- suppressWarnings(as.numeric(df$n_tier1_leaves))
      is_t1 <- is_t1 | (is.na(mts) & !is.na(nt1) & nt1 >= 1)
    }
  } else if (has_nt1) {
    nt1 <- suppressWarnings(as.numeric(df$n_tier1_leaves))
    is_t1 <- !is.na(nt1) & (nt1 >= 1)
  } else if ("tier_assignment" %in% names(df)) {
    is_t1 <- grepl("^tier_?1", as.character(df$tier_assignment), ignore.case = TRUE)
  } else if ("assigned_tier" %in% names(df)) {
    is_t1 <- grepl("^tier_?1", as.character(df$assigned_tier), ignore.case = TRUE)
  } else {
    is_t1 <- rep(TRUE, n)
  }
  is_t1[is.na(is_t1)] <- FALSE
  if ("is_contour_outlier" %in% names(df)) {
    out_flag <- as.logical(df$is_contour_outlier)
    out_flag[is.na(out_flag)] <- FALSE
    if (sum(is_t1 & !out_flag) >= 2L) is_t1 <- is_t1 & !out_flag
  }
  is_t1
}

.discover_numeric_sorted_cols <- function(nms, pattern) {
  hits <- grep(pattern, nms, value = TRUE)
  if (length(hits) <= 1L) return(hits)
  nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", hits)))
  hits[order(nums, hits)]
}

select_morphometric_subspace <- function(X, group_vec = NULL, min_p = 2L, max_p = NULL) {
  X <- as.matrix(X)
  M <- ncol(X); N <- nrow(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(M))
  if (M <= min_p) return(colnames(X))
  pc_idx <- grep("^(Sym_)?PC[0-9]+$", colnames(X))
  if (length(pc_idx) >= 3L) {
    v_obs <- apply(X[, pc_idx, drop = FALSE], 2, stats::var, na.rm = TRUE)
    v_prop <- v_obs / max(sum(v_obs, na.rm = TRUE), 1e-12)
    k_m <- length(pc_idx)
    bs_prop <- vapply(seq_len(k_m), function(j) sum(1 / (j:k_m)) / k_m, numeric(1))
    horn_thr <- (1 + sqrt(k_m / max(N, k_m + 1)))^2 / k_m
    keep_pc <- which(is.finite(v_prop) & (v_prop >= pmin(bs_prop * 0.55, horn_thr * 0.55)))
    if (length(keep_pc) < min_p) keep_pc <- seq_len(min(min_p, k_m))
    X <- X[, c(pc_idx[keep_pc], setdiff(seq_len(M), pc_idx)), drop = FALSE]
  }
  valid_g <- if (!is.null(group_vec)) stats::na.omit(as.character(group_vec)) else character(0)
  valid_g <- valid_g[nzchar(valid_g) & !(valid_g %in% c("Unassigned_Packera", "Unknown"))]
  g_tab <- table(valid_g); g_tab_multi <- g_tab[g_tab >= 2L]
  N_min <- if (length(g_tab_multi) >= 2L) min(as.integer(g_tab_multi)) else N
  P_max <- if (N_min >= 6L) max(min_p, floor(N_min / 3)) else ncol(X)
  if (!is.null(max_p) && length(max_p) == 1L && is.finite(as.numeric(max_p)) && as.numeric(max_p) >= 1) {
    P_max <- min(P_max, as.integer(max_p))
  }
  colnames(X)[seq_len(min(ncol(X), P_max))]
}

.prepare_anchor_feature_matrix <- function(df, feature_matrix = NULL, max_p = NULL, group_vec = NULL) {
  n <- nrow(df)
  if (is.null(feature_matrix)) {
    auto_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
    if (length(auto_cols) == 0L) auto_cols <- .discover_numeric_sorted_cols(names(df), "^Sym_PC[0-9]+$")
    if (length(auto_cols) == 0L) auto_cols <- .discover_numeric_sorted_cols(names(df), "^[AD][0-9]+$")
    if (length(auto_cols) == 0L) return(NULL)
    feature_matrix <- df[, auto_cols, drop = FALSE]
  }
  X <- as.matrix(feature_matrix)
  storage.mode(X) <- "double"
  if (nrow(X) != n || ncol(X) == 0L) return(NULL)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))

  primary_cols <- grep("^(Sym_)?PC[0-9]+$|^[AD]_?[0-9]+$", colnames(X))
  if (length(primary_cols) == 0L) primary_cols <- seq_len(ncol(X))
  has_feat <- apply(X[, primary_cols, drop = FALSE], 1, function(r) all(is.finite(r)))
  if (sum(has_feat) < 2L) {
    has_feat <- apply(X, 1, function(r) any(is.finite(r)))
  }

  keep_cols <- which(vapply(seq_len(ncol(X)), function(j) {
    v <- X[has_feat, j]
    v <- v[is.finite(v)]
    length(v) >= 3L && stats::var(v) > 1e-10
  }, logical(1)))
  if (length(keep_cols) == 0L) return(list(X = NULL, has_feat = has_feat))

  X_use <- X[, keep_cols, drop = FALSE]
  for (j in seq_len(ncol(X_use))) {
    bad_j <- !is.finite(X_use[, j])
    if (any(bad_j)) {
      med_j <- stats::median(X_use[!bad_j, j], na.rm = TRUE)
      X_use[bad_j, j] <- if (is.finite(med_j)) med_j else 0.0
    }
  }
  if (!any(grepl("^(Sym_)?PC[0-9]+$", colnames(X_use))) && all(grepl("^[ABCD]_?[0-9]+$", colnames(X_use)))) {
    h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(X_use))))
    if (!any(is.na(h_nums))) X_use <- X_use[, order(h_nums, colnames(X_use)), drop = FALSE]
    raw_cap <- if (!is.null(max_p) && length(max_p) == 1L && is.finite(as.numeric(max_p)) && as.numeric(max_p) >= 1) {
      as.integer(max_p)
    } else {
      max(2L, min(ncol(X_use), floor(max(sum(has_feat), 4L) / 2L), 10L))
    }
    X_use <- X_use[, seq_len(min(ncol(X_use), raw_cap)), drop = FALSE]
  } else {
    sel_cols <- select_morphometric_subspace(
      X_use[has_feat, , drop = FALSE],
      group_vec = if (!is.null(group_vec)) group_vec[has_feat] else NULL,
      min_p = min(2L, ncol(X_use)),
      max_p = max_p
    )
    X_use <- X_use[, sel_cols, drop = FALSE]
  }
  list(X = X_use, has_feat = has_feat)
}

.compute_robust_taxon_medoid <- function(X_sub, S_inv = NULL) {
  X_sub <- as.matrix(X_sub)
  m <- nrow(X_sub)
  p <- ncol(X_sub)
  if (m == 0L) return(rep(0.0, p))
  if (m == 1L) return(as.numeric(X_sub[1, ]))
  if (m <= 4L) return(as.numeric(apply(X_sub, 2, stats::median, na.rm = TRUE)))
  if (is.null(S_inv)) {
    mads <- apply(X_sub, 2, stats::mad, na.rm = TRUE)
    sds <- apply(X_sub, 2, stats::sd, na.rm = TRUE)
    sc <- ifelse(is.finite(mads) & mads > 1e-8, mads, ifelse(is.finite(sds) & sds > 1e-8, sds, 1.0))
    Z <- sweep(X_sub, 2, sc, "/")
    dmat <- as.matrix(stats::dist(Z, method = "euclidean"))
  } else {
    # Pairwise Mahalanobis distances under precision matrix S_inv
    G <- X_sub %*% S_inv %*% t(X_sub)
    dg <- diag(G)
    d2mat <- outer(dg, dg, "+") - 2.0 * G
    d2mat[d2mat < 0.0] <- 0.0
    dmat <- as.matrix(sqrt(d2mat))
  }
  row_d <- rowSums(dmat, na.rm = TRUE)
  k_core <- max(2L, ceiling(0.60 * m))
  core_idx <- order(row_d)[seq_len(min(m, k_core))]
  as.numeric(colMeans(X_sub[core_idx, , drop = FALSE]))
}

.mcd_consistency_factor <- function(alpha, p) {
  alpha <- min(0.99, max(0.50, as.numeric(alpha)))
  denom <- stats::pchisq(stats::qchisq(alpha, df = p), df = p + 2)
  if (is.finite(denom) && denom > 1e-8) (alpha / denom) else 1.0
}

.estimate_pooled_within_cov <- function(X_use, taxa_vec, det_tier, valid_pool) {
  p <- ncol(X_use)
  resid_blocks <- list()
  for (tx in unique(taxa_vec[valid_pool])) {
    if (tx == "Unassigned_Packera") next
    idx_gs <- which(valid_pool & (taxa_vec == tx) & (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")))
    idx_tx <- if (length(idx_gs) >= 3L) idx_gs else which(valid_pool & (taxa_vec == tx))
    if (length(idx_tx) >= 3L) {
      X_t <- X_use[idx_tx, , drop = FALSE]
      med_t <- .compute_robust_taxon_medoid(X_t)
      diff_t <- sweep(X_t, 2, med_t, "-")
      nz_mask <- rowSums(diff_t^2) > 1e-12
      if (sum(nz_mask) >= 2L) diff_t <- diff_t[nz_mask, , drop = FALSE]
      if (length(idx_gs) < 3L && nrow(diff_t) >= 5L) {
        r_norm <- sqrt(rowSums(diff_t^2))
        keep_r <- r_norm <= stats::quantile(r_norm, probs = 0.75, na.rm = TRUE)
        if (sum(keep_r) >= 3L) diff_t <- diff_t[keep_r, , drop = FALSE]
      }
      df_scale <- sqrt(nrow(X_t) / max(1.0, nrow(X_t) - 1.0))
      resid_blocks[[length(resid_blocks) + 1L]] <- diff_t * df_scale
    }
  }
  if (length(resid_blocks) > 0L) {
    R_all <- do.call(rbind, resid_blocks)
    if (nrow(R_all) >= (p + 2L)) {
      if (nrow(R_all) >= 25L && requireNamespace("MASS", quietly = TRUE)) {
        q_used <- min(nrow(R_all) - 1L, max(floor(0.80 * nrow(R_all)), floor((nrow(R_all) + p + 1L) / 2L)))
        rob <- tryCatch(MASS::cov.rob(R_all, method = "mcd", quantile.used = q_used), error = function(e) NULL)
        if (!is.null(rob) && !is.null(rob$cov)) {
          c_mcd <- .mcd_consistency_factor(q_used / nrow(R_all), p)
          return(c_mcd * rob$cov + diag(1e-5, nrow = p, ncol = p))
        }
      }
      return(stats::cov(R_all) + diag(1e-5, nrow = p, ncol = p))
    }
  }
  v_diag <- apply(X_use[valid_pool, , drop = FALSE], 2, stats::var, na.rm = TRUE)
  v_diag[!is.finite(v_diag) | v_diag <= 1e-8] <- 1.0
  diag(v_diag, nrow = p, ncol = p)
}

select_adaptive_reference_anchors <- function(
  df,
  feature_matrix = NULL,
  taxon_col = "species_standardized",
  min_anchors_per_taxon = 8L,
  mahal_q_gold_silver = 0.95,
  mahal_q_bronze_promo = 0.80
) {
  if (is.null(df) || nrow(df) == 0L) {
    return(list(
      active_mask = logical(0),
      is_gold_anchor = logical(0),
      is_silver_anchor = logical(0),
      is_promoted_bronze = logical(0),
      is_promoted_anchor = logical(0),
      anchor_tier_source = character(0)
    ))
  }

  n <- nrow(df)
  min_anc <- max(1L, as.integer(min_anchors_per_taxon))
  det_tier <- .standardize_determiner_tier_vec(df)
  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  is_t1_pure <- .identify_tier1_pure_mask(df)
  if (sum(is_t1_pure) == 0L) is_t1_pure <- rep(TRUE, n)

  feat_prep <- .prepare_anchor_feature_matrix(df, feature_matrix = feature_matrix, max_p = NULL, group_vec = taxa_vec)
  X_use <- if (!is.null(feat_prep)) feat_prep$X else NULL
  has_feat <- if (!is.null(feat_prep)) feat_prep$has_feat else rep(TRUE, n)
  valid_pool <- is_t1_pure & has_feat

  is_gold <- (det_tier == "Tier_1_Gold") & valid_pool
  is_silver <- (det_tier == "Tier_2_Silver") & valid_pool
  is_gs <- is_gold | is_silver
  is_promoted_bronze <- rep(FALSE, n)

  S_pool <- if (!is.null(X_use) && sum(valid_pool) >= 4L) {
    .estimate_pooled_within_cov(X_use, taxa_vec, det_tier, valid_pool)
  } else {
    NULL
  }
  S_pool_inv <- if (!is.null(S_pool)) {
    tryCatch(solve(S_pool), error = function(e) diag(1.0 / pmax(diag(S_pool), 1e-6), ncol(S_pool)))
  } else {
    NULL
  }

  unique_taxa <- unique(taxa_vec[valid_pool])
  for (tx in unique_taxa) {
    idx_gs <- which(is_gs & (taxa_vec == tx))

    # Step 1: Within-taxon robust MCD / medoid Mahalanobis pruning of Gold + Silver (chi-sq_0.95, p)
    if (!is.null(X_use) && length(idx_gs) >= 5L) {
      p_use <- min(ncol(X_use), length(idx_gs) - 2L)
      if (p_use >= 1L) {
        X_gs <- X_use[idx_gs, seq_len(p_use), drop = FALSE]
        v_ok <- apply(X_gs, 2, stats::var, na.rm = TRUE) > 1e-10
        if (any(v_ok)) {
          X_gs <- X_gs[, v_ok, drop = FALSE]
          p_eff <- ncol(X_gs)
          rob_gs <- NULL; q_mcd <- nrow(X_gs)
          if (nrow(X_gs) >= 15L && requireNamespace("MASS", quietly = TRUE)) {
            q_mcd <- min(nrow(X_gs) - 1L, max(floor(0.80 * nrow(X_gs)), floor((nrow(X_gs) + p_eff + 1L) / 2L)))
            rob_gs <- tryCatch(MASS::cov.rob(X_gs, method = "mcd", quantile.used = q_mcd), error = function(e) NULL)
          }
          center_gs <- if (!is.null(rob_gs)) rob_gs$center else .compute_robust_taxon_medoid(X_gs)
          cov_gs <- if (!is.null(rob_gs)) (.mcd_consistency_factor(q_mcd / nrow(X_gs), p_eff) * rob_gs$cov + diag(1e-5, nrow = p_eff, ncol = p_eff)) else (stats::cov(X_gs) + diag(1e-4, nrow = p_eff, ncol = p_eff))
          d2_gs <- tryCatch(stats::mahalanobis(X_gs, center = center_gs, cov = cov_gs), error = function(e) rep(0.0, nrow(X_gs)))
          drop_gs <- is.finite(d2_gs) & (d2_gs > stats::qchisq(mahal_q_gold_silver, df = p_eff))
          if (any(drop_gs) && sum(!drop_gs) >= 3L) {
            pruned_ids <- idx_gs[drop_gs]
            is_gs[pruned_ids] <- FALSE
            is_gold[pruned_ids] <- FALSE
            is_silver[pruned_ids] <- FALSE
          }
        }
      }
    }

    # Step 2: Adaptive Anchor Promotion for Gold/Silver-poor taxa (raw GS < min_anchors_per_taxon or clean n_gs < 3)
    idx_gs_clean <- which(is_gs & (taxa_vec == tx))
    n_gs <- length(idx_gs_clean)
    if (tx != "Unassigned_Packera" && (length(idx_gs) < min_anc || n_gs < min(3L, min_anc))) {
      idx_br <- which((det_tier == "Tier_3_Bronze") & valid_pool & (taxa_vec == tx))
      n_needed <- min_anc - n_gs
      if (length(idx_br) > 0L && n_needed > 0L) {
        if (is.null(X_use) || ncol(X_use) == 0L) {
          ord_q <- if ("leaf_quality_weight" %in% names(df)) {
            order(suppressWarnings(as.numeric(df$leaf_quality_weight[idx_br])), decreasing = TRUE, na.last = TRUE)
          } else {
            seq_along(idx_br)
          }
          promo_ids <- idx_br[ord_q[seq_len(min(length(idx_br), n_needed))]]
          is_promoted_bronze[promo_ids] <- TRUE
        } else {
          p_eff <- ncol(X_use)
          X_ref <- if (n_gs >= 1L) X_use[idx_gs_clean, , drop = FALSE] else X_use[idx_br, , drop = FALSE]
          medoid_tx <- .compute_robust_taxon_medoid(X_ref, S_inv = S_pool_inv)
          base_cov_tx <- if (n_gs >= (p_eff + 3L)) {
            0.75 * stats::cov(X_use[idx_gs_clean, , drop = FALSE]) + 0.25 * S_pool + diag(1e-5, nrow = p_eff, ncol = p_eff)
          } else if (!is.null(S_pool)) {
            S_pool
          } else {
            diag(pmax(apply(X_use, 2, stats::var, na.rm = TRUE), 1e-4), nrow = p_eff, ncol = p_eff)
          }
          pred_scale <- 1.0 + (1.0 / max(1.0, nrow(X_ref)))
          cov_tx <- pred_scale * base_cov_tx
          d2_br <- tryCatch(
            stats::mahalanobis(X_use[idx_br, , drop = FALSE], center = medoid_tx, cov = cov_tx),
            error = function(e) rowSums(sweep(X_use[idx_br, , drop = FALSE], 2, medoid_tx, "-")^2)
          )
          ord_br <- order(d2_br, na.last = NA)
          q80_thr <- stats::qchisq(mahal_q_bronze_promo, df = p_eff)
          elig_ord <- ord_br[d2_br[ord_br] <= q80_thr]
          promo_ids <- idx_br[elig_ord[seq_len(min(length(elig_ord), n_needed))]]

          # Minimum-viability safeguard for 0/1-anchor rare taxa: allow up to 2 total anchors if within chi-sq_0.95
          min_viable <- min(2L, min_anc)
          if ((n_gs + length(promo_ids)) < min_viable && length(ord_br) > length(promo_ids)) {
            q95_thr <- stats::qchisq(0.95, df = p_eff)
            rem_ord <- setdiff(ord_br[d2_br[ord_br] <= q95_thr], elig_ord)
            if (length(rem_ord) > 0L) {
              n_extra <- min(length(rem_ord), min_viable - (n_gs + length(promo_ids)))
              promo_ids <- c(promo_ids, idx_br[rem_ord[seq_len(n_extra)]])
            }
          }
          if (length(promo_ids) > 0L) is_promoted_bronze[promo_ids] <- TRUE
        }
      }
    }
  }

  active_mask <- is_gs | is_promoted_bronze
  anchor_tier_source <- ifelse(
    is_gold, "Tier_1_Gold_Anchor",
    ifelse(
      is_silver, "Tier_2_Silver_Anchor",
      ifelse(is_promoted_bronze, "Tier_3_Bronze_Promoted", "Passive_Projected")
    )
  )

  list(
    active_mask = as.logical(active_mask),
    is_gold_anchor = as.logical(is_gold),
    is_silver_anchor = as.logical(is_silver),
    is_promoted_bronze = as.logical(is_promoted_bronze),
    is_promoted_anchor = as.logical(is_promoted_bronze),
    anchor_tier_source = as.character(anchor_tier_source)
  )
}

compute_ledoit_wolf_cov <- function(
  S_w,
  X_centered = NULL,
  n_obs = NULL,
  min_lambda = 1e-5,
  max_lambda = 0.25
) {
  S_mat <- as.matrix(S_w)
  p <- ncol(S_mat)
  dn <- dimnames(S_mat)
  S_mat[!is.finite(S_mat)] <- 0.0
  S_mat <- 0.5 * (S_mat + t(S_mat))
  min_l <- max(0.0, as.numeric(min_lambda))
  max_l <- max(min_l, min(1.0, as.numeric(max_lambda)))

  tr_s <- sum(diag(S_mat))
  mu_s <- tr_s / max(1L, p)
  if (!is.finite(mu_s) || mu_s <= 1e-12) {
    mu_s <- 1e-4
  }
  if (p <= 1L) {
    out_1d <- matrix(max(as.numeric(S_mat[1L, 1L]), mu_s * min_l, 1e-8), 1L, 1L, dimnames = dn)
    attr(out_1d, "shrinkage_lambda") <- min_l
    return(out_1d)
  }

  target_mat <- diag(mu_s, p)
  diff_target <- S_mat - target_mat
  delta2 <- sum(diff_target * diff_target) / p

  lambda_raw <- min_l
  if (is.finite(delta2) && delta2 > 1e-16) {
    if (!is.null(X_centered)) {
      Xc <- as.matrix(X_centered)
      valid_rows <- apply(Xc, 1L, function(r) all(is.finite(r)))
      Xc <- Xc[valid_rows, , drop = FALSE]
      n_c <- nrow(Xc)
      if (n_c >= 2L && ncol(Xc) == p) {
        emp_tr <- sum(Xc * Xc) / n_c
        if (is.finite(emp_tr) && emp_tr > 1e-12 && tr_s > 1e-12) {
          Xc <- Xc * sqrt(tr_s / emp_tr)
        }
        q_i <- rowSums(Xc * Xc)
        s_norm2 <- sum(S_mat * S_mat)
        cross_i <- rowSums((Xc %*% S_mat) * Xc)
        b_bar2 <- mean(pmax(0.0, q_i * q_i - 2.0 * cross_i + s_norm2)) / (p * n_c)
        beta2 <- min(delta2, max(0.0, b_bar2))
        lambda_raw <- beta2 / delta2
      }
    } else {
      n_eff <- if (!is.null(n_obs) && is.finite(as.numeric(n_obs))) max(2.0, as.numeric(n_obs)) else max(2.0, as.numeric(p) + 2.0)
      tr_s2 <- sum(S_mat * S_mat)
      denom_oas <- (n_eff + 1.0 - 2.0 / p) * (tr_s2 - (tr_s * tr_s) / p)
      if (is.finite(denom_oas) && denom_oas > 1e-16) {
        num_oas <- (1.0 - 2.0 / p) * tr_s2 + (tr_s * tr_s)
        lambda_raw <- num_oas / denom_oas
      }
    }
  }

  lambda <- max(min_l, min(max_l, if (is.finite(lambda_raw)) lambda_raw else min_l))
  S_reg <- (1.0 - lambda) * S_mat + lambda * target_mat
  S_reg <- 0.5 * (S_reg + t(S_reg))
  dimnames(S_reg) <- dn
  attr(S_reg, "shrinkage_lambda") <- lambda
  S_reg
}

solve_shrinkage_cda_eigen <- function(S_reg, B, max_axes = NULL) {
  S_mat <- 0.5 * (as.matrix(S_reg) + t(as.matrix(S_reg)))
  B_mat <- 0.5 * (as.matrix(B) + t(as.matrix(B)))
  p <- ncol(S_mat)
  k_keep <- if (is.null(max_axes)) p else max(1L, min(as.integer(max_axes), p))

  chol_R <- tryCatch(base::chol(S_mat), error = function(e) NULL)
  if (is.null(chol_R)) {
    jitter_s <- max(1e-6, 1e-4 * mean(diag(S_mat)))
    S_mat <- S_mat + diag(jitter_s, p)
    chol_R <- tryCatch(base::chol(S_mat), error = function(e) NULL)
  }

  if (!is.null(chol_R)) {
    # Symmetric whitened generalized eigensystem: M = R^{-T} B R^{-1}
    left_whiten <- matrix(backsolve(chol_R, B_mat, transpose = TRUE), nrow = p, ncol = p)
    M_whiten <- t(matrix(backsolve(chol_R, t(left_whiten), transpose = TRUE), nrow = p, ncol = p))
    M_sym <- 0.5 * (M_whiten + t(M_whiten))
    eig_sym <- eigen(M_sym, symmetric = TRUE)
    vals <- pmax(0.0, as.numeric(eig_sym$values)[seq_len(k_keep)])
    U_sub <- eig_sym$vectors[, seq_len(k_keep), drop = FALSE]
    V_sub <- matrix(backsolve(chol_R, U_sub), nrow = p, ncol = k_keep)
  } else {
    eig_ns <- eigen(solve(S_mat, B_mat))
    real_idx <- which(abs(Im(eig_ns$values)) < 1e-6)
    if (length(real_idx) == 0L) real_idx <- seq_len(p)
    r_vals <- Re(eig_ns$values[real_idx])
    r_vecs <- Re(eig_ns$vectors[, real_idx, drop = FALSE])
    ord <- order(r_vals, decreasing = TRUE)
    take <- ord[seq_len(min(k_keep, length(ord)))]
    vals <- pmax(0.0, r_vals[take])
    V_sub <- r_vecs[, take, drop = FALSE]
  }

  for (j in seq_len(ncol(V_sub))) {
    v <- as.numeric(V_sub[, j])
    sc <- sqrt(max(as.numeric(crossprod(v, S_mat %*% v)), 1e-12))
    v <- v / sc
    dom_idx <- which.max(abs(v))
    if (length(dom_idx) == 1L && is.finite(v[dom_idx]) && v[dom_idx] < 0.0) {
      v <- -v
    }
    V_sub[, j] <- v
  }
  rownames(V_sub) <- colnames(S_mat)
  colnames(V_sub) <- paste0("Can", seq_len(ncol(V_sub)))
  list(values = vals, vectors = V_sub)
}

.normalize_determiner_id_vec <- function(df) {
  n <- nrow(df)
  raw_det <- rep("", n)
  invalid_raw <- c("na", "n/a", "none", "unknown", "null", "anonymous", "anon", "unspecified", "not recorded", "s.n.", "sn", "indet", "indet.")
  for (col in intersect(c("identifiedBy", "determiner_raw"), names(df))) {
    v <- as.character(df[[col]])
    v[is.na(v)] <- ""
    v <- trimws(v)
    need <- !nzchar(raw_det) & nzchar(v) & !(tolower(v) %in% invalid_raw)
    if (any(need)) raw_det[need] <- v[need]
  }
  out <- rep(NA_character_, n)
  has_raw <- nzchar(raw_det)
  if (!any(has_raw)) return(out)

  s <- tolower(raw_det[has_raw])
  for (iter in 1:2) {
    cand <- sub(
      "^(det\\.?|determined\\s+by|annot\\.?|annotated\\s+by|ver\\.?|verified\\s+by|rev\\.?|revised\\s+by|conf\\.?|confirmed\\s+by|identified\\s+by)[:;,\\s]+",
      "",
      s,
      perl = TRUE
    )
    cand_norm <- trimws(gsub("[[:punct:][:space:]]+", " ", cand, perl = TRUE))
    keep_strip <- nchar(cand_norm) >= 3L
    if (any(keep_strip)) s[keep_strip] <- cand[keep_strip]
  }
  s <- gsub("\\b(1[789][0-9]{2}|20[0-9]{2})\\b", " ", s, perl = TRUE)
  s <- gsub("[[:punct:][:space:]]+", " ", s, perl = TRUE)
  s <- trimws(s)
  invalid_norm <- c("", "na", "n a", "none", "unknown", "null", "anonymous", "anon", "unspecified", "not recorded", "s n", "indet", "unidentified")
  s[s %in% invalid_norm | nchar(s) < 2L] <- NA_character_
  out[has_raw] <- s
  out
}

.compute_loo_pc_likelihoods <- function(df, taxa_vec, det_tier, is_promoted, eval_mask, valid_taxa, prior_vec) {
  n <- nrow(df)
  G <- length(valid_taxa)
  if (G < 2L) return(NULL)

  pc_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
  sym_cols <- .discover_numeric_sorted_cols(names(df), "^Sym_PC[0-9]+$")
  use_cols <- if (length(pc_cols) >= 2L) {
    pc_cols
  } else if (length(sym_cols) >= 2L) {
    sym_cols
  } else if (length(pc_cols) == 1L) {
    pc_cols
  } else {
    sym_cols
  }
  if (length(use_cols) > 0L) {
    use_cols <- use_cols[vapply(use_cols, function(col) sum(is.finite(suppressWarnings(as.numeric(df[[col]])))) >= 3L, logical(1))]
  }
  if (length(use_cols) == 0L) return(NULL)

  X_pc <- do.call(cbind, lapply(use_cols, function(col) suppressWarnings(as.numeric(df[[col]]))))
  colnames(X_pc) <- use_cols
  finite_pc <- apply(X_pc, 1L, function(r) all(is.finite(r)))
  eval_pc <- eval_mask & finite_pc & (taxa_vec %in% valid_taxa)
  if (sum(eval_pc) < 3L) return(NULL)
  if (ncol(X_pc) > 2L) {
    sel_pc <- select_morphometric_subspace(X_pc[eval_pc, , drop = FALSE], group_vec = taxa_vec[eval_pc], min_p = 2L, max_p = NULL)
    X_pc <- X_pc[, sel_pc, drop = FALSE]
  }
  p <- ncol(X_pc)

  is_high_tier <- (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted
  ref_idx_by_k <- vector("list", G)
  names(ref_idx_by_k) <- valid_taxa
  taxa_with_ref <- logical(G)

  for (k in seq_len(G)) {
    tx <- valid_taxa[k]
    cand_hi <- which(eval_pc & (taxa_vec == tx) & is_high_tier)
    cand_all <- which(eval_pc & (taxa_vec == tx))
    cand_k <- if (length(cand_hi) >= 2L) cand_hi else cand_all
    if (length(cand_k) >= 5L) {
      Xk <- X_pc[cand_k, , drop = FALSE]
      med_k <- apply(Xk, 2L, stats::median)
      mad_k <- pmax(apply(Xk, 2L, stats::mad), 1e-3)
      r_d2 <- rowSums(sweep(sweep(Xk, 2L, med_k, "-"), 2L, mad_k, "/")^2)
      inlier_k <- r_d2 <= stats::qchisq(0.975, df = p)
      min_keep <- max(2L, ceiling(0.55 * length(cand_k)))
      if (sum(inlier_k) >= min_keep) {
        cand_k <- cand_k[inlier_k]
      } else {
        cand_k <- cand_k[order(r_d2)[seq_len(min_keep)]]
      }
    }
    ref_idx_by_k[[k]] <- cand_k
    taxa_with_ref[k] <- length(cand_k) >= 1L
  }
  if (sum(taxa_with_ref) < 2L) return(NULL)

  sum_w_k <- numeric(G)
  sum_wx_k <- matrix(0.0, nrow = G, ncol = p)
  in_ref_mat <- matrix(FALSE, nrow = n, ncol = G)
  res_list <- list()
  n_ref_total <- 0L
  w_ref_total <- 0.0

  for (k in which(taxa_with_ref)) {
    idx_k <- ref_idx_by_k[[k]]
    uk <- pmax(0.05, prior_vec[idx_k])
    Xk <- X_pc[idx_k, , drop = FALSE]
    sw <- sum(uk)
    swx <- colSums(Xk * uk)
    sum_w_k[k] <- sw
    sum_wx_k[k, ] <- swx
    in_ref_mat[idx_k, k] <- TRUE
    mk <- swx / sw
    diff_k <- sweep(Xk, 2L, mk, "-")
    res_list[[length(res_list) + 1L]] <- diff_k * sqrt(uk)
    n_ref_total <- n_ref_total + length(idx_k)
    w_ref_total <- w_ref_total + sw
  }

  X_res <- do.call(rbind, res_list)
  denom_w <- max(w_ref_total - sum(taxa_with_ref), 1.0)
  S_w <- crossprod(X_res) / denom_w
  if (sum(diag(S_w)) <= 1e-10) {
    S_w <- S_w + diag(1e-3, p)
  }
  S_reg <- compute_ledoit_wolf_cov(S_w, X_centered = X_res, n_obs = n_ref_total, min_lambda = 1e-4, max_lambda = 0.25)
  chol_R <- tryCatch(base::chol(S_reg), error = function(e) base::chol(S_reg + diag(1e-4, p)))

  obs_idx <- which(eval_mask & finite_pc)
  if (length(obs_idx) == 0L) return(NULL)

  d2_mat <- matrix(Inf, nrow = length(obs_idx), ncol = G)
  X_obs <- X_pc[obs_idx, , drop = FALSE]
  u_obs <- pmax(0.05, prior_vec[obs_idx])

  for (k in which(taxa_with_ref)) {
    mu_mat <- matrix(sum_wx_k[k, ], nrow = length(obs_idx), ncol = p, byrow = TRUE) / sum_w_k[k]
    is_in_k <- in_ref_mat[obs_idx, k]
    if (any(is_in_k) && length(ref_idx_by_k[[k]]) >= 2L) {
      rem_w <- sum_w_k[k] - u_obs[is_in_k]
      can_loo <- rem_w > 1e-8
      if (any(can_loo)) {
        loo_rows <- which(is_in_k)[can_loo]
        num_loo <- sweep(-X_obs[loo_rows, , drop = FALSE] * u_obs[loo_rows], 2L, sum_wx_k[k, ], "+")
        mu_mat[loo_rows, ] <- num_loo / rem_w[can_loo]
      }
    }
    diff_obs <- X_obs - mu_mat
    whiten_diff <- matrix(backsolve(chol_R, t(diff_obs), transpose = TRUE), nrow = p, ncol = nrow(diff_obs))
    d2_mat[, k] <- colSums(whiten_diff * whiten_diff)
  }

  temp_scale <- max(1.0, 0.5 * p)
  log_lik <- -0.5 * (d2_mat / temp_scale)
  row_max <- apply(log_lik, 1L, max)
  lik_sub <- exp(sweep(log_lik, 1L, row_max, "-"))
  lik_sub <- sweep(lik_sub, 1L, pmax(rowSums(lik_sub), 1e-15), "/")

  L_pc <- matrix(NA_real_, nrow = n, ncol = G, dimnames = list(NULL, valid_taxa))
  L_pc[obs_idx, ] <- lik_sub
  L_pc
}

.extract_specimen_discordance_signals <- function(df, eval_mask) {
  n <- nrow(df)
  err_cols <- intersect(
    c("dawid_skene_error_rate", "label_error_rate", "error_rate_k", "c_error", "loo_discordance", "loo_error_rate"),
    names(df)
  )
  eps_err <- rep(NA_real_, n)
  if (length(err_cols) > 0L) {
    for (ec in err_cols) {
      need <- eval_mask & !is.finite(eps_err)
      if (!any(need)) break
      v <- suppressWarnings(as.numeric(df[[ec]]))
      ok <- need & is.finite(v)
      if (any(ok)) eps_err[ok] <- pmax(0.0, pmin(1.0, v[ok]))
    }
  }

  apply_bool_flag <- function(cur_eps, flag_vec, pos_rate = 0.80, neg_rate = 0.05) {
    fv <- as.logical(flag_vec)
    ok_f <- eval_mask & !is.na(fv)
    if (!any(ok_f)) return(cur_eps)
    only_f <- ok_f & !is.finite(cur_eps)
    if (any(only_f)) {
      cur_eps[only_f] <- ifelse(fv[only_f], pos_rate, neg_rate)
    }
    both_pos <- ok_f & is.finite(cur_eps) & fv
    if (any(both_pos)) {
      cur_eps[both_pos] <- pmax(cur_eps[both_pos], 1.0 - (1.0 - cur_eps[both_pos]) * 0.35)
    }
    cur_eps
  }

  if ("is_label_corrupted" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$is_label_corrupted, pos_rate = 0.85, neg_rate = 0.05)
  }
  if ("morphometric_outlier_flag" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$morphometric_outlier_flag, pos_rate = 0.75, neg_rate = 0.05)
  }
  if ("multimodal_conflict_flag" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$multimodal_conflict_flag, pos_rate = 0.78, neg_rate = 0.05)
  }
  eps_err
}

compute_empirical_bayes_authority_weights <- function(
  df,
  mode = c("active_passive", "continuous_prior_posterior"),
  active_mask = NULL,
  prior_alpha = 5.0
) {
  mode <- match.arg(mode)
  if (is.null(df) || nrow(df) == 0L) return(numeric(0))
  n <- nrow(df)
  alpha0 <- max(1e-6, as.numeric(prior_alpha))

  det_tier <- .standardize_determiner_tier_vec(df)
  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = "species_standardized")
  det_id <- .normalize_determiner_id_vec(df)

  is_promoted <- rep(FALSE, n)
  if (!is.null(active_mask) && length(active_mask) == n) {
    act_in <- as.logical(active_mask)
    act_in[is.na(act_in)] <- FALSE
    is_promoted <- (det_tier == "Tier_3_Bronze") & act_in
  } else if ("is_promoted_bronze" %in% names(df)) {
    is_promoted <- as.logical(df$is_promoted_bronze)
    is_promoted[is.na(is_promoted)] <- FALSE
  } else if ("is_promoted_anchor" %in% names(df)) {
    is_promoted <- (det_tier == "Tier_3_Bronze") & as.logical(df$is_promoted_anchor)
    is_promoted[is.na(is_promoted)] <- FALSE
  } else if ("anchor_tier_source" %in% names(df)) {
    is_promoted <- as.character(df$anchor_tier_source) == "Tier_3_Bronze_Promoted"
    is_promoted[is.na(is_promoted)] <- FALSE
  }

  act_mask <- if (!is.null(active_mask) && length(active_mask) == n) {
    m <- as.logical(active_mask)
    m[is.na(m)] <- FALSE
    m
  } else {
    (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted
  }

  if (mode == "active_passive" && !any(act_mask)) {
    return(rep(0.0, n))
  }
  eval_mask <- if (mode == "active_passive") act_mask else rep(TRUE, n)

  tier_cell <- ifelse(
    det_tier == "Tier_1_Gold", "Tier_1_Gold",
    ifelse(
      det_tier == "Tier_2_Silver", "Tier_2_Silver",
      ifelse(is_promoted, "Tier_3_Bronze_Promoted", "Tier_3_Bronze")
    )
  )

  prior_map <- c(
    "Tier_1_Gold" = 0.95,
    "Tier_2_Silver" = 0.85,
    "Tier_3_Bronze_Promoted" = 0.65,
    "Tier_3_Bronze" = 0.35
  )
  pi_vec <- unname(prior_map[tier_cell])
  tier_ceil <- pi_vec + 0.35 * (1.0 - pi_vec)

  safe_logit <- function(p) {
    pc <- pmax(1e-5, pmin(1.0 - 1e-5, p))
    log(pc / (1.0 - pc))
  }

  blend_tier_and_determiner <- function(conc_vec, obs_mask) {
    theta_cell <- pi_vec
    cell_keys <- paste0(taxa_vec, "||", tier_cell)
    for (ck in unique(cell_keys[eval_mask])) {
      idx_all <- which(eval_mask & (cell_keys == ck))
      idx_obs <- idx_all[obs_mask[idx_all]]
      pi_t <- pi_vec[idx_all[1L]]
      n_c <- length(idx_obs)
      s_c <- if (n_c > 0L) sum(conc_vec[idx_obs]) else 0.0
      theta_cell[idx_all] <- (alpha0 * pi_t + s_c) / (alpha0 + n_c)
    }

    theta_blend <- theta_cell
    known_dets <- unique(det_id[eval_mask & !is.na(det_id)])
    if (length(known_dets) > 0L) {
      for (d_key in known_dets) {
        idx_d_all <- which(eval_mask & !is.na(det_id) & (det_id == d_key))
        idx_d_obs <- idx_d_all[obs_mask[idx_d_all]]
        n_d <- length(idx_d_obs)
        if (n_d == 0L) next
        pi_bar_d <- mean(pi_vec[idx_d_all])
        theta_d <- (alpha0 * pi_bar_d + sum(conc_vec[idx_d_obs])) / (alpha0 + n_d)
        if (length(unique(pi_vec[idx_d_all])) > 1L) {
          delta_d <- safe_logit(theta_d) - safe_logit(pi_bar_d)
          theta_d_i <- stats::plogis(safe_logit(pi_vec[idx_d_all]) + delta_d)
        } else {
          theta_d_i <- rep(theta_d, length(idx_d_all))
        }
        gamma_d <- n_d / (0.5 * alpha0 + n_d)
        theta_blend[idx_d_all] <- gamma_d * theta_d_i + (1.0 - gamma_d) * theta_cell[idx_d_all]
      }
    }
    pmax(0.05, pmin(0.995, theta_blend))
  }

  eps_err <- .extract_specimen_discordance_signals(df, eval_mask)

  is_art <- if ("is_artifact_driven" %in% names(df)) {
    av <- as.logical(df$is_artifact_driven)
    av[is.na(av)] <- FALSE
    av
  } else {
    rep(FALSE, n)
  }

  safe_std_tx <- function(v) {
    tx_v <- trimws(as.character(v))
    if (exists("standardize_packera_taxon", mode = "function")) {
      tx_v <- as.character(standardize_packera_taxon(tx_v))
    } else if (exists("standardize_taxon_name", mode = "function")) {
      std_v <- as.character(standardize_taxon_name(tx_v))
      tx_v <- ifelse(!is.na(std_v) & nzchar(std_v), std_v, tx_v)
    }
    tx_v
  }

  pred_sources <- list()
  if ("cda_predicted_taxon" %in% names(df)) {
    p_cda <- safe_std_tx(df$cda_predicted_taxon)
    c_cda <- if ("cda_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(0.82, n)
    c_cda[!is.finite(c_cda)] <- 0.82
    pred_sources[["cda"]] <- list(taxon = p_cda, conf = pmax(0.55, pmin(0.92, c_cda)), valid = rep(TRUE, n))
  }
  vis_col <- if ("vision_predicted_label" %in% names(df)) "vision_predicted_label" else if ("predicted_label" %in% names(df)) "predicted_label" else NULL
  if (!is.null(vis_col)) {
    p_vis <- safe_std_tx(df[[vis_col]])
    pred_sources[["vision"]] <- list(taxon = p_vis, conf = rep(0.82, n), valid = !is_art)
  }
  rf_col <- if ("rf_predicted_taxon" %in% names(df)) "rf_predicted_taxon" else if ("edaphic_best_fit_taxon" %in% names(df)) "edaphic_best_fit_taxon" else NULL
  if (!is.null(rf_col)) {
    p_rf <- safe_std_tx(df[[rf_col]])
    pred_sources[["rf"]] <- list(taxon = p_rf, conf = rep(0.78, n), valid = rep(TRUE, n))
  }

  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  valid_taxa <- sort(unique(taxa_vec[eval_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_tx)]))
  G <- length(valid_taxa)
  L_pc <- if (G >= 2L) .compute_loo_pc_likelihoods(df, taxa_vec, det_tier, is_promoted, eval_mask, valid_taxa, pi_vec) else NULL

  y_idx <- match(taxa_vec, valid_taxa)
  has_ds_Row <- rep(FALSE, n)
  L_mat <- NULL

  if (G >= 2L) {
    L_mat <- matrix(1.0 / G, nrow = n, ncol = G, dimnames = list(NULL, valid_taxa))
    ds_rows <- which(eval_mask & !is.na(y_idx))
    for (i in ds_rows) {
      yi <- y_idx[i]
      log_odds_vec <- numeric(0)
      wt_vec <- numeric(0)
      off_accum <- rep(0.0, G)
      off_count <- 0.0

      if (is.finite(eps_err[i])) {
        p_self_e <- max(0.005, min(0.995, 1.0 - eps_err[i]))
        lo_e <- log(p_self_e / (1.0 - p_self_e))
        log_odds_vec <- c(log_odds_vec, lo_e)
        wt_vec <- c(wt_vec, if (lo_e < 0.0) 2.2 else 1.4)
      }

      if (!is.null(L_pc) && all(is.finite(L_pc[i, ]))) {
        p_self_pc <- max(0.002, min(0.998, L_pc[i, yi]))
        odds_pc <- (G - 1.0) * p_self_pc / max(1.0 - p_self_pc, 1e-6)
        lo_pc <- log(max(odds_pc, 1e-4))
        log_odds_vec <- c(log_odds_vec, lo_pc)
        wt_vec <- c(wt_vec, if (lo_pc < 0.0) 1.5 else 1.0)
        off_v <- L_pc[i, ]
        off_v[yi] <- 0.0
        if (sum(off_v) > 1e-12) {
          off_accum <- off_accum + (off_v / sum(off_v))
          off_count <- off_count + 1.0
        }
      }

      if (length(pred_sources) > 0L) {
        for (ps in pred_sources) {
          if (!isTRUE(ps$valid[i])) next
          ptx <- ps$taxon[i]
          if (is.na(ptx) || (ptx %in% bad_tx)) next
          c_m <- ps$conf[i]
          if (ptx == valid_taxa[yi]) {
            lo_m <- log(c_m / max(1.0 - c_m, 1e-4))
            log_odds_vec <- c(log_odds_vec, lo_m)
            wt_vec <- c(wt_vec, 1.1)
          } else {
            lo_m <- log(max(1.0 - c_m, 0.02) / c_m)
            log_odds_vec <- c(log_odds_vec, lo_m)
            wt_vec <- c(wt_vec, 1.8)
            pk_idx <- match(ptx, valid_taxa)
            if (!is.na(pk_idx) && pk_idx != yi) {
              v_m <- rep((1.0 - c_m) / max(1L, G - 1L), G)
              v_m[yi] <- 0.0
              v_m[pk_idx] <- c_m
              if (sum(v_m) > 1e-12) {
                off_accum <- off_accum + (v_m / sum(v_m))
                off_count <- off_count + 1.2
              }
            }
          }
        }
      }

      if (length(log_odds_vec) > 0L) {
        has_ds_Row[i] <- TRUE
        comb_odds <- exp(sum(wt_vec * log_odds_vec) / sum(wt_vec))
        B_i <- if (comb_odds <= 1.0) {
          max(0.005, comb_odds)
        } else {
          1.0 + 0.75 * log(comb_odds) / (log(comb_odds) + 1.5)
        }
        l_self <- B_i / (B_i + G - 1.0)
        if (off_count > 0.0 && sum(off_accum) > 1e-12) {
          q_off <- off_accum / sum(off_accum)
        } else {
          q_off <- rep(1.0 / (G - 1.0), G)
          q_off[yi] <- 0.0
        }
        L_mat[i, ] <- (1.0 - l_self) * q_off
        L_mat[i, yi] <- l_self
      }
    }
  }

  w_eb <- pi_vec

  if (G >= 2L && any(has_ds_Row)) {
    # 2-Step Dawid-Skene EM Update (Dawid & Skene 1979) with Beta/Dirichlet Empirical Bayes Shrinkage
    T_mat <- matrix(0.0, nrow = n, ncol = G, dimnames = list(NULL, valid_taxa))
    ds_obs_idx <- which(has_ds_Row)
    for (i in ds_obs_idx) {
      yi <- y_idx[i]
      l_yi <- L_mat[i, yi]
      B_init <- (G - 1.0) * l_yi / max(1.0 - l_yi, 1e-8)
      r0 <- (pi_vec[i] * B_init) / max(pi_vec[i] * B_init + (1.0 - pi_vec[i]), 1e-12)
      r0 <- max(0.005, min(0.995, r0))
      off_l <- L_mat[i, ]
      off_l[yi] <- 0.0
      off_s <- sum(off_l)
      if (off_s > 1e-12) {
        T_mat[i, ] <- (1.0 - r0) * (off_l / off_s)
      } else {
        T_mat[i, ] <- (1.0 - r0) / (G - 1.0)
      }
      T_mat[i, yi] <- r0
    }

    theta_curr <- pi_vec
    for (em_step in 1:2) {
      # M-step: determiner + (species_standardized, determiner_tier) Beta-Binomial shrinkage & confusion matrix Q
      t_self_prev <- numeric(n)
      t_self_prev[ds_obs_idx] <- T_mat[cbind(ds_obs_idx, y_idx[ds_obs_idx])]
      conc_m <- pmin(t_self_prev, tier_ceil)
      theta_curr <- blend_tier_and_determiner(conc_m, has_ds_Row)

      Q_off <- matrix(1.0 / (G - 1.0), nrow = G, ncol = G)
      diag(Q_off) <- 0.0
      if (G > 2L) {
        prior_q <- alpha0 / (G - 1.0)
        for (k in seq_len(G)) {
          counts_k <- rep(prior_q, G)
          counts_k[k] <- 0.0
          for (m in seq_len(G)) {
            if (m == k) next
            m_rows <- ds_obs_idx[y_idx[ds_obs_idx] == m]
            if (length(m_rows) > 0L) {
              counts_k[m] <- counts_k[m] + sum(T_mat[m_rows, k])
            }
          }
          Q_off[k, ] <- counts_k / max(sum(counts_k), 1e-12)
        }
      }

      # E-step: update latent true-taxon posteriors T_{i, k}^{(s)}
      for (i in ds_obs_idx) {
        yi <- y_idx[i]
        th_i <- theta_curr[i]
        u_vec <- (1.0 - th_i) * Q_off[, yi] * L_mat[i, ]
        u_vec[yi] <- th_i * L_mat[i, yi]
        u_sum <- sum(u_vec)
        if (is.finite(u_sum) && u_sum > 1e-15) {
          T_mat[i, ] <- u_vec / u_sum
        }
      }
    }

    t_self_final <- numeric(n)
    t_self_final[ds_obs_idx] <- T_mat[cbind(ds_obs_idx, y_idx[ds_obs_idx])]
    conc_final <- pmin(t_self_final, tier_ceil)
    w_eb[eval_mask] <- theta_curr[eval_mask]
    w_eb[ds_obs_idx] <- pmax(0.05, pmin(0.995, 0.80 * theta_curr[ds_obs_idx] + 0.20 * conc_final[ds_obs_idx]))

    # Handle any evaluated rows with unknown taxon string that still carry explicit error rates
    rem_err <- which(eval_mask & !has_ds_Row & is.finite(eps_err))
    if (length(rem_err) > 0L) {
      fb_obs <- eval_mask & is.finite(eps_err)
      fb_conc <- pmin(pmax(0.0, 1.0 - eps_err), tier_ceil)
      fb_theta <- blend_tier_and_determiner(fb_conc, fb_obs)
      w_eb[rem_err] <- pmax(0.05, pmin(0.995, 0.80 * fb_theta[rem_err] + 0.20 * fb_conc[rem_err]))
    }
  } else {
    # Graceful single-taxon or feature-absent degradation: Beta-Binomial tier + determiner shrinkage
    eps_fb <- eps_err
    if (length(pred_sources) > 0L) {
      for (i in which(eval_mask & !is.finite(eps_fb))) {
        tx_i <- taxa_vec[i]
        if (is.na(tx_i) || (tx_i %in% bad_tx)) next
        mism <- numeric(0)
        for (ps in pred_sources) {
          if (!isTRUE(ps$valid[i])) next
          ptx <- ps$taxon[i]
          if (!is.na(ptx) && !(ptx %in% bad_tx)) {
            mism <- c(mism, if (ptx == tx_i) (1.0 - ps$conf[i]) else ps$conf[i])
          }
        }
        if (length(mism) > 0L) eps_fb[i] <- mean(mism)
      }
    }
    obs_fb <- eval_mask & is.finite(eps_fb)
    if (any(obs_fb)) {
      conc_fb <- pmin(pmax(0.0, 1.0 - eps_fb), tier_ceil)
      theta_fb <- blend_tier_and_determiner(conc_fb, obs_fb)
      w_eb[eval_mask] <- theta_fb[eval_mask]
      w_eb[obs_fb] <- pmax(0.05, pmin(0.995, 0.85 * theta_fb[obs_fb] + 0.15 * conc_fb[obs_fb]))
    }
  }

  w_eb[!is.finite(w_eb)] <- pi_vec[!is.finite(w_eb)]
  w_eb <- pmax(0.05, pmin(0.995, w_eb))

  if (mode == "active_passive") {
    w_out <- rep(0.0, n)
    if (any(act_mask)) {
      mean_act <- mean(w_eb[act_mask], na.rm = TRUE)
      w_out[act_mask] <- if (is.finite(mean_act) && mean_act > 1e-12) round(w_eb[act_mask] / mean_act, 6) else 1.0
    }
    return(w_out)
  }

  max_w <- max(w_eb, na.rm = TRUE)
  w_cont <- if (is.finite(max_w) && max_w > 1e-12) (3.0 * w_eb / max_w) else rep(1.0, n)
  w_cont[!is.finite(w_cont)] <- 1.0
  round(pmax(0.5, pmin(3.0, w_cont)), 6)
}
