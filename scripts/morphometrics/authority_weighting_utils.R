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

.prepare_anchor_feature_matrix <- function(df, feature_matrix = NULL, max_p = 5L) {
  n <- nrow(df)
  if (is.null(feature_matrix)) {
    auto_cols <- intersect(paste0("PC", seq_len(5)), names(df))
    if (length(auto_cols) == 0L) auto_cols <- intersect(paste0("Sym_PC", seq_len(5)), names(df))
    if (length(auto_cols) == 0L) auto_cols <- intersect(c(paste0("A", seq_len(12)), paste0("D", seq_len(12))), names(df))
    if (length(auto_cols) == 0L) return(NULL)
    feature_matrix <- df[, auto_cols, drop = FALSE]
  }
  X <- as.matrix(feature_matrix)
  storage.mode(X) <- "double"
  if (nrow(X) != n || ncol(X) == 0L) return(NULL)

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

  keep_cols <- keep_cols[seq_len(min(length(keep_cols), as.integer(max_p)))]
  X_use <- X[, keep_cols, drop = FALSE]
  for (j in seq_len(ncol(X_use))) {
    bad_j <- !is.finite(X_use[, j])
    if (any(bad_j)) {
      med_j <- stats::median(X_use[!bad_j, j], na.rm = TRUE)
      X_use[bad_j, j] <- if (is.finite(med_j)) med_j else 0.0
    }
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

  feat_prep <- .prepare_anchor_feature_matrix(df, feature_matrix = feature_matrix, max_p = 5L)
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

  err_cols <- intersect(
    c("dawid_skene_error_rate", "label_error_rate", "error_rate_k", "c_error", "loo_discordance", "loo_error_rate"),
    names(df)
  )
  eps_hat <- rep(NA_real_, n)
  if (length(err_cols) > 0L) {
    for (ec in err_cols) {
      need <- !is.finite(eps_hat)
      if (!any(need)) break
      v <- suppressWarnings(as.numeric(df[[ec]]))
      ok <- need & is.finite(v)
      if (any(ok)) eps_hat[ok] <- pmax(0.0, pmin(1.0, v[ok]))
    }
  }

  w_eb <- unname(prior_map[tier_cell])
  cell_keys <- paste0(taxa_vec, "||", tier_cell)
  for (ck in unique(cell_keys)) {
    idx_c <- which(cell_keys == ck)
    pi_tier <- unname(prior_map[tier_cell[idx_c[1]]])
    valid_eps <- eps_hat[idx_c]
    valid_eps <- valid_eps[is.finite(valid_eps)]
    n_c <- length(valid_eps)
    sum_rel <- if (n_c > 0L) sum(1.0 - valid_eps) else 0.0
    w_eb[idx_c] <- (alpha0 * pi_tier + sum_rel) / (alpha0 + n_c)
  }

  if (mode == "active_passive") {
    act_mask <- if (!is.null(active_mask) && length(active_mask) == n) {
      m <- as.logical(active_mask)
      m[is.na(m)] <- FALSE
      m
    } else {
      (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted
    }
    w_out <- rep(0.0, n)
    if (any(act_mask)) {
      mean_act <- mean(w_eb[act_mask], na.rm = TRUE)
      w_out[act_mask] <- if (is.finite(mean_act) && mean_act > 1e-12) round(w_eb[act_mask] / mean_act, 6) else 1.0
    }
    return(w_out)
  }

  max_w <- max(w_eb, na.rm = TRUE)
  w_cont <- if (is.finite(max_w) && max_w > 1e-12) (3.0 * w_eb / max_w) else rep(1.0, n)
  round(pmax(0.5, pmin(3.0, w_cont)), 6)
}
