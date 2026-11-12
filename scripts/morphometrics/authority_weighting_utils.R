#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/authority_weighting_utils.R
# Project: Packera Genus-Wide Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Unified Empirical Bayes / Dawid-Skene Hierarchical Authority-Tier Weighting
#   and Adaptive Reference Anchor Promotion utilities (Slenker et al. 2022;
#   Dawid & Skene 1979; Gelman et al. 2013; Kaufman & Rousseeuw 1990).
#
#   1. select_adaptive_reference_anchors():
#      - Identifies clean Basal_Rosette specimens across all active Packera taxa.
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

.identify_type_specimen_vec <- function(df) {
  if (is.null(df)) {
    return(list(type_rank = integer(0), is_type_specimen = logical(0), type_status_clean = character(0)))
  }
  if (!is.data.frame(df)) {
    raw_vec <- as.character(df)
    n <- length(raw_vec)
    explicit_bool <- rep(FALSE, n)
  } else {
    n <- nrow(df)
    if (n == 0L) {
      return(list(type_rank = integer(0), is_type_specimen = logical(0), type_status_clean = character(0)))
    }
    raw_vec <- if ("type_status" %in% names(df)) {
      as.character(df$type_status)
    } else if ("typeStatus" %in% names(df)) {
      as.character(df$typeStatus)
    } else {
      rep("", n)
    }
    explicit_bool <- if ("is_type_specimen" %in% names(df)) {
      b <- as.logical(df$is_type_specimen)
      b[is.na(b)] <- FALSE
      b
    } else {
      rep(FALSE, n)
    }
  }
  raw_vec[is.na(raw_vec)] <- ""
  raw_trim <- trimws(raw_vec)
  neg_mask <- grepl("not a type|non-type|non type|false type|no type|unverified|type locality|type species|type collection", raw_trim, ignore.case = TRUE)
  prim_mask <- !neg_mask & grepl("\\b(holotype|lectotype|neotype|epitype|isotype|isolectotype|isoneotype)s?\\b", raw_trim, ignore.case = TRUE)
  syntype_mask <- !neg_mask & !prim_mask & grepl("\\b(syntype|isosyntype)s?\\b", raw_trim, ignore.case = TRUE)
  sec_mask <- !neg_mask & !prim_mask & grepl("\\b(syntype|isosyntype|paratype|isoparatype|type|cotype)s?\\b", raw_trim, ignore.case = TRUE)

  type_rank <- integer(n)
  type_rank[sec_mask] <- 1L
  type_rank[prim_mask] <- 2L
  fb_mask <- (type_rank == 0L) & explicit_bool & !neg_mask
  if (any(fb_mask)) {
    type_rank[fb_mask] <- 1L
  }
  is_nb_type <- prim_mask | syntype_mask
  if (is.data.frame(df) && ("species_standardized" %in% names(df)) && any(is_nb_type)) {
    sp_std <- tolower(trimws(as.character(df$species_standardized)))
    sp_ep <- sub("^[a-z]+\\s+([a-z-]+).*$", "\\1", sp_std)
    has_of_binomial <- grepl("\\bof\\s+(packera|senecio|cineraria)\\s+[a-z]{3,}", raw_trim, ignore.case = TRUE)
    if (any(is_nb_type & has_of_binomial)) {
      of_binom <- tolower(sub("^.*\\bof\\s+((?:packera|senecio|cineraria)\\s+[a-z-]+).*$", "\\1", raw_trim))
      of_ep <- sub("^[a-z]+\\s+([a-z-]+).*$", "\\1", of_binom)
      stem6 <- function(s) substr(gsub("(us|a|um)$", "", s), 1L, 6L)
      dyn_equiv <- rep(FALSE, n)
      if (exists("standardize_taxon_name", mode = "function")) {
        chk_idx <- which(is_nb_type & has_of_binomial & nzchar(of_binom))
        if (length(chk_idx) > 0L) {
          resolved_of <- tolower(trimws(as.character(standardize_taxon_name(of_binom[chk_idx]))))
          dyn_equiv[chk_idx] <- !is.na(resolved_of) & nzchar(resolved_of) & (resolved_of == sp_std[chk_idx])
        }
      }
      het_mismatch <- is_nb_type & has_of_binomial & nzchar(sp_ep) & nzchar(of_ep) &
        (sp_ep != of_ep) & (stem6(sp_ep) != stem6(of_ep)) & !dyn_equiv
      is_nb_type[het_mismatch] <- FALSE
    }
  }
  clean_lbl <- ifelse(type_rank >= 1L & nzchar(raw_trim), raw_trim, ifelse(type_rank >= 1L, "TYPE", "None"))
  list(
    type_rank = as.integer(type_rank),
    is_type_specimen = as.logical(type_rank >= 1L),
    is_name_bearing_type = as.logical(is_nb_type),
    type_status_clean = as.character(clean_lbl)
  )
}

.standardize_determiner_tier_vec <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  raw <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_3_Bronze", n)
  raw[is.na(raw) | !nzchar(trimws(raw))] <- "Tier_3_Bronze"
  out <- rep("Tier_3_Bronze", n)
  out[grepl("tier_?2|silver", raw, ignore.case = TRUE)] <- "Tier_2_Silver"
  out[grepl("tier_?[01]|gold|type", raw, ignore.case = TRUE)] <- "Tier_1_Gold"
  t_info <- .identify_type_specimen_vec(df)
  if (any(t_info$is_type_specimen)) {
    out[t_info$is_type_specimen] <- "Tier_1_Gold"
  }
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
  if (chosen_col %in% c("species_standardized", "scientificName", "species_raw") && !any(grepl("^Morphotype_K|^Ecotype_", tx[!is.na(tx)]))) {
    if (exists("standardize_taxon_name", mode = "function")) {
      std_tx <- as.character(standardize_taxon_name(tx))
      tx <- ifelse(!is.na(std_tx) & nzchar(std_tx), std_tx, tx)
    } else if (exists("standardize_packera_taxon", mode = "function")) {
      tx <- as.character(standardize_packera_taxon(tx))
    }
  }
  tx[is.na(tx) | !nzchar(tx) | tolower(tx) == "unknown"] <- "Unassigned_Packera"
  tx
}

.identify_tier1_pure_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  if ("morph_tier_source" %in% names(df)) {
    mts <- trimws(as.character(df$morph_tier_source))
    is_t1 <- is.na(mts) | !nzchar(mts) | !(mts %in% c("Cauline_Only_Fallback", "No_Leaves"))
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

if (!exists("prune_collinear_dissection_scalars", mode = "function")) {
  for (.dl_p in c("scripts/morphometrics/dissected_leaf_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "dissected_leaf_utils.R"))) {
    if (!is.na(.dl_p) && nzchar(.dl_p) && file.exists(.dl_p)) { source(.dl_p); break }
  }
}
if (!exists("packera_mclapply", mode = "function")) {
  for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
    if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
  }
}

.select_morphometric_subspace_core <- function(X, group_vec = NULL, min_p = 2L, max_p = NULL, aux_df = NULL, allow_supervised_f = FALSE) {
  X <- as.matrix(X)
  M <- ncol(X); N <- nrow(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(M))
  if (M <= min_p) return(colnames(X))
  valid_g <- if (!is.null(group_vec)) stats::na.omit(as.character(group_vec)) else character(0)
  valid_g <- valid_g[nzchar(valid_g) & !(valid_g %in% c("Unassigned_Packera", "Unknown"))]
  g_tab <- table(valid_g); g_tab_multi <- g_tab[g_tab >= 2L]
  N_min <- if (length(g_tab_multi) >= 2L) min(as.integer(g_tab_multi)) else N
  P_max <- if (N_min >= 6L) max(min_p, floor(N_min / 3)) else M
  if (isTRUE(allow_supervised_f) && length(g_tab_multi) >= 2L) {
    # Pooled within-group Ledoit-Wolf CDA estimates a single regularized W across df_W = N - G anchors,
    # so a single rare endemic with small N_min must not collapse 13-species CDA down to 2 PCs.
    df_pooled <- max(1L, as.integer(N) - length(g_tab_multi))
    P_max_pooled <- max(min_p, min(16L, floor(df_pooled / 3L)))
    P_max <- max(P_max, min(M, P_max_pooled))
  }
  if (!is.null(max_p) && length(max_p) == 1L && is.finite(as.numeric(max_p)) && as.numeric(max_p) >= 1) {
    P_max <- min(P_max, as.integer(max_p))
  }

  g_fac <- NULL
  if (isTRUE(allow_supervised_f) && !is.null(group_vec) && length(group_vec) == N && length(g_tab_multi) >= 2L) {
    g_chr <- as.character(group_vec)
    g_chr[!(g_chr %in% names(g_tab_multi))] <- NA_character_
    g_fac <- droplevels(factor(g_chr))
    if (nlevels(g_fac) < 2L) g_fac <- NULL
  }

  pc_idx <- grep("^PC[0-9]+$", colnames(X))
  non_pc_idx <- setdiff(seq_len(M), pc_idx)
  kept_pc_cols <- pc_idx
  if (length(pc_idx) >= 3L) {
    v_obs <- apply(X[, pc_idx, drop = FALSE], 2, stats::var, na.rm = TRUE)
    v_prop <- v_obs / max(sum(v_obs, na.rm = TRUE), 1e-12)
    k_m <- length(pc_idx)
    bs_prop <- vapply(seq_len(k_m), function(j) sum(1 / (j:k_m)) / k_m, numeric(1))
    horn_thr <- (1 + sqrt(k_m / max(N, k_m + 1)))^2 / k_m
    keep_pc <- which(is.finite(v_prop) & (v_prop >= pmin(bs_prop * 0.55, horn_thr * 0.55)))
    if (length(keep_pc) < min_p) keep_pc <- sort(unique(c(seq_len(min(min_p, k_m)), keep_pc)))
    rem_pc <- setdiff(which(is.finite(v_prop) & v_prop >= 0.004), keep_pc)
    if (length(rem_pc) > 0L) {
      rescued <- integer(0)
      if (!is.null(aux_df) && is.data.frame(aux_df) && nrow(aux_df) == N) {
        dis_aux_cols <- intersect(c("solidity", "q25_solidity", "dissection_index"), names(aux_df))
        for (rp in rem_pc) {
          v_p <- X[, pc_idx[rp]]
          hit_r <- any(vapply(dis_aux_cols, function(ac) {
            va <- suppressWarnings(as.numeric(aux_df[[ac]]))
            ok <- is.finite(v_p) & is.finite(va)
            if (sum(ok) >= 8L && stats::sd(v_p[ok]) > 1e-8 && stats::sd(va[ok]) > 1e-8) {
              abs(suppressWarnings(stats::cor(v_p[ok], va[ok]))) >= 0.25
            } else FALSE
          }, logical(1)))
          if (isTRUE(hit_r)) rescued <- c(rescued, rp)
        }
      }
      if (!is.null(g_fac)) {
        p_rem <- vapply(rem_pc, function(rp) {
          v_p <- X[, pc_idx[rp]]
          ok <- is.finite(v_p) & !is.na(g_fac)
          if (sum(ok) < 8L || nlevels(droplevels(g_fac[ok])) < 2L) return(1.0)
          pv <- tryCatch(summary(stats::aov(v_p[ok] ~ droplevels(g_fac[ok])))[[1]][["Pr(>F)"]][1L], error = function(e) 1.0)
          if (is.finite(pv)) as.numeric(pv) else 1.0
        }, numeric(1))
        rescued <- union(rescued, rem_pc[stats::p.adjust(p_rem, method = "BH") < 0.01])
      }
      if (length(rescued) > 0L) keep_pc <- sort(unique(c(keep_pc, rescued)))
    }
    if (length(keep_pc) < min_p) keep_pc <- seq_len(min(min_p, k_m))
    kept_pc_cols <- pc_idx[keep_pc]
  }

  cand_idx <- c(kept_pc_cols, non_pc_idx)
  if (length(cand_idx) <= P_max) return(colnames(X)[cand_idx])
  if (!is.null(g_fac)) {
    core_n <- min(length(kept_pc_cols), min(2L, P_max))
    core_idx <- if (core_n > 0L) kept_pc_cols[seq_len(core_n)] else integer(0)
    pool_idx <- setdiff(cand_idx, core_idx)
    if (length(pool_idx) > 0L && length(core_idx) < P_max) {
      f_pool <- vapply(pool_idx, function(ci) {
        v_c <- X[, ci]; ok <- is.finite(v_c) & !is.na(g_fac)
        if (sum(ok) < 6L || nlevels(droplevels(g_fac[ok])) < 2L) return(0.0)
        fv <- tryCatch(summary(stats::aov(v_c[ok] ~ droplevels(g_fac[ok])))[[1]][["F value"]][1L], error = function(e) 0.0)
        if (is.finite(fv)) as.numeric(fv) else 0.0
      }, numeric(1))
      pick_pool <- pool_idx[order(f_pool, decreasing = TRUE)][seq_len(min(length(pool_idx), P_max - length(core_idx)))]
      final_idx <- cand_idx[cand_idx %in% c(core_idx, pick_pool)]
      return(colnames(X)[final_idx])
    }
  }
  colnames(X)[cand_idx[seq_len(P_max)]]
}

select_morphometric_subspace <- function(X, group_vec = NULL, min_p = 2L, max_p = NULL, aux_df = NULL, allow_supervised_f = FALSE) {
  sel <- .select_morphometric_subspace_core(X, group_vec = group_vec, min_p = min_p, max_p = max_p, aux_df = aux_df, allow_supervised_f = allow_supervised_f)
  message(sprintf("Subspace Selection (Horn Parallel Analysis / Broken-Stick & Williams & Titus 1988 P <= floor(N_min / 3)): retained P = %d / %d features.", length(sel), ncol(as.matrix(X))))
  sel
}

.prepare_anchor_feature_matrix <- function(df, feature_matrix = NULL, max_p = NULL, group_vec = NULL) {
  n <- nrow(df)
  if (is.null(feature_matrix)) {
    auto_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
    if (length(auto_cols) == 0L) auto_cols <- .discover_numeric_sorted_cols(names(df), "^[AD][0-9]+$")
    if (length(auto_cols) == 0L) return(NULL)
    feature_matrix <- df[, auto_cols, drop = FALSE]
  }
  X <- as.matrix(feature_matrix)
  storage.mode(X) <- "double"
  if (nrow(X) != n || ncol(X) == 0L) return(NULL)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))

  primary_cols <- grep("^PC[0-9]+$|^[AD]_?[0-9]+$", colnames(X))
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
  if (!any(grepl("^PC[0-9]+$", colnames(X_use))) && all(grepl("^[ABCD]_?[0-9]+$", colnames(X_use)))) {
    h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(X_use))))
    if (!any(is.na(h_nums))) X_use <- X_use[, order(h_nums, colnames(X_use)), drop = FALSE]
  }
  dyn_cap <- if (!is.null(max_p) && length(max_p) == 1L && is.finite(as.numeric(max_p)) && as.numeric(max_p) >= 1) {
    as.integer(max_p)
  } else {
    max(2L, min(ncol(X_use), floor(max(sum(has_feat), 4L) / 4L), 15L))
  }
  X_use <- X_use[, seq_len(min(ncol(X_use), dyn_cap)), drop = FALSE]
  list(X = X_use, has_feat = has_feat)
}

.compute_robust_taxon_medoid <- function(X_sub, S_inv = NULL) {
  X_sub <- as.matrix(X_sub)
  m <- nrow(X_sub)
  p <- ncol(X_sub)
  if (m == 0L) return(rep(0.0, p))
  if (m == 1L) return(as.numeric(X_sub[1, ]))
  if (m <= 4L) return(as.numeric(apply(X_sub, 2, stats::median, na.rm = TRUE)))
  if (m <= 600L) {
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
  } else {
    ref_idx <- unique(round(seq(1L, m, length.out = 400L)))
    if (is.null(S_inv)) {
      mads <- apply(X_sub, 2, stats::mad, na.rm = TRUE)
      sds <- apply(X_sub, 2, stats::sd, na.rm = TRUE)
      sc <- ifelse(is.finite(mads) & mads > 1e-8, mads, ifelse(is.finite(sds) & sds > 1e-8, sds, 1.0))
      Z <- sweep(X_sub, 2, sc, "/")
      Z_ref <- Z[ref_idx, , drop = FALSE]
      d2mat <- outer(rowSums(Z * Z), rowSums(Z_ref * Z_ref), "+") - 2.0 * tcrossprod(Z, Z_ref)
    } else {
      XS <- X_sub %*% S_inv
      X_ref <- X_sub[ref_idx, , drop = FALSE]
      d_all <- rowSums(XS * X_sub)
      d2mat <- outer(d_all, d_all[ref_idx], "+") - 2.0 * tcrossprod(XS, X_ref)
    }
    d2mat[d2mat < 0.0] <- 0.0
    row_d <- rowSums(sqrt(d2mat), na.rm = TRUE)
  }
  k_core <- max(2L, ceiling(0.60 * m))
  core_idx <- order(row_d)[seq_len(min(m, k_core))]
  as.numeric(colMeans(X_sub[core_idx, , drop = FALSE]))
}

.mcd_consistency_factor <- function(alpha, p) {
  alpha <- min(0.99, max(0.50, as.numeric(alpha)))
  denom <- stats::pchisq(stats::qchisq(alpha, df = p), df = p + 2)
  if (is.finite(denom) && denom > 1e-8) (alpha / denom) else 1.0
}

.fit_fast_mcd_cov <- function(X_mat, q_used = NULL, seed = 42L, allow_parallel = TRUE) {
  X_mat <- as.matrix(X_mat)
  n_m <- nrow(X_mat)
  p_m <- ncol(X_mat)
  if (n_m < max(12L, p_m + 2L) || !requireNamespace("MASS", quietly = TRUE)) return(NULL)
  if (is.null(q_used)) {
    q_used <- min(n_m - 1L, max(floor(0.80 * n_m), floor((n_m + p_m + 1L) / 2L)))
  }
  if (n_m <= 1000L) {
    return(tryCatch(MASS::cov.rob(X_mat, method = "mcd", quantile.used = q_used), error = function(e) NULL))
  }
  # Partitioned parallel FastMCD (Rousseeuw & Van Driessen 1999) + 2 full-sample C-steps across all n_m rows
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({ if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv) }, add = TRUE)
  set.seed(as.integer(seed))
  q_frac <- min(0.95, max(0.50, q_used / n_m))
  n_blk <- min(8L, max(2L, as.integer(ceiling(n_m / 600L))))
  blk_cap <- min(n_m, n_blk * 350L)
  pool_idx <- if (blk_cap < n_m) round(seq(1L, n_m, length.out = blk_cap)) else seq_len(n_m)
  perm_idx <- pool_idx[sample.int(length(pool_idx))]
  blocks <- split(perm_idx, rep_len(seq_len(n_blk), length(perm_idx)))
  fit_block <- function(b_idx) {
    Xb <- X_mat[b_idx, , drop = FALSE]
    nb <- nrow(Xb)
    qb <- min(nb - 1L, max(floor(q_frac * nb), floor((nb + p_m + 1L) / 2L)))
    tryCatch(MASS::cov.rob(Xb, method = "mcd", quantile.used = qb), error = function(e) NULL)
  }
  cand_fits <- if (isTRUE(allow_parallel) && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) {
    packera_mclapply(blocks, fit_block, mc.cores = resolve_packera_workers(min(16L, n_blk), n_blk))
  } else {
    lapply(blocks, fit_block)
  }
  cand_fits <- Filter(function(f) !is.null(f) && !is.null(f$cov) && !is.null(f$center), cand_fits)
  if (length(cand_fits) == 0L) return(NULL)
  jit_diag <- diag(1e-6, nrow = p_m, ncol = p_m)
  best_ld <- Inf; best_res <- NULL
  for (cf in cand_fits) {
    c_cur <- as.numeric(cf$center)
    s_cur <- as.matrix(cf$cov) + jit_diag
    ok_c <- TRUE
    for (c_step in 1:2) {
      d2 <- tryCatch(stats::mahalanobis(X_mat, center = c_cur, cov = s_cur), error = function(e) NULL)
      if (is.null(d2) || sum(is.finite(d2)) < q_used) { ok_c <- FALSE; break }
      inl <- order(d2, na.last = TRUE)[seq_len(q_used)]
      X_inl <- X_mat[inl, , drop = FALSE]
      c_cur <- colMeans(X_inl)
      s_raw <- stats::cov(X_inl)
      s_cur <- s_raw + jit_diag
    }
    if (!ok_c) next
    ld <- tryCatch(as.numeric(determinant(s_cur, logarithm = TRUE)$modulus), error = function(e) Inf)
    if (is.finite(ld) && ld < best_ld) {
      best_ld <- ld
      best_res <- list(center = c_cur, cov = s_raw)
    }
  }
  if (!is.null(best_res)) best_res else cand_fits[[1L]]
}

.estimate_pooled_within_cov <- function(X_use, taxa_vec, det_tier, valid_pool) {
  p <- ncol(X_use)
  u_pool_tx <- setdiff(unique(taxa_vec[valid_pool]), "Unassigned_Packera")
  build_resid <- function(tx) {
    idx_gold <- which(valid_pool & (taxa_vec == tx) & (det_tier == "Tier_1_Gold"))
    idx_gs <- which(valid_pool & (taxa_vec == tx) & (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")))
    idx_tx <- if (length(idx_gs) >= 3L) idx_gs else which(valid_pool & (taxa_vec == tx))
    if (length(idx_tx) < 3L) return(NULL)
    X_t <- X_use[idx_tx, , drop = FALSE]
    med_t <- if (length(idx_gold) >= 3L) {
      .compute_robust_taxon_medoid(X_use[idx_gold, , drop = FALSE])
    } else {
      .compute_robust_taxon_medoid(X_t)
    }
    diff_t <- sweep(X_t, 2, med_t, "-")
    nz_mask <- rowSums(diff_t^2) > 1e-12
    if (sum(nz_mask) >= 2L) diff_t <- diff_t[nz_mask, , drop = FALSE]
    if (length(idx_gs) < 3L && nrow(diff_t) >= 5L) {
      r_norm <- sqrt(rowSums(diff_t^2))
      keep_r <- r_norm <= stats::quantile(r_norm, probs = 0.75, na.rm = TRUE)
      if (sum(keep_r) >= 3L) diff_t <- diff_t[keep_r, , drop = FALSE]
    }
    df_scale <- sqrt(nrow(X_t) / max(1.0, nrow(X_t) - 1.0))
    diff_t * df_scale
  }
  resid_blocks <- if (length(u_pool_tx) >= 2L && sum(valid_pool) >= 1200L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) {
    packera_mclapply(u_pool_tx, build_resid, mc.cores = resolve_packera_workers(16L, length(u_pool_tx)))
  } else {
    lapply(u_pool_tx, build_resid)
  }
  resid_blocks <- Filter(Negate(is.null), resid_blocks)
  if (length(resid_blocks) > 0L) {
    R_all <- do.call(rbind, resid_blocks)
    if (nrow(R_all) >= (p + 2L)) {
      if (nrow(R_all) >= 25L && requireNamespace("MASS", quietly = TRUE)) {
        q_used <- min(nrow(R_all) - 1L, max(floor(0.80 * nrow(R_all)), floor((nrow(R_all) + p + 1L) / 2L)))
        rob <- .fit_fast_mcd_cov(R_all, q_used = q_used)
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

.compute_type_shrinkage_gamma <- function(n_eff_type, n_ref, type_prior_mult = 2.5, max_gamma = 0.35) {
  n_t <- max(0.0, as.numeric(n_eff_type))
  if (!is.finite(n_t) || n_t <= 0.0) return(0.0)
  n_r <- max(1.0, as.numeric(n_ref))
  min(as.numeric(max_gamma), (type_prior_mult * n_t) / (n_r + type_prior_mult * n_t))
}

.count_effective_type_gatherings <- function(df, row_idx) {
  k_len <- length(row_idx)
  if (k_len == 0L) return(0.0)
  if (!is.null(df) && "exsiccatae_key" %in% names(df)) {
    ek <- trimws(as.character(df$exsiccatae_key[row_idx]))
    valid_ek <- !is.na(ek) & nzchar(ek) & !(tolower(ek) %in% c("na", "none", "null", "unknown"))
    n_distinct_ek <- length(unique(ek[valid_ek]))
    n_no_ek <- sum(!valid_ek)
    return(as.numeric(max(1L, n_distinct_ek + n_no_ek)))
  }
  as.numeric(max(1L, k_len))
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
      is_type_anchor = logical(0),
      is_gold_anchor = logical(0),
      is_silver_anchor = logical(0),
      is_promoted_bronze = logical(0),
      is_confident_label_noise = logical(0),
      anchor_tier_source = character(0)
    ))
  }

  n <- nrow(df)
  min_anc <- max(1L, as.integer(min_anchors_per_taxon))
  type_info <- .identify_type_specimen_vec(df)
  type_rank <- type_info$type_rank
  det_tier <- .standardize_determiner_tier_vec(df)
  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  is_t1_pure <- .identify_tier1_pure_mask(df)
  if (sum(is_t1_pure) == 0L) is_t1_pure <- rep(TRUE, n)

  feat_prep <- .prepare_anchor_feature_matrix(df, feature_matrix = feature_matrix, max_p = NULL, group_vec = taxa_vec)
  X_use <- if (!is.null(feat_prep)) feat_prep$X else NULL
  has_feat <- if (!is.null(feat_prep)) feat_prep$has_feat else rep(TRUE, n)
  valid_pool <- is_t1_pure & has_feat

  is_type_valid <- (type_rank >= 1L) & valid_pool
  is_nb_type_valid <- as.logical(type_info$is_name_bearing_type) & valid_pool
  is_nb_type_valid[is.na(is_nb_type_valid)] <- FALSE
  is_gold <- ((det_tier == "Tier_1_Gold") | (type_rank >= 1L)) & valid_pool
  is_silver <- (det_tier == "Tier_2_Silver") & !is_gold & valid_pool
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

  q25_v <- if ("q25_solidity" %in% names(df)) suppressWarnings(as.numeric(df$q25_solidity)) else rep(NA_real_, n)
  sol_v <- if ("solidity" %in% names(df)) suppressWarnings(as.numeric(df$solidity)) else rep(NA_real_, n)
  dlf_v <- if ("dissected_leaf_fraction" %in% names(df)) suppressWarnings(as.numeric(df$dissected_leaf_fraction)) else rep(NA_real_, n)
  dis_b <- if ("is_dissected" %in% names(df)) (!is.na(df$is_dissected) & as.logical(df$is_dissected)) else rep(FALSE, n)
  eff_sol <- ifelse(is.finite(q25_v), q25_v, sol_v)
  is_dis_mode <- (is.finite(eff_sol) & eff_sol < 0.72) | (is.finite(dlf_v) & dlf_v >= 0.50) | (!is.finite(eff_sol) & dis_b)
  lqw_col <- if ("leaf_quality_weight" %in% names(df)) suppressWarnings(as.numeric(df$leaf_quality_weight)) else NULL

  unique_taxa <- unique(taxa_vec[valid_pool])
  par_taxa <- length(unique_taxa) >= 2L && n >= 800L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")
  eval_one_taxon <- function(tx) {
    pruned_out <- integer(0)
    promo_out <- integer(0)
    idx_gs <- which(is_gs & (taxa_vec == tx))
    idx_gs_clean <- idx_gs
    gs_type_local <- is_type_valid[idx_gs]
    gs_nb_type_local <- is_nb_type_valid[idx_gs]
    gs_paratype_local <- gs_type_local & !gs_nb_type_local
    gs_gold_local <- is_gold[idx_gs]
    gs_silver_local <- is_silver[idx_gs]
    n_gold_k <- sum(gs_gold_local)
    dis_gold_n <- sum(is_dis_mode[idx_gs[gs_gold_local]], na.rm = TRUE)
    dis_gold_frac <- if (n_gold_k > 0L) (dis_gold_n / n_gold_k) else if (length(idx_gs) > 0L) mean(is_dis_mode[idx_gs], na.rm = TRUE) else 0.0
    if (!is.finite(dis_gold_frac)) dis_gold_frac <- 0.0
    gold_sol_pool <- if (n_gold_k > 0L) eff_sol[idx_gs[gs_gold_local]] else eff_sol[idx_gs]
    dis_Mode_pool <- if (n_gold_k > 0L) is_dis_mode[idx_gs[gs_gold_local]] else is_dis_mode[idx_gs]
    ent_sol_vals <- gold_sol_pool[!dis_Mode_pool & is.finite(gold_sol_pool)]
    dis_sol_vals <- gold_sol_pool[dis_Mode_pool & is.finite(gold_sol_pool)]
    ent_med_sol <- if (length(ent_sol_vals) > 0L) stats::median(ent_sol_vals, na.rm = TRUE) else NA_real_
    dis_med_sol <- if (length(dis_sol_vals) > 0L) stats::median(dis_sol_vals, na.rm = TRUE) else NA_real_
    has_bimodal_dis_cluster <- (dis_gold_n >= 5L) && is.finite(dis_med_sol) &&
      ((dis_med_sol <= 0.65) || (is.finite(ent_med_sol) && ((ent_med_sol - dis_med_sol) >= 0.18)))
    min_dis_gold_frac <- if (has_bimodal_dis_cluster) 0.15 else 0.25
    allow_dis_split <- isTRUE(if (n_gold_k > 0L) (dis_gold_frac >= min_dis_gold_frac) else (dis_gold_frac >= 0.30))
    idx_type_k <- idx_gs[gs_type_local]
    n_type_k <- length(idx_type_k)
    n_eff_type_k <- if (n_type_k > 0L) .count_effective_type_gatherings(df, idx_type_k) else 0.0
    gamma_k <- if (n_type_k > 0L) .compute_type_shrinkage_gamma(n_eff_type_k, length(idx_gs)) else 0.0

    # Step 1: Within-taxon robust MCD / medoid Mahalanobis pruning of Gold + Silver (chi-sq_0.95, p)
    # Two-stage Gold/Silver seeding & heteroblasty gate: when a taxon has >= 5 Gold+Type anchors,
    # estimate initial reference scatter on Gold+Type first and screen Silver rows at min(mahal_q_gold_silver, 0.90).
    if (!is.null(X_use) && length(idx_gs) >= 5L) {
      p_use <- min(ncol(X_use), length(idx_gs) - 2L)
      if (p_use >= 1L) {
        X_gs_raw <- X_use[idx_gs, seq_len(p_use), drop = FALSE]
        v_ok <- apply(X_gs_raw, 2, stats::var, na.rm = TRUE) > 1e-10
        if (any(v_ok)) {
          col_sub <- which(v_ok)
          X_gs <- X_gs_raw[, col_sub, drop = FALSE]
          p_eff <- ncol(X_gs)
          extreme_type_local <- rep(FALSE, nrow(X_gs))
          if (any(gs_type_local) && sum(!gs_type_local) >= 5L) {
            nt_mask <- if (sum(gs_gold_local & !gs_type_local) >= 5L) (gs_gold_local & !gs_type_local) else !gs_type_local
            X_nt <- X_gs[nt_mask, , drop = FALSE]
            mu_nt <- .compute_robust_taxon_medoid(X_nt)
            cov_nt_raw <- stats::cov(X_nt)
            cov_nt_raw[!is.finite(cov_nt_raw)] <- 0.0
            reg_nt <- max(1e-4, 0.05 * max(mean(diag(cov_nt_raw), na.rm = TRUE), 1e-4))
            S_nt <- cov_nt_raw + diag(reg_nt, nrow = p_eff, ncol = p_eff)
            type_pos <- which(gs_type_local)
            d2_type_prelim <- tryCatch(
              stats::mahalanobis(X_gs[type_pos, , drop = FALSE], center = mu_nt, cov = S_nt),
              error = function(e) rep(0.0, length(type_pos))
            )
            extreme_type_local[type_pos] <- is.finite(d2_type_prelim) & (d2_type_prelim > stats::qchisq(0.99, df = p_eff))
          }
          non_ext_type_pos <- which(gs_type_local & !extreme_type_local)
          if (n_type_k > 0L && length(non_ext_type_pos) == 0L && sum(!gs_type_local) >= 5L) {
            mu_type_gs <- NULL
            gamma_k <- 0.0
          } else if (length(non_ext_type_pos) > 0L) {
            mu_type_gs <- .compute_robust_taxon_medoid(X_gs[non_ext_type_pos, , drop = FALSE])
          } else {
            mu_type_gs <- NULL
          }
          gs_dis_pos <- if (allow_dis_split) which(is_dis_mode[idx_gs]) else integer(0)
          gs_ent_pos <- if (allow_dis_split) which(!is_dis_mode[idx_gs]) else seq_len(nrow(X_gs))
          S_sub_pool <- if (!is.null(S_pool)) S_pool[col_sub, col_sub, drop = FALSE] else (stats::cov(X_gs) + diag(1e-4, nrow = p_eff, ncol = p_eff))
          if (allow_dis_split && length(gs_dis_pos) >= 3L && length(gs_ent_pos) >= 3L) {
            drop_gs <- rep(FALSE, nrow(X_gs))
            for (pos_m in list(gs_ent_pos, gs_dis_pos)) {
              Xm <- X_gs[pos_m, , drop = FALSE]
              pos_m_gold <- pos_m[gs_gold_local[pos_m] & !extreme_type_local[pos_m]]
              has_gold_med_m <- length(pos_m_gold) >= 3L
              use_gold_m <- length(pos_m_gold) >= max(5L, p_eff + 3L)
              pos_m_fit <- if (use_gold_m) pos_m_gold else pos_m[!extreme_type_local[pos_m]]
              if (length(pos_m_fit) < 3L) pos_m_fit <- pos_m
              Xm_fit <- X_gs[pos_m_fit, , drop = FALSE]
              rob_m <- NULL; q_m <- nrow(Xm_fit)
              if (nrow(Xm_fit) >= max(12L, 5L * p_eff) && requireNamespace("MASS", quietly = TRUE)) {
                q_m <- min(nrow(Xm_fit) - 1L, max(floor(0.80 * nrow(Xm_fit)), floor((nrow(Xm_fit) + p_eff + 1L) / 2L)))
                rob_m <- .fit_fast_mcd_cov(Xm_fit, q_used = q_m, allow_parallel = !par_taxa)
              }
              cent_raw_m <- if (has_gold_med_m) {
                if (!is.null(rob_m) && use_gold_m) rob_m$center else .compute_robust_taxon_medoid(X_gs[pos_m_gold, , drop = FALSE])
              } else if (!is.null(rob_m)) {
                rob_m$center
              } else {
                .compute_robust_taxon_medoid(Xm_fit)
              }
              cent_m <- cent_raw_m
              cov_shift_m <- matrix(0.0, nrow = p_eff, ncol = p_eff)
              if (!is.null(mu_type_gs) && gamma_k > 0.0) {
                sub_types <- pos_m[gs_type_local[pos_m] & !extreme_type_local[pos_m]]
                mu_sub <- if (length(sub_types) > 0L) .compute_robust_taxon_medoid(X_gs[sub_types, , drop = FALSE]) else mu_type_gs
                cent_m <- (1.0 - gamma_k) * cent_raw_m + gamma_k * mu_sub
                d_mu_m <- cent_m - cent_raw_m
                cov_shift_m <- tcrossprod(d_mu_m) + diag(sum(d_mu_m^2) / max(1.0, p_eff), nrow = p_eff, ncol = p_eff)
              }
              cov_m <- (if (!is.null(rob_m)) {
                .mcd_consistency_factor(q_m / nrow(Xm_fit), p_eff) * rob_m$cov + diag(1e-5, nrow = p_eff, ncol = p_eff)
              } else if (nrow(Xm_fit) >= (p_eff + 3L)) {
                0.65 * stats::cov(Xm_fit) + 0.35 * S_sub_pool + diag(1e-5, nrow = p_eff, ncol = p_eff)
              } else {
                S_sub_pool + diag(1e-4, nrow = p_eff, ncol = p_eff)
              }) + cov_shift_m
              d2_m <- tryCatch(stats::mahalanobis(Xm, center = cent_m, cov = cov_m), error = function(e) rep(0.0, nrow(Xm)))
              q_thr_m <- ifelse(
                gs_paratype_local[pos_m],
                stats::qchisq(0.99, df = p_eff),
                ifelse(
                  gs_silver_local[pos_m] & has_gold_med_m,
                  stats::qchisq(min(mahal_q_gold_silver, 0.90), df = p_eff),
                  stats::qchisq(mahal_q_gold_silver, df = p_eff)
                )
              )
              dr_m <- is.finite(d2_m) & (d2_m > q_thr_m)
              if (any(dr_m) && sum(!dr_m) >= 2L) drop_gs[pos_m[dr_m]] <- TRUE
            }
          } else {
            gold_core_pos <- which(gs_gold_local & !extreme_type_local)
            has_gold_medoid <- length(gold_core_pos) >= 3L
            use_gold_core <- length(gold_core_pos) >= max(5L, p_eff + 3L)
            X_gs_fit <- if (use_gold_core) {
              X_gs[gold_core_pos, , drop = FALSE]
            } else if (any(extreme_type_local) && sum(!extreme_type_local) >= 5L) {
              X_gs[!extreme_type_local, , drop = FALSE]
            } else {
              X_gs
            }
            rob_gs <- NULL; q_mcd <- nrow(X_gs_fit)
            if (nrow(X_gs_fit) >= max(12L, 5L * p_eff) && requireNamespace("MASS", quietly = TRUE)) {
              q_mcd <- min(nrow(X_gs_fit) - 1L, max(floor(0.80 * nrow(X_gs_fit)), floor((nrow(X_gs_fit) + p_eff + 1L) / 2L)))
              rob_gs <- .fit_fast_mcd_cov(X_gs_fit, q_used = q_mcd, allow_parallel = !par_taxa)
            }
            center_raw_gs <- if (has_gold_medoid) {
              if (!is.null(rob_gs) && use_gold_core) rob_gs$center else .compute_robust_taxon_medoid(X_gs[gold_core_pos, , drop = FALSE])
            } else if (!is.null(rob_gs)) {
              rob_gs$center
            } else {
              .compute_robust_taxon_medoid(X_gs_fit)
            }
            center_gs <- center_raw_gs
            cov_shift_gs <- matrix(0.0, nrow = p_eff, ncol = p_eff)
            if (!is.null(mu_type_gs) && gamma_k > 0.0) {
              center_gs <- (1.0 - gamma_k) * center_raw_gs + gamma_k * mu_type_gs
              d_mu_gs <- center_gs - center_raw_gs
              cov_shift_gs <- tcrossprod(d_mu_gs) + diag(sum(d_mu_gs^2) / max(1.0, p_eff), nrow = p_eff, ncol = p_eff)
            }
            cov_gs <- (if (!is.null(rob_gs)) {
              .mcd_consistency_factor(q_mcd / nrow(X_gs_fit), p_eff) * rob_gs$cov + diag(1e-5, nrow = p_eff, ncol = p_eff)
            } else if (nrow(X_gs_fit) >= (p_eff + 3L)) {
              0.75 * stats::cov(X_gs_fit) + 0.25 * S_sub_pool + diag(1e-4, nrow = p_eff, ncol = p_eff)
            } else {
              0.50 * stats::cov(X_gs_fit) + 0.50 * S_sub_pool + diag(1e-4, nrow = p_eff, ncol = p_eff)
            }) + cov_shift_gs
            d2_gs <- tryCatch(stats::mahalanobis(X_gs, center = center_gs, cov = cov_gs), error = function(e) rep(0.0, nrow(X_gs)))
            q_thr_gs <- ifelse(
              gs_paratype_local,
              stats::qchisq(0.99, df = p_eff),
              ifelse(
                gs_silver_local & has_gold_medoid,
                stats::qchisq(min(mahal_q_gold_silver, 0.90), df = p_eff),
                stats::qchisq(mahal_q_gold_silver, df = p_eff)
              )
            )
            drop_gs <- is.finite(d2_gs) & (d2_gs > q_thr_gs)
          }
          # Nomenclatural type pruning immunity: never prune verified name-bearing types (Holotype/Lectotype/Neotype/Epitype/Isotype/Syntype)
          if (any(gs_nb_type_local)) {
            drop_gs[gs_nb_type_local] <- FALSE
          }
          if (any(drop_gs) && sum(!drop_gs) >= 3L) {
            pruned_out <- idx_gs[drop_gs]
            idx_gs_clean <- idx_gs[!drop_gs]
          }
        }
      }
    }

    # Step 2: Adaptive Anchor Promotion for Gold/Silver-poor taxa (raw GS < min_anchors_per_taxon or clean n_gs < 3)
    n_gs <- length(idx_gs_clean)
    if (tx != "Unassigned_Packera" && (length(idx_gs) < min_anc || n_gs < min(3L, min_anc))) {
      idx_br <- which((det_tier == "Tier_3_Bronze") & valid_pool & (taxa_vec == tx))
      n_needed <- min_anc - n_gs
      if (length(idx_br) > 0L && n_needed > 0L) {
        if (is.null(X_use) || ncol(X_use) == 0L) {
          ord_q <- if (!is.null(lqw_col)) {
            order(lqw_col[idx_br], decreasing = TRUE, na.last = TRUE)
          } else {
            seq_along(idx_br)
          }
          promo_out <- idx_br[ord_q[seq_len(min(length(idx_br), n_needed))]]
        } else {
          p_eff <- ncol(X_use)
          ref_ids <- if (n_gs >= 1L) idx_gs_clean else idx_br
          X_ref <- X_use[ref_ids, , drop = FALSE]
          medoid_tx <- .compute_robust_taxon_medoid(X_ref, S_inv = S_pool_inv)
          if (n_type_k > 0L) {
            mu_type_tx <- .compute_robust_taxon_medoid(X_use[idx_gs[gs_type_local], , drop = FALSE], S_inv = S_pool_inv)
            medoid_tx <- (1.0 - gamma_k) * medoid_tx + gamma_k * mu_type_tx
          }
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
          ref_dis <- if (allow_dis_split) ref_ids[is_dis_mode[ref_ids]] else integer(0)
          ref_ent <- if (allow_dis_split) ref_ids[!is_dis_mode[ref_ids]] else ref_ids
          if (allow_dis_split && length(ref_dis) >= 2L && length(ref_ent) >= 2L) {
            med_dis <- .compute_robust_taxon_medoid(X_use[ref_dis, , drop = FALSE], S_inv = S_pool_inv)
            med_ent <- .compute_robust_taxon_medoid(X_use[ref_ent, , drop = FALSE], S_inv = S_pool_inv)
            d2_dis <- tryCatch(stats::mahalanobis(X_use[idx_br, , drop = FALSE], center = med_dis, cov = cov_tx), error = function(e) d2_br)
            d2_ent <- tryCatch(stats::mahalanobis(X_use[idx_br, , drop = FALSE], center = med_ent, cov = cov_tx), error = function(e) d2_br)
            d2_br <- pmin(d2_br, d2_dis, d2_ent, na.rm = TRUE)
          }
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
          if (length(promo_ids) > 0L) promo_out <- promo_ids
        }
      }
    }
    list(pruned = pruned_out, promo = promo_out)
  }

  tx_results <- if (length(unique_taxa) >= 2L && n >= 800L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) {
    packera_mclapply(unique_taxa, eval_one_taxon, mc.cores = resolve_packera_workers(min(16L, length(unique_taxa)), length(unique_taxa)))
  } else {
    lapply(unique_taxa, eval_one_taxon)
  }
  all_pruned <- unlist(lapply(tx_results, `[[`, "pruned"), use.names = FALSE)
  all_promo <- unlist(lapply(tx_results, `[[`, "promo"), use.names = FALSE)
  if (length(all_pruned) > 0L) {
    all_pruned <- setdiff(all_pruned, which(is_nb_type_valid))
    if (length(all_pruned) > 0L) {
      is_gs[all_pruned] <- FALSE
      is_gold[all_pruned] <- FALSE
      is_silver[all_pruned] <- FALSE
    }
  }
  if (length(all_promo) > 0L) {
    is_promoted_bronze[all_promo] <- TRUE
  }

  active_mask <- is_gs | is_promoted_bronze

  # Step 3: Confident Learning Out-of-Fold Anchor Pruning (Northcutt et al. 2021)
  is_confident_label_noise <- rep(FALSE, n)
  bad_cl_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  cl_act_idx <- which(active_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_cl_tx))
  taxa_levels <- sort(unique(taxa_vec[cl_act_idx]))
  if (length(taxa_levels) >= 2L && !is.null(X_use) && ncol(X_use) >= 1L) {
    K_cl <- length(taxa_levels)
    p_cl <- min(ncol(X_use), max(1L, floor(length(cl_act_idx) / 3L)), 10L)
    X_cl <- X_use[, seq_len(p_cl), drop = FALSE]
    idx_by_tx <- lapply(taxa_levels, function(tx_k) which(active_mask & (taxa_vec == tx_k)))
    names(idx_by_tx) <- taxa_levels
    if (sum(vapply(idx_by_tx, length, integer(1L)) >= 2L) >= 2L) {
      S_cl <- if (!is.null(S_pool)) {
        S_pool[seq_len(p_cl), seq_len(p_cl), drop = FALSE]
      } else {
        stats::cov(X_cl[cl_act_idx, , drop = FALSE]) + diag(1e-4, nrow = p_cl, ncol = p_cl)
      }
      S_cl <- 0.75 * S_cl + 0.25 * diag(pmax(diag(S_cl), 1e-4), nrow = p_cl, ncol = p_cl)
      S_cl_inv <- tryCatch(solve(S_cl), error = function(e) diag(1.0 / pmax(diag(S_cl), 1e-4), nrow = p_cl, ncol = p_cl))
      mu_cl_mat <- matrix(0.0, nrow = K_cl, ncol = p_cl, dimnames = list(taxa_levels, NULL))
      for (k in seq_len(K_cl)) {
        ik <- idx_by_tx[[k]]
        if (length(ik) > 0L) {
          gold_ik <- ik[is_gold[ik]]
          ref_ik <- if (length(gold_ik) >= 3L) gold_ik else ik
          mu_cl_mat[k, ] <- .compute_robust_taxon_medoid(X_cl[ref_ik, , drop = FALSE], S_inv = S_cl_inv)
        }
      }
      P_cl <- matrix(0.0, nrow = length(cl_act_idx), ncol = K_cl, dimnames = list(NULL, taxa_levels))
      y_cl_idx <- match(taxa_vec[cl_act_idx], taxa_levels)
      for (ii in seq_along(cl_act_idx)) {
        r_i <- cl_act_idx[ii]
        yk <- y_cl_idx[ii]
        mu_use <- mu_cl_mat
        ik <- idx_by_tx[[yk]]
        if (length(ik) >= 3L) {
          ik_loo <- setdiff(ik, r_i)
          gold_loo <- ik_loo[is_gold[ik_loo]]
          ref_loo <- if (length(gold_loo) >= 3L) gold_loo else ik_loo
          mu_use[yk, ] <- .compute_robust_taxon_medoid(X_cl[ref_loo, , drop = FALSE], S_inv = S_cl_inv)
        }
        diff_i <- sweep(mu_use, 2L, X_cl[r_i, ], "-")
        d2_i <- rowSums((diff_i %*% S_cl_inv) * diff_i)
        ll_i <- -0.5 * d2_i / max(1.0, 0.5 * p_cl)
        ll_i <- ll_i - max(ll_i, na.rm = TRUE)
        p_i <- exp(ll_i)
        P_cl[ii, ] <- p_i / max(sum(p_i, na.rm = TRUE), 1e-12)
      }
      t_j <- numeric(K_cl)
      for (k in seq_len(K_cl)) {
        m_k <- which(y_cl_idx == k)
        t_j[k] <- if (length(m_k) > 0L) mean(P_cl[m_k, k], na.rm = TRUE) else 0.5
      }
      for (k in seq_len(K_cl)) {
        m_k <- which(y_cl_idx == k)
        if (length(m_k) <= 3L) next
        cand_ii <- integer(0)
        for (ii in m_k) {
          r_i <- cl_act_idx[ii]
          if (is_nb_type_valid[r_i]) next
          p_self <- P_cl[ii, k]
          oth_k <- setdiff(seq_len(K_cl), k)
          best_oth <- oth_k[which.max(P_cl[ii, oth_k])]
          p_oth <- P_cl[ii, best_oth]
          if (is.finite(p_self) && is.finite(p_oth) &&
              (p_self < 0.45 * t_j[k]) && (p_oth >= max(0.55, t_j[best_oth]))) {
            cand_ii <- c(cand_ii, ii)
          }
        }
        if (length(cand_ii) > 0L) {
          max_drop_k <- max(0L, length(m_k) - 3L)
          if (max_drop_k > 0L) {
            ord_drop <- cand_ii[order(P_cl[cand_ii, k], decreasing = FALSE)]
            drop_rows <- cl_act_idx[ord_drop[seq_len(min(length(ord_drop), max_drop_k))]]
            is_confident_label_noise[drop_rows] <- TRUE
          }
        }
      }
      if (any(is_confident_label_noise)) {
        active_mask[is_confident_label_noise] <- FALSE
        is_gold[is_confident_label_noise] <- FALSE
        is_silver[is_confident_label_noise] <- FALSE
        is_promoted_bronze[is_confident_label_noise] <- FALSE
      }
    }
  }

  anchor_tier_source <- ifelse(
    is_confident_label_noise, "Passive_Confident_Noise_Demoted",
    ifelse(
      is_type_valid & active_mask, "Tier_0_Type_Anchor",
      ifelse(
        is_gold, "Tier_1_Gold_Anchor",
        ifelse(
          is_silver, "Tier_2_Silver_Anchor",
          ifelse(is_promoted_bronze, "Tier_3_Bronze_Promoted", "Passive_Projected")
        )
      )
    )
  )

  list(
    active_mask = as.logical(active_mask),
    is_type_anchor = as.logical(is_type_valid & active_mask),
    is_gold_anchor = as.logical(is_gold),
    is_silver_anchor = as.logical(is_silver),
    is_promoted_bronze = as.logical(is_promoted_bronze),
    is_confident_label_noise = as.logical(is_confident_label_noise),
    anchor_tier_source = as.character(anchor_tier_source)
  )
}

compute_ledoit_wolf_cov <- function(
  S_w,
  X_centered = NULL,
  n_obs = NULL,
  min_lambda = 1e-5,
  max_lambda = 0.25,
  weights = NULL,
  min_shrinkage = NULL,
  max_shrinkage = NULL
) {
  if (!is.null(min_shrinkage)) min_lambda <- min_shrinkage
  if (!is.null(max_shrinkage)) max_lambda <- max_shrinkage
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
  n_eff_used <- NA_real_
  if (is.finite(delta2) && delta2 > 1e-16) {
    if (!is.null(X_centered)) {
      Xc_in <- as.matrix(X_centered)
      valid_rows <- apply(Xc_in, 1L, function(r) all(is.finite(r)))
      Xc <- Xc_in[valid_rows, , drop = FALSE]
      n_raw <- nrow(Xc)
      w_sub <- if (!is.null(weights) && length(weights) == nrow(Xc_in)) suppressWarnings(as.numeric(weights))[valid_rows] else NULL
      use_w <- !is.null(w_sub) && all(is.finite(w_sub)) && sum(w_sub > 0.0) >= 2L
      if (use_w) {
        w_pos <- pmax(0.0, as.numeric(w_sub))
        n_eff <- (sum(w_pos)^2) / max(sum(w_pos^2), 1e-12)
        n_c <- max(2.0, n_eff)
        w_norm <- w_pos / sum(w_pos)
      } else {
        n_c <- as.numeric(n_raw)
        w_norm <- rep(1.0 / max(1L, n_raw), n_raw)
      }
      n_eff_used <- n_c
      if (n_raw >= 2L && ncol(Xc) == p) {
        emp_tr <- sum((Xc * Xc) * w_norm)
        if (is.finite(emp_tr) && emp_tr > 1e-12 && tr_s > 1e-12) {
          Xc <- Xc * sqrt(tr_s / emp_tr)
        }
        X_w <- Xc * sqrt(w_norm * n_c)
        Q_mat <- crossprod(X_w^2)
        q_i <- rowSums(Xc * Xc)
        s_norm2 <- sum(S_mat * S_mat)
        cross_i <- rowSums((Xc %*% S_mat) * Xc)
        u_w <- pmax(0.0, n_c * (w_norm^2))
        u_w <- u_w / max(sum(u_w), 1e-12)
        b_bar2 <- sum(u_w * pmax(0.0, q_i * q_i - 2.0 * cross_i + s_norm2)) / (p * n_c)
        beta2 <- min(delta2, max(0.0, b_bar2))
        lambda_raw <- beta2 / delta2
      }
    } else {
      w_ok <- !is.null(weights) && length(weights) >= 2L && sum(is.finite(weights) & weights > 0.0) >= 2L
      n_eff <- if (w_ok) {
        w_pos <- pmax(0.0, as.numeric(weights[is.finite(weights)]))
        max(2.0, (sum(w_pos)^2) / max(sum(w_pos^2), 1e-12))
      } else if (!is.null(n_obs) && is.finite(as.numeric(n_obs))) {
        max(2.0, as.numeric(n_obs))
      } else {
        max(2.0, as.numeric(p) + 2.0)
      }
      n_eff_used <- n_eff
      tr_s2 <- sum(S_mat * S_mat)
      denom_oas <- (n_eff + 1.0 - 2.0 / p) * (tr_s2 - (tr_s * tr_s) / p)
      if (is.finite(denom_oas) && denom_oas > 1e-16) {
        num_oas <- (1.0 - 2.0 / p) * tr_s2 + (tr_s * tr_s)
        lambda_raw <- num_oas / denom_oas
      } else {
        lambda_raw <- min(0.50, p / max(1.0, as.numeric(n_eff) * 2.0))
      }
    }
  }

  lambda <- max(min_l, min(max_l, if (is.finite(lambda_raw)) lambda_raw else min_l))
  S_reg <- (1.0 - lambda) * S_mat + lambda * target_mat
  S_reg <- 0.5 * (S_reg + t(S_reg))
  dimnames(S_reg) <- dn
  attr(S_reg, "shrinkage_lambda") <- lambda
  attr(S_reg, "n_eff") <- n_eff_used
  S_reg
}


build_taxon_blind_gathering_deme_folds <- function(coords_df, n_obs, deme_radius_km = 10.0) {
  if (is.null(coords_df) || !is.data.frame(coords_df) || nrow(coords_df) != n_obs || n_obs < 2L) {
    return(seq_len(n_obs))
  }
  if ("cv_deme_id" %in% names(coords_df) && any(!is.na(coords_df$cv_deme_id) & nzchar(as.character(coords_df$cv_deme_id)))) {
    cv_raw <- as.character(coords_df$cv_deme_id)
    cv_raw[is.na(cv_raw) | !nzchar(cv_raw)] <- paste0("ISO_", seq_len(n_obs))[is.na(cv_raw) | !nzchar(cv_raw)]
    return(match(cv_raw, unique(cv_raw)))
  }
  cat_v <- if ("catalogNumber" %in% names(coords_df)) trimws(as.character(coords_df$catalogNumber)) else as.character(seq_len(n_obs))
  cat_v[is.na(cat_v) | !nzchar(cat_v)] <- as.character(seq_len(n_obs))[is.na(cat_v) | !nzchar(cat_v)]
  exs_v <- if ("exsiccatae_key" %in% names(coords_df)) trimws(as.character(coords_df$exsiccatae_key)) else rep("", n_obs)
  exs_v[is.na(exs_v) | tolower(exs_v) %in% c("na", "none", "unknown")] <- ""
  grp_key <- ifelse(nzchar(exs_v), paste0("EXS__", exs_v), paste0("CAT__", cat_v))
  fold_id <- match(grp_key, unique(grp_key))

  lat_v <- if ("decimalLatitude" %in% names(coords_df)) suppressWarnings(as.numeric(coords_df$decimalLatitude)) else rep(NA_real_, n_obs)
  lon_v <- if ("decimalLongitude" %in% names(coords_df)) suppressWarnings(as.numeric(coords_df$decimalLongitude)) else rep(NA_real_, n_obs)
  ok_c <- which(is.finite(lat_v) & is.finite(lon_v) & abs(lat_v) <= 90.0 & abs(lon_v) <= 180.0)
  if (length(ok_c) >= 2L && is.finite(deme_radius_km) && deme_radius_km > 0.0) {
    # Complete-linkage (diameter-bounded) taxon-blind 10 km spatial locality grouping (prevents single-linkage chaining)
    lat_r <- lat_v[ok_c] * (pi / 180.0); lon_r <- lon_v[ok_c] * (pi / 180.0)
    clusters <- list()
    for (ii in seq_along(ok_c)) {
      placed <- FALSE
      for (cl_idx in seq_along(clusters)) {
        members <- clusters[[cl_idx]]
        dlat <- lat_r[members] - lat_r[ii]; dlon <- lon_r[members] - lon_r[ii]
        hav <- sin(dlat * 0.5)^2 + cos(lat_r[ii]) * cos(lat_r[members]) * (sin(dlon * 0.5)^2)
        dkm <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hav))))
        if (all(is.finite(dkm) & dkm <= deme_radius_km)) {
          clusters[[cl_idx]] <- c(members, ii); placed <- TRUE; break
        }
      }
      if (!placed) clusters[[length(clusters) + 1L]] <- ii
    }
    for (cl_idx in seq_along(clusters)) {
      rows_cl <- ok_c[clusters[[cl_idx]]]
      if (length(rows_cl) >= 2L) {
        old_f <- unique(fold_id[rows_cl])
        fold_id[fold_id %in% old_f] <- min(old_f)
      }
    }
  }
  match(fold_id, unique(fold_id))
}

assign_spatial_gathering_demes <- function(df, deme_radius_km = 10.0) {
  # Feature-preserving, taxon-blind spatial/gathering deme metadata assignment (ZERO feature shrinkage!)
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  if (n == 0L) return(df)
  if ("PC1" %in% names(df) && !("PC1_raw_individual" %in% names(df))) df$PC1_raw_individual <- df$PC1
  if ("PC2" %in% names(df) && !("PC2_raw_individual" %in% names(df))) df$PC2_raw_individual <- df$PC2
  fold_ids <- build_taxon_blind_gathering_deme_folds(df, n_obs = n, deme_radius_km = deme_radius_km)
  df$cv_deme_id <- sprintf("LOC_DEME_%04d", as.integer(fold_ids))
  df$population_deme_id <- df$cv_deme_id
  tab_d <- table(df$population_deme_id)
  df$deme_n_specimens <- as.integer(tab_d[df$population_deme_id])
  df
}

compute_weighted_shrinkage_cda_loocv <- function(
  X_a,
  y_a,
  w_a = NULL,
  active_mask = NULL,
  coords_df = NULL,
  X_cand = NULL,
  deme_radius_km = 10.0
) {
  # Purely morphological gathering/deme-grouped cross-validation:
  #   - Zero geographic coordinate priors (preserves Padial et al. 2010 modality independence)
  #   - Holds out entire gatherings / 10 km spatial localities together to prevent pseudoreplication leakage
  #   - Re-estimates feature subspace (when X_cand is supplied) and Ledoit-Wolf covariance strictly inside each training fold
  X_a <- as.matrix(X_a)
  storage.mode(X_a) <- "double"
  if (!is.null(X_cand)) {
    X_cand <- as.matrix(X_cand)
    storage.mode(X_cand) <- "double"
  }
  if (!is.null(active_mask) && length(active_mask) == nrow(X_a) && sum(active_mask, na.rm = TRUE) < nrow(X_a)) {
    act_idx <- which(as.logical(active_mask) & !is.na(active_mask))
    X_a <- X_a[act_idx, , drop = FALSE]
    if (!is.null(X_cand) && nrow(X_cand) == length(active_mask)) X_cand <- X_cand[act_idx, , drop = FALSE]
    y_a <- y_a[act_idx]
    if (!is.null(w_a) && length(w_a) == length(active_mask)) w_a <- w_a[act_idx]
    if (!is.null(coords_df) && is.data.frame(coords_df) && nrow(coords_df) == length(active_mask)) {
      coords_df <- coords_df[act_idx, , drop = FALSE]
    }
  } else if (!is.null(coords_df) && is.data.frame(coords_df) && nrow(coords_df) != nrow(X_a) &&
             !is.null(active_mask) && length(active_mask) == nrow(coords_df) &&
             sum(active_mask, na.rm = TRUE) == nrow(X_a)) {
    act_idx2 <- which(as.logical(active_mask) & !is.na(active_mask))
    coords_df <- coords_df[act_idx2, , drop = FALSE]
    if (!is.null(X_cand) && nrow(X_cand) == length(active_mask)) X_cand <- X_cand[act_idx2, , drop = FALSE]
  }
  y_f <- droplevels(as.factor(y_a))
  g_levels <- levels(y_f)
  K <- length(g_levels)
  N_a <- nrow(X_a)
  p <- ncol(X_a)
  if (N_a < K + 2L || K < 2L || p < 1L) return(NA_real_)
  if (is.null(w_a) || length(w_a) != N_a) {
    w_a <- rep(1.0, N_a)
  } else {
    w_a <- suppressWarnings(as.numeric(w_a))
    w_a[!is.finite(w_a) | w_a < 0.0] <- 0.0
    if (sum(w_a) <= 1e-12) w_a <- rep(1.0, N_a)
  }

  fold_vec <- build_taxon_blind_gathering_deme_folds(coords_df, n_obs = N_a, deme_radius_km = deme_radius_km)
  # Ensure holding out a spatial deme never wipes out an entire class in small cohorts; fall back to gathering/LOO for singleton-deme taxa
  for (k in seq_len(K)) {
    idx_k <- which(as.integer(y_f) == k)
    if (length(unique(fold_vec[idx_k])) < 2L && length(idx_k) >= 2L) {
      max_f <- max(fold_vec, na.rm = TRUE)
      fold_vec[idx_k] <- max_f + seq_along(idx_k)
    }
  }
  u_folds <- sort(unique(fold_vec))
  priors <- rep(1.0 / K, K)
  pred_cv <- character(N_a)

  for (fid in u_folds) {
    test_idx <- which(fold_vec == fid)
    train_idx <- which(fold_vec != fid)
    if (length(train_idx) < K) train_idx <- setdiff(seq_len(N_a), test_idx[1L])

    if (!is.null(X_cand) && nrow(X_cand) == N_a && ncol(X_cand) >= 2L && exists("select_morphometric_subspace", mode = "function")) {
      aux_tr <- if (!is.null(coords_df) && is.data.frame(coords_df)) coords_df[train_idx, , drop = FALSE] else NULL
      sel_f <- tryCatch(
        suppressMessages(select_morphometric_subspace(X_cand[train_idx, , drop = FALSE], group_vec = y_f[train_idx], aux_df = aux_tr, allow_supervised_f = TRUE)),
        error = function(e) colnames(X_a)
      )
      sel_f <- intersect(sel_f, colnames(X_cand))
      X_tr <- if (length(sel_f) >= 1L) X_cand[train_idx, sel_f, drop = FALSE] else X_a[train_idx, , drop = FALSE]
      X_te <- if (length(sel_f) >= 1L) X_cand[test_idx, sel_f, drop = FALSE] else X_a[test_idx, , drop = FALSE]
    } else {
      X_tr <- X_a[train_idx, , drop = FALSE]
      X_te <- X_a[test_idx, , drop = FALSE]
    }

    p_f <- ncol(X_tr)
    y_tr <- y_f[train_idx]
    w_tr <- w_a[train_idx]
    mu_k_f <- matrix(0.0, nrow = K, ncol = p_f)
    W_mat_f <- matrix(0.0, nrow = p_f, ncol = p_f)
    X_c_f <- matrix(0.0, nrow = nrow(X_tr), ncol = p_f)
    w_corr_f <- 0.0

    for (k in seq_len(K)) {
      ik <- which(as.integer(y_tr) == k)
      if (length(ik) == 0L) {
        # Fallback only if class had a single total anchor
        ik_all <- which(as.integer(y_f) == k)
        mu_k_f[k, ] <- colMeans(if (!is.null(X_cand) && all(colnames(X_tr) %in% colnames(X_cand))) X_cand[ik_all, colnames(X_tr), drop = FALSE] else X_a[ik_all, , drop = FALSE])
        next
      }
      wk <- w_tr[ik]; sw_k <- sum(wk)
      if (sw_k <= 1e-12) { wk <- rep(1.0, length(ik)); sw_k <- sum(wk) }
      w_corr_f <- w_corr_f + (sum(wk^2) / sw_k)
      Xk <- X_tr[ik, , drop = FALSE]
      mk <- colSums(Xk * wk) / sw_k
      mu_k_f[k, ] <- mk
      dk <- sweep(Xk, 2L, mk, "-")
      X_c_f[ik, ] <- dk
      W_mat_f <- W_mat_f + crossprod(dk, dk * wk)
    }

    denom_f <- max(sum(w_tr) - w_corr_f, 1.0)
    S_w_f <- W_mat_f / denom_f
    S_reg_f <- compute_ledoit_wolf_cov(S_w_f, X_centered = X_c_f, n_obs = nrow(X_tr), weights = w_tr)
    if (is.list(S_reg_f) && !is.null(S_reg_f$S_reg)) S_reg_f <- S_reg_f$S_reg
    S_inv_f <- tryCatch(solve(S_reg_f), error = function(e) diag(1.0 / pmax(diag(S_reg_f), 1e-6), nrow = p_f, ncol = p_f))

    for (ti in seq_along(test_idx)) {
      obs_i <- test_idx[ti]
      diff_i <- sweep(mu_k_f, 2L, X_te[ti, ], "-")
      d2_i <- rowSums((diff_i %*% S_inv_f) * diff_i)
      ll_i <- -0.5 * d2_i + log(priors)
      pred_cv[obs_i] <- g_levels[which.max(ll_i)]
    }
  }
  round(mean(pred_cv == as.character(y_f), na.rm = TRUE), 4)
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
  for (col in intersect(c("identifiedBy"), names(df))) {
    v <- as.character(df[[col]])
    v[is.na(v)] <- ""
    v <- trimws(v)
    need <- !nzchar(raw_det) & nzchar(v) & !(tolower(v) %in% invalid_raw)
    if (any(need)) raw_det[need] <- v[need]
  }
  out <- rep(NA_character_, n)
  has_raw <- nzchar(raw_det)
  if (!any(has_raw)) return(out)

  raw_sub <- raw_det[has_raw]
  u_raw <- unique(raw_sub)
  s <- tolower(u_raw)
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
  out[has_raw] <- s[match(raw_sub, u_raw)]
  out
}

.compute_loo_pc_likelihoods <- function(df, taxa_vec, det_tier, is_promoted, eval_mask, valid_taxa, prior_vec) {
  n <- nrow(df)
  G <- length(valid_taxa)
  if (G < 2L) return(NULL)

  use_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
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
    p_cap <- max(2L, min(ncol(X_pc), floor(max(sum(eval_pc), 4L) / 4L), 15L))
    X_pc <- X_pc[, seq_len(p_cap), drop = FALSE]
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

.extract_specimen_discordance_signals <- function(
  df,
  eval_mask = NULL,
  taxon_col = "scientificName",
  pc_cols = NULL,
  target_modality = c("all", "morph", "vision", "ecological")
) {
  target_modality <- match.arg(target_modality, c("all", "morph", "vision", "ecological"))
  n <- nrow(df)
  if (is.null(eval_mask) || is.character(eval_mask) || length(eval_mask) != n) {
    eval_mask <- rep(TRUE, n)
  }
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

  if ("exsiccatae_taxonomic_conflict" %in% names(df)) {
    ex_conf <- as.logical(df$exsiccatae_taxonomic_conflict)
    if ("exsiccatae_gold_taxon" %in% names(df)) {
      g_tx <- trimws(as.character(df$exsiccatae_gold_taxon))
      s_tx <- .resolve_authority_taxon_vec(df)
      if (exists("standardize_packera_taxon", mode = "function")) {
        has_g <- !is.na(g_tx) & nzchar(g_tx) & !(tolower(g_tx) %in% c("na", "none", "unknown"))
        if (any(has_g)) g_tx[has_g] <- as.character(standardize_packera_taxon(g_tx[has_g]))
      }
      dt_ex <- .standardize_determiner_tier_vec(df)
      agrees_gold <- !is.na(ex_conf) & ex_conf & !is.na(g_tx) & nzchar(g_tx) & (g_tx == s_tx)
      # Specialist sheet itself suffers zero conflict penalty; elevated Silver duplicate gets mild 0.18 uncertainty
      ex_disagree <- ifelse(agrees_gold, FALSE, ex_conf)
      eps_err <- apply_bool_flag(eps_err, ex_disagree, pos_rate = 0.85, neg_rate = 0.05)
      elev_silv <- eval_mask & agrees_gold & (dt_ex == "Tier_2_Silver")
      if (any(elev_silv)) {
        eps_err[elev_silv] <- ifelse(is.finite(eps_err[elev_silv]), pmax(eps_err[elev_silv], 0.18), 0.18)
      }
    } else {
      eps_err <- apply_bool_flag(eps_err, ex_conf, pos_rate = 0.85, neg_rate = 0.05)
    }
  }
  if (target_modality != "vision" && "is_label_corrupted" %in% names(df)) {
    eps_err <- apply_bool_flag(eps_err, df$is_label_corrupted, pos_rate = 0.85, neg_rate = 0.05)
  }
  if (target_modality != "morph" && "morphometric_outlier_flag" %in% names(df)) {
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
  prior_alpha = 5.0,
  anchor_info = NULL,
  target_modality = c("all", "morph", "vision", "ecological"),
  taxon_col = "species_standardized",
  ...
) {
  mode <- match.arg(mode)
  target_modality <- match.arg(target_modality, c("all", "morph", "vision", "ecological"))
  if (is.null(df) || nrow(df) == 0L) return(numeric(0))
  n <- nrow(df)
  if (is.null(active_mask) && !is.null(anchor_info) && is.list(anchor_info) && !is.null(anchor_info$active_mask)) {
    active_mask <- anchor_info$active_mask
  }
  alpha0 <- max(1e-6, as.numeric(prior_alpha))

  type_info <- .identify_type_specimen_vec(df)
  type_rank <- type_info$type_rank
  is_type_mask <- type_info$is_type_specimen
  det_tier <- .standardize_determiner_tier_vec(df)
  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  det_id <- .normalize_determiner_id_vec(df)

  is_promoted <- rep(FALSE, n)
  if (!is.null(active_mask) && length(active_mask) == n) {
    act_in <- as.logical(active_mask)
    act_in[is.na(act_in)] <- FALSE
    is_promoted <- (det_tier == "Tier_3_Bronze") & act_in
  } else if ("is_promoted_bronze" %in% names(df)) {
    is_promoted <- as.logical(df$is_promoted_bronze)
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
    is_type_mask | (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_promoted
  }

  if (mode == "active_passive" && !any(act_mask)) {
    return(rep(0.0, n))
  }
  eval_mask <- if (mode == "active_passive") act_mask else rep(TRUE, n)

  tier_cell <- ifelse(
    type_rank == 2L, "Tier_0_Primary_Type",
    ifelse(
      type_rank == 1L, "Tier_0_Secondary_Type",
      ifelse(
        det_tier == "Tier_1_Gold", "Tier_1_Gold",
        ifelse(
          det_tier == "Tier_2_Silver", "Tier_2_Silver",
          ifelse(is_promoted, "Tier_3_Bronze_Promoted", "Tier_3_Bronze")
        )
      )
    )
  )

  prior_map <- c(
    "Tier_0_Primary_Type" = 0.999,
    "Tier_0_Secondary_Type" = 0.985,
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
    theta_tier <- pi_vec
    for (tc in unique(tier_cell[eval_mask])) {
      idx_t_all <- which(eval_mask & (tier_cell == tc))
      idx_t_obs <- idx_t_all[obs_mask[idx_t_all]]
      pi_t <- pi_vec[idx_t_all[1L]]
      n_t <- length(idx_t_obs)
      s_t <- if (n_t > 0L) sum(conc_vec[idx_t_obs]) else 0.0
      theta_tier[idx_t_all] <- (alpha0 * pi_t + s_t) / (alpha0 + n_t)
    }

    theta_taxon_tier <- pi_vec
    taxon_tier_cell <- paste0(taxa_vec, "||", tier_cell)
    for (ck in unique(taxon_tier_cell[eval_mask])) {
      idx_all <- which(eval_mask & (taxon_tier_cell == ck))
      idx_obs <- idx_all[obs_mask[idx_all]]
      pi_t <- pi_vec[idx_all[1L]]
      n_c <- length(idx_obs)
      s_c <- if (n_c > 0L) sum(conc_vec[idx_obs]) else 0.0
      theta_taxon_tier[idx_all] <- (alpha0 * pi_t + s_c) / (alpha0 + n_c)
    }

    theta_cell_blend <- 0.50 * theta_tier + 0.50 * theta_taxon_tier
    theta_blend <- theta_cell_blend
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
        lambda_d <- n_d / (0.5 * alpha0 + n_d)
        theta_blend[idx_d_all] <- (1.0 - lambda_d) * theta_cell_blend[idx_d_all] + lambda_d * theta_d_i
      }
    }
    pmax(0.05, pmin(0.995, theta_blend))
  }

  eps_err <- .extract_specimen_discordance_signals(df, eval_mask = eval_mask, target_modality = target_modality)

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
  if ("exsiccatae_gold_taxon" %in% names(df)) {
    p_ex <- safe_std_tx(df$exsiccatae_gold_taxon)
    valid_ex <- !is.na(p_ex) & nzchar(p_ex) & !(tolower(p_ex) %in% c("na", "none", "unknown", "unassigned_packera"))
    if (any(valid_ex)) {
      pred_sources[["exsiccatae"]] <- list(taxon = p_ex, conf = rep(0.92, n), valid = valid_ex)
    }
  }
  if (target_modality != "morph" && "cda_predicted_taxon" %in% names(df)) {
    p_cda <- safe_std_tx(df$cda_predicted_taxon)
    c_cda <- if ("cda_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(0.82, n)
    c_cda[!is.finite(c_cda)] <- 0.82
    pred_sources[["cda"]] <- list(taxon = p_cda, conf = pmax(0.55, pmin(0.92, c_cda)), valid = rep(TRUE, n))
  }
  vis_col <- if ("vision_predicted_label" %in% names(df)) "vision_predicted_label" else if ("predicted_label" %in% names(df)) "predicted_label" else NULL
  if (target_modality != "vision" && !is.null(vis_col)) {
    p_vis <- safe_std_tx(df[[vis_col]])
    pred_sources[["vision"]] <- list(taxon = p_vis, conf = rep(0.82, n), valid = !is_art)
  }
  if (target_modality != "ecological") {
    if ("rf_predicted_taxon" %in% names(df)) {
      p_rf <- safe_std_tx(df$rf_predicted_taxon)
      pred_sources[["rf"]] <- list(taxon = p_rf, conf = rep(0.78, n), valid = rep(TRUE, n))
    }
    if ("edaphic_best_fit_taxon" %in% names(df)) {
      p_ed <- safe_std_tx(df$edaphic_best_fit_taxon)
      pred_sources[["edaphic"]] <- list(taxon = p_ed, conf = rep(0.75, n), valid = rep(TRUE, n))
    }
  }

  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  valid_taxa <- sort(unique(taxa_vec[eval_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_tx)]))
  G <- length(valid_taxa)
  exog_cols <- intersect(
    c("exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon", "vision_predicted_label", "is_label_corrupted", "rf_predicted_taxon", "edaphic_best_fit_taxon"),
    names(df)
  )
  has_exog_signals <- length(exog_cols) > 0L && any(vapply(exog_cols, function(cl) any(!is.na(df[[cl]]) & nzchar(trimws(as.character(df[[cl]])))), logical(1)))
  use_loo_pc <- (G >= 2L) && (target_modality != "morph")
  L_pc <- if (use_loo_pc) .compute_loo_pc_likelihoods(df, taxa_vec, det_tier, is_promoted, eval_mask, valid_taxa, pi_vec) else NULL

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

      has_exog_disc_i <- (is.finite(eps_err[i]) && eps_err[i] > 0.35)
      if (!has_exog_disc_i && length(pred_sources) > 0L) {
        for (ps_chk in pred_sources) {
          if (isTRUE(ps_chk$valid[i]) && !is.na(ps_chk$taxon[i]) && !(ps_chk$taxon[i] %in% bad_tx) && (ps_chk$taxon[i] != valid_taxa[yi])) {
            has_exog_disc_i <- TRUE
            break
          }
        }
      }

      if (!is.null(L_pc) && all(is.finite(L_pc[i, ]))) {
        p_self_pc <- max(0.002, min(0.998, L_pc[i, yi]))
        odds_pc <- (G - 1.0) * p_self_pc / max(1.0 - p_self_pc, 1e-6)
        lo_pc <- log(max(odds_pc, 1e-4))
        # Guard peripheral Tier_1_Gold ecotypes against circular morphology-only down-weighting when exogenous signals do not conflict
        if (!has_exog_disc_i && (det_tier[i] == "Tier_1_Gold") && (lo_pc < 0.0)) {
          lo_pc_use <- max(lo_pc, -0.30)
          wt_pc_use <- 0.45
        } else {
          lo_pc_use <- lo_pc
          wt_pc_use <- if (lo_pc < 0.0) 1.5 else 1.0
        }
        log_odds_vec <- c(log_odds_vec, lo_pc_use)
        wt_vec <- c(wt_vec, wt_pc_use)
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
    if ("institutionCode" %in% names(df) && (G >= 2L || isTRUE(use_loo_pc))) {
      inst_raw <- trimws(as.character(df$institutionCode))
      inst_valid <- !is.na(inst_raw) & nzchar(inst_raw) & !(toupper(inst_raw) %in% c("NA", "NULL", "UNKNOWN"))
      inst_obs <- inst_valid & has_ds_Row & is.finite(conc_final)
      if (any(inst_obs)) {
        inst_sum <- tapply(conc_final[inst_obs], inst_raw[inst_obs], sum, na.rm = TRUE)
        inst_cnt <- tapply(conc_final[inst_obs], inst_raw[inst_obs], length)
        theta_inst_map <- (4.0 * 0.85 + inst_sum) / (4.0 + inst_cnt)
        logit_fn <- function(p_v) log(pmax(1e-4, pmin(1.0 - 1e-4, p_v)) / (1.0 - pmax(1e-4, pmin(1.0 - 1e-4, p_v))))
        inv_logit_fn <- function(z_v) 1.0 / (1.0 + exp(-z_v))
        hit_inst <- which(eval_mask & inst_valid & (inst_raw %in% names(theta_inst_map)))
        if (length(hit_inst) > 0L) {
          th_inst_v <- as.numeric(theta_inst_map[inst_raw[hit_inst]])
          adj_lo <- logit_fn(theta_curr[hit_inst]) + 0.35 * (logit_fn(th_inst_v) - logit_fn(0.85))
          theta_curr[hit_inst] <- pmax(0.05, pmin(0.995, inv_logit_fn(adj_lo)))
        }
      }
    }
    w_eb[eval_mask] <- theta_curr[eval_mask]
    p_prior <- pmax(0.05, pmin(0.995, theta_curr[ds_obs_idx]))
    l_obs <- pmax(0.02, pmin(0.995, conc_final[ds_obs_idx]))
    l_alt <- pmax(0.02, 1.0 - l_obs)
    post_bayes <- (p_prior * l_obs) / pmax(1e-8, p_prior * l_obs + (1.0 - p_prior) * l_alt)
    w_eb[ds_obs_idx] <- pmax(0.05, pmin(0.995, 0.35 * p_prior + 0.65 * post_bayes))

    # Handle any evaluated rows with unknown taxon string that still carry explicit error rates
    rem_err <- which(eval_mask & !has_ds_Row & is.finite(eps_err))
    if (length(rem_err) > 0L) {
      fb_obs <- eval_mask & is.finite(eps_err)
      fb_conc <- pmin(pmax(0.0, 1.0 - eps_err), tier_ceil)
      fb_theta <- blend_tier_and_determiner(fb_conc, fb_obs)
      p_pr_r <- pmax(0.05, pmin(0.995, fb_theta[rem_err]))
      l_ob_r <- pmax(0.02, pmin(0.995, fb_conc[rem_err]))
      l_al_r <- pmax(0.02, 1.0 - l_ob_r)
      pb_r <- (p_pr_r * l_ob_r) / pmax(1e-8, p_pr_r * l_ob_r + (1.0 - p_pr_r) * l_al_r)
      w_eb[rem_err] <- pmax(0.05, pmin(0.995, 0.35 * p_pr_r + 0.65 * pb_r))
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
  if (any(is_type_mask)) {
    w_eb[is_type_mask] <- pi_vec[is_type_mask]
  }

  aw_cfg <- .load_authority_weight_config()
  is_nb_sec <- (type_rank == 1L) & as.logical(type_info$is_name_bearing_type)
  is_nb_sec[is.na(is_nb_sec)] <- FALSE
  base_prec <- ifelse(
    type_rank == 2L, aw_cfg$tier_0_type,
    ifelse(
      is_nb_sec, aw_cfg$tier_0_syntype * (w_eb / 0.985),
      ifelse(
        type_rank == 1L, (aw_cfg$tier_1_gold + 0.15) * (w_eb / 0.985),
        ifelse(
          det_tier == "Tier_1_Gold", aw_cfg$tier_1_gold * (w_eb / 0.95),
          ifelse(
            det_tier == "Tier_2_Silver", aw_cfg$tier_2_silver * (w_eb / 0.85),
            ifelse(
              is_promoted, aw_cfg$tier_3_bronze_promoted * (w_eb / 0.65),
              aw_cfg$tier_3_bronze * (w_eb / 0.35)
            )
          )
        )
      )
    )
  )
  base_prec[!is.finite(base_prec) | base_prec <= 0.0] <- 0.05

  if (mode == "active_passive") {
    w_out <- rep(0.0, n)
    if (any(act_mask)) {
      cap_res <- .apply_stratified_silver_pool_cap(
        w_vec = base_prec,
        taxa_vec = taxa_vec,
        det_tier = det_tier,
        is_type_mask = is_type_mask,
        is_promoted = is_promoted,
        act_mask = act_mask,
        silver_max_share = aw_cfg$silver_max_pool_share,
        min_type_gold = aw_cfg$min_type_gold_for_cap
      )
      w_act <- cap_res$weights[act_mask]
      mean_act <- mean(w_act, na.rm = TRUE)
      w_out[act_mask] <- if (is.finite(mean_act) && mean_act > 1e-12) round(w_act / mean_act, 6) else 1.0
    }
    return(w_out)
  }

  max_w <- max(base_prec, na.rm = TRUE)
  w_cont <- if (is.finite(max_w) && max_w > 1e-12) (3.0 * base_prec / max_w) else rep(1.0, n)
  w_cont[!is.finite(w_cont)] <- 1.0
  if (any(type_rank == 2L)) {
    w_cont[type_rank == 2L] <- 3.0
  }
  w_cont <- pmax(0.5, pmin(3.0, w_cont))
  cap_cont <- .apply_stratified_silver_pool_cap(
    w_vec = w_cont,
    taxa_vec = taxa_vec,
    det_tier = det_tier,
    is_type_mask = is_type_mask,
    is_promoted = is_promoted,
    act_mask = act_mask,
    silver_max_share = aw_cfg$silver_max_pool_share,
    min_type_gold = aw_cfg$min_type_gold_for_cap
  )
  w_final <- cap_cont$weights
  if (any(cap_cont$capped_mask)) {
    w_final[cap_cont$capped_mask] <- pmax(0.01, w_final[cap_cont$capped_mask])
    w_final[!cap_cont$capped_mask] <- pmax(0.5, pmin(3.0, w_final[!cap_cont$capped_mask]))
  }
  round(w_final, 6)
}

.AUTH_WEIGHT_CFG_ENV <- new.env(parent = emptyenv())
.load_authority_weight_config <- function(cfg_path = "config/config.yaml") {
  if (exists("cached", envir = .AUTH_WEIGHT_CFG_ENV, inherits = FALSE)) {
    return(get("cached", envir = .AUTH_WEIGHT_CFG_ENV, inherits = FALSE))
  }
  defs <- list(
    silver_max_pool_share = 0.40,
    min_type_gold_for_cap = 5L,
    tier_0_type = 3.0,
    tier_0_syntype = 2.5,
    tier_1_gold = 2.0,
    tier_2_silver = 0.75,
    tier_3_bronze_promoted = 0.50,
    tier_3_bronze = 0.35
  )
  if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    aw <- if (is.list(cfg) && is.list(cfg$thresholds) && is.list(cfg$thresholds$authority_weights)) {
      cfg$thresholds$authority_weights
    } else if (is.list(cfg) && is.list(cfg$authority_weights)) {
      cfg$authority_weights
    } else {
      NULL
    }
    if (is.list(aw)) {
      if (!is.null(aw$silver_max_pool_share)) defs$silver_max_pool_share <- min(0.85, max(0.05, as.numeric(aw$silver_max_pool_share)))
      if (!is.null(aw$min_type_gold_for_cap)) defs$min_type_gold_for_cap <- max(1L, as.integer(aw$min_type_gold_for_cap))
      if (!is.null(aw$tier_0_type)) defs$tier_0_type <- max(0.5, as.numeric(aw$tier_0_type))
      if (!is.null(aw$tier_0_syntype)) defs$tier_0_syntype <- max(0.5, as.numeric(aw$tier_0_syntype))
      if (!is.null(aw$tier_1_gold)) defs$tier_1_gold <- max(0.2, as.numeric(aw$tier_1_gold))
      if (!is.null(aw$tier_2_silver)) defs$tier_2_silver <- max(0.05, as.numeric(aw$tier_2_silver))
      if (!is.null(aw$tier_3_bronze_promoted)) defs$tier_3_bronze_promoted <- max(0.05, as.numeric(aw$tier_3_bronze_promoted))
      if (!is.null(aw$tier_3_bronze)) defs$tier_3_bronze <- max(0.01, as.numeric(aw$tier_3_bronze))
    }
  }
  assign("cached", defs, envir = .AUTH_WEIGHT_CFG_ENV)
  defs
}

.apply_stratified_silver_pool_cap <- function(
  w_vec,
  taxa_vec,
  det_tier,
  is_type_mask,
  is_promoted = NULL,
  act_mask = NULL,
  silver_max_share = 0.40,
  min_type_gold = 5L
) {
  n <- length(w_vec)
  capped_mask <- rep(FALSE, n)
  if (n == 0L || !is.finite(silver_max_share) || silver_max_share <= 0.0 || silver_max_share >= 1.0) {
    return(list(weights = w_vec, capped_mask = capped_mask))
  }
  if (is.null(act_mask)) act_mask <- rep(TRUE, n)
  if (is.null(is_promoted)) is_promoted <- rep(FALSE, n)
  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  u_tx <- unique(taxa_vec[act_mask & !is.na(taxa_vec) & !(taxa_vec %in% bad_tx)])
  ratio_cap <- silver_max_share / max(1e-6, 1.0 - silver_max_share)
  for (tx in u_tx) {
    idx_tx_act <- which(act_mask & (taxa_vec == tx))
    idx_tg <- idx_tx_act[is_type_mask[idx_tx_act] | (det_tier[idx_tx_act] == "Tier_1_Gold")]
    idx_sv <- idx_tx_act[!is_type_mask[idx_tx_act] & (det_tier[idx_tx_act] == "Tier_2_Silver")]
    if (length(idx_tg) >= as.integer(min_type_gold) && length(idx_sv) >= 1L) {
      w_tg_sum <- sum(w_vec[idx_tg], na.rm = TRUE)
      w_sv_sum <- sum(w_vec[idx_sv], na.rm = TRUE)
      max_sv_sum <- ratio_cap * w_tg_sum
      if (is.finite(w_sv_sum) && is.finite(max_sv_sum) && max_sv_sum > 0.0 && w_sv_sum > (max_sv_sum + 1e-12)) {
        w_vec[idx_sv] <- w_vec[idx_sv] * (max_sv_sum / w_sv_sum)
        capped_mask[idx_sv] <- TRUE
        idx_pr <- idx_tx_act[is_promoted[idx_tx_act]]
        if (length(idx_pr) >= 1L) {
          sv_min <- min(w_vec[idx_sv], na.rm = TRUE)
          w_vec[idx_pr] <- pmin(w_vec[idx_pr], 0.85 * sv_min)
          capped_mask[idx_pr] <- TRUE
        }
      }
    }
  }
  list(weights = w_vec, capped_mask = capped_mask)
}

compute_authority_precision_weights <- function(df, tier_col = "determiner_tier", active_mask = NULL, target_modality = "morph") {
  if (is.logical(tier_col) && is.null(active_mask)) {
    active_mask <- tier_col
    tier_col <- "determiner_tier"
  }
  compute_empirical_bayes_authority_weights(df, mode = "active_passive", active_mask = active_mask, target_modality = target_modality)
}

resolve_ecological_authority_weights <- function(df, allow_bronze_prior = FALSE) {
  n <- if (is.null(df) || !is.data.frame(df)) 0L else nrow(df)
  if (n == 0L) return(numeric(0))
  w_eb <- compute_empirical_bayes_authority_weights(
    df,
    mode = if (isTRUE(allow_bronze_prior)) "continuous_prior_posterior" else "active_passive",
    target_modality = "ecological"
  )
  if (!isTRUE(allow_bronze_prior) && "is_confident_label_noise" %in% names(df)) {
    cl_noise <- as.logical(df$is_confident_label_noise)
    cl_noise[is.na(cl_noise)] <- FALSE
    ti <- .identify_type_specimen_vec(df)
    w_eb[cl_noise & !ti$is_type_specimen] <- 0.0
  }
  round(w_eb, 4)
}

screen_gmm_point_mass_columns <- function(X_mat, ref_mask = NULL, max_tie_frac = 0.25) {
  X_mat <- as.matrix(X_mat)
  cn <- colnames(X_mat)
  if (is.null(cn) || ncol(X_mat) <= 2L) return(seq_len(ncol(X_mat)))
  ref_idx <- if (!is.null(ref_mask) && length(ref_mask) == nrow(X_mat) && sum(ref_mask, na.rm = TRUE) >= 10L) {
    which(as.logical(ref_mask))
  } else {
    seq_len(nrow(X_mat))
  }
  if (length(ref_idx) < 30L) return(cn)
  keep <- vapply(seq_len(ncol(X_mat)), function(j) {
    if (grepl("^PC[0-9]+$", cn[j])) return(TRUE)
    v <- X_mat[ref_idx, j]
    v <- round(v[is.finite(v)], 5)
    if (length(v) < 10L) return(FALSE)
    (max(tabulate(match(v, unique(v)))) / length(v)) <= max_tie_frac
  }, logical(1L))
  if (sum(keep) < 2L) return(cn)
  cn[keep]
}

.max_column_tie_count <- function(X_mat) {
  X_mat <- as.matrix(X_mat)
  if (nrow(X_mat) == 0L || ncol(X_mat) == 0L) return(0L)
  max(vapply(seq_len(ncol(X_mat)), function(j) {
    v <- round(X_mat[is.finite(X_mat[, j]), j], 5)
    if (length(v) == 0L) 0L else max(tabulate(match(v, unique(v))))
  }, integer(1L)))
}

run_bounded_clustvarsel <- function(X_mat, core_idx, max_k = 8L, timeout_sec = 20.0, seed = 42L) {
  if (!requireNamespace("clustvarsel", quietly = TRUE)) return(NULL)
  X_mat <- as.matrix(X_mat)
  sub_idx <- if (length(core_idx) > 350L) {
    core_idx[unique(round(seq(1L, length(core_idx), length.out = 350L)))]
  } else {
    core_idx
  }
  X_sub <- X_mat[sub_idx, , drop = FALSE]
  p_dim <- ncol(X_sub)
  if (nrow(X_sub) < 10L || p_dim < 3L) return(NULL)
  has_ties <- .max_column_tie_count(X_sub) > p_dim
  G_cvs <- seq_len(min(as.integer(max_k), 5L))

  run_cvs_inner <- function() {
    set.seed(as.integer(seed))
    if (has_ties) {
      clustvarsel::clustvarsel(
        X_sub, G = G_cvs, search = "headlong",
        hcModel = "VVI", emModels2 = c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
        itermax = 100, parallel = FALSE, verbose = FALSE
      )
    } else {
      clustvarsel::clustvarsel(X_sub, G = G_cvs, itermax = 100, parallel = FALSE, verbose = FALSE)
    }
  }

  cvs_res <- if (.Platform$OS.type == "unix") {
    jb <- parallel::mcparallel(tryCatch(run_cvs_inner(), error = function(e) NULL), silent = TRUE)
    col_res <- parallel::mccollect(jb, wait = FALSE, timeout = timeout_sec)
    if (is.null(col_res)) {
      tryCatch(tools::pskill(jb$pid, signal = 9L), error = function(e) NULL)
      tryCatch(parallel::mccollect(jb, wait = FALSE, timeout = 0.2), error = function(e) NULL)
      NULL
    } else {
      val <- col_res[[1L]]
      if (inherits(val, "try-error")) NULL else val
    }
  } else {
    tryCatch(run_cvs_inner(), error = function(e) NULL)
  }

  if (is.null(cvs_res) || length(cvs_res$subset) < 2L) return(NULL)
  sub_pos <- as.integer(cvs_res$subset)
  if (has_ties && !is.null(colnames(X_mat))) {
    lead_pc <- which(colnames(X_mat) == "PC1")
    if (length(lead_pc) > 0L && !(lead_pc[1L] %in% sub_pos)) {
      sub_pos <- c(lead_pc[1L], sub_pos)
    }
  }
  cvs_res$subset <- sort(unique(sub_pos))
  cvs_res
}

fit_parallel_mclust_gmm <- function(X_core, G_vec, modelNames = NULL, max_hc_subset = 1500L) {
  if (!requireNamespace("mclust", quietly = TRUE)) return(NULL)
  if (!("package:mclust" %in% search())) {
    tryCatch(suppressPackageStartupMessages(library(mclust, quietly = TRUE, warn.conflicts = FALSE)), error = function(e) NULL)
  }
  X_core <- as.matrix(X_core)
  n_obs <- nrow(X_core)
  p_dim <- ncol(X_core)
  G_int <- as.integer(G_vec)
  has_ties <- .max_column_tie_count(X_core) > p_dim
  ctrl <- mclust::emControl(itmax = c(180L, 180L), tol = c(1e-5, sqrt(.Machine$double.eps)))

  if (n_obs <= 300L && !has_ties) {
    return(mclust::Mclust(X_core, G = G_int, modelNames = modelNames, control = ctrl, verbose = FALSE))
  }

  prior_supported <- if (p_dim == 1L) {
    c("E", "V")
  } else {
    c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV")
  }
  cand_models <- if (!is.null(modelNames) && length(modelNames) > 0L) {
    as.character(modelNames)
  } else {
    as.character(mclust::mclust.options("emModelNames"))
  }
  if (p_dim == 1L) cand_models <- intersect(cand_models, c("E", "V"))
  if (has_ties) cand_models <- intersect(cand_models, prior_supported)
  if (length(cand_models) == 0L) cand_models <- if (p_dim == 1L) c("E", "V") else c("EII", "VII", "EEI", "VEI", "EVI", "VVI")

  prior_arg <- if (has_ties) mclust::priorControl() else NULL
  n_sub <- min(n_obs, max(600L, as.integer(max_hc_subset)))
  init_arg <- if (n_obs > n_sub) {
    ord1 <- order(X_core[, 1L])
    list(subset = sort(ord1[unique(round(seq(1, n_obs, length.out = n_sub)))]))
  } else {
    NULL
  }

  grid_df <- expand.grid(G = G_int, modelName = cand_models, stringsAsFactors = FALSE)
  n_tasks <- nrow(grid_df)
  fit_one_cell <- function(idx) {
    g_i <- grid_df$G[idx]
    mn_i <- grid_df$modelName[idx]
    tryCatch(
      mclust::Mclust(
        X_core, G = g_i, modelNames = mn_i,
        prior = prior_arg, control = ctrl, initialization = init_arg, verbose = FALSE
      ),
      error = function(e) NULL
    )
  }

  use_par <- exists("packera_mclapply", mode = "function") && n_tasks >= 2L
  fits <- if (use_par) {
    n_w <- if (exists("resolve_packera_workers", mode = "function")) {
      resolve_packera_workers(min(32L, n_tasks), n_tasks)
    } else {
      min(8L, n_tasks)
    }
    packera_mclapply(seq_len(n_tasks), fit_one_cell, mc.cores = n_w)
  } else {
    lapply(seq_len(n_tasks), fit_one_cell)
  }

  valid_idx <- which(vapply(fits, function(f) !is.null(f) && !is.null(f$bic) && is.finite(f$bic[1L]), logical(1L)))
  if (length(valid_idx) == 0L) {
    return(mclust::Mclust(
      X_core, G = G_int, modelNames = intersect(cand_models, prior_supported),
      prior = mclust::priorControl(), control = ctrl, initialization = init_arg, verbose = FALSE
    ))
  }

  bic_mat <- matrix(
    NA_real_,
    nrow = length(G_int),
    ncol = length(cand_models),
    dimnames = list(as.character(G_int), cand_models)
  )
  for (idx in valid_idx) {
    r_nm <- as.character(grid_df$G[idx])
    c_nm <- as.character(grid_df$modelName[idx])
    bic_mat[r_nm, c_nm] <- as.numeric(fits[[idx]]$bic[1L])
  }
  bics <- vapply(valid_idx, function(idx) as.numeric(fits[[idx]]$bic[1L]), numeric(1L))
  best_fit <- fits[[valid_idx[which.max(bics)]]]
  best_fit$BIC <- bic_mat
  best_fit
}


.canonicalize_type_basionym_vec <- function(df, fallback_taxa = NULL) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L) return(character(0))
  fb <- if (!is.null(fallback_taxa) && length(fallback_taxa) == n) {
    as.character(fallback_taxa)
  } else if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else {
    rep("Unknown", n)
  }
  raw_v <- rep(NA_character_, n)
  for (sc in intersect(c("species_raw", "scientificName", "taxon", "species_standardized"), names(df))) {
    cand <- trimws(as.character(df[[sc]]))
    hit <- is.na(raw_v) & !is.na(cand) & nzchar(cand) & !(tolower(cand) %in% c("na", "nan", "null", "none", "unknown"))
    raw_v[hit] <- cand[hit]
  }
  raw_v[is.na(raw_v) | !nzchar(raw_v)] <- fb[is.na(raw_v) | !nzchar(raw_v)]
  parse_one <- function(s, fb_one) {
    if (is.na(s) || !nzchar(s)) return(fb_one)
    s_clean <- gsub("\\s+", " ", trimws(s))
    m_infra <- regmatches(s_clean, regexpr("^[A-Z][a-z]+\\s+[a-z][a-z-]+\\s+(var\\.|subsp\\.|ssp\\.|f\\.)\\s+[a-z][a-z-]+", s_clean))
    if (length(m_infra) > 0L && nzchar(m_infra[1L])) return(m_infra[1L])
    m_hyb <- regmatches(s_clean, regexpr("^[A-Z][a-z]+\\s+[x×]\\s*[a-z][a-z-]+", s_clean))
    if (length(m_hyb) > 0L && nzchar(m_hyb[1L])) return(gsub("\\s+", " ", m_hyb[1L]))
    m_bin <- regmatches(s_clean, regexpr("^[A-Z][a-z]+\\s+[a-z][a-z-]+", s_clean))
    if (length(m_bin) > 0L && nzchar(m_bin[1L])) return(m_bin[1L])
    fb_one
  }
  vapply(seq_len(n), function(i) parse_one(raw_v[i], fb[i]), character(1))
}

.extract_terminal_epithet_stem <- function(taxon_vec) {
  vapply(as.character(taxon_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("")
    toks <- strsplit(tolower(trimws(s)), "\\s+")[[1L]]
    toks <- toks[!grepl("^(var\\.?|subsp\\.?|ssp\\.?|f\\.?|x|×)$", toks)]
    if (length(toks) < 2L) return(toks[length(toks)])
    ep <- gsub("[^a-z]", "", toks[length(toks)])
    sub("(us|a|um|is|e)$", "", ep)
  }, character(1), USE.NAMES = FALSE)
}

.is_homotypic_sensu_stricto_basionym <- function(basionym_vec, accepted_vec) {
  b_stem <- .extract_terminal_epithet_stem(basionym_vec)
  a_stem <- .extract_terminal_epithet_stem(accepted_vec)
  nzchar(b_stem) & nzchar(a_stem) & (b_stem == a_stem)
}

.filter_homotypic_type_rows <- function(df, type_rows, accepted_taxon) {
  if (length(type_rows) <= 1L || is.null(df) || nrow(df) == 0L) return(type_rows)
  bas_sub <- .canonicalize_type_basionym_vec(df[type_rows, , drop = FALSE], rep(accepted_taxon, length(type_rows)))
  if (length(unique(bas_sub)) <= 1L) return(type_rows)
  is_homo <- .is_homotypic_sensu_stricto_basionym(bas_sub, rep(accepted_taxon, length(type_rows)))
  if (any(is_homo)) type_rows[is_homo] else type_rows
}

compute_type_specimen_morph_control <- function(
  df,
  X_mat = NULL,
  S_reg = NULL,
  cda_res = NULL,
  gmm_res = NULL,
  flags = NULL,
  triage = NULL,
  reasons = NULL,
  taxon_col = "species_standardized"
) {
  n <- if (is.null(df)) 0L else nrow(df)
  if (n == 0L) {
    return(list(
      type_status = character(0),
      is_type_specimen = logical(0),
      type_rank = integer(0),
      has_empirical_type_anchor = logical(0),
      type_morph_mahal_d2 = numeric(0),
      nearest_type_morph_taxon = character(0),
      nearest_type_basionym = character(0),
      nearest_type_basionym_d2 = numeric(0),
      type_morph_concordance = numeric(0),
      is_typical_morphotype = logical(0),
      min_congener_mahal_d2 = numeric(0),
      min_congener_taxon = character(0),
      open_set_novelty_p = numeric(0),
      open_set_df = integer(0),
      type_centroids = NULL,
      basionym_type_centroids = NULL,
      flags = flags,
      triage = triage,
      reasons = reasons
    ))
  }

  t_info <- .identify_type_specimen_vec(df)
  type_rank <- t_info$type_rank
  is_type_specimen <- t_info$is_type_specimen
  type_status <- if ("type_status" %in% names(df)) {
    ts <- trimws(as.character(df$type_status))
    ifelse(!is.na(ts) & nzchar(ts), ts, t_info$type_status_clean)
  } else {
    t_info$type_status_clean
  }

  if (is.null(X_mat) && !is.null(cda_res) && !is.null(cda_res$X_proj)) {
    X_mat <- cda_res$X_proj
  }
  if (is.null(X_mat)) {
    cand_cols <- if (!is.null(cda_res) && !is.null(cda_res$selected_features)) {
      intersect(as.character(cda_res$selected_features), names(df))
    } else {
      character(0)
    }
    if (length(cand_cols) == 0L) {
      cand_cols <- .discover_numeric_sorted_cols(names(df), "^PC[0-9]+$")
    }
    if (length(cand_cols) > 0L) {
      X_mat <- as.matrix(df[, cand_cols, drop = FALSE])
    } else if (!is.null(cda_res) && !is.null(cda_res$canonical_scores)) {
      X_mat <- as.matrix(cda_res$canonical_scores)
    } else {
      X_mat <- matrix(0.0, nrow = n, ncol = 2L)
    }
  }
  X_mat <- as.matrix(X_mat)
  storage.mode(X_mat) <- "double"
  p <- ncol(X_mat)
  for (j in seq_len(p)) {
    bad_j <- !is.finite(X_mat[, j])
    if (any(bad_j)) {
      m_j <- mean(X_mat[!bad_j, j], na.rm = TRUE)
      X_mat[bad_j, j] <- if (is.finite(m_j)) m_j else 0.0
    }
  }

  if (is.null(S_reg) && !is.null(cda_res) && !is.null(cda_res$S_reg)) {
    S_reg <- cda_res$S_reg
  }
  if (is.null(S_reg) || !is.matrix(S_reg) || nrow(S_reg) != p || ncol(S_reg) != p) {
    S_raw <- if (n >= 3L) stats::cov(X_mat) else diag(1.0, nrow = p, ncol = p)
    S_raw[!is.finite(S_raw)] <- 0.0
    S_reg <- compute_ledoit_wolf_cov(S_raw + diag(1e-5, nrow = p, ncol = p), X_centered = scale(X_mat, center = TRUE, scale = FALSE), n_obs = n)
  }
  S_inv <- tryCatch(
    solve(S_reg),
    error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-5), nrow = p, ncol = p)
  )

  taxa_vec <- .resolve_authority_taxon_vec(df, taxon_col = taxon_col)
  if (all(taxa_vec == "Unassigned_Packera") && !is.null(cda_res) && !is.null(cda_res$grouping_col)) {
    taxa_vec <- as.character(cda_res$grouping_col)
  }
  det_tier <- .standardize_determiner_tier_vec(df)
  is_t1_pure <- .identify_tier1_pure_mask(df)
  if (sum(is_t1_pure) == 0L) is_t1_pure <- rep(TRUE, n)
  auth_w <- if ("authority_weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$authority_weight))
  } else if (!is.null(cda_res) && !is.null(cda_res$authority_weight)) {
    suppressWarnings(as.numeric(cda_res$authority_weight))
  } else {
    rep(1.0, n)
  }
  auth_w[!is.finite(auth_w) | auth_w <= 0.0] <- 1.0

  bad_tx <- c("", "Unknown", "Unassigned_Packera", "NA")
  active_taxa <- sort(unique(taxa_vec[!is.na(taxa_vec) & !(taxa_vec %in% bad_tx)]))
  if (length(active_taxa) == 0L) {
    active_taxa <- sort(unique(stats::na.omit(taxa_vec)))
  }
  G <- length(active_taxa)

  gmm_cl <- if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n) {
    as.integer(gmm_res$classification)
  } else if ("gmm_cluster" %in% names(df)) {
    suppressWarnings(as.integer(df$gmm_cluster))
  } else {
    rep(1L, n)
  }
  u_clusters <- sort(unique(stats::na.omit(gmm_cl)))
  cl_centroids <- if (length(u_clusters) > 0L) {
    cent_mat <- matrix(0.0, nrow = length(u_clusters), ncol = p, dimnames = list(as.character(u_clusters), colnames(X_mat)))
    for (ci in seq_along(u_clusters)) {
      c_val <- u_clusters[ci]
      c_rows <- which(!is.na(gmm_cl) & (gmm_cl == c_val))
      cent_mat[ci, ] <- colMeans(X_mat[c_rows, , drop = FALSE])
    }
    cent_mat
  } else {
    NULL
  }

  basionym_vec_all <- .canonicalize_type_basionym_vec(df, taxa_vec)
  type_centroids <- matrix(0.0, nrow = G, ncol = p, dimnames = list(active_taxa, colnames(X_mat)))
  has_emp_map <- stats::setNames(rep(FALSE, G), active_taxa)
  typical_cl_map <- stats::setNames(rep(NA_integer_, G), active_taxa)
  typical_cl_list <- stats::setNames(vector("list", G), active_taxa)
  n_eff_type_vec <- stats::setNames(numeric(G), active_taxa)
  type_anchor_rows_list <- stats::setNames(vector("list", G), active_taxa)
  basionym_entries <- list()
  is_noise_vec <- if (!is.null(df) && "is_confident_label_noise" %in% names(df)) {
    as.logical(df$is_confident_label_noise)
  } else {
    rep(FALSE, n)
  }
  is_noise_vec[is.na(is_noise_vec)] <- FALSE

  compute_pool_centroid <- function(rows_sub) {
    prim_s <- rows_sub[type_rank[rows_sub] == 2L]
    pool_s <- if (length(rows_sub) >= 3L && length(prim_s) > 0L) prim_s else rows_sub
    mu_s <- if (length(rows_sub) == 1L) {
      as.numeric(X_mat[rows_sub, ])
    } else if (length(rows_sub) == 2L) {
      w_t <- ifelse(type_rank[rows_sub] == 2L, 2.0, 1.0)
      as.numeric(colSums(X_mat[rows_sub, , drop = FALSE] * w_t) / sum(w_t))
    } else {
      .compute_robust_taxon_medoid(X_mat[pool_s, , drop = FALSE], S_inv = S_inv)
    }
    list(mu = mu_s, pool = pool_s, n_eff = .count_effective_type_gatherings(df, pool_s))
  }

  for (k in seq_len(G)) {
    tx_k <- active_taxa[k]
    tx_rows <- which(taxa_vec == tx_k)
    type_rows_k <- tx_rows[is_type_specimen[tx_rows] & is_t1_pure[tx_rows]]
    if (length(type_rows_k) == 0L) {
      type_rows_k <- tx_rows[is_type_specimen[tx_rows]]
    }
    ss_type_rows_k <- type_rows_k
    if (length(type_rows_k) > 0L) {
      has_emp_map[tx_k] <- TRUE
      bas_k <- basionym_vec_all[type_rows_k]
      u_bas_k <- unique(bas_k[!is.na(bas_k) & nzchar(bas_k)])
      if (length(u_bas_k) == 0L) u_bas_k <- tx_k
      bas_k_list <- lapply(u_bas_k, function(bn) {
        b_rows <- type_rows_k[bas_k == bn]
        if (length(b_rows) == 0L) b_rows <- type_rows_k
        cp <- compute_pool_centroid(b_rows)
        is_homo <- isTRUE(.is_homotypic_sensu_stricto_basionym(bn, tx_k)[1L])
        list(
          taxon = tx_k,
          basionym = bn,
          is_homotypic = is_homo,
          rows = b_rows,
          pool = cp$pool,
          n_eff = cp$n_eff,
          max_rank = max(type_rank[b_rows], na.rm = TRUE),
          mu = cp$mu
        )
      })
      for (be in bas_k_list) {
        basionym_entries[[length(basionym_entries) + 1L]] <- be
      }
      if (length(bas_k_list) == 1L) {
        mu_k <- bas_k_list[[1L]]$mu
        type_anchor_rows_list[[tx_k]] <- bas_k_list[[1L]]$pool
        n_eff_type_vec[tx_k] <- bas_k_list[[1L]]$n_eff
        ss_type_rows_k <- bas_k_list[[1L]]$rows
      } else {
        homo_idx <- which(vapply(bas_k_list, function(x) isTRUE(x$is_homotypic), logical(1L)))
        if (length(homo_idx) >= 1L) {
          ss_rows <- unlist(lapply(bas_k_list[homo_idx], `[[`, "rows"), use.names = FALSE)
          cp_ss <- compute_pool_centroid(ss_rows)
          mu_k <- cp_ss$mu
          type_anchor_rows_list[[tx_k]] <- cp_ss$pool
          n_eff_type_vec[tx_k] <- cp_ss$n_eff
          ss_type_rows_k <- ss_rows
        } else {
          mu_mat_b <- do.call(rbind, lapply(bas_k_list, `[[`, "mu"))
          max_pair_d2 <- 0.0
          for (bi in seq_len(nrow(mu_mat_b) - 1L)) {
            for (bj in (bi + 1L):nrow(mu_mat_b)) {
              dv_b <- mu_mat_b[bi, ] - mu_mat_b[bj, ]
              max_pair_d2 <- max(max_pair_d2, as.numeric(t(dv_b) %*% S_inv %*% dv_b))
            }
          }
          if (max_pair_d2 > stats::qchisq(0.95, df = max(1L, p))) {
            ord_b <- order(
              vapply(bas_k_list, `[[`, integer(1L), "max_rank"),
              vapply(bas_k_list, `[[`, numeric(1L), "n_eff"),
              decreasing = TRUE
            )
            best_be <- bas_k_list[[ord_b[1L]]]
            mu_k <- best_be$mu
            type_anchor_rows_list[[tx_k]] <- best_be$pool
            n_eff_type_vec[tx_k] <- best_be$n_eff
            ss_type_rows_k <- best_be$rows
          } else {
            cp_all <- compute_pool_centroid(type_rows_k)
            mu_k <- cp_all$mu
            type_anchor_rows_list[[tx_k]] <- cp_all$pool
            n_eff_type_vec[tx_k] <- cp_all$n_eff
            ss_type_rows_k <- type_rows_k
          }
        }
      }
    } else {
      has_emp_map[tx_k] <- FALSE
      type_anchor_rows_list[[tx_k]] <- integer(0)
      gold_rows_k <- tx_rows[(det_tier[tx_rows] == "Tier_1_Gold") & is_t1_pure[tx_rows] & !is_noise_vec[tx_rows]]
      if (length(gold_rows_k) == 0L) {
        gold_rows_k <- tx_rows[(det_tier[tx_rows] == "Tier_1_Gold") & !is_noise_vec[tx_rows]]
      }
      if (length(gold_rows_k) > 0L) {
        max_w_k <- max(auth_w[gold_rows_k], na.rm = TRUE)
        top_gold_k <- gold_rows_k[auth_w[gold_rows_k] >= 0.85 * max_w_k]
        n_eff_type_vec[tx_k] <- as.numeric(max(1L, length(top_gold_k)))
        mu_k <- .compute_robust_taxon_medoid(X_mat[top_gold_k, , drop = FALSE], S_inv = S_inv)
      } else {
        n_eff_type_vec[tx_k] <- as.numeric(max(1L, length(tx_rows)))
        mu_k <- .compute_robust_taxon_medoid(X_mat[tx_rows, , drop = FALSE], S_inv = S_inv)
      }
      basionym_entries[[length(basionym_entries) + 1L]] <- list(
        taxon = tx_k,
        basionym = tx_k,
        is_homotypic = TRUE,
        rows = integer(0),
        pool = integer(0),
        n_eff = n_eff_type_vec[tx_k],
        max_rank = 0L,
        mu = mu_k
      )
    }
    type_centroids[k, ] <- mu_k

    # Identify per-species typical GMM morphotype cluster(s) for taxon k
    if (length(u_clusters) > 0L) {
      valid_cl_rows <- tx_rows[!is.na(gmm_cl[tx_rows]) & !is_noise_vec[tx_rows]]
      if (length(valid_cl_rows) == 0L) {
        valid_cl_rows <- tx_rows[!is.na(gmm_cl[tx_rows])]
      }
      gold_type_cl_rows <- valid_cl_rows[is_type_specimen[valid_cl_rows] | (det_tier[valid_cl_rows] == "Tier_1_Gold")]
      if (length(gold_type_cl_rows) >= 3L) {
        d2_gt <- rowSums((sweep(X_mat[gold_type_cl_rows, , drop = FALSE], 2L, mu_k, "-") %*% S_inv) *
                         sweep(X_mat[gold_type_cl_rows, , drop = FALSE], 2L, mu_k, "-"))
        inlier_gt <- gold_type_cl_rows[is_type_specimen[gold_type_cl_rows] | (is.finite(d2_gt) & d2_gt <= stats::qchisq(0.975, df = max(1L, p)))]
        if (length(inlier_gt) >= 2L) gold_type_cl_rows <- inlier_gt
      }
      ref_k <- if (length(gold_type_cl_rows) >= 2L) gold_type_cl_rows else valid_cl_rows
      typ_k <- integer(0)
      if (length(ref_k) > 0L) {
        cl_mass <- tapply(auth_w[ref_k], gmm_cl[ref_k], sum, na.rm = TRUE)
        cl_mass[is.na(cl_mass)] <- 0.0
        tot_mass <- sum(cl_mass, na.rm = TRUE)
        if (tot_mass > 0.0) {
          max_mass <- max(cl_mass, na.rm = TRUE)
          top_cls <- as.integer(names(cl_mass[cl_mass >= max_mass - 1e-9]))
          if (length(top_cls) == 1L) {
            modal_cl <- top_cls[1L]
          } else {
            d2_top <- vapply(top_cls, function(cv) {
              dv <- cl_centroids[as.character(cv), ] - mu_k
              as.numeric(t(dv) %*% S_inv %*% dv)
            }, numeric(1L))
            modal_cl <- top_cls[which.min(d2_top)]
          }
          sec_cls <- as.integer(names(cl_mass[cl_mass >= 0.30 * tot_mass]))
          typ_k <- unique(c(modal_cl, sec_cls))
        }
      }
      type_cl_vals <- as.integer(stats::na.omit(gmm_cl[ss_type_rows_k]))
      if (length(type_cl_vals) > 0L) {
        if (length(ss_type_rows_k) < length(type_rows_k)) {
          typ_k <- unique(type_cl_vals)
        } else {
          typ_k <- unique(c(type_cl_vals, typ_k))
        }
      }
      if (length(typ_k) == 0L) {
        tx_cls <- sort(unique(stats::na.omit(gmm_cl[tx_rows])))
        cand_cls <- if (length(tx_cls) > 0L) tx_cls else u_clusters
        d2_cls <- vapply(cand_cls, function(cv) {
          dv <- cl_centroids[as.character(cv), ] - mu_k
          as.numeric(t(dv) %*% S_inv %*% dv)
        }, numeric(1L))
        typ_k <- cand_cls[which.min(d2_cls)]
      }
      typical_cl_map[tx_k] <- typ_k[1L]
      typical_cl_list[[tx_k]] <- typ_k
    }
  }

  n_bas <- length(basionym_entries)
  bas_keys <- if (n_bas > 0L) {
    raw_bk <- vapply(basionym_entries, `[[`, character(1L), "basionym")
    make.unique(raw_bk, sep = "__")
  } else {
    character(0)
  }
  basionym_type_centroids <- if (n_bas > 0L) {
    b_mat <- do.call(rbind, lapply(basionym_entries, `[[`, "mu"))
    rownames(b_mat) <- bas_keys
    colnames(b_mat) <- colnames(X_mat)
    attr(b_mat, "accepted_taxon") <- vapply(basionym_entries, `[[`, character(1L), "taxon")
    attr(b_mat, "basionym_name") <- vapply(basionym_entries, `[[`, character(1L), "basionym")
    attr(b_mat, "is_homotypic") <- vapply(basionym_entries, `[[`, logical(1L), "is_homotypic")
    b_mat
  } else {
    type_centroids
  }

  d2_bas_mat <- if (n_bas > 0L) matrix(NA_real_, nrow = n, ncol = n_bas, dimnames = list(NULL, bas_keys)) else NULL
  if (n_bas > 0L) {
    for (b_idx in seq_len(n_bas)) {
      be <- basionym_entries[[b_idx]]
      diff_b <- sweep(X_mat, 2L, be$mu, "-")
      raw_b <- rowSums((diff_b %*% S_inv) * diff_b)
      raw_b[!is.finite(raw_b) | raw_b < 0.0] <- 0.0
      scale_b <- 1.0 + (1.0 / max(1.0, be$n_eff))
      sc_b <- raw_b / scale_b
      if (length(be$pool) > 0L) sc_b[be$pool] <- raw_b[be$pool]
      d2_bas_mat[, b_idx] <- sc_b
    }
  }

  emp_type_all <- which(is_type_specimen)
  if (any(grepl("^Morphotype_K", active_taxa)) && length(emp_type_all) > 0L) {
    u_emp_bas <- unique(basionym_vec_all[emp_type_all])
    emp_bas_d2 <- matrix(NA_real_, nrow = n, ncol = length(u_emp_bas), dimnames = list(NULL, u_emp_bas))
    emp_bas_mu <- matrix(0.0, nrow = length(u_emp_bas), ncol = p, dimnames = list(u_emp_bas, colnames(X_mat)))
    for (ebi in seq_along(u_emp_bas)) {
      eb_nm <- u_emp_bas[ebi]
      eb_rows <- emp_type_all[basionym_vec_all[emp_type_all] == eb_nm]
      cp_eb <- compute_pool_centroid(eb_rows)
      emp_bas_mu[ebi, ] <- cp_eb$mu
      diff_eb <- sweep(X_mat, 2L, cp_eb$mu, "-")
      raw_eb <- rowSums((diff_eb %*% S_inv) * diff_eb)
      raw_eb[!is.finite(raw_eb) | raw_eb < 0.0] <- 0.0
      sc_eb <- raw_eb / (1.0 + (1.0 / max(1.0, cp_eb$n_eff)))
      if (length(cp_eb$pool) > 0L) sc_eb[cp_eb$pool] <- raw_eb[cp_eb$pool]
      emp_bas_d2[, ebi] <- sc_eb
    }
    basionym_type_centroids <- emp_bas_mu
    eb_near_idx <- max.col(-emp_bas_d2, ties.method = "first")
    nearest_type_basionym <- u_emp_bas[eb_near_idx]
    nearest_type_basionym_d2 <- round(emp_bas_d2[cbind(seq_len(n), eb_near_idx)], 4)
  } else if (!is.null(d2_bas_mat) && ncol(d2_bas_mat) >= 1L) {
    b_near_idx <- max.col(-d2_bas_mat, ties.method = "first")
    bas_names_clean <- vapply(basionym_entries, `[[`, character(1L), "basionym")
    nearest_type_basionym <- bas_names_clean[b_near_idx]
    nearest_type_basionym_d2 <- round(d2_bas_mat[cbind(seq_len(n), b_near_idx)], 4)
  } else {
    nearest_type_basionym <- rep(NA_character_, n)
    nearest_type_basionym_d2 <- rep(NA_real_, n)
  }

  d2_mat <- matrix(NA_real_, nrow = n, ncol = G, dimnames = list(NULL, active_taxa))
  for (k in seq_len(G)) {
    tx_k <- active_taxa[k]
    diff_k <- sweep(X_mat, 2L, type_centroids[k, ], "-")
    raw_d2 <- rowSums((diff_k %*% S_inv) * diff_k)
    raw_d2[!is.finite(raw_d2) | raw_d2 < 0.0] <- 0.0
    scale_k <- 1.0 + (1.0 / max(1.0, n_eff_type_vec[tx_k]))
    scaled_d2 <- raw_d2 / scale_k
    anchor_rows <- type_anchor_rows_list[[tx_k]]
    if (length(anchor_rows) > 0L) {
      scaled_d2[anchor_rows] <- raw_d2[anchor_rows]
    }
    d2_mat[, k] <- scaled_d2
  }

  given_idx <- match(taxa_vec, active_taxa)
  nearest_idx <- if (G >= 1L) max.col(-d2_mat, ties.method = "first") else rep(NA_integer_, n)
  nearest_type_morph_taxon <- if (G >= 1L) active_taxa[nearest_idx] else rep(NA_character_, n)

  type_morph_mahal_d2 <- rep(NA_real_, n)
  valid_g <- which(!is.na(given_idx))
  if (length(valid_g) > 0L) {
    type_morph_mahal_d2[valid_g] <- round(d2_mat[cbind(valid_g, given_idx[valid_g])], 4)
  }

  type_morph_concordance <- rep(NA_real_, n)
  if (G >= 2L && length(valid_g) > 0L) {
    log_lik <- -0.5 * d2_mat
    max_ll <- apply(log_lik, 1L, max, na.rm = TRUE)
    lik_mat <- exp(log_lik - max_ll)
    lik_mat <- lik_mat / pmax(rowSums(lik_mat, na.rm = TRUE), 1e-12)
    type_morph_concordance[valid_g] <- round(lik_mat[cbind(valid_g, given_idx[valid_g])], 4)
  } else if (G == 1L && length(valid_g) > 0L) {
    type_morph_concordance[valid_g] <- round(1.0 - stats::pchisq(type_morph_mahal_d2[valid_g], df = max(1L, p)), 4)
  }

  has_empirical_type_anchor <- ifelse(!is.na(given_idx), unname(has_emp_map[active_taxa[given_idx]]), FALSE)
  is_typical_morphotype <- rep(NA, n)
  chisq_typ_90 <- stats::qchisq(0.90, df = max(1L, p))
  if (length(valid_g) > 0L) {
    for (ii in valid_g) {
      tx_ii <- active_taxa[given_idx[ii]]
      typ_cls_ii <- typical_cl_list[[tx_ii]]
      has_cl_ii <- !is.na(gmm_cl[ii])
      has_d2_ii <- is.finite(type_morph_mahal_d2[ii])
      in_typ_cl <- has_cl_ii && (length(typ_cls_ii) > 0L) && (gmm_cl[ii] %in% typ_cls_ii)
      in_d2_90 <- has_d2_ii && (type_morph_mahal_d2[ii] <= chisq_typ_90)
      if (in_typ_cl || in_d2_90) {
        is_typical_morphotype[ii] <- TRUE
      } else if (has_cl_ii || has_d2_ii) {
        is_typical_morphotype[ii] <- FALSE
      }
    }
  }

  if (!is.null(flags) && !is.null(triage) && !is.null(reasons) && any(is_type_specimen)) {
    grp_vec <- if (!is.null(cda_res) && !is.null(cda_res$grouping_col)) as.character(cda_res$grouping_col) else taxa_vec
    pred_vec <- if (!is.null(cda_res) && !is.null(cda_res$predicted_taxon)) as.character(cda_res$predicted_taxon) else taxa_vec
    post_vec <- if (!is.null(cda_res) && !is.null(cda_res$posterior_prob)) suppressWarnings(as.numeric(cda_res$posterior_prob)) else rep(1.0, n)
    for (i in which(is_type_specimen)) {
      flags[i] <- FALSE
      r_sp <- if (!is.na(grp_vec[i]) && nzchar(grp_vec[i])) grp_vec[i] else taxa_vec[i]
      p_sp <- if (!is.na(pred_vec[i]) && nzchar(pred_vec[i])) pred_vec[i] else r_sp
      pp <- if (is.finite(post_vec[i])) post_vec[i] else 0.0
      if (r_sp != p_sp && pp >= 0.85) {
        triage[i] <- "CRITICAL"
        reasons[i] <- sprintf("Type_Specimen_Taxonomic_Overlap_with_%s", gsub(" ", "_", p_sp))
      } else {
        triage[i] <- "CLEAN"
        reasons[i] <- "Verified_Tier1_Gold_Anchor"
      }
    }
  }

  osr_res <- .compute_open_set_recognition_metrics(df, X_mat, S_inv, d2_mat, active_taxa, taxa_vec, gmm_cl)

  list(
    type_status = as.character(type_status),
    is_type_specimen = as.logical(is_type_specimen),
    type_rank = as.integer(type_rank),
    has_empirical_type_anchor = as.logical(has_empirical_type_anchor),
    type_morph_mahal_d2 = as.numeric(type_morph_mahal_d2),
    nearest_type_morph_taxon = as.character(nearest_type_morph_taxon),
    nearest_type_basionym = as.character(nearest_type_basionym),
    nearest_type_basionym_d2 = as.numeric(nearest_type_basionym_d2),
    type_morph_concordance = as.numeric(type_morph_concordance),
    is_typical_morphotype = as.logical(is_typical_morphotype),
    min_congener_mahal_d2 = as.numeric(osr_res$min_congener_mahal_d2),
    min_congener_taxon = as.character(osr_res$min_congener_taxon),
    open_set_novelty_p = as.numeric(osr_res$open_set_novelty_p),
    open_set_df = as.integer(osr_res$open_set_df),
    type_centroids = type_centroids,
    basionym_type_centroids = basionym_type_centroids,
    flags = flags,
    triage = triage,
    reasons = reasons
  )
}

.compute_open_set_recognition_metrics <- function(df, X_mat, S_inv, d2_mat, active_taxa, taxa_vec, gmm_cl) {
  n <- nrow(X_mat)
  p <- max(1L, ncol(X_mat))
  min_congener_mahal_d2 <- rep(NA_real_, n)
  min_congener_taxon <- rep(NA_character_, n)
  open_set_novelty_p <- rep(NA_real_, n)
  open_set_df <- rep(as.integer(p), n)
  if (n == 0L) {
    return(list(min_congener_mahal_d2 = min_congener_mahal_d2, min_congener_taxon = min_congener_taxon, open_set_novelty_p = open_set_novelty_p, open_set_df = integer(0)))
  }

  sp_raw <- if (!is.null(df) && "species_standardized" %in% names(df)) as.character(df$species_standardized) else as.character(taxa_vec)
  sp_taxa <- active_taxa[grepl("^Packera\\s+", active_taxa)]
  if (length(sp_taxa) >= 2L && !is.null(d2_mat) && ncol(d2_mat) >= 2L) {
    sp_cols <- match(sp_taxa, active_taxa)
    d2_sp <- d2_mat[, sp_cols, drop = FALSE]
    min_all <- apply(d2_sp, 1L, min, na.rm = TRUE)
    min_all[!is.finite(min_all) | min_all < 0.0] <- 0.0
    open_set_novelty_p <- round(1.0 - stats::pchisq(min_all, df = p), 6)
    for (i in seq_len(n)) {
      yi_sp <- if (!is.na(sp_raw[i]) && grepl("^Packera\\s+", sp_raw[i])) sp_raw[i] else taxa_vec[i]
      oth_idx <- which(sp_taxa != yi_sp)
      if (length(oth_idx) >= 1L) {
        d_sub <- d2_sp[i, oth_idx]
        if (any(is.finite(d_sub))) {
          best_o <- oth_idx[which.min(d_sub)]
          min_congener_mahal_d2[i] <- round(as.numeric(d2_sp[i, best_o]), 4)
          min_congener_taxon[i] <- sp_taxa[best_o]
        }
      }
    }
    return(list(min_congener_mahal_d2 = min_congener_mahal_d2, min_congener_taxon = min_congener_taxon, open_set_novelty_p = open_set_novelty_p, open_set_df = open_set_df))
  }

  # Single-taxon intraspecific mode: query canonical SQLite for congeneric Packera centroids
  # in run-invariant standardized harmonic + vestiture + topology space,
  # falling back strictly to same-run_id PC1..3 when harmonics are absent.
  db_osr <- tryCatch({
    if (!is.null(df) && exists("connect_packera_db", mode = "function") && requireNamespace("DBI", quietly = TRUE)) {
      con <- connect_packera_db(read_only = TRUE)
      on.exit(tryCatch(DBI::dbDisconnect(con), error = function(e) NULL), add = TRUE)
      if (is.null(con)) return(NULL)
      inv_cands <- intersect(
        c(paste0("A", 1:6), paste0("D", 1:6), "solidity", "dissection_index", "blade_cie_l_norm", "ECT_PC1", "ECT_PC2", "ect_lobe_complexity"),
        names(df)
      )
      inv_cands <- inv_cands[vapply(inv_cands, function(cl) is.numeric(df[[cl]]) && mean(is.finite(df[[cl]])) >= 0.50, logical(1))]
      raw_pc_cands <- intersect(c("PC1", "PC2", "PC3"), names(df))
      raw_pc_cands <- raw_pc_cands[vapply(raw_pc_cands, function(cl) is.numeric(df[[cl]]) && mean(is.finite(df[[cl]])) >= 0.80, logical(1))]
      active_rid <- if (exists("resolve_packera_run_id", mode = "function")) tryCatch(as.character(resolve_packera_run_id(NULL)), error = function(e) "") else ""

      eval_db_space <- function(cand_cols, require_same_run = FALSE) {
        if (length(cand_cols) < 2L) return(NULL)
        db_df <- NULL
        if (DBI::dbExistsTable(con, "otu_morphometrics")) {
          flds <- DBI::dbListFields(con, "otu_morphometrics")
          use_c <- intersect(cand_cols, flds)
          if (length(use_c) >= 2L && "species_standardized" %in% flds) {
            q_cols <- paste(c("species_standardized", use_c), collapse = ", ")
            where_parts <- c("species_standardized LIKE 'Packera %'")
            if ("authority_weight" %in% flds) where_parts <- c(where_parts, "COALESCE(authority_weight, 1.0) > 0")
            if ("morph_tier_source" %in% flds) where_parts <- c(where_parts, "COALESCE(morph_tier_source, 'Basal_Rosette') != 'Cauline_Only_Fallback'")
            if (isTRUE(require_same_run) && "run_id" %in% flds) {
              if (!nzchar(active_rid)) return(NULL)
              sql <- sprintf("SELECT %s FROM otu_morphometrics WHERE %s AND run_id = $1", q_cols, paste(where_parts, collapse = " AND "))
              db_df <- DBI::dbGetQuery(con, sql, params = list(active_rid))
            } else if (!isTRUE(require_same_run)) {
              sql <- sprintf("SELECT %s FROM otu_morphometrics WHERE %s", q_cols, paste(where_parts, collapse = " AND "))
              db_df <- DBI::dbGetQuery(con, sql)
            }
          }
        }
        if ((is.null(db_df) || length(unique(db_df$species_standardized)) < 2L) && !isTRUE(require_same_run) &&
            DBI::dbExistsTable(con, "leaf_efa_harmonics") && (DBI::dbExistsTable(con, "specimens") || DBI::dbExistsTable(con, "vouchers"))) {
          v_tbl <- if (DBI::dbExistsTable(con, "specimens")) "specimens" else "vouchers"
          flds_h <- DBI::dbListFields(con, "leaf_efa_harmonics")
          flds_v <- DBI::dbListFields(con, v_tbl)
          use_c <- intersect(cand_cols, flds_h)
          if (length(use_c) >= 2L && "species_standardized" %in% flds_v) {
            q_cols <- paste(paste0("h.", use_c), collapse = ", ")
            lz_clause <- if ("leaf_zone" %in% flds_h) " AND COALESCE(LOWER(h.leaf_zone), 'basal') != 'cauline'" else ""
            db_df <- DBI::dbGetQuery(con, sprintf(
              "SELECT v.species_standardized, %s FROM leaf_efa_harmonics h JOIN %s v ON h.catalogNumber = v.catalogNumber WHERE v.species_standardized LIKE 'Packera %%'%s",
              q_cols, v_tbl, lz_clause
            ))
          }
        }
        if (is.null(db_df) || nrow(db_df) < 6L) return(NULL)
        use_c <- intersect(cand_cols, names(db_df))
        if (length(use_c) < 2L) return(NULL)
        for (uc in use_c) db_df[[uc]] <- suppressWarnings(as.numeric(db_df[[uc]]))
        cc <- stats::complete.cases(db_df[, use_c, drop = FALSE])
        db_df <- db_df[cc & !is.na(db_df$species_standardized) & grepl("^Packera\\s+", db_df$species_standardized), , drop = FALSE]
        focal_sp <- unique(sp_raw[!is.na(sp_raw) & grepl("^Packera\\s+", sp_raw)])
        oth_sp <- setdiff(unique(db_df$species_standardized), focal_sp)
        if (length(oth_sp) < 2L) return(NULL)
        X_loc_raw <- as.matrix(df[, use_c, drop = FALSE]); storage.mode(X_loc_raw) <- "double"
        X_db_raw <- as.matrix(db_df[, use_c, drop = FALSE]); storage.mode(X_db_raw) <- "double"
        comb_sd <- vapply(seq_len(ncol(X_db_raw)), function(jc) stats::sd(c(X_loc_raw[, jc], X_db_raw[, jc]), na.rm = TRUE), numeric(1))
        keep_v <- which(is.finite(comb_sd) & comb_sd > 1e-6)
        if (length(keep_v) < 2L) return(NULL)
        use_c <- use_c[keep_v]
        sc_vec <- comb_sd[keep_v]
        X_loc_s <- sweep(X_loc_raw[, keep_v, drop = FALSE], 2L, sc_vec, "/")
        X_db_s <- sweep(X_db_raw[, keep_v, drop = FALSE], 2L, sc_vec, "/")
        X_loc_imp <- X_loc_s
        for (jc in seq_len(ncol(X_loc_imp))) {
          bad_c <- !is.finite(X_loc_imp[, jc])
          if (any(bad_c)) X_loc_imp[bad_c, jc] <- stats::median(c(X_loc_imp[!bad_c, jc], X_db_s[, jc]), na.rm = TRUE)
        }
        all_sp <- unique(c(focal_sp, oth_sp))
        cent_db <- matrix(0.0, nrow = length(all_sp), ncol = length(use_c), dimnames = list(all_sp, use_c))
        resid_list <- vector("list", length(all_sp))
        for (si in seq_along(all_sp)) {
          sp_nm <- all_sp[si]
          sub_m <- if (sp_nm %in% db_df$species_standardized) X_db_s[db_df$species_standardized == sp_nm, , drop = FALSE] else X_loc_imp
          med_s <- apply(sub_m, 2L, stats::median, na.rm = TRUE)
          cent_db[si, ] <- med_s
          resid_list[[si]] <- sweep(sub_m, 2L, med_s, "-")
        }
        R_pool <- do.call(rbind, resid_list)
        S_raw_db <- if (nrow(R_pool) >= 3L) stats::cov(R_pool) else diag(1.0, ncol(X_loc_imp))
        S_raw_db[!is.finite(S_raw_db)] <- 0.0
        S_pool_db <- compute_ledoit_wolf_cov(S_raw_db + diag(1e-5, ncol(X_loc_imp)), X_centered = R_pool, n_obs = nrow(R_pool))
        S_db_inv <- tryCatch(solve(S_pool_db), error = function(e) diag(1.0 / pmax(diag(S_pool_db), 1e-4), ncol(X_loc_imp)))
        d2_all <- matrix(NA_real_, nrow = n, ncol = length(all_sp), dimnames = list(NULL, all_sp))
        for (si in seq_along(all_sp)) {
          df_s <- sweep(X_loc_imp, 2L, cent_db[si, ], "-")
          d2_all[, si] <- pmax(0.0, rowSums((df_s %*% S_db_inv) * df_s))
        }
        df_vec <- rep(ncol(X_loc_imp), n)
        # Re-evaluate any multi-island rows lacking harmonics strictly on their observed topological subset
        topo_pos <- which(!(use_c %in% c(paste0("A", 1:24), paste0("D", 1:24))))
        incomplete_rows <- which(!stats::complete.cases(X_loc_s))
        if (length(incomplete_rows) > 0L && length(topo_pos) >= 2L) {
          S_sub <- S_pool_db[topo_pos, topo_pos, drop = FALSE]
          S_sub_inv <- tryCatch(solve(S_sub), error = function(e) diag(1.0 / pmax(diag(S_sub), 1e-4), length(topo_pos)))
          for (ri in incomplete_rows) {
            obs_r <- topo_pos[is.finite(X_loc_s[ri, topo_pos])]
            if (length(obs_r) >= 2L) {
              S_r <- S_pool_db[obs_r, obs_r, drop = FALSE]
              S_r_inv <- tryCatch(solve(S_r), error = function(e) diag(1.0 / pmax(diag(S_r), 1e-4), length(obs_r)))
              for (si in seq_along(all_sp)) {
                dv_r <- X_loc_s[ri, obs_r] - cent_db[si, obs_r]
                d2_all[ri, si] <- max(0.0, as.numeric(t(dv_r) %*% S_r_inv %*% dv_r))
              }
              df_vec[ri] <- length(obs_r)
            }
          }
        }
        oth_cols <- match(oth_sp, all_sp)
        d2_oth <- d2_all[, oth_cols, drop = FALSE]
        best_c <- max.col(-d2_oth, ties.method = "first")
        list(
          min_congener_mahal_d2 = round(d2_oth[cbind(seq_len(n), best_c)], 4),
          min_congener_taxon = oth_sp[best_c],
          open_set_novelty_p = {
            min_d2_r <- vapply(seq_len(n), function(ri) { vr <- d2_all[ri, ]; vr <- vr[is.finite(vr)]; if (length(vr) > 0L) min(vr) else NA_real_ }, numeric(1))
            round(1.0 - stats::pchisq(min_d2_r, df = df_vec), 6)
          },
          open_set_df = as.integer(df_vec)
        )
      }

      res_inv <- if (length(inv_cands) >= 3L) eval_db_space(inv_cands, require_same_run = FALSE) else NULL
      if (!is.null(res_inv)) {
        res_inv
      } else if (length(raw_pc_cands) >= 2L) {
        eval_db_space(raw_pc_cands, require_same_run = TRUE)
      } else {
        NULL
      }
    } else NULL
  }, error = function(e) NULL)

  if (!is.null(db_osr)) return(db_osr)

  core_rows <- which(!is.na(gmm_cl) & gmm_cl == 1L)
  if (length(core_rows) < 2L) core_rows <- seq_len(n)
  mu_core <- colMeans(X_mat[core_rows, , drop = FALSE])
  diff_c <- sweep(X_mat, 2L, mu_core, "-")
  d2_core <- pmax(0.0, rowSums((diff_c %*% S_inv) * diff_c))
  open_set_novelty_p <- round(1.0 - stats::pchisq(d2_core, df = p), 6)
  list(min_congener_mahal_d2 = min_congener_mahal_d2, min_congener_taxon = min_congener_taxon, open_set_novelty_p = open_set_novelty_p, open_set_df = open_set_df)
}


.canonical_relabel_gmm_clusters <- function(sub_v, cls_sub, best_k, ref_cand, core_idx = NULL) {
  n_eff_ratio <- 1.0
  w_all <- if ("authority_weight" %in% names(sub_v) && any(is.finite(suppressWarnings(as.numeric(sub_v$authority_weight))))) {
    suppressWarnings(as.numeric(sub_v$authority_weight))
  } else {
    rep(1.0, nrow(sub_v))
  }
  w_all[!is.finite(w_all) | !ref_cand] <- 0.0
  if (all(w_all <= 0.0)) w_all[ref_cand] <- 1.0
  if (!is.null(core_idx) && length(core_idx) >= 2L) {
    w_c <- w_all[core_idx]
    w_c <- w_c[is.finite(w_c) & w_c > 0.0]
    if (length(w_c) >= 2L && sum(w_c^2) > 1e-12) {
      n_eff_ratio <- min(1.0, max(0.1, ((sum(w_c)^2) / sum(w_c^2)) / length(w_c)))
    }
  }
  if (!is.null(cls_sub) && best_k > 1L) {
    u_cl <- sort(unique(stats::na.omit(as.integer(cls_sub))))
    if (length(u_cl) > 1L) {
      cl_w <- tapply(w_all, factor(cls_sub, levels = u_cl), sum, na.rm = TRUE)
      cl_w[is.na(cl_w)] <- 0.0
      t_info <- .identify_type_specimen_vec(sub_v)
      t_idx <- which(t_info$is_type_specimen & !is.na(cls_sub))
      if (length(t_idx) > 0L) {
        t_wt <- pmax(1.0, w_all[t_idx]) * ifelse(t_info$type_rank[t_idx] == 2L, 2.0, 1.0)
        t_cl_w <- tapply(t_wt, factor(cls_sub[t_idx], levels = u_cl), sum, na.rm = TRUE)
        t_cl_w[is.na(t_cl_w)] <- 0.0
        ord_cl <- order(-t_cl_w, -cl_w, u_cl)
      } else {
        ord_cl <- order(-cl_w, u_cl)
      }
      cls_sub <- as.integer(match(as.integer(cls_sub), u_cl[ord_cl]))
    }
  }
  list(cls_sub = cls_sub, n_eff_ratio = n_eff_ratio)
}

# ==============================================================================
# LEVEL 1: WITHIN-PLANT BASAL ROSETTE & CAULINE FALLBACK OTU MEDIAN RECONCILIATION
# ==============================================================================

# Reconciles Step 03 specimen-level OTU basal medians (unweighted across clean basal rosette
# leaves; cauline leaves excluded whenever >= 1 basal leaf exists) and 75th-percentile petiole/blade
# architecture traits. Short-circuits when aggregate_step03_otu_medians_fast() already populated all fields.
reconcile_step03_basal_otu_medians <- function(specimen_df, closed_leaf_df, all_harm_cols = character(0), sym_names = character(0)) {
  if (is.null(specimen_df) || is.null(closed_leaf_df) || nrow(specimen_df) == 0L || nrow(closed_leaf_df) == 0L) {
    return(specimen_df)
  }
  expected_fast_cols <- intersect(c("area_mm2", "petiole_width_mm", "blade_width_mm", "blade_cie_l_norm", "blade_cie_l_sd", "ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity"), names(closed_leaf_df))
  if (isTRUE(attr(specimen_df, "step03_fast_aggregated")) && all(expected_fast_cols %in% names(specimen_df))) {
    return(specimen_df)
  }
  pids_leaf <- if ("plant_individual_id" %in% names(closed_leaf_df)) ifelse(is.na(closed_leaf_df$plant_individual_id), 0L, as.integer(closed_leaf_df$plant_individual_id)) else rep(0L, nrow(closed_leaf_df))
  pids_spec <- if ("plant_individual_id" %in% names(specimen_df)) ifelse(is.na(specimen_df$plant_individual_id), 0L, as.integer(specimen_df$plant_individual_id)) else rep(0L, nrow(specimen_df))
  leaf_keys <- paste(as.character(closed_leaf_df$catalogNumber), pids_leaf, sep = "__")
  spec_keys <- paste(as.character(specimen_df$catalogNumber), pids_spec, sep = "__")
  idx_by_key <- split(seq_len(nrow(closed_leaf_df)), leaf_keys)
  harm_cols <- intersect(unique(c(all_harm_cols, sym_names)), names(specimen_df))
  trait_digits <- list(
    aspect_ratio = 4L, circularity = 4L, solidity = 4L, dissection_index = 4L,
    area_px = 1L, area_mm2 = 2L, centroid_size = 4L,
    log_centsize = 6L, harmonic_roughness = 6L, leaf_quality_weight = 6L,
    pct_asymmetry = 6L, ECT_PC1 = 6L, ECT_PC2 = 6L, ECT_PC3 = 6L,
    ect_lobe_complexity = 4L,
    blade_cie_l_norm = 4L,
    blade_cie_l_sd = 4L
  )
  pet_cols <- intersect(c("petiole_length_mm", "petiole_width_mm", "blade_length_mm", "blade_width_mm", "petiole_to_blade_ratio"), names(specimen_df))
  for (i in seq_len(nrow(specimen_df))) {
    rows_i <- idx_by_key[[spec_keys[i]]]
    if (length(rows_i) < 1L) next
    grp <- closed_leaf_df[rows_i, , drop = FALSE]
    typ_g <- if ("is_type_specimen" %in% names(grp)) (!is.na(grp$is_type_specimen) & as.logical(grp$is_type_specimen)) else rep(FALSE, nrow(grp))
    if ("is_contour_outlier" %in% names(grp) && any((!as.logical(grp$is_contour_outlier)) | typ_g, na.rm = TRUE)) {
      grp <- grp[(!is.na(grp$is_contour_outlier) & !as.logical(grp$is_contour_outlier)) | typ_g, , drop = FALSE]
    }
    if ("leaf_zone" %in% names(grp)) {
      lz_g <- tolower(trimws(as.character(grp$leaf_zone)))
      if (any(!is.na(lz_g) & lz_g == "basal")) {
        grp <- grp[!is.na(lz_g) & lz_g == "basal", , drop = FALSE]
      }
    }
    if (nrow(grp) < 1L) next
    for (hc in harm_cols) {
      if (hc %in% names(grp)) {
        hv <- suppressWarnings(as.numeric(grp[[hc]]))
        hv <- hv[is.finite(hv)]
        if (length(hv) > 0L) specimen_df[[hc]][i] <- stats::median(hv)
      }
    }
    for (tc in names(trait_digits)) {
      if ((tc %in% names(specimen_df)) && (tc %in% names(grp))) {
        tv <- suppressWarnings(as.numeric(grp[[tc]]))
        tv <- tv[is.finite(tv)]
        if (length(tv) > 0L) specimen_df[[tc]][i] <- round(stats::median(tv), trait_digits[[tc]])
      }
    }
    for (pc in pet_cols) {
      if (pc %in% names(grp)) {
        pv <- suppressWarnings(as.numeric(grp[[pc]]))
        pv <- pv[is.finite(pv)]
        if (length(pv) > 0L) specimen_df[[pc]][i] <- round(as.numeric(stats::quantile(pv, probs = 0.75, na.rm = TRUE, names = FALSE)), 4L)
      }
    }
  }
  specimen_df
}

# ==============================================================================
# STEP 04 MODULAR CLI, EXACT MCLUST PARAMETER, SHARED-TOPOLOGY & AUDIT HELPERS
# ==============================================================================
.parse_step04_cli_fallback <- function(raw_args, defaults) {
  opts <- defaults
  if (length(raw_args) == 0L) return(opts)
  expanded <- character(0)
  for (tok in as.character(raw_args)) {
    if (grepl("^--[A-Za-z0-9_-]+=", tok)) {
      eq_pos <- regexpr("=", tok, fixed = TRUE)[1L]
      expanded <- c(expanded, substr(tok, 1L, eq_pos - 1L), substr(tok, eq_pos + 1L, nchar(tok)))
    } else {
      expanded <- c(expanded, tok)
    }
  }
  i <- 1L; n <- length(expanded)
  while (i <= n) {
    a <- expanded[i]
    nxt <- if (i < n) expanded[i + 1L] else ""
    if (a %in% c("-v", "--vouchers") && i < n) { opts$vouchers <- nxt; i <- i + 2L }
    else if (a %in% c("-e", "--harmonics") && i < n) { opts$harmonics <- nxt; i <- i + 2L }
    else if (a %in% c("-m", "--reproductive") && i < n) { opts$reproductive <- nxt; i <- i + 2L }
    else if (a == "--outline-only") { opts$outline_only <- TRUE; i <- i + 1L }
    else if (a == "--correct-allometry" && i < n) { opts$correct_allometry <- nxt; i <- i + 2L }
    else if (a == "--allometry-report" && i < n) { opts$allometry_report <- nxt; i <- i + 2L }
    else if (a == "--allometry-pairwise" && i < n) { opts$allometry_pairwise <- nxt; i <- i + 2L }
    else if (a == "--allometry-plot" && i < n) { opts$allometry_plot <- nxt; i <- i + 2L }
    else if (a == "--use-clustvarsel") { opts$use_clustvarsel <- TRUE; i <- i + 1L }
    else if (a %in% c("-f", "--output-flags") && i < n) { opts$output_flags <- nxt; i <- i + 2L }
    else if (a %in% c("-p", "--output-plot") && i < n) { opts$output_plot <- nxt; i <- i + 2L }
    else if (a %in% c("-r", "--output-report") && i < n) { opts$output_report <- nxt; i <- i + 2L }
    else if (a %in% c("-b", "--output-batch-audit") && i < n) { opts$output_batch_audit <- nxt; i <- i + 2L }
    else if (a %in% c("-k", "--max-k") && i < n) { opts$max_k <- as.integer(nxt); i <- i + 2L }
    else if (a %in% c("-d", "--num-pcs") && i < n) { opts$num_pcs <- as.integer(nxt); i <- i + 2L }
    else if (a == "--permutations" && i < n) { opts$permutations <- as.integer(nxt); i <- i + 2L }
    else if (a == "--seed" && i < n) { opts$seed <- as.integer(nxt); i <- i + 2L }
    else if (a == "--intraspecific-mode") { opts$intraspecific_mode <- TRUE; i <- i + 1L }
    else if (a == "--focal-taxon" && i < n) { opts$focal_taxon <- nxt; i <- i + 2L }
    else { i <- i + 1L }
  }
  opts
}

.compute_exact_mclust_params <- function(K_vec, p_dim, bic_mat = NULL, default_model = "VVI", explicit_n_params = NULL) {
  K_int <- pmax(1L, as.integer(K_vec))
  p <- max(1L, as.integer(p_dim))
  if (!is.null(explicit_n_params) && length(explicit_n_params) == length(K_int) && all(is.finite(explicit_n_params))) {
    return(as.numeric(explicit_n_params))
  }
  row_models <- rep(sub("_RidgelineMerged$", "", as.character(default_model)[1L]), length(K_int))
  if (!is.null(bic_mat) && is.matrix(bic_mat) && !is.null(colnames(bic_mat))) {
    for (idx in seq_along(K_int)) {
      r_nm <- as.character(K_int[idx])
      if (r_nm %in% rownames(bic_mat)) {
        rv <- suppressWarnings(as.numeric(bic_mat[r_nm, ]))
        if (any(is.finite(rv))) row_models[idx] <- colnames(bic_mat)[which.max(ifelse(is.finite(rv), rv, -Inf))]
      }
    }
  }
  has_mclust <- requireNamespace("mclust", quietly = TRUE)
  vapply(seq_along(K_int), function(idx) {
    G <- K_int[idx]
    mn <- toupper(trimws(row_models[idx]))
    if (has_mclust) {
      np_m <- tryCatch(as.numeric(mclust::nMclustParams(modelName = mn, d = p, G = G)), error = function(e) NA_real_)
      if (is.finite(np_m) && np_m >= 1.0) return(np_m)
    }
    base_p <- (G - 1.0) + G * p
    cov_p <- switch(
      mn,
      "E" = 1.0,
      "V" = as.numeric(G),
      "EII" = 1.0,
      "VII" = as.numeric(G),
      "EEI" = as.numeric(p),
      "VEI" = as.numeric(G + p - 1L),
      "EVI" = as.numeric(1L + G * (p - 1L)),
      "VVI" = as.numeric(G * p),
      "EEE" = as.numeric(p * (p + 1L) / 2.0),
      "VVV" = as.numeric(G * p * (p + 1L) / 2.0),
      "FULL_COVARIANCE_EM" = as.numeric(G * p * (p + 1L) / 2.0),
      as.numeric(G * p)
    )
    base_p + cov_p
  }, numeric(1))
}

.project_shared_topology_rows_step04 <- function(df, cda_res, gmm_res, pca_cols) {
  if (is.null(df) || nrow(df) == 0L || length(pca_cols) == 0L) return(list(cda_res = cda_res, gmm_res = gmm_res))
  guard_cols <- pca_cols[seq_len(min(2L, length(pca_cols)))]
  if (!all(guard_cols %in% names(df))) return(list(cda_res = cda_res, gmm_res = gmm_res))
  no_efa_idx <- which(!stats::complete.cases(df[, guard_cols, drop = FALSE]))
  if (length(no_efa_idx) == 0L) return(list(cda_res = cda_res, gmm_res = gmm_res))

  cda_res$is_active[no_efa_idx] <- FALSE
  # Identify shared whole+dissected topological & architectural features (never EFA PCs)
  topo_pairs <- list(
    c("ECT_PC1_std", "ECT_PC1"),
    c("ECT_PC2_std", "ECT_PC2"),
    c("ECT_PC3_std", "ECT_PC3"),
    c("ect_lobe_complexity_std", "ect_lobe_complexity"),
    c("solidity_std", "solidity"),
    c("q25_solidity_std", "q25_solidity"),
    c("dissection_index_std", "dissection_index"),
    c("aspect_ratio_std", "aspect_ratio"),
    c("AR_inv_std", "AR_inv"),
    c("petiole_to_blade_ratio_std", "petiole_to_blade_ratio"),
    c("scape_to_basal_leaf_ratio_std", "scape_to_basal_leaf_ratio")
  )
  act_idx <- which(as.logical(cda_res$is_active))
  shared_cols <- character(0)
  for (pr in topo_pairs) {
    hit <- intersect(pr, names(df))
    if (length(hit) > 0L) {
      cn <- hit[1L]
      v_num <- suppressWarnings(as.numeric(df[[cn]]))
      if (length(act_idx) >= 3L && sum(is.finite(v_num[act_idx])) >= max(3L, ceiling(0.50 * length(act_idx))) &&
          any(is.finite(v_num[no_efa_idx])) && isTRUE(stats::sd(v_num[act_idx], na.rm = TRUE) > 1e-6)) {
        shared_cols <- c(shared_cols, cn)
      }
    }
  }
  shared_cols <- unique(shared_cols)

  if (length(shared_cols) >= 2L && length(act_idx) >= 4L) {
    X_sh_raw <- as.matrix(df[, shared_cols, drop = FALSE]); storage.mode(X_sh_raw) <- "double"
    mu_sh <- vapply(seq_len(ncol(X_sh_raw)), function(j) mean(X_sh_raw[act_idx, j], na.rm = TRUE), numeric(1))
    sd_sh <- vapply(seq_len(ncol(X_sh_raw)), function(j) { s <- stats::sd(X_sh_raw[act_idx, j], na.rm = TRUE); if (is.finite(s) && s > 1e-6) s else 1.0 }, numeric(1))
    X_sh <- sweep(sweep(X_sh_raw, 2L, mu_sh, "-"), 2L, sd_sh, "/")
    for (j in seq_len(ncol(X_sh))) {
      bad_j <- !is.finite(X_sh[, j])
      if (any(bad_j)) X_sh[bad_j, j] <- 0.0
    }
    obs_counts <- rowSums(is.finite(X_sh_raw[no_efa_idx, , drop = FALSE]))
    y_act <- as.character(cda_res$grouping_col[act_idx])
    g_lvls <- if (!is.null(cda_res$post_probs) && !is.null(colnames(cda_res$post_probs))) colnames(cda_res$post_probs) else sort(unique(stats::na.omit(y_act)))
    w_sb_sh <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight[act_idx])) else rep(1.0, length(act_idx)); w_sb_sh[!is.finite(w_sb_sh) | w_sb_sh <= 0] <- 1.0
    w_act <- (if (!is.null(cda_res$authority_weight)) as.numeric(cda_res$authority_weight[act_idx]) else rep(1.0, length(act_idx))) * w_sb_sh
    w_act[!is.finite(w_act) | w_act <= 0.0] <- 1.0
    G_sh <- length(g_lvls)
    if (G_sh >= 2L) {
      mu_g <- matrix(0.0, nrow = G_sh, ncol = ncol(X_sh), dimnames = list(g_lvls, shared_cols))
      X_cent_sh <- matrix(0.0, nrow = length(act_idx), ncol = ncol(X_sh))
      W_sh <- matrix(0.0, nrow = ncol(X_sh), ncol = ncol(X_sh))
      B_sh <- matrix(0.0, nrow = ncol(X_sh), ncol = ncol(X_sh))
      gm_sh <- colSums(X_sh[act_idx, , drop = FALSE] * w_act) / sum(w_act)
      for (ki in seq_along(g_lvls)) {
        km <- which(y_act == g_lvls[ki])
        if (length(km) > 0L) {
          sub_k <- X_sh[act_idx[km], , drop = FALSE]; wk <- w_act[km]
          mk <- colSums(sub_k * wk) / sum(wk); mu_g[ki, ] <- mk
          dm <- mk - gm_sh; B_sh <- B_sh + sum(wk) * tcrossprod(dm)
          dx <- sweep(sub_k, 2L, mk, "-"); X_cent_sh[km, ] <- dx
          W_sh <- W_sh + crossprod(dx, dx * wk)
        }
      }
      S_w_sh <- W_sh / max(sum(w_act) - G_sh, 1.0)
      S_reg_sh <- compute_ledoit_wolf_cov(S_w_sh, X_centered = X_cent_sh, n_obs = length(act_idx), weights = w_act)
      S_inv_sh <- tryCatch(solve(S_reg_sh), error = function(e) diag(1.0 / pmax(diag(S_reg_sh), 1e-4), ncol(X_sh)))
      eig_sh <- tryCatch(solve_shrinkage_cda_eigen(S_reg_sh, B_sh, max_axes = min(ncol(X_sh), G_sh - 1L)), error = function(e) NULL)
      # Align canonical axis signs with primary CDA centroids so joint ordination plots remain coherent
      if (!is.null(eig_sh) && !is.null(cda_res$centroids_full) && !is.null(rownames(cda_res$centroids_full)) && all(g_lvls %in% rownames(cda_res$centroids_full))) {
        c_sh <- sweep(mu_g, 2L, gm_sh, "-") %*% eig_sh$vectors
        k_match <- min(ncol(c_sh), ncol(cda_res$centroids_full))
        for (aj in seq_len(k_match)) {
          r_al <- suppressWarnings(stats::cor(c_sh[, aj], cda_res$centroids_full[g_lvls, aj], use = "complete.obs"))
          if (is.finite(r_al) && r_al < 0.0) eig_sh$vectors[, aj] <- -eig_sh$vectors[, aj]
        }
      }
      for (ii_pos in seq_along(no_efa_idx)) {
        ri <- no_efa_idx[ii_pos]
        if (obs_counts[ii_pos] < 2L) {
          cda_res$predicted_taxon[ri] <- NA_character_; cda_res$posterior_prob[ri] <- NA_real_
          cda_res$second_taxon[ri] <- NA_character_; cda_res$second_posterior_prob[ri] <- NA_real_
          cda_res$posterior_margin[ri] <- NA_real_; cda_res$canonical_scores[ri, ] <- NA_real_
          if (!is.null(cda_res$oof_posterior_prob)) cda_res$oof_posterior_prob[ri] <- NA_real_
          next
        }
        lp <- numeric(G_sh); d2_sh_vec <- numeric(G_sh)
        for (ki in seq_along(g_lvls)) {
          dv <- X_sh[ri, ] - mu_g[ki, ]
          d2_sh_vec[ki] <- max(0.0, as.numeric(t(dv) %*% S_inv_sh %*% dv))
          lp[ki] <- log(1.0 / G_sh) - 0.5 * d2_sh_vec[ki]
        }
        pp <- exp(lp - max(lp)); pp <- pp / max(sum(pp), 1e-12)
        if (!is.null(cda_res$post_probs) && ncol(cda_res$post_probs) == G_sh) cda_res$post_probs[ri, ] <- pp
        best_k <- which.max(pp); oth_k <- setdiff(seq_len(G_sh), best_k); sec_k <- oth_k[which.max(pp[oth_k])]
        cda_res$predicted_taxon[ri] <- g_lvls[best_k]; cda_res$posterior_prob[ri] <- pp[best_k]
        cda_res$second_taxon[ri] <- g_lvls[sec_k]; cda_res$second_posterior_prob[ri] <- pp[sec_k]
        cda_res$posterior_margin[ri] <- pp[best_k] - pp[sec_k]
        if (!is.null(cda_res$oof_posterior_prob)) cda_res$oof_posterior_prob[ri] <- pp[best_k]
        if (!is.null(cda_res$cl_class_map) && (g_lvls[best_k] %in% names(cda_res$cl_class_map))) {
          cda_res$cl_target_threshold[ri] <- unname(cda_res$cl_class_map[g_lvls[best_k]])
        }
        yi_sh <- match(as.character(cda_res$grouping_col[ri]), g_lvls)
        if (!is.na(yi_sh) && !is.null(cda_res$self_morph_mahal_d2)) {
          cda_res$self_morph_mahal_d2[ri] <- round(d2_sh_vec[yi_sh], 4)
          cda_res$self_morph_df[ri] <- as.integer(ncol(X_sh))
        }
        if (!is.null(eig_sh)) {
          z_r <- as.numeric((X_sh[ri, ] - gm_sh) %*% eig_sh$vectors)
          n_fill <- min(length(z_r), ncol(cda_res$canonical_scores))
          cda_res$canonical_scores[ri, ] <- 0.0
          cda_res$canonical_scores[ri, seq_len(n_fill)] <- z_r[seq_len(n_fill)]
        }
      }
    }
    if (!is.null(gmm_res) && !is.null(gmm_res$classification)) {
      cl_act <- gmm_res$classification[act_idx]
      u_cl <- sort(unique(stats::na.omit(cl_act)))
      u_cl <- u_cl[vapply(u_cl, function(ck) any(!is.na(cl_act) & cl_act == ck), logical(1))]
      if (length(u_cl) >= 1L) {
        cl_meds <- do.call(rbind, lapply(u_cl, function(ck) apply(X_sh[act_idx[!is.na(cl_act) & cl_act == ck], , drop = FALSE], 2L, stats::median, na.rm = TRUE)))
        for (ii_pos in seq_along(no_efa_idx)) {
          ri <- no_efa_idx[ii_pos]
          if (obs_counts[ii_pos] >= 2L) {
            d2_c <- rowSums(sweep(cl_meds, 2L, X_sh[ri, ], "-")^2)
            wc <- exp(-0.5 * (d2_c - min(d2_c))); wc <- wc / max(sum(wc), 1e-12)
            gmm_res$classification[ri] <- u_cl[which.max(wc)]
            if (!is.null(gmm_res$uncertainty)) gmm_res$uncertainty[ri] <- 1.0 - max(wc)
          }
        }
      }
    }
    return(list(cda_res = cda_res, gmm_res = gmm_res))
  }

  # Fewer than 2 shared features available: mark unprojected rather than mean-imputing PC1 = 0
  cda_res$predicted_taxon[no_efa_idx] <- NA_character_
  cda_res$posterior_prob[no_efa_idx] <- NA_real_
  cda_res$second_taxon[no_efa_idx] <- NA_character_
  cda_res$second_posterior_prob[no_efa_idx] <- NA_real_
  cda_res$posterior_margin[no_efa_idx] <- NA_real_
  if (!is.null(cda_res$canonical_scores)) cda_res$canonical_scores[no_efa_idx, ] <- NA_real_
  list(cda_res = cda_res, gmm_res = gmm_res)
}

.compute_oof_cda_and_self_ellipsoid <- function(
  X_proj,
  active_idx,
  y_a,
  w_a,
  g_levels,
  group_means,
  S_reg,
  grouping_col_all,
  predicted_taxon_all,
  posterior_prob_all,
  coords_df = NULL,
  intraspecific_mode = FALSE,
  deme_radius_km = 10.0
) {
  n_all <- nrow(X_proj)
  p <- max(1L, ncol(X_proj))
  G <- length(g_levels)
  self_d2 <- rep(NA_real_, n_all)
  self_df <- rep(as.integer(p), n_all)
  oof_prob <- as.numeric(posterior_prob_all)
  cl_thr_vec <- rep(NA_real_, n_all)
  cl_map <- stats::setNames(rep(NA_real_, G), g_levels)

  S_inv <- tryCatch(
    solve(S_reg),
    error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-5), nrow = p, ncol = p)
  )

  # 1. Leave-one-out within-taxon Mahalanobis tolerance distance (Zapata & Jimenez 2012)
  #    Depends strictly on active_idx so passive Tier_3_Bronze rows have zero influence on active anchors.
  sum_w_g <- stats::setNames(numeric(G), g_levels)
  sum_wx_g <- matrix(0.0, nrow = G, ncol = p, dimnames = list(g_levels, colnames(X_proj)))
  neff_g <- stats::setNames(rep(1.0, G), g_levels)
  n_act_g <- stats::setNames(integer(G), g_levels)
  act_pos_map <- rep(NA_integer_, n_all)
  act_pos_map[active_idx] <- seq_along(active_idx)

  for (k in seq_len(G)) {
    gk <- g_levels[k]
    mk <- which(y_a == gk)
    n_act_g[k] <- length(mk)
    if (length(mk) > 0L) {
      wk <- pmax(1e-6, as.numeric(w_a[mk]))
      sum_w_g[k] <- sum(wk)
      sum_wx_g[k, ] <- colSums(X_proj[active_idx[mk], , drop = FALSE] * wk)
      neff_g[k] <- max(1.0, (sum(wk)^2) / max(sum(wk^2), 1e-12))
    } else {
      sum_wx_g[k, ] <- group_means[k, ]
      sum_w_g[k] <- 1.0
    }
  }

  grp_chr <- as.character(grouping_col_all)
  for (i in seq_len(n_all)) {
    gi <- grp_chr[i]
    if (is.na(gi) || !(gi %in% g_levels)) next
    ki <- match(gi, g_levels)
    ap <- act_pos_map[i]
    if (!is.na(ap) && n_act_g[ki] >= 3L) {
      wi <- max(1e-6, as.numeric(w_a[ap]))
      rem_w <- sum_w_g[ki] - wi
      mu_i <- if (rem_w > 1e-8) (sum_wx_g[ki, ] - wi * X_proj[i, ]) / rem_w else group_means[ki, ]
      n_eff_use <- max(1.0, neff_g[ki] - 1.0)
    } else {
      mu_i <- group_means[ki, ]
      n_eff_use <- neff_g[ki]
    }
    dv <- X_proj[i, ] - mu_i
    raw_d2 <- max(0.0, as.numeric(t(dv) %*% S_inv %*% dv))
    self_d2[i] <- round(raw_d2 / (1.0 + 1.0 / max(1.0, n_eff_use)), 4)
  }

  # 2. Deme-blocked Out-of-Fold CDA Posteriors & Confident Learning Class Thresholds (Northcutt et al. 2021)
  if (!isTRUE(intraspecific_mode) && G >= 2L && length(active_idx) >= (G + 2L) && !all(grepl("^Morphotype_K", g_levels))) {
    N_a <- length(active_idx)
    X_a <- X_proj[active_idx, , drop = FALSE]
    y_fac <- factor(y_a, levels = g_levels)
    sub_coords <- if (!is.null(coords_df) && is.data.frame(coords_df) && nrow(coords_df) >= max(active_idx)) {
      coords_df[active_idx, , drop = FALSE]
    } else {
      NULL
    }
    fold_vec <- build_taxon_blind_gathering_deme_folds(sub_coords, n_obs = N_a, deme_radius_km = deme_radius_km)
    for (k in seq_len(G)) {
      ik <- which(as.integer(y_fac) == k)
      if (length(unique(fold_vec[ik])) < 2L && length(ik) >= 2L) {
        fold_vec[ik] <- max(fold_vec, na.rm = TRUE) + seq_along(ik)
      }
    }
    u_folds <- sort(unique(fold_vec))
    P_oof_act <- matrix(1.0 / G, nrow = N_a, ncol = G, dimnames = list(NULL, g_levels))
    priors <- rep(1.0 / G, G)

    for (fid in u_folds) {
      te <- which(fold_vec == fid)
      tr <- which(fold_vec != fid)
      if (length(tr) < G) tr <- setdiff(seq_len(N_a), te[1L])
      X_tr <- X_a[tr, , drop = FALSE]
      y_tr <- y_fac[tr]
      w_tr <- pmax(1e-6, as.numeric(w_a[tr]))
      mu_f <- matrix(0.0, nrow = G, ncol = p)
      W_f <- matrix(0.0, nrow = p, ncol = p)
      Xc_f <- matrix(0.0, nrow = nrow(X_tr), ncol = p)
      w_corr <- 0.0
      for (k in seq_len(G)) {
        ik <- which(as.integer(y_tr) == k)
        if (length(ik) == 0L) {
          mu_f[k, ] <- group_means[k, ]
          next
        }
        wk <- w_tr[ik]; sw <- sum(wk)
        w_corr <- w_corr + (sum(wk^2) / sw)
        mk <- colSums(X_tr[ik, , drop = FALSE] * wk) / sw
        mu_f[k, ] <- mk
        dk <- sweep(X_tr[ik, , drop = FALSE], 2L, mk, "-")
        Xc_f[ik, ] <- dk
        W_f <- W_f + crossprod(dk, dk * wk)
      }
      S_wf <- W_f / max(sum(w_tr) - w_corr, 1.0)
      S_rf <- compute_ledoit_wolf_cov(S_wf, X_centered = Xc_f, n_obs = nrow(X_tr), weights = w_tr)
      if (is.list(S_rf) && !is.null(S_rf$S_reg)) S_rf <- S_rf$S_reg
      S_if <- tryCatch(solve(S_rf), error = function(e) S_inv)
      for (ti in seq_along(te)) {
        obs_a <- te[ti]
        diff_i <- sweep(mu_f, 2L, X_a[obs_a, ], "-")
        d2_i <- pmax(0.0, rowSums((diff_i %*% S_if) * diff_i))
        ll_i <- log(priors) - 0.5 * d2_i
        pp_i <- exp(ll_i - max(ll_i))
        P_oof_act[obs_a, ] <- pp_i / max(sum(pp_i), 1e-12)
      }
    }

    for (k in seq_len(G)) {
      ik <- which(as.integer(y_fac) == k)
      if (length(ik) > 0L) {
        wk <- pmax(1e-6, as.numeric(w_a[ik]))
        tk_raw <- sum(wk * P_oof_act[ik, k]) / sum(wk)
        cl_map[k] <- round(max(0.52, min(0.86, tk_raw)), 4)
      } else {
        cl_map[k] <- 0.68
      }
    }

    pred_chr <- as.character(predicted_taxon_all)
    for (pos_a in seq_len(N_a)) {
      ri <- active_idx[pos_a]
      pr_tx <- pred_chr[ri]
      if (!is.na(pr_tx) && (pr_tx %in% g_levels)) {
        oof_prob[ri] <- round(P_oof_act[pos_a, pr_tx], 4)
      }
    }
    hit_pr <- match(pred_chr, g_levels)
    ok_pr <- which(!is.na(hit_pr))
    if (length(ok_pr) > 0L) {
      cl_thr_vec[ok_pr] <- cl_map[hit_pr[ok_pr]]
    }
  }

  list(
    self_morph_mahal_d2 = self_d2,
    self_morph_df = self_df,
    oof_posterior_prob = round(oof_prob, 4),
    cl_target_threshold = cl_thr_vec,
    cl_class_map = cl_map
  )
}

.format_step04_anchor_audit_reason <- function(i, df, cda_res, raw_sp, pred_sp, post_p, is_act) {
  if (is.na(pred_sp) || !nzchar(pred_sp) || !is.finite(post_p)) {
    return(list(flag = FALSE, triage = "LOW", reason = "Topology_Only_Unprojected_Multi_Island"))
  }
  ats <- if (!is.null(cda_res$anchor_tier_source) && length(cda_res$anchor_tier_source) >= i) as.character(cda_res$anchor_tier_source[i]) else ""
  det_t <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier[i]) else ""
  mts <- if ("morph_tier_source" %in% names(df)) as.character(df$morph_tier_source[i]) else "Basal_Rosette"
  is_typ <- (!is.na(ats) && ats == "Tier_0_Type_Anchor") || ("is_type_specimen" %in% names(df) && isTRUE(as.logical(df$is_type_specimen[i])))
  is_promo <- (!is.na(ats) && ats == "Tier_3_Bronze_Promoted") || (!is.null(cda_res$is_promoted_bronze) && isTRUE(as.logical(cda_res$is_promoted_bronze[i])))

  # Self-tolerance ellipsoid protection (Zapata & Jimenez 2012): do not flag a specimen for interspecific
  # morphological discordance when it sits inside its recorded taxon's 95% Mahalanobis tolerance region.
  is_morphotype_grp <- grepl("^Morphotype_K", raw_sp) || grepl("^Morphotype_K", pred_sp)
  self_d2_i <- if (!is.null(cda_res$self_morph_mahal_d2) && length(cda_res$self_morph_mahal_d2) >= i) cda_res$self_morph_mahal_d2[i] else NA_real_
  self_df_i <- if (!is.null(cda_res$self_morph_df) && length(cda_res$self_morph_df) >= i && !is.na(cda_res$self_morph_df[i])) as.integer(cda_res$self_morph_df[i]) else 3L
  in_self_95 <- !is_morphotype_grp && is.finite(self_d2_i) && (self_d2_i <= stats::qchisq(0.95, df = max(1L, self_df_i)))
  oof_p_i <- if (!is.null(cda_res$oof_posterior_prob) && length(cda_res$oof_posterior_prob) >= i && is.finite(cda_res$oof_posterior_prob[i])) cda_res$oof_posterior_prob[i] else post_p
  cl_thr_i <- if (!is.null(cda_res$cl_target_threshold) && length(cda_res$cl_target_threshold) >= i && is.finite(cda_res$cl_target_threshold[i])) cda_res$cl_target_threshold[i] else NA_real_

  if (is_act) {
    act_lbl <- if (is_typ) "Type_Anchor" else if (is_promo) "Promoted_Bronze_Anchor" else if (!is.na(det_t) && det_t == "Tier_2_Silver") "Tier2_Silver_Anchor" else "Tier1_Gold_Anchor"
    disc_pref <- if (is_typ) "Type_Anchor" else if (is_promo) "Promoted_Bronze" else if (!is.na(det_t) && det_t == "Tier_2_Silver") "Tier_2_Silver" else "Tier_1_Gold"
    if (raw_sp == pred_sp) return(list(flag = FALSE, triage = "CLEAN", reason = paste0("Verified_", act_lbl)))
    if (in_self_95) return(list(flag = FALSE, triage = "LOW", reason = paste0(disc_pref, "_Within_Self_Ellipsoid")))
    act_cut <- if (is.finite(cl_thr_i)) max(0.80, cl_thr_i) else 0.85
    if (oof_p_i >= act_cut && post_p >= 0.80) return(list(flag = TRUE, triage = "HIGH", reason = sprintf("%s_Morphological_Discordance_to_%s", disc_pref, gsub(" ", "_", pred_sp))))
    return(list(flag = FALSE, triage = "MEDIUM", reason = paste0(disc_pref, "_Borderline_Variant")))
  }
  cong_lbl <- if (!is.na(mts) && mts == "Cauline_Only_Fallback") {
    "Congruent_Cauline_Fallback_Determination"
  } else if (!is.na(det_t) && det_t %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver")) {
    "Congruent_Passive_Specialist_Determination"
  } else {
    "Congruent_Bronze_Determination"
  }
  if (raw_sp == pred_sp && post_p >= 0.75) return(list(flag = FALSE, triage = "CLEAN", reason = cong_lbl))
  if (raw_sp != pred_sp && in_self_95) return(list(flag = FALSE, triage = "LOW", reason = "Within_Self_Ellipsoid_Overlap_Protected"))
  pas_cut <- if (is.finite(cl_thr_i)) max(0.65, min(0.78, cl_thr_i)) else 0.70
  if (raw_sp != pred_sp && post_p >= pas_cut) return(list(flag = TRUE, triage = "HIGH", reason = sprintf("%s_Discordant_to_%s", gsub(" ", "_", raw_sp), gsub(" ", "_", pred_sp))))
  if (post_p < 0.60) return(list(flag = FALSE, triage = "MEDIUM", reason = "Morphological_Intermediate_Ambiguous"))
  if (raw_sp != pred_sp) return(list(flag = TRUE, triage = "MEDIUM", reason = "Candidate_Discordance_Moderate_Posterior"))
  list(flag = FALSE, triage = "LOW", reason = "Minor_Uncertainty")
}

.update_step04_allometry_pcs <- function(morpho_df, sub_df, gold_df, valid_idx, sym_names, active_cols, pca_resids, ref_mean_cs, eff_weight) {
  all_sym_idx <- which(stats::complete.cases(morpho_df[, sym_names, drop = FALSE]))
  if (length(all_sym_idx) == 0L) return(morpho_df)
  proj_pcs <- scale(as.matrix(morpho_df[all_sym_idx, sym_names[active_cols], drop = FALSE]), center = pca_resids$center, scale = FALSE) %*% pca_resids$rotation
  pc_allom_cols <- grep("^PC[0-9]+$", names(morpho_df), value = TRUE)
  if (length(pc_allom_cols) > 1L) pc_allom_cols <- pc_allom_cols[order(as.integer(gsub("[^0-9]", "", pc_allom_cols)))]
  asym_names <- c(grep("^B[0-9]+$", names(morpho_df), value = TRUE), grep("^C[0-9]+$", names(morpho_df), value = TRUE))
  has_asym <- length(asym_names) > 0L && any(vapply(asym_names, function(ac) isTRUE(stats::sd(morpho_df[[ac]][all_sym_idx], na.rm = TRUE) > 1e-6), logical(1)))
  if (!has_asym) {
    for (p in seq_along(pc_allom_cols)) morpho_df[all_sym_idx, pc_allom_cols[p]] <- if (p <= ncol(proj_pcs)) round(proj_pcs[, p], 6) else 0.0
  } else {
    full_h <- c(sym_names, asym_names)
    fit_as <- if (nlevels(gold_df$scientificName) >= 2L) stats::lm(as.matrix(gold_df[, asym_names, drop = FALSE]) ~ log_centsize + scientificName, data = gold_df, weights = eff_weight) else stats::lm(as.matrix(gold_df[, asym_names, drop = FALSE]) ~ log_centsize, data = gold_df, weights = eff_weight)
    b_as <- stats::coef(fit_as)["log_centsize", ]; b_as[is.na(b_as)] <- 0.0
    morpho_df[valid_idx, asym_names] <- round(as.matrix(sub_df[, asym_names, drop = FALSE]) - outer(sub_df$log_centsize - ref_mean_cs, b_as), 6)
    act_fh <- full_h[vapply(full_h, function(cl) isTRUE(stats::var(morpho_df[[cl]][valid_idx], na.rm = TRUE) > 1e-8), logical(1))]
    sub_pos <- match(rownames(gold_df), rownames(sub_df))
    gold_pos <- if (!any(is.na(sub_pos))) valid_idx[sub_pos] else valid_idx
    all_fh_idx <- which(stats::complete.cases(morpho_df[, act_fh, drop = FALSE]))
    if (length(act_fh) >= 2L && length(all_fh_idx) > 0L) {
      wt_fh <- stats::cov.wt(as.matrix(morpho_df[gold_pos, act_fh, drop = FALSE]), wt = eff_weight, center = TRUE)
      eig_fh <- eigen(wt_fh$cov, symmetric = TRUE)
      proj_fh <- sweep(as.matrix(morpho_df[all_fh_idx, act_fh, drop = FALSE]), 2L, wt_fh$center, "-") %*% eig_fh$vectors
      for (p in seq_along(pc_allom_cols)) morpho_df[all_fh_idx, pc_allom_cols[p]] <- if (p <= ncol(proj_fh)) round(proj_fh[, p], 6) else 0.0
    }
  }
  morpho_df
}

.update_step04_allometry_scalars <- function(morpho_df, sub_df, gold_df, valid_idx, ref_mean_cs) {
  sc_allom_cands <- intersect(c("ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "solidity", "q25_solidity", "dissection_index", "aspect_ratio", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio"), names(morpho_df))
  for (sc_nm in sc_allom_cands) {
    yv_g <- suppressWarnings(as.numeric(gold_df[[sc_nm]])); ok_g <- is.finite(yv_g) & is.finite(gold_df$log_centsize)
    if (sum(ok_g) >= 10L && stats::sd(yv_g[ok_g]) > 1e-6) {
      g_sub <- gold_df[ok_g, , drop = FALSE]; g_sub$.y_sc <- yv_g[ok_g]; g_sub$scientificName <- droplevels(factor(g_sub$scientificName))
      p_int_sc <- 1.0
      if (nlevels(g_sub$scientificName) >= 2L) {
        fit_i_sc <- tryCatch(stats::lm(.y_sc ~ log_centsize * scientificName, data = g_sub, weights = eff_weight), error = function(e) NULL)
        if (!is.null(fit_i_sc)) { aov_i <- stats::anova(fit_i_sc); idx_i <- grep(":", rownames(aov_i))[1L]; if (!is.na(idx_i) && is.finite(aov_i[["Pr(>F)"]][idx_i])) p_int_sc <- aov_i[["Pr(>F)"]][idx_i] }
      }
      if (p_int_sc > 0.05) {
        fit_s_sc <- if (nlevels(g_sub$scientificName) >= 2L) stats::lm(.y_sc ~ log_centsize + scientificName, data = g_sub, weights = eff_weight) else stats::lm(.y_sc ~ log_centsize, data = g_sub, weights = eff_weight)
        b_sc <- unname(stats::coef(fit_s_sc)["log_centsize"]); p_s_sc <- tryCatch(summary(fit_s_sc)$coefficients["log_centsize", "Pr(>|t|)"], error = function(e) 1.0)
        if (is.finite(b_sc) && is.finite(p_s_sc) && p_s_sc < 0.05) {
          yv_all <- suppressWarnings(as.numeric(morpho_df[[sc_nm]][valid_idx])); fin_v <- is.finite(yv_all)
          if (any(fin_v)) morpho_df[[sc_nm]][valid_idx[fin_v]] <- round(yv_all[fin_v] - b_sc * (sub_df$log_centsize[fin_v] - ref_mean_cs), 6)
        }
      }
    }
  }
  attr(morpho_df, "allometry_corrected_scalars") <- TRUE
  morpho_df
}
