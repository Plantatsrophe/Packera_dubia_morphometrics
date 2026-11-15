#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/06_genomic_phylogenomics.R
# Project: Plantmatics
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Standalone Step 06 Nuclear & Organellar Phylogenomics Stream.
#   1. Executes as a clean zero-overhead no-op when `data/genomics/` holds no payloads.
#   2. Computes pairwise allele-sharing/dosage distances, M/P replicate concordance,
#      PCoA, PAM clustering, LORO-CV shrinkage DAPC, convex admixture simplex,
#      122 bp microhaplotype co-ancestry, (h, Q12) hybrid classifier, tree-guided
#      f-branch, SNPs2CF/SNaQ quartets, and invariant-site SVDquartets/BPP/gdi.
#   3. Persists Tier B SQLite `genetic_lineage_assignments` (`20261123`) BEFORE
#      running post-hoc MMRR/FEEMS/Mantel cross-modal congruence.
#   4. Renders 300-DPI PNG plates `fig06a_genomic_pcoa_tree_concordance.png` and
#      `fig06b_reticulation_gdi_landscape.png`.
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("ape", quietly = TRUE)) library(ape)
  if (requireNamespace("cluster", quietly = TRUE)) library(cluster)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
  if (requireNamespace("patchwork", quietly = TRUE)) library(patchwork)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
}

.src_step06_rel <- function(rel_path) {
  cands <- c(
    rel_path,
    file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), basename(rel_path)),
    file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "..", "..", rel_path)
  )
  for (p in cands) {
    if (!is.na(p) && nzchar(p) && file.exists(p)) {
      source(p)
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

if (!exists("persist_step06_genomics_to_sqlite", mode = "function")) .src_step06_rel("scripts/core/db_utils.R")
if (!exists("save_packera_png_300", mode = "function")) .src_step06_rel("scripts/figures/theme_packera_pub.R")
if (!exists("compute_microhaplotype_coancestry", mode = "function")) .src_step06_rel("scripts/analysis/step06_reticulation_coancestry.R")
if (!exists("run_step06_landscape_genomics", mode = "function")) .src_step06_rel("scripts/analysis/step06_spatial_landscape_genomics.R")
if (!exists("load_public_packera_backbone", mode = "function")) .src_step06_rel("scripts/analysis/step06_public_backbone.R")
if (!exists("render_fig06a_genomic_pcoa_tree_png", mode = "function")) .src_step06_rel("scripts/figures/step06_figures.R")
if (!exists("render_fig06c_cytonuclear_paralog_cytotype_png", mode = "function")) .src_step06_rel("scripts/figures/step06c_figures.R")

.resolve_step06_db_env <- function(default = "") {
  for (k in c("PLANTMATICS_DB_PATH", "PACKERA_SQLITE_DB", "PACKERA_DB_PATH")) {
    v <- Sys.getenv(k, unset = "")
    if (nzchar(trimws(v))) return(trimws(v))
  }
  default
}

parse_step06_args <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  parser <- optparse::OptionParser()
  parser <- optparse::add_option(parser, c("--genomics-dir"), type = "character", default = "data/genomics")
  parser <- optparse::add_option(parser, c("--raw-genomics-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--parsed-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--genomic-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--out-dir"), type = "character", default = "outputs/genomics")
  parser <- optparse::add_option(parser, c("--tables-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--output-tables"), type = "character", default = "outputs/tables")
  parser <- optparse::add_option(parser, c("--output-figures"), type = "character", default = "outputs/figures")
  parser <- optparse::add_option(parser, c("--figure-a"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--figure-b"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--figure-c"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--run-id"), type = "character", default = Sys.getenv("RUN_ID", unset = "canonical"))
  parser <- optparse::add_option(parser, c("--db"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--db-path"), type = "character", default = .resolve_step06_db_env(""))
  parser <- optparse::add_option(parser, c("--demo"), action = "store_true", default = FALSE)
  optparse::parse_args(parser, args = cli_args)
}

.clean_tip_id <- function(x) {
  gsub("_sorted$", "", trimws(as.character(x)))
}

.json_escape_kv <- function(k, v) {
  k_c <- gsub('["\\\\]', "", as.character(k))
  sprintf('"%s": %.4f', k_c, as.numeric(v))
}

.compute_kappa_r <- function(y_true, y_pred) {
  keep <- !is.na(y_true) & nzchar(y_true) & !is.na(y_pred) & nzchar(y_pred)
  yt <- as.character(y_true[keep])
  yp <- as.character(y_pred[keep])
  if (length(yt) < 3L) {
    return(NA_real_)
  }
  lv <- sort(unique(c(yt, yp)))
  cm <- table(factor(yt, levels = lv), factor(yp, levels = lv))
  tot <- sum(cm)
  if (tot == 0) {
    return(NA_real_)
  }
  po <- sum(diag(cm)) / tot
  pe <- sum(rowSums(cm) * colSums(cm)) / (tot^2)
  if ((1.0 - pe) <= 1e-9) {
    return(ifelse(po >= 0.999, 1.0, 0.0))
  }
  max(-1.0, min(1.0, (po - pe) / (1.0 - pe)))
}

.compute_allele_sharing_distance <- function(dosage_mat) {
  n <- nrow(dosage_mat)
  p <- ncol(dosage_mat)
  d_mat <- matrix(0.0, nrow = n, ncol = n, dimnames = list(rownames(dosage_mat), rownames(dosage_mat)))
  n_shared_mat <- matrix(0L, nrow = n, ncol = n)
  min_shared <- if (p >= 200L) 150L else max(10L, as.integer(floor(0.20 * p)))

  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      vi <- dosage_mat[i, ]
      mi <- is.finite(vi)
      for (j in (i + 1L):n) {
        vj <- dosage_mat[j, ]
        shared <- mi & is.finite(vj)
        ns <- sum(shared)
        n_shared_mat[i, j] <- ns
        n_shared_mat[j, i] <- ns
        if (ns >= min_shared) {
          dval <- mean(abs(vi[shared] - vj[shared]))
        } else {
          dval <- NA_real_
        }
        d_mat[i, j] <- dval
        d_mat[j, i] <- dval
      }
    }
  }

  d_raw <- d_mat
  if (any(!is.finite(d_mat))) {
    global_mu <- mean(d_mat[upper.tri(d_mat)], na.rm = TRUE)
    if (!is.finite(global_mu)) global_mu <- 0.25
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j && !is.finite(d_mat[i, j])) {
          via_k <- d_mat[i, ] + d_mat[, j]
          via_k <- via_k[is.finite(via_k)]
          fill_v <- if (length(via_k) > 0L) min(max(via_k, na.rm = TRUE), mean(via_k)) else global_mu
          d_mat[i, j] <- fill_v
          d_mat[j, i] <- fill_v
        }
      }
    }
  }
  diag(d_mat) <- 0.0
  list(dist = d_mat, raw_dist = d_raw, shared = n_shared_mat, min_shared = min_shared)
}

.project_simplex <- function(v) {
  k <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1.0
  rho <- max(which((u - cssv / seq_len(k)) > 0))
  theta <- cssv[rho] / rho
  pmax(v - theta, 0.0)
}

.solve_admixture_simplex <- function(x_mat, centroids_mat, taxa_names) {
  n <- nrow(x_mat)
  g <- nrow(centroids_mat)
  q_mat <- matrix(0.0, nrow = n, ncol = g, dimnames = list(rownames(x_mat), taxa_names))
  if (g == 1L) {
    q_mat[, 1L] <- 1.0
    return(q_mat)
  }
  Mt <- t(centroids_mat) # m x G
  MtM <- crossprod(Mt) + diag(1e-4, g)
  L_lip <- max(eigen(MtM, symmetric = TRUE, only.values = TRUE)$values, 1e-4)
  step_sz <- 1.0 / L_lip

  for (i in seq_len(n)) {
    xi <- as.numeric(x_mat[i, ])
    # Restrict projection cleanly to top-3 closest centroids for sparsity + speed, or full simplex
    w <- rep(1.0 / g, g)
    Mtx <- as.numeric(crossprod(Mt, xi))
    for (it in seq_len(80L)) {
      grad <- as.numeric(MtM %*% w) - Mtx
      w_next <- .project_simplex(w - step_sz * grad)
      if (sum(abs(w_next - w)) < 1e-6) {
        w <- w_next
        break
      }
      w <- w_next
    }
    w[w < 0.02] <- 0.0
    if (sum(w) > 0) w <- w / sum(w) else w <- rep(1.0 / g, g)
    q_mat[i, ] <- w
  }
  q_mat
}

.fit_shrinkage_dapc <- function(X_train, y_train, X_test, lambda_shrink = 0.20) {
  classes <- sort(unique(y_train))
  g <- length(classes)
  m <- ncol(X_train)
  if (g < 2L || nrow(X_train) < 3L) {
    p_mat <- matrix(1.0 / max(1L, g), nrow = nrow(X_test), ncol = max(1L, g), dimnames = list(NULL, classes))
    return(list(
      pred = rep(classes[1L], nrow(X_test)),
      probs = p_mat,
      max_p = rep(1.0 / max(1L, g), nrow(X_test)),
      margin = rep(0.0, nrow(X_test)),
      centroids = matrix(colMeans(X_train), nrow = 1L)
    ))
  }

  centroids <- matrix(0.0, nrow = g, ncol = m, dimnames = list(classes, colnames(X_train)))
  S_w <- matrix(0.0, nrow = m, ncol = m)
  w_sum <- 0.0
  for (ci in seq_along(classes)) {
    idx <- which(y_train == classes[ci])
    sub <- X_train[idx, , drop = FALSE]
    mu <- colMeans(sub)
    centroids[ci, ] <- mu
    if (nrow(sub) >= 2L) {
      centered <- sweep(sub, 2L, mu, "-")
      S_w <- S_w + (crossprod(centered) / nrow(sub))
      w_sum <- w_sum + 1.0
    }
  }
  if (w_sum > 0) S_w <- S_w / w_sum else S_w <- stats::cov(X_train)
  if (!all(is.finite(S_w))) S_w <- diag(1.0, m)

  diag_sw <- diag(diag(S_w), nrow = m, ncol = m)
  mean_var <- max(1e-5, mean(diag(S_w)))
  Sigma_reg <- (1.0 - lambda_shrink) * S_w + lambda_shrink * diag_sw + diag(1e-3 * mean_var, m)
  Sigma_inv <- tryCatch(solve(Sigma_reg), error = function(e) solve(Sigma_reg + diag(1e-2, m)))

  log_scores <- matrix(0.0, nrow = nrow(X_test), ncol = g, dimnames = list(rownames(X_test), classes))
  # Scale temperature gently by feature dimensionality to avoid overconfident collapse on PCoA axes
  temp_scale <- max(1.0, sqrt(m))
  for (ci in seq_along(classes)) {
    diff_c <- sweep(X_test, 2L, centroids[ci, ], "-")
    quad <- rowSums((diff_c %*% Sigma_inv) * diff_c)
    log_scores[, ci] <- -0.5 * quad / temp_scale
  }

  probs <- matrix(0.0, nrow = nrow(X_test), ncol = g, dimnames = list(rownames(X_test), classes))
  for (r in seq_len(nrow(X_test))) {
    ls_r <- log_scores[r, ]
    ex <- exp(ls_r - max(ls_r))
    probs[r, ] <- ex / sum(ex)
  }

  pred <- classes[max.col(probs, ties.method = "first")]
  max_p <- apply(probs, 1L, max)
  sec_p <- apply(probs, 1L, function(v) sort(v, decreasing = TRUE)[min(2L, length(v))])
  margin <- if (g >= 2L) (max_p - sec_p) else max_p
  list(pred = pred, probs = probs, max_p = max_p, margin = margin, centroids = centroids)
}

.extract_phylo_metrics <- function(tree_files, pass_df) {
  if (length(tree_files) > 0L) {
    tree_files <- tree_files[!grepl("^public_packera_", basename(as.character(tree_files)), ignore.case = TRUE)]
  }
  out <- data.frame(
    genetic_sample_id = pass_df$genetic_sample_id,
    nearest_tree_sister_taxon = NA_character_,
    patristic_dist_to_dubia = NA_real_,
    stringsAsFactors = FALSE
  )
  summary_df <- data.frame(
    species_standardized = character(0),
    n_tips = integer(0),
    mean_intra_patristic = numeric(0),
    min_inter_patristic = numeric(0),
    nearest_sister_taxon = character(0),
    is_monophyletic = logical(0),
    stringsAsFactors = FALSE
  )
  if (length(tree_files) == 0L || !requireNamespace("ape", quietly = TRUE)) {
    return(list(sample_metrics = out, clade_summary = summary_df, tree = NULL))
  }

  tf <- tree_files[1L]
  if (!file.exists(tf)) {
    return(list(sample_metrics = out, clade_summary = summary_df, tree = NULL))
  }
  tr <- tryCatch(
    {
      if (grepl("\\.nex(us)?$", tf, ignore.case = TRUE)) ape::read.nexus(tf) else ape::read.tree(tf)
    },
    error = function(e) NULL
  )
  if (is.null(tr) || length(tr$tip.label) < 3L) {
    return(list(sample_metrics = out, clade_summary = summary_df, tree = NULL))
  }

  tr$tip.label <- .clean_tip_id(tr$tip.label)
  coph <- ape::cophenetic.phylo(tr)
  sp_map <- stats::setNames(pass_df$species_standardized, pass_df$genetic_sample_id)
  rep_map <- stats::setNames(pass_df$replicate_pair_id, pass_df$genetic_sample_id)
  qry_map <- stats::setNames(as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L, pass_df$genetic_sample_id)

  dubia_tips <- intersect(
    rownames(coph),
    pass_df$genetic_sample_id[
      !is.na(pass_df$species_standardized) &
        pass_df$species_standardized == "Packera dubia" &
        as.integer(pass_df$is_putative_novel_or_hybrid_query) == 0L
    ]
  )

  for (i in seq_len(nrow(pass_df))) {
    sid <- pass_df$genetic_sample_id[i]
    if (!(sid %in% rownames(coph))) next
    d_row <- coph[sid, ]
    if (length(dubia_tips) > 0L) {
      other_dub <- setdiff(dubia_tips, sid)
      out$patristic_dist_to_dubia[i] <- if (length(other_dub) > 0L) mean(d_row[other_dub], na.rm = TRUE) else 0.0
    }
    my_rep <- rep_map[[sid]]
    cand_tips <- setdiff(names(d_row), sid)
    if (!is.na(my_rep) && nzchar(my_rep)) {
      cand_tips <- cand_tips[is.na(rep_map[cand_tips]) | rep_map[cand_tips] != my_rep]
    }
    ref_cand <- cand_tips[!is.na(qry_map[cand_tips]) & !qry_map[cand_tips] & !is.na(sp_map[cand_tips])]
    use_cand <- if (length(ref_cand) > 0L) ref_cand else cand_tips[!is.na(sp_map[cand_tips])]
    if (length(use_cand) > 0L) {
      best_tip <- use_cand[which.min(d_row[use_cand])]
      out$nearest_tree_sister_taxon[i] <- sp_map[[best_tip]]
    }
  }

  # Clade cohesion summary per reference taxon
  present_tips <- intersect(tr$tip.label, pass_df$genetic_sample_id)
  ref_tips <- present_tips[!is.na(sp_map[present_tips]) & (!qry_map[present_tips] | grepl("^Packera sp\\. nov", sp_map[present_tips]))]
  taxa_u <- sort(unique(sp_map[ref_tips]))
  if (length(taxa_u) >= 2L) {
    sum_rows <- lapply(taxa_u, function(tx) {
      t_tx <- ref_tips[sp_map[ref_tips] == tx]
      t_oth <- setdiff(ref_tips, t_tx)
      intra <- if (length(t_tx) >= 2L) mean(coph[t_tx, t_tx][upper.tri(coph[t_tx, t_tx])], na.rm = TRUE) else 0.0
      inter_mat <- coph[t_tx, t_oth, drop = FALSE]
      min_inter <- min(inter_mat, na.rm = TRUE)
      min_idx <- which(inter_mat == min_inter, arr.ind = TRUE)[1L, 2L]
      sis_tx <- sp_map[[t_oth[min_idx]]]
      mono_ok <- if (length(t_tx) >= 2L) {
        tryCatch(ape::is.monophyletic(tr, tips = t_tx), error = function(e) NA)
      } else {
        NA
      }
      data.frame(
        species_standardized = tx,
        n_tips = length(t_tx),
        mean_intra_patristic = round(intra, 6),
        min_inter_patristic = round(min_inter, 6),
        nearest_sister_taxon = sis_tx,
        is_monophyletic = mono_ok,
        stringsAsFactors = FALSE
      )
    })
    summary_df <- dplyr::bind_rows(sum_rows)
  }

  list(sample_metrics = out, clade_summary = summary_df, tree = tr)
}

.read_admixture_q_if_present <- function(q_files, pass_ids, ref_taxa, pred_probs) {
  if (length(q_files) == 0L || !file.exists(q_files[1L])) {
    return(NULL)
  }
  q_df <- tryCatch(utils::read.table(q_files[1L], header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(q_df) || nrow(q_df) != length(pass_ids) || ncol(q_df) < 2L) {
    return(NULL)
  }
  if (is.numeric(q_df[[1L]])) {
    mat <- as.matrix(q_df)
    colnames(mat) <- paste0("Cluster_", seq_len(ncol(mat)))
    return(mat)
  }
  NULL
}

.compute_posthoc_congruence <- function(pass_df, d_geno_mat, pco_coords, db_path, run_id, demo) {
  congr_rows <- list()

  # 1. Isolation-by-Distance (Genomic vs Geodesic km) on G1/G2/G3 samples
  geo_idx <- which(is.finite(pass_df$decimalLatitude) & is.finite(pass_df$decimalLongitude))
  if (length(geo_idx) >= 4L) {
    lat_r <- pass_df$decimalLatitude[geo_idx] * pi / 180.0
    lon_r <- pass_df$decimalLongitude[geo_idx] * pi / 180.0
    ng <- length(geo_idx)
    d_geo <- matrix(0.0, nrow = ng, ncol = ng)
    for (i in seq_len(ng - 1L)) {
      for (j in (i + 1L):ng) {
        dlat <- lat_r[j] - lat_r[i]
        dlon <- lon_r[j] - lon_r[i]
        a <- sin(dlat / 2)^2 + cos(lat_r[i]) * cos(lat_r[j]) * sin(dlon / 2)^2
        km <- 2.0 * 6371.0 * asin(min(1.0, sqrt(max(0.0, a))))
        d_geo[i, j] <- km
        d_geo[j, i] <- km
      }
    }
    d_g_sub <- d_geno_mat[geo_idx, geo_idx, drop = FALSE]
    sd_g <- stats::sd(d_g_sub[upper.tri(d_g_sub)], na.rm = TRUE)
    sd_geo <- stats::sd(d_geo[upper.tri(d_geo)], na.rm = TRUE)
    has_var <- is.finite(sd_g) && sd_g > 0 && is.finite(sd_geo) && sd_geo > 0
    m_res <- if (has_var && requireNamespace("vegan", quietly = TRUE)) {
      tryCatch(vegan::mantel(stats::as.dist(d_g_sub), stats::as.dist(d_geo), permutations = 999L), error = function(e) NULL)
    } else {
      NULL
    }
    r_val <- if (!is.null(m_res)) as.numeric(m_res$statistic) else if (has_var) stats::cor(d_g_sub[upper.tri(d_g_sub)], d_geo[upper.tri(d_geo)]) else NA_real_
    p_val <- if (!is.null(m_res)) as.numeric(m_res$signif) else NA_real_
    congr_rows[[length(congr_rows) + 1L]] <- data.frame(
      comparison = "Genomic_vs_Geographic_IBD",
      n_samples = ng,
      mantel_r = round(r_val, 4),
      mantel_p = round(p_val, 4),
      procrustes_r = NA_real_,
      procrustes_p = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # 2. Genomic vs Morphometric / Ecological concordance via read-only DB lookup (strictly post-hoc)
  morph_linked_rows <- NULL
  if (!is.null(db_path) && nzchar(db_path) && file.exists(db_path)) {
    conn <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO), error = function(e) NULL)
    if (!is.null(conn)) {
      on.exit(tryCatch(DBI::dbDisconnect(conn), error = function(e) NULL), add = TRUE)
      cats <- unique(stats::na.omit(pass_df$catalogNumber))
      if (length(cats) >= 3L && db_table_exists_safe(conn, "otu_morphometrics")) {
        q_cats <- paste(vapply(cats, function(s) as.character(DBI::dbQuoteString(conn, s)), character(1)), collapse = ", ")
        sql_m <- sprintf(
          "SELECT catalogNumber, AVG(PC1) AS PC1, AVG(PC2) AS PC2, AVG(PC3) AS PC3 FROM otu_morphometrics WHERE run_id = %s AND catalogNumber IN (%s) GROUP BY catalogNumber;",
          as.character(DBI::dbQuoteString(conn, run_id)),
          q_cats
        )
        morph_linked_rows <- tryCatch(DBI::dbGetQuery(conn, sql_m), error = function(e) NULL)
      }
    }
  }

  if (!is.null(morph_linked_rows) && nrow(morph_linked_rows) >= 4L) {
    m_match <- match(pass_df$catalogNumber, morph_linked_rows$catalogNumber)
    ok_m <- which(!is.na(m_match) & is.finite(morph_linked_rows$PC1[m_match]) & is.finite(morph_linked_rows$PC2[m_match]))
    if (length(ok_m) >= 4L) {
      X_m <- as.matrix(morph_linked_rows[m_match[ok_m], c("PC1", "PC2", "PC3"), drop = FALSE])
      X_m[!is.finite(X_m)] <- 0.0
      X_g <- pco_coords[ok_m, seq_len(min(3L, ncol(pco_coords))), drop = FALSE]
      d_m <- as.matrix(stats::dist(X_m))
      d_g <- d_geno_mat[ok_m, ok_m, drop = FALSE]
      sd_gm <- stats::sd(d_g[upper.tri(d_g)], na.rm = TRUE)
      sd_mm <- stats::sd(d_m[upper.tri(d_m)], na.rm = TRUE)
      has_mvar <- is.finite(sd_gm) && sd_gm > 0 && is.finite(sd_mm) && sd_mm > 0
      m_gm <- if (has_mvar && requireNamespace("vegan", quietly = TRUE)) {
        tryCatch(vegan::mantel(stats::as.dist(d_g), stats::as.dist(d_m), permutations = 999L), error = function(e) NULL)
      } else {
        NULL
      }
      p_gm <- if (has_mvar && requireNamespace("vegan", quietly = TRUE)) {
        tryCatch(vegan::protest(X_g, X_m, permutations = 999L), error = function(e) NULL)
      } else {
        NULL
      }
      congr_rows[[length(congr_rows) + 1L]] <- data.frame(
        comparison = "Genomic_vs_Morphometrics_PCA",
        n_samples = length(ok_m),
        mantel_r = if (!is.null(m_gm)) round(as.numeric(m_gm$statistic), 4) else if (has_mvar) round(stats::cor(d_g[upper.tri(d_g)], d_m[upper.tri(d_m)]), 4) else NA_real_,
        mantel_p = if (!is.null(m_gm)) round(as.numeric(m_gm$signif), 4) else NA_real_,
        procrustes_r = if (!is.null(p_gm)) round(as.numeric(p_gm$t0), 4) else NA_real_,
        procrustes_p = if (!is.null(p_gm)) round(as.numeric(p_gm$signif), 4) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(congr_rows) == 0L) {
    return(data.frame(
      comparison = "Genomic_Only_No_Linked_Geo_Or_Morph_Pairs",
      n_samples = nrow(pass_df),
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      procrustes_r = NA_real_,
      procrustes_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(congr_rows)
}

.resolve_first_nonempty <- function(...) {
  for (v in list(...)) {
    if (!is.null(v) && nzchar(trimws(as.character(v)))) return(trimws(as.character(v)))
  }
  ""
}

run_step06_genomic_phylogenomics <- function(opts) {
  raw_geno_dir <- .resolve_first_nonempty(opts[["raw-genomics-dir"]], opts[["genomics-dir"]], "data/genomics")
  parsed_dir <- .resolve_first_nonempty(opts[["parsed-dir"]], opts[["genomic-dir"]], "outputs/reports")
  out_dir <- .resolve_first_nonempty(opts[["out-dir"]], "outputs/genomics")
  tables_dir <- .resolve_first_nonempty(opts[["tables-dir"]], opts[["output-tables"]], "")
  fig_dir <- .resolve_first_nonempty(opts[["output-figures"]], "outputs/figures")
  eff_db_in <- .resolve_first_nonempty(opts[["db"]], opts[["db-path"]], "")

  samples_csv <- file.path(parsed_dir, "parsed_genomic_samples.csv")
  dosage_tsv <- file.path(parsed_dir, "allele_dosage_matrix.tsv")
  manifest_json <- file.path(parsed_dir, "genomic_manifest.json")
  if (!file.exists(samples_csv) || !file.exists(dosage_tsv)) {
    message("[Step 06/07] No genomic files staged in ", parsed_dir, " - skipping standalone Step 06")
    return(invisible(NULL))
  }

  samples_df <- utils::read.csv(samples_csv, stringsAsFactors = FALSE)
  pass_df <- samples_df[!is.na(samples_df$geno_quality_pass) & as.integer(samples_df$geno_quality_pass) == 1L, , drop = FALSE]
  if (nrow(pass_df) < 3L) {
    message("[Step 06/07] Fewer than 3 QC-passing genomic samples - skipping standalone Step 06")
    return(invisible(NULL))
  }

  dos_df <- readr::read_tsv(dosage_tsv, show_col_types = FALSE, progress = FALSE)
  dos_df$genetic_sample_id <- as.character(dos_df$genetic_sample_id)
  common_ids <- intersect(pass_df$genetic_sample_id, dos_df$genetic_sample_id)
  pass_df <- pass_df[match(common_ids, pass_df$genetic_sample_id), , drop = FALSE]
  dos_sub <- dos_df[match(common_ids, dos_df$genetic_sample_id), , drop = FALSE]
  dos_mat <- as.matrix(dos_sub[, setdiff(names(dos_sub), "genetic_sample_id"), drop = FALSE])
  storage.mode(dos_mat) <- "numeric"
  rownames(dos_mat) <- common_ids

  dist_res <- .compute_allele_sharing_distance(dos_mat)
  d_mat <- dist_res$dist
  d_raw <- dist_res$raw_dist

  rep_conc <- rep(NA_real_, nrow(pass_df))
  rep_ids <- as.character(pass_df$replicate_pair_id)
  for (i in seq_len(nrow(pass_df))) {
    rid_i <- rep_ids[i]
    if (!is.na(rid_i) && nzchar(rid_i)) {
      mates <- setdiff(which(rep_ids == rid_i), i)
      if (length(mates) > 0L) {
        d_m <- d_raw[i, mates]
        if (any(is.finite(d_m))) rep_conc[i] <- max(0.0, min(1.0, 1.0 - mean(d_m[is.finite(d_m)])))
      }
    }
  }

  k_pco <- min(8L, max(2L, nrow(pass_df) - 1L))
  pco_fit <- stats::cmdscale(stats::as.dist(d_mat), k = k_pco, eig = TRUE)
  eig_pos <- pmax(0.0, pco_fit$eig)
  eig_pct <- if (sum(eig_pos) > 0) 100.0 * eig_pos / sum(eig_pos) else rep(100.0 / length(eig_pos), length(eig_pos))
  k_use <- min(6L, ncol(pco_fit$points), max(2L, sum(eig_pos > 1e-9)))
  pco_coords <- pco_fit$points[, seq_len(k_use), drop = FALSE]
  colnames(pco_coords) <- paste0("PCo", seq_len(k_use))
  pass_df$PCo1 <- pco_coords[, 1L]
  pass_df$PCo2 <- if (k_use >= 2L) pco_coords[, 2L] else 0.0

  best_k <- 2L; best_sil <- -1.0; best_clus <- rep(1L, nrow(pass_df))
  k_max <- min(8L, nrow(pass_df) - 1L)
  if (k_max >= 2L && requireNamespace("cluster", quietly = TRUE)) {
    for (kk in 2:k_max) {
      pam_k <- tryCatch(cluster::pam(stats::as.dist(d_mat), k = kk, diss = TRUE), error = function(e) NULL)
      if (!is.null(pam_k) && is.finite(pam_k$silinfo$avg.width) && pam_k$silinfo$avg.width > best_sil) {
        best_sil <- pam_k$silinfo$avg.width; best_k <- kk; best_clus <- pam_k$clustering
      }
    }
  }
  gen_cluster <- sprintf("Geno_Cluster_%d", as.integer(best_clus))

  is_ref_anchor <- (is.na(pass_df$is_putative_novel_or_hybrid_query) | as.integer(pass_df$is_putative_novel_or_hybrid_query) == 0L) &
    !is.na(pass_df$species_standardized) & nzchar(pass_df$species_standardized)
  ref_idx <- which(is_ref_anchor)
  ref_tbl <- table(pass_df$species_standardized[ref_idx])
  keep_taxa <- names(ref_tbl[ref_tbl >= 1L])
  if (length(keep_taxa) >= 2L) {
    ref_idx <- ref_idx[pass_df$species_standardized[ref_idx] %in% keep_taxa]
  } else {
    keep_taxa <- names(ref_tbl)
  }

  m_dapc <- min(5L, ncol(pco_coords), max(2L, length(ref_idx) - length(keep_taxa)))
  X_dapc <- pco_coords[, seq_len(m_dapc), drop = FALSE]
  y_ref <- pass_df$species_standardized[ref_idx]
  full_model <- .fit_shrinkage_dapc(X_dapc[ref_idx, , drop = FALSE], y_ref, X_dapc, lambda_shrink = 0.20)
  pred_taxon <- full_model$pred; post_prob <- full_model$max_p; margin_vec <- full_model$margin; probs_mat <- full_model$probs
  oof_pred <- rep(NA_character_, length(ref_idx))
  for (ri in seq_along(ref_idx)) {
    i_global <- ref_idx[ri]
    rep_i <- rep_ids[i_global]
    hold_mask <- if (!is.na(rep_i) && nzchar(rep_i)) (rep_ids[ref_idx] == rep_i) else (seq_along(ref_idx) == ri)
    tr_sub <- which(!hold_mask)
    if (sum(y_ref[tr_sub] == y_ref[i_global], na.rm = TRUE) >= 1L && length(unique(y_ref[tr_sub])) >= 2L) {
      cv_fit <- .fit_shrinkage_dapc(X_dapc[ref_idx[tr_sub], , drop = FALSE], y_ref[tr_sub], X_dapc[i_global, , drop = FALSE], lambda_shrink = 0.20)
      oof_pred[ri] <- cv_fit$pred[1L]; pred_taxon[i_global] <- cv_fit$pred[1L]
      post_prob[i_global] <- cv_fit$max_p[1L]; margin_vec[i_global] <- cv_fit$margin[1L]
      probs_mat[i_global, ] <- 0.0
      common_cols <- intersect(colnames(probs_mat), colnames(cv_fit$probs))
      if (length(common_cols) > 0L) probs_mat[i_global, common_cols] <- cv_fit$probs[1L, common_cols]
    }
  }
  valid_oof <- !is.na(oof_pred)
  loro_acc <- if (any(valid_oof)) mean(oof_pred[valid_oof] == y_ref[valid_oof], na.rm = TRUE) else NA_real_
  loro_kappa <- if (sum(valid_oof) >= 2L) .compute_kappa_r(y_ref[valid_oof], oof_pred[valid_oof]) else NA_real_

  tree_files <- character(0); q_files <- character(0)
  if (file.exists(manifest_json)) {
    m_txt <- gsub("[\r\n]+", " ", paste(readLines(manifest_json, warn = FALSE), collapse = " "))
    ext_paths <- function(key) {
      m <- regexec(sprintf('"%s"\\s*:\\s*\\[(.*?)\\]', key), m_txt)
      hit <- regmatches(m_txt, m)[[1L]]
      if (length(hit) < 2L) return(character(0))
      vals <- regmatches(hit[2L], gregexpr('"([^"]+)"', hit[2L]))[[1L]]
      gsub('^"|"$', "", vals)
    }
    tree_files <- ext_paths("tree_files"); q_files <- ext_paths("admixture_q_files")
  }
  if (length(tree_files) > 0L) {
    tree_files <- tree_files[!grepl("^public_packera_", basename(as.character(tree_files)), ignore.case = TRUE)]
  }
  if (length(tree_files) == 0L && dir.exists(raw_geno_dir)) {
    cand_tf <- sort(list.files(raw_geno_dir, full.names = TRUE, recursive = TRUE))
    b_tf <- basename(cand_tf)
    tree_files <- cand_tf[grepl("^RAxML_bestTree|\\.(treefile|tre|nwk|tree)$", b_tf, ignore.case = TRUE) &
                          !grepl("^public_packera_", b_tf, ignore.case = TRUE)]
  }

  q_ext <- .read_admixture_q_if_present(q_files, common_ids, keep_taxa, probs_mat)
  q_mat <- if (!is.null(q_ext)) q_ext else .solve_admixture_simplex(X_dapc, full_model$centroids, rownames(full_model$centroids))
  q_cols <- colnames(q_mat); q_max <- apply(q_mat, 1L, max)
  is_admixed <- integer(nrow(pass_df)); p1 <- character(nrow(pass_df)); p2 <- character(nrow(pass_df)); q_json <- character(nrow(pass_df))
  for (i in seq_len(nrow(pass_df))) {
    wi <- q_mat[i, ]; ord_w <- order(wi, decreasing = TRUE); nz <- which(wi >= 0.02)
    q_json[i] <- paste0("{", paste(vapply(nz[order(wi[nz], decreasing = TRUE)], function(j) .json_escape_kv(q_cols[j], wi[j]), character(1)), collapse = ", "), "}")
    p1[i] <- q_cols[ord_w[1L]]
    p2[i] <- if (length(ord_w) >= 2L && wi[ord_w[2L]] >= 0.15) q_cols[ord_w[2L]] else NA_character_
    is_admixed[i] <- if (q_max[i] < 0.80 && !is.na(p2[i])) 1L else 0L
  }
  phylo_res <- .extract_phylo_metrics(tree_files, pass_df)

  # Phase 2 Independent Modules: Co-ancestry, (h, Q12) Hybrid Classifier, f-branch, SNPs2CF/SNaQ
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pub_bb <- load_public_packera_backbone(raw_geno_dir)
  coanc_res <- compute_microhaplotype_coancestry(common_ids, dos_mat, raw_geno_dir)
  tri_res <- compute_hybrid_index_q12(pass_df, dos_mat, q_mat, p1, p2, d_mat)
  fb_res <- compute_tree_guided_fbranch(pass_df, dos_mat, raw_geno_dir, backbone = pub_bb)
  cf_res <- compute_quartet_concordance_snps2cf(pass_df, dos_mat)
  pub_eval <- evaluate_public_backbone_and_organellar(pass_df, raw_geno_dir, out_dir)
  # Ensure divergent cryptic lineages are never falsely marked admixed when classified Divergent_NonHybrid_Lineage
  is_admixed[tri_res$hybrid_class == "Divergent_NonHybrid_Lineage"] <- 0L

  # Phase 2 Invariant-site Coalescent Delimitation (SVDquartets, BPP MSci, gdi)
  py_bin <- if (file.exists(".venv/bin/python")) ".venv/bin/python" else Sys.which("python3")
  coal_script <- "scripts/analysis/step06_coalescent_delimitation.py"
  if (nzchar(py_bin) && file.exists(coal_script)) {
    c_args <- c(
      shQuote(coal_script),
      "--genomics-dir", shQuote(raw_geno_dir),
      "--out-dir", shQuote(out_dir),
      "--parsed-samples", shQuote(samples_csv),
      "--dosage-tsv", shQuote(dosage_tsv),
      "--n-boot", "250"
    )
    tryCatch(system2(py_bin, args = c_args, stdout = TRUE, stderr = TRUE), error = function(e) warning("[Step 06 Coalescent] ", e$message))
  }
  gdi_samp_csv <- file.path(out_dir, "coalescent_gdi_samples.csv")
  gdi_tax_csv <- file.path(out_dir, "coalescent_gdi_delimitation.csv")
  gdi_df <- if (file.exists(gdi_tax_csv)) utils::read.csv(gdi_tax_csv, stringsAsFactors = FALSE) else data.frame()
  gdi_s_df <- if (file.exists(gdi_samp_csv)) utils::read.csv(gdi_samp_csv, stringsAsFactors = FALSE) else data.frame()

  gdi_vec <- rep(NA_real_, nrow(pass_df)); bf_vec <- rep(NA_real_, nrow(pass_df))
  if (nrow(gdi_s_df) > 0L && "genetic_sample_id" %in% names(gdi_s_df)) {
    gm <- match(pass_df$genetic_sample_id, gdi_s_df$genetic_sample_id)
    if ("coalescent_gdi" %in% names(gdi_s_df)) gdi_vec <- suppressWarnings(as.numeric(gdi_s_df$coalescent_gdi[gm]))
    if ("bpp_introgression_bf" %in% names(gdi_s_df)) bf_vec <- suppressWarnings(as.numeric(gdi_s_df$bpp_introgression_bf[gm]))
  }

  assign_df <- data.frame(
    genetic_sample_id = pass_df$genetic_sample_id,
    catalogNumber = ifelse(is.na(pass_df$catalogNumber) | !nzchar(pass_df$catalogNumber), NA_character_, pass_df$catalogNumber),
    population_id = pass_df$population_id,
    geno_predicted_taxon = pred_taxon,
    geno_posterior_prob = round(post_prob, 4),
    geno_margin = round(margin_vec, 4),
    genetic_cluster = gen_cluster,
    q_max = round(q_max, 4),
    q_admixture_json = q_json,
    is_geno_admixed = is_admixed,
    geno_parent_1 = p1,
    geno_parent_2 = p2,
    nearest_tree_sister_taxon = phylo_res$sample_metrics$nearest_tree_sister_taxon,
    patristic_dist_to_dubia = round(phylo_res$sample_metrics$patristic_dist_to_dubia, 6),
    replicate_pair_concordance = round(rep_conc, 4),
    hybrid_index_h = tri_res$hybrid_index_h,
    interspecific_ancestry_q12 = tri_res$q12,
    hybrid_generation_class = tri_res$hybrid_class,
    coancestry_cluster = coanc_res$cluster,
    f_branch_max = fb_res$f_branch_max,
    f_branch_donor_taxon = fb_res$f_branch_donor_taxon,
    snaq_gamma_minor = cf_res$sample_gamma,
    coalescent_gdi = gdi_vec,
    bpp_introgression_bf = bf_vec,
    mmrr_beta_ibd = NA_real_,
    mmrr_beta_ibe = NA_real_,
    mmrr_p_ibe = NA_real_,
    public_backbone_clade = pub_eval$public_backbone_clade,
    public_phytop_ils_index = pub_eval$public_phytop_ils_index,
    public_phytop_ih_index = pub_eval$public_phytop_ih_index,
    stringsAsFactors = FALSE
  )

  # Persist Tier B SQLite BEFORE reading any Step 04/05 data (strict modality independence)
  if (!isTRUE(opts$demo) && exists("persist_step06_genomics_to_sqlite", mode = "function")) {
    tryCatch(
      persist_step06_genomics_to_sqlite(assign_df, run_id = opts[["run-id"]], db_path = if (nzchar(eff_db_in)) eff_db_in else NULL, demo = FALSE),
      error = function(e) warning("[Step 06] SQLite persistence warning: ", e$message)
    )
  }

  eff_db_ro <- if (nzchar(eff_db_in)) eff_db_in else .resolve_step06_db_env("data/plantmatics.sqlite")

  # Post-stream read-only landscape genomics (MMRR + FEEMS) & Mantel congruence
  land_res <- run_step06_landscape_genomics(pass_df, d_mat, db_path = eff_db_ro, run_id = opts[["run-id"]], demo = opts$demo)
  if (nrow(land_res$summary_df) > 0L) {
    assign_df$mmrr_beta_ibd <- land_res$summary_df$mmrr_beta_ibd[1L]
    assign_df$mmrr_beta_ibe <- land_res$summary_df$mmrr_beta_ibe[1L]
    assign_df$mmrr_p_ibe <- land_res$summary_df$mmrr_p_ibe[1L]
    if (!isTRUE(opts$demo) && exists("persist_step06_genomics_to_sqlite", mode = "function")) {
      tryCatch(persist_step06_genomics_to_sqlite(assign_df, run_id = opts[["run-id"]], db_path = if (nzchar(eff_db_in)) eff_db_in else NULL, demo = FALSE), error = function(e) NULL)
    }
  }
  congr_df <- .compute_posthoc_congruence(pass_df, d_mat, pco_coords, eff_db_ro, opts[["run-id"]], opts$demo)

  full_audit <- dplyr::left_join(samples_df, assign_df[, setdiff(names(assign_df), c("catalogNumber", "population_id")), drop = FALSE], by = "genetic_sample_id")
  m_pass <- match(full_audit$genetic_sample_id, pass_df$genetic_sample_id)
  for (mcol in c("public_plastome_match", "public_its_match", "waxy_lineage_match")) {
    if (!(mcol %in% names(full_audit))) full_audit[[mcol]] <- NA_character_
    need_fill <- !is.na(m_pass) & (is.na(full_audit[[mcol]]) | !nzchar(as.character(full_audit[[mcol]])))
    if (any(need_fill)) full_audit[[mcol]][need_fill] <- pub_eval[[mcol]][m_pass[need_fill]]
  }
  metrics_df <- data.frame(
    n_samples_total = nrow(samples_df), n_samples_pass = nrow(pass_df),
    n_reference_anchors = length(ref_idx),
    n_novel_queries = sum(as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L, na.rm = TRUE),
    n_loci_retained = ncol(dos_mat), optimal_k_pam = best_k, pam_silhouette = round(best_sil, 4),
    loro_cv_accuracy = round(loro_acc, 4), loro_cv_kappa = round(loro_kappa, 4),
    mean_replicate_concordance = round(mean(rep_conc, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )

  readr::write_csv(full_audit, file.path(out_dir, "genomic_assignments.csv"), na = "")
  readr::write_csv(phylo_res$clade_summary, file.path(out_dir, "genomic_clade_cohesion_summary.csv"), na = "")
  readr::write_csv(metrics_df, file.path(out_dir, "genomic_model_metrics.csv"), na = "")
  readr::write_csv(congr_df, file.path(out_dir, "genomic_multimodal_congruence_summary.csv"), na = "")
  readr::write_csv(coanc_res$coanc_df, file.path(out_dir, "coancestry_matrix_finerad.csv"), na = "")
  readr::write_csv(tri_res$triangle_df, file.path(out_dir, "hybrid_index_q12_triangle.csv"), na = "")
  readr::write_csv(fb_res$fbranch_df, file.path(out_dir, "introgressive_fbranch_summary.csv"), na = "")
  readr::write_csv(cf_res$cf_df, file.path(out_dir, "quartet_concordance_factors_snps2cf.csv"), na = "")
  readr::write_csv(land_res$summary_df, file.path(out_dir, "landscape_mmrr_feems_summary.csv"), na = "")
  if (nrow(land_res$feems_edges_df) > 0L) readr::write_csv(land_res$feems_edges_df, file.path(out_dir, "feems_migration_edges.csv"), na = "")

  if (nzchar(tables_dir) && normalizePath(tables_dir, mustWork = FALSE) != normalizePath(out_dir, mustWork = FALSE)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(full_audit, file.path(tables_dir, "genomic_assignments.csv"), na = "")
    readr::write_csv(phylo_res$clade_summary, file.path(tables_dir, "genomic_clade_cohesion_summary.csv"), na = "")
    readr::write_csv(metrics_df, file.path(tables_dir, "genomic_model_metrics.csv"), na = "")
    readr::write_csv(congr_df, file.path(tables_dir, "genomic_multimodal_congruence_summary.csv"), na = "")
    readr::write_csv(pub_eval$concordance_df, file.path(tables_dir, "public_backbone_concordance.csv"), na = "")
    readr::write_csv(pub_eval$placement_df, file.path(tables_dir, "public_organellar_placement.csv"), na = "")
  }

  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  fig_a <- .resolve_first_nonempty(opts[["figure-a"]], file.path(fig_dir, "fig06a_genomic_pcoa_tree_concordance.png"))
  fig_b <- .resolve_first_nonempty(opts[["figure-b"]], file.path(fig_dir, "fig06b_reticulation_gdi_landscape.png"))
  fig_c <- .resolve_first_nonempty(opts[["figure-c"]], file.path(fig_dir, "fig06c_cytonuclear_paralog_cytotype_audit.png"))
  pass_df$coancestry_cluster <- coanc_res$cluster
  render_fig06a_genomic_pcoa_tree_png(fig_a, pass_df, assign_df, q_mat, eig_pct, phylo_res$tree, congr_df, eff_db_ro, opts[["run-id"]])
  render_fig06b_reticulation_gdi_png(fig_b, tri_res$triangle_df, fb_res$fbranch_df, gdi_df, land_res$feems_edges_df, land_res$summary_df, organellar_placement_df = pub_eval$placement_df)
  render_fig06c_cytonuclear_paralog_cytotype_png(fig_c, parsed_dir, out_dir, full_audit, pub_concordance_df = pub_eval$concordance_df, organellar_placement_df = pub_eval$placement_df, db_path = eff_db_ro)

  message(sprintf(
    "[Step 06/07] Phylogenomics complete: %d/%d passing | LORO acc=%.3f kappa=%.3f | M/P conc=%.4f | median gdi=%.3f",
    nrow(pass_df), nrow(samples_df),
    ifelse(is.finite(loro_acc), loro_acc, NA_real_), ifelse(is.finite(loro_kappa), loro_kappa, NA_real_),
    ifelse(is.finite(mean(rep_conc, na.rm = TRUE)), mean(rep_conc, na.rm = TRUE), NA_real_),
    ifelse(any(is.finite(gdi_vec)), stats::median(gdi_vec, na.rm = TRUE), NA_real_)
  ))
  invisible(assign_df)
}

main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_step06_args(cli_args)
  run_step06_genomic_phylogenomics(opts)
}

if (sys.nframe() == 0L) {
  main()
}
