# ==============================================================================
# Module: scripts/analysis/step06_phylogenomics_utils.R
# Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Core Step 06 genomic distance, simplex admixture, shrinkage DAPC, tree
#   clade cohesion, post-hoc Mantel/Procrustes congruence (quarantining centroids
#   and unvouchered rooftops), and empirical cytonuclear discordance reconciliation.
# ==============================================================================

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
}

.clean_tip_id <- function(x) gsub("_sorted$", "", trimws(as.character(x)))

.json_escape_kv <- function(k, v) sprintf('"%s": %.4f', gsub('["\\\\]', "", as.character(k)), as.numeric(v))

.compute_kappa_r <- function(y_true, y_pred) {
  keep <- !is.na(y_true) & nzchar(y_true) & !is.na(y_pred) & nzchar(y_pred)
  yt <- as.character(y_true[keep]); yp <- as.character(y_pred[keep])
  if (length(yt) < 3L) return(NA_real_)
  lv <- sort(unique(c(yt, yp)))
  cm <- table(factor(yt, levels = lv), factor(yp, levels = lv)); tot <- sum(cm)
  if (tot == 0) return(NA_real_)
  po <- sum(diag(cm)) / tot; pe <- sum(rowSums(cm) * colSums(cm)) / (tot^2)
  if ((1.0 - pe) <= 1e-9) return(ifelse(po >= 0.999, 1.0, 0.0))
  max(-1.0, min(1.0, (po - pe) / (1.0 - pe)))
}

.compute_allele_sharing_distance <- function(dosage_mat) {
  n <- nrow(dosage_mat); p <- ncol(dosage_mat)
  d_mat <- matrix(0.0, nrow = n, ncol = n, dimnames = list(rownames(dosage_mat), rownames(dosage_mat)))
  n_shared_mat <- matrix(0L, nrow = n, ncol = n)
  min_shared <- if (p >= 200L) 150L else max(10L, as.integer(floor(0.20 * p)))
  if (n >= 2L) {
    eval_row_i <- function(i) {
      vi <- dosage_mat[i, ]; mi <- is.finite(vi); j_seq <- (i + 1L):n
      ns_v <- integer(length(j_seq)); dv_v <- rep(NA_real_, length(j_seq))
      for (k in seq_along(j_seq)) {
        vj <- dosage_mat[j_seq[k], ]; sh <- mi & is.finite(vj); ns <- sum(sh); ns_v[k] <- ns
        if (ns >= min_shared) dv_v[k] <- mean(abs(vi[sh] - vj[sh]))
      }
      list(i = i, j_seq = j_seq, ns = ns_v, d = dv_v)
    }
    i_seq <- seq_len(n - 1L)
    nw_d <- if (n >= 32L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, length(i_seq)) else 1L
    row_res <- if (nw_d > 1L) packera_mclapply(i_seq, eval_row_i, mc.cores = nw_d) else lapply(i_seq, eval_row_i)
    for (rr in row_res) {
      i <- rr$i; js <- rr$j_seq
      n_shared_mat[i, js] <- rr$ns; n_shared_mat[js, i] <- rr$ns
      d_mat[i, js] <- rr$d; d_mat[js, i] <- rr$d
    }
  }
  d_raw <- d_mat
  if (any(!is.finite(d_mat))) {
    global_mu <- mean(d_mat[upper.tri(d_mat)], na.rm = TRUE); if (!is.finite(global_mu)) global_mu <- 0.25
    for (i in seq_len(n)) for (j in seq_len(n)) if (i != j && !is.finite(d_mat[i, j])) {
      via_k <- d_mat[i, ] + d_mat[, j]; via_k <- via_k[is.finite(via_k)]
      fill_v <- if (length(via_k) > 0L) min(max(via_k, na.rm = TRUE), mean(via_k)) else global_mu
      d_mat[i, j] <- fill_v; d_mat[j, i] <- fill_v
    }
  }
  diag(d_mat) <- 0.0
  list(dist = d_mat, raw_dist = d_raw, shared = n_shared_mat, min_shared = min_shared)
}

.project_simplex <- function(v) {
  k <- length(v); u <- sort(v, decreasing = TRUE); cssv <- cumsum(u) - 1.0
  rho <- max(which((u - cssv / seq_len(k)) > 0))
  pmax(v - cssv[rho] / rho, 0.0)
}

.solve_admixture_simplex <- function(x_mat, centroids_mat, taxa_names) {
  n <- nrow(x_mat); g <- nrow(centroids_mat)
  q_mat <- matrix(0.0, nrow = n, ncol = g, dimnames = list(rownames(x_mat), taxa_names))
  if (g == 1L) { q_mat[, 1L] <- 1.0; return(q_mat) }
  Mt <- t(centroids_mat); MtM <- crossprod(Mt) + diag(1e-4, g)
  step_sz <- 1.0 / max(eigen(MtM, symmetric = TRUE, only.values = TRUE)$values, 1e-4)
  for (i in seq_len(n)) {
    w <- rep(1.0 / g, g); Mtx <- as.numeric(crossprod(Mt, as.numeric(x_mat[i, ])))
    for (it in seq_len(80L)) {
      w_next <- .project_simplex(w - step_sz * (as.numeric(MtM %*% w) - Mtx))
      if (sum(abs(w_next - w)) < 1e-6) { w <- w_next; break }
      w <- w_next
    }
    w[w < 0.02] <- 0.0
    q_mat[i, ] <- if (sum(w) > 0) w / sum(w) else rep(1.0 / g, g)
  }
  q_mat
}

.fit_shrinkage_dapc <- function(X_train, y_train, X_test, lambda_shrink = NULL) {
  classes <- sort(unique(y_train)); g <- length(classes); m <- ncol(X_train)
  if (g < 2L || nrow(X_train) < 3L) {
    p_mat <- matrix(1.0 / max(1L, g), nrow = nrow(X_test), ncol = max(1L, g), dimnames = list(NULL, classes))
    return(list(pred = rep(classes[1L], nrow(X_test)), probs = p_mat, max_p = rep(1.0 / max(1L, g), nrow(X_test)),
                margin = rep(0.0, nrow(X_test)), centroids = matrix(colMeans(X_train), nrow = 1L), lambda_shrink = 0.95))
  }
  centroids <- matrix(0.0, nrow = g, ncol = m, dimnames = list(classes, colnames(X_train)))
  S_w <- matrix(0.0, nrow = m, ncol = m); w_sum <- 0.0; wc_list <- vector("list", g)
  for (ci in seq_along(classes)) {
    sub <- X_train[which(y_train == classes[ci]), , drop = FALSE]; mu <- colMeans(sub); centroids[ci, ] <- mu
    if (nrow(sub) >= 2L) { centered <- sweep(sub, 2L, mu, "-"); wc_list[[ci]] <- centered; S_w <- S_w + crossprod(centered) / nrow(sub); w_sum <- w_sum + 1.0 }
  }
  S_w <- if (w_sum > 0) S_w / w_sum else stats::cov(X_train)
  if (!all(is.finite(S_w))) S_w <- diag(1.0, m)
  wc_nonempty <- Filter(Negate(is.null), wc_list)
  X_tr_c <- if (length(wc_nonempty) > 0L) do.call(rbind, wc_nonempty) else scale(X_train, center = TRUE, scale = FALSE)
  if (is.null(nrow(X_tr_c)) || nrow(X_tr_c) < 2L) X_tr_c <- scale(X_train, center = TRUE, scale = FALSE)
  lam_eff <- if (is.null(lambda_shrink) || !is.finite(lambda_shrink)) {
    if (exists("estimate_ledoit_wolf_lambda", mode = "function")) estimate_ledoit_wolf_lambda(X_tr_c) else 0.20
  } else as.numeric(lambda_shrink)
  Sigma_reg <- (1.0 - lam_eff) * S_w + lam_eff * diag(diag(S_w), nrow = m, ncol = m) + diag(1e-3 * max(1e-5, mean(diag(S_w))), m)
  Sigma_inv <- tryCatch(solve(Sigma_reg), error = function(e) solve(Sigma_reg + diag(1e-2, m)))
  log_scores <- matrix(0.0, nrow = nrow(X_test), ncol = g, dimnames = list(rownames(X_test), classes)); temp_scale <- max(1.0, sqrt(m))
  for (ci in seq_along(classes)) {
    diff_c <- sweep(X_test, 2L, centroids[ci, ], "-")
    log_scores[, ci] <- -0.5 * rowSums((diff_c %*% Sigma_inv) * diff_c) / temp_scale
  }
  probs <- matrix(0.0, nrow = nrow(X_test), ncol = g, dimnames = list(rownames(X_test), classes))
  for (r in seq_len(nrow(X_test))) { ex <- exp(log_scores[r, ] - max(log_scores[r, ])); probs[r, ] <- ex / sum(ex) }
  max_p <- apply(probs, 1L, max); sec_p <- apply(probs, 1L, function(v) sort(v, decreasing = TRUE)[min(2L, length(v))])
  list(pred = classes[max.col(probs, ties.method = "first")], probs = probs, max_p = max_p, margin = if (g >= 2L) (max_p - sec_p) else max_p, centroids = centroids, lambda_shrink = lam_eff)
}

.evaluate_tree_clade_cohesion <- function(tr, pass_df) {
  summary_df <- data.frame(species_standardized = character(0), n_tips = integer(0), mean_intra_patristic = numeric(0),
                           min_inter_patristic = numeric(0), nearest_sister_taxon = character(0), is_monophyletic = logical(0), stringsAsFactors = FALSE)
  if (is.null(tr) || length(tr$tip.label) < 3L || !requireNamespace("ape", quietly = TRUE)) return(summary_df)
  coph <- ape::cophenetic.phylo(tr)
  sp_map <- stats::setNames(pass_df$species_standardized, pass_df$genetic_sample_id)
  qry_map <- stats::setNames(as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L, pass_df$genetic_sample_id)
  present_tips <- intersect(tr$tip.label, pass_df$genetic_sample_id)
  ref_tips <- present_tips[!is.na(sp_map[present_tips]) & (!qry_map[present_tips] | grepl("^Packera sp\\. nov", sp_map[present_tips]))]
  taxa_u <- sort(unique(sp_map[ref_tips]))
  if (length(taxa_u) >= 2L) {
    sum_rows <- lapply(taxa_u, function(tx) {
      t_tx <- ref_tips[sp_map[ref_tips] == tx]; t_oth <- setdiff(ref_tips, t_tx)
      intra <- if (length(t_tx) >= 2L) mean(coph[t_tx, t_tx][upper.tri(coph[t_tx, t_tx])], na.rm = TRUE) else 0.0
      inter_mat <- coph[t_tx, t_oth, drop = FALSE]; min_inter <- min(inter_mat, na.rm = TRUE)
      sis_tx <- sp_map[[t_oth[which(inter_mat == min_inter, arr.ind = TRUE)[1L, 2L]]]]
      mono_ok <- if (length(t_tx) >= 2L) tryCatch(ape::is.monophyletic(tr, tips = t_tx), error = function(e) NA) else NA
      data.frame(species_standardized = tx, n_tips = length(t_tx), mean_intra_patristic = round(intra, 6),
                 min_inter_patristic = round(min_inter, 6), nearest_sister_taxon = sis_tx, is_monophyletic = mono_ok, stringsAsFactors = FALSE)
    })
    summary_df <- dplyr::bind_rows(sum_rows)
  }
  summary_df
}

.extract_phylo_metrics <- function(tree_files, pass_df) {
  if (length(tree_files) > 0L) tree_files <- tree_files[!grepl("^public_packera_", basename(as.character(tree_files)), ignore.case = TRUE)]
  out <- data.frame(genetic_sample_id = pass_df$genetic_sample_id, nearest_tree_sister_taxon = NA_character_, patristic_dist_to_dubia = NA_real_, stringsAsFactors = FALSE)
  summary_df <- data.frame(species_standardized = character(0), n_tips = integer(0), mean_intra_patristic = numeric(0),
                           min_inter_patristic = numeric(0), nearest_sister_taxon = character(0), is_monophyletic = logical(0), stringsAsFactors = FALSE)
  if (length(tree_files) == 0L || !requireNamespace("ape", quietly = TRUE) || !file.exists(tree_files[1L])) {
    return(list(sample_metrics = out, clade_summary = summary_df, tree = NULL))
  }
  tf <- tree_files[1L]
  tr <- tryCatch({ if (grepl("\\.nex(us)?$", tf, ignore.case = TRUE)) ape::read.nexus(tf) else ape::read.tree(tf) }, error = function(e) NULL)
  if (is.null(tr) || length(tr$tip.label) < 3L) return(list(sample_metrics = out, clade_summary = summary_df, tree = NULL))
  tr$tip.label <- .clean_tip_id(tr$tip.label); coph <- ape::cophenetic.phylo(tr)
  sp_map <- stats::setNames(pass_df$species_standardized, pass_df$genetic_sample_id)
  rep_map <- stats::setNames(pass_df$replicate_pair_id, pass_df$genetic_sample_id)
  qry_map <- stats::setNames(as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L, pass_df$genetic_sample_id)
  dubia_tips <- intersect(rownames(coph), pass_df$genetic_sample_id[!is.na(pass_df$species_standardized) & pass_df$species_standardized == "Packera dubia" & as.integer(pass_df$is_putative_novel_or_hybrid_query) == 0L])
  for (i in seq_len(nrow(pass_df))) {
    sid <- pass_df$genetic_sample_id[i]; if (!(sid %in% rownames(coph))) next
    d_row <- coph[sid, ]
    if (length(dubia_tips) > 0L) {
      other_dub <- setdiff(dubia_tips, sid)
      out$patristic_dist_to_dubia[i] <- if (length(other_dub) > 0L) mean(d_row[other_dub], na.rm = TRUE) else 0.0
    }
    my_rep <- rep_map[[sid]]; cand_tips <- setdiff(names(d_row), sid)
    if (!is.na(my_rep) && nzchar(my_rep)) cand_tips <- cand_tips[is.na(rep_map[cand_tips]) | rep_map[cand_tips] != my_rep]
    ref_cand <- cand_tips[!is.na(qry_map[cand_tips]) & !qry_map[cand_tips] & !is.na(sp_map[cand_tips])]
    use_cand <- if (length(ref_cand) > 0L) ref_cand else cand_tips[!is.na(sp_map[cand_tips])]
    if (length(use_cand) > 0L) out$nearest_tree_sister_taxon[i] <- sp_map[[use_cand[which.min(d_row[use_cand])]]]
  }
  list(sample_metrics = out, clade_summary = .evaluate_tree_clade_cohesion(tr, pass_df), tree = tr)
}

.read_admixture_q_if_present <- function(q_files, pass_ids, ref_taxa, pred_probs) {
  if (length(q_files) == 0L || !file.exists(q_files[1L])) return(NULL)
  q_df <- tryCatch(utils::read.table(q_files[1L], header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(q_df) || nrow(q_df) != length(pass_ids) || ncol(q_df) < 2L || !is.numeric(q_df[[1L]])) return(NULL)
  mat <- as.matrix(q_df); colnames(mat) <- paste0("Cluster_", seq_len(ncol(mat))); mat
}

.compute_posthoc_congruence <- function(pass_df, d_geno_mat, pco_coords, db_path, run_id, demo) {
  congr_rows <- list()
  lat_v <- suppressWarnings(as.numeric(pass_df$decimalLatitude)); lon_v <- suppressWarnings(as.numeric(pass_df$decimalLongitude))
  geo_idx <- if (exists(".filter_step06_valid_geo_idx", mode = "function")) {
    .filter_step06_valid_geo_idx(pass_df, lat_v, lon_v, db_path = db_path)
  } else which(is.finite(lat_v) & is.finite(lon_v) & !(abs(lat_v) < 1e-4 & abs(lon_v) < 1e-4))
  if (length(geo_idx) >= 4L) {
    lat_r <- lat_v[geo_idx] * pi / 180.0; lon_r <- lon_v[geo_idx] * pi / 180.0; ng <- length(geo_idx)
    d_geo <- matrix(0.0, nrow = ng, ncol = ng)
    for (i in seq_len(ng - 1L)) for (j in (i + 1L):ng) {
      a <- sin((lat_r[j] - lat_r[i]) / 2)^2 + cos(lat_r[i]) * cos(lat_r[j]) * sin((lon_r[j] - lon_r[i]) / 2)^2
      km <- 2.0 * 6371.0 * asin(min(1.0, sqrt(max(0.0, a)))); d_geo[i, j] <- km; d_geo[j, i] <- km
    }
    d_g_sub <- d_geno_mat[geo_idx, geo_idx, drop = FALSE]
    has_var <- is.finite(stats::sd(d_g_sub[upper.tri(d_g_sub)], na.rm = TRUE)) && stats::sd(d_g_sub[upper.tri(d_g_sub)], na.rm = TRUE) > 0 &&
      is.finite(stats::sd(d_geo[upper.tri(d_geo)], na.rm = TRUE)) && stats::sd(d_geo[upper.tri(d_geo)], na.rm = TRUE) > 0
    nw_vg <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L) else 1L
    m_res <- if (has_var && requireNamespace("vegan", quietly = TRUE)) tryCatch(vegan::mantel(stats::as.dist(d_g_sub), stats::as.dist(d_geo), permutations = 999L, parallel = nw_vg), error = function(e) NULL) else NULL
    congr_rows[[length(congr_rows) + 1L]] <- data.frame(
      comparison = "Genomic_vs_Geographic_IBD", n_samples = ng,
      mantel_r = round(if (!is.null(m_res)) as.numeric(m_res$statistic) else if (has_var) stats::cor(d_g_sub[upper.tri(d_g_sub)], d_geo[upper.tri(d_geo)]) else NA_real_, 4),
      mantel_p = round(if (!is.null(m_res)) as.numeric(m_res$signif) else NA_real_, 4),
      procrustes_r = NA_real_, procrustes_p = NA_real_, stringsAsFactors = FALSE
    )
  }
  morph_linked_rows <- NULL
  if (!is.null(db_path) && nzchar(db_path) && file.exists(db_path)) {
    conn <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO), error = function(e) NULL)
    if (!is.null(conn)) {
      on.exit(tryCatch(DBI::dbDisconnect(conn), error = function(e) NULL), add = TRUE)
      cats <- unique(stats::na.omit(pass_df$catalogNumber))
      if (length(cats) >= 3L && db_table_exists_safe(conn, "otu_morphometrics")) {
        q_cats <- paste(vapply(cats, function(s) as.character(DBI::dbQuoteString(conn, s)), character(1)), collapse = ", ")
        morph_linked_rows <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT catalogNumber, AVG(PC1) AS PC1, AVG(PC2) AS PC2, AVG(PC3) AS PC3 FROM otu_morphometrics WHERE run_id = %s AND catalogNumber IN (%s) GROUP BY catalogNumber;", as.character(DBI::dbQuoteString(conn, run_id)), q_cats)), error = function(e) NULL)
      }
    }
  }
  if (!is.null(morph_linked_rows) && nrow(morph_linked_rows) >= 4L) {
    m_match <- match(pass_df$catalogNumber, morph_linked_rows$catalogNumber)
    ok_m <- which(!is.na(m_match) & is.finite(morph_linked_rows$PC1[m_match]) & is.finite(morph_linked_rows$PC2[m_match]))
    if (length(ok_m) >= 4L) {
      X_m <- as.matrix(morph_linked_rows[m_match[ok_m], c("PC1", "PC2", "PC3"), drop = FALSE]); X_m[!is.finite(X_m)] <- 0.0
      X_g <- pco_coords[ok_m, seq_len(min(3L, ncol(pco_coords))), drop = FALSE]; d_m <- as.matrix(stats::dist(X_m)); d_g <- d_geno_mat[ok_m, ok_m, drop = FALSE]
      has_mvar <- is.finite(stats::sd(d_g[upper.tri(d_g)], na.rm = TRUE)) && stats::sd(d_g[upper.tri(d_g)], na.rm = TRUE) > 0 && is.finite(stats::sd(d_m[upper.tri(d_m)], na.rm = TRUE)) && stats::sd(d_m[upper.tri(d_m)], na.rm = TRUE) > 0
      nw_vg2 <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L) else 1L
      m_gm <- if (has_mvar && requireNamespace("vegan", quietly = TRUE)) tryCatch(vegan::mantel(stats::as.dist(d_g), stats::as.dist(d_m), permutations = 999L, parallel = nw_vg2), error = function(e) NULL) else NULL
      p_gm <- if (has_mvar && requireNamespace("vegan", quietly = TRUE)) tryCatch(vegan::protest(X_g, X_m, permutations = 999L, parallel = nw_vg2), error = function(e) NULL) else NULL
      congr_rows[[length(congr_rows) + 1L]] <- data.frame(
        comparison = "Genomic_vs_Morphometrics_PCA", n_samples = length(ok_m),
        mantel_r = if (!is.null(m_gm)) round(as.numeric(m_gm$statistic), 4) else if (has_mvar) round(stats::cor(d_g[upper.tri(d_g)], d_m[upper.tri(d_m)]), 4) else NA_real_,
        mantel_p = if (!is.null(m_gm)) round(as.numeric(m_gm$signif), 4) else NA_real_,
        procrustes_r = if (!is.null(p_gm)) round(as.numeric(p_gm$t0), 4) else NA_real_,
        procrustes_p = if (!is.null(p_gm)) round(as.numeric(p_gm$signif), 4) else NA_real_, stringsAsFactors = FALSE
      )
    }
  }
  if (length(congr_rows) == 0L) {
    return(data.frame(comparison = "Genomic_Only_No_Linked_Geo_Or_Morph_Pairs", n_samples = nrow(pass_df),
                      mantel_r = NA_real_, mantel_p = NA_real_, procrustes_r = NA_real_, procrustes_p = NA_real_, stringsAsFactors = FALSE))
  }
  dplyr::bind_rows(congr_rows)
}

reconcile_step06_cytonuclear_discordance <- function(full_audit, assign_df = NULL, placement_df = NULL) {
  if (is.null(full_audit) || !is.data.frame(full_audit) || nrow(full_audit) == 0L) return(full_audit)
  n <- nrow(full_audit)
  canon_tx <- function(x) {
    s <- tolower(trimws(gsub("_", " ", ifelse(is.na(x), "", as.character(x)))))
    s <- gsub("\u00d7", "x", s)
    s <- gsub("\\bpackera[[:space:]]+x[[:space:]]*memmingeri\\b", "packera xmemmingeri", s)
    s <- ifelse(s %in% c("packera tomentosa", "senecio tomentosus", "p. tomentosa"), "packera dubia", s)
    s <- ifelse(s == "packera millefolia", "packera millefolium", ifelse(s == "packera memmingeri", "packera xmemmingeri", s))
    s
  }
  sids <- as.character(full_audit$genetic_sample_id)
  if (!is.null(assign_df) && is.data.frame(assign_df) && nrow(assign_df) > 0L && "genetic_sample_id" %in% names(assign_df)) {
    ma <- match(sids, as.character(assign_df$genetic_sample_id)); ok_a <- !is.na(ma)
    for (col in intersect(c("geno_predicted_taxon", "geno_posterior_prob", "is_geno_admixed"), names(assign_df))) {
      if (!(col %in% names(full_audit))) full_audit[[col]] <- NA
      miss <- ok_a & is.na(full_audit[[col]])
      if (any(miss)) full_audit[[col]][miss] <- assign_df[[col]][ma[miss]]
    }
  }
  if (!("cytonuclear_discordant" %in% names(full_audit))) full_audit$cytonuclear_discordant <- 0L
  disc_v <- suppressWarnings(as.integer(full_audit$cytonuclear_discordant))
  g_tx <- if ("geno_predicted_taxon" %in% names(full_audit)) as.character(full_audit$geno_predicted_taxon) else rep(NA_character_, n)
  g_pr <- if ("geno_posterior_prob" %in% names(full_audit)) suppressWarnings(as.numeric(full_audit$geno_posterior_prob)) else rep(NA_real_, n)
  g_ad <- if ("is_geno_admixed" %in% names(full_audit)) suppressWarnings(as.integer(full_audit$is_geno_admixed)) else rep(0L, n)
  is_emp_ok <- !is.na(g_tx) & grepl("^Packera[[:space:]_]", trimws(g_tx)) & !is.na(g_pr) & g_pr >= 0.70 & (is.na(g_ad) | g_ad == 0L)
  sp_raw <- if ("species_standardized" %in% names(full_audit)) as.character(full_audit$species_standardized) else rep(NA_character_, n)
  eff_k <- canon_tx(ifelse(is_emp_ok, g_tx, sp_raw))

  pc_vec <- if ("plastid_clade" %in% names(full_audit)) trimws(as.character(full_audit$plastid_clade)) else rep("", n)
  pc_vec[is.na(pc_vec)] <- ""
  pub_cp <- if ("public_plastome_match" %in% names(full_audit)) as.character(full_audit$public_plastome_match) else rep("", n)
  pub_its <- if ("public_its_match" %in% names(full_audit)) as.character(full_audit$public_its_match) else rep("", n)
  pub_wx <- if ("waxy_lineage_match" %in% names(full_audit)) as.character(full_audit$waxy_lineage_match) else rep("", n)

  has_pl_df <- !is.null(placement_df) && is.data.frame(placement_df) && nrow(placement_df) > 0L && "locus" %in% names(placement_df)
  pl_cp <- if (has_pl_df) placement_df[placement_df$locus == "plastid", , drop = FALSE] else NULL

  for (i in which(is_emp_ok)) {
    ek <- eff_k[i]; if (!nzchar(ek)) next
    clade_conflict <- FALSE
    if (nzchar(pc_vec[i])) {
      c_idx <- which(pc_vec == pc_vec[i] & grepl("^packera ", eff_k))
      if (length(c_idx) >= 2L) {
        tb <- sort(table(eff_k[c_idx]), decreasing = TRUE)
        maj_k <- names(tb)[1L]; maj_cnt <- as.integer(tb[1L]); my_cnt <- if (ek %in% names(tb)) as.integer(tb[[ek]]) else 0L
        if (maj_k != ek && my_cnt < maj_cnt) clade_conflict <- TRUE
      }
    }
    cp_tx <- if (!is.na(pub_cp[i])) trimws(pub_cp[i]) else ""
    hq_cp <- nzchar(cp_tx) && grepl("^Packera[[:space:]_]", cp_tx)
    if (has_pl_df && !is.null(pl_cp)) {
      m_cp <- match(sids[i], as.character(pl_cp$genetic_sample_id))
      if (!is.na(m_cp)) {
        cont <- suppressWarnings(as.numeric(pl_cp$kmer_containment[m_cp]))
        n11 <- suppressWarnings(as.integer(pl_cp$n_shared_11mers[m_cp]))
        hq_cp <- is.finite(cont) && cont >= 0.80 && is.finite(n11) && n11 >= 80L
        cp_tx <- if (hq_cp) as.character(pl_cp$best_match_taxon[m_cp]) else ""
      }
    }
    cp_k <- canon_tx(cp_tx)
    cp_valid <- hq_cp && grepl("^packera ", cp_k) && !(cp_k %in% c("packera sp", "packera sp.", "packera unknown"))
    nuc_matches_cp <- cp_valid && (canon_tx(pub_its[i]) == cp_k || canon_tx(pub_wx[i]) == cp_k)
    cp_conflict <- cp_valid && (cp_k != ek) && !nuc_matches_cp
    disc_v[i] <- if (clade_conflict || cp_conflict) 1L else 0L
  }
  full_audit$cytonuclear_discordant <- disc_v
  full_audit
}

