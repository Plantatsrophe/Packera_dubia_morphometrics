# ==============================================================================
# Module: scripts/analysis/step06_reticulation_coancestry.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Modular Phase 2 Phylogenomics Helpers for Step 06:
#     1. `122 bp` RAD Microhaplotype Co-ancestry (`RADpainter` auto-detect + native R NN)
#     2. Continuous-dosage (h, Q12) Triangle Hybrid Classifier (guarding divergent
#        non-hybrid lineages from false F1 calls via AICc + BH-LRT)
#     3. Sources `step06_fbranch_quartets.R` for tree-guided f-branch and SNPs2CF/SNaQ
# ==============================================================================

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
}

if (!exists("compute_tree_guided_fbranch", mode = "function")) {
  for (.fq_p in c("scripts/analysis/step06_fbranch_quartets.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "step06_fbranch_quartets.R"))) {
    if (!is.na(.fq_p) && nzchar(.fq_p) && file.exists(.fq_p)) {
      source(.fq_p)
      break
    }
  }
}

.find_122bp_phylip <- function(raw_genomics_dir) {
  if (is.null(raw_genomics_dir) || !dir.exists(raw_genomics_dir)) return("")
  cand <- list.files(raw_genomics_dir, pattern = "\\.(phylip|phy|nex|nexus)$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  cand <- cand[!grepl("^(vcfref_|ref_packera|public_packera_)", basename(cand), ignore.case = TRUE)]
  if (length(cand) > 0L) cand[1L] else ""
}

.read_phylip_blocks_fast <- function(phy_path, block_len = 122L) {
  if (!nzchar(phy_path) || !file.exists(phy_path)) return(NULL)
  lines <- readLines(phy_path, warn = FALSE); lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 5L) return(NULL)
  n_taxa <- suppressWarnings(as.integer(strsplit(trimws(lines[1L]), "[[:space:]]+")[[1L]][1L]))
  if (!is.finite(n_taxa) || n_taxa < 4L || length(lines) < (n_taxa + 1L)) return(NULL)
  seq_lines <- lines[2:(n_taxa + 1L)]; sids <- character(n_taxa); seqs <- character(n_taxa)
  for (i in seq_len(n_taxa)) {
    sp <- strsplit(trimws(seq_lines[i]), "[[:space:]]+")[[1L]]
    if (length(sp) < 2L) return(NULL)
    sids[i] <- gsub("_sorted$", "", sp[1L]); seqs[i] <- toupper(paste(sp[-1L], collapse = ""))
  }
  n_sites <- min(nchar(seqs))
  if (n_sites < 40L) return(NULL)
  eff_blk <- if (n_sites >= block_len * 4L) as.integer(block_len) else max(20L, as.integer(floor(n_sites / 10L)))
  n_blk <- as.integer(floor(n_sites / eff_blk))
  sel_b <- unique(round(seq(1L, n_blk, length.out = min(400L, n_blk))))
  list(sample_ids = sids, seqs = seqs, block_len = eff_blk, blocks = sel_b)
}

.accum_coanc_from_dist <- function(C_mat, d_mat, pos_map, tol) {
  for (a in seq_len(nrow(d_mat))) {
    fin_b <- which(is.finite(d_mat[a, ]))
    if (length(fin_b) > 0L) {
      min_d <- min(d_mat[a, fin_b])
      donors <- fin_b[d_mat[a, fin_b] <= (min_d + tol)]
      C_mat[pos_map[a], pos_map[donors]] <- C_mat[pos_map[a], pos_map[donors]] + (1.0 / length(donors))
    }
  }
  C_mat
}

compute_microhaplotype_coancestry <- function(pass_ids, dos_mat, raw_genomics_dir) {
  n <- length(pass_ids)
  C_mat <- matrix(0.0, nrow = n, ncol = n, dimnames = list(pass_ids, pass_ids))

  phy <- .read_phylip_blocks_fast(.find_122bp_phylip(raw_genomics_dir))
  used_mode <- "dosage_window_haplotypes"
  if (!is.null(phy) && length(intersect(pass_ids, phy$sample_ids)) >= 4L) {
    used_mode <- "122bp_phylip_microhaplotypes"
    idx_map <- match(pass_ids, phy$sample_ids); valid_pos <- which(!is.na(idx_map)); m <- length(valid_pos)
    raw_seqs <- phy$seqs[idx_map[valid_pos]]
    raw_all <- do.call(rbind, lapply(raw_seqs, function(sq) as.integer(charToRaw(sq))))
    tol_phy <- 1.0 / phy$block_len + 1e-9
    for (b in phy$blocks) {
      st <- (b - 1L) * phy$block_len + 1L; en <- b * phy$block_len
      if (length(unique(substr(raw_seqs, st, en))) < 2L) next
      R_blk <- raw_all[, st:en, drop = FALSE]
      V <- (R_blk != 0x4eL) & (R_blk != 0x2dL)
      n_sh <- tcrossprod(matrix(as.numeric(V), nrow = m, ncol = ncol(R_blk)))
      n_same <- matrix(0.0, nrow = m, ncol = m)
      for (uc in unique(as.integer(R_blk[V]))) {
        Iu <- matrix(as.numeric(R_blk == uc), nrow = m, ncol = ncol(R_blk))
        n_same <- n_same + tcrossprod(Iu)
      }
      d_blk <- (n_sh - n_same) / pmax(1.0, n_sh)
      d_blk[n_sh < 10.0] <- Inf; diag(d_blk) <- Inf
      C_mat <- .accum_coanc_from_dist(C_mat, d_blk, valid_pos, tol_phy)
    }
  }

  if (all(rowSums(C_mat) == 0)) {
    n_loci <- ncol(dos_mat); step <- max(5L, min(20L, as.integer(floor(n_loci / 25L))))
    for (st in seq(1L, max(1L, n_loci - step + 1L), by = step)) {
      sub_m <- dos_mat[, st:min(n_loci, st + step - 1L), drop = FALSE]
      ok_m <- is.finite(sub_m); n_sh <- tcrossprod(matrix(as.numeric(ok_m), nrow = n))
      fin_vals <- sub_m[ok_m]
      if (length(fin_vals) > 0L && all(fin_vals %in% c(0.0, 0.5, 1.0))) {
        I0 <- matrix(as.numeric(ok_m & (sub_m == 0.0)), nrow = n)
        Ih <- matrix(as.numeric(ok_m & (sub_m == 0.5)), nrow = n)
        I1 <- matrix(as.numeric(ok_m & (sub_m == 1.0)), nrow = n)
        abs_sum <- 0.5 * (tcrossprod(I0, Ih) + tcrossprod(Ih, I0) + tcrossprod(Ih, I1) + tcrossprod(I1, Ih)) +
          (tcrossprod(I0, I1) + tcrossprod(I1, I0))
        d_win <- abs_sum / pmax(1.0, n_sh)
      } else {
        M0 <- ifelse(ok_m, sub_m, 0.0); tM0 <- t(M0); tok <- t(ok_m)
        d_win <- matrix(Inf, nrow = n, ncol = n)
        for (i in seq_len(n)) {
          ok_i <- tok & ok_m[i, ]
          d_win[i, ] <- colSums(abs(tM0 - M0[i, ]) * ok_i) / pmax(1L, colSums(ok_i))
        }
      }
      d_win[n_sh < 2.0] <- Inf; diag(d_win) <- Inf
      C_mat <- .accum_coanc_from_dist(C_mat, d_win, seq_len(n), 0.05)
    }
  }

  rs <- rowSums(C_mat); C_norm <- C_mat
  for (i in seq_len(n)) if (rs[i] > 0) C_norm[i, ] <- C_mat[i, ] / rs[i]
  C_sym <- 0.5 * (C_norm + t(C_norm)); max_c <- max(C_sym, na.rm = TRUE)
  d_coanc <- if (is.finite(max_c) && max_c > 0) {
    d_mat_c <- 1.0 - C_sym / max_c; d_mat_c[d_mat_c < 0.0] <- 0.0; diag(d_mat_c) <- 0.0; stats::as.dist(d_mat_c)
  } else stats::dist(dos_mat)
  cl_vec <- sprintf("Coanc_Cluster_%d", as.integer(stats::cutree(stats::hclust(d_coanc, method = "ward.D2"), k = max(2L, min(6L, n - 1L)))))
  df_out <- cbind(genetic_sample_id = pass_ids, coancestry_cluster = cl_vec, as.data.frame(round(C_norm, 5), stringsAsFactors = FALSE), stringsAsFactors = FALSE)
  list(coanc_df = df_out, cluster = cl_vec, mode = used_mode)
}

.eval_hybrid_multinomial_ll <- function(n11, n12, n22, h, q12, eps = 0.015) {
  p_raw <- pmax(0.0, c(1.0 - h - 0.5 * q12, q12, h - 0.5 * q12))
  p_raw <- if (sum(p_raw) > 0.0) p_raw / sum(p_raw) else rep(1.0 / 3.0, 3L)
  p_err <- (1.0 - 3.0 * eps) * p_raw + eps; p_err <- p_err / sum(p_err)
  lgamma(n11 + n12 + n22 + 1.0) - lgamma(n11 + 1.0) - lgamma(n12 + 1.0) - lgamma(n22 + 1.0) +
    n11 * log(p_err[1L]) + n12 * log(p_err[2L]) + n22 * log(p_err[3L])
}

.fit_advanced_introgressant_mle <- function(n11, n12, n22, h_val, q12_val, eps = 0.015) {
  best_ll <- -Inf; best_h <- 0.5; best_q <- 0.5
  for (hc in seq(0.02, 0.98, by = 0.04)) {
    q_max <- 2.0 * min(hc, 1.0 - hc)
    for (qc in (if (q_max >= 0.04) seq(0.0, q_max, by = 0.04) else 0.0)) {
      ll <- .eval_hybrid_multinomial_ll(n11, n12, n22, hc, qc, eps = eps)
      if (ll > best_ll) { best_ll <- ll; best_h <- hc; best_q <- qc }
    }
  }
  n_tot <- max(1.0, n11 + n12 + n22); h_cnt <- (n22 + 0.5 * n12) / n_tot; q_cnt <- n12 / n_tot
  h_cands <- sort(unique(unlist(lapply(unique(pmin(0.99, pmax(0.01, c(best_h, h_val, h_cnt)))), function(ch) seq(max(0.01, ch - 0.05), min(0.99, ch + 0.05), by = 0.01)))))
  for (hf in h_cands) {
    qf_max <- 2.0 * min(hf, 1.0 - hf)
    q_cands <- sort(unique(c(0.0, qf_max, unlist(lapply(unique(pmin(qf_max, pmax(0.0, c(best_q, q12_val, q_cnt)))), function(cq) seq(max(0.0, cq - 0.05), min(qf_max, cq + 0.05), by = 0.01))))))
    for (qf in q_cands) {
      ll <- .eval_hybrid_multinomial_ll(n11, n12, n22, hf, qf, eps = eps)
      if (ll > best_ll) { best_ll <- ll; best_h <- hf; best_q <- qf }
    }
  }
  list(ll = best_ll, h = best_h, q12 = best_q)
}

compute_hybrid_index_q12 <- function(pass_df, dos_mat, q_mat = NULL, p1_vec = NULL, p2_vec = NULL, d_mat = NULL, assign_df = NULL) {
  if (is.matrix(pass_df) && is.data.frame(dos_mat)) { tmp_m <- pass_df; pass_df <- dos_mat; dos_mat <- tmp_m }
  if (is.null(assign_df) && is.data.frame(q_mat)) { assign_df <- q_mat; q_mat <- NULL }
  n <- if (is.data.frame(pass_df)) nrow(pass_df) else 0L
  sids <- if (n > 0L && "genetic_sample_id" %in% names(pass_df)) as.character(pass_df$genetic_sample_id) else character(0)
  species_vec <- if (n > 0L && "species_standardized" %in% names(pass_df)) as.character(pass_df$species_standardized) else rep(NA_character_, n)
  if (n > 0L && !is.null(rownames(dos_mat)) && length(sids) == n && all(sids %in% rownames(dos_mat)) && !identical(rownames(dos_mat), sids)) {
    dos_mat <- dos_mat[match(sids, rownames(dos_mat)), , drop = FALSE]
  }

  empty_fallback <- function(cls_val, p1_val = NA_character_, p2_val = NA_character_) {
    tri_df <- data.frame(
      genetic_sample_id = sids, species_standardized = species_vec, parent_1 = rep(p1_val, n), parent_2 = rep(p2_val, n),
      parent_1_taxon = rep(p1_val, n), parent_2_taxon = rep(p2_val, n), n_diagnostic_loci = integer(n),
      hybrid_index_h = rep(NA_real_, n), interspecific_ancestry_q12 = rep(NA_real_, n), hybrid_generation_class = rep(cls_val, n),
      best_model_aicc = rep(NA_real_, n), akaike_weight = rep(NA_real_, n), delta_aicc_vs_pure = rep(NA_real_, n), lrt_q_bh = rep(NA_real_, n), stringsAsFactors = FALSE
    )
    list(triangle_df = tri_df, hybrid_index_h = tri_df$hybrid_index_h, q12 = tri_df$interspecific_ancestry_q12, hybrid_class = tri_df$hybrid_generation_class)
  }
  if (n == 0L || is.null(dos_mat) || ncol(dos_mat) == 0L) return(empty_fallback(character(0), character(0), character(0)))

  is_query <- if ("is_putative_novel_or_hybrid_query" %in% names(pass_df)) !is.na(pass_df$is_putative_novel_or_hybrid_query) & as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L else rep(FALSE, n)
  valid_sp <- unique(species_vec[!is.na(species_vec) & nzchar(species_vec)])
  q_has_taxa <- !is.null(q_mat) && !is.null(colnames(q_mat)) && sum(colnames(q_mat) %in% valid_sp) >= 2L
  ref_taxa <- if (q_has_taxa) unique(colnames(q_mat)[colnames(q_mat) %in% valid_sp]) else {
    tx_ref <- sort(unique(species_vec[!is_query & !is.na(species_vec) & nzchar(species_vec)])); if (length(tx_ref) < 2L) sort(valid_sp) else tx_ref
  }
  if (length(ref_taxa) < 2L) return(empty_fallback("Undetermined", if (length(ref_taxa) == 1L) ref_taxa[1L] else NA_character_, NA_character_))

  q_cols <- if (q_has_taxa) colnames(q_mat) else ref_taxa
  taxa_freq <- stats::setNames(lapply(ref_taxa, function(tx) {
    idx <- which(!is_query & species_vec == tx); if (length(idx) == 0L) idx <- which(species_vec == tx)
    if (length(idx) == 0L) rep(NA_real_, ncol(dos_mat)) else colMeans(dos_mat[idx, , drop = FALSE], na.rm = TRUE)
  }), ref_taxa)

  a_row <- rep(NA_integer_, n)
  if (!is.null(assign_df) && is.data.frame(assign_df) && nrow(assign_df) > 0L) {
    a_row <- if ("genetic_sample_id" %in% names(assign_df)) match(sids, as.character(assign_df$genetic_sample_id)) else if (nrow(assign_df) == n) seq_len(n) else a_row
  }
  get_a_col <- function(cands, default_vec) {
    if (is.null(assign_df) || !is.data.frame(assign_df)) return(default_vec)
    hit <- intersect(cands, names(assign_df)); if (length(hit) == 0L) default_vec else assign_df[[hit[1L]]][a_row]
  }
  if (is.null(p1_vec)) p1_vec <- as.character(get_a_col(c("geno_parent_1", "parent_1", "parent_1_taxon"), rep(NA_character_, n)))
  if (is.null(p2_vec)) p2_vec <- as.character(get_a_col(c("geno_parent_2", "parent_2", "parent_2_taxon"), rep(NA_character_, n)))
  geno_sp_vec <- as.character(get_a_col(c("geno_predicted_taxon", "predicted_taxon"), rep(NA_character_, n)))
  geno_prob_vec <- suppressWarnings(as.numeric(get_a_col(c("geno_posterior_prob", "posterior_prob"), rep(NA_real_, n))))
  geno_adm_vec <- vapply(get_a_col(c("is_geno_admixed", "geno_admixed"), rep(FALSE, n)), function(v) {
    length(v) > 0L && !is.na(v) && (isTRUE(v) || identical(suppressWarnings(as.integer(v)), 1L) || tolower(as.character(v)) %in% c("true", "t", "1", "yes"))
  }, logical(1))

  med_ref_d <- 0.10; mad_ref_d <- 0.05; div_dist_cut <- 0.225
  if (!is.null(d_mat)) {
    ref_min_dists <- vapply(which(!is_query), function(i) {
      mates <- setdiff(which(!is_query), i); if (length(mates) == 0L) NA_real_ else min(d_mat[i, mates], na.rm = TRUE)
    }, numeric(1))
    ref_fin <- ref_min_dists[is.finite(ref_min_dists)]
    if (length(ref_fin) > 0L) med_ref_d <- stats::median(ref_fin)
    if (length(ref_fin) > 2L) mad_ref_d <- max(0.02, stats::mad(ref_fin))
    div_dist_cut <- med_ref_d + 2.5 * mad_ref_d
  }

  model_names <- c("Pure_Parent_1", "Pure_Parent_2", "F1_Hybrid", "F2_Hybrid", "BC1_Parent_1", "BC1_Parent_2", "Advanced_Introgressant")
  k_vec <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
  h_out <- rep(NA_real_, n); q12_out <- rep(NA_real_, n); cls_out <- rep("Undetermined", n); n_diag_out <- integer(n)
  parent_a <- character(n); parent_b <- character(n); best_aicc_out <- rep(NA_real_, n); best_weight_out <- rep(NA_real_, n)
  delta_aicc_out <- rep(NA_real_, n); raw_p_out <- rep(NA_real_, n); cand_best_idx <- integer(n); pure_best_idx <- integer(n); is_div_vec <- logical(n)

  for (i in seq_len(n)) {
    y <- dos_mat[i, ]
    if (q_has_taxa) {
      ordered_ref <- intersect(q_cols[order(q_mat[i, ], decreasing = TRUE)], ref_taxa)
      cand_taxa <- c(ordered_ref, setdiff(ref_taxa, ordered_ref))
    } else {
      d_ref <- vapply(ref_taxa, function(tx) { tf <- taxa_freq[[tx]]; ok_t <- is.finite(y) & is.finite(tf); if (any(ok_t)) mean(abs(y[ok_t] - tf[ok_t])) else Inf }, numeric(1))
      cand_taxa <- ref_taxa[order(d_ref, decreasing = FALSE)]
    }
    pa <- if (!is.na(p1_vec[i]) && nzchar(p1_vec[i]) && p1_vec[i] %in% ref_taxa) p1_vec[i] else if (!is_query[i] && !is.na(species_vec[i]) && species_vec[i] %in% ref_taxa) species_vec[i] else cand_taxa[1L]
    pb_pool <- setdiff(cand_taxa, pa)
    pb <- if (!is.na(p2_vec[i]) && nzchar(p2_vec[i]) && p2_vec[i] %in% ref_taxa && p2_vec[i] != pa) p2_vec[i] else if (length(pb_pool) >= 1L) pb_pool[1L] else pa
    parent_a[i] <- pa; parent_b[i] <- pb; fa <- taxa_freq[[pa]]; fb <- taxa_freq[[pb]]
    if (is.null(fa) || is.null(fb) || identical(pa, pb)) { cls_out[i] <- "Undetermined"; next }
    ok <- is.finite(y) & is.finite(fa) & is.finite(fb); delta_p <- abs(fb - fa); diag_m <- ok & (delta_p >= 0.45); n_diag <- sum(diag_m); n_diag_out[i] <- n_diag
    if (n_diag == 0L) { cls_out[i] <- "Undetermined"; next }

    conserved_m <- ok & (delta_p <= 0.15)
    priv_div <- if (sum(conserved_m) >= 10L) mean(abs(y[conserved_m] - 0.5 * (fa[conserved_m] + fb[conserved_m]))) else 0.0
    nearest_ref_d <- if (!is.null(d_mat) && any(!is_query)) min(d_mat[i, setdiff(which(!is_query), i)], na.rm = TRUE) else NA_real_
    x_pol <- ifelse(fb[diag_m] >= fa[diag_m], y[diag_m], 1.0 - y[diag_m]); h_val <- mean(x_pol); q12_val <- mean(1.0 - abs(2.0 * x_pol - 1.0))
    h_out[i] <- round(h_val, 4); q12_out[i] <- round(q12_val, 4)

    is_div_assign <- !is.na(geno_sp_vec[i]) && nzchar(geno_sp_vec[i]) && !(geno_sp_vec[i] %in% c(pa, pb)) && is.finite(geno_prob_vec[i]) && geno_prob_vec[i] >= 0.70 && !isTRUE(geno_adm_vec[i]) && q12_val < 0.20
    is_div_dist <- !is.null(d_mat) && ((is_query[i] && is.finite(nearest_ref_d) && nearest_ref_d > div_dist_cut && q12_val < 0.72) || (priv_div >= 0.14 && q12_val < 0.65) || (n_diag < 8L && is.finite(nearest_ref_d) && nearest_ref_d > med_ref_d + 1.5 * mad_ref_d))
    is_div_vec[i] <- is_div_assign || is_div_dist

    if (n_diag < 6L) {
      cls_out[i] <- if (is_div_vec[i]) "Divergent_NonHybrid_Lineage" else if (h_val <= 0.10 && q12_val < 0.20) "Pure_Parent_1" else if (h_val >= 0.90 && q12_val < 0.20) "Pure_Parent_2" else if (h_val >= 0.35 && h_val <= 0.65 && q12_val >= 0.70) "F1_Hybrid" else if (h_val >= 0.30 && h_val <= 0.70 && q12_val >= 0.35) "F2_Hybrid" else if (h_val < 0.35 && q12_val >= 0.25) "BC1_Parent_1" else if (h_val > 0.65 && q12_val >= 0.25) "BC1_Parent_2" else if (q12_val < 0.25 && h_val > 0.10 && h_val < 0.90) "Advanced_Introgressant" else "Pure_Parent_1"
      next
    }
    n11 <- sum(x_pol < 0.25); n12 <- sum(x_pol >= 0.25 & x_pol <= 0.75); n22 <- sum(x_pol > 0.75)
    adv_fit <- .fit_advanced_introgressant_mle(n11, n12, n22, h_val, q12_val, eps = 0.015)
    ll_vec <- c(.eval_hybrid_multinomial_ll(n11, n12, n22, 0.0, 0.0), .eval_hybrid_multinomial_ll(n11, n12, n22, 1.0, 0.0), .eval_hybrid_multinomial_ll(n11, n12, n22, 0.5, 1.0), .eval_hybrid_multinomial_ll(n11, n12, n22, 0.5, 0.5), .eval_hybrid_multinomial_ll(n11, n12, n22, 0.25, 0.5), .eval_hybrid_multinomial_ll(n11, n12, n22, 0.75, 0.5), adv_fit$ll)
    aicc_vec <- vapply(seq_along(ll_vec), function(m) -2.0 * ll_vec[m] + 2.0 * k_vec[m] + (2.0 * k_vec[m] * (k_vec[m] + 1.0)) / max(1.0, n_diag - k_vec[m] - 1.0), numeric(1))
    w_vec <- exp(-0.5 * (aicc_vec - min(aicc_vec))); w_vec <- w_vec / sum(w_vec); best_idx <- which.min(aicc_vec)
    best_aicc_out[i] <- aicc_vec[best_idx]; best_weight_out[i] <- w_vec[best_idx]; delta_aicc_out[i] <- min(aicc_vec[1:2]) - min(aicc_vec[3:7])
    cand_best_idx[i] <- best_idx; pure_best_idx[i] <- if (aicc_vec[1L] <= aicc_vec[2L]) 1L else 2L
    raw_p_out[i] <- stats::pchisq(2.0 * max(0.0, max(ll_vec[3:7]) - max(ll_vec[1:2])), df = 2.0, lower.tail = FALSE)
  }

  lrt_q_bh <- rep(NA_real_, n); fin_p <- which(is.finite(raw_p_out))
  if (length(fin_p) > 0L) lrt_q_bh[fin_p] <- stats::p.adjust(raw_p_out[fin_p], method = "BH")
  for (i in seq_len(n)) {
    if (cand_best_idx[i] == 0L) next
    cls_out[i] <- if (is_div_vec[i]) "Divergent_NonHybrid_Lineage" else if (cand_best_idx[i] >= 3L && is.finite(delta_aicc_out[i]) && delta_aicc_out[i] >= 2.0 && is.finite(lrt_q_bh[i]) && lrt_q_bh[i] < 0.05) model_names[cand_best_idx[i]] else model_names[pure_best_idx[i]]
  }

  tri_df <- data.frame(
    genetic_sample_id = sids, species_standardized = species_vec, parent_1 = parent_a, parent_2 = parent_b,
    parent_1_taxon = parent_a, parent_2_taxon = parent_b, n_diagnostic_loci = n_diag_out, hybrid_index_h = h_out,
    interspecific_ancestry_q12 = q12_out, hybrid_generation_class = cls_out, best_model_aicc = round(best_aicc_out, 4),
    akaike_weight = round(best_weight_out, 4), delta_aicc_vs_pure = round(delta_aicc_out, 4), lrt_q_bh = lrt_q_bh, stringsAsFactors = FALSE
  )
  list(triangle_df = tri_df, hybrid_index_h = h_out, q12 = q12_out, hybrid_class = cls_out)
}
