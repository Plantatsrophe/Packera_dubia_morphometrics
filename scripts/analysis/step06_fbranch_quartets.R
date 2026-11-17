# ==============================================================================
# Module: scripts/analysis/step06_fbranch_quartets.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Tree-guided Patterson's D / f4-ratio / f-branch (f_b) jackknife introgression
#   and unlinked SNP quartet concordance factor (SNPs2CF / SNaQ gamma) helpers.
# ==============================================================================

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
}

compute_tree_guided_fbranch <- function(pass_df, dos_mat, raw_genomics_dir, backbone = NULL) {
  taxa_vec <- as.character(pass_df$species_standardized)
  is_q <- !is.na(pass_df$is_putative_novel_or_hybrid_query) & as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L
  taxa_eff <- ifelse(is_q & !is.na(pass_df$population_id) & nzchar(pass_df$population_id), paste0(taxa_vec, "_", pass_df$population_id), taxa_vec)
  taxa_eff[is.na(taxa_eff) | !nzchar(taxa_eff)] <- pass_df$genetic_sample_id[is.na(taxa_eff) | !nzchar(taxa_eff)]
  uniq_taxa <- sort(unique(taxa_eff))

  if (length(uniq_taxa) < 4L || ncol(dos_mat) < 20L) {
    empty_df <- data.frame(
      recipient_taxon = character(0), donor_taxon = character(0), sister_p1 = character(0),
      outgroup_o = character(0), d_stat = numeric(0), jackknife_z = numeric(0),
      f_branch_fb = numeric(0), polarization_source = character(0), stringsAsFactors = FALSE
    )
    n <- nrow(pass_df)
    return(list(fbranch_df = empty_df, f_branch_max = rep(0.0, n), f_branch_donor_taxon = rep(NA_character_, n), polarization_source = "insufficient_taxa"))
  }

  p_mat <- matrix(NA_real_, nrow = length(uniq_taxa), ncol = ncol(dos_mat), dimnames = list(uniq_taxa, NULL))
  for (ti in seq_along(uniq_taxa)) p_mat[ti, ] <- colMeans(dos_mat[which(taxa_eff == uniq_taxa[ti]), , drop = FALSE], na.rm = TRUE)
  p_dist <- as.matrix(stats::dist(ifelse(is.finite(p_mat), p_mat, 0.5)))
  bb_guide <- if (exists("select_backbone_guided_outgroup", mode = "function")) select_backbone_guided_outgroup(uniq_taxa, p_dist, backbone = backbone) else NULL
  if (!is.null(bb_guide) && !is.null(bb_guide$outgroup) && (bb_guide$outgroup %in% uniq_taxa) &&
      length(bb_guide$tree_order) >= 3L && setequal(bb_guide$tree_order, setdiff(uniq_taxa, bb_guide$outgroup))) {
    outgroup <- bb_guide$outgroup; ingroups <- setdiff(uniq_taxa, outgroup); tree_order <- bb_guide$tree_order
    pol_src <- if (!is.null(bb_guide$source) && nzchar(bb_guide$source)) bb_guide$source else "public_nuclear_backbone"
  } else {
    outgroup <- uniq_taxa[which.max(rowMeans(p_dist))]; ingroups <- setdiff(uniq_taxa, outgroup)
    hc_sp <- stats::hclust(stats::as.dist(p_dist[ingroups, ingroups, drop = FALSE]), method = "average")
    tree_order <- hc_sp$labels[hc_sp$order]; pol_src <- "dosage_upgma_fallback"
  }

  n_loci <- ncol(dos_mat); n_blocks <- min(200L, max(10L, as.integer(floor(n_loci / 5L))))
  blk_id <- cut(seq_len(n_loci), breaks = n_blocks, labels = FALSE); po <- p_mat[outgroup, ]
  fb_rows <- list()
  if (length(tree_order) >= 3L) {
    for (p_idx in seq_len(length(tree_order) - 1L)) {
      p1 <- tree_order[p_idx]; p2 <- tree_order[p_idx + 1L]; f_a <- p_mat[p1, ]; f_b <- p_mat[p2, ]
      for (p3 in setdiff(tree_order, c(p1, p2))) {
        f_c <- p_mat[p3, ]; ok <- is.finite(f_a) & is.finite(f_b) & is.finite(f_c) & is.finite(po)
        if (sum(ok) < 15L) next
        fa <- f_a[ok]; fb <- f_b[ok]; fc <- f_c[ok]; fd <- po[ok]
        abba <- (1.0 - fa) * fb * fc * (1.0 - fd); baba <- fa * (1.0 - fb) * fc * (1.0 - fd)
        den_d <- abba + baba; sum_den <- sum(den_d)
        if (sum_den <= 1e-8) next
        d_full <- (sum(abba) - sum(baba)) / sum_den
        den_f <- sum((1.0 - fa) * fc * fc * (1.0 - fd) - fa * (1.0 - fc) * fc * (1.0 - fd))
        fb_val <- if (abs(den_f) > 1e-8) max(0.0, min(1.0, (sum(abba) - sum(baba)) / den_f)) else 0.0
        b_sub <- blk_id[ok]
        blk_mat <- rowsum(cbind(abba - baba, den_d), group = b_sub, reorder = FALSE)
        rem_num <- (sum(abba) - sum(baba)) - blk_mat[, 1L]
        rem_den <- sum_den - blk_mat[, 2L]
        jk_d <- ifelse(rem_den > 1e-8, rem_num / rem_den, d_full)
        se_d <- sqrt(max(1e-10, ((length(jk_d) - 1.0) / length(jk_d)) * sum((jk_d - mean(jk_d))^2)))
        z_sc <- d_full / se_d
        rec_tx <- if (d_full >= 0) p2 else p1; sis_tx <- if (d_full >= 0) p1 else p2
        fb_rows[[length(fb_rows) + 1L]] <- data.frame(
          recipient_taxon = rec_tx, donor_taxon = p3, sister_p1 = sis_tx, outgroup_o = outgroup,
          d_stat = round(abs(d_full), 4), jackknife_z = round(abs(z_sc), 3),
          f_branch_fb = round(if (abs(z_sc) >= 2.0) fb_val else fb_val * 0.25, 4),
          polarization_source = pol_src, stringsAsFactors = FALSE
        )
      }
    }
  }

  fbranch_df <- if (length(fb_rows) > 0L) do.call(rbind, fb_rows) else data.frame(
    recipient_taxon = character(0), donor_taxon = character(0), sister_p1 = character(0),
    outgroup_o = character(0), d_stat = numeric(0), jackknife_z = numeric(0),
    f_branch_fb = numeric(0), polarization_source = character(0), stringsAsFactors = FALSE
  )
  f_max <- rep(0.0, nrow(pass_df)); f_don <- rep(NA_character_, nrow(pass_df))
  for (i in seq_len(nrow(pass_df))) {
    hits <- fbranch_df[fbranch_df$recipient_taxon == taxa_eff[i], , drop = FALSE]
    if (nrow(hits) > 0L) {
      top <- which.max(hits$f_branch_fb); f_max[i] <- hits$f_branch_fb[top]; f_don[i] <- hits$donor_taxon[top]
    }
  }
  list(fbranch_df = fbranch_df, f_branch_max = round(f_max, 4), f_branch_donor_taxon = f_don, polarization_source = pol_src)
}

.empty_cf_df <- function() {
  data.frame(
    quartet = character(0), major_topology = character(0), cf_12_34 = numeric(0), cf_13_24 = numeric(0),
    cf_14_23 = numeric(0), n_inform_loci = numeric(0), snaq_gamma_minor = numeric(0),
    g_test_p = numeric(0), delta_aicc_reticulation = numeric(0), stringsAsFactors = FALSE
  )
}

compute_quartet_concordance_snps2cf <- function(pass_df, dos_mat, max_quartets = 120L, seed = 202611L) {
  if (is.matrix(pass_df) && is.data.frame(dos_mat)) {
    tmp_m <- pass_df; pass_df <- dos_mat; dos_mat <- tmp_m
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))

  n <- if (is.data.frame(pass_df)) nrow(pass_df) else 0L
  if (n == 0L || is.null(dos_mat) || ncol(dos_mat) < 12L) return(list(cf_df = .empty_cf_df(), sample_gamma = rep(0.0, n)))
  sids <- if ("genetic_sample_id" %in% names(pass_df)) as.character(pass_df$genetic_sample_id) else as.character(seq_len(n))
  if (!is.null(rownames(dos_mat)) && all(sids %in% rownames(dos_mat)) && !identical(rownames(dos_mat), sids)) {
    dos_mat <- dos_mat[match(sids, rownames(dos_mat)), , drop = FALSE]
  }
  taxa_vec <- if ("species_standardized" %in% names(pass_df)) as.character(pass_df$species_standardized) else rep("", n)
  taxa_vec[is.na(taxa_vec) | !nzchar(taxa_vec)] <- sids[is.na(taxa_vec) | !nzchar(taxa_vec)]
  uniq_taxa <- sort(unique(taxa_vec))
  if (length(uniq_taxa) < 4L) return(list(cf_df = .empty_cf_df(), sample_gamma = rep(0.0, n)))

  tax_groups <- lapply(uniq_taxa, function(tx) which(taxa_vec == tx)); names(tax_groups) <- uniq_taxa
  comb_q <- utils::combn(seq_along(uniq_taxa), 4L, simplify = FALSE); max_q <- max(1L, as.integer(max_quartets))
  if (length(comb_q) > max_q) comb_q <- comb_q[sort(sample.int(length(comb_q), max_q, replace = FALSE))]

  rows <- vector("list", length(comb_q)); gamma_by_taxon <- stats::setNames(numeric(length(uniq_taxa)), uniq_taxa)
  has_mscq <- requireNamespace("MSCquartets", quietly = TRUE)
  for (qi in seq_along(comb_q)) {
    t4 <- uniq_taxa[comb_q[[qi]]]
    grid <- expand.grid(s1 = tax_groups[[t4[1L]]], s2 = tax_groups[[t4[2L]]], s3 = tax_groups[[t4[3L]]], s4 = tax_groups[[t4[4L]]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    if (nrow(grid) == 0L) next
    if (nrow(grid) > 16L) grid <- grid[sort(sample.int(nrow(grid), 16L, replace = FALSE)), , drop = FALSE]

    any_valid_sites <- FALSE; c1_vec <- numeric(0); c2_vec <- numeric(0); c3_vec <- numeric(0)
    for (r_idx in seq_len(nrow(grid))) {
      sub <- dos_mat[as.integer(grid[r_idx, ]), , drop = FALSE]
      ok <- colSums(is.finite(sub)) == 4L
      if (sum(ok) < 8L) next
      any_valid_sites <- TRUE; bin <- round(sub[, ok, drop = FALSE]); inf <- which(colSums(bin) == 2L)
      if (length(inf) == 0L) next
      b_inf <- bin[, inf, drop = FALSE]
      c1_vec <- c(c1_vec, sum(b_inf[1L, ] == b_inf[2L, ]))
      c2_vec <- c(c2_vec, sum(b_inf[1L, ] == b_inf[3L, ]))
      c3_vec <- c(c3_vec, sum(b_inf[1L, ] == b_inf[4L, ]))
    }
    if (!any_valid_sites) next
    if (length(c1_vec) == 0L) {
      rows[[qi]] <- data.frame(quartet = paste(t4, collapse = ";"), major_topology = paste0(t4[1L], ",", t4[2L], "|", t4[3L], ",", t4[4L]), cf_12_34 = 1.0, cf_13_24 = 0.0, cf_14_23 = 0.0, n_inform_loci = 0.0, snaq_gamma_minor = 0.0, g_test_p = 1.0, delta_aicc_reticulation = 0.0, stringsAsFactors = FALSE)
      next
    }
    c1 <- mean(c1_vec); c2 <- mean(c2_vec); c3 <- mean(c3_vec); n_inf <- sum(c(c1, c2, c3))
    if (n_inf <= 0.0) next
    cf_vec <- c(c1, c2, c3) / n_inf; s <- sort(c(c1, c2, c3), decreasing = TRUE)
    top_n <- c(paste0(t4[1L], ",", t4[2L], "|", t4[3L], ",", t4[4L]), paste0(t4[1L], ",", t4[3L], "|", t4[2L], ",", t4[4L]), paste0(t4[1L], ",", t4[4L], "|", t4[2L], ",", t4[3L]))[which.max(cf_vec)]
    p0 <- c(s[1L], 0.5 * (s[2L] + s[3L]), 0.5 * (s[2L] + s[3L])) / n_inf; p1 <- s / n_inf
    pos0 <- (s > 0.0) & (p0 > 0.0); pos1 <- (s > 0.0) & (p1 > 0.0)
    const_ll <- lgamma(n_inf + 1.0) - sum(lgamma(s + 1.0))
    ll_0 <- const_ll + sum(s[pos0] * log(p0[pos0])); ll_1 <- const_ll + sum(s[pos1] * log(p1[pos1]))
    g_test_p <- stats::pchisq(max(0.0, 2.0 * (ll_1 - ll_0)), df = 1.0, lower.tail = FALSE)
    if (has_mscq) tryCatch(MSCquartets::quartetTreeTest(c(c1, c2, c3)), error = function(e) NULL)
    delta_aicc_ret <- (-2.0 * ll_0 + 2.0 + 4.0 / max(1.0, n_inf - 2.0)) - (-2.0 * ll_1 + 4.0 + 12.0 / max(1.0, n_inf - 3.0))
    gamma_minor <- if (is.finite(delta_aicc_ret) && delta_aicc_ret >= 2.0 && is.finite(g_test_p) && g_test_p < 0.05) round(s[2L] / max(1e-9, s[1L] + s[2L]), 4) else 0.0
    for (tx_q in t4) gamma_by_taxon[tx_q] <- max(gamma_by_taxon[tx_q], gamma_minor)
    rows[[qi]] <- data.frame(
      quartet = paste(t4, collapse = ";"), major_topology = top_n,
      cf_12_34 = round(cf_vec[1L], 4), cf_13_24 = round(cf_vec[2L], 4), cf_14_23 = round(cf_vec[3L], 4),
      n_inform_loci = round(n_inf, 2), snaq_gamma_minor = gamma_minor, g_test_p = g_test_p,
      delta_aicc_reticulation = round(delta_aicc_ret, 4), stringsAsFactors = FALSE
    )
  }

  cf_df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(cf_df) || nrow(cf_df) == 0L) cf_df <- .empty_cf_df()
  sample_gamma <- round(as.numeric(gamma_by_taxon[taxa_vec]), 4); sample_gamma[!is.finite(sample_gamma)] <- 0.0
  list(cf_df = cf_df, sample_gamma = sample_gamma)
}

