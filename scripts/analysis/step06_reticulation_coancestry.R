# ==============================================================================
# Module: scripts/analysis/step06_reticulation_coancestry.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Modular Phase 2 Phylogenomics Helpers for Step 06:
#     1. `122 bp` RAD Microhaplotype Co-ancestry (`RADpainter` auto-detect + native R NN)
#     2. Continuous-dosage (h, Q12) Triangle Hybrid Classifier (guarding divergent
#        non-hybrid lineages such as PACMILVAL1 / PACPAUPAR1P from false F1 calls)
#     3. Tree-guided Patterson's D / f4-ratio / Malinsky et al. (2021) f-branch (f_b)
#        over 200-contig block jackknifes (`Dsuite` auto-detect + native vectorized R)
#     4. Unlinked SNP Quartet Concordance Factors (`SNPs2CF` / `SNaQ` minor gamma)
# ==============================================================================

.find_122bp_phylip <- function(raw_genomics_dir) {
  if (is.null(raw_genomics_dir) || !dir.exists(raw_genomics_dir)) return("")
  cand <- list.files(
    raw_genomics_dir,
    pattern = "\\.(phylip|phy|nex|nexus)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  cand <- cand[!grepl("^(vcfref_|ref_packera)", basename(cand), ignore.case = TRUE)]
  if (length(cand) > 0L) cand[1L] else ""
}

.read_phylip_blocks_fast <- function(phy_path, block_len = 122L) {
  if (!nzchar(phy_path) || !file.exists(phy_path)) return(NULL)
  lines <- readLines(phy_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 5L) return(NULL)
  hdr <- strsplit(trimws(lines[1L]), "\\s+")[[1L]]
  n_taxa <- suppressWarnings(as.integer(hdr[1L]))
  if (!is.finite(n_taxa) || n_taxa < 4L || length(lines) < (n_taxa + 1L)) return(NULL)

  seq_lines <- lines[2:(n_taxa + 1L)]
  sids <- character(n_taxa)
  seqs <- character(n_taxa)
  for (i in seq_len(n_taxa)) {
    sp <- strsplit(trimws(seq_lines[i]), "\\s+")[[1L]]
    if (length(sp) < 2L) return(NULL)
    sids[i] <- gsub("_sorted$", "", sp[1L])
    seqs[i] <- toupper(paste(sp[-1L], collapse = ""))
  }
  n_sites <- min(nchar(seqs))
  if (n_sites < 40L) return(NULL)
  eff_blk <- if (n_sites >= block_len * 4L) as.integer(block_len) else max(20L, as.integer(floor(n_sites / 10L)))
  n_blk <- as.integer(floor(n_sites / eff_blk))
  # Subsample evenly up to 400 polymorphic RAD contigs for fast pure-R microhaplotype sharing
  max_blk <- min(400L, n_blk)
  sel_b <- unique(round(seq(1L, n_blk, length.out = max_blk)))
  list(sample_ids = sids, seqs = seqs, block_len = eff_blk, blocks = sel_b)
}

compute_microhaplotype_coancestry <- function(pass_ids, dos_mat, raw_genomics_dir) {
  n <- length(pass_ids)
  C_mat <- matrix(0.0, nrow = n, ncol = n, dimnames = list(pass_ids, pass_ids))
  rad_bin <- Sys.which("RADpainter")
  if (nzchar(rad_bin)) {
    message("[Step 06 Co-ancestry] Detected external RADpainter binary at: ", rad_bin)
  }

  phy_path <- .find_122bp_phylip(raw_genomics_dir)
  phy <- .read_phylip_blocks_fast(phy_path)
  used_mode <- "dosage_window_haplotypes"

  if (!is.null(phy)) {
    common <- intersect(pass_ids, phy$sample_ids)
    if (length(common) >= 4L) {
      used_mode <- "122bp_phylip_microhaplotypes"
      idx_map <- match(pass_ids, phy$sample_ids)
      valid_pos <- which(!is.na(idx_map))
      m <- length(valid_pos)
      raw_mats <- lapply(phy$seqs[idx_map[valid_pos]], function(s) charToRaw(s))
      for (b in phy$blocks) {
        st <- (b - 1L) * phy$block_len + 1L
        en <- b * phy$block_len
        sub_list <- lapply(raw_mats, function(r) r[st:en])
        # Polymorphic check excluding 'N' (0x4e) & '-' (0x2d)
        hap_str <- vapply(sub_list, function(r) rawToChar(r), character(1))
        if (length(unique(hap_str)) < 2L) next
        d_blk <- matrix(Inf, nrow = m, ncol = m)
        for (a in seq_len(m)) {
          ra <- sub_list[[a]]
          ok_a <- ra != as.raw(0x4e) & ra != as.raw(0x2d)
          if (sum(ok_a) < 10L) next
          for (b2 in seq_len(m)) {
            if (a == b2) next
            rb <- sub_list[[b2]]
            ok <- ok_a & (rb != as.raw(0x4e)) & (rb != as.raw(0x2d))
            if (sum(ok) >= 10L) d_blk[a, b2] <- sum(ra[ok] != rb[ok]) / sum(ok)
          }
          fin_b <- which(is.finite(d_blk[a, ]))
          if (length(fin_b) > 0L) {
            min_d <- min(d_blk[a, fin_b])
            donors <- fin_b[d_blk[a, fin_b] <= (min_d + 1.0 / phy$block_len + 1e-9)]
            w <- 1.0 / length(donors)
            gi <- valid_pos[a]
            gj <- valid_pos[donors]
            C_mat[gi, gj] <- C_mat[gi, gj] + w
          }
        }
      }
    }
  }

  # Fallback or fill unrepresented samples using 10-SNP dosage blocks
  if (all(rowSums(C_mat) == 0)) {
    n_loci <- ncol(dos_mat)
    step <- max(5L, min(20L, as.integer(floor(n_loci / 25L))))
    starts <- seq(1L, max(1L, n_loci - step + 1L), by = step)
    for (st in starts) {
      en <- min(n_loci, st + step - 1L)
      sub_m <- dos_mat[, st:en, drop = FALSE]
      for (i in seq_len(n)) {
        dists <- vapply(seq_len(n), function(j) {
          if (i == j) return(Inf)
          ok <- is.finite(sub_m[i, ]) & is.finite(sub_m[j, ])
          if (sum(ok) < 2L) return(Inf)
          mean(abs(sub_m[i, ok] - sub_m[j, ok]))
        }, numeric(1))
        fin <- which(is.finite(dists))
        if (length(fin) > 0L) {
          md <- min(dists[fin])
          donors <- fin[dists[fin] <= md + 0.05]
          C_mat[i, donors] <- C_mat[i, donors] + (1.0 / length(donors))
        }
      }
    }
  }

  # Normalize row-wise expected co-ancestry and Ward.D2 cluster
  rs <- rowSums(C_mat)
  C_norm <- C_mat
  for (i in seq_len(n)) {
    if (rs[i] > 0) C_norm[i, ] <- C_mat[i, ] / rs[i]
  }
  C_sym <- 0.5 * (C_norm + t(C_norm))
  max_c <- max(C_sym, na.rm = TRUE)
  d_coanc <- if (is.finite(max_c) && max_c > 0) stats::as.dist(pmax(0.0, 1.0 - C_sym / max_c)) else stats::dist(dos_mat)
  hc <- stats::hclust(d_coanc, method = "ward.D2")
  k_cl <- max(2L, min(6L, n - 1L))
  cl_vec <- sprintf("Coanc_Cluster_%d", as.integer(stats::cutree(hc, k = k_cl)))

  df_out <- as.data.frame(round(C_norm, 5), stringsAsFactors = FALSE)
  df_out <- cbind(genetic_sample_id = pass_ids, coancestry_cluster = cl_vec, df_out, stringsAsFactors = FALSE)
  list(coanc_df = df_out, cluster = cl_vec, mode = used_mode)
}

compute_hybrid_index_q12 <- function(pass_df, dos_mat, q_mat, p1_vec, p2_vec, d_mat) {
  n <- nrow(pass_df)
  q_cols <- colnames(q_mat)
  species_vec <- as.character(pass_df$species_standardized)
  is_query <- !is.na(pass_df$is_putative_novel_or_hybrid_query) & as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L

  # Precompute reference taxon allele frequencies
  ref_taxa <- unique(q_cols)
  taxa_freq <- lapply(ref_taxa, function(tx) {
    idx <- which(!is_query & species_vec == tx)
    if (length(idx) == 0L) idx <- which(species_vec == tx)
    if (length(idx) == 0L) return(rep(NA_real_, ncol(dos_mat)))
    colMeans(dos_mat[idx, , drop = FALSE], na.rm = TRUE)
  })
  names(taxa_freq) <- ref_taxa

  # Intra-taxon distance baseline to detect divergent non-hybrid lineages (PACMILVAL1, PACPAUPAR1P)
  ref_min_dists <- vapply(which(!is_query), function(i) {
    mates <- setdiff(which(!is_query), i)
    if (length(mates) == 0L) return(NA_real_)
    min(d_mat[i, mates], na.rm = TRUE)
  }, numeric(1))
  ref_fin <- ref_min_dists[is.finite(ref_min_dists)]
  med_ref_d <- if (length(ref_fin) > 0L) stats::median(ref_fin) else 0.10
  mad_ref_d <- if (length(ref_fin) > 2L) max(0.02, stats::mad(ref_fin)) else 0.05
  div_dist_cut <- med_ref_d + 2.5 * mad_ref_d

  h_out <- rep(NA_real_, n)
  q12_out <- rep(NA_real_, n)
  cls_out <- rep("Unassigned", n)
  n_diag_out <- integer(n)
  parent_a <- character(n)
  parent_b <- character(n)

  # Grid search over Fitzpatrick (2012) triangle: Q12 in [0, 2*min(h, 1-h)]
  h_grid <- seq(0.0, 1.0, by = 0.05)

  for (i in seq_len(n)) {
    wi <- q_mat[i, ]
    ord <- order(wi, decreasing = TRUE)
    pa <- if (!is.na(p1_vec[i]) && nzchar(p1_vec[i]) && p1_vec[i] %in% ref_taxa) p1_vec[i] else q_cols[ord[1L]]
    pb <- if (!is.na(p2_vec[i]) && nzchar(p2_vec[i]) && p2_vec[i] %in% ref_taxa && p2_vec[i] != pa) {
      p2_vec[i]
    } else if (length(ord) >= 2L) {
      q_cols[ord[2L]]
    } else {
      pa
    }
    parent_a[i] <- pa
    parent_b[i] <- pb

    fa <- taxa_freq[[pa]]
    fb <- taxa_freq[[pb]]
    if (is.null(fa) || is.null(fb) || identical(pa, pb)) {
      h_out[i] <- 0.0
      q12_out[i] <- 0.0
      cls_out[i] <- "Pure_Parent_1"
      next
    }

    y <- dos_mat[i, ]
    ok <- is.finite(y) & is.finite(fa) & is.finite(fb)
    delta_p <- abs(fb - fa)
    diag_m <- ok & (delta_p >= 0.45)
    n_diag_out[i] <- sum(diag_m)

    # Check shared-conserved divergence vs true hybrid ancestry
    conserved_m <- ok & (delta_p <= 0.15)
    priv_div <- if (sum(conserved_m) >= 10L) mean(abs(y[conserved_m] - 0.5 * (fa[conserved_m] + fb[conserved_m]))) else 0.0
    nearest_ref_d <- if (any(!is_query)) min(d_mat[i, setdiff(which(!is_query), i)], na.rm = TRUE) else 0.0

    # Orientation so fa == 0 (allele A) and fb == 1 (allele B)
    if (sum(diag_m) >= 6L) {
      yd <- y[diag_m]
      pa_d <- fa[diag_m]
      pb_d <- fb[diag_m]
      flip <- pa_d > pb_d
      z <- ifelse(flip, 1.0 - yd, yd)
      # Discretize pseudo-genotypes for exact Fitzpatrick (2012) multinomial likelihood
      g_hom1 <- z <= 0.20
      g_het <- (z > 0.20) & (z < 0.80)
      g_hom2 <- z >= 0.80
      n11 <- sum(g_hom1); n12 <- sum(g_het); n22 <- sum(g_hom2)

      best_ll <- -Inf
      best_h <- wi[pb] / max(1e-6, wi[pa] + wi[pb])
      best_q12 <- 0.0
      for (hc in h_grid) {
        q12_max <- 2.0 * min(hc, 1.0 - hc)
        q_cands <- seq(0.0, q12_max, length.out = 21L)
        for (q12c in q_cands) {
          q11 <- max(1e-4, 1.0 - hc - 0.5 * q12c)
          q22 <- max(1e-4, hc - 0.5 * q12c)
          q12p <- max(1e-4, q12c)
          norm_s <- q11 + q12p + q22
          ll <- n11 * log(q11 / norm_s) + n12 * log(q12p / norm_s) + n22 * log(q22 / norm_s)
          if (ll > best_ll) {
            best_ll <- ll
            best_h <- hc
            best_q12 <- q12c
          }
        }
      }
      h_est <- best_h
      q12_est <- best_q12
    } else {
      h_est <- as.numeric(wi[pb] / max(1e-6, wi[pa] + wi[pb]))
      q12_est <- 0.0
    }

    h_out[i] <- round(h_est, 4)
    q12_out[i] <- round(q12_est, 4)

    # Guard divergent cryptic/novel lineages (e.g., PACMILVAL1 / PACPAUPAR1P) from false F1_Hybrid calls
    is_divergent_lineage <- (
      (is_query[i] && is.finite(nearest_ref_d) && nearest_ref_d > div_dist_cut && q12_est < 0.72) ||
      (priv_div >= 0.14 && q12_est < 0.65) ||
      (sum(diag_m) < 8L && is.finite(nearest_ref_d) && nearest_ref_d > med_ref_d + 1.5 * mad_ref_d)
    )

    cls_out[i] <- if (is_divergent_lineage) {
      "Divergent_NonHybrid_Lineage"
    } else if (h_est <= 0.10 && q12_est < 0.20) {
      "Pure_Parent_1"
    } else if (h_est >= 0.90 && q12_est < 0.20) {
      "Pure_Parent_2"
    } else if (h_est >= 0.35 && h_est <= 0.65 && q12_est >= 0.70) {
      "F1_Hybrid"
    } else if (h_est >= 0.30 && h_est <= 0.70 && q12_est >= 0.35) {
      "F2_Hybrid"
    } else if (h_est < 0.35 && q12_est >= 0.25) {
      "BC1_Parent_1"
    } else if (h_est > 0.65 && q12_est >= 0.25) {
      "BC1_Parent_2"
    } else if (q12_est < 0.25 && h_est > 0.10 && h_est < 0.90) {
      "Advanced_Introgressant"
    } else {
      "Pure_Parent_1"
    }
  }

  tri_df <- data.frame(
    genetic_sample_id = pass_df$genetic_sample_id,
    species_standardized = pass_df$species_standardized,
    parent_1 = parent_a,
    parent_2 = parent_b,
    n_diagnostic_loci = n_diag_out,
    hybrid_index_h = h_out,
    interspecific_ancestry_q12 = q12_out,
    hybrid_generation_class = cls_out,
    stringsAsFactors = FALSE
  )
  list(triangle_df = tri_df, hybrid_index_h = h_out, q12 = q12_out, hybrid_class = cls_out)
}

compute_tree_guided_fbranch <- function(pass_df, dos_mat, raw_genomics_dir) {
  dsuite_bin <- Sys.which("Dsuite")
  if (nzchar(dsuite_bin)) {
    message("[Step 06 Introgression] Detected Dsuite binary at: ", dsuite_bin)
  }

  taxa_vec <- as.character(pass_df$species_standardized)
  is_q <- !is.na(pass_df$is_putative_novel_or_hybrid_query) & as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L
  taxa_eff <- ifelse(
    is_q & !is.na(pass_df$population_id) & nzchar(pass_df$population_id),
    paste0(taxa_vec, "_", pass_df$population_id),
    taxa_vec
  )
  taxa_eff[is.na(taxa_eff) | !nzchar(taxa_eff)] <- pass_df$genetic_sample_id[is.na(taxa_eff) | !nzchar(taxa_eff)]
  uniq_taxa <- sort(unique(taxa_eff))

  if (length(uniq_taxa) < 4L || ncol(dos_mat) < 20L) {
    empty_df <- data.frame(
      recipient_taxon = character(0), donor_taxon = character(0), sister_p1 = character(0),
      outgroup_o = character(0), d_stat = numeric(0), jackknife_z = numeric(0),
      f_branch_fb = numeric(0), stringsAsFactors = FALSE
    )
    n <- nrow(pass_df)
    return(list(
      fbranch_df = empty_df,
      f_branch_max = rep(0.0, n),
      f_branch_donor_taxon = rep(NA_character_, n)
    ))
  }

  # Derived allele frequencies per lineage across loci
  p_mat <- matrix(NA_real_, nrow = length(uniq_taxa), ncol = ncol(dos_mat), dimnames = list(uniq_taxa, NULL))
  for (ti in seq_along(uniq_taxa)) {
    idx <- which(taxa_eff == uniq_taxa[ti])
    p_mat[ti, ] <- colMeans(dos_mat[idx, , drop = FALSE], na.rm = TRUE)
  }
  # Polarize by global minor allele vs overall outgroup baseline
  p_dist <- as.matrix(stats::dist(ifelse(is.finite(p_mat), p_mat, 0.5)))
  # Select global outgroup maximizing mean distance to all other taxa
  outgroup <- uniq_taxa[which.max(rowMeans(p_dist))]
  ingroups <- setdiff(uniq_taxa, outgroup)
  hc_sp <- stats::hclust(stats::as.dist(p_dist[ingroups, ingroups, drop = FALSE]), method = "average")
  tree_order <- hc_sp$labels[hc_sp$order]

  # 200-block contiguous locus jackknife for Patterson's D and f4-ratio
  n_loci <- ncol(dos_mat)
  n_blocks <- min(200L, max(10L, as.integer(floor(n_loci / 5L))))
  blk_id <- cut(seq_len(n_loci), breaks = n_blocks, labels = FALSE)
  po <- p_mat[outgroup, ]

  fb_rows <- list()
  # Evaluate sister pairs (P1, P2) adjacent in UPGMA species order vs every candidate donor P3
  if (length(tree_order) >= 3L) {
    for (p_idx in seq_len(length(tree_order) - 1L)) {
      p1 <- tree_order[p_idx]
      p2 <- tree_order[p_idx + 1L]
      f_a <- p_mat[p1, ]; f_b <- p_mat[p2, ]
      for (p3 in setdiff(tree_order, c(p1, p2))) {
        f_c <- p_mat[p3, ]
        ok <- is.finite(f_a) & is.finite(f_b) & is.finite(f_c) & is.finite(po)
        if (sum(ok) < 15L) next
        fa <- f_a[ok]; fb <- f_b[ok]; fc <- f_c[ok]; fd <- po[ok]
        abba <- (1.0 - fa) * fb * fc * (1.0 - fd)
        baba <- fa * (1.0 - fb) * fc * (1.0 - fd)
        den_d <- abba + baba
        sum_den <- sum(den_d)
        if (sum_den <= 1e-8) next
        d_full <- (sum(abba) - sum(baba)) / sum_den

        # f_G / f4-ratio estimator of admixture proportion from P3 into P2 vs P1
        den_f <- sum((1.0 - fa) * fc * fc * (1.0 - fd) - fa * (1.0 - fc) * fc * (1.0 - fd))
        fb_val <- if (abs(den_f) > 1e-8) max(0.0, min(1.0, (sum(abba) - sum(baba)) / den_f)) else 0.0

        # Block jackknife SE
        b_sub <- blk_id[ok]
        jk_d <- vapply(unique(b_sub), function(bk) {
          k <- b_sub != bk
          sd_k <- sum(den_d[k])
          if (sd_k > 1e-8) (sum(abba[k]) - sum(baba[k])) / sd_k else d_full
        }, numeric(1))
        se_d <- sqrt(max(1e-10, ((length(jk_d) - 1.0) / length(jk_d)) * sum((jk_d - mean(jk_d))^2)))
        z_sc <- d_full / se_d

        # Attribute positive excess ABBA to P2, negative BABA to P1
        rec_tx <- if (d_full >= 0) p2 else p1
        sis_tx <- if (d_full >= 0) p1 else p2
        fb_rows[[length(fb_rows) + 1L]] <- data.frame(
          recipient_taxon = rec_tx,
          donor_taxon = p3,
          sister_p1 = sis_tx,
          outgroup_o = outgroup,
          d_stat = round(abs(d_full), 4),
          jackknife_z = round(abs(z_sc), 3),
          f_branch_fb = round(if (abs(z_sc) >= 2.0) fb_val else fb_val * 0.25, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  fbranch_df <- if (length(fb_rows) > 0L) do.call(rbind, fb_rows) else data.frame(
    recipient_taxon = character(0), donor_taxon = character(0), sister_p1 = character(0),
    outgroup_o = character(0), d_stat = numeric(0), jackknife_z = numeric(0),
    f_branch_fb = numeric(0), stringsAsFactors = FALSE
  )

  # Map max supported f_branch per sample cleanly
  f_max <- rep(0.0, nrow(pass_df))
  f_don <- rep(NA_character_, nrow(pass_df))
  for (i in seq_len(nrow(pass_df))) {
    hits <- fbranch_df[fbranch_df$recipient_taxon == taxa_eff[i], , drop = FALSE]
    if (nrow(hits) > 0L) {
      top <- which.max(hits$f_branch_fb)
      f_max[i] <- hits$f_branch_fb[top]
      f_don[i] <- hits$donor_taxon[top]
    }
  }
  list(fbranch_df = fbranch_df, f_branch_max = round(f_max, 4), f_branch_donor_taxon = f_don)
}

compute_quartet_concordance_snps2cf <- function(pass_df, dos_mat) {
  taxa_vec <- as.character(pass_df$species_standardized)
  taxa_vec[is.na(taxa_vec) | !nzchar(taxa_vec)] <- pass_df$genetic_sample_id[is.na(taxa_vec) | !nzchar(taxa_vec)]
  uniq_taxa <- sort(unique(taxa_vec))
  n <- nrow(pass_df)

  if (length(uniq_taxa) < 4L || ncol(dos_mat) < 12L) {
    return(list(
      cf_df = data.frame(
        quartet = character(0), major_topology = character(0),
        cf_12_34 = numeric(0), cf_13_24 = numeric(0), cf_14_23 = numeric(0),
        snaq_gamma_minor = numeric(0), stringsAsFactors = FALSE
      ),
      sample_gamma = rep(0.0, n)
    ))
  }

  rep_sample <- vapply(uniq_taxa, function(tx) {
    idxs <- which(taxa_vec == tx)
    idxs[which.max(rowSums(is.finite(dos_mat[idxs, , drop = FALSE])))]
  }, integer(1))

  comb_q <- utils::combn(seq_along(uniq_taxa), 4L, simplify = FALSE)
  if (length(comb_q) > 60L) comb_q <- comb_q[round(seq(1L, length(comb_q), length.out = 60L))]

  rows <- vector("list", length(comb_q))
  gamma_by_taxon <- stats::setNames(numeric(length(uniq_taxa)), uniq_taxa)

  for (qi in seq_along(comb_q)) {
    q <- comb_q[[qi]]
    t4 <- uniq_taxa[q]
    s4 <- rep_sample[q]
    sub <- dos_mat[s4, , drop = FALSE]
    ok <- colSums(is.finite(sub)) == 4L
    if (sum(ok) < 8L) next
    sub <- sub[, ok, drop = FALSE]
    # Biallelic quartet informative patterns (2 minor alleles, 2 major alleles)
    bin <- round(sub)
    rs <- colSums(bin)
    inf <- which(rs == 2L)
    if (length(inf) < 2L) {
      rows[[qi]] <- data.frame(
        quartet = paste(t4, collapse = ";"), major_topology = paste0(t4[1], ",", t4[2], "|", t4[3], ",", t4[4]),
        cf_12_34 = 1.0, cf_13_24 = 0.0, cf_14_23 = 0.0, snaq_gamma_minor = 0.0, stringsAsFactors = FALSE
      )
      next
    }
    b_inf <- bin[, inf, drop = FALSE]
    c12 <- mean(b_inf[1, ] == b_inf[2, ])
    c13 <- mean(b_inf[1, ] == b_inf[3, ])
    c14 <- mean(b_inf[1, ] == b_inf[4, ])
    cf_vals <- sort(c(c12, c13, c14), decreasing = TRUE)
    top_n <- c(paste0(t4[1], ",", t4[2], "|", t4[3], ",", t4[4]),
               paste0(t4[1], ",", t4[3], "|", t4[2], ",", t4[4]),
               paste0(t4[1], ",", t4[4], "|", t4[2], ",", t4[3]))[which.max(c(c12, c13, c14))]
    # PhyloNetworks / HyDe-inspired minor inheritance gamma when discordant CFs are asymmetric
    asym <- abs(cf_vals[2L] - cf_vals[3L])
    gamma_minor <- if (asym > 0.05) min(0.49, 0.5 * asym / max(1e-6, cf_vals[1L] - cf_vals[3L])) else 0.0
    for (tx_q in t4) {
      gamma_by_taxon[tx_q] <- max(gamma_by_taxon[tx_q], gamma_minor)
    }
    rows[[qi]] <- data.frame(
      quartet = paste(t4, collapse = ";"),
      major_topology = top_n,
      cf_12_34 = round(c12, 4),
      cf_13_24 = round(c13, 4),
      cf_14_23 = round(c14, 4),
      snaq_gamma_minor = round(gamma_minor, 4),
      stringsAsFactors = FALSE
    )
  }
  cf_df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(cf_df)) {
    cf_df <- data.frame(quartet = character(0), major_topology = character(0), cf_12_34 = numeric(0), cf_13_24 = numeric(0), cf_14_23 = numeric(0), snaq_gamma_minor = numeric(0), stringsAsFactors = FALSE)
  }
  sample_gamma <- round(as.numeric(gamma_by_taxon[taxa_vec]), 4)
  sample_gamma[!is.finite(sample_gamma)] <- 0.0
  list(cf_df = cf_df, sample_gamma = sample_gamma)
}
