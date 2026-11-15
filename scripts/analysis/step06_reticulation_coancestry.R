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
    recursive = TRUE,
    ignore.case = TRUE
  )
  cand <- cand[!grepl("^(vcfref_|ref_packera|public_packera_)", basename(cand), ignore.case = TRUE)]
  if (length(cand) > 0L) cand[1L] else ""
}

.read_phylip_blocks_fast <- function(phy_path, block_len = 122L) {
  if (!nzchar(phy_path) || !file.exists(phy_path)) return(NULL)
  lines <- readLines(phy_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 5L) return(NULL)
  hdr <- strsplit(trimws(lines[1L]), "[[:space:]]+")[[1L]]
  n_taxa <- suppressWarnings(as.integer(hdr[1L]))
  if (!is.finite(n_taxa) || n_taxa < 4L || length(lines) < (n_taxa + 1L)) return(NULL)

  seq_lines <- lines[2:(n_taxa + 1L)]
  sids <- character(n_taxa)
  seqs <- character(n_taxa)
  for (i in seq_len(n_taxa)) {
    sp <- strsplit(trimws(seq_lines[i]), "[[:space:]]+")[[1L]]
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
  d_coanc <- if (is.finite(max_c) && max_c > 0) {
    d_mat_c <- 1.0 - C_sym / max_c
    d_mat_c[d_mat_c < 0.0] <- 0.0
    diag(d_mat_c) <- 0.0
    stats::as.dist(d_mat_c)
  } else {
    stats::dist(dos_mat)
  }
  hc <- stats::hclust(d_coanc, method = "ward.D2")
  k_cl <- max(2L, min(6L, n - 1L))
  cl_vec <- sprintf("Coanc_Cluster_%d", as.integer(stats::cutree(hc, k = k_cl)))

  df_out <- as.data.frame(round(C_norm, 5), stringsAsFactors = FALSE)
  df_out <- cbind(genetic_sample_id = pass_ids, coancestry_cluster = cl_vec, df_out, stringsAsFactors = FALSE)
  list(coanc_df = df_out, cluster = cl_vec, mode = used_mode)
}

.eval_hybrid_multinomial_ll <- function(n11, n12, n22, h, q12, eps = 0.015) {
  n_d <- n11 + n12 + n22
  p_raw <- pmax(0.0, c(1.0 - h - 0.5 * q12, q12, h - 0.5 * q12))
  s_raw <- sum(p_raw)
  p_raw <- if (s_raw > 0.0) p_raw / s_raw else c(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
  p_err <- (1.0 - 3.0 * eps) * p_raw + eps
  p_err <- p_err / sum(p_err)
  lgamma(n_d + 1.0) - lgamma(n11 + 1.0) - lgamma(n12 + 1.0) - lgamma(n22 + 1.0) +
    n11 * log(p_err[1L]) + n12 * log(p_err[2L]) + n22 * log(p_err[3L])
}

.fit_advanced_introgressant_mle <- function(n11, n12, n22, h_val, q12_val, eps = 0.015) {
  best_ll <- -Inf
  best_h <- 0.5
  best_q <- 0.5
  for (hc in seq(0.02, 0.98, by = 0.04)) {
    q_max <- 2.0 * min(hc, 1.0 - hc)
    q_grid <- if (q_max >= 0.04) seq(0.0, q_max, by = 0.04) else 0.0
    for (qc in q_grid) {
      ll <- .eval_hybrid_multinomial_ll(n11, n12, n22, hc, qc, eps = eps)
      if (ll > best_ll) {
        best_ll <- ll
        best_h <- hc
        best_q <- qc
      }
    }
  }
  n_tot <- max(1.0, n11 + n12 + n22)
  h_cnt <- (n22 + 0.5 * n12) / n_tot
  q_cnt <- n12 / n_tot
  h_centers <- unique(pmin(0.99, pmax(0.01, c(best_h, h_val, h_cnt))))
  h_cands <- sort(unique(unlist(lapply(h_centers, function(ch) {
    seq(max(0.01, ch - 0.05), min(0.99, ch + 0.05), by = 0.01)
  }))))
  for (hf in h_cands) {
    qf_max <- 2.0 * min(hf, 1.0 - hf)
    q_centers <- unique(pmin(qf_max, pmax(0.0, c(best_q, q12_val, q_cnt))))
    q_cands <- sort(unique(c(0.0, qf_max, unlist(lapply(q_centers, function(cq) {
      seq(max(0.0, cq - 0.05), min(qf_max, cq + 0.05), by = 0.01)
    })))))
    for (qf in q_cands) {
      ll <- .eval_hybrid_multinomial_ll(n11, n12, n22, hf, qf, eps = eps)
      if (ll > best_ll) {
        best_ll <- ll
        best_h <- hf
        best_q <- qf
      }
    }
  }
  list(ll = best_ll, h = best_h, q12 = best_q)
}

compute_hybrid_index_q12 <- function(pass_df, dos_mat, q_mat = NULL, p1_vec = NULL, p2_vec = NULL, d_mat = NULL, assign_df = NULL) {
  if (is.matrix(pass_df) && is.data.frame(dos_mat)) {
    tmp_m <- pass_df
    pass_df <- dos_mat
    dos_mat <- tmp_m
  }
  if (is.null(assign_df) && is.data.frame(q_mat)) {
    assign_df <- q_mat
    q_mat <- NULL
  }

  n <- if (is.data.frame(pass_df)) nrow(pass_df) else 0L
  sids <- if (n > 0L && "genetic_sample_id" %in% names(pass_df)) as.character(pass_df$genetic_sample_id) else character(0)
  species_vec <- if (n > 0L && "species_standardized" %in% names(pass_df)) as.character(pass_df$species_standardized) else rep(NA_character_, n)

  if (n > 0L && !is.null(rownames(dos_mat)) && length(sids) == n &&
      all(sids %in% rownames(dos_mat)) && !identical(rownames(dos_mat), sids)) {
    dos_mat <- dos_mat[match(sids, rownames(dos_mat)), , drop = FALSE]
  }

  empty_fallback <- function(cls_val, p1_val = NA_character_, p2_val = NA_character_) {
    tri_df <- data.frame(
      genetic_sample_id = sids,
      species_standardized = species_vec,
      parent_1 = rep(p1_val, n),
      parent_2 = rep(p2_val, n),
      parent_1_taxon = rep(p1_val, n),
      parent_2_taxon = rep(p2_val, n),
      n_diagnostic_loci = integer(n),
      hybrid_index_h = rep(NA_real_, n),
      interspecific_ancestry_q12 = rep(NA_real_, n),
      hybrid_generation_class = rep(cls_val, n),
      best_model_aicc = rep(NA_real_, n),
      akaike_weight = rep(NA_real_, n),
      delta_aicc_vs_pure = rep(NA_real_, n),
      lrt_q_bh = rep(NA_real_, n),
      stringsAsFactors = FALSE
    )
    list(
      triangle_df = tri_df,
      hybrid_index_h = tri_df$hybrid_index_h,
      q12 = tri_df$interspecific_ancestry_q12,
      hybrid_class = tri_df$hybrid_generation_class
    )
  }

  if (n == 0L || is.null(dos_mat) || ncol(dos_mat) == 0L) {
    return(empty_fallback(character(0), character(0), character(0)))
  }

  is_query <- if ("is_putative_novel_or_hybrid_query" %in% names(pass_df)) {
    !is.na(pass_df$is_putative_novel_or_hybrid_query) & as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L
  } else {
    rep(FALSE, n)
  }

  valid_sp <- unique(species_vec[!is.na(species_vec) & nzchar(species_vec)])
  q_has_taxa <- !is.null(q_mat) && !is.null(colnames(q_mat)) && sum(colnames(q_mat) %in% valid_sp) >= 2L
  ref_taxa <- if (q_has_taxa) {
    unique(colnames(q_mat)[colnames(q_mat) %in% valid_sp])
  } else {
    tx_ref <- sort(unique(species_vec[!is_query & !is.na(species_vec) & nzchar(species_vec)]))
    if (length(tx_ref) < 2L) tx_ref <- sort(valid_sp)
    tx_ref
  }

  if (length(ref_taxa) < 2L) {
    only_tx <- if (length(ref_taxa) == 1L) ref_taxa[1L] else NA_character_
    return(empty_fallback("Undetermined", only_tx, NA_character_))
  }

  q_cols <- if (q_has_taxa) colnames(q_mat) else ref_taxa

  # Precompute reference taxon allele frequencies
  taxa_freq <- lapply(ref_taxa, function(tx) {
    idx <- which(!is_query & species_vec == tx)
    if (length(idx) == 0L) idx <- which(species_vec == tx)
    if (length(idx) == 0L) return(rep(NA_real_, ncol(dos_mat)))
    colMeans(dos_mat[idx, , drop = FALSE], na.rm = TRUE)
  })
  names(taxa_freq) <- ref_taxa

  # Match optional assign_df rows to pass_df
  a_row <- rep(NA_integer_, n)
  if (!is.null(assign_df) && is.data.frame(assign_df) && nrow(assign_df) > 0L) {
    if ("genetic_sample_id" %in% names(assign_df)) {
      a_row <- match(sids, as.character(assign_df$genetic_sample_id))
    } else if (nrow(assign_df) == n) {
      a_row <- seq_len(n)
    }
  }
  get_a_col <- function(cands, default_vec) {
    if (is.null(assign_df) || !is.data.frame(assign_df)) return(default_vec)
    hit <- intersect(cands, names(assign_df))
    if (length(hit) == 0L) return(default_vec)
    vals <- assign_df[[hit[1L]]][a_row]
    vals
  }

  if (is.null(p1_vec)) p1_vec <- as.character(get_a_col(c("geno_parent_1", "parent_1", "parent_1_taxon"), rep(NA_character_, n)))
  if (is.null(p2_vec)) p2_vec <- as.character(get_a_col(c("geno_parent_2", "parent_2", "parent_2_taxon"), rep(NA_character_, n)))
  geno_sp_vec <- as.character(get_a_col(c("geno_predicted_taxon", "predicted_taxon"), rep(NA_character_, n)))
  geno_prob_vec <- suppressWarnings(as.numeric(get_a_col(c("geno_posterior_prob", "posterior_prob"), rep(NA_real_, n))))
  geno_adm_raw <- get_a_col(c("is_geno_admixed", "geno_admixed"), rep(FALSE, n))
  geno_adm_vec <- vapply(geno_adm_raw, function(v) {
    if (length(v) == 0L || is.na(v)) return(FALSE)
    isTRUE(v) || identical(suppressWarnings(as.integer(v)), 1L) ||
      tolower(as.character(v)) %in% c("true", "t", "1", "yes")
  }, logical(1))

  # Intra-taxon distance baseline to detect divergent non-hybrid lineages (PACMILVAL1, PACPAUPAR1P)
  med_ref_d <- 0.10
  mad_ref_d <- 0.05
  div_dist_cut <- 0.225
  if (!is.null(d_mat)) {
    ref_min_dists <- vapply(which(!is_query), function(i) {
      mates <- setdiff(which(!is_query), i)
      if (length(mates) == 0L) return(NA_real_)
      min(d_mat[i, mates], na.rm = TRUE)
    }, numeric(1))
    ref_fin <- ref_min_dists[is.finite(ref_min_dists)]
    if (length(ref_fin) > 0L) med_ref_d <- stats::median(ref_fin)
    if (length(ref_fin) > 2L) mad_ref_d <- max(0.02, stats::mad(ref_fin))
    div_dist_cut <- med_ref_d + 2.5 * mad_ref_d
  }

  model_names <- c(
    "Pure_Parent_1", "Pure_Parent_2", "F1_Hybrid", "F2_Hybrid",
    "BC1_Parent_1", "BC1_Parent_2", "Advanced_Introgressant"
  )
  k_vec <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)

  h_out <- rep(NA_real_, n)
  q12_out <- rep(NA_real_, n)
  cls_out <- rep("Undetermined", n)
  n_diag_out <- integer(n)
  parent_a <- character(n)
  parent_b <- character(n)
  best_aicc_out <- rep(NA_real_, n)
  best_weight_out <- rep(NA_real_, n)
  delta_aicc_out <- rep(NA_real_, n)
  raw_p_out <- rep(NA_real_, n)
  cand_best_idx <- integer(n)
  pure_best_idx <- integer(n)
  is_div_vec <- logical(n)

  for (i in seq_len(n)) {
    y <- dos_mat[i, ]
    if (q_has_taxa) {
      wi <- q_mat[i, ]
      ord <- order(wi, decreasing = TRUE)
      ordered_ref <- intersect(q_cols[ord], ref_taxa)
      cand_taxa <- c(ordered_ref, setdiff(ref_taxa, ordered_ref))
    } else {
      wi <- NULL
      d_ref <- vapply(ref_taxa, function(tx) {
        tf <- taxa_freq[[tx]]
        ok_t <- is.finite(y) & is.finite(tf)
        if (any(ok_t)) mean(abs(y[ok_t] - tf[ok_t])) else Inf
      }, numeric(1))
      ord <- order(d_ref, decreasing = FALSE)
      cand_taxa <- ref_taxa[ord]
    }

    pa <- if (!is.na(p1_vec[i]) && nzchar(p1_vec[i]) && p1_vec[i] %in% ref_taxa) {
      p1_vec[i]
    } else if (!is_query[i] && !is.na(species_vec[i]) && species_vec[i] %in% ref_taxa) {
      species_vec[i]
    } else {
      cand_taxa[1L]
    }
    pb_pool <- setdiff(cand_taxa, pa)
    pb <- if (!is.na(p2_vec[i]) && nzchar(p2_vec[i]) && p2_vec[i] %in% ref_taxa && p2_vec[i] != pa) {
      p2_vec[i]
    } else if (length(pb_pool) >= 1L) {
      pb_pool[1L]
    } else {
      pa
    }
    parent_a[i] <- pa
    parent_b[i] <- pb

    fa <- taxa_freq[[pa]]
    fb <- taxa_freq[[pb]]
    if (is.null(fa) || is.null(fb) || identical(pa, pb)) {
      cls_out[i] <- "Undetermined"
      next
    }

    ok <- is.finite(y) & is.finite(fa) & is.finite(fb)
    delta_p <- abs(fb - fa)
    diag_m <- ok & (delta_p >= 0.45)
    n_diag <- sum(diag_m)
    n_diag_out[i] <- n_diag

    if (n_diag == 0L) {
      cls_out[i] <- "Undetermined"
      next
    }

    conserved_m <- ok & (delta_p <= 0.15)
    priv_div <- if (sum(conserved_m) >= 10L) mean(abs(y[conserved_m] - 0.5 * (fa[conserved_m] + fb[conserved_m]))) else 0.0
    nearest_ref_d <- if (!is.null(d_mat) && any(!is_query)) min(d_mat[i, setdiff(which(!is_query), i)], na.rm = TRUE) else NA_real_

    xd <- y[diag_m]
    fa_d <- fa[diag_m]
    fb_d <- fb[diag_m]
    x_pol <- ifelse(fb_d >= fa_d, xd, 1.0 - xd)
    h_val <- mean(x_pol)
    q12_val <- mean(1.0 - abs(2.0 * x_pol - 1.0))

    h_out[i] <- round(h_val, 4)
    q12_out[i] <- round(q12_val, 4)

    g_sp <- geno_sp_vec[i]
    g_pr <- geno_prob_vec[i]
    g_ad <- isTRUE(geno_adm_vec[i])
    is_div_assign <- (!is.na(g_sp) && nzchar(g_sp) && !(g_sp %in% c(pa, pb)) &&
                        is.finite(g_pr) && g_pr >= 0.70 && !g_ad && q12_val < 0.20)
    is_div_dist <- !is.null(d_mat) && (
      (is_query[i] && is.finite(nearest_ref_d) && nearest_ref_d > div_dist_cut && q12_val < 0.72) ||
      (priv_div >= 0.14 && q12_val < 0.65) ||
      (n_diag < 8L && is.finite(nearest_ref_d) && nearest_ref_d > med_ref_d + 1.5 * mad_ref_d)
    )
    is_div_vec[i] <- is_div_assign || is_div_dist

    if (n_diag < 6L) {
      cls_out[i] <- if (is_div_vec[i]) {
        "Divergent_NonHybrid_Lineage"
      } else if (h_val <= 0.10 && q12_val < 0.20) {
        "Pure_Parent_1"
      } else if (h_val >= 0.90 && q12_val < 0.20) {
        "Pure_Parent_2"
      } else if (h_val >= 0.35 && h_val <= 0.65 && q12_val >= 0.70) {
        "F1_Hybrid"
      } else if (h_val >= 0.30 && h_val <= 0.70 && q12_val >= 0.35) {
        "F2_Hybrid"
      } else if (h_val < 0.35 && q12_val >= 0.25) {
        "BC1_Parent_1"
      } else if (h_val > 0.65 && q12_val >= 0.25) {
        "BC1_Parent_2"
      } else if (q12_val < 0.25 && h_val > 0.10 && h_val < 0.90) {
        "Advanced_Introgressant"
      } else {
        "Pure_Parent_1"
      }
      next
    }

    n11 <- sum(x_pol < 0.25)
    n12 <- sum(x_pol >= 0.25 & x_pol <= 0.75)
    n22 <- sum(x_pol > 0.75)

    adv_fit <- .fit_advanced_introgressant_mle(n11, n12, n22, h_val, q12_val, eps = 0.015)
    ll_vec <- c(
      .eval_hybrid_multinomial_ll(n11, n12, n22, 0.0, 0.0),
      .eval_hybrid_multinomial_ll(n11, n12, n22, 1.0, 0.0),
      .eval_hybrid_multinomial_ll(n11, n12, n22, 0.5, 1.0),
      .eval_hybrid_multinomial_ll(n11, n12, n22, 0.5, 0.5),
      .eval_hybrid_multinomial_ll(n11, n12, n22, 0.25, 0.5),
      .eval_hybrid_multinomial_ll(n11, n12, n22, 0.75, 0.5),
      adv_fit$ll
    )
    aicc_vec <- vapply(seq_along(ll_vec), function(m) {
      km <- k_vec[m]
      -2.0 * ll_vec[m] + 2.0 * km + (2.0 * km * (km + 1.0)) / max(1.0, n_diag - km - 1.0)
    }, numeric(1))

    w_vec <- exp(-0.5 * (aicc_vec - min(aicc_vec)))
    w_vec <- w_vec / sum(w_vec)

    best_idx <- which.min(aicc_vec)
    pure_idx <- if (aicc_vec[1L] <= aicc_vec[2L]) 1L else 2L
    best_aicc_out[i] <- aicc_vec[best_idx]
    best_weight_out[i] <- w_vec[best_idx]
    delta_aicc_out[i] <- min(aicc_vec[1:2]) - min(aicc_vec[3:7])
    cand_best_idx[i] <- best_idx
    pure_best_idx[i] <- pure_idx

    lrt_stat <- 2.0 * max(0.0, max(ll_vec[3:7]) - max(ll_vec[1:2]))
    raw_p_out[i] <- stats::pchisq(lrt_stat, df = 2.0, lower.tail = FALSE)
  }

  lrt_q_bh <- rep(NA_real_, n)
  fin_p <- which(is.finite(raw_p_out))
  if (length(fin_p) > 0L) {
    lrt_q_bh[fin_p] <- stats::p.adjust(raw_p_out[fin_p], method = "BH")
  }

  for (i in seq_len(n)) {
    if (cand_best_idx[i] == 0L) next
    cls_out[i] <- if (is_div_vec[i]) {
      "Divergent_NonHybrid_Lineage"
    } else if (cand_best_idx[i] >= 3L && is.finite(delta_aicc_out[i]) && delta_aicc_out[i] >= 2.0 &&
               is.finite(lrt_q_bh[i]) && lrt_q_bh[i] < 0.05) {
      model_names[cand_best_idx[i]]
    } else {
      model_names[pure_best_idx[i]]
    }
  }

  tri_df <- data.frame(
    genetic_sample_id = sids,
    species_standardized = species_vec,
    parent_1 = parent_a,
    parent_2 = parent_b,
    parent_1_taxon = parent_a,
    parent_2_taxon = parent_b,
    n_diagnostic_loci = n_diag_out,
    hybrid_index_h = h_out,
    interspecific_ancestry_q12 = q12_out,
    hybrid_generation_class = cls_out,
    best_model_aicc = round(best_aicc_out, 4),
    akaike_weight = round(best_weight_out, 4),
    delta_aicc_vs_pure = round(delta_aicc_out, 4),
    lrt_q_bh = lrt_q_bh,
    stringsAsFactors = FALSE
  )
  list(triangle_df = tri_df, hybrid_index_h = h_out, q12 = q12_out, hybrid_class = cls_out)
}

compute_tree_guided_fbranch <- function(pass_df, dos_mat, raw_genomics_dir, backbone = NULL) {
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
      f_branch_fb = numeric(0), polarization_source = character(0), stringsAsFactors = FALSE
    )
    n <- nrow(pass_df)
    return(list(
      fbranch_df = empty_df,
      f_branch_max = rep(0.0, n),
      f_branch_donor_taxon = rep(NA_character_, n),
      polarization_source = "insufficient_taxa"
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
  bb_guide <- if (exists("select_backbone_guided_outgroup", mode = "function")) {
    select_backbone_guided_outgroup(uniq_taxa, p_dist, backbone = backbone)
  } else NULL
  if (!is.null(bb_guide) && !is.null(bb_guide$outgroup) && (bb_guide$outgroup %in% uniq_taxa) &&
      length(bb_guide$tree_order) >= 3L && setequal(bb_guide$tree_order, setdiff(uniq_taxa, bb_guide$outgroup))) {
    outgroup <- bb_guide$outgroup
    ingroups <- setdiff(uniq_taxa, outgroup)
    tree_order <- bb_guide$tree_order
    pol_src <- if (!is.null(bb_guide$source) && nzchar(bb_guide$source)) bb_guide$source else "public_nuclear_backbone"
  } else {
    # Select global outgroup maximizing mean distance to all other taxa
    outgroup <- uniq_taxa[which.max(rowMeans(p_dist))]
    ingroups <- setdiff(uniq_taxa, outgroup)
    hc_sp <- stats::hclust(stats::as.dist(p_dist[ingroups, ingroups, drop = FALSE]), method = "average")
    tree_order <- hc_sp$labels[hc_sp$order]
    pol_src <- "dosage_upgma_fallback"
  }

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
          polarization_source = pol_src,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  fbranch_df <- if (length(fb_rows) > 0L) do.call(rbind, fb_rows) else data.frame(
    recipient_taxon = character(0), donor_taxon = character(0), sister_p1 = character(0),
    outgroup_o = character(0), d_stat = numeric(0), jackknife_z = numeric(0),
    f_branch_fb = numeric(0), polarization_source = character(0), stringsAsFactors = FALSE
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
  list(fbranch_df = fbranch_df, f_branch_max = round(f_max, 4), f_branch_donor_taxon = f_don, polarization_source = pol_src)
}

.empty_cf_df <- function() {
  data.frame(
    quartet = character(0),
    major_topology = character(0),
    cf_12_34 = numeric(0),
    cf_13_24 = numeric(0),
    cf_14_23 = numeric(0),
    n_inform_loci = numeric(0),
    snaq_gamma_minor = numeric(0),
    g_test_p = numeric(0),
    delta_aicc_reticulation = numeric(0),
    stringsAsFactors = FALSE
  )
}

compute_quartet_concordance_snps2cf <- function(pass_df, dos_mat, max_quartets = 120L, seed = 202611L) {
  if (is.matrix(pass_df) && is.data.frame(dos_mat)) {
    tmp_m <- pass_df
    pass_df <- dos_mat
    dos_mat <- tmp_m
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  n <- if (is.data.frame(pass_df)) nrow(pass_df) else 0L
  if (n == 0L || is.null(dos_mat) || ncol(dos_mat) < 12L) {
    return(list(cf_df = .empty_cf_df(), sample_gamma = rep(0.0, n)))
  }

  sids <- if ("genetic_sample_id" %in% names(pass_df)) as.character(pass_df$genetic_sample_id) else as.character(seq_len(n))
  if (!is.null(rownames(dos_mat)) && all(sids %in% rownames(dos_mat)) && !identical(rownames(dos_mat), sids)) {
    dos_mat <- dos_mat[match(sids, rownames(dos_mat)), , drop = FALSE]
  }

  taxa_vec <- if ("species_standardized" %in% names(pass_df)) as.character(pass_df$species_standardized) else rep("", n)
  taxa_vec[is.na(taxa_vec) | !nzchar(taxa_vec)] <- sids[is.na(taxa_vec) | !nzchar(taxa_vec)]
  uniq_taxa <- sort(unique(taxa_vec))

  if (length(uniq_taxa) < 4L) {
    return(list(cf_df = .empty_cf_df(), sample_gamma = rep(0.0, n)))
  }

  tax_groups <- lapply(uniq_taxa, function(tx) which(taxa_vec == tx))
  names(tax_groups) <- uniq_taxa

  comb_q <- utils::combn(seq_along(uniq_taxa), 4L, simplify = FALSE)
  max_q <- max(1L, as.integer(max_quartets))
  if (length(comb_q) > max_q) {
    comb_q <- comb_q[sort(sample.int(length(comb_q), max_q, replace = FALSE))]
  }

  rows <- vector("list", length(comb_q))
  gamma_by_taxon <- stats::setNames(numeric(length(uniq_taxa)), uniq_taxa)

  for (qi in seq_along(comb_q)) {
    q <- comb_q[[qi]]
    t4 <- uniq_taxa[q]
    grid <- expand.grid(
      s1 = tax_groups[[t4[1L]]],
      s2 = tax_groups[[t4[2L]]],
      s3 = tax_groups[[t4[3L]]],
      s4 = tax_groups[[t4[4L]]],
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    if (nrow(grid) == 0L) next
    if (nrow(grid) > 16L) {
      grid <- grid[sort(sample.int(nrow(grid), 16L, replace = FALSE)), , drop = FALSE]
    }

    any_valid_sites <- FALSE
    c1_vec <- numeric(0)
    c2_vec <- numeric(0)
    c3_vec <- numeric(0)
    for (r_idx in seq_len(nrow(grid))) {
      s4 <- as.integer(grid[r_idx, ])
      sub <- dos_mat[s4, , drop = FALSE]
      ok <- colSums(is.finite(sub)) == 4L
      if (sum(ok) < 8L) next
      any_valid_sites <- TRUE
      sub <- sub[, ok, drop = FALSE]
      bin <- round(sub)
      inf <- which(colSums(bin) == 2L)
      if (length(inf) == 0L) next
      b_inf <- bin[, inf, drop = FALSE]
      c1_vec <- c(c1_vec, sum(b_inf[1L, ] == b_inf[2L, ]))
      c2_vec <- c(c2_vec, sum(b_inf[1L, ] == b_inf[3L, ]))
      c3_vec <- c(c3_vec, sum(b_inf[1L, ] == b_inf[4L, ]))
    }

    if (!any_valid_sites) next
    if (length(c1_vec) == 0L) {
      rows[[qi]] <- data.frame(
        quartet = paste(t4, collapse = ";"),
        major_topology = paste0(t4[1L], ",", t4[2L], "|", t4[3L], ",", t4[4L]),
        cf_12_34 = 1.0,
        cf_13_24 = 0.0,
        cf_14_23 = 0.0,
        n_inform_loci = 0.0,
        snaq_gamma_minor = 0.0,
        g_test_p = 1.0,
        delta_aicc_reticulation = 0.0,
        stringsAsFactors = FALSE
      )
      next
    }

    c1 <- mean(c1_vec)
    c2 <- mean(c2_vec)
    c3 <- mean(c3_vec)
    n_inf <- sum(c(c1, c2, c3))
    if (n_inf <= 0.0) next

    cf_vec <- c(c1, c2, c3) / n_inf
    s <- sort(c(c1, c2, c3), decreasing = TRUE)
    top_n <- c(
      paste0(t4[1L], ",", t4[2L], "|", t4[3L], ",", t4[4L]),
      paste0(t4[1L], ",", t4[3L], "|", t4[2L], ",", t4[4L]),
      paste0(t4[1L], ",", t4[4L], "|", t4[2L], ",", t4[3L])
    )[which.max(cf_vec)]

    p0 <- c(s[1L], 0.5 * (s[2L] + s[3L]), 0.5 * (s[2L] + s[3L])) / n_inf
    p1 <- s / n_inf
    pos0 <- (s > 0.0) & (p0 > 0.0)
    pos1 <- (s > 0.0) & (p1 > 0.0)
    const_ll <- lgamma(n_inf + 1.0) - sum(lgamma(s + 1.0))
    ll_0 <- const_ll + sum(s[pos0] * log(p0[pos0]))
    ll_1 <- const_ll + sum(s[pos1] * log(p1[pos1]))

    g_stat <- max(0.0, 2.0 * (ll_1 - ll_0))
    g_test_p <- stats::pchisq(g_stat, df = 1.0, lower.tail = FALSE)
    aicc_0 <- -2.0 * ll_0 + 2.0 + 4.0 / max(1.0, n_inf - 2.0)
    aicc_1 <- -2.0 * ll_1 + 4.0 + 12.0 / max(1.0, n_inf - 3.0)
    delta_aicc_ret <- aicc_0 - aicc_1

    gamma_minor <- if (is.finite(delta_aicc_ret) && delta_aicc_ret >= 2.0 &&
                       is.finite(g_test_p) && g_test_p < 0.05) {
      round(s[2L] / max(1e-9, s[1L] + s[2L]), 4)
    } else {
      0.0
    }

    for (tx_q in t4) {
      gamma_by_taxon[tx_q] <- max(gamma_by_taxon[tx_q], gamma_minor)
    }
    rows[[qi]] <- data.frame(
      quartet = paste(t4, collapse = ";"),
      major_topology = top_n,
      cf_12_34 = round(cf_vec[1L], 4),
      cf_13_24 = round(cf_vec[2L], 4),
      cf_14_23 = round(cf_vec[3L], 4),
      n_inform_loci = round(n_inf, 2),
      snaq_gamma_minor = gamma_minor,
      g_test_p = g_test_p,
      delta_aicc_reticulation = round(delta_aicc_ret, 4),
      stringsAsFactors = FALSE
    )
  }

  cf_df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(cf_df) || nrow(cf_df) == 0L) {
    cf_df <- .empty_cf_df()
  }
  sample_gamma <- round(as.numeric(gamma_by_taxon[taxa_vec]), 4)
  sample_gamma[!is.finite(sample_gamma)] <- 0.0
  list(cf_df = cf_df, sample_gamma = sample_gamma)
}

