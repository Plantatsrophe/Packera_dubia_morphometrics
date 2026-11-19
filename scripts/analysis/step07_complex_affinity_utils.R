if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.with_preserved_rng_seed <- function(seed, expr_fn) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({ if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
  if (!is.null(seed)) set.seed(as.integer(seed)); expr_fn()
}

.std_complex_taxon <- function(x) {
  if (exists("standardize_packera_taxon", mode = "function")) { res <- standardize_packera_taxon(x); return(ifelse(is.na(res) | res == "Unknown", as.character(x), res)) }
  if (exists("standardize_taxon_name", mode = "function")) { res <- standardize_taxon_name(x); return(ifelse(is.na(res), as.character(x), res)) }
  as.character(x)
}

.parse_step07_complex_semicolon <- function(df, col_nm) {
  if (!(col_nm %in% names(df))) return(character(0))
  raw <- trimws(as.character(stats::na.omit(df[[col_nm]]))); raw <- raw[nzchar(raw)]; if (length(raw) == 0L) return(character(0))
  toks <- unique(trimws(unlist(strsplit(raw, ";", fixed = TRUE)))); intersect(toks[nzchar(toks) & !(tolower(toks) %in% c("na", "none", "null"))], names(df))
}

.sort_step07_complex_pcs <- function(cols, pfx) if (length(cols) <= 1L) cols else cols[order(suppressWarnings(as.integer(sub(paste0("^", pfx), "", cols))), cols)]

.complex_dna_tier0_flags <- function(df, is_type) {
  n <- nrow(df); f0 <- list(dna_ok = rep(FALSE, n), dna_bad = rep(FALSE, n))
  if (n == 0L || !any(c("geno_predicted_taxon", "geno_linkage_tier", "linkage_tier") %in% names(df))) return(f0)
  if (exists("resolve_tier0_anchor_overrides", mode = "function")) {
    t0 <- resolve_tier0_anchor_overrides(df)
    return(list(dna_ok = t0$is_tier0_dna_anchor & !t0$evict_from_given_anchor & !is_type, dna_bad = t0$evict_from_given_anchor & !is_type))
  }
  bfn <- function(v) if (is.logical(v)) !is.na(v) & v else if (is.numeric(v) || is.integer(v)) !is.na(v) & (v != 0) else tolower(trimws(as.character(v))) %in% c("1", "true", "t", "yes", "y")
  tier <- if ("geno_linkage_tier" %in% names(df)) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% names(df)) as.character(df$linkage_tier) else rep("", n)
  qc <- if ("geno_quality_pass" %in% names(df)) bfn(df$geno_quality_pass) else rep(FALSE, n)
  adm <- if ("is_geno_admixed" %in% names(df)) bfn(df$is_geno_admixed) else rep(FALSE, n)
  nov <- if ("is_putative_novel_or_hybrid_query" %in% names(df)) bfn(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  pred <- if ("geno_predicted_taxon" %in% names(df)) trimws(as.character(df$geno_predicted_taxon)) else rep("", n); pred[is.na(pred)] <- ""
  sp <- if ("species_standardized" %in% names(df)) trimws(as.character(df$species_standardized)) else rep("", n); sp[is.na(sp)] <- ""
  prob <- if ("geno_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$geno_posterior_prob)) else rep(NA_real_, n)
  margin <- if ("geno_margin" %in% names(df)) suppressWarnings(as.numeric(df$geno_margin)) else rep(NA_real_, n)
  q12_col <- intersect("interspecific_ancestry_q12", names(df))[1L]
  q12 <- if (!is.na(q12_col)) suppressWarnings(as.numeric(df[[q12_col]])) else rep(NA_real_, n)
  is_t0 <- !is_type & (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc & !adm & !nov & is.finite(prob) & (prob >= 0.85) & (is.na(margin) | margin >= 0.40) & nzchar(pred) & (pred != "Unknown")
  ev <- !is_type & (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc & ((nzchar(pred) & (pred != "Unknown") & (pred != sp) & is.finite(prob) & (prob >= 0.85)) | (is.finite(q12) & (q12 >= 0.35)) | (adm & !is.finite(q12)))
  list(dna_ok = is_t0 & !ev & !is_type, dna_bad = ev & !is_type)
}

.resolve_complex_specimen_weights <- function(df) {
  n <- nrow(df); if (n == 0L) return(numeric(0))
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else list(is_type_specimen = rep(FALSE, n))
  is_type <- as.logical(type_info$is_type_specimen); is_type[is.na(is_type)] <- FALSE
  if ("is_type_specimen" %in% names(df)) { itc <- as.logical(df$is_type_specimen); itc[is.na(itc)] <- FALSE; is_type <- is_type | itc }
  dna <- .complex_dna_tier0_flags(df, is_type)
  if ("authority_weight" %in% names(df) && any(is.finite(aw <- suppressWarnings(as.numeric(df$authority_weight))) & aw > 0)) {
    aw[!is.finite(aw) | aw <= 0] <- 0.50; wt <- pmin(4.0, pmax(0.15, aw))
  } else {
    tier <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n)
    prom <- if ("is_promoted_bronze" %in% names(df)) (!is.na(df$is_promoted_bronze) & as.logical(df$is_promoted_bronze)) else rep(FALSE, n)
    wt <- ifelse(is_type | grepl("type", tier, ignore.case = TRUE), 4.0, ifelse(grepl("tier_?1|gold", tier, ignore.case = TRUE), 3.0, ifelse(grepl("tier_?2|silver", tier, ignore.case = TRUE), 1.5, ifelse(prom, 1.0, 0.5))))
    noise <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
    wt[noise & !is_type] <- 0.05
  }
  wt[dna$dna_ok] <- 4.0; wt[dna$dna_bad] <- 0.15; wt
}

if (!exists(".is_step07_centroid_coord", mode = "function")) {
  .is_step07_centroid_coord <- function(df) {
    n <- nrow(df); if (n == 0L) return(logical(0))
    if (!("is_centroid_or_inst_coord" %in% names(df))) return(rep(FALSE, n))
    raw_c <- df$is_centroid_or_inst_coord
    c_flag <- if (is.logical(raw_c)) !is.na(raw_c) & raw_c
    else if (is.numeric(raw_c) || is.integer(raw_c)) !is.na(raw_c) & (raw_c != 0)
    else tolower(trimws(as.character(raw_c))) %in% c("1", "true", "t", "yes", "y")
    is_t <- if ("is_type_specimen" %in% names(df)) (!is.na(df$is_type_specimen) & as.logical(df$is_type_specimen)) else rep(FALSE, n)
    if (exists(".identify_type_specimen_vec", mode = "function")) {
      ti <- .identify_type_specimen_vec(df)
      if (!is.null(ti$is_type_specimen) && length(ti$is_type_specimen) == n) {
        it2 <- as.logical(ti$is_type_specimen); it2[is.na(it2)] <- FALSE; is_t <- is_t | it2
      }
    }
    c_flag & !is_t
  }
}

# 3. Dynamic Morphometric & Ecological Feature Discovery
resolve_step07_complex_morph_features <- function(morph_df, anchor_df = morph_df, n_taxa = 2L) {
  max_p <- max(2L, min(15L, nrow(anchor_df) - as.integer(n_taxa) - 1L))
  ok_num <- function(cc) {
    if (!(cc %in% names(anchor_df)) || !(cc %in% names(morph_df)) || is.logical(anchor_df[[cc]]) || is.character(anchor_df[[cc]])) return(FALSE)
    v <- suppressWarnings(as.numeric(anchor_df[[cc]])); v <- v[is.finite(v)]; length(v) >= 3L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
  }
  sel <- .parse_step07_complex_semicolon(morph_df, "selected_morph_features"); sel <- sel[vapply(sel, ok_num, logical(1))]
  if (length(sel) >= 2L) {
    is_pc <- grepl("^PC[0-9]+$", sel); ord_pc <- .sort_step07_complex_pcs(grep("^PC[0-9]+$", sel, value = TRUE), "PC")
    return(utils::head(unique(c(ord_pc, sel[!is_pc])), max_p))
  }
  pc_cols <- .sort_step07_complex_pcs(grep("^PC[0-9]+$", names(anchor_df), value = TRUE), "PC")
  pc_cols <- pc_cols[vapply(pc_cols, ok_num, logical(1))]; can_cols <- intersect(c("can1", "can2", "can3"), names(anchor_df))
  can_cols <- can_cols[vapply(can_cols, ok_num, logical(1))]; cands <- if (length(pc_cols) >= 2L) pc_cols else unique(c(can_cols, pc_cols))
  if (length(cands) <= 2L) return(cands)
  X_c <- as.matrix(anchor_df[, cands, drop = FALSE]); storage.mode(X_c) <- "double"
  ev <- apply(X_c, 2L, stats::var, na.rm = TRUE); ev[!is.finite(ev) | ev < 0] <- 0.0
  m <- length(ev); ev_prop <- ev / max(sum(ev), 1e-12); bs_prop <- rev(cumsum(1.0 / (m:1L))) / m
  pa_q95 <- .with_preserved_rng_seed(42L, function() {
    n_r <- max(4L, nrow(X_c))
    null_mat <- replicate(20L, sort(eigen(stats::cor(matrix(stats::rnorm(n_r * m), nrow = n_r, ncol = m)), symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE) / m)
    apply(null_mat, 1L, stats::quantile, probs = 0.95, na.rm = TRUE, names = FALSE)})
  keep_k <- max(2L, min(max_p, sum(cumprod(as.integer(ev_prop >= bs_prop | ev_prop >= pa_q95)) == 1L)))
  cands[seq_len(min(length(cands), keep_k))]
}

resolve_step07_complex_eco_predictors <- function(df) {
  excl <- c(
    "edaphic_z_score", "climatic_z_score", "type_edaphic_z_score", "edaphic_schoener_d", "schoener_d", "edaphic_niche_overlap_d",
    "sgmc_lithology_class", "Dominant_Lithology_Class", "edaphic_best_fit_taxon",
    "rf_predicted_taxon", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors")
  ok_env <- function(cands) {
    cands <- setdiff(unique(cands), excl)
    cands[vapply(cands, function(cc) {
      if (!(cc %in% names(df)) || is.logical(df[[cc]]) || is.character(df[[cc]])) return(FALSE)
      v <- suppressWarnings(as.numeric(df[[cc]])); v <- v[is.finite(v)]; length(v) >= 4L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
    }, logical(1))]
  }
  meta_cols <- ok_env(c(
    .parse_step07_complex_semicolon(df, "primary_env_predictors"), .parse_step07_complex_semicolon(df, "primary_edaph_predictors"),
    .parse_step07_complex_semicolon(df, "primary_clim_predictors")
  ))
  if (length(meta_cols) > 0L) return(meta_cols)
  cands <- ok_env(grep("^(bio[0-9]+_|fine_|ilr_texture_|ssurgo_|solus_|topo_)", names(df), value = TRUE))
  if (!any(grepl("^fine_", cands))) cands <- unique(c(cands, ok_env(grep("^soil_", names(df), value = TRUE))))
  if (length(cands) <= 1L) return(cands)
  E <- as.matrix(df[, cands, drop = FALSE]); storage.mode(E) <- "double"
  cands <- cands[order(-colSums(is.finite(E)), cands)]; E <- E[, cands, drop = FALSE]; kept <- cands[1L]
  for (cc in cands[-1L]) {
    cors <- vapply(kept, function(kc) {
      r <- suppressWarnings(stats::cor(E[, cc], E[, kc], use = "pairwise.complete.obs"))
      if (is.finite(r)) abs(r) else 0.0
    }, numeric(1))
    if (all(cors < 0.70)) kept <- c(kept, cc)
  }
  kept
}

# 4. Authority-Weighted Morphometric, Ecological, and Circular Phenological Affinity Matrices
compute_ecological_affinity <- function(morph_df, taxa, anchor_mask = NULL) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); n_taxa <- length(taxa); n <- nrow(morph_df)
  S_env <- matrix(NA_real_, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa)); diag(S_env) <- 0.0
  D_env <- matrix(NA_real_, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa)); diag(D_env) <- 0.0
  eco_cols <- resolve_step07_complex_eco_predictors(morph_df)
  if (length(eco_cols) == 0L) {
    return(list(S_env = S_env, D_env = D_env, has_env = FALSE, eco_predictors = character(0), specimen_schoener_d = rep(NA_real_, n)))
  }
  E_raw <- as.matrix(morph_df[, eco_cols, drop = FALSE]); storage.mode(E_raw) <- "double"
  is_cent <- .is_step07_centroid_coord(morph_df)
  has_row_env <- (rowSums(is.finite(E_raw)) > 0L) & !is_cent
  if (sum(has_row_env) < 4L) {
    return(list(S_env = S_env, D_env = D_env, has_env = FALSE, eco_predictors = eco_cols, specimen_schoener_d = rep(NA_real_, n)))
  }
  if (!is.logical(anchor_mask) || length(anchor_mask) != n) anchor_mask <- select_taxon_anchors(morph_df)$anchor_mask
  noise_m <- if ("is_confident_label_noise" %in% names(morph_df)) (!is.na(morph_df$is_confident_label_noise) & as.logical(morph_df$is_confident_label_noise)) else rep(FALSE, n)
  wt_all <- .resolve_complex_specimen_weights(morph_df); E <- E_raw; ref_env_mask <- anchor_mask & has_row_env & !noise_m
  if (sum(ref_env_mask) < 4L) ref_env_mask <- has_row_env & !noise_m
  for (j in seq_len(ncol(E))) {
    vj <- E_raw[ref_env_mask & is.finite(E_raw[, j]), j]
    if (length(vj) == 0L) vj <- E_raw[is.finite(E_raw[, j]), j]
    med_j <- if (length(vj) > 0L) stats::median(vj) else 0.0
    sd_j <- if (length(vj) >= 2L && is.finite(stats::sd(vj)) && stats::sd(vj) > 1e-8) stats::sd(vj) else 1.0
    for (tx in taxa) {
      idx_t <- which(morph_df$species_standardized == tx & has_row_env & !is.finite(E[, j]))
      if (length(idx_t) > 0L) {
        vt <- E_raw[morph_df$species_standardized == tx & ref_env_mask & is.finite(E_raw[, j]), j]
        E[idx_t, j] <- if (length(vt) > 0L) stats::median(vt) else med_j
      }
    }
    E[has_row_env, j] <- (E[has_row_env, j] - med_j) / sd_j
  }
  env_idx_all <- which(has_row_env); Z <- matrix(NA_real_, nrow = n, ncol = 2L)
  if (ncol(E) >= 2L) {
    fit_idx <- which(ref_env_mask); pca_e <- stats::prcomp(E[fit_idx, , drop = FALSE], center = TRUE, scale. = FALSE)
    Z[env_idx_all, ] <- stats::predict(pca_e, newdata = E[env_idx_all, , drop = FALSE])[, 1:2, drop = FALSE]
  } else {
    Z[env_idx_all, 1L] <- E[env_idx_all, 1L]; Z[env_idx_all, 2L] <- 0.0
  }
  z1_rng <- range(Z[env_idx_all, 1L], na.rm = TRUE); z2_rng <- range(Z[env_idx_all, 2L], na.rm = TRUE)
  grid_pts <- as.matrix(expand.grid(
    seq(z1_rng[1L] - 1.0, z1_rng[2L] + 1.0, length.out = 25L), seq(z2_rng[1L] - 1.0, z2_rng[2L] + 1.0, length.out = 25L)
  ))
  has_tx_env <- stats::setNames(logical(n_taxa), taxa); dens_by_tx <- stats::setNames(vector("list", n_taxa), taxa)
  for (tx in taxa) {
    idx <- which(morph_df$species_standardized == tx & anchor_mask & has_row_env & !noise_m)
    if (length(idx) < 2L) idx <- which(morph_df$species_standardized == tx & has_row_env & !noise_m)
    if (length(idx) == 0L) next
    has_tx_env[[tx]] <- TRUE; Zt <- Z[idx, , drop = FALSE]; wt_t <- wt_all[idx]; wt_t <- wt_t / max(sum(wt_t), 1e-12)
    bw <- max(0.35, 1.06 * mean(apply(Zt, 2L, function(v) if (length(v) >= 2L) max(stats::sd(v), 0.25) else 0.5)) * (length(idx)^(-0.2)))
    K_mat <- exp(-0.5 * (outer(grid_pts[, 1L], Zt[, 1L], "-")^2 + outer(grid_pts[, 2L], Zt[, 2L], "-")^2) / (bw^2))
    p <- as.numeric(K_mat %*% wt_t); dens_by_tx[[tx]] <- p / max(sum(p), 1e-12)
  }
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    if (has_tx_env[[taxa[i]]] && has_tx_env[[taxa[j]]]) {
      S_env[i, j] <- max(0.0, min(1.0, 1.0 - 0.5 * sum(abs(dens_by_tx[[i]] - dens_by_tx[[j]])))); D_env[i, j] <- 1.0 - S_env[i, j]
    }
  }
  spec_d <- rep(NA_real_, n)
  for (tx in taxa) {
    idx <- which(morph_df$species_standardized == tx & has_row_env); ref_i <- which(morph_df$species_standardized == tx & ref_env_mask)
    if (length(ref_i) == 0L) ref_i <- idx
    if (length(idx) > 0L) {
      med_tx <- apply(E[ref_i, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
      spec_d[idx] <- exp(-0.5 * rowMeans(sweep(E[idx, , drop = FALSE], 2L, med_tx, "-")^2))
    }
  }
  list(S_env = S_env, D_env = D_env, has_env = sum(has_tx_env) >= 2L, eco_predictors = eco_cols, specimen_schoener_d = spec_d)
}

compute_morphometric_affinity <- function(morph_df, anchor_df, taxa) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); anchor_df <- as.data.frame(anchor_df, stringsAsFactors = FALSE); n_taxa <- length(taxa); feat_cols <- resolve_step07_complex_morph_features(morph_df, anchor_df, n_taxa)
  if (length(feat_cols) == 0L) {
    fb_c <- intersect(c("PC1", "PC2", "PC3", "can1", "can2", "can3"), names(anchor_df))
    feat_cols <- fb_c[vapply(fb_c, function(cc) any(is.finite(suppressWarnings(as.numeric(anchor_df[[cc]])))), logical(1))]
    if (length(feat_cols) == 0L) feat_cols <- fb_c
  }
  for (rc in setdiff(intersect(c("involucre_height_mm", "involucre_width_mm", "AR_inv_std", "W_inv_std"), names(anchor_df)), feat_cols)) {
    vals <- suppressWarnings(as.numeric(anchor_df[[rc]]))
    if (nrow(anchor_df) >= 2L && sum(is.finite(vals) & vals > 0) >= 0.5 * nrow(anchor_df) && isTRUE(stats::sd(vals, na.rm = TRUE) > 1e-6)) {
      col_nm <- paste0("log_", rc); anchor_df[[col_nm]] <- log10(pmax(suppressWarnings(as.numeric(anchor_df[[rc]])), 1e-3))
      morph_df[[col_nm]] <- log10(pmax(suppressWarnings(as.numeric(morph_df[[rc]])), 1e-3)); feat_cols <- c(feat_cols, col_nm)
    }
  }
  X_anchor <- as.matrix(anchor_df[, feat_cols, drop = FALSE]); storage.mode(X_anchor) <- "double"
  X_all <- as.matrix(morph_df[, feat_cols, drop = FALSE]); storage.mode(X_all) <- "double"
  for (j in seq_len(ncol(X_anchor))) {
    med_j <- stats::median(X_anchor[, j], na.rm = TRUE); if (!is.finite(med_j)) med_j <- 0.0
    X_anchor[!is.finite(X_anchor[, j]), j] <- med_j; X_all[!is.finite(X_all[, j]), j] <- med_j
  }
  p <- ncol(X_anchor); centroids <- matrix(0.0, nrow = n_taxa, ncol = p, dimnames = list(taxa, feat_cols))
  anc_type <- if (exists(".identify_type_specimen_vec", mode = "function")) as.logical(.identify_type_specimen_vec(anchor_df)$is_type_specimen) else rep(FALSE, nrow(anchor_df))
  anc_type[is.na(anc_type)] <- FALSE; wt_anc <- .resolve_complex_specimen_weights(anchor_df)
  wt_all <- .resolve_complex_specimen_weights(morph_df); S_pooled <- matrix(0.0, nrow = p, ncol = p); neff_total <- 0.0
  for (tx in taxa) {
    tx_rows <- which(anchor_df$species_standardized == tx)
    if (length(tx_rows) > 0L) {
      sub_x <- X_anchor[tx_rows, , drop = FALSE]; w_tx <- wt_anc[tx_rows]
    } else {
      all_rows <- which(morph_df$species_standardized == tx); sub_x <- X_all[all_rows, , drop = FALSE]; w_tx <- wt_all[all_rows]
    }
    w_norm <- w_tx / max(sum(w_tx), 1e-12); mu_base <- as.numeric(colSums(sub_x * w_norm)); type_rows_tx <- if (exists(".filter_homotypic_type_rows", mode = "function")) .filter_homotypic_type_rows(anchor_df, tx_rows[anc_type[tx_rows]], tx) else tx_rows[anc_type[tx_rows]]
    if (length(type_rows_tx) >= 1L) {
      n_eff_t <- if (exists(".count_effective_type_gatherings", mode = "function")) .count_effective_type_gatherings(anchor_df, type_rows_tx) else length(type_rows_tx)
      gamma_k <- if (exists(".compute_type_shrinkage_gamma", mode = "function")) .compute_type_shrinkage_gamma(n_eff_t, nrow(sub_x)) else min(0.35, (2.5 * n_eff_t) / (max(1.0, nrow(sub_x)) + 2.5 * n_eff_t))
      mu_base <- (1.0 - gamma_k) * mu_base + gamma_k * colMeans(X_anchor[type_rows_tx, , drop = FALSE], na.rm = TRUE)
    }
    centroids[tx, ] <- mu_base
    if (nrow(sub_x) >= 2L) {
      neff_k <- (sum(w_tx)^2) / max(sum(w_tx^2), 1e-12)
      cw <- tryCatch(stats::cov.wt(sub_x, wt = w_norm, method = "unbiased")$cov, error = function(e) stats::cov(sub_x))
      cw[!is.finite(cw)] <- 0.0; df_k <- max(1.0, neff_k - 1.0); S_pooled <- S_pooled + df_k * cw; neff_total <- neff_total + df_k
    }
  }
  S_pooled <- if (neff_total > 0.0) S_pooled / neff_total else stats::cov(X_anchor); S_pooled[!is.finite(S_pooled)] <- 0.0
  lam_lw <- min(0.50, max(0.02, p / max(2.0 * neff_total, 4.0))); v_bar <- max(mean(diag(S_pooled), na.rm = TRUE), 1e-4)
  S_reg <- (1.0 - lam_lw) * S_pooled + lam_lw * diag(v_bar, nrow = p, ncol = p)
  S_inv <- tryCatch(solve(S_reg), error = function(e) MASS::ginv(S_reg))
  D2_mahal <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    diff_ij <- centroids[i, ] - centroids[j, ]; D2_mahal[i, j] <- max(0.0, as.numeric(t(diff_ij) %*% S_inv %*% diff_ij))
  }
  off_diag_d2 <- D2_mahal[upper.tri(D2_mahal)]; med_d2 <- stats::median(off_diag_d2[off_diag_d2 > 0], na.rm = TRUE)
  scale_d2 <- if (is.finite(med_d2) && med_d2 > 1e-8) min(med_d2, 2.0 * stats::qchisq(0.99, df = max(2L, p))) else 1.0
  K_mahal <- exp(-D2_mahal / (2.0 * max(scale_d2, 1.0))); diag(K_mahal) <- 0.0
  noise_m <- if ("is_confident_label_noise" %in% names(morph_df)) (!is.na(morph_df$is_confident_label_noise) & as.logical(morph_df$is_confident_label_noise)) else rep(FALSE, nrow(morph_df))
  sub_call <- if ("taxonomic_treatment_subtype" %in% names(morph_df)) as.character(morph_df$taxonomic_treatment_subtype) else rep("", nrow(morph_df))
  misid_m <- !is.na(sub_call) & (sub_call == "Misidentification_Reassign")
  C_morph <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  has_type_nn <- "nearest_type_morph_taxon" %in% names(morph_df)
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]; idx_i <- which(morph_df$species_standardized == tx_i & !noise_m & !misid_m)
    if (length(idx_i) == 0L) idx_i <- which(morph_df$species_standardized == tx_i)
    if (length(idx_i) == 0L) next
    X_i <- X_all[idx_i, , drop = FALSE]; w_i <- wt_all[idx_i]; w_i <- w_i / max(sum(w_i), 1e-12)
    cda_pred_i <- if ("cda_predicted_taxon" %in% names(morph_df)) as.character(morph_df$cda_predicted_taxon[idx_i]) else rep(tx_i, length(idx_i))
    type_nn_i <- if (has_type_nn) as.character(morph_df$nearest_type_morph_taxon[idx_i]) else rep(NA_character_, length(idx_i))
    log_likes <- matrix(0.0, nrow = length(idx_i), ncol = n_taxa)
    for (j in seq_len(n_taxa)) {
      diff_mat <- sweep(X_i, 2L, centroids[j, ], "-"); log_likes[, j] <- -0.5 * pmax(0.0, rowSums((diff_mat %*% S_inv) * diff_mat))
    }
    probs_i <- exp(sweep(log_likes, 1L, apply(log_likes, 1L, max), "-")); probs_i <- sweep(probs_i, 1L, pmax(rowSums(probs_i), 1e-12), "/")
    for (j in seq_len(n_taxa)) if (i != j) {
      hit_ij <- (!is.na(cda_pred_i) & cda_pred_i == taxa[j]) |
                (!is.na(type_nn_i) & type_nn_i == taxa[j] & probs_i[, j] >= 0.10) |
                (probs_i[, j] >= 0.15)
      hit_ij[is.na(hit_ij)] <- FALSE; C_morph[i, j] <- sum(w_i * as.numeric(hit_ij))
    }
  }
  C_morph_sym <- 0.5 * (C_morph + t(C_morph)); diag(C_morph_sym) <- 0.0; S_morph <- 0.60 * K_mahal + 0.40 * C_morph_sym; diag(S_morph) <- 0.0
  list(S_morph = S_morph, K_mahal = K_mahal, C_morph_sym = C_morph_sym, D2_mahal = D2_mahal, morph_features = feat_cols, p_morph = p)
}

compute_hopkins_phenological_affinity <- function(morph_df, taxa, anchor_mask = NULL) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); n_taxa <- length(taxa); n <- nrow(morph_df)
  S_pheno <- matrix(NA_real_, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa)); diag(S_pheno) <- 0.0
  if (!("doy" %in% names(morph_df)) && !("hopkins_adjusted_doy" %in% names(morph_df))) {
    return(list(S_pheno = S_pheno, has_pheno = FALSE, delta_doy_resid = rep(NA_real_, n), pheno_fdr_p = numeric(0)))
  }
  doy_vec <- if ("doy" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$doy)) else rep(NA_real_, n)
  hop_in <- if ("hopkins_adjusted_doy" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$hopkins_adjusted_doy)) else rep(NA_real_, n)
  has_valid_doy <- (is.finite(doy_vec) & doy_vec >= 1 & doy_vec <= 366) | (is.finite(hop_in) & hop_in >= 1 & hop_in <= 366)
  is_anth <- if (exists("resolve_active_anthesis_mask", mode = "function")) {
    resolve_active_anthesis_mask(morph_df)
  } else {
    st_col <- intersect(c("voucher_phenological_state", "phenological_Stage", "phenological_state"), names(morph_df))[1L]
    if (!is.na(st_col)) {
      st <- tolower(trimws(as.character(morph_df[[st_col]]))); unk <- is.na(st) | !nzchar(st) | st %in% c("na", "unknown")
      ok <- unk | grepl("anthesis|flower|floret|bloom", st)
      ok[!unk & grepl("bud_only|fruiting|fruit|senescent|sterile|vegetative|pre[_-]?anthesis|post[_-]?anthesis", st)] <- FALSE; ok
    } else rep(TRUE, n)
  }
  anthesis_mask <- has_valid_doy & is_anth
  if (sum(anthesis_mask) < 4L) {
    return(list(S_pheno = S_pheno, has_pheno = FALSE, delta_doy_resid = rep(NA_real_, n), pheno_fdr_p = numeric(0)))
  }
  if (!is.logical(anchor_mask) || length(anchor_mask) != n) anchor_mask <- select_taxon_anchors(morph_df)$anchor_mask
  noise_m <- if ("is_confident_label_noise" %in% names(morph_df)) (!is.na(morph_df$is_confident_label_noise) & as.logical(morph_df$is_confident_label_noise)) else rep(FALSE, n)
  wt_all <- .resolve_complex_specimen_weights(morph_df)
  if (!any(is.finite(hop_in)) && exists("compute_hopkins_adjusted_phenology", mode = "function")) {
    h_df <- tryCatch(compute_hopkins_adjusted_phenology(morph_df), error = function(e) NULL)
    if (!is.null(h_df) && "hopkins_adjusted_doy" %in% names(h_df)) hop_in <- suppressWarnings(as.numeric(h_df$hopkins_adjusted_doy))
  }
  doy_eff <- ifelse(is.finite(doy_vec) & doy_vec >= 1 & doy_vec <= 366, doy_vec, hop_in)
  is_cent <- .is_step07_centroid_coord(morph_df)
  lat_vec <- if ("decimalLatitude" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$decimalLatitude)) else rep(NA_real_, n)
  elev_col <- intersect(c("elevation_m", "elevation", "minimumElevationInMeters"), names(morph_df))
  elev_vec <- if (length(elev_col) > 0L) suppressWarnings(as.numeric(morph_df[[elev_col[1L]]])) else rep(NA_real_, n)
  lat_vec[is_cent] <- NA_real_; elev_vec[is_cent] <- NA_real_
  ref_lat <- if (any(is.finite(lat_vec[anthesis_mask & !noise_m]))) stats::median(lat_vec[anthesis_mask & !noise_m], na.rm = TRUE) else 35.0
  ref_elev <- if (any(is.finite(elev_vec[anthesis_mask & !noise_m]))) stats::median(elev_vec[anthesis_mask & !noise_m], na.rm = TRUE) else 150.0
  shift_fb <- ifelse(is.finite(lat_vec), 4.0 * (lat_vec - ref_lat), 0.0) + ifelse(is.finite(elev_vec), (4.0 / 121.92) * (elev_vec - ref_elev), 0.0)
  hop_eff <- ifelse(is.finite(hop_in) & hop_in >= 1 & hop_in <= 366, hop_in, ((doy_eff - shift_fb - 1.0) %% 365.25) + 1.0)
  theta_adj <- 2.0 * pi * (hop_eff / 365.25)
  alpha_eb <- min(1.0, max(0.0, (sum(anthesis_mask & anchor_mask & !noise_m) - 5.0) / 15.0))
  mu_theta <- rep(NA_real_, n_taxa); names(mu_theta) <- taxa; has_tx_pheno <- stats::setNames(logical(n_taxa), taxa)
  delta_doy_resid <- rep(NA_real_, n)
  for (tx in taxa) {
    idx_ref <- which(morph_df$species_standardized == tx & anthesis_mask & anchor_mask & !noise_m)
    if (length(idx_ref) == 0L) idx_ref <- which(morph_df$species_standardized == tx & anthesis_mask & !noise_m)
    idx_all_tx <- which(morph_df$species_standardized == tx & anthesis_mask)
    if (length(idx_ref) > 0L) {
      has_tx_pheno[[tx]] <- TRUE; w_r <- wt_all[idx_ref]; w_r <- w_r / max(sum(w_r), 1e-12)
      mu_tx <- atan2(sum(w_r * sin(theta_adj[idx_ref])), sum(w_r * cos(theta_adj[idx_ref]))); mu_theta[tx] <- mu_tx
      diff_rad <- atan2(sin(theta_adj[idx_all_tx] - mu_tx), cos(theta_adj[idx_all_tx] - mu_tx))
      delta_doy_resid[idx_all_tx] <- diff_rad * (365.25 / (2.0 * pi))
    }
  }
  .circular_watson_williams_or_perm_p <- function(th1, th2) {
    th1 <- th1[is.finite(th1)]; th2 <- th2[is.finite(th2)]; n1 <- length(th1); n2 <- length(th2); N <- n1 + n2
    if (n1 < 2L || n2 < 2L) return(1.0)
    R1 <- sqrt(sum(sin(th1))^2 + sum(cos(th1))^2); R2 <- sqrt(sum(sin(th2))^2 + sum(cos(th2))^2); th_all <- c(th1, th2)
    R_tot <- sqrt(sum(sin(th_all))^2 + sum(cos(th_all))^2); rbar <- (R1 + R2) / N; denom_ww <- max(N - R1 - R2, 1e-8)
    num_ww <- max(0.0, R1 + R2 - R_tot); k_corr <- if (rbar < 0.95) (1.0 + 3.0 / max(8.0 * max(rbar, 0.15), 1.0)) else 1.0
    F_obs <- k_corr * (N - 2.0) * num_ww / denom_ww; p_ww <- stats::pf(F_obs, df1 = 1.0, df2 = max(1.0, N - 2.0), lower.tail = FALSE)
    if (N <= 60L) {
      return(.with_preserved_rng_seed(42L + N, function() {
        B_perm <- 199L; s_all <- sin(th_all); c_all <- cos(th_all); ge_cnt <- 0L
        for (b in seq_len(B_perm)) {
          idx1 <- sample.int(N, size = n1, replace = FALSE); R1_b <- sqrt(sum(s_all[idx1])^2 + sum(c_all[idx1])^2)
          R2_b <- sqrt(sum(s_all[-idx1])^2 + sum(c_all[-idx1])^2)
          if ((R1_b + R2_b - R_tot) >= (num_ww - 1e-12)) ge_cnt <- ge_cnt + 1L
        }
        (ge_cnt + 1.0) / (B_perm + 1.0)
      }))
    }
    pmin(1.0, pmax(0.0, as.numeric(p_ww)))
  }
  pair_ij <- list()
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j && has_tx_pheno[[taxa[i]]] && has_tx_pheno[[taxa[j]]]) {
    arc_rad <- abs(atan2(sin(mu_theta[i] - mu_theta[j]), cos(mu_theta[i] - mu_theta[j]))); delta_circ <- arc_rad * (365.25 / (2.0 * pi))
    S_pheno[i, j] <- exp(-(delta_circ^2) / (2.0 * (21.0^2)))
    if (i < j) pair_ij[[length(pair_ij) + 1L]] <- c(i, j)
  }
  eval_ph_p <- function(ij) {
    .circular_watson_williams_or_perm_p(theta_adj[morph_df$species_standardized == taxa[ij[1L]] & anthesis_mask & !noise_m],
                                        theta_adj[morph_df$species_standardized == taxa[ij[2L]] & anthesis_mask & !noise_m])
  }
  nw_ph <- if (length(pair_ij) >= 2L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, length(pair_ij)) else 1L
  raw_p <- if (length(pair_ij) > 0L) unlist(if (nw_ph > 1L) packera_mclapply(pair_ij, eval_ph_p, mc.cores = nw_ph) else lapply(pair_ij, eval_ph_p), use.names = FALSE) else numeric(0)
  pheno_fdr_p <- if (length(raw_p) > 0L) stats::p.adjust(raw_p, method = "BH") else numeric(0)
  list(S_pheno = S_pheno, has_pheno = sum(has_tx_pheno) >= 2L, delta_doy_resid = delta_doy_resid, pheno_fdr_p = pheno_fdr_p, alpha_eb = alpha_eb)
}
