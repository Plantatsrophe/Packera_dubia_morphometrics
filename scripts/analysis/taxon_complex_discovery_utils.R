# ==============================================================================
# Module: scripts/analysis/taxon_complex_discovery_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline (NCU)
# Description: Unified Step 06 Unsupervised Quantitative Taxon-Complex Discovery
#   & Multimodal Affinity Engine: OOF Cohen's kappa weights (Dawid & Skene 1979),
#   authority-weighted Kish-N_eff Ledoit-Wolf morphometric affinity, anchor-clean
#   ecological KDE Schoener's D, Hopkins-adjusted circular phenological affinity
#   (Hopkins 1918; Pewsey et al. 2013; Mardia & Jupp 2000), and Leiden/Louvain
#   graph modularity + metric Ward.D2 hierarchical complex delimitation
#   (Murtagh & Legendre 2014; Padial et al. 2010; Schlick-Steiner et al. 2010).
# ==============================================================================
.with_preserved_rng_seed <- function(seed, expr_fn) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(as.integer(seed))
  expr_fn()
}
.std_complex_taxon <- function(x) {
  if (exists("standardize_packera_taxon", mode = "function")) {
    res <- standardize_packera_taxon(x)
    return(ifelse(is.na(res) | res == "Unknown", as.character(x), res))
  }
  if (exists("standardize_taxon_name", mode = "function")) {
    res <- standardize_taxon_name(x)
    return(ifelse(is.na(res), as.character(x), res))
  }
  as.character(x)
}
.parse_step06_complex_semicolon <- function(df, col_nm) {
  if (!(col_nm %in% names(df))) return(character(0))
  raw <- trimws(as.character(stats::na.omit(df[[col_nm]]))); raw <- raw[nzchar(raw)]
  if (length(raw) == 0L) return(character(0))
  toks <- unique(trimws(unlist(strsplit(raw, ";", fixed = TRUE))))
  intersect(toks[nzchar(toks) & !(tolower(toks) %in% c("na", "none", "null"))], names(df))
}
.sort_step06_complex_pcs <- function(cols, pfx) {
  if (length(cols) <= 1L) return(cols)
  cols[order(suppressWarnings(as.integer(sub(paste0("^", pfx), "", cols))), cols)]
}
.resolve_complex_specimen_weights <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(numeric(0))
  if ("authority_weight" %in% names(df)) {
    aw <- suppressWarnings(as.numeric(df$authority_weight))
    if (any(is.finite(aw) & aw > 0)) {
      aw[!is.finite(aw) | aw <= 0] <- 0.50
      return(pmin(4.0, pmax(0.15, aw)))
    }
  }
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else list(is_type_specimen = rep(FALSE, n))
  is_type <- as.logical(type_info$is_type_specimen); is_type[is.na(is_type)] <- FALSE
  tier <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier) else rep("Tier_1_Gold", n)
  prom <- if ("is_promoted_bronze" %in% names(df)) (!is.na(df$is_promoted_bronze) & as.logical(df$is_promoted_bronze)) else rep(FALSE, n)
  wt <- ifelse(
    is_type | grepl("type", tier, ignore.case = TRUE), 4.0, ifelse(
      grepl("tier_?1|gold", tier, ignore.case = TRUE), 3.0, ifelse(
        grepl("tier_?2|silver", tier, ignore.case = TRUE), 1.5, ifelse(prom, 1.0, 0.5))))
  noise <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
  wt[noise & !is_type] <- 0.05; wt
}
# 1. Clean Basal_Rosette Anchor Selection (Excluding Confident Noise, Passive & Misidentifications)
select_taxon_anchors <- function(morph_df) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); n <- nrow(morph_df)
  if (n == 0L) return(list(anchor_df = morph_df, anchor_mask = logical(0), n_gold_silver = integer(0), n_type_specimens = integer(0)))
  is_t1_pure <- rep(TRUE, n)
  if ("morph_tier_source" %in% names(morph_df)) {
    mts <- as.character(morph_df$morph_tier_source)
    if (any(!is.na(mts) & nzchar(mts))) is_t1_pure <- !is.na(mts) & (mts == "Basal_Rosette")
  } else if ("basal_leaf_count" %in% names(morph_df)) {
    blc <- suppressWarnings(as.integer(morph_df$basal_leaf_count))
    if (any(!is.na(blc))) is_t1_pure <- !is.na(blc) & (blc >= 1L)
  }
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(morph_df) else list(type_rank = rep(0L, n), is_type_specimen = rep(FALSE, n))
  is_type <- as.logical(type_info$is_type_specimen); is_type[is.na(is_type)] <- FALSE
  if ("is_type_specimen" %in% names(morph_df)) {
    it_col <- as.logical(morph_df$is_type_specimen); it_col[is.na(it_col)] <- FALSE; is_type <- is_type | it_col
  }
  noise <- if ("is_confident_label_noise" %in% names(morph_df)) (!is.na(morph_df$is_confident_label_noise) & as.logical(morph_df$is_confident_label_noise)) else rep(FALSE, n)
  if ("anchor_status" %in% names(morph_df)) {
    ast <- tolower(trimws(as.character(morph_df$anchor_status))); noise <- noise | (!is.na(ast) & grepl("confident_noise|noise_demoted", ast))
  }
  is_pass <- if ("is_passive" %in% names(morph_df)) (!is.na(morph_df$is_passive) & as.logical(morph_df$is_passive)) else rep(FALSE, n)
  mo_out <- if ("morphometric_outlier_flag" %in% names(morph_df)) (!is.na(morph_df$morphometric_outlier_flag) & as.logical(morph_df$morphometric_outlier_flag)) else rep(FALSE, n)
  sub_call <- if ("taxonomic_treatment_subtype" %in% names(morph_df)) as.character(morph_df$taxonomic_treatment_subtype) else rep("", n)
  is_misid <- !is.na(sub_call) & (sub_call == "Misidentification_Reassign")
  prom <- if ("is_promoted_bronze" %in% names(morph_df)) (!is.na(morph_df$is_promoted_bronze) & as.logical(morph_df$is_promoted_bronze)) else rep(FALSE, n)
  det_tier <- if ("determiner_tier" %in% names(morph_df)) {
    ifelse(is.na(morph_df$determiner_tier) | !nzchar(trimws(as.character(morph_df$determiner_tier))), "Tier_3_Bronze", as.character(morph_df$determiner_tier))
  } else rep("Tier_1_Gold", n)
  is_gs <- (grepl("tier_?0|tier_?1|tier_?2|gold|silver|type", det_tier, ignore.case = TRUE) | prom | is_type) & !noise & (!is_pass | is_type) & (!is_misid | is_type)
  clean_non_gs <- !noise & !is_misid & !mo_out
  # Extreme damaged nomenclatural types (> chisq_0.99) excluded from covariance estimation while counted in n_type_specimens
  d2_col <- if ("mahalanobis_d2" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$mahalanobis_d2)) else rep(NA_real_, n)
  p_est <- max(2L, length(.parse_step06_complex_semicolon(morph_df, "selected_morph_features")))
  type_extreme <- is_type & is.finite(d2_col) & (d2_col > stats::qchisq(0.99, df = p_est))
  taxa <- sort(unique(as.character(morph_df$species_standardized))); anchor_mask <- rep(FALSE, n)
  n_gs_vec <- stats::setNames(integer(length(taxa)), taxa); n_type_vec <- stats::setNames(integer(length(taxa)), taxa)
  for (tx in taxa) {
    tx_idx <- which(morph_df$species_standardized == tx); type_idx <- tx_idx[is_type[tx_idx]]; n_type_vec[tx] <- length(type_idx)
    gs_pure_idx <- tx_idx[is_t1_pure[tx_idx] & is_gs[tx_idx] & !mo_out[tx_idx]]
    if (length(gs_pure_idx) < 3L) gs_pure_idx <- tx_idx[is_t1_pure[tx_idx] & is_gs[tx_idx]]
    n_gs_vec[tx] <- sum(is_t1_pure[tx_idx] & is_gs[tx_idx]); clean_pure_idx <- tx_idx[is_t1_pure[tx_idx] & clean_non_gs[tx_idx]]
    clean_any_idx <- tx_idx[clean_non_gs[tx_idx]]
    if (length(gs_pure_idx) >= 5L) anchor_mask[gs_pure_idx] <- TRUE
    else if (length(clean_pure_idx) >= 3L) anchor_mask[clean_pure_idx] <- TRUE
    else if (length(clean_any_idx) >= 2L) anchor_mask[clean_any_idx] <- TRUE
    else anchor_mask[tx_idx] <- TRUE
    intact_types <- type_idx[!type_extreme[type_idx]]
    if (length(intact_types) > 0L) anchor_mask[intact_types] <- TRUE
  }
  list(anchor_df = morph_df[anchor_mask, , drop = FALSE], anchor_mask = anchor_mask, n_gold_silver = n_gs_vec, n_type_specimens = n_type_vec)
}
# 2. Three-Stream Empirical Cohen's Kappa Reliability Weights (Dawid & Skene 1979)
compute_cohens_kappa_stat <- function(y_true, y_pred) {
  valid <- !is.na(y_true) & nzchar(y_true) & !is.na(y_pred) & nzchar(y_pred) & !(y_pred %in% c("Unknown", "Unresolved", "NA"))
  if (sum(valid) < 2L) return(NA_real_)
  yt <- as.character(y_true[valid]); yp <- as.character(y_pred[valid]); lvls <- sort(unique(c(yt, yp)))
  if (length(lvls) < 2L || length(unique(yt)) < 2L) return(max(0.05, min(0.95, mean(yt == yp))))
  tab <- table(factor(yt, levels = lvls), factor(yp, levels = lvls)); n_tot <- sum(tab)
  p_o <- sum(diag(tab)) / n_tot; p_e <- sum((rowSums(tab) / n_tot) * (colSums(tab) / n_tot))
  if (p_e >= 1.0 - 1e-12) return(1.0)
  (p_o - p_e) / (1.0 - p_e)
}
compute_modality_kappa_weights <- function(morph_df, anchor_mask = NULL, ...) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE)
  if (is.logical(anchor_mask) && length(anchor_mask) == nrow(morph_df)) {
    anch_df <- morph_df[anchor_mask, , drop = FALSE]
  } else {
    dots <- list(...)
    if (length(dots) >= 1L && is.logical(dots[[1L]]) && length(dots[[1L]]) == nrow(morph_df)) {
      anch_df <- morph_df[dots[[1L]], , drop = FALSE]
    } else {
      anch_df <- morph_df
    }
  }
  # Use genuine base_species & base_cda_predicted_taxon when in subcluster mode to prevent self-agreement inflation
  yt_morph <- if ("base_species" %in% names(anch_df)) as.character(anch_df$base_species) else as.character(anch_df$species_standardized)
  yp_morph <- if ("base_cda_predicted_taxon" %in% names(anch_df)) as.character(anch_df$base_cda_predicted_taxon) else if ("cda_predicted_taxon" %in% names(anch_df)) as.character(anch_df$cda_predicted_taxon) else rep(NA_character_, nrow(anch_df))
  kappa_morph <- if (any(!is.na(yp_morph) & nzchar(yp_morph))) compute_cohens_kappa_stat(yt_morph, yp_morph) else NA_real_
  w_morph_raw <- if (is.finite(kappa_morph)) max(0.05, min(1.0, kappa_morph)) else if ("w_morph_kappa" %in% names(anch_df) && any(is.finite(as.numeric(anch_df$w_morph_kappa)))) stats::median(as.numeric(anch_df$w_morph_kappa), na.rm = TRUE) else 0.75
  ed_col <- intersect(c("edaphic_best_fit_taxon", "rf_predicted_taxon"), names(anch_df))[1L]
  has_eco_obs <- !is.na(ed_col) && any(!is.na(anch_df[[ed_col]]) & nzchar(trimws(as.character(anch_df[[ed_col]]))))
  if (!has_eco_obs) {
    eco_p <- resolve_step06_complex_eco_predictors(morph_df); has_eco_obs <- length(eco_p) > 0L
  }
  kappa_eco <- if (!is.na(ed_col) && any(!is.na(anch_df[[ed_col]]))) compute_cohens_kappa_stat(yt_morph, anch_df[[ed_col]]) else NA_real_
  w_eco_raw <- if (is.finite(kappa_eco)) max(0.05, min(1.0, kappa_eco)) else if (has_eco_obs && "w_edaphic_kappa" %in% names(anch_df) && any(is.finite(as.numeric(anch_df$w_edaphic_kappa)))) stats::median(as.numeric(anch_df$w_edaphic_kappa), na.rm = TRUE) else if (has_eco_obs) 0.55 else NA_real_
  has_doy <- "doy" %in% names(anch_df) && sum(is.finite(suppressWarnings(as.numeric(anch_df$doy)))) >= 4L
  dd_col <- intersect(c("delta_doy_resid", "delta_doy"), names(anch_df))[1L]
  w_pheno_raw <- if (!is.na(dd_col) && sum(is.finite(suppressWarnings(as.numeric(anch_df[[dd_col]])))) >= 4L) {
    ddv <- abs(suppressWarnings(as.numeric(anch_df[[dd_col]]))); ddv <- ddv[is.finite(ddv)]
    max(0.05, min(0.95, (mean(ddv < 18.0) - 0.25) / 0.75))
  } else if (has_doy && "w_pheno_kappa" %in% names(anch_df) && any(is.finite(as.numeric(anch_df$w_pheno_kappa)))) {
    stats::median(as.numeric(anch_df$w_pheno_kappa), na.rm = TRUE)
  } else if (has_doy) {
    0.50
  } else NA_real_
  present_w <- c(morph = w_morph_raw, edaphic = w_eco_raw, pheno = w_pheno_raw); avail <- is.finite(present_w)
  norm_w <- rep(0.0, length(present_w)); names(norm_w) <- names(present_w)
  if (any(avail)) norm_w[avail] <- present_w[avail] / sum(present_w[avail]) else norm_w["morph"] <- 1.0
  list(
    w_morph_raw = w_morph_raw, w_edaphic_raw = w_eco_raw, w_pheno_raw = w_pheno_raw, w_morph = unname(norm_w["morph"]),
    w_edaphic = unname(norm_w["edaphic"]), w_pheno = unname(norm_w["pheno"]))
}
# 3. Dynamic Morphometric & Ecological Feature Discovery
resolve_step06_complex_morph_features <- function(morph_df, anchor_df = morph_df, n_taxa = 2L) {
  max_p <- max(2L, min(15L, nrow(anchor_df) - as.integer(n_taxa) - 1L))
  ok_num <- function(cc) {
    if (!(cc %in% names(anchor_df)) || !(cc %in% names(morph_df)) || is.logical(anchor_df[[cc]]) || is.character(anchor_df[[cc]])) return(FALSE)
    v <- suppressWarnings(as.numeric(anchor_df[[cc]])); v <- v[is.finite(v)]; length(v) >= 3L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
  }
  sel <- .parse_step06_complex_semicolon(morph_df, "selected_morph_features"); sel <- sel[vapply(sel, ok_num, logical(1))]
  if (length(sel) >= 2L) {
    is_pc <- grepl("^PC[0-9]+$", sel); ord_pc <- .sort_step06_complex_pcs(grep("^PC[0-9]+$", sel, value = TRUE), "PC")
    return(utils::head(unique(c(ord_pc, sel[!is_pc])), max_p))
  }
  pc_cols <- .sort_step06_complex_pcs(grep("^PC[0-9]+$", names(anchor_df), value = TRUE), "PC")
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
resolve_step06_complex_eco_predictors <- function(df) {
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
    .parse_step06_complex_semicolon(df, "primary_env_predictors"), .parse_step06_complex_semicolon(df, "primary_edaph_predictors"),
    .parse_step06_complex_semicolon(df, "primary_clim_predictors")
  ))
  if (length(meta_cols) > 0L) return(meta_cols)
  cands <- ok_env(grep("^(bio[0-9]+_|fine_|ilr_texture_|ssurgo_|polaris_|solus_|topo_)", names(df), value = TRUE))
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
  eco_cols <- resolve_step06_complex_eco_predictors(morph_df)
  if (length(eco_cols) == 0L) {
    return(list(S_env = S_env, D_env = D_env, has_env = FALSE, eco_predictors = character(0), specimen_schoener_d = rep(NA_real_, n)))
  }
  E_raw <- as.matrix(morph_df[, eco_cols, drop = FALSE]); storage.mode(E_raw) <- "double"; has_row_env <- rowSums(is.finite(E_raw)) > 0L
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
    # Impute partial missing cells ONLY on georeferenced rows using taxon-anchor median
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
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); anchor_df <- as.data.frame(anchor_df, stringsAsFactors = FALSE); n_taxa <- length(taxa); feat_cols <- resolve_step06_complex_morph_features(morph_df, anchor_df, n_taxa)
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
  # Kish N_eff-scaled Ledoit-Wolf linear shrinkage toward isotropic target
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
  # Exclude confident label noise & curatorial Misidentification_Reassign sheets from C_morph biological confusion
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
  # Use canonical Step 06 anthesis resolver (properly matching Peak_Anthesis, Early_Anthesis, Late_Anthesis, is_verified_anthesis)
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
  mo_m <- if ("morphometric_outlier_flag" %in% names(morph_df)) (!is.na(morph_df$morphometric_outlier_flag) & as.logical(morph_df$morphometric_outlier_flag)) else rep(FALSE, n)
  wt_all <- .resolve_complex_specimen_weights(morph_df)
  # Honor upstream Step 05/06 hopkins_adjusted_doy or delegate to Step 06 compute_hopkins_adjusted_phenology()
  if (!any(is.finite(hop_in)) && exists("compute_hopkins_adjusted_phenology", mode = "function")) {
    h_df <- tryCatch(compute_hopkins_adjusted_phenology(morph_df), error = function(e) NULL)
    if (!is.null(h_df) && "hopkins_adjusted_doy" %in% names(h_df)) hop_in <- suppressWarnings(as.numeric(h_df$hopkins_adjusted_doy))
  }
  doy_eff <- ifelse(is.finite(doy_vec) & doy_vec >= 1 & doy_vec <= 366, doy_vec, hop_in)
  lat_vec <- if ("decimalLatitude" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$decimalLatitude)) else rep(NA_real_, n)
  elev_col <- intersect(c("elevation_m", "elevation", "minimumElevationInMeters"), names(morph_df))
  elev_vec <- if (length(elev_col) > 0L) suppressWarnings(as.numeric(morph_df[[elev_col[1L]]])) else rep(NA_real_, n)
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
  raw_p <- numeric(0)
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    if (has_tx_pheno[[taxa[i]]] && has_tx_pheno[[taxa[j]]]) {
      arc_rad <- abs(atan2(sin(mu_theta[i] - mu_theta[j]), cos(mu_theta[i] - mu_theta[j]))); delta_circ <- arc_rad * (365.25 / (2.0 * pi))
      S_pheno[i, j] <- exp(-(delta_circ^2) / (2.0 * (21.0^2)))
      if (i < j) {
        th_i <- theta_adj[morph_df$species_standardized == taxa[i] & anthesis_mask & !noise_m]
        th_j <- theta_adj[morph_df$species_standardized == taxa[j] & anthesis_mask & !noise_m]
        raw_p <- c(raw_p, .circular_watson_williams_or_perm_p(th_i, th_j))
      }
    }
  }
  pheno_fdr_p <- if (length(raw_p) > 0L) stats::p.adjust(raw_p, method = "BH") else numeric(0)
  list(S_pheno = S_pheno, has_pheno = sum(has_tx_pheno) >= 2L, delta_doy_resid = delta_doy_resid, pheno_fdr_p = pheno_fdr_p, alpha_eb = alpha_eb)
}
compute_multimodal_complex_affinity <- function(morph_res, eco_res, pheno_res, kappa_wts, taxa) {
  n_taxa <- length(taxa); S_m <- morph_res$S_morph; S_e <- eco_res$S_env; S_p <- pheno_res$S_pheno
  wm <- if (is.finite(kappa_wts$w_morph_raw)) kappa_wts$w_morph_raw else 0.75
  we <- if (isTRUE(eco_res$has_env) && is.finite(kappa_wts$w_edaphic_raw)) kappa_wts$w_edaphic_raw else 0.0
  wp <- if (isTRUE(pheno_res$has_pheno) && is.finite(kappa_wts$w_pheno_raw)) kappa_wts$w_pheno_raw else 0.0
  W <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    sm_ij <- max(0.0, min(1.0, S_m[i, j]))
    # Morphological proximity gate prevents distantly related sympatric clades on shared substrates from merging
    g_ij <- sqrt(sm_ij); num <- wm * sm_ij; den <- wm
    if (we > 0.0 && is.finite(S_e[i, j])) {
      num <- num + g_ij * we * max(0.0, min(1.0, S_e[i, j])); den <- den + g_ij * we
    }
    if (wp > 0.0 && is.finite(S_p[i, j])) {
      num <- num + g_ij * wp * max(0.0, min(1.0, S_p[i, j])); den <- den + g_ij * wp
    }
    W[i, j] <- if (den > 1e-12) num / den else sm_ij
  }
  W <- 0.5 * (W + t(W)); diag(W) <- 0.0; round(W, 6)
}
# 5. Row-Wise Canonical Triage Subtype Resolution (Zero Fabricated Constants)
evaluate_hierarchical_complex_delimitation <- function(morph_df, pheno_res = NULL, eco_res = NULL, ...) {
  if (!is.list(pheno_res) || !("delta_doy_resid" %in% names(pheno_res))) {
    dots <- list(...)
    if (length(dots) >= 1L && is.list(dots[[1L]]) && ("delta_doy_resid" %in% names(dots[[1L]]))) {
      pheno_res <- dots[[1L]]
      if (length(dots) >= 2L) eco_res <- dots[[2L]]
    }
  }
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); n <- nrow(morph_df)
  if (n == 0L) return(character(0))
  pre_sub <- if ("taxonomic_treatment_subtype" %in% names(morph_df)) trimws(as.character(morph_df$taxonomic_treatment_subtype)) else rep(NA_character_, n)
  has_pre <- !is.na(pre_sub) & nzchar(pre_sub) & !(tolower(pre_sub) %in% c("na", "null", "none"))
  if (all(has_pre)) return(pre_sub)
  # Delegate to canonical Step 06 apply_taxonomic_decision_matrix when loaded, else run strict zero-imputation fallback
  if (exists("apply_taxonomic_decision_matrix", mode = "function")) {
    df_tmp <- morph_df
    if (!("delta_doy_resid" %in% names(df_tmp)) && !is.null(pheno_res) && !is.null(pheno_res$delta_doy_resid)) {
      df_tmp$delta_doy_resid <- pheno_res$delta_doy_resid
    }
    if (!("edaphic_schoener_d" %in% names(df_tmp)) && !is.null(eco_res) && !is.null(eco_res$specimen_schoener_d)) {
      df_tmp$edaphic_schoener_d <- eco_res$specimen_schoener_d
    }
    res_s06 <- tryCatch(apply_taxonomic_decision_matrix(df_tmp), error = function(e) NULL)
    if (!is.null(res_s06) && "taxonomic_treatment_subtype" %in% names(res_s06)) {
      out <- as.character(res_s06$taxonomic_treatment_subtype); out[has_pre] <- pre_sub[has_pre]
      return(out)
    }
  }
  sp_given <- as.character(morph_df$species_standardized)
  cda_pred <- if ("cda_predicted_taxon" %in% names(morph_df)) as.character(morph_df$cda_predicted_taxon) else rep(NA_character_, n)
  cda_prob <- if ("cda_posterior_prob" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_posterior_prob)) else rep(NA_real_, n)
  cda_sec <- if ("cda_second_posterior_prob" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_second_posterior_prob)) else rep(NA_real_, n)
  cda_margin <- if ("cda_posterior_margin" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_posterior_margin)) else ifelse(is.finite(cda_prob) & is.finite(cda_sec), abs(cda_prob - cda_sec), NA_real_)
  ed_pred <- if ("edaphic_best_fit_taxon" %in% names(morph_df)) as.character(morph_df$edaphic_best_fit_taxon) else if ("rf_predicted_taxon" %in% names(morph_df)) as.character(morph_df$rf_predicted_taxon) else rep(NA_character_, n)
  m_flag <- if ("morphometric_outlier_flag" %in% names(morph_df)) (!is.na(morph_df$morphometric_outlier_flag) & as.logical(morph_df$morphometric_outlier_flag)) else rep(FALSE, n)
  gmm_cl <- if ("gmm_cluster" %in% names(morph_df)) suppressWarnings(as.integer(morph_df$gmm_cluster)) else rep(NA_integer_, n)
  gmm_bf <- if ("gmm_2_delta_bic" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$gmm_2_delta_bic)) else rep(NA_real_, n)
  d_edaph <- if ("edaphic_schoener_d" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$edaphic_schoener_d)) else if ("schoener_d" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$schoener_d)) else rep(NA_real_, n)
  if (!is.null(eco_res) && !is.null(eco_res$specimen_schoener_d)) d_edaph <- ifelse(is.finite(d_edaph), d_edaph, eco_res$specimen_schoener_d)
  d_doy <- if (!is.null(pheno_res) && !is.null(pheno_res$delta_doy_resid)) pheno_res$delta_doy_resid else if ("delta_doy_resid" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$delta_doy_resid)) else rep(NA_real_, n)
  is_t1_pure <- rep(TRUE, n)
  if ("morph_tier_source" %in% names(morph_df)) {
    mts <- as.character(morph_df$morph_tier_source)
    if (any(!is.na(mts) & nzchar(mts))) is_t1_pure <- !is.na(mts) & (mts == "Basal_Rosette")
  } else if ("basal_leaf_count" %in% names(morph_df)) {
    blc <- suppressWarnings(as.integer(morph_df$basal_leaf_count))
    if (any(!is.na(blc))) is_t1_pure <- !is.na(blc) & (blc >= 1L)
  }
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(morph_df) else list(is_type_specimen = rep(FALSE, n))
  is_type_spec <- as.logical(type_info$is_type_specimen); is_type_spec[is.na(is_type_spec)] <- FALSE
  p_morph <- max(2L, length(.parse_step06_complex_semicolon(morph_df, "selected_morph_features")), 3L)
  min_cong_d2 <- if ("min_congener_mahal_d2" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$min_congener_mahal_d2)) else rep(NA_real_, n)
  min_cong_tx <- if ("min_congener_taxon" %in% names(morph_df)) as.character(morph_df$min_congener_taxon) else rep(NA_character_, n)
  morph_disagrees <- !is.na(cda_pred) & nzchar(cda_pred) & (cda_pred != sp_given) & is.finite(cda_prob) & (cda_prob >= 0.50)
  edaph_disagrees <- !is.na(ed_pred) & nzchar(ed_pred) & (ed_pred != sp_given) & is.finite(d_edaph) & (d_edaph < 0.45)
  both_alt_agree <- morph_disagrees & edaph_disagrees & (cda_pred == ed_pred)
  inter_modal_mosaic <- morph_disagrees & edaph_disagrees & (cda_pred != ed_pred)
  congener_osr_reassign <- is.finite(min_cong_d2) & (min_cong_d2 <= stats::qchisq(0.90, df = p_morph)) & !is.na(min_cong_tx) & nzchar(min_cong_tx) & (min_cong_tx != sp_given) & m_flag
  single_strong_reassign <- !inter_modal_mosaic & ((morph_disagrees & is.finite(cda_prob) & cda_prob >= 0.80) | congener_osr_reassign)
  is_stage1_misid <- !is_type_spec & (both_alt_agree | single_strong_reassign)
  is_stage2_hybrid <- !is_type_spec & !is_stage1_misid & (inter_modal_mosaic | (is.finite(cda_margin) & cda_margin < 0.25 & (morph_disagrees | edaph_disagrees) & (!is.finite(d_edaph) | d_edaph >= 0.30)))
  has_evol_barrier <- (is.finite(d_doy) & abs(d_doy) >= 18.0) | (is.finite(d_edaph) & d_edaph < 0.30)
  coh_pass <- if ("cryptic_cluster_cohesion_pass" %in% names(morph_df)) { cp <- as.logical(morph_df$cryptic_cluster_cohesion_pass); cp[is.na(cp)] <- FALSE; cp } else rep(TRUE, n)
  is_stage3_cryptic <- !is_type_spec & !is_stage1_misid & !is_stage2_hybrid & is_t1_pure & coh_pass & (is.finite(cda_prob) & cda_prob < 0.70) &
    (m_flag | (!is.na(gmm_cl) & gmm_cl > length(unique(sp_given)) & is.finite(gmm_bf) & gmm_bf >= 6.0)) & has_evol_barrier
  is_stage4_ecotype <- !is_type_spec & !is_stage1_misid & !is_stage2_hybrid & !is_stage3_cryptic &
    ((is.finite(d_edaph) & d_edaph < 0.30) | (is.finite(d_doy) & abs(d_doy) >= 14.0)) & (m_flag | (!is.na(gmm_cl) & gmm_cl > 1L))
  is_stage4_plastic <- !is_type_spec & !is_stage1_misid & !is_stage2_hybrid & !is_stage3_cryptic & !is_stage4_ecotype & m_flag
  calc_sub <- ifelse(is_stage1_misid, "Misidentification_Reassign", ifelse(is_stage2_hybrid, "Hybrid_Introgressant",
      ifelse(is_stage3_cryptic, "Putative_Cryptic_Species", ifelse(is_stage4_ecotype, "Subspecies_Ecotype",
          ifelse(is_stage4_plastic, "Ecophenotypic_Plasticity", "Confirmed_Species")))))
  ifelse(has_pre, pre_sub, calc_sub)
}
# 6. Leiden/Louvain Graph Modularity & Metric Ward.D2 Complex Partitioning
relabel_clusters_sequential <- function(cl_vec, taxa) match(cl_vec, unique(cl_vec[order(taxa)]))
partition_taxon_complexes <- function(W, taxa, n_vouchers_vec, n_gs_vec, leiden_resolutions = c(0.6, 0.8, 1.0, 1.2), max_complexes = 8L) {
  n_taxa <- length(taxa); W <- as.matrix(W); diag(W) <- 0.0; nv_map <- stats::setNames(as.integer(n_vouchers_vec[taxa]), taxa)
  # Metric chord/kernel dissimilarity sqrt(1 - W_ij) satisfies Euclidean geometry for Ward.D2 (Murtagh & Legendre 2014)
  D_mat <- matrix(sqrt(pmax(1.0 - pmin(pmax(as.numeric(W), 0.0), 1.0), 0.0)), nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa)); diag(D_mat) <- 0.0
  dist_obj <- stats::as.dist(D_mat); off_vals <- W[upper.tri(W)]
  tau <- min(0.03, stats::quantile(off_vals[off_vals > 0], probs = 0.50, na.rm = TRUE))
  if (!is.finite(tau) || tau <= 0) tau <- 0.01
  W_adj <- W; W_adj[W_adj < tau] <- 0.0
  for (i in seq_len(n_taxa)) {
    if (sum(W_adj[i, ] > 0) == 0L && max(W[i, ]) > 0) {
      best_j <- which.max(W[i, ]); W_adj[i, best_j] <- W[i, best_j]; W_adj[best_j, i] <- W[best_j, i]
    }
  }
  k_cap <- min(as.integer(max_complexes), max(1L, n_taxa - 1L))
  leiden_comm <- .with_preserved_rng_seed(42L, function() {
    comm <- rep(1L, n_taxa); best_mod <- -Inf
    if (requireNamespace("igraph", quietly = TRUE) && sum(W_adj > 0) > 0) {
      g <- igraph::graph_from_adjacency_matrix(W_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
      for (gm in leiden_resolutions) {
        cl_try <- tryCatch(
          suppressWarnings(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution = gm, n_iterations = 10)),
          error = function(e) tryCatch(suppressWarnings(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution_parameter = gm, n_iterations = 10)), error = function(e2) NULL)
        )
        if (!is.null(cl_try)) {
          memb <- as.integer(igraph::membership(cl_try)); k_m <- length(unique(memb))
          q_val <- tryCatch(igraph::modularity(g, memb, weights = igraph::E(g)$weight), error = function(e) -1.0)
          if (k_m >= 2L && k_m <= k_cap && is.finite(q_val) && q_val > best_mod) { best_mod <- q_val; comm <- memb }
        }
      }
      if (length(unique(comm)) < 2L) {
        lv_try <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$weight), error = function(e) NULL)
        if (!is.null(lv_try)) comm <- as.integer(igraph::membership(lv_try))
      }
    }
    relabel_clusters_sequential(comm, taxa)})
  hc <- stats::hclust(dist_obj, method = "ward.D2"); ward_cluster <- rep(1L, n_taxa); best_sil_mean <- -Inf; best_sil_vec <- rep(0.0, n_taxa)
  if (n_taxa >= 3L && k_cap >= 2L && requireNamespace("cluster", quietly = TRUE)) {
    k_candidates <- 2L:k_cap
    sil_scores <- vapply(k_candidates, function(kv) mean(cluster::silhouette(stats::cutree(hc, k = kv), dist_obj)[, "sil_width"], na.rm = TRUE), numeric(1))
    best_idx <- which.max(sil_scores); best_k <- k_candidates[best_idx]; best_sil_mean <- sil_scores[best_idx]
    ward_cluster <- relabel_clusters_sequential(as.integer(stats::cutree(hc, k = best_k)), taxa)
    best_sil_vec <- as.numeric(cluster::silhouette(ward_cluster, dist_obj)[, "sil_width"])
  }
  leiden_sil_mean <- -Inf; leiden_sil_vec <- best_sil_vec
  if (length(unique(leiden_comm)) >= 2L && length(unique(leiden_comm)) < n_taxa && requireNamespace("cluster", quietly = TRUE)) {
    l_sil <- cluster::silhouette(leiden_comm, dist_obj)
    leiden_sil_mean <- mean(l_sil[, "sil_width"], na.rm = TRUE); leiden_sil_vec <- as.numeric(l_sil[, "sil_width"])
  }
  if (leiden_sil_mean > best_sil_mean + 0.02) { consensus_id <- leiden_comm; sil_width_vec <- leiden_sil_vec }
  else { consensus_id <- ward_cluster; sil_width_vec <- best_sil_vec }
  consensus_id <- relabel_clusters_sequential(consensus_id, taxa); complex_labels <- character(n_taxa); deg_vec <- rowSums(W)
  for (cid in sort(unique(consensus_id))) {
    m_idx <- which(consensus_id == cid); ord <- m_idx[order(-deg_vec[m_idx], -nv_map[m_idx], taxa[m_idx])]
    top_ep <- gsub("\\s*\\[K([0-9]+)\\]", "_K\\1", sub("^(?:Packera|P\\.)\\s+", "", taxa[ord[seq_len(min(2L, length(ord)))]]))
    complex_labels[m_idx] <- paste0("Complex_", cid, "_", paste(top_ep, collapse = "_"))
  }
  nn1 <- character(n_taxa); na1 <- numeric(n_taxa); nn2 <- character(n_taxa); na2 <- numeric(n_taxa); within_mean_aff <- numeric(n_taxa)
  for (i in seq_len(n_taxa)) {
    w_row <- W[i, ]; w_row[i] <- -1.0; ord_n <- order(w_row, decreasing = TRUE); nn1[i] <- if (length(ord_n) >= 1L) taxa[ord_n[1L]] else ""
    na1[i] <- if (length(ord_n) >= 1L) round(max(0.0, W[i, ord_n[1L]]), 4) else 0.0; nn2[i] <- if (length(ord_n) >= 2L) taxa[ord_n[2L]] else ""
    na2[i] <- if (length(ord_n) >= 2L) round(max(0.0, W[i, ord_n[2L]]), 4) else 0.0
    peers_i <- which(consensus_id == consensus_id[i] & seq_len(n_taxa) != i)
    within_mean_aff[i] <- if (length(peers_i) > 0L) round(mean(W[i, peers_i], na.rm = TRUE), 4) else 0.0
  }
  # Canonical single-column contract (zero twin/duplicate columns)
  complexes_df <- data.frame(
    species_standardized = taxa, n_vouchers = as.integer(nv_map[taxa]), n_gold_silver_anchors = as.integer(n_gs_vec[taxa]),
    consensus_complex_id = as.integer(consensus_id), complex_label = complex_labels, leiden_complex_id = as.integer(leiden_comm),
    ward_complex_id = as.integer(ward_cluster), silhouette_width = round(sil_width_vec, 4), within_complex_mean_affinity = within_mean_aff,
    nearest_neighbor_taxon = nn1, nearest_neighbor_affinity = na1, second_neighbor_taxon = nn2, second_neighbor_affinity = na2,
    stringsAsFactors = FALSE)
  list(complexes_df = complexes_df, hc = hc, dist_obj = dist_obj)
}
# 7. Taxon x Morphotype Subcluster Builder (Preserving Base-Species Predictions for Unbiased Kappa)
build_taxon_morphotype_subclusters <- function(morph_df, min_vouchers = 5L) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); raw_sp <- as.character(morph_df$species_standardized); u_sp <- sort(unique(raw_sp[!is.na(raw_sp) & nzchar(raw_sp) & raw_sp != "Unknown"]))
  n_sp <- length(u_sp); has_gmm <- "gmm_cluster" %in% names(morph_df) && length(unique(stats::na.omit(morph_df$gmm_cluster))) >= 2L
  within_sp_gmm <- if (has_gmm) any(vapply(u_sp, function(sp) {
    ks <- stats::na.omit(as.integer(morph_df$gmm_cluster[raw_sp == sp])); length(unique(ks)) >= 2L && sum(table(ks) >= min_vouchers) >= 2L
  }, logical(1))) else FALSE
  if (n_sp > 3L || (n_sp > 2L && (!within_sp_gmm || nrow(morph_df) > 100L))) {
    return(list(morph_df = morph_df, subcluster_mode = FALSE))
  }
  if (nrow(morph_df) < 2L * min_vouchers && !within_sp_gmm && n_sp < 2L) {
    return(list(morph_df = morph_df, subcluster_mode = FALSE))
  }
  k_vec <- if (has_gmm) ifelse(is.na(morph_df$gmm_cluster), 1L, as.integer(morph_df$gmm_cluster)) else rep(1L, nrow(morph_df))
  short_sp <- sub("^Packera\\s+", "P. ", raw_sp); cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
  if (sum(table(cand_lbl) >= min_vouchers) < 3L && nrow(morph_df) >= 3L * min_vouchers) {
    fcols <- resolve_step06_complex_morph_features(morph_df, morph_df, max(2L, n_sp))
    if (length(fcols) == 0L) fcols <- intersect(c("can1", "can2", "can3"), names(morph_df))
    if (length(fcols) >= 1L) {
      X_m <- scale(as.matrix(morph_df[, fcols, drop = FALSE])); X_m[!is.finite(X_m)] <- 0.0
      k_vec <- .with_preserved_rng_seed(42L, function() {
        kv <- k_vec
        for (sp in u_sp) {
          idx_s <- which(raw_sp == sp)
          k_s <- if (n_sp == 1L) min(3L, max(2L, floor(length(idx_s) / min_vouchers))) else min(2L, max(1L, floor(length(idx_s) / min_vouchers)))
          if (k_s >= 2L && length(unique(kv[idx_s])) < k_s && length(idx_s) > k_s) {
            km <- tryCatch(stats::kmeans(X_m[idx_s, , drop = FALSE], centers = k_s, nstart = 10L), error = function(e) NULL)
            if (!is.null(km)) kv[idx_s] <- as.integer(km$cluster)
          }
        }
        kv})
      cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
    }
  }
  tab_c <- table(cand_lbl); valid_c <- names(tab_c)[tab_c >= min_vouchers]
  if (length(valid_c) < 2L) {
    if (n_sp >= 2L) return(list(morph_df = morph_df, subcluster_mode = FALSE))
    return(list(morph_df = morph_df, subcluster_mode = FALSE))
  }
  keep_m <- cand_lbl %in% valid_c; morph_df <- morph_df[keep_m, , drop = FALSE]; morph_df$base_species <- raw_sp[keep_m]
  morph_df$base_cda_predicted_taxon <- if ("cda_predicted_taxon" %in% names(morph_df)) .std_complex_taxon(morph_df$cda_predicted_taxon) else raw_sp[keep_m]
  morph_df$species_standardized <- cand_lbl[keep_m]
  if ("cda_predicted_taxon" %in% names(morph_df)) {
    cp <- morph_df$base_cda_predicted_taxon; cp_short <- sub("^Packera\\s+", "P. ", cp)
    morph_df$cda_predicted_taxon <- vapply(seq_len(nrow(morph_df)), function(i) {
      if (!is.na(cp[i]) && cp[i] %in% valid_c) return(cp[i])
      if (!is.na(cp[i]) && cp[i] == morph_df$base_species[i]) return(morph_df$species_standardized[i])
      m_alt <- valid_c[startsWith(valid_c, paste0(cp_short[i], " ["))]
      if (length(m_alt) > 0L) m_alt[1L] else NA_character_
    }, character(1))
  }
  message(sprintf("[Subcluster Mode] Built %d Taxon x Morphotype clusters (>= %d vouchers): %s", length(valid_c), min_vouchers, paste(valid_c, collapse = ", ")))
  list(morph_df = morph_df, subcluster_mode = TRUE)
}
# 8. Unified Step 06 Complex Discovery Orchestrator
.atomic_write_complex_csv <- function(df, path) {
  if (is.null(path) || !nzchar(trimws(as.character(path)))) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  utils::write.csv(df, file = tmp, row.names = FALSE); file.rename(tmp, path); invisible(path)
}
.atomic_write_complex_lines <- function(lines_vec, path) {
  if (is.null(path) || !nzchar(trimws(as.character(path)))) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  writeLines(as.character(lines_vec), con = tmp); file.rename(tmp, path); invisible(path)
}
run_step06_complex_discovery <- function(full_df, opts = list()) {
  if (is.null(full_df) || nrow(full_df) == 0L) return(NULL)
  full_df <- as.data.frame(full_df, stringsAsFactors = FALSE)
  valid_sp <- !is.na(full_df$species_standardized) & nzchar(trimws(as.character(full_df$species_standardized))) & (as.character(full_df$species_standardized) != "Unknown")
  work_df <- full_df[valid_sp, , drop = FALSE]
  if (nrow(work_df) < 2L) return(NULL)
  min_aff <- if (!is.null(opts$min_affinity) && is.finite(suppressWarnings(as.numeric(opts$min_affinity)))) as.numeric(opts$min_affinity) else 0.03
  min_sub <- if (!is.null(opts$min_subcluster_vouchers) && is.finite(suppressWarnings(as.integer(opts$min_subcluster_vouchers)))) as.integer(opts$min_subcluster_vouchers) else 5L
  sub_res <- build_taxon_morphotype_subclusters(work_df, min_vouchers = min_sub); morph_df <- sub_res$morph_df
  taxa <- sort(unique(as.character(morph_df$species_standardized))); n_taxa <- length(taxa)
  if (n_taxa < 2L) {
    message("[Step 06 Complex Discovery] Single taxon without >=2 qualifying morphotype subclusters; skipping multi-complex graph partitioning.")
    return(NULL)
  }
  n_vouchers_vec <- stats::setNames(as.integer(table(morph_df$species_standardized)[taxa]), taxa); anc <- select_taxon_anchors(morph_df)
  kappa_wts <- compute_modality_kappa_weights(morph_df, anc$anchor_mask)
  message(sprintf(
    "[Step 06 Complex Discovery] Evaluating %d vouchers across %d Packera groups (%d anchors; kappa weights: morph=%.3f, edaphic=%.3f, pheno=%.3f).",
    nrow(morph_df), n_taxa, nrow(anc$anchor_df), kappa_wts$w_morph, kappa_wts$w_edaphic, kappa_wts$w_pheno
  ))
  morph_res <- compute_morphometric_affinity(morph_df, anc$anchor_df, taxa)
  pheno_res <- compute_hopkins_phenological_affinity(morph_df, taxa, anc$anchor_mask)
  eco_res <- compute_ecological_affinity(morph_df, taxa, anc$anchor_mask)
  specimen_subtypes <- evaluate_hierarchical_complex_delimitation(morph_df, pheno_res, eco_res)
  W <- compute_multimodal_complex_affinity(morph_res, eco_res, pheno_res, kappa_wts, taxa)
  part_res <- partition_taxon_complexes(W = W, taxa = taxa, n_vouchers_vec = n_vouchers_vec, n_gs_vec = anc$n_gold_silver)
  cdf <- part_res$complexes_df; hop_ov <- rep(NA_real_, n_taxa); eco_ov <- rep(NA_real_, n_taxa); int_cong <- logical(n_taxa)
  n_misid <- integer(n_taxa); n_hyb <- integer(n_taxa); n_crypt <- integer(n_taxa)
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]; idx_tx <- which(morph_df$species_standardized == tx_i)
    n_misid[i] <- sum(specimen_subtypes[idx_tx] == "Misidentification_Reassign", na.rm = TRUE)
    n_hyb[i] <- sum(specimen_subtypes[idx_tx] == "Hybrid_Introgressant", na.rm = TRUE)
    n_crypt[i] <- sum(specimen_subtypes[idx_tx] == "Putative_Cryptic_Species", na.rm = TRUE); nn_tx <- cdf$nearest_neighbor_taxon[i]
    if (nzchar(nn_tx) && nn_tx %in% taxa) {
      if (isTRUE(pheno_res$has_pheno) && is.finite(pheno_res$S_pheno[tx_i, nn_tx])) hop_ov[i] <- round(pheno_res$S_pheno[tx_i, nn_tx], 4)
      if (isTRUE(eco_res$has_env) && is.finite(eco_res$S_env[tx_i, nn_tx])) eco_ov[i] <- round(eco_res$S_env[tx_i, nn_tx], 4)
      n_mod <- as.integer(morph_res$S_morph[tx_i, nn_tx] >= min_aff) +
               as.integer(is.finite(hop_ov[i]) && hop_ov[i] >= 0.25) +
               as.integer(is.finite(eco_ov[i]) && eco_ov[i] >= 0.25)
      int_cong[i] <- n_mod >= 2L
    }
  }
  part_res$complexes_df$n_type_specimens <- as.integer(anc$n_type_specimens[taxa]); part_res$complexes_df$hopkins_pheno_overlap <- hop_ov
  part_res$complexes_df$eco_schoener_overlap <- eco_ov; part_res$complexes_df$integrative_congruence <- int_cong
  part_res$complexes_df$w_morph_kappa <- round(kappa_wts$w_morph, 4); part_res$complexes_df$w_edaphic_kappa <- round(kappa_wts$w_edaphic, 4)
  part_res$complexes_df$w_pheno_kappa <- round(kappa_wts$w_pheno, 4); part_res$complexes_df$n_misidentification_reassign <- n_misid
  part_res$complexes_df$n_hybrid_introgressant <- n_hyb; part_res$complexes_df$n_putative_cryptic <- n_crypt
  # Reorder columns into canonical 24-column sequence (zero duplicate/twin columns)
  canon_order <- c(
    "species_standardized", "n_vouchers", "n_gold_silver_anchors", "n_type_specimens",
    "consensus_complex_id", "complex_label", "leiden_complex_id", "ward_complex_id", "silhouette_width", "within_complex_mean_affinity",
    "nearest_neighbor_taxon", "nearest_neighbor_affinity", "second_neighbor_taxon", "second_neighbor_affinity",
    "hopkins_pheno_overlap", "eco_schoener_overlap", "integrative_congruence", "w_morph_kappa", "w_edaphic_kappa", "w_pheno_kappa",
    "n_misidentification_reassign", "n_hybrid_introgressant", "n_putative_cryptic")
  part_res$complexes_df <- part_res$complexes_df[, canon_order, drop = FALSE]
  if (!is.null(opts$output_complexes) && nzchar(opts$output_complexes)) {
    .atomic_write_complex_csv(part_res$complexes_df, opts$output_complexes)
  }
  aff_df <- data.frame(species_standardized = taxa, as.data.frame(W, check.names = FALSE), stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(opts$output_affinity) && nzchar(opts$output_affinity)) {
    .atomic_write_complex_csv(aff_df, opts$output_affinity)
  }
  focal_raw <- if (!is.null(opts$focal_taxon)) trimws(as.character(opts$focal_taxon)) else ""; selected_neighbors <- character(0)
  if (nzchar(focal_raw)) {
    focal_std <- .std_complex_taxon(focal_raw); focal_short <- sub("^Packera\\s+", "P. ", focal_std)
    focal_nodes <- taxa[taxa == focal_std | startsWith(taxa, paste0(focal_short, " ["))]
    if (length(focal_nodes) > 0L) {
      cdf_now <- part_res$complexes_df; focal_cids <- unique(cdf_now$consensus_complex_id[cdf_now$species_standardized %in% focal_nodes])
      same_complex <- cdf_now$species_standardized[cdf_now$consensus_complex_id %in% focal_cids]
      aff_mask <- apply(morph_res$C_morph_sym[focal_nodes, , drop = FALSE] >= min_aff, 2L, any)
      selected_neighbors <- sort(setdiff(union(same_complex, taxa[aff_mask]), focal_nodes))
    }
    if (!is.null(opts$output_neighbors) && nzchar(opts$output_neighbors)) {
      .atomic_write_complex_lines(selected_neighbors, opts$output_neighbors)
    }
  } else if (!is.null(opts$output_neighbors) && nzchar(opts$output_neighbors)) {
    .atomic_write_complex_lines(part_res$complexes_df$species_standardized, opts$output_neighbors)
  }
  if (!is.null(opts$output_complex_plot) && nzchar(opts$output_complex_plot)) {
    if (!exists("render_step06_complex_discovery_pdf", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R")
    }
    if (exists("render_step06_complex_discovery_pdf", mode = "function")) {
      tryCatch(render_step06_complex_discovery_pdf(W = W, part_res = part_res, out_pdf = opts$output_complex_plot), error = function(e) warning("Step 06 complex network figure warning: ", e$message))
    }
  }
  if (!isTRUE(opts$demo) && exists("persist_step06_complexes_to_sqlite", mode = "function")) {
    tryCatch(persist_step06_complexes_to_sqlite(part_res$complexes_df, vouchers_path = opts$vouchers, demo = isTRUE(opts$demo)), error = function(e) warning("SQLite taxon_complexes persistence note: ", e$message))
  }
  list(
    complexes_df = part_res$complexes_df, affinity_matrix = W, part_res = part_res, morph_res = morph_res, eco_res = eco_res,
    pheno_res = pheno_res, selected_neighbors = selected_neighbors, subcluster_mode = isTRUE(sub_res$subcluster_mode))
}
