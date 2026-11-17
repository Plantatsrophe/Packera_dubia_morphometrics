# Module: scripts/analysis/taxon_complex_discovery_utils.R (Packera Step 07 Complex Discovery)
# Dawid-Skene OOF kappa weights, Ledoit-Wolf morphometric + Schoener D eco + Hopkins
# circular phenology affinities, Leiden modularity & Ward.D2 complex delimitation (Schlick-Steiner et al. 2010).
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
for (.ca_p in c("scripts/analysis/step07_complex_affinity_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step07_complex_affinity_utils.R"))) {
  if (!is.na(.ca_p) && nzchar(.ca_p) && file.exists(.ca_p)) { source(.ca_p); break }
}
# 1. Clean Basal_Rosette Anchor Selection (Excluding Confident Noise, Passive, DNA-Contradicted & Misidentifications)
select_taxon_anchors <- function(morph_df) {
  morph_df <- as.data.frame(morph_df, stringsAsFactors = FALSE); n <- nrow(morph_df)
  if (n == 0L) return(list(anchor_df = morph_df, anchor_mask = logical(0), n_gold_silver = integer(0), n_type_specimens = integer(0)))
  is_t1_pure <- rep(TRUE, n)
  if ("morph_tier_source" %in% names(morph_df) && any(!is.na(mts <- as.character(morph_df$morph_tier_source)) & nzchar(mts))) { is_t1_pure <- !is.na(mts) & (mts == "Basal_Rosette") }
  else if ("basal_leaf_count" %in% names(morph_df) && any(!is.na(blc <- suppressWarnings(as.integer(morph_df$basal_leaf_count))))) { is_t1_pure <- !is.na(blc) & (blc >= 1L) }
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(morph_df) else list(type_rank = rep(0L, n), is_type_specimen = rep(FALSE, n))
  is_type <- as.logical(type_info$is_type_specimen); is_type[is.na(is_type)] <- FALSE
  if ("is_type_specimen" %in% names(morph_df)) { it_col <- as.logical(morph_df$is_type_specimen); it_col[is.na(it_col)] <- FALSE; is_type <- is_type | it_col }
  dna <- .complex_dna_tier0_flags(morph_df, is_type)
  noise <- if ("is_confident_label_noise" %in% names(morph_df)) (!is.na(morph_df$is_confident_label_noise) & as.logical(morph_df$is_confident_label_noise)) else rep(FALSE, n)
  if ("anchor_status" %in% names(morph_df)) { ast <- tolower(trimws(as.character(morph_df$anchor_status))); noise <- noise | (!is.na(ast) & grepl("confident_noise|noise_demoted", ast)) }
  is_pass <- if ("is_passive" %in% names(morph_df)) (!is.na(morph_df$is_passive) & as.logical(morph_df$is_passive)) else rep(FALSE, n)
  mo_out <- if ("morphometric_outlier_flag" %in% names(morph_df)) (!is.na(morph_df$morphometric_outlier_flag) & as.logical(morph_df$morphometric_outlier_flag)) else rep(FALSE, n)
  sub_call <- if ("taxonomic_treatment_subtype" %in% names(morph_df)) as.character(morph_df$taxonomic_treatment_subtype) else rep("", n)
  is_misid <- !is.na(sub_call) & (sub_call == "Misidentification_Reassign")
  prom <- if ("is_promoted_bronze" %in% names(morph_df)) (!is.na(morph_df$is_promoted_bronze) & as.logical(morph_df$is_promoted_bronze)) else rep(FALSE, n)
  det_tier <- if ("determiner_tier" %in% names(morph_df)) ifelse(is.na(morph_df$determiner_tier) | !nzchar(trimws(as.character(morph_df$determiner_tier))), "Tier_3_Bronze", as.character(morph_df$determiner_tier)) else rep("Tier_1_Gold", n)
  is_gs <- ((grepl("tier_?0|tier_?1|tier_?2|gold|silver|type", det_tier, ignore.case = TRUE) | prom | is_type) & !noise & (!is_pass | is_type) & (!is_misid | is_type) & !dna$dna_bad) | dna$dna_ok
  clean_non_gs <- !noise & !is_misid & !mo_out & !dna$dna_bad
  d2_col <- if ("mahalanobis_d2" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$mahalanobis_d2)) else rep(NA_real_, n)
  p_est <- max(2L, length(.parse_step07_complex_semicolon(morph_df, "selected_morph_features")))
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
    intact_types <- type_idx[!type_extreme[type_idx]]; if (length(intact_types) > 0L) anchor_mask[intact_types] <- TRUE
    dna_pure <- tx_idx[dna$dna_ok[tx_idx] & is_t1_pure[tx_idx] & !mo_out[tx_idx]]; if (length(dna_pure) > 0L) anchor_mask[dna_pure] <- TRUE
    if (any(anchor_mask[tx_idx] & !dna$dna_bad[tx_idx])) anchor_mask[tx_idx[dna$dna_bad[tx_idx]]] <- FALSE
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
    eco_p <- resolve_step07_complex_eco_predictors(morph_df); has_eco_obs <- length(eco_p) > 0L
  }
  kappa_eco <- if (!is.na(ed_col) && any(!is.na(anch_df[[ed_col]]))) compute_cohens_kappa_stat(yt_morph, anch_df[[ed_col]]) else NA_real_
  w_eco_raw <- if (is.finite(kappa_eco)) max(0.05, min(1.0, kappa_eco)) else if (has_eco_obs && "w_edaphic_kappa" %in% names(anch_df) && any(is.finite(as.numeric(anch_df$w_edaphic_kappa)))) stats::median(as.numeric(anch_df$w_edaphic_kappa), na.rm = TRUE) else if (has_eco_obs) 0.55 else NA_real_
  doy_col_k <- intersect(c("hopkins_adjusted_doy", "doy"), names(anch_df))[1L]; dd_col <- intersect(c("delta_doy_resid", "delta_doy"), names(anch_df))[1L]
  doy_v_k <- if (!is.na(doy_col_k)) suppressWarnings(as.numeric(anch_df[[doy_col_k]])) else numeric(0)
  anth_k <- if ("is_active_anthesis" %in% names(anch_df)) {
    as.logical(anch_df$is_active_anthesis)
  } else if (exists("resolve_active_anthesis_mask", mode = "function")) {
    resolve_active_anthesis_mask(anch_df)
  } else {
    st_col <- intersect(c("voucher_phenological_state", "phenological_Stage", "phenological_state", "phenophase", "head_state", "reproductive_trait_status", "capitulum_phenology_state"), names(anch_df))[1L]
    if (!is.na(st_col)) {
      st <- tolower(trimws(as.character(anch_df[[st_col]]))); unk <- is.na(st) | !nzchar(st) | st %in% c("na", "unknown")
      ok <- unk | grepl("anthesis|flower|floret|bloom", st)
      ok[!unk & grepl("bud_only|bud|fruiting|fruit|senescent|sterile|vegetative|pre[_-]?anthesis|post[_-]?anthesis|post_fruiting|non[_-]?anthesis", st)] <- FALSE; ok
    } else rep(TRUE, nrow(anch_df))
  }
  anth_k[is.na(anth_k)] <- TRUE
  has_doy <- length(doy_v_k) == nrow(anch_df) && sum(anth_k & is.finite(doy_v_k) & doy_v_k >= 1 & doy_v_k <= 366) >= 4L
  ok_p <- if (has_doy) which(anth_k & is.finite(doy_v_k) & doy_v_k >= 1 & doy_v_k <= 366 & !is.na(yt_morph) & nzchar(yt_morph) & yt_morph != "Unknown") else integer(0)
  tx_tab_p <- if (length(ok_p) >= 4L) table(yt_morph[ok_p]) else integer(0)
  tx_p <- names(tx_tab_p[tx_tab_p >= 2L])
  ok_p <- ok_p[yt_morph[ok_p] %in% tx_p]
  w_pheno_raw <- if ("w_pheno_kappa" %in% names(anch_df) && any(is.finite(suppressWarnings(as.numeric(anch_df$w_pheno_kappa))))) {
    max(0.05, min(0.95, stats::median(suppressWarnings(as.numeric(anch_df$w_pheno_kappa)), na.rm = TRUE)))
  } else if (length(tx_p) >= 2L && length(ok_p) >= 4L) {
    th_p <- 2.0 * pi * (doy_v_k[ok_p] / 365.25); yt_p <- yt_morph[ok_p]
    yp_p <- vapply(seq_along(ok_p), function(ii) {
      best_d <- Inf; best_tx <- tx_p[1L]
      for (tx_k in tx_p) {
        jj <- which(yt_p == tx_k & seq_along(ok_p) != ii)
        if (length(jj) == 0L) jj <- which(yt_p == tx_k)
        mu_k <- atan2(mean(sin(th_p[jj])), mean(cos(th_p[jj])))
        r_k <- max(1e-4, min(0.9995, sqrt(mean(sin(th_p[jj]))^2 + mean(cos(th_p[jj]))^2)))
        sd_rad_k <- max(2.0 * pi * (6.0 / 365.25), sqrt(-2.0 * log(r_k)))
        d_z <- abs(atan2(sin(th_p[ii] - mu_k), cos(th_p[ii] - mu_k))) / sd_rad_k
        if (d_z < best_d) { best_d <- d_z; best_tx <- tx_k }
      }
      best_tx
    }, character(1))
    kp_p <- compute_cohens_kappa_stat(yt_p, yp_p); if (is.finite(kp_p)) max(0.05, min(0.95, kp_p)) else 0.15
  } else if (!is.na(dd_col) && sum(anth_k & is.finite(suppressWarnings(as.numeric(anch_df[[dd_col]])))) >= 4L) { 0.20
  } else if (has_doy) { 0.20 } else NA_real_
  present_w <- c(morph = w_morph_raw, edaphic = w_eco_raw, pheno = w_pheno_raw); avail <- is.finite(present_w)
  norm_w <- rep(0.0, length(present_w)); names(norm_w) <- names(present_w)
  if (any(avail)) norm_w[avail] <- present_w[avail] / sum(present_w[avail]) else norm_w["morph"] <- 1.0
  list(
    w_morph_raw = w_morph_raw, w_edaphic_raw = w_eco_raw, w_pheno_raw = w_pheno_raw, w_morph = unname(norm_w["morph"]),
    w_edaphic = unname(norm_w["edaphic"]), w_pheno = unname(norm_w["pheno"]))
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
# 5. Row-Wise Canonical Triage Subtype Resolution (Padial et al. 2010; Zero Fabricated Constants)
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
  # Delegate to canonical Step 07 apply_taxonomic_decision_matrix when loaded, else run strict zero-imputation fallback
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
  p_morph <- max(2L, length(.parse_step07_complex_semicolon(morph_df, "selected_morph_features")), 3L)
  min_cong_d2 <- if ("min_congener_mahal_d2" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$min_congener_mahal_d2)) else rep(NA_real_, n)
  min_cong_tx <- if ("min_congener_taxon" %in% names(morph_df)) as.character(morph_df$min_congener_taxon) else rep(NA_character_, n)
  morph_disagrees <- !is.na(cda_pred) & nzchar(cda_pred) & (cda_pred != sp_given) & is.finite(cda_prob) & (cda_prob >= 0.50)
  edaph_disagrees <- !is.na(ed_pred) & nzchar(ed_pred) & (ed_pred != sp_given) & is.finite(d_edaph) & (d_edaph < 0.45)
  both_alt_agree <- morph_disagrees & edaph_disagrees & (cda_pred == ed_pred)
  can_hyb_d <- length(unique(sp_given[!is.na(sp_given) & grepl("^Packera[[:space:]]+[a-z][a-z-]+$", sp_given)])) >= 2L
  congener_osr_reassign <- is.finite(min_cong_d2) & (min_cong_d2 <= stats::qchisq(0.90, df = p_morph)) & !is.na(min_cong_tx) & nzchar(min_cong_tx) & (min_cong_tx != sp_given) & m_flag
  single_strong_reassign <- (morph_disagrees & is.finite(cda_prob) & cda_prob >= 0.72) | congener_osr_reassign
  is_stage1_misid <- !is_type_spec & (both_alt_agree | single_strong_reassign)
  is_stage2_hybrid <- can_hyb_d & !is_type_spec & !is_stage1_misid & morph_disagrees & is.finite(cda_prob) & (cda_prob >= 0.35 & cda_prob <= 0.65) & is.finite(cda_margin) & (cda_margin <= 0.22) & (!is.finite(cda_sec) | (cda_sec >= 0.25 & (cda_prob + cda_sec) >= 0.75)) & (!is.finite(d_edaph) | d_edaph >= 0.30)
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
    top_ep <- gsub("[[:space:]]*\\[K([0-9]+)\\]", "_K\\1", sub("^(?:Packera|P\\.)[[:space:]]+", "", taxa[ord[seq_len(min(2L, length(ord)))]]))
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
  if (n_sp == 1L && !within_sp_gmm) {
    eco_col <- if ("ecotype_regime" %in% names(morph_df) && any(!is.na(morph_df$ecotype_regime) & nzchar(trimws(morph_df$ecotype_regime)))) {
      trimws(as.character(morph_df$ecotype_regime))
    } else if ("analysis_group" %in% names(morph_df)) {
      trimws(as.character(morph_df$analysis_group))
    } else {
      character(0)
    }
    eco_ok <- length(eco_col) == nrow(morph_df) && sum(table(eco_col[!is.na(eco_col) & nzchar(eco_col)]) >= min_vouchers) >= 2L
    if (!eco_ok) return(list(morph_df = morph_df, subcluster_mode = FALSE))
    cand_lbl <- ifelse(!is.na(eco_col) & nzchar(eco_col), eco_col, raw_sp)
  } else {
    k_vec <- if (has_gmm) ifelse(is.na(morph_df$gmm_cluster), 1L, as.integer(morph_df$gmm_cluster)) else rep(1L, nrow(morph_df))
    short_sp <- sub("^Packera[[:space:]]+", "P. ", raw_sp); cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
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
    cp <- morph_df$base_cda_predicted_taxon; cp_short <- sub("^Packera[[:space:]]+", "P. ", cp)
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
# 8. Unified Step 07 Complex Discovery Orchestrator
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
run_step07_complex_discovery <- function(full_df, opts = list()) {
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
    message("[Step 07 Complex Discovery] Single taxon without >=2 qualifying morphotype subclusters; skipping multi-complex graph partitioning.")
    return(NULL)
  }
  n_vouchers_vec <- stats::setNames(as.integer(table(morph_df$species_standardized)[taxa]), taxa); anc <- select_taxon_anchors(morph_df)
  kappa_wts <- compute_modality_kappa_weights(morph_df, anc$anchor_mask)
  message(sprintf(
    "[Step 07 Complex Discovery] Evaluating %d vouchers across %d Packera groups (%d anchors; kappa weights: morph=%.3f, edaphic=%.3f, pheno=%.3f).",
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
    focal_std <- .std_complex_taxon(focal_raw); focal_short <- sub("^Packera[[:space:]]+", "P. ", focal_std)
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
    if (!exists("render_step07_complex_discovery_png", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R")
    }
    if (exists("render_step07_complex_discovery_png", mode = "function")) {
      tryCatch(render_step07_complex_discovery_png(W = W, part_res = part_res, out_path = opts$output_complex_plot), error = function(e) warning("Step 07 complex network figure warning: ", e$message))
    }
  }
  if (!isTRUE(opts$demo) && exists("persist_step07_complexes_to_sqlite", mode = "function")) {
    tryCatch(persist_step07_complexes_to_sqlite(part_res$complexes_df, vouchers_path = opts$vouchers, demo = isTRUE(opts$demo)), error = function(e) warning("SQLite taxon_complexes persistence note: ", e$message))
  }
  list(
    complexes_df = part_res$complexes_df, affinity_matrix = W, part_res = part_res, morph_res = morph_res, eco_res = eco_res,
    pheno_res = pheno_res, selected_neighbors = selected_neighbors, subcluster_mode = isTRUE(sub_res$subcluster_mode))
}
