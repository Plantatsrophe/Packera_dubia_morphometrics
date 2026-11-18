if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.step07_safe_bool <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & (x != 0))
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

.has_step07_geno_row_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  if (exists("has_active_genomic_row_mask", mode = "function")) return(has_active_genomic_row_mask(df))
  nms <- names(df); if (!("geno_predicted_taxon" %in% nms)) return(rep(FALSE, n))
  pred <- as.character(df$geno_predicted_taxon)
  pred_ok <- !is.na(pred) & nzchar(trimws(pred)) & (pred != "Unknown")
  tier_col <- if ("geno_linkage_tier" %in% nms) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% nms) as.character(df$linkage_tier) else rep("G1_Direct_Voucher", n)
  tier_ok <- tier_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G2_Conditional_Deme_Hull")
  qc_ok <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(TRUE, n)
  prob_ok <- if ("geno_posterior_prob" %in% nms) is.finite(suppressWarnings(as.numeric(df$geno_posterior_prob))) else rep(TRUE, n)
  pred_ok & tier_ok & qc_ok & prob_ok
}

.step07_type_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  is_t <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df)$is_type_specimen else rep(FALSE, n)
  is_t[is.na(is_t)] <- FALSE
  if ("is_type_specimen" %in% names(df)) {
    tv <- suppressWarnings(as.logical(df$is_type_specimen)); tv[is.na(tv)] <- FALSE; is_t <- is_t | tv
  }
  if ("determiner_tier" %in% names(df)) is_t <- is_t | (trimws(as.character(df$determiner_tier)) == "Tier_0_Type_Anchor")
  is_t
}

.is_step07_centroid_coord <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  if (!("is_centroid_or_inst_coord" %in% names(df))) return(rep(FALSE, n))
  .step07_safe_bool(df$is_centroid_or_inst_coord) & !.step07_type_mask(df)
}

.circ_median_doy_step07 <- function(dv) {
  v <- dv[is.finite(dv) & dv >= 1.0 & dv <= 366.0]; if (length(v) == 0L) return(NA_real_)
  th <- 2.0 * pi * (v / 365.25); mu <- atan2(mean(sin(th)), mean(cos(th))); if (mu < 0.0) mu <- mu + 2.0 * pi
  off <- atan2(sin(th - mu), cos(th - mu)) * (365.25 / (2.0 * pi))
  ((mu * (365.25 / (2.0 * pi)) + stats::median(off) - 1.0) %% 365.25) + 1.0
}

resolve_active_anthesis_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  ver_vec <- if ("is_verified_anthesis" %in% names(df)) suppressWarnings(as.logical(df$is_verified_anthesis)) else rep(NA, n)
  na_cnt <- if ("n_anthesis" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis)) else if ("n_anthesis_heads" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis_heads)) else rep(NA_real_, n)
  state_cols <- intersect(c("voucher_phenological_state", "phenological_state", "phenophase", "head_state"), names(df))
  is_anthesis <- if (length(state_cols) > 0L) {
    st_raw <- tolower(trimws(as.character(df[[state_cols[1L]]])))
    unk_mask <- is.na(st_raw) | !nzchar(st_raw) | st_raw == "unknown"
    ia <- (!unk_mask) & (st_raw %in% c("anthesis", "flowering", "active_anthesis", "open_flower") | grepl("anthesis|flower", st_raw))
    ia[!is.na(st_raw) & (st_raw %in% c("bud", "fruit", "sterile", "vegetative", "fruiting", "senescent") | grepl("pre[_-]?anthesis|post[_-]?anthesis|non[_-]?anthesis", st_raw))] <- FALSE
    if (any(unk_mask)) {
      has_ver <- unk_mask & !is.na(ver_vec)
      if (any(has_ver)) ia[has_ver] <- ver_vec[has_ver]
      rem_unk <- unk_mask & is.na(ver_vec) & is.finite(na_cnt)
      if (any(rem_unk)) ia[rem_unk] <- na_cnt[rem_unk] >= 1.0
    }
    ia
  } else if ("is_verified_anthesis" %in% names(df)) {
    ia <- ver_vec
    if (any(is.na(ia)) && any(is.finite(na_cnt))) ia[is.na(ia) & is.finite(na_cnt)] <- na_cnt[is.na(ia) & is.finite(na_cnt)] >= 1.0
    if (all(is.na(ia)) && ("is_active_anthesis" %in% names(df))) { ia_fb <- suppressWarnings(as.logical(df$is_active_anthesis)); ia_fb[is.na(ia_fb)] <- FALSE; ia_fb } else { ia[is.na(ia)] <- FALSE; ia }
  } else if ("is_active_anthesis" %in% names(df)) {
    ia <- suppressWarnings(as.logical(df$is_active_anthesis)); ia[is.na(ia)] <- FALSE; ia
  } else if ("is_anthesis" %in% names(df)) {
    ia <- as.logical(df$is_anthesis); ia[is.na(ia)] <- FALSE; ia
  } else if ("n_anthesis" %in% names(df)) {
    !is.na(na_cnt) & is.finite(na_cnt) & na_cnt >= 1.0
  } else if ("n_anthesis_heads" %in% names(df)) {
    nh <- suppressWarnings(as.numeric(df$n_anthesis_heads)); !is.na(nh) & is.finite(nh) & nh >= 1.0
  } else rep(TRUE, n)
  defs <- .load_step07_pheno_config_defaults()
  doy_v <- .num_col(df, "doy", n)
  if (all(!is.finite(doy_v)) && ("hopkins_adjusted_doy" %in% names(df))) doy_v <- .num_col(df, "hopkins_adjusted_doy", n)
  has_d <- is.finite(doy_v) & (doy_v >= 1.0) & (doy_v <= 366.0)
  if (any(has_d)) {
    is_anthesis[has_d & (doy_v < defs$min_flowering_doy)] <- FALSE
    coh_med <- .circ_median_doy_step07(doy_v[is_anthesis & has_d])
    sp_col <- if ("species_standardized" %in% names(df)) trimws(as.character(df$species_standardized)) else rep("", n)
    med_by_row <- rep(coh_med, n); ok_sp <- !is.na(sp_col) & nzchar(sp_col) & sp_col != "Unknown"
    for (tx in unique(sp_col[ok_sp])) {
      m_tx <- ok_sp & (sp_col == tx) & is_anthesis & has_d
      if (sum(m_tx) >= 5L && is.finite(tx_med <- .circ_median_doy_step07(doy_v[m_tx]))) med_by_row[ok_sp & (sp_col == tx)] <- tx_med
    }
    is_vernal <- is.finite(med_by_row) & (med_by_row < 170.0)
    is_anthesis[is_vernal & has_d & ((doy_v < defs$min_flowering_doy) | (doy_v > defs$max_flowering_doy) | (doy_v > 240.0))] <- FALSE
  }
  is_anthesis
}

.num_col <- function(df, col_nm, n = nrow(df), default = NA_real_) if (col_nm %in% names(df)) suppressWarnings(as.numeric(df[[col_nm]])) else rep(default, n)
.chr_col <- function(df, col_nm, n = nrow(df), default = NA_character_) if (col_nm %in% names(df)) as.character(df[[col_nm]]) else rep(default, n)
.lgl_col <- function(df, col_nm, n = nrow(df)) if (col_nm %in% names(df)) (!is.na(df[[col_nm]]) & as.logical(df[[col_nm]])) else rep(FALSE, n)

.tier0_geno_ctx <- function(df) {
  n <- nrow(df); nms <- names(df)
  tier_g <- if ("geno_linkage_tier" %in% nms) trimws(as.character(df$geno_linkage_tier)) else if ("linkage_tier" %in% nms) trimws(as.character(df$linkage_tier)) else rep("", n)
  qc_g <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(FALSE, n)
  adm_g <- if ("is_geno_admixed" %in% nms) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
  qry_g <- if ("is_putative_novel_or_hybrid_query" %in% nms) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  pred_g <- trimws(.chr_col(df, "geno_predicted_taxon", n, "")); pred_g[is.na(pred_g)] <- ""
  sp_g <- trimws(.chr_col(df, "species_standardized", n, "")); sp_g[is.na(sp_g)] <- ""
  prob_g <- .num_col(df, "geno_posterior_prob", n); marg_g <- .num_col(df, "geno_margin", n)
  g1_ok <- (tier_g %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc_g & nzchar(pred_g) & (pred_g != "Unknown") & is.finite(prob_g) & (prob_g >= 0.85)
  list(g1_ok = g1_ok, adm_g = adm_g, qry_g = qry_g, pred_g = pred_g, sp_g = sp_g, marg_g = marg_g)
}

resolve_step07_ecological_authority_weights <- function(df, allow_bronze_prior = TRUE) {
  n <- nrow(df); if (n == 0L) return(numeric(0))
  w_out <- resolve_ecological_authority_weights(df, allow_bronze_prior = allow_bronze_prior)
  if (isTRUE(allow_bronze_prior) && length(w_out) == n) w_out <- round(pmax(0.01, w_out), 4)
  if (any(.has_step07_geno_row_mask(df))) {
    is_type_m <- .step07_type_mask(df); gx <- .tier0_geno_ctx(df)
    t0_conc <- !is_type_m & gx$g1_ok & !gx$adm_g & !gx$qry_g & (is.na(gx$marg_g) | gx$marg_g >= 0.40) & (gx$pred_g == gx$sp_g)
    w_out[t0_conc] <- 4.0; w_out[!is_type_m & gx$g1_ok & (gx$pred_g != gx$sp_g)] <- 0.0
  }
  w_out
}

.load_step07_pheno_config_defaults <- function(cfg_path = "config/config.yaml") {
  defs <- list(tau_ecotype = 14.0, tau_cryptic = 18.0, min_flowering_doy = 30.0, max_flowering_doy = 220.0)
  if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    if (is.list(cfg) && is.list(cfg$synthesis)) {
      if (!is.null(cfg$synthesis$ecotype_ddoy_cutoff_days)) defs$tau_ecotype <- as.numeric(cfg$synthesis$ecotype_ddoy_cutoff_days)
      if (!is.null(cfg$synthesis$cryptic_ddoy_cutoff_days)) defs$tau_cryptic <- as.numeric(cfg$synthesis$cryptic_ddoy_cutoff_days)
    }
    ph_c <- if (is.list(cfg) && is.list(cfg$macroecology) && is.list(cfg$macroecology$phenology)) cfg$macroecology$phenology else if (is.list(cfg) && is.list(cfg$phenology)) cfg$phenology else NULL
    if (is.list(ph_c)) {
      mn <- suppressWarnings(as.numeric(ph_c$min_flowering_doy)); mx <- suppressWarnings(as.numeric(ph_c$max_flowering_doy))
      if (length(mn) == 1L && is.finite(mn)) defs$min_flowering_doy <- mn; if (length(mx) == 1L && is.finite(mx)) defs$max_flowering_doy <- mx
    }
  }
  defs
}

.circular_sd_days_from_residuals <- function(ddoy_days_vec) {
  v <- ddoy_days_vec[is.finite(ddoy_days_vec)]; if (length(v) < 2L) return(NA_real_)
  phi <- 2.0 * pi * (v / 365.25); r_clamped <- max(1e-6, min(0.999999, sqrt(mean(sin(phi))^2 + mean(cos(phi))^2)))
  (365.25 / (2.0 * pi)) * sqrt(-2.0 * log(r_clamped))
}

.is_step07_extralimital_non_type <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0)); nms <- names(df)
  is_type_m <- .step07_type_mask(df)
  u_stat <- if ("usda_county_status" %in% nms) tolower(trimws(ifelse(is.na(df$usda_county_status), "", as.character(df$usda_county_status)))) else rep("", n)
  ext_u <- u_stat %in% c("out_of_range", "natureserve_exotic_or_waif")
  exs_c <- if ("exsiccatae_taxonomic_conflict" %in% nms) {
    raw_ec <- df$exsiccatae_taxonomic_conflict
    if (is.logical(raw_ec)) (!is.na(raw_ec) & raw_ec) else (!is.na(suppressWarnings(as.integer(raw_ec))) & suppressWarnings(as.integer(raw_ec)) != 0L)
  } else rep(FALSE, n)
  pas_aud <- if ("anchor_tier_source" %in% nms) (trimws(ifelse(is.na(df$anchor_tier_source), "", as.character(df$anchor_tier_source))) == "Passive_Extralimital_Audit") else rep(FALSE, n)
  flag <- !is_type_m & (ext_u | exs_c | pas_aud)
  if (any(flag) && ("species_standardized" %in% nms)) {
    sp_v <- trimws(ifelse(is.na(df$species_standardized), "", as.character(df$species_standardized)))
    noise_m <- if ("is_confident_label_noise" %in% nms) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
    for (tx in unique(sp_v[flag & nzchar(sp_v)])) {
      tx_idx <- which(sp_v == tx & !noise_m); if (sum(!flag[tx_idx]) < 2L) flag[sp_v == tx] <- FALSE
    }
  }
  flag
}

.resolve_step07_anchor_mask <- function(df) {
  n <- nrow(df); if (n == 0L) return(logical(0))
  tier_col <- if ("determiner_tier" %in% names(df)) {
    if (exists(".standardize_determiner_tier_vec", mode = "function")) .standardize_determiner_tier_vec(df)
    else ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier))
  } else rep("Tier_1_Gold", n)
  is_type_m <- .step07_type_mask(df); promo_mask <- .lgl_col(df, "is_promoted_bronze", n) | .lgl_col(df, "is_type_specimen", n) | is_type_m
  noise_mask <- .lgl_col(df, "is_confident_label_noise", n); passive_mask <- .lgl_col(df, "is_passive", n)
  if ("authority_weight" %in% names(df)) { aw_v <- .num_col(df, "authority_weight", n); if (any(is.finite(aw_v) & aw_v > 0.0)) passive_mask <- passive_mask | (is.finite(aw_v) & aw_v <= 0.0) }
  anchor_mask <- (((tier_col %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze_Promoted")) | promo_mask) & !noise_mask & !passive_mask & !.is_step07_extralimital_non_type(df)) | (is_type_m & !noise_mask)
  if (any(.has_step07_geno_row_mask(df))) {
    gx <- .tier0_geno_ctx(df); g1_base <- gx$g1_ok & !gx$adm_g
    t0_contra <- !is_type_m & g1_base & (gx$pred_g != gx$sp_g); t0_conc <- g1_base & !gx$qry_g & (gx$pred_g == gx$sp_g) & !noise_mask
    anchor_mask <- (anchor_mask | t0_conc) & !t0_contra
    if (sum(anchor_mask, na.rm = TRUE) < 4L) {
      clean_fb <- !noise_mask & !passive_mask & !t0_contra
      anchor_mask <- if (sum(clean_fb, na.rm = TRUE) >= 1L) clean_fb else if (sum(!noise_mask & !t0_contra, na.rm = TRUE) >= 1L) (!noise_mask & !t0_contra) else if (sum(!noise_mask, na.rm = TRUE) >= 1L) !noise_mask else rep(TRUE, n)
    }
  } else if (sum(anchor_mask, na.rm = TRUE) < 4L) {
    clean_fb <- !noise_mask & !passive_mask
    anchor_mask <- if (sum(clean_fb, na.rm = TRUE) >= 1L) clean_fb else if (sum(!noise_mask, na.rm = TRUE) >= 1L) !noise_mask else rep(TRUE, n)
  }
  anchor_mask
}

.compute_taxon_pheno_thresholds <- function(df) {
  n <- nrow(df); defs <- .load_step07_pheno_config_defaults()
  if (n == 0L) return(list(tau_ecotype = numeric(0), tau_cryptic = numeric(0), sigma_by_taxon = numeric(0)))
  sp_vec <- .chr_col(df, "species_standardized", n, "Unknown")
  is_anth <- if ("is_active_anthesis" %in% names(df)) as.logical(df$is_active_anthesis) else resolve_active_anthesis_mask(df); is_anth[is.na(is_anth)] <- TRUE
  ddoy_vec <- if ("delta_doy_resid" %in% names(df)) .num_col(df, "delta_doy_resid", n) else .num_col(df, "delta_doy", n)
  valid_ph <- is_anth & is.finite(ddoy_vec); mo_m <- .lgl_col(df, "morphometric_outlier_flag", n); noise_m <- .lgl_col(df, "is_confident_label_noise", n)
  clean_anch <- .resolve_step07_anchor_mask(df) & !mo_m & !noise_m
  coh_idx <- which(clean_anch & valid_ph); if (length(coh_idx) < 5L) coh_idx <- which(!mo_m & !noise_m & valid_ph); if (length(coh_idx) < 5L) coh_idx <- which(!noise_m & valid_ph)
  sigma_cohort <- if (length(coh_idx) >= 5L) .circular_sd_days_from_residuals(ddoy_vec[coh_idx]) else NA_real_
  u_taxa <- unique(sp_vec[!is.na(sp_vec) & nzchar(sp_vec)])
  tau_eco_map <- stats::setNames(rep(defs$tau_ecotype, length(u_taxa)), u_taxa); tau_cry_map <- stats::setNames(rep(defs$tau_cryptic, length(u_taxa)), u_taxa); sigma_map <- stats::setNames(rep(sigma_cohort, length(u_taxa)), u_taxa)
  if (is.finite(sigma_cohort)) { tau_eco_map[] <- min(28.0, max(10.0, 2.0 * sigma_cohort)); tau_cry_map[] <- min(35.0, max(14.0, 2.5 * sigma_cohort)) }
  for (tx in u_taxa) {
    idx_t <- which(sp_vec == tx); tx_pool <- idx_t[clean_anch[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) < 5L) tx_pool <- idx_t[!mo_m[idx_t] & !noise_m[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) < 5L) tx_pool <- idx_t[!noise_m[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) >= 5L && is.finite(sig_tx <- .circular_sd_days_from_residuals(ddoy_vec[tx_pool]))) {
      sigma_map[[tx]] <- sig_tx; tau_eco_map[[tx]] <- min(28.0, max(10.0, 2.0 * sig_tx)); tau_cry_map[[tx]] <- min(35.0, max(14.0, 2.5 * sig_tx))
    }
  }
  fallback_eco <- if (is.finite(sigma_cohort)) min(28.0, max(10.0, 2.0 * sigma_cohort)) else defs$tau_ecotype
  fallback_cry <- if (is.finite(sigma_cohort)) min(35.0, max(14.0, 2.5 * sigma_cohort)) else defs$tau_cryptic
  list(tau_ecotype = as.numeric(ifelse(sp_vec %in% names(tau_eco_map), tau_eco_map[sp_vec], fallback_eco)),
       tau_cryptic = as.numeric(ifelse(sp_vec %in% names(tau_cry_map), tau_cry_map[sp_vec], fallback_cry)), sigma_by_taxon = sigma_map)
}

compute_hopkins_adjusted_phenology <- function(df) {
  n <- nrow(df); if (n == 0L) return(df)
  is_anthesis <- resolve_active_anthesis_mask(df); df$is_active_anthesis <- is_anthesis
  doy_raw <- .num_col(df, "doy", n); has_valid_doy <- is.finite(doy_raw) & (doy_raw >= 1.0) & (doy_raw <= 366.0)
  anchor_mask <- .resolve_step07_anchor_mask(df); mo_m <- .lgl_col(df, "morphometric_outlier_flag", n); noise_m <- .lgl_col(df, "is_confident_label_noise", n)
  clean_anth_doy <- is_anthesis & has_valid_doy & anchor_mask & !mo_m & !noise_m
  doy_pool <- doy_raw[clean_anth_doy]; if (length(doy_pool) == 0L) doy_pool <- doy_raw[is_anthesis & has_valid_doy & !noise_m]; if (length(doy_pool) == 0L) doy_pool <- doy_raw[has_valid_doy]
  if (length(doy_pool) > 0L) {
    th_pool <- 2.0 * pi * (doy_pool / 365.25); ref_theta <- atan2(mean(sin(th_pool)), mean(cos(th_pool)))
    if (ref_theta < 0.0) ref_theta <- ref_theta + 2.0 * pi; ref_doy <- ref_theta * (365.25 / (2.0 * pi))
  } else { ref_doy <- 120.0; ref_theta <- 2.0 * pi * (ref_doy / 365.25) }
  doy_vec <- ifelse(has_valid_doy, doy_raw, ref_doy)

  is_cent_h <- .is_step07_centroid_coord(df); lat_vec <- .num_col(df, "decimalLatitude", n)
  elev_col <- intersect(c("elevation_m", "elevation", "minimumElevationInMeters"), names(df))
  elev_vec <- if (length(elev_col) > 0L) .num_col(df, elev_col[1L], n) else rep(NA_real_, n)
  mat_col <- intersect(c("bio1_temp_mean", "bio1", "mat_c"), names(df))
  mat_vec <- if (length(mat_col) > 0L) .num_col(df, mat_col[1L], n) else rep(NA_real_, n)
  lat_vec[is_cent_h] <- NA_real_; elev_vec[is_cent_h] <- NA_real_; mat_vec[is_cent_h] <- NA_real_
  pick_pool <- function(v) { p <- v[clean_anth_doy & is.finite(v)]; if (length(p) == 0L) p <- v[is_anthesis & is.finite(v)]; if (length(p) == 0L) p <- v[is.finite(v)]; p }
  ref_lat <- if (length(lp <- pick_pool(lat_vec)) > 0L) stats::median(lp, na.rm = TRUE) else 35.0
  ref_elev <- if (length(ep <- pick_pool(elev_vec)) > 0L) stats::median(ep, na.rm = TRUE) else 150.0
  ref_mat <- if (length(mp <- pick_pool(mat_vec)) > 0L) stats::median(mp, na.rm = TRUE) else 14.5
  lat_shift <- ifelse(is.finite(lat_vec), 4.0 * (lat_vec - ref_lat), 0.0); elev_shift <- ifelse(is.finite(elev_vec), (4.0 / 121.92) * (elev_vec - ref_elev), 0.0); mat_shift <- ifelse(is.finite(mat_vec), -1.5 * (mat_vec - ref_mat), 0.0)
  has_geo_cov <- is.finite(lat_vec) | is.finite(elev_vec); has_mat_cov <- is.finite(mat_vec); has_any_cov <- has_geo_cov | has_mat_cov
  theo_shift <- ifelse(has_geo_cov & has_mat_cov, 0.70 * (lat_shift + elev_shift) + 0.30 * mat_shift, ifelse(has_geo_cov, lat_shift + elev_shift, ifelse(has_mat_cov, mat_shift, 0.0)))
  theta_obs <- 2.0 * pi * (doy_vec / 365.25); theta_theo <- ref_theta + 2.0 * pi * (theo_shift / 365.25); sin_theo <- sin(theta_theo); cos_theo <- cos(theta_theo)
  fit_mask <- is_anthesis & has_valid_doy & has_any_cov & anchor_mask & !mo_m & !noise_m
  if (sum(fit_mask) < 5L) fit_mask <- is_anthesis & has_valid_doy & has_any_cov & anchor_mask & !noise_m
  if (sum(fit_mask) < 5L) fit_mask <- is_anthesis & has_valid_doy & has_any_cov & !noise_m
  cov_raw <- list(decimalLatitude = lat_vec, elevation = elev_vec, bio1_temp_mean = mat_vec)
  cov_terms <- names(cov_raw)[vapply(cov_raw, function(v) { vf <- v[fit_mask & is.finite(v)]; length(vf) >= 3L && is.finite(stats::sd(vf)) && stats::sd(vf) > 1e-6 }, logical(1))]
  alpha_eb <- if (length(cov_terms) > 0L) min(1.0, max(0.0, (sum(fit_mask) - 5.0) / 15.0)) else 0.0; theta_blend <- theta_theo; theta_blend_ref <- ref_theta
  if (alpha_eb > 0.0 && length(cov_terms) > 0L) {
    w_hop <- resolve_step07_ecological_authority_weights(df, allow_bronze_prior = TRUE); if (length(w_hop) != n) w_hop <- rep(1.0, n)
    w_hop[!is.finite(w_hop) | w_hop < 0.0] <- 0.25; if (!any(w_hop[fit_mask] > 0.0)) w_hop[fit_mask] <- 1.0
    hop_df <- data.frame(sin_doy = sin(theta_obs), cos_doy = cos(theta_obs), decimalLatitude = ifelse(is.finite(lat_vec), lat_vec, ref_lat),
                         elevation = ifelse(is.finite(elev_vec), elev_vec, ref_elev), bio1_temp_mean = ifelse(is.finite(mat_vec), mat_vec, ref_mat), .w_hop = w_hop)
    ref_df <- data.frame(decimalLatitude = ref_lat, elevation = ref_elev, bio1_temp_mean = ref_mat); rhs_str <- paste(cov_terms, collapse = " + ")
    lm_sin <- tryCatch(stats::lm(stats::as.formula(paste("sin_doy ~", rhs_str)), data = hop_df[fit_mask, , drop = FALSE], weights = .w_hop), error = function(e) NULL)
    lm_cos <- tryCatch(stats::lm(stats::as.formula(paste("cos_doy ~", rhs_str)), data = hop_df[fit_mask, , drop = FALSE], weights = .w_hop), error = function(e) NULL)
    if (!is.null(lm_sin) && !is.null(lm_cos)) {
      pred_s <- tryCatch(as.numeric(stats::predict(lm_sin, newdata = hop_df)), error = function(e) rep(NA_real_, n)); pred_c <- tryCatch(as.numeric(stats::predict(lm_cos, newdata = hop_df)), error = function(e) rep(NA_real_, n))
      ref_s <- tryCatch(as.numeric(stats::predict(lm_sin, newdata = ref_df)), error = function(e) NA_real_); ref_c <- tryCatch(as.numeric(stats::predict(lm_cos, newdata = ref_df)), error = function(e) NA_real_)
      if (all(is.finite(pred_s)) && all(is.finite(pred_c)) && is.finite(ref_s) && is.finite(ref_c)) {
        th_e <- atan2(pred_s, pred_c); th_er <- atan2(ref_s, ref_c)
        theta_blend <- atan2(alpha_eb * sin(th_e) + (1.0 - alpha_eb) * sin_theo, alpha_eb * cos(th_e) + (1.0 - alpha_eb) * cos_theo)
        theta_blend_ref <- atan2(alpha_eb * sin(th_er) + (1.0 - alpha_eb) * sin(ref_theta), alpha_eb * cos(th_er) + (1.0 - alpha_eb) * cos(ref_theta))
      }
    }
  }

  delta_theta_env <- ifelse(has_any_cov, atan2(sin(theta_blend - theta_blend_ref), cos(theta_blend - theta_blend_ref)), 0.0)
  theta_adj <- theta_obs - delta_theta_env; shift_days <- delta_theta_env * (365.25 / (2.0 * pi))
  adj_doy_calc <- ifelse(has_valid_doy, round(((doy_vec - shift_days - 1.0) %% 365.25) + 1.0, 2), NA_real_)
  adj_in <- .num_col(df, "hopkins_adjusted_doy", n); has_valid_adj <- is.finite(adj_in) & (adj_in >= 1.0) & (adj_in <= 366.0)
  has_pheno_doy <- has_valid_doy | has_valid_adj
  df$hopkins_adjusted_doy <- ifelse(has_valid_adj, round(adj_in, 2), adj_doy_calc)
  df$sin_doy <- ifelse(has_valid_doy, round(sin(theta_obs), 4), NA_real_); df$cos_doy <- ifelse(has_valid_doy, round(cos(theta_obs), 4), NA_real_)

  theta_adj_use <- ifelse(is.finite(df$hopkins_adjusted_doy), 2.0 * pi * (df$hopkins_adjusted_doy / 365.25), theta_adj)
  circ_resid_days <- numeric(n); sp_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unassigned_Packera", n)
  t0_hop <- if (any(.has_step07_geno_row_mask(df))) {
    if (exists("get_tier0_rectified_species_and_evictions", mode = "function")) get_tier0_rectified_species_and_evictions(df)
    else if (exists("resolve_tier0_anchor_overrides", mode = "function")) resolve_tier0_anchor_overrides(df) else NULL
  } else NULL
  for (tx in unique(sp_vec[!is.na(sp_vec)])) {
      idx_tx <- which(sp_vec == tx)
      if (!is.null(t0_hop)) {
        ev_tx <- t0_hop$evict_from_given_anchor & (!t0_hop$is_tier0_dna_anchor | t0_hop$rectified_taxon != tx)
        rect_tx <- which(t0_hop$is_tier0_dna_anchor & !is.na(t0_hop$rectified_taxon) & (t0_hop$rectified_taxon == tx) & is_anthesis & has_pheno_doy & !noise_m)
        anth_tx <- unique(c(idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & anchor_mask[idx_tx] & !mo_m[idx_tx] & !noise_m[idx_tx] & !ev_tx[idx_tx]], rect_tx))
        if (length(anth_tx) < 2L) anth_tx <- unique(c(idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & !noise_m[idx_tx] & !ev_tx[idx_tx]], rect_tx))
        if (length(anth_tx) == 0L) anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx]]
        w_h <- if (length(anth_tx) > 0L) ifelse(t0_hop$is_tier0_dna_anchor[anth_tx] & (t0_hop$rectified_taxon[anth_tx] == tx), ifelse(is.finite(t0_hop$tier0_weight[anth_tx]), t0_hop$tier0_weight[anth_tx], 4.0), 1.0) else numeric(0)
        theta_tx_c <- if (length(anth_tx) > 0L) atan2(sum(w_h * sin(theta_adj_use[anth_tx]), na.rm = TRUE), sum(w_h * cos(theta_adj_use[anth_tx]), na.rm = TRUE)) else ref_theta
      } else {
        anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & anchor_mask[idx_tx] & !mo_m[idx_tx] & !noise_m[idx_tx]]
        if (length(anth_tx) < 2L) anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & !noise_m[idx_tx]]
        if (length(anth_tx) == 0L) anth_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx]]
        theta_tx_c <- if (length(anth_tx) > 0L) atan2(mean(sin(theta_adj_use[anth_tx]), na.rm = TRUE), mean(cos(theta_adj_use[anth_tx]), na.rm = TRUE)) else ref_theta
      }
      d_th <- theta_adj_use[idx_tx] - theta_tx_c; circ_resid_days[idx_tx] <- atan2(sin(d_th), cos(d_th)) * (365.25 / (2.0 * pi))
    }
  circ_resid_gated <- ifelse(is_anthesis & has_pheno_doy, round(circ_resid_days, 2), NA_real_)
  ddoy_in <- if ("delta_doy_resid" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy_resid)) else if ("delta_doy" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy)) else rep(NA_real_, n)
  ddoy_use <- ifelse(is_anthesis & (has_pheno_doy | is.finite(ddoy_in)) & is.finite(ddoy_in), ddoy_in, circ_resid_gated)
  ddoy_use[!is_anthesis] <- NA_real_; df$delta_doy_resid <- ddoy_use; df$delta_doy <- ddoy_use
  df
}

compute_cohens_kappa <- function(y_true, y_pred, default_kappa = 0.70) {
  valid <- !is.na(y_true) & !is.na(y_pred) & nzchar(as.character(y_true)) & nzchar(as.character(y_pred))
  if (sum(valid) < 4L) return(max(0.05, default_kappa))
  yt <- as.character(y_true[valid]); yp <- as.character(y_pred[valid]); lvls <- sort(unique(c(yt, yp)))
  if (length(lvls) < 2L || length(unique(yt)) < 2L) return(max(0.05, min(0.95, mean(yt == yp))))
  tab <- table(factor(yt, levels = lvls), factor(yp, levels = lvls)); n_tot <- sum(tab)
  p_o <- sum(diag(tab)) / n_tot; p_e <- sum((rowSums(tab) / n_tot) * (colSums(tab) / n_tot))
  if (!is.finite(p_e) || (1.0 - p_e) <= 1e-8) return(max(0.05, p_o))
  kv <- (p_o - p_e) / (1.0 - p_e); if (!is.finite(kv)) return(max(0.05, default_kappa))
  max(0.05, min(0.99, round(kv, 4)))
}

.compute_oof_class_precision_map <- function(y_true_vec, y_pred_vec, prior_alpha = 0.5) {
  ok <- !is.na(y_true_vec) & !is.na(y_pred_vec) & nzchar(y_true_vec) & nzchar(y_pred_vec) & (y_true_vec != "Unknown") & (y_pred_vec != "Unknown")
  if (sum(ok) < 6L) return(numeric(0))
  yt <- as.character(y_true_vec[ok]); yp <- as.character(y_pred_vec[ok]); true_taxa <- sort(unique(yt)); all_taxa <- sort(unique(c(true_taxa, yp)))
  if (length(true_taxa) < 2L || length(all_taxa) < 2L) return(numeric(0))
  K_cols <- length(all_taxa); C_mat <- matrix(0.0, nrow = length(true_taxa), ncol = K_cols, dimnames = list(true_taxa, all_taxa))
  for (k_tx in true_taxa) { idx_k <- which(yt == k_tx); n_k <- length(idx_k); for (j_tx in all_taxa) C_mat[k_tx, j_tx] <- (sum(yp[idx_k] == j_tx) + prior_alpha / K_cols) / (n_k + prior_alpha) }
  prec_out <- stats::setNames(rep(1.0, K_cols), all_taxa)
  for (j_tx in all_taxa) {
    cs <- sum(C_mat[, j_tx]); dv <- if (j_tx %in% rownames(C_mat)) C_mat[j_tx, j_tx] else (prior_alpha / K_cols)
    prec_out[[j_tx]] <- min(1.0, max(0.05, round(if (is.finite(cs) && cs > 1e-12) dv / cs else 1.0, 4)))
  }
  prec_out
}

.resolve_step07_anthesis_flag <- function(df, n = nrow(df)) {
  if ("is_active_anthesis" %in% names(df)) { v <- as.logical(df$is_active_anthesis); v[is.na(v)] <- TRUE; v } else resolve_active_anthesis_mask(df)
}

compute_modality_reliability_weights <- function(df) {
  n <- nrow(df); if (n == 0L) return(c(morph = 0.80, edaphic = 0.65, pheno = 0.55))
  sub_a <- df[.resolve_step07_anchor_mask(df), , drop = FALSE]; na_s <- nrow(sub_a); y_true <- as.character(sub_a$species_standardized)
  n_true_taxa <- length(unique(y_true[!is.na(y_true) & nzchar(y_true) & y_true != "Unknown"]))
  m_pure <- if ("morph_tier_source" %in% names(sub_a)) (is.na(sub_a$morph_tier_source) | sub_a$morph_tier_source == "Basal_Rosette") else rep(TRUE, na_s)
  y_morph <- .chr_col(sub_a, "cda_predicted_taxon", na_s)
  def_morph <- if ("loocv_accuracy" %in% names(df) && any(is.finite(as.numeric(df$loocv_accuracy)))) stats::median(as.numeric(df$loocv_accuracy), na.rm = TRUE) else 0.80
  y_edaph <- if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else .chr_col(sub_a, "rf_predicted_taxon", na_s)
  cv_med <- if ("spatial_block_cv_acc" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$spatial_block_cv_acc))))) stats::median(suppressWarnings(as.numeric(df$spatial_block_cv_acc)), na.rm = TRUE) else 0.65
  if (n_true_taxa < 2L) {
    out_m_a <- if ("morphometric_outlier_flag" %in% names(sub_a)) as.logical(sub_a$morphometric_outlier_flag[m_pure]) else rep(FALSE, sum(m_pure)); out_m_a[is.na(out_m_a)] <- FALSE
    w_morph <- round(min(0.80, max(0.30, (if (length(out_m_a) > 0L) mean(!out_m_a) else 0.80) * 0.65)), 4)
    eco_reg_a <- if ("ecotype_regime" %in% names(sub_a)) as.character(sub_a$ecotype_regime) else .chr_col(sub_a, "analysis_group", na_s)
    rf_reg_a <- if ("rf_predicted_taxon" %in% names(sub_a)) as.character(sub_a$rf_predicted_taxon) else .chr_col(sub_a, "edaphic_best_fit_taxon", na_s)
    u_eco_reg <- unique(eco_reg_a[!is.na(eco_reg_a) & nzchar(eco_reg_a) & eco_reg_a != "Unknown"])
    w_edaph <- if (length(u_eco_reg) >= 2L && any(rf_reg_a %in% u_eco_reg, na.rm = TRUE)) max(0.15, min(0.85, compute_cohens_kappa(eco_reg_a, rf_reg_a, default_kappa = 0.50))) else round(min(0.75, max(0.20, 2.0 * (cv_med - 0.50))), 4)
  } else {
    w_morph <- max(0.05, compute_cohens_kappa(y_true[m_pure], y_morph[m_pure], default_kappa = def_morph))
    w_edaph <- max(0.05, compute_cohens_kappa(y_true, y_edaph, default_kappa = cv_med))
  }
  anth_a <- .resolve_step07_anthesis_flag(sub_a, na_s); hop_col_a <- intersect(c("hopkins_adjusted_doy", "doy"), names(sub_a))[1L]
  doy_a <- if (!is.na(hop_col_a)) .num_col(sub_a, hop_col_a, na_s) else rep(NA_real_, na_s)
  ddoy_a <- if ("delta_doy_resid" %in% names(sub_a)) .num_col(sub_a, "delta_doy_resid", na_s) else .num_col(sub_a, "delta_doy", na_s)
  ok_ph_cls <- which(anth_a & is.finite(doy_a) & doy_a >= 1.0 & doy_a <= 366.0 & !is.na(y_true) & nzchar(y_true) & y_true != "Unknown")
  tx_ph_valid <- if (length(ok_ph_cls) >= 4L) names(which(table(y_true[ok_ph_cls]) >= 2L)) else character(0)
  ok_ph_cls <- ok_ph_cls[y_true[ok_ph_cls] %in% tx_ph_valid]
  w_pheno <- if (length(tx_ph_valid) >= 2L && length(ok_ph_cls) >= 4L) {
    th_ph <- 2.0 * pi * (doy_a[ok_ph_cls] / 365.25); yt_ph <- y_true[ok_ph_cls]; yp_ph <- character(length(ok_ph_cls))
    for (ii in seq_along(ok_ph_cls)) {
      best_d <- Inf; best_tx <- tx_ph_valid[1L]
      for (tx_k in tx_ph_valid) {
        jj <- which(yt_ph == tx_k & seq_along(ok_ph_cls) != ii); if (length(jj) == 0L) jj <- which(yt_ph == tx_k)
        mu_k <- atan2(mean(sin(th_ph[jj])), mean(cos(th_ph[jj]))); r_k <- max(1e-4, min(0.9995, sqrt(mean(sin(th_ph[jj]))^2 + mean(cos(th_ph[jj]))^2)))
        d_z <- abs(atan2(sin(th_ph[ii] - mu_k), cos(th_ph[ii] - mu_k))) / max(2.0 * pi * (6.0 / 365.25), sqrt(-2.0 * log(r_k)))
        if (d_z < best_d) { best_d <- d_z; best_tx <- tx_k }
      }
      yp_ph[ii] <- best_tx
    }
    max(0.05, min(0.95, compute_cohens_kappa(yt_ph, yp_ph, default_kappa = 0.15)))
  } else if (sum(is.finite(ddoy_a)) >= 4L) 0.20 else 0.15
  w_geno <- if (exists("compute_genomic_kappa_weight", mode = "function")) compute_genomic_kappa_weight(df) else NA_real_
  if (is.finite(w_geno)) c(morph = w_morph, edaphic = w_edaph, pheno = w_pheno, geno = w_geno) else c(morph = w_morph, edaphic = w_edaph, pheno = w_pheno)
}

compute_probabilistic_consensus_scores <- function(df) {
  n <- nrow(df); if (n == 0L) return(df)
  prop_sum_attr <- attr(df, "genomic_propagation_summary"); w_vec <- compute_modality_reliability_weights(df)
  df$w_morph_kappa <- w_vec[["morph"]]; df$w_edaphic_kappa <- w_vec[["edaphic"]]; df$w_pheno_kappa <- w_vec[["pheno"]]
  has_geno_row <- .has_step07_geno_row_mask(df); w_geno_base <- if ("geno" %in% names(w_vec) && is.finite(w_vec[["geno"]])) w_vec[["geno"]] else NA_real_
  sub_a <- df[.resolve_step07_anchor_mask(df), , drop = FALSE]
  y_true_a <- .chr_col(sub_a, "species_standardized", nrow(sub_a))
  m_pure_a <- if ("morph_tier_source" %in% names(sub_a)) (is.na(sub_a$morph_tier_source) | sub_a$morph_tier_source == "Basal_Rosette") else rep(TRUE, nrow(sub_a))
  y_morph_a <- .chr_col(sub_a, "cda_predicted_taxon", nrow(sub_a))
  y_edaph_a <- if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else .chr_col(sub_a, "rf_predicted_taxon", nrow(sub_a))
  prec_morph <- .compute_oof_class_precision_map(y_true_a[m_pure_a], y_morph_a[m_pure_a]); prec_edaph <- .compute_oof_class_precision_map(y_true_a, y_edaph_a)
  lookup_prec <- function(pmap, pred_tx) if (length(pmap) == 0L || is.na(pred_tx) || !nzchar(pred_tx) || !(pred_tx %in% names(pmap))) 1.0 else as.numeric(pmap[[pred_tx]])

  p_conf_out <- numeric(n); p_score_out <- numeric(n); n_mod_out <- integer(n); df_nms_c <- names(df)
  w_morph_sx <- numeric(n); w_edaph_sx <- numeric(n); w_pheno_sx <- numeric(n); w_geno_sx <- numeric(n)
  given_sp_cv <- as.character(df$species_standardized)
  morph_sp_cv <- .chr_col(df, "cda_predicted_taxon", n); morph_sp_cv[!is.na(morph_sp_cv) & !nzchar(morph_sp_cv)] <- NA_character_
  cda_p_cv <- .num_col(df, "cda_posterior_prob", n); is_mo_cv <- .lgl_col(df, "morphometric_outlier_flag", n); has_geo_cv <- is.finite(.num_col(df, "decimalLatitude", n))
  edaph_sp_cv <- if ("edaphic_best_fit_taxon" %in% df_nms_c) as.character(df$edaphic_best_fit_taxon) else .chr_col(df, "rf_predicted_taxon", n)
  edaph_sp_cv[!is.na(edaph_sp_cv) & !nzchar(edaph_sp_cv)] <- NA_character_
  rf_p_cv <- .num_col(df, "rf_predicted_prob", n); rf_g_cv <- .num_col(df, "rf_given_taxon_prob", n)
  ez_cv <- .num_col(df, "edaphic_z_score", n); tez_cv <- .num_col(df, "type_edaphic_z_score", n); ddoy_r_cv <- .num_col(df, "delta_doy_resid", n)
  ded_cv <- if ("edaphic_niche_overlap_d" %in% df_nms_c) .num_col(df, "edaphic_niche_overlap_d", n) else .num_col(df, "edaphic_schoener_d", n)
  is_anth_cv <- if ("is_active_anthesis" %in% df_nms_c) .lgl_col(df, "is_active_anthesis", n) else rep(TRUE, n)
  w_morph_base <- w_vec[["morph"]]; w_edaph_base <- w_vec[["edaphic"]]; w_pheno_base <- w_vec[["pheno"]]
  is_single_sp_cons <- length(unique(given_sp_cv[!is.na(given_sp_cv) & nzchar(given_sp_cv) & given_sp_cv != "Unknown"])) < 2L
  eco_reg_cv <- if ("ecotype_regime" %in% df_nms_c) as.character(df$ecotype_regime) else .chr_col(df, "analysis_group", n)
  sd2_cv <- .num_col(df, "self_morph_mahal_d2", n); sdf_cv <- .num_col(df, "self_morph_df", n, 6.0); sdf_cv[!is.finite(sdf_cv) | sdf_cv < 1.0] <- 6.0
  geno_sp_cv <- .chr_col(df, "geno_predicted_taxon", n); geno_p_cv <- .num_col(df, "geno_posterior_prob", n); geno_qmax_cv <- .num_col(df, "q_max", n); deme_w_cv <- .num_col(df, "geno_deme_weight_scale", n)
  geno_tier_cv <- if ("geno_linkage_tier" %in% df_nms_c) as.character(df$geno_linkage_tier) else .chr_col(df, "linkage_tier", n)
  geno_adm_cv <- if ("is_geno_admixed" %in% df_nms_c) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
  geno_qry_cv <- if ("is_putative_novel_or_hybrid_query" %in% df_nms_c) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  geno_link_scale <- ifelse(geno_tier_cv == "G1_Direct_Voucher", 1.00, ifelse(geno_tier_cv == "G1_Exsiccatae_Isoduplicate", 0.95,
                     ifelse(geno_tier_cv == "G2_Population_Voucher_Proxy", 0.90, ifelse(geno_tier_cv == "G2_Conditional_Deme_Hull", ifelse(is.finite(deme_w_cv) & deme_w_cv > 0.0, pmin(0.75, pmax(0.45, deme_w_cv)), 0.75), 0.00))))
  geno_link_scale[is.na(geno_link_scale)] <- 0.0; run_has_geno <- any(has_geno_row) && is.finite(w_geno_base)
  df$w_geno_kappa <- if (run_has_geno) ifelse(has_geno_row, round(w_geno_base * geno_link_scale, 4), 0.0) else NA_real_
  t0_abstain_cv <- if ("tier0_abstention_flag" %in% df_nms_c) suppressWarnings(as.integer(df$tier0_abstention_flag)) else rep(NA_integer_, n)
  t0_cal_tax_cv <- trimws(.chr_col(df, "tier0_calibrated_taxon", n)); t0_cal_prob_cv <- .num_col(df, "tier0_calibrated_prob", n)
  lopo_raw <- if (!is.null(prop_sum_attr) && !is.null(prop_sum_attr$lopo_cv_kappa)) suppressWarnings(as.numeric(prop_sum_attr$lopo_cv_kappa)) else NA_real_
  lopo_cv_kappa <- min(1.0, max(0.0, if (is.finite(lopo_raw)) lopo_raw else 0.60))
  known_species_pool <- unique(given_sp_cv[!is.na(given_sp_cv) & nzchar(trimws(given_sp_cv)) & given_sp_cv != "Unknown"])

  for (i in seq_len(n)) {
    given_sp <- given_sp_cv[i]; morph_sp <- morph_sp_cv[i]; cda_p <- cda_p_cv[i]; is_mo <- is_mo_cv[i]
    has_geo <- has_geo_cv[i]; edaph_sp <- edaph_sp_cv[i]; rf_p <- rf_p_cv[i]; rf_g <- rf_g_cv[i]
    ez_i <- ez_cv[i]; tez_i <- tez_cv[i]; ded_i <- ded_cv[i]; ddoy_r <- ddoy_r_cv[i]; is_anth <- is_anth_cv[i]
    morph_is_subcluster <- !is.na(morph_sp) && grepl("^Morphotype_K", morph_sp); edaph_is_ecotype <- !is.na(edaph_sp) && grepl("^Ecotype_", edaph_sp)

    p_m <- if (is_single_sp_cons || morph_is_subcluster) {
      if (is_mo) 0.72 else if (is.finite(sd2_cv[i])) min(0.85, max(0.05, stats::pchisq(sd2_cv[i], df = sdf_cv[i]) * 0.45)) else if (!is.na(morph_sp) || "gmm_cluster" %in% df_nms_c) 0.08 else NA_real_
    } else if (!is.na(morph_sp)) {
      if (morph_sp != given_sp) max(if (is.finite(cda_p)) cda_p else 0.75, if (is_mo) 0.80 else 0.55) else if (is_mo) 0.65 else if (is.finite(cda_p)) max(0.0, 1.0 - cda_p) else 0.10
    } else NA_real_
    p_e <- if (is_single_sp_cons || edaph_is_ecotype) {
      if (has_geo && (is.finite(ez_i) || is.finite(tez_i))) { z_u <- if (is.finite(ez_i)) ez_i else tez_i; min(0.92, max(0.04, stats::pchisq(z_u^2, df = 1) * 0.55)) }
      else if (edaph_is_ecotype && !is.na(eco_reg_cv[i]) && nzchar(eco_reg_cv[i])) { if (edaph_sp == eco_reg_cv[i]) 0.08 else 0.35 } else NA_real_
    } else if (!is.na(edaph_sp)) {
      if (edaph_sp != given_sp) { if (is.finite(rf_p)) max(0.55, min(0.98, rf_p)) else if (is.finite(ez_i) && ez_i > 2.0) min(0.95, 0.55 + 0.10 * (ez_i - 2.0)) else if (is.finite(ded_i) && ded_i < 0.35) 0.80 else 0.65 }
      else { if (is.finite(rf_g)) max(0.05, min(0.85, 1.0 - rf_g)) else if (is.finite(rf_p)) max(0.05, min(0.85, 1.0 - rf_p)) else if (is.finite(ez_i) && ez_i > 2.5) min(0.75, 0.25 * ez_i) else 0.10 }
    } else if (has_geo && (is.finite(ez_i) || is.finite(tez_i))) { z_u <- if (is.finite(ez_i)) ez_i else tez_i; min(0.95, max(0.05, stats::pchisq(z_u^2, df = 1) * 0.85)) } else NA_real_
    p_p <- if (is_anth && is.finite(ddoy_r)) min(1.0, abs(ddoy_r) / 30.0) else NA_real_
    p_g <- if (has_geno_row[i] && is.finite(w_geno_base) && geno_link_scale[i] > 0.0) {
      g_sp_i <- geno_sp_cv[i]; g_p_i <- geno_p_cv[i]; g_qm_i <- geno_qmax_cv[i]
      if (geno_adm_cv[i] || geno_qry_cv[i]) { if (is.finite(g_qm_i)) max(if (is_single_sp_cons) 0.35 else 0.40, min(0.92, 1.0 - g_qm_i)) else if (is_single_sp_cons) 0.72 else 0.65 }
      else if (!is.na(g_sp_i) && g_sp_i != given_sp) { if (is.finite(g_p_i)) max(0.55, min(0.99, g_p_i)) else if (is_single_sp_cons) 0.08 else 0.80 }
      else if (is_single_sp_cons) { if (is.finite(g_qm_i) && g_qm_i < 0.75) 0.55 else 0.08 }
      else { if (is.finite(g_p_i)) max(0.02, min(0.65, 1.0 - g_p_i)) else 0.08 }
    } else NA_real_
    has_t0_soft_i <- (geno_link_scale[i] == 0.0) && !is.na(t0_abstain_cv[i]) && (t0_abstain_cv[i] == 0L) &&
      is.finite(t0_cal_prob_cv[i]) && (t0_cal_prob_cv[i] >= 0.80) && !is.na(t0_cal_tax_cv[i]) && nzchar(t0_cal_tax_cv[i]) && (t0_cal_tax_cv[i] != "Sensu_Lato_DNA_Required")
    w_t0_i <- if (has_t0_soft_i && lopo_cv_kappa > 0.0) 0.50 * lopo_cv_kappa else NA_real_
    p_t0_i <- if (is.finite(w_t0_i)) {
      s_t0 <- if ((grepl(" ", t0_cal_tax_cv[i]) || (t0_cal_tax_cv[i] %in% known_species_pool)) && t0_cal_tax_cv[i] != given_sp) (1.0 - t0_cal_prob_cv[i]) else t0_cal_prob_cv[i]
      min(1.0, max(0.0, 1.0 - s_t0))
    } else NA_real_
    w_i_vec <- c(morph = w_morph_base * (0.5 + 0.5 * lookup_prec(prec_morph, morph_sp)), edaphic = w_edaph_base * (0.5 + 0.5 * lookup_prec(prec_edaph, edaph_sp)),
                 pheno = w_pheno_base, geno = if (has_geno_row[i] && is.finite(w_geno_base)) w_geno_base * geno_link_scale[i] else NA_real_, t0_soft = w_t0_i)
    p_vec <- c(morph = p_m, edaphic = p_e, pheno = p_p, geno = p_g, t0_soft = p_t0_i)
    avail <- !is.na(p_vec) & is.finite(p_vec) & !is.na(w_i_vec) & is.finite(w_i_vec); n_mod_out[i] <- sum(avail)
    p_conf_i <- if (any(avail)) sum(w_i_vec[avail] * p_vec[avail]) / sum(w_i_vec[avail]) else 0.0
    w4_vec <- pmax(0.0, w_i_vec[c("morph", "edaphic", "pheno", "geno")]); av4 <- avail[c("morph", "edaphic", "pheno", "geno")]; w4_sx <- rep(0.0, 4L)
    if (any(av4) && sum(w4_vec[av4]) > 0.0) { w4_rnd <- round(w4_vec[av4] / sum(w4_vec[av4]), 6L); im <- which.max(w4_rnd); w4_rnd[im] <- round(w4_rnd[im] + (1.0 - sum(w4_rnd)), 6L); w4_sx[av4] <- w4_rnd }
    w_morph_sx[i] <- w4_sx[1L]; w_edaph_sx[i] <- w4_sx[2L]; w_pheno_sx[i] <- w4_sx[3L]; w_geno_sx[i] <- w4_sx[4L]
    p_conf_out[i] <- round(min(1.0, max(0.0, p_conf_i)), 4); p_score_out[i] <- round(1.0 - p_conf_out[i], 4)
  }
  df$w_morph_simplex <- w_morph_sx; df$w_edaphic_simplex <- w_edaph_sx; df$w_pheno_simplex <- w_pheno_sx; df$w_geno_simplex <- w_geno_sx
  df$n_available_modalities <- n_mod_out; df$probabilistic_consensus_conflict <- p_conf_out; df$probabilistic_consensus_score <- p_score_out
  if (!is.null(prop_sum_attr)) attr(df, "genomic_propagation_summary") <- prop_sum_attr; df
}
