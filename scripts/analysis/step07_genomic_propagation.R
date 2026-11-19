# Step 07 Genomic Propagation, Tier_0 Calibration, Deme Priors & gdi Governance
for (.dp_p in c("scripts/analysis/step07_deme_propagation_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step07_deme_propagation_utils.R"))) {
  if (!exists(".eval_mcl_cluster_gdi_prop", mode = "function") && !is.na(.dp_p) && nzchar(.dp_p) && file.exists(.dp_p)) { source(.dp_p); break }
}

# Step A: Propagate G1_Direct_Voucher across duplicate exsiccatae sheets
propagate_exsiccatae_isoduplicates <- function(df, ga_df = NULL) {
  df <- .init_propagation_cols(df, reset = FALSE)
  n <- nrow(df); if (n < 2L) return(df)
  tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), ""); qc <- .pick_bool_prop(df, "geno_quality_pass", FALSE)
  adm <- .pick_bool_prop(df, "is_geno_admixed", FALSE); pred <- trimws(.pick_chr_prop(df, "geno_predicted_taxon", "")); prob <- .pick_num_prop(df, "geno_posterior_prob")

  if (!is.null(ga_df) && is.data.frame(ga_df) && nrow(ga_df) > 0L && "catalogNumber" %in% names(df) && "catalogNumber" %in% names(ga_df)) {
    ga_keep <- (.pick_chr_prop(ga_df, c("geno_linkage_tier", "linkage_tier"), "") == "G1_Direct_Voucher") & .pick_bool_prop(ga_df, c("geno_quality_pass", "qc_pass"), FALSE)
    if (any(ga_keep)) {
      ga_sub <- ga_df[ga_keep, , drop = FALSE]; m_idx <- match(as.character(df$catalogNumber), as.character(ga_sub$catalogNumber))
      fill_idx <- which(!nzchar(tier) & !is.na(m_idx))
      if (length(fill_idx) > 0L) for (col in intersect(c("geno_predicted_taxon", "geno_posterior_prob", "geno_quality_pass", "is_geno_admixed", "genetic_sample_id"), names(ga_sub))) if (!(col %in% names(df))) df[[col]] <- rep(NA, n)
      for (r in fill_idx) {
        g <- m_idx[r]; tier[r] <- "G1_Direct_Voucher"; qc[r] <- TRUE
        adm[r] <- .pick_bool_prop(ga_sub[g, , drop = FALSE], "is_geno_admixed", FALSE)[1L]
        pred[r] <- trimws(.pick_chr_prop(ga_sub[g, , drop = FALSE], "geno_predicted_taxon", "")[1L])
        prob[r] <- .pick_num_prop(ga_sub[g, , drop = FALSE], "geno_posterior_prob")[1L]
        for (col in setdiff(intersect(names(ga_sub), names(df)), "catalogNumber")) {
          if (is.na(df[[col]][r]) || (is.character(df[[col]][r]) && !nzchar(trimws(df[[col]][r])))) df[[col]][r] <- ga_sub[[col]][g]
        }
      }
    }
  }

  ex_k <- trimws(.pick_chr_prop(df, "exsiccatae_key", "")); sur_c <- trimws(.pick_chr_prop(df, "collector_surname", ""))
  sur <- .extract_collector_surname(ifelse(nzchar(sur_c), sur_c, .pick_chr_prop(df, "recordedBy", "")))
  recn <- .norm_record_number(.pick_chr_prop(df, "recordNumber", "")); evdt <- tolower(substr(trimws(.pick_chr_prop(df, "eventDate", "")), 1L, 10L))
  gkey <- ifelse(nzchar(ex_k), ex_k, ifelse(nzchar(sur) & nzchar(recn), paste(sur, recn, evdt, sep = "|"), ""))
  nonempty_k <- gkey[nzchar(gkey)]; if (length(nonempty_k) < 2L) return(df)
  ktab <- table(nonempty_k); is_donor <- (tier == "G1_Direct_Voucher") & qc & !adm & nzchar(pred) & (pred != "Unknown") & nzchar(gkey)
  if (!any(is_donor)) return(df)

  is_type <- .is_type_prop(df); cat_vec <- if ("catalogNumber" %in% names(df)) as.character(df$catalogNumber) else as.character(seq_len(n))
  n_pl <- if ("n_plant_individuals" %in% names(df)) suppressWarnings(as.integer(df$n_plant_individuals)) else rep(1L, n)
  is_multi <- (as.integer(table(cat_vec)[cat_vec]) > 1L) | (!is.na(n_pl) & n_pl > 1L); is_mo <- .pick_bool_prop(df, c("morphometric_outlier_flag", "is_morph_outlier"), FALSE)
  pool <- .compute_morph_pool_prop(df)
  copy_cols <- intersect(c(
    "genetic_sample_id", "geno_predicted_taxon", "geno_posterior_prob", "geno_margin", "q_max", "q_admixture_json",
    "genetic_cluster", "is_geno_admixed", "is_putative_novel_or_hybrid_query", "geno_quality_pass", "geno_parent_1", "geno_parent_2",
    "nearest_tree_sister_taxon", "patristic_dist_to_dubia", "replicate_pair_concordance", "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl",
    "plastid_clade", "its_ribotype", "cytonuclear_discordant", "hybrid_index_h", "interspecific_ancestry_q12", "hybrid_generation_class",
    "coancestry_cluster", "f_branch_max", "f_branch_donor_taxon", "snaq_gamma_minor", "coalescent_gdi",
    "bpp_introgression_bf", "mmrr_beta_ibd", "mmrr_beta_ibe", "mmrr_p_ibe", "public_backbone_clade", "public_phytop_ils_index",
    "public_phytop_ih_index", "public_plastome_match", "public_its_match", "waxy_lineage_match"
  ), names(df))

  for (k in unique(gkey[is_donor])) {
    if (is.na(ktab[k]) || ktab[k] < 2L) next
    d_idx <- which(is_donor & gkey == k); if (length(unique(pred[d_idx])) != 1L) next
    d <- d_idx[which.max(ifelse(is.finite(prob[d_idx]), prob[d_idx], -Inf))]
    for (j in setdiff(which(gkey == k), d_idx)) {
      if ((tier[j] == "G1_Direct_Voucher") || is_type[j] || (is_multi[j] && is_mo[j])) next
      if (!is.null(pool) && length(pool$pc_cols) >= 2L && pool$is_basal[d] && pool$is_basal[j] && all(is.finite(pool$X[d, ])) && all(is.finite(pool$X[j, ]))) {
        dv <- pool$X[j, ] - pool$X[d, ]
        if (is.finite(d2 <- as.numeric(crossprod(dv, pool$S_inv %*% dv))) && d2 > stats::qchisq(0.99, df = length(pool$pc_cols))) next
      }
      for (col in copy_cols) df[[col]][j] <- df[[col]][d]
      if (!("geno_linkage_tier" %in% names(df))) df$geno_linkage_tier <- rep(NA_character_, n)
      df$geno_linkage_tier[j] <- "G1_Exsiccatae_Isoduplicate"
      if ("linkage_tier" %in% names(df)) df$linkage_tier[j] <- "G1_Exsiccatae_Isoduplicate"
      df$genomic_propagation_tier[j] <- "G1_Exsiccatae_Isoduplicate"
    }
  }
  df
}

# Step B: Resolve Tier_0 DNA-verified anchors and anchor evictions
resolve_tier0_anchor_overrides <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(list(is_tier0_dna_anchor = logical(0), evict_from_given_anchor = logical(0), rectified_taxon = character(0), tier0_weight = numeric(0)))
  is_type <- .is_type_prop(df); tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), "")
  qc <- .pick_bool_prop(df, "geno_quality_pass", FALSE); adm <- .pick_bool_prop(df, "is_geno_admixed", FALSE)
  nov <- .pick_bool_prop(df, "is_putative_novel_or_hybrid_query", FALSE); pred <- trimws(.pick_chr_prop(df, "geno_predicted_taxon", ""))
  sp <- trimws(.pick_chr_prop(df, "species_standardized", "")); prob <- .pick_num_prop(df, "geno_posterior_prob")
  margin <- if ("geno_margin" %in% names(df)) suppressWarnings(as.numeric(df$geno_margin)) else rep(NA_real_, n)
  q12 <- .pick_num_prop(df, "interspecific_ancestry_q12")
  is_t0 <- !is_type & (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc & !adm & !nov &
    is.finite(prob) & (prob >= 0.85) & (is.na(margin) | margin >= 0.40) & nzchar(pred) & (pred != "Unknown")
  evict <- !is_type & (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc &
    ((nzchar(pred) & (pred != "Unknown") & (pred != sp) & is.finite(prob) & (prob >= 0.85)) | (is.finite(q12) & (q12 >= 0.35)) | (adm & !is.finite(q12)))
  list(is_tier0_dna_anchor = is_t0, evict_from_given_anchor = evict, rectified_taxon = ifelse(is_t0, pred, sp), tier0_weight = ifelse(is_t0, 4.0, NA_real_))
}

get_tier0_rectified_species_and_evictions <- function(df) resolve_tier0_anchor_overrides(df)

# Step C: Conditional spatial deme priors (conflict-only vetoes + soft reliability attenuation)
propagate_conditional_deme_priors <- function(df, ga_df = NULL, morph_X = NULL, morph_S_inv = NULL) {
  df <- .init_propagation_cols(df, reset = FALSE)
  n <- nrow(df); if (n == 0L) return(df)
  anc_df <- .extract_deme_anchors_prop(df, ga_df = ga_df); if (nrow(anc_df) == 0L) return(df)
  if (is.null(morph_X) || is.null(morph_S_inv)) {
    pool <- .compute_morph_pool_prop(df); if (is.null(pool)) return(df)
    morph_X <- pool$X; morph_S_inv <- pool$S_inv; is_basal <- pool$is_basal
  } else is_basal <- .is_basal_mask_prop(df)
  p <- ncol(morph_X); if (p < 2L) return(df)

  t0 <- resolve_tier0_anchor_overrides(df); mu_taxon <- .build_taxon_centroids(df, morph_X, is_basal, t0)
  for (tx in unique(anc_df$taxon)) {
    if (is.null(mu_taxon[[tx]])) {
      rws <- anc_df$df_row[anc_df$taxon == tx & !is.na(anc_df$df_row)]
      rws <- rws[is_basal[rws] & stats::complete.cases(morph_X[rws, , drop = FALSE])]
      if (length(rws) > 0L) mu_taxon[[tx]] <- colMeans(morph_X[rws, , drop = FALSE])
    }
  }
  tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), ""); sp <- trimws(.pick_chr_prop(df, "species_standardized", ""))
  lat <- .pick_num_prop(df, "decimalLatitude"); lon <- .pick_num_prop(df, "decimalLongitude"); az <- .pick_num_prop(df, "is_azonal_lithology")
  ph <- .pick_num_prop(df, "fine_soil_ph"); ed_z <- .pick_num_prop(df, c("type_edaphic_z_score", "edaphic_z_score"))
  unc <- .pick_num_prop(df, "coordinateUncertaintyInMeters"); doy <- .pick_num_prop(df, "hopkins_adjusted_doy")
  rstat <- tolower(.pick_chr_prop(df, c("reproductive_trait_status", "capitulum_phenology_state"), ""))
  is_anth <- if ("is_active_anthesis" %in% names(df)) .safe_bool_prop(df$is_active_anthesis) else (is.finite(doy) & !grepl("vegetative|bud_only|post_fruiting", rstat))
  tau_eco_v <- if ("tau_ecotype" %in% names(df)) suppressWarnings(as.numeric(df$tau_ecotype)) else if (!is.null(attr(df, "tau_ecotype"))) suppressWarnings(as.numeric(attr(df, "tau_ecotype"))) else rep(NA_real_, n)
  if (length(tau_eco_v) == 1L) tau_eco_v <- rep(tau_eco_v, n)

  pure_ev <- t0$evict_from_given_anchor & !t0$is_tier0_dna_anchor; mu_doy <- list(); ph_ref <- list()
  for (tx in unique(c(names(mu_taxon), anc_df$taxon))) {
    a_idx <- which(t0$rectified_taxon == tx & !pure_ev & is_anth & is.finite(doy))
    if (length(a_idx) >= 2L) { ang <- doy[a_idx] * 2.0 * pi / 365.25; mu_doy[[tx]] <- ((atan2(mean(sin(ang)), mean(cos(ang))) * 365.25 / (2.0 * pi)) %% 365.25) }
    p_idx <- which(t0$rectified_taxon == tx & !pure_ev & is.finite(ph))
    if (length(p_idx) >= 2L) ph_ref[[tx]] <- c(mean(ph[p_idx]), max(stats::sd(ph[p_idx]), 0.4))
  }

  is_cent_raw <- .pick_bool_prop(df, "is_centroid_or_inst_coord", FALSE)
  cands <- which(!(tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy")) & !.is_type_prop(df) & !is_cent_raw & is.finite(lat) & is.finite(lon))
  df_geo_verified <- .resolve_verified_sympatry_rows_prop(df, t0, lat, lon, pure_ev); crit_d2 <- stats::qchisq(0.95, df = p)
  qc_num <- is.numeric(df$geno_quality_pass) || is.integer(df$geno_quality_pass); adm_num <- is.numeric(df$is_geno_admixed) || is.integer(df$is_geno_admixed)
  rec_ploidy <- .extract_known_ploidy_vec(df); is_multi_tx <- .pick_num_prop(df, "is_multi_cytotype_taxon")
  multi_sp_set <- unique(sp[!is.na(is_multi_tx) & is_multi_tx == 1L & nzchar(sp)]); row_typed_pl <- rec_ploidy
  has_dr <- !is.na(anc_df$df_row) & !is.na(anc_df$ploidy)
  if (any(has_dr)) { dr <- anc_df$df_row[has_dr]; m <- is.na(row_typed_pl[dr]); if (any(m)) row_typed_pl[dr[m]] <- anc_df$ploidy[has_dr][m] }
  df_typed_ll <- which(!is.na(row_typed_pl) & !is_cent_raw & is.finite(lat) & is.finite(lon))

  for (i in cands) {
    dists <- .haversine_km_vec(lat[i], lon[i], anc_df$lat, anc_df$lon); in_2km <- which(dists <= 2.0)
    if (length(in_2km) == 0L || (is.finite(unc[i]) && unc[i] > 2000.0)) next
    min_d <- min(dists[in_2km]); soft_mult <- if ((is.finite(unc[i]) && unc[i] > 1000.0) || (is.na(unc[i]) && min_d > 1.0)) 0.85 else 1.0
    np <- in_2km[dists[in_2km] <= min_d + 1e-6]; near_a <- np[1L]
    if (length(np) > 1L && is_basal[i] && all(is.finite(morph_X[i, ]))) {
      cd2 <- vapply(np, function(a) if (is.null(mu_taxon[[anc_df$taxon[a]]])) Inf else { dv0 <- morph_X[i, ] - mu_taxon[[anc_df$taxon[a]]]; as.numeric(crossprod(dv0, morph_S_inv %*% dv0)) }, numeric(1L))
      if (any(is.finite(cd2))) near_a <- np[which.min(cd2)]
    }
    cand_sp <- anc_df$taxon[near_a]
    if (is.finite(az[i]) && is.finite(anc_df$azonal[near_a]) && as.integer(az[i]) != as.integer(anc_df$azonal[near_a])) next
    if (!is.na(anc_df$ploidy[near_a]) && !is.na(rec_ploidy[i]) && anc_df$ploidy[near_a] != rec_ploidy[i]) next
    if ((!is.na(is_multi_tx[i]) && is_multi_tx[i] == 1L) || (cand_sp %in% multi_sp_set)) {
      df_2km_pl <- if (length(df_typed_ll) > 0L) df_typed_ll[.haversine_km_vec(lat[i], lon[i], lat[df_typed_ll], lon[df_typed_ll]) <= 2.0] else integer(0)
      deme_pls <- unique(stats::na.omit(c(anc_df$ploidy[in_2km], row_typed_pl[df_2km_pl])))
      if (length(deme_pls) > 1L || (length(deme_pls) == 1L && !is.na(rec_ploidy[i]) && rec_ploidy[i] != deme_pls[1L])) next
      if (length(deme_pls) == 0L || is.na(rec_ploidy[i])) soft_mult <- soft_mult * 0.85
    }
    dr_near <- anc_df$df_row[near_a]; ph_donor <- if (!is.na(dr_near) && is.finite(ph[dr_near])) ph[dr_near] else anc_df$ph[near_a]
    ed_gate <- .eval_deme_edaphic_gate_prop(cand_sp, sp[i], ed_z[i], ph[i], ph_donor, ph_ref[[cand_sp]])
    if (!ed_gate$pass) next
    soft_mult <- soft_mult * ed_gate$mult
    if (!is_basal[i] || !all(is.finite(morph_X[i, ])) || is.null(mu_taxon[[cand_sp]])) next
    dv <- morph_X[i, ] - mu_taxon[[cand_sp]]; d2_cand <- as.numeric(crossprod(dv, morph_S_inv %*% dv))
    df$deme_mahal_d2[i] <- d2_cand; if (!is.finite(d2_cand) || d2_cand > crit_d2) next
    doy_donor <- if (!is.na(dr_near) && is_anth[dr_near] && is.finite(doy[dr_near])) doy[dr_near] else anc_df$doy[near_a]
    ph_gate <- .eval_deme_anthesis_gate_prop(is_anth[i], doy[i], doy_donor, mu_doy[[cand_sp]], tau_eco_v[i])
    if (!ph_gate$pass) next
    soft_mult <- soft_mult * ph_gate$mult
    df_2km_ver <- if (length(df_geo_verified) > 0L) t0$rectified_taxon[df_geo_verified[.haversine_km_vec(lat[i], lon[i], lat[df_geo_verified], lon[df_geo_verified]) <= 2.0]] else character(0)
    alt_taxa <- setdiff(unique(c(anc_df$taxon[dists <= 5.0], df_2km_ver)), c(cand_sp, "", "Unknown"))
    if (length(alt_taxa) > 0L) {
      alt_d2 <- vapply(alt_taxa, function(atx) if (is.null(mu_taxon[[atx]])) NA_real_ else { av <- morph_X[i, ] - mu_taxon[[atx]]; as.numeric(crossprod(av, morph_S_inv %*% av)) }, numeric(1L))
      if (!any(is.finite(alt_d2)) || (min(alt_d2, na.rm = TRUE) - d2_cand) < 4.0) next
    }
    if (!("geno_linkage_tier" %in% names(df))) df$geno_linkage_tier <- rep(NA_character_, n)
    if (!("geno_predicted_taxon" %in% names(df))) df$geno_predicted_taxon <- rep(NA_character_, n)
    if (!("geno_quality_pass" %in% names(df))) df$geno_quality_pass <- rep(NA_integer_, n)
    if (!("geno_posterior_prob" %in% names(df))) df$geno_posterior_prob <- rep(NA_real_, n)
    w_scale_i <- round(pmin(0.75, pmax(0.45, 0.75 * soft_mult)), 4L)
    df$geno_deme_weight_scale[i] <- w_scale_i; df$geno_linkage_tier[i] <- "G2_Conditional_Deme_Hull"
    if ("linkage_tier" %in% names(df)) df$linkage_tier[i] <- "G2_Conditional_Deme_Hull"
    df$genomic_propagation_tier[i] <- "G2_Conditional_Deme_Hull"; df$geno_predicted_taxon[i] <- cand_sp
    df$geno_quality_pass[i] <- if (qc_num || !("geno_quality_pass" %in% names(df))) 1L else TRUE
    if ("is_geno_admixed" %in% names(df)) df$is_geno_admixed[i] <- if (adm_num) 0L else FALSE
    base_prob <- pmin(0.90, pmax(0.70, anc_df$prob[near_a] * (1.0 - 0.08 * (min_d / 2.0))))
    df$geno_posterior_prob[i] <- round(pmin(0.90, pmax(0.55, base_prob * (w_scale_i / 0.75))), 4L)
    if ("coalescent_gdi" %in% names(df)) df$coalescent_gdi[i] <- NA_real_
  }
  df
}

# Step D: Tier_0 LOPO-CV Regularized GDA Projection & Abstention Gate
fit_tier0_calibrated_projection <- function(df, opts = list()) {
  df <- .init_propagation_cols(df, reset = FALSE); n <- nrow(df)
  empty_sum <- list(active = FALSE, lopo_cv_kappa = NA_real_, n_calibrated = 0L, n_abstained = 0L, confusion_matrix = NULL)
  t0 <- resolve_tier0_anchor_overrides(df); is_basal <- .is_basal_mask_prop(df)
  cand_feats <- intersect(c(paste0("PC", 1:6), paste0("ECT_PC", 1:3), "ect_lobe_complexity", "circularity", "aspect_ratio", "AR_inv_std", "tooth_density", "sinus_depth_ratio"), names(df))
  if (n == 0L || length(cand_feats) < 2L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }

  X_raw <- matrix(NA_real_, nrow = n, ncol = length(cand_feats))
  for (j in seq_along(cand_feats)) X_raw[, j] <- suppressWarnings(as.numeric(df[[cand_feats[j]]]))
  t0_raw <- which(t0$is_tier0_dna_anchor & is_basal)
  if (length(t0_raw) < 4L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }
  keep_f <- vapply(seq_len(ncol(X_raw)), function(j) { v <- X_raw[t0_raw, j]; mean(is.finite(v)) >= 0.60 && isTRUE(stats::sd(v, na.rm = TRUE) > 1e-8) }, logical(1L))
  if (sum(keep_f) < 2L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }
  X_use <- X_raw[, keep_f, drop = FALSE]; d <- ncol(X_use); row_fin <- apply(X_use, 1L, function(r) mean(is.finite(r)) >= 0.60)
  y_all <- trimws(t0$rectified_taxon); cls_tab <- table(y_all[t0_raw[row_fin[t0_raw]]])
  classes <- sort(names(cls_tab[cls_tab >= 2L])); train_idx <- which(t0$is_tier0_dna_anchor & is_basal & row_fin & (y_all %in% classes))
  if (length(classes) < 2L || length(train_idx) < 4L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }

  pop_col <- .pick_chr_prop(df, c("population_deme_id", "coancestry_cluster"), "")
  lat <- .pick_num_prop(df, "decimalLatitude"); lon <- .pick_num_prop(df, "decimalLongitude")
  pop_vec <- ifelse(nzchar(pop_col), pop_col, ifelse(is.finite(lat) & is.finite(lon), sprintf("%.1f_%.1f", round(lat, 1), round(lon, 1)), as.character(seq_len(n))))
  oof_logits <- matrix(NA_real_, nrow = length(train_idx), ncol = length(classes))
  u_pops <- unique(pop_vec[train_idx])
  eval_pop <- function(pop) {
    va_pos <- which(pop_vec[train_idx] == pop); tr_pos <- setdiff(seq_along(train_idx), va_pos)
    if (length(unique(y_all[train_idx[tr_pos]])) == length(classes)) {
      list(vp = va_pos, lg = .fit_gda_predict(X_use[train_idx[tr_pos], , drop = FALSE], y_all[train_idx[tr_pos]], X_use[train_idx[va_pos], , drop = FALSE], classes)$logits)
    } else {
      m_vp <- matrix(NA_real_, nrow = length(va_pos), ncol = length(classes))
      for (k in seq_along(va_pos)) {
        tr1 <- setdiff(seq_along(train_idx), va_pos[k])
        if (length(unique(y_all[train_idx[tr1]])) == length(classes)) m_vp[k, ] <- .fit_gda_predict(X_use[train_idx[tr1], , drop = FALSE], y_all[train_idx[tr1]], X_use[train_idx[va_pos[k]], , drop = FALSE], classes)$logits[1L, ]
      }
      list(vp = va_pos, lg = m_vp)
    }
  }
  nw_pop <- if (length(u_pops) >= 4L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, length(u_pops)) else 1L
  pop_res <- if (nw_pop > 1L) packera_mclapply(u_pops, eval_pop, mc.cores = nw_pop) else lapply(u_pops, eval_pop)
  for (pr in pop_res) oof_logits[pr$vp, ] <- pr$lg

  oof_ok <- stats::complete.cases(oof_logits); y_tr <- y_all[train_idx]; T_opt <- 1.0; kappa <- NA_real_
  cm <- matrix(0L, nrow = length(classes), ncol = length(classes), dimnames = list(classes, classes))
  if (sum(oof_ok) >= 4L) {
    y_idx <- match(y_tr[oof_ok], classes)
    T_opt <- max(1.0, min(5.0, stats::optimize(function(Tv) { pr <- .softmax_T(oof_logits[oof_ok, , drop = FALSE], Tv); -mean(log(pmax(pr[cbind(seq_len(sum(oof_ok)), y_idx)], 1e-12))) }, c(1.0, 5.0))$minimum))
    oof_pred <- classes[max.col(.softmax_T(oof_logits[oof_ok, , drop = FALSE], T_opt), ties.method = "first")]
    cm <- as.matrix(table(factor(y_tr[oof_ok], levels = classes), factor(oof_pred, levels = classes)))
    n_tot <- sum(cm); pe <- sum(rowSums(cm) * colSums(cm)) / max(1.0, n_tot^2)
    kappa <- if (n_tot == 0L) NA_real_ else if (abs(1.0 - pe) < 1e-12) 1.0 else (sum(diag(cm)) / n_tot - pe) / (1.0 - pe)
  }

  proj_idx <- which(is_basal & row_fin & !t0$is_tier0_dna_anchor)
  if (length(proj_idx) > 0L) {
    full_fit <- .fit_gda_predict(X_use[train_idx, , drop = FALSE], y_tr, X_use[proj_idx, , drop = FALSE], classes)
    probs <- .softmax_T(full_fit$logits, T_opt); crit_osr <- stats::qchisq(0.975, df = d); rs <- rowSums(cm)
    for (r in seq_along(proj_idx)) {
      i <- proj_idx[r]; p_row <- probs[r, ]; ord <- order(p_row, decreasing = TRUE)
      c1 <- ord[1L]; c2 <- ord[2L]; p1 <- p_row[c1]; p2 <- p_row[c2]
      pair_err <- if ((rs[c1] + rs[c2]) > 0L) (cm[c1, c2] + cm[c2, c1]) / (rs[c1] + rs[c2]) else 0.0
      min_d2 <- min(full_fit$d2[r, ])
      abstain <- (p1 < 0.70) || ((p1 - p2) < 0.15) || (pair_err > 0.20 && p1 < 0.90) || (!is.finite(min_d2) || min_d2 > crit_osr)
      df$tier0_calibrated_prob[i] <- round(p1, 4L); df$tier0_abstention_flag[i] <- if (abstain) 1L else 0L
      df$tier0_calibrated_taxon[i] <- if (abstain) "Sensu_Lato_DNA_Required" else classes[c1]
      if (!abstain && is.na(df$genomic_propagation_tier[i])) df$genomic_propagation_tier[i] <- "Tier0_Projected_Calibrated"
    }
  }
  attr(df, "tier0_projection_summary") <- list(
    active = TRUE, lopo_cv_kappa = kappa, confusion_matrix = cm,
    n_calibrated = sum(!is.na(df$tier0_calibrated_taxon) & df$tier0_abstention_flag == 0L, na.rm = TRUE),
    n_abstained = sum(df$tier0_abstention_flag == 1L, na.rm = TRUE)
  )
  df
}

# Step E: Propagate coalescent gdi to cohesive clusters/lineages (conflict-only cytotype gating)
propagate_cluster_gdi_governance <- function(df, gdi_df = NULL, morph_X = NULL, morph_S_inv = NULL) {
  df <- .init_propagation_cols(df, reset = FALSE); n <- nrow(df); if (n == 0L) return(df)
  tx_gdi <- list(); cl_gdi <- list(); cl_anchor_pl <- list(); cl_conflict <- list(); tx_anchor_pls <- list()
  tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), ""); qc <- .pick_bool_prop(df, "geno_quality_pass", FALSE)
  adm <- .pick_bool_prop(df, "is_geno_admixed", FALSE); pred <- trimws(.pick_chr_prop(df, "geno_predicted_taxon", ""))
  g_col <- .pick_num_prop(df, "coalescent_gdi"); cl_col <- trimws(.pick_chr_prop(df, c("coancestry_cluster"), ""))
  gmm_v <- .pick_num_prop(df, "gmm_cluster"); eco_v <- trimws(.pick_chr_prop(df, "ecotype_regime", ""))
  sp <- trimws(.pick_chr_prop(df, "species_standardized", "")); rec_ploidy <- .extract_known_ploidy_vec(df)
  lat <- .pick_num_prop(df, "decimalLatitude"); lon <- .pick_num_prop(df, "decimalLongitude")
  mcl_col <- ifelse(is.finite(gmm_v) & nzchar(sp) & sp != "Unknown", paste(sp, as.integer(gmm_v), eco_v, sep = "||"), "")
  src_ok <- (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy")) & qc & !adm & is.finite(g_col) & nzchar(pred) & (pred != "Unknown")
  for (i in which(src_ok)) {
    tx_gdi[[pred[i]]] <- c(tx_gdi[[pred[i]]], g_col[i])
    if (nzchar(cl_col[i])) cl_gdi[[cl_col[i]]] <- c(cl_gdi[[cl_col[i]]], g_col[i])
    if (!is.na(rec_ploidy[i])) tx_anchor_pls[[pred[i]]] <- unique(c(tx_anchor_pls[[pred[i]]], rec_ploidy[i]))
  }
  for (ck in unique(cl_col[nzchar(cl_col)])) {
    pls <- unique(rec_ploidy[cl_col == ck & !is.na(rec_ploidy)])
    if (length(pls) == 1L) cl_anchor_pl[[ck]] <- pls[1L] else if (length(pls) > 1L) { cl_conflict[[ck]] <- TRUE; cl_gdi[[ck]] <- NULL }
  }
  if (!is.null(gdi_df) && is.data.frame(gdi_df) && nrow(gdi_df) > 0L) {
    gt <- trimws(.pick_chr_prop(gdi_df, c("taxon", "species_standardized", "lineage"), "")); gv <- .pick_num_prop(gdi_df, c("coalescent_gdi", "gdi"))
    gcl <- trimws(.pick_chr_prop(gdi_df, c("coancestry_cluster"), ""))
    for (r in which(is.finite(gv))) {
      if (nzchar(gt[r])) tx_gdi[[gt[r]]] <- c(tx_gdi[[gt[r]]], gv[r])
      if (nzchar(gcl[r]) && !isTRUE(cl_conflict[[gcl[r]]])) cl_gdi[[gcl[r]]] <- c(cl_gdi[[gcl[r]]], gv[r])
    }
  }
  if (length(tx_gdi) == 0L && length(cl_gdi) == 0L && !any(src_ok & nzchar(mcl_col))) return(df)
  tx_med <- lapply(tx_gdi, stats::median); cl_med <- lapply(cl_gdi, stats::median)
  if (is.null(morph_X) || is.null(morph_S_inv)) {
    pool <- .compute_morph_pool_prop(df)
    if (!is.null(pool)) { morph_X <- pool$X; morph_S_inv <- pool$S_inv; is_basal <- pool$is_basal } else is_basal <- rep(FALSE, n)
  } else is_basal <- .is_basal_mask_prop(df)
  t0 <- resolve_tier0_anchor_overrides(df)
  mu_taxon <- if (!is.null(morph_X)) .build_taxon_centroids(df, morph_X, is_basal, t0) else list()
  crit_d2 <- if (!is.null(morph_X)) stats::qchisq(0.95, df = ncol(morph_X)) else NA_real_
  cal_tx <- trimws(.pick_chr_prop(df, "tier0_calibrated_taxon", "")); abst <- df$tier0_abstention_flag
  is_cent_raw <- .pick_bool_prop(df, "is_centroid_or_inst_coord", FALSE)
  pure_ev <- t0$evict_from_given_anchor & !t0$is_tier0_dna_anchor; typed_ll <- which(!is.na(rec_ploidy) & !is_cent_raw & is.finite(lat) & is.finite(lon))
  for (i in seq_len(n)) {
    if (adm[i] || pure_ev[i] || (!is.na(abst[i]) && abst[i] == 1L)) next
    if (nzchar(cl_col[i]) && isTRUE(cl_conflict[[cl_col[i]]])) next
    if (nzchar(cl_col[i]) && !is.na(rec_ploidy[i]) && !is.null(cl_anchor_pl[[cl_col[i]]]) && rec_ploidy[i] != cl_anchor_pl[[cl_col[i]]]) next
    if (!is_cent_raw[i] && is.finite(lat[i]) && is.finite(lon[i]) && length(typed_ll) > 0L) {
      deme_pls <- unique(rec_ploidy[typed_ll[.haversine_km_vec(lat[i], lon[i], lat[typed_ll], lon[typed_ll]) <= 2.0]])
      if (length(deme_pls) > 1L || (length(deme_pls) == 1L && !is.na(rec_ploidy[i]) && rec_ploidy[i] != deme_pls[1L])) next
    }
    mcl_res <- .eval_mcl_cluster_gdi_prop(i, mcl_col[i], mcl_col, src_ok, g_col, pred, rec_ploidy, lat, lon, is_cent = is_cent_raw)
    if (mcl_res$conflict) next
    cl_ok <- nzchar(cl_col[i]) && !is.null(cl_med[[cl_col[i]]])
    gval <- if (cl_ok) cl_med[[cl_col[i]]] else if (is.finite(mcl_res$gdi)) mcl_res$gdi else {
      cand <- ""
      if ((df$genomic_propagation_tier[i] %in% c("G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull") ||
        tier[i] %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy")) && nzchar(pred[i]) && pred[i] != "Unknown") {
        cand <- pred[i]
      } else if (!is.na(abst[i]) && abst[i] == 0L && nzchar(cal_tx[i]) && cal_tx[i] != "Sensu_Lato_DNA_Required") {
        cand <- cal_tx[i]
      } else if (nzchar(sp[i]) && !is.null(tx_med[[sp[i]]]) && is_basal[i] && !is.null(morph_X) && all(is.finite(morph_X[i, ])) && !is.null(mu_taxon[[sp[i]]])) {
        dv <- morph_X[i, ] - mu_taxon[[sp[i]]]
        if (is.finite(d2 <- as.numeric(crossprod(dv, morph_S_inv %*% dv))) && d2 <= crit_d2) cand <- sp[i]
      }
      if (nzchar(cand) && !is.null(tx_med[[cand]])) {
        t_pls <- tx_anchor_pls[[cand]]
        if (!is.na(rec_ploidy[i]) && length(t_pls) == 1L && rec_ploidy[i] != t_pls[1L]) NA_real_ else tx_med[[cand]]
      } else NA_real_
    }
    if (is.finite(gval)) {
      df$cluster_propagated_gdi[i] <- round(gval, 4L)
      if (is.na(df$genomic_propagation_tier[i]) && !(tier[i] %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy"))) df$genomic_propagation_tier[i] <- "Cluster_GDI_Governed"
    }
  }
  df
}

# Top-level orchestrator
run_step07_genomic_propagation <- function(df, opts = list()) {
  df <- .init_propagation_cols(df, reset = FALSE)
  ga_df <- .resolve_prop_csv(df, opts, "genomic_assignments_df", c("genomic_assignments_df", "ga_df"), "genomic_assignments.csv")
  gdi_df <- .resolve_prop_csv(df, opts, "genomic_gdi_delimitation", c("genomic_gdi_delimitation", "gdi_df"), "coalescent_gdi_delimitation.csv")
  df_act <- if (nrow(df) > 0L) any(.pick_bool_prop(df, "geno_quality_pass", FALSE) & (.pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), "") %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy", "G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull"))) else FALSE
  ga_act <- (!is.null(ga_df) && is.data.frame(ga_df) && nrow(ga_df) > 0L) && any(.pick_bool_prop(ga_df, c("geno_quality_pass", "qc_pass"), FALSE) & (.pick_chr_prop(ga_df, c("geno_linkage_tier", "linkage_tier"), "") %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy", "G3_Georeferenced_Unvouchered", "G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull")))
  if (nrow(df) == 0L || (!df_act && !ga_act)) {
    attr(df, "genomic_propagation_summary") <- list(active = FALSE, n_g1_exsiccatae = 0L, n_g2_deme_hull = 0L, n_tier0_anchors = 0L, n_tier0_evicted = 0L, n_tier0_calibrated = 0L, n_tier0_abstained = 0L, n_cluster_gdi = 0L, lopo_cv_kappa = NA_real_, table = data.frame())
    return(df)
  }

  df <- .init_propagation_cols(df, reset = TRUE); df <- propagate_exsiccatae_isoduplicates(df, ga_df = ga_df)
  pool <- .compute_morph_pool_prop(df); df <- propagate_conditional_deme_priors(df, ga_df = ga_df, morph_X = pool$X, morph_S_inv = pool$S_inv)
  t0 <- resolve_tier0_anchor_overrides(df); df <- fit_tier0_calibrated_projection(df, opts = opts)
  df <- propagate_cluster_gdi_governance(df, gdi_df = gdi_df, morph_X = pool$X, morph_S_inv = pool$S_inv)

  proj_sum <- attr(df, "tier0_projection_summary"); sp <- trimws(.pick_chr_prop(df, "species_standardized", "Unknown")); sp[!nzchar(sp)] <- "Unknown"
  tx_levels <- sort(unique(c(sp, t0$rectified_taxon[nzchar(t0$rectified_taxon)]))); cnt_tx <- function(m) vapply(tx_levels, function(x) sum(m(x), na.rm = TRUE), integer(1L))
  pred_v <- .pick_chr_prop(df, "geno_predicted_taxon", "")
  tbl <- data.frame(
    taxon = tx_levels, n_total = cnt_tx(function(x) sp == x),
    n_tier0_anchors = cnt_tx(function(x) t0$is_tier0_dna_anchor & t0$rectified_taxon == x),
    n_tier0_evicted = cnt_tx(function(x) t0$evict_from_given_anchor & sp == x),
    n_g1_exsiccatae = cnt_tx(function(x) df$genomic_propagation_tier == "G1_Exsiccatae_Isoduplicate" & t0$rectified_taxon == x),
    n_g2_deme_hull = cnt_tx(function(x) df$genomic_propagation_tier == "G2_Conditional_Deme_Hull" & pred_v == x),
    n_tier0_calibrated = cnt_tx(function(x) df$tier0_calibrated_taxon == x & df$tier0_abstention_flag == 0L),
    n_tier0_abstained = cnt_tx(function(x) sp == x & df$tier0_abstention_flag == 1L),
    n_cluster_gdi = cnt_tx(function(x) t0$rectified_taxon == x & is.finite(df$cluster_propagated_gdi)), stringsAsFactors = FALSE
  )
  attr(df, "genomic_propagation_summary") <- list(
    active = TRUE, n_g1_exsiccatae = sum(df$genomic_propagation_tier == "G1_Exsiccatae_Isoduplicate", na.rm = TRUE),
    n_g2_deme_hull = sum(df$genomic_propagation_tier == "G2_Conditional_Deme_Hull", na.rm = TRUE),
    n_tier0_anchors = sum(t0$is_tier0_dna_anchor, na.rm = TRUE), n_tier0_evicted = sum(t0$evict_from_given_anchor, na.rm = TRUE),
    n_tier0_calibrated = sum(!is.na(df$tier0_calibrated_taxon) & df$tier0_abstention_flag == 0L, na.rm = TRUE),
    n_tier0_abstained = sum(df$tier0_abstention_flag == 1L, na.rm = TRUE), n_cluster_gdi = sum(is.finite(df$cluster_propagated_gdi)),
    lopo_cv_kappa = if (!is.null(proj_sum)) proj_sum$lopo_cv_kappa else NA_real_, table = tbl
  )
  df
}
