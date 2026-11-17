if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.evaluate_step07_specimen_stages <- function(df, mv, bg_ctx, cl_ctx, org_sig, hybrid_eval_fn, gmm_two_delta_bic = 12.0, is_intra = FALSE) {
  n <- nrow(df); df_nms_d <- names(df)
  status_calls <- character(n); rec_taxa <- character(n); triage_prios <- character(n)
  triage_acts <- character(n); rationales <- character(n); treat_broad <- character(n); treat_sub <- character(n)
  cryptic_cluster_cohesion_pass <- logical(n); cryptic_non_clinal_pass <- logical(n)
  resurrection_basionym <- rep(NA_character_, n); hyb_cand_pass <- logical(n)
  hyb_cand_p1 <- rep(NA_character_, n); hyb_cand_p2 <- rep(NA_character_, n)
  hyb_cand_d1_km <- rep(Inf, n); hyb_cand_d2_km <- rep(Inf, n); hyb_cand_ddoy <- rep(Inf, n)
  hyb_cand_t_ab <- rep(NA_real_, n); hyb_cand_d2_perp <- rep(Inf, n); is_marginal_ext_v <- logical(n)
  outside_self_95_v <- logical(n); tier_auth_ok_v <- logical(n)
  cyto_bar_lbl_attr <- attr(df, "cytotype_barrier_labels"); n_known_taxa <- max(1L, length(bg_ctx$active_taxa))
  env_db_path <- env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")

  set_i <- function(i, st, rc, pr, ac, rt, tb, ts) {
    status_calls[i] <<- st; rec_taxa[i] <<- rc; triage_prios[i] <<- pr; triage_acts[i] <<- ac
    rationales[i] <<- rt; treat_broad[i] <<- tb; treat_sub[i] <<- ts
  }

  for (i in seq_len(n)) {
    given_sp <- if (!is.na(mv$sp_all_vec[i]) && nzchar(mv$sp_all_vec[i])) mv$sp_all_vec[i] else "Unknown"
    raw_sp <- if (!is.na(mv$raw_sp_dv[i])) mv$raw_sp_dv[i] else given_sp
    tier <- if (!is.na(mv$tier_dv[i])) mv$tier_dv[i] else "Tier_3_Bronze"; det_t <- tier
    has_obs_doy_i <- mv$valid_doy_all[i]
    sp_med_doy <- if (given_sp %in% names(mv$med_doy)) mv$med_doy[[given_sp]] else mv$glob_circ$med_doy
    sp_q75_doy <- if (given_sp %in% names(mv$q75_doy)) mv$q75_doy[[given_sp]] else mv$glob_circ$q75_doy
    sp_mu_theta <- if (given_sp %in% names(mv$mu_theta_hop)) mv$mu_theta_hop[[given_sp]] else mv$glob_circ$mu_theta
    sp_q75_arc <- if (given_sp %in% names(mv$q75_arc_hop)) mv$q75_arc_hop[[given_sp]] else mv$glob_circ$q75_arc
    hop_doy_i <- if (has_obs_doy_i && is.finite(mv$hop_doy_all[i])) mv$hop_doy_all[i] else sp_med_doy; doy <- hop_doy_i
    hop_arc_i <- if (has_obs_doy_i && is.finite(mv$hop_doy_all[i])) {
      th_i <- 2.0 * pi * (mv$hop_doy_all[i] / 365.25); atan2(sin(th_i - sp_mu_theta), cos(th_i - sp_mu_theta)) * (365.25 / (2.0 * pi))
    } else NA_real_

    morph_sp <- if (!is.na(mv$cda_pred_dv[i]) && nzchar(mv$cda_pred_dv[i])) mv$cda_pred_dv[i] else NULL
    cda_sec_taxon <- if (!is.na(mv$cda_sec1_dv[i]) && nzchar(mv$cda_sec1_dv[i])) mv$cda_sec1_dv[i] else if (!is.na(mv$cda_sec2_dv[i]) && nzchar(mv$cda_sec2_dv[i])) mv$cda_sec2_dv[i] else NULL
    cda_prob <- mv$cda_prob_dv[i]; cda_sec <- mv$cda_sec_p_dv[i]
    cda_marg <- if ("cda_posterior_margin" %in% df_nms_d && is.finite(mv$cda_marg_dv[i])) as.numeric(mv$cda_marg_dv[i]) else if (!is.null(morph_sp) && is.finite(cda_prob)) max(0.0, 2.0 * cda_prob - 1.0) else NA_real_
    is_mo_flag <- if (!is.na(mv$mo_flag_dv[i])) isTRUE(mv$mo_flag_dv[i]) else if (!is.na(mv$misid_flag_dv[i])) isTRUE(mv$misid_flag_dv[i]) else FALSE
    triage_cat <- if (!is.na(mv$triage_cat_dv[i])) mv$triage_cat_dv[i] else "Clean_MultiModal_Consensus"

    edaph_sp <- if (!is.na(mv$edaph_sp_dv[i]) && nzchar(mv$edaph_sp_dv[i])) mv$edaph_sp_dv[i] else if (!is.na(mv$rf_pred_dv[i]) && nzchar(mv$rf_pred_dv[i])) mv$rf_pred_dv[i] else given_sp
    rf_sp <- if (!is.na(mv$rf_pred_dv[i]) && nzchar(mv$rf_pred_dv[i])) mv$rf_pred_dv[i] else edaph_sp
    rf_prob <- mv$rf_prob_dv[i]; edaph_z_i <- mv$edaph_z_vec[i]; te_z_i <- mv$te_z_dv[i]
    edaph_z <- if (is.finite(edaph_z_i)) edaph_z_i else te_z_i
    d_edaph <- if (!is.na(mv$d_ed1_dv[i])) mv$d_ed1_dv[i] else if (!is.na(mv$d_ed2_dv[i])) mv$d_ed2_dv[i] else if (!is.na(mv$d_ed3_dv[i])) mv$d_ed3_dv[i] else if (is_intra) {
      if (is.finite(edaph_z)) exp(-0.5 * (edaph_z / 1.8)^2) else 0.75
    } else if (edaph_sp != given_sp && is.finite(edaph_z) && edaph_z > 2.0) 0.22 else 0.75
    ph_str <- if (is.finite(mv$ph_num_all[i])) sprintf("%.2f", mv$ph_num_all[i]) else "NA"
    sand_str <- if (is.finite(mv$sand_num_all[i])) sprintf("%.1f%%", mv$sand_num_all[i]) else "NA"
    gmm_unc <- if (!is.na(mv$gmm_unc_dv[i])) mv$gmm_unc_dv[i] else 0.05
    g_cl <- if (!is.na(mv$gmm_cl_dv[i])) mv$gmm_cl_dv[i] else 1L
    ddoy <- if (is.finite(mv$ddoy1_dv[i])) mv$ddoy1_dv[i] else if (is.finite(mv$ddoy2_dv[i])) mv$ddoy2_dv[i] else 0.0
    is_anth <- mv$is_anth_dv[i]; bf_val <- if (is.finite(mv$bf1_dv[i])) mv$bf1_dv[i] else if (is.finite(mv$bf2_dv[i])) mv$bf2_dv[i] else gmm_two_delta_bic
    tau_ecotype_i <- mv$tau_ecotype_vec[i]; tau_cryptic_i <- mv$tau_cryptic_vec[i]

    morph_is_taxon <- .is_valid_packera_taxon_step07(morph_sp, bg_ctx$active_taxa)
    edaph_is_taxon <- .is_valid_packera_taxon_step07(edaph_sp, bg_ctx$active_taxa)
    rf_is_taxon <- .is_valid_packera_taxon_step07(rf_sp, bg_ctx$active_taxa)
    sec_is_taxon <- .is_valid_packera_taxon_step07(cda_sec_taxon, bg_ctx$active_taxa)
    m_tier_src <- mv$m_tier_src_dv[i]; n_t1_i <- mv$n_t1_dv[i]
    is_tier1_pure_morph <- (is.na(m_tier_src) || !nzchar(m_tier_src) || m_tier_src == "Basal_Rosette") && (is.na(n_t1_i) || n_t1_i >= 1L)
    is_tier1_pure_i <- is_tier1_pure_morph

    is_type_specimen_i <- isTRUE(mv$is_type_dv[i])
    type_stat_i <- if (!is.na(mv$type_stat_dv[i]) && nzchar(trimws(mv$type_stat_dv[i])) && tolower(trimws(mv$type_stat_dv[i])) != "none") trimws(mv$type_stat_dv[i]) else "TYPE"
    typ_m_i <- !is.na(mv$is_typ_morph_dv[i]) && isTRUE(mv$is_typ_morph_dv[i]); tm_d2_i <- mv$tm_d2_dv[i]
    nt_m_i <- if (.is_valid_packera_taxon_step07(mv$nt_morph_dv[i], bg_ctx$active_taxa)) mv$nt_morph_dv[i] else NA_character_
    nt_e_i <- if (.is_valid_packera_taxon_step07(mv$nt_niche_dv[i], bg_ctx$active_taxa)) mv$nt_niche_dv[i] else NA_character_
    clim_z_i <- mv$clim_z_dv[i]; rf_given_p_i <- mv$rf_given_p_dv[i]; u_stat_i <- mv$usda_stat_dv[i]
    has_geo_ll_i <- is.finite(bg_ctx$lat_vec_geo[i]) && is.finite(bg_ctx$lon_vec_geo[i])
    is_extralimital_st <- u_stat_i %in% c("out_of_range", "natureserve_exotic_or_waif")
    is_exs_conflict_i <- isTRUE(mv$exs_conf_dv[i]) && !is.na(mv$exs_gold_tx_dv[i]) && nzchar(mv$exs_gold_tx_dv[i]) && (mv$exs_gold_tx_dv[i] != given_sp)
    is_marginal_range_extension_i <- is_extralimital_st && is.finite(rf_given_p_i) && (rf_given_p_i >= 0.35) && (!is.finite(clim_z_i) || clim_z_i <= 3.0)
    if (is_marginal_range_extension_i) is_marginal_ext_v[i] <- TRUE
    has_obs_edaph_i <- (!is.na(mv$edaph_sp_dv[i]) && nzchar(trimws(mv$edaph_sp_dv[i]))) || (!is.na(mv$rf_pred_dv[i]) && nzchar(trimws(mv$rf_pred_dv[i]))) ||
      is.finite(edaph_z_i) || is.finite(te_z_i) || !is.na(mv$d_ed1_dv[i]) || !is.na(mv$d_ed2_dv[i]) || !is.na(mv$d_ed3_dv[i])
    has_obs_edaph_chk <- (!is.na(mv$edaph_sp_dv[i]) && nzchar(trimws(mv$edaph_sp_dv[i]))) || (!is.na(mv$rf_pred_dv[i]) && nzchar(trimws(mv$rf_pred_dv[i]))) ||
      is.finite(edaph_z_i) || !is.na(mv$d_ed1_dv[i]) || !is.na(mv$d_ed2_dv[i]) || !is.na(mv$d_ed3_dv[i])
    in_type_edaphic_env_i <- if (!has_geo_ll_i && is_extralimital_st && !has_obs_edaph_chk) FALSE else if (is.finite(te_z_i)) (te_z_i <= 2.0) else (!is.finite(edaph_z_i) || edaph_z_i <= 2.0)
    bas_i <- if (!is.na(mv$nt_bas_dv[i]) && nzchar(trimws(mv$nt_bas_dv[i]))) trimws(mv$nt_bas_dv[i]) else NA_character_
    is_heterotypic_bas_match_i <- !is.na(bas_i) && exists(".is_homotypic_sensu_stricto_basionym", mode = "function") && !isTRUE(.is_homotypic_sensu_stricto_basionym(bas_i, given_sp)[1L])
    is_sensu_stricto_locked_i <- typ_m_i && in_type_edaphic_env_i && !is_heterotypic_bas_match_i && !is_extralimital_st && !is_exs_conflict_i

    cl_key_i <- paste(given_sp, g_cl, sep = "||")
    cryptic_non_clinal_pass[i] <- if (g_cl > 1L && (cl_key_i %in% names(cl_ctx$cryptic_non_clinal_map))) isTRUE(cl_ctx$cryptic_non_clinal_map[[cl_key_i]]) else FALSE
    clinal_r2_i <- if (g_cl > 1L && (cl_key_i %in% names(cl_ctx$clinal_r2_map))) cl_ctx$clinal_r2_map[[cl_key_i]] else NA_real_
    cryptic_cluster_cohesion_pass[i] <- if (g_cl > 1L && (cl_key_i %in% names(cl_ctx$cryptic_cluster_pass_map))) isTRUE(cl_ctx$cryptic_cluster_pass_map[[cl_key_i]]) else FALSE
    osr_chisq95_i <- mv$osr_chisq95_vec[i]
    osr_cryptic_ok_i <- if (is.finite(mv$min_cong_d2_dv[i])) ((mv$min_cong_d2_dv[i] > osr_chisq95_i) && (!is.finite(mv$osr_p_dv[i]) || mv$osr_p_dv[i] < 0.05)) else if (is.finite(mv$osr_p_dv[i])) (mv$osr_p_dv[i] < 0.05) else (!is_intra || !nzchar(env_db_path))

    edaph_strong_support <- has_obs_edaph_i && ((is.finite(rf_prob) && rf_prob >= 0.55) || (!is.finite(rf_prob) && ((is.finite(d_edaph) && d_edaph < 0.35) || (is.finite(edaph_z) && edaph_z > 2.0))))
    eff_cda_prob <- if (is.finite(mv$cda_oof_p_dv[i]) && is.finite(cda_prob)) min(cda_prob, mv$cda_oof_p_dv[i]) else if (is.finite(mv$cda_oof_p_dv[i])) mv$cda_oof_p_dv[i] else cda_prob
    morph_edaph_alt_agree <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && ((edaph_is_taxon && edaph_sp == morph_sp && edaph_strong_support) || (!is.na(nt_e_i) && nt_e_i == morph_sp && is.finite(te_z_i) && te_z_i > 2.0)) && (is.na(eff_cda_prob) || eff_cda_prob >= 0.55))
    type_morph_alt_agree <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && !is.na(nt_m_i) && nt_m_i == morph_sp && (!is.na(eff_cda_prob) && eff_cda_prob >= 0.85))
    both_alt_agree <- morph_edaph_alt_agree || type_morph_alt_agree
    type_synonymy_agree <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && is.finite(cda_prob) && cda_prob >= 0.80 && ((edaph_is_taxon && edaph_sp == morph_sp && edaph_strong_support) || (!is.na(nt_e_i) && nt_e_i == morph_sp) || (!is.na(nt_m_i) && nt_m_i == morph_sp && cda_prob >= 0.85)))

    if (is_type_specimen_i) {
      alt_cands <- stats::na.omit(c(
        if (morph_is_taxon && morph_sp != given_sp) morph_sp else NA_character_,
        if (edaph_is_taxon && edaph_sp != given_sp && edaph_strong_support) edaph_sp else NA_character_,
        if (!is.na(nt_m_i) && nt_m_i != given_sp) nt_m_i else NA_character_, if (!is.na(nt_e_i) && nt_e_i != given_sp) nt_e_i else NA_character_
      )); alt_cands <- alt_cands[nzchar(alt_cands)]
      if (type_synonymy_agree) {
        set_i(i, "Type_Synonymy_Candidate", given_sp, "CRITICAL", "Evaluate_Heterotypic_Synonymy",
              sprintf("Rule A Nomenclatural Type Control (%s): Protected from Misidentification_Reassign; cross-modal type-to-taxon affinity with %s (evaluate heterotypic synonymy)", type_stat_i, morph_sp), "Confirmed Species", "Confirmed_Species")
      } else if (length(alt_cands) > 0L) {
        set_i(i, "Confirmed_Type_Control", given_sp, "LOW", "Accept_Current_Determination",
              sprintf("Rule A Nomenclatural Type Control (%s): Pinned as %s sensu stricto baseline; partial type-to-taxon affinity noted with %s", type_stat_i, given_sp, names(sort(table(alt_cands), decreasing = TRUE))[1L]), "Confirmed Species", "Confirmed_Species")
      } else {
        set_i(i, "Confirmed_Type_Control", given_sp, "RESOLVED", "Accept_Current_Determination",
              sprintf("Rule A Nomenclatural Type Control (%s): Immutable %s sensu stricto baseline concordant across morphological, phenological, and locus classicus streams", type_stat_i, given_sp), "Confirmed Species", "Confirmed_Species")
      }
      next
    }

    self_df_i <- if (is.finite(mv$self_df_raw[i]) && mv$self_df_raw[i] >= 1.0) as.integer(round(mv$self_df_raw[i])) else mv$p_type_morph
    self_d2_eff <- if (is.finite(mv$self_d2_dv[i])) mv$self_d2_dv[i] else if (is.finite(tm_d2_i)) tm_d2_i else if (!is.null(bg_ctx$morph_X) && !is.null(bg_ctx$morph_centroids) && !is.null(bg_ctx$self_morph_S_inv) && (given_sp %in% rownames(bg_ctx$morph_centroids)) && all(is.finite(bg_ctx$morph_X[i, ]))) {
      dv_s <- bg_ctx$morph_X[i, ] - bg_ctx$morph_centroids[given_sp, ]; self_df_i <- max(1L, ncol(bg_ctx$morph_X)); max(0.0, as.numeric(t(dv_s) %*% bg_ctx$self_morph_S_inv %*% dv_s))
    } else NA_real_
    tm_df_i <- self_df_i; tm_chisq95_i <- stats::qchisq(0.95, df = tm_df_i)
    outside_self_morph_95 <- !is.finite(self_d2_eff) || (self_d2_eff > stats::qchisq(0.95, df = self_df_i)); outside_self_95_v[i] <- outside_self_morph_95

    cl_floor_i <- if (is.finite(mv$cl_thr_dv[i])) min(0.82, max(0.58, mv$cl_thr_dv[i])) else 0.0
    cda_cut_base <- if (det_t == "Tier_1_Gold") 0.78 else if (det_t == "Tier_2_Silver") 0.70 else 0.62
    cda_cut <- max(cda_cut_base, if (det_t == "Tier_3_Bronze") min(0.72, cl_floor_i) else cl_floor_i)

    min_cong_tx_i <- if (!is.na(mv$min_cong_tx_dv[i]) && nzchar(trimws(mv$min_cong_tx_dv[i]))) trimws(mv$min_cong_tx_dv[i]) else NA_character_
    min_cong_d2_i <- mv$min_cong_d2_dv[i]; osr_df_i <- mv$osr_df_vec[i]; osr_chisq90 <- mv$osr_chisq90_vec[i]
    osr_self_outlier <- (is.finite(tm_d2_i) && (tm_d2_i > tm_chisq95_i)) || (is.finite(self_d2_eff) && (self_d2_eff > stats::qchisq(0.95, df = self_df_i)))
    ref_self_d2_norm <- if (is.finite(tm_d2_i) && tm_d2_i > 0) tm_d2_i / max(1.0, tm_df_i) else if (is.finite(self_d2_eff) && self_d2_eff > 0) self_d2_eff / max(1.0, self_df_i) else NA_real_
    osr_strictly_closer <- if (is.finite(ref_self_d2_norm)) (min_cong_d2_i / max(1.0, osr_df_i)) < (0.70 * ref_self_d2_norm) else is_mo_flag
    osr_cda_support <- (!is.null(morph_sp) && morph_is_taxon && morph_sp == min_cong_tx_i && is.finite(eff_cda_prob) && eff_cda_prob >= 0.50) ||
      (!is.na(nt_m_i) && nzchar(nt_m_i) && nt_m_i == min_cong_tx_i) || ((is_intra || !morph_is_taxon) && is.finite(min_cong_d2_i) && (min_cong_d2_i <= mv$osr_chisq85_vec[i]))
    geo_osr_res <- if (!is.na(min_cong_tx_i)) .eval_taxon_geo_plausibility_ctx(i, min_cong_tx_i, bg_ctx, is_intra) else list(plausible = FALSE, weight = 0.0, dist_km = Inf)
    osr_geo_ok <- isTRUE(geo_osr_res$plausible) && (!is_intra || (if (bg_ctx$has_ext_geo_ref) (is.finite(geo_osr_res$dist_km) && geo_osr_res$dist_km <= 250.0) else (!is.finite(geo_osr_res$dist_km) || geo_osr_res$dist_km <= 250.0)))
    osr_d2_ceil <- if (is_intra) stats::qchisq(0.50, df = max(1.0, osr_df_i)) else osr_chisq90
    osr_gold_eco_lock <- (det_t == "Tier_1_Gold") && has_obs_edaph_i && in_type_edaphic_env_i && !is_extralimital_st && !is_exs_conflict_i &&
      !(is.finite(clim_z_i) && clim_z_i > 3.5 && (!is.finite(rf_given_p_i) || rf_given_p_i < 0.25)) &&
      (is.na(mv$exs_gold_tx_dv[i]) || mv$exs_gold_tx_dv[i] != min_cong_tx_i) && (!edaph_is_taxon || edaph_sp == given_sp)
    congener_osr_reassign <- !is_type_specimen_i && is_tier1_pure_i && !osr_gold_eco_lock && is.finite(min_cong_d2_i) && !is.na(min_cong_tx_i) && nzchar(min_cong_tx_i) &&
      (min_cong_tx_i != given_sp) && (min_cong_d2_i <= osr_d2_ceil) && osr_self_outlier && osr_strictly_closer && osr_cda_support && osr_geo_ok &&
      (!is_intra || ((is_mo_flag || (is.finite(edaph_z_i) && edaph_z_i > 1.75)) && (det_t != "Tier_1_Gold" || (is.finite(edaph_z_i) && edaph_z_i > 2.0)))) &&
      (!isTRUE(typ_m_i) || (is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i))

    morph_concordant_self <- (morph_is_taxon && morph_sp == given_sp && is.finite(cda_prob) && cda_prob >= 0.65) || (isTRUE(typ_m_i) && is.finite(tm_d2_i) && tm_d2_i <= tm_chisq95_i)
    vote_morph <- if (morph_concordant_self) given_sp else if (morph_is_taxon && morph_sp != given_sp && is.finite(eff_cda_prob) && eff_cda_prob >= 0.55) morph_sp else if (!is.na(nt_m_i) && nzchar(nt_m_i) && is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i) nt_m_i else if (morph_is_taxon && is.finite(eff_cda_prob) && eff_cda_prob >= 0.55) morph_sp else NA_character_
    rf_cand_sp <- if (rf_is_taxon && rf_sp != given_sp) rf_sp else if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else if (edaph_is_taxon) edaph_sp else if (rf_is_taxon) rf_sp else NA_character_
    rf_vote <- if (has_obs_edaph_i && !is.na(rf_cand_sp) && ((is.finite(rf_prob) && rf_prob >= 0.55) || (!is.finite(rf_prob) && ((is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z_i) && edaph_z_i > 2.0))))) rf_cand_sp else NA_character_
    nt_e_vote <- if (!is.na(nt_e_i) && nzchar(nt_e_i) && is.finite(te_z_i) && te_z_i > 2.0) nt_e_i else NA_character_
    vote_edaph <- if (!is.na(nt_e_vote) && !is.na(rf_vote)) {
      if (nt_e_vote == rf_vote) nt_e_vote else if (is.finite(te_z_i) && te_z_i > 2.5 && (!is.finite(rf_prob) || rf_prob < 0.70)) nt_e_vote else if (is.finite(rf_prob) && rf_prob >= 0.70) rf_vote else NA_character_
    } else if (!is.na(nt_e_vote)) nt_e_vote else rf_vote
    deviates_from_given_type <- (is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i) || (is.finite(mv$tm_conc_dv[i]) && mv$tm_conc_dv[i] < 0.25) || (is.finite(te_z_i) && te_z_i > 2.0)
    type_control_reassign <- !is_type_specimen_i && !morph_concordant_self && outside_self_morph_95 && (!isTRUE(typ_m_i) || (is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i)) &&
      !is.na(vote_morph) && !is.na(vote_edaph) && (vote_morph == vote_edaph) && (vote_morph != given_sp) && deviates_from_given_type
    type_control_target_sp <- if (type_control_reassign) vote_morph else NULL

    has_geno_i <- mv$has_geno_dv[i]; geno_sp_i <- if (has_geno_i && !is.na(mv$geno_pred_dv[i]) && nzchar(mv$geno_pred_dv[i])) mv$geno_pred_dv[i] else NA_character_
    geno_is_taxon <- has_geno_i && !is.na(geno_sp_i) && .is_valid_packera_taxon_step07(geno_sp_i, bg_ctx$active_taxa)
    geno_prob_i <- mv$geno_prob_dv[i]; geno_marg_i <- if (is.finite(mv$geno_marg_dv[i])) mv$geno_marg_dv[i] else if (is.finite(geno_prob_i)) max(0.0, 2.0 * geno_prob_i - 1.0) else 0.0
    geno_tier_i <- mv$geno_tier_dv[i]; geno_g1_direct <- !is.na(geno_tier_i) && (geno_tier_i %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate"))
    geno_unilat_tier_ok <- !is.na(geno_tier_i) && (geno_tier_i %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy"))
    geno_gdi_i <- if (has_geno_i) mv$geno_gdi_dv[i] else NA_real_; eff_gdi_i <- mv$eff_gdi_dv[i]
    is_cluster_prop_gdi_only_i <- mv$is_cluster_prop_gdi_only[i]; tier0_abstain_i <- mv$tier0_abstain_dv[i]
    geno_h_i <- if (has_geno_i) mv$geno_h_dv[i] else NA_real_; geno_q12_i <- if (has_geno_i) mv$geno_q12_dv[i] else NA_real_
    geno_divergent_nonhybrid_i <- mv$geno_divergent_nonhybrid_v[i]; geno_contemp_hybrid_i <- mv$geno_contemp_hybrid_v[i]; geno_ancient_introg_i <- mv$geno_ancient_introg_v[i]
    geno_gdi_oversplit_i <- mv$geno_gdi_oversplit_v[i]; geno_gdi_cryptic_support_i <- mv$geno_gdi_cryptic_support_v[i]
    geno_strong_i <- geno_is_taxon && !mv$geno_adm_dv[i] && !mv$geno_qry_dv[i] && is.finite(geno_prob_i) && (geno_prob_i >= 0.85) && (geno_marg_i >= 0.40)

    cand_alt_sp <- if (geno_unilat_tier_ok && geno_strong_i && (geno_sp_i != given_sp) && (geno_g1_direct || !morph_edaph_alt_agree || (!is.null(morph_sp) && morph_sp == geno_sp_i))) geno_sp_i
      else if (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp) morph_sp else if (!is.null(type_control_target_sp)) type_control_target_sp
      else if (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp) min_cong_tx_i else if (geno_unilat_tier_ok && geno_strong_i && (geno_sp_i != given_sp)) geno_sp_i else NA_character_
    geno_supports_given <- geno_is_taxon && (geno_sp_i == given_sp) && !mv$geno_adm_dv[i] && !mv$geno_qry_dv[i] && is.finite(geno_prob_i) && (geno_prob_i >= 0.65)
    geno_supports_alt <- geno_is_taxon && !is.na(cand_alt_sp) && (geno_sp_i == cand_alt_sp) && !mv$geno_adm_dv[i] && is.finite(geno_prob_i) && (geno_prob_i >= 0.70)
    is_single_ramet_sheet <- is.na(mv$n_sheet_ramets_vec[i]) || (mv$n_sheet_ramets_vec[i] <= 1L)
    eco_supports_given <- has_obs_edaph_i && is_single_ramet_sheet && ((edaph_is_taxon && edaph_sp == given_sp && ((is.finite(rf_prob) && rf_prob >= 0.65) || (!is.finite(rf_prob) && is.finite(d_edaph) && d_edaph >= 0.60))) || (is.finite(edaph_z) && edaph_z <= 1.5 && (!edaph_is_taxon || edaph_sp == given_sp)))
    eco_supports_alt <- has_obs_edaph_i && !is.na(cand_alt_sp) && ((((edaph_is_taxon && edaph_sp == cand_alt_sp) || (rf_is_taxon && rf_sp == cand_alt_sp)) && edaph_strong_support) || (!is.na(nt_e_i) && nt_e_i == cand_alt_sp && is.finite(te_z_i) && te_z_i > 2.0))

    n_sec_giv <- 0L; n_sec_alt <- 0L; vestiture_supports_given <- FALSE
    if (!is.null(bg_ctx$sec_organ_mat) && length(bg_ctx$sec_organ_cols) > 0L && !is.na(cand_alt_sp) && (given_sp %in% names(bg_ctx$sec_organ_med_by_tx)) && (cand_alt_sp %in% names(bg_ctx$sec_organ_med_by_tx))) {
      mg_v <- bg_ctx$sec_organ_med_by_tx[[given_sp]]; ma_v <- bg_ctx$sec_organ_med_by_tx[[cand_alt_sp]]
      for (sj in seq_along(bg_ctx$sec_organ_cols)) {
        vij <- bg_ctx$sec_organ_mat[i, sj]; mj <- bg_ctx$sec_organ_mad_pool[sj]
        if (is.finite(vij) && !is.null(mg_v) && !is.null(ma_v) && is.finite(mg_v[sj]) && is.finite(ma_v[sj]) && is.finite(mj) && abs(ma_v[sj] - mg_v[sj]) >= 1.5 * mj) {
          zg <- abs(vij - mg_v[sj]) / mj; za <- abs(vij - ma_v[sj]) / mj
          if (zg <= 1.25 && za >= 1.75) { n_sec_giv <- n_sec_giv + 1L; if (grepl("^blade_cie_l", bg_ctx$sec_organ_cols[sj])) vestiture_supports_given <- TRUE } else if (za <= 1.25 && zg >= 1.75) n_sec_alt <- n_sec_alt + 1L
        }
      }
    }
    sec_organ_supports_given <- (n_sec_giv >= 1L) && (n_sec_alt == 0L); sec_organ_supports_alt <- (n_sec_alt >= 1L) && (n_sec_giv == 0L)
    pheno_supports_given <- FALSE; pheno_supports_alt <- FALSE
    if (is_single_ramet_sheet && !is.na(cand_alt_sp) && has_obs_doy_i && is_anth && is.finite(mv$hop_doy_all[i]) &&
        (given_sp %in% names(bg_ctx$has_anthesis_doy_by_tx)) && (cand_alt_sp %in% names(bg_ctx$has_anthesis_doy_by_tx)) &&
        isTRUE(bg_ctx$has_anthesis_doy_by_tx[[given_sp]]) && isTRUE(bg_ctx$has_anthesis_doy_by_tx[[cand_alt_sp]]) &&
        (given_sp %in% names(mv$mu_theta_hop)) && (cand_alt_sp %in% names(mv$mu_theta_hop))) {
      th_g <- mv$mu_theta_hop[[given_sp]]; th_a <- mv$mu_theta_hop[[cand_alt_sp]]
      if (is.finite(.circ_arc_days_step07(th_g, th_a)) && .circ_arc_days_step07(th_g, th_a) >= 14.0) {
        th_o <- 2.0 * pi * (mv$hop_doy_all[i] / 365.25); dg <- .circ_arc_days_step07(th_o, th_g); da <- .circ_arc_days_step07(th_o, th_a)
        pheno_supports_given <- (dg <= 10.0) && (da >= 16.0); pheno_supports_alt <- (da <= 10.0) && (dg >= 16.0)
      }
    }
    if (is.finite(clim_z_i) && clim_z_i > 3.5 && (!is.finite(rf_given_p_i) || rf_given_p_i < 0.25)) pheno_supports_given <- FALSE
    exs_supports_given <- !is.na(mv$exs_gold_tx_dv[i]) && nzchar(mv$exs_gold_tx_dv[i]) && (mv$exs_gold_tx_dv[i] == given_sp)
    exs_supports_alt <- !is.na(mv$exs_gold_tx_dv[i]) && nzchar(mv$exs_gold_tx_dv[i]) && !is.na(cand_alt_sp) && (mv$exs_gold_tx_dv[i] == cand_alt_sp)
    any_non_outline_supports_alt <- eco_supports_alt || sec_organ_supports_alt || pheno_supports_alt || exs_supports_alt || geno_supports_alt || (!is.na(nt_m_i) && !is.na(cand_alt_sp) && nt_m_i == cand_alt_sp && is.finite(tm_d2_i) && tm_d2_i > mv$type_morph_chisq95)
    any_non_outline_supports_given <- eco_supports_given || sec_organ_supports_given || pheno_supports_given || exs_supports_given || geno_supports_given

    exs_conflict_reassign <- !is_type_specimen_i && is_exs_conflict_i && !exs_supports_given
    extralimital_biogeo_reassign <- !is_type_specimen_i && !is_sensu_stricto_locked_i && !exs_supports_given && !is_marginal_range_extension_i && is_extralimital_st &&
      is.finite(clim_z_i) && (clim_z_i > 3.5) && (!is.finite(rf_given_p_i) || rf_given_p_i < 0.25) && !(morph_is_taxon && !is.null(morph_sp) && morph_sp == given_sp && is.finite(eff_cda_prob) && eff_cda_prob >= 0.75)
    ungeoref_extralimital_reassign <- !is_type_specimen_i && !is_sensu_stricto_locked_i && !exs_supports_given && !is_marginal_range_extension_i && is_extralimital_st && !has_geo_ll_i &&
      ((!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && is.finite(eff_cda_prob) && eff_cda_prob >= 0.80) || (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp && is.finite(min_cong_d2_i) && min_cong_d2_i <= osr_chisq90 && (!morph_is_taxon || is.null(morph_sp) || morph_sp != given_sp)))
    portal_lag_reassign <- exs_conflict_reassign || extralimital_biogeo_reassign || ungeoref_extralimital_reassign

    in_edaphic_envelope <- (is.finite(edaph_z_i) && edaph_z_i <= 2.0) || (edaph_is_taxon && edaph_sp == given_sp)
    single_borderline_div <- ((!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp) || is_mo_flag) && (is.null(morph_sp) || !morph_is_taxon || morph_sp == given_sp || is.na(eff_cda_prob) || eff_cda_prob < cda_cut)
    dlf_i <- mv$dlf_dv[i]; asym_i <- mv$asym_dv[i]; is_sym_dissected_i <- ((is.finite(dlf_i) && dlf_i > 0.0) || mv$is_dis_dv[i]) && !(is.finite(asym_i) && asym_i > 0.22)
    is_late_foliar_wear <- has_obs_doy_i && is.finite(hop_arc_i) && is.finite(sp_q75_arc) && (hop_arc_i >= sp_q75_arc) && (g_cl == 1L) && single_borderline_div && !is_sym_dissected_i
    is_taphonomic_plasticity <- (triage_cat != "Severe_Triple_Stream_Conflict") && !both_alt_agree && !morph_edaph_alt_agree && !type_control_reassign && !congener_osr_reassign && !portal_lag_reassign && in_edaphic_envelope && is_late_foliar_wear

    non_assigned_known <- (is.na(cda_prob) || cda_prob < 0.70) && !both_alt_agree && !morph_edaph_alt_agree && !type_control_reassign && !congener_osr_reassign
    has_cytotype_barrier_i <- isTRUE(mv$cyto_barrier_pass_v[i]); has_morph_or_gmm_outlier <- is_mo_flag || (g_cl > max(2L, n_known_taxa) && bf_val >= 6.0)
    has_allochronic_barrier <- is_anth && (abs(ddoy) >= tau_cryptic_i); has_edaphic_barrier <- (is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z) && edaph_z > 2.5)
    has_eco_pheno_barrier_ms <- has_allochronic_barrier || has_edaphic_barrier || geno_gdi_cryptic_support_i
    morph_exceeds_chisq95_i <- if (is.finite(self_d2_eff)) (self_d2_eff > stats::qchisq(0.95, df = self_df_i)) else if (is.finite(tm_d2_i)) (tm_d2_i > tm_chisq95_i) else if (is.finite(mv$min_cong_d2_dv[i])) (mv$min_cong_d2_dv[i] > osr_chisq95_i) else (is_mo_flag && outside_self_morph_95)
    cyto_barrier_valid_ms <- has_cytotype_barrier_i && (has_eco_pheno_barrier_ms || (morph_exceeds_chisq95_i && outside_self_morph_95 && osr_cryptic_ok_i && cryptic_cluster_cohesion_pass[i]))
    mod_morph <- is_tier1_pure_i && ((has_morph_or_gmm_outlier && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i))) || (has_cytotype_barrier_i && morph_exceeds_chisq95_i && outside_self_morph_95 && (is_mo_flag || (g_cl > 1L && bf_val >= 6.0) || is.finite(self_d2_eff))))
    n_cryptic_mod <- as.integer(mod_morph) + as.integer(non_assigned_known) + as.integer(has_allochronic_barrier) + as.integer(has_edaphic_barrier) + as.integer(geno_gdi_cryptic_support_i) + as.integer(cyto_barrier_valid_ms)
    raw_is_crypt_pre_gdi <- !is_sensu_stricto_locked_i && !is_taphonomic_plasticity && (cryptic_cluster_cohesion_pass[i] || geno_gdi_cryptic_support_i) &&
      (osr_cryptic_ok_i || (geno_gdi_cryptic_support_i && !is.finite(mv$min_cong_d2_dv[i]))) && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i)) &&
      is_tier1_pure_i && (mod_morph || (geno_gdi_cryptic_support_i && (mv$geno_qry_dv[i] || outside_self_morph_95))) && (non_assigned_known || geno_gdi_cryptic_support_i) && (has_eco_pheno_barrier_ms || cyto_barrier_valid_ms) && (n_cryptic_mod >= 2L)
    is_crypt <- raw_is_crypt_pre_gdi && !geno_gdi_oversplit_i && !tier0_abstain_i; cryptic_oversplit_blocked_i <- raw_is_crypt_pre_gdi && geno_gdi_oversplit_i

    edaph_contradicts_morph <- has_obs_edaph_i && is_single_ramet_sheet && edaph_is_taxon && (edaph_sp == given_sp) && (rf_sp == given_sp) &&
      ((is.finite(rf_prob) && rf_prob >= 0.80) || (!is.finite(rf_prob) && is.finite(d_edaph) && d_edaph >= 0.70)) && (is.na(eff_cda_prob) || eff_cda_prob < 0.88) &&
      (!is.finite(d_edaph) || d_edaph >= 0.65) && !sec_organ_supports_alt && !exs_supports_alt
    tier_authority_gate_ok <- if (det_t == "Tier_1_Gold") outside_self_morph_95 && !exs_supports_given && !vestiture_supports_given && any_non_outline_supports_alt
      else if (det_t == "Tier_2_Silver") outside_self_morph_95 && !exs_supports_given && !vestiture_supports_given && (if (!is_single_ramet_sheet) (any_non_outline_supports_alt || (is.finite(eff_cda_prob) && eff_cda_prob >= 0.88)) else if (any_non_outline_supports_given) any_non_outline_supports_alt else (any_non_outline_supports_alt || (is.finite(eff_cda_prob) && eff_cda_prob >= 0.85)))
      else outside_self_morph_95 && !exs_supports_given && (!vestiture_supports_given || (eco_supports_alt && pheno_supports_alt))
    tier_auth_ok_v[i] <- tier_authority_gate_ok

    single_morph_reassign <- morph_is_taxon && (morph_sp != given_sp) && !is.na(eff_cda_prob) && (eff_cda_prob >= if (is_tier1_pure_i) cda_cut else max(0.85, cda_cut)) && !edaph_contradicts_morph && tier_authority_gate_ok
    eff_triple_conflict <- (triage_cat == "Severe_Triple_Stream_Conflict") && outside_self_morph_95 && (both_alt_agree || morph_edaph_alt_agree || single_morph_reassign || congener_osr_reassign) && tier_authority_gate_ok
    morph_corrupt_cand <- is_mo_flag && morph_is_taxon && (morph_sp != given_sp) && (!is.na(eff_cda_prob) && eff_cda_prob >= (if (is_tier1_pure_i) max(0.68, cda_cut) else 0.85)) && !edaph_contradicts_morph && tier_authority_gate_ok
    geno_uncontradicted_reassign <- !is_type_specimen_i && !is_sensu_stricto_locked_i && geno_unilat_tier_ok && geno_strong_i && (geno_sp_i != given_sp) && (geno_prob_i >= if (geno_g1_direct) 0.85 else 0.92) && !(morph_concordant_self && eco_supports_given)
    is_misid <- (!morph_concordant_self && ((!is_taphonomic_plasticity && !is_crypt && outside_self_morph_95 && tier_authority_gate_ok && (both_alt_agree || morph_edaph_alt_agree || eff_triple_conflict || single_morph_reassign || morph_corrupt_cand)) || (type_control_reassign && !exs_supports_given) || (congener_osr_reassign && !exs_supports_given))) || portal_lag_reassign || (geno_uncontradicted_reassign && !exs_supports_given)

    if (is_misid) {
      morph_allo_w <- (!morph_is_taxon || is.null(morph_sp) || !isTRUE(.eval_taxon_geo_plausibility_ctx(i, morph_sp, bg_ctx, is_intra)$plausible) || (morph_sp == given_sp && is_extralimital_st)); rank_pl <- NULL
      if (portal_lag_reassign) { mv$portal_lag_flag[i] <- TRUE; mv$portal_lag_gate[i] <- if (exs_conflict_reassign) "Gate1_Exsiccatae_Conflict" else if (extralimital_biogeo_reassign) "Gate3_Extralimital_Bioclimatic" else "Gate3_Ungeoref_Extralimital" }
      geno_prio_ov <- geno_uncontradicted_reassign && !is.na(geno_sp_i) && (geno_g1_direct || !morph_edaph_alt_agree || (!is.null(morph_sp) && morph_sp == geno_sp_i))
      raw_target_sp <- if (geno_prio_ov) geno_sp_i else if (exs_conflict_reassign) mv$exs_gold_tx_dv[i]
        else if (extralimital_biogeo_reassign || ungeoref_extralimital_reassign) {
          rank_pl <- .rank_sympatric_congeners_ctx(i, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec, edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob, min_cong_tx_i, nt_m_i, nt_e_i, te_z_i, geno_is_taxon, geno_sp_i, geno_prob_i, cl_ctx, bg_ctx, mv, TRUE)
          if (!is.null(rank_pl$best) && is.finite(rank_pl$max_score) && rank_pl$max_score >= 0.12) rank_pl$best else if (rf_is_taxon && rf_sp != given_sp && isTRUE(.eval_taxon_geo_plausibility_ctx(i, rf_sp, bg_ctx, is_intra)$plausible)) rf_sp else if (morph_is_taxon && !is.null(morph_sp) && morph_sp != given_sp && isTRUE(.eval_taxon_geo_plausibility_ctx(i, morph_sp, bg_ctx, is_intra)$plausible)) morph_sp else given_sp
        } else if (type_control_reassign && !is.null(type_control_target_sp)) type_control_target_sp else if (both_alt_agree && morph_is_taxon) morph_sp
        else if (congener_osr_reassign && !is.na(min_cong_tx_i)) min_cong_tx_i else if (geno_uncontradicted_reassign && !is.na(geno_sp_i)) geno_sp_i
        else if (morph_is_taxon && morph_sp != given_sp) morph_sp else if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else given_sp
      raw_geo <- .eval_taxon_geo_plausibility_ctx(i, raw_target_sp, bg_ctx, is_intra)
      geno_ex_geo <- geno_uncontradicted_reassign && !is.na(geno_sp_i) && (raw_target_sp == geno_sp_i) && is.finite(geno_prob_i) && (geno_prob_i >= 0.85)
      rg_note <- if (geno_ex_geo && !isTRUE(raw_geo$plausible)) sprintf(" [Disjunct Population Verified by Genomics (%.0f km from reference hull)]", ifelse(is.finite(raw_geo$dist_km), raw_geo$dist_km, 500.0)) else ""
      if (raw_target_sp != given_sp && !isTRUE(raw_geo$plausible) && !geno_ex_geo) {
        if (is.null(rank_pl)) rank_pl <- .rank_sympatric_congeners_ctx(i, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec, edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob, min_cong_tx_i, nt_m_i, nt_e_i, te_z_i, geno_is_taxon, geno_sp_i, geno_prob_i, cl_ctx, bg_ctx, mv, morph_allo_w)
        if (!is.null(rank_pl$best) && is.finite(rank_pl$max_score) && rank_pl$max_score >= 0.12) { rg_note <- sprintf(" [Range-Gated from allopatric %s (%.0f km outside range)]", raw_target_sp, ifelse(is.finite(raw_geo$dist_km), raw_geo$dist_km, 500.0)); raw_target_sp <- rank_pl$best } else raw_target_sp <- given_sp
      }
      target_sp <- raw_target_sp
      if (target_sp == given_sp) {
        mv$portal_lag_flag[i] <- FALSE; mv$portal_lag_gate[i] <- NA_character_
        set_i(i, "Ecophenotypic_Plasticity", given_sp, "MEDIUM", "Annotate_Taphonomic_Ecophenotype",
              sprintf("Stage 4 Plasticity (Allopatric Reassignment Blocked): Morphological outlier retained as %s because candidate congener lies outside documented geographic range", given_sp), "Variety / Ecotype", "Ecophenotypic_Plasticity"); next
      }
      p_note <- if (mv$portal_lag_flag[i]) {
        sprintf(" [Portal Metadata Lag / Extralimital Audit: %s; identifiedBy=%s; loc=%s/%s; clim_z=%s; rf_given_p=%s%s%s]", mv$portal_lag_gate[i],
                ifelse(!is.na(mv$id_by_dv[i]) && nzchar(trimws(mv$id_by_dv[i])), trimws(mv$id_by_dv[i]), "Unknown"),
                ifelse(nzchar(bg_ctx$state_vec_geo[i]), bg_ctx$state_vec_geo[i], "NA"), ifelse(nzchar(bg_ctx$cnty_raw_geo[i]), bg_ctx$cnty_raw_geo[i], "NA"),
                ifelse(is.finite(clim_z_i), sprintf("%.2f", clim_z_i), "NA"), ifelse(is.finite(rf_given_p_i), sprintf("%.3f", rf_given_p_i), "NA"),
                ifelse(exs_conflict_reassign, sprintf("; exsiccatae_gold=%s", mv$exs_gold_tx_dv[i]), ""),
                if (!is.null(rank_pl) && !is.null(rank_pl$peer_votes) && (target_sp %in% rank_pl$peer_votes$same_coll)) sprintf("; same-collector/deme peer=%s", target_sp) else "")
      } else ""
      r_msg <- if (geno_uncontradicted_reassign && target_sp == geno_sp_i) sprintf("Stage 1 Reassign (Independent Genomic Assignment): %s contradicted by %s genome (%s, posterior=%.2f, margin=%.2f); reassign to %s", given_sp, target_sp, ifelse(is.na(geno_tier_i), "Genomic", geno_tier_i), geno_prob_i, geno_marg_i, target_sp)
        else if (exs_conflict_reassign) sprintf("Stage 1 Reassign (Exsiccatae Portal Metadata Lag): %s conflicts with verified in-range Tier_1_Gold isoduplicate %s; reassign to %s", given_sp, target_sp, target_sp)
        else if (extralimital_biogeo_reassign || ungeoref_extralimital_reassign) sprintf("Stage 1 Reassign (Extralimital / Portal Metadata Lag Audit): %s fails USDA/NatureServe county range (status=%s) and bioclimatic/morphometric support; reassign to %s", given_sp, u_stat_i, target_sp)
        else if (congener_osr_reassign && target_sp == min_cong_tx_i) sprintf("Stage 1 Reassign (Open-Set Congener Match): %s deviates from %s baseline and falls inside the 90%% Mahalanobis ellipsoid of sympatric congener %s (min_congener_mahal_d2=%.2f); reassign to %s", given_sp, given_sp, target_sp, min_cong_d2_i, target_sp)
        else if (type_control_reassign && !is.null(type_control_target_sp) && target_sp == type_control_target_sp) sprintf("Stage 1 Reassign (Type-Baseline Concordance): %s deviates from recorded type control (type_mahal_d2=%s, type_edaphic_z=%s) and matches %s nomenclatural type baseline across >=2 modalities; reassign to %s", given_sp, ifelse(is.finite(tm_d2_i), sprintf("%.2f", tm_d2_i), "NA"), ifelse(is.finite(te_z_i), sprintf("%.2f", te_z_i), "NA"), target_sp, target_sp)
        else sprintf("Stage 1 Reassign: %s contradicted across modalities (consensus_conflict=%.2f); reassign to %s", given_sp, mv$p_conf_dv[i], target_sp)
      set_i(i, "Misidentification_Severe", target_sp, if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH", "Reassign_Determination", paste0(r_msg, rg_note, p_note), "Misidentification / Cryptic Species", "Misidentification_Reassign")
      next
    }

    res_bas_i <- if (g_cl > 1L && (cl_key_i %in% names(cl_ctx$cluster_resurrect_bas_map))) cl_ctx$cluster_resurrect_bas_map[[cl_key_i]] else if (is_heterotypic_bas_match_i && (!is.finite(mv$bas_d2_dv[i]) || mv$bas_d2_dv[i] <= mv$type_morph_chisq95) && (!is.finite(tm_d2_i) || is.finite(mv$bas_d2_dv[i]) || tm_d2_i <= mv$type_morph_chisq95)) bas_i else NA_character_

    if (is_intra) {
      has_ep_in <- (is_anth && abs(ddoy) >= tau_cryptic_i) || (is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z) && edaph_z > 2.5) || geno_gdi_cryptic_support_i
      cyto_v_in <- has_cytotype_barrier_i && (has_ep_in || (outside_self_morph_95 && morph_exceeds_chisq95_i && osr_cryptic_ok_i && cryptic_cluster_cohesion_pass[i]))
      in_m_ok <- is_tier1_pure_morph && ((is_mo_flag && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i))) || ((g_cl > 1L && bf_val >= 6.0) && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i))) || (geno_gdi_cryptic_support_i && outside_self_morph_95) || (cyto_v_in && outside_self_morph_95 && morph_exceeds_chisq95_i))
      n_in_m <- as.integer(in_m_ok) + as.integer(is_anth && abs(ddoy) >= tau_cryptic_i) + as.integer((is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z) && edaph_z > 2.5)) + as.integer(geno_gdi_cryptic_support_i) + as.integer(cyto_v_in)
      raw_in_c <- !is_sensu_stricto_locked_i && !is_taphonomic_plasticity && (cryptic_cluster_cohesion_pass[i] || geno_gdi_cryptic_support_i) && (osr_cryptic_ok_i || (geno_gdi_cryptic_support_i && !is.finite(mv$min_cong_d2_dv[i]))) && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i)) && in_m_ok && (is.na(cda_prob) || cda_prob < 0.70 || geno_gdi_cryptic_support_i) && (has_ep_in || cyto_v_in) && (n_in_m >= 2L)
      is_cry_in <- raw_in_c && !geno_gdi_oversplit_i && !tier0_abstain_i
      has_ep_div <- (is_anth && abs(ddoy) >= tau_ecotype_i) || (d_edaph < 0.45) || (is.finite(edaph_z) && edaph_z > 2.0) || (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60)
      if (is_cry_in) {
        p_c_n <- if (is_cluster_prop_gdi_only_i && is.finite(eff_gdi_i) && (eff_gdi_i > 0.70) && !tier0_abstain_i) sprintf(" [Candidate Member via Cluster-Propagated gdi=%.2f; confirm by targeted sequencing]", eff_gdi_i) else ""
        c_b_n <- if (has_cytotype_barrier_i) { bp <- .resolve_cytotype_barrier_pair(i, given_sp, g_cl, mv$cat_sheet_vec[i], mv$pid_cyto_dv[i], cyto_bar_lbl_attr, mv$sp_all_vec, mv$gmm_cl_dv, mv$voucher_ploidy_v, mv$taxon_cyto_str_v, mv$is_type_dv, mv$is_typ_morph_dv); sprintf(" [Postzygotic Heteroploid Barrier: cluster ploidy %s vs type-cluster %s]", bp$cl_lbl, bp$type_lbl) } else ""
        g_s <- if (is.finite(geno_gdi_i)) sprintf(", gdi=%.2f", geno_gdi_i) else ""
        if (!is.na(res_bas_i) && nzchar(res_bas_i)) { resurrection_basionym[i] <- res_bas_i; set_i(i, "Taxonomic_Resurrection_Candidate", sprintf("%s [Resurrect Basionym: %s]", given_sp, res_bas_i), if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH", "Evaluate_Basionym_Resurrection", paste0(sprintf("Stage 3 Basionym Resurrection Candidate: Congruent Basal_Rosette outlier lineage (K%d) matches heterotypic type %s with evolutionary barrier (delta_DOY=%.1fd, tau=%.1fd, D_edaph=%.2f%s)", g_cl, res_bas_i, ddoy, tau_cryptic_i, d_edaph, g_s), c_b_n, p_c_n), "Misidentification / Cryptic Species", "Putative_Cryptic_Species") }
        else set_i(i, "Putative_Cryptic_Species", sprintf("%s (Putative Cryptic Lineage K%d)", given_sp, g_cl), if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH", "Flag_Putative_Cryptic_Species", paste0(sprintf("Congruent Basal_Rosette morphometric outlier + evolutionary barrier (delta_DOY=%.1fd, tau=%.1fd, D_edaph=%.2f%s)", ddoy, tau_cryptic_i, d_edaph, g_s), c_b_n, p_c_n), "Misidentification / Cryptic Species", "Putative_Cryptic_Species")
        next
      } else if (geno_gdi_oversplit_i && raw_in_c) {
        ot_in <- if (is_cluster_prop_gdi_only_i) sprintf("[Cluster Coalescent Oversplit Block: propagated gdi=%.2f]", eff_gdi_i) else sprintf("[Coalescent Oversplit Block: gdi=%.2f < 0.20]", eff_gdi_i)
        if (has_cytotype_barrier_i) { sc <- .build_stage4_cytotype_decision(given_sp, raw_sp, g_cl, sprintf("Stage 4 Intraspecific Cytotype %s: Conspecific coalescent support retains %s at subspecific rank", ot_in, given_sp)); set_i(i, sc$status, sc$rec, sc$prio, sc$act, sc$rat, sc$broad, sc$sub); next }
        else if (is_taphonomic_plasticity || is_mo_flag || g_cl > 1L) { set_i(i, "Ecophenotypic_Plasticity", given_sp, "MEDIUM", "Annotate_Taphonomic_Ecophenotype", sprintf("Stage 4 Plasticity %s: Conspecific coalescent support retains %s within intraspecific envelope", ot_in, given_sp), "Variety / Ecotype", "Ecophenotypic_Plasticity"); next }
        else { set_i(i, "Species_Confirmed", given_sp, "RESOLVED", "Accept_Current_Determination", sprintf("Stage 5 Confirmed %s: Conspecific coalescent support retains %s", ot_in, given_sp), "Confirmed Species", "Confirmed_Species"); next }
      } else if (!is_tier1_pure_morph && is_mo_flag && !morph_concordant_self && (d_edaph < 0.30 || (is.finite(edaph_z) && edaph_z > 2.5))) {
        set_i(i, "Ecophenotypic_Plasticity", given_sp, "MEDIUM", "Annotate_Taphonomic_Ecophenotype", "Cauline_Only_Fallback morphology with multimodal discordance retained within Ecophenotypic_Plasticity envelope", "Variety / Ecotype", "Ecophenotypic_Plasticity"); next
      } else if (has_cytotype_barrier_i && !is_sensu_stricto_locked_i) {
        br_c <- if (g_cl > 1L && has_ep_div) paste0(sprintf("Secondary GMM morphotype (K%d) with cytotype barrier and edaphic/phenological divergence (delta_DOY=%.1fd, pH=%s, Sand=%s)", g_cl, ddoy, ph_str, sand_str), if (is_heterotypic_bas_match_i) sprintf(" [Available Heterotypic Basionym: %s]", bas_i) else "", if (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60) sprintf(" [Allopatric Spatial Cline R2=%.2f; Chambers & Hillis 2020]", clinal_r2_i) else "")
          else sprintf("Stage 4 Intraspecific Cytotype Differentiation: Morphologically overlapping or non-disjunct heteroploid cluster (GMM K%d, edaphic_z=%s) retained at subspecific rank", max(1L, g_cl), ifelse(is.finite(edaph_z), sprintf("%.2f", edaph_z), "NA"))
        sc <- .build_stage4_cytotype_decision(given_sp, raw_sp, g_cl, br_c); set_i(i, sc$status, sc$rec, sc$prio, sc$act, sc$rat, sc$broad, sc$sub); next
      } else if (!is_sensu_stricto_locked_i && (g_cl > 1L) && !tier0_abstain_i) {
        if (has_ep_div) {
          set_i(i, "Subspecies_Ecotype", sprintf("%s var. morphotype_K%d", given_sp, g_cl), "MEDIUM", "Accept_Subspecific_Treatment",
                paste0(sprintf("Secondary GMM morphotype (K%d) with cohesive edaphic/phenological divergence (delta_DOY=%.1fd, pH=%s, Sand=%s)", g_cl, ddoy, ph_str, sand_str),
                       if (is_heterotypic_bas_match_i) sprintf(" [Available Heterotypic Basionym: %s]", bas_i) else "",
                       if (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60) sprintf(" [Allopatric Spatial Cline R2=%.2f; Chambers & Hillis 2020]", clinal_r2_i) else ""),
                "Variety / Ecotype", "Subspecies_Ecotype"); next
        } else if (!cryptic_cluster_cohesion_pass[i] && (!is.finite(edaph_z) || edaph_z <= 2.5) && (!is_anth || abs(ddoy) < tau_cryptic_i)) {
          if (is_taphonomic_plasticity || is_mo_flag) set_i(i, "Ecophenotypic_Plasticity", given_sp, "MEDIUM", "Annotate_Taphonomic_Ecophenotype", sprintf("Stage 4 Plasticity: Secondary shape mode (GMM K%d) within %s ecological/phenological envelope (edaphic_z=%s)", g_cl, given_sp, ifelse(is.finite(edaph_z), sprintf("%.2f", edaph_z), "NA")), "Variety / Ecotype", "Ecophenotypic_Plasticity")
          else set_i(i, "Species_Confirmed", given_sp, "RESOLVED", "Accept_Current_Determination", sprintf("Stage 5 Confirmed: Secondary unweighted shape mode (GMM K%d) within %s ecological and phenological niche (continuous intraspecific foliar polymorphism / heteroblasty)", g_cl, given_sp), "Confirmed Species", "Confirmed_Species")
          next
        }
      }
    }

    has_morph_i <- isTRUE(mv$has_observed_or_sheet_morph[i]) && is_tier1_pure_i && (!is.null(morph_sp) && !is.na(morph_sp) && nzchar(morph_sp) && is.finite(cda_prob))
    geno_contemp_direct_i <- geno_g1_direct && geno_contemp_hybrid_i && is.finite(geno_q12_i) && (geno_q12_i >= 0.35)
    if (bg_ctx$can_evaluate_hybrids && !geno_divergent_nonhybrid_i && !geno_ancient_introg_i && (has_morph_i || geno_contemp_direct_i) && !both_alt_agree && !is_sensu_stricto_locked_i && (!morph_concordant_self || geno_contemp_direct_i) && (given_sp %in% bg_ctx$cohort_species)) {
      nn_mahal_sp <- .resolve_nearest_mahal_ctx(i, given_sp, bg_ctx$cohort_species, bg_ctx); explicit_p2 <- NA_character_
      if (grepl("(^|[[:space:]])[x×][[:space:]]+[a-z]", raw_sp, ignore.case = TRUE)) {
        mh <- regmatches(raw_sp, regexpr("[x×][[:space:]]+[A-Za-z-]+([[:space:]]+[a-z-]+)?", raw_sp, ignore.case = TRUE))
        if (length(mh) > 0L && nzchar(mh[1L])) {
          rtk <- trimws(sub("^[x×][[:space:]]+", "", mh[1L]))
          g_act <- if (exists("resolve_r_target_genus", mode = "function")) resolve_r_target_genus() else "Packera"
          g_esc <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", g_act)
          crhs <- if (grepl(paste0("^", g_esc, "[[:space:]]+"), rtk, ignore.case = TRUE)) standardize_packera_taxon(rtk) else standardize_packera_taxon(paste(g_act, rtk))
          if (!is.na(crhs) && (crhs %in% bg_ctx$cohort_species) && crhs != given_sp) explicit_p2 <- crhs
        }
      }
      dcands <- unique(stats::na.omit(c(explicit_p2,
        if (has_geno_i && (mv$geno_adm_dv[i] || geno_contemp_hybrid_i) && !is.na(mv$geno_p1_dv[i]) && mv$geno_p1_dv[i] != given_sp) mv$geno_p1_dv[i] else NA_character_,
        if (has_geno_i && (mv$geno_adm_dv[i] || geno_contemp_hybrid_i) && !is.na(mv$geno_p2_dv[i]) && mv$geno_p2_dv[i] != given_sp) mv$geno_p2_dv[i] else NA_character_,
        if (!is.null(org_sig)) org_sig$hq_donor_hint[i] else NA_character_, if (morph_is_taxon && morph_sp != given_sp) morph_sp else NA_character_,
        if (sec_is_taxon && cda_sec_taxon != given_sp) cda_sec_taxon else NA_character_, if (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp) min_cong_tx_i else NA_character_,
        if (!is.na(nt_m_i) && nt_m_i != given_sp) nt_m_i else NA_character_, if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else NA_character_, nn_mahal_sp)))
      dcands <- dcands[dcands %in% bg_ctx$cohort_species & dcands != given_sp]; best_g <- NULL
      if (has_morph_i) for (dc in dcands) {
        ge <- hybrid_eval_fn(i, given_sp, dc, morph_sp, cda_prob, cda_sec_taxon, cda_sec, cda_marg, gmm_unc, edaph_sp, d_edaph, edaph_z, asym_i)
        if (isTRUE(ge$pass) && (is.null(best_g) || ge$score > best_g$score)) { best_g <- ge; best_g$p2 <- dc }
      }
      if (!is.null(best_g)) {
        hyb_cand_pass[i] <- TRUE; hyb_cand_p1[i] <- given_sp; hyb_cand_p2[i] <- best_g$p2; hyb_cand_d1_km[i] <- best_g$d1_km
        hyb_cand_d2_km[i] <- best_g$d2_km; hyb_cand_ddoy[i] <- best_g$ddoy_12; hyb_cand_t_ab[i] <- best_g$t_ab; hyb_cand_d2_perp[i] <- best_g$d2_perp
      } else if (geno_contemp_direct_i) {
        p1_g <- if (!is.na(mv$geno_p1_dv[i]) && .is_valid_packera_taxon_step07(mv$geno_p1_dv[i], bg_ctx$active_taxa)) mv$geno_p1_dv[i] else given_sp
        p2_g <- if (!is.na(mv$geno_p2_dv[i]) && .is_valid_packera_taxon_step07(mv$geno_p2_dv[i], bg_ctx$active_taxa) && mv$geno_p2_dv[i] != p1_g) mv$geno_p2_dv[i] else if (length(dcands) > 0L) dcands[1L] else NA_character_
        if (!is.na(p1_g) && !is.na(p2_g) && p1_g != p2_g) {
          hyb_cand_pass[i] <- TRUE; hyb_cand_p1[i] <- p1_g; hyb_cand_p2[i] <- p2_g; hyb_cand_d1_km[i] <- 0.0; hyb_cand_d2_km[i] <- 0.0
          hyb_cand_ddoy[i] <- if (is.finite(ddoy)) abs(ddoy) else 0.0; hyb_cand_t_ab[i] <- if (is.finite(geno_h_i)) geno_h_i else 0.5; hyb_cand_d2_perp[i] <- 0.0
        }
      }
    }

    has_infra <- grepl("(^|[^A-Za-z0-9_])(var\\.|subsp\\.|ssp\\.)[[:space:]]+[a-z-]+", raw_sp, ignore.case = TRUE)
    is_cohesive_ecotype <- !is_sensu_stricto_locked_i && !is_crypt && !tier0_abstain_i && (has_infra || ((g_cl > 1L) && ((d_edaph < 0.45) || (is.finite(edaph_z) && edaph_z > 2.0) || (is_anth && abs(ddoy) >= tau_ecotype_i) || (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60))))
    if (is_crypt && !hyb_cand_pass[i]) {
      pn_m <- if (is_cluster_prop_gdi_only_i && is.finite(eff_gdi_i) && (eff_gdi_i > 0.70) && !tier0_abstain_i) sprintf(" [Candidate Member via Cluster-Propagated gdi=%.2f; confirm by targeted sequencing]", eff_gdi_i) else ""
      cn_m <- if (has_cytotype_barrier_i) { bp <- .resolve_cytotype_barrier_pair(i, given_sp, g_cl, mv$cat_sheet_vec[i], mv$pid_cyto_dv[i], cyto_bar_lbl_attr, mv$sp_all_vec, mv$gmm_cl_dv, mv$voucher_ploidy_v, mv$taxon_cyto_str_v, mv$is_type_dv, mv$is_typ_morph_dv); sprintf(" [Postzygotic Heteroploid Barrier: cluster ploidy %s vs type-cluster %s]", bp$cl_lbl, bp$type_lbl) } else ""
      gs <- if (is.finite(geno_gdi_i)) sprintf(", gdi=%.2f", geno_gdi_i) else ""
      if (!is.na(res_bas_i) && nzchar(res_bas_i)) { resurrection_basionym[i] <- res_bas_i; set_i(i, "Taxonomic_Resurrection_Candidate", sprintf("%s [Resurrect Basionym: %s]", given_sp, res_bas_i), if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH", "Evaluate_Basionym_Resurrection", paste0(sprintf("Stage 3 Basionym Resurrection Candidate (Padial et al. 2010): Basal_Rosette outlier lineage (K%d) matches heterotypic type %s (CDA_p=%.2f) + barrier (|delta_DOY|=%.1fd >= %.1fd, D_edaph=%.2f%s)", g_cl, res_bas_i, ifelse(is.na(cda_prob), 0.45, cda_prob), abs(ddoy), tau_cryptic_i, d_edaph, gs), cn_m, pn_m), "Misidentification / Cryptic Species", "Putative_Cryptic_Species") }
      else set_i(i, "Putative_Cryptic_Species", sprintf("%s (Putative Cryptic Lineage K%d)", given_sp, g_cl), if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH", "Flag_Putative_Cryptic_Species", paste0(sprintf("Stage 3 Cryptic Congruence (Padial et al. 2010): Basal_Rosette outlier unassigned to known taxa (CDA_p=%.2f) + barrier (|delta_DOY|=%.1fd >= %.1fd, D_edaph=%.2f%s)", ifelse(is.na(cda_prob), 0.45, cda_prob), abs(ddoy), tau_cryptic_i, d_edaph, gs), cn_m, pn_m), "Misidentification / Cryptic Species", "Putative_Cryptic_Species")
    } else if (has_cytotype_barrier_i && !is_sensu_stricto_locked_i && !hyb_cand_pass[i]) {
      br <- if (cryptic_oversplit_blocked_i) sprintf("Stage 4 Intraspecific Cytotype %s: Conspecific coalescent support retains %s at subspecific rank", if (is_cluster_prop_gdi_only_i) sprintf("[Cluster Coalescent Oversplit Block: propagated gdi=%.2f]", eff_gdi_i) else sprintf("[Coalescent Oversplit Block: gdi=%.2f < 0.20]", eff_gdi_i), given_sp) else sprintf("Stage 4 Intraspecific Cytotype Differentiation: Morphologically overlapping or non-disjunct heteroploid cluster (GMM K%d, D_edaph=%.2f, delta_DOY=%.1fd)", max(1L, g_cl), d_edaph, ddoy)
      sc <- .build_stage4_cytotype_decision(given_sp, raw_sp, g_cl, br); set_i(i, sc$status, sc$rec, sc$prio, sc$act, sc$rat, sc$broad, sc$sub)
    } else if (cryptic_oversplit_blocked_i && !hyb_cand_pass[i]) {
      ot <- if (is_cluster_prop_gdi_only_i) sprintf("[Cluster Coalescent Oversplit Block: propagated gdi=%.2f]", eff_gdi_i) else sprintf("[Coalescent Oversplit Block: gdi=%.2f < 0.20]", eff_gdi_i)
      if (is_taphonomic_plasticity || is_mo_flag || g_cl > 1L) set_i(i, "Ecophenotypic_Plasticity", given_sp, "MEDIUM", "Annotate_Taphonomic_Ecophenotype", sprintf("Stage 4 Plasticity %s: Conspecific coalescent support retains %s within intraspecific envelope", ot, given_sp), "Variety / Ecotype", "Ecophenotypic_Plasticity")
      else set_i(i, "Species_Confirmed", given_sp, "RESOLVED", "Accept_Current_Determination", sprintf("Stage 5 Confirmed %s: Conspecific coalescent support retains %s", ot, given_sp), "Confirmed Species", "Confirmed_Species")
    } else if (is_cohesive_ecotype && !hyb_cand_pass[i]) {
      ilbl <- if (has_infra) { mi <- regmatches(raw_sp, regexpr("(var\\.|subsp\\.|ssp\\.)[[:space:]]+[A-Za-z-]+", raw_sp, ignore.case = TRUE)); if (length(mi) > 0L && nzchar(mi[1L])) paste(given_sp, tolower(mi[1L])) else paste0(given_sp, " var. ecotype") } else if (g_cl > 1L) sprintf("%s var. morphotype_K%d", given_sp, g_cl) else paste0(given_sp, " var. ecotype")
      set_i(i, "Subspecies_Ecotype", ilbl, "LOW", "Accept_Subspecific_Treatment", paste0(sprintf("Stage 4 Ecotype: Cohesive intraspecific lineage (GMM K%d, D_edaph=%.2f, delta_DOY=%.1fd, tau=%.1fd)", g_cl, d_edaph, ddoy, tau_ecotype_i), if (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60) sprintf(" [Allopatric Spatial Cline R2=%.2f; Chambers & Hillis 2020]", clinal_r2_i) else ""), "Variety / Ecotype", "Subspecies_Ecotype")
    } else if ((!is_sensu_stricto_locked_i && !is_cohesive_ecotype && !is_crypt && (is_taphonomic_plasticity || is_mo_flag)) || hyb_cand_pass[i]) {
      dp <- if (is.finite(doy) && is.finite(sp_q75_doy) && doy >= sp_q75_doy) sprintf("DOY %.0f >= Q75 %.0f", doy, sp_q75_doy) else sprintf("DOY %s, Q75 %s", ifelse(is.finite(doy), sprintf("%.0f", doy), "NA"), ifelse(is.finite(sp_q75_doy), sprintf("%.0f", sp_q75_doy), "NA"))
      set_i(i, "Ecophenotypic_Plasticity", given_sp, "MEDIUM", "Annotate_Taphonomic_Ecophenotype", sprintf("Stage 4 Plasticity: Within-envelope taphonomic/ontogenetic variation (%s, edaphic_z=%s, pH=%s, Sand=%s)", dp, ifelse(is.finite(edaph_z_i), sprintf("%.2f", edaph_z_i), "NA"), ph_str, sand_str), "Variety / Ecotype", "Ecophenotypic_Plasticity")
    } else set_i(i, "Species_Confirmed", given_sp, "RESOLVED", "Accept_Current_Determination", "Stage 5 Confirmed: Concordant across available morphological, reproductive phenological, and edaphic streams", "Confirmed Species", "Confirmed_Species")
  }
  list(status_calls = status_calls, rec_taxa = rec_taxa, triage_prios = triage_prios, triage_acts = triage_acts, rationales = rationales,
       treat_broad = treat_broad, treat_sub = treat_sub, cryptic_cluster_cohesion_pass = cryptic_cluster_cohesion_pass,
       cryptic_non_clinal_pass = cryptic_non_clinal_pass, resurrection_basionym = resurrection_basionym, hyb_cand_pass = hyb_cand_pass,
       hyb_cand_p1 = hyb_cand_p1, hyb_cand_p2 = hyb_cand_p2, hyb_cand_d1_km = hyb_cand_d1_km, hyb_cand_d2_km = hyb_cand_d2_km,
       hyb_cand_ddoy = hyb_cand_ddoy, hyb_cand_t_ab = hyb_cand_t_ab, hyb_cand_d2_perp = hyb_cand_d2_perp,
       portal_lag_flag = mv$portal_lag_flag, portal_lag_gate = mv$portal_lag_gate, is_marginal_ext_v = is_marginal_ext_v,
       outside_self_morph_95 = outside_self_95_v, tier_authority_gate_ok = tier_auth_ok_v)
}
