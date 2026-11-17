if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.apply_step07_postprocess_and_guards <- function(df, st_res, mv, bg_ctx, cl_ctx, org_sig, opts) {
  n <- nrow(df)
  status_calls <- st_res$status_calls; rec_taxa <- st_res$rec_taxa; triage_prios <- st_res$triage_prios
  triage_acts <- st_res$triage_acts; rationales <- st_res$rationales; treat_broad <- st_res$treat_broad; treat_sub <- st_res$treat_sub
  resurrection_basionym <- st_res$resurrection_basionym
  hybrid_zone_cohesion_pass <- logical(n); hybrid_parent_1 <- rep(NA_character_, n); hybrid_parent_2 <- rep(NA_character_, n)

  if (bg_ctx$can_evaluate_hybrids && any(st_res$hyb_cand_pass)) {
    cand_all <- which(st_res$hyb_cand_pass)
    pk_cand <- vapply(cand_all, function(ci) paste(sort(c(st_res$hyb_cand_p1[ci], st_res$hyb_cand_p2[ci])), collapse = "||"), character(1))
    for (pos_c in seq_along(cand_all)) {
      ci <- cand_all[pos_c]; same_pk <- cand_all[pk_cand == pk_cand[pos_c]]
      lat_ci <- bg_ctx$lat_vec_geo[ci]; lon_ci <- bg_ctx$lon_vec_geo[ci]; sc_ci <- bg_ctx$state_cnty_geo[ci]; local_peers <- integer(0)
      for (cj in same_pk) {
        if (cj == ci) local_peers <- c(local_peers, cj)
        else if (is.finite(lat_ci) && is.finite(lon_ci) && is.finite(bg_ctx$lat_vec_geo[cj]) && is.finite(bg_ctx$lon_vec_geo[cj])) {
          if (.haversine_km_step07(lat_ci, lon_ci, bg_ctx$lat_vec_geo[cj], bg_ctx$lon_vec_geo[cj]) <= 25.0) local_peers <- c(local_peers, cj)
        } else if (nzchar(sc_ci) && nzchar(bg_ctx$state_cnty_geo[cj]) && sc_ci == bg_ctx$state_cnty_geo[cj]) local_peers <- c(local_peers, cj)
      }
      n_dm_v <- length(unique(mv$cat_sheet_vec[local_peers]))
      raw_ci <- if (!is.na(mv$raw_sp_dv[ci])) mv$raw_sp_dv[ci] else ""
      has_x <- grepl("(^|[[:space:]])[x×][[:space:]]+[a-z]", raw_ci, ignore.case = TRUE)
      geno_adm_ci <- mv$has_geno_dv[ci] && !mv$geno_divergent_nonhybrid_v[ci] && !mv$geno_ancient_introg_v[ci] &&
        (mv$geno_adm_dv[ci] || mv$geno_contemp_hybrid_v[ci]) && identical(mv$geno_tier_dv[ci], "G1_Direct_Voucher")
      if (n_dm_v >= 2L || has_x || geno_adm_ci) {
        p_a <- st_res$hyb_cand_p1[ci]; p_b <- st_res$hyb_cand_p2[ci]; t_ci <- if (!is.na(mv$tier_dv[ci])) mv$tier_dv[ci] else "Tier_3_Bronze"
        hybrid_zone_cohesion_pass[ci] <- TRUE; hybrid_parent_1[ci] <- p_a; hybrid_parent_2[ci] <- p_b
        status_calls[ci] <- "Hybrid_Intergrade_Swarm"; rec_taxa[ci] <- sprintf("%s x %s", p_a, p_b)
        triage_prios[ci] <- if (t_ci == "Tier_1_Gold") "HIGH" else "MEDIUM"; triage_acts[ci] <- "Flag_Hybrid_Swarm_Intergrade"
        hcn <- .hybrid_cytotype_annotation(mv$voucher_ploidy_v[ci], mv$taxon_cyto_str_v[ci], p_a, p_b, mv$sp_all_vec, mv$voucher_ploidy_v, mv$taxon_cyto_str_v)
        rationales[ci] <- paste0(sprintf(
          "Stage 2 Confirmed Hybrid Swarm (%s x %s): Syntopic parents verified within %.1f km and %.1f km (co-flowering |dDOY|=%.1fd), parental chord t_AB=%.2f (d2_perp=%.2f), corroborated across n=%d contact-zone vouchers",
          p_a, p_b, ifelse(is.finite(st_res$hyb_cand_d1_km[ci]), st_res$hyb_cand_d1_km[ci], 0.0), ifelse(is.finite(st_res$hyb_cand_d2_km[ci]), st_res$hyb_cand_d2_km[ci], 0.0),
          ifelse(is.finite(st_res$hyb_cand_ddoy[ci]), st_res$hyb_cand_ddoy[ci], 0.0), st_res$hyb_cand_t_ab[ci], ifelse(is.finite(st_res$hyb_cand_d2_perp[ci]), st_res$hyb_cand_d2_perp[ci], 0.0), n_dm_v
        ), hcn, if (!is.null(org_sig)) org_sig$format_stage2_maternal_note(ci, p_a, p_b) else "")
        treat_broad[ci] <- "Hybrid / Introgressant"; treat_sub[ci] <- "Hybrid_Introgressant"
      } else if (isTRUE(mv$cyto_barrier_pass_v[ci])) {
        raw_c_sp <- if (!is.na(mv$raw_sp_dv[ci])) mv$raw_sp_dv[ci] else st_res$hyb_cand_p1[ci]
        g_ci <- if (!is.na(mv$gmm_cl_dv[ci]) && mv$gmm_cl_dv[ci] >= 1L) mv$gmm_cl_dv[ci] else 1L
        st4_ci <- .build_stage4_cytotype_decision(st_res$hyb_cand_p1[ci], raw_c_sp, g_ci,
          sprintf("Stage 4 Intraspecific Cytotype Differentiation: Unreplicated sympatric intermediate toward %s (t_AB=%.2f, n_deme=1) retained as subspecific cytotype series", st_res$hyb_cand_p2[ci], st_res$hyb_cand_t_ab[ci]))
        status_calls[ci] <- st4_ci$status; rec_taxa[ci] <- st4_ci$rec; triage_prios[ci] <- st4_ci$prio
        triage_acts[ci] <- st4_ci$act; rationales[ci] <- st4_ci$rat; treat_broad[ci] <- st4_ci$broad; treat_sub[ci] <- st4_ci$sub
      } else {
        if (!isTRUE(mv$mo_flag_dv[ci])) {
          status_calls[ci] <- "Species_Confirmed"; rec_taxa[ci] <- st_res$hyb_cand_p1[ci]; triage_prios[ci] <- "LOW"
          triage_acts[ci] <- "Accept_Current_Determination"; treat_broad[ci] <- "Confirmed Species"; treat_sub[ci] <- "Confirmed_Species"
        }
        rationales[ci] <- paste0(rationales[ci], sprintf(" [Unreplicated sympatric intermediate toward %s (t_AB=%.2f, n_deme=1); retained as %s under conservative population hybrid gate]",
          st_res$hyb_cand_p2[ci], st_res$hyb_cand_t_ab[ci], st_res$hyb_cand_p1[ci]))
      }
    }
  }

  aw_sheet_vec <- resolve_step07_ecological_authority_weights(df, allow_bronze_prior = TRUE)
  for (s_idx in mv$sheet_idx_list) {
    if (length(s_idx) < 2L || any(mv$is_type_dv[s_idx]) || any(status_calls[s_idx] %in% c("Confirmed_Type_Control", "Type_Synonymy_Candidate"))) next
    if (length(unique(treat_sub[s_idx])) <= 1L && length(unique(status_calls[s_idx])) <= 1L && length(unique(rec_taxa[s_idx])) <= 1L) next
    is_conf_s <- (treat_sub[s_idx] == "Confirmed_Species") | (status_calls[s_idx] == "Species_Confirmed")
    is_basal_s <- cl_ctx$t1_pure_all[s_idx] & mv$orig_has_obs_morph_v[s_idx] & (mv$orig_n_leaves_v[s_idx] >= 1L)
    w_r <- pmax(0.25, aw_sheet_vec[s_idx]) * (1.0 + mv$orig_n_leaves_v[s_idx]) * (1.0 + 1.5 * as.numeric(is_basal_s)) * (1.0 + 0.5 * as.numeric(is_conf_s)) *
      ifelse(is.finite(mv$orig_cda_prob_v[s_idx]), pmax(0.35, mv$orig_cda_prob_v[s_idx]), 0.35)
    has_leaf_s <- mv$orig_has_obs_morph_v[s_idx] & (mv$orig_n_leaves_v[s_idx] >= 1L)
    conf_leaf_s <- is_basal_s & is_conf_s & ((!is.na(mv$orig_cda_pred_v[s_idx]) & mv$orig_cda_pred_v[s_idx] == cl_ctx$sp_clean_vec[s_idx]) | (!is.na(mv$is_typ_morph_dv[s_idx]) & mv$is_typ_morph_dv[s_idx]))
    is_mixed_r <- function(r, ref_r) {
      if (!cl_ctx$t1_pure_all[r] || !cl_ctx$t1_pure_all[ref_r]) return(FALSE)
      nbr <- if (!is.na(mv$bas_dv[r])) mv$bas_dv[r] else mv$orig_n_leaves_v[r]; nbref <- if (!is.na(mv$bas_dv[ref_r])) mv$bas_dv[ref_r] else mv$orig_n_leaves_v[ref_r]
      if (is.na(nbr) || nbr < 2L || is.na(nbref) || nbref < 2L || !mv$orig_has_obs_morph_v[r] || is.na(mv$orig_cda_pred_v[r]) || !nzchar(mv$orig_cda_pred_v[r])) return(FALSE)
      if (!is.na(mv$orig_cda_pred_v[ref_r]) && mv$orig_cda_pred_v[r] == mv$orig_cda_pred_v[ref_r]) return(FALSE)
      d2r <- if (is.finite(mv$orig_tm_d2_v[r])) mv$orig_tm_d2_v[r] else if ("mahalanobis_d2" %in% names(df) && is.finite(suppressWarnings(as.numeric(df$mahalanobis_d2[r])))) suppressWarnings(as.numeric(df$mahalanobis_d2[r])) else if ("cda_mahal_d2" %in% names(df) && is.finite(suppressWarnings(as.numeric(df$cda_mahal_d2[r])))) suppressWarnings(as.numeric(df$cda_mahal_d2[r])) else NA_real_
      mgr <- if (is.finite(mv$cda_marg_dv[r])) mv$cda_marg_dv[r] else if (is.finite(mv$orig_cda_prob_v[r])) max(0.0, 2.0 * mv$orig_cda_prob_v[r] - 1.0) else 0.0
      is.finite(mv$orig_cda_prob_v[r]) && (mv$orig_cda_prob_v[r] >= 0.90) && (mgr >= 0.50) && ((is.finite(d2r) && (d2r > mv$type_morph_chisq95)) || isTRUE(mv$mo_flag_dv[r]))
    }
    copy_ramet_from <- function(r, src_idx) {
      status_calls[r] <<- status_calls[src_idx]; treat_sub[r] <<- treat_sub[src_idx]; treat_broad[r] <<- treat_broad[src_idx]
      rec_taxa[r] <<- rec_taxa[src_idx]; hybrid_zone_cohesion_pass[r] <<- hybrid_zone_cohesion_pass[src_idx]
      hybrid_parent_1[r] <<- hybrid_parent_1[src_idx]; hybrid_parent_2[r] <<- hybrid_parent_2[src_idx]
      resurrection_basionym[r] <<- resurrection_basionym[src_idx]; triage_prios[r] <<- triage_prios[src_idx]; triage_acts[r] <<- triage_acts[src_idx]
      if (!grepl("\\[Sheet-Coherent Ramet Consensus\\]", rationales[r])) rationales[r] <<- paste0(rationales[r], " [Sheet-Coherent Ramet Consensus]")
    }
    if (any(conf_leaf_s)) {
      ref_conf_r <- s_idx[which(conf_leaf_s)[which.max(w_r[which(conf_leaf_s)])]]
      for (kp in seq_along(s_idx)) {
        r <- s_idx[kp]; if (st_res$portal_lag_flag[r] || (is_conf_s[kp] && rec_taxa[r] == cl_ctx$sp_clean_vec[r]) || isTRUE(mv$cyto_barrier_pass_v[r])) next
        if (!(treat_sub[r] %in% c("Hybrid_Introgressant", "Misidentification_Reassign", "Ecophenotypic_Plasticity") || status_calls[r] %in% c("Hybrid_Intergrade_Swarm", "Putative_Hybrid_Introgression", "Ambiguous_Needs_Expert_Review", "Misidentification_Severe"))) next
        if (!is_mixed_r(r, ref_conf_r)) {
          status_calls[r] <- "Species_Confirmed"; treat_sub[r] <- "Confirmed_Species"; treat_broad[r] <- "Confirmed Species"; rec_taxa[r] <- cl_ctx$sp_clean_vec[r]
          hybrid_zone_cohesion_pass[r] <- FALSE; hybrid_parent_1[r] <- NA_character_; hybrid_parent_2[r] <- NA_character_; resurrection_basionym[r] <- NA_character_
          triage_prios[r] <- "RESOLVED"; triage_acts[r] <- "Accept_Current_Determination"
          if (!grepl("\\[Sheet-Coherent Ramet Consensus\\]", rationales[r])) rationales[r] <- paste0(rationales[r], " [Sheet-Coherent Ramet Consensus]")
        }
      }
    } else if (any(has_leaf_s)) {
      cp_idx <- if (any(is_basal_s)) which(is_basal_s) else which(has_leaf_s); best_r <- s_idx[cp_idx[which.max(w_r[cp_idx])]]
      for (r in s_idx) {
        if (r == best_r || st_res$portal_lag_flag[r]) next
        if ((mv$geno_gdi_oversplit_v[r] && treat_sub[best_r] == "Putative_Cryptic_Species") || (mv$tier0_abstain_dv[r] && (treat_sub[best_r] %in% c("Putative_Cryptic_Species", "Subspecies_Ecotype"))) ||
            (isTRUE(mv$cyto_barrier_pass_v[r]) && (treat_sub[best_r] %in% c("Ecophenotypic_Plasticity", "Confirmed_Species")))) next
        if (!is_mixed_r(r, best_r)) copy_ramet_from(r, best_r)
      }
    }
  }

  for (si in which((status_calls %in% c("Misidentification_Reassign", "Misidentification_Severe") | treat_sub == "Misidentification_Reassign") & !is.na(rec_taxa) & !is.na(mv$sp_all_vec) & nzchar(mv$sp_all_vec) & (rec_taxa == mv$sp_all_vec))) {
    if (isTRUE(mv$mo_flag_dv[si])) { status_calls[si] <- "Ecophenotypic_Plasticity"; treat_broad[si] <- "Variety / Ecotype"; treat_sub[si] <- "Ecophenotypic_Plasticity"; triage_prios[si] <- "MEDIUM"; triage_acts[si] <- "Annotate_Taphonomic_Ecophenotype" }
    else { status_calls[si] <- "Species_Confirmed"; treat_broad[si] <- "Confirmed Species"; treat_sub[si] <- "Confirmed_Species"; triage_prios[si] <- "LOW"; triage_acts[si] <- "Accept_Current_Determination" }
    rationales[si] <- paste0(rationales[si], " [Self-reassignment guard retained given determination]")
  }
  for (gi in which(mv$has_geno_dv & !is.na(mv$geno_pred_dv) & nzchar(mv$geno_pred_dv))) {
    if (!grepl("\\[Genomic|Genomic Assignment|Rule A Nomenclatural", rationales[gi])) {
      pb <- character(0)
      if (!is.na(mv$geno_cyt_dv[gi]) && nzchar(mv$geno_cyt_dv[gi])) pb <- c(pb, sprintf("cytotype=%s", mv$geno_cyt_dv[gi]))
      if (is.finite(mv$geno_gdi_dv[gi])) pb <- c(pb, sprintf("gdi=%.2f", mv$geno_gdi_dv[gi]))
      if (is.finite(mv$geno_q12_dv[gi])) pb <- c(pb, sprintf("Q12=%.2f", mv$geno_q12_dv[gi]))
      if (!is.na(mv$geno_hyb_cls_dv[gi]) && nzchar(mv$geno_hyb_cls_dv[gi])) pb <- c(pb, mv$geno_hyb_cls_dv[gi])
      p2s <- if (length(pb) > 0L) paste0("; ", paste(pb, collapse = ", ")) else ""
      anc <- if (mv$geno_ancient_introg_v[gi]) sprintf("; Ancient Introgression / Organellar Capture (f_b=%s, cytonuc=%s%s; nuclear determination retained)", ifelse(is.finite(mv$geno_fb_max_dv[gi]), sprintf("%.3f", mv$geno_fb_max_dv[gi]), "NA"), ifelse(isTRUE(mv$geno_cytonuc_dv[gi]), "1", "0"), if (!is.null(org_sig) && nzchar(org_sig$ref_cytonuc_note[gi])) paste0(", ", org_sig$ref_cytonuc_note[gi]) else "") else ""
      qt <- if (mv$geno_divergent_nonhybrid_v[gi]) "; Divergent_NonHybrid_Lineage" else if (mv$geno_qry_dv[gi]) "; Putative_Novel_Or_Hybrid_Query" else if (mv$geno_adm_dv[gi] && !is.na(mv$geno_p1_dv[gi]) && !is.na(mv$geno_p2_dv[gi])) sprintf("; Admixed q_max=%.2f (%s x %s)", ifelse(is.finite(mv$geno_qmax_dv[gi]), mv$geno_qmax_dv[gi], 0.5), mv$geno_p1_dv[gi], mv$geno_p2_dv[gi]) else ""
      rationales[gi] <- paste0(rationales[gi], sprintf(" [Genomic %s: pred=%s (p=%.2f%s%s%s)]", ifelse(is.na(mv$geno_tier_dv[gi]), "Stream", mv$geno_tier_dv[gi]), mv$geno_pred_dv[gi], mv$geno_prob_dv[gi], p2s, anc, qt))
    }
  }
  for (ai in which(mv$tier0_abstain_dv)) if (!grepl("Tier_0 Abstention", rationales[ai], fixed = TRUE)) rationales[ai] <- paste0(rationales[ai], " [Tier_0 Abstention: Phenotypic overlap across genetically distinct lineages; DNA required]")
  for (mi in which(st_res$is_marginal_ext_v & !st_res$portal_lag_flag & (treat_sub == "Confirmed_Species"))) {
    if (!grepl("Marginal Range Extension Preserved", rationales[mi], fixed = TRUE)) rationales[mi] <- paste0(rationales[mi], sprintf(" [Marginal Range Extension Preserved: rf_given_p=%.3f, clim_z=%s]", mv$rf_given_p_dv[mi], ifelse(is.finite(mv$clim_z_dv[mi]), sprintf("%.2f", mv$clim_z_dv[mi]), "NA")))
  }

  df$authority_weight <- aw_sheet_vec; df$is_type_specimen <- as.logical(mv$is_type_dv)
  if (!("type_status" %in% names(df))) df$type_status <- mv$type_stat_dv
  df$edaphic_z_score <- mv$edaph_z_vec; df$cryptic_cluster_cohesion_pass <- as.logical(st_res$cryptic_cluster_cohesion_pass)
  df$cryptic_non_clinal_pass <- as.logical(st_res$cryptic_non_clinal_pass); df$resurrection_basionym <- resurrection_basionym
  df$hybrid_zone_cohesion_pass <- as.logical(hybrid_zone_cohesion_pass); df$hybrid_parent_1 <- hybrid_parent_1; df$hybrid_parent_2 <- hybrid_parent_2
  df$portal_metadata_lag_flag <- as.logical(st_res$portal_lag_flag); df$portal_lag_trigger_gate <- st_res$portal_lag_gate
  df$is_marginal_range_extension <- as.logical(st_res$is_marginal_ext_v)
  df$taxonomic_status_call <- status_calls; df$recommended_determination <- rec_taxa; df$synthesis_triage_priority <- triage_prios
  df$synthesis_triage_action <- triage_acts; df$synthesis_rationale <- rationales; df$taxonomic_treatment_broad <- treat_broad; df$taxonomic_treatment_subtype <- treat_sub

  if (exists("evaluate_step07_hydro_edaphic_and_basionym_guards", mode = "function")) {
    gw_conn <- NULL
    if (!is.null(opts$db_conn) && requireNamespace("DBI", quietly = TRUE) && isTRUE(tryCatch(DBI::dbIsValid(opts$db_conn), error = function(e) FALSE))) { gw_conn <- opts$db_conn }
    else if (!isTRUE(opts$demo) && trimws(env_plantmatics_or_packera("PLANTMATICS_DEMO_MODE", "PACKERA_DEMO_MODE")) != "1" && exists("connect_packera_db", mode = "function") && requireNamespace("DBI", quietly = TRUE) &&
             (nzchar(as.character(if (!is.null(mv$opts_cyto$db_path)) mv$opts_cyto$db_path else "")) || nzchar(Sys.getenv("RUN_ID", unset = "")) || nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")) ||
              (!is.null(opts$vouchers) && exists("is_sqlite_active_for_inputs", mode = "function") && is_sqlite_active_for_inputs(opts$vouchers, demo = isTRUE(opts$demo))))) {
      gw_conn <- tryCatch(connect_packera_db(db_path = mv$opts_cyto$db_path, read_only = TRUE, demo = FALSE), error = function(e) NULL)
      if (!is.null(gw_conn)) on.exit(tryCatch(if (!is.null(gw_conn) && DBI::dbIsValid(gw_conn)) DBI::dbDisconnect(gw_conn), error = function(e) NULL), add = TRUE)
    }
    pre_res_bas <- df$resurrection_basionym; df <- evaluate_step07_hydro_edaphic_and_basionym_guards(df, db_conn = gw_conn)
    nwpl_v <- if ("nwpl_indicator_code" %in% names(df)) toupper(trimws(as.character(df$nwpl_indicator_code))) else rep(NA_character_, n)
    nwi_d <- if ("nwi_wetland_dist_m" %in% names(df)) suppressWarnings(as.numeric(df$nwi_wetland_dist_m)) else rep(NA_real_, n)
    w_ek <- if ("w_edaphic_kappa" %in% names(df)) suppressWarnings(as.numeric(df$w_edaphic_kappa)) else rep(NA_real_, n)
    cda_p_v <- if ("cda_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(NA_real_, n)
    for (hi in seq_len(n)) {
      if (!is.na(pre_res_bas[hi]) && nzchar(trimws(pre_res_bas[hi])) && identical(df$basionym_synonymy_type[hi], "Homotypic") && is.na(df$resurrection_basionym[hi])) {
        if (identical(df$taxonomic_status_call[hi], "Taxonomic_Resurrection_Candidate")) {
          gk_h <- if ("gmm_cluster" %in% names(df) && is.finite(suppressWarnings(as.integer(df$gmm_cluster[hi])))) max(1L, as.integer(df$gmm_cluster[hi])) else 2L
          df$taxonomic_status_call[hi] <- "Putative_Cryptic_Species"; df$recommended_determination[hi] <- sprintf("%s (Putative Cryptic Lineage K%d)", as.character(df$species_standardized[hi]), gk_h)
          df$synthesis_triage_action[hi] <- "Flag_Putative_Cryptic_Species"; df$taxonomic_treatment_broad[hi] <- "Misidentification / Cryptic Species"; df$taxonomic_treatment_subtype[hi] <- "Putative_Cryptic_Species"
          df$synthesis_rationale[hi] <- sub("Stage 3 Basionym Resurrection Candidate", "Stage 3 Cryptic Congruence (Homotypic Basionym Guarded)", ifelse(is.na(df$synthesis_rationale[hi]), "", df$synthesis_rationale[hi]), fixed = TRUE)
        }
        crh <- ifelse(is.na(df$synthesis_rationale[hi]), "", as.character(df$synthesis_rationale[hi]))
        if (!grepl("\\[Homotypic Basionym Guard:", crh)) df$synthesis_rationale[hi] <- paste0(crh, sprintf(" [Homotypic Basionym Guard: %s shares type with %s; resurrection blocked]", pre_res_bas[hi], as.character(df$species_standardized[hi])))
      }
      is_dry <- (!is.na(df$hydro_edaphic_niche_regime[hi]) && df$hydro_edaphic_niche_regime[hi] %in% c("Ultramafic_Serpentine_Barren", "Calcareous_Carbonate_Glade", "Siliceous_Felsic_Outcrop")) || (is.finite(nwi_d[hi]) && nwi_d[hi] >= 200.0)
      crw <- ifelse(is.na(df$synthesis_rationale[hi]), "", as.character(df$synthesis_rationale[hi]))
      if (!is.na(nwpl_v[hi]) && (nwpl_v[hi] %in% c("OBL", "FACW")) && is_dry && is.finite(df$wetland_fidelity_concordance[hi]) && df$wetland_fidelity_concordance[hi] < 0.25 && !grepl("\\[Hydro-Edaphic Discordance:", crw)) {
        df$synthesis_rationale[hi] <- paste0(crw, sprintf(" [Hydro-Edaphic Discordance: taxon NWPL=%s in %s (wetland_concordance=%.2f%s%s)]", nwpl_v[hi],
          if (!is.na(df$hydro_edaphic_niche_regime[hi])) df$hydro_edaphic_niche_regime[hi] else "Dry_Upland", df$wetland_fidelity_concordance[hi],
          if (is.finite(w_ek[hi])) sprintf(", w_edaphic_kappa=%.2f", w_ek[hi]) else "", if (is.finite(cda_p_v[hi])) sprintf(", morph_cda_p=%.2f", cda_p_v[hi]) else ""))
      }
    }
  }
  if (!is.null(org_sig)) for (oi in which(nzchar(org_sig$oop_rationale_note))) {
    cro <- ifelse(is.na(df$synthesis_rationale[oi]), "", as.character(df$synthesis_rationale[oi]))
    if (!grepl("Out-of-Panel Nuclear Locus Match:", cro, fixed = TRUE)) df$synthesis_rationale[oi] <- paste0(cro, org_sig$oop_rationale_note[oi])
  }
  has_m <- if ("can1" %in% names(df)) !is.na(df$can1) else if ("cda_predicted_taxon" %in% names(df)) !is.na(df$cda_predicted_taxon) else rep(TRUE, n)
  has_e <- if ("edaphic_best_fit_taxon" %in% names(df)) !is.na(df$edaphic_best_fit_taxon) else rep(FALSE, n)
  df$is_completed_voucher <- has_m | has_e | mv$has_geno_dv
  for (cc in c("genomic_propagation_tier", "tier0_calibrated_taxon", "cytology_linkage_tier", "taxon_reported_cytotypes", "hydro_edaphic_niche_regime", "basionym_synonymy_type")) {
    df[[cc]] <- if (cc %in% names(df)) as.character(df[[cc]]) else rep(NA_character_, n)
  }
  for (nc in c("tier0_calibrated_prob", "deme_mahal_d2", "geno_deme_weight_scale", "cluster_propagated_gdi", "wetland_fidelity_concordance")) {
    df[[nc]] <- if (nc %in% names(df)) suppressWarnings(as.numeric(df[[nc]])) else rep(NA_real_, n)
  }
  for (ic in c("tier0_abstention_flag", "voucher_sporophytic_2n", "voucher_ploidy_level", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass")) {
    df[[ic]] <- if (ic %in% names(df)) suppressWarnings(as.integer(df[[ic]])) else rep(NA_integer_, n)
  }
  df
}

