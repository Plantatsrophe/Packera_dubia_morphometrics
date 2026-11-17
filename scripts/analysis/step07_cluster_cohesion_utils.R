if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.norm_collector_surname <- function(cv) {
  vapply(as.character(cv), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("")
    x <- tolower(trimws(s))
    x <- gsub("(^|[^a-z0-9_])(et[[:space:]]+al\\.?|jr\\.?|sr\\.?|ii|iii|iv|dr\\.?|prof\\.?)($|[^a-z0-9_])", " ", x)
    x <- strsplit(x, "([[:space:]]*[;&][[:space:]]*|[[:space:]]+(and|with)([[:space:]]+|$))")[[1L]][1L]
    if (grepl(",", x, fixed = TRUE)) {
      sur <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]][1L])
    } else {
      toks <- strsplit(gsub("[^a-z[:space:]-]+", " ", x), "[[:space:]]+")[[1L]]
      toks <- toks[nchar(toks) >= 2L]
      sur <- if (length(toks) > 0L) toks[length(toks)] else ""
    }
    gsub("[^a-z-]", "", sur)
  }, character(1), USE.NAMES = FALSE)
}

.county_eval_for_taxon_ctx <- function(row_i, ctx, bg_ctx) {
  sc_i <- bg_ctx$state_cnty_geo[row_i]; st_i <- bg_ctx$state_vec_geo[row_i]; cty_i <- bg_ctx$cnty_raw_geo[row_i]
  exact_cnty <- nzchar(sc_i) && (ctx %in% names(bg_ctx$usda_counties_by_tx)) && (sc_i %in% bg_ctx$usda_counties_by_tx[[ctx]])
  n_st_cnty <- if (nzchar(st_i) && (ctx %in% names(bg_ctx$usda_state_county_counts_by_tx)) && (st_i %in% names(bg_ctx$usda_state_county_counts_by_tx[[ctx]]))) {
    as.integer(bg_ctx$usda_state_county_counts_by_tx[[ctx]][[st_i]])
  } else 0L
  list(exact_county = exact_cnty, absent_comprehensive_county = nzchar(st_i) && nzchar(cty_i) && !exact_cnty && (n_st_cnty >= 5L))
}

.eval_deme_peer_votes_ctx <- function(row_i, given_sp, cl_ctx, bg_ctx, mv) {
  sur_i <- cl_ctx$coll_surname_all[row_i]; dm_i <- cl_ctx$spatial_deme_vec[row_i]
  sc_i <- bg_ctx$state_cnty_geo[row_i]; lat_i <- bg_ctx$lat_vec_geo[row_i]; lon_i <- bg_ctx$lon_vec_geo[row_i]; cat_i <- mv$cat_sheet_vec[row_i]
  has_ll_i <- is.finite(lat_i) && is.finite(lon_i); same_coll_sp <- character(0); same_deme_sp <- character(0)
  cp <- cl_ctx$clean_peer_df
  if (!is.null(cp) && nrow(cp) > 0L) {
    ok_c <- (cp$cat != cat_i) & (cp$sp != given_sp); near_c <- rep(FALSE, nrow(cp))
    if (has_ll_i) {
      ok_p <- ok_c & is.finite(cp$lat) & is.finite(cp$lon)
      if (any(ok_p)) {
        dkm_p <- .haversine_km_step07(lat_i, lon_i, cp$lat[ok_p], cp$lon[ok_p])
        near_c[which(ok_p)[is.finite(dkm_p) & dkm_p <= 10.0]] <- TRUE
      }
    }
    if (!grepl("^SHEET_", dm_i) && nzchar(dm_i)) near_c[ok_c & nzchar(cp$deme_id) & (cp$deme_id == dm_i)] <- TRUE
    if (!has_ll_i && nzchar(sc_i)) near_c[ok_c & (cp$state_cnty == sc_i)] <- TRUE
    if (any(near_c)) {
      same_deme_sp <- c(same_deme_sp, unique(cp$sp[near_c]))
      if (nzchar(sur_i) && any(sc_hit <- near_c & (cp$surname == sur_i))) same_coll_sp <- c(same_coll_sp, unique(cp$sp[sc_hit]))
    }
  }
  n_all <- length(mv$cat_sheet_vec)
  peer_j <- if (has_ll_i) {
    ok_ll <- which((seq_len(n_all) != row_i) & (mv$cat_sheet_vec != cat_i) & is.finite(bg_ctx$lat_vec_geo) & is.finite(bg_ctx$lon_vec_geo))
    if (length(ok_ll) > 0L) { dkm_j <- .haversine_km_step07(lat_i, lon_i, bg_ctx$lat_vec_geo[ok_ll], bg_ctx$lon_vec_geo[ok_ll]); ok_ll[is.finite(dkm_j) & dkm_j <= 10.0] } else integer(0)
  } else if (nzchar(sc_i)) which((seq_len(n_all) != row_i) & (mv$cat_sheet_vec != cat_i) & (bg_ctx$state_cnty_geo == sc_i)) else integer(0)
  for (pj in peer_j) {
    rf_pj <- mv$rf_pred_dv[pj]; cd_pj <- mv$cda_pred_dv[pj]
    cand_pj <- unique(stats::na.omit(c(
      if (.is_valid_packera_taxon_step07(rf_pj, bg_ctx$active_taxa) && rf_pj != given_sp) rf_pj else NA_character_,
      if (.is_valid_packera_taxon_step07(cd_pj, bg_ctx$active_taxa) && cd_pj != given_sp && isTRUE(.eval_taxon_geo_plausibility_ctx(pj, cd_pj, bg_ctx, mv$is_intra)$plausible)) cd_pj else NA_character_
    )))
    for (cpj in cand_pj) {
      same_deme_sp <- c(same_deme_sp, cpj)
      if (nzchar(sur_i) && cl_ctx$coll_surname_all[pj] == sur_i) same_coll_sp <- c(same_coll_sp, cpj)
    }
  }
  list(same_coll = unique(same_coll_sp), same_deme = unique(same_deme_sp))
}

.score_sympatric_cand_ctx <- function(row_i, ctx, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
                                      edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob, min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
                                      geno_is_taxon, geno_sp_i, geno_prob_i, bg_ctx, mv, morph_allopatric_or_weak_self = FALSE, peer_votes = NULL) {
  gw <- .eval_taxon_geo_plausibility_ctx(row_i, ctx, bg_ctx, mv$is_intra)$weight
  ce <- .county_eval_for_taxon_ctx(row_i, ctx, bg_ctx)
  cnty_mult <- if (ce$absent_comprehensive_county) 0.25 else 1.0; cnty_bonus <- if (ce$exact_county) 0.14 else 0.0
  ms <- if (!is.null(morph_sp) && morph_sp == ctx && is.finite(cda_prob)) cda_prob
  else if (!is.null(cda_sec_taxon) && cda_sec_taxon == ctx && is.finite(cda_sec)) cda_sec
  else if (!is.null(bg_ctx$morph_X) && !is.null(bg_ctx$morph_centroids) && !is.null(bg_ctx$morph_S_inv) && (ctx %in% rownames(bg_ctx$morph_centroids))) {
    xi <- bg_ctx$morph_X[row_i, ]; bm <- if (given_sp %in% rownames(bg_ctx$morph_centroids)) bg_ctx$morph_centroids[given_sp, ] else colMeans(bg_ctx$morph_centroids)
    xi[!is.finite(xi)] <- bm[!is.finite(xi)]; dv <- xi - bg_ctx$morph_centroids[ctx, ]; exp(-0.25 * max(0.0, as.numeric(t(dv) %*% bg_ctx$morph_S_inv %*% dv)))
  } else 0.15
  es <- if ((rf_is_taxon && rf_sp == ctx) || (edaph_is_taxon && edaph_sp == ctx)) { if (is.finite(rf_prob)) (if (morph_allopatric_or_weak_self) 0.55 else 0.35) * rf_prob else 0.20 } else 0.0
  ts <- (if ((!is.na(min_cong_tx_i) && min_cong_tx_i == ctx) || (!is.na(nt_m_i) && nt_m_i == ctx)) 0.20 else 0.0) +
        (if (!is.na(nt_e_i) && nt_e_i == ctx && is.finite(te_z_i) && te_z_i > 2.0) 0.15 else 0.0)
  gs <- if (geno_is_taxon && !is.na(geno_sp_i) && geno_sp_i == ctx && is.finite(geno_prob_i)) 0.85 * geno_prob_i else 0.0
  ps <- if (!is.null(peer_votes) && (ctx %in% peer_votes$same_coll)) 0.20 else if (!is.null(peer_votes) && (ctx %in% peer_votes$same_deme)) 0.10 else 0.0
  gw * cnty_mult * (0.65 * ms + es + ts + gs + cnty_bonus + ps)
}

.rank_sympatric_congeners_ctx <- function(row_i, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
                                          edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob, min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
                                          geno_is_taxon, geno_sp_i, geno_prob_i, cl_ctx, bg_ctx, mv, morph_allopatric_or_weak_self = FALSE) {
  cand_pool <- unique(stats::na.omit(c(
    if (geno_is_taxon && !is.na(geno_sp_i) && geno_sp_i != given_sp) geno_sp_i else NA_character_,
    if (rf_is_taxon && rf_sp != given_sp) rf_sp else NA_character_, if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else NA_character_,
    if (!is.na(nt_e_i) && nt_e_i != given_sp) nt_e_i else NA_character_, if (!is.null(morph_sp) && morph_sp != given_sp) morph_sp else NA_character_,
    if (!is.null(cda_sec_taxon) && cda_sec_taxon != given_sp) cda_sec_taxon else NA_character_,
    if (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp) min_cong_tx_i else NA_character_,
    if (!is.na(nt_m_i) && nt_m_i != given_sp) nt_m_i else NA_character_, setdiff(bg_ctx$active_taxa, given_sp)
  )))
  cand_pool <- cand_pool[vapply(cand_pool, .is_valid_packera_taxon_step07, logical(1), active_taxa = bg_ctx$active_taxa) & (cand_pool != given_sp)]
  plaus_cands <- cand_pool[vapply(cand_pool, function(ctx) isTRUE(.eval_taxon_geo_plausibility_ctx(row_i, ctx, bg_ctx, mv$is_intra)$plausible), logical(1))]
  if (length(plaus_cands) == 0L && nzchar(bg_ctx$state_cnty_geo[row_i])) {
    chits <- names(bg_ctx$usda_counties_by_tx)[vapply(names(bg_ctx$usda_counties_by_tx), function(ctx) ctx != given_sp && .is_valid_packera_taxon_step07(ctx, bg_ctx$active_taxa) && (bg_ctx$state_cnty_geo[row_i] %in% bg_ctx$usda_counties_by_tx[[ctx]]), logical(1))]
    plaus_cands <- chits[vapply(chits, function(ctx) isTRUE(.eval_taxon_geo_plausibility_ctx(row_i, ctx, bg_ctx, mv$is_intra)$plausible), logical(1))]
  }
  if (length(plaus_cands) == 0L) return(list(best = NULL, max_score = -Inf, peer_votes = NULL))
  pv <- .eval_deme_peer_votes_ctx(row_i, given_sp, cl_ctx, bg_ctx, mv)
  scores <- vapply(plaus_cands, function(ctx) {
    .score_sympatric_cand_ctx(row_i, ctx, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec, edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob,
                              min_cong_tx_i, nt_m_i, nt_e_i, te_z_i, geno_is_taxon, geno_sp_i, geno_prob_i, bg_ctx, mv, morph_allopatric_or_weak_self, pv)
  }, numeric(1))
  list(best = plaus_cands[which.max(scores)], max_score = max(scores), peer_votes = pv)
}

.build_step07_cluster_cohesion_context <- function(df, bg_ctx, mv, gmm_two_delta_bic, is_intra) {
  n <- nrow(df); nms <- names(df)
  sp_clean_vec <- ifelse(!is.na(mv$sp_all_vec) & nzchar(mv$sp_all_vec), mv$sp_all_vec, "Unknown")
  g_cl_all <- ifelse(!is.na(mv$gmm_cl_dv), mv$gmm_cl_dv, 1L)
  t1_pure_all <- (is.na(mv$m_tier_src_dv) | !nzchar(mv$m_tier_src_dv) | mv$m_tier_src_dv == "Basal_Rosette") & (is.na(mv$n_t1_dv) | mv$n_t1_dv >= 1L)
  edaph_z_all <- ifelse(is.finite(mv$edaph_z_vec), mv$edaph_z_vec, mv$te_z_dv)
  ed_sp_all <- ifelse(!is.na(mv$edaph_sp_dv) & nzchar(mv$edaph_sp_dv), mv$edaph_sp_dv, ifelse(!is.na(mv$rf_pred_dv) & nzchar(mv$rf_pred_dv), mv$rf_pred_dv, sp_clean_vec))
  d_edaph_all <- ifelse(!is.na(mv$d_ed1_dv), mv$d_ed1_dv, ifelse(!is.na(mv$d_ed2_dv), mv$d_ed2_dv, ifelse(!is.na(mv$d_ed3_dv), mv$d_ed3_dv,
    if (is_intra) ifelse(is.finite(edaph_z_all), exp(-0.5 * (edaph_z_all / 1.8)^2), 0.75) else ifelse(ed_sp_all != sp_clean_vec & is.finite(edaph_z_all) & edaph_z_all > 2.0, 0.22, 0.75))))
  ddoy_all <- ifelse(is.finite(mv$ddoy1_dv), mv$ddoy1_dv, ifelse(is.finite(mv$ddoy2_dv), mv$ddoy2_dv, NA_real_))
  bf_all <- ifelse(is.finite(mv$bf1_dv), mv$bf1_dv, ifelse(is.finite(mv$bf2_dv), mv$bf2_dv, gmm_two_delta_bic))
  rec_by_v <- .chr_col(df, "recordedBy", n, NA_character_); rec_num_v <- .chr_col(df, "recordNumber", n, NA_character_)
  exs_key_v <- .chr_col(df, "exsiccatae_key", n, NA_character_); yr_v <- .chr_col(df, "year", n, NA_character_)
  if (all(is.na(yr_v)) && "eventDate" %in% nms) yr_v <- substr(as.character(df$eventDate), 1L, 4L)
  lat_v_cl <- bg_ctx$lat_vec_geo; lon_v_cl <- bg_ctx$lon_vec_geo
  inst_v <- .chr_col(df, "institutionCode", n, NA_character_); deme_id_v <- .chr_col(df, "population_deme_id", n, NA_character_)
  if (all(is.na(deme_id_v))) deme_id_v <- .chr_col(df, "cv_deme_id", n, NA_character_)
  azonal_flag_all <- rep(FALSE, n)
  for (az_c in intersect(c("is_azonal_lithology", "is_rock_outcrop", "substrate_specialist_flag"), nms)) { az_b <- suppressWarnings(as.logical(df[[az_c]])); azonal_flag_all <- azonal_flag_all | (!is.na(az_b) & az_b) }
  if ("matches_type_lithology" %in% nms) { mtl <- suppressWarnings(as.logical(df$matches_type_lithology)); azonal_flag_all <- azonal_flag_all | (!is.na(mtl) & !mtl) }

  spatial_deme_vec <- ifelse(!is.na(deme_id_v) & nzchar(trimws(deme_id_v)), trimws(deme_id_v), "")
  miss_dm <- which(!nzchar(spatial_deme_vec))
  if (length(miss_dm) > 0L) {
    ok_dm <- miss_dm[is.finite(lat_v_cl[miss_dm]) & is.finite(lon_v_cl[miss_dm])]
    if (length(ok_dm) >= 1L) {
      dm_cls <- list()
      for (ii in seq_along(ok_dm)) {
        placed <- FALSE; la_i <- lat_v_cl[ok_dm[ii]]; lo_i <- lon_v_cl[ok_dm[ii]]
        for (cid in seq_along(dm_cls)) {
          mem <- dm_cls[[cid]]; dkm <- .haversine_km_step07(la_i, lo_i, lat_v_cl[ok_dm[mem]], lon_v_cl[ok_dm[mem]])
          if (all(is.finite(dkm) & dkm <= 10.0)) { dm_cls[[cid]] <- c(mem, ii); placed <- TRUE; break }
        }
        if (!placed) dm_cls[[length(dm_cls) + 1L]] <- ii
      }
      for (cid in seq_along(dm_cls)) spatial_deme_vec[ok_dm[dm_cls[[cid]]]] <- sprintf("S06_DEME_%04d", cid)
    }
    sm <- which(!nzchar(spatial_deme_vec))
    if (length(sm) > 0L) spatial_deme_vec[sm] <- ifelse(nzchar(bg_ctx$state_cnty_geo[sm]), paste0("CNTY_", bg_ctx$state_cnty_geo[sm]), paste0("SHEET_", mv$cat_sheet_vec[sm]))
  }
  coll_surname_all <- .norm_collector_surname(rec_by_v)

  idx_cp <- which(mv$anchor_all & !mv$mo_all & !mv$noise_all & !mv$ext_non_type_all & nzchar(sp_clean_vec) & (sp_clean_vec != "Unknown"))
  df_cp <- if (length(idx_cp) > 0L) data.frame(sp = sp_clean_vec[idx_cp], surname = coll_surname_all[idx_cp], deme_id = spatial_deme_vec[idx_cp],
    state_cnty = bg_ctx$state_cnty_geo[idx_cp], lat = bg_ctx$lat_vec_geo[idx_cp], lon = bg_ctx$lon_vec_geo[idx_cp], cat = mv$cat_sheet_vec[idx_cp], stringsAsFactors = FALSE) else NULL
  db_cp <- NULL
  if (!is.null(bg_ctx$db_syntopy_rows) && nrow(bg_ctx$db_syntopy_rows) > 0L) {
    dbr <- bg_ctx$db_syntopy_rows; tx_cp <- .resolve_db_syntopy_taxa(dbr); kp <- !is.na(tx_cp) & nzchar(tx_cp) & (tx_cp != "Unknown")
    if (any(kp)) {
      is_cdb <- if (exists(".is_step07_centroid_coord", mode = "function")) .is_step07_centroid_coord(dbr[kp, , drop = FALSE]) else rep(FALSE, sum(kp))
      st_cp <- .norm_state_step07(dbr$stateProvince[kp]); cn_cp <- .norm_county_step07(dbr$county[kp])
      la_cp <- suppressWarnings(as.numeric(dbr$decimalLatitude[kp])); lo_cp <- suppressWarnings(as.numeric(dbr$decimalLongitude[kp]))
      la_cp[is_cdb] <- NA_real_; lo_cp[is_cdb] <- NA_real_
      db_cp <- data.frame(sp = tx_cp[kp], surname = if ("recordedBy" %in% names(dbr)) .norm_collector_surname(dbr$recordedBy[kp]) else rep("", sum(kp)),
        deme_id = rep("", sum(kp)), state_cnty = ifelse(nzchar(st_cp) & nzchar(cn_cp), paste(st_cp, cn_cp, sep = "||"), ""),
        lat = la_cp, lon = lo_cp, cat = as.character(dbr$catalogNumber[kp]), stringsAsFactors = FALSE)
    }
  }
  res_cp <- if (!is.null(df_cp) && !is.null(db_cp)) rbind(df_cp, db_cp) else if (!is.null(df_cp)) df_cp else db_cp
  clean_peer_df <- if (!is.null(res_cp) && nrow(res_cp) > 0L) res_cp else NULL

  rn_clean <- tolower(gsub("[^a-z0-9]+", "", ifelse(is.na(rec_num_v), "", rec_num_v))); rn_clean[rn_clean %in% c("sn", "na", "none", "unknown", "0")] <- ""
  exs_clean <- tolower(trimws(ifelse(is.na(exs_key_v), "", exs_key_v))); exs_clean[exs_clean %in% c("na", "none", "unknown", "null")] <- ""
  gathering_event_key_all <- ifelse(nzchar(exs_clean), paste0("EXS::", exs_clean),
    ifelse(nzchar(coll_surname_all) & nzchar(rn_clean), paste0("COLLNUM::", coll_surname_all, "::", rn_clean),
           paste0("EVT::", ifelse(nzchar(coll_surname_all), coll_surname_all, "unk_coll"), "::", ifelse(!is.na(yr_v) & nzchar(trimws(yr_v)), trimws(yr_v), "unk_yr"), "::", spatial_deme_vec)))

  cryptic_cluster_pass_map <- list(); cryptic_non_clinal_map <- list(); clinal_r2_map <- list(); cluster_resurrect_bas_map <- list()
  u_sc <- unique(data.frame(sp = sp_clean_vec, cl = g_cl_all, stringsAsFactors = FALSE))
  for (ri in seq_len(nrow(u_sc))) {
    sp_k <- u_sc$sp[ri]; cl_k <- u_sc$cl[ri]; if (is.na(cl_k) || cl_k <= 1L) next
    mk <- paste(sp_k, cl_k, sep = "||"); idx_cl <- which(sp_clean_vec == sp_k & g_cl_all == cl_k & t1_pure_all)
    if (length(idx_cl) == 0L) { cryptic_cluster_pass_map[[mk]] <- FALSE; cryptic_non_clinal_map[[mk]] <- FALSE; next }
    has_cc <- any(!is.na(rec_by_v) & nzchar(trimws(rec_by_v))) || any(nzchar(exs_clean))
    pass_events <- if (n < 50L && !has_cc) length(idx_cl) >= 1L else {
      u_co <- unique(coll_surname_all[idx_cl][nzchar(coll_surname_all[idx_cl])]); u_in <- unique(tolower(trimws(inst_v[idx_cl][!is.na(inst_v[idx_cl]) & nzchar(trimws(inst_v[idx_cl]))])))
      req_dm <- if (length(unique(spatial_deme_vec[sp_clean_vec == sp_k & t1_pure_all])) >= 3L) 3L else max(1L, length(unique(spatial_deme_vec[sp_clean_vec == sp_k & t1_pure_all])))
      (length(unique(gathering_event_key_all[idx_cl])) >= 5L) && (length(unique(spatial_deme_vec[idx_cl])) >= req_dm) && (if ((length(u_co) > 0L) || (length(u_in) > 0L)) (length(u_co) >= 2L || length(u_in) >= 2L) else TRUE)
    }
    fin_c <- idx_cl[is.finite(mv$min_cong_d2_dv[idx_cl])]
    pass_osr_d2 <- if (length(fin_c) > 0L) stats::median(mv$min_cong_d2_dv[fin_c] - mv$osr_chisq95_vec[fin_c], na.rm = TRUE) > 0.0 else (!is_intra || !nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")))
    fin_p <- idx_cl[is.finite(mv$osr_p_dv[idx_cl]) & is.finite(mv$min_cong_d2_dv[idx_cl])]
    pass_osr <- isTRUE(pass_osr_d2 && (if (length(fin_p) > 0L) (stats::median(mv$osr_p_dv[fin_p], na.rm = TRUE) < 0.05) else TRUE))
    cy_cl <- mv$cyto_barrier_pass_v[idx_cl]
    div_v <- (is.finite(d_edaph_all[idx_cl]) & d_edaph_all[idx_cl] < 0.30) | (is.finite(edaph_z_all[idx_cl]) & edaph_z_all[idx_cl] > 2.5) | (mv$is_anth_dv[idx_cl] & is.finite(ddoy_all[idx_cl]) & abs(ddoy_all[idx_cl]) >= mv$tau_cryptic_vec[idx_cl]) | cy_cl
    ev_m <- is.finite(edaph_z_all[idx_cl]) | !is.na(mv$d_ed1_dv[idx_cl]) | !is.na(mv$d_ed2_dv[idx_cl]) | !is.na(mv$d_ed3_dv[idx_cl]) | (mv$is_anth_dv[idx_cl] & is.finite(ddoy_all[idx_cl])) | cy_cl
    div_f <- if (any(ev_m)) mean(div_v[ev_m], na.rm = TRUE) else mean(div_v, na.rm = TRUE); bf_m <- stats::median(bf_all[idx_cl], na.rm = TRUE)
    pass_div <- is.finite(div_f) && (div_f >= 0.50) && is.finite(bf_m) && (bf_m >= (if (n < 50L) 6.0 else 10.0))

    non_clinal_ok <- TRUE; r2_cl <- NA_real_; idx_core <- which(sp_clean_vec == sp_k & g_cl_all == 1L & t1_pure_all)
    if (length(idx_cl) >= 5L && length(idx_core) >= 5L) {
      ll_cl <- idx_cl[is.finite(lat_v_cl[idx_cl]) & is.finite(lon_v_cl[idx_cl])]; ll_co <- idx_core[is.finite(lat_v_cl[idx_core]) & is.finite(lon_v_cl[idx_core])]
      min_km <- if (length(ll_cl) >= 3L && length(ll_co) >= 3L) min(vapply(ll_cl, function(ic1) min(.haversine_km_step07(lat_v_cl[ic1], lon_v_cl[ic1], lat_v_cl[ll_co], lon_v_cl[ll_co])), numeric(1))) else Inf
      parap <- (is.finite(min_km) && min_km <= 50.0) || any(nzchar(bg_ctx$state_cnty_geo[idx_cl]) & (bg_ctx$state_cnty_geo[idx_cl] %in% bg_ctx$state_cnty_geo[idx_core]))
      anth_cl <- idx_cl[mv$is_anth_dv[idx_cl] & is.finite(ddoy_all[idx_cl])]
      allo_j <- (length(anth_cl) >= 3L) && (stats::median(abs(ddoy_all[anth_cl]), na.rm = TRUE) >= stats::median(mv$tau_cryptic_vec[anth_cl], na.rm = TRUE))
      if (!parap && !allo_j && !(mean(azonal_flag_all[idx_cl], na.rm = TRUE) >= 0.50) && length(ll_cl) >= 4L && length(ll_co) >= 4L) {
        comb <- c(ll_co, ll_cl)
        Xc <- if (!is.null(bg_ctx$morph_X)) bg_ctx$morph_X[comb, , drop = FALSE] else if (length(bg_ctx$morph_coord_cols) >= 1L) as.matrix(df[comb, bg_ctx$morph_coord_cols, drop = FALSE]) else NULL
        if (!is.null(Xc) && ncol(Xc) >= 1L && sum(ccm <- stats::complete.cases(Xc)) >= 8L) {
          Xo <- Xc[ccm, , drop = FALSE]; ro <- comb[ccm]; is_ok <- ro %in% idx_cl
          if (sum(is_ok) >= 3L && sum(!is_ok) >= 3L) {
            mu0 <- colMeans(Xo[!is_ok, , drop = FALSE]); dv <- colMeans(Xo[is_ok, , drop = FALSE]) - mu0; Sc <- stats::cov(Xo) + diag(1e-3, ncol(Xo))
            Sic <- tryCatch(solve(Sc), error = function(e) diag(1.0 / pmax(diag(Sc), 1e-3), ncol(Xo))); proj_c <- as.numeric(sweep(Xo, 2L, mu0, "-") %*% as.numeric(Sic %*% dv))
            if (stats::sd(proj_c) > 1e-8) {
              fit_cl <- tryCatch(stats::lm(proj_c ~ as.numeric(scale(lat_v_cl[ro])) + as.numeric(scale(lon_v_cl[ro]))), error = function(e) NULL)
              r2_cl <- if (!is.null(fit_cl)) summary(fit_cl)$r.squared else NA_real_
              has_gap <- as.numeric(stats::quantile(proj_c[is_ok], 0.10, names = FALSE)) > (as.numeric(stats::quantile(proj_c[!is_ok], 0.90, names = FALSE)) + 0.50 * stats::sd(proj_c[!is_ok]))
              if (is.finite(r2_cl) && ((r2_cl >= 0.60 && !has_gap) || r2_cl >= 0.72)) non_clinal_ok <- FALSE
            }
          }
        }
      }
    }
    bas_cl <- trimws(ifelse(is.na(mv$nt_bas_dv[idx_cl]), "", mv$nt_bas_dv[idx_cl])); hm_cl <- nzchar(bas_cl)
    if (any(hm_cl) && exists(".is_homotypic_sensu_stricto_basionym", mode = "function")) hm_cl[hm_cl] <- !vapply(bas_cl[hm_cl], function(bs) isTRUE(.is_homotypic_sensu_stricto_basionym(bs, sp_k)[1L]), logical(1))
    if (any(hm_cl)) {
      tb <- names(sort(table(bas_cl[hm_cl]), decreasing = TRUE))[1L]; tbi <- idx_cl[hm_cl & bas_cl == tb]
      d2b <- ifelse(is.finite(mv$bas_d2_dv[tbi]), mv$bas_d2_dv[tbi], mv$tm_d2_dv[tbi]); d2bf <- d2b[is.finite(d2b)]
      if ((length(tbi) / length(idx_cl)) >= 0.50 && ((length(d2bf) == 0L) || (stats::median(d2bf) <= mv$type_morph_chisq95))) cluster_resurrect_bas_map[[mk]] <- tb
    }
    cryptic_non_clinal_map[[mk]] <- isTRUE(non_clinal_ok); clinal_r2_map[[mk]] <- r2_cl; cryptic_cluster_pass_map[[mk]] <- isTRUE(pass_events && pass_osr && pass_div && non_clinal_ok)
  }
  list(sp_clean_vec = sp_clean_vec, g_cl_all = g_cl_all, t1_pure_all = t1_pure_all, edaph_z_all = edaph_z_all, ed_sp_all = ed_sp_all, d_edaph_all = d_edaph_all,
       ddoy_all = ddoy_all, bf_all = bf_all, spatial_deme_vec = spatial_deme_vec, coll_surname_all = coll_surname_all, clean_peer_df = clean_peer_df,
       azonal_flag_all = azonal_flag_all, cryptic_cluster_pass_map = cryptic_cluster_pass_map, cryptic_non_clinal_map = cryptic_non_clinal_map,
       clinal_r2_map = clinal_r2_map, cluster_resurrect_bas_map = cluster_resurrect_bas_map)
}
