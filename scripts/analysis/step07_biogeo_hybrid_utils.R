if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.haversine_km_step07 <- function(lat1, lon1, lat2_vec, lon2_vec) {
  dlat <- (lat2_vec - lat1) * (pi / 180.0); dlon <- (lon2_vec - lon1) * (pi / 180.0)
  a_hv <- sin(dlat * 0.5)^2 + cos(lat1 * (pi / 180.0)) * cos(lat2_vec * (pi / 180.0)) * sin(dlon * 0.5)^2
  6371.0088 * 2.0 * atan2(sqrt(pmax(0.0, a_hv)), sqrt(pmax(0.0, 1.0 - a_hv)))
}

.circ_arc_days_step07 <- function(theta1, theta2) abs(atan2(sin(theta1 - theta2), cos(theta1 - theta2))) * (365.25 / (2.0 * pi))

.resolve_db_syntopy_taxa <- function(db_rows) {
  if (is.null(db_rows) || nrow(db_rows) == 0L) return(character(0))
  sc_nm <- if ("species_standardized" %in% names(db_rows)) as.character(db_rows$species_standardized) else rep("", nrow(db_rows))
  sn_nm <- if ("scientificName" %in% names(db_rows)) as.character(db_rows$scientificName) else rep("", nrow(db_rows))
  sr_nm <- if ("species_raw" %in% names(db_rows)) as.character(db_rows$species_raw) else rep("", nrow(db_rows))
  raw_nm <- ifelse(!is.na(sc_nm) & nzchar(trimws(sc_nm)), sc_nm, ifelse(!is.na(sn_nm) & nzchar(trimws(sn_nm)), sn_nm, sr_nm))
  standardize_packera_taxon(raw_nm)
}

if (!exists(".parse_step07_semicolon_cols", mode = "function")) {
  .parse_step07_semicolon_cols <- function(df, col_name) {
    if (!(col_name %in% names(df))) return(character(0))
    rv <- as.character(stats::na.omit(df[[col_name]])); rv <- trimws(rv[nzchar(trimws(rv))])
    if (length(rv) == 0L) return(character(0))
    toks <- unique(trimws(unlist(strsplit(rv, ";", fixed = TRUE))))
    intersect(toks[nzchar(toks) & !(tolower(toks) %in% c("na", "none", "null"))], names(df))
  }
}

if (!exists(".sort_numbered_pc_cols", mode = "function")) {
  .sort_numbered_pc_cols <- function(cols, prefix) {
    if (length(cols) <= 1L) return(cols)
    cols[order(suppressWarnings(as.integer(sub(paste0("^", prefix), "", cols))), cols)]
  }
}

.resolve_step07_edaphic_predictors <- function(df) {
  exclude_cols <- c(
    "edaphic_z_score", "climatic_z_score", "edaphic_niche_overlap_d", "edaphic_schoener_d",
    "schoener_d", "sgmc_lithology_class", "Dominant_Lithology_Class", "edaphic_best_fit_taxon",
    "rf_predicted_taxon", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors",
    "type_edaphic_z_score", "nearest_type_niche_taxon", "matches_type_lithology"
  )
  keep_num <- function(cands) {
    cs <- setdiff(unique(cands), exclude_cols)
    cs[vapply(cs, function(cc) cc %in% names(df) && !is.logical(df[[cc]]) && any(is.finite(suppressWarnings(as.numeric(df[[cc]])))), logical(1))]
  }
  p_ed <- keep_num(.parse_step07_semicolon_cols(df, "primary_edaph_predictors")); if (length(p_ed) > 0L) return(p_ed)
  p_en <- keep_num(.parse_step07_semicolon_cols(df, "primary_env_predictors")); if (length(p_en) > 0L) return(p_en)
  p_fn <- keep_num(grep("^(fine_|ilr_texture_|ssurgo_|solus_|topo_)", names(df), value = TRUE)); if (length(p_fn) > 0L) return(p_fn)
  keep_num(grep("^(soil_|bio[0-9]+_)", names(df), value = TRUE))
}

.compute_dynamic_edaphic_z_scores <- function(df, sp_all_vec, active_taxa, ref_idx_by_taxon, clean_coh_mask = NULL) {
  n <- nrow(df)
  out_z <- if ("edaphic_z_score" %in% names(df)) suppressWarnings(as.numeric(df$edaphic_z_score)) else rep(NA_real_, n)
  out_z[!is.finite(out_z)] <- NA_real_; need_idx <- which(is.na(out_z))
  if (length(need_idx) == 0L) return(round(out_z, 4))
  pred_cols <- .resolve_step07_edaphic_predictors(df); if (length(pred_cols) == 0L) return(round(out_z, 4))
  p_env <- length(pred_cols); E_mat <- matrix(NA_real_, nrow = n, ncol = p_env, dimnames = list(NULL, pred_cols))
  for (j in seq_len(p_env)) E_mat[, j] <- suppressWarnings(as.numeric(df[[pred_cols[j]]]))
  noise_m <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
  coh_mask <- if (!is.null(clean_coh_mask) && length(clean_coh_mask) == n) (clean_coh_mask & !noise_m) else (!noise_m)
  robust_disp <- function(v) {
    v <- v[is.finite(v)]; if (length(v) < 2L) return(NA_real_)
    mv <- stats::mad(v, constant = 1.4826); sv <- stats::sd(v)
    if (is.finite(mv) && mv > 1e-6) (if (is.finite(sv)) max(mv, 0.25 * sv) else mv) else if (is.finite(sv) && sv > 1e-6) sv else NA_real_
  }
  coh_center <- numeric(p_env); coh_disp <- numeric(p_env)
  for (j in seq_len(p_env)) {
    vj <- E_mat[coh_mask & is.finite(E_mat[, j]), j]
    if (length(vj) < 3L) vj <- E_mat[!noise_m & is.finite(E_mat[, j]), j]
    if (length(vj) == 0L) vj <- E_mat[is.finite(E_mat[, j]), j]
    coh_center[j] <- if (length(vj) > 0L) stats::median(vj) else NA_real_; coh_disp[j] <- robust_disp(vj)
  }
  tx_center <- matrix(NA_real_, nrow = length(active_taxa), ncol = p_env, dimnames = list(active_taxa, pred_cols))
  tx_disp <- matrix(NA_real_, nrow = length(active_taxa), ncol = p_env, dimnames = list(active_taxa, pred_cols))
  for (tx in active_taxa) {
    r_idx <- ref_idx_by_taxon[[tx]]; clean_t <- which(sp_all_vec == tx & !noise_m)
    for (j in seq_len(p_env)) {
      v_ref <- E_mat[r_idx[is.finite(E_mat[r_idx, j])], j]; if (length(v_ref) == 0L) v_ref <- E_mat[clean_t[is.finite(E_mat[clean_t, j])], j]
      tx_center[tx, j] <- if (length(v_ref) > 0L) stats::median(v_ref) else coh_center[j]
      dt <- if (length(v_ref) >= 3L) { d <- robust_disp(v_ref); if (is.finite(d) && is.finite(coh_disp[j])) max(d, 0.15 * coh_disp[j]) else d } else NA_real_
      tx_disp[tx, j] <- if (is.finite(dt)) dt else coh_disp[j]
    }
  }
  for (i in need_idx) {
    xi <- E_mat[i, ]; if (!any(is.finite(xi))) { out_z[i] <- NA_real_; next }
    gsp <- sp_all_vec[i]; mu_i <- if (!is.na(gsp) && gsp %in% rownames(tx_center)) tx_center[gsp, ] else coh_center
    sd_i <- if (!is.na(gsp) && gsp %in% rownames(tx_disp)) tx_disp[gsp, ] else coh_disp
    ok_j <- is.finite(xi) & is.finite(mu_i) & is.finite(sd_i) & (sd_i > 1e-8)
    out_z[i] <- if (any(ok_j)) sqrt(mean(((xi[ok_j] - mu_i[ok_j]) / sd_i[ok_j])^2)) else {
      fb <- is.finite(xi) & is.finite(mu_i); if (any(fb) && all(abs(xi[fb] - mu_i[fb]) <= 1e-8)) 0.0 else NA_real_
    }
  }
  round(out_z, 4)
}

.parse_cytotype_x_vals <- function(s) {
  if (is.null(s) || length(s) == 0L) return(integer(0))
  txt <- as.character(s); txt <- txt[!is.na(txt) & nzchar(trimws(txt)) & !(toupper(trimws(txt)) %in% c("NA", "NAN", "NONE", "NULL"))]
  if (length(txt) == 0L) return(integer(0))
  out <- integer(0)
  for (item in txt) {
    m <- regmatches(item, gregexpr("([0-9]+)[[:space:]]*x($|[^A-Za-z0-9_])", item, ignore.case = TRUE))[[1L]]
    if (length(m) > 0L) {
      vals <- suppressWarnings(as.integer(sub("[[:space:]]*x.*$", "", tolower(m))))
      out <- c(out, vals[!is.na(vals) & vals >= 2L & vals <= 16L])
    }
  }
  sort(unique(out))
}

.infer_taxon_ploidy_levels <- function(tx, sp_vec, vp_vec, tx_cyto_vec) {
  if (is.null(tx) || is.na(tx) || !nzchar(trimws(tx))) return(integer(0))
  idx <- which(sp_vec == tx); if (length(idx) == 0L) return(integer(0))
  vp <- vp_vec[idx]; vp <- vp[!is.na(vp) & is.finite(vp) & vp >= 2L]
  if (length(vp) > 0L) return(as.integer(round(stats::median(vp))))
  .parse_cytotype_x_vals(tx_cyto_vec[idx])
}

.hybrid_cytotype_annotation <- function(vp_i, tax_cyto_i, p_a, p_b, sp_vec, vp_vec, tx_cyto_vec) {
  is_odd_vp <- !is.na(vp_i) && (vp_i >= 3L) && ((vp_i %% 2L) == 1L)
  has_odd_str <- !is.na(tax_cyto_i) && nzchar(trimws(tax_cyto_i)) && grepl("(^|[^0-9])[3579][[:space:]]*x($|[^A-Za-z0-9_])", tax_cyto_i, ignore.case = TRUE)
  pa_l <- .infer_taxon_ploidy_levels(p_a, sp_vec, vp_vec, tx_cyto_vec); pb_l <- .infer_taxon_ploidy_levels(p_b, sp_vec, vp_vec, tx_cyto_vec)
  both_k <- (length(pa_l) > 0L) && (length(pb_l) > 0L); both_2x <- both_k && all(pa_l == 2L) && all(pb_l == 2L)
  if (is_odd_vp || has_odd_str || (both_k && !identical(sort(unique(pa_l)), sort(unique(pb_l)))) || (both_2x && !is.na(vp_i) && vp_i != 2L)) {
    return(" [Heteroploid / Odd-Ploidy (3x/5x) or Allopolyploid Candidate]")
  }
  if (both_2x && (is.na(vp_i) || vp_i == 2L)) return(" [Homoploid Hybrid Candidate: shared 2x]")
  ""
}

.resolve_cytotype_barrier_pair <- function(i, sp_i, g_cl_i, cat_i, pid_i, barrier_attr, sp_vec, g_cl_vec, vp_vec, tx_cyto_vec, is_type_vec, is_typ_morph_vec) {
  raw_lbl <- NA_character_
  if (!is.null(barrier_attr)) {
    if (is.data.frame(barrier_attr) && nrow(barrier_attr) > 0L) {
      hit <- integer(0)
      if (all(c("catalogNumber", "plant_individual_id") %in% names(barrier_attr))) hit <- which(as.character(barrier_attr$catalogNumber) == cat_i & as.integer(barrier_attr$plant_individual_id) == pid_i)
      if (length(hit) == 0L && all(c("species_standardized", "gmm_cluster") %in% names(barrier_attr))) hit <- which(as.character(barrier_attr$species_standardized) == sp_i & as.integer(barrier_attr$gmm_cluster) == g_cl_i)
      if (length(hit) > 0L) { lcol <- intersect(c("barrier_label", "label", "cytotype_barrier_label"), names(barrier_attr)); if (length(lcol) > 0L) raw_lbl <- as.character(barrier_attr[[lcol[1L]]][hit[1L]]) }
    } else if (is.character(barrier_attr) || is.list(barrier_attr)) {
      nm_b <- names(barrier_attr); keys <- c(paste0(sp_i, "||", g_cl_i), paste0(cat_i, "__", pid_i), cat_i, paste0(sp_i, "__K", g_cl_i), paste0(sp_i, "__", g_cl_i), sp_i)
      if (!is.null(nm_b) && any(keys %in% nm_b)) {
        val <- barrier_attr[[keys[keys %in% nm_b][1L]]]
        if (is.list(val) && !is.null(val$cl_lbl) && !is.null(val$type_lbl)) return(list(cl_lbl = as.character(val$cl_lbl), type_lbl = as.character(val$type_lbl)))
        raw_lbl <- as.character(val)[1L]
      } else if (is.character(barrier_attr) && length(barrier_attr) >= i) { cv <- trimws(as.character(barrier_attr[i])); if (!is.na(cv) && nzchar(cv)) raw_lbl <- cv }
    }
  }
  if (!is.na(raw_lbl) && grepl("[[:space:]]+vs[[:space:]]+", raw_lbl)) {
    pts <- trimws(strsplit(raw_lbl, "[[:space:]]+vs[[:space:]]+")[[1L]])
    if (length(pts) >= 2L && nzchar(pts[1L]) && nzchar(pts[2L])) {
      return(list(cl_lbl = sub("^cluster ploidy[[:space:]]+", "", pts[1L], ignore.case = TRUE), type_lbl = sub("^type-cluster[[:space:]]+", "", pts[2L], ignore.case = TRUE)))
    }
  }
  cl_vp <- vp_vec[sp_vec == sp_i & g_cl_vec == g_cl_i]; cl_vp <- cl_vp[!is.na(cl_vp) & is.finite(cl_vp)]
  cl_p <- if (length(cl_vp) > 0L) as.integer(round(stats::median(cl_vp))) else if (!is.na(vp_vec[i]) && is.finite(vp_vec[i])) as.integer(round(vp_vec[i])) else NA_integer_
  type_vp <- vp_vec[sp_vec == sp_i & is_type_vec & !is.na(vp_vec) & is.finite(vp_vec)]
  if (length(type_vp) == 0L) type_vp <- vp_vec[sp_vec == sp_i & (g_cl_vec == 1L | (!is.na(is_typ_morph_vec) & is_typ_morph_vec)) & !is.na(vp_vec) & is.finite(vp_vec)]
  if (length(type_vp) == 0L) type_vp <- vp_vec[sp_vec == sp_i & (g_cl_vec != g_cl_i) & !is.na(vp_vec) & is.finite(vp_vec)]
  type_p <- if (length(type_vp) > 0L) as.integer(round(stats::median(type_vp))) else if (!is.na(cl_p) && !is.na(tx_cyto_vec[i])) { oth <- setdiff(.parse_cytotype_x_vals(tx_cyto_vec[i]), cl_p); if (length(oth) == 1L) oth[1L] else NA_integer_ } else NA_integer_
  list(cl_lbl = if (!is.na(cl_p)) sprintf("%dx", cl_p) else "heteroploid", type_lbl = if (!is.na(type_p)) sprintf("%dx", type_p) else "type-cluster")
}

.format_infra_cytotype_lbl <- function(raw_sp, given_sp, g_cl) {
  if (grepl("(^|[^A-Za-z0-9_])(var\\.|subsp\\.|ssp\\.)[[:space:]]+[a-z-]+", raw_sp, ignore.case = TRUE)) {
    m_inf <- regmatches(raw_sp, regexpr("(var\\.|subsp\\.|ssp\\.)[[:space:]]+[A-Za-z-]+", raw_sp, ignore.case = TRUE))
    if (length(m_inf) > 0L && nzchar(m_inf[1L])) return(paste(given_sp, tolower(m_inf[1L])))
  }
  sprintf("%s var. cytotype_K%d", given_sp, max(1L, as.integer(g_cl)))
}

.build_stage4_cytotype_decision <- function(given_sp, raw_sp, g_cl, base_rationale) {
  suffix <- " [Intraspecific Cytotype Series: FCM/squash verification recommended]"
  list(status = "Subspecies_Ecotype", broad = "Variety / Ecotype", sub = "Subspecies_Ecotype", prio = "MEDIUM", act = "Accept_Subspecific_Treatment",
       rec = .format_infra_cytotype_lbl(raw_sp, given_sp, g_cl),
       rat = if (grepl("Intraspecific Cytotype Series:", base_rationale, fixed = TRUE)) base_rationale else paste0(base_rationale, suffix))
}

.US_ABBR_TO_STATE_STEP07 <- c(
  AL="alabama",AK="alaska",AZ="arizona",AR="arkansas",CA="california",CO="colorado",CT="connecticut",DE="delaware",DC="district of columbia",FL="florida",
  GA="georgia",HI="hawaii",ID="idaho",IL="illinois",IN="indiana",IA="iowa",KS="kansas",KY="kentucky",LA="louisiana",ME="maine",MD="maryland",MA="massachusetts",
  MI="michigan",MN="minnesota",MS="mississippi",MO="missouri",MT="montana",NE="nebraska",NV="nevada",NH="new hampshire",NJ="new jersey",NM="new mexico",NY="new york",
  NC="north carolina",ND="north dakota",OH="ohio",OK="oklahoma",OR="oregon",PA="pennsylvania",RI="rhode island",SC="south carolina",SD="south dakota",TN="tennessee",
  TX="texas",UT="utah",VT="vermont",VA="virginia",WA="washington",WV="west virginia",WI="wisconsin",WY="wyoming",AB="alberta",BC="british columbia",MB="manitoba",
  NB="new brunswick",NL="newfoundland and labrador",NS="nova scotia",NT="northwest territories",NU="nunavut",ON="ontario",PE="prince edward island",QC="quebec",SK="saskatchewan",YT="yukon"
)

.norm_state_step07 <- function(sv) {
  s <- tolower(trimws(as.character(sv))); s[is.na(s) | s %in% c("", "na", "nan", "none", "null", "unknown")] <- ""
  up <- toupper(s); hit <- nzchar(s) & (up %in% names(.US_ABBR_TO_STATE_STEP07))
  if (any(hit)) s[hit] <- unname(.US_ABBR_TO_STATE_STEP07[up[hit]])
  s
}

.norm_county_step07 <- function(cv) {
  s <- tolower(trimws(as.character(cv))); s[is.na(s) | s %in% c("", "na", "nan", "none", "null", "unknown")] <- ""
  s <- gsub("(^|[^A-Za-z0-9_])(county|parish|co\\.|borough)($|[^A-Za-z0-9_])", " ", s)
  trimws(gsub("[^a-z0-9]+", " ", s))
}

.is_valid_packera_taxon_step07 <- function(s, active_taxa) !is.null(s) && !is.na(s) && nzchar(s) && s != "Unknown" && (s %in% active_taxa || grepl("^(Packera|P\\.)[[:space:]]+", s))

.eval_taxon_geo_plausibility_ctx <- function(row_i, cand_tx, bg_ctx, is_intra) {
  if (!.is_valid_packera_taxon_step07(cand_tx, bg_ctx$active_taxa)) return(list(plausible = FALSE, weight = 0.0, dist_km = Inf))
  st_i <- bg_ctx$state_vec_geo[row_i]; lat_i <- bg_ctx$lat_vec_geo[row_i]; lon_i <- bg_ctx$lon_vec_geo[row_i]
  has_st_i <- nzchar(st_i); has_ll_i <- is.finite(lat_i) && is.finite(lon_i)
  if (!has_st_i && !has_ll_i) return(if (is_intra && bg_ctx$has_ext_geo_ref) list(plausible = FALSE, weight = 0.10, dist_km = Inf) else list(plausible = TRUE, weight = 1.0, dist_km = NA_real_))
  known_sts <- unique(c(bg_ctx$anch_states_by_tx[[cand_tx]], bg_ctx$usda_states_by_tx[[cand_tx]]))
  c_mat <- bg_ctx$anch_coords_by_tx[[cand_tx]]
  min_km <- if (has_ll_i && !is.null(c_mat) && nrow(c_mat) > 0L) min(.haversine_km_step07(lat_i, lon_i, c_mat[, 1L], c_mat[, 2L])) else NA_real_
  if (has_st_i && length(known_sts) > 0L && (st_i %in% known_sts) && !(has_ll_i && is.finite(min_km) && min_km > 450.0)) {
    return(list(plausible = TRUE, weight = 1.0, dist_km = if (is.finite(min_km)) min_km else NA_real_))
  }
  if (has_ll_i && is.finite(min_km)) {
    max_km <- if (nrow(c_mat) >= 4L || length(bg_ctx$usda_states_by_tx[[cand_tx]]) > 0L) 450.0 else 650.0
    if (min_km <= max_km) return(list(plausible = TRUE, weight = max(0.35, if (min_km <= 200.0) 1.0 else exp(-0.5 * ((min_km - 200.0) / 250.0)^2)), dist_km = min_km))
    return(list(plausible = FALSE, weight = max(0.01, exp(-0.5 * ((min_km - 200.0) / 200.0)^2)), dist_km = min_km))
  }
  if (has_st_i && length(known_sts) >= 2L && !(st_i %in% known_sts)) return(list(plausible = FALSE, weight = 0.05, dist_km = 800.0))
  if (is_intra && bg_ctx$has_ext_geo_ref) return(list(plausible = FALSE, weight = 0.05, dist_km = Inf))
  list(plausible = TRUE, weight = 0.85, dist_km = NA_real_)
}

.resolve_nearest_mahal_ctx <- function(row_idx, given_sp, cohort_species, bg_ctx) {
  other_tx <- setdiff(cohort_species, given_sp); if (length(other_tx) == 0L) return(NA_character_)
  if (length(other_tx) == 1L || is.null(bg_ctx$morph_X) || is.null(bg_ctx$morph_centroids) || is.null(bg_ctx$morph_S_inv)) return(other_tx[1L])
  other_m <- intersect(other_tx, rownames(bg_ctx$morph_centroids)); if (length(other_m) == 0L) return(other_tx[1L])
  x_i <- bg_ctx$morph_X[row_idx, ]; base_mu <- if (given_sp %in% rownames(bg_ctx$morph_centroids)) bg_ctx$morph_centroids[given_sp, ] else colMeans(bg_ctx$morph_centroids)
  x_i[!is.finite(x_i)] <- base_mu[!is.finite(x_i)]
  d2_v <- vapply(other_m, function(otx) { dv <- x_i - bg_ctx$morph_centroids[otx, ]; max(0.0, as.numeric(t(dv) %*% bg_ctx$morph_S_inv %*% dv)) }, numeric(1))
  other_m[which.min(d2_v)]
}

.build_step07_biogeo_and_morph_context <- function(df, opts, sp_all_vec, cat_sheet_vec, ref_idx_by_taxon, active_taxa, is_intra, active_geno_ref, t0_info,
                                                   valid_doy_all, anth_all, hop_doy_all, anchor_all, mo_all, noise_all) {
  n <- nrow(df); is_cent_geo <- if (exists(".is_step07_centroid_coord", mode = "function")) .is_step07_centroid_coord(df) else rep(FALSE, n)
  state_vec_geo <- if ("stateProvince" %in% names(df)) .norm_state_step07(df$stateProvince) else rep("", n)
  cnty_raw_geo <- if ("county" %in% names(df)) .norm_county_step07(df$county) else rep("", n)
  state_cnty_geo <- ifelse(nzchar(state_vec_geo) & nzchar(cnty_raw_geo), paste(state_vec_geo, cnty_raw_geo, sep = "||"), "")
  lat_vec_geo <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, n)
  lon_vec_geo <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, n)
  lat_vec_geo[is_cent_geo] <- NA_real_; lon_vec_geo[is_cent_geo] <- NA_real_
  cohort_species <- sort(unique(trimws(as.character(sp_all_vec))[!is.na(sp_all_vec) & grepl("^Packera[[:space:]]+[a-z][a-z-]+$", trimws(as.character(sp_all_vec)))]))
  can_evaluate_hybrids <- !is_intra && (length(cohort_species) >= 2L)
  has_any_geo <- any(nzchar(state_vec_geo)) || any(is.finite(lat_vec_geo) & is.finite(lon_vec_geo))
  use_ext_db <- nzchar(Sys.getenv("RUN_ID", unset = "")) || nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")) || has_any_geo
  usda_states_by_tx <- list(); usda_counties_by_tx <- list(); usda_state_county_counts_by_tx <- list(); db_syntopy_rows <- NULL
  if (use_ext_db && trimws(env_plantmatics_or_packera("PLANTMATICS_DEMO_MODE", "PACKERA_DEMO_MODE")) != "1" && exists("connect_packera_db", mode = "function") && requireNamespace("DBI", quietly = TRUE)) {
    tryCatch({
      c_geo <- connect_packera_db(read_only = TRUE, demo = FALSE)
      on.exit(tryCatch(if (!is.null(c_geo) && DBI::dbIsValid(c_geo)) DBI::dbDisconnect(c_geo), error = function(e) NULL), add = TRUE)
      if (!is.null(c_geo) && DBI::dbExistsTable(c_geo, "usda_county_ranges")) {
        uf <- tryCatch(DBI::dbListFields(c_geo, "usda_county_ranges"), error = function(e) character(0))
        u_df <- DBI::dbGetQuery(c_geo, paste0("SELECT DISTINCT canonical_packera, state_norm", if ("county_norm" %in% uf) ", county_norm" else "",
          if ("is_state_only" %in% uf) ", is_state_only" else "", if ("usda_symbol" %in% uf) ", usda_symbol" else "", " FROM usda_county_ranges WHERE state_norm IS NOT NULL AND TRIM(state_norm) != '';"))
        if (!is.null(u_df) && nrow(u_df) > 0L) {
          u_tx <- trimws(as.character(u_df$canonical_packera)); u_st <- .norm_state_step07(u_df$state_norm); ok_u <- !is.na(u_tx) & nzchar(u_tx) & nzchar(u_st)
          if ("usda_symbol" %in% names(u_df) && any(ok_u)) {
            u_sym <- toupper(trimws(ifelse(is.na(u_df$usda_symbol), "", as.character(u_df$usda_symbol)))); hs <- ok_u & nzchar(u_sym)
            if (any(hs) && length(sh_s <- names(which(tapply(u_tx[hs], u_sym[hs], function(v) length(unique(v))) > 1L))) > 0L) {
              ok_u <- ok_u & (!(hs & (u_sym %in% sh_s)) | (nchar(u_sym) >= 4L & toupper(substr(sub("^[A-Za-z-]+[[:space:]]+", "", u_tx), 1L, 2L)) == substr(u_sym, 3L, 4L)))
            }
          }
          if (any(ok_u)) {
            usda_states_by_tx <- lapply(split(u_st[ok_u], u_tx[ok_u]), unique)
            if ("county_norm" %in% names(u_df)) {
              u_cty <- .norm_county_step07(u_df$county_norm); u_so <- if ("is_state_only" %in% names(u_df)) (!is.na(u_df$is_state_only) & suppressWarnings(as.integer(u_df$is_state_only)) != 0L) else rep(FALSE, nrow(u_df))
              if (any(ok_uc <- ok_u & nzchar(u_cty) & !u_so)) {
                usda_counties_by_tx <- lapply(split(paste(u_st[ok_uc], u_cty[ok_uc], sep = "||"), u_tx[ok_uc]), unique)
                for (tx_c in names(usda_counties_by_tx)) usda_state_county_counts_by_tx[[tx_c]] <- tapply(u_cty[ok_uc & u_tx == tx_c], u_st[ok_uc & u_tx == tx_c], function(v) length(unique(v)))
              }
            }
          }
        }
      }
      if (!is.null(c_geo) && DBI::dbExistsTable(c_geo, "natureserve_subnational_ranges")) {
        ns_df <- tryCatch(DBI::dbGetQuery(c_geo, "SELECT DISTINCT species_standardized, subnation_code FROM natureserve_subnational_ranges WHERE subnation_code IS NOT NULL AND TRIM(subnation_code) != '' AND COALESCE(is_native, 0) = 1 AND COALESCE(is_exotic, 0) = 0;"), error = function(e) NULL)
        if (!is.null(ns_df) && nrow(ns_df) > 0L) {
          ns_tx <- trimws(as.character(ns_df$species_standardized)); ns_st <- .norm_state_step07(ns_df$subnation_code); ok_ns <- !is.na(ns_tx) & nzchar(ns_tx) & nzchar(ns_st)
          for (txn in unique(ns_tx[ok_ns])) usda_states_by_tx[[txn]] <- unique(c(usda_states_by_tx[[txn]], ns_st[ok_ns & ns_tx == txn]))
        }
      }
      if ((can_evaluate_hybrids || is_intra) && (nzchar(Sys.getenv("RUN_ID", unset = "")) || nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"))) && DBI::dbExistsTable(c_geo, "vouchers")) {
        vf <- tryCatch(DBI::dbListFields(c_geo, "vouchers"), error = function(e) character(0))
        q_syn <- paste(paste0("SELECT v.catalogNumber, ", if ("species_standardized" %in% vf) "v.species_standardized, " else "", if ("recordedBy" %in% vf) "v.recordedBy, " else "",
          "v.scientificName, v.species_raw, v.determiner_tier, v.stateProvince, v.county, v.decimalLatitude, v.decimalLongitude", if ("is_centroid_or_inst_coord" %in% vf) ", v.is_centroid_or_inst_coord" else "", " FROM vouchers v"),
          "WHERE v.determiner_tier IN ('Tier_0_Type_Anchor', 'Tier_1_Gold', 'Tier_2_Silver')",
          if ("usda_county_status" %in% vf) "AND (v.determiner_tier = 'Tier_0_Type_Anchor' OR COALESCE(LOWER(TRIM(v.usda_county_status)), '') NOT IN ('out_of_range', 'natureserve_exotic_or_waif'))" else "",
          if ("exsiccatae_taxonomic_conflict" %in% vf) "AND (v.determiner_tier = 'Tier_0_Type_Anchor' OR COALESCE(v.exsiccatae_taxonomic_conflict, 0) = 0)" else "",
          if ("is_centroid_or_inst_coord" %in% vf) "AND (v.determiner_tier = 'Tier_0_Type_Anchor' OR COALESCE(v.is_centroid_or_inst_coord, 0) = 0)" else "",
          "AND NOT EXISTS (SELECT 1 FROM curation_triage ct WHERE ct.catalogNumber = v.catalogNumber AND ct.taxonomic_treatment_subtype = 'Misidentification_Reassign');")
        db_syntopy_rows <- tryCatch(DBI::dbGetQuery(c_geo, q_syn), error = function(e) NULL)
      }
    }, error = function(e) NULL)
  }

  anch_states_by_tx <- list(); anch_coords_by_tx <- list()
  for (tx in active_taxa) {
    ri <- ref_idx_by_taxon[[tx]]; anch_states_by_tx[[tx]] <- unique(c(state_vec_geo[ri][nzchar(state_vec_geo[ri])], usda_states_by_tx[[tx]]))
    ok_ll <- ri[is.finite(lat_vec_geo[ri]) & is.finite(lon_vec_geo[ri])]
    anch_coords_by_tx[[tx]] <- if (length(ok_ll) > 0L) cbind(lat = lat_vec_geo[ok_ll], lon = lon_vec_geo[ok_ll]) else NULL
  }
  db_tx_vec <- character(0); db_lat_vec <- numeric(0); db_lon_vec <- numeric(0)
  if (!is.null(db_syntopy_rows) && nrow(db_syntopy_rows) > 0L) {
    db_tx_vec <- .resolve_db_syntopy_taxa(db_syntopy_rows); is_cent_db <- if (exists(".is_step07_centroid_coord", mode = "function")) .is_step07_centroid_coord(db_syntopy_rows) else rep(FALSE, nrow(db_syntopy_rows))
    st_db <- .norm_state_step07(db_syntopy_rows$stateProvince)
    db_lat_vec <- suppressWarnings(as.numeric(db_syntopy_rows$decimalLatitude)); db_lon_vec <- suppressWarnings(as.numeric(db_syntopy_rows$decimalLongitude))
    db_lat_vec[is_cent_db] <- NA_real_; db_lon_vec[is_cent_db] <- NA_real_
    for (tx_u in unique(db_tx_vec[!is.na(db_tx_vec) & nzchar(db_tx_vec) & db_tx_vec != "Unknown"])) {
      mu <- !is.na(db_tx_vec) & (db_tx_vec == tx_u)
      anch_states_by_tx[[tx_u]] <- unique(c(anch_states_by_tx[[tx_u]], st_db[mu & nzchar(st_db)], usda_states_by_tx[[tx_u]]))
      ok_dll <- mu & is.finite(db_lat_vec) & is.finite(db_lon_vec)
      if (any(ok_dll)) { cm <- cbind(lat = db_lat_vec[ok_dll], lon = db_lon_vec[ok_dll]); anch_coords_by_tx[[tx_u]] <- if (is.null(anch_coords_by_tx[[tx_u]])) cm else rbind(anch_coords_by_tx[[tx_u]], cm) }
    }
  }

  max_mc_p <- max(2L, min(15L, n - length(active_taxa) - 1L))
  ok_fn <- function(cc) !is.logical(df[[cc]]) && sum(is.finite(suppressWarnings(as.numeric(df[[cc]])))) >= 2L
  sel_mc <- Filter(ok_fn, .parse_step07_semicolon_cols(df, "selected_morph_features"))
  can_c <- Filter(ok_fn, intersect(c("can1", "can2", "can3"), names(df)))
  raw_pc <- Filter(ok_fn, .sort_numbered_pc_cols(grep("^PC[0-9]+$", names(df), value = TRUE), "PC"))
  pc_cand <- if (length(sel_mc) >= 2L) utils::head(sel_mc, max_mc_p) else utils::head(raw_pc, max_mc_p)
  morph_coord_cols <- if (length(active_taxa) > 4L && length(pc_cand) >= 2L) pc_cand else if (length(can_c) >= 2L) can_c else if (length(pc_cand) >= 1L) pc_cand else can_c
  morph_X <- NULL; morph_centroids <- NULL; morph_S_inv <- NULL
  if (length(morph_coord_cols) >= 1L && length(active_taxa) >= 2L) {
    pd <- length(morph_coord_cols); morph_X <- matrix(NA_real_, nrow = n, ncol = pd, dimnames = list(NULL, morph_coord_cols))
    for (j in seq_len(pd)) morph_X[, j] <- suppressWarnings(as.numeric(df[[morph_coord_cols[j]]]))
    morph_centroids <- matrix(0.0, nrow = length(active_taxa), ncol = pd, dimnames = list(active_taxa, morph_coord_cols))
    S_pool <- matrix(0.0, nrow = pd, ncol = pd); df_pool <- 0L
    for (tx in active_taxa) {
      ri <- ref_idx_by_taxon[[tx]]; sub_m <- morph_X[ri, , drop = FALSE]
      u_wt <- active_geno_ref && !is.null(t0_info) && length(ri) > 0L && any(t0_info$is_tier0_dna_anchor[ri] & (t0_info$rectified_taxon[ri] == tx))
      wr <- if (u_wt) ifelse(t0_info$is_tier0_dna_anchor[ri] & (t0_info$rectified_taxon[ri] == tx), ifelse(is.finite(t0_info$tier0_weight[ri]), t0_info$tier0_weight[ri], 4.0), 1.0) else NULL
      for (j in seq_len(pd)) {
        fm <- is.finite(sub_m[, j]); fj <- sub_m[fm, j]
        morph_centroids[tx, j] <- if (length(fj) > 0L) (if (u_wt) sum(wr[fm] * fj) / sum(wr[fm]) else mean(fj)) else { fa <- morph_X[is.finite(morph_X[, j]), j]; if (length(fa) > 0L) mean(fa) else 0.0 }
        sub_m[!fm, j] <- morph_centroids[tx, j]
      }
      if (nrow(sub_m) >= 2L) {
        c_tx <- if (u_wt && length(unique(wr)) > 1L) tryCatch(stats::cov.wt(sub_m, wt = wr, method = "unbiased")$cov, error = function(e) stats::cov(sub_m)) else stats::cov(sub_m)
        c_tx[!is.finite(c_tx)] <- 0.0; S_pool <- S_pool + (nrow(sub_m) - 1L) * c_tx; df_pool <- df_pool + (nrow(sub_m) - 1L)
      }
    }
    if (df_pool > 0L) S_pool <- S_pool / df_pool else {
      X_imp <- morph_X; for (j in seq_len(pd)) { mdj <- stats::median(X_imp[, j], na.rm = TRUE); X_imp[!is.finite(X_imp[, j]), j] <- if (is.finite(mdj)) mdj else 0.0 }
      S_pool <- if (nrow(X_imp) >= 2L) stats::cov(X_imp) else diag(1.0, nrow = pd, ncol = pd)
    }
    S_pool[!is.finite(S_pool)] <- 0.0; S_reg <- S_pool + diag(max(1e-4, 0.05 * mean(diag(S_pool), na.rm = TRUE)), nrow = pd, ncol = pd)
    morph_S_inv <- tryCatch(solve(S_reg), error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-4), nrow = pd, ncol = pd))
  }
  self_morph_S_inv <- if (!is.null(morph_S_inv)) tryCatch({ S_e <- solve(morph_S_inv); diag(S_e) <- pmax(diag(S_e), 0.25); solve(S_e) }, error = function(e) morph_S_inv) else NULL

  organ_taxa <- unique(c(active_taxa, cohort_species)); has_anthesis_doy_by_tx <- stats::setNames(rep(FALSE, length(organ_taxa)), organ_taxa)
  sec_cols <- intersect(c("blade_cie_l_sheet_delta", "blade_cie_l_norm", "blade_cie_l_median", "AR_inv", "involucre_length_width_ratio", "log1p_capitulum_count", "petiole_to_blade_ratio"), names(df))
  sec_cols <- sec_cols[vapply(sec_cols, function(cc) sum(is.finite(suppressWarnings(as.numeric(df[[cc]])))) >= 4L, logical(1))]
  sec_mat <- if (length(sec_cols) > 0L) { m_t <- matrix(NA_real_, nrow = n, ncol = length(sec_cols), dimnames = list(NULL, sec_cols)); for (sj in seq_along(sec_cols)) m_t[, sj] <- suppressWarnings(as.numeric(df[[sec_cols[sj]]])); m_t } else NULL
  sec_med_tx <- list(); sec_mad <- stats::setNames(rep(NA_real_, length(sec_cols)), sec_cols)
  for (tx in organ_taxa) {
    ri <- ref_idx_by_taxon[[tx]]; if (length(ri) == 0L) next
    has_anthesis_doy_by_tx[[tx]] <- any(valid_doy_all[ri] & anth_all[ri] & is.finite(hop_doy_all[ri]))
    if (!is.null(sec_mat)) { mv <- stats::setNames(rep(NA_real_, ncol(sec_mat)), colnames(sec_mat)); for (sj in seq_len(ncol(sec_mat))) { vv <- sec_mat[ri, sj]; vv <- vv[is.finite(vv)]; if (length(vv) >= 2L) mv[sj] <- stats::median(vv) }; sec_med_tx[[tx]] <- mv }
  }
  if (!is.null(sec_mat)) for (sj in seq_len(ncol(sec_mat))) {
    vv <- sec_mat[anchor_all & !mo_all & !noise_all, sj]; vv <- vv[is.finite(vv)]; if (length(vv) < 4L) { vv <- sec_mat[, sj]; vv <- vv[is.finite(vv)] }
    md <- if (length(vv) >= 2L) stats::mad(vv, constant = 1.4826) else NA_real_; sdv <- if (length(vv) >= 2L) stats::sd(vv) else NA_real_
    sec_mad[sj] <- if (is.finite(md) && md > 1e-6) md else if (is.finite(sdv) && sdv > 1e-6) sdv else NA_real_
  }
  syntopy_pool_by_tx <- list()
  if (can_evaluate_hybrids) {
    db_cat <- if (length(db_tx_vec) > 0L) trimws(as.character(db_syntopy_rows$catalogNumber)) else character(0)
    db_sc <- if (length(db_tx_vec) > 0L) { st <- .norm_state_step07(db_syntopy_rows$stateProvince); cn <- .norm_county_step07(db_syntopy_rows$county); ifelse(nzchar(st) & nzchar(cn), paste(st, cn, sep = "||"), "") } else character(0)
    for (tx in cohort_species) {
      ri <- ref_idx_by_taxon[[tx]]; ct <- cat_sheet_vec[ri]; lt <- lat_vec_geo[ri]; ln <- lon_vec_geo[ri]; sc <- state_cnty_geo[ri]
      if (length(db_tx_vec) > 0L) {
        md <- which(!is.na(db_tx_vec) & db_tx_vec == tx & !(db_cat %in% ct))
        if (length(md) > 0L) { ct <- c(ct, db_cat[md]); lt <- c(lt, db_lat_vec[md]); ln <- c(ln, db_lon_vec[md]); sc <- c(sc, db_sc[md]) }
      }
      syntopy_pool_by_tx[[tx]] <- list(cats = ct, lats = lt, lons = ln, scs = sc)
    }
  }
  list(active_taxa = active_taxa, cohort_species = cohort_species, can_evaluate_hybrids = can_evaluate_hybrids, is_cent_geo = is_cent_geo,
       state_vec_geo = state_vec_geo, cnty_raw_geo = cnty_raw_geo, state_cnty_geo = state_cnty_geo, lat_vec_geo = lat_vec_geo, lon_vec_geo = lon_vec_geo,
       usda_states_by_tx = usda_states_by_tx, usda_counties_by_tx = usda_counties_by_tx, usda_state_county_counts_by_tx = usda_state_county_counts_by_tx,
       db_syntopy_rows = db_syntopy_rows, anch_states_by_tx = anch_states_by_tx, anch_coords_by_tx = anch_coords_by_tx,
       has_ext_geo_ref = (length(usda_states_by_tx) > 0L) || (!is.null(db_syntopy_rows) && nrow(db_syntopy_rows) > 0L),
       morph_X = morph_X, morph_centroids = morph_centroids, morph_S_inv = morph_S_inv, self_morph_S_inv = self_morph_S_inv, sel_mc = sel_mc,
       morph_coord_cols = morph_coord_cols, sec_organ_cols = sec_cols, sec_organ_mat = sec_mat, sec_organ_med_by_tx = sec_med_tx,
       sec_organ_mad_pool = sec_mad, has_anthesis_doy_by_tx = has_anthesis_doy_by_tx, syntopy_pool_by_tx = syntopy_pool_by_tx)
}

.eval_five_gate_hybrid_pair_impl <- function(row_i, p1_tx, p2_tx, morph_sp_i, cda_prob_i, cda_sec_tx_i, cda_sec_p_i, cda_marg_i,
                                            gmm_unc_i, edaph_sp_i, d_edaph_i, edaph_z_i, asym_i, bg_ctx, cat_sheet_vec, mu_theta_hop,
                                            valid_doy_all, anth_all, hop_doy_all) {
  fail_res <- list(pass = FALSE, d1_km = Inf, d2_km = Inf, ddoy_12 = Inf, t_ab = NA_real_, d2_perp = Inf, score = -Inf)
  if (!bg_ctx$can_evaluate_hybrids || !(p1_tx %in% bg_ctx$cohort_species) || !(p2_tx %in% bg_ctx$cohort_species) || p1_tx == p2_tx) return(fail_res)
  p1 <- bg_ctx$syntopy_pool_by_tx[[p1_tx]]; p2 <- bg_ctx$syntopy_pool_by_tx[[p2_tx]]; if (is.null(p1) || is.null(p2)) return(fail_res)
  self_cat <- cat_sheet_vec[row_i]; lat_i <- bg_ctx$lat_vec_geo[row_i]; lon_i <- bg_ctx$lon_vec_geo[row_i]; sc_i <- bg_ctx$state_cnty_geo[row_i]
  has_ll <- is.finite(lat_i) && is.finite(lon_i); has_sc <- nzchar(sc_i); if (!has_ll && !has_sc) return(fail_res)
  k1 <- which(p1$cats != self_cat); k2 <- which(p2$cats != self_cat); if (length(k1) == 0L || length(k2) == 0L) return(fail_res)
  same_sc1 <- has_sc && any(nzchar(p1$scs[k1]) & p1$scs[k1] == sc_i); same_sc2 <- has_sc && any(nzchar(p2$scs[k2]) & p2$scs[k2] == sc_i)
  d1_km <- Inf; d2_km <- Inf
  if (has_ll) {
    ll1 <- k1[is.finite(p1$lats[k1]) & is.finite(p1$lons[k1])]; ll2 <- k2[is.finite(p2$lats[k2]) & is.finite(p2$lons[k2])]
    if (length(ll1) > 0L) d1_km <- min(.haversine_km_step07(lat_i, lon_i, p1$lats[ll1], p1$lons[ll1]))
    if (length(ll2) > 0L) d2_km <- min(.haversine_km_step07(lat_i, lon_i, p2$lats[ll2], p2$lons[ll2]))
    syn_ok <- is.finite(d1_km) && is.finite(d2_km) && ((d1_km <= 15.0 && d2_km <= 15.0) || (same_sc1 && same_sc2 && d1_km <= 25.0 && d2_km <= 25.0))
  } else { syn_ok <- same_sc1 && same_sc2; if (syn_ok) { d1_km <- 0.0; d2_km <- 0.0 } }
  if (!syn_ok) return(fail_res)

  ddoy_12 <- 0.0
  if ((p1_tx %in% names(bg_ctx$has_anthesis_doy_by_tx)) && (p2_tx %in% names(bg_ctx$has_anthesis_doy_by_tx)) &&
      isTRUE(bg_ctx$has_anthesis_doy_by_tx[[p1_tx]]) && isTRUE(bg_ctx$has_anthesis_doy_by_tx[[p2_tx]]) &&
      (p1_tx %in% names(mu_theta_hop)) && (p2_tx %in% names(mu_theta_hop))) {
    th1 <- mu_theta_hop[[p1_tx]]; th2 <- mu_theta_hop[[p2_tx]]; ddoy_12 <- .circ_arc_days_step07(th1, th2)
    if (!is.finite(ddoy_12) || ddoy_12 > 21.0) return(fail_res)
    if (valid_doy_all[row_i] && anth_all[row_i] && is.finite(hop_doy_all[row_i])) {
      thi <- 2.0 * pi * (hop_doy_all[row_i] / 365.25); thm <- atan2(0.5 * (sin(th1) + sin(th2)), 0.5 * (cos(th1) + cos(th2)))
      if (min(.circ_arc_days_step07(thi, th1), .circ_arc_days_step07(thi, th2), .circ_arc_days_step07(thi, thm)) > 21.0) return(fail_res)
    }
  }
  if ((is.finite(d_edaph_i) && d_edaph_i < 0.25) || (is.finite(edaph_z_i) && edaph_z_i > 2.5)) return(fail_res)

  has_mc <- !is.null(bg_ctx$morph_X) && !is.null(bg_ctx$morph_centroids) && !is.null(bg_ctx$morph_S_inv) &&
            (p1_tx %in% rownames(bg_ctx$morph_centroids)) && (p2_tx %in% rownames(bg_ctx$morph_centroids)) && all(is.finite(bg_ctx$morph_X[row_i, ]))
  if (has_mc) {
    mu1 <- bg_ctx$morph_centroids[p1_tx, ]; d12 <- bg_ctx$morph_centroids[p2_tx, ] - mu1
    d2_12 <- max(0.0, as.numeric(t(d12) %*% bg_ctx$morph_S_inv %*% d12)); if (!is.finite(d2_12) || d2_12 < 2.25) return(fail_res)
    xi <- bg_ctx$morph_X[row_i, ]; t_ab <- as.numeric(t(xi - mu1) %*% bg_ctx$morph_S_inv %*% d12) / d2_12
    rp <- (xi - mu1) - t_ab * d12; d2_perp <- max(0.0, as.numeric(t(rp) %*% bg_ctx$morph_S_inv %*% rp))
    if (!is.finite(t_ab) || t_ab < 0.22 || t_ab > 0.78 || !is.finite(d2_perp) || d2_perp > stats::qchisq(0.90, df = max(1L, ncol(bg_ctx$morph_X) - 1L))) return(fail_res)
  } else {
    if (is.null(morph_sp_i) || !(morph_sp_i %in% c(p1_tx, p2_tx)) || !is.finite(cda_prob_i) || cda_prob_i < 0.32 || cda_prob_i > 0.65 || !is.finite(cda_marg_i) || cda_marg_i > 0.22) return(fail_res)
    if (is.finite(cda_sec_p_i) && (cda_sec_p_i < 0.25 || (cda_prob_i + cda_sec_p_i) < 0.75)) return(fail_res)
    t_ab <- if (morph_sp_i == p2_tx) cda_prob_i else (1.0 - cda_prob_i); d2_perp <- 0.0
  }
  if (is.finite(cda_prob_i) && cda_prob_i >= 0.68) return(fail_res)
  if (!((is.finite(cda_marg_i) && cda_marg_i <= 0.25) || (is.finite(gmm_unc_i) && gmm_unc_i >= 0.22))) return(fail_res)
  has_p2 <- (!is.null(morph_sp_i) && morph_sp_i == p2_tx) || (!is.null(cda_sec_tx_i) && cda_sec_tx_i == p2_tx && (!is.finite(cda_sec_p_i) || cda_sec_p_i >= 0.20)) ||
            (!is.na(edaph_sp_i) && edaph_sp_i == p2_tx) || (is.finite(gmm_unc_i) && gmm_unc_i >= 0.28)
  if (!has_p2) return(fail_res)
  if (!is.null(bg_ctx$sec_organ_mat) && length(bg_ctx$sec_organ_cols) > 0L) {
    m1 <- bg_ctx$sec_organ_med_by_tx[[p1_tx]]; m2 <- bg_ctx$sec_organ_med_by_tx[[p2_tx]]; n_d <- 0L; adm_o <- FALSE
    for (sj in seq_along(bg_ctx$sec_organ_cols)) {
      vij <- bg_ctx$sec_organ_mat[row_i, sj]; mj <- bg_ctx$sec_organ_mad_pool[sj]
      if (is.finite(vij) && !is.null(m1) && !is.null(m2) && is.finite(m1[sj]) && is.finite(m2[sj]) && is.finite(mj) && abs(m2[sj] - m1[sj]) >= 1.2 * mj) {
        n_d <- n_d + 1L; to <- (vij - m1[sj]) / (m2[sj] - m1[sj])
        if ((to >= 0.15 && to <= 0.85) || (t_ab >= 0.50 && to <= 0.45) || (t_ab <= 0.50 && to >= 0.55)) adm_o <- TRUE
      }
    }
    if (n_d >= 1L && !adm_o && !(is.finite(asym_i) && asym_i >= 0.18)) return(fail_res)
  }
  list(pass = TRUE, d1_km = d1_km, d2_km = d2_km, ddoy_12 = ddoy_12, t_ab = t_ab, d2_perp = d2_perp,
       score = (1.0 - 2.0 * abs(t_ab - 0.50)) - 0.08 * d2_perp - 0.02 * (if (is.finite(d2_km)) d2_km else 10.0))
}
