# ==============================================================================
# LeafMachine3 + Plantmatics: Soil Geochemistry, NWI Wetland, and Basionym Guard Helpers
# Purpose: Backfill DS-801 soil geochemistry and USFWS NWI wetland predictors from
#          SQLite environmental_cache, derive variance-guarded composite indices,
#          and compute Step 07 hydro-edaphic regimes plus homotypic synonymy guards.
# ==============================================================================

.GW_RAW_NUM_COLS <- c("soil_ca_mg_molar_ratio", "soil_ultramafic_ni_cr_index", "soil_carbonate_c_pct", "nwi_wetland_dist_m")

.is_valid_gw_db <- function(db_conn) {
  !is.null(db_conn) && requireNamespace("DBI", quietly = TRUE) &&
    isTRUE(tryCatch(DBI::dbIsValid(db_conn), error = function(e) FALSE))
}

.has_finite_variation <- function(x, min_n = 3L, min_sd = 1e-6) {
  v <- suppressWarnings(as.numeric(x))
  v <- v[is.finite(v)]
  length(v) >= min_n && is.finite(stats::sd(v)) && stats::sd(v) > min_sd
}

.pick_gw_chr <- function(df, col_nm, n) {
  if (!(col_nm %in% names(df))) return(rep(NA_character_, n))
  s <- trimws(as.character(df[[col_nm]]))
  s[!nzchar(s) | tolower(s) %in% c("na", "null", "nan")] <- NA_character_
  s
}

.pick_gw_num <- function(df, col_nm, n) {
  if (!(col_nm %in% names(df))) return(rep(NA_real_, n))
  v <- df[[col_nm]]
  suppressWarnings(as.numeric(if (is.factor(v)) as.character(v) else v))
}

.coerce_int_flag <- function(x, n) {
  if (is.null(x)) return(rep(NA_integer_, n))
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.finite(x), suppressWarnings(as.integer(round(x))), NA_integer_))
  s <- tolower(trimws(as.character(x)))
  out <- suppressWarnings(as.integer(suppressWarnings(as.numeric(s))))
  out[is.na(out) & s %in% c("true", "t", "yes", "y")] <- 1L
  out[is.na(out) & s %in% c("false", "f", "no", "n")] <- 0L
  out
}

.ensure_gw_cols <- function(df) {
  n <- nrow(df)
  for (nm in c(.GW_RAW_NUM_COLS, "ultramafic_calcareous_gradient", "wetland_seepage_index")) df[[nm]] <- .pick_gw_num(df, nm, n)
  df$is_nwi_wetland <- if ("is_nwi_wetland" %in% names(df)) .coerce_int_flag(df$is_nwi_wetland, n) else rep(NA_integer_, n)
  df$nwi_wetland_system <- .pick_gw_chr(df, "nwi_wetland_system", n)
  df
}

.apply_gw_lookup <- function(df, hit_idx, lookup_df, avail_cols) {
  ok_pos <- which(!is.na(hit_idx) & hit_idx >= 1L & hit_idx <= nrow(lookup_df))
  if (!length(ok_pos)) return(df)
  for (col_nm in intersect(avail_cols, names(lookup_df))) {
    raw_v <- lookup_df[[col_nm]][hit_idx[ok_pos]]
    if (col_nm == "nwi_wetland_system") {
      fv <- trimws(as.character(raw_v))
      fv[!nzchar(fv) | tolower(fv) %in% c("na", "null", "nan")] <- NA_character_
      need <- is.na(df[[col_nm]][ok_pos]) | !nzchar(trimws(df[[col_nm]][ok_pos]))
    } else if (col_nm == "is_nwi_wetland") {
      fv <- .coerce_int_flag(raw_v, length(ok_pos))
      need <- is.na(df[[col_nm]][ok_pos])
    } else {
      fv <- suppressWarnings(as.numeric(raw_v))
      fv[!is.finite(fv)] <- NA_real_
      need <- !is.finite(df[[col_nm]][ok_pos])
    }
    if (any(need)) df[[col_nm]][ok_pos[need]] <- fv[need]
  }
  df
}

.query_in_chunks <- function(db_conn, keys, build_sql_fn, chunk_size = 450L) {
  keys <- unique(keys[!is.na(keys) & nzchar(trimws(keys))])
  if (length(keys) == 0L) return(NULL)
  parts <- list()
  for (s_i in seq(1L, length(keys), by = chunk_size)) {
    sub_k <- keys[s_i:min(length(keys), s_i + chunk_size - 1L)]
    q_str <- paste(vapply(sub_k, function(k) as.character(DBI::dbQuoteString(db_conn, k)), character(1)), collapse = ", ")
    res <- tryCatch(DBI::dbGetQuery(db_conn, build_sql_fn(q_str)), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0L) parts[[length(parts) + 1L]] <- res
  }
  if (length(parts) == 0L) NULL else do.call(rbind, parts)
}

.backfill_gw_from_sqlite <- function(df, db_conn) {
  if (nrow(df) == 0L || !.is_valid_gw_db(db_conn)) return(df)
  if (!isTRUE(tryCatch(DBI::dbExistsTable(db_conn, "environmental_cache"), error = function(e) FALSE))) return(df)
  ec_fields <- tryCatch(DBI::dbListFields(db_conn, "environmental_cache"), error = function(e) character(0))
  avail_cols <- intersect(c(.GW_RAW_NUM_COLS, "is_nwi_wetland", "nwi_wetland_system"), ec_fields)
  if (length(avail_cols) == 0L) return(df)

  needs_any <- Reduce(`|`, lapply(avail_cols, function(cn) is.na(df[[cn]])))
  if (!any(needs_any)) return(df)

  key_vec <- if ("lat_lon_key" %in% names(df)) {
    .pick_gw_chr(df, "lat_lon_key", nrow(df))
  } else if (all(c("decimalLatitude", "decimalLongitude") %in% names(df)) && exists("format_lat_lon_key", mode = "function")) {
    format_lat_lon_key(df$decimalLatitude, df$decimalLongitude)
  } else {
    rep(NA_character_, nrow(df))
  }
  if ("lat_lon_key" %in% ec_fields && any(!is.na(key_vec) & needs_any)) {
    sel_ec <- paste(avail_cols, collapse = ", ")
    ec_df <- .query_in_chunks(db_conn, key_vec[!is.na(key_vec) & needs_any], function(q_in) {
      sprintf("SELECT lat_lon_key, %s FROM environmental_cache WHERE lat_lon_key IN (%s);", sel_ec, q_in)
    })
    if (!is.null(ec_df) && nrow(ec_df) > 0L) {
      ec_df <- ec_df[!duplicated(ec_df$lat_lon_key), , drop = FALSE]
      df <- .apply_gw_lookup(df, match(key_vec, trimws(as.character(ec_df$lat_lon_key))), ec_df, avail_cols)
    }
  }

  needs_still <- Reduce(`|`, lapply(avail_cols, function(cn) is.na(df[[cn]])))
  if (any(needs_still) && ("catalogNumber" %in% names(df)) &&
      isTRUE(tryCatch(DBI::dbExistsTable(db_conn, "vouchers"), error = function(e) FALSE))) {
    v_fields <- tryCatch(DBI::dbListFields(db_conn, "vouchers"), error = function(e) character(0))
    if (all(c("catalogNumber", "lat_lon_key") %in% v_fields) && ("lat_lon_key" %in% ec_fields)) {
      cat_vec <- .pick_gw_chr(df, "catalogNumber", nrow(df))
      sel_join <- paste(paste0("ec.", avail_cols), collapse = ", ")
      v_df <- .query_in_chunks(db_conn, cat_vec[needs_still], function(q_in) {
        sprintf("SELECT v.catalogNumber, %s FROM vouchers v INNER JOIN environmental_cache ec ON v.lat_lon_key = ec.lat_lon_key WHERE v.catalogNumber IN (%s);", sel_join, q_in)
      })
      if (!is.null(v_df) && nrow(v_df) > 0L) {
        v_df <- v_df[!duplicated(v_df$catalogNumber), , drop = FALSE]
        df <- .apply_gw_lookup(df, match(cat_vec, trimws(as.character(v_df$catalogNumber))), v_df, avail_cols)
      }
    }
  }
  df
}

enrich_geochem_wetland_predictors <- function(df, db_conn = NULL, env_dir = "data/environmental") {
  if (is.null(df) || !is.data.frame(df)) return(df)
  df <- .ensure_gw_cols(df)
  if (nrow(df) == 0L) return(df)
  df <- .backfill_gw_from_sqlite(df, db_conn)

  ca_mg <- df$soil_ca_mg_molar_ratio
  ni_cr <- df$soil_ultramafic_ni_cr_index
  carb_c <- df$soil_carbonate_c_pct
  if (.has_finite_variation(ca_mg) && (.has_finite_variation(ni_cr) || .has_finite_variation(carb_c))) {
    ok_grad <- is.finite(ca_mg) & (is.finite(ni_cr) | is.finite(carb_c)) & !is.finite(df$ultramafic_calcareous_gradient)
    if (any(ok_grad)) {
      ni_u <- ifelse(is.finite(ni_cr), pmax(ni_cr, 0.0), 0.0)
      cb_u <- ifelse(is.finite(carb_c), pmax(carb_c, 0.0), 0.0)
      grad <- round(log1p(1.0 / pmax(ca_mg, 0.05)) * (0.5 + 0.5 * pmin(ni_u / 6.0, 1.5)) - 0.5 * log1p(cb_u), 4)
      df$ultramafic_calcareous_gradient[ok_grad] <- grad[ok_grad]
    }
  }

  dist_m <- df$nwi_wetland_dist_m
  is_wet <- suppressWarnings(as.numeric(df$is_nwi_wetland))
  if (.has_finite_variation(dist_m) || .has_finite_variation(is_wet)) {
    ok_wsi <- (is.finite(dist_m) | is.finite(is_wet)) & !is.finite(df$wetland_seepage_index)
    if (any(ok_wsi)) {
      nwi_sys <- ifelse(is.na(df$nwi_wetland_system), "", df$nwi_wetland_system)
      wet_u <- ifelse(is.finite(is_wet), ifelse(is_wet > 0, 1.0, 0.0), 0.0)
      dist_u <- ifelse(is.finite(dist_m), pmax(dist_m, 0.0), ifelse(wet_u > 0, 0.0, Inf))
      decay <- ifelse(grepl("^(PEM|PFO|PSS)", nwi_sys, ignore.case = TRUE), exp(-dist_u / 60.0), 0.35 * exp(-dist_u / 60.0))
      twi_norm <- rep(NA_real_, nrow(df))
      if ("topo_twi_30m" %in% names(df)) {
        twi_v <- suppressWarnings(as.numeric(df$topo_twi_30m))
        if (.has_finite_variation(twi_v)) {
          rng <- range(twi_v[is.finite(twi_v)])
          twi_norm <- ifelse(is.finite(twi_v), pmin(1.0, pmax(0.0, (twi_v - rng[1L]) / max(rng[2L] - rng[1L], 1e-6))), NA_real_)
        }
      }
      wsi_raw <- ifelse(is.finite(twi_norm), 0.50 * wet_u + 0.35 * decay + 0.15 * twi_norm, 0.55 * wet_u + 0.45 * decay)
      df$wetland_seepage_index[ok_wsi] <- round(pmin(1.0, pmax(0.0, wsi_raw[ok_wsi])), 4)
    }
  }
  df
}

.parse_syn_json <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(character(0))
  s <- trimws(as.character(x)[1L])
  if (!nzchar(s) || tolower(s) %in% c("na", "null", "[]")) return(character(0))
  parsed <- if (requireNamespace("jsonlite", quietly = TRUE)) tryCatch(jsonlite::fromJSON(s, simplifyVector = TRUE), error = function(e) character(0)) else s
  vals <- trimws(as.character(unlist(parsed, use.names = FALSE)))
  vals[!is.na(vals) & nzchar(vals) & !(tolower(vals) %in% c("na", "null"))]
}

.norm_name_tokens <- function(vals) {
  vals <- trimws(as.character(vals))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0L) return(character(0))
  toks <- character(0)
  for (v in vals) {
    s <- gsub("\\bssp\\.", "subsp.", tolower(trimws(gsub("\\s+", " ", v))))
    m_inf <- regmatches(s, regexpr("^[a-z-]+\\s+(?:[×x]\\s+)?[a-z-]+\\s+(?:var\\.|subsp\\.|f\\.)\\s+[a-z-]+", s))
    core <- if (length(m_inf) > 0L && nzchar(m_inf[1L])) m_inf[1L] else {
      m_bin <- regmatches(s, regexpr("^[a-z-]+\\s+(?:[×x]\\s+)?[a-z-]+", s))
      if (length(m_bin) > 0L && nzchar(m_bin[1L])) m_bin[1L] else character(0)
    }
    toks <- c(toks, s, core)
  }
  unique(toks)
}

.load_backbone_profiles <- function(db_conn) {
  prof_map <- list()
  if (!.is_valid_gw_db(db_conn) || !isTRUE(tryCatch(DBI::dbExistsTable(db_conn, "taxon_backbone_profiles"), error = function(e) FALSE))) return(prof_map)
  flds <- tryCatch(DBI::dbListFields(db_conn, "taxon_backbone_profiles"), error = function(e) character(0))
  keep <- intersect(c("species_standardized", "accepted_scientific_name", "basionym_name",
                      "homotypic_synonyms_json", "heterotypic_synonyms_json", "nwpl_wetland_score", "nwpl_indicator_code"), flds)
  if (length(keep) == 0L) return(prof_map)
  tbl <- tryCatch(DBI::dbGetQuery(db_conn, sprintf("SELECT %s FROM taxon_backbone_profiles;", paste(keep, collapse = ", "))), error = function(e) NULL)
  if (is.null(tbl) || nrow(tbl) == 0L) return(prof_map)
  n_t <- nrow(tbl)
  bas_v <- .pick_gw_chr(tbl, "basionym_name", n_t); acc_v <- .pick_gw_chr(tbl, "accepted_scientific_name", n_t)
  sp_v <- .pick_gw_chr(tbl, "species_standardized", n_t); ic_v <- .pick_gw_chr(tbl, "nwpl_indicator_code", n_t)
  sc_v <- .pick_gw_num(tbl, "nwpl_wetland_score", n_t)
  for (r in seq_len(n_t)) {
    ho_r <- if ("homotypic_synonyms_json" %in% names(tbl)) .parse_syn_json(tbl$homotypic_synonyms_json[r]) else character(0)
    he_r <- if ("heterotypic_synonyms_json" %in% names(tbl)) .parse_syn_json(tbl$heterotypic_synonyms_json[r]) else character(0)
    entry <- list(homo_toks = .norm_name_tokens(c(bas_v[r], ho_r)), hetero_toks = .norm_name_tokens(he_r),
                  acc_toks = .norm_name_tokens(c(sp_v[r], acc_v[r])),
                  nwpl_score = if (is.finite(sc_v[r])) sc_v[r] else NA_real_, nwpl_code = ic_v[r])
    for (k in .norm_name_tokens(c(sp_v[r], acc_v[r]))) prof_map[[k]] <- entry
  }
  prof_map
}

evaluate_step07_hydro_edaphic_and_basionym_guards <- function(master_df, db_conn = NULL) {
  if (is.null(master_df) || !is.data.frame(master_df)) return(master_df)
  n <- nrow(master_df)
  if (!("hydro_edaphic_niche_regime" %in% names(master_df))) master_df$hydro_edaphic_niche_regime <- rep(NA_character_, n)
  if (!("wetland_fidelity_concordance" %in% names(master_df))) master_df$wetland_fidelity_concordance <- rep(NA_real_, n)
  if (!("basionym_synonymy_type" %in% names(master_df))) master_df$basionym_synonymy_type <- rep(NA_character_, n)
  if (n == 0L) return(master_df)

  master_df <- enrich_geochem_wetland_predictors(master_df, db_conn = db_conn)
  prof_map <- .load_backbone_profiles(db_conn)
  is_wet <- suppressWarnings(as.numeric(master_df$is_nwi_wetland))
  wsi <- suppressWarnings(as.numeric(master_df$wetland_seepage_index))
  ca_mg <- master_df$soil_ca_mg_molar_ratio; ni_cr <- master_df$soil_ultramafic_ni_cr_index
  carb_c <- master_df$soil_carbonate_c_pct; uc_grad <- master_df$ultramafic_calcareous_gradient
  dist_m <- master_df$nwi_wetland_dist_m; sgmc <- .pick_gw_chr(master_df, "sgmc_lithology_class", n)

  regime <- as.character(master_df$hydro_edaphic_niche_regime)
  for (i in seq_len(n)) {
    if (!(is.finite(ca_mg[i]) || is.finite(ni_cr[i]) || is.finite(carb_c[i]) || is.finite(uc_grad[i]) ||
          is.finite(is_wet[i]) || is.finite(dist_m[i]) || is.finite(wsi[i]) || !is.na(sgmc[i]))) next
    if ((is.finite(ca_mg[i]) && ca_mg[i] <= 0.7) || (is.finite(uc_grad[i]) && uc_grad[i] >= 0.9) || (is.finite(ni_cr[i]) && ni_cr[i] >= 4.0)) {
      regime[i] <- "Ultramafic_Serpentine_Barren"
    } else if ((!is.na(is_wet[i]) && is_wet[i] == 1) || (is.finite(wsi[i]) && wsi[i] >= 0.65)) {
      regime[i] <- "Palustrine_Wetland_Fen"
    } else if ((is.finite(ca_mg[i]) && ca_mg[i] >= 2.0) || (is.finite(carb_c[i]) && carb_c[i] >= 0.5) || (is.finite(uc_grad[i]) && uc_grad[i] <= -0.35)) {
      regime[i] <- "Calcareous_Carbonate_Glade"
    } else if ((is.finite(ca_mg[i]) && ca_mg[i] < 1.2 && (!is.finite(ni_cr[i]) || ni_cr[i] < 4.0)) ||
               (!is.na(sgmc[i]) && grepl("felsic|granit|quartzit", sgmc[i], ignore.case = TRUE) &&
                (!is.finite(ca_mg[i]) || (ca_mg[i] < 2.0 && (!is.finite(ni_cr[i]) || ni_cr[i] < 4.0) && (!is.finite(carb_c[i]) || carb_c[i] < 0.5))))) {
      regime[i] <- "Siliceous_Felsic_Outcrop"
    } else {
      regime[i] <- "Mesic_Upland_Matrix"
    }
  }
  master_df$hydro_edaphic_niche_regime <- regime

  sp_vec <- .pick_gw_chr(master_df, "species_standardized", n); acc_vec <- .pick_gw_chr(master_df, "accepted_scientific_name", n)
  nwpl_sc <- .pick_gw_num(master_df, "nwpl_wetland_score", n); nwpl_ic <- .pick_gw_chr(master_df, "nwpl_indicator_code", n)
  has_res_col <- "resurrection_basionym" %in% names(master_df)
  res_vec <- .pick_gw_chr(master_df, "resurrection_basionym", n); ntb_vec <- .pick_gw_chr(master_df, "nearest_type_basionym", n)
  row_bas <- if ("backbone_basionym_name" %in% names(master_df)) .pick_gw_chr(master_df, "backbone_basionym_name", n) else .pick_gw_chr(master_df, "basionym_name", n)
  row_ho <- if ("homotypic_synonyms_json" %in% names(master_df)) as.character(master_df$homotypic_synonyms_json) else rep(NA_character_, n)
  row_he <- if ("heterotypic_synonyms_json" %in% names(master_df)) as.character(master_df$heterotypic_synonyms_json) else rep(NA_character_, n)
  is_tp <- if ("is_type_specimen" %in% names(master_df)) (.coerce_int_flag(master_df$is_type_specimen, n) == 1L) else rep(FALSE, n)
  is_tp[is.na(is_tp)] <- FALSE

  syn_type <- as.character(master_df$basionym_synonymy_type)
  for (i in seq_len(n)) {
    prof <- NULL
    for (k in .norm_name_tokens(c(sp_vec[i], acc_vec[i]))) if (k %in% names(prof_map)) { prof <- prof_map[[k]]; break }
    if (!is.null(prof)) {
      if (!is.finite(nwpl_sc[i]) && is.finite(prof$nwpl_score)) nwpl_sc[i] <- prof$nwpl_score
      if (is.na(nwpl_ic[i]) && !is.na(prof$nwpl_code)) nwpl_ic[i] <- prof$nwpl_code
    }
    ho_toks <- unique(c(if (!is.null(prof)) prof$homo_toks else character(0), .norm_name_tokens(c(row_bas[i], .parse_syn_json(row_ho[i])))))
    he_toks <- unique(c(if (!is.null(prof)) prof$hetero_toks else character(0), .norm_name_tokens(.parse_syn_json(row_he[i]))))
    ac_toks <- unique(c(if (!is.null(prof)) prof$acc_toks else character(0), .norm_name_tokens(c(sp_vec[i], acc_vec[i]))))
    has_prof <- !is.null(prof) || length(ho_toks) > 0L || length(he_toks) > 0L

    if (!is.na(res_vec[i])) {
      r_toks <- .norm_name_tokens(res_vec[i])
      if (any(r_toks %in% ho_toks)) {
        syn_type[i] <- "Homotypic"; res_vec[i] <- NA_character_
      } else if (any(r_toks %in% he_toks)) {
        syn_type[i] <- "Heterotypic"
      }
    } else if (!is.na(ntb_vec[i])) {
      n_toks <- .norm_name_tokens(ntb_vec[i])
      if (any(n_toks %in% he_toks)) {
        syn_type[i] <- "Heterotypic"
      } else if (any(n_toks %in% ho_toks)) {
        syn_type[i] <- "Homotypic"
      } else if (is_tp[i] && any(n_toks %in% ac_toks)) {
        syn_type[i] <- "Accepted_Type"
      }
    } else if (is_tp[i] && has_prof) {
      syn_type[i] <- "Accepted_Type"
    }
  }

  ok_wfc <- is.finite(wsi) & is.finite(nwpl_sc)
  if (any(ok_wfc)) master_df$wetland_fidelity_concordance[ok_wfc] <- round(pmin(1.0, pmax(0.0, 1.0 - abs(wsi[ok_wfc] - nwpl_sc[ok_wfc]))), 4)
  master_df$basionym_synonymy_type <- syn_type
  if (has_res_col) master_df$resurrection_basionym <- res_vec
  if (any(is.finite(nwpl_sc)) || ("nwpl_wetland_score" %in% names(master_df))) master_df$nwpl_wetland_score <- nwpl_sc
  if (any(!is.na(nwpl_ic)) || ("nwpl_indicator_code" %in% names(master_df))) master_df$nwpl_indicator_code <- nwpl_ic
  master_df
}
