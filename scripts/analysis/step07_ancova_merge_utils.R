if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists(".parse_step07_semicolon_cols", mode = "function")) {
  .parse_step07_semicolon_cols <- function(df, col_name) {
    if (!(col_name %in% names(df))) return(character(0))
    raw_vals <- as.character(stats::na.omit(df[[col_name]]))
    raw_vals <- trimws(raw_vals[nzchar(trimws(raw_vals))])
    if (length(raw_vals) == 0L) return(character(0))
    toks <- unique(trimws(unlist(strsplit(raw_vals, ";", fixed = TRUE))))
    intersect(toks[nzchar(toks) & !(tolower(toks) %in% c("na", "none", "null"))], names(df))
  }
}

if (!exists(".sort_numbered_pc_cols", mode = "function")) {
  .sort_numbered_pc_cols <- function(cols, prefix) {
    if (length(cols) <= 1L) return(cols)
    cols[order(suppressWarnings(as.integer(sub(paste0("^", prefix), "", cols))), cols)]
  }
}

.STEP07_MORPHO_NUM_COLS <- c(
  "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "plant_height_cm", "scape_to_basal_leaf_ratio",
  "capitula_count", "involucre_height_mm", "involucre_width_mm", "AR_inv", "blade_cie_l_sheet_delta",
  "blade_cie_a_browning_delta", "blade_cie_b_yellowing_delta", "corymb_head_count_est", "blade_cie_l_norm", "blade_cie_l_sd"
)

# 5. Taxon-Adjusted Institutional Batch Audit (Marginal Drop-1 ANCOVA + Partial Eta^2)
.resolve_step07_ancova_features <- function(df) {
  is_usable_num <- function(cc) {
    if (!(cc %in% names(df)) || is.logical(df[[cc]]) || is.character(df[[cc]]) || is.factor(df[[cc]])) return(FALSE)
    v <- suppressWarnings(as.numeric(df[[cc]])); v <- v[is.finite(v)]
    length(v) >= 6L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
  }
  sel_cols <- .parse_step07_semicolon_cols(df, "selected_morph_features")
  sel_cols <- sel_cols[!grepl("^blade_cie_", sel_cols)]
  sel_cols <- sel_cols[vapply(sel_cols, is_usable_num, logical(1))]
  if (length(sel_cols) > 0L) return(sel_cols)

  pc_cols <- .sort_numbered_pc_cols(grep("^PC[0-9]+$", names(df), value = TRUE), "PC")
  pc_cols <- pc_cols[!grepl("^blade_cie_", pc_cols)]
  pc_cols <- pc_cols[vapply(pc_cols, is_usable_num, logical(1))]
  if (length(pc_cols) > 0L) return(pc_cols)
  can_cols <- intersect(c("can1", "can2", "can3"), names(df))
  can_cols[vapply(can_cols, is_usable_num, logical(1))]
}

audit_institutional_batch_effects_ancova <- function(df, min_obs = 5L) {
  if (!("institutionCode" %in% names(df))) return(NULL)
  resp_cols <- .resolve_step07_ancova_features(df)
  if (length(resp_cols) == 0L) return(NULL)
  n_all <- nrow(df)
  inst_clean <- trimws(as.character(df$institutionCode))
  sp_clean <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unassigned_Packera", n_all)
  noise_m <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n_all)
  pass_m <- if ("is_passive" %in% names(df)) (!is.na(df$is_passive) & as.logical(df$is_passive)) else rep(FALSE, n_all)
  is_basal_anc <- if ("morph_tier_source" %in% names(df) && any(!is.na(df$morph_tier_source) & nzchar(as.character(df$morph_tier_source)))) {
    (!is.na(df$morph_tier_source) & as.character(df$morph_tier_source) == "Basal_Rosette")
  } else rep(TRUE, n_all)
  aw_v <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n_all); aw_v[!is.finite(aw_v) | aw_v <= 0.0] <- 1.0
  sbw_v <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, n_all); sbw_v[!is.finite(sbw_v) | sbw_v <= 0.0] <- 1.0
  w_all <- pmax(0.05, aw_v * sbw_v)

  base_raw <- !is.na(inst_clean) & nzchar(inst_clean) & !is.na(sp_clean) & nzchar(sp_clean)
  base_clean <- base_raw & !noise_m & !pass_m & is_basal_anc
  base_valid <- if (sum(base_clean) >= 6L) base_clean else if (sum(base_raw & is_basal_anc) >= 6L) (base_raw & is_basal_anc) else base_raw
  if (sum(base_valid) < 6L) return(NULL)

  rows_out <- list()
  for (rc in resp_cols) {
    rc_num <- suppressWarnings(as.numeric(df[[rc]]))
    valid_rc <- base_valid & is.finite(rc_num)
    if (sum(valid_rc) < 6L) next
    sub_df <- df[valid_rc, , drop = FALSE]; sub_df[[rc]] <- rc_num[valid_rc]; sub_df$.w_ancova <- w_all[valid_rc]
    inst_tab <- table(sub_df$institutionCode); keep_inst <- names(inst_tab[inst_tab >= min_obs])
    if (length(keep_inst) >= 2L) sub_df <- sub_df[sub_df$institutionCode %in% keep_inst, , drop = FALSE]
    if (nrow(sub_df) < 6L) next
    sub_df$institutionCode <- droplevels(factor(sub_df$institutionCode))
    sub_df$species_standardized <- droplevels(factor(sub_df$species_standardized))
    if (nlevels(sub_df$institutionCode) < 2L) next

    is_cent_sub <- if (exists(".is_step07_centroid_coord", mode = "function")) .is_step07_centroid_coord(sub_df) else rep(FALSE, nrow(sub_df))
    has_reg <- ("regional_group" %in% names(sub_df)) && {
      rg_v <- trimws(as.character(sub_df$regional_group))
      rg_ok <- !is_cent_sub & !is.na(rg_v) & nzchar(rg_v) & !(tolower(rg_v) %in% c("na", "unknown", "unassigned"))
      sum(rg_ok) == nrow(sub_df) && length(unique(rg_v[rg_ok])) >= 2L &&
        (nrow(sub_df) - nlevels(sub_df$species_standardized) - length(unique(rg_v[rg_ok])) - nlevels(sub_df$institutionCode) >= 3L)
    }
    if (has_reg) sub_df$regional_group <- droplevels(factor(trimws(as.character(sub_df$regional_group))))

    rhs_terms <- if (nlevels(sub_df$species_standardized) >= 2L && has_reg) "species_standardized + regional_group + institutionCode"
    else if (nlevels(sub_df$species_standardized) >= 2L) "species_standardized + institutionCode" else "institutionCode"
    model_lbl <- paste0("lm(Response ~ ", rhs_terms, ")")
    fit_ancova <- stats::lm(stats::as.formula(paste(rc, "~", rhs_terms)), data = sub_df, weights = sub_df$.w_ancova)
    df_inst <- NA_integer_; ss_inst <- NA_real_; ss_res <- NA_real_; f_val <- NA_real_; p_val <- NA_real_
    d1_tab <- tryCatch(stats::drop1(fit_ancova, scope = "institutionCode", test = "F"), error = function(e) NULL)
    if (!is.null(d1_tab) && "institutionCode" %in% rownames(d1_tab) && "<none>" %in% rownames(d1_tab)) {
      df_inst <- suppressWarnings(as.integer(d1_tab["institutionCode", "Df"]))
      ss_inst <- suppressWarnings(as.numeric(d1_tab["institutionCode", "Sum of Sq"]))
      ss_res <- suppressWarnings(as.numeric(d1_tab["<none>", "RSS"]))
      f_val <- suppressWarnings(as.numeric(d1_tab["institutionCode", "F value"]))
      p_val <- suppressWarnings(as.numeric(d1_tab["institutionCode", "Pr(>F)"]))
    } else if (rhs_terms == "institutionCode") {
      aov_tab <- tryCatch(stats::anova(fit_ancova), error = function(e) NULL)
      if (!is.null(aov_tab)) {
        rn <- rownames(aov_tab); inst_idx <- which(rn == "institutionCode")[1L]; res_idx <- which(rn == "Residuals")[1L]
        if (!is.na(inst_idx) && !is.na(res_idx)) {
          df_inst <- as.integer(aov_tab$Df[inst_idx]); ss_inst <- as.numeric(aov_tab$`Sum Sq`[inst_idx])
          ss_res <- as.numeric(aov_tab$`Sum Sq`[res_idx]); f_val <- as.numeric(aov_tab$`F value`[inst_idx])
          p_val <- as.numeric(aov_tab$`Pr(>F)`[inst_idx])
        }
      }
    }
    if (is.finite(df_inst) && df_inst >= 1L && is.finite(ss_inst) && is.finite(ss_res)) {
      y_rc <- sub_df[[rc]]; y_mu <- stats::weighted.mean(y_rc, w = sub_df$.w_ancova)
      ss_tot <- sum(sub_df$.w_ancova * (y_rc - y_mu)^2)
      partial_eta_sq <- if ((ss_inst + ss_res) > 0) ss_inst / (ss_inst + ss_res) else 0.0
      rows_out[[length(rows_out) + 1L]] <- data.frame(
        Response = rc, Model = model_lbl, Term = "institutionCode", Df = df_inst,
        Sum_Sq = round(ss_inst, 4), Residual_Sum_Sq = round(ss_res, 4),
        F_value = round(f_val, 4), P_value = signif(p_val, 4),
        partial_eta_sq = round(partial_eta_sq, 4),
        R_squared = round(if (is.finite(ss_tot) && ss_tot > 0) ss_inst / ss_tot else 0.0, 4),
        Batch_Effect_Status = ifelse(partial_eta_sq < 0.05, "Pass (<0.05)", "Exceeds 0.05"),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows_out) == 0L) return(NULL)
  out_df <- do.call(rbind, rows_out)
  out_df$p_adj_bh <- signif(stats::p.adjust(out_df$P_value, method = "BH"), 4)
  out_df
}

validate_upstream_contracts <- function(vouchers_df = NULL, morpho_df = NULL, multimodal_df = NULL) {
  if (!is.null(vouchers_df) && !("catalogNumber" %in% names(vouchers_df))) warning("Upstream vouchers missing canonical column: catalogNumber")
  if (!is.null(morpho_df)) {
    if (!("plant_individual_id" %in% names(morpho_df))) morpho_df$plant_individual_id <- 0L
    if (!("morphometric_outlier_flag" %in% names(morpho_df))) morpho_df$morphometric_outlier_flag <- FALSE
    if (!("misidentification_flag" %in% names(morpho_df))) morpho_df$misidentification_flag <- FALSE
  }
  TRUE
}

# 6. Merging Evidence Streams & Coalescing Non-Georeferenced Vouchers
safe_left_join_and_coalesce <- function(base_df, right_df, numeric_coalesce_cols = character(0), append_missing = TRUE) {
  if (is.null(right_df) || nrow(right_df) == 0L || !("catalogNumber" %in% names(right_df))) return(base_df)
  right_df$catalogNumber <- as.character(right_df$catalogNumber); base_df$catalogNumber <- as.character(base_df$catalogNumber)
  right_has_pid <- ("plant_individual_id" %in% names(right_df)) && any(!is.na(right_df$plant_individual_id))
  if (right_has_pid) {
    r_pid <- ifelse(is.na(right_df$plant_individual_id), 0L, suppressWarnings(as.integer(right_df$plant_individual_id))); r_pid[is.na(r_pid)] <- 0L; right_df$plant_individual_id <- r_pid
    right_sub <- right_df[!duplicated(paste0(right_df$catalogNumber, "__", right_df$plant_individual_id)), , drop = FALSE]
    multi_cats <- unique(right_sub$catalogNumber[duplicated(right_sub$catalogNumber)])
    if (length(multi_cats) > 0L) {
      if (!("plant_individual_id" %in% names(base_df))) base_df$plant_individual_id <- NA_integer_
      for (mc in multi_cats) {
        b_rows <- which(base_df$catalogNumber == mc); r_pids <- right_sub$plant_individual_id[right_sub$catalogNumber == mc]
        if (length(b_rows) == 1L && (is.na(base_df$plant_individual_id[b_rows]) || suppressWarnings(as.integer(base_df$plant_individual_id[b_rows])) == 0L)) {
          rep_blk <- base_df[rep(b_rows, length(r_pids)), , drop = FALSE]; rep_blk$plant_individual_id <- r_pids
          base_df <- rbind(base_df[seq_len(nrow(base_df)) != b_rows, , drop = FALSE], rep_blk)
        } else if (isTRUE(append_missing) && length(b_rows) >= 1L) {
          miss_pids <- setdiff(r_pids, stats::na.omit(suppressWarnings(as.integer(base_df$plant_individual_id[b_rows]))))
          if (length(miss_pids) > 0L) { rep_ex <- base_df[rep(b_rows[1L], length(miss_pids)), , drop = FALSE]; rep_ex$plant_individual_id <- miss_pids; base_df <- rbind(base_df, rep_ex) }
        }
      }
    }
  } else right_sub <- right_df[!duplicated(right_df$catalogNumber), , drop = FALSE]
  missing_cats <- setdiff(unique(right_sub$catalogNumber), unique(base_df$catalogNumber))
  if (isTRUE(append_missing) && length(missing_cats) > 0L) {
    extra_rows <- right_sub[right_sub$catalogNumber %in% missing_cats, , drop = FALSE]
    for (col in setdiff(names(base_df), names(extra_rows))) extra_rows[[col]] <- NA
    base_df <- rbind(base_df, extra_rows[, names(base_df), drop = FALSE])
  }
  join_cols <- setdiff(names(right_sub), names(base_df)); dedup_cols <- unique(c("catalogNumber", join_cols))
  if (!right_has_pid && !("plant_individual_id" %in% names(base_df)) && length(join_cols) > 0L && !any(duplicated(base_df$catalogNumber))) {
    base_df <- dplyr::left_join(base_df, right_sub[, dedup_cols, drop = FALSE], by = "catalogNumber")
  }
  idx_r <- if (right_has_pid && ("plant_individual_id" %in% names(base_df))) {
    b_pid <- ifelse(is.na(base_df$plant_individual_id), 0L, suppressWarnings(as.integer(base_df$plant_individual_id))); b_pid[is.na(b_pid)] <- 0L
    idx_ex <- match(paste0(base_df$catalogNumber, "__", b_pid), paste0(right_sub$catalogNumber, "__", right_sub$plant_individual_id))
    r_cnt <- table(right_sub$catalogNumber); s_cats <- names(r_cnt[r_cnt == 1L])
    ifelse(!is.na(idx_ex), idx_ex, ifelse(base_df$catalogNumber %in% s_cats, match(base_df$catalogNumber, right_sub$catalogNumber), NA_integer_))
  } else match(base_df$catalogNumber, right_sub$catalogNumber)

  for (col in setdiff(names(right_sub), "catalogNumber")) {
    if (col == "plant_individual_id" && ("plant_individual_id" %in% names(base_df))) {
      cur_p <- suppressWarnings(as.integer(base_df$plant_individual_id)); new_p <- suppressWarnings(as.integer(right_sub$plant_individual_id[idx_r]))
      up_m <- !is.na(new_p) & (is.na(cur_p) | (cur_p == 0L & new_p != 0L)); if (any(up_m)) base_df$plant_individual_id[up_m] <- new_p[up_m]
    } else if (col %in% numeric_coalesce_cols && col %in% names(base_df)) {
      base_df[[col]] <- dplyr::coalesce(suppressWarnings(as.numeric(base_df[[col]])), suppressWarnings(as.numeric(right_sub[[col]][idx_r])))
    } else if (col %in% names(base_df)) {
      na_m <- (is.na(base_df[[col]]) | (is.character(base_df[[col]]) & !nzchar(trimws(base_df[[col]])))) & !is.na(idx_r)
      if (any(na_m)) base_df[[col]][na_m] <- right_sub[[col]][idx_r[na_m]]
    } else base_df[[col]] <- right_sub[[col]][idx_r]
  }
  base_df
}

resolve_gmm_bic_stats <- function(opts) {
  opt_k <- 4L; delta_bic_val <- 12.0
  cands <- unique(stats::na.omit(c(opts$gmm_bic, "data/tables/gmm_bic_model_selection.csv", opts$gmm_summary, "outputs/reports/gmm_bayes_factors_summary.csv")))
  cands <- cands[nzchar(cands) & file.exists(cands)]
  if (length(cands) > 0L) {
    gdf <- tryCatch(readr::read_csv(cands[1L], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(gdf) && nrow(gdf) > 0L && "BIC" %in% names(gdf)) {
      bi <- if ("is_optimal" %in% names(gdf) && any(gdf$is_optimal == TRUE, na.rm = TRUE)) which(gdf$is_optimal == TRUE)[1L] else which.max(gdf$BIC)
      if (length(bi) > 0L && "K" %in% names(gdf)) opt_k <- as.integer(gdf$K[bi])
      for (bc in c("BayesFactor_vs_Null", "Two_Delta_BIC", "Delta_BIC")) if (bc %in% names(gdf) && length(bi) > 0L && is.finite(as.numeric(gdf[[bc]][bi]))) { delta_bic_val <- as.numeric(gdf[[bc]][bi]); break }
    }
  }
  list(opt_k = opt_k, two_delta_bic = delta_bic_val)
}

build_triage_queue <- function(opts) {
  if (is.data.frame(opts)) return(apply_taxonomic_decision_matrix(opts))
  message("Ingesting and joining evidence streams...")
  if (!file.exists(opts$vouchers)) stop("Missing vouchers: ", opts$vouchers)
  vouchers <- readr::read_csv(opts$vouchers, show_col_types = FALSE); vouchers$catalogNumber <- as.character(vouchers$catalogNumber)
  validate_upstream_contracts(vouchers_df = vouchers)
  sp_cand <- rep("Unknown", nrow(vouchers))
  for (sc in intersect(c("species_standardized", "scientificName", "species_raw"), names(vouchers))) {
    v <- standardize_packera_taxon(vouchers[[sc]]); m <- sp_cand == "Unknown" & !is.na(v) & v != "Unknown"; sp_cand[m] <- v[m]
  }
  vouchers$species_standardized <- sp_cand
  sql_master <- if (!isTRUE(opts$demo) && exists("try_fetch_master_view_for_cohort", mode = "function")) {
    tryCatch(try_fetch_master_view_for_cohort(vouchers_df = vouchers, vouchers_path = opts$vouchers, demo = isTRUE(opts$demo), require_step04 = file.exists(opts$morphometrics), require_step05 = file.exists(opts$multimodal_flags)), error = function(e) NULL)
  } else NULL
  if (!is.null(sql_master) && nrow(sql_master) > 0L) {
    vouchers <- sql_master
    if ("plant_height_mm" %in% names(vouchers) && !("plant_height_cm" %in% names(vouchers))) vouchers$plant_height_cm <- round(suppressWarnings(as.numeric(vouchers$plant_height_mm)) / 10.0, 2)
    if (!("can3" %in% names(vouchers))) vouchers$can3 <- if ("PC3" %in% names(vouchers)) suppressWarnings(as.numeric(vouchers$PC3)) else NA_real_
    message(sprintf("Hydrated %d specimen rows from SQLite v_integrative_triage_master.", nrow(vouchers)))
  } else {
  mm_numeric_coalesce <- c(
    "doy", "decimalLatitude", "decimalLongitude", "soil_ph", "soil_sand", "soil_cec", "soil_bulk_density", "delta_doy",
    "fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "ssurgo_rock_outcrop_pct",
    "edaphic_z_score", "type_edaphic_z_score", "climatic_z_score", "spatial_block_cv_acc", "rf_predicted_prob", "rf_given_taxon_prob",
    "log1p_capitulum_count", .STEP07_MORPHO_NUM_COLS
  )
  run_dir_env <- env_plantmatics_or_packera("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR")
  def_dt_dir <- normalizePath(if (nzchar(run_dir_env)) file.path(run_dir_env, "data/tables") else "data/tables", mustWork = FALSE)
  vouch_dir <- normalizePath(dirname(opts$vouchers), mustWork = FALSE); is_isolated_vouch_dir <- !identical(vouch_dir, def_dt_dir)
  if (is_isolated_vouch_dir) {
    sib_defs <- list(
      list(k = "multimodal_flags", f = "multimodal_conflict_flags.csv", fl = isTRUE(opts$cli_has_multimodal)),
      list(k = "harmonics", f = "leaf_efa_harmonics.csv", fl = isTRUE(opts$cli_has_harmonics)),
      list(k = "morphometrics", f = "morphometric_flags.csv", fl = isTRUE(opts$cli_has_morphometrics))
    )
    for (sb in sib_defs) {
      val <- opts[[sb$k]]
      if (!sb$fl && !is.null(val) && nzchar(val) && identical(normalizePath(dirname(val), mustWork = FALSE), def_dt_dir)) {
        sp_f <- file.path(dirname(opts$vouchers), sb$f); opts[[sb$k]] <- if (file.exists(sp_f)) sp_f else ""
      }
    }
  }
    if (!is.null(opts$multimodal_flags) && nzchar(opts$multimodal_flags) && file.exists(opts$multimodal_flags)) {
      mm_df <- readr::read_csv(opts$multimodal_flags, show_col_types = FALSE); mm_df$catalogNumber <- as.character(mm_df$catalogNumber)
      if ("plant_height_mm" %in% names(mm_df) && !("plant_height_cm" %in% names(mm_df))) mm_df$plant_height_cm <- round(suppressWarnings(as.numeric(mm_df$plant_height_mm)) / 10.0, 2)
      dyn_mm_num <- grep("^(bio[0-9]+_|fine_|ilr_texture_|ssurgo_|solus_|topo_)", names(mm_df), value = TRUE)
      dyn_mm_num <- dyn_mm_num[vapply(dyn_mm_num, function(cc) !is.character(mm_df[[cc]]) && !is.logical(mm_df[[cc]]), logical(1))]
      vouchers <- safe_left_join_and_coalesce(vouchers, mm_df, numeric_coalesce_cols = unique(c(mm_numeric_coalesce, dyn_mm_num)))
      for (mc_col in intersect(c("rf_predicted_taxon", "is_verified_anthesis", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors", "nearest_type_niche_taxon", "matches_type_lithology", "type_edaphic_z_score"), names(mm_df))) {
        if (!(mc_col %in% names(vouchers))) vouchers[[mc_col]] <- NA
      }
    }
    ph_cands <- file.path(dirname(opts$vouchers), "voucher_phenological_states.csv")
    if (file.exists(ph_cands)) {
      ph_df <- tryCatch(readr::read_csv(ph_cands, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(ph_df) && "catalogNumber" %in% names(ph_df)) {
        vouchers <- safe_left_join_and_coalesce(vouchers, ph_df, numeric_coalesce_cols = c("n_anthesis", "n_bud", "n_fruiting", "n_senescent", "reproductive_stage_confidence", "total_capitula"), append_missing = TRUE)
      }
    }
    med_or_na <- function(x) { v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]; if (length(v) > 0L) stats::median(v) else NA_real_ }
    first_chr_or_na <- function(x) { v <- as.character(stats::na.omit(x)); v <- v[nzchar(trimws(v))]; if (length(v) > 0L) v[1L] else NA_character_ }
    first_int_or_na <- function(x) { v <- suppressWarnings(as.integer(stats::na.omit(x))); v <- v[!is.na(v)]; if (length(v) > 0L) v[1L] else NA_integer_ }
    any_lgl_or_na <- function(x) { v <- stats::na.omit(suppressWarnings(as.logical(x))); if (length(v) > 0L) any(v) else NA }
    if (is.null(opts$morphometrics) || !nzchar(opts$morphometrics) || !file.exists(opts$morphometrics)) {
      sibling_morph <- file.path(dirname(opts$vouchers), "morphometric_flags.csv")
      if (file.exists(sibling_morph)) opts$morphometrics <- sibling_morph else if (!is_isolated_vouch_dir && file.exists("data/tables/morphometric_flags.csv")) opts$morphometrics <- "data/tables/morphometric_flags.csv"
    }
    if (!is.null(opts$morphometrics) && nzchar(opts$morphometrics) && file.exists(opts$morphometrics)) {
      mo_df <- readr::read_csv(opts$morphometrics, show_col_types = FALSE); mo_df$catalogNumber <- as.character(mo_df$catalogNumber)
      validate_upstream_contracts(morpho_df = mo_df)
      if ("cda_predicted" %in% names(mo_df) && !("cda_predicted_taxon" %in% names(mo_df))) mo_df$cda_predicted_taxon <- as.character(mo_df$cda_predicted)
      if ("plant_height_mm" %in% names(mo_df) && !("plant_height_cm" %in% names(mo_df))) mo_df$plant_height_cm <- round(suppressWarnings(as.numeric(mo_df$plant_height_mm)) / 10.0, 2)
      for (bc in c("is_promoted_bronze", "is_bent_stem", "is_typical_morphotype", "has_empirical_type_anchor", "is_type_specimen")) if (!(bc %in% names(mo_df))) mo_df[[bc]] <- NA
      for (fc in c("morphometric_outlier_flag", "misidentification_flag")) if (!(fc %in% names(mo_df))) mo_df[[fc]] <- FALSE
      if (!("plant_individual_id" %in% names(mo_df))) mo_df$plant_individual_id <- 0L
      mo_df$plant_individual_id <- ifelse(is.na(mo_df$plant_individual_id), 0L, suppressWarnings(as.integer(mo_df$plant_individual_id))); mo_df$plant_individual_id[is.na(mo_df$plant_individual_id)] <- 0L
      if (!("can3" %in% names(mo_df))) mo_df$can3 <- if ("PC3" %in% names(mo_df)) suppressWarnings(as.numeric(mo_df$PC3)) else NA_real_
      if (!("PC1" %in% names(mo_df))) mo_df$PC1 <- if ("can1" %in% names(mo_df)) suppressWarnings(as.numeric(mo_df$can1)) else NA_real_
      if (!("PC2" %in% names(mo_df))) mo_df$PC2 <- if ("can2" %in% names(mo_df)) suppressWarnings(as.numeric(mo_df$can2)) else NA_real_
      mo_num_all <- c(.STEP07_MORPHO_NUM_COLS, "cda_posterior_prob", "cda_second_posterior_prob", "cda_posterior_margin",
                      "loocv_accuracy", "can1", "can2", "gmm_uncertainty", "ecoregion_cda_posterior_prob",
                      "authority_weight", "type_morph_mahal_d2", "nearest_type_basionym_d2", "type_morph_concordance")
      for (nc in mo_num_all) if (!(nc %in% names(mo_df))) mo_df[[nc]] <- NA_real_
      if (!("cda_second_taxon" %in% names(mo_df))) mo_df$cda_second_taxon <- if ("cda_second_predicted_taxon" %in% names(mo_df)) as.character(mo_df$cda_second_predicted_taxon) else NA_character_
      for (cc in c("cda_predicted_taxon", "ecoregion_cda_predicted", "nearest_type_morph_taxon", "nearest_type_basionym", "morph_tier_source")) if (!(cc %in% names(mo_df))) mo_df[[cc]] <- NA_character_
      for (ic in c("gmm_cluster", "leaf_count", "basal_leaf_count", "cauline_leaf_count")) if (!(ic %in% names(mo_df))) mo_df[[ic]] <- NA_integer_

      mo_agg <- mo_df %>% dplyr::group_by(catalogNumber, plant_individual_id) %>% dplyr::summarise(
        morph_tier_source = first_chr_or_na(morph_tier_source), leaf_count = first_int_or_na(leaf_count),
        basal_leaf_count = first_int_or_na(basal_leaf_count), cauline_leaf_count = first_int_or_na(cauline_leaf_count),
        cda_predicted_taxon = first_chr_or_na(cda_predicted_taxon), cda_second_taxon = first_chr_or_na(cda_second_taxon),
        cda_posterior_prob = med_or_na(cda_posterior_prob), cda_second_posterior_prob = med_or_na(cda_second_posterior_prob),
        cda_posterior_margin = med_or_na(cda_posterior_margin), loocv_accuracy = med_or_na(loocv_accuracy),
        PC1 = med_or_na(PC1), PC2 = med_or_na(PC2), can1 = med_or_na(can1), can2 = med_or_na(can2), can3 = med_or_na(can3),
        gmm_cluster = first_int_or_na(gmm_cluster), gmm_uncertainty = med_or_na(gmm_uncertainty),
        is_typical_morphotype = any_lgl_or_na(is_typical_morphotype), type_morph_mahal_d2 = med_or_na(type_morph_mahal_d2),
        nearest_type_morph_taxon = first_chr_or_na(nearest_type_morph_taxon), nearest_type_basionym = first_chr_or_na(nearest_type_basionym),
        nearest_type_basionym_d2 = med_or_na(nearest_type_basionym_d2), type_morph_concordance = med_or_na(type_morph_concordance),
        has_empirical_type_anchor = any_lgl_or_na(has_empirical_type_anchor), is_type_specimen = any_lgl_or_na(is_type_specimen),
        petiole_length_mm = med_or_na(petiole_length_mm), blade_length_mm = med_or_na(blade_length_mm),
        petiole_to_blade_ratio = med_or_na(petiole_to_blade_ratio), plant_height_cm = med_or_na(plant_height_cm),
        scape_to_basal_leaf_ratio = med_or_na(scape_to_basal_leaf_ratio), is_bent_stem = any_lgl_or_na(is_bent_stem),
        capitula_count = med_or_na(capitula_count), involucre_height_mm = med_or_na(involucre_height_mm),
        involucre_width_mm = med_or_na(involucre_width_mm), AR_inv = med_or_na(AR_inv),
        blade_cie_l_sheet_delta = med_or_na(blade_cie_l_sheet_delta), blade_cie_a_browning_delta = med_or_na(blade_cie_a_browning_delta),
        blade_cie_b_yellowing_delta = med_or_na(blade_cie_b_yellowing_delta), corymb_head_count_est = med_or_na(corymb_head_count_est),
        blade_cie_l_norm = med_or_na(blade_cie_l_norm), blade_cie_l_sd = med_or_na(blade_cie_l_sd),
        ecoregion_cda_predicted = first_chr_or_na(ecoregion_cda_predicted), ecoregion_cda_posterior_prob = med_or_na(ecoregion_cda_posterior_prob),
        is_promoted_bronze = any_lgl_or_na(is_promoted_bronze), authority_weight = med_or_na(authority_weight),
        morphometric_outlier_flag = any(morphometric_outlier_flag == TRUE, na.rm = TRUE),
        misidentification_flag = any(misidentification_flag == TRUE, na.rm = TRUE), .groups = "drop"
      )
      rem_mo_cols <- setdiff(names(mo_df), c(names(mo_agg), "species_standardized", "scientificName", "species_raw", "determiner_tier"))
      dyn_mo_num_cols <- character(0)
      if (length(rem_mo_cols) > 0L) {
        mo_key <- paste0(mo_df$catalogNumber, "__", mo_df$plant_individual_id); agg_key <- paste0(mo_agg$catalogNumber, "__", mo_agg$plant_individual_id)
        lgl_known <- c("is_dissected", "is_contour_outlier", "is_confident_label_noise", "is_passive", "is_calibrated_mm", "cryptic_cluster_cohesion_pass", "cryptic_non_clinal_pass")
        int_known <- c("efd_harmonic_opt", "efd_k95_pcs", "self_morph_df", "open_set_df", "voucher_sporophytic_2n", "voucher_ploidy_level", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass")
        for (rmc in rem_mo_cols) {
          col_v <- mo_df[[rmc]]
          if (is.logical(col_v) || rmc %in% lgl_known) mo_agg[[rmc]] <- as.logical(tapply(col_v, mo_key, any_lgl_or_na)[agg_key])
          else if (is.integer(col_v) || rmc %in% int_known) mo_agg[[rmc]] <- as.integer(tapply(col_v, mo_key, first_int_or_na)[agg_key])
          else if (is.numeric(col_v)) { mo_agg[[rmc]] <- as.numeric(tapply(col_v, mo_key, med_or_na)[agg_key]); dyn_mo_num_cols <- c(dyn_mo_num_cols, rmc) }
          else mo_agg[[rmc]] <- as.character(tapply(col_v, mo_key, first_chr_or_na)[agg_key])
        }
      }
      vouchers <- safe_left_join_and_coalesce(vouchers, mo_agg, numeric_coalesce_cols = unique(c(mo_num_all, "PC1", "PC2", "can3", dyn_mo_num_cols)))
    }
    if (!is.null(opts$linear_traits) && nzchar(opts$linear_traits) && file.exists(opts$linear_traits)) {
      lt_df <- tryCatch(readr::read_csv(opts$linear_traits, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(lt_df) && "catalogNumber" %in% names(lt_df)) vouchers <- safe_left_join_and_coalesce(vouchers, lt_df, append_missing = FALSE)
    }
    if (!is.null(opts$harmonics) && nzchar(opts$harmonics) && file.exists(opts$harmonics)) {
      hm_df <- tryCatch(readr::read_csv(opts$harmonics, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(hm_df) && "catalogNumber" %in% names(hm_df)) vouchers <- safe_left_join_and_coalesce(vouchers, hm_df, append_missing = FALSE)
    }
  }

  if (any(is.na(vouchers$species_standardized) | vouchers$species_standardized == "Unknown")) {
    for (sc in intersect(c("species_standardized", "scientificName", "species_raw", "taxon"), names(vouchers))) {
      v <- standardize_packera_taxon(vouchers[[sc]]); m <- (is.na(vouchers$species_standardized) | vouchers$species_standardized == "Unknown") & !is.na(v) & v != "Unknown"
      if (any(m)) vouchers$species_standardized[m] <- v[m]
    }
    vouchers$species_standardized[is.na(vouchers$species_standardized)] <- "Unknown"
  }
  if (!("doy" %in% names(vouchers))) vouchers$doy <- NA_real_
  vouchers$doy <- suppressWarnings(as.numeric(vouchers$doy))
  if ("dayOfYear" %in% names(vouchers)) vouchers$doy <- dplyr::coalesce(vouchers$doy, suppressWarnings(as.numeric(vouchers$dayOfYear)))
  if ("eventDate" %in% names(vouchers) && any(is.na(vouchers$doy))) {
    na_d <- which(is.na(vouchers$doy) & !is.na(vouchers$eventDate))
    if (length(na_d) > 0) {
      dt_str <- substr(trimws(as.character(vouchers$eventDate[na_d])), 1, 10)
      parsed_dt <- suppressWarnings(as.Date(dt_str, format = "%Y-%m-%d"))
      slash_idx <- which(is.na(parsed_dt) & nzchar(dt_str))
      if (length(slash_idx) > 0) parsed_dt[slash_idx] <- suppressWarnings(as.Date(dt_str[slash_idx], format = "%Y/%m/%d"))
      vouchers$doy[na_d] <- suppressWarnings(as.numeric(format(parsed_dt, "%j")))
    }
  }
  if (exists("load_step07_genomic_metrics", mode = "function")) vouchers <- load_step07_genomic_metrics(opts, vouchers)
  if (exists("hydrate_step07_organellar_placement_metrics", mode = "function")) vouchers <- hydrate_step07_organellar_placement_metrics(vouchers, opts)

  opts_cyto_q <- opts
  if (is.null(opts_cyto_q$db_path)) { db_env_q <- env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"); if (nzchar(db_env_q)) opts_cyto_q$db_path <- db_env_q }
  if (is.null(opts_cyto_q$report_out) && !is.null(opts_cyto_q$output_report) && nzchar(opts_cyto_q$output_report)) opts_cyto_q$report_out <- opts_cyto_q$output_report
  cyto_dir_q <- if (!is.null(opts_cyto_q$cytology_dir) && nzchar(opts_cyto_q$cytology_dir)) opts_cyto_q$cytology_dir else {
    env_cdir_q <- env_plantmatics_or_packera("PLANTMATICS_CYTOLOGY_DIR", "PACKERA_CYTOLOGY_DIR"); if (nzchar(env_cdir_q)) env_cdir_q else "data/cytology"
  }
  if (exists("run_step07_cytology_integration", mode = "function")) {
    vouchers <- run_step07_cytology_integration(vouchers, cytology_dir = cyto_dir_q, opts = opts_cyto_q); attr(vouchers, "step07_cytology_executed") <- TRUE
  }
  cyto_md_q <- attr(vouchers, "cytology_markdown_lines"); cyto_bar_lbl_q <- attr(vouchers, "cytotype_barrier_labels"); gigas_df_q <- attr(vouchers, "gigas_concordance_df")
  if (exists("run_step07_genomic_propagation", mode = "function")) vouchers <- run_step07_genomic_propagation(vouchers, opts)
  prop_sum <- attr(vouchers, "genomic_propagation_summary")

  gmm_stats <- resolve_gmm_bic_stats(opts)
  vouchers <- apply_taxonomic_decision_matrix(vouchers, intraspecific_mode = isTRUE(opts$intraspecific_mode), gmm_two_delta_bic = gmm_stats$two_delta_bic, opts = opts)
  if (!is.null(attr(vouchers, "genomic_propagation_summary"))) prop_sum <- attr(vouchers, "genomic_propagation_summary")
  if (!is.null(prop_sum)) attr(vouchers, "genomic_propagation_summary") <- prop_sum
  if (!is.null(attr(vouchers, "cytology_markdown_lines"))) cyto_md_q <- attr(vouchers, "cytology_markdown_lines")
  if (!is.null(attr(vouchers, "cytotype_barrier_labels"))) cyto_bar_lbl_q <- attr(vouchers, "cytotype_barrier_labels")
  if (!is.null(attr(vouchers, "gigas_concordance_df"))) gigas_df_q <- attr(vouchers, "gigas_concordance_df")

  comp_res <- if (exists("run_step07_complex_discovery", mode = "function")) {
    tryCatch(run_step07_complex_discovery(vouchers, opts), error = function(e) { warning("Step 07 complex discovery warning: ", e$message); NULL })
  } else NULL
  if (!is.null(comp_res) && !is.null(comp_res$complexes_df) && nrow(comp_res$complexes_df) > 0L && !isTRUE(comp_res$subcluster_mode)) {
    cdf_m <- comp_res$complexes_df; idx_cm <- match(as.character(vouchers$species_standardized), as.character(cdf_m$species_standardized))
    vouchers$consensus_complex_id <- as.integer(cdf_m$consensus_complex_id[idx_cm]); vouchers$complex_label <- as.character(cdf_m$complex_label[idx_cm])
  }
  if (!isTRUE(opts$demo) && exists("persist_step07_to_sqlite", mode = "function")) {
    tryCatch(persist_step07_to_sqlite(vouchers, vouchers_path = opts$vouchers, demo = isTRUE(opts$demo)), error = function(e) warning("SQLite Step 07 persistence note: ", e$message))
  }
  can_eval_hyb_attr <- attr(vouchers, "can_evaluate_hybrids"); full_complex_df <- vouchers
  focal_active_flag <- !is.null(opts$focal_taxon) && nzchar(trimws(opts$focal_taxon))
  if (focal_active_flag) {
    f_std <- standardize_packera_taxon(opts$focal_taxon)
    rec_det_v <- if ("recommended_determination" %in% names(vouchers)) as.character(vouchers$recommended_determination) else rep(NA_character_, nrow(vouchers))
    hp1_v <- if ("hybrid_parent_1" %in% names(vouchers)) as.character(vouchers$hybrid_parent_1) else rep(NA_character_, nrow(vouchers))
    hp2_v <- if ("hybrid_parent_2" %in% names(vouchers)) as.character(vouchers$hybrid_parent_2) else rep(NA_character_, nrow(vouchers))
    vouchers <- vouchers[(vouchers$species_standardized == f_std) | (!is.na(rec_det_v) & rec_det_v == f_std) | (!is.na(hp1_v) & hp1_v == f_std) | (!is.na(hp2_v) & hp2_v == f_std), , drop = FALSE]
  }
  prio_order <- c("CRITICAL" = 1, "HIGH" = 2, "MEDIUM" = 3, "LOW" = 4, "RESOLVED" = 5)
  sort_triage_by_prio <- function(d_in) {
    if (!is.data.frame(d_in) || nrow(d_in) == 0L) return(d_in)
    d_in$prio_rank <- prio_order[d_in$synthesis_triage_priority]
    d_in <- if ("coordinateUncertainty" %in% names(d_in)) d_in[order(d_in$prio_rank, -d_in$coordinateUncertainty), ] else d_in[order(d_in$prio_rank), ]
    d_in$prio_rank <- NULL; d_in
  }
  vouchers <- sort_triage_by_prio(vouchers); full_complex_df <- sort_triage_by_prio(full_complex_df)
  if (!is.null(attr_keep <- list(
    complex_discovery = comp_res, can_evaluate_hybrids = can_eval_hyb_attr,
    genomic_propagation_summary = prop_sum, cytology_markdown_lines = cyto_md_q,
    cytotype_barrier_labels = cyto_bar_lbl_q, gigas_concordance_df = gigas_df_q
  ))) for (nm in names(attr_keep)) if (!is.null(attr_keep[[nm]])) attr(vouchers, nm) <- attr_keep[[nm]]
  if (!is.null(full_complex_df) && length(unique(stats::na.omit(as.character(full_complex_df$species_standardized)))) >= 2L) attr(vouchers, "full_complex_df") <- full_complex_df

  .write_triage_atomic <- function(d, p) {
    tmp_p <- paste0(p, ".tmp")
    if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(d, file = tmp_p, na = "NA", quote = "auto", logical01 = FALSE) else readr::write_csv(d, tmp_p, na = "NA")
    file.rename(tmp_p, p)
  }
  dir.create(dirname(opts$output_queue), recursive = TRUE, showWarnings = FALSE); .write_triage_atomic(vouchers, opts$output_queue)
  if (focal_active_flag && !is.null(full_complex_df) && nrow(full_complex_df) > 0L) .write_triage_atomic(full_complex_df, file.path(dirname(opts$output_queue), "triage_queue_complex_full.csv"))
  if (exists("export_genomic_field_triage_csv", mode = "function")) tryCatch(export_genomic_field_triage_csv(opts, vouchers), error = function(e) warning("Genomic field triage export note: ", e$message))

  reports_dir <- if (!is.null(opts$output_report) && nzchar(opts$output_report)) dirname(opts$output_report)
  else if (!is.null(opts$output_complexes) && nzchar(opts$output_complexes)) dirname(opts$output_complexes)
  else if (!is.null(opts$output_html) && nzchar(opts$output_html)) dirname(opts$output_html)
  else if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) { qd <- dirname(opts$output_queue); if (basename(qd) == "tables") file.path(dirname(qd), "reports") else qd } else "outputs/reports"
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    audit_cols <- c("catalogNumber", "institutionCode", "species_standardized", "identifiedBy", "determiner_tier", "recordedBy", "recordNumber", "stateProvince", "county",
                    "usda_county_status", "exsiccatae_key", "exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon", "climatic_z_score", "rf_given_taxon_prob",
                    "cda_predicted_taxon", "cda_posterior_prob", "rf_predicted_taxon", "rf_predicted_prob", "portal_lag_trigger_gate", "portal_metadata_lag_flag",
                    "is_marginal_range_extension", "recommended_determination", "synthesis_triage_action", "taxonomic_treatment_subtype", "synthesis_rationale")
    n_v <- nrow(vouchers)
    u_aud <- if ("usda_county_status" %in% names(vouchers)) tolower(trimws(ifelse(is.na(vouchers$usda_county_status), "", as.character(vouchers$usda_county_status)))) else rep("", n_v)
    ex_aud <- if ("exsiccatae_taxonomic_conflict" %in% names(vouchers)) (!is.na(suppressWarnings(as.integer(vouchers$exsiccatae_taxonomic_conflict))) & suppressWarnings(as.integer(vouchers$exsiccatae_taxonomic_conflict)) != 0L) else rep(FALSE, n_v)
    pl_aud <- if ("portal_metadata_lag_flag" %in% names(vouchers)) (!is.na(vouchers$portal_metadata_lag_flag) & as.logical(vouchers$portal_metadata_lag_flag)) else rep(FALSE, n_v)
    aud_df <- vouchers[(u_aud %in% c("out_of_range", "natureserve_exotic_or_waif")) | ex_aud | pl_aud, , drop = FALSE]
    for (ac in audit_cols) if (!(ac %in% names(aud_df))) aud_df[[ac]] <- if (nrow(aud_df) > 0L) NA else character(0)
    aud_df <- aud_df[, audit_cols, drop = FALSE]
    if (nrow(aud_df) > 0L && ("institutionCode" %in% names(aud_df))) aud_df <- aud_df[order(ifelse(is.na(aud_df$institutionCode), "", as.character(aud_df$institutionCode)), ifelse(is.na(aud_df$catalogNumber), "", as.character(aud_df$catalogNumber))), , drop = FALSE]
    .write_triage_atomic(aud_df, file.path(reports_dir, "portal_metadata_lag_audit.csv"))
  }, error = function(e) warning("Portal metadata lag audit CSV export note: ", e$message))
  if (!is.null(prop_sum) && isTRUE(prop_sum$active) && is.data.frame(prop_sum$table) && nrow(prop_sum$table) > 0L) {
    prop_out <- file.path(reports_dir, "unsequenced_genomic_propagation_summary.csv"); tmp_p <- paste0(prop_out, ".tmp")
    utils::write.csv(prop_sum$table, tmp_p, row.names = FALSE); file.rename(tmp_p, prop_out)
  }
  vouchers
}
load_and_merge_modalities <- build_triage_queue
