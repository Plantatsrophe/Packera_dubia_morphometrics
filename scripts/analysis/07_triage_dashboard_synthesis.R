#!/usr/bin/env Rscript
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Script: 07_triage_dashboard_synthesis.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
# ==============================================================================

suppressPackageStartupMessages({
  for (.pkg_s7 in c("dplyr", "readr", "tibble", "ggplot2", "gridExtra", "optparse")) {
    if (requireNamespace(.pkg_s7, quietly = TRUE)) library(.pkg_s7, character.only = TRUE)
  }
})

.script_dir_step07 <- {
  f_cmd <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f_frm <- unlist(lapply(sys.frames(), function(fr) fr$ofile))
  cand_d <- dirname(c(f_cmd, f_frm))
  cand_d <- cand_d[!is.na(cand_d) & nzchar(cand_d)]
  if (length(cand_d) > 0L) cand_d[1L] else "scripts/analysis"
}

.step07_load_env <- environment()
.source_step07_dep <- function(repo_rel, rel_to_script, symbol_chk) {
  if (exists(symbol_chk, envir = .step07_load_env, mode = "function")) return(invisible(TRUE))
  for (cand_p in unique(c(repo_rel, file.path(.script_dir_step07, rel_to_script)))) {
    if (!is.na(cand_p) && nzchar(cand_p) && file.exists(cand_p)) {
      source(cand_p, local = .step07_load_env)
      if (exists(symbol_chk, envir = .step07_load_env, mode = "function")) return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

.source_step07_dep("scripts/morphometrics/authority_weighting_utils.R", "../morphometrics/authority_weighting_utils.R", ".identify_type_specimen_vec")
.source_step07_dep("scripts/morphometrics/pub_figures_utils.R", "../morphometrics/pub_figures_utils.R", "generate_step07_pub_figures")
.source_step07_dep("scripts/core/db_utils.R", "../core/db_utils.R", "persist_step07_to_sqlite")
.source_step07_dep("scripts/analysis/taxon_complex_discovery_utils.R", "taxon_complex_discovery_utils.R", "run_step07_complex_discovery")
.source_step07_dep("scripts/analysis/genomic_triage_utils.R", "genomic_triage_utils.R", "load_step07_genomic_metrics")
.source_step07_dep("scripts/analysis/step07_genomic_propagation.R", "step07_genomic_propagation.R", "run_step07_genomic_propagation")
.source_step07_dep("scripts/analysis/step07_cytology_integration.R", "step07_cytology_integration.R", "run_step07_cytology_integration")
.source_step07_dep("scripts/analysis/edaphic_geochem_wetland_utils.R", "edaphic_geochem_wetland_utils.R", "evaluate_step07_hydro_edaphic_and_basionym_guards")
.source_step07_dep("scripts/analysis/step07_organellar_triage_utils.R", "step07_organellar_triage_utils.R", "evaluate_step07_organellar_and_oop_signals")
.source_step07_dep("scripts/analysis/step07_hopkins_consensus_utils.R", "step07_hopkins_consensus_utils.R", "compute_probabilistic_consensus_scores")
.source_step07_dep("scripts/analysis/step07_biogeo_hybrid_utils.R", "step07_biogeo_hybrid_utils.R", ".eval_five_gate_hybrid_pair_impl")
.source_step07_dep("scripts/analysis/step07_cluster_cohesion_utils.R", "step07_cluster_cohesion_utils.R", ".build_step07_cluster_cohesion_context")
.source_step07_dep("scripts/analysis/step07_decision_stages_utils.R", "step07_decision_stages_utils.R", ".evaluate_step07_specimen_stages")
.source_step07_dep("scripts/analysis/step07_postprocess_guards_utils.R", "step07_postprocess_guards_utils.R", ".apply_step07_postprocess_and_guards")
.source_step07_dep("scripts/analysis/step07_ancova_merge_utils.R", "step07_ancova_merge_utils.R", "build_triage_queue")
.source_step07_dep("scripts/analysis/step07_reports_html_utils.R", "step07_reports_html_utils.R", "generate_interactive_triage_html") # delegates render_step07_synthesis_plate_png

parse_args_robust <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  run_dir <- env_plantmatics_or_packera("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR")
  rp <- function(rel) if (nzchar(run_dir)) file.path(run_dir, rel) else rel
  option_list <- list(
    optparse::make_option(c("-v", "--vouchers"), type = "character", default = rp("data/tables/curated_vouchers.csv")),
    optparse::make_option(c("--phenology"), type = "character", default = ""),
    optparse::make_option(c("-m", "--morphometrics"), type = "character", default = rp("data/tables/morphometric_flags.csv")),
    optparse::make_option(c("--linear-traits"), type = "character", default = ""),
    optparse::make_option(c("-c", "--multimodal-flags"), type = "character", default = rp("data/tables/multimodal_conflict_flags.csv")),
    optparse::make_option(c("--spatial-flags"), type = "character", default = ""),
    optparse::make_option(c("-g", "--gmm-summary"), type = "character", default = rp("outputs/reports/gmm_bayes_factors_summary.csv")),
    optparse::make_option(c("--gmm-bic"), type = "character", default = rp("data/tables/gmm_bic_model_selection.csv")),
    optparse::make_option(c("--cda-summary"), type = "character", default = rp("data/tables/cda_summary.csv")),
    optparse::make_option(c("-s", "--niche-summary"), type = "character", default = rp("outputs/reports/multimodal_spatial_rf_summary.csv")),
    optparse::make_option(c("-e", "--harmonics"), type = "character", default = rp("data/tables/leaf_efa_harmonics.csv")),
    optparse::make_option(c("-q", "--output-queue"), type = "character", default = rp("data/tables/triage_queue.csv")),
    optparse::make_option(c("-p", "--output-plot"), type = "character", default = rp("outputs/figures/fig07a_integrative_taxonomic_synthesis.png")),
    optparse::make_option(c("-r", "--output-report"), type = "character", default = rp("outputs/reports/Packera_Taxonomic_Revision_Summary.md")),
    optparse::make_option(c("--output-html"), type = "character", default = rp("outputs/reports/Triage_Dashboard_Interactive.html")),
    optparse::make_option(c("--min-affinity"), type = "double", default = 0.03),
    optparse::make_option(c("--output-complexes"), type = "character", default = rp("outputs/reports/empirical_taxon_complexes.csv")),
    optparse::make_option(c("--output-affinity"), type = "character", default = rp("outputs/reports/pairwise_taxon_affinity_matrix.csv")),
    optparse::make_option(c("--output-neighbors"), type = "character", default = rp("outputs/reports/focal_auto_neighbors.txt")),
    optparse::make_option(c("--output-complex-plot"), type = "character", default = rp("outputs/figures/fig07b_taxon_complex_delimitation_network.png")),
    optparse::make_option(c("--genomic-dir"), type = "character", default = ""),
    optparse::make_option(c("--focal-taxon"), type = "character", default = ""),
    optparse::make_option(c("--intraspecific-mode"), action = "store_true", default = FALSE),
    optparse::make_option(c("--demo"), action = "store_true", default = FALSE)
  )
  if (requireNamespace("optparse", quietly = TRUE)) {
    raw_cli <- cli_args
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]", option_list = option_list), args = cli_args)
    names(opts) <- gsub("-", "_", names(opts))
    if (!is.null(opts$spatial_flags) && nzchar(opts$spatial_flags)) opts$multimodal_flags <- opts$spatial_flags
    opts$cli_has_multimodal <- any(grepl("^(-c|--multimodal-flags|--spatial-flags)(=|$)", raw_cli))
    opts$cli_has_harmonics <- any(grepl("^(-e|--harmonics)(=|$)", raw_cli))
    opts$cli_has_morphometrics <- any(grepl("^(-m|--morphometrics)(=|$)", raw_cli))
    return(opts)
  }
  list(vouchers = rp("data/tables/curated_vouchers.csv"), harmonics = rp("data/tables/leaf_efa_harmonics.csv"),
       morphometrics = rp("data/tables/morphometric_flags.csv"), multimodal_flags = rp("data/tables/multimodal_conflict_flags.csv"),
       gmm_summary = rp("outputs/reports/gmm_bayes_factors_summary.csv"), gmm_bic = rp("data/tables/gmm_bic_model_selection.csv"),
       cda_summary = rp("data/tables/cda_summary.csv"), niche_summary = rp("outputs/reports/multimodal_spatial_rf_summary.csv"),
       output_queue = rp("data/tables/triage_queue.csv"), output_plot = rp("outputs/figures/fig07a_integrative_taxonomic_synthesis.png"),
       output_report = rp("outputs/reports/Packera_Taxonomic_Revision_Summary.md"),
       output_html = rp("outputs/reports/Triage_Dashboard_Interactive.html"),
       min_affinity = 0.03, output_complexes = rp("outputs/reports/empirical_taxon_complexes.csv"),
       output_affinity = rp("outputs/reports/pairwise_taxon_affinity_matrix.csv"),
       output_neighbors = rp("outputs/reports/focal_auto_neighbors.txt"),
       output_complex_plot = rp("outputs/figures/fig07b_taxon_complex_delimitation_network.png"),
       genomic_dir = "", focal_taxon = "", intraspecific_mode = FALSE, demo = FALSE)
}

apply_taxonomic_decision_matrix <- function(df, intraspecific_mode = FALSE, gmm_two_delta_bic = 12.0, opts = list()) {
  message("Evaluating vouchers across Five-Stage Hierarchical Integrative Decision Matrix...")
  if (!("morphometric_outlier_flag" %in% names(df))) df$morphometric_outlier_flag <- FALSE
  if (!("misidentification_flag" %in% names(df))) df$misidentification_flag <- FALSE
  is_intra <- isTRUE(intraspecific_mode); n <- nrow(df); if (n == 0L) return(df)

  opts_cyto <- opts
  if (is.null(opts_cyto$db_path)) { dbe <- env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"); if (nzchar(dbe)) opts_cyto$db_path <- dbe }
  if (is.null(opts_cyto$report_out) && !is.null(opts_cyto$output_report) && nzchar(opts_cyto$output_report)) opts_cyto$report_out <- opts_cyto$output_report
  cyto_dir_arg <- if (!is.null(opts_cyto$cytology_dir) && nzchar(opts_cyto$cytology_dir)) opts_cyto$cytology_dir else {
    ce <- env_plantmatics_or_packera("PLANTMATICS_CYTOLOGY_DIR", "PACKERA_CYTOLOGY_DIR"); if (nzchar(ce)) ce else "data/cytology"
  }
  has_prepop_cyto <- is.null(opts$cytological_records_df) && is.null(opts$taxon_cytology_summary_df) && (
    isTRUE(attr(df, "step07_cytology_executed")) ||
    (("cytotype_reproductive_barrier_pass" %in% names(df)) && any(!is.na(df$cytotype_reproductive_barrier_pass))) ||
    (("cytology_linkage_tier" %in% names(df)) && any(!is.na(df$cytology_linkage_tier) & nzchar(trimws(as.character(df$cytology_linkage_tier)))))
  )
  if (!has_prepop_cyto && exists("run_step07_cytology_integration", mode = "function")) {
    df <- run_step07_cytology_integration(df, cytology_dir = cyto_dir_arg, opts = opts_cyto); attr(df, "step07_cytology_executed") <- TRUE
  }
  cyto_md_attr <- attr(df, "cytology_markdown_lines"); cyto_bar_lbl_attr <- attr(df, "cytotype_barrier_labels"); gigas_df_attr <- attr(df, "gigas_concordance_df")

  if (!("genomic_propagation_tier" %in% names(df)) && any(.has_step07_geno_row_mask(df)) && exists("run_step07_genomic_propagation", mode = "function")) {
    df <- run_step07_genomic_propagation(df, opts)
  }
  prop_sum_attr <- attr(df, "genomic_propagation_summary")
  sp_all_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unknown", n)
  sp_sheet_vec <- ifelse(!is.na(sp_all_vec) & nzchar(sp_all_vec), sp_all_vec, "Unknown")
  cat_sheet_vec <- if ("catalogNumber" %in% names(df)) as.character(df$catalogNumber) else paste0("__ROW_", seq_len(n))
  cat_sheet_vec[is.na(cat_sheet_vec) | !nzchar(trimws(cat_sheet_vec))] <- paste0("__ROW_", which(is.na(cat_sheet_vec) | !nzchar(trimws(cat_sheet_vec))))
  sheet_idx_list <- split(seq_len(n), paste(cat_sheet_vec, sp_sheet_vec, sep = "||"))

  cda_pred_init <- .chr_col(df, "cda_predicted_taxon", n, NA_character_); cda_prob_init <- .num_col(df, "cda_posterior_prob", n); tm_d2_init <- .num_col(df, "type_morph_mahal_d2", n)
  bas_raw_i <- if ("basal_leaf_count" %in% names(df)) suppressWarnings(as.integer(df$basal_leaf_count)) else rep(NA_integer_, n)
  caul_raw_i <- if ("cauline_leaf_count" %in% names(df)) suppressWarnings(as.integer(df$cauline_leaf_count)) else rep(NA_integer_, n)
  bc_sum <- ifelse(is.na(bas_raw_i) & is.na(caul_raw_i), NA_integer_, ifelse(is.na(bas_raw_i), 0L, bas_raw_i) + ifelse(is.na(caul_raw_i), 0L, caul_raw_i))
  n_lvs_i <- dplyr::coalesce(if ("n_leaves" %in% names(df)) suppressWarnings(as.integer(df$n_leaves)) else rep(NA_integer_, n),
                             if ("leaf_count" %in% names(df)) suppressWarnings(as.integer(df$leaf_count)) else rep(NA_integer_, n), bc_sum, bas_raw_i)
  orig_has_obs_morph_v <- !is.na(cda_pred_init) & nzchar(trimws(cda_pred_init)) & is.finite(cda_prob_init)
  orig_n_leaves_v <- ifelse(!is.na(n_lvs_i), pmax(0L, n_lvs_i), ifelse(orig_has_obs_morph_v, 1L, 0L))
  has_observed_or_sheet_morph <- orig_has_obs_morph_v
  prop_morph_cols <- intersect(c("cda_predicted_taxon", "cda_posterior_prob", "cda_oof_posterior_prob", "cl_target_threshold", "self_morph_mahal_d2", "self_morph_df",
    "cda_posterior_margin", "cda_runner_up_taxon", "cda_second_taxon", "cda_second_predicted_taxon", "cda_second_posterior_prob", "cda_mahal_d2", "type_morph_mahal_d2",
    "nearest_type_morph_taxon", "nearest_type_basionym", "nearest_type_basionym_d2", "type_morph_concordance", "is_typical_morphotype", "min_congener_mahal_d2",
    "min_congener_taxon", "open_set_novelty_p", "open_set_df", "gmm_cluster", "gmm_uncertainty", "is_morph_outlier", "morphometric_outlier_flag", "misidentification_flag"), names(df))
  mts_init_v <- .chr_col(df, "morph_tier_source", n, "Basal_Rosette"); is_caul_init <- !is.na(mts_init_v) & (mts_init_v == "Cauline_Only_Fallback")
  is_basal_obs <- orig_has_obs_morph_v & !is_caul_init
  for (s_idx in sheet_idx_list) {
    if (length(s_idx) < 2L) next
    hb <- any(is_basal_obs[s_idx]); obs_s <- if (hb) s_idx[is_basal_obs[s_idx]] else s_idx[orig_has_obs_morph_v[s_idx]]
    miss_s <- if (hb) s_idx[!orig_has_obs_morph_v[s_idx] | is_caul_init[s_idx]] else s_idx[!orig_has_obs_morph_v[s_idx]]
    if (length(obs_s) >= 1L && length(miss_s) >= 1L) {
      di <- obs_s[order(-orig_n_leaves_v[obs_s], -cda_prob_init[obs_s])[1L]]
      for (mc in prop_morph_cols) df[[mc]][miss_s] <- df[[mc]][di]
      if ("cda_posterior_prob" %in% names(df) && is.finite(df$cda_posterior_prob[di])) df$cda_posterior_prob[miss_s] <- 0.90 * as.numeric(df$cda_posterior_prob[di])
      if ("cda_posterior_margin" %in% names(df) && is.finite(df$cda_posterior_margin[di])) df$cda_posterior_margin[miss_s] <- 0.90 * as.numeric(df$cda_posterior_margin[di])
      has_observed_or_sheet_morph[miss_s] <- TRUE
    }
  }
  df$has_observed_or_sheet_morph <- has_observed_or_sheet_morph
  df <- compute_hopkins_adjusted_phenology(df); df <- compute_probabilistic_consensus_scores(df)

  active_taxa <- sort(unique(sp_all_vec[!is.na(sp_all_vec) & nzchar(sp_all_vec) & sp_all_vec != "Unknown"]))
  anchor_all <- .resolve_step07_anchor_mask(df); mo_all <- .lgl_col(df, "morphometric_outlier_flag", n); noise_all <- .lgl_col(df, "is_confident_label_noise", n)
  doy_num_all <- .num_col(df, "doy", n); valid_raw_doy <- is.finite(doy_num_all) & doy_num_all >= 1.0 & doy_num_all <= 366.0
  hop_in_all <- .num_col(df, "hopkins_adjusted_doy", n); valid_hop <- is.finite(hop_in_all) & hop_in_all >= 1.0 & hop_in_all <= 366.0
  valid_doy_all <- valid_raw_doy | valid_hop; hop_doy_all <- ifelse(valid_hop, hop_in_all, ifelse(valid_raw_doy, doy_num_all, NA_real_))
  anth_all <- if ("is_active_anthesis" %in% names(df)) { v <- as.logical(df$is_active_anthesis); v[is.na(v)] <- TRUE; v } else rep(TRUE, n)
  ph_num_all <- if ("fine_soil_ph" %in% names(df)) .num_col(df, "fine_soil_ph", n) else rep(NA_real_, n)
  if ("soil_ph" %in% names(df)) ph_num_all <- ifelse(is.finite(ph_num_all), ph_num_all, .num_col(df, "soil_ph", n))
  sand_num_all <- if ("fine_soil_sand" %in% names(df)) .num_col(df, "fine_soil_sand", n) else rep(NA_real_, n)
  if ("soil_sand" %in% names(df)) sand_num_all <- ifelse(is.finite(sand_num_all), sand_num_all, .num_col(df, "soil_sand", n))

  .circ_hop_stats <- function(h_vals, w_vals = NULL) {
    ok_h <- is.finite(h_vals); hv <- h_vals[ok_h]; if (length(hv) == 0L) return(list(mu_theta = 2.0 * pi * (120.0 / 365.25), q75_arc = 15.0, q75_doy = 135.0, med_doy = 120.0))
    th <- 2.0 * pi * (hv / 365.25)
    mu_th <- if (!is.null(w_vals) && length(w_vals) == length(h_vals)) { wv <- w_vals[ok_h]; wv[!is.finite(wv) | wv <= 0.0] <- 1.0; atan2(sum(wv * sin(th)) / sum(wv), sum(wv * cos(th)) / sum(wv)) } else atan2(mean(sin(th)), mean(cos(th)))
    if (mu_th < 0.0) mu_th <- mu_th + 2.0 * pi
    q75_a <- as.numeric(stats::quantile(atan2(sin(th - mu_th), cos(th - mu_th)) * (365.25 / (2.0 * pi)), probs = 0.75, na.rm = TRUE, names = FALSE))
    med_d <- ((mu_th * (365.25 / (2.0 * pi)) - 1.0) %% 365.25) + 1.0
    list(mu_theta = mu_th, q75_arc = q75_a, q75_doy = ((med_d + q75_a - 1.0) %% 365.25) + 1.0, med_doy = med_d)
  }
  glob_hop_pool <- hop_doy_all[valid_doy_all & anth_all & anchor_all & !mo_all & !noise_all]
  if (length(glob_hop_pool) == 0L) glob_hop_pool <- hop_doy_all[valid_doy_all & anth_all & !noise_all]
  if (length(glob_hop_pool) == 0L) glob_hop_pool <- hop_doy_all[valid_doy_all]
  glob_circ <- .circ_hop_stats(glob_hop_pool)
  mu_theta_hop <- stats::setNames(rep(glob_circ$mu_theta, length(active_taxa)), active_taxa); q75_arc_hop <- stats::setNames(rep(glob_circ$q75_arc, length(active_taxa)), active_taxa)
  q75_doy <- stats::setNames(rep(glob_circ$q75_doy, length(active_taxa)), active_taxa); med_doy <- stats::setNames(rep(glob_circ$med_doy, length(active_taxa)), active_taxa)
  active_geno_ref <- any(.has_step07_geno_row_mask(df))
  t0_info <- if (active_geno_ref) { if (exists("get_tier0_rectified_species_and_evictions", mode = "function")) get_tier0_rectified_species_and_evictions(df) else if (exists("resolve_tier0_anchor_overrides", mode = "function")) resolve_tier0_anchor_overrides(df) else NULL } else NULL
  ext_non_type_all <- .is_step07_extralimital_non_type(df); ref_idx_by_taxon <- list()
  for (tx in active_taxa) {
    idx_t <- which(sp_all_vec == tx)
    if (!is.null(t0_info)) {
      pe <- t0_info$evict_from_given_anchor & (!t0_info$is_tier0_dna_anchor | t0_info$rectified_taxon != tx)
      t0r <- which(t0_info$is_tier0_dna_anchor & !is.na(t0_info$rectified_taxon) & (t0_info$rectified_taxon == tx) & !noise_all & !ext_non_type_all)
      ig <- unique(c(idx_t[anchor_all[idx_t] & !mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t] & !pe[idx_t]], t0r))
      ic <- unique(c(idx_t[!mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t] & !pe[idx_t]], t0r)); inn <- unique(c(idx_t[!noise_all[idx_t] & !pe[idx_t]], t0r))
      ri <- if (length(ig) >= 2L) ig else if (length(ic) >= 2L) ic else if (length(inn) >= 1L) inn else idx_t; ref_idx_by_taxon[[tx]] <- ri
      ra <- ri[valid_doy_all[ri] & anth_all[ri]]; if (length(ra) == 0L) ra <- idx_t[valid_doy_all[idx_t] & anth_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t] & !pe[idx_t]]
      if (length(ra) == 0L) ra <- idx_t[valid_doy_all[idx_t] & !pe[idx_t]]; if (length(ra) == 0L) ra <- idx_t[valid_doy_all[idx_t]]
      if (length(ra) > 0L) {
        wa <- ifelse(t0_info$is_tier0_dna_anchor[ra] & (t0_info$rectified_taxon[ra] == tx), ifelse(is.finite(t0_info$tier0_weight[ra]), t0_info$tier0_weight[ra], 4.0), 1.0)
        tc <- .circ_hop_stats(hop_doy_all[ra], w_vals = wa); mu_theta_hop[[tx]] <- tc$mu_theta; q75_arc_hop[[tx]] <- tc$q75_arc; q75_doy[[tx]] <- tc$q75_doy; med_doy[[tx]] <- tc$med_doy
      }
    } else {
      ig <- idx_t[anchor_all[idx_t] & !mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t]]; ic <- idx_t[!mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t]]; inn <- idx_t[!noise_all[idx_t]]
      ri <- if (length(ig) >= 2L) ig else if (length(ic) >= 2L) ic else if (length(inn) >= 1L) inn else idx_t; ref_idx_by_taxon[[tx]] <- ri
      dr <- hop_doy_all[ri[valid_doy_all[ri] & anth_all[ri]]]; if (length(dr) == 0L) dr <- hop_doy_all[idx_t[valid_doy_all[idx_t] & anth_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t]]]; if (length(dr) == 0L) dr <- hop_doy_all[idx_t[valid_doy_all[idx_t]]]
      if (length(dr) > 0L) { tc <- .circ_hop_stats(dr); mu_theta_hop[[tx]] <- tc$mu_theta; q75_arc_hop[[tx]] <- tc$q75_arc; q75_doy[[tx]] <- tc$q75_doy; med_doy[[tx]] <- tc$med_doy }
    }
  }
  edaph_z_vec <- .compute_dynamic_edaphic_z_scores(df, sp_all_vec, active_taxa, ref_idx_by_taxon, clean_coh_mask = (anchor_all & !mo_all & !noise_all))
  pheno_taus <- .compute_taxon_pheno_thresholds(df); bg_ctx <- .build_step07_biogeo_and_morph_context(df, opts, sp_all_vec, cat_sheet_vec, ref_idx_by_taxon, active_taxa, is_intra, active_geno_ref, t0_info, valid_doy_all, anth_all, hop_doy_all, anchor_all, mo_all, noise_all)
  can_evaluate_hybrids <- bg_ctx$can_evaluate_hybrids

  .eval_five_gate_hybrid_pair <- function(row_i, p1_tx, p2_tx, morph_sp_i, cda_prob_i, cda_sec_tx_i, cda_sec_p_i, cda_marg_i, gmm_unc_i, edaph_sp_i, d_edaph_i, edaph_z_i, asym_i) {
    .eval_five_gate_hybrid_pair_impl(row_i, p1_tx, p2_tx, morph_sp_i, cda_prob_i, cda_sec_tx_i, cda_sec_p_i, cda_marg_i, gmm_unc_i, edaph_sp_i, d_edaph_i, edaph_z_i, asym_i, bg_ctx, cat_sheet_vec, mu_theta_hop, valid_doy_all, anth_all, hop_doy_all)
  }

  has_geno_dv <- .has_step07_geno_row_mask(df); geno_pred_dv <- .chr_col(df, "geno_predicted_taxon", n, NA_character_)
  geno_gdi_dv <- .num_col(df, "coalescent_gdi", n); gdi_att <- attr(df, "genomic_gdi_delimitation")
  if (!is.null(gdi_att) && is.data.frame(gdi_att) && nrow(gdi_att) > 0L && all(c("taxon", "coalescent_gdi") %in% names(gdi_att))) {
    mgd <- has_geno_dv & !is.finite(geno_gdi_dv); if (any(mgd)) {
      gtk <- ifelse(!is.na(geno_pred_dv) & nzchar(geno_pred_dv), geno_pred_dv, sp_all_vec); m_g <- match(gtk[mgd], as.character(gdi_att$taxon))
      ok_g <- !is.na(m_g); if (any(ok_g)) geno_gdi_dv[which(mgd)[ok_g]] <- suppressWarnings(as.numeric(gdi_att$coalescent_gdi[m_g[ok_g]]))
    }
  }
  if (exists("hydrate_step07_organellar_placement_metrics", mode = "function")) df <- hydrate_step07_organellar_placement_metrics(df, opts)
  org_sig <- if (exists("evaluate_step07_organellar_and_oop_signals", mode = "function")) evaluate_step07_organellar_and_oop_signals(df) else NULL
  geno_adm_dv <- if ("is_geno_admixed" %in% names(df)) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
  geno_qry_dv <- if ("is_putative_novel_or_hybrid_query" %in% names(df)) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  geno_hyb_cls_dv <- .chr_col(df, "hybrid_generation_class", n, NA_character_); geno_q12_dv <- .num_col(df, "interspecific_ancestry_q12", n); geno_fb_max_dv <- .num_col(df, "f_branch_max", n)
  geno_cytonuc_dv <- (if ("cytonuclear_discordant" %in% names(df)) .step07_safe_bool(df$cytonuclear_discordant) else rep(FALSE, n)) | (if (!is.null(org_sig)) org_sig$ref_cytonuc_capture else FALSE)
  geno_divergent_nonhybrid_v <- has_geno_dv & ((!is.na(geno_hyb_cls_dv) & geno_hyb_cls_dv == "Divergent_NonHybrid_Lineage") | (is.finite(geno_gdi_dv) & geno_gdi_dv > 0.70 & (!is.finite(geno_q12_dv) | geno_q12_dv < 0.15)))
  geno_contemp_hybrid_v <- has_geno_dv & !geno_divergent_nonhybrid_v & ((is.finite(geno_q12_dv) & geno_q12_dv >= 0.35) | (!is.na(geno_hyb_cls_dv) & (geno_hyb_cls_dv %in% c("F1_Hybrid", "F2_Hybrid", "BC1_Parent_1", "BC1_Parent_2", "F2_Backcross", "Recent_Admixed"))) | (geno_adm_dv & (!is.finite(geno_q12_dv) | geno_q12_dv >= 0.15)))
  geno_ancient_introg_v <- has_geno_dv & !geno_contemp_hybrid_v & !geno_divergent_nonhybrid_v & (!is.finite(geno_q12_dv) | geno_q12_dv < 0.15) & ((is.finite(geno_fb_max_dv) & geno_fb_max_dv >= 0.05) | geno_cytonuc_dv | (!is.na(geno_hyb_cls_dv) & (geno_hyb_cls_dv %in% c("Advanced_Introgressant", "Ancient_Introgressant"))))
  cluster_prop_gdi_dv <- .num_col(df, "cluster_propagated_gdi", n); tier0_abstain_dv <- if ("tier0_abstention_flag" %in% names(df)) (!is.na(df$tier0_abstention_flag) & suppressWarnings(as.integer(df$tier0_abstention_flag)) == 1L) else rep(FALSE, n)
  has_dg <- has_geno_dv & is.finite(geno_gdi_dv); eff_gdi_dv <- ifelse(has_dg, geno_gdi_dv, cluster_prop_gdi_dv); is_cluster_prop_gdi_only <- !has_dg & is.finite(cluster_prop_gdi_dv)
  cyto_b_int <- if ("cytotype_reproductive_barrier_pass" %in% names(df)) suppressWarnings(as.integer(df$cytotype_reproductive_barrier_pass)) else rep(NA_integer_, n)
  type_info_d <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else list(is_type_specimen = rep(FALSE, n), type_status_clean = rep("None", n))
  is_type_dv <- type_info_d$is_type_specimen | .lgl_col(df, "is_type_specimen", n)
  exs_gold_tx_dv <- .chr_col(df, "exsiccatae_gold_taxon", n, NA_character_); if (any(ok_eg <- which(!is.na(exs_gold_tx_dv) & nzchar(trimws(exs_gold_tx_dv))))) exs_gold_tx_dv[ok_eg] <- standardize_packera_taxon(exs_gold_tx_dv[ok_eg])
  exs_conf_r <- if ("exsiccatae_taxonomic_conflict" %in% names(df)) df$exsiccatae_taxonomic_conflict else rep(0L, n)
  p_type_morph <- max(2L, if (length(bg_ctx$sel_mc) > 0L) length(bg_ctx$sel_mc) else if (length(bg_ctx$morph_coord_cols) > 0L) length(bg_ctx$morph_coord_cols) else 3L)
  osr_df_raw <- .num_col(df, "open_set_df", n); osr_df_vec <- ifelse(is.finite(osr_df_raw) & osr_df_raw >= 1.0, as.integer(round(osr_df_raw)), p_type_morph)
  osr_chisq85_vec <- stats::qchisq(0.85, df = osr_df_vec); osr_chisq90_vec <- stats::qchisq(0.90, df = osr_df_vec); osr_chisq95_vec <- stats::qchisq(0.95, df = osr_df_vec)

  mv <- list(
    is_intra = is_intra, opts_cyto = opts_cyto, sp_all_vec = sp_all_vec, cat_sheet_vec = cat_sheet_vec, sheet_idx_list = sheet_idx_list,
    raw_sp_dv = .chr_col(df, "species_raw", n, NA_character_), tier_dv = .chr_col(df, "determiner_tier", n, NA_character_),
    cda_pred_dv = .chr_col(df, "cda_predicted_taxon", n, NA_character_), cda_sec1_dv = .chr_col(df, "cda_second_taxon", n, NA_character_),
    cda_sec2_dv = .chr_col(df, "cda_second_predicted_taxon", n, NA_character_), cda_prob_dv = .num_col(df, "cda_posterior_prob", n),
    cda_sec_p_dv = .num_col(df, "cda_second_posterior_prob", n), cda_marg_dv = .num_col(df, "cda_posterior_margin", n),
    mo_flag_dv = if ("morphometric_outlier_flag" %in% names(df)) as.logical(df$morphometric_outlier_flag) else rep(NA, n),
    misid_flag_dv = if ("misidentification_flag" %in% names(df)) as.logical(df$misidentification_flag) else rep(NA, n),
    triage_cat_dv = .chr_col(df, "triage_category", n, NA_character_), edaph_sp_dv = .chr_col(df, "edaphic_best_fit_taxon", n, NA_character_),
    rf_pred_dv = .chr_col(df, "rf_predicted_taxon", n, NA_character_), rf_prob_dv = .num_col(df, "rf_predicted_prob", n),
    d_ed1_dv = .num_col(df, "edaphic_niche_overlap_d", n), d_ed2_dv = .num_col(df, "edaphic_schoener_d", n), d_ed3_dv = .num_col(df, "schoener_d", n),
    gmm_unc_dv = .num_col(df, "gmm_uncertainty", n), gmm_cl_dv = if ("gmm_cluster" %in% names(df)) suppressWarnings(as.integer(df$gmm_cluster)) else rep(NA_integer_, n),
    ddoy1_dv = .num_col(df, "delta_doy_resid", n), ddoy2_dv = .num_col(df, "delta_doy", n), is_anth_dv = if ("is_active_anthesis" %in% names(df)) .lgl_col(df, "is_active_anthesis", n) else rep(TRUE, n),
    bf1_dv = .num_col(df, "two_delta_bic", n), bf2_dv = .num_col(df, "gmm_2_delta_bic", n), m_tier_src_dv = mts_init_v, n_t1_dv = dplyr::coalesce(bas_raw_i, if ("leaf_count" %in% names(df)) suppressWarnings(as.integer(df$leaf_count)) else rep(NA_integer_, n)),
    bas_dv = bas_raw_i, dlf_dv = .num_col(df, "dissected_leaf_fraction", n), asym_dv = .num_col(df, "pct_asymmetry", n), is_dis_dv = .lgl_col(df, "is_dissected", n),
    p_conf_dv = .num_col(df, "probabilistic_consensus_conflict", n, 0.0), orig_has_obs_morph_v = orig_has_obs_morph_v, orig_n_leaves_v = orig_n_leaves_v,
    orig_cda_pred_v = cda_pred_init, orig_cda_prob_v = cda_prob_init, orig_tm_d2_v = tm_d2_init, has_observed_or_sheet_morph = has_observed_or_sheet_morph,
    anchor_all = anchor_all, mo_all = mo_all, noise_all = noise_all, ext_non_type_all = ext_non_type_all, valid_doy_all = valid_doy_all, anth_all = anth_all, hop_doy_all = hop_doy_all,
    ph_num_all = ph_num_all, sand_num_all = sand_num_all, glob_circ = glob_circ, mu_theta_hop = mu_theta_hop, q75_arc_hop = q75_arc_hop, q75_doy = q75_doy, med_doy = med_doy,
    edaph_z_vec = edaph_z_vec, tau_ecotype_vec = pheno_taus$tau_ecotype, tau_cryptic_vec = pheno_taus$tau_cryptic,
    has_geno_dv = has_geno_dv, geno_pred_dv = geno_pred_dv, geno_prob_dv = .num_col(df, "geno_posterior_prob", n), geno_marg_dv = .num_col(df, "geno_margin", n),
    geno_tier_dv = if ("geno_linkage_tier" %in% names(df)) as.character(df$geno_linkage_tier) else .chr_col(df, "linkage_tier", n, NA_character_),
    geno_adm_dv = geno_adm_dv, geno_p1_dv = .chr_col(df, "geno_parent_1", n, NA_character_), geno_p2_dv = .chr_col(df, "geno_parent_2", n, NA_character_),
    geno_qmax_dv = .num_col(df, "q_max", n), geno_qry_dv = geno_qry_dv, geno_gdi_dv = geno_gdi_dv, eff_gdi_dv = eff_gdi_dv, is_cluster_prop_gdi_only = is_cluster_prop_gdi_only,
    tier0_abstain_dv = tier0_abstain_dv, geno_h_dv = .num_col(df, "hybrid_index_h", n), geno_q12_dv = geno_q12_dv, geno_hyb_cls_dv = geno_hyb_cls_dv,
    geno_fb_max_dv = geno_fb_max_dv, geno_cytonuc_dv = geno_cytonuc_dv, geno_cyt_dv = .chr_col(df, "inferred_cytotype", n, NA_character_),
    geno_divergent_nonhybrid_v = geno_divergent_nonhybrid_v, geno_contemp_hybrid_v = geno_contemp_hybrid_v, geno_ancient_introg_v = geno_ancient_introg_v,
    geno_gdi_oversplit_v = is.finite(eff_gdi_dv) & (eff_gdi_dv < 0.20),
    geno_gdi_cryptic_support_v = (has_geno_dv & ((is.finite(geno_gdi_dv) & geno_gdi_dv > 0.70) | (!is.na(geno_hyb_cls_dv) & geno_hyb_cls_dv == "Divergent_NonHybrid_Lineage"))) | (is_cluster_prop_gdi_only & (cluster_prop_gdi_dv > 0.70) & !tier0_abstain_dv),
    cyto_barrier_pass_v = !is.na(cyto_b_int) & (cyto_b_int == 1L), voucher_ploidy_v = if ("voucher_ploidy_level" %in% names(df)) suppressWarnings(as.integer(df$voucher_ploidy_level)) else rep(NA_integer_, n),
    taxon_cyto_str_v = .chr_col(df, "taxon_reported_cytotypes", n, NA_character_), pid_cyto_dv = if ("plant_individual_id" %in% names(df)) suppressWarnings(as.integer(df$plant_individual_id)) else rep(0L, n),
    is_type_dv = is_type_dv, type_stat_dv = if ("type_status" %in% names(df)) as.character(df$type_status) else type_info_d$type_status_clean,
    is_typ_morph_dv = if ("is_typical_morphotype" %in% names(df)) as.logical(df$is_typical_morphotype) else rep(NA, n),
    tm_d2_dv = tm_d2_init, nt_morph_dv = .chr_col(df, "nearest_type_morph_taxon", n, NA_character_), nt_bas_dv = .chr_col(df, "nearest_type_basionym", n, NA_character_),
    bas_d2_dv = .num_col(df, "nearest_type_basionym_d2", n), tm_conc_dv = .num_col(df, "type_morph_concordance", n), te_z_dv = .num_col(df, "type_edaphic_z_score", n),
    nt_niche_dv = .chr_col(df, "nearest_type_niche_taxon", n, NA_character_), min_cong_d2_dv = .num_col(df, "min_congener_mahal_d2", n),
    min_cong_tx_dv = .chr_col(df, "min_congener_taxon", n, NA_character_), osr_p_dv = .num_col(df, "open_set_novelty_p", n), osr_df_vec = osr_df_vec,
    self_d2_dv = .num_col(df, "self_morph_mahal_d2", n), self_df_raw = .num_col(df, "self_morph_df", n), cda_oof_p_dv = .num_col(df, "cda_oof_posterior_prob", n),
    cl_thr_dv = .num_col(df, "cl_target_threshold", n), exs_gold_tx_dv = exs_gold_tx_dv,
    usda_stat_dv = tolower(trimws(ifelse(is.na(.chr_col(df, "usda_county_status", n, "")), "", .chr_col(df, "usda_county_status", n, "")))),
    exs_conf_dv = if (is.logical(exs_conf_r)) (!is.na(exs_conf_r) & exs_conf_r) else (!is.na(suppressWarnings(as.integer(exs_conf_r))) & suppressWarnings(as.integer(exs_conf_r)) != 0L),
    clim_z_dv = .num_col(df, "climatic_z_score", n), rf_given_p_dv = .num_col(df, "rf_given_taxon_prob", n), id_by_dv = .chr_col(df, "identifiedBy", n, NA_character_),
    portal_lag_flag = logical(n), portal_lag_gate = rep(NA_character_, n), p_type_morph = p_type_morph, type_morph_chisq95 = stats::qchisq(0.95, df = p_type_morph),
    osr_chisq85_vec = osr_chisq85_vec, osr_chisq90_vec = osr_chisq90_vec, osr_chisq95_vec = osr_chisq95_vec,
    n_sheet_ramets_vec = as.integer(table(cat_sheet_vec)[cat_sheet_vec])
  )
  cl_ctx <- .build_step07_cluster_cohesion_context(df, bg_ctx, mv, gmm_two_delta_bic = gmm_two_delta_bic, is_intra = is_intra)
  st_res <- .evaluate_step07_specimen_stages(df, mv, bg_ctx, cl_ctx, org_sig, .eval_five_gate_hybrid_pair, gmm_two_delta_bic = gmm_two_delta_bic, is_intra = is_intra)
  outside_self_morph_95 <- st_res$outside_self_morph_95; tier_authority_gate_ok <- st_res$tier_authority_gate_ok
  df <- .apply_step07_postprocess_and_guards(df, st_res, mv, bg_ctx, cl_ctx, org_sig, opts)

  df$cryptic_non_clinal_pass <- as.logical(df$cryptic_non_clinal_pass)
  df$resurrection_basionym <- as.character(df$resurrection_basionym)
  df$hybrid_zone_cohesion_pass <- as.logical(df$hybrid_zone_cohesion_pass)
  df$hybrid_parent_1 <- as.character(df$hybrid_parent_1); df$hybrid_parent_2 <- as.character(df$hybrid_parent_2)
  attr(df, "can_evaluate_hybrids") <- can_evaluate_hybrids
  attr(df, "n_resurrection_candidates") <- sum(df$taxonomic_status_call == "Taxonomic_Resurrection_Candidate", na.rm = TRUE)
  if (!is.null(prop_sum_attr)) attr(df, "genomic_propagation_summary") <- prop_sum_attr
  if (!is.null(cyto_md_attr)) attr(df, "cytology_markdown_lines") <- cyto_md_attr
  if (!is.null(cyto_bar_lbl_attr)) attr(df, "cytotype_barrier_labels") <- cyto_bar_lbl_attr
  if (!is.null(gigas_df_attr)) attr(df, "gigas_concordance_df") <- gigas_df_attr
  df
}

main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  set.seed(42); opts <- parse_args_robust(cli_args)
  message("=== Starting Step 07 Triage Dashboard & Taxonomic Synthesis ===")
  triage_df <- build_triage_queue(opts)
  render_synthesis_figure_png(triage_df, opts$output_plot, cda_summary = opts$cda_summary, focal_taxon = opts$focal_taxon, output_complex_plot = opts$output_complex_plot)
  tryCatch({
    .source_step07_dep("scripts/morphometrics/pub_figures_utils.R", "../morphometrics/pub_figures_utils.R", "generate_step07_pub_figures")
    if (exists("generate_step07_pub_figures", mode = "function")) generate_step07_pub_figures(triage_df, opts$output_plot)
  }, error = function(e) warning("Step 07 publication figures export warning: ", e$message))
  generate_taxonomic_revision_report(triage_df, opts)
  if (!is.null(opts$output_html) && nzchar(trimws(opts$output_html))) {
    generate_interactive_triage_html(triage_df, opts$output_html, focal_taxon = opts$focal_taxon, opts = opts)
  }
  message("=== Taxonomic Synthesis & Triage Dashboard Complete ===")
}

if (sys.nframe() == 0L) main()
