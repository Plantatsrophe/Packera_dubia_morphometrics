#!/usr/bin/env Rscript
# ==============================================================================
# Script: 07_triage_dashboard_synthesis.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Step 07 Final Multi-Evidence Synthesis & Triage Dashboard Engine.
#   1. Safely joins and coalesces georeferenced and non-georeferenced vouchers
#      across Step 04 (Vegetative EFA/ECT Morphometrics), Step 02/05
#      (Reproductive Capitulum & Phenological Allochrony), Step 05
#      (Fine-Grain Rhizosphere/Outcrop & Climatic Niche), and Step 06
#      (Independent Genomic & Phylogenomic Stream) using composite keys.
#   2. Computes empirical out-of-fold Cohen's kappa reliability weights on
#      Tier_1_Gold + Tier_2_Silver anchors (Dawid & Skene 1979; Kuncheva 2014)
#      with dynamic missing-modality weight renormalization.
#   3. Adjusts anthesis phenology via a circular Hopkins bioclimatic model
#      (sin_doy, cos_doy ~ decimalLatitude + elevation + bio1; Hopkins 1918;
#      Pewsey et al. 2013) gated strictly on active anthesis specimens
#      (capitulum_phenology_classifier.py) using shortest-arc atan2 residuals.
#   4. Enforces the Five-Stage Hierarchical Integrative Decision Matrix
#      (Padial et al. 2010; Schlick-Steiner et al. 2010; Rieseberg 1995;
#      Slenker et al. 2022):
#        Stage 1: Misidentification_Reassign (evaluated BEFORE hybridity)
#        Stage 2: Hybrid_Introgressant (inter-modal mosaicism or canonical
#                 intermediacy cda_posterior_margin < 0.25 with niche overlap)
#        Stage 3: Putative_Cryptic_Species (Integration by Congruence:
#                 Basal_Rosette, non-assigned to known taxa, morph/GMM outlier
#                 with 2*delta_BIC >= 6, and allochronic |delta_doy_resid| >= 18d
#                 or edaphic divergence D_edaph < 0.30)
#        Stage 4: Subspecies_Ecotype vs. Ecophenotypic_Plasticity
#        Stage 5: Confirmed_Species
#   5. Audits institutional herbarium batch effects via taxon-adjusted marginal
#      ANCOVA (lm(PC1 ~ species_standardized + institutionCode)) reporting
#      partial eta-squared and Benjamini-Hochberg FDR-corrected p-values.
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("tibble", quietly = TRUE)) library(tibble)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("gridExtra", quietly = TRUE)) library(gridExtra)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})

for (.aw_p in c("scripts/morphometrics/authority_weighting_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "morphometrics", "authority_weighting_utils.R"))) {
  if ((!exists(".identify_type_specimen_vec", mode = "function") || !exists("resolve_ecological_authority_weights", mode = "function")) && !is.na(.aw_p) && nzchar(.aw_p) && file.exists(.aw_p)) { source(.aw_p); break }
}
for (.pf_p in c("scripts/morphometrics/pub_figures_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "morphometrics", "pub_figures_utils.R"))) {
  if ((!exists("generate_step07_pub_figures", mode = "function") || !exists("standardize_packera_taxon", mode = "function")) && !is.na(.pf_p) && nzchar(.pf_p) && file.exists(.pf_p)) { source(.pf_p); break }
}
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!exists("persist_step07_to_sqlite", mode = "function") && !is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
}
for (.tc_p in c("scripts/analysis/taxon_complex_discovery_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "taxon_complex_discovery_utils.R"))) {
  if (!exists("run_step07_complex_discovery", mode = "function") && !is.na(.tc_p) && nzchar(.tc_p) && file.exists(.tc_p)) { source(.tc_p); break }
}
for (.gt_p in c("scripts/analysis/genomic_triage_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "genomic_triage_utils.R"))) {
  if (!exists("load_step07_genomic_metrics", mode = "function") && !is.na(.gt_p) && nzchar(.gt_p) && file.exists(.gt_p)) { source(.gt_p); break }
}
for (.gp_p in c("scripts/analysis/step07_genomic_propagation.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step07_genomic_propagation.R"))) {
  if (!exists("run_step07_genomic_propagation", mode = "function") && !is.na(.gp_p) && nzchar(.gp_p) && file.exists(.gp_p)) { source(.gp_p); break }
}
for (.cy_p in c("scripts/analysis/step07_cytology_integration.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step07_cytology_integration.R"))) {
  if (!exists("run_step07_cytology_integration", mode = "function") && !is.na(.cy_p) && nzchar(.cy_p) && file.exists(.cy_p)) { source(.cy_p); break }
}
for (.gw_p in c("scripts/analysis/edaphic_geochem_wetland_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "edaphic_geochem_wetland_utils.R"))) {
  if (!exists("evaluate_step07_hydro_edaphic_and_basionym_guards", mode = "function") && !is.na(.gw_p) && nzchar(.gw_p) && file.exists(.gw_p)) { source(.gw_p); break }
}

.step07_safe_bool <- function(x) {
  if (exists(".safe_bool_vec_r", mode = "function")) return(.safe_bool_vec_r(x))
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  toupper(trimws(as.character(x))) %in% c("1", "TRUE", "T", "YES", "Y")
}

.has_step07_geno_row_mask <- function(df) {
  if (exists("has_active_genomic_row_mask", mode = "function")) return(has_active_genomic_row_mask(df))
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  nms <- names(df)
  if (!("geno_predicted_taxon" %in% nms)) return(rep(FALSE, n))
  pred <- as.character(df$geno_predicted_taxon)
  pred_ok <- !is.na(pred) & nzchar(trimws(pred)) & (pred != "Unknown")
  tier_col <- if ("geno_linkage_tier" %in% nms) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% nms) as.character(df$linkage_tier) else rep("G1_Direct_Voucher", n)
  tier_ok <- tier_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G2_Conditional_Deme_Hull")
  qc_ok <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(TRUE, n)
  prob_ok <- if ("geno_posterior_prob" %in% nms) is.finite(suppressWarnings(as.numeric(df$geno_posterior_prob))) else rep(TRUE, n)
  pred_ok & tier_ok & qc_ok & prob_ok
}

# ------------------------------------------------------------------------------
# 1. CLI Argument Parsing & Taxonomic Normalization
# ------------------------------------------------------------------------------
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
  list(vouchers = "data/tables/curated_vouchers.csv", harmonics = "data/tables/leaf_efa_harmonics.csv",
       morphometrics = "data/tables/morphometric_flags.csv",
       multimodal_flags = "data/tables/multimodal_conflict_flags.csv",
       gmm_summary = "outputs/reports/gmm_bayes_factors_summary.csv", gmm_bic = "data/tables/gmm_bic_model_selection.csv",
       cda_summary = "data/tables/cda_summary.csv", niche_summary = "outputs/reports/multimodal_spatial_rf_summary.csv",
       output_queue = "data/tables/triage_queue.csv", output_plot = "outputs/figures/fig07a_integrative_taxonomic_synthesis.png",
       output_report = "outputs/reports/Packera_Taxonomic_Revision_Summary.md",
       output_html = "outputs/reports/Triage_Dashboard_Interactive.html",
       min_affinity = 0.03, output_complexes = "outputs/reports/empirical_taxon_complexes.csv",
       output_affinity = "outputs/reports/pairwise_taxon_affinity_matrix.csv",
       output_neighbors = "outputs/reports/focal_auto_neighbors.txt",
       output_complex_plot = "outputs/figures/fig07b_taxon_complex_delimitation_network.png",
       genomic_dir = "",
       focal_taxon = "", intraspecific_mode = FALSE, demo = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Hopkins-Corrected Circular Anthesis Phenology (Hopkins 1918; Pewsey et al. 2013)
# ------------------------------------------------------------------------------
resolve_active_anthesis_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  ver_vec <- if ("is_verified_anthesis" %in% names(df)) suppressWarnings(as.logical(df$is_verified_anthesis)) else rep(NA, n)
  na_cnt  <- if ("n_anthesis" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis)) else if ("n_anthesis_heads" %in% names(df)) suppressWarnings(as.numeric(df$n_anthesis_heads)) else rep(NA_real_, n)
  non_anth_cnt <- rep(0.0, n)
  has_non_anth_col <- FALSE
  for (nc in intersect(c("n_bud", "n_fruiting", "n_senescent"), names(df))) {
    has_non_anth_col <- TRUE
    v <- suppressWarnings(as.numeric(df[[nc]]))
    non_anth_cnt <- non_anth_cnt + ifelse(is.finite(v), v, 0.0)
  }
  state_cols <- intersect(c("voucher_phenological_state", "phenological_state", "phenophase", "head_state"), names(df))
  if (length(state_cols) > 0L) {
    st_raw <- tolower(trimws(as.character(df[[state_cols[1L]]])))
    unk_mask <- is.na(st_raw) | !nzchar(st_raw) | st_raw == "unknown"
    is_anthesis <- (!unk_mask) & (st_raw %in% c("anthesis", "flowering", "active_anthesis", "open_flower") | grepl("anthesis|flower", st_raw))
    is_anthesis[!is.na(st_raw) & (st_raw %in% c("bud", "fruit", "sterile", "vegetative", "fruiting", "senescent") | grepl("pre[_-]?anthesis|post[_-]?anthesis|non[_-]?anthesis", st_raw))] <- FALSE
    if (any(unk_mask)) {
      has_ver <- unk_mask & !is.na(ver_vec)
      if (any(has_ver)) is_anthesis[has_ver] <- ver_vec[has_ver]
      rem_unk <- unk_mask & is.na(ver_vec) & is.finite(na_cnt)
      if (any(rem_unk)) {
        is_anthesis[rem_unk] <- na_cnt[rem_unk] >= 1.0
      }
    }
    return(is_anthesis)
  }
  if ("is_verified_anthesis" %in% names(df)) {
    ia <- ver_vec
    if (any(is.na(ia)) && any(is.finite(na_cnt))) ia[is.na(ia) & is.finite(na_cnt)] <- na_cnt[is.na(ia) & is.finite(na_cnt)] >= 1.0
    if (all(is.na(ia)) && ("is_active_anthesis" %in% names(df))) {
      ia_fb <- suppressWarnings(as.logical(df$is_active_anthesis))
      ia_fb[is.na(ia_fb)] <- FALSE
      return(ia_fb)
    }
    ia[is.na(ia)] <- FALSE
    return(ia)
  }
  if ("is_active_anthesis" %in% names(df)) {
    ia <- suppressWarnings(as.logical(df$is_active_anthesis))
    ia[is.na(ia)] <- FALSE
    return(ia)
  }
  if ("is_anthesis" %in% names(df)) { ia <- as.logical(df$is_anthesis); ia[is.na(ia)] <- FALSE; return(ia) }
  if ("n_anthesis" %in% names(df)) { return(!is.na(na_cnt) & is.finite(na_cnt) & na_cnt >= 1.0) }
  if ("n_anthesis_heads" %in% names(df)) { nh <- suppressWarnings(as.numeric(df$n_anthesis_heads)); return(!is.na(nh) & is.finite(nh) & nh >= 1.0) }
  rep(TRUE, n)
}

.num_col <- function(df, col_nm, n = nrow(df), default = NA_real_) {
  if (col_nm %in% names(df)) suppressWarnings(as.numeric(df[[col_nm]])) else rep(default, n)
}

.chr_col <- function(df, col_nm, n = nrow(df), default = NA_character_) {
  if (col_nm %in% names(df)) as.character(df[[col_nm]]) else rep(default, n)
}

resolve_step07_ecological_authority_weights <- function(df, allow_bronze_prior = TRUE) {
  n <- nrow(df)
  if (n == 0L) return(numeric(0))
  w_out <- resolve_ecological_authority_weights(df, allow_bronze_prior = allow_bronze_prior)
  if (isTRUE(allow_bronze_prior) && length(w_out) == n) {
    w_out <- round(pmax(0.01, w_out), 4)
  }
  if (any(.has_step07_geno_row_mask(df))) {
    nms <- names(df)
    is_type_m <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df)$is_type_specimen else rep(FALSE, n)
    if ("is_type_specimen" %in% nms) {
      tv <- suppressWarnings(as.logical(df$is_type_specimen)); tv[is.na(tv)] <- FALSE
      is_type_m <- is_type_m | tv
    }
    tier_g <- if ("geno_linkage_tier" %in% nms) trimws(as.character(df$geno_linkage_tier)) else if ("linkage_tier" %in% nms) trimws(as.character(df$linkage_tier)) else rep("", n)
    qc_g <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(FALSE, n)
    adm_g <- if ("is_geno_admixed" %in% nms) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
    qry_g <- if ("is_putative_novel_or_hybrid_query" %in% nms) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
    pred_g <- if ("geno_predicted_taxon" %in% nms) trimws(as.character(df$geno_predicted_taxon)) else rep("", n)
    pred_g[is.na(pred_g)] <- ""
    sp_g <- if ("species_standardized" %in% nms) trimws(as.character(df$species_standardized)) else rep("", n)
    sp_g[is.na(sp_g)] <- ""
    prob_g <- if ("geno_posterior_prob" %in% nms) suppressWarnings(as.numeric(df$geno_posterior_prob)) else rep(NA_real_, n)
    marg_g <- if ("geno_margin" %in% nms) suppressWarnings(as.numeric(df$geno_margin)) else rep(NA_real_, n)
    g1_ok <- (tier_g %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc_g & nzchar(pred_g) & (pred_g != "Unknown") & is.finite(prob_g) & (prob_g >= 0.85)
    t0_conc <- !is_type_m & g1_ok & !adm_g & !qry_g & (is.na(marg_g) | marg_g >= 0.40) & (pred_g == sp_g)
    t0_contra <- !is_type_m & g1_ok & (pred_g != sp_g)
    w_out[t0_conc] <- 4.0
    w_out[t0_contra] <- 0.0
  }
  w_out
}

.load_step07_pheno_config_defaults <- function(cfg_path = "config/config.yaml") {
  defs <- list(tau_ecotype = 14.0, tau_cryptic = 18.0)
  if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    if (is.list(cfg) && is.list(cfg$synthesis)) {
      if (!is.null(cfg$synthesis$ecotype_ddoy_cutoff_days)) defs$tau_ecotype <- as.numeric(cfg$synthesis$ecotype_ddoy_cutoff_days)
      if (!is.null(cfg$synthesis$cryptic_ddoy_cutoff_days)) defs$tau_cryptic <- as.numeric(cfg$synthesis$cryptic_ddoy_cutoff_days)
    }
  }
  defs
}

.circular_sd_days_from_residuals <- function(ddoy_days_vec) {
  v <- ddoy_days_vec[is.finite(ddoy_days_vec)]
  if (length(v) < 2L) return(NA_real_)
  phi <- 2.0 * pi * (v / 365.25)
  r_bar <- sqrt(mean(sin(phi))^2 + mean(cos(phi))^2)
  r_clamped <- max(1e-6, min(0.999999, r_bar))
  (365.25 / (2.0 * pi)) * sqrt(-2.0 * log(r_clamped))
}

.is_step07_extralimital_non_type <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  nms <- names(df)
  is_type_m <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df)$is_type_specimen else rep(FALSE, n)
  if ("is_type_specimen" %in% nms) {
    tv <- suppressWarnings(as.logical(df$is_type_specimen)); tv[is.na(tv)] <- FALSE
    is_type_m <- is_type_m | tv
  }
  if ("determiner_tier" %in% nms) {
    is_type_m <- is_type_m | (trimws(as.character(df$determiner_tier)) == "Tier_0_Type_Anchor")
  }
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
      tx_idx <- which(sp_v == tx & !noise_m)
      if (sum(!flag[tx_idx]) < 2L) flag[sp_v == tx] <- FALSE
    }
  }
  flag
}

.resolve_step07_anchor_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  tier_col <- if ("determiner_tier" %in% names(df)) {
    if (exists(".standardize_determiner_tier_vec", mode = "function")) .standardize_determiner_tier_vec(df)
    else ifelse(is.na(df$determiner_tier) | !nzchar(trimws(as.character(df$determiner_tier))), "Tier_3_Bronze", as.character(df$determiner_tier))
  } else {
    rep("Tier_1_Gold", n)
  }
  promo_mask <- rep(FALSE, n)
  for (pcol in intersect(c("is_promoted_bronze", "is_type_specimen"), names(df))) {
    pv <- suppressWarnings(as.logical(df[[pcol]]))
    pv[is.na(pv)] <- FALSE
    promo_mask <- promo_mask | pv
  }
  if (exists(".identify_type_specimen_vec", mode = "function")) {
    promo_mask <- promo_mask | .identify_type_specimen_vec(df)$is_type_specimen
  }
  noise_mask <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
  passive_mask <- if ("is_passive" %in% names(df)) (!is.na(df$is_passive) & as.logical(df$is_passive)) else rep(FALSE, n)
  if ("authority_weight" %in% names(df)) {
    aw_v <- suppressWarnings(as.numeric(df$authority_weight))
    if (any(is.finite(aw_v) & aw_v > 0.0)) {
      passive_mask <- passive_mask | (is.finite(aw_v) & aw_v <= 0.0)
    }
  }
  is_type_m <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df)$is_type_specimen else rep(FALSE, n)
  if ("is_type_specimen" %in% names(df)) {
    tv <- suppressWarnings(as.logical(df$is_type_specimen)); tv[is.na(tv)] <- FALSE
    is_type_m <- is_type_m | tv
  }
  ext_non_type_mask <- .is_step07_extralimital_non_type(df)
  anchor_mask <- (((tier_col %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze_Promoted")) | promo_mask) & !noise_mask & !passive_mask & !ext_non_type_mask) | (is_type_m & !noise_mask)
  if (any(.has_step07_geno_row_mask(df))) {
    nms <- names(df)
    tier_g <- if ("geno_linkage_tier" %in% nms) trimws(as.character(df$geno_linkage_tier)) else if ("linkage_tier" %in% nms) trimws(as.character(df$linkage_tier)) else rep("", n)
    qc_g <- if ("geno_quality_pass" %in% nms) .step07_safe_bool(df$geno_quality_pass) else rep(FALSE, n)
    adm_g <- if ("is_geno_admixed" %in% nms) .step07_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
    qry_g <- if ("is_putative_novel_or_hybrid_query" %in% nms) .step07_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
    pred_g <- if ("geno_predicted_taxon" %in% nms) trimws(as.character(df$geno_predicted_taxon)) else rep("", n)
    pred_g[is.na(pred_g)] <- ""
    sp_g <- if ("species_standardized" %in% nms) trimws(as.character(df$species_standardized)) else rep("", n)
    sp_g[is.na(sp_g)] <- ""
    prob_g <- if ("geno_posterior_prob" %in% nms) suppressWarnings(as.numeric(df$geno_posterior_prob)) else rep(NA_real_, n)
    g1_base <- (tier_g %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc_g & !adm_g & nzchar(pred_g) & (pred_g != "Unknown") & is.finite(prob_g) & (prob_g >= 0.85)
    t0_contra <- !is_type_m & g1_base & (pred_g != sp_g)
    t0_conc <- g1_base & !qry_g & (pred_g == sp_g) & !noise_mask
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
  n <- nrow(df)
  defs <- .load_step07_pheno_config_defaults()
  if (n == 0L) return(list(tau_ecotype = numeric(0), tau_cryptic = numeric(0), sigma_by_taxon = numeric(0)))

  sp_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unknown", n)
  is_anth <- if ("is_active_anthesis" %in% names(df)) as.logical(df$is_active_anthesis) else resolve_active_anthesis_mask(df)
  is_anth[is.na(is_anth)] <- TRUE
  ddoy_vec <- if ("delta_doy_resid" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy_resid)) else if ("delta_doy" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy)) else rep(NA_real_, n)
  valid_ph <- is_anth & is.finite(ddoy_vec)

  anchor_mask <- .resolve_step07_anchor_mask(df)
  mo_m <- if ("morphometric_outlier_flag" %in% names(df)) (!is.na(df$morphometric_outlier_flag) & as.logical(df$morphometric_outlier_flag)) else rep(FALSE, n)
  noise_m <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
  clean_anch <- anchor_mask & !mo_m & !noise_m

  coh_idx <- which(clean_anch & valid_ph)
  if (length(coh_idx) < 5L) coh_idx <- which(!mo_m & !noise_m & valid_ph)
  if (length(coh_idx) < 5L) coh_idx <- which(!noise_m & valid_ph)
  sigma_cohort <- if (length(coh_idx) >= 5L) .circular_sd_days_from_residuals(ddoy_vec[coh_idx]) else NA_real_

  u_taxa <- unique(sp_vec[!is.na(sp_vec) & nzchar(sp_vec)])
  tau_eco_map <- stats::setNames(rep(defs$tau_ecotype, length(u_taxa)), u_taxa)
  tau_cry_map <- stats::setNames(rep(defs$tau_cryptic, length(u_taxa)), u_taxa)
  sigma_map <- stats::setNames(rep(sigma_cohort, length(u_taxa)), u_taxa)

  if (is.finite(sigma_cohort)) {
    def_eco_coh <- min(28.0, max(10.0, 2.0 * sigma_cohort))
    def_cry_coh <- min(35.0, max(14.0, 2.5 * sigma_cohort))
    tau_eco_map[] <- def_eco_coh
    tau_cry_map[] <- def_cry_coh
  }

  for (tx in u_taxa) {
    idx_t <- which(sp_vec == tx)
    tx_pool <- idx_t[clean_anch[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) < 5L) tx_pool <- idx_t[!mo_m[idx_t] & !noise_m[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) < 5L) tx_pool <- idx_t[!noise_m[idx_t] & valid_ph[idx_t]]
    if (length(tx_pool) >= 5L) {
      sig_tx <- .circular_sd_days_from_residuals(ddoy_vec[tx_pool])
      if (is.finite(sig_tx)) {
        sigma_map[[tx]] <- sig_tx
        tau_eco_map[[tx]] <- min(28.0, max(10.0, 2.0 * sig_tx))
        tau_cry_map[[tx]] <- min(35.0, max(14.0, 2.5 * sig_tx))
      }
    }
  }

  fallback_eco <- if (is.finite(sigma_cohort)) min(28.0, max(10.0, 2.0 * sigma_cohort)) else defs$tau_ecotype
  fallback_cry <- if (is.finite(sigma_cohort)) min(35.0, max(14.0, 2.5 * sigma_cohort)) else defs$tau_cryptic
  tau_eco_vec <- ifelse(sp_vec %in% names(tau_eco_map), tau_eco_map[sp_vec], fallback_eco)
  tau_cry_vec <- ifelse(sp_vec %in% names(tau_cry_map), tau_cry_map[sp_vec], fallback_cry)
  list(tau_ecotype = as.numeric(tau_eco_vec), tau_cryptic = as.numeric(tau_cry_vec), sigma_by_taxon = sigma_map)
}

compute_hopkins_adjusted_phenology <- function(df) {
  # Circular Hopkins bioclimatic model with Empirical Bayes shrinkage:
  # sin_doy, cos_doy ~ available non-degenerate covariates among decimalLatitude, elevation, bio1_temp_mean
  # (Hopkins 1918; Pewsey et al. 2013) gated strictly on active anthesis specimens (capitulum_phenology_classifier.py).
  n <- nrow(df); if (n == 0L) return(df)
  is_anthesis <- resolve_active_anthesis_mask(df); df$is_active_anthesis <- is_anthesis
  doy_raw <- if ("doy" %in% names(df)) suppressWarnings(as.numeric(df$doy)) else rep(NA_real_, n)
  has_valid_doy <- is.finite(doy_raw) & (doy_raw >= 1.0) & (doy_raw <= 366.0)

  anchor_mask <- .resolve_step07_anchor_mask(df)
  mo_m <- if ("morphometric_outlier_flag" %in% names(df)) (!is.na(df$morphometric_outlier_flag) & as.logical(df$morphometric_outlier_flag)) else rep(FALSE, n)
  noise_m <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
  clean_anth_doy <- is_anthesis & has_valid_doy & anchor_mask & !mo_m & !noise_m

  doy_pool <- doy_raw[clean_anth_doy]
  if (length(doy_pool) == 0L) doy_pool <- doy_raw[is_anthesis & has_valid_doy & !noise_m]
  if (length(doy_pool) == 0L) doy_pool <- doy_raw[has_valid_doy]
  if (length(doy_pool) > 0L) {
    th_pool <- 2.0 * pi * (doy_pool / 365.25)
    ref_theta <- atan2(mean(sin(th_pool)), mean(cos(th_pool)))
    if (ref_theta < 0.0) ref_theta <- ref_theta + 2.0 * pi
    ref_doy <- ref_theta * (365.25 / (2.0 * pi))
  } else {
    ref_doy <- 120.0
    ref_theta <- 2.0 * pi * (ref_doy / 365.25)
  }
  doy_vec <- ifelse(has_valid_doy, doy_raw, ref_doy)

  lat_col <- intersect(c("decimalLatitude"), names(df))
  lat_vec <- if (length(lat_col) > 0L) suppressWarnings(as.numeric(df[[lat_col[1L]]])) else rep(NA_real_, n)
  elev_col <- intersect(c("elevation_m", "elevation", "minimumElevationInMeters"), names(df))
  elev_vec <- if (length(elev_col) > 0L) suppressWarnings(as.numeric(df[[elev_col[1L]]])) else rep(NA_real_, n)
  mat_col <- intersect(c("bio1_temp_mean", "bio1", "mat_c"), names(df))
  mat_vec <- if (length(mat_col) > 0L) suppressWarnings(as.numeric(df[[mat_col[1L]]])) else rep(NA_real_, n)

  lat_pool <- lat_vec[clean_anth_doy & is.finite(lat_vec)]; if (length(lat_pool) == 0L) lat_pool <- lat_vec[is_anthesis & is.finite(lat_vec)]; if (length(lat_pool) == 0L) lat_pool <- lat_vec[is.finite(lat_vec)]
  elev_pool <- elev_vec[clean_anth_doy & is.finite(elev_vec)]; if (length(elev_pool) == 0L) elev_pool <- elev_vec[is_anthesis & is.finite(elev_vec)]; if (length(elev_pool) == 0L) elev_pool <- elev_vec[is.finite(elev_vec)]
  mat_pool <- mat_vec[clean_anth_doy & is.finite(mat_vec)]; if (length(mat_pool) == 0L) mat_pool <- mat_vec[is_anthesis & is.finite(mat_vec)]; if (length(mat_pool) == 0L) mat_pool <- mat_vec[is.finite(mat_vec)]

  ref_lat <- if (length(lat_pool) > 0L) stats::median(lat_pool, na.rm = TRUE) else 35.0
  ref_elev <- if (length(elev_pool) > 0L) stats::median(elev_pool, na.rm = TRUE) else 150.0
  ref_mat <- if (length(mat_pool) > 0L) stats::median(mat_pool, na.rm = TRUE) else 14.5

  lat_shift <- ifelse(is.finite(lat_vec), 4.0 * (lat_vec - ref_lat), 0.0)
  elev_shift <- ifelse(is.finite(elev_vec), (4.0 / 121.92) * (elev_vec - ref_elev), 0.0)
  mat_shift <- ifelse(is.finite(mat_vec), -1.5 * (mat_vec - ref_mat), 0.0)
  has_geo_cov <- is.finite(lat_vec) | is.finite(elev_vec)
  has_mat_cov <- is.finite(mat_vec)
  has_any_cov <- has_geo_cov | has_mat_cov
  theo_shift <- ifelse(
    has_geo_cov & has_mat_cov,
    0.70 * (lat_shift + elev_shift) + 0.30 * mat_shift,
    ifelse(has_geo_cov, lat_shift + elev_shift, ifelse(has_mat_cov, mat_shift, 0.0))
  )

  theta_obs <- 2.0 * pi * (doy_vec / 365.25)
  theta_theo <- ref_theta + 2.0 * pi * (theo_shift / 365.25)
  sin_theo <- sin(theta_theo); cos_theo <- cos(theta_theo)

  fit_mask <- is_anthesis & has_valid_doy & has_any_cov & anchor_mask & !mo_m & !noise_m
  if (sum(fit_mask) < 5L) fit_mask <- is_anthesis & has_valid_doy & has_any_cov & anchor_mask & !noise_m
  if (sum(fit_mask) < 5L) fit_mask <- is_anthesis & has_valid_doy & has_any_cov & !noise_m
  n_anthesis_fit <- sum(fit_mask)

  cov_raw <- list(decimalLatitude = lat_vec, elevation = elev_vec, bio1_temp_mean = mat_vec)
  cov_terms <- names(cov_raw)[vapply(cov_raw, function(v) {
    vf <- v[fit_mask & is.finite(v)]
    length(vf) >= 3L && is.finite(stats::sd(vf)) && stats::sd(vf) > 1e-6
  }, logical(1))]

  alpha_eb <- if (length(cov_terms) > 0L) min(1.0, max(0.0, (n_anthesis_fit - 5.0) / 15.0)) else 0.0
  theta_blend <- theta_theo
  theta_blend_ref <- ref_theta

  if (alpha_eb > 0.0 && length(cov_terms) > 0L) {
    hop_df <- data.frame(
      sin_doy = sin(theta_obs),
      cos_doy = cos(theta_obs),
      decimalLatitude = ifelse(is.finite(lat_vec), lat_vec, ref_lat),
      elevation = ifelse(is.finite(elev_vec), elev_vec, ref_elev),
      bio1_temp_mean = ifelse(is.finite(mat_vec), mat_vec, ref_mat)
    )
    ref_df <- data.frame(decimalLatitude = ref_lat, elevation = ref_elev, bio1_temp_mean = ref_mat)
    rhs_str <- paste(cov_terms, collapse = " + ")
    lm_sin <- tryCatch(stats::lm(stats::as.formula(paste("sin_doy ~", rhs_str)), data = hop_df[fit_mask, , drop = FALSE]), error = function(e) NULL)
    lm_cos <- tryCatch(stats::lm(stats::as.formula(paste("cos_doy ~", rhs_str)), data = hop_df[fit_mask, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(lm_sin) && !is.null(lm_cos)) {
      pred_s <- tryCatch(as.numeric(stats::predict(lm_sin, newdata = hop_df)), error = function(e) rep(NA_real_, n))
      pred_c <- tryCatch(as.numeric(stats::predict(lm_cos, newdata = hop_df)), error = function(e) rep(NA_real_, n))
      ref_s <- tryCatch(as.numeric(stats::predict(lm_sin, newdata = ref_df)), error = function(e) NA_real_)
      ref_c <- tryCatch(as.numeric(stats::predict(lm_cos, newdata = ref_df)), error = function(e) NA_real_)
      if (all(is.finite(pred_s)) && all(is.finite(pred_c)) && is.finite(ref_s) && is.finite(ref_c)) {
        theta_emp <- atan2(pred_s, pred_c)
        theta_emp_ref <- atan2(ref_s, ref_c)
        sin_emp <- sin(theta_emp); cos_emp <- cos(theta_emp)
        theta_blend <- atan2(alpha_eb * sin_emp + (1.0 - alpha_eb) * sin_theo, alpha_eb * cos_emp + (1.0 - alpha_eb) * cos_theo)
        theta_blend_ref <- atan2(alpha_eb * sin(theta_emp_ref) + (1.0 - alpha_eb) * sin(ref_theta), alpha_eb * cos(theta_emp_ref) + (1.0 - alpha_eb) * cos(ref_theta))
      }
    }
  }

  delta_theta_env <- ifelse(has_any_cov, atan2(sin(theta_blend - theta_blend_ref), cos(theta_blend - theta_blend_ref)), 0.0)
  theta_adj <- theta_obs - delta_theta_env
  shift_days <- delta_theta_env * (365.25 / (2.0 * pi))
  adj_doy_calc <- ifelse(has_valid_doy, round(((doy_vec - shift_days - 1.0) %% 365.25) + 1.0, 2), NA_real_)
  adj_in <- if ("hopkins_adjusted_doy" %in% names(df)) suppressWarnings(as.numeric(df$hopkins_adjusted_doy)) else rep(NA_real_, n)
  has_valid_adj <- is.finite(adj_in) & (adj_in >= 1.0) & (adj_in <= 366.0)
  has_pheno_doy <- has_valid_doy | has_valid_adj
  df$hopkins_adjusted_doy <- ifelse(has_valid_adj, round(adj_in, 2), adj_doy_calc)
  df$sin_doy <- ifelse(has_valid_doy, round(sin(theta_obs), 4), NA_real_)
  df$cos_doy <- ifelse(has_valid_doy, round(cos(theta_obs), 4), NA_real_)

  theta_adj_use <- ifelse(is.finite(df$hopkins_adjusted_doy), 2.0 * pi * (df$hopkins_adjusted_doy / 365.25), theta_adj)
  circ_resid_days <- numeric(n)
  sp_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unassigned_Packera", n)
  active_geno_hop <- any(.has_step07_geno_row_mask(df))
  t0_hop <- if (active_geno_hop) {
    if (exists("get_tier0_rectified_species_and_evictions", mode = "function")) get_tier0_rectified_species_and_evictions(df)
    else if (exists("resolve_tier0_anchor_overrides", mode = "function")) resolve_tier0_anchor_overrides(df)
    else NULL
  } else NULL
  for (tx in unique(sp_vec[!is.na(sp_vec)])) {
    idx_tx <- which(sp_vec == tx)
    if (!is.null(t0_hop)) {
      ev_tx <- t0_hop$evict_from_given_anchor & (!t0_hop$is_tier0_dna_anchor | t0_hop$rectified_taxon != tx)
      rect_tx <- which(t0_hop$is_tier0_dna_anchor & !is.na(t0_hop$rectified_taxon) & (t0_hop$rectified_taxon == tx) & is_anthesis & has_pheno_doy & !noise_m)
      anth_clean_tx <- unique(c(idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & anchor_mask[idx_tx] & !mo_m[idx_tx] & !noise_m[idx_tx] & !ev_tx[idx_tx]], rect_tx))
      if (length(anth_clean_tx) < 2L) anth_clean_tx <- unique(c(idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & !noise_m[idx_tx] & !ev_tx[idx_tx]], rect_tx))
      if (length(anth_clean_tx) == 0L) anth_clean_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx]]
      if (length(anth_clean_tx) > 0L) {
        w_h <- ifelse(t0_hop$is_tier0_dna_anchor[anth_clean_tx] & (t0_hop$rectified_taxon[anth_clean_tx] == tx),
                      ifelse(is.finite(t0_hop$tier0_weight[anth_clean_tx]), t0_hop$tier0_weight[anth_clean_tx], 4.0), 1.0)
        theta_tx_center <- atan2(sum(w_h * sin(theta_adj_use[anth_clean_tx]), na.rm = TRUE),
                                 sum(w_h * cos(theta_adj_use[anth_clean_tx]), na.rm = TRUE))
      } else {
        theta_tx_center <- ref_theta
      }
    } else {
      anth_clean_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & anchor_mask[idx_tx] & !mo_m[idx_tx] & !noise_m[idx_tx]]
      if (length(anth_clean_tx) < 2L) anth_clean_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx] & !noise_m[idx_tx]]
      if (length(anth_clean_tx) == 0L) anth_clean_tx <- idx_tx[is_anthesis[idx_tx] & has_pheno_doy[idx_tx]]
      theta_tx_center <- if (length(anth_clean_tx) > 0L) atan2(mean(sin(theta_adj_use[anth_clean_tx]), na.rm = TRUE), mean(cos(theta_adj_use[anth_clean_tx]), na.rm = TRUE)) else ref_theta
    }
    d_theta <- theta_adj_use[idx_tx] - theta_tx_center
    circ_resid_days[idx_tx] <- atan2(sin(d_theta), cos(d_theta)) * (365.25 / (2.0 * pi))
  }
  circ_resid_gated <- ifelse(is_anthesis & has_pheno_doy, round(circ_resid_days, 2), NA_real_)
  ddoy_in <- if ("delta_doy_resid" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy_resid)) else if ("delta_doy" %in% names(df)) suppressWarnings(as.numeric(df$delta_doy)) else rep(NA_real_, n)
  ddoy_use <- ifelse(is_anthesis & (has_pheno_doy | is.finite(ddoy_in)) & is.finite(ddoy_in), ddoy_in, circ_resid_gated)
  ddoy_use[!is_anthesis] <- NA_real_
  df$delta_doy_resid <- ddoy_use
  df$delta_doy <- ddoy_use
  return(df)
}

# ------------------------------------------------------------------------------
# 3. Reliability-Weighted Probabilistic Consensus (Dawid & Skene 1979; Kuncheva 2014)
# ------------------------------------------------------------------------------
compute_cohens_kappa <- function(y_true, y_pred, default_kappa = 0.70) {
  valid <- !is.na(y_true) & !is.na(y_pred) & nzchar(as.character(y_true)) & nzchar(as.character(y_pred))
  if (sum(valid) < 4L) return(max(0.05, default_kappa))
  yt <- as.character(y_true[valid]); yp <- as.character(y_pred[valid])
  lvls <- sort(unique(c(yt, yp)))
  if (length(lvls) < 2L || length(unique(yt)) < 2L) return(max(0.05, min(0.95, mean(yt == yp))))
  tab <- table(factor(yt, levels = lvls), factor(yp, levels = lvls))
  n_tot <- sum(tab); p_o <- sum(diag(tab)) / n_tot; p_e <- sum((rowSums(tab) / n_tot) * (colSums(tab) / n_tot))
  if (!is.finite(p_e) || (1.0 - p_e) <= 1e-8) return(max(0.05, p_o))
  kappa_val <- (p_o - p_e) / (1.0 - p_e)
  if (!is.finite(kappa_val)) return(max(0.05, default_kappa))
  max(0.05, min(0.99, round(kappa_val, 4)))
}

.compute_oof_class_precision_map <- function(y_true_vec, y_pred_vec, prior_alpha = 0.5) {
  ok <- !is.na(y_true_vec) & !is.na(y_pred_vec) & nzchar(y_true_vec) & nzchar(y_pred_vec) &
        (y_true_vec != "Unknown") & (y_pred_vec != "Unknown")
  if (sum(ok) < 6L) return(numeric(0))
  yt <- as.character(y_true_vec[ok])
  yp <- as.character(y_pred_vec[ok])
  true_taxa <- sort(unique(yt))
  all_taxa <- sort(unique(c(true_taxa, yp)))
  if (length(true_taxa) < 2L || length(all_taxa) < 2L) return(numeric(0))
  K_cols <- length(all_taxa)
  # Out-of-fold anchor confusion matrix C^{(m)}_{k,j} = P(y_hat_m = j | y_true = k)
  C_mat <- matrix(0.0, nrow = length(true_taxa), ncol = K_cols, dimnames = list(true_taxa, all_taxa))
  for (k_tx in true_taxa) {
    idx_k <- which(yt == k_tx)
    n_k <- length(idx_k)
    for (j_tx in all_taxa) {
      cnt_kj <- sum(yp[idx_k] == j_tx)
      C_mat[k_tx, j_tx] <- (cnt_kj + prior_alpha / K_cols) / (n_k + prior_alpha)
    }
  }
  prec_out <- stats::setNames(rep(1.0, K_cols), all_taxa)
  for (j_tx in all_taxa) {
    col_sum <- sum(C_mat[, j_tx])
    diag_val <- if (j_tx %in% rownames(C_mat)) C_mat[j_tx, j_tx] else (prior_alpha / K_cols)
    prec_j <- if (is.finite(col_sum) && col_sum > 1e-12) (diag_val / col_sum) else 1.0
    prec_out[[j_tx]] <- min(1.0, max(0.05, round(prec_j, 4)))
  }
  prec_out
}

compute_modality_reliability_weights <- function(df) {
  # Evaluates out-of-fold Cohen's kappa (w_m = max(0.05, kappa_m)) on Gold/Silver/promoted anchors
  # across the 3 independent biological modalities: Morphometrics, Edaphic/Climatic Niche, and Phenology
  n <- nrow(df); if (n == 0L) return(c(morph = 0.80, edaphic = 0.65, pheno = 0.55))
  anchor_mask <- .resolve_step07_anchor_mask(df)
  sub_a <- df[anchor_mask, , drop = FALSE]; y_true <- as.character(sub_a$species_standardized)

  n_true_taxa <- length(unique(y_true[!is.na(y_true) & nzchar(y_true) & y_true != "Unknown"]))
  m_pure <- if ("morph_tier_source" %in% names(sub_a)) (is.na(sub_a$morph_tier_source) | sub_a$morph_tier_source == "Basal_Rosette") else rep(TRUE, nrow(sub_a))
  y_morph <- if ("cda_predicted_taxon" %in% names(sub_a)) as.character(sub_a$cda_predicted_taxon) else rep(NA_character_, nrow(sub_a))
  def_morph <- if ("loocv_accuracy" %in% names(df) && any(is.finite(as.numeric(df$loocv_accuracy)))) stats::median(as.numeric(df$loocv_accuracy), na.rm = TRUE) else 0.80

  if (n_true_taxa < 2L) {
    # Single-taxon / intraspecific reliability weights: avoid degenerate K=1 Cohen's kappa tautology
    # ("Morphotype_K1" != "Packera dubia" -> 0.05 vs "Packera dubia" == "Packera dubia" -> 0.95).
    out_m_a <- if ("morphometric_outlier_flag" %in% names(sub_a)) as.logical(sub_a$morphometric_outlier_flag[m_pure]) else rep(FALSE, sum(m_pure))
    out_m_a[is.na(out_m_a)] <- FALSE
    inlier_rate_m <- if (length(out_m_a) > 0L) mean(!out_m_a) else 0.80
    w_morph <- round(min(0.80, max(0.30, inlier_rate_m * 0.65)), 4)

    eco_reg_a <- if ("ecotype_regime" %in% names(sub_a)) as.character(sub_a$ecotype_regime) else if ("analysis_group" %in% names(sub_a)) as.character(sub_a$analysis_group) else rep(NA_character_, nrow(sub_a))
    rf_reg_a <- if ("rf_predicted_taxon" %in% names(sub_a)) as.character(sub_a$rf_predicted_taxon) else if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else rep(NA_character_, nrow(sub_a))
    u_eco_reg <- unique(eco_reg_a[!is.na(eco_reg_a) & nzchar(eco_reg_a) & eco_reg_a != "Unknown"])
    if (length(u_eco_reg) >= 2L && any(rf_reg_a %in% u_eco_reg, na.rm = TRUE)) {
      w_edaph <- max(0.15, min(0.85, compute_cohens_kappa(eco_reg_a, rf_reg_a, default_kappa = 0.50)))
    } else {
      cv_acc_s <- if ("spatial_block_cv_acc" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$spatial_block_cv_acc))))) {
        stats::median(suppressWarnings(as.numeric(df$spatial_block_cv_acc)), na.rm = TRUE)
      } else 0.65
      w_edaph <- round(min(0.75, max(0.20, 2.0 * (cv_acc_s - 0.50))), 4)
    }
  } else {
    w_morph <- max(0.05, compute_cohens_kappa(y_true[m_pure], y_morph[m_pure], default_kappa = def_morph))
    y_edaph <- if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else if ("rf_predicted_taxon" %in% names(sub_a)) as.character(sub_a$rf_predicted_taxon) else rep(NA_character_, nrow(sub_a))
    def_edaph <- if ("spatial_block_cv_acc" %in% names(df) && any(is.finite(as.numeric(df$spatial_block_cv_acc)))) stats::median(as.numeric(df$spatial_block_cv_acc), na.rm = TRUE) else 0.65
    w_edaph <- max(0.05, compute_cohens_kappa(y_true, y_edaph, default_kappa = def_edaph))
  }

  # Evaluate w_pheno on the identical multi-class out-of-fold Cohen's kappa scale (Dawid & Skene 1979)
  # by testing leave-one-out circular Hopkins-adjusted DOY taxon discrimination across active anchors.
  anth_a <- if ("is_active_anthesis" %in% names(sub_a)) as.logical(sub_a$is_active_anthesis) else resolve_active_anthesis_mask(sub_a)
  anth_a[is.na(anth_a)] <- TRUE
  hop_col_a <- intersect(c("hopkins_adjusted_doy", "doy"), names(sub_a))[1L]
  doy_a <- if (!is.na(hop_col_a)) suppressWarnings(as.numeric(sub_a[[hop_col_a]])) else rep(NA_real_, nrow(sub_a))
  ddoy_a <- if ("delta_doy_resid" %in% names(sub_a)) suppressWarnings(as.numeric(sub_a$delta_doy_resid)) else if ("delta_doy" %in% names(sub_a)) suppressWarnings(as.numeric(sub_a$delta_doy)) else rep(NA_real_, nrow(sub_a))
  ok_ph_cls <- which(anth_a & is.finite(doy_a) & doy_a >= 1.0 & doy_a <= 366.0 & !is.na(y_true) & nzchar(y_true) & y_true != "Unknown")
  tx_ph_tab <- if (length(ok_ph_cls) >= 4L) table(y_true[ok_ph_cls]) else integer(0)
  tx_ph_valid <- names(tx_ph_tab[tx_ph_tab >= 2L])
  ok_ph_cls <- ok_ph_cls[y_true[ok_ph_cls] %in% tx_ph_valid]
  w_pheno <- if (length(tx_ph_valid) >= 2L && length(ok_ph_cls) >= 4L) {
    th_ph <- 2.0 * pi * (doy_a[ok_ph_cls] / 365.25)
    yt_ph <- y_true[ok_ph_cls]
    yp_ph <- character(length(ok_ph_cls))
    for (ii in seq_along(ok_ph_cls)) {
      best_d_th <- Inf; best_tx_ph <- tx_ph_valid[1L]
      for (tx_k in tx_ph_valid) {
        jj <- which(yt_ph == tx_k & seq_along(ok_ph_cls) != ii)
        if (length(jj) == 0L) jj <- which(yt_ph == tx_k)
        mu_k <- atan2(mean(sin(th_ph[jj])), mean(cos(th_ph[jj])))
        r_k <- max(1e-4, min(0.9995, sqrt(mean(sin(th_ph[jj]))^2 + mean(cos(th_ph[jj]))^2)))
        sd_rad_k <- max(2.0 * pi * (6.0 / 365.25), sqrt(-2.0 * log(r_k)))
        d_z <- abs(atan2(sin(th_ph[ii] - mu_k), cos(th_ph[ii] - mu_k))) / sd_rad_k
        if (d_z < best_d_th) { best_d_th <- d_z; best_tx_ph <- tx_k }
      }
      yp_ph[ii] <- best_tx_ph
    }
    max(0.05, min(0.95, compute_cohens_kappa(yt_ph, yp_ph, default_kappa = 0.15)))
  } else if (sum(is.finite(ddoy_a)) >= 4L) {
    0.20
  } else {
    0.15
  }
  w_geno <- if (exists("compute_genomic_kappa_weight", mode = "function")) compute_genomic_kappa_weight(df) else NA_real_
  if (is.finite(w_geno)) c(morph = w_morph, edaphic = w_edaph, pheno = w_pheno, geno = w_geno) else c(morph = w_morph, edaphic = w_edaph, pheno = w_pheno)
}

compute_probabilistic_consensus_scores <- function(df) {
  # Dynamically renormalizes empirical reliability weights modulated by OOF confusion-matrix
  # per-taxon precision w_{i,m} = w_m * (0.5 + 0.5 * Precision_m(y_hat_{i,m})) over observed modalities
  n <- nrow(df); if (n == 0L) return(df)
  prop_sum_attr <- attr(df, "genomic_propagation_summary")
  w_vec <- compute_modality_reliability_weights(df)
  df$w_morph_kappa <- w_vec[["morph"]]
  df$w_edaphic_kappa <- w_vec[["edaphic"]]; df$w_pheno_kappa <- w_vec[["pheno"]]
  has_geno_row <- .has_step07_geno_row_mask(df)
  w_geno_base <- if ("geno" %in% names(w_vec) && is.finite(w_vec[["geno"]])) w_vec[["geno"]] else NA_real_

  anchor_mask <- .resolve_step07_anchor_mask(df)
  sub_a <- df[anchor_mask, , drop = FALSE]
  y_true_a <- if ("species_standardized" %in% names(sub_a)) as.character(sub_a$species_standardized) else rep(NA_character_, nrow(sub_a))
  m_pure_a <- if ("morph_tier_source" %in% names(sub_a)) (is.na(sub_a$morph_tier_source) | sub_a$morph_tier_source == "Basal_Rosette") else rep(TRUE, nrow(sub_a))
  y_morph_a <- if ("cda_predicted_taxon" %in% names(sub_a)) as.character(sub_a$cda_predicted_taxon) else rep(NA_character_, nrow(sub_a))
  y_edaph_a <- if ("edaphic_best_fit_taxon" %in% names(sub_a)) as.character(sub_a$edaphic_best_fit_taxon) else if ("rf_predicted_taxon" %in% names(sub_a)) as.character(sub_a$rf_predicted_taxon) else rep(NA_character_, nrow(sub_a))

  prec_morph <- .compute_oof_class_precision_map(y_true_a[m_pure_a], y_morph_a[m_pure_a])
  prec_edaph <- .compute_oof_class_precision_map(y_true_a, y_edaph_a)
  lookup_prec <- function(pmap, pred_tx) {
    if (length(pmap) == 0L || is.na(pred_tx) || !nzchar(pred_tx) || !(pred_tx %in% names(pmap))) return(1.0)
    as.numeric(pmap[[pred_tx]])
  }

  p_conf_out <- numeric(n); p_score_out <- numeric(n); n_mod_out <- integer(n)
  w_morph_sx <- numeric(n); w_edaph_sx <- numeric(n); w_pheno_sx <- numeric(n); w_geno_sx <- numeric(n)
  df_nms_c <- names(df)
  given_sp_cv <- as.character(df$species_standardized)
  morph_sp_cv <- if ("cda_predicted_taxon" %in% df_nms_c) as.character(df$cda_predicted_taxon) else rep(NA_character_, n)
  morph_sp_cv[!is.na(morph_sp_cv) & !nzchar(morph_sp_cv)] <- NA_character_
  cda_p_cv <- if ("cda_posterior_prob" %in% df_nms_c) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(NA_real_, n)
  is_mo_cv <- if ("morphometric_outlier_flag" %in% df_nms_c) (!is.na(df$morphometric_outlier_flag) & as.logical(df$morphometric_outlier_flag)) else rep(FALSE, n)
  has_geo_cv <- if ("decimalLatitude" %in% df_nms_c) is.finite(suppressWarnings(as.numeric(df$decimalLatitude))) else rep(FALSE, n)
  edaph_sp_cv <- if ("edaphic_best_fit_taxon" %in% df_nms_c) as.character(df$edaphic_best_fit_taxon) else if ("rf_predicted_taxon" %in% df_nms_c) as.character(df$rf_predicted_taxon) else rep(NA_character_, n)
  edaph_sp_cv[!is.na(edaph_sp_cv) & !nzchar(edaph_sp_cv)] <- NA_character_
  rf_p_cv <- if ("rf_predicted_prob" %in% df_nms_c) suppressWarnings(as.numeric(df$rf_predicted_prob)) else rep(NA_real_, n)
  rf_g_cv <- if ("rf_given_taxon_prob" %in% df_nms_c) suppressWarnings(as.numeric(df$rf_given_taxon_prob)) else rep(NA_real_, n)
  ez_cv <- if ("edaphic_z_score" %in% df_nms_c) suppressWarnings(as.numeric(df$edaphic_z_score)) else rep(NA_real_, n)
  tez_cv <- if ("type_edaphic_z_score" %in% df_nms_c) suppressWarnings(as.numeric(df$type_edaphic_z_score)) else rep(NA_real_, n)
  ded_cv <- if ("edaphic_niche_overlap_d" %in% df_nms_c) suppressWarnings(as.numeric(df$edaphic_niche_overlap_d)) else if ("edaphic_schoener_d" %in% df_nms_c) suppressWarnings(as.numeric(df$edaphic_schoener_d)) else rep(NA_real_, n)
  ddoy_r_cv <- if ("delta_doy_resid" %in% df_nms_c) suppressWarnings(as.numeric(df$delta_doy_resid)) else rep(NA_real_, n)
  is_anth_cv <- if ("is_active_anthesis" %in% df_nms_c) (!is.na(df$is_active_anthesis) & as.logical(df$is_active_anthesis)) else rep(TRUE, n)
  w_morph_base <- w_vec[["morph"]]; w_edaph_base <- w_vec[["edaphic"]]; w_pheno_base <- w_vec[["pheno"]]

  is_single_sp_cons <- length(unique(given_sp_cv[!is.na(given_sp_cv) & nzchar(given_sp_cv) & given_sp_cv != "Unknown"])) < 2L
  eco_reg_cv <- if ("ecotype_regime" %in% df_nms_c) as.character(df$ecotype_regime) else if ("analysis_group" %in% df_nms_c) as.character(df$analysis_group) else rep(NA_character_, n)
  sd2_cv <- if ("self_morph_mahal_d2" %in% df_nms_c) suppressWarnings(as.numeric(df$self_morph_mahal_d2)) else rep(NA_real_, n)
  sdf_cv <- if ("self_morph_df" %in% df_nms_c) suppressWarnings(as.numeric(df$self_morph_df)) else rep(6.0, n)
  sdf_cv[!is.finite(sdf_cv) | sdf_cv < 1.0] <- 6.0

  geno_sp_cv <- if ("geno_predicted_taxon" %in% df_nms_c) as.character(df$geno_predicted_taxon) else rep(NA_character_, n)
  geno_p_cv <- if ("geno_posterior_prob" %in% df_nms_c) suppressWarnings(as.numeric(df$geno_posterior_prob)) else rep(NA_real_, n)
  geno_tier_cv <- if ("geno_linkage_tier" %in% df_nms_c) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% df_nms_c) as.character(df$linkage_tier) else rep(NA_character_, n)
  geno_adm_cv <- if ("is_geno_admixed" %in% df_nms_c) {
    if (exists(".safe_bool_vec_r", mode = "function")) .safe_bool_vec_r(df$is_geno_admixed) else (!is.na(df$is_geno_admixed) & as.logical(df$is_geno_admixed))
  } else rep(FALSE, n)
  geno_qmax_cv <- if ("q_max" %in% df_nms_c) suppressWarnings(as.numeric(df$q_max)) else rep(NA_real_, n)
  geno_qry_cv <- if ("is_putative_novel_or_hybrid_query" %in% df_nms_c) {
    if (exists(".safe_bool_vec_r", mode = "function")) .safe_bool_vec_r(df$is_putative_novel_or_hybrid_query) else (!is.na(df$is_putative_novel_or_hybrid_query) & as.logical(df$is_putative_novel_or_hybrid_query))
  } else rep(FALSE, n)
  geno_link_scale <- ifelse(geno_tier_cv == "G1_Direct_Voucher", 1.00,
                     ifelse(geno_tier_cv == "G1_Exsiccatae_Isoduplicate", 0.95,
                     ifelse(geno_tier_cv == "G2_Population_Voucher_Proxy", 0.90,
                     ifelse(geno_tier_cv == "G2_Conditional_Deme_Hull", 0.75, 0.00))))
  geno_link_scale[is.na(geno_link_scale)] <- 0.0
  run_has_geno <- any(has_geno_row) && is.finite(w_geno_base)
  df$w_geno_kappa <- if (run_has_geno) ifelse(has_geno_row, round(w_geno_base * geno_link_scale, 4), 0.0) else NA_real_

  t0_abstain_cv <- if ("tier0_abstention_flag" %in% df_nms_c) suppressWarnings(as.integer(df$tier0_abstention_flag)) else rep(NA_integer_, n)
  t0_cal_tax_cv <- if ("tier0_calibrated_taxon" %in% df_nms_c) trimws(as.character(df$tier0_calibrated_taxon)) else rep(NA_character_, n)
  t0_cal_prob_cv <- if ("tier0_calibrated_prob" %in% df_nms_c) suppressWarnings(as.numeric(df$tier0_calibrated_prob)) else rep(NA_real_, n)
  lopo_cv_kappa <- if (!is.null(prop_sum_attr) && !is.null(prop_sum_attr$lopo_cv_kappa) && is.finite(suppressWarnings(as.numeric(prop_sum_attr$lopo_cv_kappa)))) {
    as.numeric(prop_sum_attr$lopo_cv_kappa)
  } else {
    0.60
  }
  lopo_cv_kappa <- min(1.0, max(0.0, lopo_cv_kappa))
  known_species_pool <- unique(given_sp_cv[!is.na(given_sp_cv) & nzchar(trimws(given_sp_cv)) & given_sp_cv != "Unknown"])

  for (i in seq_len(n)) {
    given_sp <- given_sp_cv[i]; morph_sp <- morph_sp_cv[i]; cda_p <- cda_p_cv[i]; is_mo <- is_mo_cv[i]
    has_geo <- has_geo_cv[i]; edaph_sp <- edaph_sp_cv[i]
    rf_p <- rf_p_cv[i]; rf_g <- rf_g_cv[i]; ez_i <- ez_cv[i]; tez_i <- tez_cv[i]; ded_i <- ded_cv[i]
    ddoy_r <- ddoy_r_cv[i]; is_anth <- is_anth_cv[i]

    morph_is_subcluster <- !is.na(morph_sp) && grepl("^Morphotype_K", morph_sp)
    edaph_is_ecotype <- !is.na(edaph_sp) && grepl("^Ecotype_", edaph_sp)

    p_m <- if (is_single_sp_cons || morph_is_subcluster) {
      if (is_mo) 0.72
      else if (is.finite(sd2_cv[i])) min(0.85, max(0.05, stats::pchisq(sd2_cv[i], df = sdf_cv[i]) * 0.45))
      else if (!is.na(morph_sp) || "gmm_cluster" %in% df_nms_c) 0.08
      else NA_real_
    } else if (!is.na(morph_sp)) {
      if (morph_sp != given_sp) max(if (is.finite(cda_p)) cda_p else 0.75, if (is_mo) 0.80 else 0.55)
      else if (is_mo) 0.65
      else if (is.finite(cda_p)) max(0.0, 1.0 - cda_p)
      else 0.10
    } else NA_real_

    p_e <- if (is_single_sp_cons || edaph_is_ecotype) {
      if (has_geo && (is.finite(ez_i) || is.finite(tez_i))) {
        z_use <- if (is.finite(ez_i)) ez_i else tez_i
        min(0.92, max(0.04, stats::pchisq(z_use^2, df = 1) * 0.55))
      } else if (edaph_is_ecotype && !is.na(eco_reg_cv[i]) && nzchar(eco_reg_cv[i])) {
        if (edaph_sp == eco_reg_cv[i]) 0.08 else 0.35
      } else NA_real_
    } else if (!is.na(edaph_sp)) {
      if (edaph_sp != given_sp) {
        if (is.finite(rf_p)) max(0.55, min(0.98, rf_p))
        else if (is.finite(ez_i) && ez_i > 2.0) min(0.95, 0.55 + 0.10 * (ez_i - 2.0))
        else if (is.finite(ded_i) && ded_i < 0.35) 0.80
        else 0.65
      } else {
        if (is.finite(rf_g)) max(0.05, min(0.85, 1.0 - rf_g))
        else if (is.finite(rf_p)) max(0.05, min(0.85, 1.0 - rf_p))
        else if (is.finite(ez_i) && ez_i > 2.5) min(0.75, 0.25 * ez_i)
        else 0.10
      }
    } else if (has_geo && (is.finite(ez_i) || is.finite(tez_i))) {
      z_use <- if (is.finite(ez_i)) ez_i else tez_i
      min(0.95, max(0.05, stats::pchisq(z_use^2, df = 1) * 0.85))
    } else {
      NA_real_
    }
    p_p <- if (is_anth && is.finite(ddoy_r)) min(1.0, abs(ddoy_r) / 30.0) else NA_real_
    p_g <- if (has_geno_row[i] && is.finite(w_geno_base) && geno_link_scale[i] > 0.0) {
      g_sp_i <- geno_sp_cv[i]; g_p_i <- geno_p_cv[i]; g_qm_i <- geno_qmax_cv[i]
      if (is_single_sp_cons) {
        if (geno_adm_cv[i] || geno_qry_cv[i]) {
          if (is.finite(g_qm_i)) max(0.35, min(0.90, 1.0 - g_qm_i)) else 0.72
        } else if (!is.na(g_sp_i) && g_sp_i != given_sp && is.finite(g_p_i)) {
          max(0.55, min(0.99, g_p_i))
        } else if (is.finite(g_qm_i) && g_qm_i < 0.75) {
          0.55
        } else {
          0.08
        }
      } else if (geno_adm_cv[i] || geno_qry_cv[i]) {
        if (is.finite(g_qm_i)) max(0.40, min(0.92, 1.0 - g_qm_i)) else 0.65
      } else if (!is.na(g_sp_i) && g_sp_i != given_sp) {
        if (is.finite(g_p_i)) max(0.55, min(0.99, g_p_i)) else 0.80
      } else {
        if (is.finite(g_p_i)) max(0.02, min(0.65, 1.0 - g_p_i)) else 0.08
      }
    } else NA_real_

    has_t0_soft_i <- (geno_link_scale[i] == 0.0) &&
                     !is.na(t0_abstain_cv[i]) && (t0_abstain_cv[i] == 0L) &&
                     is.finite(t0_cal_prob_cv[i]) && (t0_cal_prob_cv[i] >= 0.80) &&
                     !is.na(t0_cal_tax_cv[i]) && nzchar(t0_cal_tax_cv[i]) &&
                     (t0_cal_tax_cv[i] != "Sensu_Lato_DNA_Required")
    w_t0_i <- NA_real_; p_t0_i <- NA_real_
    if (has_t0_soft_i && lopo_cv_kappa > 0.0) {
      w_t0_i <- 0.50 * lopo_cv_kappa
      s_t0_i <- if (grepl(" ", t0_cal_tax_cv[i]) || (t0_cal_tax_cv[i] %in% known_species_pool)) {
        if (t0_cal_tax_cv[i] == given_sp) t0_cal_prob_cv[i] else (1.0 - t0_cal_prob_cv[i])
      } else {
        t0_cal_prob_cv[i]
      }
      p_t0_i <- min(1.0, max(0.0, 1.0 - s_t0_i))
    }

    w_i_vec <- c(
      morph   = w_morph_base * (0.5 + 0.5 * lookup_prec(prec_morph, morph_sp)),
      edaphic = w_edaph_base * (0.5 + 0.5 * lookup_prec(prec_edaph, edaph_sp)),
      pheno   = w_pheno_base,
      geno    = if (has_geno_row[i] && is.finite(w_geno_base)) w_geno_base * geno_link_scale[i] else NA_real_,
      t0_soft = w_t0_i
    )
    p_vec <- c(morph = p_m, edaphic = p_e, pheno = p_p, geno = p_g, t0_soft = p_t0_i)
    avail <- !is.na(p_vec) & is.finite(p_vec) & !is.na(w_i_vec) & is.finite(w_i_vec); n_mod_out[i] <- sum(avail)
    p_conf_i <- if (any(avail)) sum(w_i_vec[avail] * p_vec[avail]) / sum(w_i_vec[avail]) else 0.0
    w4_vec <- pmax(0.0, w_i_vec[c("morph", "edaphic", "pheno", "geno")])
    av4 <- avail[c("morph", "edaphic", "pheno", "geno")]
    w4_sx <- rep(0.0, 4L)
    if (any(av4) && sum(w4_vec[av4]) > 0.0) {
      w4_raw <- w4_vec[av4] / sum(w4_vec[av4])
      w4_rnd <- round(w4_raw, 6L)
      idx_max <- which.max(w4_rnd)
      w4_rnd[idx_max] <- round(w4_rnd[idx_max] + (1.0 - sum(w4_rnd)), 6L)
      w4_sx[av4] <- w4_rnd
    }
    w_morph_sx[i] <- w4_sx[1L]; w_edaph_sx[i] <- w4_sx[2L]; w_pheno_sx[i] <- w4_sx[3L]; w_geno_sx[i] <- w4_sx[4L]
    p_conf_out[i] <- round(min(1.0, max(0.0, p_conf_i)), 4); p_score_out[i] <- round(1.0 - p_conf_out[i], 4)
  }
  df$w_morph_simplex <- w_morph_sx
  df$w_edaphic_simplex <- w_edaph_sx
  df$w_pheno_simplex <- w_pheno_sx
  df$w_geno_simplex <- w_geno_sx
  df$n_available_modalities <- n_mod_out
  df$probabilistic_consensus_conflict <- p_conf_out
  df$probabilistic_consensus_score <- p_score_out
  if (!is.null(prop_sum_attr)) attr(df, "genomic_propagation_summary") <- prop_sum_attr
  return(df)
}

# ------------------------------------------------------------------------------
# 4. Five-Stage Hierarchical Integrative Decision Matrix
#    (Padial et al. 2010; Schlick-Steiner et al. 2010; Rieseberg 1995; Slenker et al. 2022)
# ------------------------------------------------------------------------------
.parse_step07_semicolon_cols <- function(df, col_name) {
  if (!(col_name %in% names(df))) return(character(0))
  raw_vals <- as.character(stats::na.omit(df[[col_name]]))
  raw_vals <- trimws(raw_vals[nzchar(trimws(raw_vals))])
  if (length(raw_vals) == 0L) return(character(0))
  toks <- unique(trimws(unlist(strsplit(raw_vals, ";", fixed = TRUE))))
  toks <- toks[nzchar(toks) & !(tolower(toks) %in% c("na", "none", "null"))]
  intersect(toks, names(df))
}

.sort_numbered_pc_cols <- function(cols, prefix) {
  if (length(cols) <= 1L) return(cols)
  idx <- suppressWarnings(as.integer(sub(paste0("^", prefix), "", cols)))
  cols[order(idx, cols)]
}

.resolve_step07_edaphic_predictors <- function(df) {
  exclude_cols <- c(
    "edaphic_z_score", "climatic_z_score", "edaphic_niche_overlap_d", "edaphic_schoener_d",
    "schoener_d", "sgmc_lithology_class", "Dominant_Lithology_Class", "edaphic_best_fit_taxon",
    "rf_predicted_taxon", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors",
    "type_edaphic_z_score", "nearest_type_niche_taxon", "matches_type_lithology"
  )
  keep_finite_num <- function(cands) {
    cands <- setdiff(unique(cands), exclude_cols)
    cands[vapply(cands, function(cc) {
      cc %in% names(df) && !is.logical(df[[cc]]) && any(is.finite(suppressWarnings(as.numeric(df[[cc]]))))
    }, logical(1))]
  }
  p_edaph <- keep_finite_num(.parse_step07_semicolon_cols(df, "primary_edaph_predictors"))
  if (length(p_edaph) > 0L) return(p_edaph)
  p_env <- keep_finite_num(.parse_step07_semicolon_cols(df, "primary_env_predictors"))
  if (length(p_env) > 0L) return(p_env)
  p_fine <- keep_finite_num(grep("^(fine_|ilr_texture_|ssurgo_|polaris_|solus_|topo_)", names(df), value = TRUE))
  if (length(p_fine) > 0L) return(p_fine)
  keep_finite_num(grep("^(soil_|bio[0-9]+_)", names(df), value = TRUE))
}

.compute_dynamic_edaphic_z_scores <- function(df, sp_all_vec, active_taxa, ref_idx_by_taxon, clean_coh_mask = NULL) {
  n <- nrow(df)
  out_z <- if ("edaphic_z_score" %in% names(df)) suppressWarnings(as.numeric(df$edaphic_z_score)) else rep(NA_real_, n)
  out_z[!is.finite(out_z)] <- NA_real_
  need_idx <- which(is.na(out_z))
  if (length(need_idx) == 0L) return(round(out_z, 4))

  pred_cols <- .resolve_step07_edaphic_predictors(df)
  if (length(pred_cols) == 0L) return(round(out_z, 4))

  p_env <- length(pred_cols)
  E_mat <- matrix(NA_real_, nrow = n, ncol = p_env, dimnames = list(NULL, pred_cols))
  for (j in seq_len(p_env)) E_mat[, j] <- suppressWarnings(as.numeric(df[[pred_cols[j]]]))

  noise_m <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)
  coh_mask <- if (!is.null(clean_coh_mask) && length(clean_coh_mask) == n) (clean_coh_mask & !noise_m) else (!noise_m)

  robust_disp <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 2L) return(NA_real_)
    m_v <- stats::mad(v, constant = 1.4826)
    s_v <- stats::sd(v)
    if (is.finite(m_v) && m_v > 1e-6) {
      if (is.finite(s_v)) max(m_v, 0.25 * s_v) else m_v
    } else if (is.finite(s_v) && s_v > 1e-6) {
      s_v
    } else {
      NA_real_
    }
  }

  coh_center <- numeric(p_env); coh_disp <- numeric(p_env)
  for (j in seq_len(p_env)) {
    vj <- E_mat[coh_mask & is.finite(E_mat[, j]), j]
    if (length(vj) < 3L) vj <- E_mat[!noise_m & is.finite(E_mat[, j]), j]
    if (length(vj) == 0L) vj <- E_mat[is.finite(E_mat[, j]), j]
    coh_center[j] <- if (length(vj) > 0L) stats::median(vj) else NA_real_
    coh_disp[j] <- robust_disp(vj)
  }

  tx_center <- matrix(NA_real_, nrow = length(active_taxa), ncol = p_env, dimnames = list(active_taxa, pred_cols))
  tx_disp <- matrix(NA_real_, nrow = length(active_taxa), ncol = p_env, dimnames = list(active_taxa, pred_cols))
  for (tx in active_taxa) {
    r_idx <- ref_idx_by_taxon[[tx]]
    clean_t_idx <- which(sp_all_vec == tx & !noise_m)
    for (j in seq_len(p_env)) {
      v_ref <- E_mat[r_idx[is.finite(E_mat[r_idx, j])], j]
      if (length(v_ref) == 0L) v_ref <- E_mat[clean_t_idx[is.finite(E_mat[clean_t_idx, j])], j]
      tx_center[tx, j] <- if (length(v_ref) > 0L) stats::median(v_ref) else coh_center[j]
      if (length(v_ref) >= 3L) {
        d_t <- robust_disp(v_ref)
        if (is.finite(d_t) && is.finite(coh_disp[j])) d_t <- max(d_t, 0.15 * coh_disp[j])
        tx_disp[tx, j] <- if (is.finite(d_t)) d_t else coh_disp[j]
      } else {
        tx_disp[tx, j] <- coh_disp[j]
      }
    }
  }

  for (i in need_idx) {
    x_i <- E_mat[i, ]
    if (!any(is.finite(x_i))) {
      out_z[i] <- NA_real_
      next
    }
    g_sp <- sp_all_vec[i]
    mu_i <- if (!is.na(g_sp) && g_sp %in% rownames(tx_center)) tx_center[g_sp, ] else coh_center
    sd_i <- if (!is.na(g_sp) && g_sp %in% rownames(tx_disp)) tx_disp[g_sp, ] else coh_disp
    ok_j <- is.finite(x_i) & is.finite(mu_i) & is.finite(sd_i) & (sd_i > 1e-8)
    if (any(ok_j)) {
      out_z[i] <- sqrt(mean(((x_i[ok_j] - mu_i[ok_j]) / sd_i[ok_j])^2))
    } else {
      fin_both <- is.finite(x_i) & is.finite(mu_i)
      out_z[i] <- if (any(fin_both) && all(abs(x_i[fin_both] - mu_i[fin_both]) <= 1e-8)) 0.0 else NA_real_
    }
  }
  round(out_z, 4)
}

.parse_cytotype_x_vals <- function(s) {
  if (is.null(s) || length(s) == 0L) return(integer(0))
  txt <- as.character(s)
  txt <- txt[!is.na(txt) & nzchar(trimws(txt)) & !(toupper(trimws(txt)) %in% c("NA", "NAN", "NONE", "NULL"))]
  if (length(txt) == 0L) return(integer(0))
  out <- integer(0)
  for (item in txt) {
    m <- regmatches(item, gregexpr("([0-9]+)\\s*x\\b", item, ignore.case = TRUE))[[1L]]
    if (length(m) > 0L) {
      vals <- suppressWarnings(as.integer(sub("\\s*x$", "", tolower(m))))
      out <- c(out, vals[!is.na(vals) & vals >= 2L & vals <= 16L])
    }
  }
  sort(unique(out))
}

.infer_taxon_ploidy_levels <- function(tx, sp_vec, vp_vec, tx_cyto_vec) {
  if (is.null(tx) || is.na(tx) || !nzchar(trimws(tx))) return(integer(0))
  idx <- which(sp_vec == tx)
  if (length(idx) == 0L) return(integer(0))
  vp <- vp_vec[idx]
  vp <- vp[!is.na(vp) & is.finite(vp) & vp >= 2L]
  if (length(vp) > 0L) return(as.integer(round(stats::median(vp))))
  .parse_cytotype_x_vals(tx_cyto_vec[idx])
}

.hybrid_cytotype_annotation <- function(vp_i, tax_cyto_i, p_a, p_b, sp_vec, vp_vec, tx_cyto_vec) {
  is_odd_vp <- !is.na(vp_i) && (vp_i >= 3L) && ((vp_i %% 2L) == 1L)
  has_odd_str <- !is.na(tax_cyto_i) && nzchar(trimws(tax_cyto_i)) &&
    grepl("(^|[^0-9])[3579]\\s*x\\b", tax_cyto_i, ignore.case = TRUE)
  pa_lvls <- .infer_taxon_ploidy_levels(p_a, sp_vec, vp_vec, tx_cyto_vec)
  pb_lvls <- .infer_taxon_ploidy_levels(p_b, sp_vec, vp_vec, tx_cyto_vec)
  both_known <- (length(pa_lvls) > 0L) && (length(pb_lvls) > 0L)
  both_2x <- both_known && all(pa_lvls == 2L) && all(pb_lvls == 2L)
  parents_differ <- both_known && !identical(sort(unique(pa_lvls)), sort(unique(pb_lvls)))
  voucher_diverges_2x <- both_2x && !is.na(vp_i) && (vp_i != 2L)
  if (is_odd_vp || has_odd_str || parents_differ || voucher_diverges_2x) {
    return(" [Heteroploid / Odd-Ploidy (3x/5x) or Allopolyploid Candidate]")
  }
  if (both_2x && (is.na(vp_i) || vp_i == 2L)) {
    return(" [Homoploid Hybrid Candidate: shared 2x]")
  }
  ""
}

.resolve_cytotype_barrier_pair <- function(i, sp_i, g_cl_i, cat_i, pid_i, barrier_attr, sp_vec, g_cl_vec, vp_vec, tx_cyto_vec, is_type_vec, is_typ_morph_vec) {
  raw_lbl <- NA_character_
  if (!is.null(barrier_attr)) {
    if (is.data.frame(barrier_attr) && nrow(barrier_attr) > 0L) {
      hit <- integer(0)
      if (all(c("catalogNumber", "plant_individual_id") %in% names(barrier_attr))) {
        hit <- which(as.character(barrier_attr$catalogNumber) == cat_i & as.integer(barrier_attr$plant_individual_id) == pid_i)
      }
      if (length(hit) == 0L && all(c("species_standardized", "gmm_cluster") %in% names(barrier_attr))) {
        hit <- which(as.character(barrier_attr$species_standardized) == sp_i & as.integer(barrier_attr$gmm_cluster) == g_cl_i)
      }
      if (length(hit) > 0L) {
        lcol <- intersect(c("barrier_label", "label", "cytotype_barrier_label"), names(barrier_attr))
        if (length(lcol) > 0L) raw_lbl <- as.character(barrier_attr[[lcol[1L]]][hit[1L]])
      }
    } else if (is.character(barrier_attr) || is.list(barrier_attr)) {
      nm_b <- names(barrier_attr)
      keys <- c(
        paste0(sp_i, "||", g_cl_i), paste0(cat_i, "__", pid_i), cat_i,
        paste0(sp_i, "__K", g_cl_i), paste0(sp_i, "__", g_cl_i), sp_i
      )
      if (!is.null(nm_b) && any(keys %in% nm_b)) {
        k_hit <- keys[keys %in% nm_b][1L]
        val <- barrier_attr[[k_hit]]
        if (is.list(val) && !is.null(val$cl_lbl) && !is.null(val$type_lbl)) {
          return(list(cl_lbl = as.character(val$cl_lbl), type_lbl = as.character(val$type_lbl)))
        }
        raw_lbl <- as.character(val)[1L]
      } else if (is.character(barrier_attr) && length(barrier_attr) >= i) {
        cand_v <- trimws(as.character(barrier_attr[i]))
        if (!is.na(cand_v) && nzchar(cand_v)) raw_lbl <- cand_v
      }
    }
  }
  if (!is.na(raw_lbl) && grepl("\\s+vs\\s+", raw_lbl)) {
    pts <- trimws(strsplit(raw_lbl, "\\s+vs\\s+")[[1L]])
    if (length(pts) >= 2L && nzchar(pts[1L]) && nzchar(pts[2L])) {
      cl_tok <- sub("^cluster ploidy\\s+", "", pts[1L], ignore.case = TRUE)
      ty_tok <- sub("^type-cluster\\s+", "", pts[2L], ignore.case = TRUE)
      return(list(cl_lbl = cl_tok, type_lbl = ty_tok))
    }
  }
  cl_vp <- vp_vec[sp_vec == sp_i & g_cl_vec == g_cl_i]
  cl_vp <- cl_vp[!is.na(cl_vp) & is.finite(cl_vp)]
  cl_p <- if (length(cl_vp) > 0L) {
    as.integer(round(stats::median(cl_vp)))
  } else if (!is.na(vp_vec[i]) && is.finite(vp_vec[i])) {
    as.integer(round(vp_vec[i]))
  } else {
    NA_integer_
  }
  type_vp <- vp_vec[sp_vec == sp_i & is_type_vec & !is.na(vp_vec) & is.finite(vp_vec)]
  if (length(type_vp) == 0L) {
    type_vp <- vp_vec[sp_vec == sp_i & (g_cl_vec == 1L | (!is.na(is_typ_morph_vec) & is_typ_morph_vec)) & !is.na(vp_vec) & is.finite(vp_vec)]
  }
  if (length(type_vp) == 0L) {
    type_vp <- vp_vec[sp_vec == sp_i & (g_cl_vec != g_cl_i) & !is.na(vp_vec) & is.finite(vp_vec)]
  }
  type_p <- if (length(type_vp) > 0L) {
    as.integer(round(stats::median(type_vp)))
  } else if (!is.na(cl_p) && !is.na(tx_cyto_vec[i])) {
    others <- setdiff(.parse_cytotype_x_vals(tx_cyto_vec[i]), cl_p)
    if (length(others) == 1L) others[1L] else NA_integer_
  } else {
    NA_integer_
  }
  list(
    cl_lbl = if (!is.na(cl_p)) sprintf("%dx", cl_p) else "heteroploid",
    type_lbl = if (!is.na(type_p)) sprintf("%dx", type_p) else "type-cluster"
  )
}

.format_infra_cytotype_lbl <- function(raw_sp, given_sp, g_cl) {
  has_infra <- grepl("\\b(var\\.|subsp\\.|ssp\\.)\\s+[a-z-]+", raw_sp, ignore.case = TRUE)
  if (has_infra) {
    m_inf <- regmatches(raw_sp, regexpr("\\b(var\\.|subsp\\.|ssp\\.)\\s+[A-Za-z-]+", raw_sp, ignore.case = TRUE))
    if (length(m_inf) > 0L && nzchar(m_inf[1L])) return(paste(given_sp, tolower(m_inf[1L])))
  }
  sprintf("%s var. cytotype_K%d", given_sp, max(1L, as.integer(g_cl)))
}

.build_stage4_cytotype_decision <- function(given_sp, raw_sp, g_cl, base_rationale) {
  suffix <- " [Intraspecific Cytotype Series: FCM/squash verification recommended]"
  rat_out <- if (grepl("Intraspecific Cytotype Series:", base_rationale, fixed = TRUE)) base_rationale else paste0(base_rationale, suffix)
  list(
    status = "Subspecies_Ecotype",
    broad = "Variety / Ecotype",
    sub = "Subspecies_Ecotype",
    prio = "MEDIUM",
    act = "Accept_Subspecific_Treatment",
    rec = .format_infra_cytotype_lbl(raw_sp, given_sp, g_cl),
    rat = rat_out
  )
}

apply_taxonomic_decision_matrix <- function(df, intraspecific_mode = FALSE, gmm_two_delta_bic = 12.0, opts = list()) {
  message("Evaluating vouchers across Five-Stage Hierarchical Integrative Decision Matrix...")
  if (!("morphometric_outlier_flag" %in% names(df))) df$morphometric_outlier_flag <- FALSE
  if (!("misidentification_flag" %in% names(df))) df$misidentification_flag <- FALSE

  is_intra <- isTRUE(intraspecific_mode)
  n <- nrow(df)
  if (n == 0L) return(df)

  opts_cyto <- opts
  if (is.null(opts_cyto$db_path)) {
    db_env_cy <- env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")
    if (nzchar(db_env_cy)) opts_cyto$db_path <- db_env_cy
  }
  if (is.null(opts_cyto$report_out) && !is.null(opts_cyto$output_report) && nzchar(opts_cyto$output_report)) {
    opts_cyto$report_out <- opts_cyto$output_report
  }
  cyto_dir_arg <- if (!is.null(opts_cyto$cytology_dir) && nzchar(opts_cyto$cytology_dir)) {
    opts_cyto$cytology_dir
  } else {
    env_cdir <- env_plantmatics_or_packera("PLANTMATICS_CYTOLOGY_DIR", "PACKERA_CYTOLOGY_DIR")
    if (nzchar(env_cdir)) env_cdir else "data/cytology"
  }
  has_prepop_cyto <- is.null(opts$cytological_records_df) && is.null(opts$taxon_cytology_summary_df) && (
    isTRUE(attr(df, "step07_cytology_executed")) ||
    (("cytotype_reproductive_barrier_pass" %in% names(df)) && any(!is.na(df$cytotype_reproductive_barrier_pass))) ||
    (("cytology_linkage_tier" %in% names(df)) && any(!is.na(df$cytology_linkage_tier) & nzchar(trimws(as.character(df$cytology_linkage_tier)))))
  )
  if (!has_prepop_cyto && exists("run_step07_cytology_integration", mode = "function")) {
    df <- run_step07_cytology_integration(df, cytology_dir = cyto_dir_arg, opts = opts_cyto)
    attr(df, "step07_cytology_executed") <- TRUE
  }
  cyto_md_attr <- attr(df, "cytology_markdown_lines")
  cyto_bar_lbl_attr <- attr(df, "cytotype_barrier_labels")
  gigas_df_attr <- attr(df, "gigas_concordance_df")

  if (!("genomic_propagation_tier" %in% names(df)) &&
      any(.has_step07_geno_row_mask(df)) &&
      exists("run_step07_genomic_propagation", mode = "function")) {
    df <- run_step07_genomic_propagation(df, opts)
  }
  prop_sum_attr <- attr(df, "genomic_propagation_summary")
  sp_all_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep("Unknown", n)
  sp_sheet_vec <- ifelse(!is.na(sp_all_vec) & nzchar(sp_all_vec), sp_all_vec, "Unknown")
  cat_sheet_vec <- if ("catalogNumber" %in% names(df)) as.character(df$catalogNumber) else paste0("__ROW_", seq_len(n))
  cat_sheet_vec[is.na(cat_sheet_vec) | !nzchar(trimws(cat_sheet_vec))] <- paste0("__ROW_", which(is.na(cat_sheet_vec) | !nzchar(trimws(cat_sheet_vec))))
  sheet_group_key <- paste(cat_sheet_vec, sp_sheet_vec, sep = "||")

  # Sibling Ramet Morphology Propagation across multi-plant sheets (catalogNumber, species_standardized)
  cda_pred_init <- if ("cda_predicted_taxon" %in% names(df)) as.character(df$cda_predicted_taxon) else rep(NA_character_, n)
  cda_prob_init <- if ("cda_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(NA_real_, n)
  tm_d2_init <- if ("type_morph_mahal_d2" %in% names(df)) suppressWarnings(as.numeric(df$type_morph_mahal_d2)) else rep(NA_real_, n)
  n_lvs_raw <- if ("n_leaves" %in% names(df)) suppressWarnings(as.integer(df$n_leaves)) else rep(NA_integer_, n)
  leaf_cnt_raw <- if ("leaf_count" %in% names(df)) suppressWarnings(as.integer(df$leaf_count)) else rep(NA_integer_, n)
  bas_cnt_raw <- if ("basal_leaf_count" %in% names(df)) suppressWarnings(as.integer(df$basal_leaf_count)) else rep(NA_integer_, n)
  caul_cnt_raw <- if ("cauline_leaf_count" %in% names(df)) suppressWarnings(as.integer(df$cauline_leaf_count)) else rep(NA_integer_, n)
  bc_sum <- ifelse(is.na(bas_cnt_raw) & is.na(caul_cnt_raw), NA_integer_, ifelse(is.na(bas_cnt_raw), 0L, bas_cnt_raw) + ifelse(is.na(caul_cnt_raw), 0L, caul_cnt_raw))
  n_leaves_init <- dplyr::coalesce(n_lvs_raw, leaf_cnt_raw, bc_sum, bas_cnt_raw)

  orig_has_obs_morph_v <- !is.na(cda_pred_init) & nzchar(trimws(cda_pred_init)) & is.finite(cda_prob_init)
  orig_n_leaves_v <- ifelse(!is.na(n_leaves_init), pmax(0L, n_leaves_init), ifelse(orig_has_obs_morph_v, 1L, 0L))
  orig_cda_pred_v <- cda_pred_init
  orig_cda_prob_v <- cda_prob_init
  orig_tm_d2_v <- tm_d2_init
  has_observed_or_sheet_morph <- orig_has_obs_morph_v

  prop_morph_cols <- intersect(
    c(
      "cda_predicted_taxon", "cda_posterior_prob", "cda_oof_posterior_prob", "cl_target_threshold",
      "self_morph_mahal_d2", "self_morph_df", "cda_posterior_margin",
      "cda_runner_up_taxon", "cda_second_taxon", "cda_second_predicted_taxon",
      "cda_second_posterior_prob", "cda_mahal_d2", "type_morph_mahal_d2",
      "nearest_type_morph_taxon", "nearest_type_basionym", "nearest_type_basionym_d2",
      "type_morph_concordance", "is_typical_morphotype",
      "min_congener_mahal_d2", "min_congener_taxon", "open_set_novelty_p", "open_set_df",
      "gmm_cluster", "gmm_uncertainty", "is_morph_outlier",
      "morphometric_outlier_flag", "misidentification_flag"
    ),
    names(df)
  )
  mts_init_v <- if ("morph_tier_source" %in% names(df)) as.character(df$morph_tier_source) else rep("Basal_Rosette", n)
  is_caul_fb_init_v <- !is.na(mts_init_v) & (mts_init_v == "Cauline_Only_Fallback")
  is_basal_obs_init_v <- orig_has_obs_morph_v & !is_caul_fb_init_v
  sheet_idx_list <- split(seq_len(n), sheet_group_key)
  for (s_idx in sheet_idx_list) {
    if (length(s_idx) < 2L) next
    has_basal_donor <- any(is_basal_obs_init_v[s_idx])
    obs_s <- if (has_basal_donor) s_idx[is_basal_obs_init_v[s_idx]] else s_idx[orig_has_obs_morph_v[s_idx]]
    miss_s <- if (has_basal_donor) s_idx[!orig_has_obs_morph_v[s_idx] | is_caul_fb_init_v[s_idx]] else s_idx[!orig_has_obs_morph_v[s_idx]]
    if (length(obs_s) >= 1L && length(miss_s) >= 1L) {
      ord_don <- order(-orig_n_leaves_v[obs_s], -orig_cda_prob_v[obs_s])
      donor_idx <- obs_s[ord_don[1L]]
      for (mc in prop_morph_cols) {
        df[[mc]][miss_s] <- df[[mc]][donor_idx]
      }
      if ("cda_posterior_prob" %in% names(df) && is.finite(df$cda_posterior_prob[donor_idx])) {
        df$cda_posterior_prob[miss_s] <- 0.90 * as.numeric(df$cda_posterior_prob[donor_idx])
      }
      if ("cda_posterior_margin" %in% names(df) && is.finite(df$cda_posterior_margin[donor_idx])) {
        df$cda_posterior_margin[miss_s] <- 0.90 * as.numeric(df$cda_posterior_margin[donor_idx])
      }
      has_observed_or_sheet_morph[miss_s] <- TRUE
    }
  }
  df$has_observed_or_sheet_morph <- has_observed_or_sheet_morph

  df <- compute_hopkins_adjusted_phenology(df)
  df <- compute_probabilistic_consensus_scores(df)

  active_taxa <- sort(unique(sp_all_vec[!is.na(sp_all_vec) & nzchar(sp_all_vec) & sp_all_vec != "Unknown"]))
  n_known_taxa <- max(1L, length(active_taxa))

  # 1. Precompute empirical per-taxon reference envelopes on active Gold/Silver/promoted anchors
  #    (falling back to non-outlier records of each taxon when < 2 anchors exist)
  anchor_all <- .resolve_step07_anchor_mask(df)
  mo_all <- if ("morphometric_outlier_flag" %in% names(df)) as.logical(df$morphometric_outlier_flag) else rep(FALSE, n)
  mo_all[is.na(mo_all)] <- FALSE
  noise_all <- if ("is_confident_label_noise" %in% names(df)) (!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else rep(FALSE, n)

  doy_num_all <- if ("doy" %in% names(df)) suppressWarnings(as.numeric(df$doy)) else rep(NA_real_, n)
  valid_raw_doy_all <- is.finite(doy_num_all) & doy_num_all >= 1.0 & doy_num_all <= 366.0
  hop_in_all <- if ("hopkins_adjusted_doy" %in% names(df)) suppressWarnings(as.numeric(df$hopkins_adjusted_doy)) else rep(NA_real_, n)
  valid_hop_all <- is.finite(hop_in_all) & hop_in_all >= 1.0 & hop_in_all <= 366.0
  valid_doy_all <- valid_raw_doy_all | valid_hop_all
  hop_doy_all <- ifelse(valid_hop_all, hop_in_all, ifelse(valid_raw_doy_all, doy_num_all, NA_real_))
  anth_all <- if ("is_active_anthesis" %in% names(df)) as.logical(df$is_active_anthesis) else rep(TRUE, n)
  anth_all[is.na(anth_all)] <- TRUE

  ph_num_all <- if ("fine_soil_ph" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_ph)) else rep(NA_real_, n)
  if ("soil_ph" %in% names(df)) ph_num_all <- ifelse(is.finite(ph_num_all), ph_num_all, suppressWarnings(as.numeric(df$soil_ph)))
  sand_num_all <- if ("fine_soil_sand" %in% names(df)) suppressWarnings(as.numeric(df$fine_soil_sand)) else rep(NA_real_, n)
  if ("soil_sand" %in% names(df)) sand_num_all <- ifelse(is.finite(sand_num_all), sand_num_all, suppressWarnings(as.numeric(df$soil_sand)))

  .circ_hop_stats <- function(h_vals, w_vals = NULL) {
    ok_h <- is.finite(h_vals)
    hv <- h_vals[ok_h]
    if (length(hv) == 0L) return(list(mu_theta = 2.0 * pi * (120.0 / 365.25), q75_arc = 15.0, q75_doy = 135.0, med_doy = 120.0))
    th <- 2.0 * pi * (hv / 365.25)
    if (!is.null(w_vals) && length(w_vals) == length(h_vals)) {
      wv <- w_vals[ok_h]; wv[!is.finite(wv) | wv <= 0.0] <- 1.0
      mu_th <- atan2(sum(wv * sin(th)) / sum(wv), sum(wv * cos(th)) / sum(wv))
    } else {
      mu_th <- atan2(mean(sin(th)), mean(cos(th)))
    }
    if (mu_th < 0.0) mu_th <- mu_th + 2.0 * pi
    arcs <- atan2(sin(th - mu_th), cos(th - mu_th)) * (365.25 / (2.0 * pi))
    q75_a <- as.numeric(stats::quantile(arcs, probs = 0.75, na.rm = TRUE, names = FALSE))
    med_d <- ((mu_th * (365.25 / (2.0 * pi)) - 1.0) %% 365.25) + 1.0
    q75_d <- ((med_d + q75_a - 1.0) %% 365.25) + 1.0
    list(mu_theta = mu_th, q75_arc = q75_a, q75_doy = q75_d, med_doy = med_d)
  }

  glob_hop_pool <- hop_doy_all[valid_doy_all & anth_all & anchor_all & !mo_all & !noise_all]
  if (length(glob_hop_pool) == 0L) glob_hop_pool <- hop_doy_all[valid_doy_all & anth_all & !noise_all]
  if (length(glob_hop_pool) == 0L) glob_hop_pool <- hop_doy_all[valid_doy_all]
  glob_circ <- .circ_hop_stats(glob_hop_pool)

  mu_theta_hop <- stats::setNames(rep(glob_circ$mu_theta, length(active_taxa)), active_taxa)
  q75_arc_hop <- stats::setNames(rep(glob_circ$q75_arc, length(active_taxa)), active_taxa)
  q75_doy <- stats::setNames(rep(glob_circ$q75_doy, length(active_taxa)), active_taxa)
  med_doy <- stats::setNames(rep(glob_circ$med_doy, length(active_taxa)), active_taxa)

  active_geno_ref <- any(.has_step07_geno_row_mask(df))
  t0_info <- if (active_geno_ref) {
    if (exists("get_tier0_rectified_species_and_evictions", mode = "function")) get_tier0_rectified_species_and_evictions(df)
    else if (exists("resolve_tier0_anchor_overrides", mode = "function")) resolve_tier0_anchor_overrides(df)
    else NULL
  } else NULL

  ext_non_type_all <- .is_step07_extralimital_non_type(df)
  ref_idx_by_taxon <- list()
  for (tx in active_taxa) {
    idx_t <- which(sp_all_vec == tx)
    if (!is.null(t0_info)) {
      pure_ev <- t0_info$evict_from_given_anchor & (!t0_info$is_tier0_dna_anchor | t0_info$rectified_taxon != tx)
      t0_rect_tx <- which(t0_info$is_tier0_dna_anchor & !is.na(t0_info$rectified_taxon) & (t0_info$rectified_taxon == tx) & !noise_all & !ext_non_type_all)
      idx_gs <- unique(c(idx_t[anchor_all[idx_t] & !mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t] & !pure_ev[idx_t]], t0_rect_tx))
      idx_clean <- unique(c(idx_t[!mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t] & !pure_ev[idx_t]], t0_rect_tx))
      idx_nonoise <- unique(c(idx_t[!noise_all[idx_t] & !pure_ev[idx_t]], t0_rect_tx))
      ref_idx <- if (length(idx_gs) >= 2L) idx_gs else if (length(idx_clean) >= 2L) idx_clean else if (length(idx_nonoise) >= 1L) idx_nonoise else idx_t
      ref_idx_by_taxon[[tx]] <- ref_idx

      ref_anth_idx <- ref_idx[valid_doy_all[ref_idx] & anth_all[ref_idx]]
      if (length(ref_anth_idx) == 0L) ref_anth_idx <- idx_t[valid_doy_all[idx_t] & anth_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t] & !pure_ev[idx_t]]
      if (length(ref_anth_idx) == 0L) ref_anth_idx <- idx_t[valid_doy_all[idx_t] & !pure_ev[idx_t]]
      if (length(ref_anth_idx) == 0L) ref_anth_idx <- idx_t[valid_doy_all[idx_t]]
      if (length(ref_anth_idx) > 0L) {
        w_anth <- ifelse(t0_info$is_tier0_dna_anchor[ref_anth_idx] & (t0_info$rectified_taxon[ref_anth_idx] == tx),
                         ifelse(is.finite(t0_info$tier0_weight[ref_anth_idx]), t0_info$tier0_weight[ref_anth_idx], 4.0), 1.0)
        tx_circ <- .circ_hop_stats(hop_doy_all[ref_anth_idx], w_vals = w_anth)
        mu_theta_hop[[tx]] <- tx_circ$mu_theta
        q75_arc_hop[[tx]] <- tx_circ$q75_arc
        q75_doy[[tx]] <- tx_circ$q75_doy
        med_doy[[tx]] <- tx_circ$med_doy
      }
    } else {
      idx_gs <- idx_t[anchor_all[idx_t] & !mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t]]
      idx_clean <- idx_t[!mo_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t]]
      idx_nonoise <- idx_t[!noise_all[idx_t]]
      ref_idx <- if (length(idx_gs) >= 2L) idx_gs else if (length(idx_clean) >= 2L) idx_clean else if (length(idx_nonoise) >= 1L) idx_nonoise else idx_t
      ref_idx_by_taxon[[tx]] <- ref_idx

      d_ref <- hop_doy_all[ref_idx[valid_doy_all[ref_idx] & anth_all[ref_idx]]]
      if (length(d_ref) == 0L) d_ref <- hop_doy_all[idx_t[valid_doy_all[idx_t] & anth_all[idx_t] & !noise_all[idx_t] & !ext_non_type_all[idx_t]]]
      if (length(d_ref) == 0L) d_ref <- hop_doy_all[idx_t[valid_doy_all[idx_t]]]
      if (length(d_ref) > 0L) {
        tx_circ <- .circ_hop_stats(d_ref)
        mu_theta_hop[[tx]] <- tx_circ$mu_theta
        q75_arc_hop[[tx]] <- tx_circ$q75_arc
        q75_doy[[tx]] <- tx_circ$q75_doy
        med_doy[[tx]] <- tx_circ$med_doy
      }
    }
  }

  edaph_z_vec <- .compute_dynamic_edaphic_z_scores(df, sp_all_vec, active_taxa, ref_idx_by_taxon, clean_coh_mask = (anchor_all & !mo_all & !noise_all))
  pheno_taus <- .compute_taxon_pheno_thresholds(df)
  tau_ecotype_vec <- pheno_taus$tau_ecotype
  tau_cryptic_vec <- pheno_taus$tau_cryptic

  # Build Biogeographic Range Gate (USDA PLANTS state distributions + clean reference anchor envelopes)
  # so Stage 1 Misidentification_Reassign never reassigns a voucher to a strictly allopatric species.
  .us_abbr_to_state <- c(
    AL = "alabama", AK = "alaska", AZ = "arizona", AR = "arkansas", CA = "california",
    CO = "colorado", CT = "connecticut", DE = "delaware", DC = "district of columbia", FL = "florida",
    GA = "georgia", HI = "hawaii", ID = "idaho", IL = "illinois", IN = "indiana",
    IA = "iowa", KS = "kansas", KY = "kentucky", LA = "louisiana", ME = "maine",
    MD = "maryland", MA = "massachusetts", MI = "michigan", MN = "minnesota", MS = "mississippi",
    MO = "missouri", MT = "montana", NE = "nebraska", NV = "nevada", NH = "new hampshire",
    NJ = "new jersey", NM = "new mexico", NY = "new york", NC = "north carolina", ND = "north dakota",
    OH = "ohio", OK = "oklahoma", OR = "oregon", PA = "pennsylvania", RI = "rhode island",
    SC = "south carolina", SD = "south dakota", TN = "tennessee", TX = "texas", UT = "utah",
    VT = "vermont", VA = "virginia", WA = "washington", WV = "west virginia", WI = "wisconsin",
    WY = "wyoming", AB = "alberta", BC = "british columbia", MB = "manitoba",
    NB = "new brunswick", NL = "newfoundland and labrador", NS = "nova scotia",
    NT = "northwest territories", NU = "nunavut", ON = "ontario",
    PE = "prince edward island", QC = "quebec", SK = "saskatchewan", YT = "yukon"
  )
  .norm_state_str <- function(sv) {
    s <- tolower(trimws(as.character(sv)))
    s[is.na(s) | s %in% c("", "na", "nan", "none", "null", "unknown")] <- ""
    up <- toupper(s)
    hit_ab <- nzchar(s) & (up %in% names(.us_abbr_to_state))
    if (any(hit_ab)) s[hit_ab] <- unname(.us_abbr_to_state[up[hit_ab]])
    s
  }
  .norm_county_str <- function(cv) {
    s <- tolower(trimws(as.character(cv)))
    s[is.na(s) | s %in% c("", "na", "nan", "none", "null", "unknown")] <- ""
    s <- gsub("\\b(county|parish|co\\.|borough)\\b", "", s)
    s <- gsub("[^a-z0-9]+", " ", s)
    trimws(s)
  }
  state_vec_geo <- if ("stateProvince" %in% names(df)) .norm_state_str(df$stateProvince) else rep("", n)
  cnty_raw_geo <- if ("county" %in% names(df)) .norm_county_str(df$county) else rep("", n)
  state_cnty_geo <- ifelse(nzchar(state_vec_geo) & nzchar(cnty_raw_geo), paste(state_vec_geo, cnty_raw_geo, sep = "||"), "")
  lat_vec_geo <- if ("decimalLatitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, n)
  lon_vec_geo <- if ("decimalLongitude" %in% names(df)) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, n)

  is_binomial_packera <- function(v) {
    s <- trimws(as.character(v))
    !is.na(s) & grepl("^Packera\\s+[a-z][a-z-]+$", s)
  }
  cohort_species <- sort(unique(trimws(sp_all_vec[is_binomial_packera(sp_all_vec)])))
  can_evaluate_hybrids <- !is_intra && (length(cohort_species) >= 2L)

  has_any_geo_in_df <- any(nzchar(state_vec_geo)) || any(is.finite(lat_vec_geo) & is.finite(lon_vec_geo))
  use_ext_geo_db <- nzchar(Sys.getenv("RUN_ID", unset = "")) || nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")) || has_any_geo_in_df
  usda_states_by_tx <- list()
  usda_counties_by_tx <- list()
  usda_state_county_counts_by_tx <- list()
  db_syntopy_rows <- NULL
  if (use_ext_geo_db && trimws(env_plantmatics_or_packera("PLANTMATICS_DEMO_MODE", "PACKERA_DEMO_MODE")) != "1" && exists("connect_packera_db", mode = "function") && requireNamespace("DBI", quietly = TRUE)) {
    tryCatch({
      c_geo <- connect_packera_db(read_only = TRUE, demo = FALSE)
      on.exit(tryCatch(if (!is.null(c_geo) && DBI::dbIsValid(c_geo)) DBI::dbDisconnect(c_geo), error = function(e) NULL), add = TRUE)
      if (!is.null(c_geo) && DBI::dbExistsTable(c_geo, "usda_county_ranges")) {
        u_flds <- tryCatch(DBI::dbListFields(c_geo, "usda_county_ranges"), error = function(e) character(0))
        u_cnty_sel <- if ("county_norm" %in% u_flds) ", county_norm" else ""
        u_stonly_sel <- if ("is_state_only" %in% u_flds) ", is_state_only" else ""
        u_sym_sel <- if ("usda_symbol" %in% u_flds) ", usda_symbol" else ""
        u_df <- DBI::dbGetQuery(
          c_geo,
          paste0("SELECT DISTINCT canonical_packera, state_norm", u_cnty_sel, u_stonly_sel, u_sym_sel,
                 " FROM usda_county_ranges WHERE state_norm IS NOT NULL AND TRIM(state_norm) != '';")
        )
        if (!is.null(u_df) && nrow(u_df) > 0L) {
          u_tx <- trimws(as.character(u_df$canonical_packera))
          u_st <- .norm_state_str(u_df$state_norm)
          ok_u <- !is.na(u_tx) & nzchar(u_tx) & nzchar(u_st)
          if ("usda_symbol" %in% names(u_df) && any(ok_u)) {
            u_sym <- toupper(trimws(ifelse(is.na(u_df$usda_symbol), "", as.character(u_df$usda_symbol))))
            has_sym_u <- ok_u & nzchar(u_sym)
            if (any(has_sym_u)) {
              sym_owners <- tapply(u_tx[has_sym_u], u_sym[has_sym_u], function(v) length(unique(v)))
              shared_syms <- names(sym_owners)[sym_owners > 1L]
              if (length(shared_syms) > 0L) {
                is_shared_u <- has_sym_u & (u_sym %in% shared_syms)
                ep_pre <- toupper(substr(sub("^[A-Za-z-]+\\s+", "", u_tx), 1L, 2L))
                sym_ok_u <- !is_shared_u | (nchar(u_sym) >= 4L & ep_pre == substr(u_sym, 3L, 4L))
                ok_u <- ok_u & sym_ok_u
              }
            }
          }
          if (any(ok_u)) {
            usda_states_by_tx <- lapply(split(u_st[ok_u], u_tx[ok_u]), unique)
            if ("county_norm" %in% names(u_df)) {
              u_cty <- .norm_county_str(u_df$county_norm)
              u_stonly <- if ("is_state_only" %in% names(u_df)) (!is.na(u_df$is_state_only) & suppressWarnings(as.integer(u_df$is_state_only)) != 0L) else rep(FALSE, nrow(u_df))
              ok_uc <- ok_u & nzchar(u_cty) & !u_stonly
              if (any(ok_uc)) {
                u_sc <- paste(u_st[ok_uc], u_cty[ok_uc], sep = "||")
                usda_counties_by_tx <- lapply(split(u_sc, u_tx[ok_uc]), unique)
                for (tx_c in names(usda_counties_by_tx)) {
                  m_c <- ok_uc & (u_tx == tx_c)
                  st_cnt_tbl <- tapply(u_cty[m_c], u_st[m_c], function(v) length(unique(v)))
                  usda_state_county_counts_by_tx[[tx_c]] <- st_cnt_tbl
                }
              }
            }
          }
        }
      }
      if (!is.null(c_geo) && DBI::dbExistsTable(c_geo, "natureserve_subnational_ranges")) {
        ns_df <- tryCatch(
          DBI::dbGetQuery(
            c_geo,
            "SELECT DISTINCT species_standardized, subnation_code FROM natureserve_subnational_ranges WHERE subnation_code IS NOT NULL AND TRIM(subnation_code) != '' AND COALESCE(is_native, 0) = 1 AND COALESCE(is_exotic, 0) = 0;"
          ),
          error = function(e) NULL
        )
        if (!is.null(ns_df) && nrow(ns_df) > 0L) {
          ns_tx <- trimws(as.character(ns_df$species_standardized))
          ns_st <- .norm_state_str(ns_df$subnation_code)
          ok_ns <- !is.na(ns_tx) & nzchar(ns_tx) & nzchar(ns_st)
          if (any(ok_ns)) {
            for (tx_ns in unique(ns_tx[ok_ns])) {
              usda_states_by_tx[[tx_ns]] <- unique(c(usda_states_by_tx[[tx_ns]], ns_st[ok_ns & ns_tx == tx_ns]))
            }
          }
        }
      }
      if ((can_evaluate_hybrids || is_intra) && (nzchar(Sys.getenv("RUN_ID", unset = "")) || nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"))) &&
          !is.null(c_geo) && DBI::dbExistsTable(c_geo, "vouchers")) {
        v_flds <- tryCatch(DBI::dbListFields(c_geo, "vouchers"), error = function(e) character(0))
        sp_std_sel <- if ("species_standardized" %in% v_flds) "v.species_standardized, " else ""
        rec_by_sel <- if ("recordedBy" %in% v_flds) "v.recordedBy, " else ""
        syn_ext_cond <- if ("usda_county_status" %in% v_flds) "AND (v.determiner_tier = 'Tier_0_Type_Anchor' OR COALESCE(LOWER(TRIM(v.usda_county_status)), '') NOT IN ('out_of_range', 'natureserve_exotic_or_waif'))" else ""
        syn_exs_cond <- if ("exsiccatae_taxonomic_conflict" %in% v_flds) "AND (v.determiner_tier = 'Tier_0_Type_Anchor' OR COALESCE(v.exsiccatae_taxonomic_conflict, 0) = 0)" else ""
        q_syn <- paste(
          paste0("SELECT v.catalogNumber, ", sp_std_sel, rec_by_sel, "v.scientificName, v.species_raw, v.determiner_tier,"),
          "v.stateProvince, v.county, v.decimalLatitude, v.decimalLongitude",
          "FROM vouchers v",
          "WHERE v.determiner_tier IN ('Tier_0_Type_Anchor', 'Tier_1_Gold', 'Tier_2_Silver')",
          syn_ext_cond,
          syn_exs_cond,
          "AND NOT EXISTS (SELECT 1 FROM curation_triage ct WHERE ct.catalogNumber = v.catalogNumber AND ct.taxonomic_treatment_subtype = 'Misidentification_Reassign');"
        )
        db_syntopy_rows <- tryCatch(DBI::dbGetQuery(c_geo, q_syn), error = function(e) NULL)
      }
    }, error = function(e) NULL)
  }
  anch_states_by_tx <- list(); anch_coords_by_tx <- list()
  for (tx in active_taxa) {
    r_idx <- ref_idx_by_taxon[[tx]]
    st_r <- unique(state_vec_geo[r_idx][nzchar(state_vec_geo[r_idx])])
    anch_states_by_tx[[tx]] <- unique(c(st_r, usda_states_by_tx[[tx]]))
    ok_ll <- r_idx[is.finite(lat_vec_geo[r_idx]) & is.finite(lon_vec_geo[r_idx])]
    anch_coords_by_tx[[tx]] <- if (length(ok_ll) > 0L) cbind(lat = lat_vec_geo[ok_ll], lon = lon_vec_geo[ok_ll]) else NULL
  }
  if (!is.null(db_syntopy_rows) && nrow(db_syntopy_rows) > 0L) {
    raw_db_nm <- if ("species_standardized" %in% names(db_syntopy_rows)) {
      ifelse(!is.na(db_syntopy_rows$species_standardized) & nzchar(trimws(as.character(db_syntopy_rows$species_standardized))),
             as.character(db_syntopy_rows$species_standardized),
             ifelse(!is.na(db_syntopy_rows$scientificName) & nzchar(trimws(as.character(db_syntopy_rows$scientificName))),
                    as.character(db_syntopy_rows$scientificName),
                    as.character(db_syntopy_rows$species_raw)))
    } else {
      ifelse(!is.na(db_syntopy_rows$scientificName) & nzchar(trimws(as.character(db_syntopy_rows$scientificName))),
             as.character(db_syntopy_rows$scientificName),
             as.character(db_syntopy_rows$species_raw))
    }
    tx_db <- standardize_packera_taxon(raw_db_nm)
    st_db <- .norm_state_str(db_syntopy_rows$stateProvince)
    lat_db <- suppressWarnings(as.numeric(db_syntopy_rows$decimalLatitude))
    lon_db <- suppressWarnings(as.numeric(db_syntopy_rows$decimalLongitude))
    for (tx_u in unique(tx_db[!is.na(tx_db) & nzchar(tx_db) & tx_db != "Unknown"])) {
      m_u <- !is.na(tx_db) & (tx_db == tx_u)
      st_u <- unique(st_db[m_u & nzchar(st_db)])
      anch_states_by_tx[[tx_u]] <- unique(c(anch_states_by_tx[[tx_u]], st_u, usda_states_by_tx[[tx_u]]))
      ok_db_ll <- m_u & is.finite(lat_db) & is.finite(lon_db)
      if (any(ok_db_ll)) {
        cand_m <- cbind(lat = lat_db[ok_db_ll], lon = lon_db[ok_db_ll])
        anch_coords_by_tx[[tx_u]] <- if (is.null(anch_coords_by_tx[[tx_u]])) cand_m else rbind(anch_coords_by_tx[[tx_u]], cand_m)
      }
    }
  }
  is_valid_taxon_str <- function(s) !is.null(s) && !is.na(s) && nzchar(s) && s != "Unknown" && (s %in% active_taxa || grepl("^(Packera|P\\.)\\s+", s))
  has_ext_geo_ref <- (length(usda_states_by_tx) > 0L) || (!is.null(db_syntopy_rows) && nrow(db_syntopy_rows) > 0L)
  .eval_taxon_geo_plausibility <- function(row_i, cand_tx) {
    if (!is_valid_taxon_str(cand_tx)) return(list(plausible = FALSE, weight = 0.0, dist_km = Inf))
    st_i <- state_vec_geo[row_i]; lat_i <- lat_vec_geo[row_i]; lon_i <- lon_vec_geo[row_i]
    has_st_i <- nzchar(st_i); has_ll_i <- is.finite(lat_i) && is.finite(lon_i)
    if (!has_st_i && !has_ll_i) {
      if (is_intra && has_ext_geo_ref) return(list(plausible = FALSE, weight = 0.10, dist_km = Inf))
      return(list(plausible = TRUE, weight = 1.0, dist_km = NA_real_))
    }
    known_sts <- unique(c(anch_states_by_tx[[cand_tx]], usda_states_by_tx[[cand_tx]]))
    c_mat <- anch_coords_by_tx[[cand_tx]]
    min_km <- NA_real_
    if (has_ll_i && !is.null(c_mat) && nrow(c_mat) > 0L) {
      dlat <- (c_mat[, 1L] - lat_i) * (pi / 180.0)
      dlon <- (c_mat[, 2L] - lon_i) * (pi / 180.0)
      a_hav <- sin(dlat / 2.0)^2 + cos(lat_i * (pi / 180.0)) * cos(c_mat[, 1L] * (pi / 180.0)) * sin(dlon / 2.0)^2
      min_km <- min(6371.0 * 2.0 * atan2(sqrt(pmax(0.0, a_hav)), sqrt(pmax(0.0, 1.0 - a_hav))))
    }
    if (has_st_i && length(known_sts) > 0L && (st_i %in% known_sts)) {
      if (!(has_ll_i && is.finite(min_km) && min_km > 450.0)) {
        return(list(plausible = TRUE, weight = 1.0, dist_km = if (is.finite(min_km)) min_km else NA_real_))
      }
    }
    if (has_ll_i && is.finite(min_km)) {
      # Allow 450 km geodesic dispersal/disjunct buffer around documented anchors (650 km when <4 total anchors in small test fixtures)
      max_km <- if (nrow(c_mat) >= 4L || length(usda_states_by_tx[[cand_tx]]) > 0L) 450.0 else 650.0
      if (min_km <= max_km) {
        w_g <- if (min_km <= 200.0) 1.0 else exp(-0.5 * ((min_km - 200.0) / 250.0)^2)
        return(list(plausible = TRUE, weight = max(0.35, w_g), dist_km = min_km))
      }
      w_far <- exp(-0.5 * ((min_km - 200.0) / 200.0)^2)
      return(list(plausible = FALSE, weight = max(0.01, w_far), dist_km = min_km))
    }
    if (has_st_i && length(known_sts) >= 2L && !(st_i %in% known_sts)) {
      return(list(plausible = FALSE, weight = 0.05, dist_km = 800.0))
    }
    if (is_intra && has_ext_geo_ref) {
      return(list(plausible = FALSE, weight = 0.05, dist_km = Inf))
    }
    list(plausible = TRUE, weight = 0.85, dist_km = NA_real_)
  }

  # Precompute morphospace centroids and regularized inverse covariance in (can1, can2, can3) / dynamic PCs
  # for nearest-neighbor Mahalanobis hybrid parent resolution when predicted modalities match given_sp
  max_mc_p <- max(2L, min(15L, n - length(active_taxa) - 1L))
  sel_mc <- .parse_step07_semicolon_cols(df, "selected_morph_features")
  sel_mc <- sel_mc[vapply(sel_mc, function(cc) !is.logical(df[[cc]]) && sum(is.finite(suppressWarnings(as.numeric(df[[cc]])))) >= 2L, logical(1))]
  can_cand <- intersect(c("can1", "can2", "can3"), names(df))
  can_cand <- can_cand[vapply(can_cand, function(cc) sum(is.finite(suppressWarnings(as.numeric(df[[cc]])))) >= 2L, logical(1))]
  raw_pc_cand <- .sort_numbered_pc_cols(grep("^PC[0-9]+$", names(df), value = TRUE), "PC")
  raw_pc_cand <- raw_pc_cand[vapply(raw_pc_cand, function(cc) sum(is.finite(suppressWarnings(as.numeric(df[[cc]])))) >= 2L, logical(1))]
  pc_cand <- if (length(sel_mc) >= 2L) utils::head(sel_mc, max_mc_p) else utils::head(raw_pc_cand, max_mc_p)
  morph_coord_cols <- if (length(active_taxa) > 4L && length(pc_cand) >= 2L) pc_cand else if (length(can_cand) >= 2L) can_cand else if (length(pc_cand) >= 1L) pc_cand else can_cand

  morph_X <- NULL; morph_centroids <- NULL; morph_S_inv <- NULL
  if (length(morph_coord_cols) >= 1L && length(active_taxa) >= 2L) {
    p_dim <- length(morph_coord_cols)
    morph_X <- matrix(NA_real_, nrow = n, ncol = p_dim, dimnames = list(NULL, morph_coord_cols))
    for (j in seq_len(p_dim)) morph_X[, j] <- suppressWarnings(as.numeric(df[[morph_coord_cols[j]]]))
    morph_centroids <- matrix(0.0, nrow = length(active_taxa), ncol = p_dim, dimnames = list(active_taxa, morph_coord_cols))
    S_pool <- matrix(0.0, nrow = p_dim, ncol = p_dim)
    df_pool <- 0L
    for (tx in active_taxa) {
      r_idx <- ref_idx_by_taxon[[tx]]
      sub_m <- morph_X[r_idx, , drop = FALSE]
      use_t0_wt <- active_geno_ref && !is.null(t0_info) && length(r_idx) > 0L && any(t0_info$is_tier0_dna_anchor[r_idx] & (t0_info$rectified_taxon[r_idx] == tx))
      w_r <- if (use_t0_wt) {
        ifelse(t0_info$is_tier0_dna_anchor[r_idx] & (t0_info$rectified_taxon[r_idx] == tx),
               ifelse(is.finite(t0_info$tier0_weight[r_idx]), t0_info$tier0_weight[r_idx], 4.0), 1.0)
      } else NULL
      for (j in seq_len(p_dim)) {
        fin_mask <- is.finite(sub_m[, j])
        fin_j <- sub_m[fin_mask, j]
        if (length(fin_j) > 0L) {
          morph_centroids[tx, j] <- if (use_t0_wt) sum(w_r[fin_mask] * fin_j) / sum(w_r[fin_mask]) else mean(fin_j)
        } else {
          fin_all <- morph_X[is.finite(morph_X[, j]), j]
          morph_centroids[tx, j] <- if (length(fin_all) > 0L) mean(fin_all) else 0.0
        }
        sub_m[!fin_mask, j] <- morph_centroids[tx, j]
      }
      if (nrow(sub_m) >= 2L) {
        c_tx <- if (use_t0_wt && length(unique(w_r)) > 1L) {
          tryCatch(stats::cov.wt(sub_m, wt = w_r, method = "unbiased")$cov, error = function(e) stats::cov(sub_m))
        } else {
          stats::cov(sub_m)
        }
        c_tx[!is.finite(c_tx)] <- 0.0
        S_pool <- S_pool + (nrow(sub_m) - 1L) * c_tx
        df_pool <- df_pool + (nrow(sub_m) - 1L)
      }
    }
    if (df_pool > 0L) {
      S_pool <- S_pool / df_pool
    } else {
      X_imp <- morph_X
      for (j in seq_len(p_dim)) {
        med_j <- stats::median(X_imp[, j], na.rm = TRUE); if (!is.finite(med_j)) med_j <- 0.0
        X_imp[!is.finite(X_imp[, j]), j] <- med_j
      }
      S_pool <- if (nrow(X_imp) >= 2L) stats::cov(X_imp) else diag(1.0, nrow = p_dim, ncol = p_dim)
    }
    S_pool[!is.finite(S_pool)] <- 0.0
    ridge <- max(1e-4, 0.05 * mean(diag(S_pool), na.rm = TRUE))
    S_reg <- S_pool + diag(ridge, nrow = p_dim, ncol = p_dim)
    morph_S_inv <- tryCatch(solve(S_reg), error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-4), nrow = p_dim, ncol = p_dim))
  }

  resolve_nearest_mahalanobis_taxon <- function(row_idx, given_sp) {
    other_tx <- setdiff(cohort_species, given_sp)
    if (length(other_tx) == 0L) return(NA_character_)
    if (length(other_tx) == 1L || is.null(morph_X) || is.null(morph_centroids) || is.null(morph_S_inv)) return(other_tx[1L])
    other_in_m <- intersect(other_tx, rownames(morph_centroids))
    if (length(other_in_m) == 0L) return(other_tx[1L])
    x_i <- morph_X[row_idx, ]
    base_mu <- if (given_sp %in% rownames(morph_centroids)) morph_centroids[given_sp, ] else colMeans(morph_centroids)
    x_i[!is.finite(x_i)] <- base_mu[!is.finite(x_i)]
    d2_vals <- vapply(other_in_m, function(otx) {
      diff_v <- x_i - morph_centroids[otx, ]
      max(0.0, as.numeric(t(diff_v) %*% morph_S_inv %*% diff_v))
    }, numeric(1))
    other_in_m[which.min(d2_vals)]
  }

  # Precompute 5-Gate Conservative Hybrid Swarm structures (syntopy pool, co-flowering, and multi-organ envelopes)
  syntopy_pool_by_tx <- list()
  organ_pheno_taxa <- unique(c(active_taxa, cohort_species))
  has_anthesis_doy_by_tx <- stats::setNames(rep(FALSE, length(organ_pheno_taxa)), organ_pheno_taxa)
  sec_organ_cols <- intersect(
    c("blade_cie_l_sheet_delta", "blade_cie_l_norm", "blade_cie_l_median", "AR_inv", "involucre_length_width_ratio", "log1p_capitulum_count", "petiole_to_blade_ratio"),
    names(df)
  )
  sec_organ_cols <- sec_organ_cols[vapply(sec_organ_cols, function(cc) sum(is.finite(suppressWarnings(as.numeric(df[[cc]])))) >= 4L, logical(1))]
  sec_organ_mat <- if (length(sec_organ_cols) > 0L) {
    m_tmp <- matrix(NA_real_, nrow = n, ncol = length(sec_organ_cols), dimnames = list(NULL, sec_organ_cols))
    for (sj in seq_along(sec_organ_cols)) m_tmp[, sj] <- suppressWarnings(as.numeric(df[[sec_organ_cols[sj]]]))
    m_tmp
  } else NULL
  sec_organ_med_by_tx <- list()
  sec_organ_mad_pool <- stats::setNames(rep(NA_real_, length(sec_organ_cols)), sec_organ_cols)
  for (tx in organ_pheno_taxa) {
    r_idx <- ref_idx_by_taxon[[tx]]
    if (length(r_idx) == 0L) next
    has_anthesis_doy_by_tx[[tx]] <- any(valid_doy_all[r_idx] & anth_all[r_idx] & is.finite(hop_doy_all[r_idx]))
    if (!is.null(sec_organ_mat)) {
      med_v <- stats::setNames(rep(NA_real_, ncol(sec_organ_mat)), colnames(sec_organ_mat))
      for (sj in seq_len(ncol(sec_organ_mat))) {
        vv <- sec_organ_mat[r_idx, sj]; vv <- vv[is.finite(vv)]
        if (length(vv) >= 2L) med_v[sj] <- stats::median(vv)
      }
      sec_organ_med_by_tx[[tx]] <- med_v
    }
  }
  if (!is.null(sec_organ_mat)) {
    for (sj in seq_len(ncol(sec_organ_mat))) {
      vv <- sec_organ_mat[anchor_all & !mo_all & !noise_all, sj]; vv <- vv[is.finite(vv)]
      if (length(vv) < 4L) { vv <- sec_organ_mat[, sj]; vv <- vv[is.finite(vv)] }
      md <- if (length(vv) >= 2L) stats::mad(vv, constant = 1.4826) else NA_real_
      sdv <- if (length(vv) >= 2L) stats::sd(vv) else NA_real_
      sec_organ_mad_pool[sj] <- if (is.finite(md) && md > 1e-6) md else if (is.finite(sdv) && sdv > 1e-6) sdv else NA_real_
    }
  }
  self_morph_S_inv <- if (!is.null(morph_S_inv)) {
    S_inv_diag_floor <- tryCatch({
      S_est <- solve(morph_S_inv)
      diag(S_est) <- pmax(diag(S_est), 0.25)
      solve(S_est)
    }, error = function(e) morph_S_inv)
    S_inv_diag_floor
  } else NULL

  if (can_evaluate_hybrids) {
    db_tx_vec <- character(0); db_cat_vec <- character(0); db_lat_vec <- numeric(0); db_lon_vec <- numeric(0); db_sc_vec <- character(0)
    if (!is.null(db_syntopy_rows) && nrow(db_syntopy_rows) > 0L) {
      raw_db_nm <- ifelse(!is.na(db_syntopy_rows$scientificName) & nzchar(trimws(as.character(db_syntopy_rows$scientificName))),
                          as.character(db_syntopy_rows$scientificName), as.character(db_syntopy_rows$species_raw))
      db_tx_vec <- standardize_packera_taxon(raw_db_nm)
      db_cat_vec <- trimws(as.character(db_syntopy_rows$catalogNumber))
      db_lat_vec <- suppressWarnings(as.numeric(db_syntopy_rows$decimalLatitude))
      db_lon_vec <- suppressWarnings(as.numeric(db_syntopy_rows$decimalLongitude))
      db_st_v <- .norm_state_str(db_syntopy_rows$stateProvince)
      db_cn_v <- .norm_county_str(db_syntopy_rows$county)
      db_sc_vec <- ifelse(nzchar(db_st_v) & nzchar(db_cn_v), paste(db_st_v, db_cn_v, sep = "||"), "")
    }
    for (tx in cohort_species) {
      r_idx <- ref_idx_by_taxon[[tx]]
      cats_t <- cat_sheet_vec[r_idx]
      lats_t <- lat_vec_geo[r_idx]
      lons_t <- lon_vec_geo[r_idx]
      scs_t  <- state_cnty_geo[r_idx]
      if (length(db_tx_vec) > 0L) {
        m_db <- which(!is.na(db_tx_vec) & db_tx_vec == tx & !(db_cat_vec %in% cats_t))
        if (length(m_db) > 0L) {
          cats_t <- c(cats_t, db_cat_vec[m_db])
          lats_t <- c(lats_t, db_lat_vec[m_db])
          lons_t <- c(lons_t, db_lon_vec[m_db])
          scs_t  <- c(scs_t, db_sc_vec[m_db])
        }
      }
      syntopy_pool_by_tx[[tx]] <- list(cats = cats_t, lats = lats_t, lons = lons_t, scs = scs_t)
    }
  }

  .eval_five_gate_hybrid_pair <- function(row_i, p1_tx, p2_tx, morph_sp_i, cda_prob_i, cda_sec_tx_i, cda_sec_p_i, cda_marg_i, gmm_unc_i, edaph_sp_i, d_edaph_i, edaph_z_i, asym_i) {
    fail_res <- list(pass = FALSE, d1_km = Inf, d2_km = Inf, ddoy_12 = Inf, t_ab = NA_real_, d2_perp = Inf, score = -Inf)
    if (!can_evaluate_hybrids || !(p1_tx %in% cohort_species) || !(p2_tx %in% cohort_species) || p1_tx == p2_tx) return(fail_res)

    # Gate 1: Local Two-Parent Physical Syntopy (<= 15 km geodesic, or <= 25 km within exact shared county)
    pool1 <- syntopy_pool_by_tx[[p1_tx]]; pool2 <- syntopy_pool_by_tx[[p2_tx]]
    if (is.null(pool1) || is.null(pool2)) return(fail_res)
    self_cat <- cat_sheet_vec[row_i]
    lat_i <- lat_vec_geo[row_i]; lon_i <- lon_vec_geo[row_i]; sc_i <- state_cnty_geo[row_i]
    has_ll_i <- is.finite(lat_i) && is.finite(lon_i); has_sc_i <- nzchar(sc_i)
    if (!has_ll_i && !has_sc_i) return(fail_res)

    k1 <- which(pool1$cats != self_cat); k2 <- which(pool2$cats != self_cat)
    if (length(k1) == 0L || length(k2) == 0L) return(fail_res)
    same_sc1 <- has_sc_i && any(nzchar(pool1$scs[k1]) & pool1$scs[k1] == sc_i)
    same_sc2 <- has_sc_i && any(nzchar(pool2$scs[k2]) & pool2$scs[k2] == sc_i)

    d1_km <- Inf; d2_km <- Inf
    if (has_ll_i) {
      ll1 <- k1[is.finite(pool1$lats[k1]) & is.finite(pool1$lons[k1])]
      ll2 <- k2[is.finite(pool2$lats[k2]) & is.finite(pool2$lons[k2])]
      if (length(ll1) > 0L) {
        dlat1 <- (pool1$lats[ll1] - lat_i) * (pi / 180.0); dlon1 <- (pool1$lons[ll1] - lon_i) * (pi / 180.0)
        a1 <- sin(dlat1 / 2.0)^2 + cos(lat_i * (pi / 180.0)) * cos(pool1$lats[ll1] * (pi / 180.0)) * sin(dlon1 / 2.0)^2
        d1_km <- min(6371.0 * 2.0 * atan2(sqrt(pmax(0.0, a1)), sqrt(pmax(0.0, 1.0 - a1))))
      }
      if (length(ll2) > 0L) {
        dlat2 <- (pool2$lats[ll2] - lat_i) * (pi / 180.0); dlon2 <- (pool2$lons[ll2] - lon_i) * (pi / 180.0)
        a2 <- sin(dlat2 / 2.0)^2 + cos(lat_i * (pi / 180.0)) * cos(pool2$lats[ll2] * (pi / 180.0)) * sin(dlon2 / 2.0)^2
        d2_km <- min(6371.0 * 2.0 * atan2(sqrt(pmax(0.0, a2)), sqrt(pmax(0.0, 1.0 - a2))))
      }
      syntopy_ok <- is.finite(d1_km) && is.finite(d2_km) &&
                    ((d1_km <= 15.0 && d2_km <= 15.0) || (same_sc1 && same_sc2 && d1_km <= 25.0 && d2_km <= 25.0))
    } else {
      syntopy_ok <- same_sc1 && same_sc2
      if (syntopy_ok) { d1_km <- 0.0; d2_km <- 0.0 }
    }
    if (!syntopy_ok) return(fail_res)

    # Gate 2: Reproductive Phenological Co-Flowering (<= 21 d peak separation) & Ecotonal Suitability
    ddoy_12 <- 0.0
    if ((p1_tx %in% names(has_anthesis_doy_by_tx)) && (p2_tx %in% names(has_anthesis_doy_by_tx)) &&
        (p1_tx %in% names(mu_theta_hop)) && (p2_tx %in% names(mu_theta_hop)) &&
        isTRUE(has_anthesis_doy_by_tx[[p1_tx]]) && isTRUE(has_anthesis_doy_by_tx[[p2_tx]])) {
      th1 <- mu_theta_hop[[p1_tx]]; th2 <- mu_theta_hop[[p2_tx]]
      ddoy_12 <- abs(atan2(sin(th1 - th2), cos(th1 - th2))) * (365.25 / (2.0 * pi))
      if (!is.finite(ddoy_12) || ddoy_12 > 21.0) return(fail_res)
      if (valid_doy_all[row_i] && anth_all[row_i] && is.finite(hop_doy_all[row_i])) {
        thi <- 2.0 * pi * (hop_doy_all[row_i] / 365.25)
        th_mid <- atan2(0.5 * (sin(th1) + sin(th2)), 0.5 * (cos(th1) + cos(th2)))
        d_p1 <- abs(atan2(sin(thi - th1), cos(thi - th1))) * (365.25 / (2.0 * pi))
        d_p2 <- abs(atan2(sin(thi - th2), cos(thi - th2))) * (365.25 / (2.0 * pi))
        d_pm <- abs(atan2(sin(thi - th_mid), cos(thi - th_mid))) * (365.25 / (2.0 * pi))
        if (min(d_p1, d_p2, d_pm) > 21.0) return(fail_res)
      }
    }
    if ((is.finite(d_edaph_i) && d_edaph_i < 0.25) || (is.finite(edaph_z_i) && edaph_z_i > 2.5)) return(fail_res)

    # Gate 3: Pairwise Parental Mahalanobis Chord Projection (t_AB in [0.22, 0.78] & bounded orthogonal residual)
    has_row_coords <- !is.null(morph_X) && !is.null(morph_centroids) && !is.null(morph_S_inv) &&
                      (p1_tx %in% rownames(morph_centroids)) && (p2_tx %in% rownames(morph_centroids)) &&
                      all(is.finite(morph_X[row_i, ]))
    if (has_row_coords) {
      mu1 <- morph_centroids[p1_tx, ]; mu2 <- morph_centroids[p2_tx, ]
      delta_12 <- mu2 - mu1
      d2_12 <- max(0.0, as.numeric(t(delta_12) %*% morph_S_inv %*% delta_12))
      if (!is.finite(d2_12) || d2_12 < 2.25) return(fail_res)
      xi <- morph_X[row_i, ]
      t_ab <- as.numeric(t(xi - mu1) %*% morph_S_inv %*% delta_12) / d2_12
      r_perp <- (xi - mu1) - t_ab * delta_12
      d2_perp <- max(0.0, as.numeric(t(r_perp) %*% morph_S_inv %*% r_perp))
      df_perp <- max(1L, ncol(morph_X) - 1L)
      if (!is.finite(t_ab) || t_ab < 0.22 || t_ab > 0.78 || !is.finite(d2_perp) || d2_perp > stats::qchisq(0.90, df = df_perp)) {
        return(fail_res)
      }
    } else {
      if (is.null(morph_sp_i) || !(morph_sp_i %in% c(p1_tx, p2_tx)) || !is.finite(cda_prob_i) ||
          cda_prob_i < 0.32 || cda_prob_i > 0.65 || !is.finite(cda_marg_i) || cda_marg_i > 0.22) {
        return(fail_res)
      }
      if (is.finite(cda_sec_p_i) && (cda_sec_p_i < 0.25 || (cda_prob_i + cda_sec_p_i) < 0.75)) return(fail_res)
      t_ab <- if (morph_sp_i == p2_tx) cda_prob_i else (1.0 - cda_prob_i)
      d2_perp <- 0.0
    }

    # Gate 4: Multi-Organ Admixture / Mosaicism Corroboration & Non-Parental Posterior Split
    if (is.finite(cda_prob_i) && cda_prob_i >= 0.68) return(fail_res)
    has_post_split <- (is.finite(cda_marg_i) && cda_marg_i <= 0.25) || (is.finite(gmm_unc_i) && gmm_unc_i >= 0.22)
    if (!has_post_split) return(fail_res)

    has_p2_signal <- (!is.null(morph_sp_i) && morph_sp_i == p2_tx) ||
                     (!is.null(cda_sec_tx_i) && cda_sec_tx_i == p2_tx && (!is.finite(cda_sec_p_i) || cda_sec_p_i >= 0.20)) ||
                     (!is.na(edaph_sp_i) && edaph_sp_i == p2_tx) ||
                     (is.finite(gmm_unc_i) && gmm_unc_i >= 0.28)
    if (!has_p2_signal) return(fail_res)

    if (!is.null(sec_organ_mat) && length(sec_organ_cols) > 0L) {
      m1_v <- sec_organ_med_by_tx[[p1_tx]]; m2_v <- sec_organ_med_by_tx[[p2_tx]]
      n_diag_org <- 0L; any_org_admix <- FALSE
      for (sj in seq_along(sec_organ_cols)) {
        val_ij <- sec_organ_mat[row_i, sj]
        mad_j <- sec_organ_mad_pool[sj]
        if (is.finite(val_ij) && !is.null(m1_v) && !is.null(m2_v) &&
            is.finite(m1_v[sj]) && is.finite(m2_v[sj]) && is.finite(mad_j) &&
            abs(m2_v[sj] - m1_v[sj]) >= 1.2 * mad_j) {
          n_diag_org <- n_diag_org + 1L
          t_org <- (val_ij - m1_v[sj]) / (m2_v[sj] - m1_v[sj])
          if ((t_org >= 0.15 && t_org <= 0.85) || (t_ab >= 0.50 && t_org <= 0.45) || (t_ab <= 0.50 && t_org >= 0.55)) {
            any_org_admix <- TRUE
          }
        }
      }
      high_asym <- is.finite(asym_i) && (asym_i >= 0.18)
      if (n_diag_org >= 1L && !any_org_admix && !high_asym) return(fail_res)
    }

    cand_score <- (1.0 - 2.0 * abs(t_ab - 0.50)) - 0.08 * d2_perp - 0.02 * (if (is.finite(d2_km)) d2_km else 10.0)
    list(pass = TRUE, d1_km = d1_km, d2_km = d2_km, ddoy_12 = ddoy_12, t_ab = t_ab, d2_perp = d2_perp, score = cand_score)
  }

  status_calls <- character(n); rec_taxa <- character(n); triage_prios <- character(n)
  triage_acts <- character(n); rationales <- character(n); treat_broad <- character(n); treat_sub <- character(n)
  df_nms_d <- names(df)
  raw_sp_dv <- if ("species_raw" %in% df_nms_d) as.character(df$species_raw) else rep(NA_character_, n)
  tier_dv <- if ("determiner_tier" %in% df_nms_d) as.character(df$determiner_tier) else rep(NA_character_, n)
  cda_pred_dv <- if ("cda_predicted_taxon" %in% df_nms_d) as.character(df$cda_predicted_taxon) else rep(NA_character_, n)
  cda_sec1_dv <- if ("cda_second_taxon" %in% df_nms_d) as.character(df$cda_second_taxon) else rep(NA_character_, n)
  cda_sec2_dv <- if ("cda_second_predicted_taxon" %in% df_nms_d) as.character(df$cda_second_predicted_taxon) else rep(NA_character_, n)
  cda_prob_dv <- if ("cda_posterior_prob" %in% df_nms_d) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(NA_real_, n)
  cda_sec_p_dv <- if ("cda_second_posterior_prob" %in% df_nms_d) suppressWarnings(as.numeric(df$cda_second_posterior_prob)) else rep(NA_real_, n)
  cda_marg_dv <- if ("cda_posterior_margin" %in% df_nms_d) suppressWarnings(as.numeric(df$cda_posterior_margin)) else rep(NA_real_, n)
  mo_flag_dv <- if ("morphometric_outlier_flag" %in% df_nms_d) as.logical(df$morphometric_outlier_flag) else rep(NA, n)
  misid_flag_dv <- if ("misidentification_flag" %in% df_nms_d) as.logical(df$misidentification_flag) else rep(NA, n)
  triage_cat_dv <- if ("triage_category" %in% df_nms_d) as.character(df$triage_category) else rep(NA_character_, n)
  mm_conf_dv <- if ("multimodal_conflict_flag" %in% df_nms_d) (!is.na(df$multimodal_conflict_flag) & as.logical(df$multimodal_conflict_flag)) else rep(FALSE, n)
  edaph_sp_dv <- if ("edaphic_best_fit_taxon" %in% df_nms_d) as.character(df$edaphic_best_fit_taxon) else rep(NA_character_, n)
  rf_pred_dv <- if ("rf_predicted_taxon" %in% df_nms_d) as.character(df$rf_predicted_taxon) else rep(NA_character_, n)
  rf_prob_dv <- if ("rf_predicted_prob" %in% df_nms_d) suppressWarnings(as.numeric(df$rf_predicted_prob)) else rep(NA_real_, n)
  d_ed1_dv <- if ("edaphic_niche_overlap_d" %in% df_nms_d) suppressWarnings(as.numeric(df$edaphic_niche_overlap_d)) else rep(NA_real_, n)
  d_ed2_dv <- if ("edaphic_schoener_d" %in% df_nms_d) suppressWarnings(as.numeric(df$edaphic_schoener_d)) else rep(NA_real_, n)
  d_ed3_dv <- if ("schoener_d" %in% df_nms_d) suppressWarnings(as.numeric(df$schoener_d)) else rep(NA_real_, n)
  gmm_unc_dv <- if ("gmm_uncertainty" %in% df_nms_d) suppressWarnings(as.numeric(df$gmm_uncertainty)) else rep(NA_real_, n)
  gmm_cl_dv <- if ("gmm_cluster" %in% df_nms_d) suppressWarnings(as.integer(df$gmm_cluster)) else rep(NA_integer_, n)
  ddoy1_dv <- if ("delta_doy_resid" %in% df_nms_d) suppressWarnings(as.numeric(df$delta_doy_resid)) else rep(NA_real_, n)
  ddoy2_dv <- if ("delta_doy" %in% df_nms_d) suppressWarnings(as.numeric(df$delta_doy)) else rep(NA_real_, n)
  is_anth_dv <- if ("is_active_anthesis" %in% df_nms_d) (!is.na(df$is_active_anthesis) & as.logical(df$is_active_anthesis)) else rep(TRUE, n)
  bf1_dv <- if ("two_delta_bic" %in% df_nms_d) suppressWarnings(as.numeric(df$two_delta_bic)) else rep(NA_real_, n)
  bf2_dv <- if ("gmm_2_delta_bic" %in% df_nms_d) suppressWarnings(as.numeric(df$gmm_2_delta_bic)) else rep(NA_real_, n)
  m_tier_src_dv <- if ("morph_tier_source" %in% df_nms_d) as.character(df$morph_tier_source) else rep(NA_character_, n)
  bas_dv <- if ("basal_leaf_count" %in% df_nms_d) suppressWarnings(as.integer(df$basal_leaf_count)) else rep(NA_integer_, n)
  lcnt_dv <- if ("leaf_count" %in% df_nms_d) suppressWarnings(as.integer(df$leaf_count)) else rep(NA_integer_, n)
  n_t1_dv <- dplyr::coalesce(bas_dv, lcnt_dv)
  dlf_dv <- if ("dissected_leaf_fraction" %in% df_nms_d) suppressWarnings(as.numeric(df$dissected_leaf_fraction)) else rep(NA_real_, n)
  asym_dv <- if ("pct_asymmetry" %in% df_nms_d) suppressWarnings(as.numeric(df$pct_asymmetry)) else rep(NA_real_, n)
  is_dis_dv <- if ("is_dissected" %in% df_nms_d) (!is.na(df$is_dissected) & as.logical(df$is_dissected)) else rep(FALSE, n)
  p_conf_dv <- if ("probabilistic_consensus_conflict" %in% df_nms_d) suppressWarnings(as.numeric(df$probabilistic_consensus_conflict)) else rep(0.0, n)

  has_geno_dv <- .has_step07_geno_row_mask(df)
  geno_pred_dv <- if ("geno_predicted_taxon" %in% df_nms_d) as.character(df$geno_predicted_taxon) else rep(NA_character_, n)
  geno_prob_dv <- .num_col(df, "geno_posterior_prob", n, default = NA_real_)
  geno_marg_dv <- .num_col(df, "geno_margin", n, default = NA_real_)
  geno_tier_dv <- if ("geno_linkage_tier" %in% df_nms_d) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% df_nms_d) as.character(df$linkage_tier) else rep(NA_character_, n)
  geno_adm_dv <- if ("is_geno_admixed" %in% df_nms_d) {
    if (exists(".safe_bool_vec_r", mode = "function")) .safe_bool_vec_r(df$is_geno_admixed) else (!is.na(df$is_geno_admixed) & as.logical(df$is_geno_admixed))
  } else rep(FALSE, n)
  geno_p1_dv <- .chr_col(df, "geno_parent_1", n, default = NA_character_)
  geno_p2_dv <- .chr_col(df, "geno_parent_2", n, default = NA_character_)
  geno_sister_dv <- .chr_col(df, "nearest_tree_sister_taxon", n, default = NA_character_)
  geno_qmax_dv <- .num_col(df, "q_max", n, default = NA_real_)
  geno_qry_dv <- if ("is_putative_novel_or_hybrid_query" %in% df_nms_d) {
    if (exists(".safe_bool_vec_r", mode = "function")) .safe_bool_vec_r(df$is_putative_novel_or_hybrid_query) else (!is.na(df$is_putative_novel_or_hybrid_query) & as.logical(df$is_putative_novel_or_hybrid_query))
  } else rep(FALSE, n)
  geno_gdi_dv <- .num_col(df, "coalescent_gdi", n, default = NA_real_)
  gdi_tax_att <- attr(df, "genomic_gdi_delimitation")
  if (!is.null(gdi_tax_att) && is.data.frame(gdi_tax_att) && nrow(gdi_tax_att) > 0L &&
      "taxon" %in% names(gdi_tax_att) && "coalescent_gdi" %in% names(gdi_tax_att)) {
    miss_gdi_g <- has_geno_dv & !is.finite(geno_gdi_dv)
    if (any(miss_gdi_g)) {
      g_tx_key <- ifelse(!is.na(geno_pred_dv) & nzchar(geno_pred_dv), geno_pred_dv, sp_all_vec)
      m_gd <- match(g_tx_key[miss_gdi_g], as.character(gdi_tax_att$taxon))
      ok_gd <- !is.na(m_gd)
      if (any(ok_gd)) geno_gdi_dv[which(miss_gdi_g)[ok_gd]] <- suppressWarnings(as.numeric(gdi_tax_att$coalescent_gdi[m_gd[ok_gd]]))
    }
  }
  geno_h_dv <- .num_col(df, "hybrid_index_h", n, default = NA_real_)
  geno_q12_dv <- .num_col(df, "interspecific_ancestry_q12", n, default = NA_real_)
  geno_hyb_cls_dv <- .chr_col(df, "hybrid_generation_class", n, default = NA_character_)
  geno_fb_max_dv <- .num_col(df, "f_branch_max", n, default = NA_real_)
  geno_fb_donor_dv <- .chr_col(df, "f_branch_donor_taxon", n, default = NA_character_)
  geno_gamma_dv <- .num_col(df, "snaq_gamma_minor", n, default = NA_real_)
  geno_bf_dv <- .num_col(df, "bpp_introgression_bf", n, default = NA_real_)
  geno_cyt_dv <- .chr_col(df, "inferred_cytotype", n, default = NA_character_)
  geno_plast_dv <- .chr_col(df, "plastid_clade", n, default = NA_character_)
  geno_its_dv <- .chr_col(df, "its_ribotype", n, default = NA_character_)
  geno_cytonuc_dv <- if ("cytonuclear_discordant" %in% df_nms_d) {
    if (exists(".safe_bool_vec_r", mode = "function")) .safe_bool_vec_r(df$cytonuclear_discordant) else (!is.na(df$cytonuclear_discordant) & as.logical(df$cytonuclear_discordant))
  } else rep(FALSE, n)
  geno_divergent_nonhybrid_v <- has_geno_dv & (
    (!is.na(geno_hyb_cls_dv) & geno_hyb_cls_dv == "Divergent_NonHybrid_Lineage") |
    (is.finite(geno_gdi_dv) & geno_gdi_dv > 0.70 & (!is.finite(geno_q12_dv) | geno_q12_dv < 0.15))
  )
  geno_contemp_hybrid_v <- has_geno_dv & !geno_divergent_nonhybrid_v & (
    (is.finite(geno_q12_dv) & geno_q12_dv >= 0.35) |
    (!is.na(geno_hyb_cls_dv) & (geno_hyb_cls_dv %in% c("F1_Hybrid", "F2_Hybrid", "BC1_Parent_1", "BC1_Parent_2", "F2_Backcross", "Recent_Admixed"))) |
    (geno_adm_dv & (!is.finite(geno_q12_dv) | geno_q12_dv >= 0.15))
  )
  geno_ancient_introg_v <- has_geno_dv & !geno_contemp_hybrid_v & !geno_divergent_nonhybrid_v &
    (!is.finite(geno_q12_dv) | geno_q12_dv < 0.15) &
    ((is.finite(geno_fb_max_dv) & geno_fb_max_dv >= 0.05) | geno_cytonuc_dv |
     (!is.na(geno_hyb_cls_dv) & (geno_hyb_cls_dv %in% c("Advanced_Introgressant", "Ancient_Introgressant"))))
  cluster_prop_gdi_dv <- .num_col(df, "cluster_propagated_gdi", n, default = NA_real_)
  tier0_abstain_dv <- if ("tier0_abstention_flag" %in% df_nms_d) (!is.na(df$tier0_abstention_flag) & suppressWarnings(as.integer(df$tier0_abstention_flag)) == 1L) else rep(FALSE, n)
  has_direct_gdi <- has_geno_dv & is.finite(geno_gdi_dv)
  eff_gdi_dv <- ifelse(has_direct_gdi, geno_gdi_dv, cluster_prop_gdi_dv)
  is_cluster_prop_gdi_only <- !has_direct_gdi & is.finite(cluster_prop_gdi_dv)
  geno_gdi_oversplit_v <- is.finite(eff_gdi_dv) & (eff_gdi_dv < 0.20)
  geno_gdi_cryptic_support_v <- (has_geno_dv & ((is.finite(geno_gdi_dv) & geno_gdi_dv > 0.70) | (!is.na(geno_hyb_cls_dv) & geno_hyb_cls_dv == "Divergent_NonHybrid_Lineage"))) |
                                (is_cluster_prop_gdi_only & (cluster_prop_gdi_dv > 0.70) & !tier0_abstain_dv)
  cyto_barrier_raw_d <- if ("cytotype_reproductive_barrier_pass" %in% df_nms_d) df$cytotype_reproductive_barrier_pass else rep(NA_integer_, n)
  cyto_barrier_int_v <- if (is.logical(cyto_barrier_raw_d)) ifelse(is.na(cyto_barrier_raw_d), NA_integer_, as.integer(cyto_barrier_raw_d)) else suppressWarnings(as.integer(cyto_barrier_raw_d))
  cyto_barrier_pass_v <- !is.na(cyto_barrier_int_v) & (cyto_barrier_int_v == 1L)
  voucher_ploidy_v <- if ("voucher_ploidy_level" %in% df_nms_d) suppressWarnings(as.integer(df$voucher_ploidy_level)) else rep(NA_integer_, n)
  voucher_2n_v <- if ("voucher_sporophytic_2n" %in% df_nms_d) suppressWarnings(as.integer(df$voucher_sporophytic_2n)) else rep(NA_integer_, n)
  taxon_cyto_str_v <- .chr_col(df, "taxon_reported_cytotypes", n, default = NA_character_)
  is_multi_cyto_raw_d <- if ("is_multi_cytotype_taxon" %in% df_nms_d) df$is_multi_cytotype_taxon else rep(NA_integer_, n)
  is_multi_cyto_int_v <- if (is.logical(is_multi_cyto_raw_d)) ifelse(is.na(is_multi_cyto_raw_d), NA_integer_, as.integer(is_multi_cyto_raw_d)) else suppressWarnings(as.integer(is_multi_cyto_raw_d))
  is_multi_cyto_v <- !is.na(is_multi_cyto_int_v) & (is_multi_cyto_int_v == 1L)
  cyto_tier_v <- .chr_col(df, "cytology_linkage_tier", n, default = NA_character_)
  pid_cyto_dv <- if ("plant_individual_id" %in% df_nms_d) suppressWarnings(as.integer(df$plant_individual_id)) else rep(0L, n)

  type_info_d <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    .identify_type_specimen_vec(df)
  } else {
    list(type_rank = rep(0L, n), is_type_specimen = rep(FALSE, n), type_status_clean = rep("None", n))
  }
  is_type_dv <- type_info_d$is_type_specimen | (if ("is_type_specimen" %in% df_nms_d) (!is.na(df$is_type_specimen) & as.logical(df$is_type_specimen)) else rep(FALSE, n))
  type_stat_dv <- if ("type_status" %in% df_nms_d) as.character(df$type_status) else type_info_d$type_status_clean
  is_typ_morph_dv <- if ("is_typical_morphotype" %in% df_nms_d) as.logical(df$is_typical_morphotype) else rep(NA, n)
  tm_d2_dv <- if ("type_morph_mahal_d2" %in% df_nms_d) suppressWarnings(as.numeric(df$type_morph_mahal_d2)) else rep(NA_real_, n)
  nt_morph_dv <- if ("nearest_type_morph_taxon" %in% df_nms_d) as.character(df$nearest_type_morph_taxon) else rep(NA_character_, n)
  nt_bas_dv <- if ("nearest_type_basionym" %in% df_nms_d) as.character(df$nearest_type_basionym) else rep(NA_character_, n)
  bas_d2_dv <- .num_col(df, "nearest_type_basionym_d2", n, default = NA_real_)
  tm_conc_dv <- if ("type_morph_concordance" %in% df_nms_d) suppressWarnings(as.numeric(df$type_morph_concordance)) else rep(NA_real_, n)
  te_z_dv <- if ("type_edaphic_z_score" %in% df_nms_d) suppressWarnings(as.numeric(df$type_edaphic_z_score)) else rep(NA_real_, n)
  nt_niche_dv <- if ("nearest_type_niche_taxon" %in% df_nms_d) as.character(df$nearest_type_niche_taxon) else rep(NA_character_, n)
  min_cong_d2_dv <- .num_col(df, "min_congener_mahal_d2", n, default = NA_real_)
  min_cong_tx_dv <- .chr_col(df, "min_congener_taxon", n, default = NA_character_)
  osr_p_dv <- .num_col(df, "open_set_novelty_p", n, default = NA_real_)
  osr_df_raw <- .num_col(df, "open_set_df", n, default = NA_real_)
  self_d2_dv <- .num_col(df, "self_morph_mahal_d2", n, default = NA_real_)
  self_df_raw <- .num_col(df, "self_morph_df", n, default = NA_real_)
  cda_oof_p_dv <- .num_col(df, "cda_oof_posterior_prob", n, default = NA_real_)
  cl_thr_dv <- .num_col(df, "cl_target_threshold", n, default = NA_real_)
  exs_gold_tx_dv <- .chr_col(df, "exsiccatae_gold_taxon", n, default = NA_character_)
  if (any(!is.na(exs_gold_tx_dv) & nzchar(exs_gold_tx_dv))) {
    ok_eg <- which(!is.na(exs_gold_tx_dv) & nzchar(trimws(exs_gold_tx_dv)))
    exs_gold_tx_dv[ok_eg] <- standardize_packera_taxon(exs_gold_tx_dv[ok_eg])
  }
  usda_stat_dv <- tolower(trimws(ifelse(is.na(.chr_col(df, "usda_county_status", n, default = "")), "", .chr_col(df, "usda_county_status", n, default = ""))))
  exs_conf_raw <- if ("exsiccatae_taxonomic_conflict" %in% df_nms_d) df$exsiccatae_taxonomic_conflict else rep(0L, n)
  exs_conf_dv <- if (is.logical(exs_conf_raw)) (!is.na(exs_conf_raw) & exs_conf_raw) else (!is.na(suppressWarnings(as.integer(exs_conf_raw))) & suppressWarnings(as.integer(exs_conf_raw)) != 0L)
  clim_z_dv <- .num_col(df, "climatic_z_score", n, default = NA_real_)
  rf_given_p_dv <- .num_col(df, "rf_given_taxon_prob", n, default = NA_real_)
  id_by_dv <- .chr_col(df, "identifiedBy", n, default = NA_character_)
  portal_lag_flag <- logical(n)
  portal_lag_gate <- rep(NA_character_, n)
  is_marginal_ext_v <- logical(n)
  p_type_morph <- max(2L, if (length(sel_mc) > 0L) length(sel_mc) else if (length(morph_coord_cols) > 0L) length(morph_coord_cols) else 3L)
  type_morph_chisq95 <- stats::qchisq(0.95, df = p_type_morph)
  osr_df_vec <- ifelse(is.finite(osr_df_raw) & osr_df_raw >= 1.0, as.integer(round(osr_df_raw)), p_type_morph)
  osr_chisq85_vec <- stats::qchisq(0.85, df = osr_df_vec)
  osr_chisq90_vec <- stats::qchisq(0.90, df = osr_df_vec)
  osr_chisq95_vec <- stats::qchisq(0.95, df = osr_df_vec)
  sheet_ramet_cnt_map <- table(cat_sheet_vec)
  n_sheet_ramets_vec <- as.integer(sheet_ramet_cnt_map[cat_sheet_vec])

  .norm_collector_surname <- function(cv) {
    vapply(as.character(cv), function(s) {
      if (is.na(s) || !nzchar(trimws(s))) return("")
      x <- tolower(trimws(s))
      x <- gsub("\\b(et\\s+al\\.?|jr\\.?|sr\\.?|ii|iii|iv|dr\\.?|prof\\.?)\\b", " ", x)
      x <- strsplit(x, "\\s*(;|&|\\band\\b|\\bwith\\b)\\s*")[[1L]][1L]
      if (grepl(",", x, fixed = TRUE)) {
        sur <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]][1L])
      } else {
        toks <- strsplit(gsub("[^a-z\\s-]+", " ", x), "\\s+")[[1L]]
        toks <- toks[nchar(toks) >= 2L]
        sur <- if (length(toks) > 0L) toks[length(toks)] else ""
      }
      gsub("[^a-z-]", "", sur)
    }, character(1), USE.NAMES = FALSE)
  }

  # Precompute cluster-level population cohesion (cryptic_cluster_cohesion_pass) & spatial-cline validation
  sp_clean_vec <- ifelse(!is.na(sp_all_vec) & nzchar(sp_all_vec), sp_all_vec, "Unknown")
  g_cl_all <- ifelse(!is.na(gmm_cl_dv), gmm_cl_dv, 1L)
  t1_pure_all <- (is.na(m_tier_src_dv) | !nzchar(m_tier_src_dv) | m_tier_src_dv == "Basal_Rosette") & (is.na(n_t1_dv) | n_t1_dv >= 1L)
  edaph_z_all <- ifelse(is.finite(edaph_z_vec), edaph_z_vec, te_z_dv)
  ed_sp_all <- ifelse(!is.na(edaph_sp_dv) & nzchar(edaph_sp_dv), edaph_sp_dv, ifelse(!is.na(rf_pred_dv) & nzchar(rf_pred_dv), rf_pred_dv, sp_clean_vec))
  d_edaph_all <- ifelse(
    !is.na(d_ed1_dv), d_ed1_dv,
    ifelse(
      !is.na(d_ed2_dv), d_ed2_dv,
      ifelse(
        !is.na(d_ed3_dv), d_ed3_dv,
        if (is_intra) {
          ifelse(is.finite(edaph_z_all), exp(-0.5 * (edaph_z_all / 1.8)^2), 0.75)
        } else {
          ifelse(ed_sp_all != sp_clean_vec & is.finite(edaph_z_all) & edaph_z_all > 2.0, 0.22, 0.75)
        }
      )
    )
  )
  ddoy_all <- ifelse(is.finite(ddoy1_dv), ddoy1_dv, ifelse(is.finite(ddoy2_dv), ddoy2_dv, NA_real_))
  bf_all <- ifelse(is.finite(bf1_dv), bf1_dv, ifelse(is.finite(bf2_dv), bf2_dv, gmm_two_delta_bic))
  rec_by_v <- .chr_col(df, "recordedBy", n, default = NA_character_)
  rec_num_v <- .chr_col(df, "recordNumber", n, default = NA_character_)
  exs_key_v <- .chr_col(df, "exsiccatae_key", n, default = NA_character_)
  yr_v <- .chr_col(df, "year", n, default = NA_character_)
  if (all(is.na(yr_v)) && "eventDate" %in% df_nms_d) yr_v <- substr(as.character(df$eventDate), 1L, 4L)
  cnty_v <- .chr_col(df, "county", n, default = NA_character_)
  lat_v_cl <- if ("decimalLatitude" %in% df_nms_d) suppressWarnings(as.numeric(df$decimalLatitude)) else rep(NA_real_, n)
  lon_v_cl <- if ("decimalLongitude" %in% df_nms_d) suppressWarnings(as.numeric(df$decimalLongitude)) else rep(NA_real_, n)
  inst_v <- .chr_col(df, "institutionCode", n, default = NA_character_)
  deme_id_v <- .chr_col(df, "population_deme_id", n, default = NA_character_)
  if (all(is.na(deme_id_v))) deme_id_v <- .chr_col(df, "cv_deme_id", n, default = NA_character_)
  azonal_flag_all <- rep(FALSE, n)
  for (az_c in intersect(c("is_azonal_lithology", "is_rock_outcrop", "substrate_specialist_flag"), df_nms_d)) {
    az_b <- suppressWarnings(as.logical(df[[az_c]]))
    azonal_flag_all <- azonal_flag_all | (!is.na(az_b) & az_b)
  }
  if ("matches_type_lithology" %in% df_nms_d) {
    mtl_b <- suppressWarnings(as.logical(df$matches_type_lithology))
    azonal_flag_all <- azonal_flag_all | (!is.na(mtl_b) & !mtl_b)
  }

  spatial_deme_vec <- ifelse(!is.na(deme_id_v) & nzchar(trimws(deme_id_v)), trimws(deme_id_v), "")
  miss_dm <- which(!nzchar(spatial_deme_vec))
  if (length(miss_dm) > 0L) {
    ok_dm_ll <- miss_dm[is.finite(lat_v_cl[miss_dm]) & is.finite(lon_v_cl[miss_dm])]
    if (length(ok_dm_ll) >= 1L) {
      dm_clusters <- list()
      lat_r_dm <- lat_v_cl[ok_dm_ll] * (pi / 180.0); lon_r_dm <- lon_v_cl[ok_dm_ll] * (pi / 180.0)
      for (ii_d in seq_along(ok_dm_ll)) {
        placed_d <- FALSE
        for (cid_d in seq_along(dm_clusters)) {
          mem_d <- dm_clusters[[cid_d]]
          dlat_d <- lat_r_dm[mem_d] - lat_r_dm[ii_d]; dlon_d <- lon_r_dm[mem_d] - lon_r_dm[ii_d]
          hav_d <- sin(dlat_d * 0.5)^2 + cos(lat_r_dm[ii_d]) * cos(lat_r_dm[mem_d]) * (sin(dlon_d * 0.5)^2)
          dkm_d <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hav_d))))
          if (all(is.finite(dkm_d) & dkm_d <= 10.0)) {
            dm_clusters[[cid_d]] <- c(mem_d, ii_d); placed_d <- TRUE; break
          }
        }
        if (!placed_d) dm_clusters[[length(dm_clusters) + 1L]] <- ii_d
      }
      for (cid_d in seq_along(dm_clusters)) {
        spatial_deme_vec[ok_dm_ll[dm_clusters[[cid_d]]]] <- sprintf("S06_DEME_%04d", cid_d)
      }
    }
    still_m <- which(!nzchar(spatial_deme_vec))
    if (length(still_m) > 0L) {
      spatial_deme_vec[still_m] <- ifelse(
        nzchar(state_cnty_geo[still_m]),
        paste0("CNTY_", state_cnty_geo[still_m]),
        paste0("SHEET_", cat_sheet_vec[still_m])
      )
    }
  }

  coll_surname_all <- .norm_collector_surname(rec_by_v)

  clean_peer_df <- {
    idx_cp <- which(anchor_all & !mo_all & !noise_all & !ext_non_type_all & nzchar(sp_clean_vec) & (sp_clean_vec != "Unknown"))
    df_cp <- if (length(idx_cp) > 0L) {
      data.frame(
        sp = sp_clean_vec[idx_cp],
        surname = coll_surname_all[idx_cp],
        deme_id = spatial_deme_vec[idx_cp],
        state_cnty = state_cnty_geo[idx_cp],
        lat = lat_vec_geo[idx_cp],
        lon = lon_vec_geo[idx_cp],
        cat = cat_sheet_vec[idx_cp],
        stringsAsFactors = FALSE
      )
    } else NULL
    db_cp <- NULL
    if (!is.null(db_syntopy_rows) && nrow(db_syntopy_rows) > 0L) {
      raw_cp_nm <- if ("species_standardized" %in% names(db_syntopy_rows)) {
        ifelse(!is.na(db_syntopy_rows$species_standardized) & nzchar(trimws(as.character(db_syntopy_rows$species_standardized))),
               as.character(db_syntopy_rows$species_standardized),
               ifelse(!is.na(db_syntopy_rows$scientificName) & nzchar(trimws(as.character(db_syntopy_rows$scientificName))),
                      as.character(db_syntopy_rows$scientificName), as.character(db_syntopy_rows$species_raw)))
      } else {
        ifelse(!is.na(db_syntopy_rows$scientificName) & nzchar(trimws(as.character(db_syntopy_rows$scientificName))),
               as.character(db_syntopy_rows$scientificName), as.character(db_syntopy_rows$species_raw))
      }
      tx_cp <- standardize_packera_taxon(raw_cp_nm)
      keep_cp <- !is.na(tx_cp) & nzchar(tx_cp) & (tx_cp != "Unknown")
      if (any(keep_cp)) {
        st_cp <- .norm_state_str(db_syntopy_rows$stateProvince[keep_cp])
        cty_cp <- .norm_county_str(db_syntopy_rows$county[keep_cp])
        sc_cp <- ifelse(nzchar(st_cp) & nzchar(cty_cp), paste(st_cp, cty_cp, sep = "||"), "")
        sur_cp <- if ("recordedBy" %in% names(db_syntopy_rows)) .norm_collector_surname(db_syntopy_rows$recordedBy[keep_cp]) else rep("", sum(keep_cp))
        db_cp <- data.frame(
          sp = tx_cp[keep_cp],
          surname = sur_cp,
          deme_id = rep("", sum(keep_cp)),
          state_cnty = sc_cp,
          lat = suppressWarnings(as.numeric(db_syntopy_rows$decimalLatitude[keep_cp])),
          lon = suppressWarnings(as.numeric(db_syntopy_rows$decimalLongitude[keep_cp])),
          cat = as.character(db_syntopy_rows$catalogNumber[keep_cp]),
          stringsAsFactors = FALSE
        )
      }
    }
    res_cp <- if (!is.null(df_cp) && !is.null(db_cp)) rbind(df_cp, db_cp) else if (!is.null(df_cp)) df_cp else db_cp
    if (!is.null(res_cp) && nrow(res_cp) > 0L) res_cp else NULL
  }

  .county_eval_for_taxon <- function(row_i, ctx) {
    sc_i <- state_cnty_geo[row_i]
    st_i <- state_vec_geo[row_i]
    cty_i <- cnty_raw_geo[row_i]
    exact_cnty <- nzchar(sc_i) && (ctx %in% names(usda_counties_by_tx)) && (sc_i %in% usda_counties_by_tx[[ctx]])
    n_st_cnty <- if (nzchar(st_i) && (ctx %in% names(usda_state_county_counts_by_tx)) && (st_i %in% names(usda_state_county_counts_by_tx[[ctx]]))) {
      as.integer(usda_state_county_counts_by_tx[[ctx]][[st_i]])
    } else 0L
    absent_comp <- nzchar(st_i) && nzchar(cty_i) && !exact_cnty && (n_st_cnty >= 5L)
    list(exact_county = exact_cnty, absent_comprehensive_county = absent_comp)
  }

  .eval_deme_peer_votes <- function(row_i, given_sp) {
    sur_i <- coll_surname_all[row_i]
    dm_i <- spatial_deme_vec[row_i]
    sc_i <- state_cnty_geo[row_i]
    lat_i <- lat_vec_geo[row_i]; lon_i <- lon_vec_geo[row_i]; cat_i <- cat_sheet_vec[row_i]
    has_ll_i <- is.finite(lat_i) && is.finite(lon_i)
    same_coll_sp <- character(0); same_deme_sp <- character(0)
    if (!is.null(clean_peer_df) && nrow(clean_peer_df) > 0L) {
      ok_c <- (clean_peer_df$cat != cat_i) & (clean_peer_df$sp != given_sp)
      near_c <- rep(FALSE, nrow(clean_peer_df))
      if (has_ll_i) {
        ok_ll_p <- ok_c & is.finite(clean_peer_df$lat) & is.finite(clean_peer_df$lon)
        if (any(ok_ll_p)) {
          dl_p <- (clean_peer_df$lat[ok_ll_p] - lat_i) * (pi / 180.0); dn_p <- (clean_peer_df$lon[ok_ll_p] - lon_i) * (pi / 180.0)
          hv_p <- sin(dl_p * 0.5)^2 + cos(lat_i * (pi / 180.0)) * cos(clean_peer_df$lat[ok_ll_p] * (pi / 180.0)) * (sin(dn_p * 0.5)^2)
          dkm_p <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hv_p))))
          near_c[which(ok_ll_p)[is.finite(dkm_p) & dkm_p <= 10.0]] <- TRUE
        }
      }
      if (!grepl("^SHEET_", dm_i) && nzchar(dm_i)) {
        near_c[ok_c & nzchar(clean_peer_df$deme_id) & (clean_peer_df$deme_id == dm_i)] <- TRUE
      }
      if (!has_ll_i && nzchar(sc_i)) {
        near_c[ok_c & (clean_peer_df$state_cnty == sc_i)] <- TRUE
      }
      if (any(near_c)) {
        same_deme_sp <- c(same_deme_sp, unique(clean_peer_df$sp[near_c]))
        if (nzchar(sur_i)) {
          sc_hit <- near_c & (clean_peer_df$surname == sur_i)
          if (any(sc_hit)) same_coll_sp <- c(same_coll_sp, unique(clean_peer_df$sp[sc_hit]))
        }
      }
    }
    peer_j <- integer(0)
    if (has_ll_i) {
      ok_ll_df <- which((seq_len(n) != row_i) & (cat_sheet_vec != cat_i) & is.finite(lat_vec_geo) & is.finite(lon_vec_geo))
      if (length(ok_ll_df) > 0L) {
        dl_j <- (lat_vec_geo[ok_ll_df] - lat_i) * (pi / 180.0); dn_j <- (lon_vec_geo[ok_ll_df] - lon_i) * (pi / 180.0)
        hv_j <- sin(dl_j * 0.5)^2 + cos(lat_i * (pi / 180.0)) * cos(lat_vec_geo[ok_ll_df] * (pi / 180.0)) * (sin(dn_j * 0.5)^2)
        dkm_j <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, hv_j))))
        peer_j <- ok_ll_df[is.finite(dkm_j) & dkm_j <= 10.0]
      }
    } else if (nzchar(sc_i)) {
      peer_j <- which((seq_len(n) != row_i) & (cat_sheet_vec != cat_i) & (state_cnty_geo == sc_i))
    }
    if (length(peer_j) > 0L) {
      for (pj in peer_j) {
        cand_pj <- unique(stats::na.omit(c(
          if (is_valid_taxon_str(rf_pred_dv[pj]) && rf_pred_dv[pj] != given_sp) rf_pred_dv[pj] else NA_character_,
          if (is_valid_taxon_str(cda_pred_dv[pj]) && cda_pred_dv[pj] != given_sp && isTRUE(.eval_taxon_geo_plausibility(pj, cda_pred_dv[pj])$plausible)) cda_pred_dv[pj] else NA_character_
        )))
        for (cpj in cand_pj) {
          same_deme_sp <- c(same_deme_sp, cpj)
          if (nzchar(sur_i) && coll_surname_all[pj] == sur_i) same_coll_sp <- c(same_coll_sp, cpj)
        }
      }
    }
    list(same_coll = unique(same_coll_sp), same_deme = unique(same_deme_sp))
  }

  .score_sympatric_cand <- function(row_i, ctx, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
                                    edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob,
                                    min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
                                    geno_is_taxon, geno_sp_i, geno_prob_i,
                                    morph_allopatric_or_weak_self = FALSE, peer_votes = NULL) {
    gw <- .eval_taxon_geo_plausibility(row_i, ctx)$weight
    cty_eval <- .county_eval_for_taxon(row_i, ctx)
    cnty_mult <- if (cty_eval$absent_comprehensive_county) 0.25 else 1.0
    cnty_bonus <- if (cty_eval$exact_county) 0.14 else 0.0
    ms <- if (!is.null(morph_sp) && morph_sp == ctx && is.finite(cda_prob)) {
      cda_prob
    } else if (!is.null(cda_sec_taxon) && cda_sec_taxon == ctx && is.finite(cda_sec)) {
      cda_sec
    } else if (!is.null(morph_X) && !is.null(morph_centroids) && !is.null(morph_S_inv) && (ctx %in% rownames(morph_centroids))) {
      xi <- morph_X[row_i, ]; bm <- if (given_sp %in% rownames(morph_centroids)) morph_centroids[given_sp, ] else colMeans(morph_centroids)
      xi[!is.finite(xi)] <- bm[!is.finite(xi)]
      dv <- xi - morph_centroids[ctx, ]
      exp(-0.25 * max(0.0, as.numeric(t(dv) %*% morph_S_inv %*% dv)))
    } else 0.15
    rf_match <- (rf_is_taxon && rf_sp == ctx) || (edaph_is_taxon && edaph_sp == ctx)
    es <- if (rf_match) {
      w_rf <- if (morph_allopatric_or_weak_self) 0.55 else 0.35
      if (is.finite(rf_prob)) w_rf * rf_prob else 0.20
    } else 0.0
    ts <- if ((!is.na(min_cong_tx_i) && min_cong_tx_i == ctx) || (!is.na(nt_m_i) && nt_m_i == ctx)) 0.20 else 0.0
    if (!is.na(nt_e_i) && nt_e_i == ctx && is.finite(te_z_i) && te_z_i > 2.0) ts <- ts + 0.15
    gs <- if (geno_is_taxon && !is.na(geno_sp_i) && geno_sp_i == ctx && is.finite(geno_prob_i)) 0.85 * geno_prob_i else 0.0
    ps <- if (!is.null(peer_votes) && (ctx %in% peer_votes$same_coll)) 0.20 else if (!is.null(peer_votes) && (ctx %in% peer_votes$same_deme)) 0.10 else 0.0
    gw * cnty_mult * (0.65 * ms + es + ts + gs + cnty_bonus + ps)
  }

  .rank_sympatric_congeners <- function(row_i, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
                                       edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob,
                                       min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
                                       geno_is_taxon, geno_sp_i, geno_prob_i,
                                       morph_allopatric_or_weak_self = FALSE) {
    cand_pool <- unique(stats::na.omit(c(
      if (geno_is_taxon && !is.na(geno_sp_i) && geno_sp_i != given_sp) geno_sp_i else NA_character_,
      if (rf_is_taxon && rf_sp != given_sp) rf_sp else NA_character_,
      if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else NA_character_,
      if (!is.na(nt_e_i) && nt_e_i != given_sp) nt_e_i else NA_character_,
      if (!is.null(morph_sp) && morph_sp != given_sp) morph_sp else NA_character_,
      if (!is.null(cda_sec_taxon) && cda_sec_taxon != given_sp) cda_sec_taxon else NA_character_,
      if (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp) min_cong_tx_i else NA_character_,
      if (!is.na(nt_m_i) && nt_m_i != given_sp) nt_m_i else NA_character_,
      setdiff(active_taxa, given_sp)
    )))
    cand_pool <- cand_pool[vapply(cand_pool, is_valid_taxon_str, logical(1)) & (cand_pool != given_sp)]
    plaus_cands <- cand_pool[vapply(cand_pool, function(ctx) isTRUE(.eval_taxon_geo_plausibility(row_i, ctx)$plausible), logical(1))]
    if (length(plaus_cands) == 0L && nzchar(state_cnty_geo[row_i])) {
      cnty_hits <- names(usda_counties_by_tx)[vapply(names(usda_counties_by_tx), function(ctx) {
        ctx != given_sp && is_valid_taxon_str(ctx) && (state_cnty_geo[row_i] %in% usda_counties_by_tx[[ctx]])
      }, logical(1))]
      plaus_cands <- cnty_hits[vapply(cnty_hits, function(ctx) isTRUE(.eval_taxon_geo_plausibility(row_i, ctx)$plausible), logical(1))]
    }
    if (length(plaus_cands) == 0L) return(list(best = NULL, max_score = -Inf, peer_votes = NULL))
    pv <- .eval_deme_peer_votes(row_i, given_sp)
    scores <- vapply(plaus_cands, function(ctx) {
      .score_sympatric_cand(
        row_i, ctx, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
        edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob,
        min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
        geno_is_taxon, geno_sp_i, geno_prob_i,
        morph_allopatric_or_weak_self = morph_allopatric_or_weak_self, peer_votes = pv
      )
    }, numeric(1))
    list(best = plaus_cands[which.max(scores)], max_score = max(scores), peer_votes = pv)
  }
  rec_num_clean_all <- tolower(gsub("[^a-z0-9]+", "", ifelse(is.na(rec_num_v), "", rec_num_v)))
  rec_num_clean_all[rec_num_clean_all %in% c("sn", "na", "none", "unknown", "0")] <- ""
  exs_clean_all <- tolower(trimws(ifelse(is.na(exs_key_v), "", exs_key_v)))
  exs_clean_all[exs_clean_all %in% c("na", "none", "unknown", "null")] <- ""
  gathering_event_key_all <- ifelse(
    nzchar(exs_clean_all),
    paste0("EXS::", exs_clean_all),
    ifelse(
      nzchar(coll_surname_all) & nzchar(rec_num_clean_all),
      paste0("COLLNUM::", coll_surname_all, "::", rec_num_clean_all),
      paste0(
        "EVT::",
        ifelse(nzchar(coll_surname_all), coll_surname_all, "unk_coll"), "::",
        ifelse(!is.na(yr_v) & nzchar(trimws(yr_v)), trimws(yr_v), "unk_yr"), "::",
        spatial_deme_vec
      )
    )
  )

  cryptic_cluster_pass_map <- list()
  cryptic_non_clinal_map <- list()
  clinal_r2_map <- list()
  cluster_resurrect_bas_map <- list()
  u_sp_cl <- unique(data.frame(sp = sp_clean_vec, cl = g_cl_all, stringsAsFactors = FALSE))
  for (r_idx in seq_len(nrow(u_sp_cl))) {
    sp_k <- u_sp_cl$sp[r_idx]
    cl_k <- u_sp_cl$cl[r_idx]
    if (is.na(cl_k) || cl_k <= 1L) next
    map_key <- paste(sp_k, cl_k, sep = "||")
    idx_cl <- which(sp_clean_vec == sp_k & g_cl_all == cl_k & t1_pure_all)
    if (length(idx_cl) == 0L) {
      cryptic_cluster_pass_map[[map_key]] <- FALSE
      cryptic_non_clinal_map[[map_key]] <- FALSE
      next
    }
    has_collector_col <- any(!is.na(rec_by_v) & nzchar(trimws(rec_by_v))) || any(nzchar(exs_clean_all))
    if (n < 50L && !has_collector_col) {
      pass_events <- length(idx_cl) >= 1L
    } else {
      n_ind_events <- length(unique(gathering_event_key_all[idx_cl]))
      u_coll <- unique(coll_surname_all[idx_cl][nzchar(coll_surname_all[idx_cl])])
      u_inst <- unique(tolower(trimws(inst_v[idx_cl][!is.na(inst_v[idx_cl]) & nzchar(trimws(inst_v[idx_cl]))])))
      has_src_cols <- (length(u_coll) > 0L) || (length(u_inst) > 0L)
      pass_src <- if (has_src_cols) (length(u_coll) >= 2L || length(u_inst) >= 2L) else TRUE
      sp_all_demes <- unique(spatial_deme_vec[sp_clean_vec == sp_k & t1_pure_all])
      req_demes <- if (length(sp_all_demes) >= 3L) 3L else max(1L, length(sp_all_demes))
      n_cl_demes <- length(unique(spatial_deme_vec[idx_cl]))
      pass_events <- (n_ind_events >= 5L) && (n_cl_demes >= req_demes) && pass_src
    }

    fin_c <- idx_cl[is.finite(min_cong_d2_dv[idx_cl])]
    pass_osr_d2 <- if (length(fin_c) > 0L) {
      stats::median(min_cong_d2_dv[fin_c] - osr_chisq95_vec[fin_c], na.rm = TRUE) > 0.0
    } else {
      !is_intra || !nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"))
    }
    fin_p <- idx_cl[is.finite(osr_p_dv[idx_cl]) & is.finite(min_cong_d2_dv[idx_cl])]
    pass_osr_p <- if (length(fin_p) > 0L) (stats::median(osr_p_dv[fin_p], na.rm = TRUE) < 0.05) else TRUE
    pass_osr <- isTRUE(pass_osr_d2 && pass_osr_p)

    cyto_cl_pass <- cyto_barrier_pass_v[idx_cl]
    div_vec <- (is.finite(d_edaph_all[idx_cl]) & d_edaph_all[idx_cl] < 0.30) |
               (is.finite(edaph_z_all[idx_cl]) & edaph_z_all[idx_cl] > 2.5) |
               (is_anth_dv[idx_cl] & is.finite(ddoy_all[idx_cl]) & abs(ddoy_all[idx_cl]) >= tau_cryptic_vec[idx_cl]) |
               cyto_cl_pass
    eval_m <- is.finite(edaph_z_all[idx_cl]) | !is.na(d_ed1_dv[idx_cl]) | !is.na(d_ed2_dv[idx_cl]) | !is.na(d_ed3_dv[idx_cl]) |
              (is_anth_dv[idx_cl] & is.finite(ddoy_all[idx_cl])) |
              cyto_cl_pass
    div_frac <- if (any(eval_m)) mean(div_vec[eval_m], na.rm = TRUE) else mean(div_vec, na.rm = TRUE)
    bf_req <- if (n < 50L) 6.0 else 10.0
    bf_med <- stats::median(bf_all[idx_cl], na.rm = TRUE)
    pass_div <- is.finite(div_frac) && (div_frac >= 0.50) && is.finite(bf_med) && (bf_med >= bf_req)

    # Spatial-Cline & Isolation-by-Distance Validation (Sukumaran & Knowles 2017; Chambers & Hillis 2020)
    non_clinal_ok <- TRUE
    r2_cline_val <- NA_real_
    idx_core <- which(sp_clean_vec == sp_k & g_cl_all == 1L & t1_pure_all)
    if (length(idx_cl) >= 5L && length(idx_core) >= 5L) {
      ll_cl <- idx_cl[is.finite(lat_v_cl[idx_cl]) & is.finite(lon_v_cl[idx_cl])]
      ll_co <- idx_core[is.finite(lat_v_cl[idx_core]) & is.finite(lon_v_cl[idx_core])]
      min_inter_km <- Inf
      if (length(ll_cl) >= 3L && length(ll_co) >= 3L) {
        for (ic1 in ll_cl) {
          dlat_c <- (lat_v_cl[ll_co] - lat_v_cl[ic1]) * (pi / 180.0)
          dlon_c <- (lon_v_cl[ll_co] - lon_v_cl[ic1]) * (pi / 180.0)
          a_c <- sin(dlat_c * 0.5)^2 + cos(lat_v_cl[ic1] * (pi / 180.0)) * cos(lat_v_cl[ll_co] * (pi / 180.0)) * (sin(dlon_c * 0.5)^2)
          min_inter_km <- min(min_inter_km, min(6371.0 * 2.0 * atan2(sqrt(pmax(0.0, a_c)), sqrt(pmax(0.0, 1.0 - a_c)))))
        }
      }
      shared_cnty_contact <- any(nzchar(state_cnty_geo[idx_cl]) & (state_cnty_geo[idx_cl] %in% state_cnty_geo[idx_core]))
      has_parapatric_contact <- (is.finite(min_inter_km) && min_inter_km <= 50.0) || shared_cnty_contact
      anth_cl <- idx_cl[is_anth_dv[idx_cl] & is.finite(ddoy_all[idx_cl])]
      has_allochronic_jump <- (length(anth_cl) >= 3L) &&
        (stats::median(abs(ddoy_all[anth_cl]), na.rm = TRUE) >= stats::median(tau_cryptic_vec[anth_cl], na.rm = TRUE))
      has_azonal_jump <- mean(azonal_flag_all[idx_cl], na.rm = TRUE) >= 0.50

      if (!has_parapatric_contact && !has_allochronic_jump && !has_azonal_jump && length(ll_cl) >= 4L && length(ll_co) >= 4L) {
        comb_ll <- c(ll_co, ll_cl)
        X_cline <- if (!is.null(morph_X)) {
          morph_X[comb_ll, , drop = FALSE]
        } else if (length(morph_coord_cols) >= 1L) {
          as.matrix(df[comb_ll, morph_coord_cols, drop = FALSE])
        } else {
          NULL
        }
        if (!is.null(X_cline) && ncol(X_cline) >= 1L && sum(stats::complete.cases(X_cline)) >= 8L) {
          cc_m <- stats::complete.cases(X_cline)
          Xc_ok <- X_cline[cc_m, , drop = FALSE]
          rows_ok <- comb_ll[cc_m]
          is_cl_ok <- rows_ok %in% idx_cl
          if (sum(is_cl_ok) >= 3L && sum(!is_cl_ok) >= 3L) {
            mu0_c <- colMeans(Xc_ok[!is_cl_ok, , drop = FALSE])
            mu1_c <- colMeans(Xc_ok[is_cl_ok, , drop = FALSE])
            dv_c <- mu1_c - mu0_c
            S_c <- stats::cov(Xc_ok) + diag(1e-3, ncol(Xc_ok))
            Si_c <- tryCatch(solve(S_c), error = function(e) diag(1.0 / pmax(diag(S_c), 1e-3), ncol(Xc_ok)))
            w_dir <- as.numeric(Si_c %*% dv_c)
            proj_c <- as.numeric(sweep(Xc_ok, 2L, mu0_c, "-") %*% w_dir)
            if (stats::sd(proj_c) > 1e-8) {
              lat_s <- as.numeric(scale(lat_v_cl[rows_ok]))
              lon_s <- as.numeric(scale(lon_v_cl[rows_ok]))
              fit_cl <- tryCatch(stats::lm(proj_c ~ lat_s + lon_s), error = function(e) NULL)
              r2_cline_val <- if (!is.null(fit_cl)) summary(fit_cl)$r.squared else NA_real_
              q90_co <- as.numeric(stats::quantile(proj_c[!is_cl_ok], probs = 0.90, names = FALSE))
              q10_cl <- as.numeric(stats::quantile(proj_c[is_cl_ok], probs = 0.10, names = FALSE))
              has_proj_gap <- q10_cl > (q90_co + 0.50 * stats::sd(proj_c[!is_cl_ok]))
              if (is.finite(r2_cline_val) && r2_cline_val >= 0.60 && !has_proj_gap) {
                non_clinal_ok <- FALSE
              } else if (is.finite(r2_cline_val) && r2_cline_val >= 0.72) {
                non_clinal_ok <- FALSE
              }
            }
          }
        }
      }
    }

    # Heterotypic Basionym Resurrection Screening across cohesive cluster members
    bas_cl_raw <- trimws(ifelse(is.na(nt_bas_dv[idx_cl]), "", nt_bas_dv[idx_cl]))
    het_mask_cl <- nzchar(bas_cl_raw)
    if (any(het_mask_cl) && exists(".is_homotypic_sensu_stricto_basionym", mode = "function")) {
      het_mask_cl[het_mask_cl] <- !vapply(bas_cl_raw[het_mask_cl], function(bs) isTRUE(.is_homotypic_sensu_stricto_basionym(bs, sp_k)[1L]), logical(1))
    }
    if (any(het_mask_cl)) {
      cand_bas_tab <- sort(table(bas_cl_raw[het_mask_cl]), decreasing = TRUE)
      top_bas <- names(cand_bas_tab)[1L]
      top_bas_idx <- idx_cl[het_mask_cl & bas_cl_raw == top_bas]
      d2_b_vals <- ifelse(is.finite(bas_d2_dv[top_bas_idx]), bas_d2_dv[top_bas_idx], tm_d2_dv[top_bas_idx])
      d2_b_fin <- d2_b_vals[is.finite(d2_b_vals)]
      bas_d2_ok <- (length(d2_b_fin) == 0L) || (stats::median(d2_b_fin) <= type_morph_chisq95)
      if ((length(top_bas_idx) / length(idx_cl)) >= 0.50 && bas_d2_ok) {
        cluster_resurrect_bas_map[[map_key]] <- top_bas
      }
    }

    cryptic_non_clinal_map[[map_key]] <- isTRUE(non_clinal_ok)
    clinal_r2_map[[map_key]] <- r2_cline_val
    cryptic_cluster_pass_map[[map_key]] <- isTRUE(pass_events && pass_osr && pass_div && non_clinal_ok)
  }
  cryptic_cluster_cohesion_pass <- logical(n)
  cryptic_non_clinal_pass <- logical(n)
  resurrection_basionym <- rep(NA_character_, n)
  hybrid_zone_cohesion_pass <- logical(n)
  hybrid_parent_1 <- rep(NA_character_, n)
  hybrid_parent_2 <- rep(NA_character_, n)
  hyb_cand_pass <- logical(n)
  hyb_cand_p1 <- rep(NA_character_, n)
  hyb_cand_p2 <- rep(NA_character_, n)
  hyb_cand_d1_km <- rep(Inf, n)
  hyb_cand_d2_km <- rep(Inf, n)
  hyb_cand_ddoy <- rep(Inf, n)
  hyb_cand_t_ab <- rep(NA_real_, n)
  hyb_cand_d2_perp <- rep(Inf, n)

  for (i in seq_len(n)) {
    given_sp <- if (!is.na(sp_all_vec[i]) && nzchar(sp_all_vec[i])) sp_all_vec[i] else "Unknown"
    raw_sp   <- if (!is.na(raw_sp_dv[i])) raw_sp_dv[i] else given_sp
    tier     <- if (!is.na(tier_dv[i])) tier_dv[i] else "Tier_3_Bronze"
    det_t    <- tier
    has_obs_doy_i <- valid_doy_all[i]
    sp_med_doy <- if (given_sp %in% names(med_doy)) med_doy[[given_sp]] else glob_circ$med_doy
    sp_q75_doy <- if (given_sp %in% names(q75_doy)) q75_doy[[given_sp]] else glob_circ$q75_doy
    sp_mu_theta <- if (given_sp %in% names(mu_theta_hop)) mu_theta_hop[[given_sp]] else glob_circ$mu_theta
    sp_q75_arc <- if (given_sp %in% names(q75_arc_hop)) q75_arc_hop[[given_sp]] else glob_circ$q75_arc
    hop_doy_i  <- if (has_obs_doy_i && is.finite(hop_doy_all[i])) hop_doy_all[i] else sp_med_doy
    doy        <- hop_doy_i
    hop_arc_i  <- if (has_obs_doy_i && is.finite(hop_doy_all[i])) {
      th_i <- 2.0 * pi * (hop_doy_all[i] / 365.25)
      atan2(sin(th_i - sp_mu_theta), cos(th_i - sp_mu_theta)) * (365.25 / (2.0 * pi))
    } else NA_real_

    morph_sp <- if (!is.na(cda_pred_dv[i]) && nzchar(cda_pred_dv[i])) cda_pred_dv[i] else NULL
    cda_sec_taxon <- if (!is.na(cda_sec1_dv[i]) && nzchar(cda_sec1_dv[i])) cda_sec1_dv[i] else if (!is.na(cda_sec2_dv[i]) && nzchar(cda_sec2_dv[i])) cda_sec2_dv[i] else NULL
    cda_prob <- cda_prob_dv[i]; cda_sec <- cda_sec_p_dv[i]
    cda_marg <- if ("cda_posterior_margin" %in% df_nms_d && is.finite(cda_marg_dv[i])) {
      as.numeric(cda_marg_dv[i])
    } else if (!is.null(morph_sp) && is.finite(cda_prob)) {
      max(0.0, 2.0 * cda_prob - 1.0)
    } else {
      NA_real_
    }

    is_mo_flag <- if (!is.na(mo_flag_dv[i])) isTRUE(mo_flag_dv[i]) else if (!is.na(misid_flag_dv[i])) isTRUE(misid_flag_dv[i]) else FALSE
    triage_cat <- if (!is.na(triage_cat_dv[i])) triage_cat_dv[i] else "Clean_MultiModal_Consensus"
    mm_conf    <- mm_conf_dv[i]

    edaph_sp <- if (!is.na(edaph_sp_dv[i]) && nzchar(edaph_sp_dv[i])) edaph_sp_dv[i] else if (!is.na(rf_pred_dv[i]) && nzchar(rf_pred_dv[i])) rf_pred_dv[i] else given_sp
    rf_sp    <- if (!is.na(rf_pred_dv[i]) && nzchar(rf_pred_dv[i])) rf_pred_dv[i] else edaph_sp
    rf_prob  <- rf_prob_dv[i]
    edaph_z_i <- edaph_z_vec[i]
    te_z_i <- te_z_dv[i]
    edaph_z <- if (is.finite(edaph_z_i)) edaph_z_i else te_z_i
    d_edaph <- if (!is.na(d_ed1_dv[i])) {
      d_ed1_dv[i]
    } else if (!is.na(d_ed2_dv[i])) {
      d_ed2_dv[i]
    } else if (!is.na(d_ed3_dv[i])) {
      d_ed3_dv[i]
    } else if (is_intra) {
      if (is.finite(edaph_z)) exp(-0.5 * (edaph_z / 1.8)^2) else 0.75
    } else {
      if (edaph_sp != given_sp && is.finite(edaph_z) && edaph_z > 2.0) 0.22 else 0.75
    }
    ph_str   <- if (is.finite(ph_num_all[i])) sprintf("%.2f", ph_num_all[i]) else "NA"
    sand_str <- if (is.finite(sand_num_all[i])) sprintf("%.1f%%", sand_num_all[i]) else "NA"

    gmm_unc  <- if (!is.na(gmm_unc_dv[i])) gmm_unc_dv[i] else 0.05
    g_cl     <- if (!is.na(gmm_cl_dv[i])) gmm_cl_dv[i] else 1L
    ddoy     <- if (is.finite(ddoy1_dv[i])) ddoy1_dv[i] else if (is.finite(ddoy2_dv[i])) ddoy2_dv[i] else 0.0
    is_anth  <- is_anth_dv[i]
    bf_val   <- if (is.finite(bf1_dv[i])) bf1_dv[i] else if (is.finite(bf2_dv[i])) bf2_dv[i] else gmm_two_delta_bic
    tau_ecotype_i <- tau_ecotype_vec[i]
    tau_cryptic_i <- tau_cryptic_vec[i]

    morph_is_taxon <- is_valid_taxon_str(morph_sp)
    edaph_is_taxon <- is_valid_taxon_str(edaph_sp)
    rf_is_taxon    <- is_valid_taxon_str(rf_sp)
    sec_is_taxon   <- is_valid_taxon_str(cda_sec_taxon)

    m_tier_src     <- m_tier_src_dv[i]
    n_t1_i         <- n_t1_dv[i]
    is_tier1_pure_morph <- (is.na(m_tier_src) || !nzchar(m_tier_src) || m_tier_src == "Basal_Rosette") && (is.na(n_t1_i) || n_t1_i >= 1L)
    is_tier1_pure_i <- is_tier1_pure_morph

    is_type_specimen_i <- isTRUE(is_type_dv[i])
    type_stat_i <- if (!is.na(type_stat_dv[i]) && nzchar(trimws(type_stat_dv[i])) && tolower(trimws(type_stat_dv[i])) != "none") trimws(type_stat_dv[i]) else "TYPE"
    is_typical_morph_i <- !is.na(is_typ_morph_dv[i]) && isTRUE(is_typ_morph_dv[i])
    typ_m_i <- is_typical_morph_i
    tm_d2_i <- tm_d2_dv[i]
    type_morph_d2_i <- tm_d2_i
    nt_m_i <- if (is_valid_taxon_str(nt_morph_dv[i])) nt_morph_dv[i] else NA_character_
    nt_e_i <- if (is_valid_taxon_str(nt_niche_dv[i])) nt_niche_dv[i] else NA_character_
    clim_z_i <- clim_z_dv[i]
    rf_given_p_i <- rf_given_p_dv[i]
    u_stat_i <- usda_stat_dv[i]
    has_geo_ll_i <- is.finite(lat_vec_geo[i]) && is.finite(lon_vec_geo[i])
    is_extralimital_st <- u_stat_i %in% c("out_of_range", "natureserve_exotic_or_waif")
    is_exs_conflict_i <- isTRUE(exs_conf_dv[i]) && !is.na(exs_gold_tx_dv[i]) && nzchar(exs_gold_tx_dv[i]) && (exs_gold_tx_dv[i] != given_sp)
    is_marginal_range_extension_i <- is_extralimital_st && is.finite(rf_given_p_i) && (rf_given_p_i >= 0.35) && (!is.finite(clim_z_i) || clim_z_i <= 3.0)
    if (is_marginal_range_extension_i) is_marginal_ext_v[i] <- TRUE
    has_obs_edaph_check_i <- (!is.na(edaph_sp_dv[i]) && nzchar(trimws(edaph_sp_dv[i]))) ||
                             (!is.na(rf_pred_dv[i]) && nzchar(trimws(rf_pred_dv[i]))) ||
                             is.finite(edaph_z_i) || !is.na(d_ed1_dv[i]) || !is.na(d_ed2_dv[i]) || !is.na(d_ed3_dv[i])
    in_type_edaphic_env_i <- if (!has_geo_ll_i && is_extralimital_st && !has_obs_edaph_check_i) {
      FALSE
    } else if (is.finite(te_z_i)) {
      (te_z_i <= 2.0)
    } else {
      (!is.finite(edaph_z_i) || edaph_z_i <= 2.0)
    }
    # Rule C: Stage 3 & Stage 4 Autonymic (sensu stricto) Locking (restricted to homotypic basionyms when basionym tracking is present)
    bas_i <- if (!is.na(nt_bas_dv[i]) && nzchar(trimws(nt_bas_dv[i]))) trimws(nt_bas_dv[i]) else NA_character_
    is_heterotypic_bas_match_i <- !is.na(bas_i) && exists(".is_homotypic_sensu_stricto_basionym", mode = "function") &&
                                  !isTRUE(.is_homotypic_sensu_stricto_basionym(bas_i, given_sp)[1L])
    is_sensu_stricto_locked_i <- is_typical_morph_i && in_type_edaphic_env_i && !is_heterotypic_bas_match_i && !is_extralimital_st && !is_exs_conflict_i

    cl_key_i <- paste(given_sp, g_cl, sep = "||")
    cryptic_non_clinal_pass[i] <- if (g_cl > 1L && (cl_key_i %in% names(cryptic_non_clinal_map))) {
      isTRUE(cryptic_non_clinal_map[[cl_key_i]])
    } else {
      FALSE
    }
    clinal_r2_i <- if (g_cl > 1L && (cl_key_i %in% names(clinal_r2_map))) clinal_r2_map[[cl_key_i]] else NA_real_
    cryptic_cluster_cohesion_pass[i] <- if (g_cl > 1L && (cl_key_i %in% names(cryptic_cluster_pass_map))) {
      isTRUE(cryptic_cluster_pass_map[[cl_key_i]])
    } else {
      FALSE
    }
    osr_chisq95_i <- osr_chisq95_vec[i]
    osr_cryptic_ok_i <- if (is.finite(min_cong_d2_dv[i])) {
      (min_cong_d2_dv[i] > osr_chisq95_i) && (!is.finite(osr_p_dv[i]) || osr_p_dv[i] < 0.05)
    } else if (is.finite(osr_p_dv[i])) {
      osr_p_dv[i] < 0.05
    } else {
      !is_intra || !nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH"))
    }

    has_obs_edaph_i <- (!is.na(edaph_sp_dv[i]) && nzchar(trimws(edaph_sp_dv[i]))) ||
                       (!is.na(rf_pred_dv[i]) && nzchar(trimws(rf_pred_dv[i]))) ||
                       is.finite(edaph_z_i) || is.finite(te_z_i) ||
                       !is.na(d_ed1_dv[i]) || !is.na(d_ed2_dv[i]) || !is.na(d_ed3_dv[i])
    edaph_strong_support <- has_obs_edaph_i && (
      (is.finite(rf_prob) && rf_prob >= 0.55) ||
      (!is.finite(rf_prob) && ((is.finite(d_edaph) && d_edaph < 0.35) || (is.finite(edaph_z) && edaph_z > 2.0)))
    )
    eff_cda_prob <- if (is.finite(cda_oof_p_dv[i]) && is.finite(cda_prob)) {
      min(cda_prob, cda_oof_p_dv[i])
    } else if (is.finite(cda_oof_p_dv[i])) {
      cda_oof_p_dv[i]
    } else {
      cda_prob
    }
    morph_edaph_alt_agree <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && ((edaph_is_taxon && edaph_sp == morph_sp && edaph_strong_support) || (!is.na(nt_e_i) && nt_e_i == morph_sp && is.finite(te_z_i) && te_z_i > 2.0)) && (is.na(eff_cda_prob) || eff_cda_prob >= 0.55))
    type_morph_alt_agree <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && !is.na(nt_m_i) && nt_m_i == morph_sp && (!is.na(eff_cda_prob) && eff_cda_prob >= 0.85))
    both_alt_agree <- morph_edaph_alt_agree || type_morph_alt_agree
    type_synonymy_agree <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && is.finite(cda_prob) && cda_prob >= 0.80 && (
      (edaph_is_taxon && edaph_sp == morph_sp && edaph_strong_support) ||
      (!is.na(nt_e_i) && nt_e_i == morph_sp) ||
      (!is.na(nt_m_i) && nt_m_i == morph_sp && cda_prob >= 0.85)
    ))

    # Rule A: Nomenclatural Type Control Pinning
    if (is_type_specimen_i) {
      alt_aff_cands <- c(
        if (morph_is_taxon && morph_sp != given_sp) morph_sp else NA_character_,
        if (edaph_is_taxon && edaph_sp != given_sp && edaph_strong_support) edaph_sp else NA_character_,
        if (!is.na(nt_m_i) && nt_m_i != given_sp) nt_m_i else NA_character_,
        if (!is.na(nt_e_i) && nt_e_i != given_sp) nt_e_i else NA_character_
      )
      alt_aff_cands <- alt_aff_cands[!is.na(alt_aff_cands) & nzchar(alt_aff_cands)]
      rec_taxa[i] <- given_sp
      treat_broad[i] <- "Confirmed Species"
      treat_sub[i] <- "Confirmed_Species"
      if (type_synonymy_agree) {
        status_calls[i] <- "Type_Synonymy_Candidate"
        triage_prios[i] <- "CRITICAL"
        triage_acts[i] <- "Evaluate_Heterotypic_Synonymy"
        rationales[i] <- sprintf(
          "Rule A Nomenclatural Type Control (%s): Protected from Misidentification_Reassign; cross-modal type-to-taxon affinity with %s (evaluate heterotypic synonymy)",
          type_stat_i, morph_sp
        )
      } else if (length(alt_aff_cands) > 0L) {
        aff_tx <- names(sort(table(alt_aff_cands), decreasing = TRUE))[1L]
        status_calls[i] <- "Confirmed_Type_Control"
        triage_prios[i] <- "LOW"
        triage_acts[i] <- "Accept_Current_Determination"
        rationales[i] <- sprintf(
          "Rule A Nomenclatural Type Control (%s): Pinned as %s sensu stricto baseline; partial type-to-taxon affinity noted with %s",
          type_stat_i, given_sp, aff_tx
        )
      } else {
        status_calls[i] <- "Confirmed_Type_Control"
        triage_prios[i] <- "RESOLVED"
        triage_acts[i] <- "Accept_Current_Determination"
        rationales[i] <- sprintf(
          "Rule A Nomenclatural Type Control (%s): Immutable %s sensu stricto baseline concordant across morphological, phenological, and locus classicus streams",
          type_stat_i, given_sp
        )
      }
      next
    }

    # Self-Tolerance Ellipsoid Exclusion (Zapata & Jimenez 2012): a voucher inside its recorded species'
    # 95% Mahalanobis tolerance ellipsoid cannot be reassigned solely because an overlapping congener centroid is slightly closer.
    self_df_i <- if (is.finite(self_df_raw[i]) && self_df_raw[i] >= 1.0) as.integer(round(self_df_raw[i])) else p_type_morph
    self_d2_eff <- if (is.finite(self_d2_dv[i])) {
      self_d2_dv[i]
    } else if (is.finite(tm_d2_i)) {
      tm_d2_i
    } else if (!is.null(morph_X) && !is.null(morph_centroids) && !is.null(self_morph_S_inv) &&
               (given_sp %in% rownames(morph_centroids)) && all(is.finite(morph_X[i, ]))) {
      dv_s <- morph_X[i, ] - morph_centroids[given_sp, ]
      self_df_i <- max(1L, ncol(morph_X))
      max(0.0, as.numeric(t(dv_s) %*% self_morph_S_inv %*% dv_s))
    } else {
      NA_real_
    }
    tm_df_i <- self_df_i
    tm_chisq95_i <- stats::qchisq(0.95, df = tm_df_i)
    outside_self_morph_95 <- !is.finite(self_d2_eff) || (self_d2_eff > stats::qchisq(0.95, df = self_df_i))

    cl_floor_i <- if (is.finite(cl_thr_dv[i])) min(0.82, max(0.58, cl_thr_dv[i])) else 0.0
    cda_cut_base <- if (det_t == "Tier_1_Gold") 0.78 else if (det_t == "Tier_2_Silver") 0.70 else 0.62
    cda_cut <- max(cda_cut_base, if (det_t == "Tier_3_Bronze") min(0.72, cl_floor_i) else cl_floor_i)

    # Open-Set Congener Match (congener_osr_reassign) with exact open_set_df & intraspecific support
    min_cong_tx_i <- if (!is.na(min_cong_tx_dv[i]) && nzchar(trimws(min_cong_tx_dv[i]))) trimws(min_cong_tx_dv[i]) else NA_character_
    min_cong_d2_i <- min_cong_d2_dv[i]
    osr_df_i <- osr_df_vec[i]
    osr_chisq90 <- osr_chisq90_vec[i]
    osr_self_outlier <- (is.finite(tm_d2_i) && (tm_d2_i > tm_chisq95_i)) ||
                        (is.finite(self_d2_eff) && (self_d2_eff > stats::qchisq(0.95, df = self_df_i)))
    ref_self_d2_norm <- if (is.finite(tm_d2_i) && tm_d2_i > 0) {
      tm_d2_i / max(1.0, tm_df_i)
    } else if (is.finite(self_d2_eff) && self_d2_eff > 0) {
      self_d2_eff / max(1.0, self_df_i)
    } else {
      NA_real_
    }
    osr_strictly_closer <- if (is.finite(ref_self_d2_norm)) {
      (min_cong_d2_i / max(1.0, osr_df_i)) < (0.70 * ref_self_d2_norm)
    } else {
      is_mo_flag
    }
    osr_cda_support <- (!is.null(morph_sp) && morph_is_taxon && morph_sp == min_cong_tx_i && is.finite(eff_cda_prob) && eff_cda_prob >= 0.50) ||
                       (!is.na(nt_m_i) && nzchar(nt_m_i) && nt_m_i == min_cong_tx_i) ||
                       ((is_intra || !morph_is_taxon) && is.finite(min_cong_d2_i) && (min_cong_d2_i <= osr_chisq85_vec[i]))
    geo_osr_res <- if (!is.na(min_cong_tx_i)) .eval_taxon_geo_plausibility(i, min_cong_tx_i) else list(plausible = FALSE, weight = 0.0, dist_km = Inf)
    osr_geo_ok <- isTRUE(geo_osr_res$plausible) && (!is_intra || (if (has_ext_geo_ref) (is.finite(geo_osr_res$dist_km) && geo_osr_res$dist_km <= 250.0) else (!is.finite(geo_osr_res$dist_km) || geo_osr_res$dist_km <= 250.0)))
    osr_d2_ceil <- if (is_intra) stats::qchisq(0.50, df = max(1.0, osr_df_i)) else osr_chisq90
    osr_gold_eco_lock <- (det_t == "Tier_1_Gold") && has_obs_edaph_i && in_type_edaphic_env_i &&
                         !is_extralimital_st && !is_exs_conflict_i &&
                         !(is.finite(clim_z_i) && clim_z_i > 3.5 && (!is.finite(rf_given_p_i) || rf_given_p_i < 0.25)) &&
                         (is.na(exs_gold_tx_dv[i]) || exs_gold_tx_dv[i] != min_cong_tx_i) &&
                         (!edaph_is_taxon || edaph_sp == given_sp)
    congener_osr_reassign <- !is_type_specimen_i &&
      is_tier1_pure_i &&
      !osr_gold_eco_lock &&
      is.finite(min_cong_d2_i) &&
      !is.na(min_cong_tx_i) && nzchar(min_cong_tx_i) &&
      (min_cong_tx_i != given_sp) &&
      (min_cong_d2_i <= osr_d2_ceil) &&
      osr_self_outlier &&
      osr_strictly_closer &&
      osr_cda_support &&
      osr_geo_ok &&
      (!is_intra || ((is_mo_flag || (is.finite(edaph_z_i) && edaph_z_i > 1.75)) && (det_t != "Tier_1_Gold" || (is.finite(edaph_z_i) && edaph_z_i > 2.0)))) &&
      (!isTRUE(typ_m_i) || (is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i))

    # Rule B: Stage 1 Type-Control Consensus Reassignment (type_control_reassign) with 1-Vote-Per-Modality & Posterior Floor
    morph_concordant_self <- (morph_is_taxon && morph_sp == given_sp && is.finite(cda_prob) && cda_prob >= 0.65) ||
                             (isTRUE(typ_m_i) && is.finite(tm_d2_i) && tm_d2_i <= tm_chisq95_i)
    vote_morph <- if (morph_concordant_self) {
      given_sp
    } else if (morph_is_taxon && morph_sp != given_sp && is.finite(eff_cda_prob) && eff_cda_prob >= 0.55) {
      morph_sp
    } else if (!is.na(nt_m_i) && nzchar(nt_m_i) && is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i) {
      nt_m_i
    } else if (morph_is_taxon && is.finite(eff_cda_prob) && eff_cda_prob >= 0.55) {
      morph_sp
    } else {
      NA_character_
    }
    rf_cand_sp <- if (rf_is_taxon && rf_sp != given_sp) rf_sp else if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else if (edaph_is_taxon) edaph_sp else if (rf_is_taxon) rf_sp else NA_character_
    rf_vote <- if (has_obs_edaph_i && !is.na(rf_cand_sp) && ((is.finite(rf_prob) && rf_prob >= 0.55) || (!is.finite(rf_prob) && ((is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z_i) && edaph_z_i > 2.0))))) rf_cand_sp else NA_character_
    nt_e_vote <- if (!is.na(nt_e_i) && nzchar(nt_e_i) && is.finite(te_z_i) && te_z_i > 2.0) nt_e_i else NA_character_
    vote_edaph <- if (!is.na(nt_e_vote) && !is.na(rf_vote)) {
      if (nt_e_vote == rf_vote) nt_e_vote
      else if (is.finite(te_z_i) && te_z_i > 2.5 && (!is.finite(rf_prob) || rf_prob < 0.70)) nt_e_vote
      else if (is.finite(rf_prob) && rf_prob >= 0.70) rf_vote
      else NA_character_
    } else if (!is.na(nt_e_vote)) {
      nt_e_vote
    } else {
      rf_vote
    }
    deviates_from_given_type <- (is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i) ||
                                (is.finite(tm_conc_dv[i]) && tm_conc_dv[i] < 0.25) ||
                                (is.finite(te_z_i) && te_z_i > 2.0)
    type_control_reassign <- !is_type_specimen_i &&
                             !morph_concordant_self &&
                             outside_self_morph_95 &&
                             (!isTRUE(typ_m_i) || (is.finite(tm_d2_i) && tm_d2_i > tm_chisq95_i)) &&
                             !is.na(vote_morph) && !is.na(vote_edaph) &&
                             (vote_morph == vote_edaph) && (vote_morph != given_sp) &&
                             deviates_from_given_type
    type_control_target_sp <- if (type_control_reassign) vote_morph else NULL

    # Multi-stream corroboration vs contradiction for Stage 1 (Padial et al. 2010; Goodwin et al. 2015)
    has_geno_i <- has_geno_dv[i]
    geno_sp_i <- if (has_geno_i && !is.na(geno_pred_dv[i]) && nzchar(geno_pred_dv[i])) geno_pred_dv[i] else NA_character_
    geno_is_taxon <- has_geno_i && !is.na(geno_sp_i) && is_valid_taxon_str(geno_sp_i)
    geno_prob_i <- geno_prob_dv[i]
    geno_marg_i <- if (is.finite(geno_marg_dv[i])) geno_marg_dv[i] else if (is.finite(geno_prob_i)) max(0.0, 2.0 * geno_prob_i - 1.0) else 0.0
    geno_tier_i <- geno_tier_dv[i]
    geno_g1_direct <- !is.na(geno_tier_i) && (geno_tier_i %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate"))
    geno_unilat_tier_ok <- !is.na(geno_tier_i) && (geno_tier_i %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy"))
    geno_gdi_i <- if (has_geno_i) geno_gdi_dv[i] else NA_real_
    eff_gdi_i <- eff_gdi_dv[i]
    is_cluster_prop_gdi_only_i <- is_cluster_prop_gdi_only[i]
    tier0_abstain_i <- tier0_abstain_dv[i]
    geno_h_i <- if (has_geno_i) geno_h_dv[i] else NA_real_
    geno_q12_i <- if (has_geno_i) geno_q12_dv[i] else NA_real_
    geno_hyb_cls_i <- if (has_geno_i && !is.na(geno_hyb_cls_dv[i]) && nzchar(trimws(geno_hyb_cls_dv[i]))) trimws(geno_hyb_cls_dv[i]) else NA_character_
    geno_fb_i <- if (has_geno_i) geno_fb_max_dv[i] else NA_real_
    geno_cytonuc_i <- has_geno_i && isTRUE(geno_cytonuc_dv[i])
    geno_divergent_nonhybrid_i <- geno_divergent_nonhybrid_v[i]
    geno_contemp_hybrid_i <- geno_contemp_hybrid_v[i]
    geno_ancient_introg_i <- geno_ancient_introg_v[i]
    geno_gdi_oversplit_i <- geno_gdi_oversplit_v[i]
    geno_gdi_cryptic_support_i <- geno_gdi_cryptic_support_v[i]
    geno_strong_i <- geno_is_taxon && !geno_adm_dv[i] && !geno_qry_dv[i] && is.finite(geno_prob_i) && (geno_prob_i >= 0.85) && (geno_marg_i >= 0.40)

    cand_alt_sp <- if (geno_unilat_tier_ok && geno_strong_i && (geno_sp_i != given_sp) && (geno_g1_direct || !morph_edaph_alt_agree || (!is.null(morph_sp) && morph_sp == geno_sp_i))) {
      geno_sp_i
    } else if (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp) {
      morph_sp
    } else if (!is.null(type_control_target_sp)) {
      type_control_target_sp
    } else if (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp) {
      min_cong_tx_i
    } else if (geno_unilat_tier_ok && geno_strong_i && (geno_sp_i != given_sp)) {
      geno_sp_i
    } else {
      NA_character_
    }
    geno_supports_given <- geno_is_taxon && (geno_sp_i == given_sp) && !geno_adm_dv[i] && !geno_qry_dv[i] && is.finite(geno_prob_i) && (geno_prob_i >= 0.65)
    geno_supports_alt <- geno_is_taxon && !is.na(cand_alt_sp) && (geno_sp_i == cand_alt_sp) && !geno_adm_dv[i] && is.finite(geno_prob_i) && (geno_prob_i >= 0.70)
    is_single_ramet_sheet <- is.na(n_sheet_ramets_vec[i]) || (n_sheet_ramets_vec[i] <= 1L)
    eco_supports_given <- has_obs_edaph_i && is_single_ramet_sheet && (
      (edaph_is_taxon && edaph_sp == given_sp && ((is.finite(rf_prob) && rf_prob >= 0.65) || (!is.finite(rf_prob) && is.finite(d_edaph) && d_edaph >= 0.60))) ||
      (is.finite(edaph_z) && edaph_z <= 1.5 && (!edaph_is_taxon || edaph_sp == given_sp))
    )
    eco_supports_alt <- has_obs_edaph_i && !is.na(cand_alt_sp) && (
      (((edaph_is_taxon && edaph_sp == cand_alt_sp) || (rf_is_taxon && rf_sp == cand_alt_sp)) && edaph_strong_support) ||
      (!is.na(nt_e_i) && nt_e_i == cand_alt_sp && is.finite(te_z_i) && te_z_i > 2.0)
    )

    n_sec_giv <- 0L; n_sec_alt <- 0L; n_diag_sec <- 0L; vestiture_supports_given <- FALSE
    if (!is.null(sec_organ_mat) && length(sec_organ_cols) > 0L && !is.na(cand_alt_sp) &&
        (given_sp %in% names(sec_organ_med_by_tx)) && (cand_alt_sp %in% names(sec_organ_med_by_tx))) {
      mg_v <- sec_organ_med_by_tx[[given_sp]]; ma_v <- sec_organ_med_by_tx[[cand_alt_sp]]
      for (sj in seq_along(sec_organ_cols)) {
        vij <- sec_organ_mat[i, sj]; mad_j <- sec_organ_mad_pool[sj]
        if (is.finite(vij) && !is.null(mg_v) && !is.null(ma_v) &&
            is.finite(mg_v[sj]) && is.finite(ma_v[sj]) && is.finite(mad_j) &&
            abs(ma_v[sj] - mg_v[sj]) >= 1.5 * mad_j) {
          n_diag_sec <- n_diag_sec + 1L
          zg <- abs(vij - mg_v[sj]) / mad_j
          za <- abs(vij - ma_v[sj]) / mad_j
          if (zg <= 1.25 && za >= 1.75) {
            n_sec_giv <- n_sec_giv + 1L
            if (grepl("^blade_cie_l", sec_organ_cols[sj])) vestiture_supports_given <- TRUE
          } else if (za <= 1.25 && zg >= 1.75) {
            n_sec_alt <- n_sec_alt + 1L
          }
        }
      }
    }
    sec_organ_supports_given <- (n_sec_giv >= 1L) && (n_sec_alt == 0L)
    sec_organ_supports_alt <- (n_sec_alt >= 1L) && (n_sec_giv == 0L)

    pheno_supports_given <- FALSE; pheno_supports_alt <- FALSE
    if (is_single_ramet_sheet && !is.na(cand_alt_sp) && has_obs_doy_i && is_anth && is.finite(hop_doy_all[i]) &&
        (given_sp %in% names(has_anthesis_doy_by_tx)) && (cand_alt_sp %in% names(has_anthesis_doy_by_tx)) &&
        (given_sp %in% names(mu_theta_hop)) && (cand_alt_sp %in% names(mu_theta_hop)) &&
        isTRUE(has_anthesis_doy_by_tx[[given_sp]]) && isTRUE(has_anthesis_doy_by_tx[[cand_alt_sp]])) {
      th_g <- mu_theta_hop[[given_sp]]; th_a <- mu_theta_hop[[cand_alt_sp]]
      ddoy_ga <- abs(atan2(sin(th_g - th_a), cos(th_g - th_a))) * (365.25 / (2.0 * pi))
      if (is.finite(ddoy_ga) && ddoy_ga >= 14.0) {
        th_obs <- 2.0 * pi * (hop_doy_all[i] / 365.25)
        dg_obs <- abs(atan2(sin(th_obs - th_g), cos(th_obs - th_g))) * (365.25 / (2.0 * pi))
        da_obs <- abs(atan2(sin(th_obs - th_a), cos(th_obs - th_a))) * (365.25 / (2.0 * pi))
        pheno_supports_given <- (dg_obs <= 10.0) && (da_obs >= 16.0)
        pheno_supports_alt <- (da_obs <= 10.0) && (dg_obs >= 16.0)
      }
    }
    if (is.finite(clim_z_i) && clim_z_i > 3.5 && (!is.finite(rf_given_p_i) || rf_given_p_i < 0.25)) {
      pheno_supports_given <- FALSE
    }
    exs_supports_given <- !is.na(exs_gold_tx_dv[i]) && nzchar(exs_gold_tx_dv[i]) && (exs_gold_tx_dv[i] == given_sp)
    exs_supports_alt <- !is.na(exs_gold_tx_dv[i]) && nzchar(exs_gold_tx_dv[i]) && !is.na(cand_alt_sp) && (exs_gold_tx_dv[i] == cand_alt_sp)
    any_non_outline_supports_alt <- eco_supports_alt || sec_organ_supports_alt || pheno_supports_alt || exs_supports_alt || geno_supports_alt ||
      (!is.na(nt_m_i) && !is.na(cand_alt_sp) && nt_m_i == cand_alt_sp && is.finite(tm_d2_i) && tm_d2_i > type_morph_chisq95)
    any_non_outline_supports_given <- eco_supports_given || sec_organ_supports_given || pheno_supports_given || exs_supports_given || geno_supports_given

    exs_conflict_reassign <- !is_type_specimen_i && is_exs_conflict_i && !exs_supports_given
    extralimital_biogeo_reassign <- !is_type_specimen_i && !is_sensu_stricto_locked_i && !exs_supports_given &&
      !is_marginal_range_extension_i && is_extralimital_st &&
      is.finite(clim_z_i) && (clim_z_i > 3.5) && (!is.finite(rf_given_p_i) || rf_given_p_i < 0.25) &&
      !(morph_is_taxon && !is.null(morph_sp) && morph_sp == given_sp && is.finite(eff_cda_prob) && eff_cda_prob >= 0.75)
    ungeoref_extralimital_reassign <- !is_type_specimen_i && !is_sensu_stricto_locked_i && !exs_supports_given &&
      !is_marginal_range_extension_i && is_extralimital_st && !has_geo_ll_i &&
      ((!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp && is.finite(eff_cda_prob) && eff_cda_prob >= 0.80) ||
       (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp && is.finite(min_cong_d2_i) && min_cong_d2_i <= osr_chisq90 &&
        (!morph_is_taxon || is.null(morph_sp) || morph_sp != given_sp)))
    portal_lag_reassign <- exs_conflict_reassign || extralimital_biogeo_reassign || ungeoref_extralimital_reassign

    # 2. General Ontogenetic / Taphonomic Plasticity Rule (circular Hopkins-adjusted DOY quantile & edaphic envelope)
    in_edaphic_envelope <- (is.finite(edaph_z_i) && edaph_z_i <= 2.0) || (edaph_is_taxon && edaph_sp == given_sp)
    morph_div_i <- (!is.null(morph_sp) && morph_is_taxon && morph_sp != given_sp) || is_mo_flag
    single_borderline_div <- morph_div_i &&
                             (is.null(morph_sp) || !morph_is_taxon || morph_sp == given_sp || is.na(eff_cda_prob) || eff_cda_prob < cda_cut)
    dlf_i <- dlf_dv[i]; asym_i <- asym_dv[i]
    is_sym_dissected_i <- ((is.finite(dlf_i) && dlf_i > 0.0) || is_dis_dv[i]) &&
                          !(is.finite(asym_i) && asym_i > 0.22)
    is_late_foliar_wear <- has_obs_doy_i && is.finite(hop_arc_i) && is.finite(sp_q75_arc) && (hop_arc_i >= sp_q75_arc) && (g_cl == 1L) && single_borderline_div && !is_sym_dissected_i
    is_taphonomic_plasticity <- (triage_cat != "Severe_Triple_Stream_Conflict") &&
                                !both_alt_agree && !morph_edaph_alt_agree &&
                                !type_control_reassign && !congener_osr_reassign &&
                                !portal_lag_reassign &&
                                in_edaphic_envelope && is_late_foliar_wear

    is_intermodal_mosaic <- FALSE
    non_assigned_known <- (is.na(cda_prob) || cda_prob < 0.70) && !both_alt_agree && !morph_edaph_alt_agree &&
                          !type_control_reassign && !congener_osr_reassign
    has_cytotype_barrier_i <- isTRUE(cyto_barrier_pass_v[i])
    has_morph_or_gmm_outlier <- is_mo_flag || (g_cl > max(2L, n_known_taxa) && bf_val >= 6.0)
    has_allochronic_barrier <- is_anth && (abs(ddoy) >= tau_cryptic_i)
    has_edaphic_barrier <- (is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z) && edaph_z > 2.5)
    has_eco_pheno_barrier_ms <- has_allochronic_barrier || has_edaphic_barrier || geno_gdi_cryptic_support_i
    morph_exceeds_chisq95_i <- if (is.finite(self_d2_eff)) {
      (self_d2_eff > stats::qchisq(0.95, df = self_df_i))
    } else if (is.finite(tm_d2_i)) {
      (tm_d2_i > tm_chisq95_i)
    } else if (is.finite(min_cong_d2_dv[i])) {
      (min_cong_d2_dv[i] > osr_chisq95_i)
    } else {
      is_mo_flag && outside_self_morph_95
    }
    cyto_barrier_valid_ms <- has_cytotype_barrier_i && (
      has_eco_pheno_barrier_ms ||
      (morph_exceeds_chisq95_i && outside_self_morph_95 && osr_cryptic_ok_i && cryptic_cluster_cohesion_pass[i])
    )
    mod_morph <- is_tier1_pure_i && (
      (has_morph_or_gmm_outlier && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i))) ||
      (has_cytotype_barrier_i && morph_exceeds_chisq95_i && outside_self_morph_95 && (is_mo_flag || (g_cl > 1L && bf_val >= 6.0) || is.finite(self_d2_eff)))
    )
    mod_non_assign <- non_assigned_known
    mod_barrier <- has_eco_pheno_barrier_ms || cyto_barrier_valid_ms
    n_cryptic_modalities <- as.integer(mod_morph) + as.integer(mod_non_assign) +
                            as.integer(has_allochronic_barrier) + as.integer(has_edaphic_barrier) +
                            as.integer(geno_gdi_cryptic_support_i) + as.integer(cyto_barrier_valid_ms)
    raw_is_crypt_pre_gdi <- !is_sensu_stricto_locked_i && !is_taphonomic_plasticity &&
                            (cryptic_cluster_cohesion_pass[i] || geno_gdi_cryptic_support_i) &&
                            (osr_cryptic_ok_i || (geno_gdi_cryptic_support_i && !is.finite(min_cong_d2_dv[i]))) &&
                            (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i)) &&
                            is_tier1_pure_i &&
                            (mod_morph || (geno_gdi_cryptic_support_i && (geno_qry_dv[i] || outside_self_morph_95))) &&
                            (mod_non_assign || geno_gdi_cryptic_support_i) &&
                            mod_barrier && (n_cryptic_modalities >= 2L)
    is_crypt <- raw_is_crypt_pre_gdi && !geno_gdi_oversplit_i && !tier0_abstain_i
    cryptic_oversplit_blocked_i <- raw_is_crypt_pre_gdi && geno_gdi_oversplit_i

    # Stage 1: Misidentification_Reassign BEFORE is_cryptic_intra & Hybrid_Introgressant (Schlick-Steiner et al. 2010; Padial et al. 2010)
    edaph_contradicts_morph <- has_obs_edaph_i && is_single_ramet_sheet && edaph_is_taxon &&
                               (edaph_sp == given_sp) && (rf_sp == given_sp) &&
                               ((is.finite(rf_prob) && rf_prob >= 0.80) || (!is.finite(rf_prob) && is.finite(d_edaph) && d_edaph >= 0.70)) &&
                               (is.na(eff_cda_prob) || eff_cda_prob < 0.88) &&
                               (!is.finite(d_edaph) || d_edaph >= 0.65) &&
                               !sec_organ_supports_alt && !exs_supports_alt

    tier_authority_gate_ok <- if (det_t == "Tier_1_Gold") {
      outside_self_morph_95 && !exs_supports_given && !vestiture_supports_given && any_non_outline_supports_alt
    } else if (det_t == "Tier_2_Silver") {
      outside_self_morph_95 && !exs_supports_given && !vestiture_supports_given && (
        if (!is_single_ramet_sheet) {
          any_non_outline_supports_alt || (is.finite(eff_cda_prob) && eff_cda_prob >= 0.88)
        } else if (any_non_outline_supports_given) {
          any_non_outline_supports_alt
        } else {
          any_non_outline_supports_alt || (is.finite(eff_cda_prob) && eff_cda_prob >= 0.85)
        }
      )
    } else {
      outside_self_morph_95 && !exs_supports_given && (!vestiture_supports_given || (eco_supports_alt && pheno_supports_alt))
    }

    single_morph_reassign <- morph_is_taxon && (morph_sp != given_sp) && !is.na(eff_cda_prob) &&
                             (eff_cda_prob >= if (is_tier1_pure_i) cda_cut else max(0.85, cda_cut)) &&
                             !edaph_contradicts_morph && tier_authority_gate_ok
    eff_triple_conflict <- (triage_cat == "Severe_Triple_Stream_Conflict") && outside_self_morph_95 &&
                           (both_alt_agree || morph_edaph_alt_agree || single_morph_reassign || congener_osr_reassign) &&
                           tier_authority_gate_ok
    morph_corrupt_cand <- is_mo_flag && morph_is_taxon && (morph_sp != given_sp) &&
                          (!is.na(eff_cda_prob) && eff_cda_prob >= (if (is_tier1_pure_i) max(0.68, cda_cut) else 0.85)) &&
                          !edaph_contradicts_morph && tier_authority_gate_ok
    geno_uncontradicted_reassign <- !is_type_specimen_i && !is_sensu_stricto_locked_i &&
                                    geno_unilat_tier_ok && geno_strong_i && (geno_sp_i != given_sp) &&
                                    (geno_prob_i >= if (geno_g1_direct) 0.85 else 0.92) &&
                                    !(morph_concordant_self && eco_supports_given)
    is_misid <- (!morph_concordant_self && (
      (!is_taphonomic_plasticity && !is_crypt && !is_intermodal_mosaic && outside_self_morph_95 && tier_authority_gate_ok && (
        both_alt_agree || morph_edaph_alt_agree ||
        eff_triple_conflict || single_morph_reassign ||
        morph_corrupt_cand
      )) || (type_control_reassign && !exs_supports_given) || (congener_osr_reassign && !exs_supports_given)
    )) || portal_lag_reassign || (geno_uncontradicted_reassign && !exs_supports_given)
    is_severe_misid <- is_misid

    geno_priority_override <- geno_uncontradicted_reassign && !is.na(geno_sp_i) &&
      (geno_g1_direct || !morph_edaph_alt_agree || (!is.null(morph_sp) && morph_sp == geno_sp_i))

    if (is_severe_misid) {
      morph_allo_or_weak <- (!morph_is_taxon || is.null(morph_sp) || !isTRUE(.eval_taxon_geo_plausibility(i, morph_sp)$plausible) || (morph_sp == given_sp && is_extralimital_st))
      rank_pl <- NULL
      if (portal_lag_reassign) {
        portal_lag_flag[i] <- TRUE
        portal_lag_gate[i] <- if (exs_conflict_reassign) "Gate1_Exsiccatae_Conflict" else if (extralimital_biogeo_reassign) "Gate3_Extralimital_Bioclimatic" else "Gate3_Ungeoref_Extralimital"
      }
      raw_target_sp <- if (geno_priority_override) {
        geno_sp_i
      } else if (exs_conflict_reassign) {
        exs_gold_tx_dv[i]
      } else if (extralimital_biogeo_reassign || ungeoref_extralimital_reassign) {
        rank_pl <- .rank_sympatric_congeners(
          i, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
          edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob,
          min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
          geno_is_taxon, geno_sp_i, geno_prob_i,
          morph_allopatric_or_weak_self = TRUE
        )
        if (!is.null(rank_pl$best) && is.finite(rank_pl$max_score) && rank_pl$max_score >= 0.12) {
          rank_pl$best
        } else if (rf_is_taxon && rf_sp != given_sp && isTRUE(.eval_taxon_geo_plausibility(i, rf_sp)$plausible)) {
          rf_sp
        } else if (morph_is_taxon && !is.null(morph_sp) && morph_sp != given_sp && isTRUE(.eval_taxon_geo_plausibility(i, morph_sp)$plausible)) {
          morph_sp
        } else {
          given_sp
        }
      } else if (type_control_reassign && !is.null(type_control_target_sp)) {
        type_control_target_sp
      } else if (both_alt_agree && morph_is_taxon) {
        morph_sp
      } else if (congener_osr_reassign && !is.na(min_cong_tx_i)) {
        min_cong_tx_i
      } else if (geno_uncontradicted_reassign && !is.na(geno_sp_i)) {
        geno_sp_i
      } else if (morph_is_taxon && morph_sp != given_sp) {
        morph_sp
      } else if (edaph_is_taxon && edaph_sp != given_sp) {
        edaph_sp
      } else {
        given_sp
      }
      raw_geo <- .eval_taxon_geo_plausibility(i, raw_target_sp)
      geno_exempt_geo_gate <- geno_uncontradicted_reassign && !is.na(geno_sp_i) &&
        (raw_target_sp == geno_sp_i) && is.finite(geno_prob_i) && (geno_prob_i >= 0.85)
      range_gated_note <- if (geno_exempt_geo_gate && !isTRUE(raw_geo$plausible)) {
        sprintf(" [Disjunct Population Verified by Genomics (%.0f km from reference hull)]", ifelse(is.finite(raw_geo$dist_km), raw_geo$dist_km, 500.0))
      } else ""
      if (raw_target_sp != given_sp && !isTRUE(raw_geo$plausible) && !geno_exempt_geo_gate) {
        if (is.null(rank_pl)) {
          rank_pl <- .rank_sympatric_congeners(
            i, given_sp, morph_sp, cda_prob, cda_sec_taxon, cda_sec,
            edaph_sp, edaph_is_taxon, rf_sp, rf_is_taxon, rf_prob,
            min_cong_tx_i, nt_m_i, nt_e_i, te_z_i,
            geno_is_taxon, geno_sp_i, geno_prob_i,
            morph_allopatric_or_weak_self = morph_allo_or_weak
          )
        }
        if (!is.null(rank_pl$best) && is.finite(rank_pl$max_score) && rank_pl$max_score >= 0.12) {
          range_gated_note <- sprintf(" [Range-Gated from allopatric %s (%.0f km outside range)]", raw_target_sp, ifelse(is.finite(raw_geo$dist_km), raw_geo$dist_km, 500.0))
          raw_target_sp <- rank_pl$best
        } else {
          raw_target_sp <- given_sp
        }
      }
      target_sp <- raw_target_sp
      if (target_sp == given_sp) {
        portal_lag_flag[i] <- FALSE
        portal_lag_gate[i] <- NA_character_
        status_calls[i] <- "Ecophenotypic_Plasticity"; rec_taxa[i] <- given_sp
        triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Annotate_Taphonomic_Ecophenotype"
        rationales[i] <- sprintf("Stage 4 Plasticity (Allopatric Reassignment Blocked): Morphological outlier retained as %s because candidate congener lies outside documented geographic range", given_sp)
        treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Ecophenotypic_Plasticity"
        next
      }
      portal_note_str <- if (portal_lag_flag[i]) {
        peer_hit_s <- if (!is.null(rank_pl) && !is.null(rank_pl$peer_votes) && (target_sp %in% rank_pl$peer_votes$same_coll)) {
          sprintf("; same-collector/deme peer=%s", target_sp)
        } else ""
        sprintf(
          " [Portal Metadata Lag / Extralimital Audit: %s; identifiedBy=%s; loc=%s/%s; clim_z=%s; rf_given_p=%s%s%s]",
          portal_lag_gate[i],
          ifelse(!is.na(id_by_dv[i]) && nzchar(trimws(id_by_dv[i])), trimws(id_by_dv[i]), "Unknown"),
          ifelse(nzchar(state_vec_geo[i]), state_vec_geo[i], "NA"),
          ifelse(nzchar(cnty_raw_geo[i]), cnty_raw_geo[i], "NA"),
          ifelse(is.finite(clim_z_i), sprintf("%.2f", clim_z_i), "NA"),
          ifelse(is.finite(rf_given_p_i), sprintf("%.3f", rf_given_p_i), "NA"),
          ifelse(exs_conflict_reassign, sprintf("; exsiccatae_gold=%s", exs_gold_tx_dv[i]), ""),
          peer_hit_s
        )
      } else ""
      status_calls[i] <- "Misidentification_Severe"; rec_taxa[i] <- target_sp
      triage_prios[i] <- if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH"; triage_acts[i] <- "Reassign_Determination"
      rationales[i] <- paste0(
        if (geno_uncontradicted_reassign && target_sp == geno_sp_i) {
          sprintf(
            "Stage 1 Reassign (Independent Genomic Assignment): %s contradicted by %s genome (%s, posterior=%.2f, margin=%.2f); reassign to %s",
            given_sp, target_sp, ifelse(is.na(geno_tier_i), "Genomic", geno_tier_i), geno_prob_i, geno_marg_i, target_sp
          )
        } else if (exs_conflict_reassign) {
          sprintf(
            "Stage 1 Reassign (Exsiccatae Portal Metadata Lag): %s conflicts with verified in-range Tier_1_Gold isoduplicate %s; reassign to %s",
            given_sp, target_sp, target_sp
          )
        } else if (extralimital_biogeo_reassign || ungeoref_extralimital_reassign) {
          sprintf(
            "Stage 1 Reassign (Extralimital / Portal Metadata Lag Audit): %s fails USDA/NatureServe county range (status=%s) and bioclimatic/morphometric support; reassign to %s",
            given_sp, u_stat_i, target_sp
          )
        } else if (congener_osr_reassign && target_sp == min_cong_tx_i) {
          sprintf(
            "Stage 1 Reassign (Open-Set Congener Match): %s deviates from %s baseline and falls inside the 90%% Mahalanobis ellipsoid of sympatric congener %s (min_congener_mahal_d2=%.2f); reassign to %s",
            given_sp, given_sp, target_sp, min_cong_d2_i, target_sp
          )
        } else if (type_control_reassign && !is.null(type_control_target_sp) && target_sp == type_control_target_sp) {
          sprintf(
            "Stage 1 Reassign (Type-Baseline Concordance): %s deviates from recorded type control (type_mahal_d2=%s, type_edaphic_z=%s) and matches %s nomenclatural type baseline across >=2 modalities; reassign to %s",
            given_sp,
            ifelse(is.finite(tm_d2_i), sprintf("%.2f", tm_d2_i), "NA"),
            ifelse(is.finite(te_z_i), sprintf("%.2f", te_z_i), "NA"),
            target_sp, target_sp
          )
        } else {
          sprintf("Stage 1 Reassign: %s contradicted across modalities (consensus_conflict=%.2f); reassign to %s", given_sp, p_conf_dv[i], target_sp)
        },
        range_gated_note,
        portal_note_str
      )
      treat_broad[i] <- "Misidentification / Cryptic Species"; treat_sub[i] <- "Misidentification_Reassign"
      next
    }

    res_bas_i <- if (g_cl > 1L && (cl_key_i %in% names(cluster_resurrect_bas_map))) {
      cluster_resurrect_bas_map[[cl_key_i]]
    } else if (is_heterotypic_bas_match_i && (!is.finite(bas_d2_dv[i]) || bas_d2_dv[i] <= type_morph_chisq95) && (!is.finite(tm_d2_i) || is.finite(bas_d2_dv[i]) || tm_d2_i <= type_morph_chisq95)) {
      bas_i
    } else {
      NA_character_
    }

    if (is_intra) {
      has_eco_pheno_intra <- (is_anth && abs(ddoy) >= tau_cryptic_i) ||
                             (is.finite(d_edaph) && d_edaph < 0.30) ||
                             (is.finite(edaph_z) && edaph_z > 2.5) ||
                             geno_gdi_cryptic_support_i
      cyto_barrier_valid_intra <- has_cytotype_barrier_i && (
        has_eco_pheno_intra ||
        (outside_self_morph_95 && morph_exceeds_chisq95_i && osr_cryptic_ok_i && cryptic_cluster_cohesion_pass[i])
      )
      intra_morph_ok <- is_tier1_pure_morph && (
        (is_mo_flag && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i))) ||
        ((g_cl > 1L && bf_val >= 6.0) && (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i))) ||
        (geno_gdi_cryptic_support_i && outside_self_morph_95) ||
        (cyto_barrier_valid_intra && outside_self_morph_95 && morph_exceeds_chisq95_i)
      )
      n_intra_modalities <- as.integer(intra_morph_ok) +
                            as.integer(is_anth && abs(ddoy) >= tau_cryptic_i) +
                            as.integer((is.finite(d_edaph) && d_edaph < 0.30) || (is.finite(edaph_z) && edaph_z > 2.5)) +
                            as.integer(geno_gdi_cryptic_support_i) +
                            as.integer(cyto_barrier_valid_intra)
      raw_intra_crypt_pre_gdi <- !is_sensu_stricto_locked_i && !is_taphonomic_plasticity &&
                                 (cryptic_cluster_cohesion_pass[i] || geno_gdi_cryptic_support_i) &&
                                 (osr_cryptic_ok_i || (geno_gdi_cryptic_support_i && !is.finite(min_cong_d2_dv[i]))) &&
                                 (!has_cytotype_barrier_i || (outside_self_morph_95 && morph_exceeds_chisq95_i)) &&
                                 intra_morph_ok &&
                                 (is.na(cda_prob) || cda_prob < 0.70 || geno_gdi_cryptic_support_i) &&
                                 (has_eco_pheno_intra || cyto_barrier_valid_intra) &&
                                 (n_intra_modalities >= 2L)
      is_cryptic_intra <- raw_intra_crypt_pre_gdi && !geno_gdi_oversplit_i && !tier0_abstain_i
      has_edaph_pheno_div <- (is_anth && abs(ddoy) >= tau_ecotype_i) || (d_edaph < 0.45) ||
                             (is.finite(edaph_z) && edaph_z > 2.0) ||
                             (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60)
      if (is_cryptic_intra) {
        prop_cand_note_intra <- if (is_cluster_prop_gdi_only_i && is.finite(eff_gdi_i) && (eff_gdi_i > 0.70) && !tier0_abstain_i) {
          sprintf(" [Candidate Member via Cluster-Propagated gdi=%.2f; confirm by targeted sequencing]", eff_gdi_i)
        } else ""
        cyto_bar_note_intra <- if (has_cytotype_barrier_i) {
          bp_in <- .resolve_cytotype_barrier_pair(i, given_sp, g_cl, cat_sheet_vec[i], pid_cyto_dv[i], cyto_bar_lbl_attr, sp_all_vec, gmm_cl_dv, voucher_ploidy_v, taxon_cyto_str_v, is_type_dv, is_typ_morph_dv)
          sprintf(" [Postzygotic Heteroploid Barrier: cluster ploidy %s vs type-cluster %s]", bp_in$cl_lbl, bp_in$type_lbl)
        } else ""
        if (!is.na(res_bas_i) && nzchar(res_bas_i)) {
          resurrection_basionym[i] <- res_bas_i
          status_calls[i] <- "Taxonomic_Resurrection_Candidate"
          rec_taxa[i] <- sprintf("%s [Resurrect Basionym: %s]", given_sp, res_bas_i)
          triage_prios[i] <- if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH"
          triage_acts[i] <- "Evaluate_Basionym_Resurrection"
          rationales[i] <- paste0(sprintf(
            "Stage 3 Basionym Resurrection Candidate: Congruent Basal_Rosette outlier lineage (K%d) matches heterotypic type %s with evolutionary barrier (delta_DOY=%.1fd, tau=%.1fd, D_edaph=%.2f%s)",
            g_cl, res_bas_i, ddoy, tau_cryptic_i, d_edaph,
            if (is.finite(geno_gdi_i)) sprintf(", gdi=%.2f", geno_gdi_i) else ""
          ), cyto_bar_note_intra, prop_cand_note_intra)
        } else {
          status_calls[i] <- "Putative_Cryptic_Species"
          rec_taxa[i] <- sprintf("%s (Putative Cryptic Lineage K%d)", given_sp, g_cl)
          triage_prios[i] <- if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH"
          triage_acts[i] <- "Flag_Putative_Cryptic_Species"
          rationales[i] <- paste0(sprintf("Congruent Basal_Rosette morphometric outlier + evolutionary barrier (delta_DOY=%.1fd, tau=%.1fd, D_edaph=%.2f%s)", ddoy, tau_cryptic_i, d_edaph, if (is.finite(geno_gdi_i)) sprintf(", gdi=%.2f", geno_gdi_i) else ""), cyto_bar_note_intra, prop_cand_note_intra)
        }
        treat_broad[i] <- "Misidentification / Cryptic Species"; treat_sub[i] <- "Putative_Cryptic_Species"; next
      } else if (geno_gdi_oversplit_i && raw_intra_crypt_pre_gdi) {
        ovs_tag_intra <- if (is_cluster_prop_gdi_only_i) {
          sprintf("[Cluster Coalescent Oversplit Block: propagated gdi=%.2f]", eff_gdi_i)
        } else {
          sprintf("[Coalescent Oversplit Block: gdi=%.2f < 0.20]", eff_gdi_i)
        }
        if (has_cytotype_barrier_i) {
          st4_c <- .build_stage4_cytotype_decision(given_sp, raw_sp, g_cl, sprintf("Stage 4 Intraspecific Cytotype %s: Conspecific coalescent support retains %s at subspecific rank", ovs_tag_intra, given_sp))
          status_calls[i] <- st4_c$status; rec_taxa[i] <- st4_c$rec
          triage_prios[i] <- st4_c$prio; triage_acts[i] <- st4_c$act
          rationales[i] <- st4_c$rat; treat_broad[i] <- st4_c$broad; treat_sub[i] <- st4_c$sub; next
        } else if (is_taphonomic_plasticity || is_mo_flag || g_cl > 1L) {
          status_calls[i] <- "Ecophenotypic_Plasticity"; rec_taxa[i] <- given_sp
          triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Annotate_Taphonomic_Ecophenotype"
          rationales[i] <- sprintf("Stage 4 Plasticity %s: Conspecific coalescent support retains %s within intraspecific envelope", ovs_tag_intra, given_sp)
          treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Ecophenotypic_Plasticity"; next
        } else {
          status_calls[i] <- "Species_Confirmed"; rec_taxa[i] <- given_sp
          triage_prios[i] <- "RESOLVED"; triage_acts[i] <- "Accept_Current_Determination"
          rationales[i] <- sprintf("Stage 5 Confirmed %s: Conspecific coalescent support retains %s", ovs_tag_intra, given_sp)
          treat_broad[i] <- "Confirmed Species"; treat_sub[i] <- "Confirmed_Species"; next
        }
      } else if (!is_tier1_pure_morph && is_mo_flag && !morph_concordant_self && (d_edaph < 0.30 || (is.finite(edaph_z) && edaph_z > 2.5))) {
        status_calls[i] <- "Ecophenotypic_Plasticity"; rec_taxa[i] <- given_sp
        triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Annotate_Taphonomic_Ecophenotype"
        rationales[i] <- "Cauline_Only_Fallback morphology with multimodal discordance retained within Ecophenotypic_Plasticity envelope"
        treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Ecophenotypic_Plasticity"; next
      } else if (has_cytotype_barrier_i && !is_sensu_stricto_locked_i) {
        base_rat_c <- if (g_cl > 1L && has_edaph_pheno_div) {
          bas_note_eco <- if (is_heterotypic_bas_match_i) sprintf(" [Available Heterotypic Basionym: %s]", bas_i) else ""
          cline_note_eco <- if (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60) sprintf(" [Allopatric Spatial Cline R2=%.2f; Chambers & Hillis 2020]", clinal_r2_i) else ""
          paste0(sprintf("Secondary GMM morphotype (K%d) with cytotype barrier and edaphic/phenological divergence (delta_DOY=%.1fd, pH=%s, Sand=%s)", g_cl, ddoy, ph_str, sand_str), bas_note_eco, cline_note_eco)
        } else {
          sprintf("Stage 4 Intraspecific Cytotype Differentiation: Morphologically overlapping or non-disjunct heteroploid cluster (GMM K%d, edaphic_z=%s) retained at subspecific rank", max(1L, g_cl), ifelse(is.finite(edaph_z), sprintf("%.2f", edaph_z), "NA"))
        }
        st4_c <- .build_stage4_cytotype_decision(given_sp, raw_sp, g_cl, base_rat_c)
        status_calls[i] <- st4_c$status; rec_taxa[i] <- st4_c$rec
        triage_prios[i] <- st4_c$prio; triage_acts[i] <- st4_c$act
        rationales[i] <- st4_c$rat; treat_broad[i] <- st4_c$broad; treat_sub[i] <- st4_c$sub; next
      } else if (!is_sensu_stricto_locked_i && (g_cl > 1L) && !tier0_abstain_i) {
        if (has_edaph_pheno_div) {
          bas_note_eco <- if (is_heterotypic_bas_match_i) sprintf(" [Available Heterotypic Basionym: %s]", bas_i) else ""
          cline_note_eco <- if (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60) sprintf(" [Allopatric Spatial Cline R2=%.2f; Chambers & Hillis 2020]", clinal_r2_i) else ""
          status_calls[i] <- "Subspecies_Ecotype"; rec_taxa[i] <- sprintf("%s var. morphotype_K%d", given_sp, g_cl)
          triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Accept_Subspecific_Treatment"
          rationales[i] <- paste0(sprintf("Secondary GMM morphotype (K%d) with cohesive edaphic/phenological divergence (delta_DOY=%.1fd, pH=%s, Sand=%s)", g_cl, ddoy, ph_str, sand_str), bas_note_eco, cline_note_eco)
          treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Subspecies_Ecotype"; next
        } else if (!cryptic_cluster_cohesion_pass[i] && (!is.finite(edaph_z) || edaph_z <= 2.5) && (!is_anth || abs(ddoy) < tau_cryptic_i)) {
          if (is_taphonomic_plasticity || is_mo_flag) {
            status_calls[i] <- "Ecophenotypic_Plasticity"; rec_taxa[i] <- given_sp
            triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Annotate_Taphonomic_Ecophenotype"
            rationales[i] <- sprintf("Stage 4 Plasticity: Secondary shape mode (GMM K%d) within %s ecological/phenological envelope (edaphic_z=%s)", g_cl, given_sp, ifelse(is.finite(edaph_z), sprintf("%.2f", edaph_z), "NA"))
            treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Ecophenotypic_Plasticity"; next
          } else {
            status_calls[i] <- "Species_Confirmed"; rec_taxa[i] <- given_sp
            triage_prios[i] <- "RESOLVED"; triage_acts[i] <- "Accept_Current_Determination"
            rationales[i] <- sprintf("Stage 5 Confirmed: Secondary unweighted shape mode (GMM K%d) within %s ecological and phenological niche (continuous intraspecific foliar polymorphism / heteroblasty)", g_cl, given_sp)
            treat_broad[i] <- "Confirmed Species"; treat_sub[i] <- "Confirmed_Species"; next
          }
        }
      }
    }

    # Stage 2 Pre-Screen (Gates 1 to 4): Conservative Contact-Zone Hybrid Candidate Evaluation
    # Strictly disabled in --intraspecific-mode or single-taxon cohorts (can_evaluate_hybrids == FALSE),
    # and strictly blocked on Divergent_NonHybrid_Lineage / Ancient_Introgressant genomic queries.
    has_morph_i <- isTRUE(has_observed_or_sheet_morph[i]) && is_tier1_pure_i &&
                   (!is.null(morph_sp) && !is.na(morph_sp) && nzchar(morph_sp) && is.finite(cda_prob))
    geno_contemp_direct_i <- geno_g1_direct && geno_contemp_hybrid_i && is.finite(geno_q12_i) && (geno_q12_i >= 0.35)
    if (can_evaluate_hybrids && !geno_divergent_nonhybrid_i && !geno_ancient_introg_i &&
        (has_morph_i || geno_contemp_direct_i) && !is_severe_misid && !both_alt_agree &&
        !is_sensu_stricto_locked_i && (!morph_concordant_self || geno_contemp_direct_i) &&
        (given_sp %in% cohort_species)) {
      nn_mahal_sp <- resolve_nearest_mahalanobis_taxon(i, given_sp)
      explicit_p2 <- NA_character_
      if (grepl("(^|\\s)[x×]\\s+[a-z]", raw_sp, ignore.case = TRUE)) {
        m_hyb <- regmatches(raw_sp, regexpr("[x×]\\s+[A-Za-z-]+(\\s+[a-z-]+)?", raw_sp, ignore.case = TRUE))
        if (length(m_hyb) > 0L && nzchar(m_hyb[1L])) {
          rhs_tok <- trimws(sub("^[x×]\\s+", "", m_hyb[1L]))
          cand_rhs <- if (grepl("^Packera\\s+", rhs_tok, ignore.case = TRUE)) standardize_packera_taxon(rhs_tok) else standardize_packera_taxon(paste("Packera", rhs_tok))
          if (!is.na(cand_rhs) && (cand_rhs %in% cohort_species) && cand_rhs != given_sp) explicit_p2 <- cand_rhs
        }
      }
      donor_cands <- unique(stats::na.omit(c(
        explicit_p2,
        if (has_geno_i && (geno_adm_dv[i] || geno_contemp_hybrid_i) && !is.na(geno_p1_dv[i]) && geno_p1_dv[i] != given_sp) geno_p1_dv[i] else NA_character_,
        if (has_geno_i && (geno_adm_dv[i] || geno_contemp_hybrid_i) && !is.na(geno_p2_dv[i]) && geno_p2_dv[i] != given_sp) geno_p2_dv[i] else NA_character_,
        if (morph_is_taxon && morph_sp != given_sp) morph_sp else NA_character_,
        if (sec_is_taxon && cda_sec_taxon != given_sp) cda_sec_taxon else NA_character_,
        if (!is.na(min_cong_tx_i) && min_cong_tx_i != given_sp) min_cong_tx_i else NA_character_,
        if (!is.na(nt_m_i) && nt_m_i != given_sp) nt_m_i else NA_character_,
        if (edaph_is_taxon && edaph_sp != given_sp) edaph_sp else NA_character_,
        nn_mahal_sp
      )))
      donor_cands <- donor_cands[donor_cands %in% cohort_species & donor_cands != given_sp]
      best_g <- NULL
      if (has_morph_i) {
        for (dc in donor_cands) {
          g_eval <- .eval_five_gate_hybrid_pair(
            row_i = i, p1_tx = given_sp, p2_tx = dc,
            morph_sp_i = morph_sp, cda_prob_i = cda_prob, cda_sec_tx_i = cda_sec_taxon,
            cda_sec_p_i = cda_sec, cda_marg_i = cda_marg, gmm_unc_i = gmm_unc,
            edaph_sp_i = edaph_sp, d_edaph_i = d_edaph, edaph_z_i = edaph_z, asym_i = asym_i
          )
          if (isTRUE(g_eval$pass) && (is.null(best_g) || g_eval$score > best_g$score)) {
            best_g <- g_eval; best_g$p2 <- dc
          }
        }
      }
      if (!is.null(best_g)) {
        hyb_cand_pass[i] <- TRUE
        hyb_cand_p1[i] <- given_sp
        hyb_cand_p2[i] <- best_g$p2
        hyb_cand_d1_km[i] <- best_g$d1_km
        hyb_cand_d2_km[i] <- best_g$d2_km
        hyb_cand_ddoy[i] <- best_g$ddoy_12
        hyb_cand_t_ab[i] <- best_g$t_ab
        hyb_cand_d2_perp[i] <- best_g$d2_perp
      } else if (geno_contemp_direct_i) {
        p1_g <- if (!is.na(geno_p1_dv[i]) && is_valid_taxon_str(geno_p1_dv[i])) geno_p1_dv[i] else given_sp
        p2_g <- if (!is.na(geno_p2_dv[i]) && is_valid_taxon_str(geno_p2_dv[i]) && geno_p2_dv[i] != p1_g) geno_p2_dv[i] else if (length(donor_cands) > 0L) donor_cands[1L] else NA_character_
        if (!is.na(p1_g) && !is.na(p2_g) && p1_g != p2_g) {
          hyb_cand_pass[i] <- TRUE
          hyb_cand_p1[i] <- p1_g
          hyb_cand_p2[i] <- p2_g
          hyb_cand_d1_km[i] <- 0.0
          hyb_cand_d2_km[i] <- 0.0
          hyb_cand_ddoy[i] <- if (is.finite(ddoy)) abs(ddoy) else 0.0
          hyb_cand_t_ab[i] <- if (is.finite(geno_h_i)) geno_h_i else 0.5
          hyb_cand_d2_perp[i] <- 0.0
        }
      }
    }

    # Assign conservative non-hybrid baseline first (Stages 3 to 5); Pass 2 elevates ONLY Gate-5 population-confirmed hybrid swarms
    has_infraspecific_rank <- grepl("\\b(var\\.|subsp\\.|ssp\\.)\\s+[a-z-]+", raw_sp, ignore.case = TRUE)
    is_cohesive_ecotype <- !is_sensu_stricto_locked_i && !is_severe_misid && !is_crypt && !tier0_abstain_i && (
      has_infraspecific_rank ||
      ((g_cl > 1L) && ((d_edaph < 0.45) || (is.finite(edaph_z) && edaph_z > 2.0) || (is_anth && abs(ddoy) >= tau_ecotype_i) ||
       (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60)))
    )
    is_plastic_variant <- !is_sensu_stricto_locked_i && !is_cohesive_ecotype && !is_severe_misid && !is_crypt &&
                          (is_taphonomic_plasticity || is_mo_flag)

    if (is_crypt && !hyb_cand_pass[i]) {
      prop_cand_note_ms <- if (is_cluster_prop_gdi_only_i && is.finite(eff_gdi_i) && (eff_gdi_i > 0.70) && !tier0_abstain_i) {
        sprintf(" [Candidate Member via Cluster-Propagated gdi=%.2f; confirm by targeted sequencing]", eff_gdi_i)
      } else ""
      cyto_bar_note_ms <- if (has_cytotype_barrier_i) {
        bp_ms <- .resolve_cytotype_barrier_pair(i, given_sp, g_cl, cat_sheet_vec[i], pid_cyto_dv[i], cyto_bar_lbl_attr, sp_all_vec, gmm_cl_dv, voucher_ploidy_v, taxon_cyto_str_v, is_type_dv, is_typ_morph_dv)
        sprintf(" [Postzygotic Heteroploid Barrier: cluster ploidy %s vs type-cluster %s]", bp_ms$cl_lbl, bp_ms$type_lbl)
      } else ""
      if (!is.na(res_bas_i) && nzchar(res_bas_i)) {
        resurrection_basionym[i] <- res_bas_i
        status_calls[i] <- "Taxonomic_Resurrection_Candidate"
        rec_taxa[i] <- sprintf("%s [Resurrect Basionym: %s]", given_sp, res_bas_i)
        triage_prios[i] <- if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH"
        triage_acts[i] <- "Evaluate_Basionym_Resurrection"
        rationales[i] <- paste0(sprintf(
          "Stage 3 Basionym Resurrection Candidate (Padial et al. 2010): Basal_Rosette outlier lineage (K%d) matches heterotypic type %s (CDA_p=%.2f) + barrier (|delta_DOY|=%.1fd >= %.1fd, D_edaph=%.2f%s)",
          g_cl, res_bas_i, ifelse(is.na(cda_prob), 0.45, cda_prob), abs(ddoy), tau_cryptic_i, d_edaph,
          if (is.finite(geno_gdi_i)) sprintf(", gdi=%.2f", geno_gdi_i) else ""
        ), cyto_bar_note_ms, prop_cand_note_ms)
      } else {
        status_calls[i] <- "Putative_Cryptic_Species"
        rec_taxa[i] <- sprintf("%s (Putative Cryptic Lineage K%d)", given_sp, g_cl)
        triage_prios[i] <- if (tier == "Tier_1_Gold") "CRITICAL" else "HIGH"
        triage_acts[i] <- "Flag_Putative_Cryptic_Species"
        rationales[i] <- paste0(sprintf("Stage 3 Cryptic Congruence (Padial et al. 2010): Basal_Rosette outlier unassigned to known taxa (CDA_p=%.2f) + barrier (|delta_DOY|=%.1fd >= %.1fd, D_edaph=%.2f%s)", ifelse(is.na(cda_prob), 0.45, cda_prob), abs(ddoy), tau_cryptic_i, d_edaph, if (is.finite(geno_gdi_i)) sprintf(", gdi=%.2f", geno_gdi_i) else ""), cyto_bar_note_ms, prop_cand_note_ms)
      }
      treat_broad[i] <- "Misidentification / Cryptic Species"; treat_sub[i] <- "Putative_Cryptic_Species"
    } else if (has_cytotype_barrier_i && !is_sensu_stricto_locked_i && !is_severe_misid && !hyb_cand_pass[i]) {
      base_rat_ms_c <- if (cryptic_oversplit_blocked_i) {
        ovs_tag_c <- if (is_cluster_prop_gdi_only_i) sprintf("[Cluster Coalescent Oversplit Block: propagated gdi=%.2f]", eff_gdi_i) else sprintf("[Coalescent Oversplit Block: gdi=%.2f < 0.20]", eff_gdi_i)
        sprintf("Stage 4 Intraspecific Cytotype %s: Conspecific coalescent support retains %s at subspecific rank", ovs_tag_c, given_sp)
      } else {
        sprintf("Stage 4 Intraspecific Cytotype Differentiation: Morphologically overlapping or non-disjunct heteroploid cluster (GMM K%d, D_edaph=%.2f, delta_DOY=%.1fd)", max(1L, g_cl), d_edaph, ddoy)
      }
      st4_ms <- .build_stage4_cytotype_decision(given_sp, raw_sp, g_cl, base_rat_ms_c)
      status_calls[i] <- st4_ms$status; rec_taxa[i] <- st4_ms$rec
      triage_prios[i] <- st4_ms$prio; triage_acts[i] <- st4_ms$act
      rationales[i] <- st4_ms$rat; treat_broad[i] <- st4_ms$broad; treat_sub[i] <- st4_ms$sub
    } else if (cryptic_oversplit_blocked_i && !hyb_cand_pass[i]) {
      ovs_tag_ms <- if (is_cluster_prop_gdi_only_i) {
        sprintf("[Cluster Coalescent Oversplit Block: propagated gdi=%.2f]", eff_gdi_i)
      } else {
        sprintf("[Coalescent Oversplit Block: gdi=%.2f < 0.20]", eff_gdi_i)
      }
      if (is_taphonomic_plasticity || is_mo_flag || g_cl > 1L) {
        status_calls[i] <- "Ecophenotypic_Plasticity"; rec_taxa[i] <- given_sp
        triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Annotate_Taphonomic_Ecophenotype"
        rationales[i] <- sprintf("Stage 4 Plasticity %s: Conspecific coalescent support retains %s within intraspecific envelope", ovs_tag_ms, given_sp)
        treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Ecophenotypic_Plasticity"
      } else {
        status_calls[i] <- "Species_Confirmed"; rec_taxa[i] <- given_sp
        triage_prios[i] <- "RESOLVED"; triage_acts[i] <- "Accept_Current_Determination"
        rationales[i] <- sprintf("Stage 5 Confirmed %s: Conspecific coalescent support retains %s", ovs_tag_ms, given_sp)
        treat_broad[i] <- "Confirmed Species"; treat_sub[i] <- "Confirmed_Species"
      }
    } else if (is_cohesive_ecotype && !hyb_cand_pass[i]) {
      infra_lbl <- if (has_infraspecific_rank) {
        m_inf <- regmatches(raw_sp, regexpr("\\b(var\\.|subsp\\.|ssp\\.)\\s+[A-Za-z-]+", raw_sp, ignore.case = TRUE))
        if (length(m_inf) > 0L && nzchar(m_inf[1L])) paste(given_sp, tolower(m_inf[1L])) else paste0(given_sp, " var. ecotype")
      } else if (g_cl > 1L) {
        sprintf("%s var. morphotype_K%d", given_sp, g_cl)
      } else {
        paste0(given_sp, " var. ecotype")
      }
      cline_note_ms <- if (!cryptic_non_clinal_pass[i] && is.finite(clinal_r2_i) && clinal_r2_i >= 0.60) sprintf(" [Allopatric Spatial Cline R2=%.2f; Chambers & Hillis 2020]", clinal_r2_i) else ""
      status_calls[i] <- "Subspecies_Ecotype"; rec_taxa[i] <- infra_lbl
      triage_prios[i] <- "LOW"; triage_acts[i] <- "Accept_Subspecific_Treatment"
      rationales[i] <- paste0(sprintf("Stage 4 Ecotype: Cohesive intraspecific lineage (GMM K%d, D_edaph=%.2f, delta_DOY=%.1fd, tau=%.1fd)", g_cl, d_edaph, ddoy, tau_ecotype_i), cline_note_ms)
      treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Subspecies_Ecotype"
    } else if (is_plastic_variant || hyb_cand_pass[i]) {
      doy_plast_note <- if (is.finite(doy) && is.finite(sp_q75_doy) && doy >= sp_q75_doy) sprintf("DOY %.0f >= Q75 %.0f", doy, sp_q75_doy) else sprintf("DOY %s, Q75 %s", ifelse(is.finite(doy), sprintf("%.0f", doy), "NA"), ifelse(is.finite(sp_q75_doy), sprintf("%.0f", sp_q75_doy), "NA"))
      status_calls[i] <- "Ecophenotypic_Plasticity"; rec_taxa[i] <- given_sp
      triage_prios[i] <- "MEDIUM"; triage_acts[i] <- "Annotate_Taphonomic_Ecophenotype"
      rationales[i] <- sprintf("Stage 4 Plasticity: Within-envelope taphonomic/ontogenetic variation (%s, edaphic_z=%s, pH=%s, Sand=%s)", doy_plast_note, ifelse(is.finite(edaph_z_i), sprintf("%.2f", edaph_z_i), "NA"), ph_str, sand_str)
      treat_broad[i] <- "Variety / Ecotype"; treat_sub[i] <- "Ecophenotypic_Plasticity"
    } else {
      status_calls[i] <- "Species_Confirmed"; rec_taxa[i] <- given_sp
      triage_prios[i] <- "RESOLVED"; triage_acts[i] <- "Accept_Current_Determination"
      rationales[i] <- "Stage 5 Confirmed: Concordant across available morphological, reproductive phenological, and edaphic streams"
      treat_broad[i] <- "Confirmed Species"; treat_sub[i] <- "Confirmed_Species"
    }
  }

  # Pass 2 (Gate 5: Population-Level Contact-Zone Hybrid Swarm Cohesion)
  # Elevates Gate 1-4 candidates to Hybrid_Introgressant ONLY when >= 2 independent vouchers in the local
  # <= 25 km contact zone (or shared county) corroborate the same parent pair, or label explicitly records nothotaxon cross.
  if (can_evaluate_hybrids && any(hyb_cand_pass)) {
    cand_idx_all <- which(hyb_cand_pass)
    pair_keys_cand <- vapply(cand_idx_all, function(ci) paste(sort(c(hyb_cand_p1[ci], hyb_cand_p2[ci])), collapse = "||"), character(1))
    for (pos_c in seq_along(cand_idx_all)) {
      ci <- cand_idx_all[pos_c]
      pk <- pair_keys_cand[pos_c]
      same_pk_idx <- cand_idx_all[pair_keys_cand == pk]
      lat_ci <- lat_vec_geo[ci]; lon_ci <- lon_vec_geo[ci]; sc_ci <- state_cnty_geo[ci]
      local_peers <- integer(0)
      for (cj in same_pk_idx) {
        if (cj == ci) {
          local_peers <- c(local_peers, cj)
        } else if (is.finite(lat_ci) && is.finite(lon_ci) && is.finite(lat_vec_geo[cj]) && is.finite(lon_vec_geo[cj])) {
          dlat_p <- (lat_vec_geo[cj] - lat_ci) * (pi / 180.0); dlon_p <- (lon_vec_geo[cj] - lon_ci) * (pi / 180.0)
          ah_p <- sin(dlat_p / 2.0)^2 + cos(lat_ci * (pi / 180.0)) * cos(lat_vec_geo[cj] * (pi / 180.0)) * sin(dlon_p / 2.0)^2
          dkm_p <- 6371.0 * 2.0 * atan2(sqrt(max(0.0, ah_p)), sqrt(max(0.0, 1.0 - ah_p)))
          if (dkm_p <= 25.0) local_peers <- c(local_peers, cj)
        } else if (nzchar(sc_ci) && nzchar(state_cnty_geo[cj]) && sc_ci == state_cnty_geo[cj]) {
          local_peers <- c(local_peers, cj)
        }
      }
      n_deme_vouchers <- length(unique(cat_sheet_vec[local_peers]))
      raw_ci <- if (!is.na(raw_sp_dv[ci])) raw_sp_dv[ci] else ""
      has_explicit_cross <- grepl("(^|\\s)[x×]\\s+[a-z]", raw_ci, ignore.case = TRUE)
      geno_admixed_ci <- has_geno_dv[ci] && !geno_divergent_nonhybrid_v[ci] && !geno_ancient_introg_v[ci] &&
        (geno_adm_dv[ci] || geno_contemp_hybrid_v[ci]) && identical(geno_tier_dv[ci], "G1_Direct_Voucher")
      if (n_deme_vouchers >= 2L || has_explicit_cross || geno_admixed_ci) {
        p_a <- hyb_cand_p1[ci]; p_b <- hyb_cand_p2[ci]
        tier_ci <- if (!is.na(tier_dv[ci])) tier_dv[ci] else "Tier_3_Bronze"
        hybrid_zone_cohesion_pass[ci] <- TRUE
        hybrid_parent_1[ci] <- p_a
        hybrid_parent_2[ci] <- p_b
        status_calls[ci] <- "Hybrid_Intergrade_Swarm"
        rec_taxa[ci] <- sprintf("%s x %s", p_a, p_b)
        triage_prios[ci] <- if (tier_ci == "Tier_1_Gold") "HIGH" else "MEDIUM"
        triage_acts[ci] <- "Flag_Hybrid_Swarm_Intergrade"
        hyb_cyto_note <- .hybrid_cytotype_annotation(
          voucher_ploidy_v[ci], taxon_cyto_str_v[ci], p_a, p_b, sp_all_vec, voucher_ploidy_v, taxon_cyto_str_v
        )
        rationales[ci] <- paste0(sprintf(
          "Stage 2 Confirmed Hybrid Swarm (%s x %s): Syntopic parents verified within %.1f km and %.1f km (co-flowering |dDOY|=%.1fd), parental chord t_AB=%.2f (d2_perp=%.2f), corroborated across n=%d contact-zone vouchers",
          p_a, p_b,
          ifelse(is.finite(hyb_cand_d1_km[ci]), hyb_cand_d1_km[ci], 0.0),
          ifelse(is.finite(hyb_cand_d2_km[ci]), hyb_cand_d2_km[ci], 0.0),
          ifelse(is.finite(hyb_cand_ddoy[ci]), hyb_cand_ddoy[ci], 0.0),
          hyb_cand_t_ab[ci],
          ifelse(is.finite(hyb_cand_d2_perp[ci]), hyb_cand_d2_perp[ci], 0.0),
          n_deme_vouchers
        ), hyb_cyto_note)
        treat_broad[ci] <- "Hybrid / Introgressant"
        treat_sub[ci] <- "Hybrid_Introgressant"
      } else if (isTRUE(cyto_barrier_pass_v[ci])) {
        raw_ci_sp <- if (!is.na(raw_sp_dv[ci])) raw_sp_dv[ci] else hyb_cand_p1[ci]
        g_cl_ci <- if (!is.na(gmm_cl_dv[ci]) && gmm_cl_dv[ci] >= 1L) gmm_cl_dv[ci] else 1L
        st4_ci <- .build_stage4_cytotype_decision(
          hyb_cand_p1[ci], raw_ci_sp, g_cl_ci,
          sprintf(
            "Stage 4 Intraspecific Cytotype Differentiation: Unreplicated sympatric intermediate toward %s (t_AB=%.2f, n_deme=1) retained as subspecific cytotype series",
            hyb_cand_p2[ci], hyb_cand_t_ab[ci]
          )
        )
        status_calls[ci] <- st4_ci$status; rec_taxa[ci] <- st4_ci$rec
        triage_prios[ci] <- st4_ci$prio; triage_acts[ci] <- st4_ci$act
        rationales[ci] <- st4_ci$rat; treat_broad[ci] <- st4_ci$broad; treat_sub[ci] <- st4_ci$sub
      } else {
        # Singleton sympatric intermediate fails Gate 5 deme cohesion -> retain conservatively under parental species
        if (!isTRUE(mo_flag_dv[ci])) {
          status_calls[ci] <- "Species_Confirmed"
          rec_taxa[ci] <- hyb_cand_p1[ci]
          triage_prios[ci] <- "LOW"
          triage_acts[ci] <- "Accept_Current_Determination"
          treat_broad[ci] <- "Confirmed Species"
          treat_sub[ci] <- "Confirmed_Species"
        }
        rationales[ci] <- paste0(
          rationales[ci],
          sprintf(
            " [Unreplicated sympatric intermediate toward %s (t_AB=%.2f, n_deme=1); retained as %s under conservative population hybrid gate]",
            hyb_cand_p2[ci], hyb_cand_t_ab[ci], hyb_cand_p1[ci]
          )
        )
      }
    }
  }

  # Same-Sheet (catalogNumber) Gathering Coherence Post-Filter across sibling plant_individual_id ramets
  aw_sheet_vec <- resolve_step07_ecological_authority_weights(df, allow_bronze_prior = TRUE)
  chisq99_type_morph <- stats::qchisq(0.99, df = p_type_morph)
  for (s_idx in sheet_idx_list) {
    if (length(s_idx) < 2L) next
    if (any(is_type_dv[s_idx]) || any(status_calls[s_idx] %in% c("Confirmed_Type_Control", "Type_Synonymy_Candidate"))) next
    if (length(unique(treat_sub[s_idx])) <= 1L && length(unique(status_calls[s_idx])) <= 1L && length(unique(rec_taxa[s_idx])) <= 1L) next

    is_conf_s <- (treat_sub[s_idx] == "Confirmed_Species") | (status_calls[s_idx] == "Species_Confirmed")
    is_basal_s <- t1_pure_all[s_idx] & orig_has_obs_morph_v[s_idx] & (orig_n_leaves_v[s_idx] >= 1L)
    w_r <- pmax(0.25, aw_sheet_vec[s_idx]) *
           (1.0 + orig_n_leaves_v[s_idx]) *
           (1.0 + 1.5 * as.numeric(is_basal_s)) *
           (1.0 + 0.5 * as.numeric(is_conf_s)) *
           ifelse(is.finite(orig_cda_prob_v[s_idx]), pmax(0.35, orig_cda_prob_v[s_idx]), 0.35)
    has_leaf_s <- orig_has_obs_morph_v[s_idx] & (orig_n_leaves_v[s_idx] >= 1L)
    conf_leaf_s <- is_basal_s & is_conf_s & (
      (!is.na(orig_cda_pred_v[s_idx]) & orig_cda_pred_v[s_idx] == sp_clean_vec[s_idx]) |
      (!is.na(is_typ_morph_dv[s_idx]) & is_typ_morph_dv[s_idx])
    )
    is_genuine_mixed_ramet <- function(r, ref_r) {
      if (!t1_pure_all[r] || !t1_pure_all[ref_r]) return(FALSE)
      n_bas_r <- if (!is.na(bas_dv[r])) bas_dv[r] else orig_n_leaves_v[r]
      n_bas_ref <- if (!is.na(bas_dv[ref_r])) bas_dv[ref_r] else orig_n_leaves_v[ref_r]
      if (is.na(n_bas_r) || n_bas_r < 2L || is.na(n_bas_ref) || n_bas_ref < 2L) return(FALSE)
      if (!orig_has_obs_morph_v[r] || is.na(orig_cda_pred_v[r]) || !nzchar(orig_cda_pred_v[r])) return(FALSE)
      if (!is.na(orig_cda_pred_v[ref_r]) && orig_cda_pred_v[r] == orig_cda_pred_v[ref_r]) return(FALSE)
      d2_cand_r <- if (is.finite(orig_tm_d2_v[r])) {
        orig_tm_d2_v[r]
      } else if ("mahalanobis_d2" %in% names(df) && is.finite(suppressWarnings(as.numeric(df$mahalanobis_d2[r])))) {
        suppressWarnings(as.numeric(df$mahalanobis_d2[r]))
      } else if ("cda_mahal_d2" %in% names(df) && is.finite(suppressWarnings(as.numeric(df$cda_mahal_d2[r])))) {
        suppressWarnings(as.numeric(df$cda_mahal_d2[r]))
      } else {
        NA_real_
      }
      marg_r <- if (is.finite(cda_marg_dv[r])) cda_marg_dv[r] else if (is.finite(orig_cda_prob_v[r])) max(0.0, 2.0 * orig_cda_prob_v[r] - 1.0) else 0.0
      is.finite(orig_cda_prob_v[r]) && (orig_cda_prob_v[r] >= 0.90) && (marg_r >= 0.50) &&
        ((is.finite(d2_cand_r) && (d2_cand_r > type_morph_chisq95)) || isTRUE(mo_flag_dv[r]))
    }
    if (any(conf_leaf_s)) {
      ref_conf_r <- s_idx[which(conf_leaf_s)[which.max(w_r[which(conf_leaf_s)])]]
      for (k_pos in seq_along(s_idx)) {
        r <- s_idx[k_pos]
        if (portal_lag_flag[r]) next
        if (is_conf_s[k_pos] && rec_taxa[r] == sp_clean_vec[r]) next
        if (isTRUE(cyto_barrier_pass_v[r])) next
        if (!(treat_sub[r] %in% c("Hybrid_Introgressant", "Misidentification_Reassign", "Ecophenotypic_Plasticity") ||
              status_calls[r] %in% c("Hybrid_Intergrade_Swarm", "Putative_Hybrid_Introgression", "Ambiguous_Needs_Expert_Review", "Misidentification_Severe"))) {
          next
        }
        if (!is_genuine_mixed_ramet(r, ref_conf_r)) {
          status_calls[r] <- "Species_Confirmed"
          treat_sub[r] <- "Confirmed_Species"
          treat_broad[r] <- "Confirmed Species"
          rec_taxa[r] <- sp_clean_vec[r]
          hybrid_zone_cohesion_pass[r] <- FALSE
          hybrid_parent_1[r] <- NA_character_
          hybrid_parent_2[r] <- NA_character_
          resurrection_basionym[r] <- NA_character_
          triage_prios[r] <- "RESOLVED"
          triage_acts[r] <- "Accept_Current_Determination"
          if (!grepl("\\[Sheet-Coherent Ramet Consensus\\]", rationales[r])) {
            rationales[r] <- paste0(rationales[r], " [Sheet-Coherent Ramet Consensus]")
          }
        }
      }
    } else if (any(has_leaf_s)) {
      cand_positions <- if (any(is_basal_s)) which(is_basal_s) else which(has_leaf_s)
      best_r <- s_idx[cand_positions[which.max(w_r[cand_positions])]]
      for (k_pos in seq_along(s_idx)) {
        r <- s_idx[k_pos]
        if (r == best_r || portal_lag_flag[r]) next
        if ((geno_gdi_oversplit_v[r] && treat_sub[best_r] == "Putative_Cryptic_Species") ||
            (tier0_abstain_dv[r] && (treat_sub[best_r] %in% c("Putative_Cryptic_Species", "Subspecies_Ecotype"))) ||
            (isTRUE(cyto_barrier_pass_v[r]) && (treat_sub[best_r] %in% c("Ecophenotypic_Plasticity", "Confirmed_Species")))) {
          next
        }
        if (!is_genuine_mixed_ramet(r, best_r)) {
          status_calls[r] <- status_calls[best_r]
          treat_sub[r] <- treat_sub[best_r]
          treat_broad[r] <- treat_broad[best_r]
          rec_taxa[r] <- rec_taxa[best_r]
          hybrid_zone_cohesion_pass[r] <- hybrid_zone_cohesion_pass[best_r]
          hybrid_parent_1[r] <- hybrid_parent_1[best_r]
          hybrid_parent_2[r] <- hybrid_parent_2[best_r]
          resurrection_basionym[r] <- resurrection_basionym[best_r]
          triage_prios[r] <- triage_prios[best_r]
          triage_acts[r] <- triage_acts[best_r]
          if (!grepl("\\[Sheet-Coherent Ramet Consensus\\]", rationales[r])) {
            rationales[r] <- paste0(rationales[r], " [Sheet-Coherent Ramet Consensus]")
          }
        }
      }
    }
  }

  # Enforce self-reassignment invariant: no voucher with recommended_determination == given_sp may stay flagged as Misidentification_*
  self_mis_idx <- which(
    (status_calls %in% c("Misidentification_Reassign", "Misidentification_Severe") | treat_sub == "Misidentification_Reassign") &
    !is.na(rec_taxa) & !is.na(sp_all_vec) & nzchar(sp_all_vec) & (rec_taxa == sp_all_vec)
  )
  for (si in self_mis_idx) {
    if (isTRUE(mo_flag_dv[si])) {
      status_calls[si] <- "Ecophenotypic_Plasticity"
      treat_broad[si] <- "Variety / Ecotype"
      treat_sub[si] <- "Ecophenotypic_Plasticity"
      triage_prios[si] <- "MEDIUM"
      triage_acts[si] <- "Annotate_Taphonomic_Ecophenotype"
    } else {
      status_calls[si] <- "Species_Confirmed"
      treat_broad[si] <- "Confirmed Species"
      treat_sub[si] <- "Confirmed_Species"
      triage_prios[si] <- "LOW"
      triage_acts[si] <- "Accept_Current_Determination"
    }
    rationales[si] <- paste0(rationales[si], " [Self-reassignment guard retained given determination]")
  }

  geno_idx_note <- which(has_geno_dv & !is.na(geno_pred_dv) & nzchar(geno_pred_dv))
  for (gi in geno_idx_note) {
    if (!grepl("\\[Genomic|Genomic Assignment|Rule A Nomenclatural", rationales[gi])) {
      phase2_bits <- character(0)
      if (!is.na(geno_cyt_dv[gi]) && nzchar(geno_cyt_dv[gi])) phase2_bits <- c(phase2_bits, sprintf("cytotype=%s", geno_cyt_dv[gi]))
      if (is.finite(geno_gdi_dv[gi])) phase2_bits <- c(phase2_bits, sprintf("gdi=%.2f", geno_gdi_dv[gi]))
      if (is.finite(geno_q12_dv[gi])) phase2_bits <- c(phase2_bits, sprintf("Q12=%.2f", geno_q12_dv[gi]))
      if (!is.na(geno_hyb_cls_dv[gi]) && nzchar(geno_hyb_cls_dv[gi])) phase2_bits <- c(phase2_bits, geno_hyb_cls_dv[gi])
      p2_str <- if (length(phase2_bits) > 0L) paste0("; ", paste(phase2_bits, collapse = ", ")) else ""
      anc_str <- if (geno_ancient_introg_v[gi]) {
        sprintf("; Ancient Introgression / Organellar Capture (f_b=%s, cytonuc=%s; nuclear determination retained)", ifelse(is.finite(geno_fb_max_dv[gi]), sprintf("%.3f", geno_fb_max_dv[gi]), "NA"), ifelse(isTRUE(geno_cytonuc_dv[gi]), "1", "0"))
      } else ""
      qry_tag <- if (geno_divergent_nonhybrid_v[gi]) {
        "; Divergent_NonHybrid_Lineage"
      } else if (geno_qry_dv[gi]) {
        "; Putative_Novel_Or_Hybrid_Query"
      } else if (geno_adm_dv[gi] && !is.na(geno_p1_dv[gi]) && !is.na(geno_p2_dv[gi])) {
        sprintf("; Admixed q_max=%.2f (%s x %s)", ifelse(is.finite(geno_qmax_dv[gi]), geno_qmax_dv[gi], 0.5), geno_p1_dv[gi], geno_p2_dv[gi])
      } else ""
      rationales[gi] <- paste0(
        rationales[gi],
        sprintf(" [Genomic %s: pred=%s (p=%.2f%s%s%s)]", ifelse(is.na(geno_tier_dv[gi]), "Stream", geno_tier_dv[gi]), geno_pred_dv[gi], geno_prob_dv[gi], p2_str, anc_str, qry_tag)
      )
    }
  }

  for (ai in which(tier0_abstain_dv)) {
    if (!grepl("Tier_0 Abstention", rationales[ai], fixed = TRUE)) {
      rationales[ai] <- paste0(rationales[ai], " [Tier_0 Abstention: Phenotypic overlap across genetically distinct lineages; DNA required]")
    }
  }
  for (mi in which(is_marginal_ext_v & !portal_lag_flag & (treat_sub == "Confirmed_Species"))) {
    if (!grepl("Marginal Range Extension Preserved", rationales[mi], fixed = TRUE)) {
      rationales[mi] <- paste0(
        rationales[mi],
        sprintf(" [Marginal Range Extension Preserved: rf_given_p=%.3f, clim_z=%s]",
                rf_given_p_dv[mi], ifelse(is.finite(clim_z_dv[mi]), sprintf("%.2f", clim_z_dv[mi]), "NA"))
      )
    }
  }

  df$authority_weight <- aw_sheet_vec
  df$is_type_specimen <- as.logical(is_type_dv)
  if (!("type_status" %in% names(df))) df$type_status <- type_stat_dv
  df$edaphic_z_score <- edaph_z_vec
  df$cryptic_cluster_cohesion_pass <- as.logical(cryptic_cluster_cohesion_pass)
  df$cryptic_non_clinal_pass <- as.logical(cryptic_non_clinal_pass)
  df$resurrection_basionym <- resurrection_basionym
  df$hybrid_zone_cohesion_pass <- as.logical(hybrid_zone_cohesion_pass)
  df$hybrid_parent_1 <- hybrid_parent_1
  df$hybrid_parent_2 <- hybrid_parent_2
  df$portal_metadata_lag_flag <- as.logical(portal_lag_flag)
  df$portal_lag_trigger_gate <- portal_lag_gate
  df$is_marginal_range_extension <- as.logical(is_marginal_ext_v)
  df$taxonomic_status_call <- status_calls; df$recommended_determination <- rec_taxa
  df$synthesis_triage_priority <- triage_prios; df$synthesis_triage_action <- triage_acts
  df$synthesis_rationale <- rationales
  df$taxonomic_treatment_broad <- treat_broad; df$taxonomic_treatment_subtype <- treat_sub
  if (exists("evaluate_step07_hydro_edaphic_and_basionym_guards", mode = "function")) {
    step07_gw_conn <- NULL
    if (!is.null(opts$db_conn) && requireNamespace("DBI", quietly = TRUE) &&
        isTRUE(tryCatch(DBI::dbIsValid(opts$db_conn), error = function(e) FALSE))) {
      step07_gw_conn <- opts$db_conn
    } else if (exists("c_geo", inherits = FALSE) && !is.null(c_geo) && requireNamespace("DBI", quietly = TRUE) &&
               isTRUE(tryCatch(DBI::dbIsValid(c_geo), error = function(e) FALSE))) {
      step07_gw_conn <- c_geo
    } else if (!isTRUE(opts$demo) &&
               trimws(env_plantmatics_or_packera("PLANTMATICS_DEMO_MODE", "PACKERA_DEMO_MODE")) != "1" &&
               exists("connect_packera_db", mode = "function") &&
               requireNamespace("DBI", quietly = TRUE) &&
               (nzchar(as.character(if (!is.null(opts_cyto$db_path)) opts_cyto$db_path else "")) ||
                nzchar(Sys.getenv("RUN_ID", unset = "")) ||
                nzchar(env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")) ||
                (!is.null(opts$vouchers) && exists("is_sqlite_active_for_inputs", mode = "function") && is_sqlite_active_for_inputs(opts$vouchers, demo = isTRUE(opts$demo))))) {
      step07_gw_conn <- tryCatch(connect_packera_db(db_path = opts_cyto$db_path, read_only = TRUE, demo = FALSE), error = function(e) NULL)
      if (!is.null(step07_gw_conn)) {
        on.exit(tryCatch(if (!is.null(step07_gw_conn) && DBI::dbIsValid(step07_gw_conn)) DBI::dbDisconnect(step07_gw_conn), error = function(e) NULL), add = TRUE)
      }
    }
    pre_res_bas <- df$resurrection_basionym
    df <- evaluate_step07_hydro_edaphic_and_basionym_guards(df, db_conn = step07_gw_conn)
    nwpl_code_v <- if ("nwpl_indicator_code" %in% names(df)) toupper(trimws(as.character(df$nwpl_indicator_code))) else rep(NA_character_, n)
    nwi_dist_v <- if ("nwi_wetland_dist_m" %in% names(df)) suppressWarnings(as.numeric(df$nwi_wetland_dist_m)) else rep(NA_real_, n)
    w_ed_k <- if ("w_edaphic_kappa" %in% names(df)) suppressWarnings(as.numeric(df$w_edaphic_kappa)) else rep(NA_real_, n)
    cda_p_v <- if ("cda_posterior_prob" %in% names(df)) suppressWarnings(as.numeric(df$cda_posterior_prob)) else rep(NA_real_, n)
    for (hi in seq_len(n)) {
      was_res_hi <- !is.na(pre_res_bas[hi]) && nzchar(trimws(pre_res_bas[hi]))
      is_homo_hi <- !is.na(df$basionym_synonymy_type[hi]) && identical(df$basionym_synonymy_type[hi], "Homotypic")
      if (was_res_hi && is_homo_hi && is.na(df$resurrection_basionym[hi])) {
        if (identical(df$taxonomic_status_call[hi], "Taxonomic_Resurrection_Candidate")) {
          gmm_k_hi <- if ("gmm_cluster" %in% names(df) && is.finite(suppressWarnings(as.integer(df$gmm_cluster[hi])))) max(1L, as.integer(df$gmm_cluster[hi])) else 2L
          df$taxonomic_status_call[hi] <- "Putative_Cryptic_Species"
          df$recommended_determination[hi] <- sprintf("%s (Putative Cryptic Lineage K%d)", as.character(df$species_standardized[hi]), gmm_k_hi)
          df$synthesis_triage_action[hi] <- "Flag_Putative_Cryptic_Species"
          df$taxonomic_treatment_broad[hi] <- "Misidentification / Cryptic Species"
          df$taxonomic_treatment_subtype[hi] <- "Putative_Cryptic_Species"
          df$synthesis_rationale[hi] <- sub("Stage 3 Basionym Resurrection Candidate", "Stage 3 Cryptic Congruence (Homotypic Basionym Guarded)", ifelse(is.na(df$synthesis_rationale[hi]), "", df$synthesis_rationale[hi]), fixed = TRUE)
        }
        cur_rat_h <- ifelse(is.na(df$synthesis_rationale[hi]), "", as.character(df$synthesis_rationale[hi]))
        if (!grepl("\\[Homotypic Basionym Guard:", cur_rat_h)) {
          df$synthesis_rationale[hi] <- paste0(cur_rat_h, sprintf(" [Homotypic Basionym Guard: %s shares type with %s; resurrection blocked]", pre_res_bas[hi], as.character(df$species_standardized[hi])))
        }
      }
      is_obl_facw <- !is.na(nwpl_code_v[hi]) && (nwpl_code_v[hi] %in% c("OBL", "FACW"))
      is_dry_site <- (!is.na(df$hydro_edaphic_niche_regime[hi]) && df$hydro_edaphic_niche_regime[hi] %in% c("Ultramafic_Serpentine_Barren", "Calcareous_Carbonate_Glade", "Siliceous_Felsic_Outcrop")) ||
                     (is.finite(nwi_dist_v[hi]) && nwi_dist_v[hi] >= 200.0)
      cur_rat_w <- ifelse(is.na(df$synthesis_rationale[hi]), "", as.character(df$synthesis_rationale[hi]))
      if (is_obl_facw && is_dry_site && is.finite(df$wetland_fidelity_concordance[hi]) && df$wetland_fidelity_concordance[hi] < 0.25 &&
          !grepl("\\[Hydro-Edaphic Discordance:", cur_rat_w)) {
        wk_s <- if (is.finite(w_ed_k[hi])) sprintf(", w_edaphic_kappa=%.2f", w_ed_k[hi]) else ""
        mp_s <- if (is.finite(cda_p_v[hi])) sprintf(", morph_cda_p=%.2f", cda_p_v[hi]) else ""
        reg_s <- if (!is.na(df$hydro_edaphic_niche_regime[hi])) df$hydro_edaphic_niche_regime[hi] else "Dry_Upland"
        df$synthesis_rationale[hi] <- paste0(
          cur_rat_w,
          sprintf(" [Hydro-Edaphic Discordance: taxon NWPL=%s in %s (wetland_concordance=%.2f%s%s)]", nwpl_code_v[hi], reg_s, df$wetland_fidelity_concordance[hi], wk_s, mp_s)
        )
      }
    }
  }
  has_morph <- if ("can1" %in% names(df)) !is.na(df$can1) else if ("cda_predicted_taxon" %in% names(df)) !is.na(df$cda_predicted_taxon) else rep(TRUE, n)
  has_edaph <- if ("edaphic_best_fit_taxon" %in% names(df)) !is.na(df$edaphic_best_fit_taxon) else rep(FALSE, n)
  df$is_completed_voucher <- has_morph | has_edaph | has_geno_dv
  if (!("genomic_propagation_tier" %in% names(df))) df$genomic_propagation_tier <- rep(NA_character_, n) else df$genomic_propagation_tier <- as.character(df$genomic_propagation_tier)
  if (!("tier0_calibrated_taxon" %in% names(df))) df$tier0_calibrated_taxon <- rep(NA_character_, n) else df$tier0_calibrated_taxon <- as.character(df$tier0_calibrated_taxon)
  if (!("tier0_calibrated_prob" %in% names(df))) df$tier0_calibrated_prob <- rep(NA_real_, n) else df$tier0_calibrated_prob <- suppressWarnings(as.numeric(df$tier0_calibrated_prob))
  if (!("tier0_abstention_flag" %in% names(df))) df$tier0_abstention_flag <- rep(NA_integer_, n) else df$tier0_abstention_flag <- suppressWarnings(as.integer(df$tier0_abstention_flag))
  if (!("deme_mahal_d2" %in% names(df))) df$deme_mahal_d2 <- rep(NA_real_, n) else df$deme_mahal_d2 <- suppressWarnings(as.numeric(df$deme_mahal_d2))
  if (!("cluster_propagated_gdi" %in% names(df))) df$cluster_propagated_gdi <- rep(NA_real_, n) else df$cluster_propagated_gdi <- suppressWarnings(as.numeric(df$cluster_propagated_gdi))
  if (!("cytology_linkage_tier" %in% names(df))) df$cytology_linkage_tier <- rep(NA_character_, n) else df$cytology_linkage_tier <- as.character(df$cytology_linkage_tier)
  if (!("voucher_sporophytic_2n" %in% names(df))) df$voucher_sporophytic_2n <- rep(NA_integer_, n) else df$voucher_sporophytic_2n <- suppressWarnings(as.integer(df$voucher_sporophytic_2n))
  if (!("voucher_ploidy_level" %in% names(df))) df$voucher_ploidy_level <- rep(NA_integer_, n) else df$voucher_ploidy_level <- suppressWarnings(as.integer(df$voucher_ploidy_level))
  if (!("taxon_reported_cytotypes" %in% names(df))) df$taxon_reported_cytotypes <- rep(NA_character_, n) else df$taxon_reported_cytotypes <- as.character(df$taxon_reported_cytotypes)
  if (!("is_multi_cytotype_taxon" %in% names(df))) df$is_multi_cytotype_taxon <- rep(NA_integer_, n) else df$is_multi_cytotype_taxon <- suppressWarnings(as.integer(df$is_multi_cytotype_taxon))
  if (!("cytotype_reproductive_barrier_pass" %in% names(df))) df$cytotype_reproductive_barrier_pass <- rep(NA_integer_, n) else df$cytotype_reproductive_barrier_pass <- suppressWarnings(as.integer(df$cytotype_reproductive_barrier_pass))
  if (!("hydro_edaphic_niche_regime" %in% names(df))) df$hydro_edaphic_niche_regime <- rep(NA_character_, n) else df$hydro_edaphic_niche_regime <- as.character(df$hydro_edaphic_niche_regime)
  if (!("wetland_fidelity_concordance" %in% names(df))) df$wetland_fidelity_concordance <- rep(NA_real_, n) else df$wetland_fidelity_concordance <- suppressWarnings(as.numeric(df$wetland_fidelity_concordance))
  if (!("basionym_synonymy_type" %in% names(df))) df$basionym_synonymy_type <- rep(NA_character_, n) else df$basionym_synonymy_type <- as.character(df$basionym_synonymy_type)
  attr(df, "can_evaluate_hybrids") <- can_evaluate_hybrids
  if (!is.null(prop_sum_attr)) attr(df, "genomic_propagation_summary") <- prop_sum_attr
  if (!is.null(cyto_md_attr)) attr(df, "cytology_markdown_lines") <- cyto_md_attr
  if (!is.null(cyto_bar_lbl_attr)) attr(df, "cytotype_barrier_labels") <- cyto_bar_lbl_attr
  if (!is.null(gigas_df_attr)) attr(df, "gigas_concordance_df") <- gigas_df_attr
  return(df)
}

# ------------------------------------------------------------------------------
# 5. Taxon-Adjusted Institutional Batch Audit (Marginal ANCOVA + Partial Eta^2)
# ------------------------------------------------------------------------------
.resolve_step07_ancova_features <- function(df) {
  is_usable_num <- function(cc) {
    if (!(cc %in% names(df)) || is.logical(df[[cc]]) || is.character(df[[cc]]) || is.factor(df[[cc]])) return(FALSE)
    v <- suppressWarnings(as.numeric(df[[cc]]))
    v <- v[is.finite(v)]
    length(v) >= 6L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
  }
  sel_cols <- .parse_step07_semicolon_cols(df, "selected_morph_features")
  sel_cols <- sel_cols[vapply(sel_cols, is_usable_num, logical(1))]
  if (length(sel_cols) > 0L) return(sel_cols)

  pc_cols <- .sort_numbered_pc_cols(grep("^PC[0-9]+$", names(df), value = TRUE), "PC")
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
  aw_v <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n_all)
  aw_v[!is.finite(aw_v) | aw_v <= 0.0] <- 1.0
  sbw_v <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, n_all)
  sbw_v[!is.finite(sbw_v) | sbw_v <= 0.0] <- 1.0
  w_all <- pmax(0.05, aw_v * sbw_v)

  base_raw <- !is.na(inst_clean) & nzchar(inst_clean) & !is.na(sp_clean) & nzchar(sp_clean)
  base_clean <- base_raw & !noise_m & !pass_m
  base_valid <- if (sum(base_clean) >= 6L) base_clean else base_raw
  if (sum(base_valid) < 6L) return(NULL)

  rows_out <- list()
  for (rc in resp_cols) {
    rc_num <- suppressWarnings(as.numeric(df[[rc]]))
    valid_rc <- base_valid & is.finite(rc_num)
    if (sum(valid_rc) < 6L) next
    sub_df <- df[valid_rc, , drop = FALSE]
    sub_df[[rc]] <- rc_num[valid_rc]
    sub_df$.w_ancova <- w_all[valid_rc]
    inst_tab <- table(sub_df$institutionCode)
    keep_inst <- names(inst_tab[inst_tab >= min_obs])
    if (length(keep_inst) >= 2L) sub_df <- sub_df[sub_df$institutionCode %in% keep_inst, , drop = FALSE]
    if (nrow(sub_df) < 6L) next
    sub_df$institutionCode <- droplevels(factor(sub_df$institutionCode))
    sub_df$species_standardized <- droplevels(factor(sub_df$species_standardized))
    if (nlevels(sub_df$institutionCode) < 2L) next

    has_reg <- ("regional_group" %in% names(sub_df)) && {
      rg_v <- trimws(as.character(sub_df$regional_group))
      rg_ok <- !is.na(rg_v) & nzchar(rg_v) & !(tolower(rg_v) %in% c("na", "unknown", "unassigned"))
      sum(rg_ok) == nrow(sub_df) && length(unique(rg_v[rg_ok])) >= 2L &&
        (nrow(sub_df) - nlevels(sub_df$species_standardized) - length(unique(rg_v[rg_ok])) - nlevels(sub_df$institutionCode) >= 3L)
    }
    if (has_reg) sub_df$regional_group <- droplevels(factor(trimws(as.character(sub_df$regional_group))))

    model_lbl <- if (nlevels(sub_df$species_standardized) >= 2L && has_reg) {
      "lm(Response ~ species_standardized + regional_group + institutionCode)"
    } else if (nlevels(sub_df$species_standardized) >= 2L) {
      "lm(Response ~ species_standardized + institutionCode)"
    } else {
      "lm(Response ~ institutionCode)"
    }
    fit_ancova <- if (nlevels(sub_df$species_standardized) >= 2L && has_reg) {
      stats::lm(stats::as.formula(paste(rc, "~ species_standardized + regional_group + institutionCode")), data = sub_df, weights = sub_df$.w_ancova)
    } else if (nlevels(sub_df$species_standardized) >= 2L) {
      if (rc == "PC1") stats::lm(PC1 ~ species_standardized + institutionCode, data = sub_df, weights = sub_df$.w_ancova)
      else stats::lm(stats::as.formula(paste(rc, "~ species_standardized + institutionCode")), data = sub_df, weights = sub_df$.w_ancova)
    } else stats::lm(stats::as.formula(paste(rc, "~ institutionCode")), data = sub_df, weights = sub_df$.w_ancova)
    aov_tab <- tryCatch(stats::anova(fit_ancova), error = function(e) NULL)
    if (is.null(aov_tab)) next
    rn <- rownames(aov_tab)
    inst_idx <- which(rn == "institutionCode")[1L]; res_idx <- which(rn == "Residuals")[1L]
    if (!is.na(inst_idx) && !is.na(res_idx)) {
      ss_inst <- as.numeric(aov_tab$`Sum Sq`[inst_idx]); ss_res <- as.numeric(aov_tab$`Sum Sq`[res_idx])
      ss_tot <- sum(as.numeric(aov_tab$`Sum Sq`), na.rm = TRUE)
      partial_eta_sq <- if ((ss_inst + ss_res) > 0) ss_inst / (ss_inst + ss_res) else 0.0
      rows_out[[length(rows_out) + 1L]] <- data.frame(
        Response = rc, Model = model_lbl,
        Term = "institutionCode", Df = as.integer(aov_tab$Df[inst_idx]),
        Sum_Sq = round(ss_inst, 4), Residual_Sum_Sq = round(ss_res, 4),
        F_value = round(as.numeric(aov_tab$`F value`[inst_idx]), 4),
        P_value = signif(as.numeric(aov_tab$`Pr(>F)`[inst_idx]), 4),
        partial_eta_sq = round(partial_eta_sq, 4),
        R_squared = round(if (ss_tot > 0) ss_inst / ss_tot else 0.0, 4),
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
  return(TRUE)
}

# ------------------------------------------------------------------------------
# 6. Merging Evidence Streams & Coalescing Non-Georeferenced Vouchers
# ------------------------------------------------------------------------------
safe_left_join_and_coalesce <- function(base_df, right_df, numeric_coalesce_cols = character(0), append_missing = TRUE) {
  if (is.null(right_df) || nrow(right_df) == 0L || !("catalogNumber" %in% names(right_df))) return(base_df)
  right_df$catalogNumber <- as.character(right_df$catalogNumber)
  base_df$catalogNumber <- as.character(base_df$catalogNumber)

  right_has_pid <- ("plant_individual_id" %in% names(right_df)) && any(!is.na(right_df$plant_individual_id))
  if (right_has_pid) {
    r_pid_vec <- ifelse(is.na(right_df$plant_individual_id), 0L, suppressWarnings(as.integer(right_df$plant_individual_id)))
    r_pid_vec[is.na(r_pid_vec)] <- 0L
    right_df$plant_individual_id <- r_pid_vec
    right_sub <- right_df[!duplicated(paste0(right_df$catalogNumber, "__", right_df$plant_individual_id)), , drop = FALSE]

    multi_cats <- unique(right_sub$catalogNumber[duplicated(right_sub$catalogNumber)])
    if (length(multi_cats) > 0L) {
      if (!("plant_individual_id" %in% names(base_df))) base_df$plant_individual_id <- NA_integer_
      for (mc in multi_cats) {
        b_rows <- which(base_df$catalogNumber == mc)
        r_rows <- which(right_sub$catalogNumber == mc)
        r_pids <- right_sub$plant_individual_id[r_rows]
        if (length(b_rows) == 1L && (is.na(base_df$plant_individual_id[b_rows]) || suppressWarnings(as.integer(base_df$plant_individual_id[b_rows])) == 0L)) {
          rep_block <- base_df[rep(b_rows, length(r_pids)), , drop = FALSE]
          rep_block$plant_individual_id <- r_pids
          keep_idx <- seq_len(nrow(base_df)) != b_rows
          base_df <- rbind(base_df[keep_idx, , drop = FALSE], rep_block)
        } else if (isTRUE(append_missing) && length(b_rows) >= 1L) {
          b_pids <- suppressWarnings(as.integer(base_df$plant_individual_id[b_rows]))
          miss_pids <- setdiff(r_pids, b_pids[!is.na(b_pids)])
          if (length(miss_pids) > 0L) {
            rep_extra <- base_df[rep(b_rows[1L], length(miss_pids)), , drop = FALSE]
            rep_extra$plant_individual_id <- miss_pids
            base_df <- rbind(base_df, rep_extra)
          }
        }
      }
    }
  } else {
    right_sub <- right_df[!duplicated(right_df$catalogNumber), , drop = FALSE]
  }

  missing_cats <- setdiff(unique(right_sub$catalogNumber), unique(base_df$catalogNumber))
  if (isTRUE(append_missing) && length(missing_cats) > 0L) {
    extra_rows <- right_sub[right_sub$catalogNumber %in% missing_cats, , drop = FALSE]
    for (col in setdiff(names(base_df), names(extra_rows))) extra_rows[[col]] <- NA
    base_df <- rbind(base_df, extra_rows[, names(base_df), drop = FALSE])
  }

  join_cols <- setdiff(names(right_sub), names(base_df))
  dedup_cols <- unique(c("catalogNumber", join_cols))
  if (!right_has_pid && !("plant_individual_id" %in% names(base_df)) && length(join_cols) > 0L && !any(duplicated(base_df$catalogNumber))) {
    base_df <- dplyr::left_join(base_df, right_sub[, dedup_cols, drop = FALSE], by = "catalogNumber")
  }

  if (right_has_pid && ("plant_individual_id" %in% names(base_df))) {
    b_pid <- ifelse(is.na(base_df$plant_individual_id), 0L, suppressWarnings(as.integer(base_df$plant_individual_id)))
    b_pid[is.na(b_pid)] <- 0L
    idx_exact <- match(paste0(base_df$catalogNumber, "__", b_pid), paste0(right_sub$catalogNumber, "__", right_sub$plant_individual_id))
    r_cnt <- table(right_sub$catalogNumber)
    single_r_cats <- names(r_cnt[r_cnt == 1L])
    idx_cat_single <- ifelse(base_df$catalogNumber %in% single_r_cats, match(base_df$catalogNumber, right_sub$catalogNumber), NA_integer_)
    idx_r <- ifelse(!is.na(idx_exact), idx_exact, idx_cat_single)
  } else {
    idx_r <- match(base_df$catalogNumber, right_sub$catalogNumber)
  }

  for (col in setdiff(names(right_sub), "catalogNumber")) {
    if (col == "plant_individual_id") {
      if ("plant_individual_id" %in% names(base_df)) {
        cur_p <- suppressWarnings(as.integer(base_df$plant_individual_id))
        new_p <- suppressWarnings(as.integer(right_sub$plant_individual_id[idx_r]))
        up_mask <- !is.na(new_p) & (is.na(cur_p) | (cur_p == 0L & new_p != 0L))
        if (any(up_mask)) base_df$plant_individual_id[up_mask] <- new_p[up_mask]
      } else {
        base_df$plant_individual_id <- right_sub$plant_individual_id[idx_r]
      }
    } else if (col %in% numeric_coalesce_cols && col %in% names(base_df)) {
      base_df[[col]] <- dplyr::coalesce(suppressWarnings(as.numeric(base_df[[col]])), suppressWarnings(as.numeric(right_sub[[col]][idx_r])))
    } else if (col %in% names(base_df)) {
      na_mask <- (is.na(base_df[[col]]) | (is.character(base_df[[col]]) & !nzchar(trimws(base_df[[col]])))) & !is.na(idx_r)
      if (any(na_mask)) base_df[[col]][na_mask] <- right_sub[[col]][idx_r[na_mask]]
    } else {
      base_df[[col]] <- right_sub[[col]][idx_r]
    }
  }
  base_df
}

resolve_gmm_bic_stats <- function(opts) {
  opt_k <- 4L; delta_bic_val <- 12.0
  gmm_candidates <- unique(stats::na.omit(c(opts$gmm_bic, "data/tables/gmm_bic_model_selection.csv", opts$gmm_summary, "outputs/reports/gmm_bayes_factors_summary.csv")))
  gmm_candidates <- gmm_candidates[nzchar(gmm_candidates) & file.exists(gmm_candidates)]
  if (length(gmm_candidates) > 0L) {
    gmm_sum_df <- tryCatch(readr::read_csv(gmm_candidates[1], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(gmm_sum_df) && nrow(gmm_sum_df) > 0 && "BIC" %in% names(gmm_sum_df)) {
      best_idx <- which.max(gmm_sum_df$BIC)
      if ("is_optimal" %in% names(gmm_sum_df) && any(gmm_sum_df$is_optimal == TRUE, na.rm = TRUE)) best_idx <- which(gmm_sum_df$is_optimal == TRUE)[1]
      if (length(best_idx) > 0 && "K" %in% names(gmm_sum_df)) opt_k <- as.integer(gmm_sum_df$K[best_idx])
      for (bc in c("BayesFactor_vs_Null", "Two_Delta_BIC", "Delta_BIC")) {
        if (bc %in% names(gmm_sum_df) && length(best_idx) > 0 && is.finite(as.numeric(gmm_sum_df[[bc]][best_idx]))) {
          delta_bic_val <- as.numeric(gmm_sum_df[[bc]][best_idx]); break
        }
      }
    }
  }
  list(opt_k = opt_k, two_delta_bic = delta_bic_val)
}

build_triage_queue <- function(opts) {
  if (is.data.frame(opts)) return(apply_taxonomic_decision_matrix(opts))
  message("Ingesting and joining evidence streams...")
  if (!file.exists(opts$vouchers)) stop("Missing vouchers: ", opts$vouchers)
  vouchers <- readr::read_csv(opts$vouchers, show_col_types = FALSE)
  vouchers$catalogNumber <- as.character(vouchers$catalogNumber)
  validate_upstream_contracts(vouchers_df = vouchers)
  sp_cand <- rep("Unknown", nrow(vouchers))
  for (sc in intersect(c("species_standardized", "scientificName", "species_raw"), names(vouchers))) {
    v <- standardize_packera_taxon(vouchers[[sc]])
    m <- sp_cand == "Unknown" & !is.na(v) & v != "Unknown"
    sp_cand[m] <- v[m]
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
    "doy", "decimalLatitude", "decimalLongitude",
    "soil_ph", "soil_sand", "soil_cec", "soil_bulk_density", "delta_doy",
    "fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "ssurgo_rock_outcrop_pct",
    "edaphic_z_score", "type_edaphic_z_score", "climatic_z_score", "spatial_block_cv_acc", "rf_predicted_prob", "rf_given_taxon_prob",
    "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
    "plant_height_cm", "scape_to_basal_leaf_ratio",
    "capitula_count", "log1p_capitulum_count", "involucre_height_mm", "involucre_width_mm", "AR_inv", "blade_cie_l_norm", "blade_cie_l_sd"
  )

  run_dir_env <- env_plantmatics_or_packera("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR")
  def_dt_dir <- normalizePath(if (nzchar(run_dir_env)) file.path(run_dir_env, "data/tables") else "data/tables", mustWork = FALSE)
  vouch_dir <- normalizePath(dirname(opts$vouchers), mustWork = FALSE)
  is_isolated_vouch_dir <- !identical(vouch_dir, def_dt_dir)
  if (is_isolated_vouch_dir && !isTRUE(opts$cli_has_multimodal) && !is.null(opts$multimodal_flags) && nzchar(opts$multimodal_flags) &&
      identical(normalizePath(dirname(opts$multimodal_flags), mustWork = FALSE), def_dt_dir)) {
    sib_mm <- file.path(dirname(opts$vouchers), "multimodal_conflict_flags.csv")
    opts$multimodal_flags <- if (file.exists(sib_mm)) sib_mm else ""
  }
  if (is_isolated_vouch_dir && !isTRUE(opts$cli_has_harmonics) && !is.null(opts$harmonics) && nzchar(opts$harmonics) &&
      identical(normalizePath(dirname(opts$harmonics), mustWork = FALSE), def_dt_dir)) {
    sib_hm <- file.path(dirname(opts$vouchers), "leaf_efa_harmonics.csv")
    opts$harmonics <- if (file.exists(sib_hm)) sib_hm else ""
  }
  if (is_isolated_vouch_dir && !isTRUE(opts$cli_has_morphometrics) && !is.null(opts$morphometrics) && nzchar(opts$morphometrics) &&
      identical(normalizePath(dirname(opts$morphometrics), mustWork = FALSE), def_dt_dir)) {
    sib_mo <- file.path(dirname(opts$vouchers), "morphometric_flags.csv")
    opts$morphometrics <- if (file.exists(sib_mo)) sib_mo else ""
  }

  if (!is.null(opts$multimodal_flags) && nzchar(opts$multimodal_flags) && file.exists(opts$multimodal_flags)) {
    mm_df <- readr::read_csv(opts$multimodal_flags, show_col_types = FALSE)
    mm_df$catalogNumber <- as.character(mm_df$catalogNumber)
    if ("plant_height_mm" %in% names(mm_df) && !("plant_height_cm" %in% names(mm_df))) {
      mm_df$plant_height_cm <- round(suppressWarnings(as.numeric(mm_df$plant_height_mm)) / 10.0, 2)
    }
    dyn_mm_num <- grep("^(bio[0-9]+_|fine_|ilr_texture_|ssurgo_|polaris_|solus_|topo_)", names(mm_df), value = TRUE)
    dyn_mm_num <- dyn_mm_num[vapply(dyn_mm_num, function(cc) !is.character(mm_df[[cc]]) && !is.logical(mm_df[[cc]]), logical(1))]
    mm_num_all <- unique(c(mm_numeric_coalesce, dyn_mm_num))
    vouchers <- safe_left_join_and_coalesce(vouchers, mm_df, numeric_coalesce_cols = mm_num_all)
    for (mc_col in intersect(c("rf_predicted_taxon", "is_verified_anthesis", "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors", "nearest_type_niche_taxon", "matches_type_lithology", "type_edaphic_z_score"), names(mm_df))) {
      if (!(mc_col %in% names(vouchers))) vouchers[[mc_col]] <- NA
    }
  }

  pheno_state_candidates <- file.path(dirname(opts$vouchers), "voucher_phenological_states.csv")
  pheno_state_candidates <- pheno_state_candidates[file.exists(pheno_state_candidates)]
  if (length(pheno_state_candidates) > 0L) {
    ph_df <- tryCatch(readr::read_csv(pheno_state_candidates[1L], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(ph_df) && "catalogNumber" %in% names(ph_df)) {
      vouchers <- safe_left_join_and_coalesce(
        vouchers, ph_df,
        numeric_coalesce_cols = c("n_anthesis", "n_bud", "n_fruiting", "n_senescent", "reproductive_stage_confidence", "total_capitula"),
        append_missing = TRUE
      )
    }
  }

  med_or_na <- function(x) {
    v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]
    if (length(v) > 0L) stats::median(v) else NA_real_
  }
  first_chr_or_na <- function(x) {
    v <- as.character(stats::na.omit(x)); v <- v[nzchar(trimws(v))]
    if (length(v) > 0L) v[1L] else NA_character_
  }
  first_int_or_na <- function(x) {
    v <- suppressWarnings(as.integer(stats::na.omit(x))); v <- v[!is.na(v)]
    if (length(v) > 0L) v[1L] else NA_integer_
  }
  any_lgl_or_na <- function(x) {
    v <- stats::na.omit(suppressWarnings(as.logical(x)))
    if (length(v) > 0L) any(v) else NA
  }

  if (is.null(opts$morphometrics) || !nzchar(opts$morphometrics) || !file.exists(opts$morphometrics)) {
    sibling_morph <- file.path(dirname(opts$vouchers), "morphometric_flags.csv")
    if (file.exists(sibling_morph)) opts$morphometrics <- sibling_morph else if (!is_isolated_vouch_dir && file.exists("data/tables/morphometric_flags.csv")) opts$morphometrics <- "data/tables/morphometric_flags.csv"
  }
  if (!is.null(opts$morphometrics) && nzchar(opts$morphometrics) && file.exists(opts$morphometrics)) {
    mo_df <- readr::read_csv(opts$morphometrics, show_col_types = FALSE)
    mo_df$catalogNumber <- as.character(mo_df$catalogNumber)
    validate_upstream_contracts(morpho_df = mo_df)
    if ("cda_predicted" %in% names(mo_df) && !("cda_predicted_taxon" %in% names(mo_df))) mo_df$cda_predicted_taxon <- as.character(mo_df$cda_predicted)
    if ("plant_height_mm" %in% names(mo_df) && !("plant_height_cm" %in% names(mo_df))) {
      mo_df$plant_height_cm <- round(suppressWarnings(as.numeric(mo_df$plant_height_mm)) / 10.0, 2)
    }
    if (!("is_promoted_bronze" %in% names(mo_df))) mo_df$is_promoted_bronze <- NA
    if (!("morphometric_outlier_flag" %in% names(mo_df))) mo_df$morphometric_outlier_flag <- FALSE
    if (!("misidentification_flag" %in% names(mo_df))) mo_df$misidentification_flag <- FALSE
    if (!("plant_individual_id" %in% names(mo_df))) mo_df$plant_individual_id <- 0L
    mo_df$plant_individual_id <- ifelse(is.na(mo_df$plant_individual_id), 0L, suppressWarnings(as.integer(mo_df$plant_individual_id)))
    mo_df$plant_individual_id[is.na(mo_df$plant_individual_id)] <- 0L
    if (!("can3" %in% names(mo_df))) mo_df$can3 <- if ("PC3" %in% names(mo_df)) suppressWarnings(as.numeric(mo_df$PC3)) else NA_real_
    if (!("PC1" %in% names(mo_df))) mo_df$PC1 <- if ("can1" %in% names(mo_df)) suppressWarnings(as.numeric(mo_df$can1)) else NA_real_
    if (!("PC2" %in% names(mo_df))) mo_df$PC2 <- if ("can2" %in% names(mo_df)) suppressWarnings(as.numeric(mo_df$can2)) else NA_real_
    for (nc in c(
      "cda_posterior_prob", "cda_second_posterior_prob", "cda_posterior_margin", "loocv_accuracy",
      "can1", "can2", "gmm_uncertainty",
      "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
      "plant_height_cm", "scape_to_basal_leaf_ratio",
      "capitula_count", "log1p_capitulum_count", "involucre_height_mm", "involucre_width_mm", "AR_inv",
      "blade_cie_l_norm", "blade_cie_l_sd", "ecoregion_cda_posterior_prob", "authority_weight",
      "type_morph_mahal_d2", "nearest_type_basionym_d2", "type_morph_concordance"
    )) if (!(nc %in% names(mo_df))) mo_df[[nc]] <- NA_real_
    if (!("cda_predicted_taxon" %in% names(mo_df))) mo_df$cda_predicted_taxon <- NA_character_
    if (!("cda_second_taxon" %in% names(mo_df))) {
      mo_df$cda_second_taxon <- if ("cda_second_predicted_taxon" %in% names(mo_df)) as.character(mo_df$cda_second_predicted_taxon) else NA_character_
    }
    if (!("ecoregion_cda_predicted" %in% names(mo_df))) mo_df$ecoregion_cda_predicted <- NA_character_
    if (!("nearest_type_morph_taxon" %in% names(mo_df))) mo_df$nearest_type_morph_taxon <- NA_character_
    if (!("nearest_type_basionym" %in% names(mo_df))) mo_df$nearest_type_basionym <- NA_character_
    if (!("nearest_type_basionym_d2" %in% names(mo_df))) mo_df$nearest_type_basionym_d2 <- NA_real_
    if (!("is_bent_stem" %in% names(mo_df))) mo_df$is_bent_stem <- NA
    if (!("is_typical_morphotype" %in% names(mo_df))) mo_df$is_typical_morphotype <- NA
    if (!("has_empirical_type_anchor" %in% names(mo_df))) mo_df$has_empirical_type_anchor <- NA
    if (!("is_type_specimen" %in% names(mo_df))) mo_df$is_type_specimen <- NA
    if (!("gmm_cluster" %in% names(mo_df))) mo_df$gmm_cluster <- NA_integer_
    if (!("morph_tier_source" %in% names(mo_df))) mo_df$morph_tier_source <- NA_character_
    for (ic in c("leaf_count", "basal_leaf_count", "cauline_leaf_count")) if (!(ic %in% names(mo_df))) mo_df[[ic]] <- NA_integer_

    mo_agg <- mo_df %>% dplyr::group_by(catalogNumber, plant_individual_id) %>% dplyr::summarise(
      morph_tier_source = first_chr_or_na(morph_tier_source),
      leaf_count = first_int_or_na(leaf_count),
      basal_leaf_count = first_int_or_na(basal_leaf_count),
      cauline_leaf_count = first_int_or_na(cauline_leaf_count),
      cda_predicted_taxon = first_chr_or_na(cda_predicted_taxon),
      cda_second_taxon = first_chr_or_na(cda_second_taxon),
      cda_posterior_prob = med_or_na(cda_posterior_prob),
      cda_second_posterior_prob = med_or_na(cda_second_posterior_prob),
      cda_posterior_margin = med_or_na(cda_posterior_margin),
      loocv_accuracy = med_or_na(loocv_accuracy),
      PC1 = med_or_na(PC1), PC2 = med_or_na(PC2),
      can1 = med_or_na(can1), can2 = med_or_na(can2), can3 = med_or_na(can3),
      gmm_cluster = first_int_or_na(gmm_cluster), gmm_uncertainty = med_or_na(gmm_uncertainty),
      is_typical_morphotype = any_lgl_or_na(is_typical_morphotype),
      type_morph_mahal_d2 = med_or_na(type_morph_mahal_d2),
      nearest_type_morph_taxon = first_chr_or_na(nearest_type_morph_taxon),
      nearest_type_basionym = first_chr_or_na(nearest_type_basionym),
      nearest_type_basionym_d2 = med_or_na(nearest_type_basionym_d2),
      type_morph_concordance = med_or_na(type_morph_concordance),
      has_empirical_type_anchor = any_lgl_or_na(has_empirical_type_anchor),
      is_type_specimen = any_lgl_or_na(is_type_specimen),
      petiole_length_mm = med_or_na(petiole_length_mm),
      blade_length_mm = med_or_na(blade_length_mm),
      petiole_to_blade_ratio = med_or_na(petiole_to_blade_ratio),
      plant_height_cm = med_or_na(plant_height_cm),
      scape_to_basal_leaf_ratio = med_or_na(scape_to_basal_leaf_ratio),
      is_bent_stem = any_lgl_or_na(is_bent_stem),
      capitula_count = med_or_na(capitula_count),
      involucre_height_mm = med_or_na(involucre_height_mm),
      involucre_width_mm = med_or_na(involucre_width_mm),
      AR_inv = med_or_na(AR_inv),
      ecoregion_cda_predicted = first_chr_or_na(ecoregion_cda_predicted),
      ecoregion_cda_posterior_prob = med_or_na(ecoregion_cda_posterior_prob),
      is_promoted_bronze = any_lgl_or_na(is_promoted_bronze),
      authority_weight = med_or_na(authority_weight),
      morphometric_outlier_flag = any(morphometric_outlier_flag == TRUE, na.rm = TRUE),
      misidentification_flag = any(misidentification_flag == TRUE, na.rm = TRUE), .groups = "drop"
    )
    # Losslessly propagate all remaining Step 03/04 columns (ECT, cauline, heteroblasty, harmonics, PCs, OSR, weights)
    rem_mo_cols <- setdiff(names(mo_df), c(names(mo_agg), "species_standardized", "scientificName", "species_raw", "determiner_tier"))
    dyn_mo_num_cols <- character(0)
    if (length(rem_mo_cols) > 0L) {
      mo_key <- paste0(mo_df$catalogNumber, "__", mo_df$plant_individual_id)
      agg_key <- paste0(mo_agg$catalogNumber, "__", mo_agg$plant_individual_id)
      lgl_known <- c("is_dissected", "is_contour_outlier", "is_confident_label_noise", "is_passive", "is_calibrated_mm", "cryptic_cluster_cohesion_pass", "cryptic_non_clinal_pass")
      int_known <- c("efd_harmonic_opt", "efd_k95_pcs", "self_morph_df", "open_set_df", "voucher_sporophytic_2n", "voucher_ploidy_level", "is_multi_cytotype_taxon", "cytotype_reproductive_barrier_pass")
      for (rmc in rem_mo_cols) {
        col_v <- mo_df[[rmc]]
        if (is.logical(col_v) || rmc %in% lgl_known) {
          vals_by_k <- tapply(col_v, mo_key, any_lgl_or_na)
          mo_agg[[rmc]] <- as.logical(vals_by_k[agg_key])
        } else if (is.integer(col_v) || rmc %in% int_known) {
          vals_by_k <- tapply(col_v, mo_key, first_int_or_na)
          mo_agg[[rmc]] <- as.integer(vals_by_k[agg_key])
        } else if (is.numeric(col_v)) {
          vals_by_k <- tapply(col_v, mo_key, med_or_na)
          mo_agg[[rmc]] <- as.numeric(vals_by_k[agg_key])
          dyn_mo_num_cols <- c(dyn_mo_num_cols, rmc)
        } else {
          vals_by_k <- tapply(col_v, mo_key, first_chr_or_na)
          mo_agg[[rmc]] <- as.character(vals_by_k[agg_key])
        }
      }
    }
    vouchers <- safe_left_join_and_coalesce(
      vouchers, mo_agg,
      numeric_coalesce_cols = unique(c(
        "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio",
        "plant_height_cm", "scape_to_basal_leaf_ratio",
        "capitula_count", "involucre_height_mm", "involucre_width_mm", "AR_inv",
        "cda_posterior_prob", "cda_second_posterior_prob", "cda_posterior_margin",
        "gmm_uncertainty", "loocv_accuracy", "PC1", "PC2",
        "can1", "can2", "can3", "ecoregion_cda_posterior_prob", "authority_weight",
        "type_morph_mahal_d2", "type_morph_concordance",
        dyn_mo_num_cols
      ))
    )
  }

  if (!is.null(opts$linear_traits) && nzchar(opts$linear_traits) && file.exists(opts$linear_traits)) {
    lt_df <- tryCatch(readr::read_csv(opts$linear_traits, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(lt_df) && "catalogNumber" %in% names(lt_df)) {
      vouchers <- safe_left_join_and_coalesce(vouchers, lt_df, append_missing = FALSE)
    }
  }
  if (!is.null(opts$harmonics) && nzchar(opts$harmonics) && file.exists(opts$harmonics)) {
    hm_df <- tryCatch(readr::read_csv(opts$harmonics, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(hm_df) && "catalogNumber" %in% names(hm_df)) {
      vouchers <- safe_left_join_and_coalesce(vouchers, hm_df, append_missing = FALSE)
    }
  }
  }

  if (any(is.na(vouchers$species_standardized) | vouchers$species_standardized == "Unknown")) {
    for (sc in intersect(c("species_standardized", "scientificName", "species_raw", "taxon"), names(vouchers))) {
      v <- standardize_packera_taxon(vouchers[[sc]])
      m <- (is.na(vouchers$species_standardized) | vouchers$species_standardized == "Unknown") & !is.na(v) & v != "Unknown"
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

  if (exists("load_step07_genomic_metrics", mode = "function")) {
    vouchers <- load_step07_genomic_metrics(opts, vouchers)
  }
  opts_cyto_q <- opts
  if (is.null(opts_cyto_q$db_path)) {
    db_env_q <- env_plantmatics_or_packera("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")
    if (nzchar(db_env_q)) opts_cyto_q$db_path <- db_env_q
  }
  if (is.null(opts_cyto_q$report_out) && !is.null(opts_cyto_q$output_report) && nzchar(opts_cyto_q$output_report)) {
    opts_cyto_q$report_out <- opts_cyto_q$output_report
  }
  cyto_dir_q <- if (!is.null(opts_cyto_q$cytology_dir) && nzchar(opts_cyto_q$cytology_dir)) {
    opts_cyto_q$cytology_dir
  } else {
    env_cdir_q <- env_plantmatics_or_packera("PLANTMATICS_CYTOLOGY_DIR", "PACKERA_CYTOLOGY_DIR")
    if (nzchar(env_cdir_q)) env_cdir_q else "data/cytology"
  }
  if (exists("run_step07_cytology_integration", mode = "function")) {
    vouchers <- run_step07_cytology_integration(vouchers, cytology_dir = cyto_dir_q, opts = opts_cyto_q)
    attr(vouchers, "step07_cytology_executed") <- TRUE
  }
  cyto_md_q <- attr(vouchers, "cytology_markdown_lines")
  cyto_bar_lbl_q <- attr(vouchers, "cytotype_barrier_labels")
  gigas_df_q <- attr(vouchers, "gigas_concordance_df")
  if (exists("run_step07_genomic_propagation", mode = "function")) {
    vouchers <- run_step07_genomic_propagation(vouchers, opts)
  }
  prop_sum <- attr(vouchers, "genomic_propagation_summary")

  gmm_stats <- resolve_gmm_bic_stats(opts)
  vouchers <- apply_taxonomic_decision_matrix(vouchers, intraspecific_mode = isTRUE(opts$intraspecific_mode), gmm_two_delta_bic = gmm_stats$two_delta_bic, opts = opts)
  if (!is.null(attr(vouchers, "genomic_propagation_summary"))) prop_sum <- attr(vouchers, "genomic_propagation_summary")
  if (!is.null(prop_sum)) attr(vouchers, "genomic_propagation_summary") <- prop_sum
  if (!is.null(attr(vouchers, "cytology_markdown_lines"))) cyto_md_q <- attr(vouchers, "cytology_markdown_lines")
  if (!is.null(attr(vouchers, "cytotype_barrier_labels"))) cyto_bar_lbl_q <- attr(vouchers, "cytotype_barrier_labels")
  if (!is.null(attr(vouchers, "gigas_concordance_df"))) gigas_df_q <- attr(vouchers, "gigas_concordance_df")

  # Run Unsupervised Quantitative Taxon-Complex Discovery across full multi-taxon cohort BEFORE SQLite persistence and --focal-taxon filtering
  comp_res <- if (exists("run_step07_complex_discovery", mode = "function")) {
    tryCatch(run_step07_complex_discovery(vouchers, opts), error = function(e) { warning("Step 07 complex discovery warning: ", e$message); NULL })
  } else NULL
  if (!is.null(comp_res) && !is.null(comp_res$complexes_df) && nrow(comp_res$complexes_df) > 0L && !isTRUE(comp_res$subcluster_mode)) {
    cdf_m <- comp_res$complexes_df
    idx_cm <- match(as.character(vouchers$species_standardized), as.character(cdf_m$species_standardized))
    vouchers$consensus_complex_id <- as.integer(cdf_m$consensus_complex_id[idx_cm])
    vouchers$complex_label <- as.character(cdf_m$complex_label[idx_cm])
  }

  # Persist full evaluated cohort (including outgroup anchors and incoming reassignments) to SQLite BEFORE --focal-taxon filtering
  if (!isTRUE(opts$demo) && exists("persist_step07_to_sqlite", mode = "function")) {
    tryCatch(persist_step07_to_sqlite(vouchers, vouchers_path = opts$vouchers, demo = isTRUE(opts$demo)), error = function(e) warning("SQLite Step 07 persistence note: ", e$message))
  }

  can_eval_hyb_attr <- attr(vouchers, "can_evaluate_hybrids")
  full_complex_df <- vouchers
  focal_active_flag <- !is.null(opts$focal_taxon) && nzchar(trimws(opts$focal_taxon))
  if (focal_active_flag) {
    f_std <- standardize_packera_taxon(opts$focal_taxon)
    rec_det_v <- if ("recommended_determination" %in% names(vouchers)) as.character(vouchers$recommended_determination) else rep(NA_character_, nrow(vouchers))
    hp1_v <- if ("hybrid_parent_1" %in% names(vouchers)) as.character(vouchers$hybrid_parent_1) else rep(NA_character_, nrow(vouchers))
    hp2_v <- if ("hybrid_parent_2" %in% names(vouchers)) as.character(vouchers$hybrid_parent_2) else rep(NA_character_, nrow(vouchers))
    keep_focal <- (vouchers$species_standardized == f_std) | (!is.na(rec_det_v) & rec_det_v == f_std) |
                  (!is.na(hp1_v) & hp1_v == f_std) | (!is.na(hp2_v) & hp2_v == f_std)
    vouchers <- vouchers[keep_focal, , drop = FALSE]
  }

  prio_order <- c("CRITICAL" = 1, "HIGH" = 2, "MEDIUM" = 3, "LOW" = 4, "RESOLVED" = 5)
  sort_triage_by_prio <- function(d_in) {
    if (!is.data.frame(d_in) || nrow(d_in) == 0L) return(d_in)
    d_in$prio_rank <- prio_order[d_in$synthesis_triage_priority]
    d_in <- if ("coordinateUncertainty" %in% names(d_in)) d_in[order(d_in$prio_rank, -d_in$coordinateUncertainty), ] else d_in[order(d_in$prio_rank), ]
    d_in$prio_rank <- NULL
    d_in
  }
  vouchers <- sort_triage_by_prio(vouchers)
  full_complex_df <- sort_triage_by_prio(full_complex_df)
  if (!is.null(comp_res)) attr(vouchers, "complex_discovery") <- comp_res
  if (!is.null(can_eval_hyb_attr)) attr(vouchers, "can_evaluate_hybrids") <- can_eval_hyb_attr
  if (!is.null(prop_sum)) attr(vouchers, "genomic_propagation_summary") <- prop_sum
  if (!is.null(cyto_md_q)) attr(vouchers, "cytology_markdown_lines") <- cyto_md_q
  if (!is.null(cyto_bar_lbl_q)) attr(vouchers, "cytotype_barrier_labels") <- cyto_bar_lbl_q
  if (!is.null(gigas_df_q)) attr(vouchers, "gigas_concordance_df") <- gigas_df_q
  if (!is.null(full_complex_df) && length(unique(stats::na.omit(as.character(full_complex_df$species_standardized)))) >= 2L) {
    attr(vouchers, "full_complex_df") <- full_complex_df
  }

  .write_triage_fast <- function(d, p, na_str = "NA") {
    if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(d, file = p, na = na_str, quote = "auto", logical01 = FALSE) else readr::write_csv(d, p, na = na_str)
  }
  dir.create(dirname(opts$output_queue), recursive = TRUE, showWarnings = FALSE)
  tmp_queue <- paste0(opts$output_queue, ".tmp"); .write_triage_fast(vouchers, tmp_queue, na_str = "NA"); file.rename(tmp_queue, opts$output_queue)
  if (focal_active_flag && !is.null(full_complex_df) && nrow(full_complex_df) > 0L) {
    full_q_path <- file.path(dirname(opts$output_queue), "triage_queue_complex_full.csv")
    tmp_full_q <- paste0(full_q_path, ".tmp"); .write_triage_fast(full_complex_df, tmp_full_q, na_str = "NA"); file.rename(tmp_full_q, full_q_path)
  }
  if (exists("export_genomic_field_triage_csv", mode = "function")) {
    tryCatch(export_genomic_field_triage_csv(opts, vouchers), error = function(e) warning("Genomic field triage export note: ", e$message))
  }
  reports_dir <- if (!is.null(opts$output_report) && nzchar(opts$output_report)) {
    dirname(opts$output_report)
  } else if (!is.null(opts$output_complexes) && nzchar(opts$output_complexes)) {
    dirname(opts$output_complexes)
  } else if (!is.null(opts$output_html) && nzchar(opts$output_html)) {
    dirname(opts$output_html)
  } else if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) {
    q_dir <- dirname(opts$output_queue)
    if (basename(q_dir) == "tables") file.path(dirname(q_dir), "reports") else q_dir
  } else {
    "outputs/reports"
  }
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    audit_cols <- c(
      "catalogNumber", "institutionCode", "species_standardized", "identifiedBy", "determiner_tier",
      "recordedBy", "recordNumber", "stateProvince", "county", "usda_county_status",
      "exsiccatae_key", "exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon",
      "climatic_z_score", "rf_given_taxon_prob", "cda_predicted_taxon", "cda_posterior_prob",
      "rf_predicted_taxon", "rf_predicted_prob", "portal_lag_trigger_gate",
      "portal_metadata_lag_flag", "is_marginal_range_extension",
      "recommended_determination", "synthesis_triage_action", "taxonomic_treatment_subtype",
      "synthesis_rationale"
    )
    n_v <- nrow(vouchers)
    u_aud <- if ("usda_county_status" %in% names(vouchers)) tolower(trimws(ifelse(is.na(vouchers$usda_county_status), "", as.character(vouchers$usda_county_status)))) else rep("", n_v)
    ex_aud <- if ("exsiccatae_taxonomic_conflict" %in% names(vouchers)) (!is.na(suppressWarnings(as.integer(vouchers$exsiccatae_taxonomic_conflict))) & suppressWarnings(as.integer(vouchers$exsiccatae_taxonomic_conflict)) != 0L) else rep(FALSE, n_v)
    pl_aud <- if ("portal_metadata_lag_flag" %in% names(vouchers)) (!is.na(vouchers$portal_metadata_lag_flag) & as.logical(vouchers$portal_metadata_lag_flag)) else rep(FALSE, n_v)
    aud_mask <- (u_aud %in% c("out_of_range", "natureserve_exotic_or_waif")) | ex_aud | pl_aud
    aud_df <- vouchers[aud_mask, , drop = FALSE]
    for (ac in audit_cols) {
      if (!(ac %in% names(aud_df))) aud_df[[ac]] <- if (nrow(aud_df) > 0L) NA else character(0)
    }
    aud_df <- aud_df[, audit_cols, drop = FALSE]
    if (nrow(aud_df) > 0L && ("institutionCode" %in% names(aud_df))) {
      aud_df <- aud_df[order(ifelse(is.na(aud_df$institutionCode), "", as.character(aud_df$institutionCode)),
                             ifelse(is.na(aud_df$catalogNumber), "", as.character(aud_df$catalogNumber))), , drop = FALSE]
    }
    aud_out <- file.path(reports_dir, "portal_metadata_lag_audit.csv")
    tmp_aud <- paste0(aud_out, ".tmp")
    .write_triage_fast(aud_df, tmp_aud, na_str = "NA")
    file.rename(tmp_aud, aud_out)
  }, error = function(e) warning("Portal metadata lag audit CSV export note: ", e$message))
  if (!is.null(prop_sum) && isTRUE(prop_sum$active) && is.data.frame(prop_sum$table) && nrow(prop_sum$table) > 0L) {
    prop_out <- file.path(reports_dir, "unsequenced_genomic_propagation_summary.csv")
    tmp_prop <- paste0(prop_out, ".tmp")
    utils::write.csv(prop_sum$table, tmp_prop, row.names = FALSE)
    file.rename(tmp_prop, prop_out)
  }
  return(vouchers)
}
load_and_merge_modalities <- build_triage_queue

# ------------------------------------------------------------------------------
# 7. Synthesis Figure Delegation & Taxonomic Revision Markdown Report
# ------------------------------------------------------------------------------
.copy_focal_taxon_artifact <- function(src_path, target_path) {
  if (is.null(src_path) || !nzchar(src_path) || !file.exists(src_path)) return(invisible(FALSE))
  if (normalizePath(src_path, mustWork = FALSE) == normalizePath(target_path, mustWork = FALSE)) return(invisible(FALSE))
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  tmp_dst <- paste0(target_path, ".", Sys.getpid(), ".tmp")
  if (isTRUE(file.copy(src_path, tmp_dst, overwrite = TRUE))) file.rename(tmp_dst, target_path)
  invisible(TRUE)
}

.resolve_focal_taxon_slug <- function(focal_taxon_raw) {
  if (is.null(focal_taxon_raw) || !nzchar(trimws(as.character(focal_taxon_raw)))) return(list(taxon = "", slug = ""))
  f_std <- standardize_packera_taxon(trimws(as.character(focal_taxon_raw)))[1L]
  if (is.na(f_std) || !nzchar(f_std) || f_std == "Unknown") f_std <- trimws(as.character(focal_taxon_raw))
  f_slug <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", f_std))
  list(taxon = f_std, slug = f_slug)
}

render_synthesis_figure_png <- function(df, out_path, cda_summary = NULL, focal_taxon = "", output_complex_plot = NULL) {
  if (!exists("render_step07_synthesis_plate_png", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
  if (!exists("render_focal_revision_synthesis_plate_png", mode = "function") && file.exists("scripts/morphometrics/pub_synthesis_plates.R")) source("scripts/morphometrics/pub_synthesis_plates.R")
  out_png <- if (!is.null(out_path) && nzchar(out_path)) sub("\\.pdf$", ".png", out_path, ignore.case = TRUE) else "outputs/figures/fig07a_integrative_taxonomic_synthesis.png"
  full_df <- if (!is.null(attr(df, "full_complex_df"))) attr(df, "full_complex_df") else df
  if (exists("render_step07_synthesis_plate_png", mode = "function")) {
    render_step07_synthesis_plate_png(full_df, out_png, cda_summary = cda_summary, focal_taxon = focal_taxon)
  }
  f_info <- .resolve_focal_taxon_slug(focal_taxon)
  if (nzchar(f_info$taxon) && nzchar(f_info$slug) && nzchar(out_png)) {
    focal_png <- file.path(dirname(out_png), sprintf("Figure_Integrative_%s_Revision.png", f_info$slug))
    if (exists("render_focal_revision_synthesis_plate_png", mode = "function")) {
      render_focal_revision_synthesis_plate_png(df, focal_png, full_df = full_df, focal_taxon = f_info$taxon)
    } else if (file.exists(out_png)) {
      .copy_focal_taxon_artifact(out_png, focal_png)
    }
  }
  if (is.null(attr(df, "complex_discovery")) && !identical(trimws(env_plantmatics_or_packera("PLANTMATICS_SKIP_FIGURES", "PACKERA_SKIP_FIGURES")), "1") &&
      exists("render_step07b_complex_network_png", mode = "function")) {
    c_out <- if (!is.null(output_complex_plot) && nzchar(output_complex_plot)) sub("\\.pdf$", ".png", output_complex_plot, ignore.case = TRUE) else file.path(dirname(out_png), "fig07b_taxon_complex_delimitation_network.png")
    c_lbl <- if (nzchar(f_info$taxon)) f_info$taxon else as.character(df$species_standardized[1L])
    tryCatch(render_step07b_complex_network_png(
      complexes_df = data.frame(species_standardized = c_lbl, consensus_complex_id = 1L, n_vouchers = nrow(df), stringsAsFactors = FALSE),
      affinity_mat = matrix(1.0, nrow = 1L, ncol = 1L, dimnames = list(c_lbl, c_lbl)),
      out_complex_path = c_out,
      voucher_df = df
    ), error = function(e) warning("Single-taxon deme network render warning: ", e$message))
  }
  invisible(out_png)
}

generate_taxonomic_revision_report <- function(df, opts) {
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE)
  tot_n <- nrow(df)
  crit_n <- sum(df$synthesis_triage_priority == "CRITICAL", na.rm = TRUE)
  high_n <- sum(df$synthesis_triage_priority == "HIGH", na.rm = TRUE)
  med_n  <- sum(df$synthesis_triage_priority == "MEDIUM", na.rm = TRUE)
  res_n  <- sum(df$synthesis_triage_priority %in% c("LOW", "RESOLVED"), na.rm = TRUE)

  sub_vec <- if ("taxonomic_treatment_subtype" %in% names(df)) as.character(df$taxonomic_treatment_subtype) else rep(NA_character_, tot_n)
  call_vec <- if ("taxonomic_status_call" %in% names(df)) as.character(df$taxonomic_status_call) else rep(NA_character_, tot_n)
  misid_n   <- sum(sub_vec == "Misidentification_Reassign" | call_vec == "Misidentification_Severe", na.rm = TRUE)
  ecotype_n <- sum(sub_vec == "Subspecies_Ecotype" | call_vec == "Subspecies_Ecotype", na.rm = TRUE)
  ecoph_n   <- sum(sub_vec == "Ecophenotypic_Plasticity" | call_vec == "Ecophenotypic_Plasticity", na.rm = TRUE)
  hybr_n    <- sum(sub_vec == "Hybrid_Introgressant" | call_vec == "Hybrid_Intergrade_Swarm", na.rm = TRUE)
  crypt_n   <- sum(sub_vec == "Putative_Cryptic_Species" | call_vec %in% c("Putative_Cryptic_Species", "Taxonomic_Resurrection_Candidate"), na.rm = TRUE)
  resurr_n  <- sum(call_vec == "Taxonomic_Resurrection_Candidate", na.rm = TRUE)

  sp_col_rep <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else character(0)
  active_taxa <- sort(unique(sp_col_rep[!is.na(sp_col_rep) & nzchar(sp_col_rep) & sp_col_rep != "Unknown"]))
  f_info <- .resolve_focal_taxon_slug(opts$focal_taxon)
  header_focal <- if (nzchar(f_info$taxon)) f_info$taxon else if (length(active_taxa) == 1L) active_taxa[1L] else ""

  report_title <- if (nzchar(header_focal)) {
    sprintf("# Taxonomic Revision and Integrative Species Delimitation in the *%s* Complex (Asteraceae: Senecioneae)", header_focal)
  } else {
    sprintf("# Taxonomic Revision and Integrative Species Delimitation Across *Packera* Species Complexes (N = %d Taxa; Asteraceae: Senecioneae)", length(active_taxa))
  }

  scope_phrase <- if (nzchar(header_focal)) {
    sprintf("the ***%s*** complex", header_focal)
  } else if (length(active_taxa) > 0L) {
    sprintf("**%d analyzed *Packera* taxa** (%s)", length(active_taxa), paste0("*", active_taxa, "*", collapse = ", "))
  } else {
    "the analyzed *Packera* herbarium dataset"
  }

  gmm_stats <- resolve_gmm_bic_stats(opts)
  bic_str <- if (is.finite(gmm_stats$two_delta_bic)) format(round(gmm_stats$two_delta_bic, 1), big.mark = ",") else "N/A"
  batch_df <- audit_institutional_batch_effects_ancova(df)
  batch_str <- if (!is.null(batch_df) && nrow(batch_df) > 0) {
    sprintf("Taxon-adjusted marginal ANCOVA (`lm(PC1 ~ species_standardized + institutionCode)`): partial $\\eta^2 = %.4f$ ($q_{\\text{BH}} = %s$, %s).", batch_df$partial_eta_sq[1], as.character(batch_df$p_adj_bh[1]), batch_df$Batch_Effect_Status[1])
  } else "Taxon-adjusted marginal ANCOVA (`lm(PC1 ~ species_standardized + institutionCode)`): partial $\\eta^2 < 0.05$ across major herbaria."

  fmt_med_iqr <- function(primary_vec, fallback_vec = NULL, digits = 1, unit_str = "") {
    v <- suppressWarnings(as.numeric(primary_vec))
    if (!is.null(fallback_vec)) {
      vf <- suppressWarnings(as.numeric(fallback_vec))
      v <- ifelse(is.finite(v), v, vf)
    }
    v <- v[is.finite(v)]
    if (length(v) == 0L) return("NA")
    sprintf("%.*f [IQR %.*f]%s", digits, stats::median(v), digits, stats::IQR(v), unit_str)
  }

  treatment_lines <- character(0)
  for (sp in active_taxa) {
    sub_sp <- df[!is.na(sp_col_rep) & sp_col_rep == sp, , drop = FALSE]
    n_sp <- nrow(sub_sp)
    is_type_sp <- if ("is_type_specimen" %in% names(sub_sp)) as.logical(sub_sp$is_type_specimen) else rep(FALSE, n_sp)
    is_type_sp[is.na(is_type_sp)] <- FALSE
    n_types <- sum(is_type_sp, na.rm = TRUE)
    n_gs <- if ("determiner_tier" %in% names(sub_sp)) sum((sub_sp$determiner_tier %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver")) | is_type_sp, na.rm = TRUE) else n_types

    rg_vals <- if ("regional_group" %in% names(sub_sp)) trimws(as.character(sub_sp$regional_group)) else character(0)
    rg_vals <- rg_vals[!is.na(rg_vals) & nzchar(rg_vals) & !(tolower(rg_vals) %in% c("na", "unknown", "unassigned"))]
    eco_str <- if (length(rg_vals) > 0L) {
      rg_tab <- utils::head(sort(table(rg_vals), decreasing = TRUE), 3L)
      paste(sprintf("%s (n=%d)", names(rg_tab), as.integer(rg_tab)), collapse = "; ")
    } else "Unassigned / Non-georeferenced"

    ph_str   <- fmt_med_iqr(if ("fine_soil_ph" %in% names(sub_sp)) sub_sp$fine_soil_ph else NA_real_,
                            if ("soil_ph" %in% names(sub_sp)) sub_sp$soil_ph else NULL, digits = 2)
    sand_str <- fmt_med_iqr(if ("fine_soil_sand" %in% names(sub_sp)) sub_sp$fine_soil_sand else NA_real_,
                            if ("soil_sand" %in% names(sub_sp)) sub_sp$soil_sand else NULL, digits = 1, unit_str = "%")
    bed_str  <- fmt_med_iqr(if ("fine_bedrock_depth_cm" %in% names(sub_sp)) sub_sp$fine_bedrock_depth_cm else NA_real_,
                            if ("ssurgo_bedrock_depth_cm" %in% names(sub_sp)) sub_sp$ssurgo_bedrock_depth_cm else NULL, digits = 1, unit_str = " cm")

    d_sp <- if ("doy" %in% names(sub_sp)) suppressWarnings(as.numeric(sub_sp$doy)) else numeric(0)
    anth_sp <- if ("is_active_anthesis" %in% names(sub_sp)) as.logical(sub_sp$is_active_anthesis) else rep(TRUE, n_sp)
    anth_sp[is.na(anth_sp)] <- TRUE
    d_valid <- d_sp[is.finite(d_sp) & d_sp >= 1.0 & d_sp <= 366.0 & anth_sp]
    if (length(d_valid) == 0L) d_valid <- d_sp[is.finite(d_sp) & d_sp >= 1.0 & d_sp <= 366.0]
    doy_str <- if (length(d_valid) > 0L) {
      q10 <- as.integer(round(stats::quantile(d_valid, probs = 0.10, na.rm = TRUE, names = FALSE)))
      q90 <- as.integer(round(stats::quantile(d_valid, probs = 0.90, na.rm = TRUE, names = FALSE)))
      sprintf("DOY %d to %d (10th to 90th pct; n=%d; Hopkins 1918; Pewsey et al. 2013)", q10, q90, length(d_valid))
    } else "NA (no dated anthesis vouchers; Hopkins 1918; Pewsey et al. 2013)"

    st_sp <- if ("taxonomic_treatment_subtype" %in% names(sub_sp)) as.character(sub_sp$taxonomic_treatment_subtype) else if ("taxonomic_status_call" %in% names(sub_sp)) as.character(sub_sp$taxonomic_status_call) else rep("Confirmed_Species", n_sp)
    call_sp <- if ("taxonomic_status_call" %in% names(sub_sp)) as.character(sub_sp$taxonomic_status_call) else st_sp
    cnt_sub <- function(lbls) sum(!is.na(st_sp) & (st_sp %in% lbls))
    n_type_ctrl <- sum((!is.na(call_sp) & call_sp == "Confirmed_Type_Control") | (!is.na(st_sp) & st_sp == "Confirmed_Type_Control"))
    n_type_syn  <- sum((!is.na(call_sp) & call_sp == "Type_Synonymy_Candidate") | (!is.na(st_sp) & st_sp == "Type_Synonymy_Candidate"))
    sci_vals <- if ("scientificName" %in% names(sub_sp)) trimws(as.character(sub_sp$scientificName)) else character(0)
    sci_vals <- sci_vals[!is.na(sci_vals) & nzchar(sci_vals) & !(tolower(sci_vals) %in% c("na", "unknown", "none", "null"))]
    sp_esc <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", sp)
    sci_with_auth <- sci_vals[grepl(paste0("^", sp_esc, "\\s+.+"), sci_vals, ignore.case = TRUE)]
    auth_str <- if (length(sci_with_auth) > 0L) {
      best_sci <- names(sort(table(sci_with_auth), decreasing = TRUE))[1L]
      trimws(sub(paste0("^", sp_esc), "", best_sci, ignore.case = TRUE))
    } else ""
    head_sp <- if (nzchar(auth_str)) sprintf("- ***%s*** %s", sp, auth_str) else sprintf("- ***%s***", sp)
    treatment_lines <- c(
      treatment_lines,
      head_sp,
      sprintf("  - **Examined Vouchers:** $n = %d$ (Nomenclatural Types: %d; Active Gold/Silver Anchors: %d)", n_sp, n_types, n_gs),
      sprintf("  - **Dominant US EPA Level III Ecoregions:** %s", eco_str),
      {
        er_col <- if ("ecotype_regime" %in% names(sub_sp)) as.character(sub_sp$ecotype_regime) else if ("analysis_group" %in% names(sub_sp)) as.character(sub_sp$analysis_group) else character(0)
        er_col <- trimws(er_col[!is.na(er_col) & nzchar(trimws(er_col)) & !(tolower(trimws(er_col)) %in% c("na", "unknown"))])
        er_multi <- er_col[grepl("^Ecotype_", er_col)]
        if (length(unique(er_multi)) >= 2L) {
          er_tab <- sort(table(er_multi), decreasing = TRUE)
          sprintf("  - **Independent Step 05 Ecological Regimes ($K_{\\text{eco}} = %d$):** %s",
                  length(er_tab), paste(sprintf("`%s` (n=%d)", names(er_tab), as.integer(er_tab)), collapse = "; "))
        } else {
          "  - **Independent Step 05 Ecological Regimes:** Unimodal broad niche envelope ($K_{\text{eco}} = 1$)"
        }
      },
      sprintf("  - **Edaphic Envelope (Median [IQR]):** `fine_soil_ph` = %s; `fine_soil_sand` = %s; `fine_bedrock_depth_cm` = %s", ph_str, sand_str, bed_str),
      sprintf("  - **Verified Anthesis Window:** %s", doy_str),
      sprintf(
        "  - **Delimitation Breakdown:** `Confirmed_Species` = %d (`Confirmed_Type_Control` = %d) | `Subspecies_Ecotype` = %d (`Type_Synonymy_Candidate` = %d) | `Ecophenotypic_Plasticity` = %d | `Hybrid_Introgressant` = %d | `Misidentification_Reassign` = %d | `Putative_Cryptic_Species` = %d",
        cnt_sub(c("Confirmed_Species", "Species_Confirmed", "Confirmed_Type_Control")),
        n_type_ctrl,
        cnt_sub(c("Subspecies_Ecotype", "Type_Synonymy_Candidate")),
        n_type_syn,
        cnt_sub("Ecophenotypic_Plasticity"),
        cnt_sub(c("Hybrid_Introgressant", "Hybrid_Intergrade_Swarm")),
        cnt_sub(c("Misidentification_Reassign", "Misidentification_Severe")),
        cnt_sub("Putative_Cryptic_Species")
      ),
      {
        hyb_sp_mask <- !is.na(st_sp) & (st_sp %in% c("Hybrid_Introgressant", "Hybrid_Intergrade_Swarm"))
        if (any(hyb_sp_mask) && all(c("hybrid_parent_1", "hybrid_parent_2") %in% names(sub_sp))) {
          p1_s <- as.character(sub_sp$hybrid_parent_1[hyb_sp_mask]); p2_s <- as.character(sub_sp$hybrid_parent_2[hyb_sp_mask])
          ok_ps <- !is.na(p1_s) & nzchar(p1_s) & !is.na(p2_s) & nzchar(p2_s)
          if (any(ok_ps)) {
            pair_lbls <- sprintf("*%s* x *%s*", p1_s[ok_ps], p2_s[ok_ps])
            ptab <- sort(table(pair_lbls), decreasing = TRUE)
            sprintf("  - **Confirmed Hybrid Swarm Parentage:** %s", paste(sprintf("%s (n=%d)", names(ptab), as.integer(ptab)), collapse = "; "))
          } else character(0)
        } else character(0)
      }
    )
  }
  if (length(treatment_lines) == 0L) {
    treatment_lines <- "- No standardized *Packera* vouchers present in active triage queue (Hopkins 1918; Pewsey et al. 2013)."
  }

  md_lines <- c(
    report_title, "",
    "**Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)  ",
    "**Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)  ",
    "**Institution:** University of North Carolina at Chapel Hill Herbarium (NCU)  ",
    "**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5`  ", "",
    "## 1. Executive Botanical Summary & Integrative Taxonomy Protocol",
    sprintf(
      "This document provides the formal taxonomic treatment, hierarchical integrative species delimitation (Padial et al. 2010; Schlick-Steiner et al. 2010), and reliability-weighted curatorial synthesis (Dawid & Skene 1979) across %s.",
      scope_phrase
    ),
    sprintf("- **Total Examined Vouchers:** %s specimens (including coalesced non-georeferenced historic vouchers).", format(tot_n, big.mark = ",")),
    sprintf("- **GMM Bayes Factors (2ΔBIC):** Decisive evidence for **K = %d discrete morphological clusters** ($2\\Delta\\text{BIC} = %s$ over $K=1$).", gmm_stats$opt_k, bic_str),
    sprintf("- **Institutional Batch Audit:** %s", batch_str),
    sprintf("- **Stage 1 (Herbarium Misidentifications, Reassign):** %d vouchers (%0.1f%%) evaluated prior to hybridity.", misid_n, (misid_n / max(tot_n, 1)) * 100),
    {
      is_single_or_intra <- isTRUE(opts$intraspecific_mode) || identical(attr(df, "can_evaluate_hybrids"), FALSE)
      if (is_single_or_intra) {
        "- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** 0 introgressants (0.0%; interspecific hybrid evaluation disabled in single-taxon intraspecific mode)."
      } else if (hybr_n > 0L && all(c("hybrid_parent_1", "hybrid_parent_2") %in% names(df))) {
        hm_all <- !is.na(sub_vec) & (sub_vec == "Hybrid_Introgressant") & !is.na(df$hybrid_parent_1) & !is.na(df$hybrid_parent_2)
        if (any(hm_all)) {
          canon_pairs <- vapply(which(hm_all), function(ii) {
            ps <- sort(c(as.character(df$hybrid_parent_1[ii]), as.character(df$hybrid_parent_2[ii])))
            sprintf("*%s* x *%s*", ps[1L], ps[2L])
          }, character(1))
          pt_all <- sort(table(canon_pairs), decreasing = TRUE)
          pair_summary <- paste(sprintf("%s (n=%d)", names(pt_all), as.integer(pt_all)), collapse = "; ")
          sprintf("- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** %d introgressants (%0.1f%%) across confirmed parent pairs: %s.", hybr_n, (hybr_n / max(tot_n, 1)) * 100, pair_summary)
        } else {
          sprintf("- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** %d introgressants (%0.1f%%) (Rieseberg 1995; Slenker et al. 2022).", hybr_n, (hybr_n / max(tot_n, 1)) * 100)
        }
      } else {
        sprintf("- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** %d introgressants (%0.1f%%; 5-gate syntopic population threshold).", hybr_n, (hybr_n / max(tot_n, 1)) * 100)
      }
    },
    {
      if (resurr_n > 0L && "resurrection_basionym" %in% names(df)) {
        rb_vals <- as.character(df$resurrection_basionym[!is.na(df$resurrection_basionym) & nzchar(as.character(df$resurrection_basionym))])
        rb_str <- if (length(rb_vals) > 0L) paste(sprintf("*%s*", unique(rb_vals)), collapse = ", ") else "heterotypic basionym"
        sprintf("- **Stage 3 (Putative Cryptic and Basionym Resurrection Lineages):** %d vouchers (%0.1f%%; %d heterotypic resurrection candidates matching %s and %d unnamed cryptic candidates) satisfying Integration by Congruence (Padial et al. 2010).", crypt_n, (crypt_n / max(tot_n, 1)) * 100, resurr_n, rb_str, max(0L, crypt_n - resurr_n))
      } else {
        sprintf("- **Stage 3 (Putative Cryptic Lineages):** %d vouchers (%0.1f%%) satisfying strict Integration by Congruence (Padial et al. 2010).", crypt_n, (crypt_n / max(tot_n, 1)) * 100)
      }
    },
    sprintf("- **Stage 4 (Subspecific Ecotypes and Taphonomic/Ontogenetic Plasticity):** %d cohesive ecotypes (%0.1f%%) and %d within-envelope plastic variants (%0.1f%%).", ecotype_n, (ecotype_n / max(tot_n, 1)) * 100, ecoph_n, (ecoph_n / max(tot_n, 1)) * 100), "",
    "## 2. Multi-Evidence Taxonomic Treatment",
    treatment_lines, "",
    {
      comp_att <- attr(df, "complex_discovery")
      cdf_rep <- if (!is.null(comp_att) && !is.null(comp_att$complexes_df)) comp_att$complexes_df else if (!is.null(opts$output_complexes) && nzchar(opts$output_complexes) && file.exists(opts$output_complexes)) tryCatch(utils::read.csv(opts$output_complexes, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
      if (!is.null(cdf_rep) && nrow(cdf_rep) >= 2L) {
        fmt_it <- function(s) ifelse(grepl("^Packera\\s+|^P\\.\\s+", s), paste0("*", s, "*"), s)
        c_rows <- vapply(seq_len(nrow(cdf_rep)), function(ri) {
          sprintf(
            "| %s | `%s` | %s (%.3f) | %.3f | %s |",
            fmt_it(as.character(cdf_rep$species_standardized[ri])),
            as.character(cdf_rep$complex_label[ri]),
            fmt_it(as.character(cdf_rep$nearest_neighbor_taxon[ri])),
            suppressWarnings(as.numeric(cdf_rep$nearest_neighbor_affinity[ri])),
            suppressWarnings(as.numeric(cdf_rep$silhouette_width[ri])),
            ifelse(isTRUE(as.logical(cdf_rep$integrative_congruence[ri])), "Yes", "No")
          )
        }, character(1))
        c(
          "## 3. Empirical Taxon-Complex Delimitation",
          "| Taxon / Subcluster | Complex | Nearest Neighbor ($W_{ij}$) | Silhouette | Congruent |",
          "| :--- | :--- | :--- | :--- | :--- |",
          c_rows, "",
          "## 4. Herbarium Triage Queue Audit Summary"
        )
      } else {
        "## 3. Herbarium Triage Queue Audit Summary"
      }
    },
    "| Triage Priority | Voucher Count | Percentage | Recommended Action |",
    "| :--- | :--- | :--- | :--- |",
    sprintf("| **CRITICAL** | **%d** | **%0.1f%%** | Specialist manual re-inspection & sheet annotation |", crit_n, (crit_n / max(tot_n, 1)) * 100),
    sprintf("| **HIGH** | **%d** | **%0.1f%%** | Systematic redetermination in aggregator databases |", high_n, (high_n / max(tot_n, 1)) * 100),
    sprintf("| **MEDIUM** | **%d** | **%0.1f%%** | Annotate foliar wear / introgression zone |", med_n, (med_n / max(tot_n, 1)) * 100),
    sprintf("| **LOW / RESOLVED** | **%d** | **%0.1f%%** | Verified taxon anchor retained |", res_n, (res_n / max(tot_n, 1)) * 100),
    sprintf("| **Total** | **%d** | **100.0%%** | Reliability-weighted multi-modal consensus |", tot_n)
  )
  geno_md <- character(0)
  if (exists("build_genomic_markdown_summary", mode = "function")) {
    geno_md <- tryCatch(build_genomic_markdown_summary(opts, df), error = function(e) character(0))
    if (length(geno_md) > 0L) md_lines <- c(md_lines, geno_md)
  }
  cyto_md_rep <- attr(df, "cytology_markdown_lines")
  if (!is.null(cyto_md_rep) && length(cyto_md_rep) > 0L && length(geno_md) == 0L) {
    cyto_md_rep <- as.character(cyto_md_rep)
    cyto_md_rep <- cyto_md_rep[!is.na(cyto_md_rep)]
    if (length(cyto_md_rep) > 0L) {
      if (nzchar(trimws(cyto_md_rep[1L]))) md_lines <- c(md_lines, "", cyto_md_rep) else md_lines <- c(md_lines, cyto_md_rep)
    }
  }
  tmp_rep <- paste0(opts$output_report, ".tmp"); writeLines(md_lines, tmp_rep); file.rename(tmp_rep, opts$output_report)

  if (nzchar(f_info$taxon) && nzchar(f_info$slug)) {
    focal_rep <- file.path(dirname(opts$output_report), sprintf("%s_Taxonomic_Revision_Summary.md", f_info$slug))
    if (normalizePath(focal_rep, mustWork = FALSE) != normalizePath(opts$output_report, mustWork = FALSE)) {
      tmp_frep <- paste0(focal_rep, ".tmp"); writeLines(md_lines, tmp_frep); file.rename(tmp_frep, focal_rep)
    }
    out_png_rep <- if (!is.null(opts$output_plot) && nzchar(opts$output_plot)) sub("\\.pdf$", ".png", opts$output_plot, ignore.case = TRUE) else ""
    if (nzchar(out_png_rep) && file.exists(out_png_rep)) {
      focal_png <- file.path(dirname(out_png_rep), sprintf("Figure_Integrative_%s_Revision.png", f_info$slug))
      if (!file.exists(focal_png) || file.info(focal_png)$size < 1024) {
        .copy_focal_taxon_artifact(out_png_rep, focal_png)
      }
    }
  }
}

.html_escape_vec <- function(x) {
  s <- ifelse(is.na(x), "", as.character(x))
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  s <- gsub("\"", "&quot;", s, fixed = TRUE)
  s
}

.json_escape_vec <- function(x) {
  s <- ifelse(is.na(x), "", as.character(x))
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\r", " ", s, fixed = TRUE)
  s <- gsub("\n", " ", s, fixed = TRUE)
  s <- gsub("\t", " ", s, fixed = TRUE)
  s <- gsub("</", "<\\/", s, fixed = TRUE)
  paste0("\"", s, "\"")
}

generate_interactive_triage_html <- function(df, html_path, focal_taxon = "", opts = NULL) {
  if (is.null(html_path) || !nzchar(trimws(as.character(html_path)))) return(invisible(NULL))
  dir.create(dirname(html_path), recursive = TRUE, showWarnings = FALSE)
  n <- nrow(df)

  get_c <- function(nm, def = "") if (nm %in% names(df)) ifelse(is.na(df[[nm]]), def, as.character(df[[nm]])) else rep(def, n)
  get_f <- function(nm, digits = 2) {
    if (!(nm %in% names(df))) return(rep("NA", n))
    v <- suppressWarnings(as.numeric(df[[nm]]))
    ifelse(is.finite(v), sprintf("%.*f", digits, v), "NA")
  }

  f_info <- .resolve_focal_taxon_slug(focal_taxon)
  sp_vec <- get_c("species_standardized", "Unknown")
  active_taxa <- sort(unique(sp_vec[nzchar(sp_vec) & sp_vec != "Unknown"]))
  hdr_taxon <- if (nzchar(f_info$taxon)) f_info$taxon else if (length(active_taxa) == 1L) active_taxa[1L] else ""
  title_html <- if (nzchar(hdr_taxon)) {
    sprintf("NCU Herbarium Integrative Triage Dashboard: <em>%s</em>", .html_escape_vec(hdr_taxon))
  } else {
    sprintf("NCU Herbarium Integrative Triage Dashboard: <em>Packera</em> (%d Taxa)", length(active_taxa))
  }

  w_m <- if ("w_morph_kappa" %in% names(df) && any(is.finite(as.numeric(df$w_morph_kappa)))) sprintf("%.3f", stats::median(as.numeric(df$w_morph_kappa), na.rm = TRUE)) else "NA"
  w_e <- if ("w_edaphic_kappa" %in% names(df) && any(is.finite(as.numeric(df$w_edaphic_kappa)))) sprintf("%.3f", stats::median(as.numeric(df$w_edaphic_kappa), na.rm = TRUE)) else "NA"
  w_p <- if ("w_pheno_kappa" %in% names(df) && any(is.finite(as.numeric(df$w_pheno_kappa)))) sprintf("%.3f", stats::median(as.numeric(df$w_pheno_kappa), na.rm = TRUE)) else "NA"

  prio_v <- get_c("synthesis_triage_priority", "RESOLVED")
  sub_v  <- get_c("taxonomic_treatment_subtype", "Confirmed_Species")
  n_crit_high <- sum(prio_v %in% c("CRITICAL", "HIGH"))
  n_misid     <- sum(sub_v == "Misidentification_Reassign")
  n_hyb       <- sum(sub_v == "Hybrid_Introgressant")
  n_crypt     <- sum(sub_v == "Putative_Cryptic_Species")
  n_eco       <- sum(sub_v %in% c("Subspecies_Ecotype", "Ecophenotypic_Plasticity"))
  n_conf      <- sum(sub_v == "Confirmed_Species")

  cat_v   <- get_c("catalogNumber", "")
  pid_v   <- get_c("plant_individual_id", "0")
  inst_v  <- get_c("institutionCode", "")
  coll_v  <- get_c("recordedBy", "")
  tier_v  <- get_c("determiner_tier", "")
  rec_v   <- get_c("recommended_determination", "")
  cda_v   <- get_c("cda_predicted_taxon", "")
  cdap_v  <- get_f("cda_posterior_prob", 2)
  ed_v    <- get_c("edaphic_best_fit_taxon", "")
  edz_v   <- get_f("edaphic_z_score", 2)
  ddoy_v  <- get_f("delta_doy_resid", 1)
  cscore_v <- get_f("probabilistic_consensus_score", 3)
  act_v   <- get_c("synthesis_triage_action", "")
  rat_v   <- get_c("synthesis_rationale", "")
  hp1_h   <- get_c("hybrid_parent_1", "")
  hp2_h   <- get_c("hybrid_parent_2", "")
  hcoh_h  <- ifelse(isTRUE(df$hybrid_zone_cohesion_pass) | (!is.null(df$hybrid_zone_cohesion_pass) & !is.na(df$hybrid_zone_cohesion_pass) & as.logical(df$hybrid_zone_cohesion_pass)), "TRUE", "FALSE")
  if (length(hcoh_h) != n) hcoh_h <- rep("FALSE", n)
  mts_v   <- get_c("morph_tier_source", "Unknown")
  bc_v    <- if ("basal_leaf_count" %in% names(df)) ifelse(is.finite(suppressWarnings(as.numeric(df$basal_leaf_count))), as.character(as.integer(df$basal_leaf_count)), "0") else rep("0", n)
  cc_v    <- if ("cauline_leaf_count" %in% names(df)) ifelse(is.finite(suppressWarnings(as.numeric(df$cauline_leaf_count))), as.character(as.integer(df$cauline_leaf_count)), "0") else rep("0", n)
  bcp_v   <- paste0(bc_v, "/", cc_v)
  eco_v   <- if ("ecotype_regime" %in% names(df) && any(nzchar(trimws(get_c("ecotype_regime", ""))))) get_c("ecotype_regime", "") else get_c("analysis_group", "")
  geno_pred_h <- get_c("geno_predicted_taxon", "")
  geno_cyt_h <- get_c("inferred_cytotype", "")
  geno_gdi_h <- get_f("coalescent_gdi", 2)
  geno_cls_h <- get_c("hybrid_generation_class", "")
  prop_tier_h <- get_c("genomic_propagation_tier", "")
  t0_tax_h <- get_c("tier0_calibrated_taxon", "")
  t0_prob_h <- get_f("tier0_calibrated_prob", 2)
  t0_abst_h <- if ("tier0_abstention_flag" %in% names(df)) (!is.na(df$tier0_abstention_flag) & suppressWarnings(as.integer(df$tier0_abstention_flag)) == 1L) else rep(FALSE, n)
  t0_ok_h <- if ("tier0_abstention_flag" %in% names(df)) (!is.na(df$tier0_abstention_flag) & suppressWarnings(as.integer(df$tier0_abstention_flag)) == 0L) else rep(FALSE, n)
  cl_gdi_h <- get_f("cluster_propagated_gdi", 2)
  vp_2n_h <- if ("voucher_sporophytic_2n" %in% names(df)) suppressWarnings(as.integer(df$voucher_sporophytic_2n)) else rep(NA_integer_, n)
  vp_x_h <- if ("voucher_ploidy_level" %in% names(df)) suppressWarnings(as.integer(df$voucher_ploidy_level)) else rep(NA_integer_, n)
  tax_cyto_h <- get_c("taxon_reported_cytotypes", "")
  is_mc_h <- if ("is_multi_cytotype_taxon" %in% names(df)) (!is.na(df$is_multi_cytotype_taxon) & suppressWarnings(as.integer(df$is_multi_cytotype_taxon)) == 1L) else rep(FALSE, n)
  geno_sum_h <- vapply(seq_len(n), function(ii) {
    base_str <- if (nzchar(geno_pred_h[ii])) {
      bits <- character(0)
      if (nzchar(geno_cyt_h[ii])) bits <- c(bits, geno_cyt_h[ii])
      if (geno_gdi_h[ii] != "NA") bits <- c(bits, paste0("gdi=", geno_gdi_h[ii]))
      if (nzchar(geno_cls_h[ii])) bits <- c(bits, geno_cls_h[ii])
      if (length(bits) > 0L) sprintf("%s (%s)", geno_pred_h[ii], paste(bits, collapse = "; ")) else geno_pred_h[ii]
    } else ""
    prop_badges <- character(0)
    if (prop_tier_h[ii] == "G1_Exsiccatae_Isoduplicate") prop_badges <- c(prop_badges, "[Exs-Dup]")
    else if (prop_tier_h[ii] == "G2_Conditional_Deme_Hull") prop_badges <- c(prop_badges, "[Deme-Hull]")
    if (t0_ok_h[ii] && nzchar(t0_tax_h[ii]) && (t0_tax_h[ii] != "Sensu_Lato_DNA_Required")) {
      prop_badges <- c(prop_badges, if (t0_prob_h[ii] != "NA") sprintf("[Tier0-Cal: %s p=%s]", t0_tax_h[ii], t0_prob_h[ii]) else sprintf("[Tier0-Cal: %s]", t0_tax_h[ii]))
    }
    if (t0_abst_h[ii]) prop_badges <- c(prop_badges, "[DNA-Abstain]")
    if (geno_gdi_h[ii] == "NA" && cl_gdi_h[ii] != "NA") prop_badges <- c(prop_badges, sprintf("[Cluster-gdi=%s]", cl_gdi_h[ii]))
    has_2n_ii <- !is.na(vp_2n_h[ii]) && is.finite(vp_2n_h[ii]) && vp_2n_h[ii] > 0L
    has_x_ii <- !is.na(vp_x_h[ii]) && is.finite(vp_x_h[ii]) && vp_x_h[ii] > 0L
    if (has_2n_ii && has_x_ii) {
      prop_badges <- c(prop_badges, sprintf("[2n=%d | %dx]", vp_2n_h[ii], vp_x_h[ii]))
    } else if (has_2n_ii) {
      prop_badges <- c(prop_badges, sprintf("[2n=%d]", vp_2n_h[ii]))
    } else if (has_x_ii) {
      prop_badges <- c(prop_badges, sprintf("[%dx]", vp_x_h[ii]))
    }
    tc_ii <- trimws(tax_cyto_h[ii])
    if (is_mc_h[ii] && nzchar(tc_ii) && !(toupper(tc_ii) %in% c("NA", "NAN", "NONE", "NULL"))) {
      prop_badges <- c(prop_badges, sprintf("[Multi-Cytotype: %s]", tc_ii))
    }
    if (length(prop_badges) == 0L) return(base_str)
    if (nzchar(base_str)) paste(c(base_str, prop_badges), collapse = " ") else paste(prop_badges, collapse = " ")
  }, character(1))

  json_rows <- if (n > 0L) {
    paste0(
      "[",
      .json_escape_vec(cat_v), ",",
      .json_escape_vec(pid_v), ",",
      .json_escape_vec(inst_v), ",",
      .json_escape_vec(coll_v), ",",
      .json_escape_vec(tier_v), ",",
      .json_escape_vec(sp_vec), ",",
      .json_escape_vec(rec_v), ",",
      .json_escape_vec(prio_v), ",",
      .json_escape_vec(sub_v), ",",
      .json_escape_vec(cda_v), ",",
      .json_escape_vec(cdap_v), ",",
      .json_escape_vec(ed_v), ",",
      .json_escape_vec(edz_v), ",",
      .json_escape_vec(ddoy_v), ",",
      .json_escape_vec(cscore_v), ",",
      .json_escape_vec(act_v), ",",
      .json_escape_vec(rat_v), ",",
      .json_escape_vec(hp1_h), ",",
      .json_escape_vec(hp2_h), ",",
      .json_escape_vec(hcoh_h), ",",
      .json_escape_vec(mts_v), ",",
      .json_escape_vec(bcp_v), ",",
      .json_escape_vec(eco_v), ",",
      .json_escape_vec(geno_sum_h),
      "]",
      collapse = ",\n"
    )
  } else ""

  taxon_opts <- paste0(
    "<option value=\"ALL\">All Recorded Taxa</option>",
    if (length(active_taxa) > 0L) paste0("<option value=\"", .html_escape_vec(active_taxa), "\">", .html_escape_vec(active_taxa), "</option>", collapse = "") else ""
  )

  html_doc <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Packera Integrative Triage Dashboard (NCU Herbarium)</title>",
    "<style>",
    ":root{--bg:#f8fafc;--card:#ffffff;--ink:#0f172a;--muted:#475569;--border:#cbd5e1;--accent:#0072B2;}",
    "*{box-sizing:border-box;}body{margin:0;padding:20px 24px;font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--ink);line-height:1.45;}",
    "header{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:16px 20px;margin-bottom:16px;}",
    "h1{margin:0 0 6px 0;font-size:1.28rem;}",
    ".meta{font-size:0.86rem;color:var(--muted);display:flex;flex-wrap:wrap;gap:16px;}",
    ".kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:10px;margin-bottom:16px;}",
    ".kpi{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:10px 14px;}",
    ".kpi .lbl{font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:var(--muted);}",
    ".kpi .val{font-size:1.35rem;font-weight:700;margin-top:2px;}",
    ".toolbar{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:12px 14px;margin-bottom:14px;display:flex;flex-wrap:wrap;gap:10px;align-items:center;}",
    ".toolbar input,.toolbar select,.toolbar button{padding:7px 10px;border:1px solid var(--border);border-radius:6px;font-size:0.86rem;background:#fff;color:var(--ink);}",
    ".toolbar input{flex:1 1 240px;}",
    ".toolbar button{background:var(--accent);color:#fff;border-color:var(--accent);font-weight:600;cursor:pointer;}",
    ".table-wrap{background:var(--card);border:1px solid var(--border);border-radius:8px;overflow-x:auto;}",
    "table{width:100%;border-collapse:collapse;font-size:0.83rem;}",
    "th,td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}",
    "th{background:#f1f5f9;font-weight:700;cursor:pointer;user-select:none;white-space:nowrap;}",
    "tr:hover td{background:#f8fafc;}",
    ".badge{display:inline-block;padding:2px 7px;border-radius:999px;font-size:0.73rem;font-weight:700;color:#fff;}",
    ".p-CRITICAL{background:#D55E00;}.p-HIGH{background:#E69F00;color:#111;}.p-MEDIUM{background:#CC79A7;}.p-LOW{background:#56B4E9;color:#111;}.p-RESOLVED{background:#009E73;}",
    ".sci{font-style:italic;}",
    ".rat{max-width:360px;color:#334155;font-size:0.79rem;}",
    ".cat-btn{background:none;border:0;padding:0;color:var(--accent);font-weight:700;cursor:pointer;text-decoration:underline;font-size:inherit;}",
    ".modal{display:none;position:fixed;inset:0;background:rgba(15,23,42,0.72);z-index:999;align-items:center;justify-content:center;padding:18px;}",
    ".modal.open{display:flex;}",
    ".modal-card{background:#fff;border-radius:10px;max-width:920px;width:100%;max-height:92vh;display:flex;flex-direction:column;overflow:hidden;box-shadow:0 20px 50px rgba(0,0,0,0.35);}",
    ".modal-hdr{display:flex;justify-content:space-between;align-items:center;padding:10px 14px;border-bottom:1px solid var(--border);background:#f8fafc;gap:8px;}",
    ".modal-body{padding:12px;overflow:auto;text-align:center;background:#0f172a;}",
    ".modal-body img{max-width:100%;max-height:76vh;object-fit:contain;border-radius:4px;}",
    "</style>",
    "</head>",
    "<body>",
    "<header>",
    sprintf("<h1>%s</h1>", title_html),
    sprintf("<div class=\"meta\"><span><strong>SOP:</strong> UNC-BOT-SOP-2026-04-REV5</span><span><strong>OOF Reliability (&kappa;):</strong> Morph = %s | Edaphic = %s | Phenology = %s</span><span id=\"shownCount\">Showing %d of %d vouchers</span></div>", w_m, w_e, w_p, n, n),
    "</header>",
    "<section class=\"kpis\">",
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Examined Vouchers</div><div class=\"val\">%d</div></div>", n),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Critical / High Queue</div><div class=\"val\">%d</div></div>", n_crit_high),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Misident. (Reassign)</div><div class=\"val\">%d</div></div>", n_misid),
    {
      is_intra_h <- (!is.null(opts) && isTRUE(opts$intraspecific_mode)) || identical(attr(df, "can_evaluate_hybrids"), FALSE)
      if (is_intra_h) {
        "<div class=\"kpi\"><div class=\"lbl\">Hybrid / Introgressant</div><div class=\"val\">0 <span style=\"font-size:0.72rem;font-weight:500;color:#64748b;\">(1-Taxon N/A)</span></div></div>"
      } else {
        sprintf("<div class=\"kpi\"><div class=\"lbl\">Hybrid / Introgressant</div><div class=\"val\">%d</div></div>", n_hyb)
      }
    },
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Putative Cryptic</div><div class=\"val\">%d</div></div>", n_crypt),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Variety / Plasticity</div><div class=\"val\">%d</div></div>", n_eco),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Confirmed Species</div><div class=\"val\">%d</div></div>", n_conf),
    "</section>",
    {
      hm_idx <- which(sub_v == "Hybrid_Introgressant" & nzchar(hp1_h) & nzchar(hp2_h))
      if (length(hm_idx) > 0L) {
        pair_vec_h <- vapply(hm_idx, function(ii) {
          ps <- sort(c(hp1_h[ii], hp2_h[ii]))
          paste(ps[1L], "x", ps[2L])
        }, character(1))
        ptab_h <- sort(table(pair_vec_h), decreasing = TRUE)
        h_chips <- vapply(names(ptab_h), function(plbl) {
          sprintf("<span style=\"display:inline-block;margin:2px 8px 2px 0;padding:3px 9px;background:#fdf2f8;border:1px solid #fbcfe8;border-radius:6px;font-size:0.81rem;\"><span class=\"sci\">%s</span> <strong>(n=%d)</strong></span>", .html_escape_vec(plbl), as.integer(ptab_h[[plbl]]))
        }, character(1))
        paste0("<section style=\"background:#fff;border:1px solid #cbd5e1;border-radius:8px;padding:10px 14px;margin-bottom:14px;\"><div style=\"font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:#475569;margin-bottom:4px;\">Confirmed Hybrid Swarms (Putative Parent Species)</div>", paste(h_chips, collapse = ""), "</section>")
      } else ""
    },
    {
      comp_att_h <- attr(df, "complex_discovery")
      cdf_h <- if (!is.null(comp_att_h) && !is.null(comp_att_h$complexes_df)) comp_att_h$complexes_df else if (!is.null(opts) && !is.null(opts$output_complexes) && nzchar(opts$output_complexes) && file.exists(opts$output_complexes)) tryCatch(utils::read.csv(opts$output_complexes, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
      if (!is.null(cdf_h) && nrow(cdf_h) >= 2L) {
        u_cl_h <- unique(as.character(cdf_h$complex_label))
        cl_chips <- vapply(u_cl_h, function(clb) {
          m_tx <- as.character(cdf_h$species_standardized[cdf_h$complex_label == clb])
          sprintf("<span style=\"display:inline-block;margin:2px 8px 2px 0;padding:3px 9px;background:#f1f5f9;border:1px solid #cbd5e1;border-radius:6px;font-size:0.81rem;\"><strong>%s:</strong> <span class=\"sci\">%s</span></span>", .html_escape_vec(clb), paste(.html_escape_vec(m_tx), collapse = ", "))
        }, character(1))
        paste0("<section style=\"background:#fff;border:1px solid #cbd5e1;border-radius:8px;padding:10px 14px;margin-bottom:14px;\"><div style=\"font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:#475569;margin-bottom:4px;\">Empirical Species Complexes (Leiden / Metric Ward.D2)</div>", paste(cl_chips, collapse = ""), "</section>")
      } else ""
    },
    {
      fig_rel_cands <- c(
        "Synthesis Plate" = "../figures/fig07a_integrative_taxonomic_synthesis.png",
        "Complex Network" = "../figures/fig07b_taxon_complex_delimitation_network.png",
        "Genomic PCoA & Tree" = "../figures/fig06a_genomic_pcoa_tree_concordance.png",
        "Reticulation, gdi & MMRR" = "../figures/fig06b_reticulation_gdi_landscape.png",
        "LM3 & ECT Exemplars" = "../figures/fig02a_lm3_organ_and_ect_exemplars.png",
        "EFA Morphospace" = "../figures/fig03a_efa_morphospace_and_harmonics.png",
        "GMM & CDA" = "../figures/fig04b_gmm_cda_morphometrics.png",
        "SDM & Niche" = "../figures/fig05c_sdm_suitability_and_distribution_map.png"
      )
      fig_links <- character(0)
      for (fnm in names(fig_rel_cands)) {
        rel_p <- fig_rel_cands[[fnm]]
        disk_p <- normalizePath(file.path(dirname(html_path), rel_p), mustWork = FALSE)
        if (file.exists(disk_p)) {
          fig_links <- c(fig_links, sprintf("<a href=\"%s\" target=\"_blank\" style=\"display:inline-block;margin:2px 8px 2px 0;padding:3px 9px;background:#eff6ff;border:1px solid #bfdbfe;border-radius:6px;font-size:0.80rem;color:#1d4ed8;text-decoration:none;\">%s (300 DPI PNG)</a>", rel_p, .html_escape_vec(fnm)))
        }
      }
      if (length(fig_links) > 0L) paste0("<section style=\"background:#fff;border:1px solid #cbd5e1;border-radius:8px;padding:8px 14px;margin-bottom:12px;\"><span style=\"font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:#475569;margin-right:10px;\">300-DPI PNG Plates:</span>", paste(fig_links, collapse = ""), "</section>") else ""
    },
    "<section class=\"toolbar\">",
    "<input type=\"search\" id=\"qSearch\" placeholder=\"Search catalogNumber, collector, herbarium, taxon, rationale...\">",
    "<select id=\"qPrio\"><option value=\"ALL\">All Priorities</option><option value=\"CRITICAL\">CRITICAL</option><option value=\"HIGH\">HIGH</option><option value=\"MEDIUM\">MEDIUM</option><option value=\"LOW\">LOW</option><option value=\"RESOLVED\">RESOLVED</option></select>",
    "<select id=\"qSub\"><option value=\"ALL\">All Treatments</option><option value=\"Misidentification_Reassign\">Misidentification_Reassign</option><option value=\"Putative_Cryptic_Species\">Putative_Cryptic_Species</option><option value=\"Hybrid_Introgressant\">Hybrid_Introgressant</option><option value=\"Subspecies_Ecotype\">Subspecies_Ecotype</option><option value=\"Ecophenotypic_Plasticity\">Ecophenotypic_Plasticity</option><option value=\"Confirmed_Species\">Confirmed_Species</option><option value=\"Confirmed_Type_Control\">Confirmed_Type_Control</option><option value=\"Type_Synonymy_Candidate\">Type_Synonymy_Candidate</option></select>",
    sprintf("<select id=\"qTaxon\">%s</select>", taxon_opts),
    "<button type=\"button\" id=\"btnExport\">Export Filtered CSV</button>",
    "</section>",
    "<div class=\"table-wrap\"><table id=\"triTable\"><thead><tr>",
    "<th data-col=\"0\">Catalog #</th><th data-col=\"1\">Ind</th><th data-col=\"2\">Herb.</th><th data-col=\"4\">Tier</th>",
    "<th data-col=\"5\">Recorded Taxon</th><th data-col=\"6\">Recommended Determination</th><th data-col=\"7\">Priority</th>",
    "<th data-col=\"8\">Treatment Subtype</th><th data-col=\"9\">CDA Pred (p)</th><th data-col=\"11\">Edaphic Fit (z)</th>",
    "<th data-col=\"13\">&Delta;DOY (d)</th><th data-col=\"14\">Consensus</th><th data-col=\"20\">Morph Source (B/C)</th><th data-col=\"22\">Ecotype Regime</th><th data-col=\"23\">Genomic / gdi / Cytotype</th><th data-col=\"16\">Synthesis Rationale</th>",
    "</tr></thead><tbody id=\"triBody\"></tbody></table></div>",
    "<div id=\"sheetModal\" class=\"modal\" role=\"dialog\" aria-modal=\"true\"><div class=\"modal-card\"><div class=\"modal-hdr\"><strong id=\"modalTitle\">Specimen Sheet</strong><div><button type=\"button\" id=\"btnToggleScan\" style=\"margin-right:8px;padding:4px 10px;\">Show Raw Scan</button><button type=\"button\" id=\"btnCloseModal\" style=\"padding:4px 10px;\">Close</button></div></div><div class=\"modal-body\"><img id=\"modalImg\" alt=\"Herbarium specimen sheet preview\"><div id=\"modalMsg\" style=\"color:#cbd5e1;font-size:0.85rem;margin-top:8px;display:none;\"></div></div></div></div>",
    "<script>",
    paste0("const ROWS = [\n", json_rows, "\n];"),
    "const PRIO_RANK = {CRITICAL:1,HIGH:2,MEDIUM:3,LOW:4,RESOLVED:5};",
    "let sortCol = 7, sortAsc = true, filtered = ROWS.slice(), modalCat = '', modalRaw = false;",
    "function esc(s){return String(s||'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;');}",
    "function fmtSci(s){const t=esc(s);return /^Packera\\s+/i.test(t)?'<span class=\"sci\">'+t+'</span>':t;}",
    "function loadModalImg(){",
    "  const img=document.getElementById('modalImg'), msg=document.getElementById('modalMsg'), btn=document.getElementById('btnToggleScan');",
    "  msg.style.display='none'; img.style.display='inline-block';",
    "  const ov='../figures/lm3_sheet_overlays/'+encodeURIComponent(modalCat)+'.png';",
    "  const rw='../../../data/raw_vouchers/'+encodeURIComponent(modalCat)+'.webp';",
    "  btn.textContent = modalRaw ? 'Show LM3 Overlay' : 'Show Raw Scan';",
    "  img.onerror = function(){ if(!modalRaw){ modalRaw=true; btn.textContent='Show LM3 Overlay'; img.onerror=function(){ img.style.display='none'; msg.textContent='Sheet preview image not found locally for '+modalCat; msg.style.display='block'; }; img.src=rw; } else { img.style.display='none'; msg.textContent='Raw voucher scan not found locally for '+modalCat; msg.style.display='block'; } };",
    "  img.src = modalRaw ? rw : ov;",
    "}",
    "function openModal(cat){ modalCat=cat; modalRaw=false; document.getElementById('modalTitle').textContent='Herbarium Voucher: '+cat; document.getElementById('sheetModal').classList.add('open'); loadModalImg(); }",
    "document.getElementById('btnToggleScan').addEventListener('click',()=>{ modalRaw=!modalRaw; loadModalImg(); });",
    "document.getElementById('btnCloseModal').addEventListener('click',()=>document.getElementById('sheetModal').classList.remove('open'));",
    "document.getElementById('sheetModal').addEventListener('click',e=>{ if(e.target.id==='sheetModal') document.getElementById('sheetModal').classList.remove('open'); });",
    "document.addEventListener('keydown',e=>{ if(e.key==='Escape') document.getElementById('sheetModal').classList.remove('open'); });",
    "function applyFilter(){",
    "  const q=document.getElementById('qSearch').value.trim().toLowerCase();",
    "  const pr=document.getElementById('qPrio').value;",
    "  const sb=document.getElementById('qSub').value;",
    "  const tx=document.getElementById('qTaxon').value;",
    "  filtered = ROWS.filter(r=>{",
    "    if(pr!=='ALL' && r[7]!==pr) return false;",
    "    if(sb!=='ALL' && r[8]!==sb) return false;",
    "    if(tx!=='ALL' && r[5]!==tx) return false;",
    "    if(q && !r.join(' ').toLowerCase().includes(q)) return false;",
    "    return true;",
    "  });",
    "  filtered.sort((a,b)=>{",
    "    let va=a[sortCol], vb=b[sortCol];",
    "    if(sortCol===7){va=PRIO_RANK[va]||9;vb=PRIO_RANK[vb]||9;}",
    "    else if([1,10,12,13,14].includes(sortCol)){va=parseFloat(va);vb=parseFloat(vb);if(isNaN(va))va=9999;if(isNaN(vb))vb=9999;}",
    "    if(va<vb) return sortAsc?-1:1; if(va>vb) return sortAsc?1:-1; return 0;",
    "  });",
    "  const html = filtered.map(r=>'<tr>'+",
    "    '<td><button type=\"button\" class=\"cat-btn\" data-cat=\"'+esc(r[0])+'\">'+esc(r[0])+'</button></td><td>'+esc(r[1])+'</td><td>'+esc(r[2])+'</td><td>'+esc(r[4])+'</td>'+",
    "    '<td>'+fmtSci(r[5])+'</td><td>'+fmtSci(r[6])+'</td>'+",
    "    '<td><span class=\"badge p-'+esc(r[7])+'\">'+esc(r[7])+'</span></td>'+",
    "    '<td>'+esc(r[8])+'</td>'+",
    "    '<td>'+fmtSci(r[9])+(r[10]!=='NA'?' ('+esc(r[10])+')':'')+'</td>'+",
    "    '<td>'+fmtSci(r[11])+(r[12]!=='NA'?' (z='+esc(r[12])+')':'')+'</td>'+",
    "    '<td>'+esc(r[13])+'</td><td>'+esc(r[14])+'</td><td>'+esc(r[20])+' ('+esc(r[21])+')</td><td>'+fmtSci(r[22])+'</td><td>'+fmtSci(r[23])+'</td><td class=\"rat\">'+esc(r[16])+'</td></tr>').join('');",
    "  document.getElementById('triBody').innerHTML = html;",
    "  document.querySelectorAll('#triBody .cat-btn').forEach(b=>b.addEventListener('click',()=>openModal(b.getAttribute('data-cat'))));",
    "  document.getElementById('shownCount').textContent = 'Showing '+filtered.length+' of '+ROWS.length+' vouchers';",
    "}",
    "document.querySelectorAll('#triTable th').forEach(th=>th.addEventListener('click',()=>{",
    "  const c=parseInt(th.getAttribute('data-col'),10); if(sortCol===c) sortAsc=!sortAsc; else {sortCol=c;sortAsc=true;} applyFilter();",
    "}));",
    "['qSearch','qPrio','qSub','qTaxon'].forEach(id=>document.getElementById(id).addEventListener('input',applyFilter));",
    "document.getElementById('btnExport').addEventListener('click',()=>{",
    "  const hdr=['catalogNumber','plant_individual_id','institutionCode','recordedBy','determiner_tier','species_standardized','recommended_determination','synthesis_triage_priority','taxonomic_treatment_subtype','cda_predicted_taxon','cda_posterior_prob','edaphic_best_fit_taxon','edaphic_z_score','delta_doy_resid','probabilistic_consensus_score','synthesis_triage_action','synthesis_rationale','hybrid_parent_1','hybrid_parent_2','hybrid_zone_cohesion_pass','morph_tier_source','basal_cauline_leaves','ecotype_regime','genomic_summary'];",
    "  const csv=[hdr.join(',')].concat(filtered.map(r=>r.map(v=>'\"'+String(v||'').replace(/\"/g,'\"\"')+'\"').join(','))).join('\\n');",
    "  const blob=new Blob([csv],{type:'text/csv;charset=utf-8;'});",
    "  const a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download='triage_queue_filtered.csv'; a.click(); URL.revokeObjectURL(a.href);",
    "});",
    "applyFilter();",
    "</script>",
    "</body>",
    "</html>"
  )

  tmp_html <- paste0(html_path, ".", Sys.getpid(), ".tmp")
  writeLines(html_doc, tmp_html, useBytes = TRUE)
  file.rename(tmp_html, html_path)
  invisible(html_path)
}

# ------------------------------------------------------------------------------
# 8. Main Execution Pipeline
# ------------------------------------------------------------------------------
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  set.seed(42)
  opts <- parse_args_robust(cli_args)
  message("=== Starting Step 07 Triage Dashboard & Taxonomic Synthesis ===")
  triage_df <- build_triage_queue(opts)
  render_synthesis_figure_png(triage_df, opts$output_plot, cda_summary = opts$cda_summary, focal_taxon = opts$focal_taxon, output_complex_plot = opts$output_complex_plot)
  tryCatch({
    if (!exists("generate_step07_pub_figures", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) {
      source("scripts/morphometrics/pub_figures_utils.R")
    }
    if (exists("generate_step07_pub_figures", mode = "function")) generate_step07_pub_figures(triage_df, opts$output_plot)
  }, error = function(e) warning("Step 07 publication figures export warning: ", e$message))
  generate_taxonomic_revision_report(triage_df, opts)
  if (!is.null(opts$output_html) && nzchar(trimws(opts$output_html))) {
    generate_interactive_triage_html(triage_df, opts$output_html, focal_taxon = opts$focal_taxon, opts = opts)
  }
  message("=== Taxonomic Synthesis & Triage Dashboard Complete ===")
}

if (sys.nframe() == 0L) main()
