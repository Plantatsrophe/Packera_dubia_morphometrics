#!/usr/bin/env python3
"""Batched warm-R and static hygiene audit suite for the Step 04 refactor."""

from __future__ import annotations

import json, math, os, re, shutil, sqlite3, subprocess, sys, tempfile, unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db_schema import init_db_schema

SKIP_R_REASON = "Skipping batched R Step 04 audit hygiene tests under fast audit or missing Rscript."
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _is_finite(val: object) -> bool:
    if val is None or isinstance(val, bool):
        return False
    try:
        return math.isfinite(float(val))
    except (TypeError, ValueError):
        return False


def _seed_step04_inputs(tmp_p: Path) -> dict[str, Path]:
    (tmp_p / "custom_step04.yaml").write_text(
        "morphometrics:\n  correct_allometry: \"true\"\n  allometry_r2_threshold: 0.15\n  num_pcs: 6\n  max_k: 5\n  random_seed: 77\n",
        encoding="utf-8",
    )
    (tmp_p / "extracted_leaf_manifest.csv").write_text(
        "catalogNumber,plant_individual_id,leaf_id,leaf_zone,is_calibrated_mm,scale_source,petiole_length_mm,blade_length_mm,petiole_to_blade_ratio,blade_cie_l_norm,blade_cie_l_sd,blade_cie_a_browning_delta\n"
        "CAT_CAL,0,1,basal,TRUE,ruler_calibrated,10.0,30.0,0.33,0.45,0.05,1.2\n"
        "CAT_CAL,0,2,basal,TRUE,ruler_calibrated,14.0,40.0,0.35,0.47,0.07,1.6\n"
        "CAT_UNBOXED,0,1,basal,FALSE,tick_lattice_unboxed,20.0,50.0,0.40,0.38,0.08,4.0\n"
        "CAT_UNBOXED,0,2,basal,FALSE,tick_lattice_unboxed,24.0,60.0,0.40,0.40,0.10,5.0\n"
        "CAT_MP,0,1,basal,TRUE,mp_estimated,15.0,30.0,0.50,0.41,0.06,3.2\n",
        encoding="utf-8",
    )
    (tmp_p / "voucher_reproductive_metrics.csv").write_text(
        "catalogNumber,plant_individual_id,capitula_count,corymb_head_count_est,AR_inv\n"
        "REP_CORYMB,0,3,11,1.45\nREP_SOL,0,6,2,1.30\nREP_STERILE,0,0,0,\n",
        encoding="utf-8",
    )
    s04_dir = tmp_p / "s04_run"; s04_dir.mkdir(parents=True, exist_ok=True)
    db_path = tmp_p / "step04_audit.sqlite"; conn = sqlite3.connect(str(db_path)); init_db_schema(conn)
    v_rows, p_rows = [], []
    v_lines = ["catalogNumber,scientificName,species_standardized,determiner_tier,institutionCode,regional_group,decimalLatitude,decimalLongitude,is_primary_duplicate,pixels_per_mm,is_calibrated_mm,scale_source"]
    for idx in range(1, 25):
        cat, sp = f"S04_{idx:02d}", ("Packera aurea" if idx <= 12 else "Packera dubia")
        inst, rg = ("NCU" if idx % 2 else "DUKE"), ("Blue_Ridge" if ((idx - 1) // 6) % 2 == 0 else "Piedmont")
        lat, lon = round(35.2 + 0.08 * idx, 4), round(-82.5 + 0.12 * idx, 4)
        v_rows.append((cat, sp, sp, "Tier_1_Gold", inst, rg, lat, lon, 1, 10.0, 1)); p_rows.append((cat, 0))
        v_lines.append(f"{cat},{sp},{sp},Tier_1_Gold,{inst},{rg},{lat},{lon},1,10.0,1,ruler_calibrated")
    conn.executemany(
        "INSERT INTO vouchers (catalogNumber,scientificName,species_standardized,determiner_tier,institutionCode,regional_group,decimalLatitude,decimalLongitude,is_primary_duplicate,pixels_per_mm,is_calibrated_mm) VALUES (?,?,?,?,?,?,?,?,?,?,?);",
        v_rows,
    )
    conn.executemany("INSERT INTO plant_individuals (catalogNumber,plant_individual_id) VALUES (?,?);", p_rows)
    conn.commit(); conn.close()
    v_csv = s04_dir / "curated_vouchers.csv"; v_csv.write_text("\n".join(v_lines) + "\n", encoding="utf-8")
    return {"db_path": db_path, "v_csv": v_csv, "s04_dir": s04_dir}


_WARM_R_SCRIPT = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; tmp_root <- args[2L]; db_path <- args[3L]; v_csv <- args[4L]; s04_dir <- args[5L]
Sys.setenv(PLANTMATICS_SKIP_FIGURES = "1", PACKERA_SKIP_FIGURES = "1")
suppressPackageStartupMessages({ library(jsonlite) })
for (f in c("scripts/core/db_utils.R", "scripts/morphometrics/pub_figures_utils.R",
            "scripts/morphometrics/authority_weighting_utils.R", "scripts/morphometrics/plasticity_mitigation_utils.R",
            "scripts/morphometrics/gmm_icl_bootstrap_utils.R", "scripts/morphometrics/pub_ecoregion_manova_utils.R",
            "scripts/morphometrics/step04_allometry_merger_utils.R", "scripts/morphometrics/04_gmm_morphotools.R")) source(f)
res <- list()

# 1. Allometry heterogeneous-slope guard
n_a <- 32L; lc_v <- rep(seq(-1.5, 1.5, length.out = 16L), 2L); sp_a <- rep(c("Packera aurea", "Packera dubia"), each = 16L)
mk_allom_df <- function(s_a, s_d) {
  set.seed(101L); sl <- ifelse(sp_a == "Packera aurea", s_a, s_d)
  d <- data.frame(
    catalogNumber = sprintf("AL_%02d", seq_len(n_a)), plant_individual_id = 0L, species_standardized = sp_a, scientificName = sp_a,
    determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", solidity = 0.85, is_dissected = FALSE,
    is_calibrated_mm = TRUE, scale_source = "ruler_calibrated", decimalLatitude = seq(34.5, 38.5, length.out = n_a),
    decimalLongitude = seq(-83.5, -77.5, length.out = n_a), log_centsize = 3.0 + lc_v,
    A1 = sl * lc_v + rnorm(n_a, 0, 0.02), D1 = 0.8 * sl * lc_v + rnorm(n_a, 0, 0.02),
    A2 = 0.6 * sl * lc_v + rnorm(n_a, 0, 0.02), D2 = 0.5 * sl * lc_v + rnorm(n_a, 0, 0.02),
    A3 = ifelse(sp_a == "Packera aurea", 0.25, -0.25) + rnorm(n_a, 0, 0.03), D3 = rnorm(n_a, 0, 0.03), A4 = rnorm(n_a, 0, 0.03), D4 = rnorm(n_a, 0, 0.03), stringsAsFactors = FALSE
  )
  d$PC1 <- d$A1 + 0.5 * d$D1; d$PC2 <- d$A2 - 0.5 * d$D2; d
}
df_het <- mk_allom_df(0.90, 0.20); df_par <- mk_allom_df(0.55, 0.55)
r_ha <- test_and_correct_allometry(df_het, correct_mode = "auto", permutations = 299L, seed = 42L, r2_threshold = 0.10, output_report = file.path(tmp_root, "ha.csv"), output_plot = "", output_pairwise = "")
r_ht <- test_and_correct_allometry(df_het, correct_mode = "true", permutations = 299L, seed = 42L, r2_threshold = 0.10, output_report = file.path(tmp_root, "ht.csv"), output_plot = "", output_pairwise = "")
r_pa <- test_and_correct_allometry(df_par, correct_mode = "auto", permutations = 299L, seed = 42L, r2_threshold = 0.10, output_report = file.path(tmp_root, "pa.csv"), output_plot = "", output_pairwise = "")
res$allom <- list(het_auto_corr = isTRUE(r_ha$corrected), het_auto_r2 = as.numeric(r_ha$r2_size), het_auto_p_int = as.numeric(r_ha$p_interaction), het_auto_a1_diff = max(abs(r_ha$df$A1 - df_het$A1)), het_true_corr = isTRUE(r_ht$corrected), par_auto_corr = isTRUE(r_pa$corrected), par_auto_r2 = as.numeric(r_pa$r2_size), par_auto_p_int = as.numeric(r_pa$p_interaction))

# 2. Calibration gating, photometry/petiole manifest merge, and corymb head merge
cal_in <- data.frame(is_calibrated_mm = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, NA), scale_source = c("ruler_calibrated", "ruler_calibrated", "mp_estimated", "predicted_from_megapixels", "uncalibrated", "tick_lattice_unboxed", "tick_lattice_unboxed"), stringsAsFactors = FALSE)
agg_m <- merge_manifest_petiole_traits_by_individual(data.frame(catalogNumber = c("CAT_CAL", "CAT_UNBOXED", "CAT_MP"), plant_individual_id = 0L, stringsAsFactors = FALSE), file.path(tmp_root, "extracted_leaf_manifest.csv"))
agg_r <- merge_reproductive_traits_by_individual(data.frame(catalogNumber = c("REP_CORYMB", "REP_SOL", "REP_STERILE"), plant_individual_id = 0L, stringsAsFactors = FALSE), file.path(tmp_root, "voucher_reproductive_metrics.csv"))
res$calib <- list(cal_mask = as.logical(is_step04_calibrated_mm(cal_in)), man_pet_na = is.na(agg_m$petiole_length_mm), man_bld_na = is.na(agg_m$blade_length_mm), man_pbr = as.numeric(agg_m$petiole_to_blade_ratio), man_brown = as.numeric(agg_m$blade_cie_a_browning_delta), repro_log1p = as.numeric(agg_r$log1p_capitulum_count), repro_is_na = is.na(agg_r$log1p_capitulum_count))

# 3. Ecoregion centroid quarantine gate
set.seed(202L); n_e <- 28L; rg_e <- c(rep("Blue_Ridge", 11L), rep("Piedmont", 11L), rep("Middle_Atlantic_Coastal_Plain", 6L))
cent_e <- c(rep(FALSE, 8L), rep(TRUE, 3L), rep(FALSE, 8L), rep(TRUE, 3L), rep(FALSE, 2L), rep(TRUE, 4L))
df_eco <- data.frame(catalogNumber = sprintf("EC_%02d", seq_len(n_e)), plant_individual_id = 0L, species_standardized = rep(c("Packera aurea", "Packera dubia"), length.out = n_e), determiner_tier = "Tier_1_Gold", authority_weight = 2.0, morph_tier_source = "Basal_Rosette", regional_group = rg_e, is_centroid_or_inst_coord = cent_e, is_type_specimen = FALSE, decimalLatitude = c(rep(35.6, 11L), rep(35.9, 11L), rep(35.2, 6L)), decimalLongitude = c(rep(-82.5, 11L), rep(-79.0, 11L), rep(-77.5, 6L)), PC1 = c(rnorm(11L, -2.5, 0.25), rnorm(11L, 2.5, 0.25), rnorm(6L, 0.0, 0.25)), PC2 = c(rnorm(22L, 0.0, 0.25), rnorm(6L, 2.0, 0.25)), PC3 = rnorm(n_e, 0.0, 0.25), stringsAsFactors = FALSE)
eco_res <- run_supervised_ecoregion_cda(df_eco, c("PC1", "PC2", "PC3"), output_summary_csv = file.path(tmp_root, "eco_sum.csv"), permutations = 99L, seed = 42L)
man_res <- run_multivariate_manova_audit(df_eco, c("PC1", "PC2", "PC3"), output_group_csv = file.path(tmp_root, "man_grp.csv"), output_pairwise_csv = file.path(tmp_root, "man_pw.csv"), include_gmm_clusters = FALSE)
rg_row <- man_res$group_table[man_res$group_table$grouping_factor == "regional_group", , drop = FALSE]
res$eco <- list(cent_active_n = sum(eco_res$is_active[cent_e]), cent_can1_fin = all(is.finite(eco_res$ecoregion_can1[cent_e])), active_ecos = as.character(eco_res$active_ecoregions), sum_max_active = max(eco_res$summary_table$n_active_anchors), man_rg_groups = as.integer(rg_row$n_groups[1L]), man_rg_vouchers = as.integer(rg_row$n_active_vouchers[1L]))

# 4. Deme CV fold-internal median imputation & outlier trimming
set.seed(303L); n_cv <- 32L; y_cv <- rep(c(rep("Packera aurea", 4L), rep("Packera dubia", 4L)), 4L)
coords_cv <- data.frame(catalogNumber = sprintf("CV_%02d", seq_len(n_cv)), cv_deme_id = rep(c("D1", "D2", "D3", "D4"), each = 8L), stringsAsFactors = FALSE)
X_cln <- cbind(PC1 = ifelse(y_cv == "Packera aurea", -2.5, 2.5) + rnorm(n_cv, 0, 0.25), PC2 = rnorm(n_cv, 0, 0.35), PC3 = rnorm(n_cv, 0, 0.35), V_aux = rnorm(n_cv, 0, 0.35))
X_cln[c(2L, 3L), c("PC2", "PC3", "V_aux")] <- NA_real_; X_cln[c(10L, 11L), "PC3"] <- NA_real_
X_cand_cln <- cbind(X_cln, V_noise = rnorm(n_cv, 0, 0.35)); X_cand_cln[2L, "V_noise"] <- NA_real_
X_spk <- X_cln; X_spk[1L, c("PC2", "PC3", "V_aux")] <- 9e6; X_cand_spk <- X_cand_cln; X_cand_spk[1L, c("PC2", "PC3", "V_aux", "V_noise")] <- 9e6
X_pc1_spk <- X_cln; X_pc1_spk[1L, "PC1"] <- 9e6
fv_cv <- build_taxon_blind_gathering_deme_folds(coords_cv, n_obs = n_cv, deme_radius_km = 10.0); tr1 <- which(fv_cv != fv_cv[1L])
res$cv <- list(
  acc_clean = as.numeric(compute_weighted_shrinkage_cda_loocv(X_a = X_cln, y_a = y_cv, coords_df = coords_cv, X_cand = X_cand_cln, allow_supervised_f = TRUE, trim_fold_outliers = TRUE)),
  acc_spiked = as.numeric(compute_weighted_shrinkage_cda_loocv(X_a = X_spk, y_a = y_cv, coords_df = coords_cv, X_cand = X_cand_spk, allow_supervised_f = TRUE, trim_fold_outliers = TRUE)),
  acc_pc1_trim = as.numeric(compute_weighted_shrinkage_cda_loocv(X_a = X_pc1_spk, y_a = y_cv, coords_df = coords_cv, allow_supervised_f = TRUE, trim_fold_outliers = TRUE)),
  acc_pc1_notrim = as.numeric(compute_weighted_shrinkage_cda_loocv(X_a = X_pc1_spk, y_a = y_cv, coords_df = coords_cv, allow_supervised_f = TRUE, trim_fold_outliers = FALSE)),
  acc_unsup = as.numeric(compute_weighted_shrinkage_cda_loocv(X_a = X_cln, y_a = y_cv, coords_df = coords_cv, X_cand = X_cand_cln, allow_supervised_f = FALSE, trim_fold_outliers = TRUE)),
  med_max_diff = max(abs(apply(X_cand_spk[tr1, , drop = FALSE], 2L, stats::median, na.rm = TRUE) - apply(X_cand_cln[tr1, , drop = FALSE], 2L, stats::median, na.rm = TRUE)))
)

# 5. YAML config loader + n=96 skewed 2-cluster GMM & sequential bootstrap step-down
cfg_def <- load_step04_morphometrics_config("config/config.yaml"); cfg_cus <- load_step04_morphometrics_config(file.path(tmp_root, "custom_step04.yaml"))
set.seed(4L)
df_gmm <- data.frame(catalogNumber = sprintf("GM_%02d", 1:96), plant_individual_id = 0L, species_standardized = "Packera dubia", determiner_tier = "Tier_1_Gold", morph_tier_source = "Basal_Rosette", solidity = 0.82, is_dissected = FALSE, exsiccatae_key = paste0("EX_", 1:96), PC1 = c(rnorm(48L, -4.0, 0.85), 4.0 + (rgamma(48L, shape = 2.0, scale = 0.8) - 1.6)), PC2 = rnorm(96L, 0.0, 0.85), PC3 = rnorm(96L, 0.0, 0.85), stringsAsFactors = FALSE)
gmm_skew <- run_gmm_cluster_analysis(df_gmm, c("PC1", "PC2", "PC3"), max_k = 4L, use_clustvarsel = FALSE, intraspecific_mode = TRUE, seed = 42L)
X_g3 <- as.matrix(df_gmm[, c("PC1", "PC2", "PC3")]); suppressPackageStartupMessages(library(mclust))
f_g3 <- mclust::Mclust(X_g3, G = 1:3, verbose = FALSE); bm3 <- as.matrix(f_g3$BIC); bt3 <- data.frame(K = 1:3, BIC = round(apply(bm3, 1L, max, na.rm = TRUE), 2))
for (mn in colnames(bm3)) bt3[[paste0("BIC_", mn)]] <- round(as.numeric(bm3[, mn]), 2)
pr3 <- predict(f_g3, newdata = X_g3)
g_seq <- apply_gmm_icl_bootstrap_gate(X_g3, X_g3, bt3, best_k = 3L, best_name = f_g3$modelName, cls_sub = as.integer(pr3$classification), unc_sub = as.numeric(pr3$uncertainty), is_intra_gmm = TRUE, seed = 42L)
res$gmm_cfg <- list(cfg_def = cfg_def, cfg_cus = cfg_cus, has_pair_fn = exists("run_gmm_pair_bootstrap_lrt", mode = "function"), has_gate_fn = exists("apply_gmm_icl_bootstrap_gate", mode = "function"), best_k = as.integer(gmm_skew$best_k), boot_p2 = as.numeric(gmm_skew$bic_table$Boot_LRT_P_prev[gmm_skew$bic_table$K == 2L]), seq_k = as.integer(g_seq$best_k), seq_name = as.character(g_seq$best_name), seq_p_prev = as.numeric(g_seq$bic_table$Boot_LRT_P_prev[g_seq$bic_table$K %in% c(2L, 3L)]))

# 6. Ecoregion-controlled marginal drop-1 batch ANCOVA & aliased safety
set.seed(404L); rg_u <- c(rep("Blue_Ridge", 20L), rep("Piedmont", 20L)); inst_u <- c(rep("NCU", 16L), rep("DUKE", 4L), rep("NCU", 4L), rep("DUKE", 16L)); sp_u <- rep(c("Packera aurea", "Packera dubia"), 20L)
df_u <- data.frame(catalogNumber = sprintf("BU_%02d", 1:40), species_standardized = sp_u, institutionCode = inst_u, regional_group = rg_u, authority_weight = 2.0, is_passive = FALSE, morph_tier_source = "Basal_Rosette", PC1 = ifelse(rg_u == "Blue_Ridge", -1.3, 1.3) + ifelse(sp_u == "Packera aurea", -0.2, 0.2) + rnorm(40L, 0, 0.15), PC2 = rnorm(40L, 0, 0.2), stringsAsFactors = FALSE)
bu <- audit_institutional_batch_effects(df_u, output_csv = file.path(tmp_root, "bu.csv"), min_obs = 5L, selected_features = c("PC1", "PC2"))
bu_pc1 <- bu$anova_table[bu$anova_table$Response == "PC1" & bu$anova_table$Term == "institutionCode", , drop = FALSE]
df_al <- data.frame(catalogNumber = sprintf("BA_%02d", 1:32), species_standardized = rep(c("Packera aurea", "Packera dubia"), 16L), institutionCode = rep(c("NCU", "LSU"), each = 16L), regional_group = rep(c("Blue_Ridge", "South_Central_Plains"), each = 16L), authority_weight = 2.0, is_passive = FALSE, morph_tier_source = "Basal_Rosette", PC1 = c(rnorm(16L, -1.3, 0.2), rnorm(16L, 1.3, 0.2)), PC2 = rnorm(32L, 0, 0.25), stringsAsFactors = FALSE)
ba <- audit_institutional_batch_effects(df_al, output_csv = file.path(tmp_root, "ba.csv"), min_obs = 5L, selected_features = c("PC1", "PC2"))
ba_inst <- ba$anova_table[ba$anova_table$Term == "institutionCode", , drop = FALSE]
res$batch <- list(uneven_eta = as.numeric(bu_pc1$partial_eta_sq[1L]), uneven_status = as.character(bu_pc1$Batch_Effect_Status[1L]), alias_rows = nrow(ba_inst), alias_responses = as.character(ba_inst$Response), alias_eta_fin = all(is.finite(as.numeric(ba_inst$partial_eta_sq))))

# 7. SQLite-only Step 04 rerun idempotency after deleting harmonics.csv
set.seed(505L); n_id <- 24L; lc_id <- seq(-1.0, 1.0, length.out = n_id); sp_id <- rep(c("Packera aurea", "Packera dubia"), each = 12L)
cs_id <- round(10^(3.0 + lc_id), 4); lcs_id <- round(log10(cs_id), 6)
A1_v <- round(0.55 * lc_id + ifelse(sp_id == "Packera aurea", 0.10, -0.10) + rnorm(n_id, 0, 0.015), 6); D1_v <- round(0.45 * lc_id + rnorm(n_id, 0, 0.015), 6)
A2_v <- round(0.35 * lc_id + ifelse(sp_id == "Packera aurea", 0.30, -0.30) + rnorm(n_id, 0, 0.015), 6); D2_v <- round(0.25 * lc_id + rnorm(n_id, 0, 0.015), 6)
A3_v <- round(ifelse(sp_id == "Packera aurea", 0.45, -0.45) + rnorm(n_id, 0, 0.02), 6); D3_v <- round(rnorm(n_id, 0, 0.02), 6); A4_v <- round(rnorm(n_id, 0, 0.02), 6); D4_v <- round(rnorm(n_id, 0, 0.02), 6)
sol_v <- round(0.84 + 0.03 * lc_id + ifelse(sp_id == "Packera aurea", 0.02, -0.02), 4)
pc_m <- round(stats::prcomp(cbind(A1_v, D1_v, A2_v, D2_v, A3_v, D3_v, A4_v, D4_v))$x[, 1:5], 6)
spec_df <- data.frame(catalogNumber = sprintf("S04_%02d", 1:n_id), plant_individual_id = 0L, species_standardized = sp_id, scientificName = sp_id, determiner_tier = "Tier_1_Gold", institutionCode = rep(c("NCU", "DUKE"), length.out = n_id), regional_group = rep(c(rep("Blue_Ridge", 6L), rep("Piedmont", 6L)), 2L), decimalLatitude = round(35.2 + 0.08 * seq_len(n_id), 4), decimalLongitude = round(-82.5 + 0.12 * seq_len(n_id), 4), pixels_per_mm = 10.0, is_calibrated_mm = TRUE, scale_source = "ruler_calibrated", centroid_size = cs_id, log_centsize = lcs_id, solidity = sol_v, dissection_index = 0.15, aspect_ratio = 1.6, morph_tier_source = "Basal_Rosette", A1 = A1_v, D1 = D1_v, A2 = A2_v, D2 = D2_v, A3 = A3_v, D3 = D3_v, A4 = A4_v, D4 = D4_v, PC1 = pc_m[, 1L], PC2 = pc_m[, 2L], PC3 = pc_m[, 3L], PC4 = pc_m[, 4L], PC5 = pc_m[, 5L], stringsAsFactors = FALSE)
leaf_df <- spec_df; leaf_df$leaf_id <- 1L; leaf_df$leaf_zone <- "basal"; leaf_df$B1 <- 0; leaf_df$C1 <- 0; leaf_df$B2 <- 0; leaf_df$C2 <- 0
harm_csv <- file.path(s04_dir, "leaf_efa_harmonics.csv"); write.csv(spec_df, harm_csv, row.names = FALSE)
Sys.setenv(PLANTMATICS_DB_PATH = db_path, PACKERA_DB_PATH = db_path, PLANTMATICS_RUN_DIR = s04_dir, PACKERA_RUN_DIR = s04_dir, RUN_ID = "run_s04_idemp")
ok_p3 <- persist_step03_to_sqlite(efa_df = leaf_df, specimen_df = spec_df, vouchers_path = v_csv, run_id = "run_s04_idemp")
mk_opts <- function(tag) {
  d <- file.path(s04_dir, tag); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  list(vouchers = v_csv, harmonics = harm_csv, reproductive = "", outline_only = FALSE, correct_allometry = "auto", allometry_report = file.path(d, "allom.csv"), allometry_pairwise = "", allometry_plot = file.path(d, "fig04a.png"), use_clustvarsel = FALSE, output_flags = file.path(d, "flags.csv"), output_plot = file.path(d, "fig04b.png"), output_report = file.path(d, "gmm.csv"), output_batch_audit = file.path(d, "batch.csv"), max_k = 4L, num_pcs = 5L, permutations = 199L, seed = 42L, intraspecific_mode = FALSE, focal_taxon = "")
}
ins_allom <- NULL
out1 <- run_gmm_morphotools_pipeline(mk_opts("r1")); unlink(harm_csv); harm_gone <- !file.exists(harm_csv)
out2 <- run_gmm_morphotools_pipeline(mk_opts("r2"))
Sys.unsetenv(c("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH", "PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR", "RUN_ID", "PLANTMATICS_SKIP_FIGURES", "PACKERA_SKIP_FIGURES"))
Sys.setenv(PLANTMATICS_DEMO_MODE = "1")
ins_allom <- test_and_correct_allometry(spec_df[1:6, , drop = FALSE], output_report = file.path(tmp_root, "insuff_allom.csv"), output_plot = file.path(tmp_root, "insuff_04a.png"))
spec_ph <- spec_df; spec_ph$blade_cie_l_norm <- round(stats::rnorm(n_id, 0.70, 0.08), 4); spec_ph$blade_cie_l_sheet_delta <- round(rep(c(-32.0, -8.0, -28.0, -12.0), length.out = n_id), 4)
harm_ph <- file.path(s04_dir, "harm_ph.csv"); write.csv(spec_ph, harm_ph, row.names = FALSE); o_ph <- mk_opts("r_ph"); o_ph$harmonics <- harm_ph
out_ph <- run_gmm_morphotools_pipeline(o_ph); bic_sel <- read.csv(file.path(s04_dir, "r_ph", "gmm_bic_model_selection.csv"), stringsAsFactors = FALSE)
Sys.unsetenv("PLANTMATICS_DEMO_MODE")
out1 <- out1[order(out1$catalogNumber, out1$plant_individual_id), , drop = FALSE]; out2 <- out2[order(out2$catalogNumber, out2$plant_individual_id), , drop = FALSE]
al_r1 <- read.csv(file.path(s04_dir, "r1", "allom.csv"), stringsAsFactors = FALSE)
res$w3 <- list(insuff_csv = file.exists(file.path(tmp_root, "insuff_allom.csv")), insuff_png = file.exists(file.path(tmp_root, "insuff_04a.png")), r1_png = file.exists(file.path(s04_dir, "r_ph", "fig04a.png")), bic_cols = all(c("optimal_K_raw", "optimal_K_merged", "selected_model") %in% names(bic_sel)), has_sheet_delta = "blade_cie_l_sheet_delta_std" %in% names(out_ph), sel_feats = as.character(out_ph$selected_morph_features[1L]))
res$rerun <- list(ok_p3 = isTRUE(ok_p3), harm_gone = harm_gone, n1 = nrow(out1), n2 = nrow(out2), r2_size = as.numeric(al_r1$R2[grepl("^log_centsize$", al_r1$Term)][1L]), slopes_homog = isTRUE(as.logical(al_r1$Slopes_Homogeneous[grepl(":", al_r1$Term)][1L])), a1_changed = max(abs(out1$A1 - spec_df$A1[match(out1$catalogNumber, spec_df$catalogNumber)])), diff_pc1 = max(abs(out1$PC1 - out2$PC1)), diff_a1 = max(abs(out1$A1 - out2$A1)), diff_sol = max(abs(out1$solidity - out2$solidity)), diff_prob = max(abs(out1$cda_posterior_prob - out2$cda_posterior_prob)))
jsonlite::write_json(res, out_json, auto_unbox = TRUE, na = "null", digits = 10)
"""


class TestStep04AuditHygiene(unittest.TestCase):
    """Contract and hygiene audit suite for the refactored Step 04 morphometrics stack."""

    results: dict = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            cls.results = {}
            return
        cls._tmp_dir = tempfile.TemporaryDirectory(prefix="step04_hyg_")
        tmp_p = Path(cls._tmp_dir.name)
        fx = _seed_step04_inputs(tmp_p)
        r_path, out_json = tmp_p / "run_step04_hygiene.R", tmp_p / "results.json"
        r_path.write_text(_WARM_R_SCRIPT, encoding="utf-8")
        proc = subprocess.run(
            ["Rscript", str(r_path), str(out_json), str(tmp_p), str(fx["db_path"]), str(fx["v_csv"]), str(fx["s04_dir"])],
            cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False, timeout=180,
        )
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(f"Batched Step 04 R failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls.results = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_allometry_heterogeneous_slopes_auto_guard(self) -> None:
        al = self.results["allom"]
        self.assertGreater(float(al["het_auto_r2"]), 0.25)
        self.assertLess(float(al["het_auto_p_int"]), 0.01)
        self.assertFalse(bool(al["het_auto_corr"]))
        self.assertEqual(float(al["het_auto_a1_diff"]), 0.0)
        self.assertTrue(bool(al["het_true_corr"]))
        self.assertGreater(float(al["par_auto_p_int"]), 0.20)
        self.assertTrue(bool(al["par_auto_corr"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_calibration_photometry_and_corymb_handshakes(self) -> None:
        cb = self.results["calib"]
        self.assertEqual(cb["cal_mask"], [True, False, False, False, False, False, False])
        self.assertEqual(cb["man_pet_na"], [False, True, True])
        self.assertEqual(cb["man_bld_na"], [False, True, True])
        self.assertTrue(all(_is_finite(v) for v in cb["man_pbr"] + cb["man_brown"]))
        self.assertAlmostEqual(float(cb["man_brown"][1]), 4.5, places=4)
        self.assertAlmostEqual(float(cb["repro_log1p"][0]), math.log1p(11.0), places=4)
        self.assertAlmostEqual(float(cb["repro_log1p"][1]), math.log1p(6.0), places=4)
        self.assertEqual(cb["repro_is_na"], [False, False, True])

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_ecoregion_centroid_quarantine_gate(self) -> None:
        ec = self.results["eco"]
        self.assertEqual(int(ec["cent_active_n"]), 0)
        self.assertTrue(bool(ec["cent_can1_fin"]))
        self.assertEqual(sorted(ec["active_ecos"]), ["Blue_Ridge", "Piedmont"])
        self.assertLessEqual(int(ec["sum_max_active"]), 8)
        self.assertEqual(int(ec["man_rg_groups"]), 2)
        self.assertLessEqual(int(ec["man_rg_vouchers"]), 16)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_deme_cv_zero_imputation_leakage_and_no_mass_lda(self) -> None:
        cv = self.results["cv"]
        self.assertEqual(float(cv["med_max_diff"]), 0.0)
        for k in ("acc_clean", "acc_spiked", "acc_pc1_trim", "acc_unsup"):
            self.assertTrue(_is_finite(cv[k]))
            self.assertGreaterEqual(float(cv[k]), 0.85)
        self.assertGreater(float(cv["acc_pc1_trim"]), float(cv["acc_pc1_notrim"]))
        for r_file in sorted((PROJECT_ROOT / "scripts" / "morphometrics").glob("*.R")):
            txt = "\n".join(ln for ln in r_file.read_text(encoding="utf-8").splitlines() if not ln.lstrip().startswith("#"))
            self.assertIsNone(re.search(r"MASS::lda|\blda\s*\(", txt), f"Forbidden MASS::lda in {r_file.name}")

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_locked_sequential_bootstrap_lrt_and_yaml_config(self) -> None:
        gc = self.results["gmm_cfg"]
        self.assertEqual(gc["cfg_def"]["correct_allometry"], "auto")
        self.assertAlmostEqual(float(gc["cfg_def"]["allometry_r2_threshold"]), 0.10, places=6)
        self.assertEqual((int(gc["cfg_def"]["num_pcs"]), int(gc["cfg_def"]["max_k"]), int(gc["cfg_def"]["random_seed"])), (5, 8, 42))
        self.assertEqual(gc["cfg_cus"]["correct_allometry"], "true")
        self.assertAlmostEqual(float(gc["cfg_cus"]["allometry_r2_threshold"]), 0.15, places=6)
        self.assertEqual((int(gc["cfg_cus"]["num_pcs"]), int(gc["cfg_cus"]["max_k"]), int(gc["cfg_cus"]["random_seed"])), (6, 5, 77))
        self.assertTrue(bool(gc["has_pair_fn"]) and bool(gc["has_gate_fn"]))
        self.assertEqual(int(gc["best_k"]), 2)
        self.assertTrue(_is_finite(gc["boot_p2"]) and float(gc["boot_p2"]) < 0.05)
        self.assertEqual(int(gc["seq_k"]), 2)
        self.assertIn("_SeqBootStepped", str(gc["seq_name"]))
        self.assertTrue(all(_is_finite(p) for p in gc["seq_p_prev"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_ecoregion_controlled_batch_ancova(self) -> None:
        bt = self.results["batch"]
        self.assertLess(float(bt["uneven_eta"]), 0.05)
        self.assertIn("Pass", str(bt["uneven_status"]))
        self.assertGreaterEqual(int(bt["alias_rows"]), 2)
        self.assertEqual(sorted(bt["alias_responses"]), ["PC1", "PC2"])
        self.assertTrue(bool(bt["alias_eta_fin"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_sqlite_only_step04_rerun_idempotency(self) -> None:
        rr = self.results["rerun"]
        self.assertTrue(bool(rr["ok_p3"]) and bool(rr["harm_gone"]))
        self.assertEqual((int(rr["n1"]), int(rr["n2"])), (24, 24))
        self.assertGreaterEqual(float(rr["r2_size"]), 0.10)
        self.assertTrue(bool(rr["slopes_homog"]))
        self.assertGreater(float(rr["a1_changed"]), 1e-4)
        for dk in ("diff_pc1", "diff_a1", "diff_sol", "diff_prob"):
            self.assertLess(float(rr[dk]), 1e-6, f"Non-idempotent SQLite rerun delta on {dk}: {rr[dk]}")
        w3 = self.results["w3"]
        self.assertTrue(bool(w3["insuff_csv"]) and bool(w3["insuff_png"]) and bool(w3["r1_png"]))
        self.assertTrue(bool(w3["bic_cols"]) and bool(w3["has_sheet_delta"]))
        self.assertNotIn("blade_cie_l_norm_std", str(w3["sel_feats"]).split(";"))

    def test_pruned_dead_paths_and_line_ceilings(self) -> None:
        ceilings = {
            PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R": 500,
            PROJECT_ROOT / "scripts" / "morphometrics" / "step04_allometry_merger_utils.R": 380,
            PROJECT_ROOT / "scripts" / "morphometrics" / "gmm_icl_bootstrap_utils.R": 280,
            PROJECT_ROOT / "scripts" / "morphometrics" / "pub_synthesis_plates.R": 120,
            PROJECT_ROOT / "scripts" / "figures" / "step04_figures.R": 500,
            Path(__file__).resolve(): 380,
        }
        for path, cap in ceilings.items():
            self.assertTrue(path.is_file(), f"Missing expected Step 04 file: {path}")
            n_lines = len(path.read_text(encoding="utf-8").splitlines())
            self.assertLess(n_lines, cap, f"{path.name}: {n_lines} lines >= {cap}")
        self.assertNotIn("2.0 * iqr_d", (PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R").read_text(encoding="utf-8"))
        self.assertNotIn("asym_names", (PROJECT_ROOT / "scripts" / "morphometrics" / "step04_allometry_merger_utils.R").read_text(encoding="utf-8"))
        self.assertNotIn("unlink(", (PROJECT_ROOT / "scripts" / "figures" / "step04_figures.R").read_text(encoding="utf-8"))
        dead_tokens = ("_" + "archive", "MASS::lda", "render_allometry_" + "diagnostic_pdf", "generate_cda_" + "biplot_pdf")
        contract_fns = (
            "load_step04_morphometrics_config", "is_step04_calibrated_mm", "merge_manifest_petiole_traits_by_individual",
            "merge_reproductive_traits_by_individual", "test_and_correct_allometry", "compute_weighted_shrinkage_cda_loocv",
            "apply_gmm_icl_bootstrap_gate", "run_gmm_pair_bootstrap_lrt", "run_supervised_ecoregion_cda",
            "run_multivariate_manova_audit", "run_gmm_morphotools_pipeline",
        )
        fn_defs: dict[str, list[str]] = {fn: [] for fn in contract_fns}
        r_paths = list(sorted((PROJECT_ROOT / "scripts" / "morphometrics").glob("*.R"))) + [PROJECT_ROOT / "scripts" / "figures" / "step04_figures.R"]
        for r_file in r_paths:
            txt = r_file.read_text(encoding="utf-8")
            for tok in dead_tokens:
                self.assertNotIn(tok, txt, f"Dead token {tok} found in {r_file.name}")
            self.assertIsNone(re.search(r"outputs/figures/[^\"'\s]+\.pdf", txt), f"PDF figure path in {r_file.name}")
            for fn in contract_fns:
                if re.search(rf"^{re.escape(fn)}\s*<-\s*function\(", txt, flags=re.MULTILINE):
                    fn_defs[fn].append(r_file.name)
        for fn, owners in fn_defs.items():
            self.assertEqual(len(owners), 1, f"Expected single definition for {fn}, found {owners}")


if __name__ == "__main__":
    unittest.main()
