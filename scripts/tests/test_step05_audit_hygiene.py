#!/usr/bin/env python3
"""Batched warm-R and static hygiene audit suite for the Step 05 macroecology stack."""

from __future__ import annotations

import json, math, os, re, shutil, sqlite3, subprocess, sys, tempfile, unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db_schema import init_db_schema

SKIP_R_REASON = "Skipping batched R Step 05 audit hygiene tests under fast audit or missing Rscript."
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _is_finite(val: object) -> bool:
    """Return True if val is a non-boolean finite scalar."""
    if val is None or isinstance(val, bool):
        return False
    try:
        return math.isfinite(float(val))
    except (TypeError, ValueError):
        return False


def _seed_step05_db(tmp_p: Path) -> Path:
    """Create an isolated SQLite DB seeded with 8 clean anthesis vouchers and empty temp subdirs."""
    db_path = tmp_p / "step05_audit.sqlite"
    conn = sqlite3.connect(str(db_path))
    init_db_schema(conn)
    v_rows, p_rows, ph_rows = [], [], []
    for idx in range(1, 9):
        cat = f"S05_{idx:02d}"
        lat, lon = round(32.0 + (idx - 1) * (10.0 / 7.0), 4), round(-84.0 + (idx - 1) * (8.0 / 7.0), 4)
        doy = round(90.0 + 4.0 * (lat - 32.0), 1)
        v_rows.append((cat, "NCU", "Packera dubia", "Packera dubia", "Tier_1_Gold", lat, lon, 100.0, doy, 1, 0))
        p_rows.append((cat, 0))
        ph_rows.append((cat, 0, "Peak_Anthesis", 3, 5, 1))
    conn.executemany(
        "INSERT INTO vouchers (catalogNumber,institutionCode,scientificName,species_standardized,determiner_tier,decimalLatitude,decimalLongitude,coordinateUncertainty,doy,is_primary_duplicate,is_centroid_or_inst_coord) VALUES (?,?,?,?,?,?,?,?,?,?,?);",
        v_rows,
    )
    conn.executemany("INSERT INTO plant_individuals (catalogNumber,plant_individual_id) VALUES (?,?);", p_rows)
    conn.executemany(
        "INSERT INTO voucher_reproductive_phenology (catalogNumber,plant_individual_id,voucher_phenological_state,n_anthesis,total_capitula,is_verified_anthesis) VALUES (?,?,?,?,?,?);",
        ph_rows,
    )
    conn.commit()
    conn.close()
    for sub in ("sdm_env", "sdm_out", "cache_env", "cache_ssurgo"):
        (tmp_p / sub).mkdir(parents=True, exist_ok=True)
    return db_path


_WARM_R_SCRIPT = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; tmp_root <- args[2L]; db_path <- args[3L]
Sys.setenv(PLANTMATICS_SKIP_FIGURES = "1", PACKERA_SKIP_FIGURES = "1")
suppressPackageStartupMessages({ library(jsonlite) })
for (f in c("scripts/morphometrics/pub_figures_utils.R", "scripts/morphometrics/authority_weighting_utils.R", "scripts/core/db_utils.R",
            "scripts/analysis/step05_config_geo_utils.R", "scripts/analysis/step05_layer_hygiene_utils.R", "scripts/analysis/step05_env_outlier_utils.R",
            "scripts/analysis/step05_circular_hopkins_models.R", "scripts/analysis/edaphic_geochem_wetland_utils.R", "scripts/analysis/edaphic_fine_grain_utils.R",
            "scripts/analysis/step05_phenology_ecotypes.R", "scripts/analysis/05_multimodal_spatial_rf.R", "scripts/analysis/step05_sdm_maxent.R")) source(f)
res <- list()

# 1. Post-selection ecotype quarantine, Warren p-value guard, and fold re-clustering
set.seed(101L); n_ec <- 28L
df_ec <- data.frame(
  catalogNumber = sprintf("EC_%02d", seq_len(n_ec)), plant_individual_id = 0L, species_standardized = "Packera dubia",
  determiner_tier = "Tier_1_Gold", is_type_specimen = FALSE, is_centroid_or_inst_coord = c(rep(FALSE, 24L), TRUE, TRUE, FALSE, FALSE),
  is_confident_label_noise = c(rep(FALSE, 26L), TRUE, TRUE), is_training_eligible_coord = TRUE, coordinateUncertaintyInMeters = 100.0,
  latitude = c(seq(33.8, 34.3, length.out = 12L), seq(35.6, 36.1, length.out = 12L), 39.0, 39.2, 31.0, 31.2),
  longitude = c(seq(-77.9, -77.4, length.out = 12L), seq(-80.3, -79.8, length.out = 12L), -88.0, -88.2, -75.0, -75.2),
  bio1_temp_mean = c(rnorm(12L, 18.6, 0.2), rnorm(12L, 14.2, 0.2), -10.0, 35.0, -8.0, 34.0),
  bio12_precip_annual = c(rnorm(12L, 1390.0, 15.0), rnorm(12L, 1110.0, 15.0), 100.0, 3900.0, 150.0, 3800.0),
  fine_soil_ph = c(rnorm(12L, 7.0, 0.12), rnorm(12L, 4.8, 0.12), 2.5, 9.5, 2.8, 9.2),
  fine_bedrock_depth_cm = c(rnorm(12L, 165.0, 5.0), rnorm(12L, 18.0, 3.0), 200.0, 1.0, 200.0, 2.0),
  ssurgo_rock_outcrop_pct = c(runif(12L, 0.0, 1.5), rnorm(12L, 58.0, 4.0), 98.0, 0.0, 95.0, 0.0),
  topo_slope_deg = c(rnorm(12L, 1.4, 0.2), rnorm(12L, 16.5, 1.2), 45.0, 0.1, 44.0, 0.2),
  topo_twi_30m = c(rnorm(12L, 9.2, 0.3), rnorm(12L, 3.1, 0.3), 20.0, 0.5, 19.0, 0.6), stringsAsFactors = FALSE
)
df_ec$decimalLatitude <- df_ec$latitude; df_ec$decimalLongitude <- df_ec$longitude
res_ec <- resolve_analysis_group(df_ec, intraspecific_mode = TRUE, allow_unimodal_null = FALSE)
w_ec <- run_warren_niche_identity_tests(res_ec, n_perm = 29L, intraspecific_mode = TRUE)$summary
df_2sp <- df_ec[1:24, , drop = FALSE]; df_2sp$species_standardized <- rep(c("Packera aurea", "Packera dubia"), each = 12L)
w_2sp <- run_warren_niche_identity_tests(df_2sp, n_perm = 19L, intraspecific_mode = FALSE)$summary

set.seed(42L)
mk_sep <- function(n_per) {
  X <- rbind(cbind(rnorm(n_per, -4, 0.5), rnorm(n_per, 0, 0.5), rnorm(n_per, -3, 0.5), rnorm(n_per, 0, 0.5), matrix(rnorm(n_per * 2L, 0, 0.5), n_per, 2L)),
             cbind(rnorm(n_per, 0, 0.5), rnorm(n_per, 4, 0.5), rnorm(n_per, 0, 0.5), rnorm(n_per, 3, 0.5), matrix(rnorm(n_per * 2L, 0, 0.5), n_per, 2L)),
             cbind(rnorm(n_per, 4, 0.5), rnorm(n_per, -4, 0.5), rnorm(n_per, 3, 0.5), rnorm(n_per, -3, 0.5), matrix(rnorm(n_per * 2L, 0, 0.5), n_per, 2L)))
  colnames(X) <- paste0("V", 1:6); list(X = X, y = rep(c("Ecotype_1", "Ecotype_2", "Ecotype_3"), each = n_per))
}
sep_tr <- mk_sep(20L); sep_te <- mk_sep(10L)
raw_sep <- compute_fold_reclustered_ecotype_ref(sep_tr$X, sep_tr$y, sep_te$X, seed = 42L)
pred_sep <- if (is.list(raw_sep) && !is.data.frame(raw_sep)) raw_sep[[intersect(c("pred", "predictions", "labels", "y_pred"), names(raw_sep))[1L]]] else raw_sep
mk_uni <- function(n_u) {
  X <- matrix(rnorm(n_u * 6L), nrow = n_u, ncol = 6L); colnames(X) <- paste0("V", 1:6); ang <- atan2(X[, 2L], X[, 1L])
  list(X = X, y = ifelse(ang < -pi / 3.0, "S1", ifelse(ang < pi / 3.0, "S2", "S3")))
}
uni_tr <- mk_uni(90L); uni_te <- mk_uni(45L)
raw_uni <- compute_fold_reclustered_ecotype_ref(uni_tr$X, uni_tr$y, uni_te$X, seed = 42L)
pred_uni <- if (is.list(raw_uni) && !is.data.frame(raw_uni)) raw_uni[[intersect(c("pred", "predictions", "labels", "y_pred"), names(raw_uni))[1L]]] else raw_uni
g_ec <- as.character(res_ec$analysis_group); u_cln <- unique(g_ec[1:24])
res$ecotype <- list(
  n_groups = length(unique(g_ec)), all_ecotype_prefix = all(grepl("^Ecotype_", g_ec)),
  grp1_uniform = length(unique(g_ec[1:12])) == 1L, grp2_uniform = length(unique(g_ec[13:24])) == 1L, grps_distinct = g_ec[1L] != g_ec[13L],
  quarantined_passive = all(!is.na(g_ec[25:28]) & (g_ec[25:28] %in% u_cln)), attr_post_sel = isTRUE(attr(res_ec, "is_post_selection_ecotype")),
  ec_pvals = as.list(w_ec[1L, c("p_equiv_D_clim", "p_bg_D_clim", "p_equiv_D_edaph", "p_bg_D_edaph", "p_value_D", "p_value_I", "p_bg_D_joint"), drop = FALSE]),
  ec_dvals = as.list(w_ec[1L, c("D_clim", "D_edaph", "D_joint", "Schoener_D"), drop = FALSE]),
  sp2_p_d = as.numeric(w_2sp$p_value_D[1L]), sp2_p_bg = as.numeric(w_2sp$p_bg_D_joint[1L]),
  sep_agree = mean(as.character(pred_sep) == sep_te$y), uni_agree = mean(as.character(pred_uni) == uni_te$y)
)

# 2. SQLite phenology ingestion + herbarium rooftop quarantine
old_p_db <- Sys.getenv("PLANTMATICS_DB_PATH", unset = NA_character_); old_k_db <- Sys.getenv("PACKERA_DB_PATH", unset = NA_character_)
Sys.setenv(PLANTMATICS_DB_PATH = db_path, PACKERA_DB_PATH = db_path)
ph_sql <- fetch_step02_phenology_from_sqlite(sprintf("S05_%02d", 1:8), vouchers_path = file.path(tmp_root, "nonexistent.csv"), demo = FALSE)
if (is.na(old_p_db)) Sys.unsetenv("PLANTMATICS_DB_PATH") else Sys.setenv(PLANTMATICS_DB_PATH = old_p_db)
if (is.na(old_k_db)) Sys.unsetenv("PACKERA_DB_PATH") else Sys.setenv(PACKERA_DB_PATH = old_k_db)
lat_8 <- seq(32.0, 42.0, length.out = 8L)
df_ph <- data.frame(
  catalogNumber = c(sprintf("S05_%02d", 1:8), "ROOF_01", "ROOF_02"), plant_individual_id = 0L, species_standardized = "Packera dubia",
  determiner_tier = "Tier_1_Gold", is_type_specimen = FALSE, is_centroid_or_inst_coord = c(rep(FALSE, 8L), TRUE, TRUE),
  decimalLatitude = c(lat_8, 35.9108, 36.0014), decimalLongitude = c(seq(-84.0, -76.0, length.out = 8L), -79.0505, -78.9382),
  doy = c(round(90.0 + 4.0 * (lat_8 - 32.0)), 128.0, 92.0), voucher_phenological_state = "Peak_Anthesis", n_anthesis = 3L, total_capitula = 5L, stringsAsFactors = FALSE
)
ph_res <- compute_phenological_anomalies(df_ph, output_csv = NULL, summary_csv = NULL, demo = FALSE)
res$pheno_sql <- list(sql_n = nrow(ph_sql), sql_states = as.character(ph_sql$voucher_phenological_state), sql_n_anthesis = as.integer(ph_sql$n_anthesis), baseline_n = as.integer(ph_res$baseline$n))

# 3. Weighted circular phenology replaces linear ANOVA/Tukey
df_c <- data.frame(
  catalogNumber = sprintf("CW_%02d", 1:16), plant_individual_id = 0L, species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 8L),
  determiner_tier = "Tier_1_Gold", is_type_specimen = FALSE, is_centroid_or_inst_coord = FALSE, decimalLatitude = rep(35.0, 16L), decimalLongitude = rep(-80.0, 16L),
  doy = c(360, 361, 362, 364, 1, 2, 4, 6, 48, 50, 51, 52, 54, 55, 57, 58), voucher_phenological_state = "Peak_Anthesis", n_anthesis = 3L, total_capitula = 5L, stringsAsFactors = FALSE
)
c_sum <- compute_phenological_anomalies(df_c, output_csv = NULL, summary_csv = NULL, demo = FALSE)$summary
df_sw <- data.frame(
  catalogNumber = sprintf("SW_%02d", 1:12), species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 6L),
  determiner_tier = "Tier_1_Gold", authority_weight = 2.0, spatial_balance_weight = 1.0,
  hopkins_adjusted_doy = c(360, 362, 364, 2, 4, 6, 48, 50, 52, 54, 56, 58),
  delta_doy = c(-8.25, -6.25, -4.25, -1.0, 1.0, 3.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0), is_verified_anthesis = TRUE, stringsAsFactors = FALSE
)
ph_sw_test <- run_circular_multigroup_phase_test(doy_to_rad(df_sw$hopkins_adjusted_doy), df_sw$species_standardized, B = 199L, seed = 42L)
sw_raw <- summarize_weighted_circular_phenology(df_sw, c("Packera aurea", "Packera dubia"), ph_sw_test, is_unif = FALSE, B = 199L, seed = 42L)
sw_sum <- if (is.data.frame(sw_raw)) sw_raw else if (is.list(sw_raw) && "summary" %in% names(sw_raw)) sw_raw$summary else do.call(rbind, Filter(is.data.frame, sw_raw))
sw_aurea_m <- if (is.data.frame(sw_sum) && all(c("Analysis_Section", "Taxon_or_Comparison", "Mean_Delta_DOY") %in% names(sw_sum))) as.numeric(sw_sum$Mean_Delta_DOY[sw_sum$Analysis_Section == "Taxon_Distribution" & sw_sum$Taxon_or_Comparison == "Packera aurea"][1L]) else NA_real_
sw_pw_d <- if (is.data.frame(sw_sum) && all(c("Analysis_Section", "Mean_Delta_DOY") %in% names(sw_sum))) as.numeric(sw_sum$Mean_Delta_DOY[sw_sum$Analysis_Section == "Circular_Pairwise_PostHoc"][1L]) else NA_real_
sw_p <- if (is.data.frame(sw_sum) && "P_Value" %in% names(sw_sum)) suppressWarnings(as.numeric(sw_sum$P_Value[sw_sum$Analysis_Section == "Circular_Pairwise_PostHoc"][1L])) else NA_real_
res$circ <- list(
  sections = as.character(c_sum$Analysis_Section), phase_p = as.numeric(c_sum$Circ_Perm_P[c_sum$Analysis_Section == "Circular_Phase_Test"][1L]),
  aurea_mean = as.numeric(c_sum$Mean_Delta_DOY[c_sum$Analysis_Section == "Taxon_Distribution" & c_sum$Taxon_or_Comparison == "Packera aurea"][1L]),
  pw_diff = as.numeric(c_sum$Mean_Delta_DOY[c_sum$Analysis_Section == "Circular_Pairwise_PostHoc"][1L]),
  sw_rows = if (is.data.frame(sw_sum)) nrow(sw_sum) else 0L, sw_aurea_m = sw_aurea_m, sw_pw_d = sw_pw_d, sw_p = sw_p
)

# 4. Held-out OOF-only Spatial Block CV & SDM metrics
set.seed(77L)
df_cv <- data.frame(
  catalogNumber = sprintf("CV_%02d", 1:13), species_standardized = c(rep("Packera aurea", 12L), "Packera dubia"), analysis_group = c(rep("Packera aurea", 12L), "Packera dubia"),
  determiner_tier = "Tier_1_Gold", is_type_specimen = FALSE, is_confident_label_noise = FALSE, is_centroid_or_inst_coord = FALSE,
  is_training_eligible_coord = TRUE, is_edaphic_eligible_coord = TRUE, coordinateUncertaintyInMeters = 100.0,
  latitude = c(seq(34.0, 38.0, length.out = 12L), 42.5), longitude = c(seq(-84.0, -79.0, length.out = 12L), -72.0),
  bio1_temp_mean = c(rnorm(12L, 14.5, 0.35), 8.0), bio12_precip_annual = c(rnorm(12L, 1320.0, 20.0), 920.0),
  fine_soil_ph = c(rnorm(12L, 5.4, 0.15), 7.4), stringsAsFactors = FALSE
)
srf_cv <- compute_spatial_rf_mems(df_cv, env_vars = c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_ph"), n_trees = 100L)
oof_ok <- srf_cv$rf_prediction_source == "Spatial_Block_OOF_Anchor"
oof_acc <- if (any(oof_ok)) mean(as.character(srf_cv$rf_predicted_taxon[oof_ok]) == df_cv$species_standardized[oof_ok]) else NA_real_

set.seed(88L)
df_sdm <- data.frame(
  catalogNumber = sprintf("SDM_%02d", 1:16), species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 8L),
  determiner_tier = "Tier_1_Gold", is_type_specimen = FALSE, is_confident_label_noise = FALSE, is_centroid_or_inst_coord = FALSE,
  is_training_eligible_coord = TRUE, coordinateUncertaintyInMeters = 100.0,
  longitude = c(seq(-83.5, -79.5, length.out = 8L), seq(-82.5, -77.5, length.out = 8L)), latitude = c(seq(35.0, 38.5, length.out = 8L), seq(33.5, 36.5, length.out = 8L)),
  bio1_temp_mean = c(rnorm(8L, 12.0, 0.5), rnorm(8L, 17.5, 0.5)), bio12_precip_annual = c(rnorm(8L, 1450.0, 35.0), rnorm(8L, 1150.0, 35.0)),
  fine_soil_ph = c(rnorm(8L, 5.8, 0.2), rnorm(8L, 4.9, 0.2)), fine_soil_sand = c(rnorm(8L, 35.0, 3.0), rnorm(8L, 68.0, 3.0)), stringsAsFactors = FALSE
)
sdm_perf <- fit_and_project_maxent_sdm_maps(df = df_sdm, srf_res = NULL, env_dir = file.path(tmp_root, "sdm_env"), out_dir = file.path(tmp_root, "sdm_out"),
                                            kept_vars = c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_ph", "fine_soil_sand"), focal_taxon = "Packera aurea", intraspecific_mode = FALSE, demo = FALSE)$performance
res$cv_sdm <- list(b_src = as.character(srf_cv$rf_prediction_source[13L]), b_pred = as.character(srf_cv$rf_predicted_taxon[13L]), a_oof_n = sum(oof_ok[1:12]),
                   cv_acc = as.numeric(srf_cv$Spatial_Block_CV_Accuracy), oof_acc = as.numeric(oof_acc), sdm_n = nrow(sdm_perf),
                   auc_mx = as.numeric(sdm_perf$AUC_Block_MaxEnt), auc_rf = as.numeric(sdm_perf$AUC_Block_SpatialRF), auc_ens = as.numeric(sdm_perf$AUC_Block_Ensemble))

# 5. Haversine km geometry & 125 km accessible M buffering
D_km <- compute_haversine_dist_km_matrix(rbind(c(-80.0, 45.0), c(-79.0, 45.0), c(-80.0, 46.0)))
prep_g <- list(
  scores_glob = data.frame(PC1 = rep(0.0, 5L), PC2 = rep(0.0, 5L)),
  scores_M_bg = data.frame(PC1 = c(seq(0.1, 2.0, length.out = 20L), 777.0, rep(-999.0, 10L)), PC2 = c(seq(0.1, 2.0, length.out = 20L), 777.0, rep(-999.0, 10L))),
  coords_occ = cbind(longitude = rep(-80.0, 5L), latitude = rep(60.0, 5L)),
  coords_M_bg = cbind(longitude = c(rep(-80.0, 20L), -78.7, rep(-80.0, 10L)), latitude = c(seq(60.05, 60.40, length.out = 20L), 60.0, seq(55.0, 55.9, length.out = 10L)))
)
g1_out <- .resolve_lineage_accessible_glob1(prep_g, 1:5)
res$geo <- list(sym_err = max(abs(D_km - t(D_km))), diag_err = max(abs(diag(D_km))), d_ew = as.numeric(D_km[1L, 2L]), d_ns = as.numeric(D_km[1L, 3L]),
                incl_72km = any(abs(g1_out[, 1L] - 777.0) < 1e-6), excl_far = !any(abs(g1_out[, 1L] + 999.0) < 1e-6))

# 6. BIO2 column fix, YAML config loader, and coastal cache hit idempotency
df_mc <- data.frame(catalogNumber = sprintf("MC_%02d", 1:6), bio2_mean_diurnal_range = c(5.5, 5.5, 5.5, 15.5, 15.5, 15.5),
                    bio4_temp_seasonality = rep(800.0, 6L), bio6_min_temp_coldest_month = rep(2.0, 6L), stringsAsFactors = FALSE)
res_mc <- compute_derived_habitat_predictors(df_mc)
cfg_s05 <- load_step05_macroecology_config("config/config.yaml")
ss_calls <- 0L
mock_ssurgo <- function(d_in, cache_dir = NULL, batch_size = 75, demo = FALSE) {
  ss_calls <<- ss_calls + 1L
  d_in$musym <- "BeB"; d_in$muname <- "Beach sand, coastal"; d_in$taxsuborder <- "Psamments"; d_in$taxgreatgroup <- "Quartzipsamments"
  d_in$is_rock_outcrop <- FALSE; d_in$outcrop_class <- "Non-outcrop matrix"; d_in$fine_soil_ph <- NA_real_; d_in$ssurgo_ph <- NA_real_
  d_in$ilr_texture_1 <- 0.1; d_in$ilr_texture_2 <- -0.1; d_in
}
df_coast <- data.frame(catalogNumber = "COAST_CACHE_01", latitude = 34.2104, longitude = -77.8868, coordinateUncertaintyInMeters = 500.0, stringsAsFactors = FALSE)
c_path <- file.path(tmp_root, "fine_grain_cache_test.csv")
load_and_update_fine_grain_cache(df_coast, env_dir = file.path(tmp_root, "cache_env"), ssurgo_cache_dir = file.path(tmp_root, "cache_ssurgo"), demo = FALSE, query_ssurgo_fn = mock_ssurgo, cache_path = c_path)
calls_1 <- ss_calls
load_and_update_fine_grain_cache(df_coast, env_dir = file.path(tmp_root, "cache_env"), ssurgo_cache_dir = file.path(tmp_root, "cache_ssurgo"), demo = FALSE, query_ssurgo_fn = mock_ssurgo, cache_path = c_path)
calls_2 <- ss_calls - calls_1
res$cfg_cache <- list(
  has_mc = "maritime_continentality_idx" %in% names(res_mc),
  mc_vals = if ("maritime_continentality_idx" %in% names(res_mc)) as.numeric(res_mc$maritime_continentality_idx) else numeric(0),
  cfg = cfg_s05, calls_1 = as.integer(calls_1), calls_2 = as.integer(calls_2)
)

jsonlite::write_json(res, out_json, auto_unbox = TRUE, na = "null", digits = 10)
"""


class TestStep05AuditHygiene(unittest.TestCase):
    """Contract and hygiene audit suite for the refactored Step 05 macroecology stack."""

    results: dict = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            cls.results = {}
            return
        cls._tmp_dir = tempfile.TemporaryDirectory(prefix="step05_hyg_")
        tmp_p = Path(cls._tmp_dir.name)
        db_path = _seed_step05_db(tmp_p)
        r_path, out_json = tmp_p / "run_step05_hygiene.R", tmp_p / "results.json"
        r_path.write_text(_WARM_R_SCRIPT, encoding="utf-8")
        proc = subprocess.run(
            ["Rscript", str(r_path), str(out_json), str(tmp_p), str(db_path)],
            cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False, timeout=180,
        )
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(f"Batched Step 05 R failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls.results = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_post_selection_ecotype_quarantine_and_warren_pvalue_guard(self) -> None:
        ec = self.results["ecotype"]
        self.assertEqual(int(ec["n_groups"]), 2)
        for flag in ("all_ecotype_prefix", "grp1_uniform", "grp2_uniform", "grps_distinct", "quarantined_passive", "attr_post_sel"):
            self.assertTrue(bool(ec[flag]), f"Ecotype quarantine check failed: {flag}")
        for pk, pv in ec["ec_pvals"].items():
            self.assertIsNone(pv, f"Expected NA/None on post-selection {pk}, got {pv}")
        for dk, dv in ec["ec_dvals"].items():
            self.assertTrue(_is_finite(dv) and 0.0 <= float(dv) <= 1.0, f"Invalid {dk}: {dv}")
        self.assertTrue(_is_finite(ec["sp2_p_d"]) and 0.0 < float(ec["sp2_p_d"]) <= 1.0)
        self.assertTrue(_is_finite(ec["sp2_p_bg"]) and 0.0 < float(ec["sp2_p_bg"]) <= 1.0)
        self.assertGreaterEqual(float(ec["sep_agree"]), 0.80)
        self.assertLess(float(ec["uni_agree"]), 0.72)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_sqlite_phenology_ingestion_and_rooftop_quarantine(self) -> None:
        ps = self.results["pheno_sql"]
        self.assertEqual(int(ps["sql_n"]), 8)
        self.assertEqual(ps["sql_states"], ["Peak_Anthesis"] * 8)
        self.assertEqual([int(x) for x in ps["sql_n_anthesis"]], [3] * 8)
        self.assertEqual(int(ps["baseline_n"]), 8)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_weighted_circular_phenology_replaces_linear_anova_tukey(self) -> None:
        cr = self.results["circ"]
        self.assertIn("Circular_Phase_Test", cr["sections"])
        self.assertIn("Circular_Pairwise_PostHoc", cr["sections"])
        self.assertNotIn("One_Way_ANOVA", cr["sections"])
        self.assertNotIn("Tukey_HSD_PostHoc", cr["sections"])
        self.assertTrue(_is_finite(cr["phase_p"]) and float(cr["phase_p"]) < 0.01)
        self.assertTrue(_is_finite(cr["aurea_mean"]) and abs(float(cr["aurea_mean"])) < 40.0)
        self.assertTrue(_is_finite(cr["pw_diff"]) and 40.0 <= abs(float(cr["pw_diff"])) <= 65.0)
        self.assertGreaterEqual(int(cr["sw_rows"]), 2)
        if _is_finite(cr["sw_aurea_m"]):
            self.assertLessEqual(abs(float(cr["sw_aurea_m"])), 10.0)
        if _is_finite(cr["sw_pw_d"]):
            self.assertTrue(40.0 <= abs(float(cr["sw_pw_d"])) <= 65.0)
        if _is_finite(cr["sw_p"]):
            self.assertLess(float(cr["sw_p"]), 0.01)
        pheno_txt = (PROJECT_ROOT / "scripts" / "analysis" / "step05_phenology_ecotypes.R").read_text(encoding="utf-8")
        self.assertNotIn("TukeyHSD", pheno_txt)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_heldout_oof_only_cv_and_sdm_metrics(self) -> None:
        cs = self.results["cv_sdm"]
        self.assertEqual(cs["b_src"], "Fallback_OOB_Anchor")
        self.assertNotEqual(cs["b_pred"], "Packera dubia")
        self.assertGreaterEqual(int(cs["a_oof_n"]), 1)
        self.assertTrue(_is_finite(cs["cv_acc"]) and _is_finite(cs["oof_acc"]))
        self.assertAlmostEqual(float(cs["cv_acc"]), float(cs["oof_acc"]), places=5)
        self.assertGreaterEqual(int(cs["sdm_n"]), 1)
        for ak in ("auc_mx", "auc_rf", "auc_ens"):
            vals = cs[ak] if isinstance(cs[ak], list) else [cs[ak]]
            self.assertTrue(len(vals) >= 1 and all(_is_finite(v) and 0.0 <= float(v) <= 1.0 for v in vals), f"Invalid {ak}: {vals}")
        sdm_txt = (PROJECT_ROOT / "scripts" / "analysis" / "step05_sdm_maxent.R").read_text(encoding="utf-8")
        sdm_fit = sdm_txt[sdm_txt.index("fit_and_project_maxent_sdm_maps <-"):]
        for tok in ("ok_oof_mx", "ok_oof_rf", "ok_oof_ens"):
            self.assertGreaterEqual(sdm_fit.count(tok), 2, f"Missing {tok} gating in step05_sdm_maxent.R")

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_haversine_km_geometry_in_mems_moran_and_accessible_m(self) -> None:
        go = self.results["geo"]
        self.assertLess(float(go["sym_err"]), 1e-10)
        self.assertLess(float(go["diag_err"]), 1e-10)
        self.assertTrue(77.0 <= float(go["d_ew"]) <= 80.0, f"Unexpected E-W 1 deg km at 45N: {go['d_ew']}")
        self.assertTrue(110.0 <= float(go["d_ns"]) <= 112.5, f"Unexpected N-S 1 deg km: {go['d_ns']}")
        self.assertTrue(bool(go["incl_72km"]))
        self.assertTrue(bool(go["excl_far"]))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_bio2_column_fix_yaml_config_and_cache_hit_idempotency(self) -> None:
        cc = self.results["cfg_cache"]
        self.assertTrue(bool(cc["has_mc"]))
        vals = [float(v) for v in cc["mc_vals"]]
        self.assertEqual(len(vals), 6)
        self.assertTrue(all(_is_finite(v) for v in vals))
        self.assertGreater(abs(vals[0] - vals[3]), 1.0)
        cfg = cc["cfg"]
        self.assertAlmostEqual(float(cfg["min_flowering_doy"]), 30.0, places=6)
        self.assertAlmostEqual(float(cfg["cor_threshold"]), 0.70, places=6)
        self.assertAlmostEqual(float(cfg["vif_threshold"]), 5.0, places=6)
        self.assertAlmostEqual(float(cfg["cumulative_importance_cutoff"]), 0.85, places=6)
        self.assertEqual((int(cc["calls_1"]), int(cc["calls_2"])), (1, 0))

    def test_step05_figure_hygiene_and_line_ceilings(self) -> None:
        ceilings = {
            PROJECT_ROOT / "scripts" / "figures" / "step05_figures.R": 500,
            PROJECT_ROOT / "scripts" / "figures" / "step05_sdm_figures.R": 320,
            PROJECT_ROOT / "scripts" / "analysis" / "step05_config_geo_utils.R": 320,
            Path(__file__).resolve(): 380,
        }
        for path, cap in ceilings.items():
            self.assertTrue(path.is_file(), f"Missing expected Step 05 file: {path}")
            n_lines = len(path.read_text(encoding="utf-8").splitlines())
            self.assertLess(n_lines, cap, f"{path.name}: {n_lines} lines >= {cap}")
        fig_txt = (PROJECT_ROOT / "scripts" / "figures" / "step05_figures.R").read_text(encoding="utf-8")
        self.assertNotIn("unlink(", fig_txt)
        self.assertNotIn("POLARIS", fig_txt)
        sdm_fig_txt = (PROJECT_ROOT / "scripts" / "figures" / "step05_sdm_figures.R").read_text(encoding="utf-8")
        self.assertNotIn("unlink(", sdm_fig_txt)
        rf_txt = (PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R").read_text(encoding="utf-8")
        m_clim = re.search(r"clim_cand_cv\s*<-.*?(?=\n\s*[A-Za-z_][A-Za-z0-9_]*\s*<-|\n\s*\})", rf_txt, flags=re.DOTALL)
        self.assertIsNotNone(m_clim, "Missing clim_cand_cv definition in 05_multimodal_spatial_rf.R")
        self.assertNotIn("sin_doy", m_clim.group(0))
        self.assertNotIn("cos_doy", m_clim.group(0))


if __name__ == "__main__":
    unittest.main()
