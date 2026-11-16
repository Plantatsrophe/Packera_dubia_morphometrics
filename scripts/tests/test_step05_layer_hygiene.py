"""Warm-R and Python contract tests for Step 05 layer hygiene and MEM decoupling."""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.data_prep import download_geochem_wetland_layers as dl_geochem

SKIP_R_REASON = "Skipping batched R Step 05 layer hygiene tests under fast audit or missing Rscript."
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

_WARM_R_SCRIPT = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]
tmp_root <- args[2L]

suppressPackageStartupMessages({
  library(jsonlite)
  library(terra)
  if (requireNamespace("ranger", quietly = TRUE)) library(ranger)
})

source("scripts/core/db_utils.R")
source("scripts/morphometrics/authority_weighting_utils.R")
source("scripts/morphometrics/pub_figures_utils.R")
source("scripts/analysis/step05_layer_hygiene_utils.R")
source("scripts/analysis/edaphic_geochem_wetland_utils.R")
source("scripts/analysis/edaphic_fine_grain_utils.R")
source("scripts/analysis/step05_phenology_ecotypes.R")
source("scripts/analysis/05_multimodal_spatial_rf.R")

results <- list()

# 1. DS-801 ultramafic_calcareous_gradient stays in audit df, excluded from X_env & RF
set.seed(101L)
n1 <- 20L
sp1 <- rep(c("Packera serpenticola", "Packera aurea"), each = 10L)
df_x <- data.frame(
  catalogNumber = sprintf("DS_%02d", seq_len(n1)),
  species_standardized = sp1,
  determiner_tier = "Tier_1_Gold",
  latitude = c(seq(35.0, 35.4, length.out = 10L), seq(36.0, 36.4, length.out = 10L)),
  longitude = c(seq(-83.8, -83.4, length.out = 10L), seq(-81.6, -81.2, length.out = 10L)),
  coordinateUncertaintyInMeters = 150.0,
  fine_soil_ph = c(rnorm(10L, 6.4, 0.2), rnorm(10L, 5.4, 0.2)),
  ilr_texture_1 = rnorm(n1, 0.1, 0.3),
  ilr_texture_2 = rnorm(n1, -0.1, 0.2),
  topo_slope_deg = rnorm(n1, 8.0, 2.0),
  bio1_temp_mean = c(rnorm(10L, 14.0, 0.4), rnorm(10L, 12.2, 0.4)),
  bio12_precip_annual = rnorm(n1, 1350.0, 45.0),
  soil_ca_mg_molar_ratio = c(rnorm(10L, 0.12, 0.02), rnorm(10L, 3.6, 0.3)),
  soil_ultramafic_ni_cr_index = c(rnorm(10L, 2.5, 0.2), rnorm(10L, 0.2, 0.05)),
  soil_carbonate_c_pct = c(rnorm(10L, 0.04, 0.01), rnorm(10L, 2.1, 0.2)),
  stringsAsFactors = FALSE
)
df_comp <- compute_derived_habitat_predictors(df_x)
cand_ds <- discover_candidate_eco_predictors(df_comp, include_pheno = FALSE)
srf_ds <- compute_spatial_rf_mems(df_comp, c(cand_ds, "ultramafic_calcareous_gradient"), n_trees = 60)
df_intra_ds <- df_comp
df_intra_ds$species_standardized <- "Packera dubia"
rag_ds <- resolve_analysis_group(df_intra_ds, intraspecific_mode = TRUE, primary_env_vars = c(cand_ds, "ultramafic_calcareous_gradient"))
results$ds801 <- list(
  ucg_finite_count = sum(is.finite(df_comp$ultramafic_calcareous_gradient)),
  ucg_in_cand = "ultramafic_calcareous_gradient" %in% cand_ds,
  ucg_in_kept = "ultramafic_calcareous_gradient" %in% srf_ds$collinearity_screening$kept_vars,
  ucg_in_imp = "ultramafic_calcareous_gradient" %in% srf_ds$importance$Predictor,
  rag_ok = nrow(rag_ds) == n1
)

# 2. Bedrock buffer radius capped at 300 m while general radii keep [100, 1000] m
radii_chk <- as.numeric(compute_uncertainty_buffer_radius(c(NA, 25, 250, 600, 3500)))
bd_env1 <- file.path(tmp_root, "bd_env1")
bd_env2 <- file.path(tmp_root, "bd_env2")
dir.create(file.path(bd_env1, "solus_100m"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(bd_env2, "solus_100m"), recursive = TRUE, showWarnings = FALSE)
r_bd1 <- terra::rast(xmin = -80.10, xmax = -79.90, ymin = 35.90, ymax = 36.10, nrows = 600, ncols = 600, crs = "EPSG:4326")
terra::values(r_bd1) <- 180.0
r_bd2 <- r_bd1
c_near <- terra::cellFromXY(r_bd1, cbind(-80.0010, 36.0000))
c_far <- terra::cellFromXY(r_bd1, cbind(-80.0061, 36.0000))
r_bd1[c_near] <- 15.0
r_bd1[c_far] <- 10.0
r_bd2[c_far] <- 10.0
terra::writeRaster(r_bd1, file.path(bd_env1, "solus_100m", "anylithicdpt_cm_p.tif"), overwrite = TRUE)
terra::writeRaster(r_bd2, file.path(bd_env2, "solus_100m", "anylithicdpt_cm_p.tif"), overwrite = TRUE)
df_bd_in <- data.frame(catalogNumber = "BD_01", latitude = 36.0000, longitude = -80.0000, coordinateUncertaintyInMeters = 900.0, stringsAsFactors = FALSE)
df_bd_out1 <- extract_uncertainty_buffered_predictors(df_bd_in, env_dir = bd_env1, demo = FALSE)
df_bd_out2 <- extract_uncertainty_buffered_predictors(df_bd_in, env_dir = bd_env2, demo = FALSE)
results$bedrock <- list(
  radii = radii_chk,
  buffer_radius_m = as.numeric(df_bd_out1$buffer_radius_m[1L]),
  bd_with_90m_and_550m = as.numeric(df_bd_out1$min_bedrock_depth_buffer_cm[1L]),
  bd_with_550m_only = as.numeric(df_bd_out2$min_bedrock_depth_buffer_cm[1L])
)

# 3. 30 m microtopography and NWI seepage gated at <= 250 m without cache refresh loops
micro_env <- file.path(tmp_root, "micro_env")
dir.create(file.path(micro_env, "dem_3dep_30m"), recursive = TRUE, showWarnings = FALSE)
r_dem <- terra::rast(xmin = -80.10, xmax = -79.90, ymin = 35.90, ymax = 36.10, nrows = 300, ncols = 300, crs = "EPSG:4326")
xy_dem <- terra::xyFromCell(r_bd <- r_dem, seq_len(terra::ncell(r_dem)))
terra::values(r_dem) <- 200.0 + 35.0 * sin((xy_dem[, 1] + 80.0) * 60.0) + 25.0 * cos((xy_dem[, 2] - 36.0) * 60.0)
terra::writeRaster(r_dem, file.path(micro_env, "dem_3dep_30m", "test_3dep_tile.tif"), overwrite = TRUE)

df_m_in <- data.frame(
  catalogNumber = sprintf("MIC_%d", 1:5),
  latitude = c(35.98, 35.99, 36.00, 36.01, 36.02),
  longitude = c(-80.02, -80.01, -80.00, -79.99, -79.98),
  coordinateUncertaintyInMeters = c(120.0, 250.0, 260.0, 800.0, NA_real_),
  soil_ph = c(5.2, 5.4, 5.6, 5.8, 6.0),
  soil_sand = c(40.0, 42.0, 44.0, 46.0, 48.0),
  soil_clay = c(20.0, 22.0, 24.0, 26.0, 28.0),
  soil_silt = c(40.0, 36.0, 32.0, 28.0, 24.0),
  is_nwi_wetland = rep(1L, 5L),
  nwi_wetland_system = rep("PEM1A", 5L),
  nwi_wetland_dist_m = c(5.0, 10.0, 8.0, 12.0, 15.0),
  stringsAsFactors = FALSE
)
df_topo_direct <- extract_micro_topography_3dep(df_m_in, env_dir = micro_env, demo = FALSE)
df_wsi_direct <- enrich_geochem_wetland_predictors(df_topo_direct, env_dir = micro_env)
cache_csv <- file.path(micro_env, "fine_grain_edaphic_cache.csv")
df_c1 <- load_and_update_fine_grain_cache(df_m_in, env_dir = micro_env, demo = FALSE, cache_path = cache_csv)
df_c1 <- enrich_geochem_wetland_predictors(df_c1, env_dir = micro_env)
msg_c2 <- capture.output({
  df_c2 <- load_and_update_fine_grain_cache(df_m_in, env_dir = micro_env, demo = FALSE, cache_path = cache_csv)
}, type = "message")
results$micro_gate <- list(
  direct_tpi_fin = is.finite(df_topo_direct$topo_tpi_30m),
  direct_hli_fin = is.finite(df_topo_direct$topo_hli_30m),
  direct_twi_fin = is.finite(df_topo_direct$topo_twi_30m),
  direct_elev_fin = is.finite(df_topo_direct$elevation_m),
  direct_slope_fin = is.finite(df_topo_direct$topo_slope_deg),
  wsi_fin = is.finite(df_wsi_direct$wetland_seepage_index),
  nwi_flag = as.integer(df_wsi_direct$is_nwi_wetland),
  cache_tpi_fin = is.finite(df_c1$topo_tpi_30m),
  cache_c2_tpi_fin = is.finite(df_c2$topo_tpi_30m),
  cache_c2_zero_refresh = any(grepl("0 coordinates requiring refresh", msg_c2, fixed = TRUE))
)

# 4. Partial 30 m POLARIS raster override disabled
pol_env <- file.path(tmp_root, "pol_env")
dir.create(file.path(pol_env, "polaris_30m"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(pol_env, "solus_100m"), recursive = TRUE, showWarnings = FALSE)
r_pol <- terra::rast(xmin = -80.10, xmax = -79.90, ymin = 35.90, ymax = 36.10, nrows = 100, ncols = 100, crs = "EPSG:4326")
terra::values(r_pol) <- 8.4
r_sol <- r_pol
terra::values(r_sol) <- 510.0
terra::writeRaster(r_pol, file.path(pol_env, "polaris_30m", "ph_0_5.tif"), overwrite = TRUE)
terra::writeRaster(r_sol, file.path(pol_env, "solus_100m", "ph1to1h2o_0_cm_p.tif"), overwrite = TRUE)
df_pol_in <- data.frame(
  catalogNumber = "POL_01", latitude = 36.0, longitude = -80.0,
  coordinateUncertaintyInMeters = 100.0, polaris_ph_30m = 8.4, stringsAsFactors = FALSE
)
df_pol_out <- extract_polaris_and_solus_predictors(df_pol_in, env_dir = pol_env, demo = FALSE)
results$polaris <- list(
  polaris_ph_is_na = is.na(df_pol_out$polaris_ph_30m[1L]),
  fine_soil_ph = as.numeric(df_pol_out$fine_soil_ph[1L])
)

# 5. MEMs decoupled from rf_predicted_taxon, rf_predicted_prob, OOB_Accuracy, and importance
set.seed(42L)
n_cl <- 20L
n_tot <- 2L * n_cl
sp_m <- rep(c("Packera dubia", "Packera anonyma"), each = n_cl)
df_mem_a <- data.frame(
  catalogNumber = sprintf("MEM_%02d", seq_len(n_tot)),
  species_standardized = sp_m,
  analysis_group = sp_m,
  determiner_tier = "Tier_1_Gold",
  spatial_balance_weight = rep(1.0, n_tot),
  latitude = c(seq(34.0, 34.4, length.out = n_cl), seq(39.0, 39.4, length.out = n_cl)),
  longitude = c(seq(-81.0, -80.6, length.out = n_cl), seq(-80.9, -80.5, length.out = n_cl)),
  coordinateUncertaintyInMeters = 150.0,
  fine_soil_ph = rnorm(n_tot, 5.8, 0.5),
  ilr_texture_1 = rnorm(n_tot, 0.0, 0.6),
  bio1_temp_mean = rnorm(n_tot, 14.0, 1.0),
  bio12_precip_annual = rnorm(n_tot, 1250.0, 90.0),
  stringsAsFactors = FALSE
)
df_mem_b <- df_mem_a
idx_b <- sp_m == "Packera anonyma"
df_mem_b$latitude[idx_b] <- df_mem_b$latitude[idx_b] + 8.0
df_mem_b$longitude[idx_b] <- df_mem_b$longitude[idx_b] + 6.0
idx_a <- which(!idx_b)
perm_a <- rev(idx_a)
df_mem_b$latitude[idx_a] <- df_mem_b$latitude[perm_a]
df_mem_b$longitude[idx_a] <- df_mem_b$longitude[perm_a]
env_m <- c("fine_soil_ph", "ilr_texture_1", "bio1_temp_mean", "bio12_precip_annual")
srf_a <- compute_spatial_rf_mems(df_mem_a, env_m, n_trees = 100)
folds_anc_fixed <- as.integer(srf_a$spatial_block_id[srf_a$anchor_idx])
orig_cv_sp <- getFromNamespace("cv_spatial", "blockCV")
srf_b <- tryCatch({
  utils::assignInNamespace("cv_spatial", function(...) list(folds_ids = folds_anc_fixed), ns = "blockCV")
  compute_spatial_rf_mems(df_mem_b, env_m, n_trees = 100)
}, finally = {
  utils::assignInNamespace("cv_spatial", orig_cv_sp, ns = "blockCV")
})
imp_ns <- ranger::importance(srf_a$rf_nonspatial)
env_rows <- srf_a$importance[srf_a$importance$Predictor %in% names(imp_ns), , drop = FALSE]
max_imp_diff <- max(abs(env_rows$MeanDecreaseAccuracy - as.numeric(imp_ns[env_rows$Predictor])))
results$mem_decouple <- list(
  same_folds = identical(srf_a$spatial_block_id, srf_b$spatial_block_id),
  mems_changed = max(abs(srf_a$mems - srf_b$mems)) > 1e-6,
  same_pred_taxon = identical(srf_a$rf_predicted_taxon, srf_b$rf_predicted_taxon),
  max_prob_delta = max(abs(srf_a$rf_predicted_prob - srf_b$rf_predicted_prob)),
  oob_equals_nonspatial = isTRUE(all.equal(srf_a$OOB_Accuracy, srf_a$NonSpatial_OOB_Accuracy, tolerance = 1e-12)),
  oob_nonspatial = as.numeric(srf_a$NonSpatial_OOB_Accuracy),
  oob_spatial_diag = as.numeric(srf_a$Spatial_MEM_OOB_Accuracy),
  max_env_imp_diff_vs_nonspatial = as.numeric(max_imp_diff)
)

writeLines(toJSON(results, auto_unbox = TRUE, null = "null"), out_json)
cat("STEP05_LAYER_HYGIENE_R_OK\n")
"""


@unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
class TestStep05LayerHygiene(unittest.TestCase):
    """Verify Step 05 DS-801 exclusion, 300 m bedrock cap, 250 m gate, POLARIS retirement, and MEM decoupling."""

    results: dict = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            cls.results = {}
            return
        with tempfile.TemporaryDirectory(prefix="step05_hygiene_") as tmp_dir:
            tmp_path = Path(tmp_dir)
            r_path = tmp_path / "run_hygiene.R"
            out_json = tmp_path / "results.json"
            r_path.write_text(_WARM_R_SCRIPT, encoding="utf-8")
            proc = subprocess.run(
                ["Rscript", str(r_path), str(out_json), str(tmp_path)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
                timeout=180,
            )
            if proc.returncode != 0 or not out_json.exists():
                raise RuntimeError(f"Warm R hygiene check failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.results = json.loads(out_json.read_text(encoding="utf-8"))

    def test_ds801_10km_cap_and_x_env_exclusion(self) -> None:
        """DS-801 IDW caps at 10.0 km by default and ultramafic_calcareous_gradient stays audit-only."""
        sites = [{"latitude": 35.0, "longitude": -83.62, "ca_pct": 1.0, "mg_pct": 14.0, "ni_ppm": 1800.0, "cr_ppm": 2200.0, "c_inorg_pct": 0.02}]
        far_res = dl_geochem.compute_idw_soil_geochem(35.225, -83.62, sites)
        near_res = dl_geochem.compute_idw_soil_geochem(35.036, -83.62, sites)
        self.assertIsNone(far_res["soil_ca_mg_molar_ratio"])
        self.assertIsNone(far_res["soil_ultramafic_ni_cr_index"])
        self.assertIsNone(far_res["soil_carbonate_c_pct"])
        self.assertIsNotNone(near_res["soil_ca_mg_molar_ratio"])
        self.assertTrue(math.isfinite(float(near_res["soil_ca_mg_molar_ratio"])))

        ds = self.results["ds801"]
        self.assertEqual(int(ds["ucg_finite_count"]), 20)
        self.assertFalse(bool(ds["ucg_in_cand"]))
        self.assertFalse(bool(ds["ucg_in_kept"]))
        self.assertFalse(bool(ds["ucg_in_imp"]))
        self.assertTrue(bool(ds["rag_ok"]))

    def test_min_bedrock_buffer_300m_radius_cap(self) -> None:
        """Bedrock min depth buffer caps at 300 m while compute_uncertainty_buffer_radius retains [100, 1000] m."""
        bd = self.results["bedrock"]
        self.assertEqual([float(x) for x in bd["radii"]], [250.0, 100.0, 250.0, 600.0, 1000.0])
        self.assertEqual(float(bd["buffer_radius_m"]), 900.0)
        self.assertAlmostEqual(float(bd["bd_with_90m_and_550m"]), 15.0, places=1)
        self.assertAlmostEqual(float(bd["bd_with_550m_only"]), 180.0, places=1)

    def test_micro_topography_and_nwi_decay_250m_uncertainty_gate(self) -> None:
        """30 m TPI/HLI/TWI and 60 m NWI decay require <= 250 m without triggering cache refresh loops."""
        mg = self.results["micro_gate"]
        self.assertEqual(mg["direct_tpi_fin"], [True, True, False, False, False])
        self.assertEqual(mg["direct_hli_fin"], [True, True, False, False, False])
        self.assertEqual(mg["direct_twi_fin"], [True, True, False, False, False])
        self.assertEqual(mg["wsi_fin"], [True, True, False, False, False])
        self.assertEqual(mg["direct_elev_fin"], [True, True, True, True, True])
        self.assertEqual(mg["direct_slope_fin"], [True, True, True, True, True])
        self.assertEqual(mg["nwi_flag"], [1, 1, 1, 1, 1])
        self.assertEqual(mg["cache_tpi_fin"], [True, True, False, False, False])
        self.assertEqual(mg["cache_c2_tpi_fin"], [True, True, False, False, False])
        self.assertTrue(bool(mg["cache_c2_zero_refresh"]))

    def test_partial_polaris_override_disabled(self) -> None:
        """Local partial polaris_30m tiles never override SOLUS100 and polaris_ph_30m stays NA."""
        pol = self.results["polaris"]
        self.assertTrue(bool(pol["polaris_ph_is_na"]))
        self.assertAlmostEqual(float(pol["fine_soil_ph"]), 5.1, places=2)

    def test_mems_decoupled_from_rf_predicted_taxon_and_env_importance(self) -> None:
        """Shifting coordinates under identical fold splits leaves RF predictions invariant and importance non-spatial."""
        md = self.results["mem_decouple"]
        self.assertTrue(bool(md["same_folds"]))
        self.assertTrue(bool(md["mems_changed"]))
        self.assertTrue(bool(md["same_pred_taxon"]))
        self.assertEqual(float(md["max_prob_delta"]), 0.0)
        self.assertTrue(bool(md["oob_equals_nonspatial"]))
        self.assertLess(float(md["max_env_imp_diff_vs_nonspatial"]), 1e-12)
        self.assertGreater(float(md["oob_spatial_diag"]), float(md["oob_nonspatial"]) + 0.25)


if __name__ == "__main__":
    unittest.main()

