#!/usr/bin/env python3
"""Unit and warm-R parity suite for multi-stage thread scaling (--threads / PLANTMATICS_MAX_WORKERS)."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_REASON = "Skipping warm-R thread-scaling parity tests under fast audit or missing Rscript."
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


_WARM_R_PARITY_SCRIPT = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]
tmp_root <- args[2L]
Sys.setenv(PLANTMATICS_SKIP_FIGURES = "1", PACKERA_SKIP_FIGURES = "1")

for (f in c(
  "scripts/core/db_utils.R",
  "scripts/morphometrics/fourier_parallel_utils.R",
  "scripts/morphometrics/step03_config_scale_utils.R",
  "scripts/morphometrics/gmm_icl_bootstrap_utils.R",
  "scripts/analysis/step05_circular_hopkins_models.R",
  "scripts/analysis/step05_config_geo_utils.R",
  "scripts/analysis/05_multimodal_spatial_rf.R",
  "scripts/analysis/step05_sdm_maxent.R"
)) {
  source(f)
}

res <- list()

# 1. Worker cap lift to 58 and nested fork guard inside packera_mclapply
Sys.setenv(PLANTMATICS_MAX_WORKERS = "58", PACKERA_MAX_WORKERS = "58")
Sys.unsetenv("PLANTMATICS_IN_FORK_WORKER")
w58 <- resolve_packera_workers(64L)
in_fork_vals <- unlist(packera_mclapply(1:4, function(i) resolve_packera_workers(64L), mc.cores = 2L))
w58_after <- resolve_packera_workers(64L)
res$workers <- list(w58 = w58, in_fork_vals = as.integer(in_fork_vals), w58_after = w58_after)

# 2. Circular phenology permutation + pairwise summary bit parity (1 vs 8 workers)
set.seed(11L)
act_tx <- c("Packera aurea", "Packera dubia", "Packera anonyma")
th_vec <- c(stats::rnorm(12L, 1.8, 0.18), stats::rnorm(12L, 2.7, 0.18), stats::rnorm(12L, 2.2, 0.18))
grp_vec <- rep(act_tx, each = 12L)
df_ph <- data.frame(
  catalogNumber = sprintf("PH_%02d", seq_along(th_vec)),
  species_standardized = grp_vec,
  is_active = TRUE,
  hopkins_adjusted_doy = round(th_vec * 365.25 / (2.0 * pi), 3),
  delta_doy = round((th_vec - mean(th_vec)) * 365.25 / (2.0 * pi), 3),
  authority_weight = 2.0,
  spatial_balance_weight = 1.0,
  stringsAsFactors = FALSE
)
Sys.setenv(PLANTMATICS_MAX_WORKERS = "1", PACKERA_MAX_WORKERS = "1")
c1 <- run_circular_multigroup_phase_test(th_vec, grp_vec, B = 199L, seed = 42L)
s1 <- summarize_weighted_circular_phenology(df_ph, act_tx, c1, B = 199L, seed = 42L)
Sys.setenv(PLANTMATICS_MAX_WORKERS = "8", PACKERA_MAX_WORKERS = "8")
c8 <- run_circular_multigroup_phase_test(th_vec, grp_vec, B = 199L, seed = 42L)
s8 <- summarize_weighted_circular_phenology(df_ph, act_tx, c8, B = 199L, seed = 42L)
res$circ <- list(
  wf_diff = abs(c1$watson_f - c8$watson_f),
  pp_diff = abs(c1$perm_p - c8$perm_p),
  sum_identical = identical(s1, s8)
)

# 3. Step 03 EFA parallel & Step 04 ICL + bootstrap LRT parity (1 vs 4 workers)
t_grid <- seq(0, 2 * pi, length.out = 65L)[-65L]
coo_list <- lapply(seq_len(28L), function(idx) {
  rx <- if (idx <= 14L) 20 + (idx %% 3L) else 12 + (idx %% 3L)
  ry <- if (idx <= 14L) 10 else 22
  cbind(x = -rx * sin(t_grid), y = ry * cos(t_grid))
})
names(coo_list) <- sprintf("LEAF_%02d", seq_len(28L))
fac_df <- data.frame(shape_id = names(coo_list), is_calibrated_mm = TRUE, stringsAsFactors = FALSE)
cal_stub <- function(obj, nb_h = 8L) 8L
Sys.setenv(PLANTMATICS_MAX_WORKERS = "1", PACKERA_MAX_WORKERS = "1")
ef1 <- run_parallel_momocs_efa(coo_list, fac_df, nb_h = 8L, calibrate_fn = cal_stub)
Sys.setenv(PLANTMATICS_MAX_WORKERS = "4", PACKERA_MAX_WORKERS = "4")
ef4 <- run_parallel_momocs_efa(coo_list, fac_df, nb_h = 8L, calibrate_fn = cal_stub)
res$efa_max_diff <- max(abs(ef1$harm_mat - ef4$harm_mat))

set.seed(42L)
X_gmm <- rbind(matrix(rnorm(30L * 3L, -2.0, 0.7), ncol = 3L), matrix(rnorm(30L * 3L, 2.0, 0.7), ncol = 3L))
bt <- data.frame(K = 1:3, BIC = c(-500, -450, -470))
Sys.setenv(PLANTMATICS_MAX_WORKERS = "1", PACKERA_MAX_WORKERS = "1")
icl1 <- compute_gmm_icl_profile(X_gmm, bt, seed = 42L)
lrt1 <- run_gmm_pair_bootstrap_lrt(X_gmm, g_null = 1L, g_alt = 2L, B = 29L, seed = 42L)
Sys.setenv(PLANTMATICS_MAX_WORKERS = "4", PACKERA_MAX_WORKERS = "4")
icl4 <- compute_gmm_icl_profile(X_gmm, bt, seed = 42L)
lrt4 <- run_gmm_pair_bootstrap_lrt(X_gmm, g_null = 1L, g_alt = 2L, B = 29L, seed = 42L)
res$gmm <- list(icl_diff = max(abs(icl1$ICL - icl4$ICL)), lrt_diff = abs(lrt1$lrt_stat - lrt4$lrt_stat), p_diff = abs(lrt1$p_value - lrt4$p_value))

# 4. Step 05 SDM CV + terra raster projection safety across fork (1 vs 4 workers)
if (requireNamespace("terra", quietly = TRUE)) {
  set.seed(77L)
  n_s <- 36L
  df_sdm <- data.frame(
    species_standardized = rep(c("Packera aurea", "Packera dubia"), each = n_s / 2L),
    longitude = c(runif(n_s / 2L, -83, -78), runif(n_s / 2L, -81, -76)),
    latitude = c(runif(n_s / 2L, 35, 39), runif(n_s / 2L, 33, 37)),
    bio1_mean_temp = c(rnorm(n_s / 2L, 12, 1), rnorm(n_s / 2L, 16, 1)),
    bio12_annual_precip = c(rnorm(n_s / 2L, 1200, 50), rnorm(n_s / 2L, 1400, 50)),
    authority_weight = 2.0,
    sdm_prior_weight_01 = 1.0,
    coord_uncertainty_m = 50.0,
    determiner_tier = "Tier_1_Gold",
    stringsAsFactors = FALSE
  )
  r1 <- terra::rast(nrows = 12L, ncols = 12L, xmin = -84, xmax = -75, ymin = 32, ymax = 40, vals = seq(10, 18, length.out = 144L))
  r2 <- terra::rast(nrows = 12L, ncols = 12L, xmin = -84, xmax = -75, ymin = 32, ymax = 40, vals = seq(1100, 1500, length.out = 144L))
  r_stk <- c(r1, r2); names(r_stk) <- c("bio1_mean_temp", "bio12_annual_precip")
  .build_regional_sdm_raster_stack <<- function(...) r_stk
  .ensure_step05_fig_modules <<- function() invisible(NULL)
  render_step05c_sdm_suitability_png <<- function(...) invisible(NULL)
  kv <- c("bio1_mean_temp", "bio12_annual_precip")
  srf_stub <- list(collinearity_screening = list(kept_vars = kv), importance = data.frame(Predictor = kv, MeanDecreaseAccuracy = c(0.5, 0.4)), spatial_autocorr_range_km = 55.0)
  Sys.setenv(PLANTMATICS_MAX_WORKERS = "1", PACKERA_MAX_WORKERS = "1")
  m1 <- fit_and_project_maxent_sdm_maps(df_sdm, srf_res = srf_stub, env_dir = tmp_root, out_dir = file.path(tmp_root, "sdm1"), kept_vars = kv, focal_taxon = "Packera aurea", demo = FALSE)
  Sys.setenv(PLANTMATICS_MAX_WORKERS = "4", PACKERA_MAX_WORKERS = "4")
  m4 <- fit_and_project_maxent_sdm_maps(df_sdm, srf_res = srf_stub, env_dir = tmp_root, out_dir = file.path(tmp_root, "sdm4"), kept_vars = kv, focal_taxon = "Packera aurea", demo = FALSE)
  v1 <- as.numeric(terra::values(m1$focal$rast_ens))
  v4 <- as.numeric(terra::values(m4$focal$rast_ens))
  res$sdm <- list(
    n_rows = nrow(m1$performance),
    auc_diff = max(abs(m1$performance$AUC_Block_Ensemble - m4$performance$AUC_Block_Ensemble)),
    rast_fin = all(is.finite(v1)) && all(is.finite(v4)),
    rast_max_diff = max(abs(v1 - v4))
  )
}

jsonlite::write_json(res, out_json, auto_unbox = TRUE, digits = 12)
"""


class TestParallelThreadScaling(unittest.TestCase):
    """Verify --threads / PLANTMATICS_MAX_WORKERS scales safely across stages with deterministic results."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._r_res: dict = {}
        if SKIP_R_COND:
            return
        cls._tmp_dir = tempfile.TemporaryDirectory(prefix="par_scale_")
        tmp_p = Path(cls._tmp_dir.name)
        r_path, out_json = tmp_p / "run_par_parity.R", tmp_p / "res.json"
        r_path.write_text(_WARM_R_PARITY_SCRIPT, encoding="utf-8")
        proc = run_warm_r_cmd(
            ["Rscript", str(r_path), str(out_json), str(tmp_p)],
            stage_key="step05", cwd=PROJECT_ROOT, skip_figures=True,
        )
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(f"Batched R thread-scaling test failed (rc={proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls._r_res = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    def test_python_resolve_max_workers_and_fork_guard(self) -> None:
        from scripts.core.config import resolve_max_workers
        with mock.patch.dict(os.environ, {"PLANTMATICS_MAX_WORKERS": "58", "PACKERA_MAX_WORKERS": "58"}, clear=False):
            os.environ.pop("PLANTMATICS_IN_FORK_WORKER", None)
            self.assertEqual(resolve_max_workers(), 58)
            self.assertEqual(resolve_max_workers(cap=32), 32)
            with mock.patch.dict(os.environ, {"PLANTMATICS_IN_FORK_WORKER": "1"}, clear=False):
                self.assertEqual(resolve_max_workers(), 1)
        with mock.patch.dict(os.environ, {"PLANTMATICS_MAX_WORKERS": "", "PACKERA_MAX_WORKERS": "44"}, clear=False):
            os.environ.pop("PLANTMATICS_IN_FORK_WORKER", None)
            self.assertEqual(resolve_max_workers(), 44)

    def test_eastern_na_monograph_species_jobs_parallel_dispatch(self) -> None:
        from scripts.pipeline.resolve_eastern_na_packera import build_parser, main as monograph_main
        args = build_parser().parse_args(["--dry-run", "-j", "58", "--species-jobs", "2"])
        self.assertEqual((args.threads, args.species_jobs), (58, 2))
        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp).resolve(); (ws / "run_pipeline.sh").write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            db_p = ws / "data" / "plantmatics.sqlite"
            calls = []
            def _rec_run(cmd, *a, **kw):
                calls.append((list(cmd), kw.get("env")))
                return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="000000\n", stderr="")
            with mock.patch("scripts.pipeline.resolve_eastern_na_packera.discover_qualifying_eastern_species", return_value=[("Packera aurea", 20), ("Packera dubia", 20)]), \
                 mock.patch("subprocess.run", side_effect=_rec_run):
                rc = monograph_main(["--skip-harvest", "--skip-discovery", "-j", "58", "--species-jobs", "2", "--db", str(db_p)], workspace_root=ws)
            self.assertEqual(rc, 0)
            sp_calls = [c for c in calls if any("Packera_aurea" in x or "Packera_dubia" in x for x in c[0])]
            self.assertEqual(len(sp_calls), 6)
            for cmd, env_d in sp_calls:
                self.assertIsNotNone(env_d)
                self.assertEqual(env_d.get("PLANTMATICS_NO_LATEST_SYMLINK"), "1")
                if "--threads" in cmd:
                    self.assertEqual(cmd[cmd.index("--threads") + 1], "29")

    def test_coalescent_models_and_svdquartets_bit_exact_across_workers(self) -> None:
        from scripts.analysis.step06_coalescent_alignment_utils import compute_svdquartets_support
        from scripts.analysis.step06_coalescent_models import compute_model_corrected_divergence_blocks
        rng = np.random.default_rng(7)
        a1 = rng.integers(0, 4, size=(6, 160), dtype=np.int32)
        a2 = a1.copy()
        blks = [(i * 10, (i + 1) * 10) for i in range(16)]
        sids = [f"S{i}" for i in range(6)]
        taxa = ["TaxA", "TaxA", "TaxB", "TaxC", "TaxD", "TaxE"]
        with mock.patch.dict(os.environ, {"PLANTMATICS_MAX_WORKERS": "1", "PACKERA_MAX_WORKERS": "1"}):
            r1 = compute_model_corrected_divergence_blocks(a1, a2, blks)
            q1 = compute_svdquartets_support(a1, a2, taxa, sids, n_boot=25, seed=42)
        with mock.patch.dict(os.environ, {"PLANTMATICS_MAX_WORKERS": "8", "PACKERA_MAX_WORKERS": "8"}):
            r8 = compute_model_corrected_divergence_blocks(a1, a2, blks)
            q8 = compute_svdquartets_support(a1, a2, taxa, sids, n_boot=25, seed=42)
        for arr1, arr8 in zip(r1, r8):
            np.testing.assert_array_equal(arr1, arr8)
        self.assertTrue(q1.equals(q8))

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_r_worker_cap_fork_guard_and_stage_parity(self) -> None:
        wk = self._r_res["workers"]
        self.assertEqual((int(wk["w58"]), int(wk["w58_after"])), (58, 58))
        self.assertEqual(wk["in_fork_vals"], [1, 1, 1, 1])
        cr = self._r_res["circ"]
        self.assertLessEqual(float(cr["wf_diff"]), 1e-12)
        self.assertLessEqual(float(cr["pp_diff"]), 1e-12)
        self.assertTrue(bool(cr["sum_identical"]))
        self.assertLessEqual(float(self._r_res["efa_max_diff"]), 1e-12)
        gm = self._r_res["gmm"]
        self.assertLessEqual(float(gm["icl_diff"]), 1e-12)
        self.assertLessEqual(float(gm["lrt_diff"]), 1e-12)
        self.assertLessEqual(float(gm["p_diff"]), 1e-12)
        if "sdm" in self._r_res:
            sdm = self._r_res["sdm"]
            self.assertEqual(int(sdm["n_rows"]), 2)
            self.assertLessEqual(float(sdm["auc_diff"]), 1e-12)
            self.assertTrue(bool(sdm["rast_fin"]))
            self.assertLessEqual(float(sdm["rast_max_diff"]), 1e-12)


if __name__ == "__main__":
    unittest.main()
