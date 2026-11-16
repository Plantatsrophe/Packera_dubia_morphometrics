"""Batched warm-R verification tests for Step 04 and Step 05 statistical model selection."""

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

SKIP_R_REASON = "Skipping batched R Step 04/05 model selection tests under fast audit or missing Rscript."
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _is_finite(val: object) -> bool:
    """Return True when a decoded JSON scalar represents a finite real number."""
    if val is None or isinstance(val, bool):
        return False
    try:
        return math.isfinite(float(val))
    except (TypeError, ValueError):
        return False


_WARM_R_BENCHMARK = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]

suppressPackageStartupMessages({
  if (requireNamespace("jsonlite", quietly = TRUE)) library(jsonlite)
})

source("scripts/morphometrics/04_gmm_morphotools.R")
source("scripts/analysis/edaphic_fine_grain_utils.R")
source("scripts/analysis/step05_sdm_maxent.R")
source("scripts/analysis/step05_phenology_ecotypes.R")
if (!exists("select_circular_hopkins_baseline_model", mode = "function")) {
  source("scripts/analysis/step05_circular_hopkins_models.R")
}

results <- list()

# ------------------------------------------------------------------------------
# 1. Step 04 GMM ICL + Parametric Bootstrap K=1 vs K=2 LRT Gate
# ------------------------------------------------------------------------------
set.seed(202611L)
n_g <- 90L
z1 <- stats::rnorm(n_g)
df_gmm_a <- data.frame(
  catalogNumber = sprintf("G_A_%03d", seq_len(n_g)),
  plant_individual_id = rep(1L, n_g),
  species_standardized = rep("Packera dubia", n_g),
  scientificName = rep("Packera dubia", n_g),
  determiner_tier = rep("Tier_1_Gold", n_g),
  morph_tier_source = rep("Basal_Rosette", n_g),
  solidity = rep(0.75, n_g),
  is_dissected = rep(FALSE, n_g),
  exsiccatae_key = paste0("EX_", seq_len(n_g)),
  PC1 = exp(0.55 * z1) - exp((0.55^2) / 2.0),
  PC2 = stats::rnorm(n_g),
  PC3 = stats::rnorm(n_g),
  stringsAsFactors = FALSE
)
res_gmm_a <- run_gmm_cluster_analysis(
  df_gmm_a, c("PC1", "PC2", "PC3"),
  max_k = 4L, use_clustvarsel = FALSE, intraspecific_mode = TRUE, seed = 202611L
)

df_gmm_b <- data.frame(
  catalogNumber = sprintf("G_B_%03d", seq_len(n_g)),
  plant_individual_id = rep(1L, n_g),
  species_standardized = rep("Packera dubia", n_g),
  scientificName = rep("Packera dubia", n_g),
  determiner_tier = rep("Tier_1_Gold", n_g),
  morph_tier_source = rep("Basal_Rosette", n_g),
  solidity = rep(0.75, n_g),
  is_dissected = rep(FALSE, n_g),
  exsiccatae_key = paste0("EX_B_", seq_len(n_g)),
  PC1 = c(stats::rnorm(45L, mean = -3.8, sd = 1.0), stats::rnorm(45L, mean = 3.8, sd = 1.0)),
  PC2 = stats::rnorm(n_g),
  PC3 = stats::rnorm(n_g),
  stringsAsFactors = FALSE
)
res_gmm_b <- run_gmm_cluster_analysis(
  df_gmm_b, c("PC1", "PC2", "PC3"),
  max_k = 4L, use_clustvarsel = FALSE, intraspecific_mode = TRUE, seed = 202611L
)

results$step04_gmm <- list(
  case_a = list(best_k = as.integer(res_gmm_a$best_k), bic_table = as.data.frame(res_gmm_a$bic_table)),
  case_b = list(best_k = as.integer(res_gmm_b$best_k), bic_table = as.data.frame(res_gmm_b$bic_table))
)

# ------------------------------------------------------------------------------
# 2. Step 05 Fold-Internal Predictor Screening & Zero-Leakage Imputation
# ------------------------------------------------------------------------------
set.seed(202611L)
n_cv <- 48L
block_id <- rep(1:4, each = 12L)
y_cv <- factor(unlist(lapply(1:4, function(b) sample(c(rep("A", 6L), rep("B", 6L))))))

clim_10 <- c(
  "bio1_temp_mean", "bio2_mean_diurnal_range", "bio4_temp_seasonality",
  "bio5_max_temp_warmest_month", "bio6_min_temp_coldest_month", "bio10_mean_temp_warmest_qtr",
  "bio11_mean_temp_coldest_qtr", "bio12_precip_annual", "bio14_precip_driest_month", "bio15_precip_seasonality"
)
edaph_14 <- c(
  "fine_ph_0_15cm", "fine_clay_pct", "fine_sand_pct", "fine_silt_pct",
  "fine_soc_g_kg", "fine_cec_cmol_kg", "fine_bdod_kg_dm3", "fine_bedrock_depth_cm",
  "polaris_ph_0_15cm", "solus_ph_0_15cm", "solus_clay_0_15cm", "ssurgo_aws_150cm",
  "lithic_outcrop_index", "ultramafic_calcareous_gradient"
)
topo_4 <- c("topo_slope_deg", "topo_tpi_30m", "topo_hli_30m", "elevation_m")
all_28 <- c(clim_10, edaph_14, topo_4)

z_c <- stats::rnorm(n_cv); z_e <- stats::rnorm(n_cv); z_t <- stats::rnorm(n_cv)
X_list <- list()
for (nm in clim_10) X_list[[nm]] <- 0.85 * z_c + 0.35 * stats::rnorm(n_cv)
for (nm in edaph_14) X_list[[nm]] <- 0.85 * z_e + 0.35 * stats::rnorm(n_cv)
for (nm in topo_4) X_list[[nm]] <- 0.85 * z_t + 0.35 * stats::rnorm(n_cv)
X_base_df <- as.data.frame(X_list, stringsAsFactors = FALSE, check.names = FALSE)
for (nm in all_28) {
  na_pos <- sample.int(n_cv, size = 3L)
  X_base_df[[nm]][na_pos] <- NA_real_
}

oof_pred <- rep(NA_character_, n_cv)
fold_p <- integer(4L)
fold_limits <- integer(4L)
max_abs_tr_med <- 0.0
max_med_diff <- 0.0

for (f in 1:4) {
  tr_i <- which(block_id != f)
  te_i <- which(block_id == f)
  fold_limits[f] <- as.integer(floor(length(tr_i) / 5.0))
  X_fold_df <- X_base_df
  X_fold_df[te_i[1L], all_28] <- 9e6
  scr_f <- screen_predictors_within_cv_fold(
    X_fold_df[tr_i, , drop = FALSE],
    cand_vars = all_28, clim_vars = clim_10, edaph_vars = edaph_14, topo_vars = topo_4,
    y_tr = y_cv[tr_i], cor.threshold = 0.70, vif.threshold = 5.0, intraspecific_mode = FALSE
  )
  scr_clean <- screen_predictors_within_cv_fold(
    X_base_df[tr_i, , drop = FALSE],
    cand_vars = all_28, clim_vars = clim_10, edaph_vars = edaph_14, topo_vars = topo_4,
    y_tr = y_cv[tr_i], cor.threshold = 0.70, vif.threshold = 5.0, intraspecific_mode = FALSE
  )
  fold_p[f] <- length(scr_f$selected_vars)
  max_abs_tr_med <- max(max_abs_tr_med, max(abs(scr_f$train_medians), na.rm = TRUE))
  max_med_diff <- max(max_med_diff, max(abs(scr_f$train_medians - scr_clean$train_medians), na.rm = TRUE))

  X_te_imp <- apply_fold_imputation(X_fold_df[te_i, , drop = FALSE], scr_f$selected_vars, scr_f$train_medians)
  fit_r <- ranger::ranger(
    y ~ ., data = data.frame(y = y_cv[tr_i], scr_f$X_tr_imp, check.names = FALSE),
    num.trees = 200L, probability = TRUE, num.threads = 1L, seed = 202611L + f
  )
  pr_m <- stats::predict(fit_r, data = X_te_imp, num.threads = 1L)$predictions
  oof_pred[te_i] <- colnames(pr_m)[max.col(pr_m, ties.method = "first")]
}

tab_cv <- table(factor(y_cv, levels = c("A", "B")), factor(oof_pred, levels = c("A", "B")))
p_obs <- sum(diag(tab_cv)) / n_cv
p_exp <- sum(rowSums(tab_cv) * colSums(tab_cv)) / (n_cv^2)
oof_kappa <- if ((1.0 - p_exp) > 1e-12) (p_obs - p_exp) / (1.0 - p_exp) else 0.0

results$step05_cv <- list(
  fold_p = as.integer(fold_p),
  fold_limits = as.integer(fold_limits),
  max_abs_train_median = as.numeric(max_abs_tr_med),
  max_median_diff_vs_clean = as.numeric(max_med_diff),
  oof_kappa = as.numeric(oof_kappa)
)

# ------------------------------------------------------------------------------
# 3. Step 05 ENMeval-Style IPP MaxEnt AICc Feature/Regularization Selection
# ------------------------------------------------------------------------------
set.seed(202611L)
x1_pool <- stats::rnorm(2000L)
p_idx <- sample.int(2000L, size = 18L, replace = FALSE, prob = exp(1.8 * x1_pool))
X_pres_mx <- data.frame(X1 = x1_pool[p_idx], X2 = stats::rnorm(18L), X3 = stats::rnorm(18L))
X_bg_mx <- data.frame(X1 = stats::rnorm(200L), X2 = stats::rnorm(200L), X3 = stats::rnorm(200L))
p_mx <- c(rep(1L, 18L), rep(0L, 200L))
X_mx <- rbind(X_pres_mx, X_bg_mx)
mx_fit <- .fit_single_maxent_model(p_mx, X_mx)

X_fb_mx <- data.frame(
  X1 = c(stats::rnorm(5L, mean = 1.5, sd = 1.0), stats::rnorm(60L)),
  X2 = stats::rnorm(65L),
  X3 = stats::rnorm(65L)
)
mx_fb <- .fit_single_maxent_model(c(rep(1L, 5L), rep(0L, 60L)), X_fb_mx)
pred_fb <- as.numeric(mx_fb$predict_cloglog(X_fb_mx[1:10, , drop = FALSE]))

results$step05_maxent <- list(
  selected_feature_class = as.character(mx_fit$selected_feature_class),
  n_nonzero_betas = as.integer(mx_fit$n_nonzero_betas),
  maxent_aicc = as.numeric(mx_fit$maxent_aicc),
  pred_fb = pred_fb
)

# ------------------------------------------------------------------------------
# 4. Step 05 Circular von Mises Hopkins AICc Selection & Wrap-Around Phase Tests
# ------------------------------------------------------------------------------
.rvonmises_bf <- function(n, mu = 0.0, kappa = 12.0) {
  if (kappa <= 1e-8) return(stats::runif(n, -pi, pi))
  tau <- 1.0 + sqrt(1.0 + 4.0 * kappa^2)
  rho <- (tau - sqrt(2.0 * tau)) / (2.0 * kappa)
  r <- (1.0 + rho^2) / (2.0 * rho)
  out <- numeric(n)
  for (i in seq_len(n)) {
    repeat {
      z <- cos(pi * stats::runif(1L))
      f <- (1.0 + r * z) / (r + z)
      c_v <- kappa * (r - f)
      u2 <- stats::runif(1L)
      if (c_v * (2.0 - c_v) - u2 > 0.0 || log(c_v / u2) + 1.0 - c_v >= 0.0) {
        ang <- sign(stats::runif(1L) - 0.5) * acos(pmin(1.0, pmax(-1.0, f)))
        out[i] <- atan2(sin(mu + ang), cos(mu + ang))
        break
      }
    }
  }
  out
}

set.seed(202611L)
lat_h <- seq(30.0, 42.0, length.out = 45L)
lat_hc <- lat_h - mean(lat_h)
eps_d <- .rvonmises_bf(45L, mu = 0.0, kappa = 12.0) * (365.25 / (2.0 * pi))
eps_d <- eps_d - as.numeric(stats::coef(stats::lm(eps_d ~ lat_hc))[2L]) * lat_hc
doy_h <- as.numeric(((round(364.0 + 1.8 * lat_hc + eps_d) - 1L) %% 365L) + 1L)
doy_h[11L] <- 358.0
doy_h[43L] <- 8.0
df_hop_a <- data.frame(decimalLatitude = lat_h, doy = doy_h, stringsAsFactors = FALSE)

sel_hop_a <- select_circular_hopkins_baseline_model(df_hop_a)
base_hop_a <- model_latitudinal_spring_baseline(df_hop_a)

mu0_rad <- doy_to_rad(363.0)
th_g1 <- wrap_rad(mu0_rad + .rvonmises_bf(24L, 0.0, 12.0))
th_g2 <- wrap_rad(mu0_rad + doy_to_rad(35.0) + .rvonmises_bf(24L, 0.0, 12.0))
phase_shift_res <- run_circular_multigroup_phase_test(
  c(th_g1, th_g2), c(rep("Early", 24L), rep("Late", 24L)), B = 999L, seed = 202611L
)

th_unif <- (seq(0.05, 2.0 * pi - 0.15, length.out = 30L) + stats::rnorm(30L, 0.0, 0.03)) %% (2.0 * pi)
ray_unif <- run_rayleigh_uniformity_test(th_unif)
df_unif <- data.frame(
  decimalLatitude = seq(30.0, 42.0, length.out = 30L),
  doy = as.numeric(((round(th_unif * (365.25 / (2.0 * pi))) - 1L) %% 365L) + 1L),
  stringsAsFactors = FALSE
)
sel_unif <- select_circular_hopkins_baseline_model(df_unif)

results$step05_hopkins <- list(
  sel_model = as.character(sel_hop_a$selected_model),
  base_sel_model = as.character(base_hop_a$selected_model),
  delta_aicc_vs_hopkins1918 = as.numeric(sel_hop_a$delta_aicc_vs_hopkins1918),
  max_abs_delta_doy = as.numeric(max(abs(base_hop_a$df$delta_doy), na.rm = TRUE)),
  shift_perm_p = as.numeric(phase_shift_res$perm_p),
  rayleigh_p_direct = as.numeric(ray_unif$p_value),
  rayleigh_p_model = as.numeric(sel_unif$rayleigh_p)
)

jsonlite::write_json(results, out_json, auto_unbox = TRUE, na = "null")
"""


@unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
class TestStep0405ModelSelection(unittest.TestCase):
    """Verify Step 04 ICL/bootstrap LRT and Step 05 CV, MaxEnt AICc, and circular Hopkins gates."""

    results: dict = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            cls.results = {}
            return
        with tempfile.TemporaryDirectory(prefix="step0405_ms_") as tmp_dir:
            tmp_path = Path(tmp_dir)
            r_script_path = tmp_path / "run_step04_05_benchmarks.R"
            out_json_path = tmp_path / "results.json"
            r_script_path.write_text(_WARM_R_BENCHMARK, encoding="utf-8")
            proc = subprocess.run(
                ["Rscript", str(r_script_path), str(out_json_path)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
                timeout=180,
            )
            if proc.returncode != 0 or not out_json_path.exists():
                raise RuntimeError(
                    f"Warm R benchmark failed (rc={proc.returncode}):\n"
                    f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
                )
            cls.results = json.loads(out_json_path.read_text(encoding="utf-8"))

    def test_step04_gmm_icl_and_bootstrap_lrt_gate(self) -> None:
        """Skewed unimodal null collapses to K=1 while separated K=2 clusters pass ICL and LRT."""
        gmm = self.results["step04_gmm"]
        case_a = gmm["case_a"]
        self.assertEqual(int(case_a["best_k"]), 1)
        for row in case_a["bic_table"]:
            self.assertTrue(_is_finite(row.get("ICL")))
            self.assertTrue(_is_finite(row.get("Delta_ICL_vs_Null")))
        k2_a = [r for r in case_a["bic_table"] if int(r["K"]) == 2][0]
        self.assertTrue(_is_finite(k2_a.get("Boot_LRT_P_K1_vs_K2")))

        case_b = gmm["case_b"]
        k2_b = [r for r in case_b["bic_table"] if int(r["K"]) == 2][0]
        self.assertGreater(float(k2_b["Delta_ICL_vs_Null"]), 0.0)
        self.assertLess(float(k2_b["Boot_LRT_P_K1_vs_K2"]), 0.05)
        self.assertEqual(int(case_b["best_k"]), 2)

    def test_step05_fold_internal_predictor_screening_zero_leakage(self) -> None:
        """Fold-internal screening respects EPV cap, blocks test-outlier leakage, and yields near-null kappa."""
        cv = self.results["step05_cv"]
        for fp, lim in zip(cv["fold_p"], cv["fold_limits"]):
            self.assertLessEqual(int(fp), int(lim))
        self.assertLess(float(cv["max_abs_train_median"]), 100.0)
        self.assertEqual(float(cv["max_median_diff_vs_clean"]), 0.0)
        self.assertLess(float(cv["oof_kappa"]), 0.20)

    def test_step05_ipp_maxent_aicc_grid_tuning(self) -> None:
        """ENMeval-style IPP AICc tuning favors parsimonious L/LQ models and preserves small-N fallback."""
        mx = self.results["step05_maxent"]
        self.assertIn(mx["selected_feature_class"], ("L", "LQ"))
        self.assertLess(int(mx["n_nonzero_betas"]), 18 - 1)
        self.assertTrue(_is_finite(mx["maxent_aicc"]))
        self.assertEqual(len(mx["pred_fb"]), 10)
        for pr in mx["pred_fb"]:
            self.assertTrue(_is_finite(pr))
            self.assertGreater(float(pr), 0.0)
            self.assertLess(float(pr), 1.0)

    def test_step05_circular_hopkins_aicc_and_wraparound(self) -> None:
        """Circular von Mises AICc picks M2 empirical slope, bounds wrap-around residuals, and gates uniform phases."""
        hop = self.results["step05_hopkins"]
        self.assertEqual(hop["sel_model"], "M2_Empirical_Lat")
        self.assertEqual(hop["base_sel_model"], "M2_Empirical_Lat")
        self.assertLess(float(hop["delta_aicc_vs_hopkins1918"]), -4.0)
        self.assertLessEqual(float(hop["max_abs_delta_doy"]), 182.625)
        self.assertLess(float(hop["shift_perm_p"]), 0.01)
        self.assertGreater(float(hop["rayleigh_p_direct"]), 0.05)
        self.assertGreater(float(hop["rayleigh_p_model"]), 0.05)


if __name__ == "__main__":
    unittest.main()
