"""Unit and batched-R tests for Step 06 evolutionary model selection and hypothesis testing."""

from __future__ import annotations

import importlib.util
import math
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tests._r_session import run_warm_r_cmd

SKIP_R_REASON = "Skipping batched R Step 06 model selection tests under fast audit or missing Rscript."
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _load_step06_ingest_module():
    return importlib.import_module("scripts.data_prep.06_ingest_genetic_data")


def _hky_transition_matrix(pi: np.ndarray, kappa: float, branch_len: float) -> np.ndarray:
    """Compute exact reversible HKY85 transition matrix P(t) = exp(t * Q)."""
    q = np.zeros((4, 4), dtype=np.float64)
    ts_pairs = {(0, 2), (2, 0), (1, 3), (3, 1)}
    for i in range(4):
        for j in range(4):
            if i != j:
                q[i, j] = (kappa if (i, j) in ts_pairs else 1.0) * pi[j]
        q[i, i] = -float(np.sum(q[i, :]))
    rate = -float(np.sum(pi * np.diag(q)))
    q /= max(1e-12, rate)
    sq_pi = np.sqrt(pi)
    inv_sq = 1.0 / sq_pi
    sym = (sq_pi[:, None] * q) * inv_sq[None, :]
    evals, evecs = np.linalg.eigh(sym)
    exp_evals = np.exp(evals * branch_len)
    p = (inv_sq[:, None] * evecs) @ (exp_evals[:, None] * evecs.T) * sq_pi[None, :]
    p = np.clip(p, 0.0, None)
    return p / np.sum(p, axis=1, keepdims=True)


def _simulate_hky_alignment(
    n_sites: int = 600,
    kappa: float = 6.0,
    pi_vec: Tuple[float, float, float, float] = (0.38, 0.12, 0.12, 0.38),
    seed: int = 202611,
) -> Tuple[List[str], Dict[str, str]]:
    """Simulate 6 taxa along an ultrametric tree under HKY85."""
    rng = np.random.default_rng(seed)
    pi = np.asarray(pi_vec, dtype=np.float64)
    bases = np.array(["A", "C", "G", "T"])

    def _evolve_step(states: np.ndarray, t: float) -> np.ndarray:
        p = _hky_transition_matrix(pi, kappa, t)
        cdf = np.cumsum(p, axis=1)
        u = rng.random(states.shape[0])
        out = np.empty_like(states)
        for b in range(4):
            mask = states == b
            if np.any(mask):
                out[mask] = np.searchsorted(cdf[b], u[mask], side="right").clip(0, 3)
        return out

    root = rng.choice(4, size=n_sites, p=pi)
    node_l = _evolve_step(root, 0.05)
    node_p = _evolve_step(root, 0.05)
    node_l1 = _evolve_step(node_l, 0.05)
    node_l2 = _evolve_step(node_l, 0.05)
    tips = {
        "T1": _evolve_step(node_l1, 0.04),
        "T2": _evolve_step(node_l1, 0.04),
        "T3": _evolve_step(node_l2, 0.04),
        "T4": _evolve_step(node_l2, 0.04),
        "T5": _evolve_step(node_p, 0.09),
        "T6": _evolve_step(node_p, 0.09),
    }
    sample_ids = ["T1", "T2", "T3", "T4", "T5", "T6"]
    seqs = {sid: "".join(bases[tips[sid]]) for sid in sample_ids}
    return sample_ids, seqs


class TestStep06ModelSelection(unittest.TestCase):
    """Test Python and R model selection and hypothesis testing across Step 06."""

    def test_substitution_model_selection_bic_and_asc(self):
        """BIC favors HKY85/TN93 over JC69/F81 and enables +ASC strictly on SNP-only data."""
        from scripts.analysis.step06_coalescent_models import select_nucleotide_substitution_model

        sample_ids, full_seqs = _simulate_hky_alignment(n_sites=600, kappa=6.0, seed=202611)
        mat = np.array([list(full_seqs[sid]) for sid in sample_ids])
        var_cols = [j for j in range(mat.shape[1]) if len(set(mat[:, j])) > 1]
        self.assertGreaterEqual(len(var_cols), 120)
        snp_seqs = {sid: "".join(mat[i, var_cols]) for i, sid in enumerate(sample_ids)}

        prev_skip = os.environ.get("PLANTMATICS_SKIP_IQTREE2")
        os.environ["PLANTMATICS_SKIP_IQTREE2"] = "1"
        try:
            res_full = select_nucleotide_substitution_model(full_seqs, sample_ids, alignment_name="full_test.phylip")
            res_snp = select_nucleotide_substitution_model(snp_seqs, sample_ids, alignment_name="snp_test.phylip")
        finally:
            if prev_skip is None:
                os.environ.pop("PLANTMATICS_SKIP_IQTREE2", None)
            else:
                os.environ["PLANTMATICS_SKIP_IQTREE2"] = prev_skip

        self.assertIn(res_full["base_model"], ("HKY85", "TN93"))
        self.assertNotIn(res_full["best_model"], ("JC69", "F81"))
        self.assertFalse(bool(res_full["asc_applied"]))
        self.assertGreater(float(res_full["ti_tv_kappa"]), 2.0)

        self.assertIn(res_snp["base_model"], ("HKY85", "TN93"))
        self.assertTrue(bool(res_snp["asc_applied"]))
        self.assertIn("+ASC", str(res_snp["best_model"]))
        self.assertGreater(len(res_snp["table_rows"]), 0)

    def test_em_nquire_ploidy_model_selection(self):
        """EM Gaussian + uniform-noise BIC mixture separates 2x diploids from 4x autotetraploids."""
        ingest_mod = _load_step06_ingest_module()
        infer_fn = ingest_mod._infer_sample_cytotype

        rng = np.random.default_rng(202611)
        n_loci = 100
        dosage = np.full(n_loci, 0.50, dtype=np.float64)

        dp_2x = rng.integers(40, 66, size=n_loci)
        is_noise_2x = rng.random(n_loci) < 0.15
        ratios_2x = np.where(
            is_noise_2x,
            rng.uniform(0.10, 0.90, size=n_loci),
            np.clip(rng.normal(0.50, 0.045, size=n_loci), 0.12, 0.88),
        )
        alt_2x = np.clip(np.round(ratios_2x * dp_2x), 1, dp_2x - 1).astype(np.float64)
        ref_2x = (dp_2x - alt_2x).astype(np.float64)
        lbl_2x, delta_2x = infer_fn(ref_2x, alt_2x, dosage, user_ploidy=None)
        self.assertEqual(lbl_2x, "2x")
        self.assertIsNotNone(delta_2x)

        dp_4x = rng.integers(40, 66, size=n_loci)
        is_noise_4x = rng.random(n_loci) < 0.15
        centers_4x = rng.choice([0.25, 0.50, 0.75], size=n_loci)
        ratios_4x = np.where(
            is_noise_4x,
            rng.uniform(0.10, 0.90, size=n_loci),
            np.clip(rng.normal(centers_4x, 0.045), 0.12, 0.88),
        )
        alt_4x = np.clip(np.round(ratios_4x * dp_4x), 1, dp_4x - 1).astype(np.float64)
        ref_4x = (dp_4x - alt_4x).astype(np.float64)
        lbl_4x, delta_4x = infer_fn(ref_4x, alt_4x, dosage, user_ploidy=None)
        self.assertEqual(lbl_4x, "4x")
        self.assertIsNotNone(delta_4x)
        self.assertGreaterEqual(float(delta_4x), 4.0)

    def test_coalescent_introgression_bayes_factor(self):
        """Coalescent mixture BF rejects introgression under pure ILS and supports bimodal secondary contact."""
        from scripts.analysis.step06_coalescent_models import compute_coalescent_introgression_bf

        rng = np.random.default_rng(202611)
        n_blocks = 80
        call_len = 200.0
        call_vec = np.full(n_blocks, call_len, dtype=np.float64)

        t_ils = rng.exponential(scale=0.015, size=n_blocks)
        k_ils = rng.poisson(lam=call_len * t_ils)
        d_ils = k_ils.astype(np.float64) / call_len
        res_ils = compute_coalescent_introgression_bf(d_ils, call_vec, bpp_out_path=None)
        self.assertIn(res_ils["bf_method"], ("bic_coalescent_mixture", "gauss_laguerre_bic"))
        self.assertLess(float(res_ils["bf_10"]), 1.0)
        self.assertTrue(math.isfinite(float(res_ils["bic_m0"])))

        n_in = int(round(0.30 * n_blocks))
        tau_vec = np.concatenate([np.full(n_in, 0.001), np.full(n_blocks - n_in, 0.070)])
        t_sec = tau_vec + rng.exponential(scale=0.005, size=n_blocks)
        k_sec = rng.poisson(lam=call_len * t_sec)
        d_sec = k_sec.astype(np.float64) / call_len
        res_sec = compute_coalescent_introgression_bf(d_sec, call_vec, bpp_out_path=None)
        self.assertIn(res_sec["bf_method"], ("bic_coalescent_mixture", "gauss_laguerre_bic"))
        self.assertGreater(float(res_sec["bf_10"]), 3.0)
        self.assertLess(float(res_sec["tau_introgr"]), float(res_sec["tau_sp"]))
        self.assertGreater(float(res_sec["w_introgr"]), 0.05)

    @unittest.skipIf(SKIP_R_COND, SKIP_R_REASON)
    def test_batched_r_step06_model_selection(self):
        """Verify R admixture K CV, Ledoit-Wolf shrinkage, hybrid AICc/BH-LRT, SNPs2CF G-test, and MMRR/FEEMS CV."""
        with tempfile.TemporaryDirectory(prefix="step06_ms_r_") as tmp_dir:
            tmp_posix = Path(tmp_dir).as_posix()
            r_script = f"""
grDevices::pdf(NULL)
Sys.setenv(PLANTMATICS_FORCE_NMF_CE = "1")
source("scripts/core/db_utils.R")
source("scripts/analysis/step06_model_selection.R")
source("scripts/analysis/step06_reticulation_coancestry.R")
source("scripts/analysis/step06_spatial_landscape_genomics.R")

set.seed(202611L)

# 1. Cross-entropy K selection: K=1 panmixia vs K=3 strong demes
p_pan <- stats::runif(120L, 0.25, 0.75)
dos_pan <- matrix(stats::rbinom(24L * 120L, 2L, rep(p_pan, each = 24L)) / 2.0, nrow = 24L, ncol = 120L)
rownames(dos_pan) <- sprintf("PAN_%02d", seq_len(24L))
res_pan <- select_admixture_k_cross_entropy(
  dos_pan, d_mat = as.matrix(stats::dist(dos_pan)),
  k_max = 5L, n_folds = 5L, seed = 202611L,
  out_csv = "{tmp_posix}/ce_pan.csv"
)
stopifnot(res_pan$best_k == 1L, file.exists("{tmp_posix}/ce_pan.csv"))

p_mat <- matrix(0.08, nrow = 3L, ncol = 120L)
p_mat[1L, 1:40] <- 0.92; p_mat[2L, 41:80] <- 0.92; p_mat[3L, 81:120] <- 0.92
dos_k3 <- do.call(rbind, lapply(1:3, function(g) {{
  matrix(stats::rbinom(8L * 120L, 2L, rep(p_mat[g, ], each = 8L)) / 2.0, nrow = 8L, ncol = 120L)
}}))
rownames(dos_k3) <- sprintf("K3_%02d", seq_len(24L))
res_k3 <- select_admixture_k_cross_entropy(
  dos_k3, d_mat = as.matrix(stats::dist(dos_k3)),
  k_max = 5L, n_folds = 5L, seed = 202611L,
  out_csv = "{tmp_posix}/admixture_k_model_selection.csv"
)
stopifnot(res_k3$best_k == 3L, is.finite(res_k3$selected_k_cross_entropy))
stopifnot(res_k3$selected_k_cross_entropy < res_k3$table_df$masked_cross_entropy[1L])

# 2. Analytical Ledoit-Wolf shrinkage monotonic in sample size under anisotropic covariance
sc_lw <- seq(4.0, 0.3, length.out = 40L)
X_small <- sweep(matrix(stats::rnorm(6L * 40L), nrow = 6L, ncol = 40L), 2L, sc_lw, "*")
X_large <- sweep(matrix(stats::rnorm(60L * 40L), nrow = 60L, ncol = 40L), 2L, sc_lw, "*")
lam_small <- estimate_ledoit_wolf_lambda(X_small)
lam_large <- estimate_ledoit_wolf_lambda(X_large)
stopifnot(is.finite(lam_small), is.finite(lam_large), lam_small >= 0.01, lam_large >= 0.01, lam_small > lam_large)

# 3. 7-model multinomial AICc + BH-LRT on hybrid triangle
m_hyb <- 40L
p1_ref <- matrix(rep(c(rep(0.0, 30L), rep(0.5, 10L)), each = 4L), nrow = 4L, ncol = m_hyb)
p2_ref <- matrix(rep(c(rep(1.0, 30L), rep(0.5, 10L)), each = 4L), nrow = 4L, ncol = m_hyb)
q_noisy <- c(0.5, 0.5, rep(0.0, 28L), rep(0.5, 10L))
q_f1 <- c(rep(0.5, 30L), rep(0.5, 10L))
q_div <- c(rep(0.0, 30L), rep(0.5, 10L))
dos_hyb <- rbind(p1_ref, p2_ref, q_noisy, q_f1, q_div)
hyb_ids <- c(sprintf("P1_%d", 1:4), sprintf("P2_%d", 1:4), "QRY_NOISY_P1", "QRY_F1", "QRY_DIV")
rownames(dos_hyb) <- hyb_ids
pass_hyb <- data.frame(
  genetic_sample_id = hyb_ids,
  species_standardized = c(rep("Parent_1", 4L), rep("Parent_2", 4L), rep("Query_Sp", 3L)),
  is_putative_novel_or_hybrid_query = c(rep(0L, 8L), rep(1L, 3L)),
  stringsAsFactors = FALSE
)
assign_hyb <- data.frame(
  genetic_sample_id = hyb_ids,
  geno_parent_1 = "Parent_1",
  geno_parent_2 = "Parent_2",
  geno_predicted_taxon = c(rep("Parent_1", 4L), rep("Parent_2", 4L), "Parent_1", "Parent_1", "Parent_3"),
  geno_posterior_prob = c(rep(0.99, 8L), 0.95, 0.55, 0.96),
  is_geno_admixed = c(rep(FALSE, 8L), FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
hyb_res <- compute_hybrid_index_q12(dos_hyb, pass_hyb, assign_hyb)
tri <- hyb_res$triangle_df
stopifnot(all(c("best_model_aicc", "akaike_weight", "delta_aicc_vs_pure", "lrt_q_bh") %in% names(tri)))
r_noisy <- tri[tri$genetic_sample_id == "QRY_NOISY_P1", , drop = FALSE]
r_f1 <- tri[tri$genetic_sample_id == "QRY_F1", , drop = FALSE]
r_div <- tri[tri$genetic_sample_id == "QRY_DIV", , drop = FALSE]
stopifnot(r_noisy$hybrid_generation_class == "Pure_Parent_1", r_noisy$lrt_q_bh >= 0.05)
stopifnot(r_f1$hybrid_generation_class == "F1_Hybrid", r_f1$lrt_q_bh < 0.05, r_f1$delta_aicc_vs_pure >= 2.0)
stopifnot(r_div$hybrid_generation_class == "Divergent_NonHybrid_Lineage")

# 4. Symmetric ILS quartet concordance zeroes snaq_gamma_minor
dos_q <- cbind(
  matrix(rep(c(1, 1, 0, 0), 60L), nrow = 4L, ncol = 60L),
  matrix(rep(c(1, 0, 1, 0), 20L), nrow = 4L, ncol = 20L),
  matrix(rep(c(1, 0, 0, 1), 20L), nrow = 4L, ncol = 20L)
)
pass_q <- data.frame(
  genetic_sample_id = c("S_A", "S_B", "S_C", "S_D"),
  species_standardized = c("TaxonA", "TaxonB", "TaxonC", "TaxonD"),
  stringsAsFactors = FALSE
)
rownames(dos_q) <- pass_q$genetic_sample_id
cf_res <- compute_quartet_concordance_snps2cf(pass_q, dos_q)
stopifnot(nrow(cf_res$cf_df) == 1L, cf_res$cf_df$snaq_gamma_minor[1L] == 0.0, cf_res$cf_df$g_test_p[1L] >= 0.05)

# 5. Nested MMRR selection and leave-one-node-out FEEMS CV
ls_pass <- data.frame(
  genetic_sample_id = sprintf("LS_%02d", 1:10),
  species_standardized = rep("P_dubia", 10L),
  population_id = sprintf("POP_%d", rep(1:5, each = 2L)),
  decimalLatitude = c(35.2, 35.4, 36.1, 36.3, 37.5, 37.7, 38.8, 39.0, 40.1, 40.3),
  decimalLongitude = c(-83.5, -83.2, -81.8, -81.5, -80.2, -79.9, -78.9, -78.6, -77.3, -77.0),
  stringsAsFactors = FALSE
)
coord_mat <- cbind(ls_pass$decimalLatitude, ls_pass$decimalLongitude)
d_geno_ls <- as.matrix(stats::dist(coord_mat))
d_geno_ls <- d_geno_ls / max(1e-6, max(d_geno_ls))
jitter_m <- matrix(stats::rnorm(100L, 0.0, 0.05), nrow = 10L, ncol = 10L)
d_geno_ls <- pmax(d_geno_ls + 0.5 * (jitter_m + t(jitter_m)), 0.0)
diag(d_geno_ls) <- 0.0
dimnames(d_geno_ls) <- list(ls_pass$genetic_sample_id, ls_pass$genetic_sample_id)
ls_out_dir <- "{tmp_posix}/landscape_out"
ls_res <- run_step06_landscape_genomics(
  pass_df = ls_pass, d_geno_mat = d_geno_ls,
  db_path = "", run_id = "test_step06_ms", demo = FALSE,
  n_perm = 99L, out_dir = ls_out_dir
)
mmrr_csv <- file.path(ls_out_dir, "mmrr_model_selection.csv")
stopifnot(file.exists(mmrr_csv))
mmrr_df <- utils::read.csv(mmrr_csv, stringsAsFactors = FALSE)
stopifnot(nrow(mmrr_df) >= 4L, "IBD_only" %in% mmrr_df$model_name, "IBD_plus_IBE" %in% mmrr_df$model_name)
stopifnot("feems_cv_lambda" %in% names(ls_res$summary_df), is.finite(ls_res$summary_df$feems_cv_lambda[1L]))
cat("ALL_STEP06_MS_R_CHECKS_OK\\n")
"""
            proc = run_warm_r_cmd(
                ["Rscript", "-e", r_script],
                stage_key="step06",
                cwd=PROJECT_ROOT,
                skip_figures=True,
            )
            self.assertEqual(
                proc.returncode,
                0,
                f"Batched R Step 06 model selection failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
            )
            self.assertIn("ALL_STEP06_MS_R_CHECKS_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
