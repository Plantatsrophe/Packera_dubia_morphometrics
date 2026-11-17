#!/usr/bin/env python3
"""Regression & code-hygiene suite for the refactored Step 07 synthesis stack."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

STEP07_CEILINGS: Dict[Path, int] = {
    PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R": 500,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_hopkins_consensus_utils.R": 480,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_ancova_merge_utils.R": 490,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_reports_html_utils.R": 490,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_biogeo_hybrid_utils.R": 490,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_cluster_cohesion_utils.R": 460,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_decision_stages_utils.R": 490,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_postprocess_guards_utils.R": 460,
    PROJECT_ROOT / "scripts" / "analysis" / "taxon_complex_discovery_utils.R": 480,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_complex_affinity_utils.R": 450,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_cytology_integration.R": 380,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_deme_propagation_utils.R": 260,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_genomic_propagation.R": 400,
    PROJECT_ROOT / "scripts" / "analysis" / "step07_organellar_triage_utils.R": 260,
    PROJECT_ROOT / "scripts" / "analysis" / "genomic_triage_utils.R": 500,
    PROJECT_ROOT / "scripts" / "figures" / "step07_figures.R": 420,
    PROJECT_ROOT / "scripts" / "figures" / "step07b_network_figures.R": 320,
}

_R_STR_RE = re.compile(r'"([^"\\]|\\.)*"|\'([^\'\\]|\\.)*\'', re.DOTALL)
_NON_POSIX_ESC_RE = re.compile(r"\\+[sdwb](?![A-Za-z])")

_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; png_a <- args[2L]; png_b <- args[3L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/morphometrics/pub_figures_utils.R")
source("scripts/analysis/07_triage_dashboard_synthesis.R")
source("scripts/figures/step07_figures.R")

# 1. Marginal drop1 ANCOVA invariance + Cauline/CIE exclusions
set.seed(11L)
df_anc <- data.frame(
  catalogNumber = sprintf("A_%02d", 1:27),
  species_standardized = c(rep("Packera aurea", 8L), rep("Packera dubia", 8L), rep("Packera anonyma", 8L), rep("Packera aurea", 3L)),
  institutionCode = c(rep("NCU", 6L), rep("DUKE", 2L), rep("DUKE", 6L), rep("US", 2L), rep("US", 6L), rep("NCU", 2L), rep("NCU", 3L)),
  regional_group = rep(c("Piedmont", "Mountains", "Coastal"), length.out = 27L),
  morph_tier_source = c(rep("Basal_Rosette", 24L), rep("Cauline_Only_Fallback", 3L)),
  selected_morph_features = rep("PC1;PC2;blade_cie_l_norm", 27L),
  PC1 = c(seq(-2, 2, length.out = 24L) + c(rep(0, 8L), rep(1.5, 8L), rep(3.0, 8L)), rep(500.0, 3L)),
  PC2 = rnorm(27L, 0, 0.4),
  blade_cie_l_norm = runif(27L, 40, 80),
  authority_weight = c(rep(2.0, 24L), rep(1.0, 3L)),
  spatial_balance_weight = rep(1.0, 27L),
  stringsAsFactors = FALSE
)
anc_1 <- audit_institutional_batch_effects_ancova(df_anc, min_obs = 4L)
df_anc_p <- df_anc; df_anc_p$PC1[df_anc_p$morph_tier_source == "Cauline_Only_Fallback"] <- -9999.0
anc_2 <- audit_institutional_batch_effects_ancova(df_anc_p, min_obs = 4L)
sub_c <- df_anc[1:24, ]
sub_c$w <- sub_c$authority_weight * sub_c$spatial_balance_weight
d1_f <- stats::drop1(stats::lm(PC1 ~ species_standardized + regional_group + institutionCode, data = sub_c, weights = w), scope = "institutionCode", test = "F")
d1_s <- stats::drop1(stats::lm(PC1 ~ institutionCode + species_standardized + regional_group, data = sub_c, weights = w), scope = "institutionCode", test = "F")
eta_f <- as.numeric(d1_f["institutionCode", "Sum of Sq"] / (d1_f["institutionCode", "Sum of Sq"] + d1_f["<none>", "RSS"]))
eta_s <- as.numeric(d1_s["institutionCode", "Sum of Sq"] / (d1_s["institutionCode", "Sum of Sq"] + d1_s["<none>", "RSS"]))

# 2a. G2_Conditional_Deme_Hull centroid quarantine
df_gp <- data.frame(
  catalogNumber = c("G1_A", "G1_B", "U_CENT", "U_CLEAN", "G1_CDON", "U_DON_C"),
  species_standardized = rep("Packera dubia", 6L),
  determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze", "Tier_3_Bronze", "Tier_1_Gold", "Tier_3_Bronze"),
  decimalLatitude = c(35.000, 35.002, 35.001, 35.001, 36.200, 36.201),
  decimalLongitude = c(-80.000, -80.002, -80.001, -80.001, -81.200, -81.201),
  coordinateUncertainty = rep(100.0, 6L),
  is_centroid_or_inst_coord = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE),
  morph_tier_source = rep("Basal_Rosette", 6L),
  PC1 = rep(0.1, 6L), PC2 = rep(0.2, 6L), doy = rep(120, 6L), hopkins_adjusted_doy = rep(120, 6L),
  geno_predicted_taxon = c("Packera dubia", "Packera dubia", NA, NA, "Packera dubia", NA),
  geno_posterior_prob = c(0.98, 0.97, NA, NA, 0.99, NA),
  geno_linkage_tier = c("G1_Direct_Voucher", "G1_Direct_Voucher", NA, NA, "G1_Direct_Voucher", NA),
  geno_quality_pass = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
  is_geno_admixed = rep(FALSE, 6L),
  stringsAsFactors = FALSE
)
gp_out <- run_step07_genomic_propagation(df_gp, opts = list())
gp_tiers <- stats::setNames(as.character(gp_out$genomic_propagation_tier), gp_out$catalogNumber)

# 2b. Stage 2 Gate 1/5 syntopy & Stage 3 spatial cline / Hopkins centroid quarantine
mk_syn <- function(cent_b) {
  data.frame(
    catalogNumber = c(paste0("PA_", 1:4), paste0("PB_", 1:4), "HY_1", "HY_2"),
    species_standardized = c(rep("Packera aurea", 4L), rep("Packera dubia", 4L), rep("Packera aurea", 2L)),
    determiner_tier = rep("Tier_1_Gold", 10L),
    is_centroid_or_inst_coord = c(rep(FALSE, 4L), rep(cent_b, 4L), FALSE, FALSE),
    decimalLatitude = c(rep(35.0, 4L), rep(35.002, 4L), 35.001, 35.003),
    decimalLongitude = c(rep(-80.0, 4L), rep(-80.002, 4L), -80.001, -80.003),
    recordedBy = paste0("Coll", 1:10), recordNumber = as.character(101:110),
    institutionCode = rep("NCU", 10L), phenological_state = rep("anthesis", 10L),
    doy = c(rep(115, 4L), rep(120, 4L), 118, 118),
    can1 = c(-3.0, -2.9, -3.1, -3.0, 3.0, 2.9, 3.1, 3.0, 0.0, 0.05),
    can2 = c(0.1, -0.1, 0.05, -0.05, 0.1, -0.1, 0.05, -0.05, 0.0, 0.02),
    cda_predicted_taxon = c(rep("Packera aurea", 4L), rep("Packera dubia", 4L), "Packera dubia", "Packera dubia"),
    cda_posterior_prob = c(rep(0.95, 8L), 0.51, 0.52),
    cda_posterior_margin = c(rep(0.90, 8L), 0.02, 0.04),
    gmm_uncertainty = c(rep(0.02, 8L), 0.35, 0.38),
    gmm_cluster = c(rep(1L, 8L), 1L, 1L),
    morph_tier_source = rep("Basal_Rosette", 10L),
    basal_leaf_count = rep(3L, 10L),
    stringsAsFactors = FALSE
  )
}
syn_cent <- apply_taxonomic_decision_matrix(mk_syn(TRUE), intraspecific_mode = FALSE, opts = list(demo = TRUE))
syn_real <- apply_taxonomic_decision_matrix(mk_syn(FALSE), intraspecific_mode = FALSE, opts = list(demo = TRUE))

mk_cline <- function(shift_cent) {
  lat_b <- c(35.0, 35.2, 35.4, 35.1, 35.3, 35.0, 35.2, 35.4, 35.1, 35.3)
  lon_b <- c(-80.0, -80.2, -80.4, -80.1, -80.3, -80.5, -80.6, -80.7, -80.8, -80.9)
  lat_c <- if (shift_cent) c(44.0, 45.0, 46.0, 47.0, 48.0, 49.0) else c(35.1, 35.2, 35.3, 35.4, 35.5, 35.6)
  lon_c <- if (shift_cent) c(-72.0, -71.0, -70.0, -69.0, -68.0, -67.0) else rep(-80.2, 6L)
  data.frame(
    catalogNumber = sprintf("CL_%02d", 1:16),
    species_standardized = rep("Packera dubia", 16L),
    determiner_tier = rep("Tier_1_Gold", 16L),
    is_centroid_or_inst_coord = c(rep(FALSE, 10L), rep(TRUE, 6L)),
    decimalLatitude = c(lat_b, lat_c), decimalLongitude = c(lon_b, lon_c),
    recordedBy = paste0("Collector", 1:16), recordNumber = as.character(1:16),
    institutionCode = rep(c("NCU", "DUKE", "US"), length.out = 16L),
    phenological_state = rep("anthesis", 16L),
    doy = c(rep(115, 5L), rep(160, 5L), rep(120, 6L)),
    bio1_temp_mean = c(rep(15.0, 10L), seq(5.0, 15.0, length.out = 6L)),
    can1 = c(rep(0.0, 5L) + c(-0.1, 0.1, 0, -0.05, 0.05), rep(5.0, 5L) + c(-0.1, 0.1, 0, -0.05, 0.05), seq(0.0, 5.0, length.out = 6L)),
    can2 = rep(0.0, 16L),
    gmm_cluster = c(rep(1L, 5L), rep(2L, 5L), rep(2L, 3L), rep(1L, 3L)),
    morph_tier_source = rep("Basal_Rosette", 16L),
    basal_leaf_count = rep(3L, 16L),
    stringsAsFactors = FALSE
  )
}
cl1 <- apply_taxonomic_decision_matrix(mk_cline(FALSE), intraspecific_mode = TRUE, opts = list(demo = TRUE))
cl2 <- apply_taxonomic_decision_matrix(mk_cline(TRUE), intraspecific_mode = TRUE, opts = list(demo = TRUE))

# 3. Active-anthesis gating + LOO circular phenology kappa
df_ph_mix <- data.frame(
  species_standardized = c(rep("Packera A", 5L), rep("Packera B", 5L), "Packera C", rep("Packera A", 4L), rep("Packera B", 4L)),
  phenological_state = c(rep("anthesis", 11L), rep("fruit", 8L)),
  doy = c(105, 107, 109, 111, 113, 210, 212, 214, 216, 218, 160, rep(300, 8L)),
  cda_predicted_taxon = c(rep("Packera A", 5L), rep("Packera B", 5L), "Packera C", rep("Packera A", 4L), rep("Packera B", 4L)),
  PC1 = c(rep(-1.0, 5L), rep(1.0, 5L), 0.0, rep(-1.0, 4L), rep(1.0, 4L)),
  stringsAsFactors = FALSE
)
kw_mix <- compute_modality_kappa_weights(df_ph_mix)
df_ph_frt <- data.frame(
  species_standardized = c(rep("Packera A", 4L), rep("Packera B", 4L)),
  phenological_state = rep("fruit", 8L),
  doy = c(90, 92, 94, 96, 240, 242, 244, 246),
  cda_predicted_taxon = c(rep("Packera A", 4L), rep("Packera B", 4L)),
  PC1 = c(rep(-1.0, 4L), rep(1.0, 4L)),
  stringsAsFactors = FALSE
)
kw_frt <- compute_modality_kappa_weights(df_ph_frt)

# 4a. Weighted Hopkins lm sensitivity: Tier_1_Gold vs Tier_3_Bronze perturbation
mk_hop <- function(b_doy, g_doy_high) {
  data.frame(
    catalogNumber = sprintf("H_%02d", 1:16),
    species_standardized = rep("Packera dubia", 16L),
    determiner_tier = c(rep("Tier_1_Gold", 12L), rep("Tier_3_Bronze", 4L)),
    authority_weight = c(rep(2.0, 12L), rep(0.05, 4L)),
    decimalLatitude = c(seq(32.0, 38.0, length.out = 12L), rep(42.0, 4L)),
    elevation_m = rep(150.0, 16L),
    bio1_temp_mean = rep(14.5, 16L),
    phenological_state = rep("anthesis", 16L),
    doy = c(100 + 4 * (seq(32.0, 38.0, length.out = 11L) - 35.0), g_doy_high, rep(b_doy, 4L)),
    stringsAsFactors = FALSE
  )
}
h_base <- compute_hopkins_adjusted_phenology(mk_hop(128, 112))
h_bpert <- compute_hopkins_adjusted_phenology(mk_hop(240, 112))
h_gpert <- compute_hopkins_adjusted_phenology(mk_hop(128, 224))
diff_bronze <- max(abs(h_base$hopkins_adjusted_doy[1:11] - h_bpert$hopkins_adjusted_doy[1:11]))
diff_gold <- max(abs(h_base$hopkins_adjusted_doy[1:11] - h_gpert$hopkins_adjusted_doy[1:11]))

# 4b. Cytology gigas calibration gate
df_g <- data.frame(
  species_standardized = rep("Packera dubia", 8L),
  area_mm2 = c(100, 110, 105, 115, 250, 260, 255, 270),
  petiole_width_mm = c(1.0, 1.1, 1.0, 1.2, 2.2, 2.3, 2.4, 2.5),
  log1p_capitulum_count = c(1.5, 1.6, 1.5, 1.7, 2.3, 2.4, 2.5, 2.6),
  scale_source = rep("tick_lattice_unboxed", 8L),
  is_calibrated_mm = rep(FALSE, 8L),
  stringsAsFactors = FALSE
)
g_uncal <- .compute_gigas_audit(df_g, c(2L, 2L, 2L, 2L, 4L, 4L, 4L, 4L), df_g$species_standardized)
df_g$scale_source <- "ruler_calibrated"; df_g$is_calibrated_mm <- TRUE
g_cal <- .compute_gigas_audit(df_g, c(2L, 2L, 2L, 2L, 4L, 4L, 4L, 4L), df_g$species_standardized)

# 4c. Deterministic render_step07a_integrative_synthesis_png
syn_fig_df <- syn_real
syn_fig_df$PC1 <- syn_fig_df$can1; syn_fig_df$PC2 <- syn_fig_df$can2
syn_fig_df$w_morph_kappa <- 0.85; syn_fig_df$w_edaphic_kappa <- 0.65; syn_fig_df$w_pheno_kappa <- 0.55
render_step07a_integrative_synthesis_png(syn_fig_df, out_plot_path = png_a)
render_step07a_integrative_synthesis_png(syn_fig_df, out_plot_path = png_b)

write_json(list(
  anc_responses = as.character(anc_1$Response),
  anc_identical_cauline = isTRUE(all.equal(anc_1, anc_2, tolerance = 1e-8)),
  anc_eta = as.numeric(anc_1$partial_eta_sq[anc_1$Response == "PC1"][1L]),
  eta_f = eta_f, eta_s = eta_s,
  gp_tiers = as.list(gp_tiers),
  syn_cent_hyb = any(syn_cent$hybrid_zone_cohesion_pass),
  syn_real_hyb = any(syn_real$hybrid_zone_cohesion_pass),
  cl1_pass = cl1$cryptic_non_clinal_pass[6:10],
  cl2_pass = cl2$cryptic_non_clinal_pass[6:10],
  cl_hop_diff = max(abs(cl1$hopkins_adjusted_doy[1:10] - cl2$hopkins_adjusted_doy[1:10])),
  kw_mix_pheno = kw_mix$w_pheno_raw,
  kw_frt_pheno = kw_frt$w_pheno_raw,
  diff_bronze = diff_bronze, diff_gold = diff_gold,
  g_uncal_traits = I(as.character(g_uncal$trait)),
  g_cal_traits = I(as.character(g_cal$trait))
), out_json, auto_unbox = TRUE, na = "null")
"""


def _iter_r_string_literals(src: str) -> List[str]:
    return [m.group(0)[1:-1] for m in _R_STR_RE.finditer(src)]


class TestStep07AuditHygiene(unittest.TestCase):
    """Static hygiene and single warm-R functional invariants for Step 07."""

    r_res: Dict[str, Any] = {}
    png_bytes: Tuple[bytes, bytes] = (b"", b"")

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            return
        with tempfile.TemporaryDirectory(prefix="step07_audit_r_") as td:
            tmp = Path(td)
            r_p, out_j = tmp / "probe.R", tmp / "out.json"
            pa, pb = tmp / "fig_a.png", tmp / "fig_b.png"
            r_p.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = subprocess.run(
                ["Rscript", str(r_p), str(out_j), str(pa), str(pb)],
                cwd=PROJECT_ROOT, capture_output=True, text=True, check=False, timeout=120,
            )
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.r_res = json.loads(out_j.read_text(encoding="utf-8"))
            ca = pa.parent / "fig07a_integrative_taxonomic_synthesis.png"
            cls.png_bytes = (pa.read_bytes() if pa.exists() else ca.read_bytes(), pb.read_bytes())

    def test_step07_static_hygiene_ceilings_and_posix_regex(self) -> None:
        bad_tokens = ("chisq99_type_morph", "non_anth_cnt", "geno_sNMF_entropy", "bisse_net_div_rate")
        for path, ceil in STEP07_CEILINGS.items():
            self.assertTrue(path.exists(), f"Missing Step 07 module: {path}")
            txt = path.read_text(encoding="utf-8")
            self.assertLess(len(txt.splitlines()), ceil, f"{path.name} exceeded {ceil} lines")
            for lit in _iter_r_string_literals(txt):
                self.assertIsNone(_NON_POSIX_ESC_RE.search(lit), f"Non-POSIX escape in {path.name}: {lit!r}")
                self.assertNotIn("polaris_", lit, f"Legacy polaris_ string in {path.name}: {lit!r}")
            for tok in bad_tokens:
                self.assertNotIn(tok, txt, f"Retired token {tok} found in {path.name}")
            if path.name != "step07_organellar_triage_utils.R":
                self.assertNotIn("interclass_het_Q12", txt, f"Retired token interclass_het_Q12 found in {path.name}")
        for fig_p in (
            PROJECT_ROOT / "scripts" / "figures" / "step07_figures.R",
            PROJECT_ROOT / "scripts" / "figures" / "step07b_network_figures.R",
        ):
            self.assertNotIn("rnorm(", fig_p.read_text(encoding="utf-8"), f"Random jitter in {fig_p.name}")

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_marginal_drop1_ancova_and_exclusions(self) -> None:
        self.assertIn("PC1", self.r_res["anc_responses"])
        self.assertNotIn("blade_cie_l_norm", self.r_res["anc_responses"])
        self.assertTrue(self.r_res["anc_identical_cauline"])
        self.assertAlmostEqual(float(self.r_res["eta_f"]), float(self.r_res["eta_s"]), places=8)
        self.assertAlmostEqual(float(self.r_res["anc_eta"]), float(self.r_res["eta_f"]), places=3)

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_centroid_quarantine_across_propagation_syntopy_clines_hopkins(self) -> None:
        gp = self.r_res["gp_tiers"]
        self.assertNotEqual(gp.get("U_CENT"), "G2_Conditional_Deme_Hull")
        self.assertNotEqual(gp.get("U_DON_C"), "G2_Conditional_Deme_Hull")
        self.assertEqual(gp.get("U_CLEAN"), "G2_Conditional_Deme_Hull")
        self.assertFalse(bool(self.r_res["syn_cent_hyb"]))
        self.assertTrue(bool(self.r_res["syn_real_hyb"]))
        self.assertEqual(self.r_res["cl1_pass"], self.r_res["cl2_pass"])
        self.assertAlmostEqual(float(self.r_res["cl_hop_diff"]), 0.0, places=5)

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_active_anthesis_gating_and_loo_phenology_kappa(self) -> None:
        self.assertIsNotNone(self.r_res["kw_mix_pheno"])
        self.assertGreaterEqual(float(self.r_res["kw_mix_pheno"]), 0.85)
        self.assertIsNone(self.r_res["kw_frt_pheno"])

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_weighted_hopkins_gigas_mm_gate_and_deterministic_fig07a(self) -> None:
        self.assertLess(float(self.r_res["diff_bronze"]), float(self.r_res["diff_gold"]))
        self.assertEqual(self.r_res["g_uncal_traits"], ["log1p_capitulum_count"])
        self.assertEqual(sorted(self.r_res["g_cal_traits"]), ["area_mm2", "log1p_capitulum_count", "petiole_width_mm"])
        ba, bb = self.png_bytes
        self.assertTrue(ba.startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertEqual(hashlib.md5(ba).hexdigest(), hashlib.md5(bb).hexdigest())


if __name__ == "__main__":
    unittest.main()
