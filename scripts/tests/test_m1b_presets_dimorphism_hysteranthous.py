#!/usr/bin/env python3
"""
Milestone M1B integration & unit test suite:
  1. All 6 GrowthFormSpec presets (including woody_branch 35 mm rosette-merge disable
     in organ_crown_metrics.py and non_laminate_bypass planar gate bypass).
  2. Family-agnostic reproductive_phenology_engine.py (generic_lm3 vs. asteraceae_radiate).
  3. Consolidated single warm Rscript subprocess (< 15s) verifying:
     - is_active_genus_binomial() & abbreviate_binomial() trinomial support (var./subsp./nothosubsp./nothovar.)
     - resolve_r_growth_form_preset() & Step 03 cauline_herb / graminoid_linear preset behavior
     - dimorphism_blocking_utils.R + plasticity_mitigation_utils.R preventing false K=2 Mclust split
       on bimodal staminate vs. pistillate Salix leaves while preserving true cryptic K=2 splits
     - pool_hysteranthous_deme_phenology_r() + compute_hopkins_adjusted_phenology() on co-located
       (<= 25 km) vs. distant Salix demes and non-disjoint Packera aurea controls.
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple
from unittest import mock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import resolve_growth_form_preset
from scripts.tests._r_session import rscript_available, run_warm_r_cmd
from scripts.vision.capitulum_phenology_classifier import VOUCHER_ANTHESIS, VOUCHER_FRUIT
from scripts.vision.growth_form_presets import GROWTH_FORM_PRESETS, GrowthFormSpec, get_growth_form_spec
from scripts.vision.organ_crown_metrics import cluster_plant_individuals
from scripts.vision.reproductive_phenology_engine import classify_voucher_reproductive_phenology, resolve_phenology_plugin


@dataclass
class _DummyLeafInst:
    bbox: Tuple[int, int, int, int]
    plant_individual_id: int = -1


class TestM1BPresetsDimorphismHysteranthous(unittest.TestCase):
    """Test suite for M1B growth-form presets, generic phenology, dimorphism blocking, and hysteranthous pooling."""

    def test_growth_form_presets_and_woody_branch_rosette_merge_disable(self) -> None:
        """Verify all 6 GrowthFormSpec presets and woody_branch disabling 35 mm rosette absorption."""
        expected = {
            "rosette_herb",
            "cauline_herb",
            "woody_branch",
            "compound_leaf",
            "graminoid_linear",
            "non_laminate_bypass",
        }
        self.assertEqual(set(GROWTH_FORM_PRESETS.keys()), expected)
        for name in expected:
            spec = get_growth_form_spec(name)
            self.assertIsInstance(spec, GrowthFormSpec)
            self.assertEqual(spec.name, name)

        self.assertEqual(get_growth_form_spec("rosette_herb").consolidate_rosette_mm, 35.0)
        self.assertEqual(get_growth_form_spec("cauline_herb").consolidate_rosette_mm, 20.0)
        self.assertTrue(get_growth_form_spec("cauline_herb").allow_cauline_primary)
        self.assertIsNone(get_growth_form_spec("woody_branch").consolidate_rosette_mm)
        self.assertTrue(get_growth_form_spec("compound_leaf").treat_leaflet_as_leaf_whole)
        self.assertTrue(get_growth_form_spec("graminoid_linear").downweight_efa_outline)
        self.assertTrue(get_growth_form_spec("non_laminate_bypass").bypass_planar_leaf_gate)
        with self.assertRaises(ValueError):
            get_growth_form_spec("unknown_preset_xyz")

        with mock.patch.dict(os.environ, {"PLANTMATICS_GROWTH_FORM_PRESET": "woody_branch"}):
            self.assertEqual(resolve_growth_form_preset(), "woody_branch")
            self.assertEqual(get_growth_form_spec().name, "woody_branch")

        # Two close leaves (cluster 0) + 1 detached leaf 32 mm away (320 px at 10 px/mm, > 0.12*2000 horizontal gap)
        def _mk_leaves():
            return [
                _DummyLeafInst(bbox=(980, 180, 1020, 220)),
                _DummyLeafInst(bbox=(980, 220, 1020, 260)),
                _DummyLeafInst(bbox=(980, 540, 1020, 580)),
            ]

        ros_clustered = cluster_plant_individuals(
            _mk_leaves(), sheet_width=2000, sheet_height=2000, pixels_per_mm=10.0, growth_form_preset="rosette_herb"
        )
        self.assertEqual(len({x.plant_individual_id for x in ros_clustered}), 1)

        woody_clustered = cluster_plant_individuals(
            _mk_leaves(), sheet_width=2000, sheet_height=2000, pixels_per_mm=10.0, growth_form_preset="woody_branch"
        )
        self.assertEqual(len({x.plant_individual_id for x in woody_clustered}), 2)

    def test_reproductive_phenology_engine_generic_vs_asteraceae(self) -> None:
        """Verify generic_lm3 handles non-Asteraceae organs while Asteraceae delegates to capitulum classifier."""
        self.assertEqual(resolve_phenology_plugin("auto", "Asteraceae"), "asteraceae_radiate")
        self.assertEqual(resolve_phenology_plugin("auto", "Salicaceae"), "generic_lm3")

        salix_res = classify_voucher_reproductive_phenology(
            [{"label": "catkin", "confidence": 0.94}, {"label": "leaf_whole", "confidence": 0.99}],
            phenology_plugin="auto",
            target_family="Salicaceae",
            catalog_number="SALIX-01",
        )
        self.assertEqual(salix_res["phenology_plugin"], "generic_lm3")
        self.assertEqual(salix_res["voucher_phenological_state"], VOUCHER_ANTHESIS)
        self.assertTrue(salix_res["is_verified_anthesis"])
        self.assertEqual(salix_res["n_anthesis"], 1)

        quercus_res = classify_voucher_reproductive_phenology(
            [{"label": "acorn", "confidence": 0.91}],
            phenology_plugin="generic_lm3",
            target_family="Fagaceae",
            catalog_number="QUERCUS-01",
        )
        self.assertEqual(quercus_res["voucher_phenological_state"], VOUCHER_FRUIT)
        self.assertFalse(quercus_res["is_verified_anthesis"])

        with mock.patch(
            "scripts.vision.reproductive_phenology_engine.classify_asteraceae_voucher_crops",
            return_value={"voucher_phenological_state": VOUCHER_ANTHESIS, "phenology_plugin": "asteraceae_radiate"},
        ) as mock_ast:
            ast_res = classify_voucher_reproductive_phenology(
                [{"label": "flower_many", "confidence": 0.9}],
                phenology_plugin="auto",
                target_family="Asteraceae",
                catalog_number="PACKERA-01",
            )
        mock_ast.assert_called_once()
        self.assertEqual(ast_res["phenology_plugin"], "asteraceae_radiate")

    def test_r_trinomials_presets_dimorphism_and_hysteranthous_warm_suite(self) -> None:
        """Run all M1B R checks in one single warm Rscript invocation (< 15s)."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            c_dir = tmp_dir / "contours"
            c_dir.mkdir(parents=True, exist_ok=True)
            t = np.linspace(0.0, 2.0 * np.pi, 60, endpoint=False)
            for idx, (cat, a_r, b_r) in enumerate((
                ("CAUL1", 60.0, 22.0), ("CAUL2", 62.0, 24.0),
                ("CAUL3", 80.0, 12.0), ("CAUL4", 82.0, 13.0),
            )):
                pts = np.column_stack([200.0 + a_r * np.cos(t), 200.0 + b_r * np.sin(t)])
                csv_p = c_dir / f"{cat}_leaf1.csv"
                lines = ["catalogNumber,leaf_id,point_index,x,y"] + [
                    f"{cat},1,{i},{x:.4f},{y:.4f}" for i, (x, y) in enumerate(pts)
                ]
                csv_p.write_text("\n".join(lines) + "\n", encoding="utf-8")

            r_code = rf"""
suppressPackageStartupMessages({{
  library(mclust)
}})
source("scripts/core/db_utils.R")
source("scripts/morphometrics/dimorphism_blocking_utils.R")
source("scripts/morphometrics/plasticity_mitigation_utils.R")
source("scripts/morphometrics/03_fourier_extractor.R")
source("scripts/morphometrics/04_gmm_morphotools.R")
source("scripts/analysis/step07_postprocess_guards_utils.R")
source("scripts/analysis/step07_hopkins_consensus_utils.R")

# 1. Trinomial binomials, abbreviations, and growth-form preset resolver
stopifnot(all(is_active_genus_binomial(
  c("Pinus ponderosa", "Pinus ponderosa var. scopulorum", "Pinus ponderosa subsp. ponderosa",
    "Pinus x attenuradiata", "Pinus x attenuradiata nothosubsp. foo", "Pinus ponderosa nothovar. bar"),
  genus = "Pinus"
)))
stopifnot(!any(is_active_genus_binomial(c("Pinus", "Pinus sp.", "Salix nigra", ""), genus = "Pinus")))
stopifnot(identical(
  abbreviate_binomial(c("Pinus ponderosa var. scopulorum", "Packera aurea", NA_character_)),
  c("P. ponderosa var. scopulorum", "P. aurea", NA_character_)
))
Sys.unsetenv("PLANTMATICS_GROWTH_FORM_PRESET")
stopifnot(identical(resolve_r_growth_form_preset(), "rosette_herb"))
Sys.setenv(PLANTMATICS_GROWTH_FORM_PRESET = "woody_branch")
stopifnot(identical(resolve_r_growth_form_preset(), "woody_branch"))
Sys.unsetenv("PLANTMATICS_GROWTH_FORM_PRESET")

# 2. Step 03 cauline_herb retention + graminoid_linear outline harmonic weight cap
tmp_r <- {repr(str(tmp_dir))}
c_r <- {repr(str(c_dir))}
v_df <- data.frame(
  catalogNumber = c("CAUL1", "CAUL2", "CAUL3", "CAUL4"),
  species_standardized = c("Salix nigra", "Salix nigra", "Salix exigua", "Salix exigua"),
  scientificName = c("Salix nigra", "Salix nigra", "Salix exigua", "Salix exigua"),
  determiner_tier = rep("Tier_1_Gold", 4L),
  stringsAsFactors = FALSE
)
v_p <- file.path(tmp_r, "vouchers.csv"); write.csv(v_df, v_p, row.names = FALSE)
m_df <- data.frame(
  catalogNumber = c("CAUL1", "CAUL2", "CAUL3", "CAUL4"),
  plant_individual_id = rep(0L, 4L), leaf_id = rep(1L, 4L),
  efa_qc_pass = rep(TRUE, 4L), qc_flag_reason = rep("PASS", 4L),
  pixels_per_mm = rep(10.0, 4L), is_calibrated_mm = rep(TRUE, 4L),
  scale_source = rep("ruler_calibrated", 4L), solidity = rep(0.95, 4L),
  archetype_score = rep(0.95, 4L), leaf_zone = rep("cauline", 4L),
  dist_from_crown_norm = rep(0.55, 4L), blade_width_mm = c(8.0, 8.2, 2.1, 2.3),
  contour_path = file.path(c_r, paste0(c("CAUL1", "CAUL2", "CAUL3", "CAUL4"), "_leaf1.csv")),
  stringsAsFactors = FALSE
)
m_p <- file.path(tmp_r, "manifest.csv"); write.csv(m_df, m_p, row.names = FALSE)
Sys.setenv(PLANTMATICS_TARGET_GENUS = "Salix", PLANTMATICS_SKIP_FIGURES = "1")
run_s03 <- function(preset, tag) {{
  Sys.setenv(PLANTMATICS_GROWTH_FORM_PRESET = preset)
  out_s <- file.path(tmp_r, paste0("efa_", tag, ".csv"))
  opts <- list(
    input = c_r, contours_dir = c_r, masks_dir = file.path(tmp_r, "masks"),
    manifest = m_p, vouchers = v_p, output = out_s,
    output_individual = file.path(tmp_r, paste0("efa_ind_", tag, ".csv")),
    figures_dir = "", config = "config/config.yaml", harmonics = 6L,
    num_pcs = 5L, symmetrize_efa_harmonics = TRUE, seed = 42L
  )
  run_fourier_extraction(opts)
  read.csv(out_s, stringsAsFactors = FALSE)
}}
res_ros <- run_s03("rosette_herb", "ros")
res_caul <- run_s03("cauline_herb", "caul")
res_gram <- run_s03("graminoid_linear", "gram")
stopifnot(all(res_ros$morph_tier_source == "Cauline_Only_Fallback"), all(res_ros$is_passive))
stopifnot(all(res_caul$morph_tier_source == "Cauline_Only_Fallback"), all(!res_caul$is_passive), all(res_caul$authority_weight > 0))
stopifnot(all(!res_gram$is_passive), all(res_gram$outline_harmonic_weight_cap <= 0.35 + 1e-6))
Sys.unsetenv("PLANTMATICS_GROWTH_FORM_PRESET")

# 3. Dimorphism blocking prevents false K=2 Mclust ICL split on Salix sex while preserving cryptic K=2
set.seed(123L)
mk_salix_df <- function(pc1_vals, sex_vals, sp = "Salix nigra") {{
  n_s <- length(pc1_vals)
  data.frame(
    catalogNumber = paste0("SAL", seq_len(n_s)),
    plant_individual_id = rep(0L, n_s),
    species_standardized = rep(sp, n_s),
    determiner_tier = rep("Tier_1_Gold", n_s),
    morph_tier_source = rep("Basal_Rosette", n_s),
    authority_weight = rep(2.0, n_s),
    spatial_balance_weight = rep(1.0, n_s),
    solidity = rep(0.92, n_s),
    is_dissected = rep(FALSE, n_s),
    sex = sex_vals,
    PC1 = pc1_vals,
    PC2 = stats::rnorm(n_s, 0.0, 0.30),
    PC3 = stats::rnorm(n_s, 0.0, 0.30),
    stringsAsFactors = FALSE
  )
}}
df_dio <- mk_salix_df(
  c(stats::rnorm(12L, -3.0, 0.35), stats::rnorm(12L, 3.0, 0.35)),
  c(rep("staminate", 12L), rep("pistillate", 12L))
)
raw_mc <- mclust::Mclust(as.matrix(df_dio[, c("PC1", "PC2", "PC3")]), G = 1:4, verbose = FALSE)
stopifnot(raw_mc$G >= 2L)
df_dio_res <- residualize_ontogenetic_phenophase_plasticity(df_dio, target_cols = c("PC1", "PC2", "PC3"), intraspecific_mode = TRUE)
gmm_dio <- run_gmm_cluster_analysis(df_dio_res, c("PC1", "PC2", "PC3"), max_k = 4L, intraspecific_mode = TRUE, seed = 42L)
stopifnot(gmm_dio$best_k == 1L)

df_cryp <- mk_salix_df(
  c(stats::rnorm(12L, -3.2, 0.30), stats::rnorm(12L, 3.2, 0.30)),
  rep(c("staminate", "pistillate"), 12L),
  sp = "Salix cryptica"
)
df_cryp_res <- residualize_ontogenetic_phenophase_plasticity(df_cryp, target_cols = c("PC1", "PC2", "PC3"), intraspecific_mode = TRUE)
gmm_cryp <- run_gmm_cluster_analysis(df_cryp_res, c("PC1", "PC2", "PC3"), max_k = 4L, intraspecific_mode = TRUE, seed = 42L)
stopifnot(gmm_cryp$best_k == 2L)

df_phase <- data.frame(
  catalogNumber = paste0("PH", 1:14),
  species_standardized = rep("Viola sororia", 14L),
  occurrenceRemarks = c(rep("chasmogamous spring flower", 7L), rep("cleistogamous late shoot", 7L)),
  stringsAsFactors = FALSE
)
covs_ph <- build_sex_and_shoot_phase_covariates_r(df_phase, seq_len(14L))
stopifnot("phase_cleistogamous_ind" %in% names(covs_ph), abs(mean(covs_ph$phase_cleistogamous_ind)) < 1e-8)

# 4. Hysteranthous deme-level phenology pooling (<= 25 km vs. distant & non-disjoint control)
h_df <- data.frame(
  catalogNumber = c("SD_A1", "SD_A2", "SD_L1", "SD_L2", "SD_DIST", "PA_A1", "PA_A2", "PA_L1"),
  plant_individual_id = rep(0L, 8L),
  species_standardized = c(rep("Salix discolor", 5L), rep("Packera aurea", 3L)),
  determiner_tier = rep("Tier_1_Gold", 8L),
  is_verified_anthesis = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
  PC1 = c(NA_real_, NA_real_, 0.4, 0.5, 0.6, 0.1, 0.2, 0.3),
  doy = c(96, 104, 200, 205, 210, 120, 125, 130),
  decimalLatitude = c(35.91, 35.92, 35.93, 35.94, 42.30, 35.91, 35.92, 35.93),
  decimalLongitude = c(-79.05, -79.06, -79.04, -79.05, -72.10, -79.05, -79.06, -79.04),
  stateProvince = c(rep("North Carolina", 4L), "Massachusetts", rep("North Carolina", 3L)),
  county = c(rep("Orange", 4L), "Hampshire", rep("Orange", 3L)),
  stringsAsFactors = FALSE
)
h_pool <- compute_hopkins_adjusted_phenology(h_df)
stopifnot(identical(h_pool$phenology_hysteranthous_pooled, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)))
stopifnot(all(h_pool$doy[3:4] >= 95 & h_pool$doy[3:4] <= 105))
stopifnot(all(is.finite(h_pool$hopkins_adjusted_doy[3:4])), all(h_pool$is_active_anthesis[3:4]))
cat("M1B_R_WARM_SUITE_OK\n")
"""
            r_path = tmp_dir / "test_m1b_warm.R"
            r_path.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(
                ["Rscript", str(r_path)],
                stage_key="step04",
                cwd=PROJECT_ROOT,
                skip_figures=True,
            )
            self.assertEqual(proc.returncode, 0, f"M1B warm R suite failed:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("M1B_R_WARM_SUITE_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()
