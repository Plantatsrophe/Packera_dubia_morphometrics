"""Warm-R contract tests for Step 07 conflict-only deme and cluster genomic propagation."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path
from typing import Any, Dict

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _run_warm_r_suite() -> Dict[str, Any]:
    r_code = textwrap.dedent(
        r"""
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        suppressPackageStartupMessages({
          source("scripts/analysis/step07_genomic_propagation.R")
          source("scripts/analysis/07_triage_dashboard_synthesis.R")
        })
        bg <- data.frame(
          catalogNumber = c(sprintf("BG_DUB_%02d", 1:6), sprintf("BG_AUR_%02d", 1:6)),
          plant_individual_id = 0L,
          species_standardized = c(rep("Packera dubia", 6L), rep("Packera aurea", 6L)),
          determiner_tier = "Tier_1_Gold",
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          decimalLatitude = c(seq(35.00, 35.05, length.out = 6L), seq(34.50, 34.55, length.out = 6L)),
          decimalLongitude = c(rep(-83.50, 6L), rep(-84.00, 6L)),
          coordinateUncertaintyInMeters = 200,
          is_azonal_lithology = 0L,
          fine_soil_ph = c(rep(6.5, 6L), rep(5.2, 6L)),
          type_edaphic_z_score = 0.3,
          hopkins_adjusted_doy = c(rep(120, 6L), rep(150, 6L)),
          is_active_anthesis = TRUE,
          PC1 = c(1.4, 1.6, 1.3, 1.7, 1.5, 1.5, -1.4, -1.6, -1.3, -1.7, -1.5, -1.5),
          PC2 = c(-0.3, 0.3, -0.2, 0.2, -0.1, 0.1, -0.3, 0.3, -0.2, 0.2, -0.1, 0.1),
          cda_predicted_taxon = c(rep("Packera dubia", 6L), rep("Packera aurea", 6L)),
          cda_posterior_prob = 0.92,
          edaphic_best_fit_taxon = c(rep("Packera dubia", 6L), rep("Packera aurea", 6L)),
          rf_predicted_prob = 0.90,
          geno_linkage_tier = c("G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 4L), "G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 4L)),
          geno_quality_pass = c(TRUE, TRUE, rep(NA, 4L), TRUE, TRUE, rep(NA, 4L)),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c("Packera dubia", "Packera dubia", rep(NA_character_, 4L), "Packera aurea", "Packera aurea", rep(NA_character_, 4L)),
          geno_posterior_prob = c(0.95, 0.94, rep(NA_real_, 4L), 0.95, 0.94, rep(NA_real_, 4L)),
          geno_margin = 0.85,
          voucher_ploidy_level = NA_integer_,
          is_multi_cytotype_taxon = 0L,
          gmm_cluster = 1L,
          coancestry_cluster = "",
          coalescent_gdi = NA_real_,
          stringsAsFactors = FALSE
        )

        # Scenario 1: Multi-cytotype missing ploidy vs explicit 2x/4x conflict
        s1_extra <- bg[1:5, ]
        s1_extra$catalogNumber <- c("S1_ANC_NA", "S1_REC_NA", "S1_ANC_2X", "S1_REC_NA_2", "S1_REC_4X")
        s1_extra$determiner_tier <- c("Tier_1_Gold", "Tier_3_Bronze", "Tier_1_Gold", "Tier_3_Bronze", "Tier_3_Bronze")
        s1_extra$decimalLatitude <- c(35.400, 35.403, 35.800, 35.803, 35.802)
        s1_extra$decimalLongitude <- c(-82.800, -82.800, -82.800, -82.800, -82.800)
        s1_extra$PC1 <- c(1.5, 1.52, 1.5, 1.52, 1.52)
        s1_extra$PC2 <- c(0.0, 0.02, 0.0, 0.02, 0.02)
        s1_extra$geno_linkage_tier <- c("G1_Direct_Voucher", NA_character_, "G1_Direct_Voucher", NA_character_, NA_character_)
        s1_extra$geno_quality_pass <- c(TRUE, NA, TRUE, NA, NA)
        s1_extra$geno_predicted_taxon <- c("Packera dubia", NA_character_, "Packera dubia", NA_character_, NA_character_)
        s1_extra$geno_posterior_prob <- c(0.94, NA_real_, 0.94, NA_real_, NA_real_)
        s1_extra$voucher_ploidy_level <- c(NA_integer_, NA_integer_, 2L, NA_integer_, 4L)
        s1_extra$is_multi_cytotype_taxon <- 1L
        # Move S1_REC_4X to a separate site with its own 2x anchor so Site B remains pure 2x+NA
        s1_anc_c <- s1_extra[3L, ]; s1_anc_c$catalogNumber <- "S1_ANC_C"; s1_anc_c$decimalLatitude <- 36.200
        s1_extra$decimalLatitude[5L] <- 36.202
        df_s1 <- rbind(bg, s1_extra, s1_anc_c)
        res_s1 <- propagate_conditional_deme_priors(df_s1)
        sc_s1 <- compute_probabilistic_consensus_scores(res_s1)
        s1_ratio <- function(rec_id, anc_id) {
          sc_s1$w_geno_kappa[sc_s1$catalogNumber == rec_id] / sc_s1$w_geno_kappa[sc_s1$catalogNumber == anc_id]
        }

        # Scenario 2: Allopatric 4x 500 km away in gmm_cluster==2 vs mixed coancestry/2km conflict
        s2_extra <- bg[1:11, ]
        s2_extra$catalogNumber <- c(
          "S2_LOC_ANC1", "S2_LOC_ANC2", "S2_LOC_REC_2X", "S2_LOC_REC_NA", "S2_DISTANT_4X",
          "S2_CLMIX_2X", "S2_CLMIX_4X", "S2_CLMIX_REC",
          "S2_DM_ANC", "S2_DM_4X", "S2_DM_REC"
        )
        s2_extra$determiner_tier <- c("Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze", "Tier_3_Bronze", "Tier_1_Gold",
                                      "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze",
                                      "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze")
        s2_extra$decimalLatitude <- c(35.500, 35.505, 35.502, 35.503, 40.000, 36.200, 36.350, 36.500, 36.800, 36.804, 36.802)
        s2_extra$decimalLongitude <- c(rep(-83.000, 5L), rep(-82.000, 3L), rep(-81.000, 3L))
        s2_extra$PC1 <- 1.5; s2_extra$PC2 <- 0.0
        s2_extra$geno_linkage_tier <- c("G1_Direct_Voucher", "G1_Direct_Voucher", NA, NA, NA,
                                        "G1_Direct_Voucher", "G1_Direct_Voucher", NA,
                                        "G1_Direct_Voucher", NA, NA)
        s2_extra$geno_quality_pass <- ifelse(is.na(s2_extra$geno_linkage_tier), NA, TRUE)
        s2_extra$geno_predicted_taxon <- ifelse(is.na(s2_extra$geno_linkage_tier), NA_character_, "Packera dubia")
        s2_extra$geno_posterior_prob <- ifelse(is.na(s2_extra$geno_linkage_tier), NA_real_, 0.94)
        s2_extra$voucher_ploidy_level <- c(2L, 2L, 2L, NA_integer_, 4L, 2L, 4L, NA_integer_, 2L, 4L, NA_integer_)
        s2_extra$is_multi_cytotype_taxon <- 1L
        s2_extra$gmm_cluster <- c(2L, 2L, 2L, 2L, 2L, 4L, 4L, 4L, 5L, 5L, 5L)
        s2_extra$coancestry_cluster <- c("CL_LOC", "CL_LOC", "", "", "CL_FAR", "CL_MIX", "CL_MIX", "CL_MIX", "CL_DM", "", "")
        s2_extra$coalescent_gdi <- c(0.82, 0.82, NA, NA, NA, 0.80, 0.80, NA, 0.80, NA, NA)
        df_s2 <- rbind(bg, s2_extra)
        res_s2 <- run_step07_genomic_propagation(df_s2)

        # Scenario 3: Uncertainty & buffer-zone soft decay vs hard vetoes
        s3_extra <- bg[1:8, ]
        s3_extra$catalogNumber <- c("S3_ANC", "S3_UNC_1600", "S3_UNC_NA_FAR", "S3_UNC_2500", "S3_BUF_EDZ", "S3_BUF_DOY", "S3_VETO_AZONAL", "S3_VETO_DOY50")
        s3_extra$determiner_tier <- c("Tier_1_Gold", rep("Tier_3_Bronze", 7L))
        s3_extra$decimalLatitude <- c(35.2000, 35.2072, 35.2144, 35.2045, 35.1955, 35.2000, 35.2045, 35.1960)
        s3_extra$decimalLongitude <- c(-83.000, -83.000, -83.000, -83.000, -83.000, -83.0055, -83.004, -83.004)
        s3_extra$coordinateUncertaintyInMeters <- c(200, 1600, NA_real_, 2500, 200, 200, 200, 200)
        s3_extra$is_azonal_lithology <- c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L)
        s3_extra$fine_soil_ph <- c(6.5, 6.5, 6.5, 6.5, 5.5, 6.5, 6.5, 6.5)
        s3_extra$type_edaphic_z_score <- c(0.3, 0.3, 0.3, 0.3, 2.15, 0.3, 0.3, 0.3)
        s3_extra$hopkins_adjusted_doy <- c(120, 120, 120, 120, 120, 143, 120, 170)
        s3_extra$PC1 <- 1.5; s3_extra$PC2 <- 0.0
        s3_extra$geno_linkage_tier <- c("G1_Direct_Voucher", rep(NA_character_, 7L))
        s3_extra$geno_quality_pass <- c(TRUE, rep(NA, 7L))
        s3_extra$geno_predicted_taxon <- c("Packera dubia", rep(NA_character_, 7L))
        s3_extra$geno_posterior_prob <- c(0.95, rep(NA_real_, 7L))
        res_s3 <- propagate_conditional_deme_priors(rbind(bg, s3_extra))

        # Scenario 4: Unsequenced Bronze/Silver misid vs verified G1/Gold micro-sympatry guard
        s4_extra <- bg[1:9, ]
        s4_extra$catalogNumber <- c(
          "S4_ANC_A", "S4_MISID_BRZ", "S4_MISID_SLV", "S4_REC_A",
          "S4_ANC_B", "S4_VER_G1_AUR", "S4_REC_B",
          "S4_ANC_C", "S4_VER_GLD_AUR"
        )
        s4_rec_c <- bg[1L, ]; s4_rec_c$catalogNumber <- "S4_REC_C"; s4_rec_c$determiner_tier <- "Tier_3_Bronze"
        s4_extra <- rbind(s4_extra, s4_rec_c)
        s4_extra$species_standardized <- c(
          "Packera dubia", "Packera aurea", "Packera aurea", "Packera dubia",
          "Packera dubia", "Packera aurea", "Packera dubia",
          "Packera dubia", "Packera aurea", "Packera dubia"
        )
        s4_extra$determiner_tier <- c(
          "Tier_1_Gold", "Tier_3_Bronze", "Tier_2_Silver", "Tier_3_Bronze",
          "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze",
          "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze"
        )
        s4_extra$decimalLatitude <- c(35.300, 35.3027, 35.2973, 35.301, 36.100, 36.1027, 36.101, 36.900, 36.9027, 36.901)
        s4_extra$decimalLongitude <- c(-83.200, -83.200, -83.200, -83.201, -82.500, -82.500, -82.501, -81.800, -81.800, -81.801)
        s4_extra$PC1 <- c(1.5, -1.5, -1.5, 0.55, 1.5, -1.5, 0.55, 1.5, -1.5, 0.55)
        s4_extra$PC2 <- 0.0
        s4_extra$geno_linkage_tier <- c("G1_Direct_Voucher", NA, NA, NA, "G1_Direct_Voucher", "G1_Direct_Voucher", NA, "G1_Direct_Voucher", NA, NA)
        s4_extra$geno_quality_pass <- ifelse(is.na(s4_extra$geno_linkage_tier), NA, TRUE)
        s4_extra$geno_predicted_taxon <- c("Packera dubia", NA, NA, NA, "Packera dubia", "Packera aurea", NA, "Packera dubia", NA, NA)
        s4_extra$geno_posterior_prob <- ifelse(is.na(s4_extra$geno_linkage_tier), NA_real_, 0.95)
        res_s4 <- propagate_conditional_deme_priors(rbind(bg, s4_extra))

        pick_row <- function(tbl, cid) {
          k <- which(tbl$catalogNumber == cid)
          list(
            tier = if (is.na(tbl$genomic_propagation_tier[k])) NULL else tbl$genomic_propagation_tier[k],
            w = if (is.na(tbl$geno_deme_weight_scale[k])) NULL else as.numeric(tbl$geno_deme_weight_scale[k]),
            gdi = if (is.na(tbl$cluster_propagated_gdi[k])) NULL else as.numeric(tbl$cluster_propagated_gdi[k])
          )
        }
        out <- list(
          s1 = list(
            rec_na = pick_row(res_s1, "S1_REC_NA"),
            rec_na2 = pick_row(res_s1, "S1_REC_NA_2"),
            rec_4x = pick_row(res_s1, "S1_REC_4X"),
            ratio_na = as.numeric(s1_ratio("S1_REC_NA", "S1_ANC_NA")),
            ratio_na2 = as.numeric(s1_ratio("S1_REC_NA_2", "S1_ANC_2X"))
          ),
          s2 = list(
            loc_2x = pick_row(res_s2, "S2_LOC_REC_2X"),
            loc_na = pick_row(res_s2, "S2_LOC_REC_NA"),
            clmix = pick_row(res_s2, "S2_CLMIX_REC"),
            dmmix = pick_row(res_s2, "S2_DM_REC")
          ),
          s3 = list(
            unc_1600 = pick_row(res_s3, "S3_UNC_1600"),
            unc_na_far = pick_row(res_s3, "S3_UNC_NA_FAR"),
            unc_2500 = pick_row(res_s3, "S3_UNC_2500"),
            buf_edz = pick_row(res_s3, "S3_BUF_EDZ"),
            buf_doy = pick_row(res_s3, "S3_BUF_DOY"),
            veto_azonal = pick_row(res_s3, "S3_VETO_AZONAL"),
            veto_doy50 = pick_row(res_s3, "S3_VETO_DOY50")
          ),
          s4 = list(
            rec_a_bronze = pick_row(res_s4, "S4_REC_A"),
            rec_b_g1 = pick_row(res_s4, "S4_REC_B"),
            rec_c_gold = pick_row(res_s4, "S4_REC_C")
          )
        )
        cat("<<<JSON_START>>>\n", jsonlite::toJSON(out, auto_unbox = TRUE, null = "null"), "\n<<<JSON_END>>>\n", sep = "")
        """
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        r_path = Path(tmpdir) / "run_deme_prop_tests.R"
        r_path.write_text(r_code, encoding="utf-8")
        proc = subprocess.run(["Rscript", str(r_path)], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(f"Warm R suite failed (exit={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    raw = proc.stdout.split("<<<JSON_START>>>", 1)[1].split("<<<JSON_END>>>", 1)[0].strip()
    return json.loads(raw)


@unittest.skipIf(SKIP_R_COND, "Skipping warm-R deme propagation suite under fast audit or missing Rscript.")
class TestStep07DemePropagation(unittest.TestCase):
    """Contract tests for Step 07 conflict-only deme and cluster genomic propagation."""

    payload: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _run_warm_r_suite()

    def test_multi_cytotype_missing_ploidy_propagates_with_soft_discount(self) -> None:
        s1 = self.payload["s1"]
        for key, rkey in (("rec_na", "ratio_na"), ("rec_na2", "ratio_na2")):
            self.assertEqual(s1[key]["tier"], "G2_Conditional_Deme_Hull")
            self.assertAlmostEqual(s1[key]["w"], 0.6375, places=4)
            self.assertAlmostEqual(s1[rkey], 0.6375, places=4)
        self.assertIsNone(s1["rec_4x"]["tier"])

    def test_allopatric_polyploid_in_same_gmm_cluster_does_not_block_local_deme_or_gdi(self) -> None:
        s2 = self.payload["s2"]
        for key in ("loc_2x", "loc_na"):
            self.assertEqual(s2[key]["tier"], "G2_Conditional_Deme_Hull")
            self.assertAlmostEqual(s2[key]["gdi"], 0.82, places=4)
        self.assertIsNone(s2["clmix"]["gdi"])
        self.assertNotIn(s2["dmmix"]["tier"], ("G2_Conditional_Deme_Hull", "Cluster_GDI_Governed"))
        self.assertIsNone(s2["dmmix"]["w"])
        self.assertIsNone(s2["dmmix"]["gdi"])

    def test_uncertainty_and_buffer_envelope_soft_decay(self) -> None:
        s3 = self.payload["s3"]
        for key in ("unc_1600", "unc_na_far"):
            self.assertEqual(s3[key]["tier"], "G2_Conditional_Deme_Hull")
            self.assertAlmostEqual(s3[key]["w"], 0.6375, places=4)
        for key in ("buf_edz", "buf_doy"):
            self.assertEqual(s3[key]["tier"], "G2_Conditional_Deme_Hull")
            self.assertAlmostEqual(s3[key]["w"], 0.66, places=4)
        for key in ("unc_2500", "veto_azonal", "veto_doy50"):
            self.assertIsNone(s3[key]["tier"])

    def test_unsequenced_bronze_misid_cannot_poison_micro_sympatry_guard(self) -> None:
        s4 = self.payload["s4"]
        self.assertEqual(s4["rec_a_bronze"]["tier"], "G2_Conditional_Deme_Hull")
        self.assertIsNone(s4["rec_b_g1"]["tier"])
        self.assertIsNone(s4["rec_c_gold"]["tier"])


if __name__ == "__main__":
    unittest.main()
