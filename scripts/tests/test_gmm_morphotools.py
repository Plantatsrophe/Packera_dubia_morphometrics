#!/usr/bin/env python3
"""
===============================================================================
Script: test_gmm_morphotools.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 04 (scripts/morphometrics/04_gmm_morphotools.R):
    - Klingenberg (2016) Ontogenetic Allometry & Slope Homogeneity
    - Unsupervised GMM Clustering & Kass & Raftery (1995) Bayes Factors
    - Full-Rank Canonical Discriminant Analysis (CDA) & Dual Ecoregion Track
    - Parametric Omnibus & Pairwise Hotelling T2 MANOVA with BH FDR
    - Dawid-Skene / Empirical Bayes Authority Weighting & Adaptive Anchor Promotion
===============================================================================
"""

import importlib
import re
import shutil
import subprocess
import tempfile
import unittest

from scripts.tests._r_session import rscript_available, run_warm_r_cmd
import zlib
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from scripts.core.config import PipelineConfig

PROJECT_ROOT = Path(__file__).resolve().parents[2]


class TestGMMAndMorphoTools(unittest.TestCase):
    """Test suite for GMM morphometrics and MorphoTools2 CDA architecture."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.r_script_path = PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R"
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.r_unit_proc = None

        if rscript_available():
            grp_out = cls._class_tmp_path / "manova_group_significance.csv"
            pw_out = cls._class_tmp_path / "manova_pairwise_posthoc.csv"
            r_batch_code = f"""
            source("scripts/morphometrics/pub_figures_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/morphometrics/04_gmm_morphotools.R")

            # Block 1: EPA Level III Ecoregion spatial helper agreement (offline + on-disk county_norm GeoJSON schema)
            local({{
              res <- infer_epa_level3_ecoregion(
                lat = c(35.9132, 34.2257, 32.5404, 31.6035, NA, NA, NA),
                lon = c(-79.0558, -77.9447, -82.9037, -94.6555, NA, NA, NA),
                state_province = c("NC", "NC", "GA", "TX", "North Carolina", "North Carolina", "NC"),
                county = c("", "", "", "", "Orange County", "Wake County", "")
              )
              off_b <- .build_offline_epa_l3_sf_r()
              eco_gj <- file.path("{cls._class_tmp_path.as_posix()}", "test_eco_l3.geojson")
              cnty_gj <- file.path("{cls._class_tmp_path.as_posix()}", "test_us_counties.geojson")
              sf::st_write(off_b$eco_sf, eco_gj, driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
              sf::st_write(off_b$cnty_sf[, c("county_key", "state_name", "county_norm", "centroid_lon", "centroid_lat")], cnty_gj, driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
              res_disk <- infer_epa_level3_ecoregion(
                lat = c(NA, NA),
                lon = c(NA, NA),
                state_province = c("North Carolina", "North Carolina"),
                county = c("Orange County", "Wake County"),
                ecoregion_path = eco_gj,
                county_path = cnty_gj
              )
              stopifnot(identical(res_disk, c("Piedmont", NA_character_)))
              cat("EPA_L3:", paste(ifelse(is.na(res), "NA", res), collapse = "|"), "\\n", sep = "")
            }})

            # Block 2: Direct run_multivariate_manova_audit single-species & collinearity safety
            local({{
              set.seed(17)
              df <- data.frame(
                catalogNumber = sprintf("V%02d", 1:18),
                species_standardized = rep("Packera dubia", 18),
                determiner_tier = rep("Tier_1_Gold", 18),
                morph_tier_source = rep("Basal_Rosette", 18),
                regional_group = c(rep("Piedmont", 10), rep("Blue_Ridge", 2), rep(NA_character_, 6)),
                PC1 = c(rnorm(9, -3.0, 0.3), rnorm(9, 3.0, 0.3)),
                PC2 = c(rnorm(9, 1.5, 0.3), rnorm(9, -1.5, 0.3)),
                stringsAsFactors = FALSE
              )
              df$PC3_collinear <- df$PC1 + df$PC2
              df$PC4_const <- 0.0
              gmm_res <- list(classification = c(rep(1L, 9), rep(2L, 9)))
              res <- run_multivariate_manova_audit(
                df,
                c("PC1", "PC2", "PC3_collinear", "PC4_const"),
                gmm_res,
                output_group_csv = "{grp_out.as_posix()}",
                output_pairwise_csv = "{pw_out.as_posix()}"
              )
              stopifnot(nrow(res$group_table) == 1L)
              stopifnot(res$group_table$grouping_factor[1] == "gmm_cluster")
              cat("MANOVA_DIRECT_OK\\n")
            }})

            # Block 3: Beta-Binomial Empirical Bayes active_passive & continuous_prior_posterior modes
            local({{
              df_t <- data.frame(
                species_standardized = rep("Packera aurea", 8),
                determiner_tier = c("Tier_1_Gold", "Tier_1_Gold", "Tier_2_Silver", "Tier_2_Silver",
                                    "Tier_3_Bronze", "Tier_3_Bronze", "Tier_3_Bronze", "Tier_3_Bronze"),
                is_promoted_bronze = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
                c_error = c(0.02, 0.04, 0.12, 0.14, 0.25, 0.30, 0.55, 0.65),
                stringsAsFactors = FALSE
              )
              w_ap <- compute_empirical_bayes_authority_weights(df_t, mode = "active_passive")
              w_cp <- compute_empirical_bayes_authority_weights(df_t, mode = "continuous_prior_posterior")
              stopifnot(all(w_ap[7:8] == 0.0), abs(mean(w_ap[1:6]) - 1.0) < 1e-5)
              stopifnot(w_ap[1] > w_ap[3], w_ap[3] > w_ap[5], w_ap[5] > 0.0)
              stopifnot(abs(max(w_cp) - 3.0) < 1e-5, min(w_cp) >= 0.5, w_cp[1] > w_cp[3], w_cp[3] > w_cp[5], w_cp[5] > w_cp[7])
              cat("EB_OK\\n")
            }})

            # Block 4: 2-step Dawid-Skene EM determiner shrinkage, no double pruning, and K_full (>3 axes) CDA
            local({{
              set.seed(20261030)
              n_sp <- 10
              df_eb <- data.frame(
                catalogNumber = sprintf("EB_%02d", 1:(2 * n_sp)),
                plant_individual_id = 1L,
                species_standardized = rep(c("Packera aurea", "Packera obovata"), each = n_sp),
                determiner_tier = rep(c(rep("Tier_2_Silver", 8), "Tier_1_Gold", "Tier_3_Bronze"), 2),
                identifiedBy = c(
                  rep("Dr. Reliable Specialist", 4),
                  rep("Noisy Generalist", 4),
                  "T.M. Barkley",
                  "Field Student",
                  rep("Dr. Reliable Specialist", 4),
                  rep("Noisy Generalist", 4),
                  "T.M. Barkley",
                  "Field Student"
                ),
                is_promoted_bronze = FALSE,
                PC1 = c(
                  rnorm(4, -3.0, 0.15),
                  rnorm(4, +3.0, 0.15),
                  -3.1, -2.9,
                  rnorm(4, +3.0, 0.15),
                  rnorm(4, -3.0, 0.15),
                  +3.1, +2.9
                ),
                PC2 = rnorm(2 * n_sp, 0.0, 0.15),
                vision_predicted_label = c(
                  rep("Packera aurea", 4),
                  rep("Packera obovata", 4),
                  "Packera aurea", "Packera aurea",
                  rep("Packera obovata", 4),
                  rep("Packera aurea", 4),
                  "Packera obovata", "Packera obovata"
                ),
                c_error = c(rep(0.05, 4), rep(0.92, 4), 0.03, 0.10, rep(0.05, 4), rep(0.92, 4), 0.03, 0.10),
                stringsAsFactors = FALSE
              )

              w_ap <- compute_empirical_bayes_authority_weights(df_eb, mode = "active_passive")
              stopifnot(length(w_ap) == nrow(df_eb))
              stopifnot(w_ap[10] == 0.0, w_ap[20] == 0.0)
              stopifnot(abs(mean(w_ap[w_ap > 0]) - 1.0) < 1e-5)
              stopifnot(mean(w_ap[df_eb$identifiedBy == "Dr. Reliable Specialist"]) >
                        mean(w_ap[df_eb$identifiedBy == "Noisy Generalist"]))

              w_cp <- compute_empirical_bayes_authority_weights(df_eb, mode = "continuous_prior_posterior")
              stopifnot(all(w_cp >= 0.5 - 1e-9 & w_cp <= 3.0 + 1e-9))
              stopifnot(mean(w_cp[df_eb$identifiedBy == "Dr. Reliable Specialist"]) >
                        mean(w_cp[df_eb$identifiedBy == "Noisy Generalist"]))

              set.seed(104)
              taxa5 <- c("Packera spA", "Packera spB", "Packera spC", "Packera spD", "Packera spE")
              centers5 <- rbind(
                c(-4.0,  0.0,  0.0,  0.0),
                c( 4.0,  0.0,  0.0,  0.0),
                c( 0.0,  4.0,  0.0,  0.0),
                c( 0.0, -2.0,  3.0, -4.5),
                c( 0.0, -2.0,  3.0, +4.5)
              )
              rows5 <- list()
              idx <- 1L
              for (g_i in seq_along(taxa5)) {{
                for (j in 1:14) {{
                  rows5[[idx]] <- data.frame(
                    catalogNumber = sprintf("G%d_%02d", g_i, j),
                    plant_individual_id = 1L,
                    species_standardized = taxa5[g_i],
                    determiner_tier = "Tier_1_Gold",
                    morph_tier_source = "Basal_Rosette",
                    PC1 = centers5[g_i, 1] + rnorm(1, 0, 0.18),
                    PC2 = centers5[g_i, 2] + rnorm(1, 0, 0.18),
                    PC3 = centers5[g_i, 3] + rnorm(1, 0, 0.18),
                    PC4 = centers5[g_i, 4] + rnorm(1, 0, 0.18),
                    stringsAsFactors = FALSE
                  )
                  idx <- idx + 1L
                }}
              }}
              df5 <- dplyr::bind_rows(rows5)
              anc5 <- select_adaptive_reference_anchors(
                df5,
                feature_matrix = as.matrix(df5[, c("PC1", "PC2", "PC3", "PC4")]),
                taxon_col = "species_standardized",
                mahal_q_gold_silver = 0.95
              )
              cda5 <- run_canonical_discriminant_analysis(
                df5,
                feature_cols = c("PC1", "PC2", "PC3", "PC4"),
                target_taxa = taxa5
              )
              stopifnot(identical(as.logical(cda5$is_active), as.logical(anc5$active_mask)))
              stopifnot(all(c("Can1", "Can2", "Can3") %in% colnames(cda5$canonical_scores)))
              stopifnot(cda5$num_canonical_axes == 4L)
              spD_idx <- which(df5$species_standardized == "Packera spD")
              spE_idx <- which(df5$species_standardized == "Packera spE")
              stopifnot(mean(cda5$predicted_taxon[spD_idx] == "Packera spD") >= 0.90)
              stopifnot(mean(cda5$predicted_taxon[spE_idx] == "Packera spE") >= 0.90)
              stopifnot(mean(cda5$posterior_prob[c(spD_idx, spE_idx)]) > 0.85)
              cat("PHASE4_DAWID_SKENE_AND_KFULL_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step04_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = run_warm_r_cmd(
                ["Rscript", str(r_batch_file)],
                stage_key="step04",
                cwd=PROJECT_ROOT,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def setUp(self) -> None:
        self.r_script_path = PROJECT_ROOT / "scripts" / "morphometrics" / "04_gmm_morphotools.R"

    @staticmethod
    def _extract_pdf_text_streams(pdf_path: Path) -> str:
        """Decompress all FlateDecode streams in an R-generated vector PDF for text verification."""
        raw = pdf_path.read_bytes()
        chunks = [raw.decode("latin1", errors="ignore")]
        for m in re.finditer(rb"stream\r?\n(.*?)endstream", raw, flags=re.DOTALL):
            payload = m.group(1).rstrip(b"\r\n")
            for cand in (payload, m.group(1)):
                try:
                    txt = zlib.decompress(cand).decode("latin1", errors="ignore")
                    txt = re.sub(r"\)\s*-?\d+(?:\.\d+)?\s*\(", "", txt)
                    txt = txt.replace(r"\(", "(").replace(r"\)", ")")
                    chunks.append(txt)
                    break
                except zlib.error:
                    continue
        return "\n".join(chunks)

    def test_batched_r_helpers_epa_manova_dawid_skene_and_k_full_cda(self):
        """Verify batched R unit checks for EPA Level III helper, single-species MANOVA collinearity, Dawid-Skene EM, and K_full (>3 axes) CDA."""
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 04 R helper tests failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        out = self.r_unit_proc.stdout

        # 1. EPA Level III Ecoregion spatial helper
        epa_line = [ln for ln in out.splitlines() if ln.startswith("EPA_L3:")][0]
        got_epa = epa_line.split("EPA_L3:", 1)[1].strip().split("|")
        self.assertEqual(
            got_epa,
            [
                "Piedmont",
                "Middle_Atlantic_Coastal_Plain",
                "Southeastern_Plains",
                "South_Central_Plains",
                "Piedmont",
                "NA",
                "NA",
            ],
        )

        # 2. Direct single-species MANOVA factor gating & collinearity safety
        self.assertIn("MANOVA_DIRECT_OK", out)
        grp_out = self._class_tmp_path / "manova_group_significance.csv"
        pw_out = self._class_tmp_path / "manova_pairwise_posthoc.csv"
        self.assertTrue(grp_out.exists() and pw_out.exists())
        g_df = pd.read_csv(grp_out)
        p_df = pd.read_csv(pw_out)
        self.assertEqual(list(g_df["grouping_factor"]), ["gmm_cluster"])
        self.assertEqual(int(g_df["n_features"].iloc[0]), 2)
        self.assertEqual(len(p_df), 1)
        self.assertEqual(p_df["contrast"].iloc[0], "Morphotype_K1 vs Morphotype_K2")
        self.assertGreater(float(p_df["hotelling_t2"].iloc[0]), 20.0)

        # 3. Empirical Bayes & Dawid-Skene + K_full (>3 axes) CDA
        self.assertIn("EB_OK", out)
        self.assertIn("PHASE4_DAWID_SKENE_AND_KFULL_OK", out)

    def test_dual_intraspecific_allometry_morphotype_and_supervised_epa_ecoregion_cda(self):
        """Verify single-taxon (--intraspecific-mode) allometry, partial AR_inv & Cauline_Only_Fallback retention, Track 1 Morphotype_K* CDA, and Track 2 EPA Level III Ecoregion CDA."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(202610)
        rng_ar = np.random.default_rng(777)
        eco_specs = [
            ("Piedmont", 35.9132, -79.0558, "North Carolina", "Orange", -3.2, -1.6),
            ("Middle_Atlantic_Coastal_Plain", 34.2257, -77.9447, "North Carolina", "New Hanover", 3.2, -1.4),
            ("South_Central_Plains", 31.6035, -94.6555, "Texas", "Nacogdoches", 0.2, 3.0),
        ]
        rows = []
        vouch_rows = []
        idx_counter = 0
        for eco_name, lat, lon, st_name, co_name, c1, c2 in eco_specs:
            for j in range(12):
                cat = f"DUAL_{idx_counter:03d}"
                idx_counter += 1
                tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                is_t2_fb = j == 11
                cs = float(24.0 + 1.1 * j)
                log_cs = float(np.log(cs))
                rec = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 0 if is_t2_fb else 3,
                    "morph_tier_source": "Cauline_Only_Fallback" if is_t2_fb else "Basal_Rosette",
                    "species_raw": "Packera anonyma",
                    "scientificName": "Packera anonyma",
                    "determiner_tier": tier,
                    "PC1": float(c1 + rng.normal(0.0, 0.22)),
                    "PC2": float(c2 + rng.normal(0.0, 0.22)),
                    "PC3": float(rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)),
                    "PC5": float(rng.normal(0.0, 0.12)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.9 + 0.25 * (c1 > 0) + rng.normal(0.0, 0.05)),
                    "solidity": float(0.88 + rng.normal(0.0, 0.02)),
                    "archetype_score": 0.94,
                    "AR_inv": np.nan if j % 5 == 0 else float(1.15 + 0.08 * c1 + rng_ar.normal(0.0, 0.015)),
                }
                for k in range(1, 13):
                    rec[f"A{k}"] = float((c1 * 0.06) / k + rng.normal(0.0, 0.003))
                    rec[f"B{k}"] = 0.0
                    rec[f"C{k}"] = 0.0
                    rec[f"D{k}"] = float((c2 * 0.06) / k + rng.normal(0.0, 0.003))
                rows.append(rec)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": "Packera anonyma",
                    "determiner_tier": tier,
                    "institutionCode": "NCU" if j % 2 == 0 else "TEX",
                    "stateProvince": st_name,
                    "county": co_name,
                    "decimalLatitude": lat if j < 10 else np.nan,
                    "decimalLongitude": lon if j < 10 else np.nan,
                    "latitude": lat if j < 10 else np.nan,
                    "longitude": lon if j < 10 else np.nan,
                    "regional_group": eco_name,
                    "doy": 118 + j,
                })

        for u_i in range(4):
            cat = f"DUAL_NA_{u_i:02d}"
            rec_na = dict(rows[u_i])
            rec_na.update({
                "catalogNumber": cat,
                "PC1": float(-2.9 + 0.2 * u_i),
                "PC2": float(-1.5 + 0.1 * u_i),
            })
            rows.append(rec_na)
            vouch_rows.append({
                "catalogNumber": cat,
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold",
                "institutionCode": "NCU",
                "stateProvince": "North Carolina",
                "county": "Wake",
                "decimalLatitude": np.nan,
                "decimalLongitude": np.nan,
                "latitude": np.nan,
                "longitude": np.nan,
                "regional_group": np.nan,
                "doy": 122,
            })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            fig_dir = tmp_p / "figures"
            rep_dir = tmp_p / "reports"
            fig_dir.mkdir(parents=True, exist_ok=True)
            rep_dir.mkdir(parents=True, exist_ok=True)

            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = rep_dir / "gmm_bayes_factors_summary.csv"
            allom_csv = rep_dir / "allometry.csv"
            allom_pdf = fig_dir / "fig04a_foliar_allometry_audit.png"
            eco_sum_csv = rep_dir / "ecoregion_supervised_cda_summary.csv"
            cda_pdf = fig_dir / "fig04b_gmm_cda_morphometrics.png"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd(
                [
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(rep_dir / "batch_audit.csv"),
                    "--allometry-report", str(allom_csv),
                    "--allometry-pairwise", str(rep_dir / "pairwise.csv"),
                    "--allometry-plot", str(allom_pdf),
                    "--focal-taxon", "Packera anonyma",
                    "--intraspecific-mode",
                    "--permutations", "29",
                ],
                stage_key="step04",
                cwd=PROJECT_ROOT,
                skip_figures=False,
            )
            self.assertEqual(proc.returncode, 0, f"Dual/Intraspecific Step 04 failed:\n{proc.stderr}\n{proc.stdout}")

            # 0. Single-taxon allometry audit CSV & PDF + zero dropped Tier_1 / Cauline_Only_Fallback / partial-AR_inv vouchers
            self.assertTrue(allom_csv.exists(), "Allometry audit CSV must be exported in single-taxon mode")
            self.assertTrue(allom_pdf.exists() and allom_pdf.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n")
            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), 40, "All 40 vouchers (including Cauline_Only_Fallback and partial AR_inv) must be retained")
            self.assertFalse(flags_df["gmm_cluster"].isna().any(), "No voucher should have NA gmm_cluster when AR_inv is partially missing")
            self.assertFalse(flags_df["cda_predicted_taxon"].isna().any(), "No voucher should have NA cda_predicted_taxon")

            # 1. Track 1: Unsupervised GMM BIC table evaluates G = 1:eff_k with Two_Delta_BIC vs G = 1
            self.assertTrue(bic_csv.exists())
            bic_df = pd.read_csv(bic_csv)
            self.assertEqual(int(bic_df["K"].iloc[0]), 1)
            self.assertAlmostEqual(float(bic_df["Two_Delta_BIC"].iloc[0]), 0.0, places=5)
            self.assertAlmostEqual(
                float(bic_df["Two_Delta_BIC"].iloc[2]),
                round(float(bic_df["BIC"].iloc[2] - bic_df["BIC"].iloc[0]), 2),
                places=2,
            )

            # 2. Track 1: cda_predicted_taxon is strictly Morphotype_K1..Kn (never regional_group)
            pred_morphs = set(flags_df["cda_predicted_taxon"].dropna().astype(str))
            self.assertGreaterEqual(len(pred_morphs), 2)
            for m_lbl in pred_morphs:
                self.assertRegex(m_lbl, r"^Morphotype_K[0-9]+$")

            # 3. Track 2: Supervised EPA Level III Ecoregion CDA columns in morphometric_flags.csv
            for req_col in ("ecoregion_cda_predicted", "ecoregion_cda_posterior_prob", "ecoregion_can1", "ecoregion_can2"):
                self.assertIn(req_col, flags_df.columns)
                self.assertFalse(flags_df[req_col].isna().any(), f"{req_col} must be populated for active + NA-projected rows")
            pred_ecos = set(flags_df["ecoregion_cda_predicted"].dropna().astype(str))
            self.assertEqual(pred_ecos, {"Piedmont", "Middle_Atlantic_Coastal_Plain", "South_Central_Plains"})
            na_rows = flags_df[flags_df["catalogNumber"].str.startswith("DUAL_NA_")]
            self.assertEqual(len(na_rows), 4)
            self.assertTrue(na_rows["regional_group"].isna().all())
            self.assertTrue((na_rows["ecoregion_cda_predicted"] == "Piedmont").all())

            # 4. Track 2 summary report & fig04c_supervised_ecoregion_cda.{pdf,png}
            self.assertTrue(eco_sum_csv.exists(), "ecoregion_supervised_cda_summary.csv must be generated")
            sum_df = pd.read_csv(eco_sum_csv)
            self.assertEqual(set(sum_df["ecoregion"]), {"Piedmont", "Middle_Atlantic_Coastal_Plain", "South_Central_Plains"})
            for s_col in (
                "n_total_vouchers",
                "n_active_anchors",
                "loocv_accuracy",
                "permanova_r2",
                "permanova_p_value",
                "ari_morphotype_vs_ecoregion",
                "fisher_exact_p_value",
            ):
                self.assertIn(s_col, sum_df.columns)
            self.assertGreater(float(sum_df["loocv_accuracy"].iloc[0]), 0.85)
            self.assertGreater(float(sum_df["permanova_r2"].iloc[0]), 0.50)
            self.assertGreater(float(sum_df["ari_morphotype_vs_ecoregion"].iloc[0]), 0.50)
            self.assertLess(float(sum_df["fisher_exact_p_value"].iloc[0]), 0.01)

            for stem in (
                "fig04a_foliar_allometry_audit",
                "fig04b_gmm_cda_morphometrics",
                "fig04c_ecoregion_and_batch_audit",
            ):
                png_f = fig_dir / f"{stem}.png"
                self.assertTrue(png_f.exists() and png_f.stat().st_size > 10000, f"Missing/empty {png_f}")
            self.assertEqual(list(fig_dir.glob("*.pdf")), [])

    def test_tier3_bronze_invariance_and_multivariate_manova_audit(self):
        """Verify Tier_3_Bronze zero-influence invariance across allometry, GMM, CDA, and MANOVA, plus batch audit R^2 < 0.05 and CDA subtitles."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(314)
        specs = [
            ("Packera anonyma", "Piedmont", 35.9132, -79.0558, "NC", -5.5, -1.5, 11),
            ("Packera dubia", "Middle_Atlantic_Coastal_Plain", 34.2257, -77.9447, "NC", 5.5, -1.5, 11),
            ("Packera paupercula", "South_Central_Plains", 31.6035, -94.6555, "TX", 0.0, 5.0, 11),
        ]
        rows = []
        vouch_rows = []
        for s_idx, (sp, eco_name, lat0, lon0, st_code, c1, c2, n_v) in enumerate(specs):
            for j in range(n_v):
                cat = f"MANOVA_{s_idx}_{j:02d}"
                det_tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 9 else "Tier_3_Bronze")
                cs = float(28.0 + 0.8 * j)
                log_cs = float(np.log(cs))
                pc1 = float(c1 + 0.25 * (log_cs - 3.4) + rng.normal(0.0, 0.25))
                pc2 = float(c2 + rng.normal(0.0, 0.25))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": det_tier,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": float(0.4 * s_idx + rng.normal(0.0, 0.15)),
                    "PC4": float(rng.normal(0.0, 0.12)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(2.0 + 0.3 * s_idx + rng.normal(0.0, 0.05)),
                    "solidity": float(0.89 + 0.02 * s_idx),
                    "archetype_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((-1.25 - 0.06 * pc1 + 0.03 * log_cs) / k + rng.normal(0.0, 0.003))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((0.46 + 0.05 * pc2 + 0.02 * log_cs) / k + rng.normal(0.0, 0.003))
                rows.append(row)
                vouch_rows.append(
                    {
                        "catalogNumber": cat,
                        "species_raw": sp,
                        "determiner_tier": det_tier,
                        "institutionCode": "NCU" if j % 2 == 0 else "BRIT",
                        "stateProvince": st_code,
                        "county": "",
                        "regional_group": eco_name,
                        "decimalLatitude": float(lat0 + rng.normal(0.0, 0.03)),
                        "decimalLongitude": float(lon0 + rng.normal(0.0, 0.03)),
                        "doy": 115 + j,
                    }
                )

        # Add 2 Blue_Ridge Gold vouchers (< 3 threshold -> excluded from regional_group MANOVA) and 2 NA-ecoregion vouchers
        for extra_idx, (rg_val, lat_v, lon_v) in enumerate(
            [
                ("Blue_Ridge", 35.5951, -82.5515),
                ("Blue_Ridge", 35.6100, -82.5400),
                ("", np.nan, np.nan),
                ("", np.nan, np.nan),
            ]
        ):
            cat = f"MANOVA_EXTRA_{extra_idx:02d}"
            row = dict(rows[0])
            row["catalogNumber"] = cat
            row["determiner_tier"] = "Tier_1_Gold"
            row["PC1"] = float(-5.4 + 0.1 * extra_idx)
            row["PC2"] = float(-1.4 - 0.1 * extra_idx)
            rows.append(row)
            vouch_rows.append(
                {
                    "catalogNumber": cat,
                    "species_raw": "Packera anonyma",
                    "determiner_tier": "Tier_1_Gold",
                    "institutionCode": "NCU",
                    "stateProvince": "NC",
                    "county": "",
                    "regional_group": rg_val,
                    "decimalLatitude": lat_v,
                    "decimalLongitude": lon_v,
                    "doy": 120,
                }
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)

            def _run_case(sub_name: str, perturb_bronze: bool = False):
                case_dir = tmp_p / sub_name
                fig_dir = case_dir / "outputs" / "figures"
                rep_dir = case_dir / "outputs" / "reports"
                tbl_dir = case_dir / "data" / "tables"
                for d in (fig_dir, rep_dir, tbl_dir):
                    d.mkdir(parents=True, exist_ok=True)
                df_h = pd.DataFrame(rows).copy()
                df_v = pd.DataFrame(vouch_rows).copy()
                if perturb_bronze:
                    b_mask = df_h["determiner_tier"] == "Tier_3_Bronze"
                    df_h.loc[b_mask, ["PC1", "PC2", "PC3"]] += 45.0
                    df_h.loc[b_mask, "log_centsize"] += 2.5
                    for k in range(1, 13):
                        df_h.loc[b_mask, f"A{k}"] += 0.5
                        df_h.loc[b_mask, f"D{k}"] -= 0.5
                    df_v.loc[df_v["determiner_tier"] == "Tier_3_Bronze", "regional_group"] = "Blue_Ridge"
                harm_csv = tbl_dir / "leaf_efa_harmonics.csv"
                vouch_csv = tbl_dir / "curated_vouchers.csv"
                flags_csv = tbl_dir / "morphometric_flags.csv"
                bic_csv = rep_dir / "gmm_bayes_factors_summary.csv"
                allom_csv = rep_dir / "allometry.csv"
                batch_csv = rep_dir / "batch_audit.csv"
                cda_pdf = fig_dir / "fig04b_gmm_cda_morphometrics.png"
                df_h.to_csv(harm_csv, index=False)
                df_v.to_csv(vouch_csv, index=False)
                proc = run_warm_r_cmd(
                    [
                        "Rscript",
                        str(self.r_script_path),
                        "--vouchers", str(vouch_csv),
                        "--harmonics", str(harm_csv),
                        "--output-flags", str(flags_csv),
                        "--output-report", str(bic_csv),
                        "--output-plot", str(cda_pdf),
                        "--output-batch-audit", str(batch_csv),
                        "--allometry-report", str(allom_csv),
                        "--allometry-pairwise", str(rep_dir / "pairwise.csv"),
                        "--allometry-plot", str(fig_dir / "fig04a_foliar_allometry_audit.png"),
                        "--correct-allometry", "true",
                        "--permutations", "19",
                    ],
                    stage_key="step04",
                    cwd=PROJECT_ROOT,
                    skip_figures=perturb_bronze,
                )
                self.assertEqual(proc.returncode, 0, f"Step 04 run failed ({sub_name}):\n{proc.stderr}\n{proc.stdout}")
                return tbl_dir, rep_dir, fig_dir

            tbl_base, rep_base, fig_base = _run_case("base", perturb_bronze=False)
            flags1 = pd.read_csv(tbl_base / "morphometric_flags.csv")
            bic1 = pd.read_csv(rep_base / "gmm_bayes_factors_summary.csv")
            allom1 = pd.read_csv(rep_base / "allometry.csv")
            batch1 = pd.read_csv(rep_base / "batch_audit.csv")

            # 1. Verify Empirical Bayes active/passive authority weights
            bronze_rows = flags1[flags1["determiner_tier"] == "Tier_3_Bronze"]
            active_rows = flags1[flags1["authority_weight"] > 0.0]
            gold_rows = active_rows[active_rows["determiner_tier"] == "Tier_1_Gold"]
            silver_rows = active_rows[active_rows["determiner_tier"] == "Tier_2_Silver"]
            self.assertTrue((bronze_rows["authority_weight"] == 0.0).all())
            self.assertTrue((bronze_rows["is_passive"] == True).all())
            self.assertAlmostEqual(float(active_rows["authority_weight"].mean()), 1.0, places=4)
            self.assertGreater(float(gold_rows["authority_weight"].mean()), float(silver_rows["authority_weight"].mean()))

            # 2. Verify institutional batch effect report R^2 < 0.05 and allometry ANOVA schema
            for resp_pc in ("PC1", "PC2"):
                b_row = batch1[(batch1["Response"] == resp_pc) & (batch1["Term"] == "institutionCode")]
                self.assertEqual(len(b_row), 1)
                self.assertLess(float(b_row["R_squared"].iloc[0]), 0.05)
                self.assertIn("Pass", str(b_row["Batch_Effect_Status"].iloc[0]))
            for req_term in ("log_centsize", "scientificName", "Residual"):
                self.assertIn(req_term, allom1["Term"].tolist())

            # 3. Verify omnibus & pairwise MANOVA tables and CDA PDF subtitles
            grp_df = pd.read_csv(rep_base / "manova_group_significance.csv")
            pw_df = pd.read_csv(rep_base / "manova_pairwise_posthoc.csv")
            req_grp_cols = {
                "grouping_factor", "pillai_trace", "wilks_lambda", "approx_F",
                "num_df", "den_df", "partial_eta2", "p_value", "q_value",
            }
            self.assertTrue(req_grp_cols.issubset(set(grp_df.columns)))
            self.assertEqual(set(grp_df["grouping_factor"]), {"species_standardized", "gmm_cluster", "regional_group"})
            for _, r in grp_df.iterrows():
                self.assertGreater(float(r["pillai_trace"]), 0.20)
                self.assertGreater(float(r["wilks_lambda"]), 0.0)
                self.assertLess(float(r["wilks_lambda"]), 0.90)
                self.assertGreater(float(r["approx_F"]), 3.0)
                self.assertGreater(float(r["partial_eta2"]), 0.10)
                self.assertLess(float(r["p_value"]), 0.01)
                self.assertGreaterEqual(float(r["q_value"]) + 1e-12, float(r["p_value"]))

            rg_row = grp_df[grp_df["grouping_factor"] == "regional_group"].iloc[0]
            self.assertEqual(int(rg_row["n_groups"]), 3)

            sp_pw = pw_df[pw_df["grouping_factor"] == "species_standardized"]
            rg_pw = pw_df[pw_df["grouping_factor"] == "regional_group"]
            cl_pw = pw_df[pw_df["grouping_factor"] == "gmm_cluster"]
            self.assertEqual(len(sp_pw), 3)
            self.assertEqual(len(rg_pw), 3)
            self.assertGreaterEqual(len(cl_pw), 1)
            for m_lbl in set(cl_pw["group_1"]).union(set(cl_pw["group_2"])):
                self.assertRegex(str(m_lbl), r"^Morphotype_K[0-9]+$")
            self.assertNotIn("Blue_Ridge", set(rg_pw["group_1"]).union(set(rg_pw["group_2"])))

            for _, pr in pw_df.iterrows():
                self.assertGreater(float(pr["hotelling_t2"]), 0.0)
                self.assertAlmostEqual(float(pr["pillai_trace"] + pr["wilks_lambda"]), 1.0, places=4)
                self.assertAlmostEqual(float(pr["partial_eta2"]), float(pr["pillai_trace"]), places=4)
                self.assertGreaterEqual(float(pr["q_value"]) + 1e-12, float(pr["p_value"]))

            self.assertTrue((fig_base / "fig04b_gmm_cda_morphometrics.png").exists())
            self.assertTrue((fig_base / "fig04c_ecoregion_and_batch_audit.png").exists())
            self.assertEqual(list(fig_base.glob("*.pdf")), [])
            sp_row = grp_df[grp_df["grouping_factor"] == "species_standardized"].iloc[0]
            self.assertGreater(float(sp_row["approx_F"]), 0.0)
            self.assertGreater(float(rg_row["approx_F"]), 0.0)

            # 4. Perturb Tier_3_Bronze drastically and assert strict invariance across allometry, GMM BIC, CDA, and MANOVA
            tbl_pert, rep_pert, _ = _run_case("perturbed_bronze", perturb_bronze=True)
            flags2 = pd.read_csv(tbl_pert / "morphometric_flags.csv")
            bic2 = pd.read_csv(rep_pert / "gmm_bayes_factors_summary.csv")
            allom2 = pd.read_csv(rep_pert / "allometry.csv")

            np.testing.assert_allclose(
                allom1.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                allom2.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                atol=1e-6,
                err_msg="Tier_3_Bronze perturbation must have zero influence on reference allometry ANOVA.",
            )
            np.testing.assert_allclose(
                bic1.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                bic2.select_dtypes(include=[np.number]).fillna(0.0).to_numpy(),
                atol=1e-6,
                err_msg="Tier_3_Bronze perturbation must have zero influence on reference GMM BIC table.",
            )
            act1 = flags1[flags1["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])].sort_values("catalogNumber")
            act2 = flags2[flags2["determiner_tier"].isin(["Tier_1_Gold", "Tier_2_Silver"])].sort_values("catalogNumber")
            for col in ["can1", "can2", "cda_posterior_prob", "PC1", "PC2"]:
                np.testing.assert_allclose(
                    act1[col].to_numpy(),
                    act2[col].to_numpy(),
                    atol=1e-6,
                    err_msg=f"Tier_3_Bronze perturbation must not alter reference {col}.",
                )
            pd.testing.assert_frame_equal(grp_df, pd.read_csv(rep_pert / "manova_group_significance.csv"))
            pd.testing.assert_frame_equal(pw_df, pd.read_csv(rep_pert / "manova_pairwise_posthoc.csv"))

    def test_dynamic_non_default_taxa_cda_second_taxon_and_multi_plant_petioles(self):
        """Verify dynamic non-default Packera species discovery, singleton passive projection, cda_second_taxon export, and multi-plant (catalogNumber, plant_individual_id) petiole aggregation."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(2026)
        non_default_specs = [
            ("Packera schweinitziana", "Senecio schweinitzianus", (-4.5, -2.5, 0.5)),
            ("Packera serpenticola", "Packera eburnea", (4.5, -2.0, -0.5)),
            ("Packera millefolium", "Senecio millefolium", (-1.5, 4.0, -2.0)),
            ("Packera aurea", "Senecio aureus", (2.5, 3.5, 2.5)),
        ]
        rows = []
        vouch_rows = []
        for idx, (canon_sp, raw_syn, (c1, c2, c3)) in enumerate(non_default_specs):
            for j in range(10):
                is_multi = idx in (0, 1) and j == 0
                cat = "MULTI_SHEET_01" if is_multi else f"ENA_{idx}_{j:02d}"
                pid = (idx + 1) if is_multi else 0
                tier = "Tier_1_Gold" if j < 6 else ("Tier_2_Silver" if j < 8 else "Tier_3_Bronze")
                cs = float(28.0 + 1.2 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": pid,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": raw_syn if j % 2 == 0 else canon_sp,
                    "scientificName": canon_sp,
                    "determiner_tier": tier,
                    "PC1": float(c1 + rng.normal(0.0, 0.15)),
                    "PC2": float(c2 + rng.normal(0.0, 0.15)),
                    "PC3": float(c3 + rng.normal(0.0, 0.12)),
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.8 + 0.25 * idx + rng.normal(0.0, 0.04)),
                    "solidity": float(0.42 if canon_sp == "Packera millefolium" else 0.89),
                    "archetype_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((c1 * 0.06) / k + rng.normal(0.0, 0.003))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((c2 * 0.06) / k + rng.normal(0.0, 0.003))
                rows.append(row)
                if not any(v["catalogNumber"] == cat for v in vouch_rows):
                    vouch_rows.append({
                        "catalogNumber": cat,
                        "species_raw": row["species_raw"],
                        "determiner_tier": tier,
                        "institutionCode": "NCU" if len(vouch_rows) % 2 == 0 else "NY",
                        "stateProvince": "North Carolina",
                        "latitude": 35.6,
                        "longitude": -82.5,
                        "doy": 132,
                    })

        # Add a singleton rare outgroup row (n = 1) to verify passive projection without singular W
        singleton_row = dict(rows[2])
        singleton_row.update({
            "catalogNumber": "ENA_SINGLETON_01",
            "plant_individual_id": 0,
            "species_raw": "Senecio glabellus",
            "scientificName": "Packera glabella",
            "determiner_tier": "Tier_1_Gold",
            "PC1": 0.0,
            "PC2": 0.0,
            "PC3": 0.0,
        })
        rows.append(singleton_row)
        vouch_rows.append({
            "catalogNumber": "ENA_SINGLETON_01",
            "species_raw": "Senecio glabellus",
            "determiner_tier": "Tier_1_Gold",
            "institutionCode": "NCU",
            "stateProvince": "North Carolina",
            "latitude": 35.6,
            "longitude": -82.5,
            "doy": 132,
        })

        man_rows = [
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 1, "leaf_id": 1, "petiole_length_mm": 32.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": round(32.0 / 28.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 1, "leaf_id": 2, "petiole_length_mm": 32.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": round(32.0 / 28.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 2, "leaf_id": 1, "petiole_length_mm": 12.0, "blade_length_mm": 30.0, "petiole_to_blade_ratio": round(12.0 / 30.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 2, "leaf_id": 2, "petiole_length_mm": 12.0, "blade_length_mm": 30.0, "petiole_to_blade_ratio": round(12.0 / 30.0, 4)},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            man_csv = tmp_p / "extracted_leaf_manifest.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            allom_csv = tmp_p / "allometry_significance_audit.csv"
            allom_pdf = tmp_p / "fig04a_foliar_allometry_audit.png"
            cda_pdf = tmp_p / "fig04b_gmm_cda_morphometrics.png"

            df_in = pd.DataFrame(rows)
            df_in.to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)
            pd.DataFrame(man_rows).to_csv(man_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(allom_csv),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
                    "--allometry-plot", str(allom_pdf),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 failed on non-default Packera taxa:\n{proc.stderr}\n{proc.stdout}")

            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), len(df_in), "All vouchers including multi-plant individuals and singleton outgroup must be retained.")

            expected_four = {
                "Packera schweinitziana",
                "Packera serpenticola",
                "Packera millefolium",
                "Packera aurea",
            }
            pred_taxa = set(flags_df["cda_predicted_taxon"].dropna().unique())
            self.assertEqual(pred_taxa, expected_four)
            self.assertTrue(
                pred_taxa.isdisjoint({"Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"}),
                f"Default 4-species names must not appear when absent from dataset: {pred_taxa}",
            )

            # Verify cda_second_taxon export
            self.assertIn("cda_second_taxon", flags_df.columns)
            self.assertTrue(flags_df["cda_second_taxon"].notna().all())
            self.assertTrue((flags_df["cda_second_taxon"] != flags_df["cda_predicted_taxon"]).all())
            self.assertTrue(set(flags_df["cda_second_taxon"].unique()).issubset(expected_four))

            # Verify high classification accuracy on the 4 well-separated species
            core_rows = flags_df[flags_df["species_standardized"].isin(expected_four)]
            acc = float((core_rows["species_standardized"] == core_rows["cda_predicted_taxon"]).mean())
            self.assertGreaterEqual(acc, 0.95, f"Expected >=95% CDA classification accuracy on 4 distinct taxa, got {acc:.3f}")

            # Verify singleton Packera glabella row was projected passively with finite 3D canonical scores
            sing_out = flags_df[flags_df["catalogNumber"] == "ENA_SINGLETON_01"]
            self.assertEqual(len(sing_out), 1)
            self.assertEqual(sing_out["species_standardized"].iloc[0], "Packera glabella")
            self.assertTrue(bool(sing_out["is_passive"].iloc[0]))
            for col in ("can1", "can2", "can3", "cda_posterior_prob"):
                self.assertTrue(np.isfinite(float(sing_out[col].iloc[0])), f"Singleton {col} must be finite.")

            # Verify multi-plant sheet individuals (pid=1, pid=2) are never collapsed and retain distinct petioles
            multi_rows = flags_df[flags_df["catalogNumber"] == "MULTI_SHEET_01"].sort_values("plant_individual_id")
            self.assertEqual(len(multi_rows), 2, "Multi-plant sheet individuals (pid=1, pid=2) must not be collapsed.")
            self.assertAlmostEqual(float(multi_rows.iloc[0]["petiole_length_mm"]), 32.0, places=2)
            self.assertAlmostEqual(float(multi_rows.iloc[1]["petiole_length_mm"]), 12.0, places=2)
            self.assertAlmostEqual(float(multi_rows.iloc[0]["petiole_to_blade_ratio"]), round(32.0 / 28.0, 4), places=3)
            self.assertAlmostEqual(float(multi_rows.iloc[1]["petiole_to_blade_ratio"]), round(12.0 / 30.0, 4), places=3)

    def test_zero_gold_species_adaptive_anchor_promotion_and_empirical_bayes_weights(self):
        """Verify multi-species GMM and CDA succeed via Adaptive Anchor Promotion when species have 0 Tier_1_Gold vouchers."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(202604)
        rows = []
        vouch_rows = []

        species_specs = [
            ("Packera dubia", (-5.0, -1.8, 0.4), ["Tier_1_Gold"] * 5 + ["Tier_2_Silver"] * 4 + ["Tier_3_Bronze"] * 3),
            ("Packera anonyma", (5.0, -1.5, -0.4), ["Tier_2_Silver"] * 3 + ["Tier_3_Bronze"] * 7),
            ("Packera serpenticola", (0.0, 4.8, 0.8), ["Tier_3_Bronze"] * 10),
        ]

        for s_idx, (sp, (c1, c2, c3), tiers) in enumerate(species_specs):
            for j, tier in enumerate(tiers):
                cat = f"ADAPT_{s_idx}_{j:02d}"
                pc1 = float(c1 + rng.normal(0.0, 0.16))
                pc2 = float(c2 + rng.normal(0.0, 0.16))
                pc3 = float(c3 + rng.normal(0.0, 0.12))
                if sp == "Packera dubia" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = 5.05, -1.45, -0.38
                elif sp == "Packera anonyma" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = -18.0, 16.0, 9.0
                elif sp == "Packera serpenticola" and j == len(tiers) - 1:
                    pc1, pc2, pc3 = 19.0, -17.0, -8.0

                cs = float(26.0 + 0.9 * j)
                log_cs = float(np.log(cs))
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": tier,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": pc3,
                    "PC4": float(rng.normal(0.0, 0.10)),
                    "PC5": float(rng.normal(0.0, 0.10)),
                    "foliar_variance": 0.02,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": log_cs,
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.9 + 0.2 * s_idx + rng.normal(0.0, 0.04)),
                    "solidity": 0.90,
                    "archetype_score": 0.94,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((pc1 * 0.06) / k + rng.normal(0.0, 0.002))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((pc2 * 0.06) / k + rng.normal(0.0, 0.002))
                rows.append(row)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "determiner_tier": tier,
                    "institutionCode": "NCU" if j % 2 == 0 else "DUKE",
                    "stateProvince": "NC",
                    "county": "Orange",
                    "decimalLatitude": 35.91,
                    "decimalLongitude": -79.05,
                    "regional_group": "Piedmont",
                    "doy": 115 + j,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            cda_pdf = tmp_p / "fig04b_gmm_cda_morphometrics.png"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(bic_csv),
                    "--output-plot", str(cda_pdf),
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(tmp_p / "allometry.csv"),
                    "--allometry-pairwise", str(tmp_p / "pairwise.csv"),
                    "--allometry-plot", str(tmp_p / "fig04a_foliar_allometry_audit.png"),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 failed with 0-Gold species:\n{proc.stderr}\n{proc.stdout}")

            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), len(rows))
            self.assertEqual(
                set(flags_df["cda_predicted_taxon"].dropna().unique()),
                {"Packera dubia", "Packera anonyma", "Packera serpenticola"},
            )

            anon_df = flags_df[flags_df["species_standardized"] == "Packera anonyma"]
            anon_promo = anon_df[(anon_df["determiner_tier"] == "Tier_3_Bronze") & (anon_df["authority_weight"] > 0.0)]
            self.assertGreaterEqual(len(anon_promo), 3, "Expected clean Bronze anchors promoted for 0-Gold Packera anonyma")
            anon_outlier = flags_df[flags_df["catalogNumber"] == "ADAPT_1_09"].iloc[0]
            self.assertEqual(float(anon_outlier["authority_weight"]), 0.0)
            self.assertTrue(bool(anon_outlier["is_passive"]))

            serp_df = flags_df[flags_df["species_standardized"] == "Packera serpenticola"]
            serp_promo = serp_df[serp_df["authority_weight"] > 0.0]
            self.assertGreaterEqual(len(serp_promo), 5, "Expected clean Bronze anchors promoted for 0-Gold/0-Silver species")
            serp_outlier = flags_df[flags_df["catalogNumber"] == "ADAPT_2_09"].iloc[0]
            self.assertEqual(float(serp_outlier["authority_weight"]), 0.0)
            self.assertTrue(bool(serp_outlier["is_passive"]))

            dub_misid = flags_df[flags_df["catalogNumber"] == "ADAPT_0_11"].iloc[0]
            self.assertEqual(str(dub_misid["cda_predicted_taxon"]), "Packera anonyma")
            self.assertEqual(str(dub_misid["discordance_reason"]), "Packera_dubia_Discordant_to_Packera_anonyma")

    def test_dynamic_subspace_beyond_five_pcs_and_multipc_bh_batch_audit(self):
        """Verify Task 1C dynamic subspace selection (>5 PCs when informative), selected_morph_features export, and multi-feature BH institutional batch audit."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        rng = np.random.default_rng(20260431)
        rows = []
        vouch_rows = []
        # 2 taxa x 24 Gold anchors each (N_min = 24 -> Williams & Titus ceiling = 8 >= 6)
        taxa = ["Packera cymbalaria", "Packera hyperborea"]
        for s_idx, sp in enumerate(taxa):
            for j in range(24):
                cat = f"HI_DIM_{s_idx}_{j:02d}"
                cs = float(18.0 + 0.3 * j)
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": "Tier_1_Gold",
                    "foliar_variance": 0.015,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": float(np.log(cs)),
                    "is_calibrated_mm": True,
                    "aspect_ratio": float(1.55 + rng.normal(0.0, 0.02)),
                    "solidity": 0.91,
                    "archetype_score": 0.95,
                }
                # 6 orthogonal high-variance signal PCs (PC1..PC6) + 2 low-noise PCs (PC7..PC8)
                for p in range(1, 7):
                    val = float(rng.normal(2.5 if (s_idx == 0 and p == 1) else (-2.5 if (s_idx == 1 and p == 1) else 0.0), 4.0))
                    row[f"PC{p}"] = val
                for p in range(7, 9):
                    val = float(rng.normal(0.0, 0.03))
                    row[f"PC{p}"] = val
                for k in range(1, 13):
                    row[f"A{k}"] = float(0.05 / k + rng.normal(0.0, 0.002))
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float(0.04 / k + rng.normal(0.0, 0.002))
                rows.append(row)
                vouch_rows.append({
                    "catalogNumber": cat,
                    "species_raw": sp,
                    "determiner_tier": "Tier_1_Gold",
                    "institutionCode": "NCU" if j % 2 == 0 else "ALA",
                    "stateProvince": "Alaska",
                    "decimalLatitude": 63.4,
                    "decimalLongitude": -148.9,
                    "regional_group": "Alaska_Range",
                    "doy": 210,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            batch_csv = tmp_p / "batch_audit.csv"
            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(tmp_p / "bic.csv"),
                    "--output-plot", str(tmp_p / "fig04b_gmm_cda_morphometrics.png"),
                    "--output-batch-audit", str(batch_csv),
                    "--allometry-report", str(tmp_p / "allom.csv"),
                    "--allometry-pairwise", str(tmp_p / "pw.csv"),
                    "--allometry-plot", str(tmp_p / "fig04a_foliar_allometry_audit.png"),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 dynamic >5 PC test failed:\n{proc.stderr}\n{proc.stdout}")
            flags_df = pd.read_csv(flags_csv)
            batch_df = pd.read_csv(batch_csv)

            self.assertIn("selected_morph_features", flags_df.columns)
            sel_feats = [t for t in str(flags_df["selected_morph_features"].iloc[0]).split(";") if t]
            self.assertGreaterEqual(len(sel_feats), 6, f"Expected >= 6 dynamically retained morphometric features without max_p=5 cap, got {sel_feats}")

            inst_rows = batch_df[batch_df["Term"] == "institutionCode"]
            self.assertGreaterEqual(len(inst_rows), len(sel_feats))
            self.assertIn("p_value_bh", batch_df.columns)
            self.assertTrue(inst_rows["p_value_bh"].notna().all())

    def test_dissection_collinearity_pruning_lobing_pc_rescue_bimodal_anchors_and_fast_permanova(self):
        """Verify |r|<0.75 scalar pruning, PC3+ lobing rescue, bimodal heteroblastic anchor preservation, and fast PERMANOVA exactness."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        r_code = """
        options(rgl.useNULL = TRUE)
        source("scripts/morphometrics/04_gmm_morphotools.R")
        set.seed(123)

        # 1. |r| < 0.75 collinearity pruning on active anchors
        n <- 60
        s <- runif(n, 0.35, 0.95)
        df_sc <- data.frame(
          solidity_std = as.numeric(scale(s)),
          q25_solidity_std = as.numeric(scale(s + rnorm(n, 0, 0.02))),
          dissection_index_std = as.numeric(scale(1 / sqrt(s) + rnorm(n, 0, 0.35))),
          within_rosette_efa_dispersion_std = as.numeric(scale(rnorm(n)))
        )
        kept_sc <- prune_collinear_dissection_scalars(
          df_sc,
          base_cols = character(0),
          cand_scalars = c("solidity_std", "q25_solidity_std", "dissection_index_std", "within_rosette_efa_dispersion_std"),
          active_mask = rep(TRUE, n),
          group_vec = rep(c("SpA", "SpB"), each = 30),
          max_r = 0.75
        )
        stopifnot(!all(c("solidity_std", "q25_solidity_std") %in% kept_sc))
        stopifnot("within_rosette_efa_dispersion_std" %in% kept_sc)

        # 2. Subspace rescue of lobing PC4 within floor(N_min / 3)
        n_sub <- 72
        tx <- rep(c("Packera anonyma", "Packera dubia", "Packera obovata"), each = 24)
        lobing_sig <- ifelse(tx == "Packera anonyma", 1.8, -0.9) + rnorm(n_sub, 0, 0.35)
        pc_mat <- cbind(
          PC1 = rnorm(n_sub, 0, 4.0),
          PC2 = rnorm(n_sub, 0, 3.0),
          PC3 = rnorm(n_sub, 0, 0.55),
          PC4 = 0.35 * lobing_sig,
          PC5 = rnorm(n_sub, 0, 0.25)
        )
        aux_df <- data.frame(
          dissection_index = 1.6 + 0.4 * lobing_sig,
          solidity = 0.75 - 0.12 * lobing_sig
        )
        sel_pcs <- select_morphometric_subspace(
          pc_mat,
          group_vec = tx,
          aux_df = aux_df
        )
        stopifnot("PC4" %in% sel_pcs)
        stopifnot(length(sel_pcs) <= floor(24 / 3))

        # 3. Heteroblasty-aware bimodal Gold+Silver anchor preservation
        df_bi <- rbind(
          data.frame(
            catalogNumber = paste0("E", 1:26), plant_individual_id = 0L,
            species_standardized = "Packera paupercula", determiner_tier = "Tier_1_Gold",
            morph_tier_source = "Basal_Rosette", solidity = rnorm(26, 0.90, 0.02),
            q25_solidity = rnorm(26, 0.88, 0.02),
            PC1 = rnorm(26, -1.4, 0.2), PC2 = rnorm(26, 0.3, 0.2), PC3 = rnorm(26, -0.4, 0.15)
          ),
          data.frame(
            catalogNumber = paste0("D", 1:6), plant_individual_id = 0L,
            species_standardized = "Packera paupercula", determiner_tier = "Tier_1_Gold",
            morph_tier_source = "Basal_Rosette", solidity = rnorm(6, 0.46, 0.03),
            q25_solidity = rnorm(6, 0.42, 0.03),
            PC1 = rnorm(6, 2.9, 0.2), PC2 = rnorm(6, -1.5, 0.2), PC3 = rnorm(6, 3.1, 0.2)
          )
        )
        anc_res <- select_adaptive_reference_anchors(
          df_bi,
          feature_matrix = as.matrix(df_bi[, c("PC1", "PC2", "PC3")]),
          taxon_col = "species_standardized"
        )
        stopifnot(sum(anc_res$active_mask[27:32]) >= 5L)

        # 4. Exact equivalence of O(N*P) fast_euclidean_permanova_trace vs vegan::adonis2
        Xp <- matrix(rnorm(90 * 4), nrow = 90, ncol = 4)
        gp <- factor(rep(c("E1", "E2", "E3"), each = 30))
        Xp[31:60, 1] <- Xp[31:60, 1] + 1.4
        fp <- fast_euclidean_permanova_trace(Xp, gp, permutations = 99L, seed = 7L)
        vp <- vegan::adonis2(Xp ~ gp, permutations = 99L, method = "euclidean")
        stopifnot(abs(fp$R2 - vp$R2[1]) < 1e-8)
        stopifnot(abs(fp$F - vp$F[1]) < 1e-8)
        cat("STEP04_DISSECTION_HELPERS_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step04.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Step 04 dissection R helper test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP04_DISSECTION_HELPERS_OK", proc.stdout)

    def test_nomenclatural_type_mahalanobis_immunity_eb_pinning_and_typical_morphotype_control(self):
        """Verify select_adaptive_reference_anchors() never prunes an intact Holotype/Lectotype during 95% chi^2 FastMCD screening, pins EB weights, and exports nearest_type_morph_taxon and is_typical_morphotype."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        r_code = """
        options(rgl.useNULL = TRUE)
        source("scripts/morphometrics/04_gmm_morphotools.R")
        set.seed(202611)

        # 1. Verify .identify_type_specimen_vec() ranks primary (2L), secondary (1L), and non-types (0L)
        df_t <- data.frame(
          type_status = c("Holotype", "Isolectotype", "Paratype", "Syntype", "Not a type", NA_character_),
          stringsAsFactors = FALSE
        )
        t_id <- .identify_type_specimen_vec(df_t)
        stopifnot(identical(as.integer(t_id$type_rank), c(2L, 2L, 1L, 1L, 0L, 0L)))
        stopifnot(identical(as.logical(t_id$is_type_specimen), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)))

        # 2. Verify select_adaptive_reference_anchors() never prunes an extreme Holotype during 95% chi^2 MCD screening,
        #    while pruning an equally extreme non-type Tier_1_Gold outlier
        df_anc <- data.frame(
          catalogNumber = c("HOLO_EXTREME", "GOLD_EXTREME_OUTLIER", paste0("GOLD_CORE_", 1:18), paste0("OBO_GOLD_", 1:12)),
          plant_individual_id = 0L,
          species_standardized = c(rep("Packera dubia", 20), rep("Packera obovata", 12)),
          determiner_tier = "Tier_1_Gold",
          type_status = c("Holotype", "None", rep("None", 18), "Lectotype", rep("None", 11)),
          morph_tier_source = "Basal_Rosette",
          PC1 = c(-4.8, 14.5, rnorm(18, -1.8, 0.18), 4.2, rnorm(11, 3.8, 0.18)),
          PC2 = c(3.9, -12.0, rnorm(18, 0.4, 0.18), -1.1, rnorm(11, -0.9, 0.18)),
          PC3 = c(2.6, 9.5, rnorm(18, 0.1, 0.15), -0.4, rnorm(11, -0.3, 0.15)),
          stringsAsFactors = FALSE
        )
        X_anc <- as.matrix(df_anc[, c("PC1", "PC2", "PC3")])
        anc_out <- select_adaptive_reference_anchors(df_anc, feature_matrix = X_anc, taxon_col = "species_standardized")
        # Holotype at row 1 MUST remain active and marked Tier_0_Type_Anchor; non-type extreme outlier at row 2 MUST be pruned
        stopifnot(isTRUE(anc_out$active_mask[1]))
        stopifnot(identical(anc_out$anchor_tier_source[1], "Tier_0_Type_Anchor"))
        stopifnot(!isTRUE(anc_out$active_mask[2]))

        # 3. Verify compute_empirical_bayes_authority_weights() pins primary types to 3.0 in continuous mode,
        #    gives type anchors >= 1.4x higher active weight than ordinary Gold anchors, and supports target_modality="morph"
        w_act <- compute_empirical_bayes_authority_weights(
          df_anc, mode = "active_passive", active_mask = anc_out$active_mask, target_modality = "morph"
        )
        w_cont <- compute_empirical_bayes_authority_weights(
          df_anc, mode = "continuous_prior_posterior", active_mask = anc_out$active_mask, target_modality = "morph"
        )
        stopifnot(w_act[1] > 1.35 * stats::median(w_act[3:20]))
        stopifnot(abs(w_cont[1] - 3.0) < 1e-6)

        # 3b. Verify Kish's N_eff in compute_ledoit_wolf_cov(..., weights = w) and two-sample scaling in compute_type_specimen_morph_control()
        X_lw <- scale(X_anc[anc_out$active_mask, , drop = FALSE], center = TRUE, scale = FALSE)
        S_emp_lw <- stats::cov(X_lw)
        lw_unw <- compute_ledoit_wolf_cov(S_emp_lw, X_centered = X_lw, n_obs = nrow(X_lw))
        w_skew <- c(100.0, rep(0.05, nrow(X_lw) - 1L))
        lw_skew <- compute_ledoit_wolf_cov(S_emp_lw, X_centered = X_lw, n_obs = nrow(X_lw), weights = w_skew)
        stopifnot(attr(lw_skew, "n_eff") < attr(lw_unw, "n_eff"))
        stopifnot(attr(lw_skew, "shrinkage_lambda") >= attr(lw_unw, "shrinkage_lambda"))

        # Verify two-sample (1 + 1/n_eff_type)^(-1) scaling in compute_type_specimen_morph_control()
        df_ts_1 <- df_anc[c(1, 3:10, 21:28), ]
        X_ts_1 <- as.matrix(df_ts_1[, c("PC1", "PC2", "PC3")])
        tc_1 <- compute_type_specimen_morph_control(df_ts_1, X_ts_1, S_reg = diag(3), taxon_col = "species_standardized")
        raw_d2_row2 <- sum((X_ts_1[2, ] - X_ts_1[1, ])^2)
        # Single Holotype (n_eff_type = 1) -> scaled by 1 / (1 + 1/1) = 0.5
        stopifnot(abs(tc_1$type_morph_mahal_d2[2] - round(0.5 * raw_d2_row2, 4)) < 1e-4)
        stopifnot(abs(tc_1$type_morph_mahal_d2[1]) < 1e-8)

        # 4. Verify audit_misidentifications() exports type_morph_mahal_d2, nearest_type_morph_taxon,
        #    type_morph_concordance, is_typical_morphotype, and protects type specimens from routine misidentification
        df_anc$gmm_cluster <- c(2L, 3L, rep(2L, 18), rep(1L, 12))
        df_anc$gmm_uncertainty <- 0.04
        cda_res <- run_canonical_discriminant_analysis(
          df_anc, feature_cols = c("PC1", "PC2", "PC3"), target_taxa = c("Packera dubia", "Packera obovata")
        )
        aud_df <- audit_misidentifications(
          df_anc, cda_res, list(classification = df_anc$gmm_cluster, uncertainty = df_anc$gmm_uncertainty)
        )
        req_cols <- c("is_type_specimen", "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance", "is_typical_morphotype", "has_empirical_type_anchor")
        stopifnot(all(req_cols %in% names(aud_df)))
        # P. dubia Holotype is in gmm_cluster 2 -> cluster 2 is typical morphotype for P. dubia
        stopifnot(isTRUE(aud_df$is_typical_morphotype[1]))
        stopifnot(!isTRUE(aud_df$is_typical_morphotype[2]))
        stopifnot(!isTRUE(aud_df$morphometric_outlier_flag[1]))
        stopifnot(identical(aud_df$nearest_type_morph_taxon[21], "Packera obovata"))
        cat("STEP04_TYPE_CONTROL_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step04.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Step 04 type control R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP04_TYPE_CONTROL_OK", proc.stdout)

    def test_morphometric_flags_exports_osr_and_confident_learning_columns(self):
        """
        Verify that `morphometric_flags.csv` exports `min_congener_mahal_d2`,
        `min_congener_taxon`, `open_set_novelty_p`, and `is_confident_label_noise`.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        rng = np.random.default_rng(20261108)
        rows = []
        vouch_rows = []
        taxa = ["Packera dubia", "Packera obovata"]
        for s_idx, sp in enumerate(taxa):
            center_pc1 = -3.0 if s_idx == 0 else 3.0
            center_pc2 = 1.2 if s_idx == 0 else -1.2
            for j in range(12):
                cat = f"OSR_CL_{s_idx}_{j:02d}"
                is_swapped_noise = (j == 11)
                pc1 = float(rng.normal(-center_pc1 if is_swapped_noise else center_pc1, 0.25))
                pc2 = float(rng.normal(-center_pc2 if is_swapped_noise else center_pc2, 0.25))
                pc3 = float(rng.normal(0.0, 0.20))
                cs = float(20.0 + 0.4 * j)
                tier = "Tier_1_Gold" if j < 10 else ("Tier_2_Silver" if j == 11 else "Tier_3_Bronze")
                t_status = "Holotype" if j == 0 else "None"
                row = {
                    "catalogNumber": cat,
                    "plant_individual_id": 0,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "species_raw": sp,
                    "scientificName": sp,
                    "determiner_tier": tier,
                    "type_status": t_status,
                    "foliar_variance": 0.015,
                    "leaf_count": 3,
                    "centroid_size": cs,
                    "log_centsize": float(np.log(cs)),
                    "is_calibrated_mm": True,
                    "aspect_ratio": 1.55,
                    "solidity": 0.90,
                    "archetype_score": 0.94,
                    "PC1": pc1,
                    "PC2": pc2,
                    "PC3": pc3,
                }
                for k in range(1, 13):
                    row[f"A{k}"] = float((pc1 * 0.04) / k)
                    row[f"B{k}"] = 0.0
                    row[f"C{k}"] = 0.0
                    row[f"D{k}"] = float((pc2 * 0.04) / k)
                rows.append(row)
                vouch_rows.append(
                    {
                        "catalogNumber": cat,
                        "species_raw": sp,
                        "scientificName": sp,
                        "determiner_tier": tier,
                        "type_status": t_status,
                        "identifiedBy": "D. K. Trock" if tier == "Tier_1_Gold" else "A. Curator",
                        "institutionCode": "NCU" if j % 2 == 0 else "DUKE",
                        "stateProvince": "NC",
                        "county": "Orange",
                        "decimalLatitude": 35.91 + 0.01 * j,
                        "decimalLongitude": -79.05 - 0.01 * j,
                        "regional_group": "Piedmont",
                        "doy": 115 + j,
                    }
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = run_warm_r_cmd([
                    "Rscript",
                    str(self.r_script_path),
                    "--vouchers", str(vouch_csv),
                    "--harmonics", str(harm_csv),
                    "--output-flags", str(flags_csv),
                    "--output-report", str(tmp_p / "bic.csv"),
                    "--output-plot", str(tmp_p / "fig04b_gmm_cda_morphometrics.png"),
                    "--output-batch-audit", str(tmp_p / "batch_audit.csv"),
                    "--allometry-report", str(tmp_p / "allom.csv"),
                    "--allometry-pairwise", str(tmp_p / "pw.csv"),
                    "--allometry-plot", str(tmp_p / "fig04a_foliar_allometry_audit.png"),
                    "--permutations", "19",
                ], stage_key='step04', cwd=PROJECT_ROOT, skip_figures=True)
            self.assertEqual(proc.returncode, 0, f"Step 04 OSR/Confident Learning test failed:\n{proc.stderr}\n{proc.stdout}")
            flags_df = pd.read_csv(flags_csv)
            for req_col in (
                "min_congener_mahal_d2",
                "min_congener_taxon",
                "open_set_novelty_p",
                "is_confident_label_noise",
            ):
                self.assertIn(req_col, flags_df.columns, f"Missing {req_col} in morphometric_flags.csv")
            self.assertTrue(flags_df["min_congener_mahal_d2"].notna().all())
            self.assertTrue(flags_df["min_congener_taxon"].notna().all())
            self.assertTrue(flags_df["open_set_novelty_p"].notna().all())
            self.assertTrue(flags_df["is_confident_label_noise"].notna().all())

    def test_spatial_gathering_demes_and_per_species_typical_morphotype(self):
        """
        Verify that:
        1. `assign_spatial_gathering_demes()` groups conspecific individuals within
           a 10 km Haversine deme without pre-shrinking feature coordinates (and
           dead `apply_hierarchical_population_deme_shrinkage` is purged).
        2. `compute_type_specimen_morph_control()` evaluates `is_typical_morphotype`
           per-species rather than requiring global `gmm_cluster == 1L`.
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        source("scripts/morphometrics/authority_weighting_utils.R")
        source("scripts/morphometrics/plasticity_mitigation_utils.R")
        stopifnot(!exists("apply_hierarchical_population_deme_shrinkage", mode = "function"))
        stopifnot(exists("assign_spatial_gathering_demes", mode = "function"))

        # 1. Test 10 km conspecific deme assignment without pre-modeling coordinate shrinkage
        df_deme <- data.frame(
          catalogNumber = c("SHEET-001", "SHEET-001", "SHEET-002", "SHEET-003", "SHEET-ISOLATE"),
          plant_individual_id = c(0L, 1L, 0L, 0L, 0L),
          species_standardized = rep("Packera dubia", 5L),
          determiner_tier = rep("Tier_1_Gold", 5L),
          authority_weight = rep(0.95, 5L),
          basal_leaf_count = c(2L, 2L, 4L, 4L, 3L),
          decimalLatitude = c(35.910, 35.910, 35.914, 35.918, 37.500),
          decimalLongitude = c(-79.050, -79.050, -79.054, -79.058, -81.200),
          PC1 = c(-2.60, -1.40, -2.05, -1.95, -2.10),
          PC2 = c(0.60, -0.40, 0.12, 0.08, 0.10),
          stringsAsFactors = FALSE
        )
        deme_df <- assign_spatial_gathering_demes(df_deme, deme_radius_km = 10.0)
        stopifnot("population_deme_id" %in% names(deme_df))
        stopifnot("deme_n_specimens" %in% names(deme_df))
        # First 4 rows are within 10 km of each other and share a deme; raw PC1 coordinates stay untouched
        stopifnot(length(unique(deme_df$population_deme_id[1:4])) == 1L)
        stopifnot(deme_df$population_deme_id[5L] != deme_df$population_deme_id[1L])
        stopifnot(all(deme_df$deme_n_specimens[1:4] == 4L))
        stopifnot(max(abs(deme_df$PC1 - df_deme$PC1)) == 0.0)

        # 2. Test per-species evaluation of is_typical_morphotype in compute_type_specimen_morph_control()
        df_typ <- data.frame(
          catalogNumber = c("DUB-TYPE", "DUB-GOLD-1", "DUB-GOLD-2", "OBO-TYPE", "OBO-GOLD-1", "OBO-GOLD-2"),
          plant_individual_id = rep(0L, 6L),
          species_standardized = c(rep("Packera dubia", 3L), rep("Packera obovata", 3L)),
          determiner_tier = rep("Tier_1_Gold", 6L),
          authority_weight = rep(0.98, 6L),
          type_status = c("Holotype", "None", "None", "Holotype", "None", "None"),
          is_type_specimen = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
          gmm_cluster = c(1L, 1L, 1L, 2L, 2L, 2L),
          PC1 = c(-2.20, -2.15, -2.25, 2.30, 2.25, 2.35),
          PC2 = c(0.10, 0.12, 0.08, -0.10, -0.12, -0.08),
          stringsAsFactors = FALSE
        )
        ctrl <- compute_type_specimen_morph_control(
          df = df_typ,
          X_mat = as.matrix(df_typ[, c("PC1", "PC2")]),
          taxon_col = "species_standardized"
        )
        stopifnot(all(ctrl$is_typical_morphotype[1:3]))
        stopifnot(all(ctrl$is_typical_morphotype[4:6]))
        cat("DEME_ASSIGNMENT_AND_PER_SPECIES_TYPICAL_OK\\n")
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_deme_and_typical.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Deme assignment / typical morphotype R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("DEME_ASSIGNMENT_AND_PER_SPECIES_TYPICAL_OK", proc.stdout)

    def test_heteroblasty_leaf_weights_pure_basal_and_cauline_fallback(self) -> None:
        """
        Verify reconcile_step03_basal_otu_medians (and absence of legacy weighted_otu_stat / compute_heteroblasty_leaf_weights):
          1. When a specimen has >= 1 basal leaf (leaf_zone == "basal"), cauline leaves are strictly
             excluded from primary OTU unweighted basal medians so harmonics and solidity are uncontaminated.
          2. When a specimen has only cauline leaves (0 basal leaves), unweighted cauline medians populate
             valid fallback OTU coordinates.
        """
        r_code = """
        source("scripts/morphometrics/authority_weighting_utils.R")
        stopifnot(!exists("compute_heteroblasty_leaf_weights", mode = "function"))
        stopifnot(!exists("weighted_otu_stat", mode = "function"))
        stopifnot(exists("reconcile_step03_basal_otu_medians", mode = "function"))

        # 1. Mixed basal + cauline specimen -> unweighted basal median (0.91, 1.01) uncontaminated by cauline leaf
        df_mixed <- data.frame(
          catalogNumber = rep("SPEC_MIX", 3L),
          plant_individual_id = rep(0L, 3L),
          leaf_id = 1:3,
          leaf_zone = c("basal", "basal", "cauline"),
          area_mm2 = c(400.0, 100.0, 90.0),
          solidity = c(0.92, 0.90, 0.42),
          dissection_index = c(1.08, 1.10, 2.40),
          blade_length_mm = c(40.0, 38.0, 19.0),
          A1 = c(1.00, 1.02, 0.30),
          D1 = c(0.50, 0.52, 0.95),
          stringsAsFactors = FALSE
        )
        spec_init <- data.frame(
          catalogNumber = "SPEC_MIX",
          plant_individual_id = 0L,
          solidity = 0.0,
          dissection_index = 0.0,
          A1 = 0.0,
          D1 = 0.0,
          stringsAsFactors = FALSE
        )
        spec_agg <- reconcile_step03_basal_otu_medians(
          spec_init,
          df_mixed,
          all_harm_cols = c("A1", "D1"),
          sym_names = c("A1", "D1")
        )
        stopifnot(abs(spec_agg$solidity[1L] - 0.91) < 1e-6)
        stopifnot(abs(spec_agg$A1[1L] - 1.01) < 1e-6)

        # 2. Cauline-only specimen -> unweighted cauline fallback median (0.50, 0.42)
        df_caul_only <- data.frame(
          catalogNumber = rep("SPEC_CAUL", 2L),
          plant_individual_id = rep(0L, 2L),
          leaf_id = 1:2,
          leaf_zone = c("cauline", "cauline"),
          area_mm2 = c(120.0, 100.0),
          solidity = c(0.52, 0.48),
          A1 = c(0.40, 0.44),
          D1 = c(0.80, 0.84),
          stringsAsFactors = FALSE
        )
        spec_caul_init <- data.frame(
          catalogNumber = "SPEC_CAUL",
          plant_individual_id = 0L,
          solidity = 0.0,
          A1 = 0.0,
          D1 = 0.0,
          stringsAsFactors = FALSE
        )
        spec_caul <- reconcile_step03_basal_otu_medians(
          spec_caul_init,
          df_caul_only,
          all_harm_cols = c("A1", "D1"),
          sym_names = c("A1", "D1")
        )
        stopifnot(abs(spec_caul$solidity[1L] - 0.50) < 1e-6)
        stopifnot(abs(spec_caul$A1[1L] - 0.42) < 1e-6)

        # 3. Scale-invariant foliar architecture scalar filtering retains dimensionless traits without ruler-missingness bias
        df_arch <- data.frame(
          PC1 = c(-2.1, -1.8, -1.5, 0.2, 1.4, 1.9, 2.3),
          authority_weight = rep(1.0, 7L),
          petiole_to_blade_ratio = c(0.35, 0.42, 0.50, 0.65, 0.72, 0.85, 0.95),
          solidity = c(0.91, 0.89, 0.87, 0.78, 0.72, 0.64, 0.58),
          petiole_to_blade_ratio_std = as.numeric(scale(c(0.35, 0.42, 0.50, 0.65, 0.72, 0.85, 0.95))),
          solidity_std = as.numeric(scale(c(0.91, 0.89, 0.87, 0.78, 0.72, 0.64, 0.58))),
          stringsAsFactors = FALSE
        )
        kept_arch <- filter_biologically_variable_dissection_scalars(
          df_arch,
          c("petiole_to_blade_ratio_std", "solidity_std"),
          active_mask = (df_arch$authority_weight > 0)
        )
        stopifnot("petiole_to_blade_ratio_std" %in% kept_arch)
        stopifnot("solidity_std" %in% kept_arch)
        cat("HETEROBLASTY_PURE_BASAL_AND_CAULINE_FALLBACK_OK\\n")
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_heteroblasty_weights.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_file)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Heteroblasty leaf weights R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("HETEROBLASTY_PURE_BASAL_AND_CAULINE_FALLBACK_OK", proc.stdout)



    def test_step04_plasticity_mitigation_module_and_sympatric_unimodal_collapse(self) -> None:
        """
        Verify plasticity_mitigation_utils.R contracts and end-to-end R behavior:
          1. Static verification that plasticity_mitigation_utils.R defines phenophase residualization,
             within-deme Q_ST covariance whitening, sPCA global test, and Hennig/Ray-Lindsay ridgeline merging,
             and is wired into 04_gmm_morphotools.R without forcing G >= 2.
          2. When Rscript is available, verify that a skewed, phenophase-confounded sympatric panmictic cloud
             collapses to K = 1 (Morphotype_K1) while a geographically segregated bimodal cohort retains K = 2.
        """
        pm_path = PROJECT_ROOT / "scripts" / "morphometrics" / "plasticity_mitigation_utils.R"
        self.assertTrue(pm_path.is_file(), "plasticity_mitigation_utils.R must exist")
        pm_src = pm_path.read_text(encoding="utf-8")
        for fn_tok in (
            "merge_voucher_phenophase_covariates",
            "residualize_ontogenetic_phenophase_plasticity",
            "merge_plastic_gmm_components_ridgeline",
        ):
            self.assertIn(fn_tok, pm_src, f"Missing {fn_tok} in plasticity_mitigation_utils.R")
        self.assertNotIn("compute_spca_global_test", pm_src)

        r04_src = self.r_script_path.read_text(encoding="utf-8")
        self.assertIn("plasticity_mitigation_utils.R", r04_src)
        self.assertNotIn("kmeans_K2_fallback", r04_src)

        if not rscript_available():
            return

        r_code = """
        source("scripts/morphometrics/pub_figures_utils.R")
        source("scripts/morphometrics/authority_weighting_utils.R")
        source("scripts/morphometrics/plasticity_mitigation_utils.R")
        source("scripts/morphometrics/04_gmm_morphotools.R")
        set.seed(20261109)

        # Case A: Sympatric skewed plastic cloud across 4 shared localities -> must collapse to K = 1 (Morphotype_K1)
        n_sym <- 48L
        site_lats <- rep(c(31.5, 33.8, 35.9, 38.4), length.out = n_sym)
        site_lons <- rep(c(-93.2, -84.1, -79.0, -76.1), length.out = n_sym)
        ecos <- rep(c("South_Central_Plains", "Piedmont", "Piedmont", "Middle_Atlantic_Coastal_Plain"), length.out = n_sym)
        # Skewed log-normal microhabitat plasticity completely randomized within every site
        plast_shock <- rexp(n_sym, rate = 1.0) - 1.0
        df_sym <- data.frame(
          catalogNumber = sprintf("SYM_%03d", seq_len(n_sym)),
          plant_individual_id = 0L,
          species_standardized = "Packera dubia",
          determiner_tier = "Tier_1_Gold",
          morph_tier_source = "Basal_Rosette",
          decimalLatitude = site_lats + rnorm(n_sym, 0, 0.02),
          decimalLongitude = site_lons + rnorm(n_sym, 0, 0.02),
          regional_group = ecos,
          solidity = 0.88,
          is_dissected = FALSE,
          PC1 = 1.6 * plast_shock + rnorm(n_sym, 0, 0.25),
          PC2 = -0.8 * plast_shock + rnorm(n_sym, 0, 0.25),
          PC3 = rnorm(n_sym, 0, 0.20),
          stringsAsFactors = FALSE
        )
        gmm_sym <- run_gmm_cluster_analysis(df_sym, c("PC1", "PC2", "PC3"), max_k = 5L, intraspecific_mode = TRUE)
        cda_sym <- run_canonical_discriminant_analysis(df_sym, c("PC1", "PC2", "PC3"), "Packera dubia", gmm_res = gmm_sym, intraspecific_mode = TRUE)
        stopifnot(gmm_sym$best_k == 1L)
        stopifnot(all(cda_sym$predicted_taxon == "Morphotype_K1"))

        # Case B: Geographically segregated true morphotypes -> must retain K = 2
        df_geo <- df_sym
        west_m <- df_geo$decimalLongitude < -88.0
        df_geo$PC1 <- ifelse(west_m, -3.2, 2.2) + rnorm(n_sym, 0, 0.30)
        df_geo$PC2 <- ifelse(west_m, 1.8, -1.2) + rnorm(n_sym, 0, 0.30)
        gmm_geo <- run_gmm_cluster_analysis(df_geo, c("PC1", "PC2", "PC3"), max_k = 5L, intraspecific_mode = TRUE)
        stopifnot(gmm_geo$best_k >= 2L)

        # Case C: Purely morphological ridgeline valley check preserves modality independence
        rl_unimodal <- .compute_ridgeline_valley_ratio(0.5, c(0, 0), diag(2), 0.5, c(0.6, 0.0), diag(2))
        rl_bimodal  <- .compute_ridgeline_valley_ratio(0.5, c(-3, 0), diag(2) * 0.25, 0.5, c(3, 0), diag(2) * 0.25)
        stopifnot(rl_unimodal$n_modes == 1L || rl_unimodal$ridge_ratio >= 0.45)
        stopifnot(rl_bimodal$n_modes >= 2L && rl_bimodal$ridge_ratio < 0.10)
        cat("PLASTICITY_MITIGATION_R_OK\n")
        """
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step04_plast.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(["Rscript", str(r_tmp)], stage_key='step04', cwd=PROJECT_ROOT)
        self.assertEqual(proc.returncode, 0, f"Plasticity mitigation R test failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("PLASTICITY_MITIGATION_R_OK", proc.stdout)


    def test_step04_zero_trust_multi_island_cauline_osr_and_no_fake_loocv(self) -> None:
        """
        Zero-trust regression suite for Step 04 fixes:
          1. Pure-R CLI fallback (.parse_step04_cli_fallback) parses all flags without optparse.
          2. Partial-NA reproductive & petiole trait backfilling populates NA rows even when other rows are finite.
          3. Whole + multi-island dissected vouchers (PC1 = NA, finite ECT_PC1..3 + solidity + dissection_index)
             are retained and classified accurately via shared topology + architecture space, and cauline
             heteroblasty features enter Step 04 at >=50% anchor coverage.
          4. Single-group (g = 1) CDA sets loocv_accuracy = NA (not 1.0), and post-GMM intraspecific figures
             never print fabricated "LOOCV = 85.0%".
          5. Cross-run SQLite open-set recognition uses run-invariant harmonics + topology (immune to PC sign flips).
        """
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            r_script = tmp_p / "zero_trust_step04.R"
            r_script.write_text(
                f"""
                options(rgl.useNULL = TRUE)
                source("scripts/morphometrics/04_gmm_morphotools.R")

                # 1. Pure-R CLI fallback
                defs <- list(vouchers = "a.csv", max_k = 8L, seed = 42L, outline_only = FALSE, intraspecific_mode = FALSE, correct_allometry = "auto")
                parsed <- .parse_step04_cli_fallback(c("--vouchers=b.csv", "--max-k", "6", "--seed", "99", "--outline-only", "--intraspecific-mode", "--correct-allometry", "false"), defs)
                stopifnot(parsed$vouchers == "b.csv", parsed$max_k == 6L, parsed$seed == 99L, isTRUE(parsed$outline_only), isTRUE(parsed$intraspecific_mode), parsed$correct_allometry == "false")

                # 2. Exact mclust parameter counts (EII vs VVI vs VVV)
                np_eii <- .compute_exact_mclust_params(c(1L, 3L), p_dim = 4L, default_model = "EII")
                np_vvi <- .compute_exact_mclust_params(c(1L, 3L), p_dim = 4L, default_model = "VVI")
                np_vvv <- .compute_exact_mclust_params(c(1L, 3L), p_dim = 4L, default_model = "VVV")
                stopifnot(np_eii[2] < np_vvi[2], np_vvi[2] < np_vvv[2])

                # 3. Partial-NA reproductive & petiole backfilling
                repro_csv <- file.path({repr(str(tmp_p))}, "repro_partial.csv")
                man_csv <- file.path({repr(str(tmp_p))}, "extracted_leaf_manifest.csv")
                write.csv(data.frame(catalogNumber = c("C1", "C2"), plant_individual_id = c(0L, 0L), AR_inv = c(1.10, 1.45)), repro_csv, row.names = FALSE)
                write.csv(data.frame(catalogNumber = c("C1", "C2"), plant_individual_id = c(0L, 0L), leaf_zone = c("basal", "basal"), petiole_length_mm = c(15.0, 28.0), blade_length_mm = c(30.0, 28.0), petiole_to_blade_ratio = c(0.5, 1.0)), man_csv, row.names = FALSE)
                agg_part <- data.frame(catalogNumber = c("C1", "C2"), plant_individual_id = c(0L, 0L), AR_inv = c(1.10, NA_real_), petiole_length_mm = c(15.0, NA_real_), blade_length_mm = c(30.0, NA_real_), petiole_to_blade_ratio = c(0.5, NA_real_), stringsAsFactors = FALSE)
                agg_part <- merge_reproductive_traits_by_individual(agg_part, repro_csv)
                agg_part <- merge_manifest_petiole_traits_by_individual(agg_part, man_csv)
                stopifnot(abs(agg_part$AR_inv[2] - 1.45) < 1e-6, abs(agg_part$petiole_length_mm[2] - 28.0) < 1e-6)

                # 4. Whole + multi-island (PC1 = NA) joint classification & 60% cauline heteroblasty retention
                set.seed(77)
                n_sp <- 12L
                mk_sp <- function(sp, idx0, pc1_mu, ect1_mu, sol_mu, di_mu, caul_sol_mu, n_multi = 0L) {{
                  rows <- lapply(seq_len(n_sp), function(j) {{
                    is_mi <- j > (n_sp - n_multi)
                    has_caul <- j <= 8L # 67% cauline coverage (< 85%, >= 50%)
                    data.frame(
                      catalogNumber = sprintf("%s_%02d", idx0, j),
                      plant_individual_id = 0L,
                      basal_leaf_count = 3L,
                      cauline_leaf_count = if (has_caul) 2L else 0L,
                      morph_tier_source = "Basal_Rosette",
                      species_raw = sp,
                      scientificName = sp,
                      determiner_tier = if (j <= 8L) "Tier_1_Gold" else "Tier_3_Bronze",
                      PC1 = if (is_mi) NA_real_ else rnorm(1, pc1_mu, 0.20),
                      PC2 = if (is_mi) NA_real_ else rnorm(1, 0.5 * pc1_mu, 0.20),
                      PC3 = if (is_mi) NA_real_ else rnorm(1, 0.0, 0.15),
                      ECT_PC1 = rnorm(1, ect1_mu, 0.18),
                      ECT_PC2 = rnorm(1, -0.4 * ect1_mu, 0.18),
                      ECT_PC3 = rnorm(1, 0.0, 0.12),
                      ect_lobe_complexity = max(0.05, 1.2 + 0.6 * ect1_mu + rnorm(1, 0, 0.08)),
                      solidity = rnorm(1, sol_mu, 0.02),
                      q25_solidity = rnorm(1, sol_mu - 0.03, 0.02),
                      dissection_index = rnorm(1, di_mu, 0.08),
                      cauline_solidity_med = if (has_caul) rnorm(1, caul_sol_mu, 0.02) else NA_real_,
                      basal_vs_cauline_solidity_delta = if (has_caul) rnorm(1, sol_mu - caul_sol_mu, 0.02) else NA_real_,
                      aspect_ratio = rnorm(1, 1.8, 0.05),
                      centroid_size = 24.0 + 1.2 * j,
                      log_centsize = log10(24.0 + 1.2 * j),
                      is_calibrated_mm = TRUE,
                      scale_source = "ruler_calibrated",
                      A1 = if (is_mi) NA_real_ else rnorm(1, -1.1 + 0.05 * pc1_mu, 0.004),
                      B1 = if (is_mi) NA_real_ else 0.0,
                      C1 = if (is_mi) NA_real_ else 0.0,
                      D1 = if (is_mi) NA_real_ else rnorm(1, 0.45 - 0.03 * pc1_mu, 0.004),
                      A2 = if (is_mi) NA_real_ else rnorm(1, 0.03 * pc1_mu, 0.004),
                      B2 = if (is_mi) NA_real_ else 0.0,
                      C2 = if (is_mi) NA_real_ else 0.0,
                      D2 = if (is_mi) NA_real_ else rnorm(1, -0.02 * pc1_mu, 0.004),
                      stringsAsFactors = FALSE
                    )
                  }})
                  do.call(rbind, rows)
                }}
                df_h <- rbind(
                  mk_sp("Packera dubia", "DUB", -3.0, -2.2, 0.90, 1.15, 0.55, n_multi = 0L),
                  mk_sp("Packera millefolium", "MIL", 3.0, 2.8, 0.36, 3.20, 0.32, n_multi = 2L)
                )
                df_v <- data.frame(
                  catalogNumber = df_h$catalogNumber,
                  species_raw = df_h$species_raw,
                  determiner_tier = df_h$determiner_tier,
                  institutionCode = rep(c("NCU", "DUKE"), length.out = nrow(df_h)),
                  stateProvince = "NC",
                  decimalLatitude = 35.8,
                  decimalLongitude = -79.1,
                  doy = 125,
                  stringsAsFactors = FALSE
                )
                h_path <- file.path({repr(str(tmp_p))}, "harm_mi.csv")
                v_path <- file.path({repr(str(tmp_p))}, "vouch_mi.csv")
                f_path <- file.path({repr(str(tmp_p))}, "flags_mi.csv")
                p_path <- file.path({repr(str(tmp_p))}, "fig04b_gmm_cda_morphometrics.png")
                write.csv(df_h, h_path, row.names = FALSE)
                write.csv(df_v, v_path, row.names = FALSE)

                opts_mi <- parse_args_robust()
                opts_mi$vouchers <- v_path; opts_mi$harmonics <- h_path; opts_mi$output_flags <- f_path
                opts_mi$output_plot <- p_path; opts_mi$output_report <- file.path({repr(str(tmp_p))}, "bic_mi.csv")
                opts_mi$output_batch_audit <- file.path({repr(str(tmp_p))}, "batch_mi.csv")
                opts_mi$allometry_report <- file.path({repr(str(tmp_p))}, "allom_mi.csv"); opts_mi$allometry_pairwise <- file.path({repr(str(tmp_p))}, "pw_mi.csv")
                opts_mi$allometry_plot <- file.path({repr(str(tmp_p))}, "fig04a_foliar_allometry_audit.png")
                opts_mi$permutations <- 19L; opts_mi$seed <- 123L
                out_mi <- run_gmm_morphotools_pipeline(opts_mi)
                stopifnot(nrow(out_mi) == 24L)
                mi_rows <- out_mi[out_mi$catalogNumber %in% c("MIL_11", "MIL_12"), , drop = FALSE]
                stopifnot(nrow(mi_rows) == 2L, all(is.na(mi_rows$PC1)))
                stopifnot(all(mi_rows$cda_predicted_taxon == "Packera millefolium"))
                stopifnot(all(mi_rows$cda_posterior_prob > 0.75))
                stopifnot("cauline_solidity_med_std" %in% names(out_mi) || "basal_vs_cauline_solidity_delta_std" %in% names(out_mi))

                # 5. Single-group (g = 1) CDA sets loocv_accuracy = NA (not 1.0) & post-GMM plot omits fake 85.0%
                df_one <- df_h[1:12, , drop = FALSE]
                gmm_one <- list(classification = rep(1L, 12L), uncertainty = rep(0.0, 12L), best_k = 1L, model_name = "EII", selected_features = c("PC1", "PC2"))
                cda_one <- run_canonical_discriminant_analysis(df_one, c("PC1", "PC2"), "Packera dubia", gmm_res = gmm_one, intraspecific_mode = TRUE)
                stopifnot(is.na(cda_one$loocv_accuracy))
                gmm_two <- list(classification = rep(1:2, each = 6L), uncertainty = rep(0.05, 12L), best_k = 2L, model_name = "VVI", selected_features = c("PC1", "PC2"), bic_table = data.frame(K = 1:2, BIC = c(-40, -28), Two_Delta_BIC = c(0, 12)))
                cda_two <- run_canonical_discriminant_analysis(df_one, c("PC1", "PC2"), "Packera dubia", gmm_res = gmm_two, intraspecific_mode = TRUE)
                stopifnot(is.na(cda_two$loocv_accuracy))
                aud_two <- audit_misidentifications(df_one, cda_two, gmm_two, selected_features = c("PC1", "PC2"))
                intra_png <- file.path({repr(str(tmp_p))}, "intra_no_fake_loocv.png")
                render_cda_biplot_4panel_pdf(aud_two, cda_two, gmm_two, intra_png, c("PC1", "PC2"))
                stopifnot(file.exists(intra_png))

                # 6. Cross-run SQLite open-set recognition invariance under PC sign flips
                if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {{
                  db_p <- file.path({repr(str(tmp_p))}, "osr_cross_run.sqlite")
                  con <- DBI::dbConnect(RSQLite::SQLite(), db_p)
                  db_cong <- rbind(
                    data.frame(run_id = "old_run_A", catalogNumber = paste0("AN_", 1:8), species_standardized = "Packera anonyma", authority_weight = 1.0, morph_tier_source = "Basal_Rosette",
                               PC1 = rnorm(8, 9.0, 0.1), PC2 = rnorm(8, -9.0, 0.1),
                               A1 = rnorm(8, -0.92, 0.01), A2 = rnorm(8, 0.12, 0.005), D1 = rnorm(8, 0.35, 0.01), D2 = rnorm(8, -0.08, 0.005), solidity = rnorm(8, 0.88, 0.01), dissection_index = rnorm(8, 1.25, 0.03)),
                    data.frame(run_id = "old_run_A", catalogNumber = paste0("AU_", 1:8), species_standardized = "Packera aurea", authority_weight = 1.0, morph_tier_source = "Basal_Rosette",
                               PC1 = rnorm(8, -9.0, 0.1), PC2 = rnorm(8, 9.0, 0.1),
                               A1 = rnorm(8, -1.35, 0.01), A2 = rnorm(8, -0.18, 0.005), D1 = rnorm(8, 0.68, 0.01), D2 = rnorm(8, 0.14, 0.005), solidity = rnorm(8, 0.82, 0.01), dissection_index = rnorm(8, 1.45, 0.03))
                  )
                  DBI::dbWriteTable(con, "otu_morphometrics", db_cong, overwrite = TRUE)
                  DBI::dbDisconnect(con)
                  Sys.setenv(PACKERA_DB_PATH = db_p)
                  df_q <- data.frame(
                    catalogNumber = paste0("Q_", 1:6), species_standardized = "Packera dubia",
                    PC1 = rnorm(6, -9.0, 0.1), PC2 = rnorm(6, 9.0, 0.1), # deliberately flipped PC signs matching aurea if PCs were leaked!
                    A1 = rnorm(6, -0.92, 0.01), A2 = rnorm(6, 0.12, 0.005), D1 = rnorm(6, 0.35, 0.01), D2 = rnorm(6, -0.08, 0.005), # true harmonics match anonyma!
                    solidity = rnorm(6, 0.88, 0.01), dissection_index = rnorm(6, 1.25, 0.03),
                    stringsAsFactors = FALSE
                  )
                  osr_res <- compute_type_specimen_morph_control(df_q, X_mat = as.matrix(df_q[, c("PC1", "PC2")]), S_reg = diag(2), taxon_col = "species_standardized")
                  Sys.unsetenv("PACKERA_DB_PATH")
                  stopifnot(all(osr_res$min_congener_taxon == "Packera anonyma"))
                }}
                cat("STEP04_ZERO_TRUST_OK\n")
                """,
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_script)], stage_key='step04', cwd=PROJECT_ROOT)
            self.assertEqual(proc.returncode, 0, f"Step 04 zero-trust R test failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertIn("STEP04_ZERO_TRUST_OK", proc.stdout)
            self.assertTrue((tmp_p / "intra_no_fake_loocv.png").exists())
            self.assertEqual(list(tmp_p.glob("*.pdf")), [])


if __name__ == "__main__":
    unittest.main()
