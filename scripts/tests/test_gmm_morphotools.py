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

        if shutil.which("Rscript"):
            grp_out = cls._class_tmp_path / "manova_group_significance.csv"
            pw_out = cls._class_tmp_path / "manova_pairwise_posthoc.csv"
            r_batch_code = f"""
            source("scripts/morphometrics/pub_figures_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/morphometrics/04_gmm_morphotools.R")

            # Block 1: EPA Level III Ecoregion spatial helper agreement
            local({{
              res <- infer_epa_level3_ecoregion(
                lat = c(35.9132, 34.2257, 32.5404, 31.6035, NA, NA, NA),
                lon = c(-79.0558, -77.9447, -82.9037, -94.6555, NA, NA, NA),
                state_province = c("NC", "NC", "GA", "TX", "North Carolina", "North Carolina", "NC"),
                county = c("", "", "", "", "Orange County", "Wake County", "")
              )
              cat("EPA_L3:", paste(ifelse(is.na(res), "NA", res), collapse = "|"), "\\n", sep = "")
            }})

            # Block 2: Direct run_multivariate_manova_audit single-species & collinearity safety
            local({{
              set.seed(17)
              df <- data.frame(
                catalogNumber = sprintf("V%02d", 1:18),
                species_standardized = rep("Packera dubia", 18),
                determiner_tier = rep("Tier_1_Gold", 18),
                morph_tier_source = rep("Tier_1_Pure", 18),
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
                    morph_tier_source = "Tier_1_Pure",
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
            cls.r_unit_proc = subprocess.run(
                ["Rscript", str(r_batch_file)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
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

    def test_step04_static_contracts_and_line_limit(self):
        """Verify 04_gmm_morphotools.R (<= 650 lines) retains all allometry, subspace, GMM, CDA, and batch-audit contracts."""
        self.assertTrue(self.r_script_path.exists(), "04_gmm_morphotools.R must exist.")
        content = self.r_script_path.read_text(encoding="utf-8")
        lines = content.splitlines()
        self.assertLessEqual(len(lines), 650, f"Script exceeds line limit: {len(lines)}")

        for required_token in (
            "--reproductive",
            "--outline-only",
            "AR_inv",
            "AR_inv_std",
            "--output-batch-audit",
            "audit_institutional_batch_effects",
            "institutionCode",
            "0.05",
            "1 row = 1 voucher",
            "group_by(catalogNumber)",
            "--correct-allometry",
            "test_and_correct_allometry",
            "generate_allometry_diagnostic_plot",
            "[INFO] Foliar allometry detected",
            "[INFO] Foliar allometry negligible",
            "evaluate_ontogenetic_allometry",
            "is_calibrated_mm",
            "Tier_1_Gold",
            "lm(Y ~ log_centsize + scientificName",
            "Y_adj <- Y - outer(log_centsize - mean(log_centsize, na.rm = TRUE), beta_within)",
            'grepl("^tier_?[12]", efa_df$assigned_tier, ignore.case = TRUE)',
            "RRPP::lm.rrpp",
            "RRPP::pairwise",
            "allometry_pairwise_trajectories.csv",
            "clustvarsel::clustvarsel",
            "is_tier1_pure",
            "morph_tier_source",
            "n_tier1_leaves",
            "n_tier2_leaves",
            '"Tier_1_Pure"',
            "PC1 ~ log_centsize * scientificName",
            "select_morphometric_subspace",
            "floor(N_min / 3",
            "rep(1/g, g)",
            "Kass & Raftery",
            "Unassigned_Packera",
        ):
            self.assertIn(required_token, content, f"Missing expected Step 04 contract token: {required_token}")

        self.assertNotIn("Glabrescent_P_dubia_Misidentified_as_P_anonyma", content)
        self.assertNotIn(
            'target_taxa <- c("Packera anonyma", "Packera dubia", "Packera plattensis", "Packera paupercula")',
            content,
        )

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

    def test_step02_and_preflight_dynamic_non_default_packera_support(self):
        """Verify Step 02 standardization/figures, config.yaml, and run_pipeline.sh support dynamic genus-wide Packera species."""
        seg_mod = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        std_fn = seg_mod._standardize_packera_species_py
        expected_map = {
            "Packera schweinitziana": "Packera schweinitziana",
            "Packera serpenticola": "Packera serpenticola",
            "Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve": "Packera millefolium",
            "Packera aurea (L.) Á.Löve & D.Löve": "Packera aurea",
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve": "Packera streptanthifolia",
            "Packera cymbalaria (Pursh) W.A.Weber & Á.Löve": "Packera cymbalaria",
            "Packera sanguisorbae (DC.) C.Jeffrey": "Packera sanguisorbae",
        }
        for raw_name, canon_name in expected_map.items():
            self.assertEqual(std_fn(raw_name), canon_name, f"Step 02 failed to standardize {raw_name!r}")

        with patch(
            "scripts.core.harvester.species.name_backbone",
            return_value={
                "matchType": "EXACT",
                "status": "SYNONYM",
                "species": "Packera schweinitziana",
                "canonicalName": "Senecio robbinsii",
            },
        ), patch("scripts.core.harvester.species.name_usage", return_value={"results": []}):
            self.assertEqual(std_fn("Senecio robbinsii"), "Packera schweinitziana")

        taxa_8 = [
            "Packera aurea",
            "Packera cana",
            "Packera cymbalaria",
            "Packera millefolium",
            "Packera sanguisorbae",
            "Packera schweinitziana",
            "Packera serpenticola",
            "Packera streptanthifolia",
        ]
        fig_rows = []
        for idx, sp in enumerate(taxa_8):
            for k in range(4):
                fig_rows.append({
                    "catalogNumber": f"FIG_{idx}_{k}",
                    "plant_individual_id": 0,
                    "leaf_id": 1,
                    "assigned_tier": "Tier_1_Direct",
                    "scientificName": sp,
                    "solidity": 0.35 + 0.05 * k if sp == "Packera millefolium" else 0.80 + 0.03 * k,
                    "ucs_score": 0.92,
                    "aspect_ratio": 1.4 + 0.3 * k,
                    "pixels_per_mm": 10.0,
                    "is_calibrated_mm": True,
                })
        with tempfile.TemporaryDirectory() as tmpdir:
            saved = seg_mod.generate_step02_publication_figures(
                extracted_df=pd.DataFrame(fig_rows),
                failed_df=pd.DataFrame(),
                repro_df=pd.DataFrame(),
                pheno_df=pd.DataFrame(),
                figures_dir=Path(tmpdir) / "figures",
            )
            self.assertTrue(Path(saved["fig02a_pdf"]).exists())
            self.assertTrue(Path(saved["fig02b_pdf"]).exists())

        cfg = PipelineConfig.from_yaml()
        self.assertEqual(cfg.taxa.target_species, [])
        self.assertEqual(cfg.taxa.outgroups, [])

        sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        self.assertIn("base_outgroups = present_taxa", sh_src)
        self.assertNotIn("mod05.TARGET_TAXA", sh_src)

    def test_dual_intraspecific_allometry_morphotype_and_supervised_epa_ecoregion_cda(self):
        """Verify single-taxon (--intraspecific-mode) allometry, partial AR_inv & Tier_2_Fallback retention, Track 1 Morphotype_K* CDA, and Track 2 EPA Level III Ecoregion CDA."""
        if not shutil.which("Rscript"):
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
                    "assigned_tier": "Tier_2_Reflected" if is_t2_fb else "Tier_1_Direct",
                    "n_tier1_leaves": 0 if is_t2_fb else 3,
                    "n_tier2_leaves": 2 if is_t2_fb else 0,
                    "morph_tier_source": "Tier_2_Fallback" if is_t2_fb else "Tier_1_Pure",
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
                    "ucs_score": 0.94,
                    "capitulum_aspect_ratio": np.nan if j % 5 == 0 else float(1.15 + 0.08 * c1 + rng_ar.normal(0.0, 0.015)),
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
            allom_pdf = fig_dir / "allometry.pdf"
            eco_sum_csv = rep_dir / "ecoregion_supervised_cda_summary.csv"
            cda_pdf = fig_dir / "cda_passive_projection.pdf"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = subprocess.run(
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
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"Dual/Intraspecific Step 04 failed:\n{proc.stderr}\n{proc.stdout}")

            # 0. Single-taxon allometry audit CSV & PDF + zero dropped Tier_1 / Tier_2_Fallback / partial-AR_inv vouchers
            self.assertTrue(allom_csv.exists(), "Allometry audit CSV must be exported in single-taxon mode")
            self.assertTrue(allom_pdf.exists() and allom_pdf.read_bytes()[:5] == b"%PDF-")
            flags_df = pd.read_csv(flags_csv)
            self.assertEqual(len(flags_df), 40, "All 40 vouchers (including Tier_2_Fallback and partial AR_inv) must be retained")
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
                "fig04a_gmm_bic_and_uncertainty",
                "fig04b_cda_biplot_and_allometry",
                "fig04_gmm_cda_clustering",
                "fig04c_supervised_ecoregion_cda",
            ):
                pdf_f = fig_dir / f"{stem}.pdf"
                png_f = fig_dir / f"{stem}.png"
                self.assertTrue(pdf_f.exists() and pdf_f.stat().st_size > 3000, f"Missing/empty {pdf_f}")
                self.assertTrue(png_f.exists() and png_f.stat().st_size > 10000, f"Missing/empty {png_f}")

    def test_tier3_bronze_invariance_and_multivariate_manova_audit(self):
        """Verify Tier_3_Bronze zero-influence invariance across allometry, GMM, CDA, and MANOVA, plus batch audit R^2 < 0.05 and CDA subtitles."""
        if not shutil.which("Rscript"):
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
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
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
                    "ucs_score": 0.94,
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
                cda_pdf = fig_dir / "cda_passive_projection.pdf"
                df_h.to_csv(harm_csv, index=False)
                df_v.to_csv(vouch_csv, index=False)
                proc = subprocess.run(
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
                        "--allometry-plot", str(fig_dir / "allometry.pdf"),
                        "--correct-allometry", "true",
                        "--permutations", "19",
                    ],
                    cwd=str(PROJECT_ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
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

            fig04_txt = self._extract_pdf_text_streams(fig_base / "fig04_gmm_cda_clustering.pdf")
            fig04c_txt = self._extract_pdf_text_streams(fig_base / "fig04c_supervised_ecoregion_cda.pdf")
            sp_row = grp_df[grp_df["grouping_factor"] == "species_standardized"].iloc[0]
            self.assertIn("Omnibus MANOVA: Pillai's Trace V =", fig04_txt)
            self.assertIn(f"= {float(sp_row['approx_F']):.2f}, p =", fig04_txt)
            self.assertIn("Omnibus MANOVA: Pillai's Trace V =", fig04c_txt)
            self.assertIn(f"= {float(rg_row['approx_F']):.2f}, p =", fig04c_txt)

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
        if not shutil.which("Rscript"):
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
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
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
                    "ucs_score": 0.94,
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
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 1, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "petiole_length_mm": 32.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": round(32.0 / 28.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 1, "leaf_id": 2, "assigned_tier": "Tier_1_Direct", "petiole_length_mm": 32.0, "blade_length_mm": 28.0, "petiole_to_blade_ratio": round(32.0 / 28.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 2, "leaf_id": 1, "assigned_tier": "Tier_1_Direct", "petiole_length_mm": 12.0, "blade_length_mm": 30.0, "petiole_to_blade_ratio": round(12.0 / 30.0, 4)},
            {"catalogNumber": "MULTI_SHEET_01", "plant_individual_id": 2, "leaf_id": 2, "assigned_tier": "Tier_1_Direct", "petiole_length_mm": 12.0, "blade_length_mm": 30.0, "petiole_to_blade_ratio": round(12.0 / 30.0, 4)},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            harm_csv = tmp_p / "leaf_efa_harmonics.csv"
            vouch_csv = tmp_p / "curated_vouchers.csv"
            man_csv = tmp_p / "extracted_leaf_manifest.csv"
            flags_csv = tmp_p / "morphometric_flags.csv"
            bic_csv = tmp_p / "gmm_bayes_factors_summary.csv"
            allom_csv = tmp_p / "allometry_significance_audit.csv"
            allom_pdf = tmp_p / "shape_allometry_regression.pdf"
            cda_pdf = tmp_p / "cda_passive_projection.pdf"

            df_in = pd.DataFrame(rows)
            df_in.to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)
            pd.DataFrame(man_rows).to_csv(man_csv, index=False)

            proc = subprocess.run(
                [
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
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
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
        if not shutil.which("Rscript"):
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
                    "assigned_tier": "Tier_1_Direct",
                    "n_tier1_leaves": 3,
                    "n_tier2_leaves": 0,
                    "morph_tier_source": "Tier_1_Pure",
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
                    "ucs_score": 0.94,
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
            cda_pdf = tmp_p / "cda_passive_projection.pdf"

            pd.DataFrame(rows).to_csv(harm_csv, index=False)
            pd.DataFrame(vouch_rows).to_csv(vouch_csv, index=False)

            proc = subprocess.run(
                [
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
                    "--allometry-plot", str(tmp_p / "allometry.pdf"),
                    "--permutations", "19",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
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


if __name__ == "__main__":
    unittest.main()
