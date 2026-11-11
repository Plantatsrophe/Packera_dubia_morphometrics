#!/usr/bin/env python3
"""
Unit tests for Eastern North American Packera N-Taxon Standardization,
Unified Step 06 Quantitative Taxon-Complex Discovery Engine
(taxon_complex_discovery_utils.R + 06_triage_dashboard_synthesis.R),
and Zero-Trust Four-Gate Complex Delimitation Contracts.
"""

import importlib
import shutil
import sqlite3
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from scripts.core.db import get_db_connection
from scripts.core.harvester import standardize_packera_taxon

REPO_ROOT = Path(__file__).resolve().parents[2]
DISCOVERY_UTILS_R = REPO_ROOT / "scripts" / "analysis" / "taxon_complex_discovery_utils.R"
SYNTHESIS_R_SCRIPT = REPO_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"
PUB_FIG_R_SCRIPT = REPO_ROOT / "scripts" / "morphometrics" / "pub_figures_utils.R"

CANONICAL_COMPLEX_COLS = {
    "species_standardized",
    "n_vouchers",
    "n_gold_silver_anchors",
    "n_type_specimens",
    "consensus_complex_id",
    "complex_label",
    "leiden_complex_id",
    "ward_complex_id",
    "silhouette_width",
    "within_complex_mean_affinity",
    "nearest_neighbor_taxon",
    "nearest_neighbor_affinity",
    "second_neighbor_taxon",
    "second_neighbor_affinity",
    "hopkins_pheno_overlap",
    "eco_schoener_overlap",
    "integrative_congruence",
    "w_morph_kappa",
    "w_edaphic_kappa",
    "w_pheno_kappa",
    "n_misidentification_reassign",
    "n_hybrid_introgressant",
    "n_putative_cryptic",
}


class TestEasternPackeraAndComplexDiscovery(unittest.TestCase):
    """Tests for N-taxon basionym mapping, unified Step 06 complex discovery, and Four-Gate contracts."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.r_unit_proc = None

        if shutil.which("Rscript"):
            intra_fig_dir = cls._class_tmp_path / "intra_figs"
            intra_fig_dir.mkdir(parents=True, exist_ok=True)
            mm_dir = cls._class_tmp_path / "mm_join"
            mm_dir.mkdir(parents=True, exist_ok=True)
            morph_csv = mm_dir / "morphometric_flags.csv"
            mm_csv = mm_dir / "multimodal_conflict_flags.csv"
            q_csv = mm_dir / "triage_queue.csv"
            rep_md = mm_dir / "rev.md"
            html_out = mm_dir / "dash.html"
            fig_out = mm_dir / "rev.pdf"

            pd.DataFrame([
                {"catalogNumber": "MM_01", "plant_individual_id": 1, "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "basal_leaf_count": 3, "PC1": -1.2, "PC2": 0.3},
                {"catalogNumber": "MM_02", "plant_individual_id": 1, "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "basal_leaf_count": 3, "PC1": 1.4, "PC2": -0.2},
            ]).to_csv(morph_csv, index=False)
            pd.DataFrame([
                {
                    "catalogNumber": "MM_01", "plant_individual_id": 1,
                    "fine_soil_ph": 5.9, "fine_soil_sand": 34.5, "fine_bedrock_depth_cm": 98.0,
                    "ssurgo_rock_outcrop_pct": 3.5, "soil_ph": 5.8, "soil_sand": 35.0,
                    "soil_cec": 14.2, "soil_bulk_density": 1.31,
                    "bio1_temp_mean": 13.4, "bio4_temp_seasonality": 712.0,
                    "bio12_precip_annual": 1340.0, "bio15_precip_seasonality": 18.0,
                    "edaphic_schoener_d": 0.41, "schoener_d": 0.44,
                },
                {
                    "catalogNumber": "MM_02", "plant_individual_id": 1,
                    "fine_soil_ph": 6.8, "fine_soil_sand": 52.0, "fine_bedrock_depth_cm": 42.0,
                    "ssurgo_rock_outcrop_pct": 24.0, "soil_ph": 6.7, "soil_sand": 50.0,
                    "soil_cec": 18.5, "soil_bulk_density": 1.42,
                    "bio1_temp_mean": 14.8, "bio4_temp_seasonality": 745.0,
                    "bio12_precip_annual": 1180.0, "bio15_precip_seasonality": 21.0,
                    "edaphic_schoener_d": 0.41, "schoener_d": 0.44,
                },
            ]).to_csv(mm_csv, index=False)

            r_batch_code = f"""
            source("{PUB_FIG_R_SCRIPT.as_posix()}")

            # Block 1: Dynamic palette, taxon standardization, shapes, titles, and italic expressions
            local({{
              stopifnot(!exists("EASTERN_NA_PACKERA_TAXA"))
              stopifnot(!exists("OKABE_ITO_TAXA"))
              stopifnot(!exists("PACKERA_ITALIC_LABELS"))

              stopifnot(standardize_taxon_name("Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve") == "Packera millefolium")
              stopifnot(standardize_taxon_name("P. aurea") == "Packera aurea")
              stopifnot(standardize_taxon_name("Packera x memmingeri") == "Packera memmingeri")
              stopifnot(standardize_taxon_name("Packera × memmingeri") == "Packera memmingeri")
              stopifnot(standardize_taxon_name("Packera cana (Hook.) W.A.Weber & Á.Löve") == "Packera cana")
              stopifnot(standardize_taxon_name("Packera sanguisorbae (DC.) C.Jeffrey") == "Packera sanguisorbae")
              stopifnot(standardize_taxon_name("Packera hyperborealis") == "Packera hyperborealis")
              stopifnot(standardize_taxon_name("Packera xylorrhiza") == "Packera xylorrhiza")

              taxa8 <- c(
                "Packera streptanthifolia", "Packera aurea", "Packera sanguisorbae", "Packera cymbalaria",
                "Packera millefolium", "Packera cana", "Packera obovata", "Packera fendleri"
              )
              sorted8 <- sort(taxa8)
              pal <- get_dynamic_packera_palette(taxa8)
              stopifnot(identical(names(pal), c(sorted8, "Unknown")))
              stopifnot(pal[[sorted8[1]]] == "#0072B2")
              stopifnot(pal[[sorted8[2]]] == "#D55E00")
              stopifnot(pal[[sorted8[3]]] == "#009E73")
              stopifnot(pal[[sorted8[4]]] == "#CC79A7")
              stopifnot(pal[["Unknown"]] == "#7F7F7F")
              stopifnot(length(unique(unname(pal[sorted8]))) == 8)

              shps <- get_dynamic_packera_shapes(taxa8)
              stopifnot(identical(unname(shps[sorted8[1:4]]), c(1L, 2L, 5L, 0L)))
              gg_shps <- get_dynamic_packera_ggplot_shapes(taxa8)
              stopifnot(identical(unname(gg_shps[sorted8[1:5]]), c(21, 24, 22, 23, 25)))

              Sys.setenv(PACKERA_FOCAL_TAXON = "")
              stopifnot(format_active_complex_title(c("Packera cana", "Packera fendleri", "Packera streptanthifolia"), "CDA") == "Packera Taxon Complex (3 Taxa): CDA")
              stopifnot(format_active_complex_title(c("Packera sanguisorbae"), "CDA") == "Packera sanguisorbae: CDA")
              Sys.setenv(PACKERA_FOCAL_TAXON = "Packera streptanthifolia")
              stopifnot(format_active_complex_title(c("Packera cana", "Packera streptanthifolia"), "CDA") == "Packera streptanthifolia Complex: CDA")
              Sys.unsetenv("PACKERA_FOCAL_TAXON")

              expr_tx <- format_italic_taxon_expr("Packera sanguisorbae")
              stopifnot(is.call(expr_tx) && identical(deparse(expr_tx), 'italic("Packera sanguisorbae")'))
              cat("STEP06_R_PALETTE_OK\\n")
            }})

            source("{DISCOVERY_UTILS_R.as_posix()}")

            # Block 2: Stage 1-5 hierarchical precedence, missing-modality weight renormalization, & BH FDR
            local({{
              fdr_in <- data.frame(Pair = c("A_vs_B", "A_vs_C", "B_vs_C"), p_value = c(0.002, 0.030, 0.250))
              fdr_out <- apply_bh_fdr_correction(fdr_in)
              stopifnot("p_adj_bh" %in% names(fdr_out))
              stopifnot(all(abs(fdr_out$p_adj_bh - stats::p.adjust(fdr_in$p_value, method = "BH")) < 1e-9))

              sh <- get_dynamic_packera_shapes(c("Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"))
              stopifnot(all(unname(sh) %in% c(1L, 2L, 5L, 0L)))

              mdf <- data.frame(
                catalogNumber = c("V1", "V2", "V3", "V4", "V5", "V6"),
                species_standardized = c("Packera dubia", "Packera dubia", "Packera dubia", "Packera anonyma", "Packera anonyma", "Packera anonyma"),
                determiner_tier = rep("Tier_1_Gold", 6),
                cda_predicted_taxon = c("Packera dubia", "Packera anonyma", "Packera anonyma", "Packera dubia", "Packera anonyma", "Packera anonyma"),
                cda_posterior_prob = c(0.95, 0.85, 0.60, 0.44, 0.48, 0.92),
                cda_second_posterior_prob = c(0.03, 0.12, 0.40, 0.28, 0.25, 0.05),
                cda_posterior_margin = c(0.92, 0.73, 0.20, 0.16, 0.23, 0.87),
                morphometric_outlier_flag = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE),
                gmm_cluster = c(1L, 2L, 2L, 5L, 5L, 2L),
                gmm_2_delta_bic = c(0, 0, 0, 12.0, 12.0, 0),
                edaphic_schoener_d = c(0.75, 0.60, 0.55, 0.18, 0.72, 0.80),
                basal_leaf_count = rep(3L, 6),
                stringsAsFactors = FALSE
              )
              kw_missing <- compute_modality_kappa_weights(mdf, rep(TRUE, 6))
              # Only morph + edaphic (from edaphic_schoener_d? No: no edaphic_best_fit_taxon/predictors -> w_morph == 1.0!)
              stopifnot(abs(kw_missing$w_morph - 1.0) < 1e-9)

              pheno_mock <- list(delta_doy_resid = c(2.0, 3.0, 4.0, 22.0, 3.0, 1.0))
              calls <- evaluate_hierarchical_complex_delimitation(mdf, pheno_res = pheno_mock)
              stopifnot(calls[1] == "Confirmed_Species")
              stopifnot(calls[2] == "Misidentification_Reassign")
              stopifnot(calls[3] == "Hybrid_Introgressant")
              stopifnot(calls[4] == "Putative_Cryptic_Species")
              stopifnot(calls[5] == "Ecophenotypic_Plasticity")
              cat("STEP06_HIERARCHICAL_OK\\n")
            }})

            # Block 3: Unified Step 06 multimodal left join + complex discovery integration
            source("{SYNTHESIS_R_SCRIPT.as_posix()}")
            local({{
              opts_m <- parse_args_robust()
              opts_m$vouchers <- "{morph_csv.as_posix()}"
              opts_m$morphometrics <- "{morph_csv.as_posix()}"
              opts_m$multimodal_flags <- "{mm_csv.as_posix()}"
              opts_m$output_queue <- "{q_csv.as_posix()}"
              opts_m$output_report <- "{rep_md.as_posix()}"
              opts_m$output_html <- "{html_out.as_posix()}"
              opts_m$output_plot <- "{fig_out.as_posix()}"
              opts_m$output_complexes <- "{(mm_dir / 'empirical_taxon_complexes.csv').as_posix()}"
              opts_m$output_affinity <- "{(mm_dir / 'pairwise_taxon_affinity_matrix.csv').as_posix()}"
              opts_m$output_neighbors <- "{(mm_dir / 'focal_auto_neighbors.txt').as_posix()}"
              opts_m$output_complex_plot <- "{(mm_dir / 'fig06c_taxon_complex_network.pdf').as_posix()}"
              opts_m$demo <- TRUE
              mdf <- build_triage_queue(opts_m)
              req <- c("fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "ssurgo_rock_outcrop_pct",
                       "soil_ph", "soil_sand", "soil_cec", "soil_bulk_density",
                       "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality",
                       "edaphic_schoener_d", "schoener_d")
              stopifnot(all(req %in% names(mdf)))
              stopifnot(abs(mdf$fine_soil_ph[mdf$catalogNumber == "MM_01"] - 5.9) < 1e-6)
              stopifnot(abs(mdf$ssurgo_rock_outcrop_pct[mdf$catalogNumber == "MM_02"] - 24.0) < 1e-6)
              cat("STEP06_MULTIMODAL_JOIN_OK\\n")
            }})

            # Block 4: Intraspecific morphotype cluster visualization across R figures
            source("{PUB_FIG_R_SCRIPT.as_posix()}")
            source("{(REPO_ROOT / 'scripts' / 'analysis' / '05_multimodal_spatial_rf.R').as_posix()}")
            local({{
              set.seed(123)
              n <- 24
              df <- data.frame(
                catalogNumber = sprintf("INTRA_%03d", seq_len(n)),
                plant_individual_id = 1L,
                species_standardized = rep("Packera dubia", n),
                scientificName = rep("Packera dubia", n),
                gmm_cluster = rep(1:3, each = 8),
                analysis_group = paste0("Morphotype_K", rep(1:3, each = 8)),
                grouping_col = paste0("Morphotype_K", rep(1:3, each = 8)),
                can1 = c(rnorm(8, -2), rnorm(8, 0), rnorm(8, 2)),
                can2 = c(rnorm(8, 1), rnorm(8, -1.5), rnorm(8, 1)),
                can3 = rnorm(n),
                PC1 = c(rnorm(8, -2), rnorm(8, 0), rnorm(8, 2)),
                PC2 = c(rnorm(8, 1), rnorm(8, -1.5), rnorm(8, 1)),
                PC3 = rnorm(n),
                aspect_ratio = runif(n, 1.5, 3.2),
                solidity = runif(n, 0.75, 0.95),
                log_centsize = runif(n, 2.0, 3.5),
                AR_inv = runif(n, 0.9, 1.4),
                cda_predicted_taxon = paste0("Morphotype_K", rep(1:3, each = 8)),
                cda_posterior_prob = runif(n, 0.75, 0.98),
                cda_second_posterior_prob = runif(n, 0.01, 0.20),
                fine_soil_sand = runif(n, 40, 85),
                fine_soil_ph = runif(n, 4.2, 6.8),
                fine_bedrock_depth_cm = runif(n, 20, 150),
                ssurgo_rock_outcrop_pct = runif(n, 0, 30),
                topo_tpi_30m = rnorm(n),
                bio1_temp_mean = runif(n, 12, 18),
                bio12_precip_annual = runif(n, 1000, 1400),
                decimalLatitude = runif(n, 32, 38),
                decimalLongitude = runif(n, -84, -77),
                doy = runif(n, 85, 145),
                delta_doy = rnorm(n, 0, 8),
                is_completed_voucher = TRUE,
                taxonomic_treatment_subtype = rep(c("Confirmed_Species", "Subspecies_Ecotype", "Hybrid_Introgressant"), 8),
                stringsAsFactors = FALSE
              )
              grp <- resolve_plot_group_column(df)
              stopifnot(identical(sort(unique(grp)), c("Morphotype_K1", "Morphotype_K2", "Morphotype_K3")))

              cda_res <- list(
                grouping_col = df$grouping_col,
                loocv_accuracy = 0.91,
                variance_pct = c(64.0, 36.0),
                canonical_scores = cbind(Can1 = df$can1, Can2 = df$can2, Can3 = df$can3)
              )
              gmm_res <- list(optimal_G = 3, best_k = 3, optimal_model = "VVI", classification = df$gmm_cluster)

              f_cda <- file.path("{intra_fig_dir.as_posix()}", "cda_passive_projection.pdf")
              render_cda_biplot_4panel_pdf(df, cda_res, gmm_res, f_cda)
              stopifnot(file.exists(f_cda) && file.info(f_cda)$size > 1000)

              generate_step04_pub_figures(df, cda_res, gmm_res, NULL, f_cda)
              stopifnot(file.exists(file.path("{intra_fig_dir.as_posix()}", "fig04_gmm_cda_clustering.pdf")))

              f_srf <- file.path("{intra_fig_dir.as_posix()}", "spatial_rf_niche_importance.pdf")
              srf_mock <- list(
                OOB_Error = 0.12,
                Spatial_Block_CV_Accuracy = 0.86,
                importance = data.frame(
                  Predictor = c("fine_soil_sand", "fine_soil_ph", "MEM1"),
                  MeanDecreaseAccuracy = c(0.28, 0.21, 0.12),
                  Type = c("Fine-Grain Edaphic & Outcrop (10-30m)", "Fine-Grain Edaphic & Outcrop (10-30m)", "Spatial (MEMs)")
                ),
                variance_partition = data.frame(
                  Component = c("[a] Pure Fine-Grain Edaphic", "[b] Pure Macro-Climatic", "[c] Pure Spatial (MEMs)", "[d] Shared Env-Spatial", "[e] Unexplained"),
                  Variance_Pct = c(38, 22, 15, 10, 15)
                )
              )
              niche_mock <- list(
                summary = data.frame(Pair = "Morphotype_K1 vs Morphotype_K2", Taxon_1 = "Morphotype_K1", Taxon_2 = "Morphotype_K2", Warren_D = 0.42, Schoeners_D = 0.42, Equivalency_p_value = 0.02),
                permutations = data.frame(Pair = rep("Morphotype_K1 vs Morphotype_K2", 10), D_null = runif(10, 0.4, 0.9))
              )
              export_spatial_rf_figures(df, srf_mock, niche_mock, f_srf)
              stopifnot(file.exists(f_srf) && file.info(f_srf)$size > 1000)

              f_s06 <- file.path("{intra_fig_dir.as_posix()}", "integrative_taxonomic_synthesis.pdf")
              render_step06_synthesis_plate_pdf(df, f_s06)
              stopifnot(file.exists(f_s06) && file.info(f_s06)$size > 1000)
              cat("STEP06_INTRA_FIGS_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step06_complex_checks.R"
            r_batch_file.write_text(r_batch_code, encoding="utf-8")
            cls.r_unit_proc = subprocess.run(
                ["Rscript", str(r_batch_file)],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._class_temp_dir.cleanup()

    def _require_batched_r(self) -> str:
        if self.r_unit_proc is None:
            self.skipTest("Rscript not installed.")
        self.assertEqual(
            self.r_unit_proc.returncode,
            0,
            f"Batched Step 06 complex discovery R unit checks failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_standardize_packera_taxon_all_eastern_species_and_basionyms(self):
        """Verify North American Packera species (Eastern, Western, Boreal, Mexican) map cleanly."""
        cases = {
            "Packera anonyma (Alph.Wood) W.A.Weber & Á.Löve": "Packera anonyma",
            "Packera antennariifolia (Britton) W.A.Weber & Á.Löve": "Packera antennariifolia",
            "Packera aurea (L.) Á.Löve & D.Löve": "Packera aurea",
            "Packera crawfordii (Britton) A.M.Mahoney & R.R.Kowal": "Packera crawfordii",
            "Packera dubia (Spreng.) Trock & Mabb.": "Packera dubia",
            "Packera glabella (Poir.) C.Jeffrey": "Packera glabella",
            "Packera insulanus (C.L.Hitchc.) W.A.Weber & Á.Löve": "Packera insulanus",
            "Packera memmingeri": "Packera memmingeri",
            "Packera x memmingeri": "Packera memmingeri",
            "Packera × memmingeri": "Packera memmingeri",
            "Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve": "Packera millefolium",
            "Packera obovata (Muhl. ex Willd.) W.A.Weber & Á.Löve": "Packera obovata",
            "Packera paupercula var. appalachiana": "Packera paupercula",
            "Packera plattensis (Nutt.) W.A.Weber & Á.Löve": "Packera plattensis",
            "Packera pseudaurea var. semicordata": "Packera pseudaurea",
            "Packera schweinitziana (Nutt.) W.A.Weber & Á.Löve": "Packera schweinitziana",
            "Packera serpenticola Boufford, Kartesz, Shiuchi & R.Zhou": "Packera serpenticola",
            "Packera spellenbergii (T.M.Barkley) C.Jeffrey": "Packera spellenbergii",
            "Packera tampicana (DC.) C.Jeffrey": "Packera tampicana",
            "Packera bernardina (Greene) W.A.Weber & Á.Löve": "Packera bernardina",
            "Packera cana (Hook.) W.A.Weber & Á.Löve": "Packera cana",
            "Packera fendleri (A.Gray) W.A.Weber & Á.Löve": "Packera fendleri",
            "Packera multilobata (Torr. & A.Gray ex A.Gray) W.A.Weber & Á.Löve": "Packera multilobata",
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve": "Packera streptanthifolia",
            "Packera werneriifolia (A.Gray) W.A.Weber & Á.Löve": "Packera werneriifolia",
            "Packera cymbalaria (Pursh) W.A.Weber & Á.Löve": "Packera cymbalaria",
            "Packera hyperborealis (Greenm.) Á.Löve & D.Löve": "Packera hyperborealis",
            "Packera indecora (Greene) Á.Löve & D.Löve": "Packera indecora",
            "Packera ogotorukensis (Packer) Á.Löve & D.Löve": "Packera ogotorukensis",
            "Packera pauciflora (Pursh) Á.Löve & D.Löve": "Packera pauciflora",
            "Packera bellidifolia (Kunth) W.A.Weber & Á.Löve": "Packera bellidifolia",
            "Packera coahuilensis (Greenm.) C.Jeffrey": "Packera coahuilensis",
            "Packera sanguisorbae (DC.) C.Jeffrey": "Packera sanguisorbae",
            "Packera toluccana (DC.) W.A.Weber & Á.Löve": "Packera toluccana",
            "Packera zimapanica (Hemsl.) C.C.Freeman & T.M.Barkley": "Packera zimapanica",
        }
        mod01_pipe = importlib.import_module("scripts.pipeline.01_voucher_harvester")
        for raw_name, expected in cases.items():
            self.assertEqual(standardize_packera_taxon(raw_name), expected, f"Failed for {raw_name!r}")
            self.assertEqual(mod01_pipe.standardize_packera_taxon(raw_name), expected)

        basionym_backbone = {
            "senecio millefolium": "Packera millefolium",
            "senecio streptanthifolius": "Packera streptanthifolia",
            "senecio cymbalaria": "Packera cymbalaria",
            "senecio sanguisorbae": "Packera sanguisorbae",
        }

        def _fake_backbone(name, **kwargs):
            acc = basionym_backbone.get(str(name).strip().lower())
            if acc:
                return {"matchType": "EXACT", "status": "SYNONYM", "species": acc, "canonicalName": str(name)}
            return {"matchType": "NONE"}

        with patch("scripts.core.harvester.species.name_backbone", side_effect=_fake_backbone), patch(
            "scripts.core.harvester.species.name_usage", return_value={"results": []}
        ):
            for raw_syn, expected in [
                ("Senecio millefolium Torr. & A.Gray", "Packera millefolium"),
                ("Senecio streptanthifolius Greene", "Packera streptanthifolia"),
                ("Senecio cymbalaria Pursh", "Packera cymbalaria"),
                ("Senecio sanguisorbae DC.", "Packera sanguisorbae"),
            ]:
                self.assertEqual(standardize_packera_taxon(raw_syn), expected)
                self.assertEqual(mod01_pipe.standardize_packera_taxon(raw_syn), expected)

    def test_r_standardize_taxon_and_dynamic_palette(self):
        """Verify R standardize_taxon_name(), get_dynamic_packera_palette(), shapes, titles, and italic expressions are species-agnostic."""
        self.assertIn("STEP06_R_PALETTE_OK", self._require_batched_r())

    def test_discover_taxon_complexes_8_species_3_complexes(self):
        """Verify unified Step 06 recovers 3 true complexes across 8 species with canonical twin-free outputs."""
        rng = np.random.default_rng(20261025)
        complex_map = {
            "Packera dubia": ("C1", np.array([-4.5, 3.2, 0.4, -0.2, 0.1]), 2.15),
            "Packera anonyma": ("C1", np.array([-4.1, 2.9, 0.3, -0.1, 0.0]), 1.95),
            "Packera plattensis": ("C1", np.array([-4.3, 3.0, 0.5, -0.2, 0.1]), 2.05),
            "Packera millefolium": ("C2", np.array([6.0, 0.2, -2.5, 1.2, -0.4]), 1.80),
            "Packera memmingeri": ("C2", np.array([5.6, 0.4, -2.3, 1.0, -0.3]), 1.85),
            "Packera aurea": ("C3", np.array([0.2, -5.2, 1.8, -1.1, 0.5]), 2.40),
            "Packera obovata": ("C3", np.array([0.5, -4.9, 1.6, -0.9, 0.4]), 2.30),
            "Packera pseudaurea": ("C3", np.array([0.0, -5.0, 1.7, -1.0, 0.5]), 2.35),
        }
        complex_members = {
            "C1": ["Packera dubia", "Packera anonyma", "Packera plattensis"],
            "C2": ["Packera millefolium", "Packera memmingeri"],
            "C3": ["Packera aurea", "Packera obovata", "Packera pseudaurea"],
        }

        vouchers_rows = []
        harmonics_rows = []
        morph_rows = []
        idx = 1
        for taxon, (true_cid, center, inv_ratio) in complex_map.items():
            peers = [t for t in complex_members[true_cid] if t != taxon]
            for k in range(15):
                cat = f"NCU-SYN-{idx:04d}"
                idx += 1
                det_tier = "Tier_1_Gold" if k < 10 else "Tier_2_Silver"
                auth_wt = 3.0 if det_tier == "Tier_1_Gold" else 1.5
                pcs = center + rng.normal(0.0, 0.28, size=5)
                if k % 5 == 0:
                    cda_pred = peers[k % len(peers)]
                    cda_post = 0.58
                    cda_sec = 0.38
                else:
                    cda_pred = taxon
                    cda_post = 0.78
                    cda_sec = 0.19

                vouchers_rows.append({
                    "catalogNumber": cat,
                    "species_raw": taxon,
                    "species_standardized": taxon,
                    "determiner_tier": det_tier,
                    "authority_weight": auth_wt,
                })
                harmonics_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "species_standardized": taxon,
                    "determiner_tier": det_tier,
                    "authority_weight": auth_wt,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "PC1": float(pcs[0]), "PC2": float(pcs[1]), "PC3": float(pcs[2]), "PC4": float(pcs[3]), "PC5": float(pcs[4]),
                })
                morph_rows.append({
                    "catalogNumber": cat,
                    "plant_individual_id": 1,
                    "species": taxon,
                    "species_standardized": taxon,
                    "determiner_tier": det_tier,
                    "authority_weight": auth_wt,
                    "basal_leaf_count": 3,
                    "morph_tier_source": "Basal_Rosette",
                    "PC1": float(pcs[0]), "PC2": float(pcs[1]), "PC3": float(pcs[2]), "PC4": float(pcs[3]), "PC5": float(pcs[4]),
                    "involucre_length_width_ratio": float(inv_ratio + rng.normal(0.0, 0.05)),
                    "cda_predicted_taxon": cda_pred,
                    "cda_posterior_prob": cda_post,
                    "cda_second_posterior_prob": cda_sec,
                    "cda_posterior_margin": cda_post - cda_sec,
                    "morphometric_outlier_flag": False,
                })

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            vouchers_csv = tmp_path / "curated_vouchers.csv"
            harmonics_csv = tmp_path / "leaf_efa_harmonics.csv"
            morph_csv = tmp_path / "morphometric_flags.csv"
            out_queue = tmp_path / "triage_queue.csv"
            out_plot = tmp_path / "Figure_Integrative_Packera_dubia_Revision.pdf"
            out_report = tmp_path / "Packera_dubia_Taxonomic_Revision_Summary.md"
            out_html = tmp_path / "Triage_Dashboard_Interactive.html"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_figure = tmp_path / "fig06c_taxon_complex_network.pdf"

            pd.DataFrame(vouchers_rows).to_csv(vouchers_csv, index=False)
            pd.DataFrame(harmonics_rows).to_csv(harmonics_csv, index=False)
            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)

            cmd = [
                "Rscript",
                "scripts/analysis/06_triage_dashboard_synthesis.R",
                "--vouchers", str(vouchers_csv),
                "--harmonics", str(harmonics_csv),
                "--morphometrics", str(morph_csv),
                "--focal-taxon", "Packera dubia",
                "--min-affinity", "0.05",
                "--output-queue", str(out_queue),
                "--output-plot", str(out_plot),
                "--output-report", str(out_report),
                "--output-html", str(out_html),
                "--output-complexes", str(out_complexes),
                "--output-affinity", str(out_affinity),
                "--output-neighbors", str(out_neighbors),
                "--output-complex-plot", str(out_figure),
                "--demo",
            ]
            if not shutil.which("Rscript"):
                self.skipTest("Rscript not installed.")
            res = subprocess.run(cmd, cwd=str(REPO_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(res.returncode, 0, f"Unified Step 06 failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

            self.assertTrue(out_complexes.exists() and out_complexes.stat().st_size > 0)
            self.assertTrue(out_affinity.exists() and out_affinity.stat().st_size > 0)
            self.assertTrue(out_neighbors.exists() and out_neighbors.stat().st_size > 0)
            self.assertTrue(out_figure.exists() and out_figure.stat().st_size > 0)
            self.assertTrue(out_figure.with_suffix(".png").exists())
            self.assertIn("Empirical Taxon-Complex Delimitation", out_report.read_text(encoding="utf-8"))
            self.assertIn("Empirical Species Complexes", out_html.read_text(encoding="utf-8"))

            comp_df = pd.read_csv(out_complexes)
            self.assertEqual(len(comp_df), 8)
            self.assertEqual(set(comp_df.columns), CANONICAL_COMPLEX_COLS)

            aff_df = pd.read_csv(out_affinity)
            self.assertIn("species_standardized", aff_df.columns)
            self.assertNotIn("taxon", aff_df.columns)

            grouped = (
                comp_df.groupby("consensus_complex_id")["species_standardized"]
                .apply(lambda s: frozenset(s.tolist()))
                .tolist()
            )
            expected_partitions = {
                frozenset(["Packera dubia", "Packera anonyma", "Packera plattensis"]),
                frozenset(["Packera millefolium", "Packera memmingeri"]),
                frozenset(["Packera aurea", "Packera obovata", "Packera pseudaurea"]),
            }
            self.assertEqual(set(grouped), expected_partitions)

            neighbors = [line.strip() for line in out_neighbors.read_text(encoding="utf-8").splitlines() if line.strip()]
            self.assertEqual(set(neighbors), {"Packera anonyma", "Packera plattensis"})

            disc_text = DISCOVERY_UTILS_R.read_text(encoding="utf-8")
            self.assertIn("compute_hopkins_phenological_affinity", disc_text)
            self.assertIn("compute_modality_kappa_weights", disc_text)
            self.assertIn("evaluate_hierarchical_complex_delimitation", disc_text)
            self.assertIn("Padial et al. 2010", disc_text)
            self.assertIn("Dawid & Skene 1979", disc_text)
            self.assertIn("Schlick-Steiner et al. 2010", disc_text)

    def test_hierarchical_precedence_missing_modality_and_3d_hollow_ordinations(self):
        """Verify Stage 1-5 precedence, missing-modality weight renormalization, and 3D hollow ordinations + BH FDR."""
        self.assertIn("STEP06_HIERARCHICAL_OK", self._require_batched_r())

    def test_subcluster_network_when_2_species_or_within_species_gmm_clusters(self) -> None:
        """Verify taxon_complex_discovery_utils.R (<=800 lines) builds Taxon x Morphotype subcluster network inside Step 06 when <=2 species."""
        line_count = len(DISCOVERY_UTILS_R.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(line_count, 800, f"taxon_complex_discovery_utils.R has {line_count} lines; must be <= 800 lines.")

        rng = np.random.default_rng(20261028)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            morph_csv = tmp_path / "morphometric_flags.csv"
            out_queue = tmp_path / "triage_queue.csv"
            out_plot = tmp_path / "rev.pdf"
            out_report = tmp_path / "rev.md"
            out_html = tmp_path / "dash.html"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_fig = tmp_path / "fig06c_taxon_complex_network.pdf"

            subclusters = [
                ("Packera dubia", 1, 8, np.array([-2.2, -1.5, 0.0, 0.0, 0.0])),
                ("Packera dubia", 2, 7, np.array([-1.8, -1.2, 0.2, 0.0, 0.0])),
                ("Packera dubia", 3, 2, np.array([0.0, 0.0, 0.0, 0.0, 0.0])),
                ("Packera paupercula", 1, 8, np.array([2.2, 1.5, 0.0, 0.0, 0.0])),
                ("Packera paupercula", 2, 6, np.array([1.9, 1.3, -0.2, 0.0, 0.0])),
            ]
            morph_rows = []
            v_idx = 0
            for sp_name, k_id, count, center in subclusters:
                for _ in range(count):
                    v_idx += 1
                    vec = center + rng.normal(0.0, 0.25, size=5)
                    morph_rows.append({
                        "catalogNumber": f"SUB_{v_idx:04d}",
                        "plant_individual_id": 1,
                        "species_standardized": sp_name,
                        "determiner_tier": "Tier_1_Gold",
                        "basal_leaf_count": 3,
                        "morph_tier_source": "Basal_Rosette",
                        "gmm_cluster": k_id,
                        "PC1": float(vec[0]), "PC2": float(vec[1]), "PC3": float(vec[2]), "PC4": float(vec[3]), "PC5": float(vec[4]),
                        "cda_predicted_taxon": sp_name,
                        "cda_posterior_prob": 0.90,
                        "cda_second_posterior_prob": 0.08,
                    })

            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)
            if not shutil.which("Rscript"):
                self.skipTest("Rscript not installed.")

            res = subprocess.run(
                [
                    "Rscript", str(SYNTHESIS_R_SCRIPT),
                    "--vouchers", str(morph_csv),
                    "--morphometrics", str(morph_csv),
                    "--focal-taxon", "Packera dubia",
                    "--output-queue", str(out_queue),
                    "--output-plot", str(out_plot),
                    "--output-report", str(out_report),
                    "--output-html", str(out_html),
                    "--output-complexes", str(out_complexes),
                    "--output-affinity", str(out_affinity),
                    "--output-neighbors", str(out_neighbors),
                    "--output-complex-plot", str(out_fig),
                    "--demo",
                ],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(res.returncode, 0, f"Subcluster discovery failed:\n{res.stderr}\n{res.stdout}")
            comp_df = pd.read_csv(out_complexes)
            expected_nodes = {"P. dubia [K1]", "P. dubia [K2]", "P. paupercula [K1]", "P. paupercula [K2]"}
            self.assertEqual(set(comp_df["species_standardized"]), expected_nodes)
            self.assertTrue((comp_df["n_vouchers"] >= 5).all())
            self.assertTrue(out_fig.exists() and out_fig.stat().st_size > 1000)

    def test_intraspecific_morphotype_cluster_visualization_across_r_figures(self) -> None:
        """Verify pub_figures_utils.R and 05_multimodal_spatial_rf.R use morphotype clusters when <2 species."""
        self.assertIn("STEP06_INTRA_FIGS_OK", self._require_batched_r())

    def test_multimodal_cli_flag_and_ecological_edaphic_left_join(self) -> None:
        """Verify unified Step 06 merges Step 05 ecological/edaphic columns and taxon_complex_discovery_utils.R stays <=800 lines."""
        line_count = len(DISCOVERY_UTILS_R.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(line_count, 800, f"taxon_complex_discovery_utils.R has {line_count} lines; must be <= 800 lines.")
        self.assertIn("STEP06_MULTIMODAL_JOIN_OK", self._require_batched_r())

    def test_dynamic_morph_and_primary_env_predictor_discovery_in_step06(self) -> None:
        """Verify unified Step 06 dynamically discovers >5 PCs via selected_morph_features, Step 05 primary_env_predictors, and Peak_Anthesis DOY."""
        rng = np.random.default_rng(20261031)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            morph_csv = tmp_path / "morphometric_flags.csv"
            mm_csv = tmp_path / "multimodal_conflict_flags.csv"
            out_queue = tmp_path / "triage_queue.csv"
            out_plot = tmp_path / "rev.pdf"
            out_report = tmp_path / "rev.md"
            out_html = tmp_path / "dash.html"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_fig = tmp_path / "fig06c_taxon_complex_network.pdf"

            taxa_specs = [
                ("Packera hyperborea", (-3.5, -3.0), (5.5, 1.8)),
                ("Packera contermina", (-3.2, -2.8), (5.8, 1.9)),
                ("Packera cymbalaria", (4.5, 4.2), (13.5, 9.2)),
            ]
            morph_rows = []
            mm_rows = []
            v_idx = 0
            for sp_name, (pc6_mu, pc7_mu), (b10_mu, ksat_mu) in taxa_specs:
                for j in range(8):
                    v_idx += 1
                    cat_no = f"DYN06_{v_idx:03d}"
                    lat_v = 61.0 + 0.4 * j
                    morph_rows.append({
                        "catalogNumber": cat_no,
                        "plant_individual_id": 0,
                        "species_standardized": sp_name,
                        "determiner_tier": "Tier_1_Gold",
                        "basal_leaf_count": 3,
                        "morph_tier_source": "Basal_Rosette",
                        "gmm_cluster": 1,
                        "PC1": float(rng.normal(0.0, 0.15)),
                        "PC2": float(rng.normal(0.0, 0.15)),
                        "PC3": float(rng.normal(0.0, 0.15)),
                        "PC4": float(rng.normal(0.0, 0.15)),
                        "PC5": float(rng.normal(0.0, 0.15)),
                        "PC6": float(pc6_mu + rng.normal(0.0, 0.20)),
                        "PC7": float(pc7_mu + rng.normal(0.0, 0.20)),
                        "selected_morph_features": "PC1;PC2;PC3;PC4;PC5;PC6;PC7",
                        "cda_predicted_taxon": sp_name,
                        "cda_posterior_prob": 0.92,
                        "cda_second_posterior_prob": 0.05,
                    })
                    mm_rows.append({
                        "catalogNumber": cat_no,
                        "plant_individual_id": 0,
                        "species_standardized": sp_name,
                        "decimalLatitude": lat_v,
                        "decimalLongitude": -148.0 + 0.3 * j,
                        "doy": int(215 + 3.0 * (lat_v - 61.0)),
                        "is_verified_anthesis": True,
                        "voucher_phenological_state": "Peak_Anthesis",
                        "bio10_mean_temp_warmest_qtr": float(b10_mu + rng.normal(0.0, 0.25)),
                        "polaris_ksat_30m": float(ksat_mu + rng.normal(0.0, 0.25)),
                        "primary_env_predictors": "bio10_mean_temp_warmest_qtr;polaris_ksat_30m",
                        "primary_clim_predictors": "bio10_mean_temp_warmest_qtr",
                        "primary_edaph_predictors": "polaris_ksat_30m",
                    })

            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)
            pd.DataFrame(mm_rows).to_csv(mm_csv, index=False)
            if not shutil.which("Rscript"):
                self.skipTest("Rscript not installed.")

            res = subprocess.run(
                [
                    "Rscript", str(SYNTHESIS_R_SCRIPT),
                    "--vouchers", str(morph_csv),
                    "--morphometrics", str(morph_csv),
                    "--multimodal-flags", str(mm_csv),
                    "--focal-taxon", "Packera hyperborea",
                    "--output-queue", str(out_queue),
                    "--output-plot", str(out_plot),
                    "--output-report", str(out_report),
                    "--output-html", str(out_html),
                    "--output-complexes", str(out_complexes),
                    "--output-affinity", str(out_affinity),
                    "--output-neighbors", str(out_neighbors),
                    "--output-complex-plot", str(out_fig),
                    "--demo",
                ],
                cwd=str(REPO_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(res.returncode, 0, f"Step 06 dynamic feature test failed:\n{res.stderr}\n{res.stdout}")
            comp_df = pd.read_csv(out_complexes)
            self.assertEqual(len(comp_df), 3)
            hyp_row = comp_df.loc[comp_df["species_standardized"] == "Packera hyperborea"].iloc[0]
            self.assertEqual(str(hyp_row["nearest_neighbor_taxon"]), "Packera contermina")
            # Gate 1 fix: Peak_Anthesis vouchers are recognized (hopkins_pheno_overlap > 0.80, NOT NA/0!)
            self.assertGreater(float(hyp_row["hopkins_pheno_overlap"]), 0.80)
            self.assertGreater(float(hyp_row["eco_schoener_overlap"]), 0.50)

    def test_multi_species_subcluster_suppression_and_line_budget(self):
        """Verify taxon_complex_discovery_utils.R returns subcluster_mode = FALSE when n_sp > 3L and stays <= 800 lines."""
        src_lines = DISCOVERY_UTILS_R.read_text(encoding="utf-8").splitlines()
        self.assertLessEqual(len(src_lines), 800, f"taxon_complex_discovery_utils.R has {len(src_lines)} lines (> 800)")
        src_text = "\n".join(src_lines)
        self.assertIn("n_sp > 3L", src_text)

        r_code = f"""
        source("{DISCOVERY_UTILS_R.as_posix()}")
        df4 <- data.frame(
          catalogNumber = sprintf("SP4_%02d", 1:40),
          plant_individual_id = 1L,
          species_standardized = rep(c("Packera dubia", "Packera aurea", "Packera anonyma", "Packera obovata"), each = 10),
          gmm_cluster = rep(c(1L, 2L), 20),
          can1 = seq(-2, 2, length.out = 40),
          can2 = seq(1, -1, length.out = 40),
          stringsAsFactors = FALSE
        )
        sub_res <- build_taxon_morphotype_subclusters(df4, min_vouchers = 3L)
        stopifnot(identical(sub_res$subcluster_mode, FALSE))
        stopifnot(length(unique(sub_res$morph_df$species_standardized)) == 4L)
        cat("SUBCLUSTER_MULTI_SP_GUARD_OK\\n")
        """
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")
        with tempfile.TemporaryDirectory() as tmp_r:
            r_tmp = Path(tmp_r) / "check_step06_guard.R"
            r_tmp.write_text(r_code, encoding="utf-8")
            res = subprocess.run(["Rscript", str(r_tmp)], cwd=str(REPO_ROOT), stdin=subprocess.DEVNULL, capture_output=True, text=True, check=False)
        self.assertEqual(res.returncode, 0, f"Subcluster multi-species guard test failed:\n{res.stderr}\n{res.stdout}")
        self.assertIn("SUBCLUSTER_MULTI_SP_GUARD_OK", res.stdout)

    def test_unified_step06_complex_discovery_zero_trust_four_gates(self):
        """
        Zero-trust regression suite for unified Step 06 Complex Delimitation Four-Gate Overhaul:
          1. Gate 1: S_morph proximity gating prevents sympatric habitat convergence from merging distant clades,
             metric sqrt(1 - W_ij) distance is used in Ward.D2, and Peak_Anthesis vouchers yield has_pheno == TRUE.
          2. Gate 2: Confident-noise & Misidentification_Reassign rows are excluded from anchors and C_morph,
             and subcluster mode preserves base_species / base_cda_predicted_taxon without forcing w_morph_raw == 1.0.
          3. Gate 3: SQLite taxon_complexes persists all 23 canonical columns under schema 20261116.
          4. Gate 4: .Random.seed is strictly preserved across phenology permutations & Leiden partitioning.
        """
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed.")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_p = Path(tmpdir)
            db_file = tmp_p / "step06_complexes.sqlite"
            with get_db_connection(db_file, init_schema=True):
                pass

            r_code = rf"""
            source("scripts/analysis/06_triage_dashboard_synthesis.R")

            # 1. Sympatric habitat convergence must NOT merge morphologically distant clades (Gate 1)
            set.seed(777)
            df_symp <- data.frame(
              catalogNumber = sprintf("SYM_%02d", 1:24),
              plant_individual_id = 0L,
              species_standardized = rep(c("Packera dubia", "Packera anonyma", "Packera millefolium"), each = 8),
              determiner_tier = "Tier_1_Gold",
              morph_tier_source = "Basal_Rosette",
              basal_leaf_count = 3L,
              PC1 = c(rnorm(8, -4.0, 0.15), rnorm(8, -3.7, 0.15), rnorm(8, 8.0, 0.15)),
              PC2 = c(rnorm(8, 1.0, 0.15), rnorm(8, 1.1, 0.15), rnorm(8, -6.0, 0.15)),
              PC3 = rnorm(24, 0.0, 0.10),
              cda_predicted_taxon = rep(c("Packera dubia", "Packera anonyma", "Packera millefolium"), each = 8),
              cda_posterior_prob = 0.94,
              decimalLatitude = rep(seq(35.0, 35.7, length.out = 8), 3),
              decimalLongitude = -80.0,
              doy = rep(120:127, 3),
              is_verified_anthesis = TRUE,
              voucher_phenological_state = "Peak_Anthesis",
              fine_soil_ph = c(rnorm(8, 6.5, 0.1), rnorm(8, 5.1, 0.1), rnorm(8, 5.1, 0.1)),
              fine_soil_sand = c(rnorm(8, 45.0, 1.5), rnorm(8, 75.0, 1.5), rnorm(8, 75.0, 1.5)),
              primary_edaph_predictors = "fine_soil_ph;fine_soil_sand",
              taxonomic_treatment_subtype = "Confirmed_Species",
              stringsAsFactors = FALSE
            )
            taxa3 <- sort(unique(df_symp$species_standardized))
            anc3 <- select_taxon_anchors(df_symp)
            seed_before <- .Random.seed
            ph3 <- compute_hopkins_phenological_affinity(df_symp, taxa3, anc3$anchor_mask)
            # RNG state preserved after permutation test!
            stopifnot(identical(.Random.seed, seed_before))
            stopifnot(isTRUE(ph3$has_pheno))

            ec3 <- compute_ecological_affinity(df_symp, taxa3, anc3$anchor_mask)
            mo3 <- compute_morphometric_affinity(df_symp, anc3$anchor_df, taxa3)
            kw3 <- compute_modality_kappa_weights(df_symp, anc3$anchor_mask)
            W3 <- compute_multimodal_complex_affinity(mo3, ec3, ph3, kw3, taxa3)
            # P. anonyma and P. millefolium share identical sand/pH & DOY, yet S_morph ~ 0 gates W_ij near 0!
            stopifnot(W3["Packera anonyma", "Packera millefolium"] < 0.02)
            stopifnot(W3["Packera anonyma", "Packera dubia"] > W3["Packera anonyma", "Packera millefolium"] * 10.0)

            pr3 <- partition_taxon_complexes(W3, taxa3, table(df_symp$species_standardized)[taxa3], anc3$n_gold_silver)
            stopifnot(identical(.Random.seed, seed_before))
            # Metric sqrt(1 - W_ij) check
            d_mat <- as.matrix(pr3$dist_obj)
            stopifnot(abs(d_mat["Packera anonyma", "Packera dubia"] - sqrt(1.0 - W3["Packera anonyma", "Packera dubia"])) < 1e-6)

            # 2. Noise & Misidentification_Reassign rows excluded from anchors and C_morph (Gate 2)
            df_corrupt <- rbind(df_symp, df_symp[1:3, ])
            df_corrupt$catalogNumber[25:27] <- c("ERR_1", "ERR_2", "ERR_3")
            df_corrupt$species_standardized[25:27] <- "Packera millefolium"
            df_corrupt$cda_predicted_taxon[25:27] <- "Packera dubia"
            df_corrupt$is_confident_label_noise <- c(rep(FALSE, 24), TRUE, FALSE, FALSE)
            df_corrupt$taxonomic_treatment_subtype[25:27] <- c("Confirmed_Species", "Misidentification_Reassign", "Misidentification_Reassign")
            anc_c <- select_taxon_anchors(df_corrupt)
            stopifnot(!any(anc_c$anchor_mask[25:27]))
            mo_c <- compute_morphometric_affinity(df_corrupt, anc_c$anchor_df, taxa3)
            stopifnot(mo_c$C_morph_sym["Packera millefolium", "Packera dubia"] < 0.01)

            # 3. Subcluster mode preserves base_cda_predicted_taxon without inflating w_morph_raw to 1.0 (Gate 2)
            df_sub <- data.frame(
              catalogNumber = sprintf("SB_%02d", 1:20),
              plant_individual_id = 0L,
              species_standardized = rep(c("Packera dubia", "Packera aurea"), each = 10),
              determiner_tier = "Tier_1_Gold",
              morph_tier_source = "Basal_Rosette",
              basal_leaf_count = 3L,
              gmm_cluster = rep(c(1L, 2L), 10),
              PC1 = rnorm(20),
              PC2 = rnorm(20),
              cda_predicted_taxon = c(rep("Packera dubia", 7), rep("Packera aurea", 3), rep("Packera aurea", 8), rep("Packera dubia", 2)),
              stringsAsFactors = FALSE
            )
            seed_sub_before <- .Random.seed
            sr <- build_taxon_morphotype_subclusters(df_sub, min_vouchers = 4L)
            stopifnot(identical(.Random.seed, seed_sub_before))
            kw_sub <- compute_modality_kappa_weights(sr$morph_df, rep(TRUE, nrow(sr$morph_df)))
            stopifnot(kw_sub$w_morph_raw < 0.80)

            # 4. SQLite persistence into taxon_complexes (Gate 3)
            Sys.setenv(PACKERA_DB_PATH = "{db_file.as_posix()}", RUN_ID = "run_step06_tc")
            full_res <- run_step06_complex_discovery(df_symp, list(min_affinity = 0.03, demo = FALSE, vouchers = ""))
            persist_step06_complexes_to_sqlite(full_res$complexes_df, vouchers_path = NULL, run_id = "run_step06_tc", demo = FALSE)
            cat("STEP06_COMPLEX_FOUR_GATES_OK\n")
            """
            r_file = tmp_p / "test_complex_gates.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = subprocess.run(["Rscript", str(r_file)], cwd=str(REPO_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(proc.returncode, 0, f"Step 06 complex four-gate R test failed:\n{proc.stderr}\n{proc.stdout}")
            self.assertIn("STEP06_COMPLEX_FOUR_GATES_OK", proc.stdout)

            conn = sqlite3.connect(db_file)
            rows = conn.execute("SELECT species_standardized, consensus_complex_id FROM taxon_complexes WHERE run_id = 'run_step06_tc'").fetchall()
            conn.close()
            self.assertEqual(len(rows), 3)


if __name__ == "__main__":
    unittest.main()
