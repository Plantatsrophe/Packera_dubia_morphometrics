#!/usr/bin/env python3
"""
Unit tests for Eastern North American Packera N-Taxon Standardization,
Bipinnatisect Dissection Gate, and Quantitative Taxon-Complex Discovery Engine
(08_discover_taxon_complexes.R).
"""

import importlib
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import cv2
import numpy as np
import pandas as pd

from scripts.core.harvester import standardize_packera_taxon
from scripts.vision.lm2_geometry_utils import (
    classify_and_route_leaf,
    evaluate_leaf_quality,
    is_botanical_dissection,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DISCOVERY_R_SCRIPT = REPO_ROOT / "scripts" / "analysis" / "08_discover_taxon_complexes.R"
PUB_FIG_R_SCRIPT = REPO_ROOT / "scripts" / "morphometrics" / "pub_figures_utils.R"


class TestEasternPackeraAndComplexDiscovery(unittest.TestCase):
    """Tests for Milestones 1-5: N-taxon basionym mapping, dissection gate, and complex discovery."""

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

            pd.DataFrame([
                {"catalogNumber": "MM_01", "plant_individual_id": 1, "species_standardized": "Packera aurea", "determiner_tier": "Tier_1_Gold", "n_tier1_leaves": 3, "Sym_PC1": -1.2, "Sym_PC2": 0.3},
                {"catalogNumber": "MM_02", "plant_individual_id": 1, "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold", "n_tier1_leaves": 3, "Sym_PC1": 1.4, "Sym_PC2": -0.2},
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
              cat("STEP08_R_PALETTE_OK\\n")
            }})

            source("{DISCOVERY_R_SCRIPT.as_posix()}")

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
                cda_posterior_prob = c(0.95, 0.62, 0.65, 0.44, 0.48, 0.92),
                cda_second_posterior_prob = c(0.03, 0.36, 0.30, 0.28, 0.25, 0.05),
                cda_posterior_margin = c(0.92, 0.26, 0.35, 0.16, 0.23, 0.87),
                morphometric_outlier_flag = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE),
                gmm_cluster = c(1L, 2L, 2L, 5L, 5L, 2L),
                gmm_2_delta_bic = c(0, 0, 0, 12.0, 12.0, 0),
                edaphic_schoener_d = c(0.75, 0.60, 0.55, 0.18, 0.72, 0.80),
                n_tier1_leaves = rep(3L, 6),
                stringsAsFactors = FALSE
              )
              kw_missing <- compute_modality_kappa_weights(mdf, NULL, rep(TRUE, 6))
              stopifnot(abs(kw_missing$w_morph - 1.0) < 1e-9 && kw_missing$w_vis == 0.0)

              vdf <- data.frame(
                catalogNumber = c("V1", "V2", "V3", "V4", "V5", "V6"),
                species_standardized = mdf$species_standardized,
                determiner_tier = rep("Tier_1_Gold", 6),
                vision_predicted_label = c("Packera dubia", "Packera anonyma", "Packera paupercula", "Packera anonyma", "Packera anonyma", "Packera anonyma"),
                c_error = c(0.04, 0.82, 0.75, 0.28, 0.25, 0.05),
                stringsAsFactors = FALSE
              )
              pheno_mock <- list(delta_doy_resid = c(2.0, 3.0, 4.0, 22.0, 3.0, 1.0))
              calls <- evaluate_hierarchical_complex_delimitation(mdf, vdf, pheno_mock)
              stopifnot(calls[1] == "Confirmed_Species")
              stopifnot(calls[2] == "Misidentification_Reassign")
              stopifnot(calls[3] == "Hybrid_Introgressant")
              stopifnot(calls[4] == "Putative_Cryptic_Species")
              stopifnot(calls[5] == "Ecophenotypic_Plasticity")
              cat("STEP08_HIERARCHICAL_OK\\n")
            }})

            # Block 3: Multimodal CLI flag & ecological/edaphic left join
            local({{
              assign("commandArgs", function(trailingOnly = TRUE) c("--morphometrics", "{morph_csv.as_posix()}", "--multimodal", "{mm_csv.as_posix()}"), envir = .GlobalEnv)
              on.exit(rm("commandArgs", envir = .GlobalEnv), add = TRUE)
              opts <- parse_discovery_args()
              stopifnot(opts$multimodal == "{mm_csv.as_posix()}")
              inp <- load_and_merge_inputs(opts)
              mdf <- inp$morph_df
              req <- c("fine_soil_ph", "fine_soil_sand", "fine_bedrock_depth_cm", "ssurgo_rock_outcrop_pct",
                       "soil_ph", "soil_sand", "soil_cec", "soil_bulk_density",
                       "bio1_temp_mean", "bio4_temp_seasonality", "bio12_precip_annual", "bio15_precip_seasonality",
                       "edaphic_schoener_d", "schoener_d")
              stopifnot(all(req %in% names(mdf)))
              stopifnot(abs(mdf$fine_soil_ph[mdf$catalogNumber == "MM_01"] - 5.9) < 1e-6)
              stopifnot(abs(mdf$ssurgo_rock_outcrop_pct[mdf$catalogNumber == "MM_02"] - 24.0) < 1e-6)
              cat("STEP08_MULTIMODAL_JOIN_OK\\n")
            }})

            # Block 4: Intraspecific morphotype cluster visualization across R figures
            source("{PUB_FIG_R_SCRIPT.as_posix()}")
            source("{(REPO_ROOT / 'scripts' / 'analysis' / '06_multimodal_spatial_rf.R').as_posix()}")
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
                Sym_PC1 = c(rnorm(8, -2), rnorm(8, 0), rnorm(8, 2)),
                Sym_PC2 = c(rnorm(8, 1), rnorm(8, -1.5), rnorm(8, 1)),
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
                latitude = runif(n, 32, 38),
                longitude = runif(n, -84, -77),
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

              f_s07 <- file.path("{intra_fig_dir.as_posix()}", "integrative_taxonomic_synthesis.pdf")
              render_step07_synthesis_plate_pdf(df, f_s07)
              stopifnot(file.exists(f_s07) && file.info(f_s07)$size > 1000)
              cat("STEP08_INTRA_FIGS_OK\\n")
            }})
            """
            r_batch_file = cls._class_tmp_path / "batched_step08_checks.R"
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
            f"Batched Step 08 R unit checks failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_standardize_packera_taxon_all_eastern_species_and_basionyms(self):
        """Verify North American Packera species (Eastern, Western, Boreal, Mexican) map cleanly."""
        cases = {
            # Eastern North American Packera species, varieties, hybrids, and authority strings
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
            # Western North American Packera species
            "Packera bernardina (Greene) W.A.Weber & Á.Löve": "Packera bernardina",
            "Packera cana (Hook.) W.A.Weber & Á.Löve": "Packera cana",
            "Packera fendleri (A.Gray) W.A.Weber & Á.Löve": "Packera fendleri",
            "Packera multilobata (Torr. & A.Gray ex A.Gray) W.A.Weber & Á.Löve": "Packera multilobata",
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve": "Packera streptanthifolia",
            "Packera werneriifolia (A.Gray) W.A.Weber & Á.Löve": "Packera werneriifolia",
            # Boreal & Arctic North American Packera species
            "Packera cymbalaria (Pursh) W.A.Weber & Á.Löve": "Packera cymbalaria",
            "Packera hyperborealis (Greenm.) Á.Löve & D.Löve": "Packera hyperborealis",
            "Packera indecora (Greene) Á.Löve & D.Löve": "Packera indecora",
            "Packera ogotorukensis (Packer) Á.Löve & D.Löve": "Packera ogotorukensis",
            "Packera pauciflora (Pursh) Á.Löve & D.Löve": "Packera pauciflora",
            # Mexican & Madrean Packera species
            "Packera bellidifolia (Kunth) W.A.Weber & Á.Löve": "Packera bellidifolia",
            "Packera coahuilensis (Greenm.) C.Jeffrey": "Packera coahuilensis",
            "Packera sanguisorbae (DC.) C.Jeffrey": "Packera sanguisorbae",
            "Packera toluccana (DC.) W.A.Weber & Á.Löve": "Packera toluccana",
            "Packera zimapanica (Hemsl.) C.C.Freeman & T.M.Barkley": "Packera zimapanica",
        }

        mod01_pipe = importlib.import_module("scripts.pipeline.01_voucher_harvester")
        mod05 = importlib.import_module("scripts.analysis.05_cleanlab_vision_xai")

        for raw_name, expected in cases.items():
            self.assertEqual(standardize_packera_taxon(raw_name), expected, f"Failed for {raw_name!r}")
            self.assertEqual(mod01_pipe.standardize_packera_taxon(raw_name), expected)
            self.assertEqual(mod05.standardize_packera_taxon(raw_name), expected)

        # Also verify dynamic GBIF backbone resolution for historical Senecio basionyms across regions
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
                self.assertEqual(mod05.standardize_packera_taxon(raw_syn), expected)

    def test_r_standardize_taxon_and_dynamic_palette(self):
        """Verify R standardize_taxon_name(), get_dynamic_packera_palette(), shapes, titles, and italic expressions are 100% species-agnostic."""
        self.assertIn("STEP08_R_PALETTE_OK", self._require_batched_r())

    def test_bipinnatisect_dissection_gate_tier1_routing(self):
        """Verify finely dissected/bipinnatisect leaf (0.18 <= solidity < 0.48, >=4 deep sinuses) routes to Tier_1_Direct."""
        mask = np.zeros((400, 400), dtype=np.uint8)
        # Central rachis along x = 200, y in [60, 340], width = 10 px
        cv2.rectangle(mask, (195, 60), (205, 340), 255, -1)
        # 5 pairs of symmetric narrow pinnae creating >= 8 deep bilateral sinuses
        # with solidity in [0.18, 0.48] and high bilateral symmetry (UCS >= 0.85)
        pinna_ys = [85, 135, 185, 235, 285]
        half_spans = [75, 95, 105, 90, 65]
        for y_c, span in zip(pinna_ys, half_spans):
            cv2.rectangle(mask, (200 - span, y_c - 6), (200 + span, y_c + 6), 255, -1)
        # Apical lobe
        cv2.rectangle(mask, (182, 45), (218, 65), 255, -1)

        cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        self.assertTrue(len(cnts) > 0)
        cnt = max(cnts, key=cv2.contourArea)
        area = cv2.contourArea(cnt)
        hull = cv2.convexHull(cnt)
        solidity = area / cv2.contourArea(hull)
        self.assertGreaterEqual(solidity, 0.18)
        self.assertLess(solidity, 0.48)

        self.assertTrue(
            is_botanical_dissection(mask, contour=cnt, solidity=solidity),
            f"Expected bipinnatisect silhouette (solidity={solidity:.3f}) to be recognized as botanical dissection.",
        )

        res = evaluate_leaf_quality(mask, contour=cnt)
        self.assertTrue(res["is_dissected"])
        self.assertEqual(res["assigned_tier"], "tier1")
        self.assertIn("Tier 1", res["tier"])

        route_res = classify_and_route_leaf(mask)
        self.assertTrue(route_res["is_dissected"])
        self.assertEqual(route_res["assigned_tier"], "tier1")

    def test_discover_taxon_complexes_8_species_3_complexes(self):
        """Verify 08_discover_taxon_complexes.R recovers 3 true complexes across 8 species with 100% accuracy."""
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
        vision_rows = []

        n_per_taxon = 15
        idx = 1
        for taxon, (true_cid, center, inv_ratio) in complex_map.items():
            peers = [t for t in complex_members[true_cid] if t != taxon]
            for k in range(n_per_taxon):
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

                if k % 4 == 0:
                    vis_pred = peers[(k // 4) % len(peers)]
                    c_err = 0.62
                else:
                    vis_pred = taxon
                    c_err = 0.12

                vouchers_rows.append(
                    {
                        "catalogNumber": cat,
                        "species_raw": taxon,
                        "species_standardized": taxon,
                        "determiner_tier": det_tier,
                        "authority_weight": auth_wt,
                    }
                )
                harmonics_rows.append(
                    {
                        "catalogNumber": cat,
                        "plant_individual_id": 1,
                        "species_standardized": taxon,
                        "determiner_tier": det_tier,
                        "authority_weight": auth_wt,
                        "n_tier1_leaves": 3,
                        "n_tier2_leaves": 0,
                        "morph_tier_source": "Tier_1_Pure",
                        "Sym_PC1": float(pcs[0]),
                        "Sym_PC2": float(pcs[1]),
                        "Sym_PC3": float(pcs[2]),
                        "Sym_PC4": float(pcs[3]),
                        "Sym_PC5": float(pcs[4]),
                        "PC1": float(pcs[0]),
                        "PC2": float(pcs[1]),
                        "PC3": float(pcs[2]),
                        "PC4": float(pcs[3]),
                        "PC5": float(pcs[4]),
                    }
                )
                morph_rows.append(
                    {
                        "catalogNumber": cat,
                        "plant_individual_id": 1,
                        "species": taxon,
                        "species_standardized": taxon,
                        "determiner_tier": det_tier,
                        "authority_weight": auth_wt,
                        "n_tier1_leaves": 3,
                        "morph_tier_source": "Tier_1_Pure",
                        "Sym_PC1": float(pcs[0]),
                        "Sym_PC2": float(pcs[1]),
                        "Sym_PC3": float(pcs[2]),
                        "Sym_PC4": float(pcs[3]),
                        "Sym_PC5": float(pcs[4]),
                        "involucre_length_width_ratio": float(inv_ratio + rng.normal(0.0, 0.05)),
                        "cda_predicted_taxon": cda_pred,
                        "cda_posterior_prob": cda_post,
                        "cda_second_posterior_prob": cda_sec,
                        "cda_posterior_margin": cda_post - cda_sec,
                        "morphometric_outlier_flag": False,
                    }
                )
                vision_rows.append(
                    {
                        "catalogNumber": cat,
                        "given_label": taxon,
                        "predicted_label": vis_pred,
                        "vision_predicted_label": vis_pred,
                        "c_error": c_err,
                        "is_label_corrupted": c_err > 0.5,
                    }
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            vouchers_csv = tmp_path / "curated_vouchers.csv"
            harmonics_csv = tmp_path / "leaf_efa_harmonics.csv"
            morph_csv = tmp_path / "morphometric_flags.csv"
            vision_csv = tmp_path / "label_noise_audit.csv"

            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_figure = tmp_path / "fig_empirical_taxon_complex_network.pdf"

            pd.DataFrame(vouchers_rows).to_csv(vouchers_csv, index=False)
            pd.DataFrame(harmonics_rows).to_csv(harmonics_csv, index=False)
            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)
            pd.DataFrame(vision_rows).to_csv(vision_csv, index=False)

            cmd = [
                "Rscript",
                "scripts/analysis/08_discover_taxon_complexes.R",
                "--morphometrics",
                str(morph_csv),
                "--harmonics",
                str(harmonics_csv),
                "--vision-audit",
                str(vision_csv),
                "--vouchers",
                str(vouchers_csv),
                "--focal-taxon",
                "Packera dubia",
                "--min-affinity",
                "0.05",
                "--output-complexes",
                str(out_complexes),
                "--output-affinity",
                str(out_affinity),
                "--output-neighbors",
                str(out_neighbors),
                "--output-figure",
                str(out_figure),
            ]
            res = subprocess.run(cmd, cwd=str(REPO_ROOT), capture_output=True, text=True, check=False)
            self.assertEqual(
                res.returncode,
                0,
                f"08_discover_taxon_complexes.R failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}",
            )

            self.assertTrue(out_complexes.exists() and out_complexes.stat().st_size > 0)
            self.assertTrue(out_affinity.exists() and out_affinity.stat().st_size > 0)
            self.assertTrue(out_neighbors.exists() and out_neighbors.stat().st_size > 0)
            self.assertTrue(out_figure.exists() and out_figure.stat().st_size > 0)
            self.assertTrue(out_figure.with_suffix(".png").exists())

            comp_df = pd.read_csv(out_complexes)
            self.assertEqual(len(comp_df), 8)
            for req_col in (
                "taxon",
                "n_Pure_Gold_Silver",
                "leiden_complex_id",
                "ward_complex_id",
                "consensus_complex_id",
                "complex_label",
                "within_complex_mean_affinity",
                "nearest_neighbor_taxon",
                "nearest_neighbor_affinity",
            ):
                self.assertIn(req_col, comp_df.columns)

            grouped = (
                comp_df.groupby("consensus_complex_id")["taxon"]
                .apply(lambda s: frozenset(s.tolist()))
                .tolist()
            )
            expected_partitions = {
                frozenset(["Packera dubia", "Packera anonyma", "Packera plattensis"]),
                frozenset(["Packera millefolium", "Packera memmingeri"]),
                frozenset(["Packera aurea", "Packera obovata", "Packera pseudaurea"]),
            }
            self.assertEqual(set(grouped), expected_partitions)

            neighbors = [
                line.strip()
                for line in out_neighbors.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertEqual(set(neighbors), {"Packera anonyma", "Packera plattensis"})

            self.assertIn("hopkins_pheno_overlap", comp_df.columns)
            self.assertIn("integrative_congruence", comp_df.columns)
            self.assertIn("kappa_morph_weight", comp_df.columns)
            self.assertIn("kappa_vis_weight", comp_df.columns)
            r08_text = DISCOVERY_R_SCRIPT.read_text(encoding="utf-8")
            self.assertIn("compute_hopkins_phenological_affinity", r08_text)
            self.assertIn("compute_modality_kappa_weights", r08_text)
            self.assertIn("evaluate_hierarchical_complex_delimitation", r08_text)
            self.assertIn("Padial et al. 2010", r08_text)
            self.assertIn("Dawid & Skene 1979", r08_text)
            self.assertIn("Schlick-Steiner et al. 2010", r08_text)

    def test_hierarchical_precedence_missing_modality_and_3d_hollow_ordinations(self):
        """Verify Stage 1-5 precedence in 08_discover_taxon_complexes.R, missing-modality weight renormalization, and 3D hollow ordinations + BH FDR in pub_figures_utils.R."""
        self.assertIn("STEP08_HIERARCHICAL_OK", self._require_batched_r())

    def test_subcluster_network_when_2_species_or_within_species_gmm_clusters(self) -> None:
        """Verify 08_discover_taxon_complexes.R (<=800 lines) builds Taxon x Morphotype subcluster network when <=2 species."""
        line_count = len(DISCOVERY_R_SCRIPT.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(
            line_count,
            800,
            f"08_discover_taxon_complexes.R has {line_count} lines; must be <= 800 lines.",
        )

        rng = np.random.default_rng(20261028)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            morph_csv = tmp_path / "morphometric_flags.csv"
            vis_csv = tmp_path / "label_noise_audit.csv"
            out_complexes = tmp_path / "empirical_taxon_complexes.csv"
            out_affinity = tmp_path / "pairwise_taxon_affinity_matrix.csv"
            out_neighbors = tmp_path / "focal_auto_neighbors.txt"
            out_fig = tmp_path / "fig_empirical_taxon_complex_network.pdf"

            subclusters = [
                ("Packera dubia", 1, 8, np.array([-2.2, -1.5, 0.0, 0.0, 0.0])),
                ("Packera dubia", 2, 7, np.array([-1.8, -1.2, 0.2, 0.0, 0.0])),
                ("Packera dubia", 3, 2, np.array([0.0, 0.0, 0.0, 0.0, 0.0])),
                ("Packera paupercula", 1, 8, np.array([2.2, 1.5, 0.0, 0.0, 0.0])),
                ("Packera paupercula", 2, 6, np.array([1.9, 1.3, -0.2, 0.0, 0.0])),
            ]
            morph_rows = []
            vis_rows = []
            v_idx = 0
            for sp_name, k_id, count, center in subclusters:
                for _ in range(count):
                    v_idx += 1
                    cat_no = f"SUB_{v_idx:04d}"
                    vec = center + rng.normal(0.0, 0.25, size=5)
                    morph_rows.append({
                        "catalogNumber": cat_no,
                        "plant_individual_id": 1,
                        "species_standardized": sp_name,
                        "determiner_tier": "Tier_1_Gold",
                        "n_tier1_leaves": 3,
                        "morph_tier_source": "Tier_1_Pure",
                        "gmm_cluster": k_id,
                        "Sym_PC1": float(vec[0]),
                        "Sym_PC2": float(vec[1]),
                        "Sym_PC3": float(vec[2]),
                        "Sym_PC4": float(vec[3]),
                        "Sym_PC5": float(vec[4]),
                        "cda_predicted_taxon": sp_name,
                        "cda_posterior_prob": 0.90,
                        "cda_second_posterior_prob": 0.08,
                    })
                    vis_rows.append({
                        "catalogNumber": cat_no,
                        "plant_individual_id": 1,
                        "species_standardized": sp_name,
                        "determiner_tier": "Tier_1_Gold",
                        "vision_predicted_label": sp_name,
                        "c_error": 0.08,
                        "is_label_corrupted": False,
                    })

            pd.DataFrame(morph_rows).to_csv(morph_csv, index=False)
            pd.DataFrame(vis_rows).to_csv(vis_csv, index=False)

            res = subprocess.run(
                [
                    "Rscript", str(DISCOVERY_R_SCRIPT),
                    "--morphometrics", str(morph_csv),
                    "--vision-audit", str(vis_csv),
                    "--focal-taxon", "Packera dubia",
                    "--output-complexes", str(out_complexes),
                    "--output-affinity", str(out_affinity),
                    "--output-neighbors", str(out_neighbors),
                    "--output-figure", str(out_fig),
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
        """Verify pub_figures_utils.R and 06_multimodal_spatial_rf.R use morphotype clusters when <2 species."""
        self.assertIn("STEP08_INTRA_FIGS_OK", self._require_batched_r())

    def test_multimodal_cli_flag_and_ecological_edaphic_left_join(self) -> None:
        """Verify 08_discover_taxon_complexes.R (<=800 lines) parses --multimodal and merges Step 06 ecological/edaphic columns onto morph_df."""
        line_count = len(DISCOVERY_R_SCRIPT.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(
            line_count,
            800,
            f"08_discover_taxon_complexes.R has {line_count} lines; must be <= 800 lines.",
        )
        self.assertIn("STEP08_MULTIMODAL_JOIN_OK", self._require_batched_r())


if __name__ == "__main__":
    unittest.main()
