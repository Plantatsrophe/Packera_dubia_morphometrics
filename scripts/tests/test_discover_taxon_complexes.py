#!/usr/bin/env python3
"""
===============================================================================
Script: test_discover_taxon_complexes.py
Project: Packera Genus-Wide Taxonomic Delimitation & Morphometrics Pipeline

Description:
    Unit & integration tests verifying genus-wide dynamic North American
    Packera generalization:
    1. Dynamic standardization across North American Packera species and
       historical Senecio basionyms (including western, boreal, and Mexican
       taxa) in Python and R with zero hardcoded species lists.
    2. Bipinnatisect/tripinnatisect dissection gate (P. millefolium, solidity ~0.35)
       passing Tier 1 without false rejection.
    3. Step 08 (08_discover_taxon_complexes.R) unsupervised discovery of 3
       distinct morphological-ecological species complexes across 8 species,
       within-species GMM morphotype subcluster resolution, and --multimodal join.
===============================================================================
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon
from scripts.vision.lm2_geometry_utils import classify_and_route_leaf


class TestGenusWidePackeraAndComplexDiscovery(unittest.TestCase):
    """Tests for dynamic North American Packera support & Step 08 complex discovery."""

    @classmethod
    def setUpClass(cls) -> None:
        cls._class_temp_dir = tempfile.TemporaryDirectory()
        cls._class_tmp_path = Path(cls._class_temp_dir.name)
        cls.intra_fig_dir = cls._class_tmp_path / "intra_figs"
        cls.intra_fig_dir.mkdir(parents=True, exist_ok=True)
        cls.mm_join_dir = cls._class_tmp_path / "mm_join"
        cls.mm_join_dir.mkdir(parents=True, exist_ok=True)
        cls.r_unit_proc = None

        if shutil.which("Rscript"):
            fig_3d = cls._class_tmp_path / "test_3d_hollow.pdf"
            fig_intra = cls.intra_fig_dir / "cda_passive_projection.pdf"
            morph_csv = cls.mm_join_dir / "morphometric_flags.csv"
            mm_csv = cls.mm_join_dir / "multimodal_conflict_flags.csv"
            out_csv = cls.mm_join_dir / "discovered_taxon_complexes.csv"
            out_fig = cls.mm_join_dir / "fig08_discovered_taxon_complexes.pdf"

            rng = np.random.default_rng(808)
            m_rows = []
            mm_rows = []
            for sp, pc1_mu, ph_mu, doy_mu in [
                ("Packera aurea", -2.5, 5.4, 118.0),
                ("Packera obovata", 2.5, 7.2, 108.0),
                ("Packera schweinitziana", 0.2, 4.8, 165.0),
            ]:
                for j in range(6):
                    cat = f"{sp.split()[1][:3].upper()}_{j:02d}"
                    m_rows.append({
                        "catalogNumber": cat,
                        "plant_individual_id": 1,
                        "species_standardized": sp,
                        "determiner_tier": "Tier_1_Gold",
                        "PC1": float(pc1_mu + rng.normal(0.0, 0.15)),
                        "PC2": float(rng.normal(0.0, 0.15)),
                        "PC3": float(rng.normal(0.0, 0.10)),
                    })
                    mm_rows.append({
                        "catalogNumber": cat,
                        "plant_individual_id": 1,
                        "species_standardized": sp,
                        "fine_soil_ph": float(ph_mu + rng.normal(0.0, 0.1)),
                        "fine_soil_sand": 45.0,
                        "fine_bedrock_depth_cm": 80.0,
                        "ssurgo_rock_outcrop_pct": 10.0,
                        "bio1_temp_mean": 13.5,
                        "bio12_precip_annual": 1150.0,
                        "doy": float(doy_mu + rng.normal(0.0, 2.0)),
                        "decimalLatitude": 35.8,
                        "decimalLongitude": -81.5,
                    })
            pd.DataFrame(m_rows).to_csv(morph_csv, index=False)
            pd.DataFrame(mm_rows).to_csv(mm_csv, index=False)

            r_batch_code = f"""
            source("scripts/morphometrics/pub_figures_utils.R")

            # Block 1: R standardize_taxon_name & dynamic Okabe-Ito/Polychrome palette (before sourcing downstream scripts)
            local({{
              stopifnot(!exists("EASTERN_NA_PACKERA_TAXA"))
              stopifnot(!exists("SENECIO_BASIONYM_TO_PACKERA"))
              stopifnot(!exists("CANONICAL_PACKERA_PALETTE"))
              stopifnot(standardize_taxon_name("Packera schweinitziana (Nutt.) W.A.Weber & A.Love") == "Packera schweinitziana")
              stopifnot(standardize_taxon_name("Packera streptanthifolia (Greene) W.A.Weber & Á.Löve") == "Packera streptanthifolia")
              stopifnot(standardize_taxon_name("Packera cymbalaria (Pursh) W.A.Weber & Á.Löve") == "Packera cymbalaria")
              stopifnot(standardize_taxon_name("P. sanguisorbae (DC.) C.Jeffrey") == "Packera sanguisorbae")
              stopifnot(standardize_taxon_name("Packera x memmingeri") == "Packera memmingeri")
              stopifnot(is.na(standardize_taxon_name("Arabidopsis thaliana")))
              pal <- get_dynamic_packera_palette(c(
                "Packera aurea", "Packera cana", "Packera cymbalaria",
                "Packera millefolium", "Packera sanguisorbae", "Packera streptanthifolia",
                "Packera schweinitziana", "Packera serpenticola", "Packeraベル",
                "Packera hartiana", "Packera pseudaurea", "Packera subnuda"
              ))
              stopifnot(length(pal) >= 12)
              stopifnot(all(grepl("^#[0-9A-Fa-f]{{6}}$", pal)))
              cat("R_PALETTE_STD_OK\\n")
            }})

            # Block 2: Hierarchical 3-level group precedence, missing-modality distance renormalization, & 3D hollow points
            local({{
              source("scripts/analysis/06_multimodal_spatial_rf.R")
              source("scripts/analysis/08_discover_taxon_complexes.R")

              df_multi <- data.frame(
                species_standardized = c("Packera aurea", "Packera aurea", "Packera obovata", "Packera obovata"),
                gmm_cluster = c(1, 2, 1, 2),
                regional_group = c("A", "B", "A", "B"),
                stringsAsFactors = FALSE
              )
              res_multi <- resolve_analysis_group(df_multi, intraspecific_mode = FALSE)
              stopifnot(identical(res_multi$analysis_group, df_multi$species_standardized))

              df_intra_gmm <- data.frame(
                species_standardized = rep("Packera dubia", 6),
                gmm_cluster = c(1, 1, 1, 2, 2, 2),
                regional_group = c("Coastal", "Piedmont", "Coastal", "Piedmont", "Coastal", "Piedmont"),
                stringsAsFactors = FALSE
              )
              res_intra_gmm <- resolve_analysis_group(df_intra_gmm, intraspecific_mode = TRUE)
              stopifnot(all(grepl("^Morphotype_K[12]$", res_intra_gmm$analysis_group)))

              df_intra_k1 <- data.frame(
                species_standardized = rep("Packera dubia", 6),
                gmm_cluster = rep(1, 6),
                fine_soil_ph = c(4.6, 4.7, 4.8, 6.7, 6.8, 6.9),
                fine_soil_sand = c(82, 80, 79, 25, 24, 22),
                fine_bedrock_depth_cm = c(180, 175, 170, 25, 22, 20),
                stringsAsFactors = FALSE
              )
              res_intra_k1 <- resolve_analysis_group(df_intra_k1, intraspecific_mode = TRUE)
              stopifnot(length(unique(res_intra_k1$analysis_group)) == 2L)
              stopifnot(all(grepl("^Ecotype_[12]_", res_intra_k1$analysis_group)))

              m_dist <- matrix(c(0, 2, 4, 2, 0, 2, 4, 2, 0), nrow = 3)
              e_dist <- matrix(c(0, NA, 6, NA, 0, NA, 6, NA, 0), nrow = 3)
              p_dist <- matrix(c(0, 1, NA, 1, 0, 3, NA, 3, 0), nrow = 3)
              g_dist <- matrix(0, nrow = 3, ncol = 3)
              comb <- combine_multimodal_distances(
                list(morph = m_dist, eco = e_dist, pheno = p_dist, geo = g_dist),
                c(morph = 0.50, eco = 0.30, pheno = 0.20, geo = 0.00)
              )
              stopifnot(all(is.finite(comb$D_joint)))
              stopifnot(abs(comb$D_joint[1, 1]) < 1e-9)
              expected_12 <- (0.50 * (2 / 4) + 0.20 * (1 / 3)) / (0.50 + 0.20)
              stopifnot(abs(comb$D_joint[1, 2] - expected_12) < 1e-6)

              pdf("{fig_3d.as_posix()}", width = 6, height = 5)
              pts_df <- data.frame(
                x = rnorm(12), y = rnorm(12), z = rnorm(12),
                grp = rep(c("Packera aurea", "Packera glabella", "Packera obovata"), each = 4)
              )
              pal3 <- get_dynamic_packera_palette(unique(pts_df$grp))
              draw_3d_morphospace_panel(
                pts_df, "x", "y", "z", "grp", pal3,
                "3D Test", "X", "Y", "Z",
                point_pch = 1, point_cex = 0.85, point_lwd = 1.2
              )
              dev.off()
              cat("PRECEDENCE_AND_3D_OK\\n")
            }})

            # Block 3: Intraspecific Morphotype_K1..Kn cluster visualization across Step 04 R figures
            local({{
              set.seed(42)
              n <- 24
              audit_df <- data.frame(
                catalogNumber = paste0("NCU", sprintf("%04d", seq_len(n))),
                species_standardized = rep("Packera dubia", n),
                determiner_tier = rep(c("Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze"), length.out = n),
                is_passive = rep(c(FALSE, FALSE, TRUE), length.out = n),
                gmm_cluster = rep(1:3, each = 8),
                gmm_uncertainty = runif(n, 0.01, 0.25),
                cda_predicted_taxon = paste0("Morphotype_K", rep(1:3, each = 8)),
                can1 = c(rnorm(8, -2.5, 0.4), rnorm(8, 0.2, 0.4), rnorm(8, 2.6, 0.4)),
                can2 = c(rnorm(8, 1.2, 0.4), rnorm(8, -1.8, 0.4), rnorm(8, 0.9, 0.4)),
                can3 = c(rnorm(8, 0.5, 0.4), rnorm(8, 0.8, 0.4), rnorm(8, -1.3, 0.4)),
                PC1 = rnorm(n),
                PC2 = rnorm(n),
                PC3 = rnorm(n),
                log_centsize = rnorm(n, 4.0, 0.2),
                stringsAsFactors = FALSE
              )
              gmm_res <- list(
                best_k = 3L,
                model_name = "VVV",
                classification = audit_df$gmm_cluster,
                uncertainty = audit_df$gmm_uncertainty,
                bic_table = data.frame(
                  K = 1:4,
                  BIC = c(-210.0, -185.0, -162.0, -170.0),
                  Two_Delta_BIC = c(-48.0, -23.0, 0.0, -8.0)
                )
              )
              cda_res <- list(
                var_explained = c(64.2, 35.8),
                loocv_accuracy = 0.917,
                loadings = matrix(
                  c(0.8, -0.3, 0.2, 0.4, 0.7, -0.2, -0.2, 0.4, 0.8),
                  nrow = 3,
                  dimnames = list(c("PC1", "PC2", "PC3"), c("Can1", "Can2", "Can3"))
                )
              )
              resolved <- resolve_plot_grouping(audit_df, prefer_clusters_single_taxon = TRUE)
              stopifnot(isTRUE(resolved$is_cluster_mode))
              stopifnot(identical(sort(unique(resolved$group_vec)), c("Morphotype_K1", "Morphotype_K2", "Morphotype_K3")))
              stopifnot(length(unique(resolved$palette[c("Morphotype_K1", "Morphotype_K2", "Morphotype_K3")])) == 3L)

              render_cda_biplot_4panel_pdf(audit_df, cda_res, gmm_res, "{fig_intra.as_posix()}", c("PC1", "PC2", "PC3"))
              generate_step04_pub_figures(audit_df, cda_res, gmm_res, NULL, "{fig_intra.as_posix()}")
              cat("INTRA_CLUSTER_FIGURES_OK\\n")
            }})

            # Block 4: Step 08 --multimodal CLI flag parsing & ecological/phenological left-join
            local({{
              assign(
                "commandArgs",
                function(trailingOnly = TRUE) c(
                  "--morphometrics", "{morph_csv.as_posix()}",
                  "--multimodal", "{mm_csv.as_posix()}",
                  "--output-csv", "{out_csv.as_posix()}",
                  "--output-fig", "{out_fig.as_posix()}"
                ),
                envir = .GlobalEnv
              )
              on.exit(rm("commandArgs", envir = .GlobalEnv), add = TRUE)
              opts <- parse_discovery_args()
              stopifnot(opts$multimodal == "{mm_csv.as_posix()}")
              res_df <- run_taxon_complex_discovery(opts)
              stopifnot(nrow(res_df) == 3L)
              stopifnot(all(is.finite(res_df$mean_soil_ph)))
              stopifnot(all(is.finite(res_df$mean_doy)))
              cat("STEP08_MULTIMODAL_JOIN_OK\\n")
            }})
            """
            cls.r_unit_proc = subprocess.run(
                ["Rscript", "-e", r_batch_code],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                timeout=90,
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
            f"Batched Step 08 R unit suite failed:\n{self.r_unit_proc.stderr}\n{self.r_unit_proc.stdout}",
        )
        return self.r_unit_proc.stdout

    def test_standardize_packera_taxon_all_eastern_species_and_basionyms(self) -> None:
        """Verify dynamic taxon standardization across North American Packera species & Senecio basionyms."""
        direct_cases = {
            "Packera anonyma (Alph.Wood) W.A.Weber & Á.Löve": "Packera anonyma",
            "Packera aurea (L.) Á.Löve & D.Löve": "Packera aurea",
            "Packera crawfordii (Britton) A.M.Mahoney & R.R.Kowal": "Packera crawfordii",
            "Packera debilis (Nutt.) W.A.Weber & Á.Löve": "Packera debilis",
            "Packera dubia (Spreng.) Trock & Mabb.": "Packera dubia",
            "Packera glabella (Poir.) C.Jeffrey": "Packera glabella",
            "Packera hartiana (A.Heller) W.A.Weber & Á.Löve": "Packera hartiana",
            "Packera indecora (Greene) Á.Löve & D.Löve": "Packera indecora",
            "Packera insulae-regalis R.R.Kowal": "Packera insulae-regalis",
            "Packera malmstenii": "Packera malmstenii",
            "Packera millefolium (Torr. & A.Gray) W.A.Weber & Á.Löve": "Packera millefolium",
            "Packera musiniensis (S.L.Welsh) Trock": "Packera musiniensis",
            "Packera obovata (Muhl. ex Willd.) W.A.Weber & Á.Löve": "Packera obovata",
            "Packera paupercula (Michx.) Á.Löve & D.Löve": "Packera paupercula",
            "Packera paupercula var. savannarum R.R.Kowal": "Packera paupercula",
            "Packera plattensis (Nutt.) W.A.Weber & Á.Löve": "Packera plattensis",
            "Packera pseudaurea (Rydb.) W.A.Weber & Á.Löve": "Packera pseudaurea",
            "Packera quaerens (Greene) W.A.Weber & Á.Löve": "Packera quaerens",
            "Packera schweinitziana (Nutt.) W.A.Weber & Á.Löve": "Packera schweinitziana",
            "Packera serpenticola Boufford, Kartesz, S.H.Shi & R.Zhou": "Packera serpenticola",
            "Packera spellenbergii (T.M.Barkley) C.Jeffrey": "Packera spellenbergii",
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve": "Packera streptanthifolia",
            "Packera cymbalaria (Pursh) W.A.Weber & Á.Löve": "Packera cymbalaria",
            "Packera cana (Hook.) W.A.Weber & Á.Löve": "Packera cana",
            "Packera fendleri (A.Gray) W.A.Weber & Á.Löve": "Packera fendleri",
            "Packera multilobata (Torr. & A.Gray) W.A.Weber & Á.Löve": "Packera multilobata",
            "Packera werneriifolia (A.Gray) W.A.Weber & Á.Löve": "Packera werneriifolia",
            "Packera crocata (Rydb.) W.A.Weber & Á.Löve": "Packera crocata",
            "Packera dimorphophylla (Greene) W.A.Weber & Á.Löve": "Packera dimorphophylla",
            "Packera neomexicana (A.Gray) W.A.Weber & Á.Löve": "Packera neomexicana",
            "Packera tampicana (DC.) C.Jeffrey": "Packera tampicana",
            "Packera coahuilensis (Greenm.) C.Jeffrey": "Packera coahuilensis",
            "Packera bellidifolia (Kunth) W.A.Weber & Á.Löve": "Packera bellidifolia",
            "Packera zimapanica (Hemsl.) C.C.Freeman & T.M.Barkley": "Packera zimapanica",
            "Packera hyperborealis (Greenm.) Á.Löve & D.Löve": "Packera hyperborealis",
            "Packera ogotorukensis (Packer) Á.Löve & D.Löve": "Packera ogotorukensis",
            "P. sanguisorbae (DC.) C.Jeffrey": "Packera sanguisorbae",
            "Packera subnuda (DC.) Trock & T.M.Barkley": "Packera subnuda",
            "Packera × memmingeri (Britton ex Small) Weakley": "Packera memmingeri",
            "Packera x memmingeri": "Packera memmingeri",
            "Arabidopsis thaliana": "Unknown",
            "Packera sp.": "Unknown",
            "": "Unknown",
        }
        for raw_name, expected_canon in direct_cases.items():
            self.assertEqual(
                standardize_packera_taxon(raw_name),
                expected_canon,
                f"Failed to standardize {raw_name!r} -> {expected_canon!r}",
            )

        gbif_synonym_mocks = {
            "Senecio tomentosus Michx.": "Packera dubia",
            "Packera tomentosa (Michx.) C.Jeffrey": "Packera dubia",
            "Senecio anonymus Alph.Wood": "Packera anonyma",
            "Senecio smallii Britton": "Packera anonyma",
            "Senecio aureus L.": "Packera aurea",
            "Senecio obovatus Muhl. ex Willd.": "Packera obovata",
            "Senecio robbinsii Oakes ex Rusby": "Packera schweinitziana",
            "Senecio schweinitzianus Nutt.": "Packera schweinitziana",
            "Senecio millefolium Torr. & A.Gray": "Packera millefolium",
            "Senecio glabellus Poir.": "Packera glabella",
            "Packera eburnea": "Packera serpenticola",
            "Senecio streptanthifolius Greene": "Packera streptanthifolia",
            "Senecio canus Hook.": "Packera cana",
            "Senecio sanguisorbae DC.": "Packera sanguisorbae",
        }

        def _fake_backbone(name: str, **kwargs):
            for k, target_sp in gbif_synonym_mocks.items():
                if name.lower().startswith(k.split()[0].lower() + " " + k.split()[1].lower()):
                    return {
                        "matchType": "EXACT",
                        "status": "SYNONYM",
                        "species": target_sp,
                        "canonicalName": " ".join(k.split()[:2]),
                    }
            return {"matchType": "NONE"}

        with patch("scripts.core.harvester.species.name_backbone", side_effect=_fake_backbone), patch(
            "scripts.core.harvester.species.name_usage", return_value={"results": []}
        ):
            for raw_syn, expected_canon in gbif_synonym_mocks.items():
                self.assertEqual(
                    standardize_packera_taxon(raw_syn),
                    expected_canon,
                    f"Dynamic GBIF synonym resolution failed for {raw_syn!r} -> {expected_canon!r}",
                )

    def test_r_standardize_taxon_and_dynamic_palette(self) -> None:
        """Verify R standardize_taxon_name() and dynamic Okabe-Ito/Polychrome palette."""
        self.assertIn("R_PALETTE_STD_OK", self._require_batched_r())

    def test_bipinnatisect_dissection_gate_tier1_routing(self) -> None:
        """Verify finely bipinnatisect leaf (P. millefolium, solidity ~0.35) passes Tier 1."""
        mask = np.zeros((500, 360), dtype=np.uint8)
        cx = 180
        cv2.line(mask, (cx, 50), (cx, 450), 255, thickness=8)
        for y in range(85, 420, 32):
            span = int(95 * np.sin(np.pi * (y - 45) / 410.0) + 25)
            cv2.line(mask, (cx - span, y - 8), (cx + span, y - 8), 255, thickness=6)
            for dx in range(22, span, 22):
                cv2.line(mask, (cx - dx, y - 22), (cx - dx, y + 6), 255, thickness=4)
                cv2.line(mask, (cx + dx, y - 22), (cx + dx, y + 6), 255, thickness=4)

        res = classify_and_route_leaf(mask)
        self.assertLess(res.solidity, 0.48, f"Expected bipinnatisect solidity < 0.48, got {res.solidity:.3f}")
        self.assertGreaterEqual(res.solidity, 0.25)
        self.assertTrue(res.is_dissected, "Bipinnatisect silhouette must trigger botanical dissection gate.")
        self.assertEqual(res.assigned_tier, "tier1", f"Expected tier1 routing, got {res.assigned_tier} ({res.rejection_reason})")

    def test_discover_taxon_complexes_8_species_3_complexes(self) -> None:
        """Verify 08_discover_taxon_complexes.R discovers 3 clean complexes from 8 synthetic species."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not available")
        rng = np.random.default_rng(2026)

        species_specs = {
            "Packera aurea": ("C1", [-3.5, 2.2, 0.4], 6.4, 30.0, 140.0, 2.0, 10.5, 1150.0, 132, 39.5, -79.5),
            "Packera pseudaurea": ("C1", [-3.2, 2.5, 0.2], 6.6, 33.0, 135.0, 3.0, 9.8, 1050.0, 136, 41.0, -88.0),
            "Packera schweinitziana": ("C1", [-3.7, 2.0, 0.5], 6.1, 35.0, 125.0, 4.0, 7.5, 1250.0, 145, 44.0, -72.5),
            "Packera dubia": ("C2", [3.6, 1.8, -0.5], 4.8, 72.0, 35.0, 38.0, 16.5, 1280.0, 108, 34.8, -79.2),
            "Packera anonyma": ("C2", [3.3, 1.5, -0.3], 5.1, 64.0, 42.0, 32.0, 15.8, 1240.0, 126, 35.6, -80.1),
            "Packera serpenticola": ("C2", [3.8, 2.0, -0.6], 6.7, 58.0, 22.0, 55.0, 14.9, 1380.0, 122, 35.1, -83.6),
            "Packera millefolium": ("C3", [0.1, -4.2, 2.8], 5.3, 68.0, 18.0, 68.0, 13.8, 1520.0, 140, 35.4, -82.8),
            "Packera glabella": ("C3", [0.4, -3.9, 2.5], 6.8, 42.0, 180.0, 0.0, 16.2, 1360.0, 105, 33.5, -88.5),
        }

        rows = []
        idx = 1
        for sp, (_, pcs, ph, sand, bdr, outc, bio1, bio12, doy, lat, lon) in species_specs.items():
            for _ in range(12):
                rows.append({
                    "catalogNumber": f"GENUS{idx:04d}",
                    "species_standardized": sp,
                    "determiner_tier": "Tier_1_Gold",
                    "PC1": pcs[0] + rng.normal(0, 0.22),
                    "PC2": pcs[1] + rng.normal(0, 0.22),
                    "PC3": pcs[2] + rng.normal(0, 0.18),
                    "PC4": rng.normal(0, 0.10),
                    "PC5": rng.normal(0, 0.10),
                    "aspect_ratio": 2.0 + 0.2 * pcs[0] + rng.normal(0, 0.05),
                    "solidity": 0.42 if pcs[1] < -2 else 0.86 + rng.normal(0, 0.02),
                    "fine_soil_ph": ph + rng.normal(0, 0.15),
                    "fine_soil_sand": sand + rng.normal(0, 3.0),
                    "fine_bedrock_depth_cm": bdr + rng.normal(0, 5.0),
                    "ssurgo_rock_outcrop_pct": max(0.0, outc + rng.normal(0, 3.0)),
                    "bio1_temp_mean": bio1 + rng.normal(0, 0.4),
                    "bio12_precip_annual": bio12 + rng.normal(0, 30.0),
                    "doy": doy + rng.normal(0, 4.0),
                    "decimalLatitude": lat + rng.normal(0, 0.4),
                    "decimalLongitude": lon + rng.normal(0, 0.4),
                })
                idx += 1

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            in_csv = tmp_path / "synthetic_genus_triage.csv"
            out_csv = tmp_path / "discovered_taxon_complexes.csv"
            out_pdf = tmp_path / "fig08_discovered_taxon_complexes.pdf"
            pd.DataFrame(rows).to_csv(in_csv, index=False)

            proc = subprocess.run(
                [
                    "Rscript",
                    "scripts/analysis/08_discover_taxon_complexes.R",
                    "--triage", str(in_csv),
                    "--output-csv", str(out_csv),
                    "--output-fig", str(out_pdf),
                    "--knn", "2",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, f"Step 08 failed:\n{proc.stdout}\n{proc.stderr}")
            self.assertTrue(out_csv.exists(), "discovered_taxon_complexes.csv was not created.")
            self.assertTrue(out_pdf.exists(), "fig08_discovered_taxon_complexes.pdf was not created.")
            self.assertTrue(out_pdf.with_suffix(".png").exists(), "Companion PNG was not created.")

            res_df = pd.read_csv(out_csv)
            self.assertEqual(len(res_df), 8)
            self.assertEqual(res_df["complex_id"].nunique(), 3, f"Expected 3 complexes, got:\n{res_df[['species', 'complex_id']]}")

            c_map = dict(zip(res_df["species"], res_df["complex_id"]))
            self.assertEqual(c_map["Packera aurea"], c_map["Packera pseudaurea"])
            self.assertEqual(c_map["Packera aurea"], c_map["Packera schweinitziana"])
            self.assertEqual(c_map["Packera dubia"], c_map["Packera anonyma"])
            self.assertEqual(c_map["Packera dubia"], c_map["Packera serpenticola"])
            self.assertEqual(c_map["Packera millefolium"], c_map["Packera glabella"])
            self.assertNotEqual(c_map["Packera aurea"], c_map["Packera dubia"])
            self.assertNotEqual(c_map["Packera dubia"], c_map["Packera millefolium"])

    def test_hierarchical_precedence_missing_modality_and_3d_hollow_ordinations(self) -> None:
        """Verify 3-level group precedence, missing-modality distance renormalization, and 3D hollow points."""
        out = self._require_batched_r()
        self.assertIn("PRECEDENCE_AND_3D_OK", out)
        fig_3d = self._class_tmp_path / "test_3d_hollow.pdf"
        self.assertTrue(fig_3d.exists() and fig_3d.stat().st_size > 1000)

    def test_intraspecific_morphotype_cluster_visualization_across_r_figures(self) -> None:
        """Verify single-taxon intraspecific Step 04 figures color/group by Morphotype_K1..Kn."""
        out = self._require_batched_r()
        self.assertIn("INTRA_CLUSTER_FIGURES_OK", out)
        for fname in (
            "cda_passive_projection.pdf",
            "fig04a_gmm_bic_and_uncertainty.pdf",
            "fig04b_cda_biplot_and_allometry.pdf",
            "fig04_gmm_cda_clustering.pdf",
        ):
            fpath = self.intra_fig_dir / fname
            self.assertTrue(fpath.exists() and fpath.stat().st_size > 1500, f"Missing or empty {fname}")

    def test_subcluster_network_when_2_species_or_within_species_gmm_clusters(self) -> None:
        """Verify Step 08 resolves Taxon x GMM Morphotype subclusters when <= 2 species are analyzed."""
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not available")
        rng = np.random.default_rng(99)
        rows = []
        sub_specs = [
            ("Packera dubia", 1, [-3.0, 1.5, 0.2], 4.6, 78.0, 160.0, 104, 33.8, -78.5),
            ("Packera dubia", 2, [-1.8, 0.9, -0.2], 5.3, 60.0, 28.0, 114, 35.3, -80.2),
            ("Packera anonyma", 1, [2.6, -1.4, 0.4], 5.5, 52.0, 35.0, 128, 35.9, -79.9),
            ("Packera anonyma", 2, [3.4, -2.1, -0.3], 6.4, 38.0, 95.0, 139, 36.4, -81.5),
        ]
        idx = 1
        for sp, k_cl, pcs, ph, sand, bdr, doy, lat, lon in sub_specs:
            for _ in range(10):
                rows.append({
                    "catalogNumber": f"SUB{idx:04d}",
                    "species_standardized": sp,
                    "determiner_tier": "Tier_1_Gold",
                    "gmm_cluster": k_cl,
                    "PC1": pcs[0] + rng.normal(0, 0.18),
                    "PC2": pcs[1] + rng.normal(0, 0.18),
                    "PC3": pcs[2] + rng.normal(0, 0.15),
                    "fine_soil_ph": ph + rng.normal(0, 0.12),
                    "fine_soil_sand": sand + rng.normal(0, 2.5),
                    "fine_bedrock_depth_cm": bdr + rng.normal(0, 4.0),
                    "doy": doy + rng.normal(0, 3.0),
                    "decimalLatitude": lat + rng.normal(0, 0.25),
                    "decimalLongitude": lon + rng.normal(0, 0.25),
                })
                idx += 1

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            in_csv = tmp_path / "two_species_subclusters.csv"
            out_csv = tmp_path / "discovered_taxon_complexes.csv"
            out_pdf = tmp_path / "fig08_discovered_taxon_complexes.pdf"
            pd.DataFrame(rows).to_csv(in_csv, index=False)

            proc = subprocess.run(
                [
                    "Rscript",
                    "scripts/analysis/08_discover_taxon_complexes.R",
                    "--triage", str(in_csv),
                    "--output-csv", str(out_csv),
                    "--output-fig", str(out_pdf),
                    "--knn", "1",
                ],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, f"Step 08 2-species subcluster run failed:\n{proc.stdout}\n{proc.stderr}")
            res_df = pd.read_csv(out_csv)
            self.assertEqual(len(res_df), 4, f"Expected 4 subcluster nodes, got:\n{res_df}")
            self.assertIn("parent_taxon", res_df.columns)
            self.assertTrue(all("Morphotype_K" in s for s in res_df["species"]))
            self.assertGreaterEqual(res_df["complex_id"].nunique(), 2)
            self.assertTrue(out_pdf.exists() and out_pdf.stat().st_size > 2000)

    def test_multimodal_cli_flag_and_ecological_edaphic_left_join(self) -> None:
        """Verify 08_discover_taxon_complexes.R (<=650 lines) parses --multimodal and left-joins Step 06 edaphic/phenology columns onto Step 04 morphometrics."""
        r_script = PROJECT_ROOT / "scripts" / "analysis" / "08_discover_taxon_complexes.R"
        line_count = len(r_script.read_text(encoding="utf-8").splitlines())
        self.assertLessEqual(line_count, 650, f"08_discover_taxon_complexes.R has {line_count} lines; must remain <= 650.")

        out = self._require_batched_r()
        self.assertIn("STEP08_MULTIMODAL_JOIN_OK", out)
        out_csv = self.mm_join_dir / "discovered_taxon_complexes.csv"
        out_fig = self.mm_join_dir / "fig08_discovered_taxon_complexes.pdf"
        self.assertTrue(out_csv.exists())
        self.assertTrue(out_fig.exists())


if __name__ == "__main__":
    unittest.main()
