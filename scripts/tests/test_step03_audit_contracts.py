#!/usr/bin/env python3
"""Batched warm-R contract suite for the Step 03 statistical and architectural refactor."""

from __future__ import annotations

import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db_schema import init_db_schema
from scripts.tests._r_session import rscript_available


def _contour_xy(rx: float, ry: float, asym: float = 0.0, tip_up: bool = True) -> np.ndarray:
    t = np.linspace(0.0, 2.0 * np.pi, 120, endpoint=False)
    if tip_up:
        x = 250.0 - rx * np.sin(t) * (1.0 + asym * np.cos(t))
        y = 250.0 + ry * np.cos(t)
    else:
        x = 250.0 - rx * np.cos(t)
        y = 250.0 + ry * np.sin(t) * (1.0 + asym * np.cos(t))
    return np.column_stack([x, y])


def _seed_sqlite_cohort(tmp_p: Path) -> Path:
    db_path = tmp_p / "step03_sql.sqlite"
    sql_dir = tmp_p / "sql_run"
    (sql_dir / "contours").mkdir(parents=True, exist_ok=True)
    (sql_dir / "tables").mkdir(parents=True, exist_ok=True)
    (sql_dir / "figures").mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    init_db_schema(conn)
    v_rows = []
    p_rows = []
    l_rows = []
    specs = [
        ("V_CAL_1", "Packera aurea", 65.0, 35.0, 0.0),
        ("V_CAL_2", "Packera aurea", 62.0, 38.0, 0.04),
        ("V_CAL_3", "Packera aurea", 60.0, 36.0, -0.03),
        ("V_CAL_4", "Packera dubia", 42.0, 72.0, 0.0),
        ("V_CAL_5", "Packera dubia", 44.0, 70.0, 0.03),
        ("V_TICK_1", "Packera dubia", 45.0, 68.0, 0.0),
    ]
    for cat, sp, rx, ry, asym in specs:
        is_tick = cat == "V_TICK_1"
        ppm = 11.4 if is_tick else 10.0
        cal = 0 if is_tick else 1
        ss = "tick_lattice_unboxed" if is_tick else "ruler_calibrated"
        v_rows.append((cat, sp, sp, "Tier_1_Gold", "", 1, ppm, cal))
        p_rows.append((cat, 0))
        leaves = [(1, "basal", 1, True, 0.50, -18.5, 3.2, 14.4), (2, "basal", 0, False, 0.60, -19.5, 3.6, 15.0), (3, "cauline", None, False, 0.15, -2.0, 0.1, 1.0)] if is_tick else [(1, "basal", 1, True, 0.45, -10.0, 1.0, 4.0)]
        for lid, lz, up, tu, pbr, ld, ad, bd in leaves:
            pts = _contour_xy(rx if tu else ry, ry if tu else rx, asym=asym, tip_up=tu)
            l_rows.append((
                cat, 0, lid, f"{cat}_p0_leaf{lid}", 0.90, 0.85, 2.0, float(np.pi * rx * ry),
                210.0 if is_tick else 180.0, 12.5, 1.2, 25.0, 14.0, pbr, ppm, cal, ss,
                up, json.dumps([[round(float(px), 2), round(float(py), 2)] for px, py in pts]),
                lz, 0.92, 0.10 if lz == "basal" else 0.45, 0.80 if lz == "basal" else 0.25,
                0.95, 0.90, 0.05, 1, 1, 1, ld, ad, bd,
            ))
    conn.executemany("INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier, type_status, is_primary_duplicate, pixels_per_mm, is_calibrated_mm) VALUES (?, ?, ?, ?, ?, ?, ?, ?);", v_rows)
    conn.executemany("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, ?);", p_rows)
    conn.executemany(
        """INSERT INTO extracted_leaves (
            catalogNumber, plant_individual_id, leaf_id, shape_id, solidity, circularity, aspect_ratio, area_px,
            area_mm2, petiole_length_mm, petiole_width_mm, blade_length_mm, blade_width_mm, petiole_to_blade_ratio,
            pixels_per_mm, is_calibrated_mm, scale_source, mask_is_upright, contour_xy_json, leaf_zone,
            leaf_zone_confidence, dist_from_crown_norm, cy_rel_plant_bbox, archetype_score, bilateral_dice, si_a,
            lm3_gates_pass, n_lamina_components, efa_qc_pass, blade_cie_l_sheet_delta, blade_cie_a_delta, blade_cie_b_delta
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);""",
        l_rows,
    )
    conn.commit()
    conn.close()
    pd.DataFrame(
        [{"catalogNumber": r[0], "scientificName": r[1], "species_standardized": r[2], "determiner_tier": r[3], "type_status": r[4], "is_primary_duplicate": 1, "pixels_per_mm": r[6], "is_calibrated_mm": r[7]} for r in v_rows]
    ).to_csv(sql_dir / "tables" / "curated_vouchers.csv", index=False)
    return db_path


def _seed_csv_pca_cohort(base_dir: Path, rep_p1: int) -> None:
    cdir = base_dir / "contours"
    cdir.mkdir(parents=True, exist_ok=True)
    v_rows = []
    m_rows = []
    items = [("P1", "Packera aurea", 70.0, 35.0), ("P2", "Packera aurea", 65.0, 38.0), ("P3", "Packera aurea", 60.0, 40.0), ("P4", "Packera dubia", 40.0, 74.0), ("P5", "Packera dubia", 42.0, 70.0), ("P6", "Packera dubia", 44.0, 66.0)]
    for cat, sp, rx, ry in items:
        v_rows.append({"catalogNumber": cat, "scientificName": sp, "species_standardized": sp, "determiner_tier": "Tier_1_Gold", "is_primary_duplicate": 1})
        n_l = rep_p1 if cat == "P1" else 1
        pts = _contour_xy(rx, ry, tip_up=True)
        for lid in range(1, n_l + 1):
            cp = cdir / f"{cat}_leaf{lid}.csv"
            pd.DataFrame({"x": pts[:, 0], "y": pts[:, 1]}).to_csv(cp, index=False)
            m_rows.append({"catalogNumber": cat, "plant_individual_id": 0, "leaf_id": lid, "leaf_zone": "basal", "leaf_zone_confidence": 0.95, "mask_is_upright": 1, "pixels_per_mm": 10.0, "is_calibrated_mm": 1, "scale_source": "ruler_calibrated", "solidity": 0.90, "archetype_score": 0.95, "contour_path": str(cp)})
    pd.DataFrame(v_rows).to_csv(base_dir / "vouchers.csv", index=False)
    pd.DataFrame(m_rows).to_csv(base_dir / "manifest.csv", index=False)


class TestStep03AuditContracts(unittest.TestCase):
    """End-to-end R/SQLite/CSV contracts for refactored Step 03 morphometrics."""

    @classmethod
    def setUpClass(cls) -> None:
        if not rscript_available():
            raise unittest.SkipTest("Rscript not installed.")
        cls._tmp_dir = tempfile.TemporaryDirectory()
        cls.tmp_p = Path(cls._tmp_dir.name)
        cls.db_path = _seed_sqlite_cohort(cls.tmp_p)
        _seed_csv_pca_cohort(cls.tmp_p / "pca1", 1)
        _seed_csv_pca_cohort(cls.tmp_p / "pca6", 6)
        (cls.tmp_p / "custom_step03.yaml").write_text(
            "morphometrics:\n  harmonics: 8\n  num_pcs: 5\n  fourier:\n    min_harmonics: 6\n    max_harmonics: 12\n    min_export_pcs: 6\n    max_export_pcs: 10\n",
            encoding="utf-8",
        )
        cls.json_out = cls.tmp_p / "contracts.json"
        r_script = cls.tmp_p / "run_contracts.R"
        r_script.write_text(
            f"""
            if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
            source("scripts/morphometrics/step03_config_scale_utils.R")
            source("scripts/morphometrics/dissected_leaf_utils.R")
            source("scripts/morphometrics/fourier_parallel_utils.R")
            source("scripts/morphometrics/plasticity_mitigation_utils.R")
            source("scripts/morphometrics/authority_weighting_utils.R")
            source("scripts/core/db_utils.R")
            source("scripts/morphometrics/03_fourier_extractor.R")

            # 1. Zone routing & proximity sheet consolidation
            df_z <- data.frame(
              catalogNumber = rep("Z_NARROW", 3L), plant_individual_id = 0L, leaf_id = 1:3,
              leaf_zone = NA_character_, petiole_to_blade_ratio = c(0.03, 0.04, 0.90),
              area_px = c(1200, 1300, 4500), dist_from_crown_norm = c(0.11, 0.12, 0.09),
              cy_rel_plant_bbox = c(0.80, 0.78, 0.85), stringsAsFactors = FALSE
            )
            z_out <- classify_and_partition_leaf_zones(df_z)

            df_sh <- data.frame(
              catalogNumber = c(rep("SH_CAL", 4L), rep("SH_UNCAL", 3L)),
              plant_individual_id = c(0L, 0L, 1L, 2L, 0L, 0L, 1L),
              leaf_id = c(1L, 2L, 3L, 4L, 1L, 2L, 3L),
              leaf_zone = rep("basal", 7L),
              dist_from_crown_norm = c(0.08, 0.09, 0.31, 0.09, 0.08, 0.09, 0.10),
              cy_rel_plant_bbox = c(0.80, 0.82, 0.35, 0.80, 0.80, 0.82, 0.80),
              petiole_to_blade_ratio = c(0.4, 0.4, 0.5, 0.4, 0.4, 0.4, 0.4),
              pixels_per_mm = c(rep(10.0, 4L), rep(NA_real_, 3L)),
              is_calibrated_mm = c(rep(TRUE, 4L), rep(FALSE, 3L)),
              scale_source = c(rep("ruler_calibrated", 4L), rep("tick_lattice_unboxed", 3L)),
              cx_px = c(500, 520, 680, 1150, NA, NA, NA),
              cy_px = c(500, 510, 500, 500, NA, NA, NA),
              cx_rel_sheet = c(NA, NA, NA, NA, 0.50, 0.52, 0.82),
              cy_rel_sheet = c(NA, NA, NA, NA, 0.50, 0.51, 0.50),
              stringsAsFactors = FALSE
            )
            sh_out <- consolidate_fragmented_sheet_individuals(df_sh)

            # 2. SQLite fetch tip-up vectorization & Step 03 SQLite run
            Sys.setenv(PLANTMATICS_DB_PATH = "{cls.db_path}", RUN_ID = "run_s03_contracts")
            sql_dir <- "{cls.tmp_p / 'sql_run'}"
            v_csv <- file.path(sql_dir, "tables", "curated_vouchers.csv")
            db_ing <- fetch_step03_contours_from_db(vouchers_path = v_csv, contours_dir = file.path(sql_dir, "contours"))
            tick_idx <- which(db_ing$meta_df$catalogNumber == "V_TICK_1")
            tick_upright <- as.logical(db_ing$meta_df$mask_is_upright[tick_idx[order(db_ing$meta_df$leaf_id[tick_idx])]])
            run_fourier_extraction(parse_args_robust(c(
              "--contours-dir", file.path(sql_dir, "contours"),
              "--vouchers", v_csv,
              "--output", file.path(sql_dir, "tables", "leaf_efa_harmonics.csv"),
              "--output-individual", file.path(sql_dir, "tables", "leaf_efa_harmonics_individual.csv"),
              "--figures-dir", file.path(sql_dir, "figures")
            )))

            # 3. Single OTU PCA invariance + custom YAML + adaptive anchor floor
            Sys.unsetenv("PLANTMATICS_DB_PATH")
            Sys.unsetenv("RUN_ID")
            cfg_y <- "{cls.tmp_p / 'custom_step03.yaml'}"
            for (tag in c("pca1", "pca6")) {{
              pdir <- file.path("{cls.tmp_p}", tag)
              run_fourier_extraction(parse_args_robust(c(
                "--config", cfg_y,
                "--contours-dir", file.path(pdir, "contours"),
                "--manifest", file.path(pdir, "manifest.csv"),
                "--vouchers", file.path(pdir, "vouchers.csv"),
                "--output", file.path(pdir, "spec.csv"),
                "--output-individual", file.path(pdir, "ind.csv")
              )))
            }}
            set.seed(17)
            mk_tx <- function(tx, n_g, n_b, cx) {{
              n_tot <- n_g + n_b
              data.frame(
                catalogNumber = paste0(tx, "_", seq_len(n_tot)),
                species_standardized = tx,
                determiner_tier = c(rep("Tier_1_Gold", n_g), rep("Tier_3_Bronze", n_b)),
                PC1 = rnorm(n_tot, cx, 0.08), PC2 = rnorm(n_tot, 0.0, 0.08),
                stringsAsFactors = FALSE
              )
            }}
            df_anc <- rbind(mk_tx("Packera small", 2L, 1L, -2.0), mk_tx("Packera med", 2L, 10L, 2.0))
            anc_res <- select_adaptive_reference_anchors(df_anc, feature_matrix = as.matrix(df_anc[, c("PC1", "PC2")]), taxon_col = "species_standardized", min_anchors_per_taxon = 5L)
            w_anc <- compute_authority_precision_weights(df_anc, active_mask = anc_res$active_mask)

            # 4. Regional herbarium batch-centering & institution invariance
            n_h <- 24L
            df_herb <- data.frame(
              catalogNumber = paste0("H_", seq_len(n_h)),
              species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 12L),
              institutionCode = rep(c("NCU", "LSU"), times = 12L),
              regional_group = rep(c("Piedmont", "South_Central_Plains"), times = 12L),
              authority_weight = rep(2.0, n_h), is_passive = FALSE,
              PC1 = rnorm(n_h, rep(c(-1.0, 1.0), times = 12L), 0.1),
              PC2 = rnorm(n_h, 0.0, 0.1), PC3 = rnorm(n_h, 0.0, 0.1),
              ECT_PC1 = rnorm(n_h, rep(c(-0.5, 0.5), times = 12L), 0.1),
              solidity = rep(c(0.82, 0.90), times = 12L),
              blade_cie_l_norm = rnorm(n_h, 0.4, 0.02),
              blade_cie_b_delta = rnorm(n_h, rep(c(5.0, 12.0), times = 12L), 0.3),
              stringsAsFactors = FALSE
            )
            eb_conf <- apply_eb_institution_batch_centering(df_herb, target_cols = c("PC1", "PC2", "PC3", "ECT_PC1", "solidity", "blade_cie_l_norm", "blade_cie_b_delta"), min_inst_n = 5L)
            df_bal <- df_herb; df_bal$regional_group <- rep(c("Piedmont", "Piedmont", "South_Central_Plains", "South_Central_Plains"), length.out = n_h)
            eb_bal <- apply_eb_institution_batch_centering(df_bal, target_cols = c("PC1", "blade_cie_b_delta"), min_inst_n = 5L)

            df_eb_w <- data.frame(
              catalogNumber = paste0("W_", 1:12),
              species_standardized = rep(c("Packera aurea", "Packera dubia"), each = 6L),
              identifiedBy = rep(c("Spec A", "Spec B"), times = 6L),
              determiner_tier = "Tier_1_Gold",
              institutionCode = rep(c("NCU", "LSU"), times = 6L),
              vision_predicted_label = rep(c("Packera aurea", "Packera dubia"), each = 6L),
              PC1 = c(rnorm(6, -2, 0.1), rnorm(6, 2, 0.1)),
              stringsAsFactors = FALSE
            )
            w_i1 <- compute_empirical_bayes_authority_weights(df_eb_w, mode = "active_passive")
            df_eb_w2 <- df_eb_w; df_eb_w2$institutionCode <- ifelse(df_eb_w2$institutionCode == "NCU", "LSU", "NCU")
            w_i2 <- compute_empirical_bayes_authority_weights(df_eb_w2, mode = "active_passive")

            res_list <- list(
              z_zones = z_out$leaf_zone,
              sh_pids = sh_out$plant_individual_id,
              sh_zones = sh_out$leaf_zone,
              tick_upright = tick_upright,
              small_bronze_active = anc_res$active_mask[df_anc$species_standardized == "Packera small" & df_anc$determiner_tier == "Tier_3_Bronze"],
              small_bronze_wt = w_anc[df_anc$species_standardized == "Packera small" & df_anc$determiner_tier == "Tier_3_Bronze"],
              med_active_n = sum(anc_res$active_mask[df_anc$species_standardized == "Packera med"]),
              eb_shape_diffs = c(
                max(abs(eb_conf$PC1 - df_herb$PC1)), max(abs(eb_conf$PC2 - df_herb$PC2)),
                max(abs(eb_conf$PC3 - df_herb$PC3)), max(abs(eb_conf$ECT_PC1 - df_herb$ECT_PC1)),
                max(abs(eb_conf$solidity - df_herb$solidity)), max(abs(eb_conf$blade_cie_b_delta - df_herb$blade_cie_b_delta)),
                max(abs(eb_bal$PC1 - df_bal$PC1))
              ),
              eb_bal_phot_diff = max(abs(eb_bal$blade_cie_b_delta - df_bal$blade_cie_b_delta)),
              pre_eta_names = names(attr(eb_conf, "pre_eb_batch_partial_eta_sq")),
              inst_weight_max_diff = max(abs(w_i1 - w_i2))
            )
            writeLines(jsonlite::toJSON(res_list, auto_unbox = TRUE), "{cls.json_out}")
            """,
            encoding="utf-8",
        )
        proc = subprocess.run(["Rscript", str(r_script)], cwd=str(PROJECT_ROOT), capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"Batched Step 03 Rscript failed:\n{proc.stderr}\n{proc.stdout}")
        cls.data = json.loads(cls.json_out.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "_tmp_dir"):
            cls._tmp_dir.cleanup()

    def test_shape_free_leaf_zone_routing_and_proximity_sheet_consolidation(self) -> None:
        self.assertEqual(self.data["z_zones"], ["basal", "basal", "basal"])
        pids = self.data["sh_pids"]
        zones = self.data["sh_zones"]
        self.assertEqual(pids[2], 0)
        self.assertEqual(zones[2], "cauline")
        self.assertEqual(pids[3], 2)
        self.assertEqual(pids[6], 1)

    def test_scale_provenance_mm_gating_tipup_vectorization_and_photometry_handoff(self) -> None:
        self.assertEqual(self.data["tick_upright"], [True, False, False])
        sql_dir = self.tmp_p / "sql_run"
        ind_df = pd.read_csv(sql_dir / "tables" / "leaf_efa_harmonics_individual.csv")
        spec_df = pd.read_csv(sql_dir / "tables" / "leaf_efa_harmonics.csv")
        mm_cols = ["centroid_size", "log_centsize", "area_mm2", "petiole_length_mm", "blade_length_mm", "petiole_width_mm", "blade_width_mm"]
        for df_chk in (ind_df, spec_df):
            t_rows = df_chk[df_chk["catalogNumber"] == "V_TICK_1"]
            self.assertGreaterEqual(len(t_rows), 1)
            self.assertTrue((t_rows["scale_source"] == "tick_lattice_unboxed").all())
            self.assertTrue((t_rows["is_calibrated_mm"].astype(str).str.upper() == "FALSE").all())
            for col in mm_cols:
                self.assertTrue(t_rows[col].isna().all(), f"Expected NA for uncalibrated {col}")
            self.assertTrue(t_rows["petiole_to_blade_ratio"].notna().all())
        t_otu = spec_df[spec_df["catalogNumber"] == "V_TICK_1"].iloc[0]
        self.assertTrue(pd.isna(t_otu["cauline_blade_length_mm"]))
        self.assertAlmostEqual(float(t_otu["blade_cie_l_sheet_delta"]), -19.0, places=4)
        self.assertAlmostEqual(float(t_otu["blade_cie_a_delta"]), 3.4, places=4)
        self.assertAlmostEqual(float(t_otu["blade_cie_b_delta"]), 14.7, places=4)

        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        om = dict(conn.execute("SELECT * FROM otu_morphometrics WHERE run_id='run_s03_contracts' AND catalogNumber='V_TICK_1';").fetchone())
        conn.close()
        self.assertEqual(om["scale_source"], "tick_lattice_unboxed")
        self.assertEqual(int(om["is_calibrated_mm"]), 0)
        for col in mm_cols:
            self.assertIsNone(om[col], f"Expected SQL NULL for uncalibrated {col}")
        self.assertAlmostEqual(float(om["blade_cie_l_sheet_delta"]), -19.0, places=4)
        self.assertAlmostEqual(float(om["blade_cie_a_delta"]), 3.4, places=4)
        self.assertAlmostEqual(float(om["blade_cie_b_delta"]), 14.7, places=4)

    def test_single_otu_pca_invariance_adaptive_anchor_floor_and_yaml_config(self) -> None:
        s1 = pd.read_csv(self.tmp_p / "pca1" / "spec.csv").sort_values("catalogNumber")
        s6 = pd.read_csv(self.tmp_p / "pca6" / "spec.csv").sort_values("catalogNumber")
        self.assertIn("PC6", s1.columns)
        self.assertIn("PC6", s6.columns)
        for pc in [f"PC{i}" for i in range(1, 7)]:
            np.testing.assert_allclose(s1[pc].to_numpy(), s6[pc].to_numpy(), atol=1e-6)
        self.assertFalse(bool(self.data["small_bronze_active"]))
        self.assertEqual(float(self.data["small_bronze_wt"]), 0.0)
        self.assertEqual(int(self.data["med_active_n"]), 5)

    def test_regional_herbarium_centering_safety_and_zero_duplicate_artifacts(self) -> None:
        for d in self.data["eb_shape_diffs"]:
            self.assertLessEqual(float(d), 1e-12)
        self.assertGreater(float(self.data["eb_bal_phot_diff"]), 0.1)
        for col in ("PC1", "PC2", "PC3", "ECT_PC1", "solidity", "blade_cie_b_delta"):
            self.assertIn(col, self.data["pre_eta_names"])
        self.assertLessEqual(float(self.data["inst_weight_max_diff"]), 1e-12)
        sql_dir = self.tmp_p / "sql_run"
        self.assertFalse((sql_dir / "tables" / "fourier_pcs.csv").exists())
        self.assertFalse((sql_dir / "tables" / "symmetric_fourier_pcs.csv").exists())
        self.assertEqual(list((sql_dir / "figures").glob("*.pdf")), [])
        self.assertFalse((PROJECT_ROOT / "Rplots.pdf").exists())


if __name__ == "__main__":
    unittest.main()

