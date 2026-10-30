#!/usr/bin/env python3
"""
===============================================================================
Test Suite: test_sql_pipeline_handshake.py
Project: Packera Species Delimitation & Morphometrics Pipeline (NCU)

Description:
    End-to-end relational contract & cross-language handshake tests for the
    Two-Tier SQLite (WAL mode) + DuckDB architecture (`scripts/core/db.py` &
    `scripts/core/db_utils.R`):
      1. DDL schema initialization, WAL mode, FOREIGN KEY / ON DELETE CASCADE
         constraints, R*Tree spatial triggers, DuckDB attachment, and --demo
         isolation.
      2. Round-trip migration of `curated_vouchers.csv`,
         `extracted_leaf_manifest.csv`, and 120-vertex contours into
         `extracted_leaves`, followed by Step 03 EFA execution with zero
         filesystem `.csv` contour files.
      3. `dinov2_feature_cache` hit/miss & SHA-256 invalidation behavior in
         `05_cleanlab_vision_xai.py` (confirming zero redundant model loads or
         forward passes on cache hit, and zero DB writes under `--demo`).
      4. Multi-plant sheet preservation `(catalogNumber, plant_individual_id)`
         and non-georeferenced historic voucher retention across
         `v_integrative_triage_master`, Step 06, and Step 07.
===============================================================================
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Tuple
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image
import torch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    DEFAULT_SQLITE_REL_PATH,
    ensure_pipeline_run,
    fetch_cached_dinov2_features,
    get_db_connection,
    get_duckdb_connection,
    ingest_existing_csvs_to_db,
    persist_vision_xai_audit,
    populate_run_active_cohort,
    query_environmental_bbox,
    resolve_db_path,
    store_dinov2_features_cache,
    upsert_environmental_cache_records,
    upsert_rows_generic,
    upsert_vouchers_records,
)


def _load_step05_module():
    script_path = PROJECT_ROOT / "scripts" / "analysis" / "05_cleanlab_vision_xai.py"
    spec = importlib.util.spec_from_file_location("step05_cleanlab_sql_test", script_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _make_closed_leaf_contour(
    n_pts: int = 120,
    rx: float = 95.0,
    ry: float = 58.0,
    taper: float = 0.14,
    serration: float = 1.8,
) -> List[Tuple[float, float]]:
    """Generate a smooth closed 120-vertex botanical lamina contour."""
    pts: List[Tuple[float, float]] = []
    for idx in range(n_pts):
        theta = 2.0 * math.pi * (idx / float(n_pts))
        r_mod = 1.0 + (serration / 100.0) * math.cos(8.0 * theta)
        x = 220.0 + rx * r_mod * math.cos(theta)
        width_scale = 1.0 + taper * math.cos(theta)
        y = 220.0 + ry * r_mod * width_scale * math.sin(theta)
        pts.append((round(x, 4), round(y, 4)))
    return pts


class _CountingFakeDINOv2(torch.nn.Module):
    """Deterministic fake DINOv2 module counting `forward_features` invocations."""

    def __init__(self) -> None:
        super().__init__()
        self.forward_calls = 0
        self.batch_sizes: List[int] = []

    def forward_features(self, x: torch.Tensor) -> Dict[str, torch.Tensor]:
        self.forward_calls += 1
        b = int(x.shape[0])
        self.batch_sizes.append(b)
        # Deterministic spatial patch tokens (B, 256, 768) derived from input pixels
        pooled_16 = torch.nn.functional.adaptive_avg_pool2d(x.float(), (16, 16))  # (B, 3, 16, 16)
        flat_rgb = pooled_16.flatten(2).transpose(1, 2)  # (B, 256, 3)
        patch_tokens = flat_rgb.repeat(1, 1, 256)  # (B, 256, 768)
        cls_token = patch_tokens.mean(dim=1)  # (B, 768)
        return {
            "x_norm_patchtokens": patch_tokens,
            "x_norm_clstoken": cls_token,
        }


class TestSQLPipelineHandshake(unittest.TestCase):
    """Integration tests verifying the central SQLite (WAL) & DuckDB architecture."""

    def test_ddl_schema_wal_foreign_keys_and_rtree_triggers(self) -> None:
        """Verify WAL mode, Tier A/B DDL, FK cascade rules, R*Tree spatial index, DuckDB, and --demo isolation."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_dir = Path(tmp_d)
            db_file = tmp_dir / "packera_test.sqlite"

            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                journal_mode = str(conn.execute("PRAGMA journal_mode;").fetchone()[0]).lower()
                fk_enabled = int(conn.execute("PRAGMA foreign_keys;").fetchone()[0])
                self.assertEqual(journal_mode, "wal")
                self.assertEqual(fk_enabled, 1)

                obj_rows = conn.execute(
                    "SELECT name, type FROM sqlite_master WHERE type IN ('table', 'view');"
                ).fetchall()
                existing_objects = {str(r[0]): str(r[1]) for r in obj_rows}

                expected_tables = {
                    "taxon_synonyms",
                    "vouchers",
                    "plant_individuals",
                    "extracted_leaves",
                    "voucher_reproductive_phenology",
                    "dinov2_feature_cache",
                    "environmental_cache",
                    "environmental_rtree",
                    "pipeline_runs",
                    "run_active_cohort",
                    "leaf_efa_harmonics",
                    "otu_morphometrics",
                    "vision_xai_audit",
                    "spatial_rf_predictions",
                    "curation_triage",
                }
                for tbl in expected_tables:
                    self.assertIn(tbl, existing_objects, f"Missing SQLite table: {tbl}")
                self.assertEqual(existing_objects.get("v_integrative_triage_master"), "view")

                # Foreign key enforcement: orphan plant_individuals and extracted_leaves must fail
                with self.assertRaises(sqlite3.IntegrityError):
                    with conn:
                        conn.execute(
                            "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('ORPHAN-1', 0);"
                        )
                with self.assertRaises(sqlite3.IntegrityError):
                    with conn:
                        conn.execute(
                            """
                            INSERT INTO extracted_leaves (
                                catalogNumber, plant_individual_id, leaf_id, assigned_tier
                            ) VALUES ('ORPHAN-1', 0, 1, 'Tier_1_Direct');
                            """
                        )

                # Valid parent -> child hierarchy + ON DELETE CASCADE
                upsert_vouchers_records(
                    conn,
                    [
                        {
                            "catalogNumber": "NCU-FK-01",
                            "institutionCode": "NCU",
                            "scientificName": "Packera aurea",
                            "species_standardized": "Packera aurea",
                            "determiner_tier": "Tier_1_Gold",
                            "decimalLatitude": 35.9132,
                            "decimalLongitude": -79.0558,
                        }
                    ],
                )
                with conn:
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU-FK-01', 0);"
                    )
                    conn.execute(
                        """
                        INSERT INTO extracted_leaves (
                            catalogNumber, plant_individual_id, leaf_id, assigned_tier, contour_xy_json, contour_n_pts
                        ) VALUES ('NCU-FK-01', 0, 1, 'Tier_1_Direct', '[[0,0],[1,0],[1,1]]', 3);
                        """
                    )
                with conn:
                    conn.execute("DELETE FROM vouchers WHERE catalogNumber = 'NCU-FK-01';")
                rem_pi = conn.execute(
                    "SELECT COUNT(*) FROM plant_individuals WHERE catalogNumber = 'NCU-FK-01';"
                ).fetchone()[0]
                rem_el = conn.execute(
                    "SELECT COUNT(*) FROM extracted_leaves WHERE catalogNumber = 'NCU-FK-01';"
                ).fetchone()[0]
                self.assertEqual(int(rem_pi), 0)
                self.assertEqual(int(rem_el), 0)

                # R*Tree spatial index trigger & bounding-box query
                upsert_vouchers_records(
                    conn,
                    [
                        {"catalogNumber": "NCU-GEO-NC", "scientificName": "Packera aurea"},
                        {"catalogNumber": "NCU-GEO-ME", "scientificName": "Packera paupercula"},
                    ],
                )
                upsert_environmental_cache_records(
                    conn,
                    [
                        {
                            "catalogNumber": "NCU-GEO-NC",
                            "latitude": 35.91320,
                            "longitude": -79.05580,
                            "fine_soil_ph": 5.6,
                            "sgmc_lithology_class": "Felsic_Crystalline",
                        },
                        {
                            "catalogNumber": "NCU-GEO-ME",
                            "latitude": 44.80160,
                            "longitude": -68.77120,
                            "fine_soil_ph": 6.2,
                            "sgmc_lithology_class": "Carbonate",
                        },
                    ],
                    demo=False,
                )
                nc_hits = query_environmental_bbox(conn, min_lon=-80.5, max_lon=-78.0, min_lat=35.0, max_lat=36.5)
                self.assertEqual(len(nc_hits), 1)
                self.assertEqual(nc_hits[0]["catalogNumber"], "NCU-GEO-NC")
                self.assertAlmostEqual(float(nc_hits[0]["fine_soil_ph"]), 5.6, places=3)

                # Verify --demo strictly blocks writes to environmental_cache and isolates resolve_db_path
                n_demo_env = upsert_environmental_cache_records(
                    conn,
                    [{"catalogNumber": "DEMO-1", "latitude": 30.0, "longitude": -80.0, "fine_soil_ph": 7.0}],
                    demo=True,
                )
                self.assertEqual(n_demo_env, 0)
                canonical_resolved = (PROJECT_ROOT / DEFAULT_SQLITE_REL_PATH).resolve()
                demo_resolved = resolve_db_path(db_path=canonical_resolved, workspace_root=PROJECT_ROOT, demo=True)
                self.assertNotEqual(str(demo_resolved), str(canonical_resolved))
            finally:
                conn.close()

            # Verify R db_utils.R spatial R*Tree bounding-box helper against the same SQLite DB
            rscript_bin = shutil.which("Rscript")
            self.assertIsNotNone(rscript_bin, "Rscript must be installed in PATH")
            r_expr = (
                "source('scripts/core/db_utils.R'); "
                f"df <- query_environmental_cache_sqlite(db_path = '{db_file.as_posix()}', "
                "bbox = c(-80.5, -78.0, 35.0, 36.5)); "
                "cat('R_BBOX_N=', nrow(df), 'CAT=', as.character(df$catalogNumber[1]), sep='')"
            )
            r_proc = subprocess.run([rscript_bin, "-e", r_expr], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True)
            self.assertIn("R_BBOX_N=1CAT=NCU-GEO-NC", r_proc.stdout)

            # Verify optional DuckDB analytical attachment if duckdb is installed in .venv
            if importlib.util.find_spec("duckdb") is not None:
                dconn = get_duckdb_connection(db_path=db_file, read_only=True)
                try:
                    cnt = dconn.execute("SELECT COUNT(*) FROM vouchers;").fetchone()[0]
                    self.assertEqual(int(cnt), 2)
                finally:
                    dconn.close()

    def test_contour_migration_and_step03_efa_zero_csv_reads(self) -> None:
        """Migrate CSVs + 120-vertex contours to SQLite, delete contour CSVs, and run Step 03 EFA purely from SQLite."""
        rscript_bin = shutil.which("Rscript")
        self.assertIsNotNone(rscript_bin, "Rscript must be installed in PATH")

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d) / "ws"
            tables_dir = ws / "data" / "tables"
            contours_dir = ws / "data" / "contours"
            reports_dir = ws / "outputs" / "reports"
            figures_dir = ws / "outputs" / "figures"
            for d in (tables_dir, contours_dir, reports_dir, figures_dir):
                d.mkdir(parents=True, exist_ok=True)

            db_path = ws / "data" / "packera_pipeline.sqlite"
            run_id = "test_efa_sql_zero_csv"

            voucher_rows: List[Dict[str, Any]] = []
            manifest_rows: List[Dict[str, Any]] = []

            # 2 taxa x 5 vouchers, with NCU-A1 having 2 distinct plant individuals (pid=0 and pid=1)
            specs: List[Tuple[str, str, List[int]]] = [
                ("NCU-A1", "Packera aurea", [0, 1]),
                ("NCU-A2", "Packera aurea", [0]),
                ("NCU-A3", "Packera aurea", [0]),
                ("NCU-A4", "Packera aurea", [0]),
                ("NCU-A5", "Packera aurea", [0]),
                ("NCU-B1", "Packera obovata", [0]),
                ("NCU-B2", "Packera obovata", [0]),
                ("NCU-B3", "Packera obovata", [0]),
                ("NCU-B4", "Packera obovata", [0]),
                ("NCU-B5", "Packera obovata", [0]),
            ]

            for v_idx, (cat, taxon, pids) in enumerate(specs):
                voucher_rows.append(
                    {
                        "catalogNumber": cat,
                        "institutionCode": "NCU" if v_idx % 2 == 0 else "US",
                        "scientificName": taxon,
                        "species_raw": taxon,
                        "species_standardized": taxon,
                        "identifiedBy": "T.M. Barkley",
                        "determiner_raw": "T.M. Barkley",
                        "determiner_tier": "Tier_1_Gold",
                        "type_status": "",
                        "county": "Orange",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": round(35.5 + 0.05 * v_idx, 5),
                        "latitude": round(35.5 + 0.05 * v_idx, 5),
                        "decimalLongitude": round(-79.5 + 0.05 * v_idx, 5),
                        "longitude": round(-79.5 + 0.05 * v_idx, 5),
                        "coordinateUncertainty": 250.0,
                        "year": 2019,
                        "month": 4,
                        "day": 18,
                        "eventDate": "2019-04-18",
                        "regional_group": "Piedmont",
                        "image_path": f"data/raw_vouchers/{cat}.jpg",
                        "exsiccatae_key": f"barkley_{v_idx}",
                        "is_primary_duplicate": "True",
                        "duplicate_count": 1,
                    }
                )
                is_aurea = taxon == "Packera aurea"
                for pid in pids:
                    for local_lid in (1, 2):
                        leaf_id = pid * 10 + local_lid
                        c_name = f"{cat}_p{pid}_leaf{leaf_id}.csv"
                        c_file = contours_dir / c_name
                        rx = (92.0 + 2.0 * local_lid + 4.0 * pid) if is_aurea else (108.0 + 2.0 * local_lid)
                        ry = (74.0 + 1.5 * local_lid - 3.0 * pid) if is_aurea else (48.0 + 1.5 * local_lid)
                        taper = 0.08 if is_aurea else -0.18
                        pts = _make_closed_leaf_contour(120, rx=rx, ry=ry, taper=taper, serration=1.2)
                        with open(c_file, "w", encoding="utf-8", newline="") as cf:
                            w = csv.writer(cf)
                            w.writerow(["x", "y"])
                            for px, py in pts:
                                w.writerow([px, py])

                        manifest_rows.append(
                            {
                                "catalogNumber": cat,
                                "plant_individual_id": pid,
                                "leaf_id": leaf_id,
                                "assigned_tier": "Tier_1_Direct",
                                "solidity": 0.91,
                                "circularity": 0.68 if is_aurea else 0.54,
                                "extent": 0.73,
                                "ucs_score": 0.94,
                                "area_px": 18500.0,
                                "pixels_per_mm": 10.0,
                                "is_calibrated_mm": "True",
                                "petiole_length_px": 280.0 if is_aurea else 95.0,
                                "petiole_length_mm": (28.0 + 5.0 * pid) if is_aurea else 9.5,
                                "blade_length_mm": 32.0 if is_aurea else 41.0,
                                "petiole_to_blade_ratio": round(((28.0 + 5.0 * pid) / 32.0) if is_aurea else (9.5 / 41.0), 4),
                                "is_dissected": "False",
                                "contour_path": str(c_file),
                            }
                        )

            vouchers_csv = tables_dir / "curated_vouchers.csv"
            with open(vouchers_csv, "w", encoding="utf-8", newline="") as vf:
                vw = csv.DictWriter(vf, fieldnames=list(voucher_rows[0].keys()))
                vw.writeheader()
                vw.writerows(voucher_rows)

            manifest_csv = tables_dir / "extracted_leaf_manifest.csv"
            with open(manifest_csv, "w", encoding="utf-8", newline="") as mf:
                mw = csv.DictWriter(mf, fieldnames=list(manifest_rows[0].keys()))
                mw.writeheader()
                mw.writerows(manifest_rows)

            mig_counts = ingest_existing_csvs_to_db(
                workspace_root=ws,
                run_id=run_id,
                db_path=db_path,
                tables_dir=tables_dir,
                contours_dir=contours_dir,
            )
            self.assertEqual(mig_counts["vouchers"], 10)
            self.assertEqual(mig_counts["plant_individuals"], 11)
            self.assertEqual(mig_counts["extracted_leaves"], 22)
            self.assertEqual(mig_counts["run_active_cohort"], 11)

            # Delete ALL .csv contour files from disk to prove Step 03 reads strictly from SQLite!
            shutil.rmtree(contours_dir)
            contours_dir.mkdir(parents=True, exist_ok=True)
            self.assertEqual(list(contours_dir.glob("*.csv")), [])

            out_otu_csv = tables_dir / "leaf_efa_harmonics.csv"
            out_ind_csv = tables_dir / "leaf_efa_harmonics_individual.csv"
            out_rep_csv = reports_dir / "tier_symmetry_validation.csv"
            out_fig_pdf = figures_dir / "tier1_vs_tier2_density_overlay.pdf"

            env = os.environ.copy()
            env["PACKERA_DB_PATH"] = str(db_path)
            env["RUN_ID"] = run_id
            env["PACKERA_RUN_DIR"] = str(ws)

            cmd = [
                rscript_bin,
                str(PROJECT_ROOT / "scripts" / "morphometrics" / "03_fourier_extractor.R"),
                "--input", str(contours_dir),
                "--contours-dir", str(contours_dir),
                "--manifest", str(manifest_csv),
                "--vouchers", str(vouchers_csv),
                "--output", str(out_otu_csv),
                "--output-individual", str(out_ind_csv),
                "--report-out", str(out_rep_csv),
                "--plot-out", str(out_fig_pdf),
                "--permutations", "19",
            ]
            proc = subprocess.run(cmd, cwd=PROJECT_ROOT, env=env, capture_output=True, text=True, check=False)
            combined_log = proc.stdout + "\n" + proc.stderr
            self.assertEqual(proc.returncode, 0, f"Step 03 failed with empty contours_dir:\n{combined_log}")
            self.assertIn("from SQLite extracted_leaves", combined_log)
            self.assertTrue(out_otu_csv.is_file())

            with open(out_otu_csv, "r", encoding="utf-8", newline="") as f:
                otu_rows = list(csv.DictReader(f))
            self.assertEqual(len(otu_rows), 11)
            otu_keys = {(r["catalogNumber"], int(float(r["plant_individual_id"]))) for r in otu_rows}
            self.assertIn(("NCU-A1", 0), otu_keys)
            self.assertIn(("NCU-A1", 1), otu_keys)
            for r in otu_rows:
                self.assertEqual(r["morph_tier_source"], "Tier_1_Pure")
                self.assertTrue(math.isfinite(float(r["Sym_PC1"])))

            # Verify SQLite persistence in leaf_efa_harmonics and otu_morphometrics
            conn = get_db_connection(db_path=db_path, read_only=True)
            try:
                n_leaf_db = conn.execute(
                    "SELECT COUNT(*) FROM leaf_efa_harmonics WHERE run_id = ?;", (run_id,)
                ).fetchone()[0]
                n_otu_db = conn.execute(
                    "SELECT COUNT(*) FROM otu_morphometrics WHERE run_id = ?;", (run_id,)
                ).fetchone()[0]
                self.assertEqual(int(n_leaf_db), 22)
                self.assertEqual(int(n_otu_db), 11)
            finally:
                conn.close()

    def test_dinov2_feature_cache_hit_miss_and_demo_isolation(self) -> None:
        """Verify Step 05 dinov2_feature_cache hit/miss, SHA-256 invalidation, and --demo isolation."""
        step05 = _load_step05_module()

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_dir = Path(tmp_d)
            db_path = tmp_dir / "vision_cache.sqlite"
            demo_db_path = tmp_dir / "demo_vision.sqlite"
            patches_dir = tmp_dir / "basal_leaves_raw"
            masks_dir = tmp_dir / "tier1_intact"
            patches_dir.mkdir(parents=True, exist_ok=True)
            masks_dir.mkdir(parents=True, exist_ok=True)

            # Seed parent vouchers + multi-plant individuals in SQLite
            conn = get_db_connection(db_path=db_path, init_schema=True)
            try:
                upsert_vouchers_records(
                    conn,
                    [{"catalogNumber": "NCU-VIS-01", "scientificName": "Packera aurea"}],
                )
                with conn:
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU-VIS-01', 0);"
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('NCU-VIS-01', 1);"
                    )
            finally:
                conn.close()

            records: List[Dict[str, Any]] = []
            for pid, rgb in ((0, (70, 125, 65)), (1, (95, 140, 80))):
                p_img = patches_dir / f"NCU-VIS-01_p{pid}_leaf1.jpg"
                p_msk = masks_dir / f"NCU-VIS-01_p{pid}_leaf1.png"
                arr = np.full((96, 96, 3), [242, 239, 230], dtype=np.uint8)
                arr[16:80, 16:80] = rgb
                Image.fromarray(arr, mode="RGB").save(p_img, quality=95)

                m_arr = np.zeros((96, 96), dtype=np.uint8)
                m_arr[16:80, 16:80] = 255
                Image.fromarray(m_arr, mode="L").save(p_msk)

                records.append(
                    {
                        "catalogNumber": "NCU-VIS-01",
                        "plant_individual_id": pid,
                        "label_idx": pid,
                        "patch_path": str(p_img),
                        "mask_path": str(p_msk),
                        "tier1_leaf_patches": [str(p_img)],
                        "tier1_leaf_masks": [str(p_msk)],
                        "tier1_leaf_weights": [1.0],
                        "surface_orientation": "adaxial",
                        "mean_luminance": 48.5 + pid,
                        "mean_saturation": 0.35,
                        "high_brightness_ratio": 0.04,
                    }
                )

            # 1. Cache MISS: first call invokes fake DINOv2 once and populates dinov2_feature_cache
            fake_model = _CountingFakeDINOv2()
            with patch.object(step05.torch.hub, "load", return_value=fake_model) as mock_hub:
                feats_1, labels_1, order_1, fg_1 = step05.extract_dinov2_embeddings(
                    records,
                    model_name="dinov2_vitb14_reg",
                    device="cpu",
                    batch_size=8,
                    demo=False,
                    return_fg_ratios=True,
                    db_path=db_path,
                    sigma=1.0,
                )
            self.assertEqual(mock_hub.call_count, 1)
            self.assertEqual(fake_model.forward_calls, 1)
            self.assertEqual(fake_model.batch_sizes, [2])
            self.assertEqual(feats_1.shape, (2, 768))
            self.assertEqual(order_1, [("NCU-VIS-01", 0), ("NCU-VIS-01", 1)])
            self.assertTrue(np.all(np.isfinite(fg_1)))

            conn_ro = get_db_connection(db_path=db_path, read_only=True)
            try:
                n_cached = conn_ro.execute("SELECT COUNT(*) FROM dinov2_feature_cache;").fetchone()[0]
                self.assertEqual(int(n_cached), 2)
                cached_map = fetch_cached_dinov2_features(
                    conn_ro,
                    [("NCU-VIS-01", 0, None), ("NCU-VIS-01", 1, None)],
                    model_name="dinov2_vitb14_reg",
                    sigma=1.0,
                )
                self.assertEqual(len(cached_map), 2)
                self.assertEqual(len(cached_map[("NCU-VIS-01", 0)]["embedding"]), 768)
            finally:
                conn_ro.close()

            # 2. Full Cache HIT: patching torch.hub.load to raise proves zero model loads & zero forward passes
            with patch.object(
                step05.torch.hub,
                "load",
                side_effect=AssertionError("DINOv2 model must NOT be loaded when all specimens hit SQLite cache!"),
            ):
                feats_2, labels_2, order_2, fg_2 = step05.extract_dinov2_embeddings(
                    records,
                    model_name="dinov2_vitb14_reg",
                    device="cpu",
                    batch_size=8,
                    demo=False,
                    return_fg_ratios=True,
                    db_path=db_path,
                    sigma=1.0,
                )
            self.assertEqual(order_2, order_1)
            np.testing.assert_array_equal(labels_2, labels_1)
            np.testing.assert_allclose(feats_2, feats_1, atol=1e-5)
            np.testing.assert_allclose(fg_2, fg_1, atol=1e-5)

            # 3. Partial SHA-256 invalidation: modifying pid=1 image re-runs strictly batch_size == 1
            p_img_pid1 = Path(records[1]["patch_path"])
            arr_mod = np.full((96, 96, 3), [242, 239, 230], dtype=np.uint8)
            arr_mod[16:80, 16:80] = (140, 90, 50)
            Image.fromarray(arr_mod, mode="RGB").save(p_img_pid1, quality=95)

            fake_model_partial = _CountingFakeDINOv2()
            with patch.object(step05.torch.hub, "load", return_value=fake_model_partial):
                feats_3, _, _, _ = step05.extract_dinov2_embeddings(
                    records,
                    model_name="dinov2_vitb14_reg",
                    device="cpu",
                    batch_size=8,
                    demo=False,
                    return_fg_ratios=True,
                    db_path=db_path,
                    sigma=1.0,
                )
            self.assertEqual(fake_model_partial.forward_calls, 1)
            self.assertEqual(fake_model_partial.batch_sizes, [1])
            np.testing.assert_allclose(feats_3[0], feats_1[0], atol=1e-5)

            # 4. --demo isolation: never persists synthetic embeddings or audit rows into SQLite
            conn_demo = get_db_connection(db_path=demo_db_path, init_schema=True)
            try:
                n_store = store_dinov2_features_cache(
                    conn_demo,
                    [{"catalogNumber": "X", "embedding": [0.1] * 768}],
                    model_name="dinov2_vitb14_reg",
                    sigma=1.0,
                    demo=True,
                )
                n_audit = persist_vision_xai_audit(conn_demo, "demo_run", [{"catalogNumber": "X"}], demo=True)
                self.assertEqual(n_store, 0)
                self.assertEqual(n_audit, 0)
            finally:
                conn_demo.close()

            with patch.object(step05.torch.hub, "load", side_effect=RuntimeError("offline")):
                _ = step05.extract_dinov2_embeddings(
                    records,
                    device="cpu",
                    demo=True,
                    db_path=demo_db_path,
                )
            conn_demo_ro = get_db_connection(db_path=demo_db_path, read_only=True)
            try:
                cnt_demo = conn_demo_ro.execute("SELECT COUNT(*) FROM dinov2_feature_cache;").fetchone()[0]
                self.assertEqual(int(cnt_demo), 0)
            finally:
                conn_demo_ro.close()

    def test_multi_plant_and_non_georeferenced_retention_across_step06_and_step07(self) -> None:
        """Verify multi-plant (catalogNumber, plant_individual_id) & non-georeferenced retention across Step 06, Step 07, and v_integrative_triage_master."""
        rscript_bin = shutil.which("Rscript")
        self.assertIsNotNone(rscript_bin, "Rscript must be installed in PATH")

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d) / "ws_master"
            tables_dir = ws / "data" / "tables"
            env_dir = ws / "data" / "environmental"
            reports_dir = ws / "outputs" / "reports"
            figures_dir = ws / "outputs" / "figures"
            for d in (tables_dir, env_dir, reports_dir, figures_dir):
                d.mkdir(parents=True, exist_ok=True)

            db_path = ws / "data" / "packera_pipeline.sqlite"
            run_id = "test_master_multi_nogeo"

            # 10 vouchers (including multi-plant MULTI-001 with pid=0 and pid=1, plus non-georeferenced HIST-NOGEO-99)
            # Gives 11 total (catalogNumber, plant_individual_id) specimen rows.
            voucher_specs: List[Dict[str, Any]] = [
                {
                    "catalogNumber": "MULTI-001",
                    "institutionCode": "NCU",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_2_Silver",
                    "decimalLatitude": 35.91000,
                    "latitude": 35.91000,
                    "decimalLongitude": -79.05000,
                    "longitude": -79.05000,
                    "coordinateUncertainty": 200.0,
                    "year": 2020,
                    "month": 4,
                    "day": 20,
                    "eventDate": "2020-04-20",
                    "doy": 111,
                    "stateProvince": "North Carolina",
                    "regional_group": "Piedmont",
                },
                {
                    "catalogNumber": "HIST-NOGEO-99",
                    "institutionCode": "NY",
                    "scientificName": "Packera obovata",
                    "species_raw": "Packera obovata",
                    "species_standardized": "Packera obovata",
                    "determiner_tier": "Tier_3_Bronze",
                    "decimalLatitude": None,
                    "latitude": None,
                    "decimalLongitude": None,
                    "longitude": None,
                    "coordinateUncertainty": None,
                    "year": 1898,
                    "month": 5,
                    "day": 4,
                    "eventDate": "1898-05-04",
                    "doy": 124,
                    "stateProvince": "Virginia",
                    "regional_group": "",
                },
            ]
            for i in range(1, 6):
                voucher_specs.append(
                    {
                        "catalogNumber": f"ANC-AUR-{i}",
                        "institutionCode": "NCU" if i % 2 == 0 else "DUKE",
                        "scientificName": "Packera aurea",
                        "species_raw": "Packera aurea",
                        "species_standardized": "Packera aurea",
                        "determiner_tier": "Tier_1_Gold",
                        "decimalLatitude": round(35.60 + 0.08 * i, 5),
                        "latitude": round(35.60 + 0.08 * i, 5),
                        "decimalLongitude": round(-79.40 + 0.08 * i, 5),
                        "longitude": round(-79.40 + 0.08 * i, 5),
                        "coordinateUncertainty": 150.0,
                        "year": 2018,
                        "month": 4,
                        "day": 15 + i,
                        "eventDate": f"2018-04-{15 + i:02d}",
                        "doy": 105 + i,
                        "stateProvince": "North Carolina",
                        "regional_group": "Piedmont",
                    }
                )
                voucher_specs.append(
                    {
                        "catalogNumber": f"ANC-OBO-{i}",
                        "institutionCode": "US" if i % 2 == 0 else "MO",
                        "scientificName": "Packera obovata",
                        "species_raw": "Packera obovata",
                        "species_standardized": "Packera obovata",
                        "determiner_tier": "Tier_1_Gold",
                        "decimalLatitude": round(37.20 + 0.08 * i, 5),
                        "latitude": round(37.20 + 0.08 * i, 5),
                        "decimalLongitude": round(-80.40 + 0.08 * i, 5),
                        "longitude": round(-80.40 + 0.08 * i, 5),
                        "coordinateUncertainty": 150.0,
                        "year": 2018,
                        "month": 4,
                        "day": 22 + i,
                        "eventDate": f"2018-04-{22 + i:02d}",
                        "doy": 112 + i,
                        "stateProvince": "Virginia",
                        "regional_group": "Ridge and Valley",
                    }
                )

            vouchers_csv = tables_dir / "curated_vouchers.csv"
            v_cols = list(voucher_specs[0].keys())
            with open(vouchers_csv, "w", encoding="utf-8", newline="") as vf:
                vw = csv.DictWriter(vf, fieldnames=v_cols)
                vw.writeheader()
                for row in voucher_specs:
                    vw.writerow({k: ("" if v is None else v) for k, v in row.items()})

            # Individual-level tuples: MULTI-001 has pid=0 (Packera aurea) and pid=1 (misidentified Packera obovata individual mounted on same sheet!)
            otu_specs: List[Tuple[str, int, str, str]] = [
                ("MULTI-001", 0, "Packera aurea", "Packera aurea"),
                ("MULTI-001", 1, "Packera aurea", "Packera obovata"),
                ("HIST-NOGEO-99", 0, "Packera obovata", "Packera obovata"),
            ]
            for i in range(1, 6):
                otu_specs.append((f"ANC-AUR-{i}", 0, "Packera aurea", "Packera aurea"))
                otu_specs.append((f"ANC-OBO-{i}", 0, "Packera obovata", "Packera obovata"))

            conn = get_db_connection(db_path=db_path, init_schema=True)
            try:
                upsert_vouchers_records(conn, voucher_specs)
                with conn:
                    for cat, pid, _, _ in otu_specs:
                        conn.execute(
                            """
                            INSERT INTO plant_individuals (
                                catalogNumber, plant_individual_id,
                                petiole_length_mm, blade_length_mm, petiole_to_blade_ratio, plant_height_mm
                            ) VALUES (?, ?, ?, ?, ?, ?);
                            """,
                            (
                                cat,
                                pid,
                                30.0 if pid == 0 else 11.0,
                                34.0 if pid == 0 else 42.0,
                                round((30.0 / 34.0) if pid == 0 else (11.0 / 42.0), 4),
                                360.0 + 25.0 * pid,
                            ),
                        )
                        conn.execute(
                            """
                            INSERT INTO voucher_reproductive_phenology (
                                catalogNumber, plant_individual_id, capitula_count, n_anthesis,
                                n_bud, n_fruiting, voucher_phenological_state
                            ) VALUES (?, ?, 12, 8, 2, 2, 'anthesis');
                            """,
                            (cat, pid),
                        )

                ensure_pipeline_run(conn, run_id=run_id, run_dir=ws)
                n_coh = populate_run_active_cohort(conn, run_id=run_id, voucher_rows=voucher_specs, clear_existing=True)
                self.assertEqual(n_coh, 13)

                otu_rows: List[Dict[str, Any]] = []
                vis_rows: List[Dict[str, Any]] = []
                env_rows: List[Dict[str, Any]] = []
                v_by_cat = {r["catalogNumber"]: r for r in voucher_specs}

                for idx_o, (cat, pid, given_sp, true_sp) in enumerate(otu_specs):
                    is_aur = true_sp == "Packera aurea"
                    second_sp = "Packera obovata" if is_aur else "Packera aurea"
                    otu_rows.append(
                        {
                            "run_id": run_id,
                            "catalogNumber": cat,
                            "plant_individual_id": pid,
                            "morph_tier_source": "Tier_1_Pure",
                            "n_tier1_leaves": 3,
                            "n_tier2_leaves": 0,
                            "PC1": (-1.6 + 0.05 * idx_o) if is_aur else (1.7 + 0.05 * idx_o),
                            "PC2": 0.12 * ((idx_o % 3) - 1),
                            "PC3": 0.05 * ((idx_o % 2) - 0.5),
                            "Sym_PC1": (-1.6 + 0.05 * idx_o) if is_aur else (1.7 + 0.05 * idx_o),
                            "Sym_PC2": 0.12 * ((idx_o % 3) - 1),
                            "can1": (-2.4 + 0.06 * idx_o) if is_aur else (2.5 + 0.06 * idx_o),
                            "can2": 0.15 * ((idx_o % 3) - 1),
                            "can3": 0.05,
                            "gmm_cluster": 1 if is_aur else 2,
                            "gmm_uncertainty": 0.04,
                            "cda_predicted_taxon": true_sp,
                            "cda_posterior_prob": 0.94,
                            "cda_second_taxon": second_sp,
                            "cda_second_posterior_prob": 0.06,
                            "cda_posterior_margin": 0.88,
                            "loocv_accuracy": 0.95,
                            "mahalanobis_d2": 1.4 if given_sp == true_sp else 14.2,
                            "morphometric_outlier_flag": 0 if given_sp == true_sp else 1,
                            "misidentification_flag": 0 if given_sp == true_sp else 1,
                        }
                    )
                    vis_rows.append(
                        {
                            "catalogNumber": cat,
                            "plant_individual_id": pid,
                            "given_label": given_sp,
                            "predicted_label": true_sp,
                            "vision_predicted_label": true_sp,
                            "vision_predicted_prob": 0.93,
                            "vision_second_taxon": second_sp,
                            "vision_second_prob": 0.07,
                            "vision_margin": 0.86,
                            "self_confidence": 0.93 if given_sp == true_sp else 0.07,
                            "c_error": 0.07 if given_sp == true_sp else 0.93,
                            "is_cleanlab_issue": given_sp != true_sp,
                            "is_label_corrupted": given_sp != true_sp,
                            "vision_ood_score": 0.22,
                            "vision_is_ood": False,
                            "foreground_attention_ratio": 0.89,
                            "is_artifact_driven": False,
                            "is_background_shortcut": False,
                            "surface_orientation": "adaxial",
                            "mean_luminance": 51.0,
                            "mean_saturation": 0.34,
                            "tomentum_taxon_discordance": False,
                        }
                    )
                    v_rec = v_by_cat[cat]
                    if v_rec["latitude"] is not None and pid == 0:
                        env_rows.append(
                            {
                                "catalogNumber": cat,
                                "latitude": v_rec["latitude"],
                                "longitude": v_rec["longitude"],
                                "fine_soil_ph": 5.4 + (0.05 * idx_o if is_aur else 1.4 + 0.04 * idx_o),
                                "fine_soil_sand": 42.0 + (1.2 * idx_o if is_aur else 18.0 + 0.9 * idx_o),
                                "fine_soil_clay": 24.0 - 0.5 * idx_o,
                                "fine_bedrock_depth_cm": 135.0 if is_aur else 48.0,
                                "ssurgo_rock_outcrop_pct": 2.0 if is_aur else 22.0,
                                "topo_elev_m": 180.0 + 5.0 * idx_o if is_aur else 460.0 + 8.0 * idx_o,
                                "topo_slope_deg": 4.5 if is_aur else 16.0,
                                "topo_tpi_30m": -0.8 if is_aur else 2.4,
                                "topo_twi_30m": 9.8 if is_aur else 5.2,
                                "sgmc_lithology_class": "Felsic_Crystalline" if is_aur else "Carbonate",
                                "musym": "WkB" if is_aur else "GlC",
                                "muname": "Wake-Wateree complex" if is_aur else "Gladeville-Rock outcrop",
                                "outcrop_class": "Granitic flatrock" if is_aur else "Limestone glade",
                                "is_rock_outcrop": False if is_aur else True,
                                "soil_ph": 5.4 + (0.05 * idx_o if is_aur else 1.4 + 0.04 * idx_o),
                                "soil_cec": 11.5 if is_aur else 22.0,
                                "soil_sand": 42.0 + (1.2 * idx_o if is_aur else 18.0 + 0.9 * idx_o),
                                "soil_bulk_density": 1.35 if is_aur else 1.28,
                                "bio1_temp_mean": 14.8 if is_aur else 11.6,
                                "bio4_temp_seasonality": 740.0 + 5.0 * idx_o,
                                "bio12_precip_annual": 1180.0 if is_aur else 980.0,
                                "bio15_precip_seasonality": 14.0 + 0.3 * idx_o,
                                "bio18_precip_warmest_qt": 345.0 if is_aur else 290.0,
                                "bio19_precip_coldest_qt": 265.0 if is_aur else 215.0,
                                "is_training_eligible_coord": True,
                            }
                        )

                upsert_rows_generic(conn, "otu_morphometrics", otu_rows, ("run_id", "catalogNumber", "plant_individual_id"))
                persist_vision_xai_audit(conn, run_id, vis_rows, demo=False)
                upsert_environmental_cache_records(conn, env_rows, demo=False)

                # Export morphometric_flags.csv & label_noise_audit.csv snapshots for Step 06 CLI
                morph_csv = tables_dir / "morphometric_flags.csv"
                audit_csv = tables_dir / "label_noise_audit.csv"
                with open(morph_csv, "w", encoding="utf-8", newline="") as mf:
                    mw = csv.DictWriter(mf, fieldnames=list(otu_rows[0].keys()))
                    mw.writeheader()
                    mw.writerows(otu_rows)
                with open(audit_csv, "w", encoding="utf-8", newline="") as af:
                    aw = csv.DictWriter(af, fieldnames=list(vis_rows[0].keys()))
                    aw.writeheader()
                    aw.writerows(vis_rows)

                # Verify v_integrative_triage_master already exposes all 13 specimen rows prior to Step 06/07
                view_rows_pre = conn.execute(
                    """
                    SELECT catalogNumber, plant_individual_id, cda_predicted_taxon,
                           vision_predicted_label, fine_soil_ph
                    FROM v_integrative_triage_master
                    WHERE run_id = ?
                    ORDER BY catalogNumber, plant_individual_id;
                    """,
                    (run_id,),
                ).fetchall()
                self.assertEqual(len(view_rows_pre), 13)
                pre_by_key = {(str(r[0]), int(r[1])): dict(r) for r in view_rows_pre}
                self.assertEqual(pre_by_key[("MULTI-001", 0)]["cda_predicted_taxon"], "Packera aurea")
                self.assertEqual(pre_by_key[("MULTI-001", 1)]["cda_predicted_taxon"], "Packera obovata")
                self.assertIn(("HIST-NOGEO-99", 0), pre_by_key)
                self.assertEqual(pre_by_key[("HIST-NOGEO-99", 0)]["vision_predicted_label"], "Packera obovata")
                self.assertIsNone(pre_by_key[("HIST-NOGEO-99", 0)]["fine_soil_ph"])
            finally:
                conn.close()

            # Execute Step 06 against SQLite environmental_cache (zero rasters on disk; cached fine-grain values in SQLite!)
            env = os.environ.copy()
            env["PACKERA_DB_PATH"] = str(db_path)
            env["RUN_ID"] = run_id
            env["PACKERA_RUN_DIR"] = str(ws)

            conflict_csv = tables_dir / "multimodal_conflict_flags.csv"
            srf_summary_csv = reports_dir / "multimodal_spatial_rf_summary.csv"
            cmd_step06 = [
                rscript_bin,
                str(PROJECT_ROOT / "scripts" / "analysis" / "06_multimodal_spatial_rf.R"),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(tables_dir / "morphometric_flags.csv"),
                "--vision-audit", str(tables_dir / "label_noise_audit.csv"),
                "--env-dir", str(env_dir),
                "--ssurgo-table", str(tables_dir / "ssurgo_edaphic_validation.csv"),
                "--contingency-table", str(reports_dir / "micro_edaphic_outcrop_contingency.csv"),
                "--pheno-anomalies", str(tables_dir / "phenological_anomalies.csv"),
                "--pheno-summary", str(reports_dir / "phenological_anomaly_summary.csv"),
                "--pheno-plot", str(figures_dir / "phenological_latitudinal_anomaly.pdf"),
                "--output-flags", str(conflict_csv),
                "--output-summary", str(srf_summary_csv),
                "--output-plot", str(figures_dir / "spatial_rf_niche_importance.pdf"),
                "--permutations", "9",
                "--n-trees", "50",
            ]
            proc06 = subprocess.run(cmd_step06, cwd=PROJECT_ROOT, env=env, capture_output=True, text=True, check=False)
            self.assertEqual(proc06.returncode, 0, f"Step 06 failed:\n{proc06.stdout}\n{proc06.stderr}")
            self.assertTrue(conflict_csv.is_file())

            # Execute Step 07 synthesis reading from SQLite v_integrative_triage_master
            # Remove intermediate CSVs (morphometric_flags, label_noise_audit, multimodal_conflict_flags)
            # to prove Step 07 hydrates all 3 upstream evidence streams directly from SQLite!
            (tables_dir / "morphometric_flags.csv").unlink()
            (tables_dir / "label_noise_audit.csv").unlink()
            conflict_csv.unlink()

            triage_queue_csv = tables_dir / "triage_queue.csv"
            final_triage_csv = reports_dir / "final_curation_triage.csv"
            rev_report_md = reports_dir / "Packera_Taxonomic_Revision_Summary.md"
            rev_fig_pdf = figures_dir / "Figure_Integrative_Packera_Revision.pdf"

            cmd_step07 = [
                rscript_bin,
                str(PROJECT_ROOT / "scripts" / "analysis" / "07_triage_dashboard_synthesis.R"),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(tables_dir / "morphometric_flags.csv"),
                "--vision-audit", str(tables_dir / "label_noise_audit.csv"),
                "--multimodal-flags", str(conflict_csv),
                "--niche-summary", str(srf_summary_csv),
                "--output-queue", str(triage_queue_csv),
                "--output-triage", str(final_triage_csv),
                "--output-figure", str(rev_fig_pdf),
                "--output-report", str(rev_report_md),
            ]
            proc07 = subprocess.run(cmd_step07, cwd=PROJECT_ROOT, env=env, capture_output=True, text=True, check=False)
            combined07 = proc07.stdout + "\n" + proc07.stderr
            self.assertEqual(proc07.returncode, 0, f"Step 07 failed:\n{combined07}")
            self.assertIn("v_integrative_triage_master", combined07)
            self.assertTrue(final_triage_csv.is_file())

            with open(final_triage_csv, "r", encoding="utf-8", newline="") as tf:
                triage_rows = list(csv.DictReader(tf))
            self.assertEqual(len(triage_rows), 13)
            by_key = {(r["catalogNumber"], int(float(r["plant_individual_id"]))): r for r in triage_rows}

            # Multi-plant sheet preservation: pid=0 confirmed Packera aurea, pid=1 reassigned to Packera obovata
            self.assertIn(("MULTI-001", 0), by_key)
            self.assertIn(("MULTI-001", 1), by_key)
            self.assertEqual(by_key[("MULTI-001", 0)]["recommended_determination"], "Packera aurea")
            self.assertEqual(by_key[("MULTI-001", 1)]["recommended_determination"], "Packera obovata")
            self.assertEqual(by_key[("MULTI-001", 1)]["taxonomic_treatment_subtype"], "Misidentification_Reassign")

            # Non-georeferenced historic voucher retention & dynamic modality renormalization
            self.assertIn(("HIST-NOGEO-99", 0), by_key)
            nogeo = by_key[("HIST-NOGEO-99", 0)]
            self.assertIn(str(nogeo["is_completed_voucher"]).strip().upper(), {"TRUE", "1"})
            self.assertTrue(str(nogeo["recommended_determination"]).startswith("Packera obovata"))

            # Verify SQLite curation_triage & v_integrative_triage_master hold all 13 rows
            conn_post = get_db_connection(db_path=db_path, read_only=True)
            try:
                n_ct = conn_post.execute(
                    "SELECT COUNT(*) FROM curation_triage WHERE run_id = ?;", (run_id,)
                ).fetchone()[0]
                self.assertEqual(int(n_ct), 13)
                nogeo_db = conn_post.execute(
                    """
                    SELECT recommended_determination, is_completed_voucher
                    FROM v_integrative_triage_master
                    WHERE run_id = ? AND catalogNumber = 'HIST-NOGEO-99' AND plant_individual_id = 0;
                    """,
                    (run_id,),
                ).fetchone()
                self.assertIsNotNone(nogeo_db)
                self.assertTrue(str(nogeo_db[0]).startswith("Packera obovata"))
                self.assertEqual(int(nogeo_db[1]), 1)
            finally:
                conn_post.close()


if __name__ == "__main__":
    unittest.main()
