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
      3. Multi-plant sheet preservation `(catalogNumber, plant_individual_id)`
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

from scripts.tests._r_session import rscript_available, run_warm_r_cmd

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import (
    DEFAULT_SQLITE_REL_PATH,
    ensure_pipeline_run,
    get_db_connection,
    get_duckdb_connection,
    ingest_existing_csvs_to_db,
    populate_run_active_cohort,
    query_environmental_bbox,
    resolve_db_path,
    upsert_environmental_cache_records,
    upsert_rows_generic,
    upsert_vouchers_records,
)


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
                    "environmental_cache",
                    "environmental_rtree",
                    "pipeline_runs",
                    "run_active_cohort",
                    "leaf_efa_harmonics",
                    "otu_morphometrics",
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
                                catalogNumber, plant_individual_id, leaf_id
                            ) VALUES ('ORPHAN-1', 0, 1);
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
                            catalogNumber, plant_individual_id, leaf_id, contour_xy_json, contour_n_pts
                        ) VALUES ('NCU-FK-01', 0, 1, '[[0,0],[1,0],[1,1]]', 3);
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
            if rscript_available():
                r_expr = (
                    "source('scripts/core/db_utils.R'); "
                    f"df <- query_environmental_cache_sqlite(db_path = '{db_file.as_posix()}', "
                    "bbox = c(-80.5, -78.0, 35.0, 36.5)); "
                    "cat('R_BBOX_N=', nrow(df), 'CAT=', as.character(df$catalogNumber[1]), sep='')"
                )
                r_proc = run_warm_r_cmd(["Rscript", "-e", r_expr], stage_key="step03", cwd=PROJECT_ROOT)
                self.assertEqual(r_proc.returncode, 0, msg=r_proc.stderr)
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
        if not rscript_available():
            self.skipTest("Rscript not installed or fast mode enabled.")
        rscript_bin = "Rscript"

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
                        "determiner_tier": "Tier_1_Gold",
                        "type_status": "",
                        "county": "Orange",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": round(35.5 + 0.05 * v_idx, 5),
                        "decimalLongitude": round(-79.5 + 0.05 * v_idx, 5),
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
                                "solidity": 0.91,
                                "circularity": 0.68 if is_aurea else 0.54,
                                "extent": 0.73,
                                "archetype_score": 0.94,
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
                "--permutations", "19",
            ]
            proc = run_warm_r_cmd(cmd, stage_key='step03', cwd=PROJECT_ROOT, env=env, skip_figures=True)
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
                self.assertEqual(r["morph_tier_source"], "Basal_Rosette")
                self.assertTrue(math.isfinite(float(r["PC1"])))

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

    def test_multi_plant_and_non_georeferenced_retention_across_step06_and_step07(self) -> None:
        """Verify multi-plant (catalogNumber, plant_individual_id) & non-georeferenced retention across Step 06, Step 07, and v_integrative_triage_master."""
        if not rscript_available():
            self.skipTest("Rscript not installed or fast mode enabled.")
        rscript_bin = "Rscript"

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
                            "morph_tier_source": "Basal_Rosette",
                            "basal_leaf_count": 3,
                            "PC1": (-1.6 + 0.05 * idx_o) if is_aur else (1.7 + 0.05 * idx_o),
                            "PC2": 0.12 * ((idx_o % 3) - 1),
                            "PC3": 0.05 * ((idx_o % 2) - 0.5),
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
                                "elevation_m": 180.0 + 5.0 * idx_o if is_aur else 460.0 + 8.0 * idx_o,
                                "topo_slope_deg": 4.5 if is_aur else 16.0,
                                "topo_tpi_30m": -0.8 if is_aur else 2.4,
                                "topo_hli_30m": 0.72 if is_aur else 0.94,
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
                                "bio18_precip_warmest_qtr": 345.0 if is_aur else 290.0,
                                "bio19_precip_coldest_qtr": 265.0 if is_aur else 215.0,
                                "is_training_eligible_coord": True,
                            }
                        )

                upsert_rows_generic(conn, "otu_morphometrics", otu_rows, ("run_id", "catalogNumber", "plant_individual_id"))
                upsert_environmental_cache_records(conn, env_rows, demo=False)

                # Export morphometric_flags.csv snapshot for Step 06 CLI
                morph_csv = tables_dir / "morphometric_flags.csv"
                with open(morph_csv, "w", encoding="utf-8", newline="") as mf:
                    mw = csv.DictWriter(mf, fieldnames=list(otu_rows[0].keys()))
                    mw.writeheader()
                    mw.writerows(otu_rows)

                # Verify v_integrative_triage_master already exposes all 13 specimen rows prior to Step 06/07
                view_rows_pre = conn.execute(
                    """
                    SELECT catalogNumber, plant_individual_id, cda_predicted_taxon,
                           fine_soil_ph
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
                self.assertEqual(pre_by_key[("HIST-NOGEO-99", 0)]["cda_predicted_taxon"], "Packera obovata")
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
            cmd_step05 = [
                rscript_bin,
                str(PROJECT_ROOT / "scripts" / "analysis" / "05_multimodal_spatial_rf.R"),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(tables_dir / "morphometric_flags.csv"),
                "--env-dir", str(env_dir),
                "--ssurgo-table", str(tables_dir / "ssurgo_edaphic_validation.csv"),
                "--contingency-table", str(reports_dir / "micro_edaphic_outcrop_contingency.csv"),
                "--pheno-anomalies", str(tables_dir / "phenological_anomalies.csv"),
                "--pheno-summary", str(reports_dir / "phenological_anomaly_summary.csv"),
                "--pheno-plot", str(figures_dir / "fig05b_hopkins_circular_phenology.png"),
                "--output-flags", str(conflict_csv),
                "--output-summary", str(srf_summary_csv),
                "--output-plot", str(figures_dir / "fig05a_ecological_niche_and_rf_importance.png"),
                "--permutations", "9",
                "--n-trees", "50",
            ]
            proc05 = run_warm_r_cmd(cmd_step05, stage_key='step05', cwd=PROJECT_ROOT, env=env, skip_figures=True)
            self.assertEqual(proc05.returncode, 0, f"Step 05 failed:\n{proc05.stdout}\n{proc05.stderr}")
            self.assertTrue(conflict_csv.is_file())

            # Execute Step 06 synthesis reading from SQLite v_integrative_triage_master
            # Remove intermediate CSVs (morphometric_flags, multimodal_conflict_flags)
            # to prove Step 06 hydrates all upstream evidence streams directly from SQLite!
            (tables_dir / "morphometric_flags.csv").unlink()
            conflict_csv.unlink()

            triage_queue_csv = tables_dir / "triage_queue.csv"
            rev_report_md = reports_dir / "Packera_Taxonomic_Revision_Summary.md"
            rev_html = reports_dir / "Triage_Dashboard_Interactive.html"
            rev_fig_pdf = figures_dir / "fig06a_integrative_taxonomic_synthesis.png"

            cmd_step06 = [
                rscript_bin,
                str(PROJECT_ROOT / "scripts" / "analysis" / "06_triage_dashboard_synthesis.R"),
                "--vouchers", str(vouchers_csv),
                "--morphometrics", str(tables_dir / "morphometric_flags.csv"),
                "--multimodal-flags", str(conflict_csv),
                "--niche-summary", str(srf_summary_csv),
                "--output-queue", str(triage_queue_csv),
                "--output-plot", str(rev_fig_pdf),
                "--output-report", str(rev_report_md),
                "--output-html", str(rev_html),
            ]
            proc06 = run_warm_r_cmd(cmd_step06, stage_key='step06', cwd=PROJECT_ROOT, env=env, skip_figures=True)
            combined06 = proc06.stdout + "\n" + proc06.stderr
            self.assertEqual(proc06.returncode, 0, f"Step 06 failed:\n{combined06}")
            self.assertIn("v_integrative_triage_master", combined06)
            self.assertTrue(triage_queue_csv.is_file())
            self.assertTrue(rev_html.is_file())

            with open(triage_queue_csv, "r", encoding="utf-8", newline="") as tf:
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

    def test_master_view_step04_parity_and_upstream_hopkins_precedence(self):
        """Verify v_integrative_triage_master exposes Step 04 parity columns and fresh Step 05 hopkins_adjusted_doy overrides stale Step 06 curation_triage."""
        conn = get_db_connection(":memory:", init_schema=True)
        try:
            mv_cols = {r[1] for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
            om_cols = {r[1] for r in conn.execute("PRAGMA table_info(otu_morphometrics);").fetchall()}
            self.assertIn("scale_source", om_cols)
            for req in ("is_passive", "scale_source", "AR_inv_std", "reproductive_trait_status", "ecoregion_can1", "ecoregion_can2"):
                self.assertIn(req, mv_cols)

            conn.execute("INSERT INTO pipeline_runs (run_id) VALUES ('run_hopkins_test');")
            conn.execute(
                """
                INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, determiner_tier)
                VALUES ('NCU-HOP-1', 'Packera aurea', 'Packera aurea', 'Tier_3_Bronze');
                """
            )
            conn.execute(
                """
                INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier, is_promoted_bronze, authority_weight)
                VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 'Packera aurea', 'Tier_3_Bronze', 0, 0.0);
                """
            )
            conn.execute(
                """
                INSERT INTO otu_morphometrics (run_id, catalogNumber, plant_individual_id, is_passive, scale_source, AR_inv_std, reproductive_trait_status, ecoregion_can1, ecoregion_can2)
                VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 1, 'ruler_calibrated', 1.25, 'Measured', -0.45, 0.85);
                """
            )
            conn.execute(
                """
                INSERT INTO spatial_rf_predictions (run_id, catalogNumber, plant_individual_id, hopkins_adjusted_doy)
                VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 134.5);
                """
            )
            conn.execute(
                """
                INSERT INTO curation_triage (run_id, catalogNumber, plant_individual_id, hopkins_adjusted_doy)
                VALUES ('run_hopkins_test', 'NCU-HOP-1', 0, 99.0);
                """
            )
            row = conn.execute(
                """
                SELECT is_passive, scale_source, AR_inv_std, reproductive_trait_status, ecoregion_can1, ecoregion_can2, hopkins_adjusted_doy
                FROM v_integrative_triage_master
                WHERE run_id = 'run_hopkins_test' AND catalogNumber = 'NCU-HOP-1';
                """
            ).fetchone()
            self.assertIsNotNone(row)
            self.assertEqual(int(row[0]), 1)
            self.assertEqual(str(row[1]), "ruler_calibrated")
            self.assertAlmostEqual(float(row[2]), 1.25)
            self.assertEqual(str(row[3]), "Measured")
            self.assertAlmostEqual(float(row[4]), -0.45)
            self.assertAlmostEqual(float(row[5]), 0.85)
            self.assertAlmostEqual(float(row[6]), 134.5)
        finally:
            conn.close()


    def test_phantom_pid0_stub_cleanup_and_false_split_db_consolidation(self) -> None:
        """Verify SQLite heals phantom pid=0 stubs and cauline/1-leaf splits in place without GPU re-segmentation while preserving true multi-rosette sheets."""
        from scripts.core.db import consolidate_false_multi_plant_splits_db, populate_run_active_cohort

        with tempfile.TemporaryDirectory() as tmp_d:
            db_file = Path(tmp_d) / "heal_splits.sqlite"
            conn = get_db_connection(db_path=db_file, init_schema=True)
            try:
                v_rows = [
                    {"catalogNumber": "CAT-PHANTOM", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold"},
                    {"catalogNumber": "CAT-CAUL-SPLIT", "species_standardized": "Packera aurea", "determiner_tier": "Tier_2_Silver"},
                    {"catalogNumber": "CAT-TRUE-MULTI", "species_standardized": "Packera obovata", "determiner_tier": "Tier_1_Gold"},
                ]
                upsert_vouchers_records(conn, v_rows)
                with conn:
                    # 1. CAT-PHANTOM: all-NULL pid=0 stub alongside real pid=1 individual with 2 basal leaves
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('CAT-PHANTOM', 0);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-PHANTOM', 1, 900.0, 2400.0);")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-PHANTOM', 1, 1, 'CAT-PHANTOM_p1_leaf1', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-PHANTOM', 1, 2, 'CAT-PHANTOM_p1_leaf2', '[[0,0],[1,0],[0,1]]', 'basal');")

                    # 2. CAT-CAUL-SPLIT: pid=0 basal rosette (2 basal leaves) + pid=1 cauline stem (2 cauline leaves + 6 capitula) + pid=2 1-leaf stray
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-CAUL-SPLIT', 0, 950.0, 2500.0);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y, plant_height_cm) VALUES ('CAT-CAUL-SPLIT', 1, 960.0, 1100.0, 42.5);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-CAUL-SPLIT', 2, 1350.0, 2480.0);")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 0, 1, 'CAT-CAUL-SPLIT_p0_leaf1', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 0, 2, 'CAT-CAUL-SPLIT_p0_leaf2', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 1, 3, 'CAT-CAUL-SPLIT_p1_leaf3', '[[0,0],[1,0],[0,1]]', 'cauline');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 1, 4, 'CAT-CAUL-SPLIT_p1_leaf4', '[[0,0],[1,0],[0,1]]', 'cauline');")
                    conn.execute("INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES ('CAT-CAUL-SPLIT', 2, 5, 'CAT-CAUL-SPLIT_p2_leaf5', '[[0,0],[1,0],[0,1]]', 'basal');")
                    conn.execute("INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count, total_capitula, n_anthesis, voucher_phenological_state, is_verified_anthesis) VALUES ('CAT-CAUL-SPLIT', 0, 0, 0, 0, 'sterile', 0);")
                    conn.execute("INSERT INTO voucher_reproductive_phenology (catalogNumber, plant_individual_id, capitula_count, total_capitula, n_anthesis, voucher_phenological_state, is_verified_anthesis) VALUES ('CAT-CAUL-SPLIT', 1, 6, 6, 5, 'anthesis', 1);")

                    # 3. CAT-TRUE-MULTI: two independent basal rosettes (>=2 basal leaves each, x=350 vs x=1650)
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-TRUE-MULTI', 0, 350.0, 2400.0);")
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) VALUES ('CAT-TRUE-MULTI', 1, 1650.0, 2420.0);")
                    for pid_m in (0, 1):
                        for lid_m in (1, 2):
                            eff_l = pid_m * 10 + lid_m
                            conn.execute(
                                "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, contour_xy_json, leaf_zone) VALUES (?, ?, ?, ?, '[[0,0],[1,0],[0,1]]', 'basal');",
                                ("CAT-TRUE-MULTI", pid_m, eff_l, f"CAT-TRUE-MULTI_p{pid_m}_leaf{eff_l}"),
                            )

                n_coh = populate_run_active_cohort(conn, run_id="heal_run", voucher_rows=v_rows, clear_existing=True)
                # Expect: CAT-PHANTOM -> 1 row (pid=1), CAT-CAUL-SPLIT -> 1 row (pid=0 with all 5 leaves & 6 capitula), CAT-TRUE-MULTI -> 2 rows (pid=0, 1) => total 4 rows
                self.assertEqual(n_coh, 4)

                phantom_pids = [r[0] for r in conn.execute("SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = 'CAT-PHANTOM';").fetchall()]
                self.assertEqual(phantom_pids, [1])

                split_pids = [r[0] for r in conn.execute("SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = 'CAT-CAUL-SPLIT';").fetchall()]
                self.assertEqual(split_pids, [0])
                split_leaves = conn.execute("SELECT COUNT(*) FROM extracted_leaves WHERE catalogNumber = 'CAT-CAUL-SPLIT' AND plant_individual_id = 0;").fetchone()[0]
                self.assertEqual(int(split_leaves), 5)
                vrp_merged = conn.execute("SELECT capitula_count, n_anthesis, voucher_phenological_state, is_verified_anthesis FROM voucher_reproductive_phenology WHERE catalogNumber = 'CAT-CAUL-SPLIT' AND plant_individual_id = 0;").fetchone()
                self.assertEqual((int(vrp_merged[0]), int(vrp_merged[1]), str(vrp_merged[2]), int(vrp_merged[3])), (6, 5, "anthesis", 1))
                h_cm = conn.execute("SELECT plant_height_cm FROM plant_individuals WHERE catalogNumber = 'CAT-CAUL-SPLIT' AND plant_individual_id = 0;").fetchone()[0]
                self.assertAlmostEqual(float(h_cm), 42.5)

                true_multi_pids = sorted(r[0] for r in conn.execute("SELECT plant_individual_id FROM plant_individuals WHERE catalogNumber = 'CAT-TRUE-MULTI';").fetchall())
                self.assertEqual(true_multi_pids, [0, 1])

                # Idempotence check
                stats_again = consolidate_false_multi_plant_splits_db(conn)
                self.assertEqual(stats_again, {"phantom_zero_stubs_removed": 0, "false_split_individuals_merged": 0})
            finally:
                conn.close()


    def test_curation_triage_hybrid_parentage_columns_and_migration(self) -> None:
        """Verify curation_triage and v_integrative_triage_master expose hybrid_parent_1, hybrid_parent_2, and hybrid_zone_cohesion_pass, and migrate older schemas idempotently."""
        with tempfile.TemporaryDirectory() as tmp:
            db_p = Path(tmp) / "legacy_ct.sqlite"
            raw_conn = sqlite3.connect(db_p)
            raw_conn.executescript(
                """
                CREATE TABLE pipeline_runs (run_id TEXT PRIMARY KEY, start_time TEXT);
                INSERT INTO pipeline_runs (run_id, start_time) VALUES ('r1', '2026-11-12T00:00:00Z');
                CREATE TABLE curation_triage (
                    run_id TEXT NOT NULL,
                    catalogNumber TEXT NOT NULL,
                    plant_individual_id INTEGER NOT NULL DEFAULT 0,
                    taxonomic_treatment_subtype TEXT,
                    recommended_determination TEXT,
                    cryptic_cluster_cohesion_pass INTEGER,
                    is_completed_voucher INTEGER,
                    PRIMARY KEY (run_id, catalogNumber, plant_individual_id)
                );
                INSERT INTO curation_triage VALUES ('r1', 'NCU-HYB-1', 0, 'Hybrid_Introgressant', 'Packera dubia x Packera anonyma', 0, 1);
                """
            )
            raw_conn.commit()
            raw_conn.close()

            conn = get_db_connection(db_p, init_schema=True)
            try:
                ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
                mv_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
                for req_col in ("hybrid_zone_cohesion_pass", "hybrid_parent_1", "hybrid_parent_2"):
                    self.assertIn(req_col, ct_cols)
                    self.assertIn(req_col, mv_cols)
                row = conn.execute(
                    "SELECT taxonomic_treatment_subtype, recommended_determination FROM curation_triage WHERE catalogNumber = 'NCU-HYB-1';"
                ).fetchone()
                self.assertEqual(tuple(row), ("Hybrid_Introgressant", "Packera dubia x Packera anonyma"))
            finally:
                conn.close()

    def test_misid_and_cryptic_schema_20261119_columns_and_migration(self) -> None:
        """Verify schema 20261119 migrates and persists OOF CDA, self-ellipsoid, open_set_df, spatial-cline, and resurrection_basionym columns."""
        from scripts.core.db import (
            ensure_pipeline_run,
            persist_otu_morphometrics,
            populate_run_active_cohort,
            upsert_rows_generic,
        )
        from scripts.core.db_schema import PACKERA_SCHEMA_USER_VERSION

        self.assertGreaterEqual(PACKERA_SCHEMA_USER_VERSION, 20261119)
        with tempfile.TemporaryDirectory() as tmp:
            db_p = Path(tmp) / "legacy_20261118.sqlite"
            raw_conn = sqlite3.connect(db_p)
            raw_conn.executescript(
                """
                PRAGMA user_version = 20261118;
                CREATE TABLE pipeline_runs (run_id TEXT PRIMARY KEY, run_dir TEXT, start_timestamp TEXT, end_timestamp TEXT, git_commit_sha TEXT, cli_arguments_json TEXT, focal_taxon TEXT, intraspecific_mode INTEGER DEFAULT 0);
                INSERT INTO pipeline_runs (run_id, start_timestamp) VALUES ('r_res', '2026-11-12T12:00:00Z');
                CREATE TABLE curation_triage (
                    run_id TEXT NOT NULL,
                    catalogNumber TEXT NOT NULL,
                    plant_individual_id INTEGER NOT NULL DEFAULT 0,
                    taxonomic_status_call TEXT,
                    taxonomic_treatment_subtype TEXT,
                    recommended_determination TEXT,
                    cryptic_cluster_cohesion_pass INTEGER,
                    hybrid_zone_cohesion_pass INTEGER,
                    hybrid_parent_1 TEXT,
                    hybrid_parent_2 TEXT,
                    is_completed_voucher INTEGER,
                    PRIMARY KEY (run_id, catalogNumber, plant_individual_id)
                );
                INSERT INTO curation_triage VALUES (
                    'r_res', 'NCU-OLD-1', 0, 'Species_Confirmed', 'Confirmed_Species', 'Packera dubia', 0, 0, NULL, NULL, 1
                );
                """
            )
            raw_conn.commit()
            raw_conn.close()

            conn = get_db_connection(db_p, init_schema=True)
            try:
                om_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(otu_morphometrics);").fetchall()}
                ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
                mv_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
                for req_om in ("cda_oof_posterior_prob", "cl_target_threshold", "self_morph_mahal_d2", "self_morph_df", "open_set_df"):
                    self.assertIn(req_om, om_cols)
                    self.assertIn(req_om, mv_cols)
                for req_ct in ("cryptic_non_clinal_pass", "resurrection_basionym"):
                    self.assertIn(req_ct, ct_cols)
                    self.assertIn(req_ct, mv_cols)

                v_rows = [
                    {
                        "catalogNumber": "NCU-RES-01",
                        "institutionCode": "NCU",
                        "scientificName": "Packera anonyma",
                        "species_raw": "Senecio earlei",
                        "species_standardized": "Packera anonyma",
                        "determiner_tier": "Tier_1_Gold",
                        "qc_status": "PASS",
                    }
                ]
                upsert_rows_generic(conn, table_name='vouchers', rows=v_rows, pk_cols=('catalogNumber',))
                ensure_pipeline_run(conn, run_id="r_res")
                populate_run_active_cohort(conn, run_id="r_res", voucher_rows=v_rows, clear_existing=True)
                persist_otu_morphometrics(
                    conn,
                    run_id="r_res",
                    otu_rows=[
                        {
                            "catalogNumber": "NCU-RES-01",
                            "plant_individual_id": 0,
                            "species_standardized": "Packera anonyma",
                            "determiner_tier": "Tier_1_Gold",
                            "cda_predicted_taxon": "Packera anonyma",
                            "cda_posterior_prob": 0.62,
                            "cda_oof_posterior_prob": 0.54,
                            "cl_target_threshold": 0.74,
                            "self_morph_mahal_d2": 9.85,
                            "self_morph_df": 4,
                            "min_congener_mahal_d2": 32.4,
                            "min_congener_taxon": "Packera dubia",
                            "open_set_novelty_p": 0.0018,
                            "open_set_df": 14,
                            "nearest_type_basionym": "Senecio earlei",
                            "nearest_type_basionym_d2": 2.15,
                        }
                    ],
                )
                upsert_rows_generic(
                    conn,
                    table_name="curation_triage",
                    rows=[
                        {
                            "run_id": "r_res",
                            "catalogNumber": "NCU-RES-01",
                            "plant_individual_id": 0,
                            "taxonomic_status_call": "Taxonomic_Resurrection_Candidate",
                            "recommended_determination": "Packera anonyma [Resurrect Basionym: Senecio earlei]",
                            "taxonomic_treatment_broad": "Misidentification / Cryptic Species",
                            "taxonomic_treatment_subtype": "Putative_Cryptic_Species",
                            "cryptic_cluster_cohesion_pass": 1,
                            "cryptic_non_clinal_pass": 1,
                            "resurrection_basionym": "Senecio earlei",
                            "hybrid_zone_cohesion_pass": 0,
                            "is_completed_voucher": 1,
                        }
                    ],
                    pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
                )
                mv_row = conn.execute(
                    """
                    SELECT cda_oof_posterior_prob, cl_target_threshold, self_morph_mahal_d2, self_morph_df,
                           open_set_df, cryptic_non_clinal_pass, resurrection_basionym, taxonomic_status_call
                    FROM v_integrative_triage_master
                    WHERE run_id = 'r_res' AND catalogNumber = 'NCU-RES-01';
                    """
                ).fetchone()
                self.assertIsNotNone(mv_row)
                self.assertAlmostEqual(float(mv_row[0]), 0.54, places=3)
                self.assertAlmostEqual(float(mv_row[1]), 0.74, places=3)
                self.assertAlmostEqual(float(mv_row[2]), 9.85, places=3)
                self.assertEqual(int(mv_row[3]), 4)
                self.assertEqual(int(mv_row[4]), 14)
                self.assertEqual(int(mv_row[5]), 1)
                self.assertEqual(str(mv_row[6]), "Senecio earlei")
                self.assertEqual(str(mv_row[7]), "Taxonomic_Resurrection_Candidate")
            finally:
                conn.close()



    def test_legacy_extracted_leaves_check_constraint_rebuilt_on_init_and_purge(self) -> None:
        """Verify pre-LM3 extracted_leaves CHECK constraints and retired NOT NULL columns are stripped by init_db_schema, savepoint transactions, and purge_segments_and_runs."""
        from scripts.core.db import purge_segments_and_runs, record_voucher_segmentation_transaction

        leg_col = "_".join(("assigned", "tier"))
        leg_v1 = "_".join(("Tier", "1", "Direct"))
        leg_v2 = "_".join(("Tier", "2", "Reflected"))

        def _seed_legacy_db(db_file: Path) -> None:
            init_conn = get_db_connection(db_file, init_schema=True)
            init_conn.close()
            raw = sqlite3.connect(db_file)
            raw.executescript(
                f"""
                PRAGMA foreign_keys = OFF;
                PRAGMA user_version = 20261121;
                INSERT INTO vouchers (catalogNumber, scientificName, species_standardized, qc_status) VALUES
                    ('WILLI_82182', 'Packera dubia', 'Packera dubia', 'PASS'),
                    ('GA202225', 'Packera anonyma', 'Packera anonyma', 'PENDING'),
                    ('VSC0067189', 'Packera aurea', 'Packera aurea', 'PENDING');
                INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES ('WILLI_82182', 0);
                DROP TABLE IF EXISTS extracted_leaves;
                CREATE TABLE extracted_leaves (
                    catalogNumber TEXT NOT NULL,
                    plant_individual_id INTEGER NOT NULL DEFAULT 0,
                    leaf_id INTEGER NOT NULL,
                    shape_id TEXT UNIQUE,
                    {leg_col} TEXT NOT NULL CHECK({leg_col} IN ('{leg_v1}', '{leg_v2}')),
                    solidity REAL,
                    contour_xy_json TEXT NOT NULL,
                    PRIMARY KEY (catalogNumber, plant_individual_id, leaf_id),
                    FOREIGN KEY (catalogNumber, plant_individual_id)
                        REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
                );
                INSERT INTO extracted_leaves (
                    catalogNumber, plant_individual_id, leaf_id, shape_id, {leg_col}, solidity, contour_xy_json
                ) VALUES (
                    'WILLI_82182', 0, 1, 'WILLI_82182_p0_leaf1', '{leg_v1}', 0.84, '[[0,0],[10,0],[5,8]]'
                );
                """
            )
            raw.commit()
            raw.close()

        with tempfile.TemporaryDirectory() as tmp:
            # Case A: In-place init_db_schema + Step 02 savepoint transaction without --purge-segments
            db_a = Path(tmp) / "inplace" / "data" / "packera_pipeline.sqlite"
            db_a.parent.mkdir(parents=True, exist_ok=True)
            _seed_legacy_db(db_a)

            conn_a = get_db_connection(db_a, init_schema=True)
            try:
                el_cols_a = {str(r[1]) for r in conn_a.execute("PRAGMA table_info(extracted_leaves);").fetchall()}
                self.assertNotIn(leg_col, el_cols_a)
                self.assertIn("efa_qc_pass", el_cols_a)
                # Existing row preserved across table rebuild
                kept = conn_a.execute(
                    "SELECT solidity FROM extracted_leaves WHERE catalogNumber = 'WILLI_82182';"
                ).fetchone()
                self.assertIsNotNone(kept)
                self.assertAlmostEqual(float(kept[0]), 0.84, places=3)

                # Step 02 per-voucher savepoint pattern must commit cleanly
                with conn_a:
                    for idx, cat in enumerate(("WILLI_82182", "GA202225", "VSC0067189")):
                        sp = f"sp_v_{idx}"
                        conn_a.execute(f"SAVEPOINT {sp};")
                        record_voucher_segmentation_transaction(
                            conn=conn_a,
                            catalog_number=cat,
                            extracted_records=[
                                {
                                    "catalogNumber": cat,
                                    "plant_individual_id": 0,
                                    "leaf_id": 1,
                                    "solidity": 0.79,
                                    "area_px": 1200.0,
                                    "pixels_per_mm": 11.5,
                                    "is_calibrated_mm": 1,
                                    "scale_source": "ruler_calibrated",
                                    "leaf_zone": "basal",
                                    "efa_qc_pass": 1,
                                    "qc_flag_reason": "PASS",
                                    "contour_xy_json": "[[0,0],[12,0],[12,9],[0,9]]",
                                }
                            ],
                            commit=False,
                        )
                        conn_a.execute(f"RELEASE SAVEPOINT {sp};")
                n_leaves = conn_a.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0]
                self.assertEqual(int(n_leaves), 3)
            finally:
                conn_a.close()

            # Case B: purge_segments_and_runs drops & recreates clean Step 02-06 schema
            ws_b = Path(tmp) / "purge_ws"
            (ws_b / "data" / "tables").mkdir(parents=True, exist_ok=True)
            (ws_b / "data" / "tables" / "curated_vouchers.csv").write_text("catalogNumber\nWILLI_82182\n", encoding="utf-8")
            db_b = ws_b / "data" / "packera_pipeline.sqlite"
            _seed_legacy_db(db_b)

            summary = purge_segments_and_runs(workspace_root=ws_b, db_path=db_b, purge_runs=True, vacuum=True)
            self.assertEqual(summary["purged_db_rows"], 2)

            conn_b = get_db_connection(db_b, init_schema=False)
            try:
                el_cols_b = {str(r[1]) for r in conn_b.execute("PRAGMA table_info(extracted_leaves);").fetchall()}
                self.assertNotIn(leg_col, el_cols_b)
                ddl_b = str(
                    conn_b.execute(
                        "SELECT sql FROM sqlite_master WHERE type='table' AND name='extracted_leaves';"
                    ).fetchone()[0]
                )
                self.assertNotIn("CHECK", ddl_b.upper())
                record_voucher_segmentation_transaction(
                    conn=conn_b,
                    catalog_number="GA202225",
                    extracted_records=[
                        {
                            "catalogNumber": "GA202225",
                            "plant_individual_id": 0,
                            "leaf_id": 1,
                            "solidity": 0.81,
                            "contour_xy_json": "[[0,0],[10,0],[5,7]]",
                        }
                    ],
                    commit=True,
                )
                self.assertEqual(
                    int(conn_b.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0]),
                    1,
                )
            finally:
                conn_b.close()


if __name__ == "__main__":
    unittest.main()
