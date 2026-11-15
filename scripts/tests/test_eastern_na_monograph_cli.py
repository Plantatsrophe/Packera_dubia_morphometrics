#!/usr/bin/env python3
"""Pure-Python unit tests for nested run IDs and Eastern NA monograph orchestration."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import List, Mapping, Sequence
from unittest import mock

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

TRIAGE_CSV_HEADERS = [
    "catalogNumber", "plant_individual_id", "institutionCode", "species_standardized", "recommended_determination",
    "taxonomic_treatment_subtype", "synthesis_triage_priority", "gmm_cluster", "cryptic_cluster_cohesion_pass", "delta_doy_resid", "synthesis_rationale",
]


def _write_csv(path: Path, headers: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=headers)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in headers})


def _seed_voucher_db(db_path: Path, species_counts: Mapping[str, int]) -> Path:
    from scripts.core.db import get_db_connection

    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = get_db_connection(db_path=db_path, init_schema=True)
    try:
        with conn:
            for taxon, count in species_counts.items():
                ep = re.sub(r"[^A-Za-z0-9]+", "", taxon.split()[-1]).upper()[:8]
                for idx in range(1, count + 1):
                    cat = f"NCU_{ep}_{idx:03d}"
                    conn.execute(
                        "INSERT INTO vouchers (catalogNumber, institutionCode, scientificName, species_raw, species_standardized, determiner_tier, county, stateProvince, decimalLatitude, decimalLongitude, regional_group, pixels_per_mm, is_calibrated_mm, qc_status) VALUES (?, 'NCU', ?, ?, ?, 'Tier_1_Gold', 'Orange', 'NC', 35.9, -79.0, 'Eastern_NA', 24.0, 1, 'PASS');",
                        (cat, taxon, taxon, taxon),
                    )
                    conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 1);", (cat,))
                    conn.execute(
                        "INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, shape_id, leaf_zone, lm3_gates_pass, efa_qc_pass, contour_xy_json) VALUES (?, 1, 1, ?, 'basal', 1, 1, '[[0,0],[10,0],[10,10],[0,10]]');",
                        (cat, f"{cat}_p1_l1"),
                    )
    finally:
        conn.close()
    return db_path


def _seed_stage_outputs(
    ws: Path,
    db_path: Path,
    run_id: str,
    triage_rows: Sequence[Mapping[str, object]],
    *,
    neighbors: Sequence[str] = (),
) -> None:
    from scripts.core.db import ensure_pipeline_run, get_db_connection

    run_dir = ws / "runs" / run_id
    tdir, rdir = run_dir / "data" / "tables", run_dir / "outputs" / "reports"
    _write_csv(tdir / "triage_queue.csv", TRIAGE_CSV_HEADERS, triage_rows)
    _write_csv(rdir / "gmm_bayes_factors_summary.csv", ["optimal_K", "two_delta_BIC"], [{"optimal_K": 2, "two_delta_BIC": 14.5}])
    _write_csv(rdir / "multimodal_spatial_rf_summary.csv", ["D_edaph"], [{"D_edaph": 0.24}])
    if neighbors:
        (rdir / "focal_auto_neighbors.txt").write_text("\n".join(neighbors) + "\n", encoding="utf-8")

    conn = get_db_connection(db_path=db_path, init_schema=True)
    try:
        ensure_pipeline_run(conn, run_id=run_id, run_dir=str(run_dir))
        with conn:
            for idx, row in enumerate(triage_rows, start=1):
                doy_val = row.get("delta_doy_resid")
                cohesive = 1 if str(row.get("cryptic_cluster_cohesion_pass", "")).upper() == "TRUE" else 0
                conn.execute(
                    "INSERT OR REPLACE INTO curation_triage (run_id, catalogNumber, plant_individual_id, recommended_determination, taxonomic_treatment_subtype, synthesis_triage_priority, cryptic_cluster_cohesion_pass, delta_doy_resid, synthesis_rationale) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);",
                    (
                        run_id,
                        str(row.get("catalogNumber", f"ROW_{idx}")),
                        int(row.get("plant_individual_id", 1)),
                        str(row.get("recommended_determination") or row.get("species_standardized", "Packera dubia")),
                        str(row.get("taxonomic_treatment_subtype", "Confirmed_Species")),
                        str(row.get("synthesis_triage_priority", "RESOLVED")),
                        cohesive,
                        float(doy_val) if doy_val not in ("", None) else None,
                        str(row.get("synthesis_rationale", "")),
                    ),
                )
            seen_sp = sorted({str(r.get("species_standardized", "")).strip() for r in triage_rows if r.get("species_standardized")})
            for c_idx, sp in enumerate(seen_sp, start=1):
                conn.execute(
                    "INSERT OR REPLACE INTO taxon_complexes (run_id, species_standardized, n_vouchers, consensus_complex_id, complex_label, nearest_neighbor_taxon, within_complex_mean_affinity, integrative_congruence) VALUES (?, ?, 16, ?, ?, ?, 0.81, 1);",
                    (run_id, sp, c_idx, f"Complex_{c_idx}", neighbors[0] if neighbors else ""),
                )
    finally:
        conn.close()


def _extract_pipeline_run_ids(mock_run: mock.MagicMock, ws: Path) -> List[str]:
    from scripts.core.db import normalize_run_id_token

    run_ids: List[str] = []
    for call in mock_run.call_args_list:
        if not call.args or not isinstance(call.args[0], (list, tuple)):
            continue
        parts = [str(x) for x in call.args[0]]
        if any("run_pipeline.sh" in p for p in parts) and "--run-id" in parts:
            idx = parts.index("--run-id")
            if idx + 1 < len(parts):
                run_ids.append(normalize_run_id_token(parts[idx + 1], workspace_root=ws))
    return run_ids


class TestEasternNAMonographCLI(unittest.TestCase):
    """Pure-Python regression suite for Eastern NA Packera monograph orchestration."""

    def test_parent_run_id_format_uniqueness_and_launcher(self) -> None:
        """Verify generate_parent_run_id uniqueness, build_parser flags, and bash launcher."""
        from scripts.pipeline.resolve_eastern_na_packera import build_parser, generate_parent_run_id

        seen = {generate_parent_run_id("eastern_na") for _ in range(10)}
        self.assertEqual(len(seen), 10)
        for rid in seen:
            self.assertRegex(rid, r"^eastern_na_\d{8}_[0-9a-f]{6}$")
        self.assertTrue((PROJECT_ROOT / "resolve_eastern_na_packera.sh").is_file())
        args = build_parser().parse_args(
            ["--dry-run", "--skip-harvest", "-j", "16", "--parent-run-id", "eastern_na_20261114_abcdef"]
        )
        self.assertTrue(args.dry_run and args.skip_harvest)
        self.assertEqual(args.parent_run_id, "eastern_na_20261114_abcdef")

    def test_nested_run_id_resolution_and_sqlite_triage_isolation(self) -> None:
        """Verify nested run ID normalization and strict SQLite curation_triage isolation."""
        from scripts.core.config import resolve_canonical_run_id
        from scripts.core.db import (
            ensure_pipeline_run,
            get_db_connection,
            normalize_run_id_token,
            query_misidentified_vouchers_from_db,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp).resolve()
            flat_dir = ws / "runs" / "dubia_03_intraspecific_clean"
            nested_id = "eastern_na_20261114_a7f3c2/Packera_dubia/03_intraspecific_clean"
            nested_dir = ws / "runs" / nested_id
            flat_dir.mkdir(parents=True, exist_ok=True)
            nested_dir.mkdir(parents=True, exist_ok=True)

            self.assertEqual(resolve_canonical_run_id(run_dir=flat_dir, explicit_run_id=None, workspace_root=ws), "dubia_03_intraspecific_clean")
            self.assertEqual(resolve_canonical_run_id(run_dir=nested_dir, explicit_run_id=None, workspace_root=ws), nested_id)
            self.assertEqual(resolve_canonical_run_id(run_dir=nested_dir, explicit_run_id=f"runs/{nested_id}", workspace_root=ws), nested_id)

            triage_id = "eastern_na_20261114_a7f3c2/Packera_dubia/02_complex_triage"
            triage_dir = ws / "runs" / triage_id
            for tok in (
                triage_id, f"runs/{triage_id}", f"{triage_id}/data/tables/triage_queue.csv",
                f"runs/{triage_id}/data/tables/triage_queue.csv", str(triage_dir), str(triage_dir / "data" / "tables" / "triage_queue.csv"),
            ):
                self.assertEqual(normalize_run_id_token(tok, workspace_root=ws), triage_id)

            conn = get_db_connection(db_path=ws / "data" / "plantmatics.sqlite", init_schema=True)
            try:
                dub_run = "eastern_na_20261114_a7f3c2/Packera_dubia/02_complex_triage"
                aur_run = "eastern_na_20261114_a7f3c2/Packera_aurea/02_complex_triage"
                ensure_pipeline_run(conn, run_id=dub_run, run_dir=str(ws / "runs" / dub_run))
                ensure_pipeline_run(conn, run_id=aur_run, run_dir=str(ws / "runs" / aur_run))
                with conn:
                    for cat, sp, rid in (("V_DUB_MISID", "Packera dubia", dub_run), ("V_AUR_MISID", "Packera aurea", aur_run)):
                        conn.execute("INSERT INTO vouchers (catalogNumber, species_standardized, qc_status) VALUES (?, ?, 'PASS');", (cat, sp))
                        conn.execute("INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 0);", (cat,))
                        conn.execute(
                            "INSERT INTO curation_triage (run_id, catalogNumber, plant_individual_id, taxonomic_treatment_broad, taxonomic_treatment_subtype) VALUES (?, ?, 0, 'Misidentified_Specimen', 'Misidentification_Reassign');",
                            (rid, cat),
                        )
                flagged, n_runs = query_misidentified_vouchers_from_db(conn, triage_sources=[dub_run], workspace_root=ws)
                self.assertEqual(n_runs, 1)
                self.assertEqual(flagged, {"V_DUB_MISID"})
                self.assertNotIn("V_AUR_MISID", flagged)
            finally:
                conn.close()

    def test_dry_run_dynamic_discovery_and_zero_disk_writes(self) -> None:
        """Verify --dry-run discovers >=15-OTU Eastern NA species cleanly without writing runs/."""
        from scripts.pipeline.resolve_eastern_na_packera import main as monograph_main
        from scripts.pipeline.resolve_focal_species_complex import _stage_run_ids, build_four_stage_commands, slugify_taxon

        self.assertTrue(callable(_stage_run_ids))
        self.assertEqual(slugify_taxon("Packera anonyma"), "Packera_anonyma")
        stages = build_four_stage_commands(run_script=PROJECT_ROOT / "run_pipeline.sh", focal_taxon="Packera dubia", prefix="dubia", geographic_scope="eastern_na")
        self.assertEqual(len(stages), 4)
        self.assertTrue(all(["--geographic-scope", "eastern_na"] in [c[i:i + 2] for i in range(len(c) - 1)] for _, _, c in stages))
        self.assertIn(["--triage-source", "dubia_01_eastern_discovery"], [stages[1][2][i:i + 2] for i in range(len(stages[1][2]) - 1)])

        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp).resolve()
            (ws / "run_pipeline.sh").write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            db_path = _seed_voucher_db(ws / "data" / "plantmatics.sqlite", {"Packera anonyma": 16, "Packera aurea": 16, "Packera obovata": 15, "Packera schweinitziana": 4})

            def _forbid_pipeline(*args, **kwargs):
                cmd = args[0] if args else []
                if isinstance(cmd, (list, tuple)) and any("run_pipeline.sh" in str(x) for x in cmd):
                    raise AssertionError(f"--dry-run must not execute run_pipeline.sh: {cmd}")
                return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="", stderr="")

            buf = io.StringIO()
            with mock.patch("subprocess.run", side_effect=_forbid_pipeline), contextlib.redirect_stdout(buf):
                rc = monograph_main(["--dry-run", "--skip-harvest", "-j", "16", "--db", str(db_path)], workspace_root=ws)
            self.assertEqual(rc, 0)
            self.assertFalse((ws / "runs").exists(), "--dry-run must not create workspace runs/")

            stdout = buf.getvalue()
            cmd_lines = [ln for ln in stdout.splitlines() if "run_pipeline.sh" in ln and "--run-id" in ln]
            self.assertEqual(len(cmd_lines), 12, f"Expected 3 genus + 9 species commands:\n{stdout}")
            self.assertTrue(all("--geographic-scope eastern_na" in ln for ln in cmd_lines))
            for stage in ("00_genus_wide/01_raw_discovery", "00_genus_wide/02_purified_genus", "00_genus_wide/03_gold_thin_audit"):
                self.assertIn(stage, stdout)
            for slug in ("Packera_anonyma", "Packera_aurea", "Packera_obovata"):
                for st in ("02_complex_triage", "03_intraspecific_clean", "04_gold_thin_audit"):
                    self.assertIn(f"{slug}/{st}", stdout)
            self.assertTrue(all("schweinitziana" not in ln.lower() for ln in cmd_lines))
            self.assertNotIn("Packera_schweinitziana", stdout)
            self.assertRegex(stdout, r"--triage-source\s+\S*00_genus_wide/01_raw_discovery")
            self.assertTrue(any("--skip-discovery" in ln for ln in cmd_lines) or all("01_eastern_discovery" not in ln for ln in cmd_lines))

    def test_resume_force_and_summarize_only_synthesis(self) -> None:
        """Verify --resume skips completed stages, --force reruns all, and --summarize-only compiles outputs."""
        from scripts.pipeline.resolve_eastern_na_packera import main as monograph_main

        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp).resolve()
            (ws / "run_pipeline.sh").write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            db_path = _seed_voucher_db(ws / "data" / "plantmatics.sqlite", {"Packera dubia": 16, "Packera aurea": 16, "Packera schweinitziana": 4})
            parent_id = "eastern_na_20261114_a7f3c2"
            runs_dir = ws / "runs"
            parent_dir = runs_dir / parent_id
            parent_dir.mkdir(parents=True, exist_ok=True)
            (runs_dir / "latest_eastern_na").symlink_to(parent_id)

            genus_rows = [
                {"catalogNumber": f"NCU_{ep}_{idx:03d}", "plant_individual_id": 1, "species_standardized": sp, "taxonomic_treatment_subtype": "Confirmed_Species", "gmm_cluster": str(cid)}
                for sp, ep, n_rows, cid in (("Packera dubia", "DUBIA", 16, 1), ("Packera aurea", "AUREA", 16, 2), ("Packera schweinitziana", "SCHWEINITZIANA", 4, 3))
                for idx in range(1, n_rows + 1)
            ]
            for st in ("01_raw_discovery", "02_purified_genus", "03_gold_thin_audit"):
                _seed_stage_outputs(ws, db_path, f"{parent_id}/00_genus_wide/{st}", genus_rows)
            for slug, sp, ep, neigh in (("Packera_dubia", "Packera dubia", "DUBIA", "Packera aurea"), ("Packera_aurea", "Packera aurea", "AUREA", "Packera dubia")):
                _seed_stage_outputs(
                    ws, db_path, f"{parent_id}/{slug}/02_complex_triage",
                    [
                        {"catalogNumber": f"NCU_{ep}_001", "plant_individual_id": 1, "species_standardized": sp, "taxonomic_treatment_subtype": "Confirmed_Type_Control", "gmm_cluster": "1"},
                        {"catalogNumber": f"NCU_{ep}_002", "plant_individual_id": 1, "species_standardized": sp, "taxonomic_treatment_subtype": "Misidentification_Reassign", "gmm_cluster": "2"},
                    ],
                    neighbors=[neigh],
                )

            def _fake_run(cmd, *args, **kwargs):
                return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="000000\n", stderr="")

            with mock.patch("subprocess.run", side_effect=_fake_run) as m_resume:
                self.assertEqual(monograph_main(["--resume", "--skip-harvest", "--db", str(db_path)], workspace_root=ws), 0)
            expected_incomplete = {
                f"{parent_id}/Packera_dubia/03_intraspecific_clean", f"{parent_id}/Packera_dubia/04_gold_thin_audit",
                f"{parent_id}/Packera_aurea/03_intraspecific_clean", f"{parent_id}/Packera_aurea/04_gold_thin_audit",
            }
            self.assertEqual(set(_extract_pipeline_run_ids(m_resume, ws)), expected_incomplete)

            with mock.patch("subprocess.run", side_effect=_fake_run) as m_force:
                self.assertEqual(monograph_main(["--parent-run-id", parent_id, "--force", "--skip-harvest", "--db", str(db_path)], workspace_root=ws), 0)
            forced_ids = set(_extract_pipeline_run_ids(m_force, ws))
            self.assertIn(f"{parent_id}/00_genus_wide/01_raw_discovery", forced_ids)
            self.assertIn(f"{parent_id}/Packera_dubia/02_complex_triage", forced_ids)
            self.assertTrue(expected_incomplete.issubset(forced_ids))

            _seed_stage_outputs(
                ws, db_path, f"{parent_id}/Packera_dubia/03_intraspecific_clean",
                [
                    {"catalogNumber": "NCU_DUBIA_001", "plant_individual_id": 1, "institutionCode": "NCU", "species_standardized": "Packera dubia", "taxonomic_treatment_subtype": "Confirmed_Type_Control", "synthesis_triage_priority": "RESOLVED", "gmm_cluster": "1", "cryptic_cluster_cohesion_pass": "FALSE", "delta_doy_resid": "2.1"},
                    {"catalogNumber": "NCU_DUBIA_003", "plant_individual_id": 1, "institutionCode": "NCU", "species_standardized": "Packera dubia", "recommended_determination": "Packera dubia", "taxonomic_treatment_subtype": "Putative_Cryptic_Species", "synthesis_triage_priority": "CRITICAL", "gmm_cluster": "2", "cryptic_cluster_cohesion_pass": "TRUE", "delta_doy_resid": "19.5", "synthesis_rationale": "Cohesive outcrop lineage"},
                ],
            )
            _seed_stage_outputs(ws, db_path, f"{parent_id}/Packera_dubia/04_gold_thin_audit",
                                [{"catalogNumber": "NCU_DUBIA_003", "plant_individual_id": 1, "species_standardized": "Packera dubia", "taxonomic_treatment_subtype": "Putative_Cryptic_Species", "gmm_cluster": "2", "cryptic_cluster_cohesion_pass": "TRUE"}])
            for st in ("03_intraspecific_clean", "04_gold_thin_audit"):
                _seed_stage_outputs(ws, db_path, f"{parent_id}/Packera_aurea/{st}",
                                    [{"catalogNumber": "NCU_AUREA_001", "plant_individual_id": 1, "species_standardized": "Packera aurea", "taxonomic_treatment_subtype": "Confirmed_Species", "gmm_cluster": "1"}])

            with mock.patch("subprocess.run", side_effect=_fake_run) as m_sum:
                self.assertEqual(monograph_main(["--parent-run-id", parent_id, "--summarize-only", "--db", str(db_path)], workspace_root=ws), 0)
            self.assertEqual(_extract_pipeline_run_ids(m_sum, ws), [])

            self.assertTrue((parent_dir / "00_genus_wide" / "Eastern_NA_Packera_Genus_Three_Stage_Comparison.csv").is_file())
            self.assertTrue((parent_dir / "00_genus_wide" / "Eastern_NA_Packera_Genus_Summary.md").is_file())
            for slug in ("Packera_dubia", "Packera_aurea"):
                self.assertTrue((parent_dir / slug / f"{slug}_Four_Stage_Comparison.csv").is_file())
                self.assertTrue((parent_dir / slug / f"{slug}_Complex_Resolution_Report.md").is_file())
            master_csv = parent_dir / "Eastern_NA_Packera_Master_Synthesis.csv"
            self.assertTrue(master_csv.is_file())
            self.assertTrue((parent_dir / "Eastern_NA_Packera_Monograph_Report.md").is_file())
            manifest_path = parent_dir / "monograph_run_manifest.json"
            self.assertTrue(manifest_path.is_file())
            self.assertIn(parent_id, json.dumps(json.loads(manifest_path.read_text(encoding="utf-8"))))

            with open(master_csv, "r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                fields = reader.fieldnames or []
                self.assertEqual(len(fields), 25)
                self.assertEqual(len(set(fields)), 25)
                for col in ("parent_run_id", "species_standardized", "median_coalescent_gdi", "public_backbone_clade", "n_tier0_abstained", "primary_verdict"):
                    self.assertIn(col, fields)
                self.assertGreaterEqual(len(list(reader)), 2)

            for md_file in parent_dir.rglob("*.md"):
                md_text = md_file.read_text(encoding="utf-8")
                self.assertNotIn("\u2013", md_text, f"En dash found in {md_file}")
                self.assertNotIn("\u2014", md_text, f"Em dash found in {md_file}")

    def test_focal_taxon_eastern_gate_and_triage_source_auto_neighbors(self) -> None:
        """Verify eastern_na --focal-taxon pruning and --auto-neighbors pinning to --triage-source."""
        import pandas as pd
        from scripts.core.db import get_db_connection, upsert_vouchers_records

        sh_src = (PROJECT_ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        fn_body = sh_src[sh_src.index("run_preflight_subset_filter() {"):sh_src.index("\nrun_preflight_subset_filter\n")]
        with tempfile.TemporaryDirectory() as tmp:
            root_data, runs_root = Path(tmp).resolve() / "root_data", Path(tmp).resolve() / "runs"
            run_dir = runs_root / "paupercula_02_complex_triage"
            (tables_dir := run_dir / "data" / "tables").mkdir(parents=True, exist_ok=True)
            (root_data / "cropped_patches" / "capitula").mkdir(parents=True, exist_ok=True)
            recs = [
                {"catalogNumber": c, "scientificName": sp, "species_standardized": sp, "stateProvince": st, "decimalLatitude": lat,
                 "decimalLongitude": lon, "type_status": ts, "determiner_tier": "Tier_1_Gold", "is_primary_duplicate": 1,
                 "duplicate_count": 1, "image_path": f"data/raw_vouchers/{c}.webp"}
                for c, sp, st, lat, lon, ts in (
                    ("PAUP_NC1", "Packera paupercula", "NC", 35.9, -79.0, ""), ("PAUP_NC2", "Packera paupercula", "NC", 35.8, -79.1, ""),
                    ("PAUP_TYPE", "Packera paupercula", "MT", 46.8, -110.0, "HOLOTYPE"), ("PAUP_MT", "Packera paupercula", "MT", 46.8, -110.0, ""),
                    ("PAUP_CO", "Packera paupercula", "CO", 39.7, -105.0, ""), ("PAUP_ND", "Packera paupercula", "ND", 46.9, -102.8, ""),
                    ("ANON_NC1", "Packera anonyma", "NC", 35.9, -79.0, ""), ("OBOV_NC1", "Packera obovata", "NC", 35.9, -79.0, ""),
                )
            ]
            with contextlib.closing(get_db_connection(db_path=root_data / "plantmatics.sqlite", init_schema=True)) as conn:
                upsert_vouchers_records(conn, recs)
            pd.DataFrame([r for r in recs if r["catalogNumber"] not in {"PAUP_MT", "PAUP_CO", "PAUP_ND"}]).to_csv(tables_dir / "curated_vouchers.csv", index=False)
            for rid, neigh in (("run_A", "Packera anonyma"), ("decoy_west", "Packera obovata")):
                _write_csv(runs_root / rid / "outputs" / "reports" / "empirical_taxon_complexes.csv", ["species_standardized", "consensus_complex_id"], [{"species_standardized": "Packera paupercula", "consensus_complex_id": 1}, {"species_standardized": neigh, "consensus_complex_id": 1}])
            (runs_root / "latest").symlink_to("decoy_west")
            env = {k: v for k, v in os.environ.items() if k not in ("PLANTMATICS_DB_PATH", "PLANTMATICS_DEMO_MODE", "PACKERA_DB_PATH", "PACKERA_DEMO_MODE")}
            env.update({"PROJECT_ROOT": str(PROJECT_ROOT), "RUN_DIR": str(run_dir), "PYTHON_MAIN": sys.executable, "GEOGRAPHIC_SCOPE": "eastern_na",
                        "INCLUDE_WESTERN": "false", "FOCAL_TAXON": "Packera paupercula", "AUTO_NEIGHBORS": "true", "TRIAGE_SOURCE": "run_A",
                        "REUSE_CANONICAL_DATA": "true", "SEED_DATA_DIR": str(root_data)})
            script = f"set -euo pipefail\n{fn_body}\nrun_preflight_subset_filter\n".replace('"${PROJECT_ROOT}/data"', f'"{root_data}"')
            res = subprocess.run(["bash", "-c", script], cwd=str(PROJECT_ROOT), capture_output=True, text=True, env=env, check=True)
            self.assertIn("Restored 3 vouchers", res.stdout)
            self.assertIn("[Preflight: Cis-Mississippi Eastern NA Scope]", res.stdout)
            self.assertEqual(set(pd.read_csv(tables_dir / "curated_vouchers.csv")["catalogNumber"]), {"PAUP_NC1", "PAUP_NC2", "PAUP_TYPE", "ANON_NC1"})
            self.assertNotIn("Packera obovata", res.stdout)

    def test_compile_complex_resolution_summary_full_queue_merged_k_and_stage2_pulls(self) -> None:
        """Verify compile_complex_resolution_summary prefers triage_queue_complex_full.csv, optimal_K_merged, dual CSVs, and Stage 2 purged misidentifications."""
        from scripts.pipeline.resolve_focal_species_complex import compile_complex_resolution_summary

        with tempfile.TemporaryDirectory() as tmp:
            ws = Path(tmp).resolve()
            stages = ("dubia_01_eastern_discovery", "dubia_02_complex_triage", "dubia_03_intraspecific_clean", "dubia_04_gold_thin_audit")
            for st in stages:
                st_dir = ws / "runs" / st
                if st == "dubia_02_complex_triage":
                    _write_csv(st_dir / "data" / "tables" / "triage_queue.csv", ["catalogNumber", "plant_individual_id", "species_standardized", "taxonomic_treatment_subtype", "gmm_cluster"],
                               [{"catalogNumber": "DUB_1", "plant_individual_id": 0, "species_standardized": "Packera dubia", "taxonomic_treatment_subtype": "Confirmed_Species", "gmm_cluster": "1"}])
                    _write_csv(st_dir / "data" / "tables" / "triage_queue_complex_full.csv", ["catalogNumber", "plant_individual_id", "species_standardized", "taxonomic_treatment_subtype", "recommended_determination", "gmm_cluster"], [
                        {"catalogNumber": "DUB_1", "plant_individual_id": 0, "species_standardized": "Packera dubia", "taxonomic_treatment_subtype": "Confirmed_Species", "recommended_determination": "Packera dubia", "gmm_cluster": "1"},
                        {"catalogNumber": "DUB_MIS2", "plant_individual_id": 0, "species_standardized": "Packera dubia", "taxonomic_treatment_subtype": "Misidentification_Reassign", "recommended_determination": "Packera anonyma", "gmm_cluster": "2"},
                        {"catalogNumber": "ANO_1", "plant_individual_id": 0, "species_standardized": "Packera anonyma", "taxonomic_treatment_subtype": "Confirmed_Species", "recommended_determination": "Packera anonyma", "gmm_cluster": "2"},
                    ])
                else:
                    _write_csv(st_dir / "data" / "tables" / "triage_queue.csv", ["catalogNumber", "plant_individual_id", "species_standardized", "taxonomic_treatment_subtype", "gmm_cluster"],
                               [{"catalogNumber": "DUB_1", "plant_individual_id": 0, "species_standardized": "Packera dubia", "taxonomic_treatment_subtype": "Confirmed_Type_Control", "gmm_cluster": "1"}])
                _write_csv(st_dir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv", ["K", "BIC", "is_optimal", "optimal_K", "optimal_K_merged"],
                           [{"K": 4, "BIC": -100.0, "is_optimal": "TRUE", "optimal_K": 4, "optimal_K_merged": 1}])

            csv_four, md_path = compile_complex_resolution_summary(focal_taxon="Packera dubia", prefix="dubia", workspace_root=ws)
            csv_cmp = ws / "outputs" / "reports" / "Packera_dubia_complex_stage_comparison.csv"
            self.assertTrue(csv_four.is_file() and csv_cmp.is_file())
            self.assertEqual(csv_four.read_text(encoding="utf-8"), csv_cmp.read_text(encoding="utf-8"))
            rows = {r["stage"]: r for r in csv.DictReader(csv_cmp.open("r", encoding="utf-8", newline=""))}
            self.assertEqual(rows["Stage_2_Complex_Triage"]["n_total_otus"], "3")
            self.assertEqual(rows["Stage_3_Intraspecific_Clean"]["gmm_optimal_k"], "1")
            md_txt = md_path.read_text(encoding="utf-8")
            self.assertIn("DUB_MIS2", md_txt)
            self.assertIn("NY00232445", md_txt)
            self.assertIn("P-MICH", md_txt)
            self.assertNotIn("\u2013", md_txt)
            self.assertNotIn("\u2014", md_txt)


if __name__ == "__main__":
    unittest.main()
