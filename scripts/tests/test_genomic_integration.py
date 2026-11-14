#!/usr/bin/env python3
"""
===============================================================================
Script: test_genomic_integration.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Batched verification suite for the independent Step 06b Genomic Evidence
    Stream: .gitignore safety, SQLite schema 20261121 -> 20261122 migration,
    1-to-1 master view cardinality under M/P duplicates, purge --runs-only,
    empty-directory no-ops, VCF/field ingestion, 4-tier voucher linking,
    LORO-CV shrinkage DAPC, 300-DPI PNG export, and Step 06 consensus isolation.
===============================================================================
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection, purge_runs_only
from scripts.core.db_schema import (
    PACKERA_SCHEMA_USER_VERSION,
    init_db_schema,
    migrate_database_schema,
)
from scripts.data_prep.ingest_genetic_data import (
    discover_genomic_inputs,
    run_genomic_ingestion,
)


class TestGenomicIntegration(unittest.TestCase):
    """End-to-end verification of the independent Step 06b genomic stream."""

    def test_gitignore_tracks_gitkeep_and_ignores_genomic_payloads(self):
        """Verify data/genomics/.gitkeep is tracked while sequence/variant/field files are ignored."""
        gitkeep = PROJECT_ROOT / "data" / "genomics" / ".gitkeep"
        self.assertTrue(gitkeep.is_file(), "data/genomics/.gitkeep must exist")

        chk_keep = subprocess.run(
            ["git", "check-ignore", "-q", "data/genomics/.gitkeep"],
            cwd=str(PROJECT_ROOT),
            check=False,
        )
        self.assertNotEqual(chk_keep.returncode, 0, "data/genomics/.gitkeep must NOT be git-ignored")

        ignored_paths = [
            "data/genomics/packera.vcf",
            "data/genomics/nested/sample.vcf.gz",
            "data/genomics/RAxML_bestTree.packera_raxml",
            "data/genomics/Packera_Field_Specimen_Database.csv",
            "data/genomics/admix.Q",
            "data/packera_sample.treefile",
        ]
        for rel in ignored_paths:
            res = subprocess.run(
                ["git", "check-ignore", "-q", rel],
                cwd=str(PROJECT_ROOT),
                check=False,
            )
            self.assertEqual(res.returncode, 0, f"Expected Git to ignore genomic payload path: {rel}")

    def test_schema_migration_master_cardinality_and_purge_runs_only(self):
        """Verify 20261121 -> 20261122 migration, 1-to-1 view cardinality, and purge --runs-only CLI."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_ws = Path(tmp)
            tmp_runs = tmp_ws / "runs" / "run_demo_old"
            tmp_runs.mkdir(parents=True, exist_ok=True)
            dummy_run_file = tmp_runs / "stale_output.csv"
            dummy_run_file.write_text("col\n1\n", encoding="utf-8")

            tmp_db = tmp_ws / "test_migrate.sqlite"
            conn = sqlite3.connect(str(tmp_db))
            try:
                init_db_schema(conn)
                # Simulate an authentic 20261121 database before Step 06b genomic tables/column
                conn.executescript(
                    """
                    DROP VIEW IF EXISTS v_integrative_triage_master;
                    DROP TABLE IF EXISTS genetic_lineage_assignments;
                    DROP TABLE IF EXISTS genomic_samples;
                    CREATE TABLE _ct_old AS SELECT run_id, catalogNumber, plant_individual_id, w_morph_kappa, w_edaphic_kappa, w_pheno_kappa, hopkins_adjusted_doy, probabilistic_consensus_score FROM curation_triage;
                    DROP TABLE curation_triage;
                    CREATE TABLE curation_triage (
                        run_id TEXT NOT NULL,
                        catalogNumber TEXT NOT NULL,
                        plant_individual_id INTEGER NOT NULL,
                        w_morph_kappa REAL,
                        w_edaphic_kappa REAL,
                        w_pheno_kappa REAL,
                        hopkins_adjusted_doy REAL,
                        probabilistic_consensus_score REAL,
                        PRIMARY KEY (run_id, catalogNumber, plant_individual_id)
                    );
                    DROP TABLE _ct_old;
                    PRAGMA user_version = 20261121;
                    INSERT INTO vouchers (catalogNumber, recordedBy, recordNumber, species_standardized, decimalLatitude, decimalLongitude, stateProvince, county, determiner_tier, is_primary_duplicate)
                    VALUES ('NCU00445952', 'PACSCH01', 'NCB1', 'Packera dubia', 35.1, -83.8, 'NC', 'Clay', 'Tier_1_Gold', 1),
                           ('NCU00470021', 'PACSER01', 'NCSP', 'Packera serpenticola', 35.09, -83.82, 'NC', 'Clay', 'Tier_1_Gold', 1);
                    INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y)
                    VALUES ('NCU00445952', 0, 100.0, 200.0),
                           ('NCU00470021', 0, 110.0, 210.0);
                    INSERT INTO extracted_leaves (catalogNumber, plant_individual_id, leaf_id, contour_xy_json)
                    VALUES ('NCU00445952', 0, 1, '[]');
                    INSERT INTO pipeline_runs (run_id, run_dir) VALUES ('run_demo_old', 'runs/run_demo_old');
                    INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier)
                    VALUES ('run_demo_old', 'NCU00445952', 0, 'Packera dubia', 'Tier_1_Gold');
                    INSERT INTO curation_triage (run_id, catalogNumber, plant_individual_id, w_morph_kappa, probabilistic_consensus_score)
                    VALUES ('run_demo_old', 'NCU00445952', 0, 0.85, 0.90);
                    """
                )
                conn.commit()
                migrate_database_schema(conn)
                uv = int(conn.execute("PRAGMA user_version;").fetchone()[0])
                self.assertEqual(uv, PACKERA_SCHEMA_USER_VERSION)
                self.assertEqual(uv, 20261122)

                ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
                self.assertIn("w_geno_kappa", ct_cols)

                # Insert M/P replicates on one voucher to confirm 1-to-1 view cardinality
                conn.execute(
                    """
                    INSERT INTO genomic_samples
                        (genetic_sample_id, catalogNumber, project_id, population_id, replicate_pair_id,
                         taxon_raw, species_standardized, is_putative_novel_or_hybrid_query, marker_platform,
                         hq_reads, missing_genotype_rate, geno_quality_pass, linkage_tier)
                    VALUES
                        ('PACSCHNCB1M', 'NCU00445952', 'PACSCH01-NCB1', 'NCB1', 'PACSCH_NCB1',
                         'Packera schweinitziana', 'Packera schweinitziana', 0, 'nextRAD',
                         120000, 0.10, 1, 'G1_Direct_Voucher'),
                        ('PACSCHNCB1P', 'NCU00445952', 'PACSCH01-NCB1', 'NCB1', 'PACSCH_NCB1',
                         'Packera schweinitziana', 'Packera schweinitziana', 0, 'nextRAD',
                         80000, 0.15, 1, 'G1_Direct_Voucher'),
                        ('PACSERNCS1M', 'NCU00470021', 'PACSER01-NCSP', 'NCS1', 'PACSER_NCS1',
                         'Packera serpenticola', 'Packera serpenticola', 0, 'nextRAD',
                         240000, 0.08, 1, 'G2_Population_Voucher_Proxy');
                    """
                )
                conn.execute(
                    """
                    INSERT INTO genetic_lineage_assignments
                        (run_id, genetic_sample_id, catalogNumber, population_id, geno_predicted_taxon,
                         geno_posterior_prob, geno_margin, genetic_cluster, q_max, is_geno_admixed)
                    VALUES
                        ('run_demo_old', 'PACSCHNCB1M', 'NCU00445952', 'NCB1', 'Packera schweinitziana', 0.99, 0.98, 'Geno_Cluster_1', 0.97, 0),
                        ('run_demo_old', 'PACSCHNCB1P', 'NCU00445952', 'NCB1', 'Packera schweinitziana', 0.95, 0.90, 'Geno_Cluster_1', 0.94, 0),
                        ('run_demo_old', 'PACSERNCS1M', 'NCU00470021', 'NCS1', 'Packera serpenticola', 0.96, 0.92, 'Geno_Cluster_2', 0.95, 0);
                    """
                )
                conn.commit()

                view_rows = conn.execute(
                    "SELECT catalogNumber, geno_predicted_taxon, geno_posterior_prob, geno_linkage_tier "
                    "FROM v_integrative_triage_master WHERE run_id = 'run_demo_old' ORDER BY catalogNumber;"
                ).fetchall()
                # Must stay 1-to-1 per (catalogNumber, plant_individual_id)
                self.assertEqual(len(view_rows), 1)
                self.assertEqual(view_rows[0][0], "NCU00445952")
                self.assertEqual(view_rows[0][1], "Packera schweinitziana")
                self.assertAlmostEqual(float(view_rows[0][2]), 0.99, places=3)
                self.assertEqual(view_rows[0][3], "G1_Direct_Voucher")
            finally:
                conn.close()

            # Verify CLI `main.py purge --runs-only --db <tmp_db>` clears Tier B while preserving Tier A
            purge_code = (
                "import sys; from pathlib import Path; "
                "import main; main.PROJECT_ROOT = Path(sys.argv[1]); "
                "sys.argv = ['main.py', 'purge', '--runs-only', '--db', sys.argv[2]]; "
                "main.main()"
            )
            env = dict(os.environ)
            env["PYTHONPATH"] = str(PROJECT_ROOT)
            proc = subprocess.run(
                [sys.executable, "-c", purge_code, str(tmp_ws), str(tmp_db)],
                cwd=str(PROJECT_ROOT),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, f"main.py purge --runs-only failed: {proc.stderr}")

            conn_chk = sqlite3.connect(str(tmp_db))
            try:
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM genomic_samples;").fetchone()[0], 3)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM vouchers;").fetchone()[0], 2)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0], 1)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM genetic_lineage_assignments;").fetchone()[0], 0)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM curation_triage;").fetchone()[0], 0)
                self.assertFalse(dummy_run_file.exists())
            finally:
                conn_chk.close()

    def test_empty_genomics_directory_is_clean_noop(self):
        """Verify empty or .gitkeep-only data/genomics triggers a zero-overhead no-op."""
        with tempfile.TemporaryDirectory() as tmp:
            geno_dir = Path(tmp) / "genomics"
            geno_dir.mkdir(parents=True, exist_ok=True)
            (geno_dir / ".gitkeep").write_text("", encoding="utf-8")
            out_dir = Path(tmp) / "parsed"

            assets = discover_genomic_inputs(geno_dir)
            self.assertFalse(assets.has_genetic_payload)

            chk = subprocess.run(
                [sys.executable, "scripts/data_prep/ingest_genetic_data.py", "--genomics-dir", str(geno_dir), "--check-only"],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("HAS_GENOMIC_INPUTS=0", chk.stdout)

            res = run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir, demo=True)
            self.assertEqual(res["status"], "skipped_empty")
            self.assertFalse(out_dir.exists())

            import argparse
            from unittest.mock import patch
            import main as root_main

            run_dir = Path(tmp) / "run_smoke"
            tables_dir = run_dir / "data" / "tables"
            figures_dir = run_dir / "outputs" / "figures"
            tables_dir.mkdir(parents=True, exist_ok=True)
            figures_dir.mkdir(parents=True, exist_ok=True)
            cli_args = argparse.Namespace(db=str(Path(tmp) / "dummy.sqlite"), demo=True)

            with patch.dict(os.environ, {"PACKERA_GENOMICS_DIR": str(geno_dir)}), patch.object(
                root_main.subprocess, "run"
            ) as mock_sub:
                res_empty = root_main._maybe_run_genomic_stream(
                    cli_args, run_dir, tables_dir, figures_dir, "Rscript"
                )
                self.assertIsNone(res_empty)
                self.assertEqual(mock_sub.call_count, 0)

            (geno_dir / "dummy.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
            with patch.dict(os.environ, {"PACKERA_GENOMICS_DIR": str(geno_dir)}), patch.object(
                root_main.subprocess, "run"
            ) as mock_sub:
                res_pop = root_main._maybe_run_genomic_stream(
                    cli_args, run_dir, tables_dir, figures_dir, "Rscript"
                )
                self.assertEqual(res_pop, run_dir / "outputs" / "genomics")
                self.assertEqual(mock_sub.call_count, 2)
                self.assertIn("ingest_genetic_data.py", " ".join(mock_sub.call_args_list[0].args[0]))
                self.assertIn("step06b_genomic_integration.R", " ".join(mock_sub.call_args_list[1].args[0]))

            with patch.dict(os.environ, {"PACKERA_GENOMICS_DIR": ""}), patch.object(
                root_main.subprocess, "run"
            ) as mock_sub:
                res_mocked_default = root_main._maybe_run_genomic_stream(
                    cli_args, run_dir, tables_dir, figures_dir, "Rscript"
                )
                self.assertIsNone(res_mocked_default)
                self.assertEqual(mock_sub.call_count, 0)

    def test_ingest_genetic_data_run5563_safeguards_paralog_filter_and_linker(self):
        """Verify smart-quote stripping, non-Packera exclusion, VAL1/PAR1 flagging, dropout QC, paralog filter, and G1-G4 linking."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            geno_dir = tmp_p / "genomics"
            out_dir = tmp_p / "parsed"
            geno_dir.mkdir(parents=True, exist_ok=True)
            db_path = tmp_p / "tmp_packera.sqlite"

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=PROJECT_ROOT)
            try:
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, species_raw, species_standardized, recordedBy, recordNumber, decimalLatitude, decimalLongitude, stateProvince, county, eventDate, determiner_tier) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                    [
                        ("NCU00445952", "Packera dubia", "Packera dubia", "B. Fuller", "PACSCH01-NCB1", 35.12, -83.85, "NC", "Clay", "2024-05-10", "Tier_1_Gold"),
                        ("NCU00445864", "Packera aurea", "Packera aurea", "B. Fuller", "PACAUR2-NCB1", 35.125, -83.855, "NC", "Clay", "2024-05-10", "Tier_1_Gold"),
                        ("NCU00470021", "Packera serpenticola", "Packera serpenticola", "B. Fuller", "PACSER-NCSP", 35.085, -83.81, "NC", "Clay", "2024-05-12", "Tier_1_Gold"),
                        ("NCU00499999", "Packera dubia", "Packera dubia", "B. Fuller", "PACTOM01-SCG1", 34.81, -79.45, "NC", "Scotland", "2024-04-20", "Tier_1_Gold"),
                    ],
                )
                conn.commit()
            finally:
                conn.close()

            with (geno_dir / "user_genetic_metadata.csv").open("w", encoding="utf-8", newline="") as fh:
                csv.writer(fh).writerows([
                    ["genetic_sample_id", "catalogNumber"],
                    ["PACTOMSCG1P", "NCU00499999"],
                ])

            # Stage Field Specimen Database CSV with smart quotes, non-Packera rows, VAL1/PAR1, and TX1 fallback
            field_rows = [
                ["Label_Printed", "Project_ID", "Collector_ID", "Flora", "State", "Family", "County", "Taxon", "Location", "Lat_Lon", "Det_by", "Collector", "Date_Collected", "Project", "Notes"],
                ["Yes", "“PACTOM-SCG1”", "F1", "NC", "NC", "Asteraceae", "Scotland", "Packera tomentosa", "Sandhills", "34.81, -79.45", "BF", "BF", "2024-04-20", "Packera", ""],
                ["Yes", "PACTOM-TX1", "F2", "TX", "TX", "Asteraceae", "Angelina", "Packera tomentosa", "Pinewoods", "31.22, -94.65", "BF", "BF", "2024-04-02", "Packera", ""],
                ["Yes", "PACMIL-VAL1", "F3", "VA", "VA", "Asteraceae", "Lee", "Packera sp. nov.", "Limestone barren", "36.661, -83.194", "BF", "BF", "2024-05-01", "Packera", "superficial resemblance to P. millefolium; possible P. paupercula parentage"],
                ["Yes", "PACPAU-PAR1", "F4", "PA", "PA", "Asteraceae", "Fayette", "Packera sp. nov. 1", "Sandstone boulders", "39.86386, -79.49601", "BF", "BF", "2024-05-18", "Packera", ""],
                ["Yes", "PYCMUT1", "F5", "NC", "NC", "Lamiaceae", "Orange", "Pycnanthemum muticum", "Piedmont", "35.9, -79.0", "BF", "BF", "2024-07-01", "Other", ""],
                ["Yes", "PLARUG1", "F6", "NC", "NC", "Plantaginaceae", "Orange", "Plantago rugelii", "Field", "35.91, -79.01", "BF", "BF", "2024-07-02", "Other", ""],
                ["Yes", "PACSER-NCSP", "F7", "NC", "NC", "Asteraceae", "Clay", "Packera serpenticola", "Buck Creek", "35.085, -83.81", "BF", "BF", "2024-05-12", "Packera", ""],
            ]
            with (geno_dir / "Packera_Field_Specimen_Database.csv").open("w", encoding="utf-8", newline="") as fh:
                csv.writer(fh).writerows(field_rows)

            samples = [
                "PACTOMSCG1M_sorted",
                "PACTOMSCG1P",
                "PACTOMTX2M",
                "PACSERNCS1M",
                "PACSCHNCB1M",
                "PACMILVAL1M",
                "PACPAUPAR1P",
                "PACPAUVAL1M",
                "PACANTL1M",
                "PACOBOJOEM",
            ]
            stats_lines = [
                "Sample\tHQ reads\tHQ reads mapped to reference\tpercent reads mapped to reference\tpercent repetitive reads\t# ./. genotypes\t# 0/0 genotypes\t# 0/1 genotypes\t# 1/1 genotypes",
                "PACTOMSCG1M\t220000\t200000\t90.9\t2.1\t0\t6\t2\t4",
                "PACTOMSCG1P\t210000\t195000\t92.8\t2.0\t0\t6\t2\t4",
                "PACTOMTX2M\t826498\t790000\t95.6\t1.9\t0\t6\t2\t4",
                "PACSERNCS1M\t246540\t220000\t89.2\t2.5\t0\t6\t2\t4",
                "PACSCHNCB1M\t180000\t165000\t91.6\t2.2\t0\t6\t2\t4",
                "PACMILVAL1M\t190000\t175000\t92.1\t2.3\t0\t6\t4\t2",
                "PACPAUPAR1P\t175000\t160000\t91.4\t2.4\t0\t6\t4\t2",
                "PACPAUVAL1M\t250000\t230000\t92.0\t2.1\t0\t6\t2\t4",
                "PACANTL1M\t160000\t140000\t87.5\t3.0\t0\t6\t2\t4",
                "PACOBOJOEM\t4123\t1200\t29.1\t5.0\t12\t0\t0\t0",
            ]
            (geno_dir / "packera_stats.txt").write_text("\n".join(stats_lines) + "\n", encoding="utf-8")

            vcf_lines = [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
            ]
            for li in range(1, 13):
                if li <= 4:
                    gts = ["0/0", "0/0", "0/0", "1/1", "1/1", "0/1", "0/1", "0/0", "1/1", "./."]
                elif li <= 8:
                    gts = ["1/1", "1/1", "1/1", "0/0", "1/1", "0/1", "0/1", "1/1", "0/0", "./."]
                else:
                    gts = ["0/0", "0/0", "0/0", "0/0", "1/1", "0/0", "0/1", "1/1", "1/1", "./."]
                vcf_lines.append(f"ctg{li}\t10\tL_NORM_{li}\tA\tG\t60\tPASS\t.\tGT\t" + "\t".join(gts))
            # Paralog locus: 100% heterozygous in all called samples (Ho = 1.0 > 0.65)
            vcf_lines.append("ctg98\t10\tL_PARALOG_HET\tA\tT\t60\tPASS\t.\tGT\t" + "\t".join(["0/1"] * 9 + ["./."]))
            # Ultra-low-MAF autapomorphy locus (still retained because Ho = 1/9 <= 0.65 and no MAF filter)
            vcf_lines.append("ctg99\t10\tL_LOWMAF\tC\tT\t60\tPASS\t.\tGT\t" + "\t".join(["0/1"] + ["0/0"] * 8 + ["./."]))
            (geno_dir / "packera.vcf").write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")

            res = run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir, db_path=db_path, demo=False)
            self.assertEqual(res["status"], "ingested")
            self.assertEqual(res["n_loci_raw"], 14)
            self.assertEqual(res["n_loci_retained"], 13)  # Paralog dropped, low-MAF kept!

            header_dos = (out_dir / "allele_dosage_matrix.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
            self.assertNotIn("L_PARALOG_HET", header_dos)
            self.assertIn("L_LOWMAF", header_dos)

            conn_chk = sqlite3.connect(str(db_path))
            try:
                rows = {
                    r[0]: r[1:]
                    for r in conn_chk.execute(
                        "SELECT genetic_sample_id, catalogNumber, taxon_raw, species_standardized, "
                        "is_putative_novel_or_hybrid_query, geno_quality_pass, linkage_tier, decimalLatitude "
                        "FROM genomic_samples;"
                    ).fetchall()
                }
            finally:
                conn_chk.close()

            self.assertNotIn("PYCMUT1", rows)
            self.assertNotIn("PLARUG1", rows)
            # Legacy PACTOM synonymy -> Packera tomentosa raw, Packera dubia standardized
            self.assertEqual(rows["PACTOMSCG1M"][1], "Packera tomentosa")
            self.assertEqual(rows["PACTOMSCG1M"][2], "Packera dubia")
            self.assertAlmostEqual(float(rows["PACTOMSCG1M"][6]), 34.81, places=3)
            # Explicit user_genetic_metadata G1 direct voucher match -> NCU00499999
            self.assertEqual(rows["PACTOMSCG1P"][0], "NCU00499999")
            self.assertEqual(rows["PACTOMSCG1P"][5], "G1_Direct_Voucher")
            # Collector-prefix disambiguated G2 population proxy match -> NCU00445952 (NOT NCU00445864)
            self.assertEqual(rows["PACSCHNCB1M"][0], "NCU00445952")
            self.assertEqual(rows["PACSCHNCB1M"][5], "G2_Population_Voucher_Proxy")
            # NCS1 -> NCSP G2 Population Proxy match -> NCU00470021
            self.assertEqual(rows["PACSERNCS1M"][0], "NCU00470021")
            self.assertEqual(rows["PACSERNCS1M"][5], "G2_Population_Voucher_Proxy")
            # TX2 -> TX1 G3 georeference fallback
            self.assertEqual(rows["PACTOMTX2M"][5], "G3_Georeferenced_Unvouchered")
            self.assertAlmostEqual(float(rows["PACTOMTX2M"][6]), 31.22, places=3)
            # Novel/hybrid focal query protection vs sympatric P. paupercula at VAL1
            self.assertEqual(rows["PACMILVAL1M"][3], 1)
            self.assertEqual(rows["PACPAUPAR1P"][3], 1)
            self.assertEqual(rows["PACPAUVAL1M"][2], "Packera paupercula")
            self.assertEqual(rows["PACPAUVAL1M"][3], 0)
            self.assertEqual(rows["PACPAUVAL1M"][5], "G3_Georeferenced_Unvouchered")
            self.assertAlmostEqual(float(rows["PACPAUVAL1M"][6]), 36.661, places=3)
            # G4 ungeoreferenced tip
            self.assertEqual(rows["PACANTL1M"][5], "G4_Ungeoreferenced_Tip")
            # High-dropout library retained in DB with geno_quality_pass = 0
            self.assertEqual(rows["PACOBOJOEM"][4], 0)

    def test_batched_r_step06b_and_step06_consensus_and_triage(self):
        """Run one batched Rscript verifying step06b outputs, 300-DPI PNG, 3-stream invariance, and 4-stream triage."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_p = Path(tmp)
            db_path = tmp_p / "synth.sqlite"
            raw_geno = tmp_p / "raw_genomics"
            parsed_geno = tmp_p / "parsed_genomics"
            empty_geno = tmp_p / "empty_genomics"
            tables_dir = tmp_p / "tables"
            base_tables_dir = tmp_p / "base_tables"
            figs_dir = tmp_p / "figures"
            reports_dir = tmp_p / "reports"
            for d in (raw_geno, parsed_geno, empty_geno, tables_dir, base_tables_dir, figs_dir, reports_dir):
                d.mkdir(parents=True, exist_ok=True)
            (empty_geno / ".gitkeep").write_text("", encoding="utf-8")

            # Seed SQLite schema 20261122 with 12 reference vouchers across 3 species
            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=PROJECT_ROOT)
            try:
                conn.execute("INSERT INTO pipeline_runs (run_id, run_dir) VALUES ('r_test', ?);", (str(tmp_p),))
                voucher_specs = [
                    ("V_DUB_1", "Packera dubia", 35.80, -79.00, "Tier_1_Gold", "PACTOM-SCG1"),
                    ("V_DUB_2", "Packera dubia", 35.85, -79.05, "Tier_1_Gold", "PACTOM-VAC1"),
                    ("V_DUB_3", "Packera dubia", 35.75, -79.10, "Tier_2_Silver", "V_DUB_3"),
                    ("V_DUB_4", "Packera dubia", 35.90, -78.95, "Tier_1_Gold", "V_DUB_4"),
                    ("V_SCH_1", "Packera schweinitziana", 36.10, -81.65, "Tier_1_Gold", "PACSCH-VT2"),
                    ("V_SCH_2", "Packera schweinitziana", 36.15, -81.70, "Tier_1_Gold", "PACSCH-NCB2"),
                    ("V_SCH_3", "Packera schweinitziana", 36.05, -81.60, "Tier_2_Silver", "V_SCH_3"),
                    ("V_SCH_4", "Packera schweinitziana", 36.20, -81.75, "Tier_1_Gold", "V_SCH_4"),
                    ("V_SER_1", "Packera serpenticola", 35.09, -83.81, "Tier_1_Gold", "PACSER-NCSP"),
                    ("V_SER_2", "Packera serpenticola", 35.08, -83.82, "Tier_1_Gold", "PACSER-NCS2"),
                    ("NCU00445952", "Packera dubia", 35.12, -83.85, "Tier_1_Gold", "PACSCH01-NCB1"),  # misidentified lookalike!
                    ("V_UNSEQ_1", "Packera dubia", 35.82, -79.02, "Tier_1_Gold", "V_UNSEQ_1"),
                ]
                for cat, sp, lat, lon, tier, recno in voucher_specs:
                    conn.execute(
                        "INSERT INTO vouchers (catalogNumber, species_standardized, decimalLatitude, decimalLongitude, stateProvince, determiner_tier, recordNumber) "
                        "VALUES (?, ?, ?, ?, 'NC', ?, ?);",
                        (cat, sp, lat, lon, tier, recno),
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, crown_x, crown_y) "
                        "VALUES (?, 0, 100.0, 200.0);",
                        (cat,),
                    )
                    conn.execute(
                        "INSERT INTO run_active_cohort (run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier) "
                        "VALUES ('r_test', ?, 0, ?, ?);",
                        (cat, sp, tier),
                    )
                conn.commit()
            finally:
                conn.close()

            # Build synthetic parsed genomics matching distinct genotypes for dubia vs schweinitziana vs serpenticola + VAL1 intermediate query
            sample_meta = [
                # (sid, cat, pid, pop, rep, taxon_raw, sp_std, is_q, hq, miss, qc, link)
                ("PACTOMSCG1M", "V_DUB_1", "PACTOM-SCG1", "SCG1", "PACTOM_SCG1", "Packera tomentosa", "Packera dubia", 0, 220000, 0.05, 1, "G1_Direct_Voucher"),
                ("PACTOMSCG1P", "V_DUB_1", "PACTOM-SCG1", "SCG1", "PACTOM_SCG1", "Packera tomentosa", "Packera dubia", 0, 210000, 0.06, 1, "G1_Direct_Voucher"),
                ("PACTOMVAC1M", "V_DUB_2", "PACTOM-VAC1", "VAC1", "PACTOM_VAC1", "Packera dubia", "Packera dubia", 0, 180000, 0.08, 1, "G1_Direct_Voucher"),
                ("PACTOMVAC1P", "V_DUB_2", "PACTOM-VAC1", "VAC1", "PACTOM_VAC1", "Packera dubia", "Packera dubia", 0, 175000, 0.09, 1, "G1_Direct_Voucher"),
                ("PACSCHVT2M", "V_SCH_1", "PACSCH-VT2", "VT2", "PACSCH_VT2", "Packera schweinitziana", "Packera schweinitziana", 0, 195000, 0.05, 1, "G1_Direct_Voucher"),
                ("PACSCHVT2P", "V_SCH_1", "PACSCH-VT2", "VT2", "PACSCH_VT2", "Packera schweinitziana", "Packera schweinitziana", 0, 190000, 0.06, 1, "G1_Direct_Voucher"),
                ("PACSCHNCB2M", "V_SCH_2", "PACSCH-NCB2", "NCB2", "PACSCH_NCB2", "Packera schweinitziana", "Packera schweinitziana", 0, 200000, 0.04, 1, "G1_Direct_Voucher"),
                ("PACSCHNCB2P", "V_SCH_2", "PACSCH-NCB2", "NCB2", "PACSCH_NCB2", "Packera schweinitziana", "Packera schweinitziana", 0, 185000, 0.07, 1, "G1_Direct_Voucher"),
                ("PACSERNCS1M", "V_SER_1", "PACSER-NCSP", "NCS1", "PACSER_NCS1", "Packera serpenticola", "Packera serpenticola", 0, 240000, 0.05, 1, "G2_Population_Voucher_Proxy"),
                ("PACSERNCS2M", "V_SER_2", "PACSER-NCS2", "NCS2", "PACSER_NCS2", "Packera serpenticola", "Packera serpenticola", 0, 230000, 0.06, 1, "G1_Direct_Voucher"),
                ("PACSCHNCB1M", "NCU00445952", "PACSCH01-NCB1", "NCB1", "PACSCH_NCB1", "Packera dubia", "Packera dubia", 0, 250000, 0.04, 1, "G1_Direct_Voucher"),
                ("PACMILVAL1M", "", "PACMIL-VAL1", "VAL1", "PACMIL_VAL1", "Packera sp. nov.", "", 1, 200000, 0.05, 1, "G3_Georeferenced_Unvouchered"),
                ("PACOBOJOEM", "", "PACOBO-JOE", "JOE", "PACOBO_JOE", "Packera obovata", "Packera obovata", 0, 4000, 0.99, 0, "G3_Georeferenced_Unvouchered"),
            ]
            with (parsed_geno / "parsed_genomic_samples.csv").open("w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["genetic_sample_id", "catalogNumber", "project_id", "population_id", "replicate_pair_id", "taxon_raw", "species_standardized", "is_putative_novel_or_hybrid_query", "hq_reads", "missing_genotype_rate", "geno_quality_pass", "linkage_tier", "marker_platform", "decimalLatitude", "decimalLongitude"])
                for row in sample_meta:
                    w.writerow([*row, "nextRAD", 35.5, -81.0])

            # Seed dosage profiles (40 loci cleanly separating 3 clades; PACSCHNCB1M matches schweinitziana; PACMILVAL1M 50/50 mix)
            rng = np.random.default_rng(42)
            loci = [f"L_{i}" for i in range(40)]
            dos_rows = []
            for row in sample_meta:
                if row[10] == 0:
                    continue
                sid, _, _, _, _, _, sp_std, is_q, *_ = row
                if sid == "PACSCHNCB1M":
                    base = np.array([1.0 if i < 20 else 0.0 for i in range(40)])  # matches schweinitziana!
                elif sp_std == "Packera dubia":
                    base = np.array([0.0] * 40)
                elif sp_std == "Packera schweinitziana":
                    base = np.array([1.0 if i < 20 else 0.0 for i in range(40)])
                elif sp_std == "Packera serpenticola":
                    base = np.array([0.0 if i < 20 else 1.0 for i in range(40)])
                else:
                    base = np.array([0.5 if i < 20 else 0.5 for i in range(40)])
                noisy = np.clip(base + rng.normal(0.0, 0.04, size=40), 0.0, 1.0)
                dos_rows.append([sid] + [f"{v:.4f}" for v in noisy])
            with (parsed_geno / "allele_dosage_matrix.tsv").open("w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["genetic_sample_id"] + loci)
                w.writerows(dos_rows)

            tree_nwk = (
                "(((PACTOMSCG1M:0.01,PACTOMSCG1P:0.01):0.02,(PACTOMVAC1M:0.01,PACTOMVAC1P:0.01):0.02):0.05,"
                "(((PACSCHVT2M:0.01,PACSCHVT2P:0.01):0.02,(PACSCHNCB2M:0.01,PACSCHNCB2P:0.01):0.02,PACSCHNCB1M:0.01):0.04,"
                "((PACSERNCS1M:0.01,PACSERNCS2M:0.01):0.04,PACMILVAL1M:0.03):0.02):0.03);"
            )
            (raw_geno / "RAxML_bestTree.packera_raxml").write_text(tree_nwk, encoding="utf-8")
            (parsed_geno / "genomic_manifest.json").write_text(
                json.dumps({"tree_files": [str(raw_geno / "RAxML_bestTree.packera_raxml")], "admixture_q_files": []}, indent=2),
                encoding="utf-8",
            )

            # Seed minimal Step 03-05 CSVs for 06_triage_dashboard_synthesis.R in both base_tables_dir and tables_dir
            v_hdr = ["catalogNumber", "plant_individual_id", "species_raw", "species_standardized", "determiner_tier", "decimalLatitude", "decimalLongitude", "stateProvince", "doy", "dayOfYear", "is_verified_anthesis", "morph_tier_source", "basal_leaf_count", "leaf_count"]
            m_hdr = ["catalogNumber", "plant_individual_id", "cda_predicted_taxon", "cda_posterior_prob", "cda_posterior_margin", "morphometric_outlier_flag", "misidentification_flag", "gmm_cluster", "gmm_uncertainty", "PC1", "PC2", "self_morph_mahal_d2", "self_morph_df"]
            f_hdr = ["catalogNumber", "plant_individual_id", "edaphic_best_fit_taxon", "rf_predicted_taxon", "rf_predicted_prob", "rf_given_taxon_prob", "edaphic_niche_overlap_d", "edaphic_z_score", "hopkins_adjusted_doy", "delta_doy_resid", "is_active_anthesis", "triage_category", "multimodal_conflict_flag", "loocv_accuracy", "spatial_block_cv_acc"]

            for tgt_dir in (base_tables_dir, tables_dir):
                with (tgt_dir / "curated_vouchers.csv").open("w", encoding="utf-8", newline="") as fh:
                    w = csv.writer(fh); w.writerow(v_hdr)
                    for cat, sp, lat, lon, tier, _ in voucher_specs:
                        w.writerow([cat, 0, sp, sp, tier, lat, lon, "NC", 125, 125, "TRUE", "Basal_Rosette", 4, 4])
                with (tgt_dir / "morphometric_flags.csv").open("w", encoding="utf-8", newline="") as fh:
                    w = csv.writer(fh); w.writerow(m_hdr)
                    for idx_v, (cat, sp, _, _, _, _) in enumerate(voucher_specs):
                        is_mis = cat == "NCU00445952"
                        pc1 = 0.8 if (sp == "Packera schweinitziana" or is_mis) else (-0.8 if sp == "Packera dubia" else 0.0)
                        pc2 = 0.9 if sp == "Packera serpenticola" else 0.1 * ((idx_v % 3) - 1)
                        d2 = 14.5 if is_mis else 1.2
                        # Give NCU00445952 ambiguous/borderline morph (pred=dubia 0.55) so ONLY strong G1 genome resolves reassignment!
                        w.writerow([cat, 0, "Packera dubia" if is_mis else sp, 0.55 if is_mis else 0.94, 0.10 if is_mis else 0.88, "FALSE", "FALSE", 1, 0.05, pc1, pc2, d2, 2])
                with (tgt_dir / "multimodal_conflict_flags.csv").open("w", encoding="utf-8", newline="") as fh:
                    w = csv.writer(fh); w.writerow(f_hdr)
                    for cat, sp, _, _, _, _ in voucher_specs:
                        w.writerow([cat, 0, sp, sp, 0.88, 0.88, 0.82, 0.4, 125.0, 2.0, "TRUE", "Concordant", "FALSE", 0.90, 0.85])
                with (tgt_dir / "leaf_efa_harmonics.csv").open("w", encoding="utf-8", newline="") as fh:
                    csv.writer(fh).writerow(["catalogNumber", "plant_individual_id", "leaf_id", "A1", "B1", "C1", "D1"])
                with (tgt_dir / "gmm_bayes_factors_summary.csv").open("w", encoding="utf-8", newline="") as fh:
                    csv.writer(fh).writerows([["species_standardized", "two_delta_bic"], ["Packera dubia", 4.0]])
                with (tgt_dir / "multimodal_spatial_rf_summary.csv").open("w", encoding="utf-8", newline="") as fh:
                    csv.writer(fh).writerows([["species_standardized", "spatial_block_cv_acc"], ["Packera dubia", 0.85]])

            # Execute single batched Rscript harness testing:
            # (1) empty no-op; (2) step06b outputs + SQLite + PNG; (3) Step 06 3-stream vs 4-stream invariants
            r_harness = tmp_p / "verify_batched.R"
            r_harness.write_text(
                f"""
                Sys.setenv(PACKERA_SQLITE_DB = "{db_path}", PACKERA_DB_PATH = "{db_path}", PACKERA_RUN_ID = "r_test", RUN_ID = "r_test")
                # 1. Empty genomic dir -> clean no-op
                R_empty <- system2("Rscript", c(
                  "scripts/analysis/step06b_genomic_integration.R",
                  "--genomic-dir", "{empty_geno}",
                  "--raw-genomics-dir", "{empty_geno}",
                  "--run-id", "r_test",
                  "--db-path", "{db_path}",
                  "--output-tables", "{base_tables_dir}",
                  "--output-figures", "{tmp_p / 'empty_figs'}"
                ))
                stopifnot(R_empty == 0L)
                stopifnot(!file.exists("{base_tables_dir / 'genomic_assignments.csv'}"))

                # 2. Populated Step 06b run
                R_pop <- system2("Rscript", c(
                  "scripts/analysis/step06b_genomic_integration.R",
                  "--genomic-dir", "{parsed_geno}",
                  "--raw-genomics-dir", "{raw_geno}",
                  "--run-id", "r_test",
                  "--db-path", "{db_path}",
                  "--output-tables", "{tables_dir}",
                  "--output-figures", "{figs_dir}"
                ))
                stopifnot(R_pop == 0L)
                stopifnot(file.exists("{tables_dir / 'genomic_assignments.csv'}"))
                stopifnot(file.exists("{tables_dir / 'genomic_model_metrics.csv'}"))
                stopifnot(file.exists("{figs_dir / 'fig06c_genomic_pcoa_tree_concordance.png'}"))
                stopifnot(file.info("{figs_dir / 'fig06c_genomic_pcoa_tree_concordance.png'}")$size > 15000)
                pdf_sidecars <- list.files("{figs_dir}", pattern = "\\\\.pdf$", full.names = TRUE, ignore.case = TRUE)
                stopifnot(length(pdf_sidecars) == 0L)

                ga <- utils::read.csv("{tables_dir / 'genomic_assignments.csv'}", stringsAsFactors = FALSE)
                stopifnot(nrow(ga) == 13L)
                ncu_row <- ga[ga$genetic_sample_id == "PACSCHNCB1M", , drop = FALSE]
                stopifnot(ncu_row$geno_predicted_taxon == "Packera schweinitziana", ncu_row$geno_posterior_prob >= 0.90)
                val_row <- ga[ga$genetic_sample_id == "PACMILVAL1M", , drop = FALSE]
                stopifnot(as.integer(val_row$is_putative_novel_or_hybrid_query) == 1L)
                stopifnot(as.integer(val_row$is_geno_admixed) == 1L)

                # 3. Run Step 06 WITHOUT genomics vs WITH genomics
                Sys.setenv(PACKERA_SKIP_FIGURES = "1")
                source("scripts/analysis/06_triage_dashboard_synthesis.R")
                opts_base <- list(
                  vouchers = "{base_tables_dir / 'curated_vouchers.csv'}",
                  harmonics = "{base_tables_dir / 'leaf_efa_harmonics.csv'}",
                  morphometrics = "{base_tables_dir / 'morphometric_flags.csv'}",
                  multimodal_flags = "{base_tables_dir / 'multimodal_conflict_flags.csv'}",
                  gmm_summary = "{base_tables_dir / 'gmm_bayes_factors_summary.csv'}",
                  niche_summary = "{base_tables_dir / 'multimodal_spatial_rf_summary.csv'}",
                  output_queue = "{base_tables_dir / 'triage_queue.csv'}",
                  output_plot = "{figs_dir / 'base_fig06a.png'}",
                  output_report = "{reports_dir / 'base_rep.md'}",
                  output_html = "",
                  min_affinity = 0.03,
                  output_complexes = "{reports_dir / 'base_comp.csv'}",
                  output_affinity = "{reports_dir / 'base_aff.csv'}",
                  output_neighbors = "{reports_dir / 'base_nb.txt'}",
                  output_complex_plot = "{figs_dir / 'base_fig06b.png'}",
                  genomic_dir = "",
                  focal_taxon = "",
                  intraspecific_mode = FALSE,
                  demo = TRUE
                )
                q_base <- build_triage_queue(opts_base)

                opts_geno <- opts_base
                opts_geno$vouchers <- "{tables_dir / 'curated_vouchers.csv'}"
                opts_geno$morphometrics <- "{tables_dir / 'morphometric_flags.csv'}"
                opts_geno$multimodal_flags <- "{tables_dir / 'multimodal_conflict_flags.csv'}"
                opts_geno$output_queue <- "{tables_dir / 'triage_queue.csv'}"
                opts_geno$output_report <- "{reports_dir / 'geno_rep.md'}"
                opts_geno$genomic_dir <- "{parsed_geno}"
                q_geno <- build_triage_queue(opts_geno)
                generate_taxonomic_revision_report(q_geno, opts_geno)

                # Unsequenced voucher MUST have bit-identical 3-stream scores!
                u_b <- q_base[q_base$catalogNumber == "V_UNSEQ_1", , drop = FALSE]
                u_g <- q_geno[q_geno$catalogNumber == "V_UNSEQ_1", , drop = FALSE]
                stopifnot(identical(u_b$probabilistic_consensus_score, u_g$probabilistic_consensus_score))
                stopifnot(identical(u_b$probabilistic_consensus_conflict, u_g$probabilistic_consensus_conflict))
                stopifnot(identical(as.integer(u_b$n_available_modalities), as.integer(u_g$n_available_modalities)))
                stopifnot(as.integer(u_g$n_available_modalities) <= 3L)

                # G1 & G2 sequenced vouchers activate 4th stream with G2 = 0.90 * G1
                d1 <- q_geno[q_geno$catalogNumber == "V_DUB_1", , drop = FALSE]
                s1 <- q_geno[q_geno$catalogNumber == "V_SER_1", , drop = FALSE]
                stopifnot(as.integer(d1$n_available_modalities) == 4L)
                stopifnot(as.integer(s1$n_available_modalities) == 4L)
                stopifnot(is.finite(d1$w_geno_kappa), d1$w_geno_kappa > 0.0)
                stopifnot(abs(s1$w_geno_kappa - round(0.90 * (d1$w_geno_kappa), 4)) <= 0.002)

                # Discordant G1 voucher NCU00445952 triggers Misidentification_Reassign -> Packera schweinitziana
                n1 <- q_geno[q_geno$catalogNumber == "NCU00445952", , drop = FALSE]
                stopifnot(n1$taxonomic_treatment_subtype == "Misidentification_Reassign")
                stopifnot(n1$recommended_determination == "Packera schweinitziana")

                stopifnot(file.exists("{tables_dir / 'genomic_field_triage.csv'}"))
                rep_txt <- paste(readLines("{reports_dir / 'geno_rep.md'}"), collapse = "\\n")
                stopifnot(grepl("Independent Genomic Evidence Stream", rep_txt))
                cat("BATCHED_R_GENOMIC_OK\\n")
                """,
                encoding="utf-8",
            )

            proc = subprocess.run(
                ["Rscript", str(r_harness)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                proc.returncode,
                0,
                f"Batched R verification failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
            )
            self.assertIn("BATCHED_R_GENOMIC_OK", proc.stdout)

            # Confirm SQLite genetic_lineage_assignments rows populated
            conn_chk = sqlite3.connect(str(db_path))
            try:
                cnt_gla = int(
                    conn_chk.execute(
                        "SELECT COUNT(*) FROM genetic_lineage_assignments WHERE run_id = 'r_test';"
                    ).fetchone()[0]
                )
                self.assertEqual(cnt_gla, 12)
            finally:
                conn_chk.close()


if __name__ == "__main__":
    unittest.main()
