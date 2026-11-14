#!/usr/bin/env python3
"""
Unit and end-to-end integration tests for Standalone Step 06 (Genomic &
Phylogenomic Evidence Stream) and Step 07 (Integrative Synthesis & Triage).

Covers:
  1. `.gitignore` tracking of `!data/genomics/.gitkeep` and flat payload ignores.
  2. SQLite WAL schema migrations (`20261121 -> 20261126`, `20261122 -> 20261126`,
     `20261123 -> 20261126`, `20261124 -> 20261126`, and `PLANTMATICS_SCHEMA_USER_VERSION - 1 -> 20261126`),
     all Phase 2 and public-reference columns in `genomic_samples`, `genetic_lineage_assignments`,
     `curation_triage`, and `v_integrative_triage_master` cardinality (`plant_individuals` preserved),
     plus `purge --runs-only` Tier A preservation vs Tier B cleanup.
  3. Empty/transient `data/genomics/` detection, flat-directory write guards, and
     `main.run_genomics_stream` zero-overhead no-op vs Phase 2 CLI orchestration.
  4. Multi-ploidy allele depth parsing (`_extract_cell_depths`), nQuire-style
     cytotype inference (`_infer_sample_cytotype`), empirical-Bayes `polyRAD`
     dosage (`_compute_eb_polyrad_dosage`), `HDplot` paralog filtering, and
     SNPsaurus `run_5563` `G1..G4` voucher linkage tiers.
  5. Off-target organellar/ITS FASTQ skimming (`skim_offtarget_organellar.py`),
     cytonuclear discordance detection, and 122 bp PHYLIP invariant-site
     coalescent delimitation (`step06_coalescent_delimitation.py`: `SVDquartets`,
     `BPP` `MSci` inputs, and empirical `coalescent_gdi`).
  6. Batched single-process R integration verifying `06_genomic_phylogenomics.R`
     (empty no-op + populated Phase 2 coancestry, `Q12` triangle classes, `f_b`,
     `SNPs2CF`, `MMRR`/Mantel, `small_n_guardrail_active`, 300-DPI `fig06a`/`fig06b`
     PNGs with zero `.pdf` outputs) and `07_triage_dashboard_synthesis.R`
     (unsequenced 3-stream invariance, Phase 2 `gdi < 0.20` oversplit block,
     `gdi > 0.70` / `Divergent_NonHybrid_Lineage` cryptic lineage support,
     `PACF1HYBM` F1 hybrid routing, `NCU00445952` assignment, and `fig07a` PNG).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib
import io
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from unittest import mock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main as root_main
from scripts.analysis import step06_coalescent_delimitation as coal_mod
from scripts.core.db import (
    get_db_connection,
    init_db_schema,
    purge_runs_only,
)
from scripts.core.db_schema import (
    PACKERA_SCHEMA_USER_VERSION,
    PLANTMATICS_SCHEMA_USER_VERSION,
    migrate_database_schema,
)
from scripts.data_prep import skim_offtarget_organellar as skim_mod
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

INGEST_MOD = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")
PUBLIC_DL_MOD = importlib.import_module("scripts.data_prep.download_public_packera_sequences")


def _build_synthetic_vcf_text(
    samples: list[str],
    locus_specs: list[tuple[str, list[tuple[str, int, int, int]]]],
) -> str:
    """Build a minimal VCF v4.2 payload with `GT:DP:AD` calls per locus."""
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for pos, (lid, cells) in enumerate(locus_specs, start=1):
        fmt_cells = [f"{gt}:{dp}:{r},{a}" for gt, dp, r, a in cells]
        lines.append(f"contig_{pos}\t{100 + pos}\t{lid}\tA\tG\t99\tPASS\t.\tGT:DP:AD\t" + "\t".join(fmt_cells))
    return "\n".join(lines) + "\n"


def _build_r1_fastq_gz_bytes(reads: list[tuple[str, str]]) -> bytes:
    """Compress synthetic FASTQ reads into `.fastq.gz` bytes."""
    buf = []
    for name, seq in reads:
        qual = "I" * len(seq)
        buf.append(f"@{name}\n{seq}\n+\n{qual}\n")
    return gzip.compress("".join(buf).encode("utf-8"))


class TestGenomicIntegration(unittest.TestCase):
    """End-to-end verification suite for Standalone Step 06 & Step 07 synthesis."""

    def test_gitignore_tracks_gitkeep_and_ignores_genomic_payloads(self) -> None:
        gitignore_path = PROJECT_ROOT / ".gitignore"
        text = gitignore_path.read_text(encoding="utf-8")
        lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]

        for required in (
            "data/genomics/*",
            "data/genomics/**",
            "!data/genomics/.gitkeep",
            "*.vcf",
            "*.vcf.gz",
            "*.fastq",
            "*.gz",
            "*.phylip",
            "*.nex",
            "*.bam",
            "*.tree",
            "*raxml*",
        ):
            self.assertIn(required, lines, f"Missing required .gitignore pattern: {required}")

        git_bin = shutil.which("git")
        if git_bin and (PROJECT_ROOT / ".git").exists():
            # Tracked .gitkeep must NOT be ignored; synthetic VCF/FASTQ/RAxML MUST be ignored.
            res_keep = subprocess.run(
                [git_bin, "check-ignore", "data/genomics/.gitkeep"],
                cwd=PROJECT_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(res_keep.returncode, 0, "data/genomics/.gitkeep must remain unignored")
            for probe in (
                "data/genomics/synthetic_test.vcf",
                "data/genomics/synthetic_reads_R1.fastq.gz",
                "data/genomics/RAxML_bestTree.synthetic",
            ):
                res_ign = subprocess.run(
                    [git_bin, "check-ignore", probe],
                    cwd=PROJECT_ROOT,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(res_ign.returncode, 0, f"Genomic payload not ignored by git: {probe}")

    def test_schema_migration_master_cardinality_and_purge_runs_only(self) -> None:
        self.assertEqual(PLANTMATICS_SCHEMA_USER_VERSION, 20261126)
        self.assertEqual(PACKERA_SCHEMA_USER_VERSION, 20261126)

        expected_gs_phase2 = {
            "hdplot_retained_loci",
            "inferred_cytotype",
            "nquire_delta_logl",
            "plastid_clade",
            "its_ribotype",
            "cytonuclear_discordant",
            "public_plastome_match",
            "public_its_match",
            "waxy_lineage_match",
        }
        expected_gla_phase2 = {
            "hybrid_index_h",
            "interspecific_ancestry_q12",
            "hybrid_generation_class",
            "coancestry_cluster",
            "f_branch_max",
            "f_branch_donor_taxon",
            "snaq_gamma_minor",
            "coalescent_gdi",
            "bpp_introgression_bf",
            "mmrr_beta_ibd",
            "mmrr_beta_ibe",
            "mmrr_p_ibe",
            "public_backbone_clade",
            "public_phytop_ils_index",
            "public_phytop_ih_index",
        }

        expected_ct_prop = {
            "genomic_propagation_tier",
            "tier0_calibrated_taxon",
            "tier0_calibrated_prob",
            "tier0_abstention_flag",
            "deme_mahal_d2",
            "cluster_propagated_gdi",
        }

        prev_schema_ver = PLANTMATICS_SCHEMA_USER_VERSION - 1
        for legacy_ver in (20261121, 20261122, 20261123, 20261124, prev_schema_ver):
            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td)
                db_path = tmp / f"legacy_{legacy_ver}.sqlite"
                conn = sqlite3.connect(str(db_path))
                init_db_schema(conn)
                if legacy_ver in (20261124, prev_schema_ver):
                    conn.executescript(
                        """
                        DROP VIEW IF EXISTS v_integrative_triage_master;
                        DROP TABLE IF EXISTS curation_triage;
                        CREATE TABLE curation_triage (
                            run_id TEXT NOT NULL,
                            catalogNumber TEXT NOT NULL,
                            plant_individual_id INTEGER NOT NULL DEFAULT 0,
                            taxonomic_status_call TEXT,
                            recommended_determination TEXT,
                            synthesis_triage_priority TEXT,
                            synthesis_triage_action TEXT,
                            synthesis_rationale TEXT,
                            taxonomic_treatment_broad TEXT,
                            taxonomic_treatment_subtype TEXT,
                            PRIMARY KEY (run_id, catalogNumber, plant_individual_id)
                        );
                        """
                    )
                else:
                    conn.executescript(
                        """
                        DROP VIEW IF EXISTS v_integrative_triage_master;
                        DROP TABLE IF EXISTS genomic_samples;
                        DROP TABLE IF EXISTS genetic_lineage_assignments;
                        """
                    )
                if legacy_ver == 20261122:
                    conn.execute(
                        """
                        CREATE TABLE genomic_samples (
                            genetic_sample_id TEXT PRIMARY KEY,
                            catalogNumber TEXT,
                            linkage_tier TEXT,
                            species_standardized TEXT,
                            geno_quality_pass INTEGER,
                            mapped_reads INTEGER
                        );
                        """
                    )
                    conn.execute(
                        """
                        CREATE TABLE genetic_lineage_assignments (
                            run_id TEXT NOT NULL,
                            genetic_sample_id TEXT NOT NULL,
                            catalogNumber TEXT,
                            geno_predicted_taxon TEXT,
                            PRIMARY KEY (run_id, genetic_sample_id)
                        );
                        """
                    )
                elif legacy_ver == 20261123:
                    conn.execute(
                        """
                        CREATE TABLE genomic_samples (
                            genetic_sample_id TEXT PRIMARY KEY,
                            catalogNumber TEXT,
                            linkage_tier TEXT,
                            species_standardized TEXT,
                            geno_quality_pass INTEGER,
                            mapped_reads INTEGER,
                            hdplot_retained_loci INTEGER,
                            inferred_cytotype TEXT,
                            nquire_delta_logl REAL,
                            plastid_clade TEXT,
                            its_ribotype TEXT,
                            cytonuclear_discordant INTEGER
                        );
                        """
                    )
                    conn.execute(
                        """
                        CREATE TABLE genetic_lineage_assignments (
                            run_id TEXT NOT NULL,
                            genetic_sample_id TEXT NOT NULL,
                            catalogNumber TEXT,
                            geno_predicted_taxon TEXT,
                            hybrid_index_h REAL,
                            interspecific_ancestry_q12 REAL,
                            hybrid_generation_class TEXT,
                            coancestry_cluster TEXT,
                            f_branch_max REAL,
                            f_branch_donor_taxon TEXT,
                            snaq_gamma_minor REAL,
                            coalescent_gdi REAL,
                            bpp_introgression_bf REAL,
                            mmrr_beta_ibd REAL,
                            mmrr_beta_ibe REAL,
                            mmrr_p_ibe REAL,
                            PRIMARY KEY (run_id, genetic_sample_id)
                        );
                        """
                    )
                conn.execute(f"PRAGMA user_version = {legacy_ver};")
                conn.commit()

                migrate_database_schema(conn)
                uv_after = conn.execute("PRAGMA user_version;").fetchone()[0]
                self.assertEqual(uv_after, 20261126)
                init_db_schema(conn)

                gs_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genomic_samples);").fetchall()}
                gla_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genetic_lineage_assignments);").fetchall()}
                ct_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(curation_triage);").fetchall()}
                view_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}

                self.assertTrue(expected_gs_phase2.issubset(gs_cols), f"Missing gs cols from {legacy_ver}: {expected_gs_phase2 - gs_cols}")
                self.assertTrue(expected_gla_phase2.issubset(gla_cols), f"Missing gla cols from {legacy_ver}: {expected_gla_phase2 - gla_cols}")
                self.assertTrue(expected_ct_prop.issubset(ct_cols), f"Missing ct cols from {legacy_ver}: {expected_ct_prop - ct_cols}")
                self.assertTrue(expected_gs_phase2.issubset(view_cols), f"View missing gs cols: {expected_gs_phase2 - view_cols}")
                self.assertTrue(expected_gla_phase2.issubset(view_cols), f"View missing gla cols: {expected_gla_phase2 - view_cols}")
                self.assertTrue(expected_ct_prop.issubset(view_cols), f"View missing ct cols: {expected_ct_prop - view_cols}")
                conn.close()

        # Verify cardinality invariance on multi-plant sheets + multiple genomic replicates, then purge --runs-only
        with tempfile.TemporaryDirectory() as td:
            ws = Path(td)
            data_dir = ws / "data"
            tables_dir = data_dir / "tables"
            runs_dir = ws / "runs" / "run_alpha"
            for d in (tables_dir, runs_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path = data_dir / "plantmatics.sqlite"

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=ws)
            try:
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, species_standardized, determiner_tier) VALUES (?, ?, ?);",
                    [
                        ("NCU001", "Packera dubia", "Tier_1_Gold"),
                        ("NCU002", "Packera schweinitziana", "Tier_1_Gold"),
                    ],
                )
                conn.executemany(
                    "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, rosette_centroid_x) VALUES (?, ?, ?);",
                    [
                        ("NCU001", 0, 100.0),
                        ("NCU001", 1, 200.0),
                        ("NCU002", 0, 150.0),
                    ],
                )
                # Two Tier A genomic replicate samples linked to NCU001 must NOT inflate master view rows!
                conn.executemany(
                    """
                    INSERT INTO genomic_samples (
                        genetic_sample_id, catalogNumber, linkage_tier, species_standardized,
                        geno_quality_pass, mapped_reads, inferred_cytotype, hdplot_retained_loci
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?);
                    """,
                    [
                        ("PACTOMNCU1M", "NCU001", "G1_Direct_Voucher", "Packera dubia", 1, 500000, "4x", 1200),
                        ("PACTOMNCU1P", "NCU001", "G2_Population_Voucher_Proxy", "Packera dubia", 1, 420000, "4x", 1180),
                    ],
                )
                conn.execute(
                    "INSERT INTO pipeline_runs (run_id, run_dir, intraspecific_mode) VALUES ('run_alpha', ?, 0);",
                    (str(runs_dir),),
                )
                conn.executemany(
                    """
                    INSERT INTO run_active_cohort (
                        run_id, catalogNumber, plant_individual_id,
                        species_standardized, determiner_tier, authority_weight
                    ) VALUES (?, ?, ?, ?, ?, ?);
                    """,
                    [
                        ("run_alpha", "NCU001", 0, "Packera dubia", "Tier_1_Gold", 2.0),
                        ("run_alpha", "NCU001", 1, "Packera dubia", "Tier_1_Gold", 2.0),
                        ("run_alpha", "NCU002", 0, "Packera schweinitziana", "Tier_1_Gold", 2.0),
                    ],
                )
                conn.executemany(
                    """
                    INSERT INTO genetic_lineage_assignments (
                        run_id, genetic_sample_id, catalogNumber, geno_predicted_taxon,
                        geno_posterior_prob, hybrid_index_h, interspecific_ancestry_q12,
                        hybrid_generation_class, coalescent_gdi
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
                    """,
                    [
                        ("run_alpha", "PACTOMNCU1M", "NCU001", "Packera dubia", 0.99, 0.04, 0.05, "Pure_Parent_1", 0.88),
                        ("run_alpha", "PACTOMNCU1P", "NCU001", "Packera dubia", 0.96, 0.06, 0.07, "Pure_Parent_1", 0.88),
                    ],
                )
                conn.execute(
                    """
                    INSERT INTO curation_triage (
                        run_id, catalogNumber, plant_individual_id,
                        recommended_determination, synthesis_triage_priority, taxonomic_treatment_subtype
                    ) VALUES ('run_alpha', 'NCU001', 0, 'Packera dubia', 'RESOLVED', 'Confirmed_Species');
                    """
                )
                conn.commit()

                master_count = conn.execute("SELECT COUNT(*) FROM v_integrative_triage_master;").fetchone()[0]
                self.assertEqual(
                    master_count,
                    3,
                    "v_integrative_triage_master MUST preserve exact plant_individuals cardinality without duplicate join inflation",
                )
            finally:
                conn.close()

            # Write dummy Step 01 curated_vouchers (preserved) + Tier B report CSV (purged)
            (tables_dir / "curated_vouchers.csv").write_text("catalogNumber\nNCU001\nNCU002\n", encoding="utf-8")
            (tables_dir / "genomic_assignments.csv").write_text("genetic_sample_id\nPACTOMNCU1M\n", encoding="utf-8")

            summary = purge_runs_only(
                workspace_root=ws,
                db_path=db_path,
                output_dir=data_dir,
                refresh_env_buffers=False,
                vacuum=True,
            )
            self.assertGreater(summary["purged_db_rows"], 0)
            self.assertTrue((tables_dir / "curated_vouchers.csv").exists())
            self.assertFalse((tables_dir / "genomic_assignments.csv").exists())

            conn_chk = get_db_connection(db_path=db_path, init_schema=True, workspace_root=ws)
            try:
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM vouchers;").fetchone()[0], 2)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM plant_individuals;").fetchone()[0], 3)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM genomic_samples;").fetchone()[0], 2)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM genetic_lineage_assignments;").fetchone()[0], 0)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM curation_triage;").fetchone()[0], 0)
                self.assertEqual(conn_chk.execute("SELECT COUNT(*) FROM pipeline_runs;").fetchone()[0], 0)
            finally:
                conn_chk.close()

    def test_empty_genomics_transient_and_flat_dir_guards(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            geno_dir.mkdir(parents=True, exist_ok=True)
            out_dir = tmp / "parsed_out"
            run_dir = tmp / "runs" / "latest"
            tables_dir = tmp / "data_tables"
            figs_dir = tmp / "figures"
            for d in (run_dir, tables_dir, figs_dir):
                d.mkdir(parents=True, exist_ok=True)

            # Populate only ignored / transient / raw FASTQ artifacts -> zero-overhead no-op!
            (geno_dir / ".gitkeep").write_text("", encoding="utf-8")
            (geno_dir / "reads.fastq.gz.crdownload").write_bytes(b"partial")
            (geno_dir / "archive.tar.gz.part").write_bytes(b"partial")
            (geno_dir / "scratch.tmp").write_text("tmp", encoding="utf-8")
            (geno_dir / ".hidden.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
            (geno_dir / "raw_sample_R1.fastq.gz").write_bytes(b"raw")
            (geno_dir / "aligned.bam").write_bytes(b"bam")

            assets = INGEST_MOD.discover_genomic_inputs(geno_dir)
            self.assertFalse(assets.has_genetic_payload)

            res = INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir)
            self.assertEqual(res["status"], "skipped_empty")
            self.assertEqual(len(list(geno_dir.iterdir())), 7, "Empty no-op must never write into data/genomics/")

            # Flat directory guard: writing parsed outputs into a subdir of genomics_dir must fail fast
            with self.assertRaises(ValueError):
                INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=geno_dir / "parsed_subdir")

            # Empty main.run_genomics_stream -> None, 0 subprocess calls
            cli_args = type(
                "Ns",
                (),
                {
                    "intraspecific_mode": False,
                    "focal_taxon": None,
                    "db_path": None,
                    "demo": False,
                    "genomics_dir": None,
                },
            )()
            with mock.patch.dict(os.environ, {"PLANTMATICS_GENOMICS_DIR": str(geno_dir)}), mock.patch("subprocess.run") as m_run:
                ret_empty = root_main.run_genomics_stream(cli_args, run_dir, tables_dir, figs_dir, "Rscript")
                self.assertIsNone(ret_empty)
                m_run.assert_not_called()

            # Populated VCF -> runs 06_ingest_genetic_data.py then standalone 06_genomic_phylogenomics.R
            (geno_dir / "packera_test.vcf").write_text(
                _build_synthetic_vcf_text(
                    ["PACTOMSCG1M", "PACSCHNCB1M"],
                    [("loc1", [("0/0", 20, 20, 0), ("1/1", 20, 0, 20)])],
                ),
                encoding="utf-8",
            )
            with mock.patch.dict(os.environ, {"PLANTMATICS_GENOMICS_DIR": str(geno_dir)}), mock.patch("subprocess.run") as m_run:
                ret_pop = root_main.run_genomics_stream(cli_args, run_dir, tables_dir, figs_dir, "Rscript")
                self.assertEqual(ret_pop, run_dir / "outputs" / "genomics")
                self.assertEqual(m_run.call_count, 2)
                cmd0 = " ".join(str(x) for x in m_run.call_args_list[0].args[0])
                cmd1 = " ".join(str(x) for x in m_run.call_args_list[1].args[0])
                self.assertIn("06_ingest_genetic_data.py", cmd0)
                self.assertIn("06_genomic_phylogenomics.R", cmd1)
                self.assertIn("fig06a_genomic_pcoa_tree_concordance.png", cmd1)
                self.assertIn("fig06b_reticulation_gdi_landscape.png", cmd1)

    def test_hdplot_eb_dosage_cytotype_and_run5563_linker(self) -> None:
        # 1. Dual GATK vs BBMap AD field depth parsing
        r_g, a_g = INGEST_MOD._extract_cell_depths(["0/1", "40", "12,28"], {"GT": 0, "DP": 1, "AD": 2}, 0.5)
        self.assertEqual((r_g, a_g), (12.0, 28.0))
        r_b, a_b = INGEST_MOD._extract_cell_depths(["0/1", "40", "18"], {"GT": 0, "DP": 1, "AD": 2}, 0.5)
        self.assertEqual((r_b, a_b), (22.0, 18.0))

        # 2. Allele-balance GMM cytotype inference (2x, 4x, 6x, and user_ploidy integer fallback)
        n_loci = 40
        ref_2x = np.full(n_loci, 20.0, dtype=np.float32)
        alt_2x = np.full(n_loci, 20.0, dtype=np.float32)
        dos_2x = np.full(n_loci, 0.5, dtype=np.float32)
        cyto_2x, _ = INGEST_MOD._infer_sample_cytotype(ref_2x, alt_2x, dos_2x, None)
        self.assertEqual(cyto_2x, "2x")

        ref_4x = np.array([30.0 if i % 2 == 0 else 10.0 for i in range(n_loci)], dtype=np.float32)
        alt_4x = np.array([10.0 if i % 2 == 0 else 30.0 for i in range(n_loci)], dtype=np.float32)
        cyto_4x, dlog_4x = INGEST_MOD._infer_sample_cytotype(ref_4x, alt_4x, dos_2x, None)
        self.assertEqual(cyto_4x, "4x")
        self.assertIsNotNone(dlog_4x)
        self.assertGreater(float(dlog_4x), 2.0)

        hex_pairs = [(52.0, 8.0)] * 30 + [(8.0, 52.0)] * 30 + [(40.0, 20.0)] * 10 + [(20.0, 40.0)] * 10
        ref_6x = np.array([p[0] for p in hex_pairs], dtype=np.float32)
        alt_6x = np.array([p[1] for p in hex_pairs], dtype=np.float32)
        dos_6x = np.full(len(hex_pairs), 0.5, dtype=np.float32)
        cyto_6x, dlog_6x = INGEST_MOD._infer_sample_cytotype(ref_6x, alt_6x, dos_6x, None)
        self.assertEqual(cyto_6x, "6x")
        self.assertIsNotNone(dlog_6x)
        self.assertGreater(float(dlog_6x), 2.0)

        cyto_fb, _ = INGEST_MOD._infer_sample_cytotype(ref_2x[:5], alt_2x[:5], dos_2x[:5], 6)
        self.assertEqual(cyto_fb, "6x")

        # 3. Empirical-Bayes polyRAD posterior allele dosage for 4x hets called as diploid 0/1
        dosage_sub = np.array([[0.5, 0.0], [0.5, 1.0], [0.5, 0.5]], dtype=np.float32)
        ref_sub = np.array([[30.0, 30.0], [10.0, 0.0], [20.0, 20.0]], dtype=np.float32)
        alt_sub = np.array([[10.0, 0.0], [30.0, 30.0], [20.0, 20.0]], dtype=np.float32)
        eb = INGEST_MOD._compute_eb_polyrad_dosage(dosage_sub, ref_sub, alt_sub, ["4x", "4x", "2x"])
        self.assertEqual(eb.shape, (3, 2))
        self.assertTrue(np.all(np.isfinite(eb)))
        self.assertTrue(np.all((eb >= 0.0) & (eb <= 1.0)))
        self.assertLess(float(eb[0, 0]), 0.40)    # 4x ref=30, alt=10 -> AAAB simplex ~0.25
        self.assertGreater(float(eb[1, 0]), 0.60) # 4x ref=10, alt=30 -> ABBB simplex ~0.75

        # 4. End-to-end HDplot paralog filter + G1..G4 linker on staged SNPsaurus fixture
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            out_dir = tmp / "parsed"
            geno_dir.mkdir(parents=True, exist_ok=True)
            db_path = tmp / "packera_test.sqlite"

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.executemany(
                    """
                    INSERT INTO vouchers (
                        catalogNumber, recordNumber, recordedBy, scientificName,
                        species_standardized, decimalLatitude, decimalLongitude,
                        stateProvince, determiner_tier
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
                    """,
                    [
                        ("NCU00445952", "PACSCH01-NCB1", "Fuller", "Packera schweinitziana", "Packera schweinitziana", 35.80, -82.90, "NC", "Tier_1_Gold"),
                        ("NCU00112233", "PACTOM-SCG1-01", "Fuller", "Packera tomentosa", "Packera dubia", 35.10, -79.20, "NC", "Tier_1_Gold"),
                    ],
                )
                conn.commit()
            finally:
                conn.close()

            samples = [
                "NCU00445952",    # G1 direct via catalogNumber
                "PACTOMSCG1M",    # G2 population proxy matching SCG1 + TOM->dubia
                "PACTOMSCG1P",    # G2 replicate pair member
                "PACMILVAL1M",    # G3 georeferenced via user_genetic_metadata (>=15 4x hets)
                "PACPAUVAL1M",    # G3 georeferenced via user_genetic_metadata
                "PACXYZ99M",      # G4 ungeoreferenced tip
            ]

            locus_specs = [
                # Excess heterozygosity H = 1.0 > 0.60 across 6 samples -> HDplot_H>0.60 filtered
                ("L_PARALOG_HET", [("0/1", 40, 20, 20) for _ in samples]),
                # Strong read-ratio bias |D| > 4.0 across 3/6 hets (H = 0.50 <= 0.60, ref=36, alt=4) -> HDplot_|D|>4.0 filtered
                ("L_PARALOG_DBIAS", [("0/1", 40, 36, 4)] * 3 + [("0/0", 40, 40, 0)] * 3),
                # Monomorphic (all 0/0) -> monomorphic_or_sparse filtered
                ("L_LOWMAF", [("0/0", 30, 30, 0) for _ in samples]),
            ]
            # 20 clean informative loci (>=15 hets so PACMILVAL1M infers 4x cytotype cleanly)
            for k in range(20):
                cells = [
                    ("1/1", 30, 0, 30),   # NCU00445952
                    ("0/0", 30, 30, 0),   # PACTOMSCG1M
                    ("0/0", 30, 30, 0),   # PACTOMSCG1P
                    ("0/1", 40, 30, 10) if k % 2 == 0 else ("0/1", 40, 10, 30),  # PACMILVAL1M (4x signature)
                    ("1/1", 30, 0, 30),   # PACPAUVAL1M
                    ("0/1", 30, 15, 15),  # PACXYZ99M
                ]
                locus_specs.append((f"L_GOOD_{k + 1}", cells))

            (geno_dir / "packera.vcf").write_text(_build_synthetic_vcf_text(samples, locus_specs), encoding="utf-8")
            (geno_dir / "packera_key.txt").write_text(
                "SampleID\tLibrary\tBarcode\n"
                + "\n".join(f"{s}\tlib1\tAA{i:02d}" for i, s in enumerate(samples))
                + "\n",
                encoding="utf-8",
            )
            (geno_dir / "packera_stats.txt").write_text(
                "Sample\tHQ reads\tHQ reads mapped to reference\tpercent reads mapped to reference\tpercent repetitive reads\t# ./. genotypes\t# 0/0 genotypes\t# 0/1 genotypes\t# 1/1 genotypes\n"
                + "\n".join(f"{s}\t400000\t360000\t90.0\t4.0\t0\t8\t4\t7" for s in samples)
                + "\n",
                encoding="utf-8",
            )
            with open(geno_dir / "user_genetic_metadata.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["genetic_sample_id", "catalogNumber", "decimalLatitude", "decimalLongitude", "taxon_raw"])
                wr.writerow(["PACMILVAL1M", "", "36.15", "-81.65", "Packera millefolium"])
                wr.writerow(["PACPAUVAL1M", "", "36.15", "-81.65", "Packera paupercula"])

            res = INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir, db_path=db_path)
            self.assertEqual(res["status"], "ingested")
            self.assertEqual(res["n_samples_total"], 6)
            self.assertEqual(res["n_loci_raw"], 23)
            self.assertEqual(res["n_loci_retained"], 20)

            # Verify HDplot TSV filtered reasons
            hd_rows = {}
            with open(out_dir / "hdplot_locus_diagnostics.tsv", "r", encoding="utf-8", newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    hd_rows[row["locus_id"]] = row
            self.assertEqual(hd_rows["L_PARALOG_HET"]["hdplot_keep"], "0")
            self.assertEqual(hd_rows["L_PARALOG_HET"]["filter_reason"], "HDplot_H>0.60")
            self.assertEqual(hd_rows["L_PARALOG_DBIAS"]["hdplot_keep"], "0")
            self.assertIn("HDplot_|D|>", hd_rows["L_PARALOG_DBIAS"]["filter_reason"])
            self.assertEqual(hd_rows["L_LOWMAF"]["hdplot_keep"], "0")
            self.assertEqual(hd_rows["L_LOWMAF"]["filter_reason"], "monomorphic_or_sparse")
            self.assertEqual(hd_rows["L_GOOD_1"]["hdplot_keep"], "1")
            self.assertEqual(hd_rows["L_GOOD_1"]["filter_reason"], "PASS")

            # Verify parsed_genomic_samples.csv tiers & SQLite persistence
            samp_rows = {}
            with open(out_dir / "parsed_genomic_samples.csv", "r", encoding="utf-8", newline="") as fh:
                for row in csv.DictReader(fh):
                    samp_rows[row["genetic_sample_id"]] = row
            self.assertEqual(samp_rows["NCU00445952"]["linkage_tier"], "G1_Direct_Voucher")
            self.assertEqual(samp_rows["NCU00445952"]["catalogNumber"], "NCU00445952")
            self.assertEqual(samp_rows["PACTOMSCG1M"]["linkage_tier"], "G2_Population_Voucher_Proxy")
            self.assertEqual(samp_rows["PACTOMSCG1M"]["catalogNumber"], "NCU00112233")
            self.assertEqual(samp_rows["PACMILVAL1M"]["linkage_tier"], "G3_Georeferenced_Unvouchered")
            self.assertEqual(samp_rows["PACMILVAL1M"]["inferred_cytotype"], "4x")
            self.assertEqual(samp_rows["PACXYZ99M"]["linkage_tier"], "G4_Ungeoreferenced_Tip")

            conn_c = sqlite3.connect(str(db_path))
            try:
                n_db = conn_c.execute("SELECT COUNT(*) FROM genomic_samples;").fetchone()[0]
                self.assertEqual(n_db, 6)
            finally:
                conn_c.close()

    def test_offtarget_organellar_skim_and_coalescent_delimitation(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            raw_dir = tmp / "raw_reads"
            out_dir = tmp / "coalescent_out"
            geno_dir.mkdir(parents=True, exist_ok=True)
            raw_dir.mkdir(parents=True, exist_ok=True)

            rbcl_seed = skim_mod._SEED_PLASTID_ANCHORS["Packera_rbcL_core"]
            its_seed = skim_mod._SEED_ITS_ANCHORS["Packera_ITS_58S_core"]

            # Give PACSCHNCB1M, PACSCHNCB2M, and PACTOMSCG1M the identical plastid haplotype,
            # while PACTOMSCG1M retains a divergent nuclear/ITS ribotype -> cytonuclear_discordant = 1
            rbcl_hap1 = rbcl_seed[:120] + ("A" if rbcl_seed[120] != "A" else "G") + rbcl_seed[121:]
            its_sch = its_seed
            its_tom = its_seed[:100] + ("C" if its_seed[100] != "C" else "T") + its_seed[101:]

            zip_path = raw_dir / "run_5563_fastqs.zip"
            with zipfile.ZipFile(zip_path, "w") as zf:
                for sid, p_seq, i_seq in (
                    ("PACSCHNCB1M", rbcl_hap1, its_sch),
                    ("PACSCHNCB2M", rbcl_hap1, its_sch),
                    ("PACTOMSCG1M", rbcl_hap1, its_tom),
                ):
                    reads = [(f"{sid}_p_{j}", p_seq) for j in range(5)] + [(f"{sid}_i_{j}", i_seq) for j in range(5)]
                    zf.writestr(f"fastqs/{sid}_R1.fastq.gz", _build_r1_fastq_gz_bytes(reads))

            skim_res = skim_mod.run_offtarget_organellar_skim(fastq_source=zip_path, genomics_dir=geno_dir)
            self.assertEqual(skim_res["status"], "skimmed")
            self.assertEqual([p for p in geno_dir.iterdir() if p.is_dir()], [], "genomics_dir must remain 100% flat")
            self.assertTrue((geno_dir / "skimmed_plastid_consensus.fasta").is_file())
            self.assertTrue((geno_dir / "skimmed_its_consensus.fasta").is_file())
            self.assertTrue((geno_dir / "cytonuclear_discordance.csv").is_file())

            cyto_rows = {}
            with open(geno_dir / "cytonuclear_discordance.csv", "r", encoding="utf-8", newline="") as fh:
                for row in csv.DictReader(fh):
                    cyto_rows[row["genetic_sample_id"]] = row
            self.assertEqual(cyto_rows["PACTOMSCG1M"]["cytonuclear_discordant"], "1")
            self.assertEqual(cyto_rows["PACSCHNCB1M"]["cytonuclear_discordant"], "0")

            # Invariant-site PHYLIP (732 bp = 6 x 122 bp contigs) coalescent delimitation
            rng = np.random.default_rng(7)
            bases = np.array(list("ACGT"))
            anc = rng.choice(bases, size=732)
            seq_b_base = anc.copy()
            for idx in range(0, 732, 14):
                seq_b_base[idx] = "G" if anc[idx] != "G" else "T"

            def _mutate(seq_arr: np.ndarray, n_mut: int, seed_off: int) -> str:
                arr = seq_arr.copy()
                r_local = np.random.default_rng(100 + seed_off)
                if n_mut > 0:
                    picks = r_local.choice(len(arr), size=n_mut, replace=False)
                    for p in picks:
                        arr[p] = "C" if arr[p] != "C" else "A"
                return "".join(arr.tolist())

            phy_records = [
                ("PACTOMSCG1M", _mutate(anc, 1, 1)),
                ("PACTOMSCG1P", _mutate(anc, 1, 2)),
                ("PACTOMVAL2M", _mutate(anc, 2, 3)),
                ("PACSCHNCB1M", _mutate(seq_b_base, 1, 4)),
                ("PACSCHNCB1P", _mutate(seq_b_base, 1, 5)),
                ("PACANOOHS1M", _mutate(anc, 0, 6)),
                ("PACANOOHS2M", _mutate(anc, 1, 7)),
            ]
            phy_path = geno_dir / "packera_invariant.phylip"
            phy_lines = [f" {len(phy_records)} {len(phy_records[0][1])}"] + [f"{nm:<16}{sq}" for nm, sq in phy_records]
            phy_path.write_text("\n".join(phy_lines) + "\n", encoding="utf-8")

            parsed_csv = tmp / "parsed_samples.csv"
            with open(parsed_csv, "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["genetic_sample_id", "species_standardized", "replicate_pair_id", "geno_quality_pass"])
                for nm, _ in phy_records:
                    sp = (
                        "Packera dubia"
                        if nm.startswith("PACTOM")
                        else ("Packera schweinitziana" if nm.startswith("PACSCH") else "Packera anonyma")
                    )
                    wr.writerow([nm, sp, nm[:-1], "1"])

            coal_res = coal_mod.run_coalescent_delimitation(
                genomics_dir=geno_dir,
                out_dir=out_dir,
                parsed_samples_csv=parsed_csv,
                block_len=122,
                n_boot=25,
            )
            self.assertEqual(coal_res["status"], "completed_alignment")
            for req in (
                "svdquartets_support.csv",
                "bpp_msci_loci.phy",
                "bpp_msci.ctl",
                "bpp_imap.txt",
                "coalescent_gdi_delimitation.csv",
                "coalescent_gdi_samples.csv",
            ):
                self.assertTrue((out_dir / req).is_file(), f"Missing coalescent output: {req}")

            gdi_map = {}
            with open(out_dir / "coalescent_gdi_delimitation.csv", "r", encoding="utf-8", newline="") as fh:
                for row in csv.DictReader(fh):
                    gdi_map[row["taxon"]] = row
            self.assertIn("Packera schweinitziana", gdi_map)
            self.assertGreater(float(gdi_map["Packera schweinitziana"]["coalescent_gdi"]), 0.70)
            self.assertIn("Packera anonyma", gdi_map)
            self.assertLess(float(gdi_map["Packera anonyma"]["coalescent_gdi"]), 0.25)

    def test_batched_r_step06_and_step07_phase2_consensus_and_triage(self) -> None:
        r_bin = shutil.which("Rscript")
        if not r_bin:
            self.skipTest("Rscript not found in PATH")

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            empty_geno = tmp / "empty_geno"
            pop_geno = tmp / "pop_geno"
            parsed_dir = tmp / "parsed"
            step06_empty_out = tmp / "step06_empty_out"
            step06_out = tmp / "step06_out"
            tables_dir = tmp / "tables"
            figs_dir = tmp / "figures"
            rep_dir = tmp / "reports"
            for d in (empty_geno, pop_geno, parsed_dir, step06_empty_out, step06_out, tables_dir, figs_dir, rep_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path = tmp / "packera_batch.sqlite"

            vouchers = [
                ("NCU00445952", 0, "Packera schweinitziana", "Packera schweinitziana", "Fuller", "NCU", "Tier_1_Gold", 35.82, -82.91, "NC", "Haywood", "2024-06-15", "Basal_Rosette", 6, 0, 1, 0.5, "Packera schweinitziana", 0.94, "Packera schweinitziana", 0.91, 166.0, 2.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.4, 0.85, 18.0, 0, -0.8, 0.2, 0.1),
                ("NCU_SCH_2", 0, "Packera schweinitziana", "Packera schweinitziana", "Fuller", "NCU", "Tier_1_Gold", 35.85, -82.88, "NC", "Haywood", "2024-06-18", "Basal_Rosette", 5, 0, 1, 0.4, "Packera schweinitziana", 0.92, "Packera schweinitziana", 0.88, 169.0, 4.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.5, 0.80, 17.5, 0, -0.7, 0.1, 0.0),
                ("NCU_DUB_1", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_1_Gold", 35.12, -79.22, "NC", "Moore", "2024-04-20", "Basal_Rosette", 7, 0, 1, 0.2, "Packera dubia", 0.95, "Packera dubia", 0.90, 110.0, -1.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.3, 0.90, 16.0, 0, 1.1, -0.1, 0.2),
                ("NCU_DUB_2", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_1_Gold", 35.14, -79.19, "NC", "Moore", "2024-04-22", "Basal_Rosette", 6, 0, 1, 0.3, "Packera dubia", 0.93, "Packera dubia", 0.89, 112.0, 1.5, "Measured", 1, 0, 1, "Basal_Rosette", 0.4, 0.88, 15.5, 0, 1.0, 0.0, 0.1),
                ("NCU_MIL_1", 0, "Packera millefolium", "Packera millefolium", "Fuller", "NCU", "Tier_1_Gold", 35.91, -82.11, "NC", "Mitchell", "2024-05-20", "Basal_Rosette", 5, 0, 1, 0.1, "Packera millefolium", 0.96, "Packera millefolium", 0.93, 140.0, 0.5, "Measured", 1, 0, 1, "Basal_Rosette", 0.2, 0.91, 19.0, 0, 0.1, 1.2, -0.3),
                ("NCU_MIL_2", 0, "Packera millefolium", "Packera millefolium", "Fuller", "NCU", "Tier_1_Gold", 35.93, -82.13, "NC", "Mitchell", "2024-05-22", "Basal_Rosette", 6, 0, 1, 0.2, "Packera millefolium", 0.95, "Packera millefolium", 0.91, 142.0, -0.5, "Measured", 1, 0, 1, "Basal_Rosette", 0.3, 0.89, 18.5, 0, 0.2, 1.1, -0.2),
                ("NCU_F1_HYB", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_2_Silver", 35.48, -81.05, "NC", "Burke", "2024-05-15", "Basal_Rosette", 4, 1, 1, 2.0, "Packera schweinitziana", 0.51, "Packera dubia", 0.55, 136.0, 12.0, "Measured", 1, 0, 1, "Basal_Rosette", 1.8, 0.35, 3.2, 0, 0.2, 0.1, 0.0),
                ("NCU_UNSEQ", 0, "Packera dubia", "Packera dubia", "Fuller", "NCU", "Tier_1_Gold", 35.10, -79.25, "NC", "Moore", "2024-04-18", "Basal_Rosette", 5, 0, 1, 0.1, "Packera dubia", 0.94, "Packera dubia", 0.92, 108.0, 0.0, "Measured", 1, 0, 1, "Basal_Rosette", 0.3, 0.85, 15.0, 0, 1.05, -0.05, 0.1),
            ]
            v_headers = [
                "catalogNumber", "plant_individual_id", "species_standardized", "species_raw",
                "recordedBy", "institutionCode", "determiner_tier", "decimalLatitude", "decimalLongitude",
                "stateProvince", "county", "eventDate", "leaf_zone", "basal_leaf_count",
                "is_morph_outlier", "gmm_cluster", "gmm_bayes_factor", "cda_predicted_taxon",
                "cda_posterior_prob", "edaphic_best_fit_taxon", "edaphic_prob", "hopkins_adjusted_doy",
                "delta_doy", "reproductive_trait_status", "is_anthesis", "is_passive", "is_basal_morph",
                "morph_tier_source", "type_edaphic_z_score", "open_set_novelty_p", "min_congener_mahal_d2",
                "cryptic_cluster_cohesion_pass", "PC1", "PC2", "PC3",
            ]
            vouchers_csv = tables_dir / "curated_vouchers.csv"
            with open(vouchers_csv, "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(v_headers)
                for row in vouchers:
                    wr.writerow(row)

            with open(tables_dir / "multimodal_conflict_flags.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["catalogNumber", "plant_individual_id", "edaphic_best_fit_taxon", "edaphic_prob", "type_edaphic_z_score", "delta_doy", "hopkins_adjusted_doy"])
                for row in vouchers:
                    wr.writerow([row[0], row[1], row[19], row[20], row[28], row[22], row[21]])

            with open(tables_dir / "morphometric_flags.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["catalogNumber", "plant_individual_id", "is_morph_outlier", "gmm_cluster", "gmm_bayes_factor", "cda_predicted_taxon", "cda_posterior_prob", "PC1", "PC2", "PC3"])
                for row in vouchers:
                    wr.writerow([row[0], row[1], row[14], row[15], row[16], row[17], row[18], row[32], row[33], row[34]])

            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                for row in vouchers:
                    conn.execute(
                        """
                        INSERT INTO vouchers (
                            catalogNumber, recordedBy, institutionCode, scientificName,
                            species_standardized, decimalLatitude, decimalLongitude,
                            stateProvince, county, eventDate, determiner_tier
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                        """,
                        (row[0], row[4], row[5], row[2], row[2], row[7], row[8], row[9], row[10], row[11], row[6]),
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, rosette_centroid_x) VALUES (?, ?, ?);",
                        (row[0], row[1], 100.0),
                    )
                conn.commit()
            finally:
                conn.close()

            s_defs = [
                ("NCU00445952", "NCU00445952", "G1_Direct_Voucher", "Packera schweinitziana", 0, 35.82, -82.91, "4x", "Clade_Sch", "ITS_Sch", 0, "SCH"),
                ("PACSCHNCB1M", "NCU_SCH_2", "G1_Direct_Voucher", "Packera schweinitziana", 0, 35.85, -82.88, "4x", "Clade_Sch", "ITS_Sch", 0, "SCH"),
                ("PACSCHNCB1P", "NCU_SCH_2", "G2_Population_Voucher_Proxy", "Packera schweinitziana", 0, 35.85, -82.88, "4x", "Clade_Sch", "ITS_Sch", 0, "SCH"),
                ("PACTOMSCG1M", "NCU_DUB_1", "G1_Direct_Voucher", "Packera dubia", 0, 35.12, -79.22, "2x", "Clade_Dub", "ITS_Dub", 0, "DUB"),
                ("PACTOMVAL2M", "NCU_DUB_2", "G1_Direct_Voucher", "Packera dubia", 0, 35.14, -79.19, "2x", "Clade_Dub", "ITS_Dub", 0, "DUB"),
                ("PACMILVAL1M", "NCU_MIL_1", "G1_Direct_Voucher", "Packera millefolium", 0, 35.91, -82.11, "4x", "Clade_Mil", "ITS_Mil", 0, "MIL"),
                ("PACMILMIT1M", "NCU_MIL_2", "G1_Direct_Voucher", "Packera millefolium", 0, 35.93, -82.13, "4x", "Clade_Mil", "ITS_Mil", 0, "MIL"),
                ("PACF1HYBM", "NCU_F1_HYB", "G1_Direct_Voucher", "Packera dubia", 1, 35.48, -81.05, "2x", "Clade_Sch", "ITS_Dub", 1, "F1_DUB_SCH"),
            ]
            loci = [f"L_{i:02d}" for i in range(30)]
            dos_rows = []
            for sid, _, _, _, _, _, _, _, _, _, _, mode in s_defs:
                vals = []
                for li in range(30):
                    if mode == "SCH":
                        v = 1.0 if li < 10 else 0.0
                    elif mode == "DUB":
                        v = 1.0 if 10 <= li < 20 else 0.0
                    elif mode == "MIL":
                        v = 1.0 if li >= 20 else 0.0
                    else:
                        v = 0.5 if li < 20 else 0.0
                    vals.append(f"{v:.4f}")
                dos_rows.append([sid] + vals)

            with open(parsed_dir / "allele_dosage_matrix.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["genetic_sample_id"] + loci)
                wr.writerows(dos_rows)

            with open(parsed_dir / "parsed_genomic_samples.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow([
                    "genetic_sample_id", "project_id", "catalogNumber", "linkage_tier",
                    "taxon_code", "taxon_raw", "species_standardized", "population_id",
                    "replicate_pair_id", "replicate_suffix", "is_putative_novel_or_hybrid_query",
                    "decimalLatitude", "decimalLongitude", "hq_reads", "mapped_reads",
                    "pct_mapped", "pct_repetitive", "n_missing_gt", "missing_genotype_rate",
                    "observed_het", "n_called_loci", "geno_quality_pass",
                    "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl",
                    "plastid_clade", "its_ribotype", "cytonuclear_discordant",
                ])
                for idx, (sid, cat, tier, sp, is_q, lat, lon, cyt, plast, its_r, cytonuc, _) in enumerate(s_defs, start=1):
                    wr.writerow([
                        sid, sid, cat, tier,
                        sid[3:6], sp, sp, f"POP_{idx:02d}" if not sid.endswith(("M", "P")) else f"POP_{sid[:-1]}",
                        sid[:-1] if sid[-1] in ("M", "P") else sid,
                        sid[-1] if sid[-1] in ("M", "P") else "",
                        str(is_q), str(lat), str(lon),
                        "450000", "410000", "91.1", "3.2", "0", "0.00", "0.22", "30", "1",
                        "30", cyt, "5.2", plast, its_r, str(cytonuc),
                    ])

            with open(parsed_dir / "genomic_manifest.json", "w", encoding="utf-8") as fh:
                json.dump({"status": "ingested", "n_samples_total": len(s_defs), "n_loci_retained": 30}, fh)
            (pop_geno / "packera_batch.vcf").write_text(
                _build_synthetic_vcf_text(
                    [r[0] for r in s_defs],
                    [("L_00", [("1/1" if r[11] == "SCH" else "0/0", 30, 10, 20) for r in s_defs])],
                ),
                encoding="utf-8",
            )
            (pop_geno / "RAxML_bestTree.batch").write_text(
                "(((NCU00445952:0.01,PACSCHNCB1M:0.01,PACSCHNCB1P:0.01):0.05,(PACTOMSCG1M:0.01,PACTOMVAL2M:0.01):0.05):0.02,PACF1HYBM:0.03,(PACMILVAL1M:0.01,PACMILMIT1M:0.01):0.08);\n",
                encoding="utf-8",
            )

            batch_r = tmp / "verify_batched_phase2.R"
            fig6a = figs_dir / "fig06a_genomic_pcoa_tree_concordance.png"
            fig6b = figs_dir / "fig06b_reticulation_gdi_landscape.png"
            fig7a_base = figs_dir / "fig07a_base_3stream.png"
            fig7a_4st = figs_dir / "fig07a_integrative_taxonomic_synthesis.png"
            q_base_csv = tables_dir / "triage_queue_base3.csv"
            q_4st_csv = tables_dir / "triage_queue.csv"
            rep_base_md = rep_dir / "report_base3.md"
            rep_4st_md = rep_dir / "Packera_Taxonomic_Revision_Summary.md"
            cplx_csv = rep_dir / "empirical_taxon_complexes.csv"
            aff_csv = rep_dir / "pairwise_taxon_affinity_matrix.csv"
            nbr_txt = rep_dir / "focal_auto_neighbors.txt"
            fig7b = figs_dir / "fig07b_taxon_complex_delimitation_network.png"
            html_out = rep_dir / "Triage_Dashboard_Interactive.html"

            r_script = f"""
Sys.setenv(
  PLANTMATICS_DB_PATH = {json.dumps(str(db_path))},
  PACKERA_SQLITE_DB = {json.dumps(str(db_path))},
  PLANTMATICS_GENOMICS_DIR = {json.dumps(str(empty_geno))},
  RUN_ID = "batch_smoke"
)

# 1. Standalone Step 06 empty no-op check via CLI
st_empty <- system2(
  {json.dumps(r_bin)},
  args = c(
    "scripts/analysis/06_genomic_phylogenomics.R",
    "--genomics-dir", {json.dumps(str(empty_geno))},
    "--parsed-dir", {json.dumps(str(tmp / "nonexistent_parsed"))},
    "--out-dir", {json.dumps(str(step06_empty_out))},
    "--tables-dir", {json.dumps(str(tmp / "empty_tables"))},
    "--output-figures", {json.dumps(str(tmp / "empty_figs"))},
    "--run-id", "batch_smoke_empty",
    "--db", {json.dumps(str(db_path))},
    "--demo"
  )
)
stopifnot(identical(as.integer(st_empty), 0L))
stopifnot(!file.exists(file.path({json.dumps(str(step06_empty_out))}, "genomic_assignments.csv")))

# 2. Standalone Step 06 populated Phase 2 execution via CLI
Sys.setenv(PLANTMATICS_GENOMICS_DIR = {json.dumps(str(pop_geno))})
st_pop <- system2(
  {json.dumps(r_bin)},
  args = c(
    "scripts/analysis/06_genomic_phylogenomics.R",
    "--genomics-dir", {json.dumps(str(pop_geno))},
    "--raw-genomics-dir", {json.dumps(str(pop_geno))},
    "--parsed-dir", {json.dumps(str(parsed_dir))},
    "--out-dir", {json.dumps(str(step06_out))},
    "--tables-dir", {json.dumps(str(tables_dir))},
    "--output-figures", {json.dumps(str(figs_dir))},
    "--figure-a", {json.dumps(str(fig6a))},
    "--figure-b", {json.dumps(str(fig6b))},
    "--run-id", "batch_smoke",
    "--db", {json.dumps(str(db_path))},
    "--demo"
  )
)
stopifnot(identical(as.integer(st_pop), 0L))

req_phase2_files <- c(
  "genomic_assignments.csv", "genomic_clade_cohesion_summary.csv",
  "genomic_model_metrics.csv", "genomic_multimodal_congruence_summary.csv",
  "coancestry_matrix_finerad.csv", "hybrid_index_q12_triangle.csv",
  "introgressive_fbranch_summary.csv", "quartet_concordance_factors_snps2cf.csv",
  "landscape_mmrr_feems_summary.csv", "feems_migration_edges.csv"
)
for (fn in req_phase2_files) {{
  fp <- file.path({json.dumps(str(step06_out))}, fn)
  if (!file.exists(fp)) stop(paste("Missing Step 06 Phase 2 file:", fp))
}}
stopifnot(file.exists({json.dumps(str(fig6a))}), file.info({json.dumps(str(fig6a))})$size > 5000)
stopifnot(file.exists({json.dumps(str(fig6b))}), file.info({json.dumps(str(fig6b))})$size > 5000)

ga <- utils::read.csv(file.path({json.dumps(str(step06_out))}, "genomic_assignments.csv"), stringsAsFactors = FALSE)
rownames(ga) <- ga$genetic_sample_id
stopifnot("NCU00445952" %in% rownames(ga))
stopifnot(identical(ga["NCU00445952", "geno_predicted_taxon"], "Packera schweinitziana"))
stopifnot(identical(ga["PACF1HYBM", "hybrid_generation_class"], "F1_Hybrid"))
stopifnot(is.finite(ga["PACF1HYBM", "interspecific_ancestry_q12"]) && ga["PACF1HYBM", "interspecific_ancestry_q12"] >= 0.80)
stopifnot(ga["PACMILVAL1M", "hybrid_generation_class"] != "F1_Hybrid")
stopifnot(is.na(ga["PACMILVAL1M", "interspecific_ancestry_q12"]) || ga["PACMILVAL1M", "interspecific_ancestry_q12"] < 0.25)

mmrr_df <- utils::read.csv(file.path({json.dumps(str(step06_out))}, "landscape_mmrr_feems_summary.csv"), stringsAsFactors = FALSE)
stopifnot(isTRUE(as.logical(mmrr_df$small_n_guardrail_active[1L])))
stopifnot(grepl("OutFLANK", mmrr_df$excluded_large_n_methods[1L]))
stopifnot(grepl("dadi", mmrr_df$excluded_large_n_methods[1L]))
stopifnot(grepl("HZAR", mmrr_df$excluded_large_n_methods[1L]))

# 3. Source Step 07 synthesis in this same R process to verify 3-stream vs 4-stream invariance + Phase 2 rules
source("scripts/analysis/07_triage_dashboard_synthesis.R")

opts_3 <- parse_args_robust(c(
  "--vouchers", {json.dumps(str(vouchers_csv))},
  "--morphometrics", {json.dumps(str(tables_dir / "morphometric_flags.csv"))},
  "--multimodal-flags", {json.dumps(str(tables_dir / "multimodal_conflict_flags.csv"))},
  "--output-queue", {json.dumps(str(q_base_csv))},
  "--output-plot", {json.dumps(str(fig7a_base))},
  "--output-report", {json.dumps(str(rep_base_md))},
  "--genomic-dir", "",
  "--demo"
))
q3 <- build_triage_queue(opts_3)

opts_4 <- parse_args_robust(c(
  "--vouchers", {json.dumps(str(vouchers_csv))},
  "--morphometrics", {json.dumps(str(tables_dir / "morphometric_flags.csv"))},
  "--multimodal-flags", {json.dumps(str(tables_dir / "multimodal_conflict_flags.csv"))},
  "--output-queue", {json.dumps(str(q_4st_csv))},
  "--output-plot", {json.dumps(str(fig7a_4st))},
  "--output-report", {json.dumps(str(rep_4st_md))},
  "--output-html", {json.dumps(str(html_out))},
  "--output-complexes", {json.dumps(str(cplx_csv))},
  "--output-affinity", {json.dumps(str(aff_csv))},
  "--output-neighbors", {json.dumps(str(nbr_txt))},
  "--output-complex-plot", {json.dumps(str(fig7b))},
  "--genomic-dir", {json.dumps(str(step06_out))},
  "--demo"
))
q4 <- build_triage_queue(opts_4)
render_synthesis_figure_png(q4, opts_4$output_plot, output_complex_plot = opts_4$output_complex_plot)
generate_taxonomic_revision_report(q4, opts_4)

u3 <- q3[q3$catalogNumber == "NCU_UNSEQ", , drop = FALSE]
u4 <- q4[q4$catalogNumber == "NCU_UNSEQ", , drop = FALSE]
stopifnot(nrow(u3) == 1L, nrow(u4) == 1L)
stopifnot(identical(as.character(u3$taxonomic_treatment_subtype), as.character(u4$taxonomic_treatment_subtype)))
stopifnot(abs(as.numeric(u3$probabilistic_consensus_score) - as.numeric(u4$probabilistic_consensus_score)) < 1e-6)

f1_row <- q4[q4$catalogNumber == "NCU_F1_HYB", , drop = FALSE]
stopifnot(nrow(f1_row) == 1L)
stopifnot(identical(as.character(f1_row$taxonomic_treatment_subtype), "Hybrid_Introgressant"))

stopifnot(file.exists(file.path({json.dumps(str(tables_dir))}, "genomic_field_triage.csv")))
stopifnot(file.exists({json.dumps(str(fig7a_4st))}), file.info({json.dumps(str(fig7a_4st))})$size > 5000)
rep_txt <- paste(readLines({json.dumps(str(rep_4st_md))}, warn = FALSE), collapse = "\\n")
stopifnot(grepl("## 5\\\\. Independent Step 06 Genomic & Phylogenomic Evidence Stream", rep_txt))
stopifnot(length(list.files({json.dumps(str(figs_dir))}, pattern = "\\\\.pdf$", ignore.case = TRUE)) == 0L)

# 4. Direct Phase 2 apply_taxonomic_decision_matrix tests
make_intra_probe <- function(prefix, gdi_val, hyb_cls, q12_val) {{
  df_k2 <- data.frame(
    catalogNumber = paste0(prefix, "_", 1:5),
    plant_individual_id = 0L,
    species_standardized = rep("Packera dubia", 5L),
    species_raw = rep("Packera dubia", 5L),
    recordedBy = c("Fuller", "Weakley", "Sorrie", "Poindexter", "Krings"),
    institutionCode = c("NCU", "DUKE", "US", "NY", "GH"),
    determiner_tier = rep("Tier_1_Gold", 5L),
    decimalLatitude = c(35.10, 35.45, 35.80, 36.15, 36.50),
    decimalLongitude = c(-79.10, -79.55, -80.05, -80.55, -81.05),
    stateProvince = rep("NC", 5L),
    county = c("Moore", "Chatham", "Orange", "Ashe", "Watauga"),
    eventDate = c("2020-04-15", "2021-04-18", "2022-04-20", "2023-04-22", "2024-04-25"),
    leaf_zone = rep("Basal_Rosette", 5L),
    basal_leaf_count = rep(6L, 5L),
    is_basal_morph = rep(TRUE, 5L),
    is_passive = rep(FALSE, 5L),
    morph_tier_source = rep("Basal_Rosette", 5L),
    is_morph_outlier = rep(TRUE, 5L),
    gmm_cluster = rep(2L, 5L),
    gmm_bayes_factor = rep(12.0, 5L),
    cda_predicted_taxon = rep("Packera dubia", 5L),
    cda_posterior_prob = rep(0.45, 5L),
    edaphic_best_fit_taxon = rep("Packera dubia", 5L),
    edaphic_prob = rep(0.20, 5L),
    type_edaphic_z_score = rep(3.2, 5L),
    hopkins_adjusted_doy = rep(155.0, 5L),
    delta_doy = rep(35.0, 5L),
    reproductive_trait_status = rep("Measured", 5L),
    is_anthesis = rep(TRUE, 5L),
    open_set_novelty_p = rep(0.01, 5L),
    open_set_df = rep(5L, 5L),
    min_congener_mahal_d2 = rep(28.0, 5L),
    cryptic_cluster_cohesion_pass = rep(TRUE, 5L),
    PC1 = c(2.4, 2.5, 2.3, 2.6, 2.4),
    PC2 = c(-1.1, -1.0, -1.2, -1.1, -0.9),
    PC3 = c(0.4, 0.3, 0.5, 0.4, 0.4),
    geno_predicted_taxon = rep("Packera dubia", 5L),
    geno_linkage_tier = rep("G1_Direct_Voucher", 5L),
    geno_quality_pass = rep(TRUE, 5L),
    geno_posterior_prob = rep(0.96, 5L),
    coalescent_gdi = rep(gdi_val, 5L),
    hybrid_generation_class = rep(hyb_cls, 5L),
    interspecific_ancestry_q12 = rep(q12_val, 5L),
    is_geno_admixed = rep(FALSE, 5L),
    stringsAsFactors = FALSE
  )
  df_k1 <- data.frame(
    catalogNumber = paste0(prefix, "_BASE_", 1:3),
    plant_individual_id = 0L,
    species_standardized = rep("Packera dubia", 3L),
    species_raw = rep("Packera dubia", 3L),
    recordedBy = c("Radford", "Bozeman", "Leonard"),
    institutionCode = c("NCU", "NCU", "DUKE"),
    determiner_tier = rep("Tier_1_Gold", 3L),
    decimalLatitude = c(35.05, 35.12, 35.18),
    decimalLongitude = c(-79.05, -79.12, -79.18),
    stateProvince = rep("NC", 3L),
    county = c("Moore", "Hoke", "Richmond"),
    eventDate = c("1960-04-10", "1968-04-12", "1972-04-14"),
    leaf_zone = rep("Basal_Rosette", 3L),
    basal_leaf_count = rep(6L, 3L),
    is_basal_morph = rep(TRUE, 3L),
    is_passive = rep(FALSE, 3L),
    morph_tier_source = rep("Basal_Rosette", 3L),
    is_morph_outlier = rep(FALSE, 3L),
    gmm_cluster = rep(1L, 3L),
    gmm_bayes_factor = rep(0.2, 3L),
    cda_predicted_taxon = rep("Packera dubia", 3L),
    cda_posterior_prob = rep(0.92, 3L),
    edaphic_best_fit_taxon = rep("Packera dubia", 3L),
    edaphic_prob = rep(0.88, 3L),
    type_edaphic_z_score = rep(0.4, 3L),
    hopkins_adjusted_doy = rep(115.0, 3L),
    delta_doy = rep(0.0, 3L),
    reproductive_trait_status = rep("Measured", 3L),
    is_anthesis = rep(TRUE, 3L),
    open_set_novelty_p = rep(0.85, 3L),
    open_set_df = rep(5L, 3L),
    min_congener_mahal_d2 = rep(16.0, 3L),
    cryptic_cluster_cohesion_pass = rep(FALSE, 3L),
    PC1 = c(-0.5, -0.4, -0.6),
    PC2 = c(0.1, 0.0, 0.2),
    PC3 = c(0.0, 0.1, -0.1),
    geno_predicted_taxon = rep("Packera dubia", 3L),
    geno_linkage_tier = rep("G1_Direct_Voucher", 3L),
    geno_quality_pass = rep(TRUE, 3L),
    geno_posterior_prob = rep(0.96, 3L),
    coalescent_gdi = rep(gdi_val, 3L),
    hybrid_generation_class = rep("Pure_Parent_1", 3L),
    interspecific_ancestry_q12 = rep(0.02, 3L),
    is_geno_admixed = rep(FALSE, 3L),
    stringsAsFactors = FALSE
  )
  rbind(df_k2, df_k1)
}}

intra_oversplit <- apply_taxonomic_decision_matrix(
  make_intra_probe("OVERSPLIT", 0.12, "Pure_Parent_1", 0.02),
  intraspecific_mode = TRUE,
  gmm_two_delta_bic = 18.0
)
stopifnot(grepl("Coalescent Oversplit Block", intra_oversplit$synthesis_rationale[1L]))
stopifnot(intra_oversplit$taxonomic_treatment_subtype[1L] %in% c("Ecophenotypic_Plasticity", "Confirmed_Species"))

intra_cryptic <- apply_taxonomic_decision_matrix(
  make_intra_probe("CRYPTIC", 0.88, "Divergent_NonHybrid_Lineage", 0.04),
  intraspecific_mode = TRUE,
  gmm_two_delta_bic = 18.0
)
stopifnot(identical(as.character(intra_cryptic$taxonomic_treatment_subtype[1L]), "Putative_Cryptic_Species"))
cat("BATCHED_PHASE2_R_OK\\n")
"""
            batch_r.write_text(r_script, encoding="utf-8")

            proc = subprocess.run(
                [r_bin, str(batch_r)],
                cwd=PROJECT_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                proc.returncode,
                0,
                f"Batched R verification failed (exit={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
            )
            self.assertIn("BATCHED_PHASE2_R_OK", proc.stdout)

    def test_download_public_packera_sequences_flat_extraction_and_sha256(self) -> None:
        """(a) Verify mocked Dryad (SHA-256 + nested zip), NCBI nuccore/SRA, and BOLD download flatly."""

        class _MockHTTPResponse(io.BytesIO):
            def __enter__(self) -> "_MockHTTPResponse":
                return self

            def __exit__(self, *args: object) -> bool:
                self.close()
                return False

        def _zip_bytes(entries: dict[str, str]) -> bytes:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
                for arcname, content in entries.items():
                    zf.writestr(arcname, content.encode("utf-8"))
            return buf.getvalue()

        tree_zip = _zip_bytes(
            {
                "nested/sub/trees/MO_Astral_lpp.tre": "((Packera_dubia:0.05,Packera_schweinitziana:0.05):0.10,Packera_aurea:0.20);\n",
                "nested/sub/trees/ASTRALPro_Astral_lpp.tre": "((Packera_dubia:0.06,Packera_schweinitziana:0.06):0.10,Packera_aurea:0.20);\n",
            }
        )
        phytop_zip = _zip_bytes(
            {
                "nested/phytop/MO_PhyTop.tsv": "taxon\tils\tih\nPackera dubia\t0.18\t0.06\nPackera schweinitziana\t0.14\t0.04\n",
            }
        )
        align_zip = _zip_bytes(
            {
                "nested/plastomes/18complete_Packera_plastomes.fasta": (
                    ">MK001|Packera_dubia\n" + ("ATGCTAGCTGATCGATCGTA" * 7) + "\n"
                ),
            }
        )
        zip_payloads = {
            "101": ("speciestree_files.zip", tree_zip),
            "102": ("PhyTop_results.zip", phytop_zip),
            "103": ("alignment_files.zip", align_zip),
        }

        versions_json = json.dumps(
            {
                "_embedded": {
                    "stash:versions": [
                        {"versionNumber": 1, "_links": {"self": {"href": "/api/v2/versions/999"}}}
                    ]
                }
            }
        ).encode("utf-8")
        files_json = json.dumps(
            {
                "_embedded": {
                    "stash:files": [
                        {
                            "id": int(fid),
                            "path": fname,
                            "size": len(zbytes),
                            "digestType": "sha-256",
                            "digest": hashlib.sha256(zbytes).hexdigest(),
                            "_links": {"stash:download": {"href": f"/api/v2/files/{fid}/download"}},
                        }
                        for fid, (fname, zbytes) in zip_payloads.items()
                    ]
                }
            }
        ).encode("utf-8")
        ncbi_nuccore_fasta = (
            ">NC_0001 Senecio tomentosus chloroplast, complete genome\n"
            + ("ACGTACGTACGTATCG" * 8)
            + "\n"
            ">MT_0002 Packera aurea internal transcribed spacer 1\n"
            + ("TGCATGCATGCATAGC" * 8)
            + "\n"
            ">JX_0003 Senecio tomentosus granule-bound starch synthase (waxy) gene\n"
            + ("GGCCAATTAAGGCCTT" * 8)
            + "\n"
        ).encode("utf-8")
        ena_tsv = (
            "run_accession\tsample_accession\tscientific_name\tlibrary_name\tfastq_ftp\tfastq_bytes\n"
            "SRR22900001\tSAMN01\tPackera dubia\tLIB1\tftp.sra.ebi.ac.uk/vol1/fastq/SRR229/001.fastq.gz\t1024000\n"
        ).encode("utf-8")
        sra_runinfo_csv = (
            "Run,BioProject,Sample,ScientificName,LibraryName,download_path,size_MB\n"
            "SRR22900002,PRJNA978568,SAMN02,Packera schweinitziana,LIB2,sra-pub-run-odp.s3.amazonaws.com/SRR22900002,12\n"
        ).encode("utf-8")
        bold_fasta = ">BOLD_01|Packera schweinitziana|ITS2\n" + ("ACGTGCATGCAT" * 11) + "\n"

        def _mock_open_url(req_or_url: object, data: object = None, timeout: object = 60) -> _MockHTTPResponse:
            _ = (data, timeout)
            url = req_or_url.full_url if hasattr(req_or_url, "full_url") else str(req_or_url)
            if "/datasets/" in url and url.endswith("/versions"):
                return _MockHTTPResponse(versions_json)
            if "/versions/999/files" in url:
                return _MockHTTPResponse(files_json)
            for fid, (_, zbytes) in zip_payloads.items():
                if f"/api/v2/files/{fid}/download" in url:
                    return _MockHTTPResponse(zbytes)
            if "esearch.fcgi" in url and "db=nuccore" in url:
                return _MockHTTPResponse(
                    json.dumps({"esearchresult": {"count": "3", "idlist": ["101", "102", "103"]}}).encode("utf-8")
                )
            if "efetch.fcgi" in url and "db=nuccore" in url:
                return _MockHTTPResponse(ncbi_nuccore_fasta)
            if "ebi.ac.uk/ena/portal/api/filereport" in url:
                return _MockHTTPResponse(ena_tsv if "PRJNA907383" in url else b"")
            if "esearch.fcgi" in url and "db=sra" in url:
                if "PRJNA978568" in url:
                    return _MockHTTPResponse(json.dumps({"esearchresult": {"idlist": ["502"]}}).encode("utf-8"))
                return _MockHTTPResponse(json.dumps({"esearchresult": {"idlist": []}}).encode("utf-8"))
            if "efetch.fcgi" in url and "db=sra" in url:
                return _MockHTTPResponse(sra_runinfo_csv)
            if "portal.boldsystems.org/api/query" in url:
                return _MockHTTPResponse(bold_fasta.encode("utf-8"))
            raise AssertionError(f"Unexpected unmocked network URL in offline test: {url}")

        def _mock_opener_open(self_opener: object, fullurl: object, data: object = None, timeout: object = 60) -> _MockHTTPResponse:
            _ = self_opener
            return _mock_open_url(fullurl, data=data, timeout=timeout)

        with tempfile.TemporaryDirectory() as td:
            geno_dir = Path(td) / "genomics"
            env_patch = {"DRYAD_API_TOKEN": "test_dryad_token", "NCBI_API_KEY": "test_ncbi_key"}
            with (
                mock.patch.dict(os.environ, env_patch),
                mock.patch("urllib.request.urlopen", side_effect=_mock_open_url),
                mock.patch("urllib.request.OpenerDirector.open", new=_mock_opener_open),
            ):
                summary = PUBLIC_DL_MOD.download_public_packera_sequences(genomics_dir=geno_dir)

            self.assertEqual(summary["dryad_status"], "downloaded")
            self.assertEqual([p for p in geno_dir.iterdir() if p.is_dir()], [], "genomics_dir must have 0 subdirectories")
            expected_files = {
                "public_packera_MO_Astral_lpp.tre",
                "public_packera_ASTRALPro_Astral_lpp.tre",
                "public_packera_PhyTop_MO.tsv",
                "public_packera_18complete_plastomes.fasta",
                "public_packera_plastomes_and_cpdna.fasta",
                "public_packera_its_ets.fasta",
                "public_packera_waxy_nuclear.fasta",
                "public_packera_sra_runs.csv",
                "public_packera_bold_barcodes.fasta",
                "public_packera_dryad_manifest.json",
                "public_packera_sources_manifest.json",
            }
            actual_files = {p.name for p in geno_dir.iterdir() if p.is_file()}
            self.assertTrue(expected_files.issubset(actual_files), f"Missing downloaded files: {expected_files - actual_files}")
            self.assertIn("Packera_dubia", (geno_dir / "public_packera_plastomes_and_cpdna.fasta").read_text(encoding="utf-8"))
            self.assertIn("Packera_dubia", (geno_dir / "public_packera_waxy_nuclear.fasta").read_text(encoding="utf-8"))
            sra_txt = (geno_dir / "public_packera_sra_runs.csv").read_text(encoding="utf-8")
            self.assertIn("SRR22900001", sra_txt)
            self.assertIn("SRR22900002", sra_txt)

    def test_public_reference_files_only_stay_noop_no_payload(self) -> None:
        """(b) Verify public_packera_* references alone keep has_genetic_payload False and Step 06 ingestion skipped_empty."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            out_dir = tmp / "parsed_out"
            geno_dir.mkdir(parents=True, exist_ok=True)

            pub_names = (
                "public_packera_MO_Astral_lpp.tre",
                "public_packera_MO_concat_best.fasta",
                "public_packera_18complete_plastomes.fasta",
                "public_packera_sra_runs.csv",
                "public_packera_dryad_manifest.json",
                "public_packera_bold_barcodes.fasta",
            )
            for fname in pub_names:
                (geno_dir / fname).write_text("placeholder\n", encoding="utf-8")

            assets = INGEST_MOD.discover_genomic_inputs(geno_dir)
            self.assertFalse(assets.has_genetic_payload)
            self.assertEqual({p.name for p in assets.public_reference_files}, set(pub_names))

            res = INGEST_MOD.run_genomic_ingestion(genomics_dir=geno_dir, out_dir=out_dir)
            self.assertEqual(res["status"], "skipped_empty")
            self.assertEqual(len(list(geno_dir.iterdir())), len(pub_names))

    def test_skim_offtarget_organellar_baits_public_refs_including_waxy(self) -> None:
        """(c) Verify skim_offtarget_organellar.py baits public plastid, ITS, and waxy FASTA refs and records matches."""

        def _seq(seed_int: int, length: int = 140) -> str:
            rng = np.random.default_rng(seed_int)
            return "".join(rng.choice(list("ACGT"), size=length).tolist())

        cp_ser, cp_aur = _seq(202601), _seq(202602)
        its_ser, its_aur = _seq(202603), _seq(202604)
        waxy_ser, waxy_aur = _seq(202605), _seq(202606)

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            raw_dir = tmp / "raw_reads"
            geno_dir.mkdir(parents=True, exist_ok=True)
            raw_dir.mkdir(parents=True, exist_ok=True)

            (geno_dir / "public_packera_18complete_plastomes.fasta").write_text(
                f">MK_CP_SER Packera serpenticola\n{cp_ser}\n>MK_CP_AUR Packera aurea\n{cp_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_its_ets.fasta").write_text(
                f">MT_ITS_SER Packera serpenticola\n{its_ser}\n>MT_ITS_AUR Packera aurea\n{its_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_waxy_nuclear.fasta").write_text(
                f">JX_WAXY_SER Packera serpenticola\n{waxy_ser}\n>JX_WAXY_AUR Packera aurea\n{waxy_aur}\n",
                encoding="utf-8",
            )

            zip_path = raw_dir / "run_5563_public_skim.zip"
            with zipfile.ZipFile(zip_path, "w") as zf:
                for sid, p_s, i_s, w_s in (
                    ("PACSERNCS1M", cp_ser, its_ser, waxy_ser),
                    ("PACAURFK1M", cp_aur, its_aur, waxy_aur),
                ):
                    reads = (
                        [(f"{sid}_cp_{j}", p_s) for j in range(4)]
                        + [(f"{sid}_its_{j}", i_s) for j in range(4)]
                        + [(f"{sid}_waxy_{j}", w_s) for j in range(4)]
                    )
                    zf.writestr(f"fastqs/{sid}_R1.fastq.gz", _build_r1_fastq_gz_bytes(reads))

            res = skim_mod.run_offtarget_organellar_skim(fastq_source=zip_path, genomics_dir=geno_dir)
            self.assertEqual(res["status"], "skimmed")
            self.assertEqual([p for p in geno_dir.iterdir() if p.is_dir()], [])
            for out_fa in (
                "skimmed_plastid_consensus.fasta",
                "skimmed_its_consensus.fasta",
                "skimmed_waxy_consensus.fasta",
            ):
                p = geno_dir / out_fa
                self.assertTrue(p.is_file(), f"Missing skimmed consensus FASTA: {out_fa}")
                txt = p.read_text(encoding="utf-8")
                self.assertIn(">PACSERNCS1M|", txt)
                self.assertIn(">PACAURFK1M|", txt)

            cyto_rows: dict[str, dict[str, str]] = {}
            with open(geno_dir / "cytonuclear_discordance.csv", "r", encoding="utf-8", newline="") as fh:
                for row in csv.DictReader(fh):
                    cyto_rows[row["genetic_sample_id"]] = row
            self.assertEqual(cyto_rows["PACSERNCS1M"]["public_plastome_match"], "Packera serpenticola")
            self.assertEqual(cyto_rows["PACSERNCS1M"]["public_its_match"], "Packera serpenticola")
            self.assertEqual(cyto_rows["PACSERNCS1M"]["waxy_lineage_match"], "Packera serpenticola")
            self.assertEqual(cyto_rows["PACAURFK1M"]["public_plastome_match"], "Packera aurea")
            self.assertEqual(cyto_rows["PACAURFK1M"]["public_its_match"], "Packera aurea")
            self.assertEqual(cyto_rows["PACAURFK1M"]["waxy_lineage_match"], "Packera aurea")

    def test_step06_public_backbone_outgroup_and_sqlite_columns(self) -> None:
        """(d) Verify Step 06 prevents focal tree shadowing, uses public backbone outgroup polarization, and persists 6 SQLite cols."""
        if not rscript_available():
            self.skipTest("Rscript unavailable or PACKERA_FAST_AUDIT=1")

        def _seq(seed_int: int, length: int = 120) -> str:
            rng = np.random.default_rng(seed_int)
            return "".join(rng.choice(list("ACGT"), size=length).tolist())

        cp_dub, cp_aur = _seq(303101), _seq(303102)
        its_dub, its_aur = _seq(303103), _seq(303104)
        waxy_dub, waxy_aur = _seq(303105), _seq(303106)

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            parsed_dir = tmp / "parsed"
            out_dir = tmp / "step06_out"
            tables_dir = tmp / "tables"
            figs_dir = tmp / "figures"
            (tmp / "data").mkdir(parents=True, exist_ok=True)
            for d in (geno_dir, parsed_dir, out_dir, tables_dir, figs_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path = tmp / "data" / "plantmatics.sqlite"

            sample_defs = [
                ("PACTOMSCG1M", "NCU_DUB_1", "Packera dubia", 35.12, -79.22, [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0.5, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0]),
                ("PACTOMVAL2M", "NCU_DUB_2", "Packera dubia", 35.14, -79.19, [1, 1, 1, 0.5, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0]),
                ("PACSCHNCB1M", "NCU_SCH_1", "Packera schweinitziana", 35.82, -82.91, [0.5, 0, 0, 0, 1, 1, 1, 1, 1, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                ("PACSCHNCB2M", "NCU_SCH_2", "Packera schweinitziana", 35.85, -82.88, [0, 0.5, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                ("PACMILVAL1M", "NCU_MIL_1", "Packera millefolium", 35.91, -82.11, [0, 0, 0, 0, 0.5, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0.5, 0, 0, 0, 0, 0, 0, 0]),
                ("PACMILMIT1M", "NCU_MIL_2", "Packera millefolium", 35.93, -82.13, [0, 0, 0, 0, 0, 0.5, 0, 0, 1, 0.5, 0, 0, 1, 1, 1, 1, 0, 0.5, 0, 0, 0, 0, 0, 0]),
                ("PACAURFK1M", "NCU_AUR_1", "Packera aurea", 35.60, -81.90, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0.5]),
                ("PACAURFK2M", "NCU_AUR_2", "Packera aurea", 35.63, -81.87, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0.5, 1]),
            ]

            run_id = "run_5563"
            conn = get_db_connection(db_path=db_path, init_schema=True, workspace_root=tmp)
            try:
                conn.execute("INSERT INTO pipeline_runs (run_id, run_dir, intraspecific_mode) VALUES (?, ?, 0);", (run_id, str(tmp)))
                for sid, cat, sp, lat, lon, _ in sample_defs:
                    conn.execute(
                        """
                        INSERT INTO vouchers (
                            catalogNumber, scientificName, species_standardized,
                            decimalLatitude, decimalLongitude, determiner_tier
                        ) VALUES (?, ?, ?, ?, ?, 'Tier_1_Gold');
                        """,
                        (cat, sp, sp, lat, lon),
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id, rosette_centroid_x) VALUES (?, 0, 100.0);",
                        (cat,),
                    )
                    conn.execute(
                        """
                        INSERT INTO run_active_cohort (
                            run_id, catalogNumber, plant_individual_id,
                            species_standardized, determiner_tier, authority_weight
                        ) VALUES (?, ?, 0, ?, 'Tier_1_Gold', 2.0);
                        """,
                        (run_id, cat, sp),
                    )
                    pub_sp = sp if sid in ("PACTOMSCG1M", "PACAURFK1M") else None
                    conn.execute(
                        """
                        INSERT INTO genomic_samples (
                            genetic_sample_id, catalogNumber, linkage_tier, species_standardized,
                            geno_quality_pass, mapped_reads, public_plastome_match,
                            public_its_match, waxy_lineage_match
                        ) VALUES (?, ?, 'G1_Direct_Voucher', ?, 1, 450000, ?, ?, ?);
                        """,
                        (sid, cat, sp, pub_sp, pub_sp, pub_sp),
                    )
                conn.commit()
            finally:
                conn.close()

            loci = [f"L_{i:02d}" for i in range(24)]
            with open(parsed_dir / "allele_dosage_matrix.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["genetic_sample_id"] + loci)
                for sid, _, _, _, _, vals in sample_defs:
                    wr.writerow([sid] + [f"{v:.4f}" for v in vals])

            with open(parsed_dir / "parsed_genomic_samples.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow([
                    "genetic_sample_id", "project_id", "catalogNumber", "linkage_tier",
                    "taxon_code", "taxon_raw", "species_standardized", "population_id",
                    "replicate_pair_id", "replicate_suffix", "is_putative_novel_or_hybrid_query",
                    "decimalLatitude", "decimalLongitude", "hq_reads", "mapped_reads",
                    "pct_mapped", "pct_repetitive", "n_missing_gt", "missing_genotype_rate",
                    "observed_het", "n_called_loci", "geno_quality_pass",
                    "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl",
                    "plastid_clade", "its_ribotype", "cytonuclear_discordant",
                ])
                for idx, (sid, cat, sp, lat, lon, _) in enumerate(sample_defs, start=1):
                    wr.writerow([
                        sid, sid, cat, "G1_Direct_Voucher", sid[3:6], sp, sp, f"POP_{idx:02d}",
                        sid[:-1], "M", "0", str(lat), str(lon), "450000", "410000",
                        "91.1", "3.0", "0", "0.00", "0.20", "24", "1", "24", "2x", "1.1",
                        f"cp_{sid[3:6]}", f"ITS_{sid[3:6]}", "0",
                    ])

            raxml_tree = geno_dir / "RAxML_bestTree.test"
            raxml_tree.write_text(
                "(((PACTOMSCG1M:0.01,PACTOMVAL2M:0.01):0.04,(PACSCHNCB1M:0.01,PACSCHNCB2M:0.01):0.04):0.05,(PACMILVAL1M:0.01,PACMILMIT1M:0.01):0.10,(PACAURFK1M:0.01,PACAURFK2M:0.01):0.25);\n",
                encoding="utf-8",
            )
            mo_tree = geno_dir / "public_packera_MO_Astral_lpp.tre"
            mo_tree.write_text(
                "(((Packera_dubia:0.04,Packera_schweinitziana:0.04):0.06,Packera_millefolium:0.10):0.15,Packera_aurea:0.25);\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_ASTRALPro_Astral_lpp.tre").write_text(
                "(((Packera_dubia:0.05,Packera_schweinitziana:0.05):0.06,Packera_millefolium:0.11):0.14,Packera_aurea:0.25);\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_PhyTop_MO.tsv").write_text(
                "taxon\tils_index\tih_index\n"
                "Packera dubia\t0.21\t0.08\n"
                "Packera schweinitziana\t0.19\t0.05\n"
                "Packera millefolium\t0.14\t0.03\n"
                "Packera aurea\t0.11\t0.02\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_18complete_plastomes.fasta").write_text(
                f">MK_CP_DUB Packera dubia\n{cp_dub}\n>MK_CP_AUR Packera aurea\n{cp_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_its_ets.fasta").write_text(
                f">MT_ITS_DUB Packera dubia\n{its_dub}\n>MT_ITS_AUR Packera aurea\n{its_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_waxy_nuclear.fasta").write_text(
                f">JX_WAXY_DUB Packera dubia\n{waxy_dub}\n>JX_WAXY_AUR Packera aurea\n{waxy_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "skimmed_plastid_consensus.fasta").write_text(
                f">PACTOMSCG1M|clade=cp_1\n{cp_dub}\n>PACAURFK1M|clade=cp_2\n{cp_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "skimmed_its_consensus.fasta").write_text(
                f">PACTOMSCG1M|ribotype=ITS_1\n{its_dub}\n>PACAURFK1M|ribotype=ITS_2\n{its_aur}\n",
                encoding="utf-8",
            )
            (geno_dir / "skimmed_waxy_consensus.fasta").write_text(
                f">PACTOMSCG1M|waxy=Waxy_1\n{waxy_dub}\n>PACAURFK1M|waxy=Waxy_2\n{waxy_aur}\n",
                encoding="utf-8",
            )
            # Place public_packera_MO_Astral_lpp.tre first so filtering proves focal RAxML_bestTree.test is not shadowed
            with open(parsed_dir / "genomic_manifest.json", "w", encoding="utf-8") as fh:
                json.dump(
                    {
                        "status": "ingested",
                        "n_samples_total": len(sample_defs),
                        "n_loci_retained": 24,
                        "tree_files": [str(mo_tree), str(raxml_tree)],
                    },
                    fh,
                )

            proc = run_warm_r_cmd(
                [
                    "Rscript",
                    str(PROJECT_ROOT / "scripts" / "analysis" / "06_genomic_phylogenomics.R"),
                    "--genomics-dir",
                    str(geno_dir),
                    "--raw-genomics-dir",
                    str(geno_dir),
                    "--parsed-dir",
                    str(parsed_dir),
                    "--out-dir",
                    str(out_dir),
                    "--tables-dir",
                    str(tables_dir),
                    "--figure-a",
                    str(figs_dir / "fig06a.png"),
                    "--figure-b",
                    str(figs_dir / "fig06b.png"),
                    "--run-id",
                    run_id,
                    "--db",
                    str(db_path),
                ],
                env={
                    "PLANTMATICS_DB_PATH": str(db_path),
                    "PLANTMATICS_GENOMICS_DIR": str(geno_dir),
                    "PYTHONPATH": str(PROJECT_ROOT),
                },
            )
            self.assertEqual(proc.returncode, 0, f"Step 06 public backbone run failed:\n{proc.stdout}\n{proc.stderr}")

            cohesion_csv = out_dir / "genomic_clade_cohesion_summary.csv"
            self.assertTrue(cohesion_csv.is_file())
            with open(cohesion_csv, "r", encoding="utf-8", newline="") as fh:
                cohesion_rows = list(csv.DictReader(fh))
            self.assertGreaterEqual(len(cohesion_rows), 4, "Focal RAxML_bestTree.test must drive clade cohesion")

            with open(out_dir / "genomic_assignments.csv", "r", encoding="utf-8", newline="") as fh:
                assign_rows = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            self.assertTrue(bool(assign_rows["PACTOMSCG1M"].get("nearest_tree_sister_taxon")))

            fb_csv = out_dir / "introgressive_fbranch_summary.csv"
            self.assertTrue(fb_csv.is_file())
            with open(fb_csv, "r", encoding="utf-8", newline="") as fh:
                fb_rows = list(csv.DictReader(fh))
            self.assertGreater(len(fb_rows), 0)
            self.assertTrue(all(r["polarization_source"] == "public_nuclear_backbone" for r in fb_rows))

            for req_csv in ("public_backbone_concordance.csv", "public_organellar_placement.csv"):
                self.assertTrue((out_dir / req_csv).is_file(), f"Missing out_dir/{req_csv}")
                self.assertTrue((tables_dir / req_csv).is_file(), f"Missing tables_dir/{req_csv}")

            conn_chk = sqlite3.connect(str(db_path))
            try:
                gla_row = conn_chk.execute(
                    """
                    SELECT public_backbone_clade, public_phytop_ils_index, public_phytop_ih_index
                    FROM genetic_lineage_assignments
                    WHERE run_id = ? AND genetic_sample_id = 'PACTOMSCG1M';
                    """,
                    (run_id,),
                ).fetchone()
                self.assertIsNotNone(gla_row)
                self.assertEqual(gla_row[0], "MO:Packera dubia+Packera schweinitziana")
                self.assertAlmostEqual(float(gla_row[1]), 0.21, places=4)
                self.assertAlmostEqual(float(gla_row[2]), 0.08, places=4)

                mv_row = conn_chk.execute(
                    """
                    SELECT public_plastome_match, public_its_match, waxy_lineage_match,
                           public_backbone_clade, public_phytop_ils_index, public_phytop_ih_index
                    FROM v_integrative_triage_master
                    WHERE run_id = ? AND catalogNumber = 'NCU_DUB_1';
                    """,
                    (run_id,),
                ).fetchone()
                self.assertIsNotNone(mv_row)
                self.assertEqual(mv_row[0], "Packera dubia")
                self.assertEqual(mv_row[1], "Packera dubia")
                self.assertEqual(mv_row[2], "Packera dubia")
                self.assertEqual(mv_row[3], "MO:Packera dubia+Packera schweinitziana")
                self.assertAlmostEqual(float(mv_row[4]), 0.21, places=4)
                self.assertAlmostEqual(float(mv_row[5]), 0.08, places=4)
            finally:
                conn_chk.close()

    def test_public_backbone_digit_prefix_and_phytop_node_tree_lookup(self) -> None:
        """(e) Verify Dryad file routing, digit-prefixed tip normalization, and PhyTop node-tree lookup."""
        self.assertEqual(
            PUBLIC_DL_MOD._route_dryad_file(Path("PhyTop_results/MO_phytop.info.tsv"), "PhyTop_results.zip"),
            "public_packera_PhyTop_MO.tsv",
        )
        self.assertEqual(
            PUBLIC_DL_MOD._route_dryad_file(Path("PhyTop_results/MO_phytop.tree"), "PhyTop_results.zip"),
            "public_packera_PhyTop_MO.tre",
        )
        self.assertEqual(
            PUBLIC_DL_MOD._route_dryad_file(Path("PhyTop_results/ASTRAL-Pro_phytop.info.tsv"), "PhyTop_results.zip"),
            "public_packera_PhyTop_ASTRALPro.tsv",
        )
        self.assertEqual(
            PUBLIC_DL_MOD._route_dryad_file(Path("RelTime_trees/S1outTest-HybPiper-df_nexus.tre"), "RelTime_trees.zip"),
            "public_packera_RelTime_S1outTest-HybPiper-df_nexus.tre",
        )
        self.assertIsNone(
            PUBLIC_DL_MOD._route_dryad_file(Path("__MACOSX/PhyTop_results/._MO_phytop.tree"), "PhyTop_results.zip"),
        )

        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir = tmp / "genomics"
            out_dir = tmp / "out"
            geno_dir.mkdir(parents=True, exist_ok=True)
            out_dir.mkdir(parents=True, exist_ok=True)

            (geno_dir / "public_packera_MO_Astral_lpp.tre").write_text(
                "(((Packera_2aurea:0.02,Packera_2obovata:0.02):0.05,Packera_dubia:0.07):0.15,Senecio_vulgaris:0.22);\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_PhyTop_MO.tre").write_text(
                "(((Packera_2aurea:0.02,Packera_2obovata:0.02)N4:0.05,Packera_dubia:0.07)N3:0.15,Senecio_vulgaris:0.22);\n",
                encoding="utf-8",
            )
            (geno_dir / "public_packera_PhyTop_MO.tsv").write_text(
                "node\tILS_index\tIH_index\nN4\t0.92\t0.04\nN3\t0.85\t0.11\n",
                encoding="utf-8",
            )

            r_script = tmp / "check_phytop_node_lookup.R"
            r_script.write_text(
                "\n".join([
                    f"source({str(PROJECT_ROOT / 'scripts' / 'analysis' / 'step06_public_backbone.R')!r})",
                    "df <- data.frame(",
                    "  genetic_sample_id = c('S1', 'S2', 'S3', 'S4'),",
                    "  species_standardized = c('Packera aurea', 'Packera obovata', 'Packera dubia', 'Senecio vulgaris'),",
                    "  stringsAsFactors = FALSE",
                    ")",
                    f"res <- evaluate_public_backbone_and_organellar(df, {str(geno_dir)!r}, {str(out_dir)!r})",
                    "cdf <- res$concordance_df",
                    "stopifnot(all(cdf$in_mo_backbone))",
                    "stopifnot(cdf$normalized_taxon[cdf$taxon == 'Senecio vulgaris'] == 'Senecio vulgaris')",
                    "stopifnot(abs(cdf$public_phytop_ils_index[cdf$taxon == 'Packera aurea'] - 0.92) < 1e-4)",
                    "stopifnot(abs(cdf$public_phytop_ih_index[cdf$taxon == 'Packera aurea'] - 0.04) < 1e-4)",
                    "stopifnot(abs(cdf$public_phytop_ils_index[cdf$taxon == 'Packera obovata'] - 0.92) < 1e-4)",
                    "stopifnot(abs(cdf$public_phytop_ils_index[cdf$taxon == 'Packera dubia'] - 0.85) < 1e-4)",
                    "stopifnot(abs(cdf$public_phytop_ih_index[cdf$taxon == 'Packera dubia'] - 0.11) < 1e-4)",
                    "cat('PHYTOP_NODE_LOOKUP_OK\\n')",
                ]) + "\n",
                encoding="utf-8",
            )
            proc = run_warm_r_cmd(["Rscript", str(r_script)], env={"PYTHONPATH": str(PROJECT_ROOT)})
            self.assertEqual(proc.returncode, 0, f"R PhyTop node lookup test failed:\n{proc.stdout}\n{proc.stderr}")
            self.assertIn("PHYTOP_NODE_LOOKUP_OK", proc.stdout)

            dotenv_file = tmp / ".env"
            dotenv_file.write_text(
                "DRYAD_CLIENT_ID=dotenv_cid\nexport DRYAD_CLIENT_SECRET='dotenv_csec'\n",
                encoding="utf-8",
            )
            with mock.patch.dict(os.environ, {"DRYAD_API_TOKEN": "keep_existing"}, clear=True):
                PUBLIC_DL_MOD._load_dotenv_defaults([dotenv_file])
                self.assertEqual(os.environ.get("DRYAD_API_TOKEN"), "keep_existing")
                self.assertEqual(os.environ.get("DRYAD_CLIENT_ID"), "dotenv_cid")
                self.assertEqual(os.environ.get("DRYAD_CLIENT_SECRET"), "dotenv_csec")


class TestStep07GenomicPropagationContract(unittest.TestCase):
    """Comprehensive verification for exsiccatae primary duplicate linking and 4-tier Step 07 genomic propagation."""

    def test_link_sample_to_voucher_exsiccatae_primary_duplicate_priority(self) -> None:
        """Verify `_link_sample_to_voucher` selects the imaged primary duplicate across herbaria."""
        vouchers = [
            {
                "catalogNumber": "DUKE001",
                "recordNumber": "4421",
                "recordedBy": "A. E. Radford",
                "scientificName": "Packera dubia",
                "species_standardized": "Packera dubia",
                "decimalLatitude": 35.80,
                "decimalLongitude": -79.05,
                "stateProvince": "North Carolina",
                "county": "Orange",
                "eventDate": "1968-05-14",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "none",
                "image_path": "",
                "exsiccatae_key": "radford|4421|1968-05-14",
                "is_primary_duplicate": 0,
            },
            {
                "catalogNumber": "US003",
                "recordNumber": "4421",
                "recordedBy": "Radford, A. E.",
                "scientificName": "Packera dubia",
                "species_standardized": "Packera dubia",
                "decimalLatitude": 35.80,
                "decimalLongitude": -79.05,
                "stateProvince": "North Carolina",
                "county": "Orange",
                "eventDate": "1968-05-14",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "none",
                "image_path": "data/raw_vouchers/US003.webp",
                "exsiccatae_key": "radford|4421|1968-05-14",
                "is_primary_duplicate": 0,
            },
            {
                "catalogNumber": "NCU002",
                "recordNumber": "4421",
                "recordedBy": "J. R. Radford",
                "scientificName": "Packera dubia",
                "species_standardized": "Packera dubia",
                "decimalLatitude": 35.80,
                "decimalLongitude": -79.05,
                "stateProvince": "North Carolina",
                "county": "Orange",
                "eventDate": "1968-05-14",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "none",
                "image_path": "data/raw_vouchers/NCU002.webp",
                "exsiccatae_key": "radford|4421|1968-05-14",
                "is_primary_duplicate": 1,
            },
        ]
        exs_index = INGEST_MOD._build_exsiccatae_index(vouchers)
        parsed = {"population_id": "OHS1", "taxon_code": "DUB", "species_standardized": "Packera dubia"}

        # 1. Pinning unimaged DUKE001 upgrades via sibling resolution to imaged primary NCU002
        res_pin = INGEST_MOD._link_sample_to_voucher(
            "PACDUBOHS1M",
            parsed,
            {"catalogNumber": "DUKE001"},
            None,
            vouchers,
            exsiccatae_index=exs_index,
        )
        self.assertEqual(res_pin["catalogNumber"], "NCU002")
        self.assertEqual(res_pin["linkage_tier"], "G1_Direct_Voucher")
        self.assertEqual(res_pin["exsiccatae_key"], "radford|4421|1968-05-14")

        # 2. Explicit exsiccatae metadata lookup resolves NCU002 over non-primary imaged US003
        res_meta = INGEST_MOD._link_sample_to_voucher(
            "PACDUBOHS1M",
            parsed,
            {"exsiccatae_key": "radford|4421|1968-05-14"},
            None,
            vouchers,
            exsiccatae_index=None,
        )
        self.assertEqual(res_meta["catalogNumber"], "NCU002")
        self.assertEqual(res_meta["linkage_tier"], "G1_Direct_Voucher")

        # 3. Record number fallback expands exsiccatae siblings and selects primary NCU002
        res_rec = INGEST_MOD._link_sample_to_voucher(
            "4421-M",
            parsed,
            None,
            None,
            vouchers,
            exsiccatae_index=exs_index,
        )
        self.assertEqual(res_rec["catalogNumber"], "NCU002")
        self.assertEqual(res_rec["linkage_tier"], "G1_Direct_Voucher")

    def test_batched_step07_genomic_propagation_four_tiers_and_invariance(self) -> None:
        """Warm-R batched contract test covering Scenarios (a) through (e) of Step 07 genomic propagation."""
        if not rscript_available():
            self.skipTest("Rscript not installed.")

        r_code = """
        source("scripts/analysis/07_triage_dashboard_synthesis.R")
        if (!exists("run_step07_genomic_propagation", mode = "function")) {
          source("scripts/analysis/step07_genomic_propagation.R")
        }

        # ------------------------------------------------------------------
        # Scenario (a): Tier 1 Exsiccatae Propagation + Mixed-Mount/Type Veto
        # ------------------------------------------------------------------
        n_bg <- 36L
        df_exs <- data.frame(
          catalogNumber = c(sprintf("BG_%02d", 1:n_bg), "EXS_DONOR", "EXS_DUP_OK", "EXS_MIXED_VETO", "EXS_TYPE_VETO", "EXS_EXP_DONOR", "EXS_EXP_MULTI_MO", "ANC_AUR_1", "ANC_AUR_2"),
          plant_individual_id = 0L,
          species_standardized = c(rep(c("Packera dubia", "Packera anonyma"), each = 18L), rep("Packera dubia", 4L), rep("Packera aurea", 4L)),
          determiner_tier = "Tier_1_Gold",
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          n_plant_individuals = c(rep(1L, n_bg + 5L), 2L, 1L, 1L),
          morphometric_outlier_flag = c(rep(FALSE, n_bg + 5L), TRUE, FALSE, FALSE),
          type_status = c(rep("none", n_bg + 3L), "ISOTYPE", rep("none", 4L)),
          is_type_specimen = c(rep(FALSE, n_bg + 3L), TRUE, rep(FALSE, 4L)),
          recordedBy = c(sprintf("Coll_%02d", 1:n_bg), "J. R. Radford", "Radford, John R.", "A. Collector", "J. R. Radford", "E. Weakley", "E. Weakley", "A. Ray", "B. Ray"),
          recordNumber = c(sprintf("%d", 1000L + seq_len(n_bg)), "4421", " 4421 ", "4421", "4421", "901", "901", "801", "802"),
          eventDate = c(rep("1990-05-01", n_bg), "1968-05-14", "1968-05-14T12:00:00Z", "1968-05-14", "1968-05-14", rep("1995-06-01", 4L)),
          exsiccatae_key = c(rep(NA_character_, n_bg + 4L), "explicit_ex_99", "explicit_ex_99", NA_character_, NA_character_),
          PC1 = c(seq(-1.5, 1.5, length.out = n_bg), 0.10, 0.14, 0.10, 0.12, 0.20, 0.22, -0.50, -0.45),
          PC2 = c(rep(c(-0.8, 0.8), 18L), 0.10, 0.12, 12.00, 0.11, 0.15, 0.16, -0.20, -0.25),
          cda_predicted_taxon = c(rep(c("Packera dubia", "Packera anonyma"), each = 18L), rep("Packera dubia", 4L), rep("Packera aurea", 4L)),
          cda_posterior_prob = 0.88,
          edaphic_best_fit_taxon = c(rep(c("Packera dubia", "Packera anonyma"), each = 18L), rep("Packera dubia", 4L), rep("Packera aurea", 4L)),
          rf_predicted_prob = 0.85,
          delta_doy_resid = 1.0,
          is_active_anthesis = TRUE,
          geno_linkage_tier = c(c("G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, n_bg - 2L)), "G1_Direct_Voucher", NA_character_, NA_character_, NA_character_, "G1_Direct_Voucher", NA_character_, "G1_Direct_Voucher", "G1_Direct_Voucher"),
          geno_quality_pass = c(c(TRUE, TRUE, rep(NA, n_bg - 2L)), TRUE, NA, NA, NA, TRUE, NA, TRUE, TRUE),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c(c("Packera dubia", "Packera anonyma", rep(NA_character_, n_bg - 2L)), "Packera serpenticola", NA_character_, NA_character_, NA_character_, "Packera aurea", NA_character_, "Packera aurea", "Packera aurea"),
          geno_posterior_prob = c(c(0.95, 0.94, rep(NA_real_, n_bg - 2L)), 0.96, NA_real_, NA_real_, NA_real_, 0.95, NA_real_, 0.95, 0.95),
          geno_margin = 0.85,
          stringsAsFactors = FALSE
        )
        # Keep EXS_MIXED_VETO in the same normalized Radford 4421 gathering
        df_exs$recordedBy[df_exs$catalogNumber == "EXS_MIXED_VETO"] <- "J. R. Radford"
        res_a <- propagate_exsiccatae_isoduplicates(df_exs)
        i_dup <- which(res_a$catalogNumber == "EXS_DUP_OK")
        i_mix <- which(res_a$catalogNumber == "EXS_MIXED_VETO")
        i_typ <- which(res_a$catalogNumber == "EXS_TYPE_VETO")
        i_mul <- which(res_a$catalogNumber == "EXS_EXP_MULTI_MO")
        stopifnot(identical(res_a$genomic_propagation_tier[i_dup], "G1_Exsiccatae_Isoduplicate"))
        stopifnot(identical(res_a$geno_linkage_tier[i_dup], "G1_Exsiccatae_Isoduplicate"))
        stopifnot(identical(res_a$geno_predicted_taxon[i_dup], "Packera serpenticola"))
        stopifnot(is.na(res_a$genomic_propagation_tier[i_mix]) && is.na(res_a$geno_predicted_taxon[i_mix]))
        stopifnot(is.na(res_a$genomic_propagation_tier[i_typ]) && is.na(res_a$geno_predicted_taxon[i_typ]))
        stopifnot(is.na(res_a$genomic_propagation_tier[i_mul]) && is.na(res_a$geno_predicted_taxon[i_mul]))

        sc_a <- compute_probabilistic_consensus_scores(res_a)
        wg_donor <- sc_a$w_geno_kappa[sc_a$catalogNumber == "EXS_DONOR"]
        wg_dup <- sc_a$w_geno_kappa[sc_a$catalogNumber == "EXS_DUP_OK"]
        stopifnot(is.finite(wg_donor) && is.finite(wg_dup) && wg_donor > 0)
        stopifnot(abs((wg_dup / wg_donor) - 0.95) < 1e-3)

        # ------------------------------------------------------------------
        # Scenario (b): Tier 2 Deme Hull + Buck Creek Sympatry Veto + Non-Override
        # ------------------------------------------------------------------
        df_deme <- data.frame(
          catalogNumber = c(
            sprintf("SERP_BG_%02d", 1:6),
            sprintf("AUR_BG_%02d", 1:6),
            "SERP_SITE1_ANC",
            "DEME_MATCH_OK",
            "DEME_VETO_EDAPH",
            "DEME_VETO_PHENO",
            "DEME_VETO_UNC",
            "SYMP_SERP_ANC",
            "SYMP_AUR_ANC",
            "DEME_VETO_SYMP"
          ),
          plant_individual_id = 0L,
          species_standardized = c(
            rep("Packera serpenticola", 6L),
            rep("Packera aurea", 6L),
            "Packera serpenticola",
            "Packera aurea",
            "Packera serpenticola",
            "Packera serpenticola",
            "Packera serpenticola",
            "Packera serpenticola",
            "Packera aurea",
            "Packera serpenticola"
          ),
          determiner_tier = c(rep("Tier_1_Gold", 13L), rep("Tier_3_Bronze", 4L), "Tier_1_Gold", "Tier_1_Gold", "Tier_3_Bronze"),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          decimalLatitude = c(
            seq(35.60, 35.65, length.out = 6L),
            seq(35.80, 35.85, length.out = 6L),
            35.100, 35.108, 35.108, 35.108, 35.108,
            35.300, 35.301, 35.302
          ),
          decimalLongitude = c(
            rep(-83.20, 6L),
            rep(-82.50, 6L),
            -83.800, -83.808, -83.808, -83.808, -83.808,
            -83.500, -83.500, -83.501
          ),
          coordinateUncertaintyInMeters = c(rep(200, 16L), 2500, 200, 200, 200),
          is_azonal_lithology = c(rep(1L, 6L), rep(0L, 6L), 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L),
          fine_soil_ph = c(rep(7.1, 6L), rep(5.5, 6L), 7.1, 7.0, 5.4, 7.1, 7.1, 7.1, 7.1, 7.1),
          type_edaphic_z_score = 0.3,
          hopkins_adjusted_doy = c(rep(125, 6L), rep(120, 6L), 125, 132, 125, 175, 125, 125, 125, 125),
          is_active_anthesis = TRUE,
          PC1 = c(
            1.4, 1.6, 1.3, 1.7, 1.5, 1.5,
            -1.4, -1.6, -1.3, -1.7, -1.5, -1.5,
            1.5, 1.55, 1.52, 1.52, 1.52, 1.5, -1.5, 0.20
          ),
          PC2 = c(
            -0.4, 0.4, -0.3, 0.3, -0.2, 0.2,
            -0.4, 0.4, -0.3, 0.3, -0.2, 0.2,
            0.0, 0.05, 0.02, 0.02, 0.02, 0.0, 0.0, 0.0
          ),
          cda_predicted_taxon = c(
            rep("Packera serpenticola", 6L),
            rep("Packera aurea", 6L),
            "Packera serpenticola",
            "Packera aurea",
            rep("Packera serpenticola", 4L),
            "Packera aurea",
            "Packera serpenticola"
          ),
          cda_posterior_prob = c(rep(0.92, 13L), 0.72, rep(0.92, 6L)),
          edaphic_best_fit_taxon = c(
            rep("Packera serpenticola", 6L),
            rep("Packera aurea", 6L),
            "Packera serpenticola",
            "Packera aurea",
            rep("Packera serpenticola", 4L),
            "Packera aurea",
            "Packera serpenticola"
          ),
          rf_predicted_prob = c(rep(0.90, 13L), 0.75, rep(0.90, 6L)),
          geno_linkage_tier = c(
            "G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 4L),
            "G1_Direct_Voucher", "G1_Direct_Voucher", rep(NA_character_, 4L),
            "G1_Direct_Voucher", rep(NA_character_, 4L),
            "G1_Direct_Voucher", "G1_Direct_Voucher", NA_character_
          ),
          geno_quality_pass = c(
            TRUE, TRUE, rep(NA, 4L),
            TRUE, TRUE, rep(NA, 4L),
            TRUE, rep(NA, 4L),
            TRUE, TRUE, NA
          ),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c(
            "Packera serpenticola", "Packera serpenticola", rep(NA_character_, 4L),
            "Packera aurea", "Packera aurea", rep(NA_character_, 4L),
            "Packera serpenticola", rep(NA_character_, 4L),
            "Packera serpenticola", "Packera aurea", NA_character_
          ),
          geno_posterior_prob = c(
            0.94, 0.93, rep(NA_real_, 4L),
            0.94, 0.93, rep(NA_real_, 4L),
            0.92, rep(NA_real_, 4L),
            0.92, 0.92, NA_real_
          ),
          geno_margin = 0.82,
          stringsAsFactors = FALSE
        )
        res_b <- propagate_conditional_deme_priors(df_deme)
        ib_ok <- which(res_b$catalogNumber == "DEME_MATCH_OK")
        stopifnot(identical(res_b$genomic_propagation_tier[ib_ok], "G2_Conditional_Deme_Hull"))
        stopifnot(identical(res_b$geno_linkage_tier[ib_ok], "G2_Conditional_Deme_Hull"))
        stopifnot(identical(res_b$geno_predicted_taxon[ib_ok], "Packera serpenticola"))
        stopifnot(is.finite(res_b$deme_mahal_d2[ib_ok]))
        stopifnot(is.finite(res_b$geno_posterior_prob[ib_ok]) && res_b$geno_posterior_prob[ib_ok] >= 0.65 && res_b$geno_posterior_prob[ib_ok] <= 0.90)
        for (vc in c("DEME_VETO_EDAPH", "DEME_VETO_PHENO", "DEME_VETO_UNC", "DEME_VETO_SYMP")) {
          iv <- which(res_b$catalogNumber == vc)
          stopifnot(is.na(res_b$genomic_propagation_tier[iv]) && is.na(res_b$geno_predicted_taxon[iv]))
        }
        sc_b <- compute_probabilistic_consensus_scores(res_b)
        stopifnot(abs((sc_b$w_geno_kappa[ib_ok] / sc_b$w_geno_kappa[sc_b$catalogNumber == "SERP_SITE1_ANC"]) - 0.75) < 1e-3)
        dec_b <- apply_taxonomic_decision_matrix(res_b)
        stopifnot(dec_b$taxonomic_treatment_subtype[dec_b$catalogNumber == "DEME_MATCH_OK"] != "Misidentification_Reassign")

        # ------------------------------------------------------------------
        # Scenario (c): Tier_0 Anchor Eviction + LOPO-CV GDA & Overlap Abstention
        # ------------------------------------------------------------------
        df_t0 <- data.frame(
          catalogNumber = c(
            sprintf("SERP_T0_%02d", 1:4),
            sprintf("DUB_T0_%02d", 1:4),
            sprintf("ANO_T0_%02d", 1:4),
            "MISID_GOLD_EVICT",
            "F1_ADMIX_EVICT",
            "UNSEQ_SERP",
            "UNSEQ_OVERLAP",
            "UNSEQ_OSR_ABST"
          ),
          plant_individual_id = 0L,
          species_standardized = c(
            rep("Packera serpenticola", 4L),
            rep("Packera dubia", 4L),
            rep("Packera anonyma", 4L),
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia",
            "Packera dubia"
          ),
          determiner_tier = c(rep("Tier_2_Silver", 12L), "Tier_1_Gold", "Tier_1_Gold", rep("Tier_3_Bronze", 3L)),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          population_deme_id = c(rep(c("P1", "P2"), 6L), "P3", "P3", "P4", "P4", "P4"),
          PC1 = c(
            5.8, 6.2, 6.0, 5.9,
            -2.1, -1.9, -2.0, -2.05,
            -1.95, -2.05, -2.0, -2.1,
            6.1, -0.5, 6.0, -2.02, 22.0
          ),
          PC2 = c(
            5.9, 6.1, 6.0, 6.05,
            0.1, -0.1, 0.0, 0.05,
            -0.05, 0.05, 0.0, -0.1,
            6.0, 0.0, 6.0, 0.0, 22.0
          ),
          hopkins_adjusted_doy = c(rep(125, 4L), rep(115, 4L), rep(116, 4L), 210, 120, 125, 115, 120),
          is_active_anthesis = TRUE,
          geno_linkage_tier = c(rep("G1_Direct_Voucher", 14L), rep(NA_character_, 3L)),
          geno_quality_pass = c(rep(TRUE, 14L), rep(NA, 3L)),
          is_geno_admixed = c(rep(FALSE, 13L), TRUE, rep(FALSE, 3L)),
          geno_predicted_taxon = c(
            rep("Packera serpenticola", 4L),
            rep("Packera dubia", 4L),
            rep("Packera anonyma", 4L),
            "Packera serpenticola",
            "Packera dubia",
            rep(NA_character_, 3L)
          ),
          geno_posterior_prob = c(rep(0.96, 13L), 0.55, rep(NA_real_, 3L)),
          geno_margin = c(rep(0.88, 13L), 0.10, rep(NA_real_, 3L)),
          interclass_het_Q12 = c(rep(0.02, 13L), 0.72, rep(NA_real_, 3L)),
          stringsAsFactors = FALSE
        )
        t0_res <- resolve_tier0_anchor_overrides(df_t0)
        i_mg <- which(df_t0$catalogNumber == "MISID_GOLD_EVICT")
        i_f1 <- which(df_t0$catalogNumber == "F1_ADMIX_EVICT")
        stopifnot(isTRUE(t0_res$evict_from_given_anchor[i_mg]) && isTRUE(t0_res$is_tier0_dna_anchor[i_mg]) && t0_res$tier0_weight[i_mg] == 4.0)
        stopifnot(isTRUE(t0_res$evict_from_given_anchor[i_f1]) && !isTRUE(t0_res$is_tier0_dna_anchor[i_f1]) && is.na(t0_res$tier0_weight[i_f1]))

        proj_c <- fit_tier0_calibrated_projection(df_t0)
        iu_s <- which(proj_c$catalogNumber == "UNSEQ_SERP")
        iu_o <- which(proj_c$catalogNumber == "UNSEQ_OVERLAP")
        iu_osr <- which(proj_c$catalogNumber == "UNSEQ_OSR_ABST")
        stopifnot(identical(proj_c$tier0_abstention_flag[iu_s], 0L))
        stopifnot(identical(proj_c$tier0_calibrated_taxon[iu_s], "Packera serpenticola"))
        stopifnot(identical(proj_c$tier0_abstention_flag[iu_o], 1L))
        stopifnot(identical(proj_c$tier0_calibrated_taxon[iu_o], "Sensu_Lato_DNA_Required"))
        stopifnot(identical(proj_c$tier0_abstention_flag[iu_osr], 1L))
        stopifnot(identical(proj_c$tier0_calibrated_taxon[iu_osr], "Sensu_Lato_DNA_Required"))

        # Confirm anchor eviction excludes MISID_GOLD_EVICT (doy=210) from Packera dubia phenology baseline
        dec_t0 <- apply_taxonomic_decision_matrix(proj_c)
        stopifnot(abs(dec_t0$delta_doy_resid[dec_t0$catalogNumber == "DUB_T0_01"]) < 3.0)
        stopifnot(grepl("[Tier_0 Abstention: Phenotypic overlap across genetically distinct lineages; DNA required]", dec_t0$synthesis_rationale[dec_t0$catalogNumber == "UNSEQ_OVERLAP"], fixed = TRUE))

        # ------------------------------------------------------------------
        # Scenario (d): Tier 4 Cluster-Wide gdi Governance (0.09 block & 0.84 support)
        # ------------------------------------------------------------------
        n_k2 <- 15L
        df_k2 <- data.frame(
          catalogNumber = c(sprintf("DUB_CORE_%02d", 1:5), sprintf("K2_ANC_%02d", 1:3), sprintf("K2_UNSEQ_%02d", 1:n_k2)),
          plant_individual_id = 0L,
          institutionCode = c(rep("NCU", 5L), rep("NCU", 3L), rep(c("NCU", "DUKE", "US"), 5L)),
          recordedBy = c(rep("CoreColl", 5L), "AncA", "AncB", "AncA", rep(c("Smith", "Jones", "Brown"), 5L)),
          recordNumber = sprintf("%d", 500L + seq_len(8L + n_k2)),
          species_standardized = "Packera dubia",
          determiner_tier = c(rep("Tier_1_Gold", 8L), rep("Tier_3_Bronze", n_k2)),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          gmm_cluster = c(rep(1L, 5L), rep(2L, 3L + n_k2)),
          morphometric_outlier_flag = c(rep(FALSE, 5L), rep(TRUE, 3L + n_k2)),
          PC1 = c(seq(-0.3, 0.3, length.out = 5L), 2.3, 2.5, 2.4, seq(2.2, 2.6, length.out = n_k2)),
          PC2 = c(rep(0.0, 5L), 0.1, -0.1, 0.0, rep(c(0.05, -0.05, 0.0), 5L)),
          cda_predicted_taxon = "Packera dubia",
          cda_posterior_prob = c(rep(0.92, 5L), rep(0.45, 3L + n_k2)),
          edaphic_best_fit_taxon = "Packera dubia",
          edaphic_niche_overlap_d = c(rep(0.85, 5L), rep(0.15, 3L + n_k2)),
          type_edaphic_z_score = c(rep(0.4, 5L), rep(3.2, 3L + n_k2)),
          decimalLatitude = c(rep(35.80, 5L), 35.81, 35.82, 35.95, rep(c(35.81, 35.95, 36.12, 36.30, 36.48), 3L)),
          decimalLongitude = c(rep(-79.00, 5L), -79.01, -79.02, -79.15, rep(c(-79.01, -79.15, -79.35, -79.55, -79.75), 3L)),
          hopkins_adjusted_doy = c(rep(115, 5L), rep(155, 3L + n_k2)),
          delta_doy_resid = c(rep(1.0, 5L), rep(40.0, 3L + n_k2)),
          is_active_anthesis = TRUE,
          min_congener_mahal_d2 = c(rep(45.0, 5L), rep(30.0, 3L + n_k2)),
          open_set_p_value = c(rep(0.001, 5L), rep(0.001, 3L + n_k2)),
          geno_linkage_tier = c(rep(NA_character_, 5L), rep("G1_Direct_Voucher", 3L), rep(NA_character_, n_k2)),
          geno_quality_pass = c(rep(NA, 5L), rep(TRUE, 3L), rep(NA, n_k2)),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c(rep(NA_character_, 5L), rep("Packera dubia", 3L), rep(NA_character_, n_k2)),
          geno_posterior_prob = c(rep(NA_real_, 5L), rep(0.95, 3L), rep(NA_real_, n_k2)),
          geno_margin = c(rep(NA_real_, 5L), rep(0.85, 3L), rep(NA_real_, n_k2)),
          coalescent_gdi = c(rep(NA_real_, 5L), rep(0.09, 3L), rep(NA_real_, n_k2)),
          stringsAsFactors = FALSE
        )
        prop_k2 <- run_step07_genomic_propagation(df_k2)
        unseq_k2 <- which(grepl("^K2_UNSEQ_", prop_k2$catalogNumber))
        stopifnot(all(is.finite(prop_k2$cluster_propagated_gdi[unseq_k2])))
        stopifnot(all(abs(prop_k2$cluster_propagated_gdi[unseq_k2] - 0.09) < 1e-4))
        dec_k2 <- apply_taxonomic_decision_matrix(prop_k2, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0)
        stopifnot(!any(dec_k2$taxonomic_treatment_subtype[unseq_k2] == "Putative_Cryptic_Species"))
        stopifnot(all(grepl("[Cluster Coalescent Oversplit Block: propagated gdi=0.09]", dec_k2$synthesis_rationale[unseq_k2], fixed = TRUE)))

        # Novel lineage with gdi=0.84: unsequenced non-abstained member supported vs overlap-abstained blocked
        df_nov <- data.frame(
          catalogNumber = c(
            sprintf("NOV_DUB_%02d", 1:4),
            sprintf("NOV_ANO_%02d", 1:4),
            sprintf("NOV_ANC_%02d", 1:4),
            "NOVEL_UNSEQ_OK",
            "NOVEL_UNSEQ_ABST"
          ),
          plant_individual_id = 0L,
          species_standardized = c(rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera dubia", 6L)),
          determiner_tier = c(rep("Tier_2_Silver", 12L), "Tier_3_Bronze", "Tier_3_Bronze"),
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          gmm_cluster = c(rep(1L, 8L), rep(3L, 6L)),
          morphometric_outlier_flag = c(rep(FALSE, 8L), rep(TRUE, 6L)),
          population_deme_id = c(rep(c("PA", "PB"), 6L), "PC", "PC"),
          PC1 = c(-3.1, -2.9, -3.0, -3.05, -2.95, -3.05, -3.0, -2.98, 4.9, 5.1, 5.0, 5.05, 5.0, -3.0),
          PC2 = c(-0.05, 0.05, 0.0, 0.02, 0.02, -0.02, 0.0, 0.04, 4.9, 5.1, 5.0, 4.95, 5.0, 0.0),
          cda_predicted_taxon = c(rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera dubia", 6L)),
          cda_posterior_prob = c(rep(0.90, 8L), rep(0.45, 6L)),
          edaphic_best_fit_taxon = "Packera dubia",
          edaphic_niche_overlap_d = c(rep(0.80, 8L), rep(0.55, 6L)),
          type_edaphic_z_score = 0.5,
          delta_doy_resid = 4.0,
          is_active_anthesis = TRUE,
          min_congener_mahal_d2 = c(rep(35.0, 8L), rep(30.0, 6L)),
          open_set_p_value = 0.001,
          geno_linkage_tier = c(rep("G1_Direct_Voucher", 12L), NA_character_, NA_character_),
          geno_quality_pass = c(rep(TRUE, 12L), NA, NA),
          is_geno_admixed = FALSE,
          geno_predicted_taxon = c(rep("Packera dubia", 4L), rep("Packera anonyma", 4L), rep("Packera sp. nov. A", 4L), NA_character_, NA_character_),
          geno_posterior_prob = c(rep(0.96, 12L), NA_real_, NA_real_),
          geno_margin = c(rep(0.88, 12L), NA_real_, NA_real_),
          coalescent_gdi = c(rep(0.50, 8L), rep(0.84, 4L), NA_real_, NA_real_),
          stringsAsFactors = FALSE
        )
        prop_nov <- run_step07_genomic_propagation(df_nov)
        i_nok <- which(prop_nov$catalogNumber == "NOVEL_UNSEQ_OK")
        i_nab <- which(prop_nov$catalogNumber == "NOVEL_UNSEQ_ABST")
        stopifnot(identical(prop_nov$tier0_abstention_flag[i_nok], 0L))
        stopifnot(identical(prop_nov$tier0_calibrated_taxon[i_nok], "Packera sp. nov. A"))
        stopifnot(abs(prop_nov$cluster_propagated_gdi[i_nok] - 0.84) < 1e-4)
        stopifnot(identical(prop_nov$tier0_abstention_flag[i_nab], 1L))
        stopifnot(is.na(prop_nov$cluster_propagated_gdi[i_nab]))

        dec_nov <- apply_taxonomic_decision_matrix(prop_nov, intraspecific_mode = TRUE, gmm_two_delta_bic = 16.0)
        stopifnot(identical(dec_nov$taxonomic_treatment_subtype[i_nok], "Putative_Cryptic_Species"))
        stopifnot(grepl("[Candidate Member via Cluster-Propagated gdi=0.84; confirm by targeted sequencing]", dec_nov$synthesis_rationale[i_nok], fixed = TRUE))
        stopifnot(dec_nov$taxonomic_treatment_subtype[i_nab] != "Putative_Cryptic_Species")
        stopifnot(grepl("[Tier_0 Abstention: Phenotypic overlap across genetically distinct lineages; DNA required]", dec_nov$synthesis_rationale[i_nab], fixed = TRUE))

        # ------------------------------------------------------------------
        # Scenario (e): Zero-Genomics Invariance
        # ------------------------------------------------------------------
        df_zero <- data.frame(
          catalogNumber = sprintf("Z_%02d", 1:6),
          plant_individual_id = 0L,
          species_standardized = rep(c("Packera dubia", "Packera anonyma"), each = 3L),
          determiner_tier = "Tier_1_Gold",
          morph_tier_source = "Basal_Rosette",
          basal_leaf_count = 3L,
          PC1 = c(-2.0, -2.1, -1.9, 2.0, 2.1, 1.9),
          PC2 = c(0.1, -0.1, 0.0, 0.1, -0.1, 0.0),
          cda_predicted_taxon = rep(c("Packera dubia", "Packera anonyma"), each = 3L),
          cda_posterior_prob = 0.90,
          edaphic_best_fit_taxon = rep(c("Packera dubia", "Packera anonyma"), each = 3L),
          rf_predicted_prob = 0.88,
          delta_doy_resid = 2.0,
          is_active_anthesis = TRUE,
          geno_linkage_tier = NA_character_,
          geno_quality_pass = NA,
          geno_predicted_taxon = NA_character_,
          geno_posterior_prob = NA_real_,
          stringsAsFactors = FALSE
        )
        zp <- run_step07_genomic_propagation(df_zero)
        stopifnot(isFALSE(attr(zp, "genomic_propagation_summary")$active))
        for (pcol in c("genomic_propagation_tier", "tier0_calibrated_taxon", "tier0_calibrated_prob", "tier0_abstention_flag", "deme_mahal_d2", "cluster_propagated_gdi")) {
          stopifnot(all(is.na(zp[[pcol]])))
        }
        r_base <- apply_taxonomic_decision_matrix(df_zero)
        r_prop <- apply_taxonomic_decision_matrix(zp)
        for (ccol in c("probabilistic_consensus_score", "authority_weight", "taxonomic_treatment_subtype", "taxonomic_status_call", "recommended_determination", "synthesis_rationale")) {
          stopifnot(identical(r_base[[ccol]], r_prop[[ccol]]))
        }
        cat("STEP07_FOUR_TIER_PROPAGATION_ALL_OK\\n")
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            r_file = Path(tmpdir) / "test_step07_genomic_prop.R"
            r_file.write_text(r_code, encoding="utf-8")
            proc = run_warm_r_cmd(
                ["Rscript", str(r_file)],
                stage_key="step07",
                cwd=PROJECT_ROOT,
                skip_figures=True,
            )
        self.assertEqual(proc.returncode, 0, f"Step 07 genomic propagation R suite failed:\n{proc.stderr}\n{proc.stdout}")
        self.assertIn("STEP07_FOUR_TIER_PROPAGATION_ALL_OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()

