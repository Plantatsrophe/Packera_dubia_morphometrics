#!/usr/bin/env python3
"""
Unit and end-to-end integration tests for Standalone Step 06 (Genomic &
Phylogenomic Evidence Stream) and Step 07 (Integrative Synthesis & Triage).

Covers:
  1. `.gitignore` tracking of `!data/genomics/.gitkeep` and flat payload ignores.
  2. SQLite WAL schema migrations (`20261121 -> 20261123` and `20261122 -> 20261123`),
     all 18 Phase 2 columns in `genomic_samples`, `genetic_lineage_assignments`,
     and `v_integrative_triage_master` cardinality (`plant_individuals` preserved),
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
import importlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
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
    migrate_database_schema,
)
from scripts.data_prep import skim_offtarget_organellar as skim_mod

INGEST_MOD = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")


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
        self.assertEqual(PACKERA_SCHEMA_USER_VERSION, 20261123)

        expected_gs_phase2 = {
            "hdplot_retained_loci",
            "inferred_cytotype",
            "nquire_delta_logl",
            "plastid_clade",
            "its_ribotype",
            "cytonuclear_discordant",
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
        }

        for legacy_ver in (20261121, 20261122):
            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td)
                db_path = tmp / f"legacy_{legacy_ver}.sqlite"
                conn = sqlite3.connect(str(db_path))
                init_db_schema(conn)
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
                conn.execute(f"PRAGMA user_version = {legacy_ver};")
                conn.commit()

                migrate_database_schema(conn)
                uv_after = conn.execute("PRAGMA user_version;").fetchone()[0]
                self.assertEqual(uv_after, 20261123)
                init_db_schema(conn)

                gs_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genomic_samples);").fetchall()}
                gla_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(genetic_lineage_assignments);").fetchall()}
                view_cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}

                self.assertTrue(expected_gs_phase2.issubset(gs_cols), f"Missing gs cols from {legacy_ver}: {expected_gs_phase2 - gs_cols}")
                self.assertTrue(expected_gla_phase2.issubset(gla_cols), f"Missing gla cols from {legacy_ver}: {expected_gla_phase2 - gla_cols}")
                self.assertTrue(expected_gs_phase2.issubset(view_cols), f"View missing gs cols: {expected_gs_phase2 - view_cols}")
                self.assertTrue(expected_gla_phase2.issubset(view_cols), f"View missing gla cols: {expected_gla_phase2 - view_cols}")
                conn.close()

        # Verify cardinality invariance on multi-plant sheets + multiple genomic replicates, then purge --runs-only
        with tempfile.TemporaryDirectory() as td:
            ws = Path(td)
            data_dir = ws / "data"
            tables_dir = data_dir / "tables"
            runs_dir = ws / "runs" / "run_alpha"
            for d in (tables_dir, runs_dir):
                d.mkdir(parents=True, exist_ok=True)
            db_path = data_dir / "packera_pipeline.sqlite"

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
                conn.execute(
                    """
                    INSERT INTO run_active_cohort (
                        run_id, catalogNumber, plant_individual_id,
                        species_standardized, determiner_tier, authority_weight
                    ) VALUES ('run_alpha', 'NCU001', 0, 'Packera dubia', 'Tier_1_Gold', 2.0);
                    """
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
            with mock.patch.dict(os.environ, {"PACKERA_GENOMICS_DIR": str(geno_dir)}), mock.patch("subprocess.run") as m_run:
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
            with mock.patch.dict(os.environ, {"PACKERA_GENOMICS_DIR": str(geno_dir)}), mock.patch("subprocess.run") as m_run:
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
                # Strong read-ratio bias |D| > 4.0 across 5 hets (ref=36, alt=4) -> HDplot_|D|>4.0 filtered
                ("L_PARALOG_DBIAS", [("0/1", 40, 36, 4)] * 5 + [("0/0", 40, 40, 0)]),
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
                for sid, cat, tier, sp, is_q, lat, lon, cyt, plast, its_r, cytonuc, _ in s_defs:
                    wr.writerow([
                        sid, sid, cat, tier,
                        sid[3:6], sp, sp, "POP1",
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
  PACKERA_DB_PATH = {json.dumps(str(db_path))},
  PACKERA_SQLITE_DB = {json.dumps(str(db_path))},
  PACKERA_GENOMICS_DIR = {json.dumps(str(empty_geno))},
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
Sys.setenv(PACKERA_GENOMICS_DIR = {json.dumps(str(pop_geno))})
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
make_intra_probe <- function(cat_id, gdi_val, hyb_cls, q12_val) {{
  data.frame(
    catalogNumber = cat_id,
    plant_individual_id = 0L,
    species_standardized = "Packera dubia",
    species_raw = "Packera dubia",
    recordedBy = "Fuller",
    institutionCode = "NCU",
    determiner_tier = "Tier_1_Gold",
    decimalLatitude = 35.2,
    decimalLongitude = -79.3,
    stateProvince = "NC",
    county = "Moore",
    leaf_zone = "Basal_Rosette",
    basal_leaf_count = 6L,
    is_basal_morph = TRUE,
    is_passive = FALSE,
    morph_tier_source = "Basal_Rosette",
    is_morph_outlier = TRUE,
    gmm_cluster = 2L,
    gmm_bayes_factor = 12.0,
    cda_predicted_taxon = "Packera dubia",
    cda_posterior_prob = 0.45,
    edaphic_best_fit_taxon = "Packera dubia",
    edaphic_prob = 0.20,
    type_edaphic_z_score = 3.2,
    hopkins_adjusted_doy = 155.0,
    delta_doy = 35.0,
    reproductive_trait_status = "Measured",
    is_anthesis = TRUE,
    open_set_novelty_p = 0.01,
    min_congener_mahal_d2 = 28.0,
    cryptic_cluster_cohesion_pass = TRUE,
    PC1 = 2.4, PC2 = -1.1, PC3 = 0.4,
    geno_predicted_taxon = "Packera dubia",
    geno_linkage_tier = "G1_Direct_Voucher",
    geno_quality_pass = TRUE,
    geno_posterior_prob = 0.96,
    coalescent_gdi = gdi_val,
    hybrid_generation_class = hyb_cls,
    interspecific_ancestry_q12 = q12_val,
    is_geno_admixed = FALSE,
    stringsAsFactors = FALSE
  )
}}

intra_oversplit <- apply_taxonomic_decision_matrix(
  make_intra_probe("OVERSPLIT_01", 0.12, "Pure_Parent_1", 0.02),
  intraspecific_mode = TRUE,
  gmm_two_delta_bic = 18.0
)
stopifnot(grepl("Coalescent Oversplit Block", intra_oversplit$synthesis_rationale[1L]))
stopifnot(intra_oversplit$taxonomic_treatment_subtype[1L] %in% c("Ecophenotypic_Plasticity", "Confirmed_Species"))

intra_cryptic <- apply_taxonomic_decision_matrix(
  make_intra_probe("CRYPTIC_01", 0.88, "Divergent_NonHybrid_Lineage", 0.04),
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


if __name__ == "__main__":
    unittest.main()
