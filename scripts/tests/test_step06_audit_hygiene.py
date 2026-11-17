#!/usr/bin/env python3
"""Audit regression suite for Step 06 ingestion, cytotypes, organellar reconciliation, and R spatial quarantine."""

from __future__ import annotations

import csv, importlib, json, math, os, shutil, sqlite3, subprocess, sys, tempfile, unittest
from pathlib import Path
from typing import Any, Dict

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection, upsert_genomic_samples
from scripts.core.db_schema import init_db_schema

ING = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")
PUB_DL = importlib.import_module("scripts.data_prep.download_public_packera_sequences")
OCU = importlib.import_module("scripts.data_prep.organellar_cytonuclear_utils")

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)

_WARM_R_CODE = r"""
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
args <- commandArgs(trailingOnly = TRUE)
out_json <- args[1L]; db_path <- args[2L]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/core/db_utils.R")
source("scripts/analysis/step06_spatial_landscape_genomics.R")
source("scripts/analysis/step06_phylogenomics_utils.R")

audit_df <- data.frame(
  genetic_sample_id = c("MISID_DUB1", "CAP_SER1", "REF_DUB1"),
  species_standardized = c("Packera aurea", "Packera serpenticola", "Packera dubia"),
  plastid_clade = c("cp_Clade_1", "cp_Clade_2", "cp_Clade_1"),
  cytonuclear_discordant = c(1L, 0L, 0L),
  public_plastome_match = c("Packera dubia", "Packera aurea", "Packera dubia"),
  public_its_match = c("Packera dubia", "Packera serpenticola", "Packera dubia"),
  waxy_lineage_match = c("Packera dubia", "Packera serpenticola", "Packera dubia"),
  geno_predicted_taxon = c("Packera dubia", "Packera serpenticola", "Packera dubia"),
  geno_posterior_prob = c(0.96, 0.95, 0.98),
  is_geno_admixed = c(0L, 0L, 0L),
  stringsAsFactors = FALSE
)
place_df <- data.frame(
  genetic_sample_id = c("MISID_DUB1", "CAP_SER1"),
  locus = c("plastid", "plastid"),
  best_match_taxon = c("Packera dubia", "Packera aurea"),
  kmer_containment = c(0.95, 0.94),
  n_shared_11mers = c(110L, 120L),
  stringsAsFactors = FALSE
)
rec_df <- reconcile_step06_cytonuclear_discordance(audit_df, placement_df = place_df)

set.seed(42L)
pass_df <- data.frame(
  genetic_sample_id = c("CLEAN_1", "CLEAN_2", "CLEAN_3", "CLEAN_4", "BAD_FLAG", "BAD_SQL", "BAD_ROOF"),
  catalogNumber = c("V_OK_1", "V_OK_2", "V_OK_3", "V_OK_4", "V_OK_5", "V_CENT", ""),
  population_id = c("P1", "P1", "P2", "P2", "P3", "P3", "P4"),
  decimalLatitude = c(35.10, 35.40, 36.10, 36.50, 35.80, 35.95, 35.9108),
  decimalLongitude = c(-79.10, -79.40, -81.10, -81.50, -80.20, -80.30, -79.0505),
  is_centroid_or_inst_coord = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
d_mat <- as.matrix(stats::dist(matrix(rnorm(7L * 8L), nrow = 7L)))
rownames(d_mat) <- colnames(d_mat) <- pass_df$genetic_sample_id
pco_c <- matrix(rnorm(7L * 3L), nrow = 7L, ncol = 3L, dimnames = list(pass_df$genetic_sample_id, paste0("PCo", 1:3)))
land <- run_step06_landscape_genomics(pass_df, d_mat, db_path = db_path, demo = TRUE, n_perm = 19L)
cong <- .compute_posthoc_congruence(pass_df, d_mat, pco_c, db_path = db_path, run_id = "test_r", demo = TRUE)

write_json(list(
  rec_disc = as.integer(rec_df$cytonuclear_discordant),
  land_n = as.integer(land$summary_df$n_georef_samples[1L]),
  cong_n = as.integer(cong$n_samples[1L])
), out_json, auto_unbox = TRUE)
"""


def _mk_vcf(samples, locus_specs):
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for pos, (lid, cells) in enumerate(locus_specs, start=100):
        fmt_cells = [f"{gt}:{dp}:{r},{a}" for gt, dp, r, a in cells]
        lines.append(f"chr1\t{pos}\t{lid}\tA\tG\t99\tPASS\t.\tGT:DP:AD\t" + "\t".join(fmt_cells))
    return "\n".join(lines) + "\n"


class TestStep06AuditHygiene(unittest.TestCase):
    """Focused Step 06 regression checks matching UNC-BOT-SOP-2026-04-REV5 audit criteria."""

    r_res: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        if SKIP_R_COND:
            return
        with tempfile.TemporaryDirectory(prefix="step06_audit_r_") as td:
            tmp = Path(td)
            db_p = tmp / "spatial_test.sqlite"
            conn = sqlite3.connect(str(db_p))
            init_db_schema(conn)
            conn.executemany(
                "INSERT INTO vouchers (catalogNumber, species_standardized, decimalLatitude, decimalLongitude, is_centroid_or_inst_coord) VALUES (?, ?, ?, ?, ?);",
                [("V_OK_1", "Packera dubia", 35.10, -79.10, 0), ("V_CENT", "Packera dubia", 35.95, -80.30, 1)],
            )
            conn.commit()
            conn.close()
            r_path, out_json = tmp / "probe.R", tmp / "out.json"
            r_path.write_text(_WARM_R_CODE, encoding="utf-8")
            proc = subprocess.run(["Rscript", str(r_path), str(out_json), str(db_p)], cwd=PROJECT_ROOT, capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                raise AssertionError(f"Warm R probe failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
            cls.r_res = json.loads(out_json.read_text(encoding="utf-8"))

    def test_taxon_stratified_eb_polyrad_dosage_and_2x_default(self) -> None:
        dos = np.array([[0.90]] * 4 + [[0.10]] * 3 + [[0.50]] + [[0.90]] + [[0.10]], dtype=np.float32)
        ref = np.full((10, 1), 3.0, dtype=np.float32)
        alt = np.full((10, 1), 1.0, dtype=np.float32)
        cytos = ["4x"] * 10
        taxa = ["Packera aurea"] * 4 + ["Packera dubia"] * 4 + ["Packera aurea", "Packera dubia"]
        eb_pan = ING._compute_eb_polyrad_dosage(dos, ref, alt, cytos, taxa_labels=None)
        eb_tax = ING._compute_eb_polyrad_dosage(dos, ref, alt, cytos, taxa_labels=taxa)
        self.assertLess(float(eb_tax[7, 0]), float(eb_pan[7, 0]) - 0.05)
        self.assertGreater(float(eb_tax[8, 0]), float(eb_pan[8, 0]) + 0.05)

        d_mix = np.array([[0.5], [0.5]], dtype=np.float32)
        r_mix = np.array([[30.0], [30.0]], dtype=np.float32)
        a_mix = np.array([[10.0], [10.0]], dtype=np.float32)
        eb_def = ING._compute_eb_polyrad_dosage(d_mix, r_mix, a_mix, ["mixed/uncertain", "mixed/uncertain"], fallback_ploidies=None)
        eb_fb = ING._compute_eb_polyrad_dosage(d_mix, r_mix, a_mix, ["mixed/uncertain", "mixed/uncertain"], fallback_ploidies=[2, 4])
        self.assertGreater(float(eb_def[0, 0]), 0.45)
        self.assertAlmostEqual(float(eb_def[0, 0]), float(eb_fb[0, 0]), places=5)
        self.assertGreater(float(eb_def[1, 0]), float(eb_fb[1, 0]) + 0.15)

    def test_nquire_em_pruned_free_8x_priors_and_signed_scores(self) -> None:
        rng = np.random.default_rng(202611)
        fit = ING._fit_nquire_em_models(np.clip(rng.normal(0.50, 0.05, size=50), 0.12, 0.88))
        self.assertNotIn("free", fit["logL_map"])
        self.assertNotIn("free", fit["bic_map"])
        self.assertNotIn("delta_logl", fit)

        self.assertEqual(ING._infer_sample_cytotype(np.array([10.0] * 5), np.array([10.0] * 5), np.full(5, 0.5, dtype=np.float32), 8), ("8x", 0.0))
        self.assertEqual(ING._infer_sample_cytotype(np.array([10.0] * 5), np.array([10.0] * 5), np.full(5, 0.5, dtype=np.float32), None, species_modal_ploidy=8, is_multi_cytotype=False), ("8x", 0.0))

        cy_2x, sc_2x = ING._infer_sample_cytotype(np.full(40, 20.0), np.full(40, 20.0), np.full(40, 0.5, dtype=np.float32), None)
        self.assertEqual(cy_2x, "2x")
        self.assertLessEqual(float(sc_2x), -4.0)

        cy_4x, sc_4x = ING._infer_sample_cytotype(np.array([30.0, 10.0] * 20), np.array([10.0, 30.0] * 20), np.full(40, 0.5, dtype=np.float32), None)
        self.assertEqual(cy_4x, "4x")
        self.assertGreaterEqual(float(sc_4x), 4.0)

        weak_4x_ratios = np.concatenate([np.full(8, 0.25), np.full(8, 0.75), np.full(4, 0.50)])
        alt_w = np.round(weak_4x_ratios * 20.0)
        ref_w = 20.0 - alt_w
        dos_w = np.full(20, 0.5, dtype=np.float32)
        cy_unp, sc_unp = ING._infer_sample_cytotype(ref_w, alt_w, dos_w, None)
        if cy_unp == "mixed/uncertain":
            self.assertGreater(float(sc_unp), 0.0)
        self.assertFalse(cy_unp == "2x" and float(sc_unp) > 0.0)
        cy_p2, sc_p2 = ING._infer_sample_cytotype(ref_w, alt_w, dos_w, 2)
        if cy_p2 == "2x":
            self.assertLessEqual(float(sc_p2), 0.0)

    def test_companion_vcf_depths_intersected_with_dosage_table(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir, db_p = tmp / "genomics", tmp / "out", tmp / "db.sqlite"
            geno_dir.mkdir()
            samples = ["S1", "S2", "S3", "S4X", "S5", "S6"]
            specs = [("L_DBIAS", [("0/1", 40, 36, 4)] * 3 + [("0/0", 40, 40, 0)] * 3)]
            for k in range(20):
                c4 = ("0/1", 40, 30, 10) if k % 2 == 0 else ("0/1", 40, 10, 30)
                specs.append((f"L_OK_{k}", [("0/0", 30, 30, 0), ("1/1", 30, 0, 30), ("0/1", 30, 15, 15), c4, ("0/0", 30, 30, 0), ("1/1", 30, 0, 30)]))
            (geno_dir / "cohort.vcf").write_text(_mk_vcf(samples, specs), encoding="utf-8")
            with open(geno_dir / "polyrad_dosage.tsv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh, delimiter="\t")
                wr.writerow(["genetic_sample_id", "L_DBIAS", *[f"L_OK_{k}" for k in range(20)]])
                for s in samples:
                    dbias_v = 0.5 if s in ("S1", "S2", "S3") else 0.0
                    vals = [
                        0.27 if (s == "S4X" and k % 2 == 0)
                        else (0.73 if s == "S4X" else (0.0 if s in ("S1", "S5") else (1.0 if s in ("S2", "S6") else 0.5)))
                        for k in range(20)
                    ]
                    wr.writerow([s, dbias_v, *vals])

            res = ING.run_genomic_ingestion(geno_dir, out_dir, db_path=db_p, min_hq_reads=0)
            self.assertEqual(res["n_loci_retained"], 20)
            with open(out_dir / "hdplot_locus_diagnostics.tsv", encoding="utf-8") as fh:
                hd = {r["locus_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(hd["L_DBIAS"]["hdplot_keep"], "0")
            self.assertIn("HDplot_|D|>", hd["L_DBIAS"]["filter_reason"])
            with open(out_dir / "parsed_genomic_samples.csv", encoding="utf-8") as fh:
                sm = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            self.assertEqual(sm["S4X"]["inferred_cytotype"], "4x")
            with open(out_dir / "allele_dosage_matrix.tsv", encoding="utf-8") as fh:
                dm = {r["genetic_sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertAlmostEqual(float(dm["S4X"]["L_OK_0"]), 0.27, places=3)

    def test_upsert_preserve_non_null_and_g3_coord_hygiene(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, out_dir, db_p = tmp / "genomics", tmp / "out", tmp / "db.sqlite"
            geno_dir.mkdir()
            conn = get_db_connection(db_path=db_p, init_schema=True, workspace_root=tmp)
            try:
                upsert_genomic_samples(conn, [{"genetic_sample_id": "S_PART", "plastid_clade": "cp_Clade_1", "hq_reads": 250000, "inferred_cytotype": "2x"}])
                upsert_genomic_samples(conn, [{"genetic_sample_id": "S_PART", "plastid_clade": None, "hq_reads": None, "inferred_cytotype": "4x"}])
                row = conn.execute("SELECT plastid_clade, hq_reads, inferred_cytotype FROM genomic_samples WHERE genetic_sample_id='S_PART';").fetchone()
                self.assertEqual(tuple(row), ("cp_Clade_1", 250000, "4x"))
            finally:
                conn.close()

            sids = ["S_ZERO", "S_FLIP", "S_NCU", "S_GOOD"]
            specs = [("L1", [("0/0", 20, 20, 0), ("1/1", 20, 0, 20), ("0/1", 20, 10, 10), ("0/1", 20, 10, 10)])]
            (geno_dir / "c.vcf").write_text(_mk_vcf(sids, specs), encoding="utf-8")
            with open(geno_dir / "user_genetic_metadata.csv", "w", encoding="utf-8", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["genetic_sample_id", "decimalLatitude", "decimalLongitude", "stateProvince", "species_standardized"])
                wr.writerow(["S_ZERO", "0.0", "0.0", "NC", "Packera dubia"])
                wr.writerow(["S_FLIP", "35.91", "79.05", "NC", "Packera dubia"])
                wr.writerow(["S_NCU", "35.9108", "-79.0505", "NC", "Packera dubia"])
                wr.writerow(["S_GOOD", "36.15", "-81.65", "NC", "Packera dubia"])
            ING.run_genomic_ingestion(geno_dir, out_dir, db_path=db_p, min_hq_reads=0)
            with open(out_dir / "parsed_genomic_samples.csv", encoding="utf-8") as fh:
                sm = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            for bad in ("S_ZERO", "S_FLIP", "S_NCU"):
                self.assertEqual(sm[bad]["linkage_tier"], "G4_Ungeoreferenced_Tip")
                self.assertEqual(sm[bad]["decimalLatitude"], "")
            self.assertEqual(sm["S_GOOD"]["linkage_tier"], "G3_Georeferenced_Unvouchered")

    def test_synonym_hygiene_ncbi_waxy_and_cytonuclear_reconciliation(self) -> None:
        for rel, bad_toks in (
            ("scripts/data_prep/06_ingest_genetic_data.py", ("KNOWN_NOVEL_OR_HYBRID_POPS", "POPULATION_ALIAS_MAP")),
            ("scripts/data_prep/step06_sample_metadata_utils.py", ("KNOWN_NOVEL_OR_HYBRID_POPS", "POPULATION_ALIAS_MAP")),
            ("scripts/data_prep/step06_synonym_utils.py", ("kf584255", "42947")),
            ("scripts/data_prep/organellar_cytonuclear_utils.py", ("kf584255", "42947")),
            ("scripts/data_prep/skim_offtarget_organellar.py", ("_IUPAC_PAIR", "get_organelle_from_reads.py")),
        ):
            txt = (PROJECT_ROOT / rel).read_text(encoding="utf-8")
            for tok in bad_toks:
                self.assertNotIn(tok, txt, f"{rel} still contains forbidden token {tok}")

        cat_unk, _ = PUB_DL._classify_ncbi_record(">JX1 Packera aurea nuclear genomic contig", "ACGT" * 40)
        cat_wx, _ = PUB_DL._classify_ncbi_record(">JX2 Packera aurea granule-bound starch synthase (GBSSI) gene", "ACGT" * 40)
        self.assertIsNone(cat_unk)
        self.assertEqual(cat_wx, "waxy")

        d_rows = OCU._detect_cytonuclear_discordance(
            {"PACAURFK1M": "C1", "PACSERNCS1M": "C2"},
            {"PACAURFK1M": "I1", "PACSERNCS1M": "I2"},
            public_plastome_matches={"PACAURFK1M": "Packera dubia", "PACSERNCS1M": "Packera aurea"},
            public_its_matches={"PACAURFK1M": "Packera dubia", "PACSERNCS1M": "Packera serpenticola"},
            waxy_lineage_matches={"PACAURFK1M": "Packera dubia", "PACSERNCS1M": "Packera serpenticola"},
            plastome_match_metrics={"PACAURFK1M": (0.95, 120), "PACSERNCS1M": (0.95, 120)},
            nuclear_taxon_map={"PACAURFK1M": "Packera dubia"},
        )
        by_id = {r["genetic_sample_id"]: r for r in d_rows}
        self.assertEqual(by_id["PACAURFK1M"]["nuclear_taxon"], "Packera dubia")
        self.assertEqual(by_id["PACAURFK1M"]["cytonuclear_discordant"], 0)
        self.assertEqual(by_id["PACSERNCS1M"]["cytonuclear_discordant"], 1)

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_r_cytonuclear_reconciliation_and_spatial_quarantine(self) -> None:
        self.assertEqual(self.r_res["rec_disc"], [0, 1, 0])
        self.assertEqual(self.r_res["land_n"], 4)
        self.assertEqual(self.r_res["cong_n"], 4)


if __name__ == "__main__":
    unittest.main()
