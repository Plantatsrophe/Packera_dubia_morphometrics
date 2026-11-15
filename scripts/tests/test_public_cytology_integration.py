#!/usr/bin/env python3
"""Integration tests for public cytology staging, Step 06 EB priors, Step 07 barriers, and R/SQLite contracts."""

from __future__ import annotations

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
import urllib.request
import zipfile
from pathlib import Path
from unittest.mock import patch

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection
from scripts.data_prep import download_public_cytology as dl_cyto

step06 = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")


class _FakeHTTPResp:
    def __init__(self, payload: bytes) -> None:
        self._payload = payload

    def read(self) -> bytes:
        return self._payload

    def __enter__(self) -> "_FakeHTTPResp":
        return self

    def __exit__(self, *args: object) -> bool:
        return False


def _zip_bytes(entries: dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, text in entries.items():
            zf.writestr(name, text)
    return buf.getvalue()


class TestPublicCytologyIntegration(unittest.TestCase):
    def test_parse_chromosome_records_and_taxon_resolution(self) -> None:
        p1 = dl_cyto.parse_packera_chromosome_record("22II")[0]
        self.assertEqual((p1["gametophytic_n"], p1["sporophytic_2n"], p1["base_x"], p1["ploidy_level"], p1["cytotype_label"]), (22, 44, 22, 2, "2x"))
        p2 = dl_cyto.parse_packera_chromosome_record("2n = 46")[0]
        self.assertEqual((p2["sporophytic_2n"], p2["base_x"], p2["ploidy_level"], p2["cytotype_label"]), (46, 23, 2, "2x"))
        p3 = dl_cyto.parse_packera_chromosome_record("22II + 22I")[0]
        self.assertEqual((p3["sporophytic_2n"], p3["ploidy_level"], p3["is_aneuploid_or_odd_ploid"], p3["cytotype_label"]), (66, 3, 1, "3x"))
        p4 = dl_cyto.parse_packera_chromosome_record("44+0-2B")[0]
        self.assertEqual((p4["sporophytic_2n"], p4["b_chromosomes"], p4["ploidy_level"]), (44, 2, 2))
        p5 = dl_cyto.parse_packera_chromosome_record("ca. 92")[0]
        self.assertEqual((p5["sporophytic_2n"], p5["base_x"], p5["ploidy_level"], p5["is_approximate"], p5["cytotype_label"]), (92, 23, 4, 1, "4x"))

        self.assertIsNone(dl_cyto.resolve_packera_cytology_taxon("Senecio vulgaris", sporophytic_2n=40))
        res = dl_cyto.resolve_packera_cytology_taxon("Senecio tomentosus", sporophytic_2n=44)
        self.assertIsNotNone(res)
        self.assertEqual((res["species_standardized"], res["basionym_name"]), ("Packera dubia", "Senecio tomentosus"))

    def test_extract_label_cytology_mentions(self) -> None:
        m1 = dl_cyto.extract_label_cytology_mentions("Chromosome voucher: n = 22 (Kowal 1975)")
        self.assertTrue(m1 and any("22" in x for x in m1))
        self.assertEqual(dl_cyto.extract_label_cytology_mentions("Flow cytometry 2C = 3.84 pg"), ["2C = 3.84 pg"])
        self.assertEqual(dl_cyto.extract_label_cytology_mentions("Occasional, n = 12 plants along ledge"), [])
        self.assertEqual(dl_cyto.extract_label_cytology_mentions("stems 2n-branched, 44 cm tall"), [])

    def test_mocked_fetch_zip_slip_sync_db_and_sqlite_linkage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            db_p, cyto_dir = root / "test.sqlite", root / "cytology"
            conn = get_db_connection(db_path=db_p, init_schema=True, demo=False)
            try:
                conn.executemany(
                    "INSERT INTO vouchers (catalogNumber, species_standardized, recordedBy, recordNumber, institutionCode, exsiccatae_key, locality) VALUES (?, ?, ?, ?, ?, ?, ?)",
                    [
                        ("NCU-CYT-1", "Packera dubia", "R. R. Kowal", "101", "NCU", "kowal_101", "Chromosome voucher: n = 22 (Kowal 1975)"),
                        ("NCU-CYT-2", "Packera dubia", "R. R. Kowal", "101", "DUKE", "kowal_101", "Granitic outcrop"),
                        ("NCU-CYT-3", "Packera aurea", "A. Gray", "55", "GH", "gray_55", "Moist meadow"),
                    ],
                )
                conn.commit()
            finally:
                conn.close()

            ccdb_payload = json.dumps([
                {"id": "c1", "original_name": "Senecio tomentosus", "resolved_name": "Packera dubia", "gametophytic": "22", "sporophytic": "44"},
                {"id": "c2", "original_name": "Senecio vulgaris", "resolved_name": "Senecio vulgaris", "sporophytic": "40"},
                {"id": "c3", "original_name": "Packera aurea", "resolved_name": "Packera aurea", "sporophytic": "88", "voucher": "A. Gray #55"},
                {"id": "c4", "original_name": "Packera cana", "resolved_name": "Packera cana", "sporophytic": "92(1)", "parsed_n": "46"},
            ]).encode("utf-8")
            icna_real = (
                "Genus,Species / Subspecies,Variety,chromosome number n=（No. of plants examined),\"chromosome number 2n=(total karyotypic length)（No. of plants examined)\",Collected Locality,Author name ＆ Published year\n"
                "Packera A. Love ＆ D. Love 1975,anonyma (Wood) A. Love as Senecio smallii,\"based on Jones, S. B. 15422 (GA) (AF459959)\",- (0),44 (80)(1),USA: NC,K. Watanabe 1995\n"
                "Senecio L.,vulgaris,-,- (0),40 (80)(1),Europe,Smith 1990\n"
            )
            icna_zip = _zip_bytes({"asteraceae_icna.csv": "Taxon,2n\nPackera obovata,44\n", "ftp_index_chromosome_number_asteraceae.csv": icna_real})
            evil_zip = _zip_bytes({"../../escape.txt": "bad"})

            def _fake_open(req: object, timeout: int = 20) -> _FakeHTTPResp:
                url = req.full_url if hasattr(req, "full_url") else str(req)
                if "countsFull" in url:
                    return _FakeHTTPResp(ccdb_payload)
                if "index_chromosome_number_asteraceae.zip" in url:
                    return _FakeHTTPResp(icna_zip)
                raise RuntimeError(f"Unexpected URL: {url}")

            class _FakeOpener:
                def open(self, req: object, timeout: int = 20) -> _FakeHTTPResp:
                    return _fake_open(req, timeout=timeout)

            with patch.object(urllib.request, "urlopen", side_effect=_fake_open), patch.object(urllib.request, "build_opener", return_value=_FakeOpener()):
                rc = dl_cyto.main(["--cytology-dir", str(cyto_dir), "--sources", "ccdb", "icna", "labels", "--db", str(db_p)])
            self.assertEqual(rc, 0)
            self.assertEqual([p for p in cyto_dir.iterdir() if p.is_dir()], [])
            for fn in ("public_packera_cytology_records.csv", "public_packera_taxon_cytology_summary.csv", "public_packera_cytology_manifest.json"):
                self.assertTrue((cyto_dir / fn).is_file())

            conn = sqlite3.connect(str(db_p))
            try:
                recs = conn.execute("SELECT catalogNumber, cytology_linkage_tier, species_standardized, basionym_name, sporophytic_2n FROM cytological_records").fetchall()
                sums = {r[0]: r for r in conn.execute("SELECT species_standardized, modal_sporophytic_2n, modal_ploidy_level FROM taxon_cytology_summary").fetchall()}
            finally:
                conn.close()
            self.assertFalse(any("vulgaris" in str(r[2]).lower() for r in recs))
            self.assertTrue(any(r[2] == "Packera dubia" and r[3] == "Senecio tomentosus" and r[4] == 44 for r in recs))
            self.assertTrue(any(r[2] == "Packera cana" and r[4] == 92 for r in recs) and not any(r[2] == "Packera cana" and r[4] == 46 for r in recs))
            self.assertTrue(any(r[2] == "Packera anonyma" and r[3] == "Senecio smallii" and r[4] == 44 for r in recs))
            by_cat = {(r[0], r[1]) for r in recs if r[0]}
            self.assertIn(("NCU-CYT-1", "C1_Direct_Cytovoucher"), by_cat)
            self.assertIn(("NCU-CYT-2", "C2_Exsiccatae_Cytoduplicate"), by_cat)
            self.assertIn(("NCU-CYT-3", "C2_Exsiccatae_Cytoduplicate"), by_cat)
            self.assertEqual((sums["Packera dubia"][1:], sums["Packera obovata"][1:]), ((44, 2), (44, 2)))

            conn = sqlite3.connect(str(db_p))
            try:
                conn.execute("DELETE FROM cytological_records"); conn.execute("DELETE FROM taxon_cytology_summary"); conn.commit()
            finally:
                conn.close()

            with patch.object(urllib.request, "urlopen", side_effect=AssertionError("HTTP forbidden in --sync-db-only")), patch.object(urllib.request, "build_opener", side_effect=AssertionError("Opener forbidden in --sync-db-only")):
                self.assertEqual(dl_cyto.main(["--cytology-dir", str(cyto_dir), "--sync-db-only", "--db", str(db_p)]), 0)
            conn = sqlite3.connect(str(db_p))
            try:
                self.assertGreater(conn.execute("SELECT COUNT(*) FROM cytological_records").fetchone()[0], 0)
                self.assertGreater(conn.execute("SELECT COUNT(*) FROM taxon_cytology_summary").fetchone()[0], 0)
            finally:
                conn.close()

            class _EvilOpener:
                def open(self, req: object, timeout: int = 20) -> _FakeHTTPResp:
                    return _FakeHTTPResp(evil_zip)

            with patch.object(urllib.request, "build_opener", return_value=_EvilOpener()):
                with self.assertRaises(ValueError):
                    dl_cyto._fetch_icna({})

    def test_run_pipeline_preflight_public_cytology(self) -> None:
        rp = PROJECT_ROOT / "run_pipeline.sh"
        self.assertEqual(subprocess.run(["bash", "-n", str(rp)], check=False).returncode, 0)
        rp_txt = rp.read_text(encoding="utf-8")
        for flag in ("--skip-cytology-download", "--force-cytology-download", "--download-public-cytology"):
            self.assertIn(flag, rp_txt)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cyto_dir, log_f, spy = root / "cytology", root / "calls.log", root / "spy_py.sh"
            cyto_dir.mkdir()
            spy.write_text(f"#!/usr/bin/env bash\necho \"$*\" >> {json.dumps(str(log_f))}\n", encoding="utf-8")
            spy.chmod(0o755)
            fn_lines = []
            in_fn = False
            for line in rp_txt.splitlines():
                if line.startswith("run_preflight_public_cytology()"):
                    in_fn = True
                if in_fn:
                    fn_lines.append(line)
                    if line.strip() == "}":
                        break
            self.assertTrue(fn_lines, "Missing run_preflight_public_cytology in run_pipeline.sh")
            harness = "\n".join(fn_lines) + "\nrun_preflight_public_cytology\n"
            env_base = {"PATH": os.environ.get("PATH", ""), "PROJECT_ROOT": str(PROJECT_ROOT), "PYTHON_MAIN": str(spy), "PLANTMATICS_CYTOLOGY_DIR": str(cyto_dir), "CUSTOM_DB_PATH": ""}
            subprocess.run(["bash", "-c", f"DEMO_FLAG=(); SKIP_CYTOLOGY_DOWNLOAD=false; FORCE_CYTOLOGY_DOWNLOAD=false; DOWNLOAD_PUBLIC_CYTOLOGY=0\n{harness}"], env=env_base, check=True)
            (cyto_dir / "public_packera_cytology_records.csv").write_text("cytology_id\nc1\n", encoding="utf-8")
            subprocess.run(["bash", "-c", f"DEMO_FLAG=(); SKIP_CYTOLOGY_DOWNLOAD=false; FORCE_CYTOLOGY_DOWNLOAD=false; DOWNLOAD_PUBLIC_CYTOLOGY=0\n{harness}"], env=env_base, check=True)
            subprocess.run(["bash", "-c", f"DEMO_FLAG=(); SKIP_CYTOLOGY_DOWNLOAD=true; FORCE_CYTOLOGY_DOWNLOAD=false; DOWNLOAD_PUBLIC_CYTOLOGY=0\n{harness}"], env=env_base, check=True)
            subprocess.run(["bash", "-c", f"DEMO_FLAG=(--demo); SKIP_CYTOLOGY_DOWNLOAD=false; FORCE_CYTOLOGY_DOWNLOAD=false; DOWNLOAD_PUBLIC_CYTOLOGY=0\n{harness}"], env=env_base, check=True)
            logged = log_f.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(logged), 2)
            self.assertNotIn("--sync-db-only", logged[0])
            self.assertIn("--sync-db-only", logged[1])

    def test_step06_cytology_priors_and_polyrad_dosage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            db_p, gen_dir, out_dir = root / "s6.sqlite", root / "genomics", root / "parsed"
            gen_dir.mkdir(); out_dir.mkdir()
            conn = get_db_connection(db_path=db_p, init_schema=True, demo=False)
            try:
                conn.execute("INSERT INTO vouchers (catalogNumber, species_standardized) VALUES ('V-MONO-2X', 'Packera dubia')")
                conn.execute("INSERT INTO taxon_cytology_summary (species_standardized, n_cytology_records, n_direct_cytovouchers, modal_sporophytic_2n, modal_base_x, modal_ploidy_level, reported_cytotypes, is_multi_cytotype_taxon) VALUES ('Packera dubia', 8, 2, 44, 22, 2, '2x', 0)")
                conn.commit()
            finally:
                conn.close()

            self.assertEqual(step06._load_db_cytology_priors(db_p)["by_taxon"]["Packera dubia"]["modal_ploidy_level"], 2)
            cy_lo, dl_lo = step06._infer_sample_cytotype(np.array([7., 6., 5., 7., 6.]), np.array([3., 4., 5., 3., 4.]), np.full(5, 0.5, dtype=np.float32), None, species_modal_ploidy=2, is_multi_cytotype=False)
            self.assertEqual((cy_lo, dl_lo), ("2x", 0.0))

            eb = step06._compute_eb_polyrad_dosage(np.array([[0.3], [0.3]], dtype=np.float32), np.array([[7.0], [7.0]], dtype=np.float32), np.array([[3.0], [3.0]], dtype=np.float32), ["2x", "mixed/uncertain"], fallback_ploidies=[2, 4])
            self.assertGreater(float(eb[0, 0]), float(eb[1, 0]) + 0.05)

            cy_4x, dl_4x = step06._infer_sample_cytotype(np.array([30., 10.] * 12), np.array([10., 30.] * 12), np.full(24, 0.5, dtype=np.float32), None, species_modal_ploidy=2, is_multi_cytotype=False)
            self.assertEqual(cy_4x, "4x")
            self.assertGreaterEqual(float(dl_4x), 4.0)

            (gen_dir / "user_genetic_metadata.csv").write_text("genetic_sample_id,catalogNumber,species_standardized\nS_MONO_2X,V-MONO-2X,Packera dubia\n", encoding="utf-8")
            vcf_lines = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">', '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS_MONO_2X\tS_REF_ALT"] + [f"chr1\t{100 + i}\t.\tA\tG\t99\tPASS\t.\tGT:AD\t0/1:6,6\t1/1:0,12" for i in range(1, 6)]
            (gen_dir / "cohort.vcf").write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")
            res = step06.run_genomic_ingestion(gen_dir, out_dir, db_path=db_p, min_hq_reads=0)
            self.assertEqual(res["status"], "ingested")
            conn = sqlite3.connect(str(db_p))
            try:
                self.assertEqual(conn.execute("SELECT inferred_cytotype FROM genomic_samples WHERE genetic_sample_id='S_MONO_2X'").fetchone()[0], "2x")
            finally:
                conn.close()

    @unittest.skipIf(os.environ.get("PLANTMATICS_FAST_AUDIT") == "1" or os.environ.get("PACKERA_FAST_AUDIT") == "1" or shutil.which("Rscript") is None, "Skipping R test in fast audit or without Rscript")
    def test_step07_cytology_integration_decision_matrix_and_propagation_r(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            csv_out = Path(tmp) / "cytotype_concordance_summary.csv"
            empty_cyto = Path(tmp) / "empty_cyto"
            empty_cyto.mkdir()
            r_code = f"""
source("scripts/analysis/step07_cytology_integration.R")
source("scripts/analysis/step07_genomic_propagation.R")
source("scripts/analysis/genomic_triage_utils.R")
source("scripts/analysis/07_triage_dashboard_synthesis.R")
make_cohort <- function(overlap = FALSE) {{
  n <- 12L
  df <- data.frame(
    catalogNumber = c(paste0("D1_", 1:5), paste0("D2_", 1:5), "AU_1", "AU_2"),
    plant_individual_id = rep(0L, n),
    species_standardized = c(rep("Packera dubia", 10), rep("Packera aurea", 2)),
    raw_species_name = c(rep("Packera dubia", 10), rep("Packera aurea", 2)),
    authority_tier = c("Tier_0_Type", rep("Tier_1_Gold", 9), "Tier_0_Type", "Tier_1_Gold"),
    determiner_tier = c("Tier_0_Type", rep("Tier_1_Gold", 9), "Tier_0_Type", "Tier_1_Gold"),
    type_status = c("holotype", rep("", 9), "holotype", ""),
    is_type_specimen = c(TRUE, rep(FALSE, 9), TRUE, FALSE),
    recordedBy = c("Smith", "Smith", "Jones", "Lee", "Brown", "Alpha", "Alpha", "Beta", "Gamma", "Delta", "Gray", "Torrey"),
    institutionCode = c("NCU", "DUKE", "GH", "NY", "US", "NCU", "GA", "TENN", "VPI", "FLAS", "GH", "NY"),
    exsiccatae_key = c("k1", "k2", "k3", "k4", "k5", "m1", "m2", "m3", "m4", "m5", "a1", "a2"),
    population_deme_id = c(rep(1L, 5), 2L, 2L, 3L, 4L, 4L, 10L, 11L),
    decimalLatitude = c(rep(35.9, 5), rep(35.95, 5), 37.5, 37.6),
    decimalLongitude = c(rep(-79.0, 5), rep(-79.05, 5), -78.5, -78.6),
    stateProvince = c(rep("NC", 10), "VA", "VA"), county = c(rep("Orange", 10), "Albemarle", "Augusta"),
    morph_tier_source = rep("Basal_Rosette", n), leaf_zone_confidence = rep(0.95, n),
    gmm_cluster = c(rep(1L, 5), rep(2L, 5), 1L, 1L), gmm_bayes_factor = c(rep(12.0, 5), rep(18.0, 5), 12.0, 12.0),
    is_typical_morphotype = c(rep(TRUE, 5), rep(FALSE, 5), TRUE, TRUE),
    self_morph_mahal_d2 = c(rep(1.5, 5), rep(ifelse(overlap, 2.0, 26.0), 5), 1.0, 1.5), self_morph_df = rep(4L, n),
    min_congener_mahal_d2 = c(rep(15.0, 5), rep(ifelse(overlap, 2.0, 30.0), 5), 1.0, 1.2), open_set_df = rep(4L, n),
    open_set_novelty_p = c(rep(0.2, 5), rep(ifelse(overlap, 0.7, 0.001), 5), 0.8, 0.7),
    type_morph_mahal_d2 = c(rep(1.0, 5), rep(ifelse(overlap, 2.0, 28.0), 5), 1.0, 1.2),
    morphometric_outlier_flag = c(rep(FALSE, 5), rep(!overlap, 5), FALSE, FALSE),
    cda_predicted_taxon = c(rep("Packera dubia", 5), rep(ifelse(overlap, "Packera dubia", "Packera aurea"), 5), "Packera aurea", "Packera aurea"),
    cda_posterior_prob = c(rep(0.95, 5), rep(ifelse(overlap, 0.72, 0.38), 5), 0.96, 0.94),
    edaphic_best_fit_taxon = c(rep("Packera dubia", 10), "Packera aurea", "Packera aurea"),
    rf_predicted_taxon = c(rep("Packera dubia", 10), "Packera aurea", "Packera aurea"),
    rf_posterior_prob = rep(0.92, n), edaphic_z_score = rep(0.4, n), type_edaphic_z_score = rep(0.4, n),
    hopkins_adjusted_doy = c(rep(120, 5), rep(124, 5), 160, 162), delta_doy_resid = c(rep(1.0, 5), rep(4.0, 5), 0.5, 1.0),
    is_active_anthesis = rep(TRUE, n), area_mm2 = c(rep(320, 5), rep(510, 5), 400, 410),
    PC1 = c(rnorm(5, 0, 0.1), rnorm(5, ifelse(overlap, 0.05, 5.0), 0.1), -4.0, -4.1),
    PC2 = c(rnorm(5, 0, 0.1), rnorm(5, ifelse(overlap, 0.05, 5.0), 0.1), -4.0, -4.1), stringsAsFactors = FALSE
  )
  rec <- data.frame(catalogNumber = c("D1_1", "D2_1", "D2_2", "AU_1"), species_standardized = c("Packera dubia", "Packera dubia", "Packera dubia", "Packera aurea"), sporophytic_2n = c(44L, 88L, 88L, 46L), ploidy_level = c(2L, 4L, 4L, 2L), base_x = c(22L, 22L, 22L, 23L), cytology_linkage_tier = rep("C1_Direct_Cytovoucher", 4), stringsAsFactors = FALSE)
  tx <- data.frame(species_standardized = c("Packera dubia", "Packera aurea"), reported_cytotypes = c("2x, 4x", "2x"), modal_ploidy_level = c(2L, 2L), is_multi_cytotype_taxon = c(1L, 0L), stringsAsFactors = FALSE)
  list(df = df, rec = rec, tx = tx)
}}
set.seed(42)
c1 <- make_cohort(FALSE)
r1 <- run_step07_cytology_integration(c1$df, opts = list(cytological_records_df = c1$rec, taxon_cytology_summary_df = c1$tx, cytology_summary_path = "{csv_out.as_posix()}"))
stopifnot(all(r1$cytotype_reproductive_barrier_pass[6:10] == 1L), file.exists("{csv_out.as_posix()}"))
d1 <- apply_taxonomic_decision_matrix(r1, intraspecific_mode = FALSE)
stopifnot(all(d1$taxonomic_treatment_subtype[6:10] == "Putative_Cryptic_Species"), all(grepl("[Postzygotic Heteroploid Barrier:", d1$synthesis_rationale[6:10], fixed = TRUE)))

c2 <- make_cohort(TRUE)
r2 <- run_step07_cytology_integration(c2$df, opts = list(cytological_records_df = c2$rec, taxon_cytology_summary_df = c2$tx))
stopifnot(all(r2$cytotype_reproductive_barrier_pass[6:10] == 1L))
d2 <- apply_taxonomic_decision_matrix(r2, intraspecific_mode = FALSE)
stopifnot(all(d2$taxonomic_treatment_subtype[6:10] == "Subspecies_Ecotype"), all(grepl("[Intraspecific Cytotype Series: FCM/squash verification recommended]", d2$synthesis_rationale[6:10], fixed = TRUE)))

df_p <- data.frame(
  catalogNumber = c("A_PURE", "A_MIX", "REF1", "REF2", "R_SAME", "R_CROSS"),
  plant_individual_id = rep(0L, 6), species_standardized = rep("Packera dubia", 6),
  authority_tier = rep("Tier_1_Gold", 6), type_status = rep("", 6), morph_tier_source = rep("Basal_Rosette", 6),
  decimalLatitude = c(36.500, 35.900, 36.498, 36.503, 36.502, 35.902),
  decimalLongitude = c(-79.050, -79.050, -79.052, -79.048, -79.050, -79.050),
  coordinateUncertaintyInMeters = rep(100, 6), is_azonal_lithology = rep(0L, 6), hopkins_adjusted_doy = rep(120, 6),
  PC1 = c(0.0, 0.0, -0.25, 0.25, 0.01, 0.01), PC2 = c(0.0, 0.0, 0.15, 0.30, -0.01, -0.01),
  geno_linkage_tier = c("G1_Direct_Voucher", "G1_Direct_Voucher", NA, NA, NA, NA),
  geno_quality_pass = c(1L, 1L, NA, NA, NA, NA), geno_predicted_taxon = c("Packera dubia", "Packera dubia", NA, NA, NA, NA),
  geno_posterior_prob = c(0.95, 0.95, NA, NA, NA, NA), voucher_ploidy_level = c(2L, 2L, 2L, 2L, 2L, 4L),
  geno_ploidy_level = c(2L, 2L, NA, NA, NA, NA), is_multi_cytotype_taxon = c(0L, 1L, 0L, 0L, 0L, 1L), stringsAsFactors = FALSE
)
p_res <- propagate_conditional_deme_priors(df_p)
stopifnot(p_res$genomic_propagation_tier[5] == "G2_Conditional_Deme_Hull", is.na(p_res$genomic_propagation_tier[6]))

df_z <- c2$df
rz <- run_step07_cytology_integration(df_z, cytology_dir = "{empty_cyto.as_posix()}")
stopifnot(all(is.na(rz$cytology_linkage_tier)), all(is.na(rz$voucher_ploidy_level)), all(is.na(rz$cytotype_reproductive_barrier_pass)))
dz_a <- apply_taxonomic_decision_matrix(df_z, intraspecific_mode = FALSE)
dz_b <- apply_taxonomic_decision_matrix(rz, intraspecific_mode = FALSE)
stopifnot(identical(dz_a$taxonomic_treatment_subtype, dz_b$taxonomic_treatment_subtype), identical(dz_a$synthesis_rationale, dz_b$synthesis_rationale))
cat("R_CYTO_PASS\\n")
"""
            proc = subprocess.run(["Rscript", "-e", r_code], cwd=str(PROJECT_ROOT), capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, msg=f"R failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
            self.assertIn("R_CYTO_PASS", proc.stdout)


if __name__ == "__main__":
    unittest.main()
