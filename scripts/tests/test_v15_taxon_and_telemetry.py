#!/usr/bin/env python3
"""
Unit and batched-R contract tests for Plantmatics v1.5:
  1. Non-Packera genus end-to-end portability (Solidago & Quercus across Python, SQLite, and R).
  2. Cache-subtracted StepETATracker & JSONL telemetry (`emit_progress` / `emit_r_progress`).
  3. Cross-OS `isal` / `gzip` FASTQ streaming, GCD cytology `base_x`, vectorized coancestry/f-branch parity,
     and Windows `latest_run.json` pointer fallback.
"""

from __future__ import annotations

import gzip, io, json, os, sqlite3, subprocess, sys, tarfile, tempfile, unittest, zipfile
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path: sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core import harvester as harv_mod, step01_biogeography_utils as bio_mod
from scripts.core.config import get_active_run_dir, initialize_run_directory, resolve_target_genus
from scripts.core.db import get_db_connection, load_taxon_synonyms
from scripts.core.progress_telemetry import StepETATracker, emit_progress
from scripts.data_prep import download_public_cytology as dl_cyto, skim_offtarget_organellar as skim_mod, step06_synonym_utils as syn_mod
from scripts.tests._r_session import rscript_available

_R_BATCH_CODE = r"""
args <- commandArgs(trailingOnly = TRUE); r_run_dir <- args[1]; out_json <- args[2]
suppressPackageStartupMessages({ library(jsonlite) })
source("scripts/core/db_utils.R")
source("scripts/analysis/step06_reticulation_coancestry.R")
source("scripts/analysis/step06_fbranch_quartets.R")
source("scripts/analysis/step07_biogeo_hybrid_utils.R")
source("scripts/analysis/step07_organellar_triage_utils.R")

# 1. Non-Packera R helpers + Step 07 retention
Sys.setenv(PLANTMATICS_TARGET_GENUS = "Solidago")
sol_genus <- resolve_r_target_genus(); sol_py <- resolve_r_python_bin()
sol_keep <- unname(is_active_genus_binomial(c("Solidago altissima", "Solidago rugosa", "Packera aurea", UNASSIGNED_TAXON_SENTINEL)))
sol_abbr <- abbreviate_binomial("Solidago altissima")
sol_unass <- unname(is_unassigned_or_unknown_taxon(c(UNASSIGNED_TAXON_SENTINEL, "Unknown", "Solidago altissima")))
sol_s07_valid <- c(.is_valid_packera_taxon_step07("Solidago altissima", character(0)), .is_valid_packera_taxon_step07("Packera aurea", character(0)))
sol_org_valid <- unname(.is_valid_packera_name_r(c("Solidago altissima", "Packera aurea")))
org_out <- evaluate_step07_organellar_and_oop_signals(data.frame(
  catalogNumber = c("S1", "S2"), species_standardized = c("Solidago altissima", "Solidago canadensis"),
  geno_predicted_taxon = c("Solidago altissima", "Solidago canadensis"), geno_posterior_prob = c(0.95, 0.95),
  public_plastome_match = c("Solidago rugosa", "Solidago canadensis"), plastid_kmer_containment = c(0.95, 0.95), plastid_shared_11mers = c(120, 120),
  public_its_match = c("Solidago altissima", "Solidago canadensis"), waxy_lineage_match = c("Solidago altissima", "Solidago canadensis"),
  cytonuclear_discordant = c(1L, 0L), stringsAsFactors = FALSE
))
Sys.setenv(PLANTMATICS_TARGET_GENUS = "Quercus")
que_genus <- resolve_r_target_genus(); que_abbr <- abbreviate_binomial("Quercus alba")
que_keep <- unname(is_active_genus_binomial(c("Quercus alba", "Quercus rubra", "Solidago altissima", "Packera dubia")))
Sys.unsetenv("PLANTMATICS_TARGET_GENUS")

# 2. R progress JSONL emission
Sys.setenv(PLANTMATICS_RUN_DIR = r_run_dir)
ok_p1 <- isTRUE(emit_r_progress(5L, "spatial_rf", 25L, 100L, unit = "vouchers", cached_skipped = 50L, status = "running"))
ok_p2 <- isTRUE(emit_r_progress(5L, "spatial_rf", 50L, 100L, unit = "vouchers", cached_skipped = 50L, status = "completed"))
Sys.unsetenv("PLANTMATICS_RUN_DIR")

# 3A. Vectorized dosage coancestry vs naive reference
set.seed(202611L); sids <- paste0("S", 1:6)
dos_mat <- matrix(sample(c(0.0, 0.5, 1.0), 6L * 50L, replace = TRUE), nrow = 6L, dimnames = list(sids, paste0("L", 1:50)))
dos_mat[1L, 3L] <- NA_real_; dos_mat[4L, 17L] <- NA_real_
empty_dir <- file.path(r_run_dir, "empty_geno"); dir.create(empty_dir, recursive = TRUE, showWarnings = FALSE)
res_cdos <- compute_microhaplotype_coancestry(sids, dos_mat, empty_dir)
C_vec_dos <- as.matrix(res_cdos$coanc_df[, sids, drop = FALSE]); rownames(C_vec_dos) <- res_cdos$coanc_df$genetic_sample_id

naive_dos_coanc <- function(pass_ids, m_dos) {
  n <- length(pass_ids); C_ref <- matrix(0.0, nrow = n, ncol = n, dimnames = list(pass_ids, pass_ids))
  n_loci <- ncol(m_dos); step <- max(5L, min(20L, as.integer(floor(n_loci / 25L))))
  for (st in seq(1L, max(1L, n_loci - step + 1L), by = step)) {
    sub_m <- m_dos[, st:min(n_loci, st + step - 1L), drop = FALSE]
    for (i in seq_len(n)) {
      dists <- vapply(seq_len(n), function(j) { if (i == j) return(Inf); ok <- is.finite(sub_m[i, ]) & is.finite(sub_m[j, ]); if (sum(ok) < 2L) Inf else mean(abs(sub_m[i, ok] - sub_m[j, ok])) }, numeric(1))
      fin <- which(is.finite(dists))
      if (length(fin) > 0L) { donors <- fin[dists[fin] <= min(dists[fin]) + 0.05]; C_ref[i, donors] <- C_ref[i, donors] + (1.0 / length(donors)) }
    }
  }
  rs <- rowSums(C_ref); for (i in seq_len(n)) if (rs[i] > 0) C_ref[i, ] <- C_ref[i, ] / rs[i]; round(C_ref, 5)
}
coanc_dos_max_diff <- max(abs(C_vec_dos - naive_dos_coanc(sids, dos_mat)))

# 3B. Vectorized 122bp PHYLIP coancestry vs naive reference
phy_dir <- file.path(r_run_dir, "phy_geno"); dir.create(phy_dir, recursive = TRUE, showWarnings = FALSE)
bases <- c("A", "C", "G", "T")
blks <- lapply(1:4, function(bi) replicate(6L, paste(sample(bases, 122L, replace = TRUE), collapse = "")))
blks[[1L]][1L] <- paste0("NNNN-", substr(blks[[1L]][1L], 6L, 122L))
phy_seqs <- do.call(paste0, blks)
writeLines(c(sprintf("6 %d", nchar(phy_seqs[1L])), sprintf("%-10s%s", sids, phy_seqs)), file.path(phy_dir, "rad_loci.phy"))
res_cphy <- compute_microhaplotype_coancestry(sids, dos_mat, phy_dir)
C_vec_phy <- as.matrix(res_cphy$coanc_df[, sids, drop = FALSE]); rownames(C_vec_phy) <- res_cphy$coanc_df$genetic_sample_id

naive_phy_coanc <- function(pass_ids, sqs) {
  n <- length(pass_ids); C_ref <- matrix(0.0, nrow = n, ncol = n, dimnames = list(pass_ids, pass_ids)); raw_mats <- lapply(sqs, charToRaw)
  for (b in 1:4) {
    st <- (b - 1L) * 122L + 1L; en <- b * 122L; sub_list <- lapply(raw_mats, function(r) r[st:en]); d_blk <- matrix(Inf, nrow = n, ncol = n)
    for (a in seq_len(n)) {
      ra <- sub_list[[a]]; ok_a <- ra != as.raw(0x4e) & ra != as.raw(0x2d)
      for (b2 in seq_len(n)) { if (a == b2) next; rb <- sub_list[[b2]]; ok <- ok_a & (rb != as.raw(0x4e)) & (rb != as.raw(0x2d)); if (sum(ok) >= 10L) d_blk[a, b2] <- sum(ra[ok] != rb[ok]) / sum(ok) }
      fin_b <- which(is.finite(d_blk[a, ]))
      if (length(fin_b) > 0L) { donors <- fin_b[d_blk[a, fin_b] <= min(d_blk[a, fin_b]) + 1.0 / 122.0 + 1e-9]; C_ref[a, donors] <- C_ref[a, donors] + (1.0 / length(donors)) }
    }
  }
  rs <- rowSums(C_ref); for (i in seq_len(n)) if (rs[i] > 0) C_ref[i, ] <- C_ref[i, ] / rs[i]; round(C_ref, 5)
}
coanc_phy_max_diff <- max(abs(C_vec_phy - naive_phy_coanc(sids, phy_seqs)))

# 3C. Vectorized f-branch jackknife vs naive block loop
fb_sids <- paste0("F", 1:8)
fb_df <- data.frame(genetic_sample_id = fb_sids, species_standardized = rep(c("Solidago altissima", "Solidago rugosa", "Solidago canadensis", "Solidago nemoralis"), each = 2L), is_putative_novel_or_hybrid_query = rep(0L, 8L), population_id = rep("P1", 8L), stringsAsFactors = FALSE)
fb_dos <- matrix(runif(8L * 60L), nrow = 8L, dimnames = list(fb_sids, paste0("FL", 1:60))); fb_dos[1L, 5L] <- NA_real_
fb_tbl <- compute_tree_guided_fbranch(fb_df, fb_dos, empty_dir)$fbranch_df
uniq_tx <- sort(unique(fb_df$species_standardized))
p_mat <- do.call(rbind, lapply(uniq_tx, function(tx) colMeans(fb_dos[fb_df$species_standardized == tx, , drop = FALSE], na.rm = TRUE))); rownames(p_mat) <- uniq_tx
blk_id <- cut(seq_len(ncol(fb_dos)), breaks = min(200L, max(10L, as.integer(floor(ncol(fb_dos) / 5L)))), labels = FALSE)
naive_diffs <- vapply(seq_len(nrow(fb_tbl)), function(ri) {
  r <- fb_tbl[ri, ]; p1 <- r$sister_p1; p2 <- r$recipient_taxon; p3 <- r$donor_taxon; og <- r$outgroup_o
  ok <- is.finite(p_mat[p1, ]) & is.finite(p_mat[p2, ]) & is.finite(p_mat[p3, ]) & is.finite(p_mat[og, ])
  fa <- p_mat[p1, ok]; fb <- p_mat[p2, ok]; fc <- p_mat[p3, ok]; fd <- p_mat[og, ok]
  abba <- (1.0 - fa) * fb * fc * (1.0 - fd); baba <- fa * (1.0 - fb) * fc * (1.0 - fd); den_d <- abba + baba
  d_full <- (sum(abba) - sum(baba)) / sum(den_d); b_sub <- blk_id[ok]
  jk_d <- vapply(unique(b_sub), function(bk) { k <- b_sub != bk; sd_k <- sum(den_d[k]); if (sd_k > 1e-8) (sum(abba[k]) - sum(baba[k])) / sd_k else d_full }, numeric(1))
  z_ref <- abs(d_full / sqrt(max(1e-10, ((length(jk_d) - 1.0) / length(jk_d)) * sum((jk_d - mean(jk_d))^2))))
  c(abs(round(abs(d_full), 4) - r$d_stat), abs(round(z_ref, 3) - r$jackknife_z))
}, numeric(2))

out <- list(
  sol_genus = sol_genus, sol_py = sol_py, sol_keep = sol_keep, sol_abbr = sol_abbr, sol_unass = sol_unass,
  sol_s07_valid = sol_s07_valid, sol_org_valid = sol_org_valid, org_cap = unname(org_out$ref_cytonuc_capture),
  que_genus = que_genus, que_keep = que_keep, que_abbr = que_abbr, ok_p1 = ok_p1, ok_p2 = ok_p2,
  coanc_dos_mode = res_cdos$mode, coanc_dos_max_diff = coanc_dos_max_diff,
  coanc_phy_mode = res_cphy$mode, coanc_phy_max_diff = coanc_phy_max_diff,
  fb_rows = nrow(fb_tbl), fb_d_max_diff = max(naive_diffs[1L, ]), fb_z_max_diff = max(naive_diffs[2L, ])
)
writeLines(toJSON(out, auto_unbox = TRUE), out_json)
"""


class TestV15TaxonAndTelemetry(unittest.TestCase):
    """Contract checks for genus portability, telemetry ETAs, fast gzip, GCD cytology, R parity, and Windows pointers."""

    r_res: dict = {}
    r_tmp_dir: tempfile.TemporaryDirectory[str] | None = None

    @classmethod
    def setUpClass(cls) -> None:
        if not rscript_available():
            return
        cls.r_tmp_dir = tempfile.TemporaryDirectory(prefix="plantmatics_v15_r_")
        tmp_p = Path(cls.r_tmp_dir.name)
        (tmp_p / "r_run").mkdir(parents=True, exist_ok=True)
        r_script, out_json = tmp_p / "batch_v15.R", tmp_p / "batch_v15.json"
        r_script.write_text(_R_BATCH_CODE, encoding="utf-8")
        proc = subprocess.run(["Rscript", "--vanilla", str(r_script), str(tmp_p / "r_run"), str(out_json)], cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=15)
        if proc.returncode != 0:
            raise AssertionError(f"Batched R check failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls.r_res = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.r_tmp_dir is not None:
            cls.r_tmp_dir.cleanup(); cls.r_tmp_dir = None

    def test_1_non_packera_genus_end_to_end(self) -> None:
        """Verify Solidago & Quercus portability across Python, SQLite taxon_synonyms, and R Step 04/05/07."""
        with patch.dict(harv_mod._DYNAMIC_SYNONYM_TO_CANONICAL, {}, clear=False), patch.dict(harv_mod._DYNAMIC_TAXON_AUTHOR_STRINGS, {}, clear=False):
            for genus, sp_a, syn_a, sp_b in (
                ("Solidago", "Solidago altissima", "Solidago canadensis var. scabra", "Solidago rugosa"),
                ("Quercus", "Quercus alba", "Quercus alba var. latiloba", "Quercus rubra"),
            ):
                with patch.dict(os.environ, {"PLANTMATICS_TARGET_GENUS": genus}):
                    self.assertEqual(resolve_target_genus(), genus)
                    syns_g = bio_mod.resolve_synonym_genera()
                    self.assertIn(genus, syns_g); self.assertNotIn("Packera", syns_g); self.assertNotIn("Senecio", syns_g)
                    harv_mod._register_dynamic_synonyms(sp_a, [syn_a], ["L."], target_genus=genus)
                    self.assertEqual(harv_mod._DYNAMIC_SYNONYM_TO_CANONICAL.get(syn_a.lower()), sp_a)
                    self.assertEqual(syn_mod.canonical_taxon_equiv_key(sp_a.replace(" ", "_")), sp_a.lower())
                    self.assertEqual(syn_mod.canonical_taxon_equiv_key("Senecio tomentosus"), "senecio tomentosus")
                    self.assertEqual(syn_mod.normalize_step06_taxon_name(sp_b), (sp_b, sp_b, False))
                    self.assertEqual(syn_mod.normalize_public_header_binomial(f"NC1|{sp_a.replace(' ', '_')}|rbcL"), sp_a.replace(" ", "_"))
                    self.assertEqual(syn_mod.normalize_public_header_binomial("Senecio_aureus"), f"{genus}_sp")
                    self.assertEqual(syn_mod.extract_species_from_seed_key(f"MK1|{sp_b.replace(' ', '_')}|ITS"), sp_b)

            with tempfile.TemporaryDirectory() as td:
                conn = get_db_connection(Path(td) / "non_packera.sqlite", init_schema=True)
                try:
                    with conn:
                        conn.executemany("INSERT INTO taxon_synonyms (raw_name_lower, canonical_taxon) VALUES (?, ?);", [("solidago canadensis var. scabra", "Solidago altissima"), ("quercus alba var. latiloba", "Quercus alba")])
                    loaded = load_taxon_synonyms(conn)
                    self.assertEqual(loaded["solidago canadensis var. scabra"], "Solidago altissima")
                    self.assertEqual(loaded["quercus alba var. latiloba"], "Quercus alba")
                finally:
                    conn.close()

        for rel in ("scripts/morphometrics/04_gmm_morphotools.R", "scripts/analysis/05_multimodal_spatial_rf.R", "scripts/analysis/step07_biogeo_hybrid_utils.R"):
            self.assertIn("is_active_genus_binomial", (PROJECT_ROOT / rel).read_text(encoding="utf-8"))

        if rscript_available():
            rr = self.r_res
            self.assertEqual(rr["sol_genus"], "Solidago"); self.assertTrue(bool(rr["sol_py"]))
            self.assertEqual(rr["sol_keep"], [True, True, False, False]); self.assertEqual(rr["sol_abbr"], "S. altissima")
            self.assertEqual(rr["sol_unass"], [True, True, False]); self.assertEqual(rr["sol_s07_valid"], [True, False])
            self.assertEqual(rr["sol_org_valid"], [True, False]); self.assertEqual(rr["org_cap"], [True, False])
            self.assertEqual(rr["que_genus"], "Quercus"); self.assertEqual(rr["que_keep"], [True, True, False, False])
            self.assertEqual(rr["que_abbr"], "Q. alba")

    def test_2_cache_subtracted_eta_and_jsonl_telemetry(self) -> None:
        """Verify cache-subtracted n_active, EWMA ETA ignoring dt=0 cached items, SQLite benchmarks, and R/Py JSONL."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td); run_dir, db_p = tmp / "run_02", tmp / "telemetry.sqlite"
            get_db_connection(db_p, init_schema=True).close()

            clk = [100.0]
            trk = StepETATracker(step=2, substage="seg", total=500, cached_skipped=400, unit="vouchers", run_dir=run_dir, db_path=db_p, clock=lambda: clk[0])
            self.assertEqual(trk.n_active, 100)
            self.assertIsNone(trk.update(done=0, status="running")["eta_sec"])
            for done_i in (1, 2, 3):
                clk[0] += 1.0; ev = trk.update(done=done_i, status="running")
            self.assertAlmostEqual(float(ev["rate_items_per_sec"]), 1.0, places=3)
            self.assertAlmostEqual(float(ev["eta_sec"]), 97.0, places=3)
            self.assertAlmostEqual(float(trk.update(done=4, dt=0.0, status="running")["rate_items_per_sec"]), 1.0, places=3)

            clk[0] += 96.0; ev_done = trk.finish()
            self.assertEqual((ev_done["status"], ev_done["eta_sec"], ev_done["progress_frac"]), ("completed", 0.0, 1.0))
            with sqlite3.connect(db_p) as c2:
                bm = c2.execute("SELECT mean_items_per_sec, n_runs FROM step_telemetry_benchmarks WHERE step_id=2 AND substage='seg';").fetchone()
            self.assertIsNotNone(bm); self.assertGreater(float(bm[0]), 0.0); self.assertEqual(int(bm[1]), 1)
            self.assertGreater(float(trk.estimate_step_eta(2, 50, substage="seg") or 0.0), 0.0)

            self.assertEqual(emit_progress(6, "genomics", 0, 0, status="skipped_cached", run_dir=run_dir)["eta_sec"], 0.0)
            self.assertEqual(trk.estimate_step_eta(6, 0), 0.0)
            lines = [json.loads(ln) for ln in (run_dir / "logs" / "progress.jsonl").read_text(encoding="utf-8").splitlines() if ln.strip()]
            self.assertGreaterEqual(len(lines), 6); self.assertEqual(lines[-1]["substage"], "genomics")

        if rscript_available():
            self.assertTrue(self.r_res["ok_p1"]); self.assertTrue(self.r_res["ok_p2"])
            assert self.r_tmp_dir is not None
            r_lines = [json.loads(ln) for ln in (Path(self.r_tmp_dir.name) / "r_run" / "logs" / "progress.jsonl").read_text(encoding="utf-8").splitlines()]
            self.assertEqual((r_lines[0]["n_active"], r_lines[0]["progress_frac"]), (50, 0.5))
            self.assertEqual((r_lines[1]["eta_sec"], r_lines[1]["progress_frac"], r_lines[1]["status"]), (0.0, 1.0, "completed"))

    def test_3_fast_gzip_cytology_gcd_coancestry_and_windows_pointer(self) -> None:
        """Verify gzip/isal streaming without subprocesses, GCD base_x=9 inference, R parity, and Windows pointer."""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td); gz_bytes = gzip.compress(b"@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")
            fq_gz, zp, tp = tmp / "S1_R1.fastq.gz", tmp / "reads.zip", tmp / "reads.tar.gz"
            fq_gz.write_bytes(gz_bytes)
            with zipfile.ZipFile(zp, "w") as zf: zf.writestr("sub/S2_R1.fastq.gz", gz_bytes)
            with tarfile.open(tp, "w:gz") as tf:
                ti = tarfile.TarInfo("sub/S3_R1.fastq.gz"); ti.size = len(gz_bytes); tf.addfile(ti, io.BytesIO(gz_bytes))
            with patch("subprocess.Popen", side_effect=AssertionError("External subprocess disallowed")):
                self.assertIn(("S1", "ACGTACGTACGTACGT"), list(skim_mod._iter_staged_fastq_reads(tmp)))
                self.assertEqual(list(skim_mod._iter_staged_fastq_reads(zp)), [("S2", "ACGTACGTACGTACGT")])
                self.assertEqual(list(skim_mod._iter_staged_fastq_reads(tp)), [("S3", "ACGTACGTACGTACGT")])
            fake_igzip = type("FakeIgzip", (), {"open": staticmethod(gzip.open)})
            with patch.object(skim_mod, "_fast_gzip", fake_igzip), patch("subprocess.Popen", side_effect=AssertionError("External subprocess disallowed")):
                self.assertEqual(list(skim_mod._iter_fastq_file_seqs(fq_gz)), ["ACGTACGTACGTACGT"])

            with patch.dict(os.environ, {"PLANTMATICS_TARGET_GENUS": "Solidago"}), patch.dict(dl_cyto._DB_SYN_CACHE, {}, clear=True):
                self.assertEqual(dl_cyto._infer_modal_base_x([18, 36, 54]), 9)
                raw_recs = [{**c, "cytology_id": f"CY_{i}", "source_db": "CCDB", "species_standardized": "Solidago altissima"} for i, val in enumerate(("2n=18", "2n=36", "2n=54", "9II + 9I")) for c in dl_cyto.parse_packera_chromosome_record(val)]
                linked, sums = dl_cyto._cross_link_and_summarize(raw_recs, {"by_cat": {}, "by_exsic": {}, "by_cn": {}, "gb_to_cat": {}})
                self.assertEqual([(r["base_x"], r["ploidy_level"], r["cytotype_label"]) for r in linked], [(9, 2, "2x"), (9, 4, "4x"), (9, 6, "6x"), (9, 3, "3x")])
                self.assertEqual((sums[0]["modal_base_x"], sums[0]["reported_cytotypes"]), (9, "2x, 3x, 4x, 6x"))

            ws = tmp / "win_ws"; ws.mkdir()
            clean_env = {k: v for k, v in os.environ.items() if k not in ("PLANTMATICS_RUN_DIR", "PACKERA_RUN_DIR", "RUN_ID")}
            with patch.dict(os.environ, clean_env, clear=True), patch.object(Path, "symlink_to", side_effect=OSError("WinError 1314")):
                rd = initialize_run_directory(run_id="win_ptr_test", workspace_root=ws, update_latest_symlink=True, export_env=False)
                ptr_file = ws / "runs" / "latest_run.json"
                self.assertTrue(ptr_file.is_file()); self.assertFalse((ws / "runs" / "latest").exists())
                self.assertEqual(json.loads(ptr_file.read_text(encoding="utf-8"))["run_id"], "win_ptr_test")
                self.assertEqual(get_active_run_dir(workspace_root=ws, fallback_to_latest=True).resolve(), rd.resolve())
                self.assertEqual(get_active_run_dir(run_id="latest", workspace_root=ws).resolve(), rd.resolve())

        if rscript_available():
            rr = self.r_res
            self.assertEqual(rr["coanc_dos_mode"], "dosage_window_haplotypes"); self.assertLess(float(rr["coanc_dos_max_diff"]), 1e-4)
            self.assertEqual(rr["coanc_phy_mode"], "122bp_phylip_microhaplotypes"); self.assertLess(float(rr["coanc_phy_max_diff"]), 1e-4)
            self.assertGreater(int(rr["fb_rows"]), 0); self.assertLess(float(rr["fb_d_max_diff"]), 1e-4); self.assertLess(float(rr["fb_z_max_diff"]), 1e-3)


if __name__ == "__main__":
    unittest.main()
