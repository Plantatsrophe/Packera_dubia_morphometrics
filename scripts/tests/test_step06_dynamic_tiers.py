#!/usr/bin/env python3
"""Consolidated warm-R & Python tests for Step 06 dynamic 4-tier routing, authority gates, and 3-plate PNG rendering."""
from __future__ import annotations

import csv, importlib, json, os, shutil, sqlite3, subprocess, sys, tempfile, textwrap, unittest
from pathlib import Path
from typing import Any, Dict, List
from unittest import mock

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db_schema import init_db_schema
from scripts.pipeline.resolve_focal_species_complex import _extract_stage_metrics

ING = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _seq(seed_int: int, length: int = 240, mut_from: str = "", mut_n: int = 0) -> str:
    rng = np.random.default_rng(seed_int)
    bases = np.array(list("ACGT"))
    if mut_from:
        arr = np.array(list(mut_from))
        idx = rng.choice(len(arr), size=min(mut_n, len(arr)), replace=False)
        for i in idx:
            arr[i] = rng.choice([b for b in "ACGT" if b != arr[i]])
        return "".join(arr.tolist())
    return "".join(rng.choice(bases, size=length).tolist())


def _mk_vcf(samples: List[str], n_loci: int = 40) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for k in range(n_loci):
        cells = []
        for idx_s, s in enumerate(samples):
            if s.startswith("PACDUB"):
                gt, r, a = ("0/1", 10, 10) if (k + idx_s) % 3 == 0 else (("0/0", 20, 0) if (k + idx_s) % 5 != 0 else ("1/1", 0, 20))
            elif s.startswith("PACANO"):
                gt, r, a = ("1/1", 0, 20) if k % 2 == 0 else ("0/1", 10, 10)
            elif s.startswith("PACOBO"):
                gt, r, a = ("0/1", 12, 12) if k % 3 == 0 else ("1/1", 0, 24)
            else:
                gt, r, a = ("1/1", 0, 20) if k % 4 != 0 else ("0/0", 20, 0)
            cells.append(f"{gt}:{r + a}:{r},{a}")
        lines.append(f"chr1\t{100 + k}\tL_{k}\tA\tG\t99\tPASS\t.\tGT:DP:AD\t" + "\t".join(cells))
    return "\n".join(lines) + "\n"


def _seed_sandbox(tmp: Path) -> Dict[str, Path]:
    geno, db_p, v_csv = tmp / "genomics", tmp / "db.sqlite", tmp / "curated_vouchers.csv"
    (geno / "snpsaurus_nextrad").mkdir(parents=True)
    (geno / "public_nuclear_backbone").mkdir(parents=True)
    (geno / "public_organellar_and_loci").mkdir(parents=True)
    (geno / "public_phytop_conflict").mkdir(parents=True)
    (geno / "public_reltime_trees").mkdir(parents=True)
    (geno / "public_catalogs").mkdir(parents=True)

    sids = ["PACDUB1M", "PACDUB2M", "PACDUB3M", "PACDUB4M", "PACANO1M", "PACANO2M", "PACANO3M", "PACOBO1M", "EREHIE1M"]
    (geno / "snpsaurus_nextrad" / "cohort.vcf").write_text(_mk_vcf(sids), encoding="utf-8")
    with open(geno / "snpsaurus_nextrad" / "user_genetic_metadata.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genetic_sample_id", "species_standardized", "decimalLatitude", "decimalLongitude", "stateProvince"])
        for s in sids:
            sp = "Packera dubia" if "DUB" in s else ("Packera anonyma" if "ANO" in s else ("Packera obovata" if "OBO" in s else "Erechtites hieraciifolius"))
            w.writerow([s, sp, "35.5", "-79.2", "NC"])

    # Dryad Hyb-Seq (includes dubia, obovata, anonyma, aurea, Erechtites; EXCLUDES crawfordii, hartiana)
    mo_tre = "(Erechtites_hieraciifolius:0.15,((Packera_dubia:0.03,Packera_anonyma:0.04)0.95:0.02,(Packera_obovata:0.03,Packera_aurea:0.03)0.92:0.02)0.98:0.05);"
    (geno / "public_nuclear_backbone" / "public_packera_MO_Astral_lpp.tre").write_text(mo_tre + "\n", encoding="utf-8")
    (geno / "public_nuclear_backbone" / "public_packera_ASTRALPro_Astral_lpp.tre").write_text(mo_tre + "\n", encoding="utf-8")
    (geno / "public_organellar_and_loci" / "public_packera_111mapping_plastomes.tre").write_text(mo_tre + "\n", encoding="utf-8")
    (geno / "public_phytop_conflict" / "public_packera_PhyTop_MO.tsv").write_text(
        "taxon\tILS_index\tIH_index\nPackera_dubia\t0.32\t0.11\nPackera_obovata\t0.41\t0.18\nPackera_anonyma\t0.28\t0.09\nPackera_aurea\t0.25\t0.05\n", encoding="utf-8"
    )
    (geno / "public_reltime_trees" / "public_packera_RelTime_MO_nexus.tre").write_text(
        f"#NEXUS\nBegin Trees;\n  tree rt = {mo_tre}\nEnd;\n", encoding="utf-8"
    )
    base_h = _seq(101, 400)
    (geno / "public_nuclear_backbone" / "public_packera_MO_concat_best.fasta").write_text(
        f">Packera_dubia\n{_seq(102, 400, base_h, 15)}\n>Packera_obovata\n{_seq(103, 400, base_h, 18)}\n"
        f">Packera_anonyma\n{_seq(104, 400, base_h, 22)}\n>Packera_aurea\n{_seq(105, 400, base_h, 20)}\n"
        f">Erechtites_hieraciifolius\n{_seq(106, 400, base_h, 90)}\n", encoding="utf-8"
    )
    with open(geno / "public_catalogs" / "public_packera_sra_runs.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["bioproject", "run_accession", "sample_accession", "scientific_name", "library_name", "specimen_voucher"])
        w.writerow(["PRJNA907383", "SRR_OBO_1", "SAMN1", "Packera obovata", "lib1", "OBO_VOUCH_1"])
        w.writerow(["PRJNA907383", "SRR_DUB_1", "SAMN2", "Packera dubia", "lib2", "DUB_VOUCH_1"])

    # GenBank / BOLD FASTAs (includes crawfordii, dubia, obovata, anonyma; EXCLUDES hartiana)
    base_its = _seq(201, 240)
    (geno / "public_organellar_and_loci" / "public_packera_its_ets.fasta").write_text(
        f">ACC_CRAW1|Packera_crawfordii Packera crawfordii voucher CRAW_VOUCH_1 internal transcribed spacer\n{_seq(202, 240, base_its, 8)}\n"
        f">ACC_CRAW2|Packera_crawfordii Packera crawfordii internal transcribed spacer\n{_seq(203, 240, base_its, 9)}\n"
        f">ACC_DUB1|Packera_dubia Packera dubia internal transcribed spacer\n{_seq(204, 240, base_its, 24)}\n"
        f">ACC_OBO1|Packera_obovata Packera obovata internal transcribed spacer\n{_seq(205, 240, base_its, 26)}\n"
        f">ACC_ANO1|Packera_anonyma Packera anonyma internal transcribed spacer\n{_seq(206, 240, base_its, 28)}\n", encoding="utf-8"
    )
    base_cp = _seq(301, 240)
    (geno / "public_organellar_and_loci" / "public_packera_plastomes_and_cpdna.fasta").write_text(
        f">CP_CRAW1|Packera_crawfordii Packera crawfordii voucher CRAW_VOUCH_1 chloroplast\n{_seq(302, 240, base_cp, 6)}\n"
        f">CP_DUB1|Packera_dubia Packera dubia chloroplast\n{_seq(303, 240, base_cp, 18)}\n"
        f">CP_OBO1|Packera_obovata Packera obovata chloroplast\n{_seq(304, 240, base_cp, 20)}\n"
        f">CP_ANO1|Packera_anonyma Packera anonyma chloroplast\n{_seq(305, 240, base_cp, 22)}\n", encoding="utf-8"
    )
    base_wx = _seq(401, 240)
    (geno / "public_organellar_and_loci" / "public_packera_waxy_nuclear.fasta").write_text(
        f">WX_CRAW1|Packera_crawfordii Packera crawfordii voucher CRAW_VOUCH_1 waxy gene\n{_seq(402, 240, base_wx, 8)}\n"
        f">WX_CRAW2|Packera_crawfordii Packera crawfordii waxy gene\n{_seq(403, 240, base_wx, 9)}\n"
        f">WX_DUB1|Packera_dubia Packera dubia waxy gene\n{_seq(404, 240, base_wx, 24)}\n"
        f">WX_ANO1|Packera_anonyma Packera anonyma waxy gene\n{_seq(405, 240, base_wx, 26)}\n", encoding="utf-8"
    )

    v_rows = [
        ("OBO_VOUCH_1", "Packera obovata", "Tier_1_Gold", 2010, 35.6, -79.1),
        ("CRAW_VOUCH_1", "Packera crawfordii", "Tier_2_Silver", 2012, 35.7, -79.2),
        ("DUB_VOUCH_1", "Packera dubia", "Tier_1_Gold", 2015, 35.5, -79.0),
        ("HART_VOUCH_1", "Packera hartiana", "Tier_1_Gold", 2018, 35.2, -111.6),
        ("ANO_VOUCH_1", "Packera anonyma", "Tier_1_Gold", 2014, 35.8, -79.4),
    ]
    with open(v_csv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["catalogNumber", "species_standardized", "determiner_tier", "year", "decimalLatitude", "decimalLongitude", "institutionCode"])
        for cat, sp, dt, yr, la, lo in v_rows:
            w.writerow([cat, sp, dt, yr, la, lo, "NCU"])

    conn = sqlite3.connect(str(db_p))
    init_db_schema(conn)
    with conn:
        for cat, sp, dt, yr, la, lo in v_rows:
            conn.execute("INSERT INTO vouchers (catalogNumber, species_standardized, scientificName, determiner_tier, decimalLatitude, decimalLongitude, institutionCode) VALUES (?, ?, ?, ?, ?, ?, 'NCU');", (cat, sp, sp, dt, la, lo))
        for sp, n2, pl in (("Packera dubia", 46, 2), ("Packera obovata", 46, 2), ("Packera crawfordii", 46, 2), ("Packera hartiana", 92, 4), ("Packera anonyma", 46, 2)):
            conn.execute("INSERT INTO taxon_cytology_summary (species_standardized, modal_sporophytic_2n, modal_ploidy_level, reported_cytotypes, is_multi_cytotype_taxon) VALUES (?, ?, ?, ?, 0);", (sp, n2, pl, f"{pl}x"))
            conn.execute("INSERT INTO cytological_records (cytology_id, source_db, original_taxon_name, species_standardized, sporophytic_2n, ploidy_level, cytology_linkage_tier) VALUES (?, 'CCDB', ?, ?, ?, ?, 'C1_Direct_Cytovoucher');", (f"CYT_{sp.replace(' ', '_')}", sp, sp, n2, pl))
    conn.close()
    return {"geno": geno, "db": db_p, "vouchers": v_csv}


class TestStep06DynamicTiers(unittest.TestCase):
    """Test Suite for Step 06 Dynamic 4-Tier Phylogenomic Routing & Figure Generation."""

    tmp_dir: tempfile.TemporaryDirectory
    paths: Dict[str, Path]
    r_summary: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp_dir = tempfile.TemporaryDirectory(prefix="step06_dyn_")
        cls.tmp = Path(cls.tmp_dir.name)
        cls.paths = _seed_sandbox(cls.tmp)
        for tag, tx, intra in (("t1_intra", "Packera dubia", True), ("t2_hyb", "Packera obovata", False), ("t3_sanger", "Packera crawfordii", False), ("t4_gap", "Packera hartiana", False)):
            ING.run_genomic_ingestion(
                cls.paths["geno"], cls.tmp / f"parsed_{tag}", db_path=cls.paths["db"], min_hq_reads=0,
                vouchers_csv=cls.paths["vouchers"], focal_taxon=tx, intraspecific_mode=intra, genomic_tier="auto",
            )
        if SKIP_R_COND:
            return
        r_script = cls.tmp / "run_batched_step06.R"
        out_json = cls.tmp / "r_summary.json"
        r_script.write_text(textwrap.dedent(f"""
        if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
        suppressPackageStartupMessages({{ library(jsonlite) }})
        source("scripts/analysis/06_genomic_phylogenomics.R")
        tmp <- {json.dumps(str(cls.tmp))}
        geno <- {json.dumps(str(cls.paths['geno']))}
        dbp <- {json.dumps(str(cls.paths['db']))}
        vcsv <- {json.dumps(str(cls.paths['vouchers']))}
        run_tag <- function(tag, tx, intra) {{
          od <- file.path(tmp, paste0("out_", tag)); td <- file.path(tmp, paste0("tbl_", tag)); fd <- file.path(tmp, paste0("fig_", tag))
          args <- c("--genomics-dir", geno, "--parsed-dir", file.path(tmp, paste0("parsed_", tag)), "--out-dir", od, "--tables-dir", td,
                    "--figure-a", file.path(fd, "fig06a_genomic_pcoa_tree_concordance.png"),
                    "--figure-b", file.path(fd, "fig06b_reticulation_gdi_landscape.png"),
                    "--figure-c", file.path(fd, "fig06c_cytonuclear_paralog_cytotype_audit.png"),
                    "--vouchers", vcsv, "--focal-taxon", tx, "--genomic-tier", "auto", "--run-id", tag, "--db", dbp, "--demo")
          if (intra) args <- c(args, "--intraspecific-mode")
          run_step06_genomic_phylogenomics(parse_step06_args(args))
        }}
        run_tag("t1_intra", "Packera dubia", TRUE)
        run_tag("t2_hyb", "Packera obovata", FALSE)
        run_tag("t3_sanger", "Packera crawfordii", FALSE)
        run_tag("t4_gap", "Packera hartiana", FALSE)
        source("scripts/analysis/genomic_triage_utils.R")
        source("scripts/analysis/07_triage_dashboard_synthesis.R")

        ga_t2 <- utils::read.csv(file.path(tmp, "out_t2_hyb", "genomic_assignments.csv"), stringsAsFactors = FALSE)
        ga_t3 <- utils::read.csv(file.path(tmp, "out_t3_sanger", "genomic_assignments.csv"), stringsAsFactors = FALSE)
        mk_v <- function(cat, sp, ga_dir, opts_g) {{
          n_a <- 8L
          df <- data.frame(
            catalogNumber = c(cat, sprintf("ANCH_%02d", seq_len(n_a))), plant_individual_id = 0L,
            species_standardized = c(sp, rep("Packera dubia", 4L), rep("Packera obovata", 4L)),
            species_raw = c(sp, rep("Packera dubia", 4L), rep("Packera obovata", 4L)),
            recordedBy = "Fuller", institutionCode = "NCU", determiner_tier = c("Tier_2_Silver", rep("Tier_1_Gold", n_a)),
            decimalLatitude = 35.5, decimalLongitude = -79.2, stateProvince = "NC", county = "Orange", eventDate = "2020-05-01",
            leaf_zone = "Basal_Rosette", basal_leaf_count = 5L, morph_tier_source = "Basal_Rosette", is_morph_outlier = FALSE,
            gmm_cluster = 1L, gmm_bayes_factor = 0.5, PC1 = c(0.1, rep(-1.2, 4L), rep(1.2, 4L)), PC2 = 0.0,
            cda_predicted_taxon = c(sp, rep("Packera dubia", 4L), rep("Packera obovata", 4L)), cda_posterior_prob = 0.90,
            cda_second_taxon = "Packera anonyma", cda_second_posterior_prob = 0.05,
            edaphic_best_fit_taxon = c(sp, rep("Packera dubia", 4L), rep("Packera obovata", 4L)), rf_predicted_prob = 0.88, type_edaphic_z_score = 0.3,
            hopkins_adjusted_doy = 120.0, delta_doy_resid = 0.0, is_active_anthesis = TRUE, open_set_df = 4L, min_congener_mahal_d2 = 16.0, stringsAsFactors = FALSE
          )
          df_h <- load_step07_genomic_metrics(list(genomic_dir = ga_dir, output_queue = file.path(ga_dir, "triage_queue.csv")), df)
          mask <- has_active_genomic_row_mask(df_h)
          attr(df_h, "step07_cytology_executed") <- TRUE
          dec <- apply_taxonomic_decision_matrix(df_h, opts = list())
          list(mask_first = mask[1L], call_first = dec$taxonomic_treatment_subtype[1L], pred_first = df_h$geno_predicted_taxon[1L])
        }}
        res_t2 <- mk_v("OBO_VOUCH_1", "Packera obovata", file.path(tmp, "out_t2_hyb"), list())
        # Force discordant Sanger voucher check on Tier 3
        df3 <- data.frame(catalogNumber = "CRAW_VOUCH_1", plant_individual_id = 0L, species_standardized = "Packera dubia",
                          species_raw = "Packera dubia", recordedBy = "Fuller", institutionCode = "NCU", determiner_tier = "Tier_2_Silver",
                          decimalLatitude = 35.7, decimalLongitude = -79.2, stateProvince = "NC", county = "Orange", eventDate = "2020-05-01",
                          leaf_zone = "Basal_Rosette", basal_leaf_count = 5L, morph_tier_source = "Basal_Rosette", is_morph_outlier = FALSE,
                          gmm_cluster = 1L, gmm_bayes_factor = 0.5, PC1 = 0.0, PC2 = 0.0,
                          cda_predicted_taxon = "Packera dubia", cda_posterior_prob = 0.90, cda_second_taxon = "Packera crawfordii", cda_second_posterior_prob = 0.05,
                          edaphic_best_fit_taxon = "Packera dubia", rf_predicted_prob = 0.88, type_edaphic_z_score = 0.3,
                          hopkins_adjusted_doy = 120.0, delta_doy_resid = 0.0, is_active_anthesis = TRUE, open_set_df = 4L, min_congener_mahal_d2 = 16.0, stringsAsFactors = FALSE)
        df3_h <- load_step07_genomic_metrics(list(genomic_dir = file.path(tmp, "out_t3_sanger"), output_queue = file.path(tmp, "out_t3_sanger", "triage_queue.csv")), df3)
        mask3 <- has_active_genomic_row_mask(df3_h)
        attr(df3_h, "step07_cytology_executed") <- TRUE
        dec3 <- apply_taxonomic_decision_matrix(df3_h, opts = list())

        writeLines(jsonlite::toJSON(list(
          t2_mask = res_t2$mask_first, t2_pred = res_t2$pred_first,
          t3_mask = mask3[1L], t3_sub = dec3$taxonomic_treatment_subtype[1L], t3_n_ga = nrow(ga_t3)
        ), auto_unbox = TRUE), {json.dumps(str(out_json))})
        """).strip() + "\n", encoding="utf-8")
        proc = subprocess.run(["Rscript", str(r_script)], cwd=PROJECT_ROOT, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise AssertionError(f"Batched Step 06 R test failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        cls.r_summary = json.loads(out_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp_dir.cleanup()

    def test_tier_autodetect_and_missing_file_stepdown(self) -> None:
        expected = {
            "t1_intra": "Tier_1_NextRAD_PopGen", "t2_hyb": "Tier_2_Public_HybSeq",
            "t3_sanger": "Tier_3_Public_MultiLocus", "t4_gap": "Tier_4_Cytology_Gap",
        }
        for tag, exp_tier in expected.items():
            mf = json.loads((self.tmp / f"parsed_{tag}" / "genomic_tier_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(mf["active_genomic_tier"], exp_tier)
        with tempfile.TemporaryDirectory() as td2:
            sub_g = Path(td2) / "geno"; shutil.copytree(self.paths["geno"], sub_g)
            shutil.rmtree(sub_g / "public_nuclear_backbone")
            r3 = ING.run_genomic_ingestion(sub_g, Path(td2) / "p3", db_path=self.paths["db"], min_hq_reads=0, vouchers_csv=self.paths["vouchers"], focal_taxon="Packera obovata", genomic_tier="auto")
            self.assertEqual(r3["active_genomic_tier"], "Tier_3_Public_MultiLocus")
            shutil.rmtree(sub_g / "public_organellar_and_loci")
            r4 = ING.run_genomic_ingestion(sub_g, Path(td2) / "p4", db_path=self.paths["db"], min_hq_reads=0, vouchers_csv=self.paths["vouchers"], focal_taxon="Packera obovata", genomic_tier="auto")
            self.assertEqual(r4["active_genomic_tier"], "Tier_4_Cytology_Gap")

    def test_stage2_complex_vs_stage3_intraspecific_scoping_and_global_cache(self) -> None:
        with tempfile.TemporaryDirectory() as td_c:
            with mock.patch.object(ING, "_parse_vcf_with_depths", side_effect=AssertionError("VCF re-parsed despite warm cache")):
                r_s2 = ING.run_genomic_ingestion(self.paths["geno"], Path(td_c) / "s2", db_path=self.paths["db"], min_hq_reads=0, vouchers_csv=self.paths["vouchers"], focal_taxon="Packera dubia", intraspecific_mode=False, genomic_tier="auto")
            self.assertEqual(r_s2["active_genomic_tier"], "Tier_1_NextRAD_PopGen")
            with open(self.tmp / "parsed_t1_intra" / "parsed_genomic_samples.csv", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            roles = {r["role_in_run"] for r in rows}
            self.assertIn("focal", roles); self.assertIn("polarizing_sister", roles)
            self.assertGreaterEqual(sum(1 for r in rows if r["role_in_run"] == "polarizing_sister"), 2)
        if not SKIP_R_COND:
            with open(self.tmp / "out_t1_intra" / "genomic_assignments.csv", encoding="utf-8") as fh:
                ga1 = list(csv.DictReader(fh))
            pol_rows = [r for r in ga1 if r.get("role_in_run") == "polarizing_sister"]
            self.assertTrue(all(r["genetic_cluster"] == "Polarizing_Sister" for r in pol_rows))
            with open(self.tmp / "out_t1_intra" / "coalescent_gdi_delimitation.csv", encoding="utf-8") as fh:
                gdi1 = list(csv.DictReader(fh))
            self.assertTrue(any(r["gdi_source"] == "nextrad_msc" for r in gdi1))

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_hybseq_g1_promotion_vs_sanger_stage1_non_mutation(self) -> None:
        with open(self.tmp / "out_t2_hyb" / "genomic_assignments.csv", encoding="utf-8") as fh:
            ga2 = list(csv.DictReader(fh))
        obo_g1 = [r for r in ga2 if r["catalogNumber"] == "OBO_VOUCH_1" and r["marker_platform"] == "HybSeq_TargetCapture"]
        self.assertTrue(obo_g1)
        self.assertEqual(obo_g1[0]["linkage_tier"], "G1_Direct_Voucher")
        self.assertEqual(obo_g1[0]["geno_predicted_taxon"], "Packera obovata")
        self.assertTrue(bool(self.r_summary["t2_mask"]))
        self.assertFalse(bool(self.r_summary["t3_mask"]))
        self.assertNotEqual(self.r_summary["t3_sub"], "Misidentification_Reassign")

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_canonical_png_plates_and_monograph_metric_hydration(self) -> None:
        self.assertEqual(list(self.tmp.rglob("*.pdf")), [])
        for tag in ("t1_intra", "t2_hyb", "t3_sanger", "t4_gap"):
            for fn in ("fig06a_genomic_pcoa_tree_concordance.png", "fig06b_reticulation_gdi_landscape.png", "fig06c_cytonuclear_paralog_cytotype_audit.png"):
                fp = self.tmp / f"fig_{tag}" / fn
                self.assertTrue(fp.is_file(), f"Missing {fp}")
                self.assertGreater(fp.stat().st_size, 20_000, f"Too small: {fp}")
        ws_m = self.tmp / "monograph_ws"; runs_m = ws_m / "runs"
        for tag, rname, sp in (("t2_hyb", "r_obo", "Packera obovata"), ("t3_sanger", "r_craw", "Packera crawfordii")):
            dst = runs_m / rname / "outputs" / "genomics"; shutil.copytree(self.tmp / f"out_{tag}", dst)
            (runs_m / rname / "data" / "tables").mkdir(parents=True, exist_ok=True)
            with open(runs_m / rname / "data" / "tables" / "triage_queue.csv", "w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh); w.writerow(["catalogNumber", "plant_individual_id", "species_standardized", "taxonomic_treatment_subtype"])
                w.writerow(["V1", "0", sp, "Confirmed_Species"])
        m2, _ = _extract_stage_metrics(stage_name="Stage_2", run_id="r_obo", focal_std="Packera obovata", runs_root=runs_m, conn=None)
        self.assertEqual(m2["active_genomic_tier"], "Tier_2_Public_HybSeq")
        self.assertIsNotNone(m2["median_coalescent_gdi"])
        self.assertIn("Packera obovata", m2["public_backbone_clade"])
        m3, _ = _extract_stage_metrics(stage_name="Stage_2", run_id="r_craw", focal_std="Packera crawfordii", runs_root=runs_m, conn=None)
        self.assertEqual(m3["active_genomic_tier"], "Tier_3_Public_MultiLocus")
        self.assertIsNone(m3["median_coalescent_gdi"])


if __name__ == "__main__":
    unittest.main()
