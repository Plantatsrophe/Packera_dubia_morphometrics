#!/usr/bin/env python3
"""Warm-R & Python QA suite for quality-gated organellar capture, maternal hybrid polarity, and OOP rationale flags."""

from __future__ import annotations

import csv, gzip, io, json, os, random, shutil, subprocess, sys, tempfile, textwrap, unittest, zipfile
from pathlib import Path
from typing import Any, Dict

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

SKIP_R_COND = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
    or shutil.which("Rscript") is None
)


def _fq_gz(reads: list[tuple[str, str]]) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for rid, seq in reads:
            gz.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n".encode("ascii"))
    return buf.getvalue()


def _run_warm_r_organellar_suite() -> Dict[str, Any]:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        out_json, r_script = tmp / "organellar_triage_results.json", tmp / "run_organellar_checks.R"
        r_script.write_text(
            textwrap.dedent(
                f"""
                if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
                suppressPackageStartupMessages({{
                  source("scripts/analysis/genomic_triage_utils.R")
                  source("scripts/analysis/step07_organellar_triage_utils.R")
                  source("scripts/analysis/07_triage_dashboard_synthesis.R")
                }})
                mk_base <- function() {{
                  n_dub <- 6L; n_ano <- 6L; n_aur <- 4L; nt <- n_dub + n_ano + n_aur
                  sp_v <- c(rep("Packera dubia", n_dub), rep("Packera anonyma", n_ano), rep("Packera aurea", n_aur))
                  data.frame(
                    catalogNumber = sprintf("ANCH_%02d", seq_len(nt)), plant_individual_id = 0L,
                    species_standardized = sp_v, species_raw = sp_v,
                    recordedBy = rep(c("Fuller", "Weakley", "Sorrie", "Radford"), length.out = nt),
                    institutionCode = "NCU", determiner_tier = "Tier_1_Gold",
                    decimalLatitude = c(seq(35.10, 35.15, length.out = n_dub), seq(35.11, 35.16, length.out = n_ano), seq(35.12, 35.15, length.out = n_aur)),
                    decimalLongitude = rep(-80.20, nt), stateProvince = "NC", county = "Montgomery", eventDate = "2023-05-01",
                    leaf_zone = "Basal_Rosette", basal_leaf_count = 5L, morph_tier_source = "Basal_Rosette", is_morph_outlier = FALSE,
                    gmm_cluster = 1L, gmm_bayes_factor = 0.5,
                    PC1 = c(seq(-1.6, -1.4, length.out = n_dub), seq(1.4, 1.6, length.out = n_ano), seq(-0.1, 0.1, length.out = n_aur)),
                    PC2 = c(seq(-0.2, 0.2, length.out = n_dub), seq(-0.2, 0.2, length.out = n_ano), seq(1.9, 2.1, length.out = n_aur)),
                    cda_predicted_taxon = sp_v, cda_posterior_prob = 0.94,
                    cda_second_taxon = c(rep("Packera anonyma", n_dub), rep("Packera dubia", n_ano), rep("Packera dubia", n_aur)),
                    cda_second_posterior_prob = 0.04, edaphic_best_fit_taxon = sp_v, rf_predicted_prob = 0.92, type_edaphic_z_score = 0.3,
                    hopkins_adjusted_doy = 120.0, delta_doy_resid = 0.0, is_active_anthesis = TRUE, open_set_df = 4L, min_congener_mahal_d2 = 16.0,
                    genetic_sample_id = sprintf("GEN_%02d", seq_len(nt)), geno_linkage_tier = "G1_Direct_Voucher", geno_quality_pass = TRUE,
                    geno_predicted_taxon = sp_v, geno_posterior_prob = 0.96, geno_margin = 0.90, is_geno_admixed = FALSE,
                    is_putative_novel_or_hybrid_query = FALSE, geno_parent_1 = NA_character_, geno_parent_2 = NA_character_,
                    interspecific_ancestry_q12 = 0.02, hybrid_generation_class = "Pure_Parent_1", f_branch_max = 0.0, cytonuclear_discordant = 0L,
                    public_plastome_match = sp_v, public_its_match = sp_v, waxy_lineage_match = sp_v, stringsAsFactors = FALSE
                  )
                }}

                b_df <- mk_base()[mk_base()$species_standardized %in% c("Packera dubia", "Packera anonyma"), , drop = FALSE]
                s1_row <- b_df[1L, ]; s1_row$catalogNumber <- "STG1_MISID"; s1_row$genetic_sample_id <- "G_STG1"
                s1_row$determiner_tier <- "Tier_2_Silver"; s1_row$PC1 <- 1.55; s1_row$PC2 <- 0.0
                s1_row$cda_predicted_taxon <- "Packera anonyma"; s1_row$edaphic_best_fit_taxon <- "Packera anonyma"
                s1_row$rf_predicted_prob <- 0.91; s1_row$type_edaphic_z_score <- 3.2; s1_row$geno_predicted_taxon <- "Packera anonyma"
                oop_row <- b_df[2L, ]; oop_row$catalogNumber <- "OOP_QRY"; oop_row$genetic_sample_id <- "G_OOP"
                oop_row$is_putative_novel_or_hybrid_query <- TRUE; oop_row$geno_predicted_taxon <- "Packera dubia"
                oop_row$geno_posterior_prob <- 0.88; oop_row$geno_margin <- 0.65

                clean_df <- rbind(b_df, s1_row, oop_row)
                clean_df$public_plastome_match <- NA_character_; clean_df$public_its_match <- NA_character_; clean_df$waxy_lineage_match <- NA_character_
                attr(clean_df, "step07_cytology_executed") <- TRUE
                conf_df <- clean_df
                conf_df$public_plastome_match[conf_df$catalogNumber == "STG1_MISID"] <- "Packera dubia"
                conf_df$public_its_match[conf_df$catalogNumber %in% c("STG1_MISID", "OOP_QRY")] <- "Packera plattensis"
                conf_df$waxy_lineage_match[conf_df$catalogNumber %in% c("STG1_MISID", "OOP_QRY")] <- "Packera plattensis"
                pop_s1 <- data.frame(
                  genetic_sample_id = c("G_STG1", "G_STG1", "G_STG1", "G_OOP", "G_OOP"),
                  locus = c("plastid", "its", "waxy", "its", "waxy"),
                  best_match_taxon = c("Packera dubia", "Packera plattensis", "Packera plattensis", "Packera plattensis", "Packera plattensis"),
                  n_shared_11mers = c(180L, 180L, 180L, 135L, 115L), kmer_containment = c(0.96, 0.96, 0.96, 0.91, 0.89), stringsAsFactors = FALSE
                )
                res_s1_clean <- apply_taxonomic_decision_matrix(clean_df, opts = list())
                res_s1_conf <- apply_taxonomic_decision_matrix(conf_df, opts = list(organellar_placement_df = pop_s1))

                m_df <- mk_base()
                cap_hq <- m_df[1L, ]; cap_hq$catalogNumber <- "CAP_HQ"; cap_hq$genetic_sample_id <- "G_CAP_HQ"
                cap_hq$PC1 <- -0.05; cap_hq$PC2 <- 0.95; cap_hq$cda_predicted_taxon <- "Packera aurea"
                cap_hq$cda_posterior_prob <- 0.52; cap_hq$cda_second_taxon <- "Packera dubia"; cap_hq$cda_second_posterior_prob <- 0.48
                cap_hq$geno_predicted_taxon <- "Packera dubia"; cap_hq$geno_posterior_prob <- 0.94
                cap_hq$interspecific_ancestry_q12 <- 0.04; cap_hq$public_plastome_match <- "Packera aurea"
                cap_low <- cap_hq; cap_low$catalogNumber <- "CAP_LOW"; cap_low$genetic_sample_id <- "G_CAP_LOW"

                hyb_hq <- m_df[1L, ]; hyb_hq$catalogNumber <- "HYB_HQ"; hyb_hq$genetic_sample_id <- "G_HYB_HQ"
                hyb_hq$species_raw <- "Packera dubia x Packera anonyma"; hyb_hq$PC1 <- 0.0; hyb_hq$PC2 <- 0.0
                hyb_hq$cda_predicted_taxon <- "Packera dubia"; hyb_hq$cda_posterior_prob <- 0.51
                hyb_hq$cda_second_taxon <- "Packera anonyma"; hyb_hq$cda_second_posterior_prob <- 0.49
                hyb_hq$geno_posterior_prob <- 0.52; hyb_hq$is_geno_admixed <- TRUE; hyb_hq$interspecific_ancestry_q12 <- 0.82
                hyb_hq$hybrid_generation_class <- "F1_Hybrid"; hyb_hq$geno_parent_1 <- "Packera dubia"; hyb_hq$geno_parent_2 <- "Packera anonyma"
                hyb_hq$public_plastome_match <- "Packera anonyma"; hyb_hq$public_its_match <- "Packera dubia"
                hyb_low <- hyb_hq; hyb_low$catalogNumber <- "HYB_LOWCP"; hyb_low$genetic_sample_id <- "G_HYB_LOWCP"
                hyb_dnr <- hyb_hq; hyb_dnr$catalogNumber <- "HYB_DONOR"; hyb_dnr$genetic_sample_id <- "G_HYB_DONOR"
                hyb_dnr$species_raw <- "Packera dubia"; hyb_dnr$geno_parent_2 <- NA_character_
                hyb_dnr$cda_predicted_taxon <- "Packera dubia"; hyb_dnr$cda_posterior_prob <- 0.52
                hyb_dnr$cda_second_taxon <- "Packera aurea"; hyb_dnr$cda_second_posterior_prob <- 0.48
                hyb_dnr$public_plastome_match <- "Packera dubia"; hyb_dnr$public_its_match <- "Packera anonyma"

                full_m <- rbind(m_df, cap_hq, cap_low, hyb_hq, hyb_low, hyb_dnr)
                attr(full_m, "step07_cytology_executed") <- TRUE
                pop_m <- data.frame(
                  genetic_sample_id = c("G_CAP_HQ", "G_CAP_LOW", "G_HYB_HQ", "G_HYB_HQ", "G_HYB_LOWCP", "G_HYB_DONOR", "G_HYB_DONOR"),
                  locus = c("plastid", "plastid", "plastid", "its", "plastid", "plastid", "its"),
                  best_match_taxon = c("Packera aurea", "Packera aurea", "Packera anonyma", "Packera dubia", "Packera anonyma", "Packera dubia", "Packera anonyma"),
                  n_shared_11mers = c(150L, 40L, 140L, 120L, 140L, 130L, 130L),
                  kmer_containment = c(0.92, 0.55, 0.91, 0.89, 0.30, 0.90, 0.90), stringsAsFactors = FALSE
                )
                res_m <- apply_taxonomic_decision_matrix(full_m, opts = list(organellar_placement_df = pop_m))
                pick_row <- function(d, cid) as.list(d[d$catalogNumber == cid, c("taxonomic_status_call", "recommended_determination", "taxonomic_treatment_subtype", "hybrid_parent_1", "hybrid_parent_2", "synthesis_rationale"), drop = FALSE][1L, ])
                out_obj <- list(
                  s1_clean_calls = as.character(res_s1_clean$taxonomic_status_call), s1_conf_calls = as.character(res_s1_conf$taxonomic_status_call),
                  s1_clean_rec = as.character(res_s1_clean$recommended_determination), s1_conf_rec = as.character(res_s1_conf$recommended_determination),
                  s1_clean_sub = as.character(res_s1_clean$taxonomic_treatment_subtype), s1_conf_sub = as.character(res_s1_conf$taxonomic_treatment_subtype),
                  stg1_row = pick_row(res_s1_conf, "STG1_MISID"), oop_row = pick_row(res_s1_conf, "OOP_QRY"),
                  cap_hq = pick_row(res_m, "CAP_HQ"), cap_low = pick_row(res_m, "CAP_LOW"),
                  hyb_hq = pick_row(res_m, "HYB_HQ"), hyb_low = pick_row(res_m, "HYB_LOWCP"), hyb_dnr = pick_row(res_m, "HYB_DONOR")
                )
                writeLines(jsonlite::toJSON(out_obj, auto_unbox = TRUE), "{out_json}")
                """
            ).strip() + "\n",
            encoding="utf-8",
        )
        proc = subprocess.run(["Rscript", str(r_script)], cwd=PROJECT_ROOT, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise AssertionError(f"Warm R suite failed ({proc.returncode}):\n{proc.stdout}\n{proc.stderr}")
        return json.loads(out_json.read_text(encoding="utf-8"))


class TestStep07OrganellarTriage(unittest.TestCase):
    """QA contracts for reference cytonuclear capture, hybrid polarity, and non-mutating OOP rationale flags."""

    r_data: Dict[str, Any] = {}

    @classmethod
    def setUpClass(cls) -> None:
        if not SKIP_R_COND:
            cls.r_data = _run_warm_r_organellar_suite()

    def test_python_small_batch_reference_cytonuclear_capture(self) -> None:
        from scripts.data_prep import organellar_cytonuclear_utils as ocu
        from scripts.data_prep import skim_offtarget_organellar as skim_mod

        rng = random.Random(202611)
        _mk = lambda n=150: "".join(rng.choice("ACGT") for _ in range(n))
        cp_aur, its_ser, wx_ser = _mk(), _mk(), _mk()
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            geno_dir, raw_dir = tmp / "genomics", tmp / "raw"
            org_dir = geno_dir / "public_organellar_and_loci"
            org_dir.mkdir(parents=True, exist_ok=True); raw_dir.mkdir(parents=True, exist_ok=True)
            (org_dir / "public_packera_18complete_plastomes.fasta").write_text(f">MK_CP_AUR Packera aurea\n{cp_aur}\n", encoding="utf-8")
            (org_dir / "public_packera_its_ets.fasta").write_text(f">MT_ITS_SER Packera serpenticola\n{its_ser}\n", encoding="utf-8")
            (org_dir / "public_packera_waxy_nuclear.fasta").write_text(f">JX_WX_SER Packera serpenticola\n{wx_ser}\n", encoding="utf-8")
            zpath = raw_dir / "ser_batch.zip"
            with zipfile.ZipFile(zpath, "w") as zf:
                for sid in ("PACSERNCS1M", "PACSERNCS2M"):
                    reads = [(f"{sid}_c{j}", cp_aur) for j in range(5)] + [(f"{sid}_i{j}", its_ser) for j in range(5)] + [(f"{sid}_w{j}", wx_ser) for j in range(5)]
                    zf.writestr(f"fq/{sid}_R1.fastq.gz", _fq_gz(reads))
            res = skim_mod.run_offtarget_organellar_skim(fastq_source=zpath, genomics_dir=geno_dir)
            self.assertEqual(res["status"], "skimmed")
            with open(org_dir / "cytonuclear_discordance.csv", "r", encoding="utf-8", newline="") as fh:
                by_id = {r["genetic_sample_id"]: r for r in csv.DictReader(fh)}
            for sid in ("PACSERNCS1M", "PACSERNCS2M"):
                self.assertEqual(by_id[sid]["cytonuclear_discordant"], "1")
                self.assertEqual(by_id[sid]["public_plastome_match"], "Packera aurea")
            sp_best, cont, n_sh = ocu.match_best_public_lineage_with_metrics(cp_aur, {"MK Packera aurea": cp_aur})
            self.assertEqual(sp_best, "Packera aurea"); self.assertGreaterEqual(cont, 0.94); self.assertGreaterEqual(n_sh, 110)

        neg = ocu._detect_cytonuclear_discordance(
            {"PACSERNCS1M": "C1", "PACSERNCS2M": "C1", "PACSERNCS3M": "C1", "EREHIENCS1M": "C1"},
            {"PACSERNCS1M": "I1", "PACSERNCS2M": "I1", "PACSERNCS3M": "I1", "EREHIENCS1M": "I2"},
            public_plastome_matches={k: "Packera aurea" for k in ("PACSERNCS1M", "PACSERNCS2M", "PACSERNCS3M", "EREHIENCS1M")},
            public_its_matches={"PACSERNCS1M": "Packera serpenticola", "PACSERNCS2M": "Packera serpenticola", "PACSERNCS3M": "Packera aurea", "EREHIENCS1M": "Erechtites"},
            waxy_lineage_matches={"PACSERNCS1M": "Packera serpenticola", "PACSERNCS2M": "Packera serpenticola", "PACSERNCS3M": "Packera serpenticola", "EREHIENCS1M": "Erechtites"},
            plastome_match_metrics={"PACSERNCS1M": (0.45, 120), "PACSERNCS2M": (0.94, 30), "PACSERNCS3M": (0.94, 120), "EREHIENCS1M": (0.94, 120)},
        )
        self.assertEqual({r["genetic_sample_id"]: r["cytonuclear_discordant"] for r in neg}, {"PACSERNCS1M": 0, "PACSERNCS2M": 0, "PACSERNCS3M": 0, "EREHIENCS1M": 0})

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_stage1_strictly_invariant_under_conflicting_organellar_and_oop_matches(self) -> None:
        self.assertEqual(self.r_data["stg1_row"]["taxonomic_treatment_subtype"], "Misidentification_Reassign")
        self.assertEqual(self.r_data["stg1_row"]["recommended_determination"], "Packera anonyma")
        self.assertEqual(self.r_data["s1_clean_calls"], self.r_data["s1_conf_calls"])
        self.assertEqual(self.r_data["s1_clean_rec"], self.r_data["s1_conf_rec"])
        self.assertEqual(self.r_data["s1_clean_sub"], self.r_data["s1_conf_sub"])

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_reference_capture_blocks_false_stage2_hybrid_and_annotates_donor(self) -> None:
        cap_hq, cap_low = self.r_data["cap_hq"], self.r_data["cap_low"]
        self.assertNotEqual(cap_hq["taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertEqual(cap_hq["recommended_determination"], "Packera dubia")
        self.assertIn("Ancient Introgression / Organellar Capture", cap_hq["synthesis_rationale"])
        self.assertIn("cpDNA=Packera aurea", cap_hq["synthesis_rationale"])
        self.assertNotIn("Ancient Introgression / Organellar Capture", cap_low["synthesis_rationale"])
        self.assertNotIn("cpDNA=Packera aurea", cap_low["synthesis_rationale"])
        self.assertEqual(self.r_data["hyb_dnr"]["taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertEqual(self.r_data["hyb_dnr"]["hybrid_parent_2"], "Packera anonyma")

    @unittest.skipIf(SKIP_R_COND, "Fast audit or Rscript unavailable")
    def test_confirmed_stage2_hybrid_maternal_polarity_and_nonmutating_oop_flag(self) -> None:
        hyb_hq, hyb_low, oop = self.r_data["hyb_hq"], self.r_data["hyb_low"], self.r_data["oop_row"]
        self.assertEqual(hyb_hq["taxonomic_treatment_subtype"], "Hybrid_Introgressant")
        self.assertIn("[Maternal cpDNA Parent: Packera anonyma (cont=0.91, 11mers=140); Putative Pollen Donor: Packera dubia, ITS/waxy corroborated]", hyb_hq["synthesis_rationale"])
        self.assertNotIn("Maternal cpDNA Parent:", hyb_low["synthesis_rationale"])
        self.assertIn("[Out-of-Panel Nuclear Locus Match: ITS+waxy=Packera plattensis (ITS_cont=0.91, waxy_cont=0.89); taxon absent from NextRAD reference panel]", oop["synthesis_rationale"])


if __name__ == "__main__":
    unittest.main()
