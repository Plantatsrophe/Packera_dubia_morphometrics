"""Targeted QA tests for v5.10.0 outlier quarantine gates (Steps 01, 05, 06)."""

from __future__ import annotations

import inspect
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from typing import Any, Callable


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _find_rscript() -> str | None:
    """Locate an Rscript executable on PATH for warm-R statistical checks."""
    return shutil.which("Rscript")


FAST_AUDIT = (
    os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
    or os.environ.get("PACKERA_FAST_AUDIT") == "1"
)
RSCRIPT_BIN = _find_rscript()


def _resolve_coord_classifier(module: Any) -> Callable[..., tuple[bool, str]]:
    """Locate the primary Step 01 coordinate quarantine classifier function."""
    for name in (
        "classify_centroid_or_institution_coord", "classify_coord_quarantine",
        "flag_centroid_or_institution_coord", "screen_centroid_or_institution_coord",
        "detect_centroid_or_institution_coord", "check_centroid_or_institution_coord",
        "evaluate_coord_hygiene", "assess_coord_hygiene",
    ):
        fn = getattr(module, name, None)
        if callable(fn):
            return fn
    expected_kw = {"locality", "uncertainty_m", "county", "state_prov", "issues", "is_type_specimen", "us_spatial_graph"}
    candidates = [
        (len(set(inspect.signature(fn).parameters.keys()) & expected_kw), fn)
        for name, fn in inspect.getmembers(module, inspect.isfunction)
        if not name.startswith("_") and fn.__module__ == module.__name__
    ]
    candidates = [c for c in candidates if c[0] >= 3]
    if candidates:
        candidates.sort(key=lambda item: item[0], reverse=True)
        return candidates[0][1]
    raise AssertionError("Unable to locate coordinate quarantine classifier in coord_hygiene_utils.")


def _call_coord_classifier(fn: Callable[..., tuple[bool, str]], lat: float, lon: float, **kwargs: Any) -> tuple[bool, str]:
    """Invoke the Step 01 classifier with signature-filtered keyword arguments."""
    params = inspect.signature(fn).parameters
    filtered = {k: v for k, v in kwargs.items() if k in params}
    if "state_prov" in kwargs and "state_prov" not in params:
        for alt in ("stateProvince", "state_province", "state"):
            if alt in params:
                filtered[alt] = kwargs["state_prov"]
                break
    result = fn(lat, lon, **filtered)
    if not (isinstance(result, tuple) and len(result) == 2):
        raise AssertionError(f"Expected 2-tuple from coordinate classifier, got {result!r}")
    return bool(result[0]), str(result[1])


_R_QUARANTINE_SCRIPT = textwrap.dedent(
    r"""
    if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
    repo_root <- normalizePath("__REPO_ROOT__", mustWork = TRUE)
    tmp_dir <- normalizePath("__TMP_DIR__", mustWork = TRUE)
    setwd(repo_root)

    source(file.path(repo_root, "scripts", "core", "db_utils.R"))
    source(file.path(repo_root, "scripts", "analysis", "edaphic_fine_grain_utils.R"))
    source(file.path(repo_root, "scripts", "analysis", "step05_env_outlier_utils.R"))
    source(file.path(repo_root, "scripts", "analysis", "05_multimodal_spatial_rf.R"))
    source(file.path(repo_root, "scripts", "analysis", "step06_rogue_sequence_qc.R"))
    source(file.path(repo_root, "scripts", "analysis", "step06_public_backbone.R"))

    on.exit(setwd(repo_root), add = TRUE)
    setwd(tmp_dir)

    # --- Step 05: flag_coordinate_uncertainty ---
    df_no_col <- data.frame(
      catalogNumber = c("V1", "V2", "V3"),
      decimalLatitude = c(35.9, 35.9, 35.9),
      decimalLongitude = c(-79.0, -79.1, -79.2),
      coordinateUncertaintyInMeters = c(100, 100, 5000),
      is_type_specimen = c(FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
    flagged_base <- flag_coordinate_uncertainty(df_no_col)
    stopifnot(identical(flagged_base$is_edaphic_eligible_coord, c(TRUE, TRUE, FALSE)))
    stopifnot(flagged_base$edaphic_coord_weight[1] > 0.0, flagged_base$edaphic_coord_weight[3] == 0.0)

    df_with_col <- df_no_col
    df_with_col$is_centroid_or_inst_coord <- c(FALSE, TRUE, FALSE)
    flagged_quar <- flag_coordinate_uncertainty(df_with_col)
    stopifnot(identical(flagged_quar$is_edaphic_eligible_coord, c(TRUE, FALSE, FALSE)))
    stopifnot(flagged_quar$edaphic_coord_weight[1] > 0.0, flagged_quar$edaphic_coord_weight[2] == 0.0)

    # --- Step 05: Multi-taxon MCD screening + azonal rescue + type exemption ---
    set.seed(20261115L)
    n_dubia <- 30L
    n_anon <- 15L
    n_norm <- n_dubia + n_anon
    n_tot <- n_norm + 3L
    env_vars <- c("bio01", "bio04", "bio12", "soil_ph_0_15", "soil_sand_0_15")

    df_srf <- data.frame(
      catalogNumber = c(paste0("DUB_", seq_len(n_dubia)), paste0("ANO_", seq_len(n_anon)),
                        "DUB_EXT_ROADSIDE", "DUB_EXT_SERP", "DUB_EXT_TYPE"),
      plant_individual_id = rep(1L, n_tot),
      species_standardized = c(rep("Packera dubia", n_dubia), rep("Packera anonyma", n_anon), rep("Packera dubia", 3L)),
      latitude = c(runif(n_norm, 34.8, 36.4), 35.3, 35.1, 35.7),
      longitude = c(runif(n_norm, -83.4, -79.2), -81.5, -83.6, -80.2),
      decimalLatitude = c(runif(n_norm, 34.8, 36.4), 35.3, 35.1, 35.7),
      decimalLongitude = c(runif(n_norm, -83.4, -79.2), -81.5, -83.6, -80.2),
      determiner_tier = c(rep("Tier_1_Gold", n_norm), "Tier_1_Gold", "Tier_1_Gold", "Tier_1_Gold"),
      authority_tier = c(rep("Tier_1_Gold", n_norm), "Tier_1_Gold", "Tier_1_Gold", "Tier_0_Type"),
      authority_weight = c(rep(2.0, n_norm), 2.0, 2.0, 3.0),
      is_passive = rep(FALSE, n_tot),
      is_type_specimen = c(rep(FALSE, n_norm), FALSE, FALSE, TRUE),
      type_status = c(rep("", n_norm), "", "", "Holotype"),
      is_centroid_or_inst_coord = rep(FALSE, n_tot),
      coord_quarantine_reason = rep("", n_tot),
      coordinateUncertaintyInMeters = rep(100.0, n_tot),
      is_edaphic_eligible_coord = rep(TRUE, n_tot),
      edaphic_coord_weight = rep(1.0, n_tot),
      is_azonal_lithology = c(rep(FALSE, n_norm), FALSE, TRUE, FALSE),
      sgmc_lithology_unit = c(rep("granite gneiss", n_norm), "granite gneiss", "ultramafic serpentine", "granite gneiss"),
      locality = c(rep("Piedmont upland forest", n_norm), "roadside", "Buck Creek serpentine barren outcrop", "Type locality"),
      habitat = c(rep("upland woodland", n_norm), "roadside", "Buck Creek serpentine barren outcrop", "granite woodland"),
      georef_qc_status = rep("valid", n_tot),
      bio01 = c(rnorm(n_dubia, 14.0, 1.0), rnorm(n_anon, 13.5, 1.0), rep(22.0, 3L)),
      bio04 = c(rnorm(n_dubia, 650.0, 40.0), rnorm(n_anon, 640.0, 40.0), rep(970.0, 3L)),
      bio12 = c(rnorm(n_dubia, 1250.0, 80.0), rnorm(n_anon, 1280.0, 80.0), rep(1890.0, 3L)),
      soil_ph_0_15 = c(rnorm(n_dubia, 5.6, 0.3), rnorm(n_anon, 5.5, 0.3), rep(8.0, 3L)),
      soil_sand_0_15 = c(rnorm(n_dubia, 45.0, 5.0), rnorm(n_anon, 46.0, 5.0), rep(85.0, 3L)),
      stringsAsFactors = FALSE
    )
    srf_res <- compute_spatial_rf_mems(df_srf, env_vars, n_trees = 50L)

    idx_road <- which(df_srf$catalogNumber == "DUB_EXT_ROADSIDE")
    idx_serp <- which(df_srf$catalogNumber == "DUB_EXT_SERP")
    idx_type <- which(df_srf$catalogNumber == "DUB_EXT_TYPE")
    stopifnot(idx_road %in% srf_res$passive_idx, !(idx_road %in% srf_res$anchor_idx))
    stopifnot(idx_serp %in% srf_res$anchor_idx, !(idx_serp %in% srf_res$passive_idx))
    stopifnot(idx_type %in% srf_res$anchor_idx, !(idx_type %in% srf_res$passive_idx))

    flag_df <- NULL
    for (nm in c("predictions", "env_outlier_audit", "env_outlier_df", "audit_df", "df")) {
      cand <- srf_res[[nm]]
      if (is.data.frame(cand) && "env_mcd_outlier_flag" %in% names(cand) && "catalogNumber" %in% names(cand)) {
        flag_df <- cand
        break
      }
    }
    if (is.null(flag_df) && !is.null(srf_res$env_mcd_outlier_flag) && length(srf_res$env_mcd_outlier_flag) == nrow(df_srf)) {
      flag_df <- df_srf
      flag_df$env_mcd_outlier_flag <- srf_res$env_mcd_outlier_flag
      if (!is.null(srf_res$env_outlier_rescued_azonal)) flag_df$env_outlier_rescued_azonal <- srf_res$env_outlier_rescued_azonal
      if (!is.null(srf_res$env_mahal_d2)) flag_df$env_mahal_d2 <- srf_res$env_mahal_d2
    }
    stopifnot(!is.null(flag_df))
    row_r <- flag_df[flag_df$catalogNumber == "DUB_EXT_ROADSIDE", , drop = FALSE]
    row_s <- flag_df[flag_df$catalogNumber == "DUB_EXT_SERP", , drop = FALSE]
    row_t <- flag_df[flag_df$catalogNumber == "DUB_EXT_TYPE", , drop = FALSE]
    stopifnot(isTRUE(as.logical(row_r$env_mcd_outlier_flag[1])))
    stopifnot(isTRUE(as.logical(row_s$env_outlier_rescued_azonal[1])), !isTRUE(as.logical(row_s$env_mcd_outlier_flag[1])))
    stopifnot(!isTRUE(as.logical(row_t$env_mcd_outlier_flag[1])))

    env_out_csv <- file.path(tmp_dir, "environmental_outlier_quarantine_audit.csv")
    write_env_outlier_quarantine_audit(
      df_srf, srf_res, env_out_csv,
      final_anchor_mask = seq_len(nrow(df_srf)) %in% srf_res$anchor_idx
    )
    stopifnot(file.exists(env_out_csv))
    env_audit <- utils::read.csv(env_out_csv, stringsAsFactors = FALSE)
    stopifnot(ncol(env_audit) == 10L, nrow(env_audit) >= 2L)
    check_env_row <- function(cat_id, want_outlier, want_rescued) {
      sub_df <- env_audit[env_audit$catalogNumber == cat_id, , drop = FALSE]
      stopifnot(nrow(sub_df) >= 1L)
      if ("env_mcd_outlier_flag" %in% names(sub_df)) stopifnot(identical(isTRUE(as.logical(sub_df$env_mcd_outlier_flag[1])), want_outlier))
      if ("env_outlier_rescued_azonal" %in% names(sub_df)) stopifnot(identical(isTRUE(as.logical(sub_df$env_outlier_rescued_azonal[1])), want_rescued))
    }
    check_env_row("DUB_EXT_ROADSIDE", TRUE, FALSE)
    check_env_row("DUB_EXT_SERP", FALSE, TRUE)
    if ("DUB_EXT_TYPE" %in% env_audit$catalogNumber) check_env_row("DUB_EXT_TYPE", FALSE, FALSE)

    # --- Step 06: Newick rogue tip screening ---
    nwk <- paste0(
      "((((((Packera_aurea:0.008,Packera_dubia:0.010):0.006,",
      "(Packera_paupercula:0.007,Packera_pseudaurea:0.012):0.005):0.004,",
      "(Packera_obovata:0.009,Packera_schweinitziana:0.011):0.006):0.005,",
      "(Packera_plattensis:0.006,Packera_tomentosa:0.014):0.007):0.005,",
      "Packera_anonyma:0.450):0.050,Jacobaea_vulgaris:0.180);"
    )
    tr <- ape::read.tree(text = nwk)
    screen_res <- screen_public_tree_rogue_tips(tr)
    tr_clean <- if (inherits(screen_res, "phylo")) screen_res else screen_res$tree
    stopifnot(inherits(tr_clean, "phylo"), length(tr_clean$tip.label) == 9L)
    clean_norm <- gsub("_", " ", tr_clean$tip.label)
    stopifnot(!("Packera anonyma" %in% clean_norm), "Jacobaea vulgaris" %in% clean_norm)

    # --- Step 06: FASTA k-mer rogue screening & query placement ---
    geno_dir <- file.path(tmp_dir, "genomics")
    out_dir <- file.path(tmp_dir, "step06_out")
    dir.create(geno_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    dna_bases <- c("A", "C", "G", "T")
    rot <- c("A" = "C", "C" = "G", "G" = "T", "T" = "A")
    mutate_k <- function(seq_vec, idxs) { seq_vec[idxs] <- unname(rot[seq_vec[idxs]]); seq_vec }
    set.seed(20261116L)
    aurea_base <- sample(dna_bases, 400L, replace = TRUE)
    glabella_base <- mutate_k(aurea_base, seq(141L, 400L, by = 2L))
    set.seed(99991L)
    fungal_vec <- sample(dna_bases, 400L, replace = TRUE)

    fasta_lines <- c(
      sprintf(">OP00000%d Packera aurea\n%s", 1:4, vapply(c(15L, 45L, 80L, 110L), function(p) paste0(mutate_k(aurea_base, c(p, p + 200L)), collapse = ""), character(1L))),
      sprintf(">OP00001%d Packera glabella\n%s", 1:3, vapply(c(20L, 55L, 95L), function(p) paste0(mutate_k(glabella_base, c(p, p + 200L)), collapse = ""), character(1L))),
      paste0(">FN999999 Packera lyrata\n", paste0(fungal_vec, collapse = "")),
      paste0(">JQ999999 Packera aurea\n", paste0(glabella_base, collapse = ""))
    )
    writeLines(fasta_lines, file.path(geno_dir, "public_packera_its_ets.fasta"))
    writeLines(c(">S1", paste0(mutate_k(glabella_base, c(35L)), collapse = "")), file.path(geno_dir, "skimmed_its_consensus.fasta"))
    rt_dir <- file.path(geno_dir, "public_reltime_trees"); dir.create(rt_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(c("#NEXUS", "Begin Taxa; Dimensions ntax=7; End;", "Begin Trees;", "  Translate",
                 "    1 Packera_aurea,", "    2 Packera_dubia,", "    3 Packera_obovata,", "    5 Packera_glabella,",
                 "    6 Packera_anonyma,", "    8 Packera_cana,", "    9 Packera_werneriifolia;",
                 "  tree TREE1 = [&R] ((((1:1.0,2[&rate=NA,divtime_95%_CI={0.5,2.0}]:1.0):1.0,3:2.0):1.0,5:3.0):1.0,(6:1.5,(8:1.0,9:1.0):0.5):2.5);", "End;"),
               file.path(rt_dir, "public_packera_RelTime_test_nexus.tre"))
    pass_df <- data.frame(
      genetic_sample_id = "S1", catalogNumber = "NCU_S1",
      species_standardized = "Packera glabella", species_raw = "Packera glabella", stringsAsFactors = FALSE
    )
    bb_res <- evaluate_public_backbone_and_organellar(pass_df, geno_dir, out_dir)
    stopifnot(length(bb_res$backbone$reltime_trees) == 1L, is.finite(ape::node.depth.edgelength(bb_res$backbone$reltime_trees[[1L]])[1L]))
    stopifnot(identical(as.character(bb_res$public_its_match[1]), "Packera glabella"))

    rogue_csv <- file.path(out_dir, "public_genbank_rogue_audit.csv")
    if (!file.exists(rogue_csv)) rogue_csv <- file.path(out_dir, "outputs", "genomics", "public_genbank_rogue_audit.csv")
    stopifnot(file.exists(rogue_csv))
    rogue_audit <- utils::read.csv(rogue_csv, stringsAsFactors = FALSE)
    stopifnot(ncol(rogue_audit) == 8L, nrow(rogue_audit) >= 2L)

    row_has_pair <- function(df_in, id_tok, reason_tok) {
      any(vapply(seq_len(nrow(df_in)), function(i) {
        vals <- as.character(unlist(df_in[i, , drop = TRUE]))
        any(grepl(id_tok, vals, fixed = TRUE)) && any(vals == reason_tok)
      }, logical(1L)))
    }
    stopifnot(row_has_pair(rogue_audit, "FN999999", "non_homologous_or_contaminant_accession"))
    stopifnot(row_has_pair(rogue_audit, "JQ999999", "mislabeled_genbank_conspecific_outlier"))
    cat("WARM_R_QUARANTINE_GATES_PASS\n")
    """
)


class OutlierQuarantineGateTests(unittest.TestCase):
    """Verify Step 01, Step 05, and Step 06 outlier quarantine gates."""

    def test_step01_centroid_and_institution_quarantine(self) -> None:
        """Verify Step 01 centroid/institution classification and SQLite round-trip."""
        from scripts.core import coord_hygiene_utils
        from scripts.core.db import get_db_connection, init_db_schema, upsert_vouchers_records

        classify_fn = _resolve_coord_classifier(coord_hygiene_utils)
        spatial_graph = {
            "counties": {"north carolina|orange": {"centroid_lat": 36.0610, "centroid_lon": -79.1205}}
        }
        cases = [
            (35.9108, -79.0505, {"locality": "Chapel Hill, Orange Co."}, (True, "herbarium_institution_rooftop")),
            (
                36.0610, -79.1205,
                {"uncertainty_m": None, "locality": "Orange County", "county": "Orange",
                 "state_prov": "North Carolina", "us_spatial_graph": spatial_graph},
                (True, "county_centroid_suspect"),
            ),
            (35.9000, -79.0500, {"locality": "", "uncertainty_m": None, "issues": set()}, (False, "")),
            (35.9108, -79.0505, {"locality": "diabase flatrock outcrop along Morgan Creek, 0.4 mi S of Hwy 54"}, (False, "")),
            (35.9108, -79.0505, {"locality": "Chapel Hill, Orange Co.", "is_type_specimen": True}, (False, "")),
        ]
        for lat, lon, kw, expected in cases:
            self.assertEqual(_call_coord_classifier(classify_fn, lat, lon, **kw), expected)

        with tempfile.TemporaryDirectory() as tmpdir:
            conn = get_db_connection(db_path=Path(tmpdir) / "quarantine_test.sqlite", init_schema=False)
            try:
                init_db_schema(conn)
                records = [
                    {"catalogNumber": "NCU_ROOFTOP_1", "species_standardized": "Packera dubia", "raw_scientificName": "Packera dubia",
                     "decimalLatitude": 35.9108, "decimalLongitude": -79.0505, "locality": "Chapel Hill, Orange Co.",
                     "georef_qc_status": "valid", "is_centroid_or_inst_coord": 1, "coord_quarantine_reason": "herbarium_institution_rooftop",
                     "is_type_specimen": 0},
                    {"catalogNumber": "NCU_CENTROID_1", "species_standardized": "Packera dubia", "raw_scientificName": "Packera dubia",
                     "decimalLatitude": 36.0610, "decimalLongitude": -79.1205, "locality": "Orange County", "county": "Orange",
                     "stateProvince": "North Carolina", "georef_qc_status": "valid", "is_centroid_or_inst_coord": 1,
                     "coord_quarantine_reason": "county_centroid_suspect", "is_type_specimen": 0},
                    {"catalogNumber": "NCU_VALID_1", "species_standardized": "Packera dubia", "raw_scientificName": "Packera dubia",
                     "decimalLatitude": 35.9000, "decimalLongitude": -79.0500, "locality": "", "georef_qc_status": "valid",
                     "is_centroid_or_inst_coord": 0, "coord_quarantine_reason": "", "is_type_specimen": 0},
                ]
                upsert_vouchers_records(conn, records)

                pi_cols = conn.execute("PRAGMA table_info(plant_individuals)").fetchall()
                col_names = [row[1] for row in pi_cols]
                insert_cols = ["catalogNumber", "plant_individual_id"] + (["is_primary_individual"] if "is_primary_individual" in col_names else [])
                insert_cols += [r[1] for r in pi_cols if r[1] not in ("catalogNumber", "plant_individual_id", "is_primary_individual") and r[3] and r[4] is None]
                type_map = {r[1]: (r[2] or "").upper() for r in pi_cols}
                placeholders, col_sql = ", ".join("?" for _ in insert_cols), ", ".join(insert_cols)
                for rec in records:
                    vals = [
                        rec["catalogNumber"] if c == "catalogNumber"
                        else (1 if c in ("plant_individual_id", "is_primary_individual")
                              else (0 if ("INT" in type_map[c] or "REAL" in type_map[c]) else ""))
                        for c in insert_cols
                    ]
                    conn.execute(f"INSERT OR REPLACE INTO plant_individuals ({col_sql}) VALUES ({placeholders})", vals)
                conn.execute(
                    "INSERT OR REPLACE INTO pipeline_runs (run_id, run_dir, start_timestamp) VALUES (?, ?, ?)",
                    ("test_run_quarantine", tmpdir, "2026-11-15T00:00:00Z"),
                )
                for rec in records:
                    conn.execute(
                        "INSERT OR REPLACE INTO run_active_cohort "
                        "(run_id, catalogNumber, plant_individual_id, species_standardized, determiner_tier, is_type_specimen, authority_weight) "
                        "VALUES (?, ?, ?, ?, ?, ?, ?)",
                        ("test_run_quarantine", rec["catalogNumber"], 1, rec["species_standardized"], "Tier_1_Gold", rec["is_type_specimen"], 2.0),
                    )
                conn.commit()

                for tbl in ("vouchers", "v_integrative_triage_master"):
                    rows = conn.execute(
                        f"SELECT catalogNumber, decimalLatitude, decimalLongitude, georef_qc_status, "
                        f"is_centroid_or_inst_coord, coord_quarantine_reason FROM {tbl} "
                        f"WHERE catalogNumber IN ('NCU_ROOFTOP_1', 'NCU_CENTROID_1', 'NCU_VALID_1')"
                    ).fetchall()
                    by_id = {r[0]: r for r in rows}
                    self.assertEqual(set(by_id.keys()), {"NCU_ROOFTOP_1", "NCU_CENTROID_1", "NCU_VALID_1"})
                    for cid, exp_flag, exp_reason in (
                        ("NCU_ROOFTOP_1", 1, "herbarium_institution_rooftop"),
                        ("NCU_CENTROID_1", 1, "county_centroid_suspect"),
                        ("NCU_VALID_1", 0, ""),
                    ):
                        _, lat_v, lon_v, qc_status, q_flag, q_reason = by_id[cid]
                        self.assertIsNotNone(lat_v)
                        self.assertIsNotNone(lon_v)
                        self.assertEqual(qc_status, "valid")
                        self.assertEqual(int(q_flag), exp_flag)
                        self.assertEqual(q_reason or "", exp_reason)
            finally:
                conn.close()

    @unittest.skipIf(
        FAST_AUDIT or RSCRIPT_BIN is None,
        "Warm-R quarantine test skipped during --fast audit or when Rscript is unavailable.",
    )
    def test_step05_and_step06_outlier_quarantine_warm_r(self) -> None:
        """Verify Step 05 MCD/azonal quarantine and Step 06 tree/FASTA rogue quarantine in one R session."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            r_script_path = tmp_path / "verify_quarantine_gates.R"
            r_code = (
                _R_QUARANTINE_SCRIPT
                .replace("__REPO_ROOT__", PROJECT_ROOT.as_posix())
                .replace("__TMP_DIR__", tmp_path.as_posix())
            )
            r_script_path.write_text(r_code, encoding="utf-8")
            proc = subprocess.run(
                [RSCRIPT_BIN, "--vanilla", str(r_script_path)],
                cwd=str(PROJECT_ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                self.fail(
                    f"Warm-R quarantine test failed (exit {proc.returncode}):\n"
                    f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
                )
            self.assertIn("WARM_R_QUARANTINE_GATES_PASS", proc.stdout)


if __name__ == "__main__":
    unittest.main()
