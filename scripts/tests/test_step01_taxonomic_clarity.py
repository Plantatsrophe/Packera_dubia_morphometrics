#!/usr/bin/env python3
"""Step 01 taxonomic clarity, determiner stratification, concept tracking, rescue, and R handshake tests."""

from __future__ import annotations

import os
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core import taxonomy_authority as ta
from scripts.core.db import get_db_connection, upsert_rows_generic, upsert_vouchers_records
from scripts.core.db_schema import PLANTMATICS_SCHEMA_USER_VERSION, init_db_schema
from scripts.core.harvester import EXPORT_COLUMNS, VoucherHarvester, export_curated_table
from scripts.core.step01_concept_authority_utils import (
    apply_concept_and_segregate_recency_caps,
    extract_identification_year,
    extract_typified_basionym,
    has_open_nomenclature_qualifier,
    is_interspecific_hybrid_formula,
    parse_voucher_taxonomic_concept,
)
from scripts.tests._r_session import rscript_available, run_warm_r_cmd

# Five new Tier A SQLite taxonomic clarity columns under schema version 20261127.
CONCEPT_COLUMNS = (
    "typified_basionym",
    "verbatim_taxon_concept",
    "concept_synonymy_status",
    "identification_year",
    "is_open_or_hybrid_concept",
)


def _occ_rec(key: int, cat: str, sci: str, state: str = "North Carolina", lat: float | None = 35.9, lon: float | None = -79.0, **extra) -> dict:
    """Build a minimal Darwin Core / GBIF occurrence record for offline harvester tests."""
    rec = {
        "key": key,
        "catalogNumber": cat,
        "institutionCode": "NCU",
        "scientificName": sci,
        "stateProvince": state,
        "county": "Orange",
        "decimalLatitude": lat,
        "decimalLongitude": lon,
        "coordinateUncertaintyInMeters": 50.0,
        "year": 2020,
        "month": 5,
        "day": 10,
        "recordedBy": f"Collector {key}",
        "recordNumber": str(key),
        "identifiedBy": "T.M. Barkley",
        "media": [{"identifier": f"https://example.org/{cat}.jpg", "type": "StillImage", "format": "image/jpeg"}],
    }
    rec.update(extra)
    return rec


class TestStep01TaxonomicClarity(unittest.TestCase):
    """Verify Step 01 taxonomic concept separation, determiner tiers, rescue passes, and R integration."""

    def test_specialist_stratification_estes_and_regional_silver(self) -> None:
        """Verify qualified L. D./Dwayne Estes -> Gold, bare J. Estes -> Bronze, regional/historical -> Silver (Gold on Holotypes)."""
        for det in ("L. D. Estes", "Dwayne Estes", "Estes, L. D.", "C. C. Freeman", "W. A. Weber"):
            tier, _, _ = ta.parse_determiner_tier(None, det, "Collector", None, "SMALL", "Meadow", None)
            self.assertEqual(tier, "Tier_1_Gold", f"{det} should map to Tier_1_Gold")

        for bare in ("J. Estes", "S. Freeman", "J. Packer"):
            tier, _, _ = ta.parse_determiner_tier(None, bare, "Collector", None, "SMALL", "Meadow", None)
            self.assertEqual(tier, "Tier_3_Bronze", f"Bare surname '{bare}' must remain Tier_3_Bronze")

        for reg in ("A. Cronquist", "M. L. Fernald", "A. Gray", "B. L. Turner", "C. Jeffrey", "J. K. Small"):
            tier, _, _ = ta.parse_determiner_tier(None, reg, "Collector", None, "SMALL", "Meadow", None)
            self.assertEqual(tier, "Tier_2_Silver", f"{reg} must map to Tier_2_Silver")
            tier_h, _, _ = ta.parse_determiner_tier("Holotype", reg, "Collector", None, "SMALL", "Meadow", None)
            self.assertEqual(tier_h, "Tier_1_Gold", f"Holotype by {reg} must promote to Tier_1_Gold")

    def test_hybrid_formulas_vs_nothospecies_and_open_nomenclature_quarantine(self) -> None:
        """Verify two-epithet hybrid formulas vs nothospecies xmemmingeri and cf./? Tier_3_Bronze quarantine."""
        for f_str, exp_norm, exp_p1 in (
            ("Packera anonyma × Packera dubia", "Packera anonyma × Packera dubia", "Packera anonyma"),
            ("Packera aurea x obovata", "Packera aurea × Packera obovata", "Packera aurea"),
            ("Packera anonyma x P. paupercula", "Packera anonyma × Packera paupercula", "Packera anonyma"),
            ("Senecio tomentosus × anonymus", "Senecio tomentosus × Senecio anonymus", "Senecio tomentosus"),
        ):
            is_h, norm, p1 = is_interspecific_hybrid_formula(f_str)
            self.assertTrue(is_h, f_str)
            self.assertEqual((norm, p1), (exp_norm, exp_p1))

        for notho in ("Packera ×memmingeri", "Packera x memmingeri"):
            self.assertEqual(is_interspecific_hybrid_formula(notho), (False, "", None))

        self.assertTrue(has_open_nomenclature_qualifier(raw_name="Packera cf. dubia"))
        self.assertTrue(has_open_nomenclature_qualifier(identification_qualifier="aff."))
        self.assertFalse(has_open_nomenclature_qualifier(raw_name="Packera ×memmingeri"))

        for q_name in ("Packera cf. dubia", "Packera dubia ?"):
            c_info = parse_voucher_taxonomic_concept({"scientificName": q_name}, None, False, {"Packera", "Senecio", "Cineraria"})
            self.assertEqual(c_info["species_standardized_hint"], "Packera dubia")
            self.assertEqual(c_info["concept_synonymy_status"], "Open_Nomenclature")
            self.assertTrue(c_info["is_open_or_hybrid_concept"])

        df_q = pd.DataFrame([
            {"catalogNumber": "Q1", "scientificName": "Packera cf. dubia", "species_standardized": "Packera dubia", "identifiedBy": "T. M. Barkley", "type_status": "", "institutionCode": "NCU", "concept_synonymy_status": "Open_Nomenclature", "is_open_or_hybrid_concept": True},
            {"catalogNumber": "H1", "scientificName": "Packera anonyma × Packera dubia", "species_standardized": "Packera anonyma", "identifiedBy": "T. M. Barkley", "type_status": "", "institutionCode": "NCU", "concept_synonymy_status": "Hybrid_Formula", "is_open_or_hybrid_concept": True},
            {"catalogNumber": "T1", "scientificName": "Packera cf. dubia", "species_standardized": "Packera dubia", "identifiedBy": "T. M. Barkley", "type_status": "Holotype", "institutionCode": "NCU", "concept_synonymy_status": "Open_Nomenclature", "is_open_or_hybrid_concept": True},
        ])
        tiered = ta.apply_corpus_authority_tiers(df_q).set_index("catalogNumber")
        self.assertEqual(tiered.loc["Q1", "determiner_tier"], "Tier_3_Bronze")
        self.assertEqual(tiered.loc["H1", "determiner_tier"], "Tier_3_Bronze")
        self.assertEqual(tiered.loc["T1", "determiner_tier"], "Tier_1_Gold")

    def test_typified_basionym_unconditional_and_country_none_type_rescue(self) -> None:
        """Verify typified_basionym priority and unconditional/country=None type rescue past max_records_per_taxon."""
        self.assertEqual(extract_typified_basionym({"typifiedName": "Senecio tomentosus", "scientificName": "Packera dubia"}, True), "Senecio tomentosus")
        self.assertEqual(extract_typified_basionym({"typeStatus": "Holotype of Senecio alabamensis Britton", "scientificName": "Packera dubia"}, True), "Senecio alabamensis")
        self.assertEqual(extract_typified_basionym({"typifiedName": "Senecio tomentosus"}, False), "")

        nontype_rec = _occ_rec(4001, "NT_01", "Packera dubia", identifiedBy="Unknown Collector")
        global_holo = _occ_rec(
            4002, "TYPE_GLOBAL_01", "Packera dubia", state="Alabama", lat=None, lon=None,
            typeStatus="Holotype of Senecio alabamensis Britton", identifiedBy="N. L. Britton",
        )
        global_holo["countryCode"] = None

        def _mock_gbif(**kw):
            if kw.get("typeStatus") == "*":
                # Only return the holotype on the country=None global rescue query.
                return {"results": [global_holo], "count": 1} if "country" not in kw else {"results": [], "count": 0}
            return {"results": [nontype_rec], "count": 1} if (kw.get("hasCoordinate") is not False and kw.get("country") == "US") else {"results": [], "count": 0}

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d)
            with patch("scripts.core.harvester.occ.search", side_effect=_mock_gbif):
                harv = VoucherHarvester(
                    taxa=["Packera dubia"], max_records_per_taxon=1, expand_synonyms=False,
                    workspace_dir=ws, db_path=ws / "p.sqlite", raw_dir=ws / "raw",
                    output_csv=ws / "curated.csv", check_sharpness=False, usda_county_filter=False,
                )
                df_h = harv.harvest()
        self.assertIn("NT_01", set(df_h["catalogNumber"]))
        self.assertIn("TYPE_GLOBAL_01", set(df_h["catalogNumber"]))
        row_t = df_h.set_index("catalogNumber").loc["TYPE_GLOBAL_01"]
        self.assertEqual(row_t["typified_basionym"], "Senecio alabamensis")
        self.assertEqual(row_t["determiner_tier"], "Tier_1_Gold")
        self.assertEqual(row_t["concept_synonymy_status"], "Heterotypic_Synonym")

    def test_autonyms_recency_cap_sqlite_migration_and_23col_csv(self) -> None:
        """Verify autonym stripping, pre-2014 P. serpenticola recency cap, schema 20261127 migration, and CSV/SQLite boundary."""
        aut = parse_voucher_taxonomic_concept({"scientificName": "Packera aurea var. aurea", "dateIdentified": "2018-06-15"}, "Packera aurea", False, {"Packera"})
        self.assertEqual((aut["infraspecific_taxon"], aut["verbatim_taxon_concept"], aut["concept_synonymy_status"], aut["identification_year"]), ("", "Packera aurea", "Accepted_Sensu_Stricto", 2018))
        var = parse_voucher_taxonomic_concept({"scientificName": "Packera aurea var. aquilonius"}, "Packera aurea", False, {"Packera"})
        self.assertEqual((var["infraspecific_taxon"], var["concept_synonymy_status"]), ("var. aquilonius", "Infraspecific_Taxon"))
        frm = parse_voucher_taxonomic_concept({"scientificName": "Packera aurea f. tristis"}, "Packera aurea", False, {"Packera"})
        self.assertEqual(frm["infraspecific_taxon"], "")
        self.assertEqual(extract_identification_year(identified_by_raw="T. M. Barkley, 1988"), 1988)

        df_cap = pd.DataFrame([
            {"catalogNumber": "PRE_NC", "species_standardized": "Packera anonyma", "determiner_tier": "Tier_1_Gold", "stateProvince": "North Carolina", "county": "Burke County", "identification_year": 1990, "year": 1990, "type_status": "", "concept_synonymy_status": "Accepted_Sensu_Stricto", "is_open_or_hybrid_concept": False},
            {"catalogNumber": "POST_NC", "species_standardized": "Packera anonyma", "determiner_tier": "Tier_1_Gold", "stateProvince": "North Carolina", "county": "Burke County", "identification_year": 2016, "year": 2016, "type_status": "", "concept_synonymy_status": "Accepted_Sensu_Stricto", "is_open_or_hybrid_concept": False},
            {"catalogNumber": "WAKE_NC", "species_standardized": "Packera anonyma", "determiner_tier": "Tier_1_Gold", "stateProvince": "North Carolina", "county": "Wake County", "identification_year": 1990, "year": 1990, "type_status": "", "concept_synonymy_status": "Accepted_Sensu_Stricto", "is_open_or_hybrid_concept": False},
            {"catalogNumber": "HET_S", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "stateProvince": "Georgia", "county": "Clarke", "identification_year": 2019, "year": 2019, "type_status": "", "concept_synonymy_status": "Heterotypic_Synonym", "is_open_or_hybrid_concept": False},
        ])
        capped = apply_concept_and_segregate_recency_caps(df_cap).set_index("catalogNumber")
        self.assertEqual(capped.loc["PRE_NC", "determiner_tier"], "Tier_2_Silver")
        self.assertEqual(capped.loc["POST_NC", "determiner_tier"], "Tier_1_Gold")
        self.assertEqual(capped.loc["WAKE_NC", "determiner_tier"], "Tier_1_Gold")
        self.assertEqual(capped.loc["HET_S", "determiner_tier"], "Tier_2_Silver")

        # Legacy 20261127 DB lacking the 5 new columns must migrate idempotently without version bump.
        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d)
            db_p = ws / "legacy_20261127.sqlite"
            c_pre = get_db_connection(db_path=db_p, init_schema=True, workspace_root=ws)
            try:
                old_cols = [r[1] for r in c_pre.execute("PRAGMA table_info(vouchers);").fetchall() if r[1] not in CONCEPT_COLUMNS]
                c_pre.execute("PRAGMA foreign_keys = OFF;")
                c_pre.execute("DROP VIEW IF EXISTS v_integrative_triage_master;")
                c_pre.execute(f"CREATE TABLE _v_pre AS SELECT {', '.join(old_cols)} FROM vouchers;")
                c_pre.execute("DROP TABLE vouchers;")
                c_pre.execute("ALTER TABLE _v_pre RENAME TO vouchers;")
                c_pre.execute("CREATE VIEW v_integrative_triage_master AS SELECT catalogNumber FROM vouchers;")
                c_pre.execute(f"PRAGMA user_version = {PLANTMATICS_SCHEMA_USER_VERSION};")
                c_pre.commit()
            finally:
                c_pre.close()

            conn2 = get_db_connection(db_path=db_p, init_schema=True, workspace_root=ws)
            try:
                self.assertEqual(conn2.execute("PRAGMA user_version;").fetchone()[0], PLANTMATICS_SCHEMA_USER_VERSION)
                v_cols = {r[1] for r in conn2.execute("PRAGMA table_info(vouchers);").fetchall()}
                mv_cols = {r[1] for r in conn2.execute("PRAGMA table_info(v_integrative_triage_master);").fetchall()}
                for c in CONCEPT_COLUMNS:
                    self.assertIn(c, v_cols)
                    self.assertIn(c, mv_cols)
            finally:
                conn2.close()

            # Verify 23 EXPORT_COLUMNS contract and SQLite-only usda_symbol/usda_county_status boundary.
            self.assertEqual(len(EXPORT_COLUMNS), 23)
            out_csv = ws / "curated.csv"
            row_exp = {col: "" for col in EXPORT_COLUMNS}
            row_exp.update({
                "catalogNumber": "EXP_1", "institutionCode": "NCU", "scientificName": "Packera dubia",
                "species_raw": "Senecio tomentosus", "species_standardized": "Packera dubia",
                "typified_basionym": "Senecio tomentosus", "verbatim_taxon_concept": "Senecio tomentosus",
                "concept_synonymy_status": "Homotypic_Synonym", "identification_year": 2021,
                "is_open_or_hybrid_concept": False, "usda_symbol": "PATO4", "usda_county_status": "in_range",
            })
            export_curated_table(pd.DataFrame([row_exp]), output_path=out_csv, db_path=db_p)
            exported_df = pd.read_csv(out_csv)
            self.assertEqual(list(exported_df.columns[: len(EXPORT_COLUMNS)]), EXPORT_COLUMNS)
            self.assertNotIn("usda_symbol", exported_df.columns)
            self.assertNotIn("usda_county_status", exported_df.columns)
            for c in CONCEPT_COLUMNS:
                self.assertIn(c, exported_df.columns)

    def test_has_coordinate_false_cis_mississippi_rescue(self) -> None:
        """Verify hasCoordinate=False pass rescues cis-Mississippi state sheets while dropping Western states."""
        nc_ungeo = _occ_rec(5001, "NC_UNG_1", "Packera dubia", state="North Carolina", lat=None, lon=None)
        or_ungeo = _occ_rec(5002, "OR_UNG_1", "Packera dubia", state="Oregon", lat=None, lon=None)

        def _mock_rescue(**kw):
            if kw.get("typeStatus") == "*":
                return {"results": [], "count": 0}
            if kw.get("hasCoordinate") is False and kw.get("country") == "US":
                return {"results": [nc_ungeo, or_ungeo], "count": 2}
            return {"results": [], "count": 0}

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d)
            with patch("scripts.core.harvester.occ.search", side_effect=_mock_rescue):
                harv = VoucherHarvester(
                    taxa=["Packera dubia"], max_records_per_taxon=5, expand_synonyms=False,
                    workspace_dir=ws, db_path=ws / "p.sqlite", raw_dir=ws / "raw",
                    output_csv=ws / "curated.csv", check_sharpness=False, usda_county_filter=False,
                    include_ungeoreferenced_state_vouchers=True,
                )
                df_r = harv.harvest()
        self.assertIn("NC_UNG_1", set(df_r["catalogNumber"]))
        self.assertNotIn("OR_UNG_1", set(df_r["catalogNumber"]))

    def test_stratify_duplicates_tier_range_rank_and_singletons(self) -> None:
        """Verify tier outranks USDA range, in_range breaks Gold ties, and incomplete dates/unknown collectors/centroids stay singletons."""
        df_dups = pd.DataFrame([
            # Duplicate 1: Tier_1_Gold out_of_range MUST outrank Tier_3_Bronze in_range
            {"catalogNumber": "D1_GOLD_OOR", "institutionCode": "NCU", "recordedBy": "A. Collector", "recordNumber": "101", "eventDate": "2010-05-01", "year": 2010, "decimalLatitude": 35.90, "decimalLongitude": -79.05, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "usda_county_status": "out_of_range"},
            {"catalogNumber": "D1_BRONZE_IN", "institutionCode": "US", "recordedBy": "A. Collector", "recordNumber": "101", "eventDate": "2010-05-01", "year": 2010, "decimalLatitude": 35.90, "decimalLongitude": -79.05, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze", "usda_county_status": "in_range"},
            # Duplicate 2: Within Tier_1_Gold, in_range outranks out_of_range
            {"catalogNumber": "D2_GOLD_OOR", "institutionCode": "DUKE", "recordedBy": "B. Collector", "recordNumber": "202", "eventDate": "2011-05-01", "year": 2011, "decimalLatitude": 35.95, "decimalLongitude": -79.15, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "usda_county_status": "out_of_range"},
            {"catalogNumber": "D2_GOLD_IN", "institutionCode": "NY", "recordedBy": "B. Collector", "recordNumber": "202", "eventDate": "2011-05-01", "year": 2011, "decimalLatitude": 35.95, "decimalLongitude": -79.15, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "usda_county_status": "in_range"},
            # Pass-2 singleton guards (< 50 m proximity must NOT merge -00-00 dates, unknown collectors, or centroids)
            {"catalogNumber": "S_UNK1", "institutionCode": "NCU", "recordedBy": "unknown", "recordNumber": "s.n.", "eventDate": "1965-05-01", "year": 1965, "decimalLatitude": 36.0000, "decimalLongitude": -79.0000, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
            {"catalogNumber": "S_UNK2", "institutionCode": "US", "recordedBy": "unknown", "recordNumber": "s.n.", "eventDate": "1965-05-01", "year": 1965, "decimalLatitude": 36.0001, "decimalLongitude": -79.0000, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
            {"catalogNumber": "S_DT1", "institutionCode": "NCU", "recordedBy": "J. Doe", "recordNumber": "s.n.", "eventDate": "1965-00-00", "year": 1965, "decimalLatitude": 36.1000, "decimalLongitude": -79.1000, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
            {"catalogNumber": "S_DT2", "institutionCode": "US", "recordedBy": "J. Doe", "recordNumber": "s.n.", "eventDate": "1965-00-00", "year": 1965, "decimalLatitude": 36.1001, "decimalLongitude": -79.1000, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
            {"catalogNumber": "S_CEN1", "institutionCode": "NCU", "recordedBy": "R. Roe", "recordNumber": "s.n.", "eventDate": "1972-04-12", "year": 1972, "decimalLatitude": 36.2000, "decimalLongitude": -79.2000, "is_centroid_or_inst_coord": True, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
            {"catalogNumber": "S_CEN2", "institutionCode": "US", "recordedBy": "R. Roe", "recordNumber": "s.n.", "eventDate": "1972-04-12", "year": 1972, "decimalLatitude": 36.2001, "decimalLongitude": -79.2000, "is_centroid_or_inst_coord": True, "scientificName": "Packera dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze"},
        ])
        strat = ta.stratify_duplicates(df_dups).set_index("catalogNumber")
        self.assertTrue(bool(strat.loc["D1_GOLD_OOR", "is_primary_duplicate"]))
        self.assertFalse(bool(strat.loc["D1_BRONZE_IN", "is_primary_duplicate"]))
        self.assertTrue(bool(strat.loc["D2_GOLD_IN", "is_primary_duplicate"]))
        self.assertFalse(bool(strat.loc["D2_GOLD_OOR", "is_primary_duplicate"]))
        for sing in ("S_UNK1", "S_UNK2", "S_DT1", "S_DT2", "S_CEN1", "S_CEN2"):
            self.assertTrue(bool(strat.loc[sing, "is_primary_duplicate"]), f"Singleton {sing} must stay unmerged")
        self.assertEqual(strat.loc[["S_UNK1", "S_UNK2", "S_DT1", "S_DT2", "S_CEN1", "S_CEN2"], "exsiccatae_key"].nunique(), 6)

    def test_r_authority_concept_and_db_roundtrip_warm_r(self) -> None:
        """Batch-verify R authority weighting, type basionym controls, and db_utils.R concept roundtrip in one warm R call."""
        if not rscript_available():
            self.skipTest("Rscript unavailable or PLANTMATICS_FAST_AUDIT=1")
        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d)
            db_p = ws / "clarity_r.sqlite"
            run_id = "run_clarity_r_01"
            conn = get_db_connection(db_path=db_p, init_schema=True, workspace_root=ws)
            try:
                conn.execute("INSERT OR IGNORE INTO pipeline_runs (run_id) VALUES (?);", (run_id,))
                conn.commit()
                upsert_vouchers_records(conn, [{
                    "catalogNumber": "RDB_01", "institutionCode": "NCU", "scientificName": "Packera dubia",
                    "species_raw": "Packera cf. dubia", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze",
                    "typified_basionym": "", "verbatim_taxon_concept": "Packera dubia",
                    "concept_synonymy_status": "Open_Nomenclature", "identification_year": 2005, "is_open_or_hybrid_concept": True,
                }])
                upsert_rows_generic(
                    conn, "plant_individuals",
                    [{"catalogNumber": "RDB_01", "plant_individual_id": 1}],
                    pk_cols=("catalogNumber", "plant_individual_id"),
                )
                upsert_rows_generic(
                    conn, "otu_morphometrics",
                    [{"run_id": run_id, "catalogNumber": "RDB_01", "plant_individual_id": 1, "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze", "EFA_PC1": 0.1, "EFA_PC2": 0.2}],
                    pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
                )
            finally:
                conn.close()

            r_code = (
                "source('scripts/morphometrics/authority_weighting_utils.R', local = FALSE)\n"
                "source('scripts/core/db_utils.R', local = FALSE)\n"
                "set.seed(42L)\n"
                "n_gold_d <- 14L\n"
                "X_d <- cbind(rnorm(n_gold_d, 0, 0.04), rnorm(n_gold_d, 0, 0.04))\n"
                "X_toy <- rbind(X_d, c(0.01, -0.01), c(8.0, 8.0), c(0.05, 0.05), c(0.02, -0.02), c(9.5, 9.5), c(-9.0, 8.5), c(4.0, 4.0), c(4.03, 3.98), c(3.97, 4.02))\n"
                "colnames(X_toy) <- c('EFA_PC1', 'EFA_PC2')\n"
                "n_tot <- nrow(X_toy)\n"
                "df_r <- data.frame(\n"
                "  catalogNumber = sprintf('R_%02d', seq_len(n_tot)),\n"
                "  species_standardized = c(rep('Packera dubia', 18L), rep('Packera aurea', 5L)),\n"
                "  determiner_tier = c(rep('Tier_1_Gold', 17L), 'Tier_3_Bronze', 'Tier_1_Gold', 'Tier_2_Silver', 'Tier_3_Bronze', 'Tier_3_Bronze', 'Tier_3_Bronze'),\n"
                "  type_status = c(rep('', n_gold_d), 'Holotype', 'Paratype', rep('', 7L)),\n"
                "  concept_synonymy_status = c(rep('Accepted_Sensu_Stricto', 16L), 'Heterotypic_Synonym', 'Accepted_Sensu_Stricto', 'Open_Nomenclature', 'Hybrid_Formula', rep('Accepted_Sensu_Stricto', 3L)),\n"
                "  is_open_or_hybrid_concept = c(rep(FALSE, 18L), TRUE, TRUE, rep(FALSE, 3L)),\n"
                "  verbatim_taxon_concept = c(rep('Packera dubia', 16L), 'Senecio alabamensis', 'Packera dubia', 'Packera cf. aurea', 'Packera aurea × Packera dubia', 'Packera ×memmingeri', 'Packera aurea', 'Packera aurea'),\n"
                "  species_raw = c(rep('Packera dubia', n_gold_d), 'Senecio tomentosus', 'Packera dubia', 'Senecio alabamensis', 'Packera dubia', 'Packera cf. aurea', 'Packera aurea x obovata', 'Packera ×memmingeri', 'Packera aurea', 'Packera aurea'),\n"
                "  typified_basionym = c(rep('', n_gold_d), 'Senecio tomentosus', 'Packera dubia', rep('', 7L)),\n"
                "  EFA_PC1 = X_toy[, 1L], EFA_PC2 = X_toy[, 2L], leaf_quality_weight = rep(1.0, n_tot), stringsAsFactors = FALSE\n"
                ")\n"
                "tiers <- .standardize_determiner_tier_vec(df_r)\n"
                "stopifnot(tiers[19L] == 'Tier_3_Bronze', tiers[20L] == 'Tier_3_Bronze', tiers[15L] == 'Tier_1_Gold')\n"
                "tv <- .identify_type_specimen_vec(df_r)\n"
                "stopifnot(isTRUE(tv$is_name_bearing_type[15L]), !isTRUE(tv$is_name_bearing_type[16L]))\n"
                "anc <- select_adaptive_reference_anchors(df_r, feature_matrix = X_toy, min_anchors_per_taxon = 5L)\n"
                "stopifnot(!anc$active_mask[19L], !anc$active_mask[20L], anc$active_mask[15L], !anc$active_mask[16L], anc$active_mask[21L])\n"
                "w_a <- compute_empirical_bayes_authority_weights(df_r, mode = 'active_passive', anchor_info = anc)\n"
                "w_c <- compute_empirical_bayes_authority_weights(df_r, mode = 'continuous_prior_posterior', anchor_info = anc)\n"
                "stopifnot(w_a[19L] == 0.0, w_a[20L] == 0.0, w_c[19L] == 0.0, w_c[20L] == 0.0)\n"
                "df_bas <- data.frame(\n"
                "  catalogNumber = c('B1', 'B2', 'B3', 'B4'), species_standardized = rep('Packera dubia', 4L),\n"
                "  type_status = c('Holotype', 'Holotype of Senecio alabamensis Britton', 'Isotype', 'Syntype'),\n"
                "  typified_basionym = c('Senecio tomentosus', '', '', ''),\n"
                "  concept_synonymy_status = c('Homotypic_Synonym', 'Heterotypic_Synonym', 'Homotypic_Synonym', 'Heterotypic_Synonym'),\n"
                "  backbone_basionym_name = c('Cineraria dubia Willd.', 'Senecio alabamensis Britton', 'Senecio tomentosus Michx.', 'Senecio anonymus Wood'),\n"
                "  species_raw = c('Packera dubia', 'Packera dubia', 'Cineraria dubia', 'Packera anonyma'), stringsAsFactors = FALSE\n"
                ")\n"
                "stopifnot(identical(.canonicalize_type_basionym_vec(df_bas)[1:2], c('Senecio tomentosus', 'Senecio alabamensis')))\n"
                "stopifnot(identical(.is_homotypic_sensu_stricto_basionym(c('cineraria dubia', 'senecio dubius', 'senecio tomentosus', 'packera tomentosa', 'Senecio alabamensis'), rep('Packera dubia', 5L)), c(TRUE, TRUE, TRUE, TRUE, FALSE)))\n"
                "stopifnot(identical(.filter_homotypic_type_rows(df_bas, 1:4, 'Packera dubia'), c(1L, 3L)))\n"
                "X_ctrl <- rbind(cbind(rnorm(6, 0, 0.05), rnorm(6, 0, 0.05)), c(0.02, -0.02), c(6.0, 6.0))\n"
                "df_ctrl <- data.frame(catalogNumber = sprintf('C%d', 1:8), species_standardized = rep('Packera dubia', 8L), determiner_tier = rep('Tier_1_Gold', 8L), type_status = c(rep('', 6L), 'Holotype', 'Holotype'), typified_basionym = c(rep('', 6L), 'Senecio tomentosus', 'Senecio alabamensis'), concept_synonymy_status = c(rep('Accepted_Sensu_Stricto', 6L), 'Homotypic_Synonym', 'Heterotypic_Synonym'), EFA_PC1 = X_ctrl[, 1L], EFA_PC2 = X_ctrl[, 2L], stringsAsFactors = FALSE)\n"
                "ctrl_res <- compute_type_specimen_morph_control(df_ctrl, X_mat = X_ctrl)\n"
                "stopifnot(all(ctrl_res$has_empirical_type_anchor), abs(ctrl_res$type_centroids['Packera dubia', 1L]) < 0.5)\n"
                "req_v_specs <- c('typified_basionym', 'verbatim_taxon_concept', 'concept_synonymy_status', 'identification_year', 'is_open_or_hybrid_concept')\n"
                "stopifnot(all(req_v_specs %in% names(.dynamic_sqlite_table_specs_r('vouchers'))), 'is_open_or_hybrid_concept' %in% PACKERA_DB_LOGICAL_COLS)\n"
                "coerced_oh <- coerce_sqlite_logical_cols(data.frame(is_open_or_hybrid_concept = c(1L, 0L)))$is_open_or_hybrid_concept\n"
                "stopifnot(identical(coerced_oh, c(TRUE, FALSE)))\n"
                "if (requireNamespace('DBI', quietly = TRUE) && requireNamespace('RSQLite', quietly = TRUE)) {\n"
                "  df_db <- fetch_step04_morphometrics_from_db(run_id = Sys.getenv('RUN_ID'))\n"
                "  stopifnot(is.data.frame(df_db), nrow(df_db) == 1L, all(req_v_specs %in% names(df_db)), is.logical(df_db$is_open_or_hybrid_concept), isTRUE(df_db$is_open_or_hybrid_concept[1L]))\n"
                "}\n"
                "cat('WARM_R_STEP01_CLARITY_OK\\n')\n"
            )
            r_env = {**os.environ, "PLANTMATICS_DB_PATH": str(db_p), "PACKERA_DB_PATH": str(db_p), "RUN_ID": run_id}
            res = run_warm_r_cmd(["Rscript", "-e", r_code], stage_key="step04", cwd=PROJECT_ROOT, env=r_env)
            self.assertEqual(res.returncode, 0, f"R failed:\nSTDOUT:{res.stdout}\nSTDERR:{res.stderr}")
            self.assertIn("WARM_R_STEP01_CLARITY_OK", res.stdout)


if __name__ == "__main__":
    unittest.main()
