#!/usr/bin/env python3
"""Offline mocked contract tests for external ENA, taxonomic backbone, geochem/wetland, and R utility sources."""

from __future__ import annotations

import importlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection
from scripts.core.usda_plants_ranges import filter_vouchers_by_usda_county_range
from scripts.data_prep import (
    download_geochem_wetland_layers as dl_geochem,
    download_public_packera_sequences as dl_pub,
    download_taxonomic_backbone as dl_backbone,
    ena_paraloss_discovery as ena_mod,
)


def _atomic_csv_writer(root_dir: Path, rel_name: str, payload, force: bool = False) -> Path:
    dest = Path(root_dir) / rel_name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if force or not dest.is_file():
        data = payload.encode("utf-8") if isinstance(payload, str) else (payload.read_bytes() if isinstance(payload, Path) else bytes(payload))
        dest.write_bytes(data)
    return dest


def _minimal_counties_geojson(path: Path) -> Path:
    payload = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "county_key": "north carolina|clay",
                    "state_name": "north carolina",
                    "county_norm": "clay",
                    "centroid_lon": -83.62,
                    "centroid_lat": 35.08,
                },
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [[[[-83.7, 35.0], [-83.5, 35.0], [-83.5, 35.2], [-83.7, 35.2], [-83.7, 35.0]]]],
                },
            },
            {
                "type": "Feature",
                "properties": {
                    "county_key": "california|alameda",
                    "state_name": "california",
                    "county_norm": "alameda",
                    "centroid_lon": -122.2,
                    "centroid_lat": 37.7,
                },
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [[[[-122.4, 37.6], [-122.0, 37.6], [-122.0, 37.9], [-122.4, 37.9], [-122.4, 37.6]]]],
                },
            },
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class TestExternalDataSourcesOffline(unittest.TestCase):
    """Offline mocked unit and R integration checks for external data provisioning."""

    def test_dynamic_ena_tax_tree_and_paraloss1272_weighting(self) -> None:
        self.assertIn("PRJNA1046046", dl_pub._SRA_BIOPROJECTS)

        def mock_http(url: str, *args, **kwargs) -> bytes:
            low = url.lower()
            if "ena/taxonomy/rest/scientific-name/packera" in low:
                return json.dumps([{"taxId": "189364", "scientificName": "Packera"}]).encode("utf-8")
            if "ena/taxonomy/rest/scientific-name/solidago" in low:
                return json.dumps([{"taxId": "59294", "scientificName": "Solidago"}]).encode("utf-8")
            if "ena/portal/api/search" in low and "tax_tree(189364)" in low:
                tsv = (
                    "study_accession\trun_accession\tsample_accession\tscientific_name\tlibrary_name\tspecimen_voucher\tfastq_ftp\tfastq_bytes\n"
                    "PRJNA1046046\tSRR2900001\tSAMN2900001\tPackera dubia\tHybSeq_1\tNCU001\tftp.sra/1.fq.gz\t12000\n"
                    "PRJNA1046046\tSRR2900002\tSAMN2900002\tPackera glabella\tHybSeq_2\tNCU002\tftp.sra/2.fq.gz\t14000\n"
                )
                return tsv.encode("utf-8")
            if "ena/portal/api/search" in low and "tax_tree(59294)" in low:
                tsv = (
                    "study_accession\trun_accession\tsample_accession\tscientific_name\tlibrary_name\tspecimen_voucher\tfastq_ftp\tfastq_bytes\n"
                    "PRJNA999001\tSRR9990001\tSAMN9990001\tSolidago altissima\tSol_1\tDUKE99\tftp.sra/s1.fq.gz\t11000\n"
                )
                return tsv.encode("utf-8")
            if "api.figshare.com" in low:
                return json.dumps({"files": [{"name": "Compositae-ParaLoss-1272.csv", "download_url": "https://mock.local/paraloss1272.csv"}]}).encode("utf-8")
            if "paraloss" in low or low.endswith((".csv", ".tsv")):
                csv_txt = (
                    "locus_id,gene_family,paralogy_class,target_length_bp\n"
                    "PL1272_0001,ParaLoss1272,single_copy,950\n"
                    "PL1272_0002,ParaLoss1272,single_copy,820\n"
                    "PL1272_0003,ParaLoss1272,multi_copy,610\n"
                )
                return csv_txt.encode("utf-8")
            return b"study_accession\trun_accession\tsample_accession\tscientific_name\tlibrary_name\tspecimen_voucher\tfastq_ftp\tfastq_bytes\n"

        self.assertEqual(ena_mod.resolve_ena_ncbi_taxid("Packera", http_get_fn=mock_http), "189364")
        self.assertEqual(ena_mod.resolve_ena_ncbi_taxid("Solidago", http_get_fn=mock_http), "59294")

        packera_rows = ena_mod.fetch_ena_tax_tree_sra_rows("Packera", ("PRJNA1046046",), http_get_fn=mock_http)
        self.assertTrue(any((r.get("bioproject") or r.get("study_accession")) == "PRJNA1046046" and "Packera dubia" in r.get("scientific_name", "") for r in packera_rows))
        self.assertTrue(any((r.get("bioproject") or r.get("study_accession")) == "PRJNA1046046" and "Packera glabella" in r.get("scientific_name", "") for r in packera_rows))

        solidago_rows = ena_mod.fetch_ena_tax_tree_sra_rows("Solidago", ("PRJNA999001",), http_get_fn=mock_http)
        self.assertTrue(any(r.get("scientific_name") == "Solidago altissima" for r in solidago_rows))

        with tempfile.TemporaryDirectory() as tmp:
            gen_dir = Path(tmp) / "genomics"
            with (
                patch.object(dl_pub, "_http_request", side_effect=mock_http),
                patch.object(dl_pub, "_ncbi_get", return_value=json.dumps({"esearchresult": {"idlist": []}}).encode("utf-8")),
            ):
                dl_pub.download_public_packera_sequences(
                    gen_dir, sources=("sra", "paraloss"), force=True, target_genus="Packera", target_family="Asteraceae",
                )

            sra_csv = gen_dir / "public_packera_sra_runs.csv"
            self.assertTrue(sra_csv.is_file())
            self.assertIn("PRJNA1046046", sra_csv.read_text(encoding="utf-8"))

            n_non_ast = ena_mod.fetch_compositae_paraloss1272_catalog(
                gen_dir, "Pinaceae", _atomic_csv_writer, http_get_fn=mock_http, force=True,
            )
            self.assertEqual(n_non_ast, 0)

            n_ast = ena_mod.fetch_compositae_paraloss1272_catalog(
                gen_dir, "Asteraceae", _atomic_csv_writer, http_get_fn=mock_http, force=True,
            )
            self.assertGreaterEqual(n_ast, 3)

            mod06 = importlib.import_module("scripts.data_prep.06_ingest_genetic_data")
            dosage = np.tile(np.array([[0.0], [0.5], [1.0], [0.0]], dtype=np.float64), (1, 3))
            depths = np.full((4, 3), 10.0, dtype=np.float64)
            keep_mask, diag = mod06._evaluate_hdplot_loci(
                dosage, dosage, depths, depths,
                ["PL1272_0001", "NOVEL_LOCUS_2", "PL1272_0003"],
                genomics_dir=gen_dir,
            )
            self.assertTrue(bool(np.all(keep_mask)))
            row_sc = diag.loc[diag["locus_id"] == "PL1272_0001"].iloc[0]
            row_nov = diag.loc[diag["locus_id"] == "NOVEL_LOCUS_2"].iloc[0]
            row_mc = diag.loc[diag["locus_id"] == "PL1272_0003"].iloc[0]
            self.assertEqual(int(row_sc["is_paraloss1272_single_copy"]), 1)
            self.assertAlmostEqual(float(row_sc["locus_paralogy_weight"]), 1.5, places=3)
            self.assertEqual(int(row_nov["is_paraloss1272_single_copy"]), 0)
            self.assertAlmostEqual(float(row_nov["locus_paralogy_weight"]), 1.0, places=3)
            self.assertAlmostEqual(float(row_mc["locus_paralogy_weight"]), 0.5, places=3)

    def test_taxonomic_backbone_and_natureserve_usda_fallback(self) -> None:
        def mock_col(url: str) -> dict:
            if "/nameusage/search" in url:
                return {
                    "last": True,
                    "result": [
                        {
                            "usage": {
                                "id": "wcvp-serp",
                                "status": "accepted",
                                "name": {"scientificName": "Packera serpenticola", "authorship": "Boufford et al.", "rank": "species", "publishedIn": "Harvard Pap. Bot. 19: 101", "publishedInYear": 2014},
                            },
                            "distributions": [{"tdwgCode": "NCU", "status": "native"}],
                        },
                        {
                            "usage": {
                                "id": "wcvp-aur",
                                "status": "accepted",
                                "name": {"scientificName": "Packera aurea", "authorship": "(L.) A. Love & D. Love", "rank": "species", "publishedIn": "Bot. Not. 128: 520", "publishedInYear": 1976},
                                "basionym": {"scientificName": "Senecio aureus", "authorship": "L."},
                            },
                            "homotypicSynonyms": ["Senecio aureus L."],
                            "distributions": [{"tdwgCode": "NCU", "status": "native"}],
                        },
                    ],
                }
            if "/taxon/wcvp-serp/info" in url:
                return {"synonyms": {"homotypic": [], "heterotypic": ["Packera anonyma var. serpenticola"]}, "distributions": [{"tdwgCode": "NCU", "status": "native"}]}
            return {}

        bhl_calls = []

        def mock_bhl(url: str) -> dict:
            bhl_calls.append(url)
            return {"Status": "ok", "Result": {"Titles": [{"PageUrl": "https://www.biodiversitylibrary.org/page/54321", "PublicationDate": "2014", "PublicationPlace": "Clay County, NC"}]}}

        def mock_ns(url: str) -> dict:
            return {
                "results": [
                    {
                        "uniqueId": "ELEMENT_GLOBAL.2.888001",
                        "scientificName": "Packera serpenticola",
                        "roundedGRank": "G2",
                        "nations": [{"nationCode": "US", "subnations": [{"subnationCode": "NC", "sRank": "S2", "native": True, "exotic": False}]}],
                    },
                    {
                        "uniqueId": "ELEMENT_GLOBAL.2.888002",
                        "scientificName": "Packera aurea",
                        "roundedGRank": "G5",
                        "nations": [{"nationCode": "CA", "subnations": [{"subnationCode": "ON", "sRank": "S5", "native": True, "exotic": False}]}],
                    },
                ]
            }

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            tax_dir = tmp_dir / "taxonomy"
            db_path = tmp_dir / "test_backbone.sqlite"
            get_db_connection(db_path=db_path, init_schema=True).close()

            # Test BHL unset fallback path cleanly skips external API calls
            profs_nokey = dl_backbone.fetch_col_powo_backbone("Packera", http_get_fn=mock_col)
            env_without_bhl = {k: v for k, v in os.environ.items() if k != "BHL_API_KEY"}
            with (
                patch.dict(os.environ, env_without_bhl, clear=True),
                patch.object(dl_backbone, "load_dotenv_defaults"),
            ):
                dl_backbone.enrich_bhl_protologues(profs_nokey, bhl_api_key=None, http_get_fn=mock_bhl)
            self.assertEqual(len(bhl_calls), 0)
            self.assertIsNone(profs_nokey["Packera serpenticola"]["protologue_url"])

            # Test full pipeline with BHL_API_KEY configured
            dl_backbone.download_taxonomic_backbone(
                taxonomy_dir=tax_dir, target_genus="Packera", db_path=db_path, force=True,
                bhl_api_key="test-bhl-key", col_http_fn=mock_col, bhl_http_fn=mock_bhl, natureserve_http_fn=mock_ns,
            )
            self.assertGreaterEqual(len(bhl_calls), 1)
            self.assertTrue((tax_dir / "public_taxon_backbone_profiles.csv").is_file())
            self.assertTrue((tax_dir / "public_natureserve_subnational_ranges.csv").is_file())
            self.assertTrue((tax_dir / "public_taxonomic_backbone_manifest.json").is_file())

            with sqlite3.connect(str(db_path)) as c2:
                prof_row = c2.execute(
                    "SELECT protologue_url, natureserve_global_rank FROM taxon_backbone_profiles WHERE species_standardized = 'Packera serpenticola'"
                ).fetchone()
                ns_row = c2.execute(
                    "SELECT country_code, subnation_code, s_rank, is_native FROM natureserve_subnational_ranges WHERE species_standardized = 'Packera serpenticola'"
                ).fetchone()
            self.assertEqual(prof_row, ("https://www.biodiversitylibrary.org/page/54321", "G2"))
            self.assertEqual(ns_row, ("US", "NC", "S2", 1))

            # Verify USDA county filter rescues Packera serpenticola in NC via NatureServe fallback
            cnty_geojson = _minimal_counties_geojson(tmp_dir / "us_counties.geojson")
            usda_cache = tmp_dir / "usda_cache"
            usda_cache.mkdir(parents=True, exist_ok=True)
            vouchers = pd.DataFrame([
                {"catalogNumber": "NCU-SERP-1", "species_standardized": "Packera serpenticola", "stateProvince": "North Carolina", "county": "Clay", "countryCode": "US", "decimalLatitude": 35.08, "decimalLongitude": -83.62, "determiner_tier": "Tier_2_Silver", "type_status": ""},
                {"catalogNumber": "NCU-SERP-2", "species_standardized": "Packera serpenticola", "stateProvince": "California", "county": "Alameda", "countryCode": "US", "decimalLatitude": 37.70, "decimalLongitude": -122.20, "determiner_tier": "Tier_2_Silver", "type_status": ""},
            ])
            filtered_df, stats = filter_vouchers_by_usda_county_range(
                vouchers, standardize_taxon_fn=lambda x: str(x), county_geojson_path=cnty_geojson,
                cache_dir=usda_cache, db_path=db_path, allow_network=False, drop_out_of_range=True,
            )
            self.assertEqual(list(filtered_df["catalogNumber"]), ["NCU-SERP-1"])
            self.assertEqual(list(filtered_df["usda_county_status"]), ["natureserve_state_in_range"])
            self.assertEqual(stats["natureserve_state_in_range"], 1)
            self.assertEqual(stats["dropped_out_of_range"], 1)

    def test_geochem_ds801_nwi_and_nwpl_sqlite_ingestion(self) -> None:
        def mock_http(url: str, *args, **kwargs) -> bytes:
            low = url.lower()
            if "sciencebase.gov" in low:
                return json.dumps({"files": [{"name": "DS801_A_horizon.csv", "url": "https://mock.local/DS801_A_horizon.csv"}]}).encode("utf-8")
            if "ds801_a_horizon.csv" in low:
                csv_txt = (
                    "site_id,latitude,longitude,horizon,ca_pct,mg_pct,ni_ppm,cr_ppm,c_inorg_pct,c_tot_pct\n"
                    "SERP1,35.0800,-83.6200,A,1.0,14.0,1800,2200,0.02,1.2\n"
                    "DOL1,36.4500,-81.4800,A,8.0,1.0,15,20,3.5,5.1\n"
                )
                return csv_txt.encode("utf-8")
            if "wetlandsmapservice" in low:
                if "-81.48" in url:
                    feat = {"attributes": {"ATTRIBUTE": "PEM1A", "WETLAND_TYPE": "Freshwater Emergent Wetland"}, "geometry": {"rings": [[[-81.4801, 36.4499], [-81.4799, 36.4499], [-81.4799, 36.4501], [-81.4801, 36.4501], [-81.4801, 36.4499]]]}}
                    return json.dumps({"features": [feat]}).encode("utf-8")
                return json.dumps({"features": []}).encode("utf-8")
            if "wetland-plants.sec.usace.army.mil" in low:
                return json.dumps([
                    {"scientificName": "Packera aurea", "region": "EMP", "indicatorStatus": "FACW"},
                    {"scientificName": "Packera serpenticola", "region": "EMP", "indicatorStatus": "FACU"},
                ]).encode("utf-8")
            return b"{}"

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            env_dir = tmp_dir / "environmental"
            db_path = tmp_dir / "geochem.sqlite"
            conn = get_db_connection(db_path=db_path, init_schema=True)
            conn.execute("INSERT INTO taxon_backbone_profiles (species_standardized, target_genus) VALUES ('Packera aurea', 'Packera')")
            conn.commit()
            conn.close()

            vouchers_csv = tmp_dir / "vouchers.csv"
            vouchers_csv.write_text("catalogNumber,decimalLatitude,decimalLongitude\nNCU-SERP,35.08000,-83.62000\nNCU-FEN,36.45000,-81.48000\n", encoding="utf-8")

            summary = dl_geochem.download_geochem_wetland_layers(
                env_dir=env_dir, vouchers_csv=vouchers_csv, db_path=db_path, genus="Packera", force=True, http_get_fn=mock_http,
            )
            self.assertEqual(summary["points_processed"], 2)
            self.assertTrue((env_dir / "usgs_ds801_soils.csv").is_file())
            self.assertTrue((env_dir / "nwi_point_cache.jsonl").is_file())
            self.assertTrue((env_dir / "usace_nwpl_ratings.csv").is_file())

            with sqlite3.connect(str(db_path)) as c2:
                serp = c2.execute("SELECT soil_ca_mg_molar_ratio, is_nwi_wetland FROM environmental_cache WHERE lat_lon_key = '35.08000_-83.62000'").fetchone()
                fen = c2.execute("SELECT soil_ca_mg_molar_ratio, is_nwi_wetland, nwi_wetland_dist_m FROM environmental_cache WHERE lat_lon_key = '36.45000_-81.48000'").fetchone()
                nwpl_row = c2.execute("SELECT nwpl_indicator_code, nwpl_wetland_score FROM taxon_backbone_profiles WHERE species_standardized = 'Packera aurea'").fetchone()

            self.assertLess(float(serp[0]), 0.7)
            self.assertEqual(int(serp[1]), 0)
            self.assertGreater(float(fen[0]), 2.5)
            self.assertEqual(int(fen[1]), 1)
            self.assertLessEqual(float(fen[2]), 30.0)
            self.assertIsNotNone(nwpl_row[1])
            self.assertAlmostEqual(float(nwpl_row[1]), 0.75, places=2)

    @unittest.skipIf(
        os.environ.get("PACKERA_FAST_AUDIT") == "1"
        or os.environ.get("PLANTMATICS_FAST_AUDIT") == "1"
        or shutil.which("Rscript") is None,
        "Skipping R hydro-edaphic and basionym guard test in fast audit mode or when Rscript is unavailable",
    )
    def test_r_hydro_edaphic_and_basionym_guards(self) -> None:
        r_util = PROJECT_ROOT / "scripts" / "analysis" / "edaphic_geochem_wetland_utils.R"
        self.assertTrue(r_util.is_file())

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            db_path = tmp_dir / "r_guards.sqlite"
            env_dir = tmp_dir / "environmental"
            env_dir.mkdir(parents=True, exist_ok=True)
            conn = get_db_connection(db_path=db_path, init_schema=True)
            conn.executemany(
                "INSERT INTO taxon_backbone_profiles (species_standardized, target_genus, basionym_name, homotypic_synonyms_json, heterotypic_synonyms_json, nwpl_indicator_code, nwpl_wetland_score) VALUES (?, 'Packera', ?, ?, ?, ?, ?)",
                [
                    ("Packera aurea", "Senecio aureus", '["Packera aurea"]', "[]", "FACW", 0.75),
                    ("Packera serpenticola", "Packera serpenticola", "[]", '["Senecio smallii var. serpentinus"]', "FACU", 0.25),
                    ("Packera anonyma", "Senecio anonymus", "[]", "[]", "UPL", 0.0),
                ],
            )
            conn.executemany(
                "INSERT INTO environmental_cache (lat_lon_key, latitude, longitude, topo_twi_30m, soil_ca_mg_molar_ratio, soil_ultramafic_ni_cr_index, soil_carbonate_c_pct, is_nwi_wetland, nwi_wetland_system, nwi_wetland_dist_m) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("36.45000_-81.48000", 36.45, -81.48, 12.0, 4.2, 0.1, 2.8, 1, "PEM1A", 5.0),
                    ("35.08000_-83.62000", 35.08, -83.62, 4.0, 0.05, 2.5, 0.02, 0, "UPL", 250.0),
                    ("36.00000_-79.00000", 36.00, -79.00, 7.0, 2.6, 0.2, 1.9, 0, "UPL", 400.0),
                ],
            )
            conn.commit()
            conn.close()

            r_expr = (
                "suppressPackageStartupMessages({library(DBI); library(RSQLite); library(jsonlite)}); "
                "args <- commandArgs(trailingOnly = TRUE); "
                "source(args[1]); "
                "con <- dbConnect(SQLite(), args[2]); "
                "on.exit(dbDisconnect(con), add = TRUE); "
                "df <- data.frame("
                "  catalogNumber = c('NCU-1', 'NCU-2', 'NCU-3'), "
                "  species_standardized = c('Packera aurea', 'Packera serpenticola', 'Packera anonyma'), "
                "  lat_lon_key = c('36.45000_-81.48000', '35.08000_-83.62000', '36.00000_-79.00000'), "
                "  taxonomic_treatment_subtype = c('Taxon_Resurrection_Candidate', 'Taxon_Resurrection_Candidate', 'Concordant_Assignment'), "
                "  resurrection_basionym = c('Senecio aureus', 'Senecio smallii var. serpentinus', NA_character_), "
                "  stringsAsFactors = FALSE"
                "); "
                "enr <- enrich_geochem_wetland_predictors(df, db_conn = con, env_dir = args[3]); "
                "res <- evaluate_step07_hydro_edaphic_and_basionym_guards(enr, db_conn = con); "
                "cat(toJSON(res[, c('catalogNumber', 'ultramafic_calcareous_gradient', 'wetland_seepage_index', 'hydro_edaphic_niche_regime', 'wetland_fidelity_concordance', 'basionym_synonymy_type', 'resurrection_basionym')], na = 'null'))"
            )
            proc = subprocess.run(
                ["Rscript", "-e", r_expr, str(r_util), str(db_path), str(env_dir)],
                capture_output=True, text=True, check=True, timeout=60,
            )
            rows = json.loads(proc.stdout.strip().splitlines()[-1])
            by_cat = {r["catalogNumber"]: r for r in rows}

            self.assertEqual(by_cat["NCU-1"]["hydro_edaphic_niche_regime"], "Palustrine_Wetland_Fen")
            self.assertEqual(by_cat["NCU-2"]["hydro_edaphic_niche_regime"], "Ultramafic_Serpentine_Barren")
            self.assertEqual(by_cat["NCU-3"]["hydro_edaphic_niche_regime"], "Calcareous_Carbonate_Glade")
            for cat in ("NCU-1", "NCU-2", "NCU-3"):
                self.assertIsNotNone(by_cat[cat]["ultramafic_calcareous_gradient"])
                self.assertIsNotNone(by_cat[cat]["wetland_seepage_index"])
                self.assertTrue(0.0 <= float(by_cat[cat]["wetland_fidelity_concordance"]) <= 1.0)

            self.assertEqual(by_cat["NCU-1"]["basionym_synonymy_type"], "Homotypic")
            self.assertIsNone(by_cat["NCU-1"]["resurrection_basionym"])
            self.assertEqual(by_cat["NCU-2"]["basionym_synonymy_type"], "Heterotypic")
            self.assertEqual(by_cat["NCU-2"]["resurrection_basionym"], "Senecio smallii var. serpentinus")


if __name__ == "__main__":
    unittest.main()
