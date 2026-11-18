#!/usr/bin/env python3
"""Offline unit tests for M1B 4-layer specialist discovery, infraspecific trinomial policy,
multi-portal iDigBio/Symbiota harvesting & deduplication, and specimen sex/shoot-phase parsing."""

from __future__ import annotations

import importlib
import json
import os
import re
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

from scripts.core import harvester as harv_mod
from scripts.core.config import resolve_infraspecific_rank_policy
from scripts.core.harvester import (
    VoucherHarvester,
    export_curated_table,
    parse_specimen_sex_and_shoot_phase,
    standardize_taxon_name,
)
from scripts.core.idigbio_symbiota_harvester import fetch_and_deduplicate_portal_records
from scripts.core.specialist_discovery_utils import (
    discover_and_sync_specialists,
    load_active_specialist_patterns,
)


class TestM1BSpecialistsAndPortals(unittest.TestCase):
    """Offline tests for specialist discovery, trinomial policy, portals, and phenology metadata."""

    def test_1_four_layer_specialist_discovery_and_guardrails(self) -> None:
        """Verify 4-layer OpenAlex/WCVP/Wikidata/Bionomia specialist sync, pre-1950 exclusion, and flat CSV."""
        def _fake_openalex(url: str) -> dict:
            if "Quercus" in url:
                return {"results": [
                    {"id": "W1", "title": "Systematics and taxonomy monograph of Quercus", "authorships": [
                        {"author": {"display_name": "K. C. Nixon", "orcid": "0000-0001-1111-2222", "works_count": 14}},
                        {"author": {"display_name": "Smith", "works_count": 50}},
                    ]},
                    {"id": "W2", "title": "Taxonomy of Quercus sect. Lobatae", "authorships": [{"author": {"display_name": "K. C. Nixon", "orcid": "0000-0001-1111-2222", "works_count": 14}}]},
                ]}
            if "Solidago" in url:
                return {"results": [{"id": "W3", "title": "Revision and taxonomy of Solidago", "authorships": [{"author": {"display_name": "J. C. Semple", "orcid": "0000-0002-3333-4444", "works_count": 28}}]}]}
            return {"results": []}

        def _fake_wikidata(url: str) -> dict:
            return {"results": {"bindings": [
                {"p428": {"value": "K.C.Nixon"}, "item": {"value": "http://www.wikidata.org/entity/Q101"}, "itemLabel": {"value": "K. C. Nixon"}, "birth": {"value": "1953-01-01T00:00:00Z"}},
                {"p428": {"value": "J.C.Semple"}, "item": {"value": "http://www.wikidata.org/entity/Q102"}, "itemLabel": {"value": "J. C. Semple"}, "birth": {"value": "1947-01-01T00:00:00Z"}},
                {"p428": {"value": "L."}, "item": {"value": "http://www.wikidata.org/entity/Q1043"}, "itemLabel": {"value": "Carl Linnaeus"}, "birth": {"value": "1707-05-23T00:00:00Z"}, "death": {"value": "1778-01-10T00:00:00Z"}},
                {"p428": {"value": "Michx."}, "item": {"value": "http://www.wikidata.org/entity/Q318"}, "itemLabel": {"value": "André Michaux"}, "death": {"value": "1802-11-16T00:00:00Z"}},
                {"p428": {"value": "Torr."}, "item": {"value": "http://www.wikidata.org/entity/Q380"}, "itemLabel": {"value": "John Torrey"}, "death": {"value": "1873-03-10T00:00:00Z"}},
            ]}}

        def _fake_bionomia(url: str) -> dict:
            if "Quercus" not in url:
                return {"records": []}
            return {"records": [
                {"canonical_name": "K. C. Nixon", "institutionCodes": ["BH", "NY", "MO", "US"], "institution_count": 4, "bionomia_det_count": 65, "variants": ["K. C. Nixon", "Nixon, K. C."]},
                {"canonical_name": "A. B. Regional", "institutionCodes": ["NCU", "DUKE"], "institution_count": 2, "bionomia_det_count": 18, "variants": ["A. B. Regional"]},
            ]}

        wcvp_quercus = {
            "q1": {"target_genus": "Quercus", "species_standardized": "Quercus welshii", "publication_year": 1985, "accepted_authorship": "K.C.Nixon"},
            "q2": {"target_genus": "Quercus", "species_standardized": "Quercus alba", "publication_year": 1753, "accepted_authorship": "L."},
            "q3": {"target_genus": "Quercus", "species_standardized": "Quercus michauxii", "publication_year": 1803, "accepted_authorship": "Michx."},
            "q4": {"target_genus": "Quercus", "species_standardized": "Quercus kelloggii", "publication_year": 1857, "accepted_authorship": "Torr."},
        }
        wcvp_solidago = {
            "s1": {"target_genus": "Solidago", "species_standardized": "Solidago kralii", "publication_year": 2003, "accepted_authorship": "J.C.Semple"},
        }

        with tempfile.TemporaryDirectory() as td:
            tmp_p = Path(td)
            tax_dir = tmp_p / "taxonomy"
            db_p = tmp_p / "spec.sqlite"

            q_rows = discover_and_sync_specialists(
                target_genus="Quercus", target_family="Fagaceae", taxonomy_dir=tax_dir, db_path=db_p,
                wcvp_profiles=wcvp_quercus, openalex_http_fn=_fake_openalex,
                wikidata_http_fn=_fake_wikidata, bionomia_http_fn=_fake_bionomia,
            )
            q_by_name = {r["canonical_name"]: r for r in q_rows}
            self.assertIn("K. C. Nixon", q_by_name)
            self.assertEqual(q_by_name["K. C. Nixon"]["tier_cap"], "Tier_1_Gold")
            self.assertIn("A. B. Regional", q_by_name)
            self.assertEqual(q_by_name["A. B. Regional"]["tier_cap"], "Tier_2_Silver")
            for blocked in ("Smith", "L.", "Carl Linnaeus", "Michx.", "André Michaux", "Torr.", "John Torrey"):
                self.assertNotIn(blocked, q_by_name)

            s_rows = discover_and_sync_specialists(
                target_genus="Solidago", target_family="Asteraceae", taxonomy_dir=tax_dir, db_path=db_p,
                wcvp_profiles=wcvp_solidago, openalex_http_fn=_fake_openalex,
                wikidata_http_fn=_fake_wikidata, bionomia_http_fn=_fake_bionomia,
            )
            s_by_name = {r["canonical_name"]: r for r in s_rows}
            self.assertIn("J. C. Semple", s_by_name)
            self.assertEqual(s_by_name["J. C. Semple"]["tier_cap"], "Tier_1_Gold")

            gold_q = load_active_specialist_patterns("Quercus", "Fagaceae", db_path=db_p, taxonomy_dir=tax_dir)
            self.assertTrue(any(re.search(p, "K. C. Nixon", re.I) for p in gold_q))
            self.assertFalse(any(re.search(p, "A. B. Regional", re.I) for p in gold_q))
            self.assertFalse(any(re.search(p, "Smith", re.I) for p in gold_q))

            self.assertTrue((tax_dir / "public_taxonomic_specialists.csv").is_file())
            self.assertEqual([p.name for p in tax_dir.iterdir() if p.is_dir()], [])

            zero_rows = discover_and_sync_specialists(
                target_genus="Raregenus", target_family="Asteraceae", taxonomy_dir=tax_dir, db_path=db_p,
                wcvp_profiles={}, openalex_http_fn=_fake_openalex,
                wikidata_http_fn=_fake_wikidata, bionomia_http_fn=_fake_bionomia,
            )
            self.assertEqual(zero_rows, [])
            self.assertEqual(load_active_specialist_patterns("Raregenus", "Asteraceae", db_path=db_p, taxonomy_dir=tax_dir), [])

    def test_2_infraspecific_rank_policy_and_sqlite_persistence(self) -> None:
        """Verify preserve_accepted_trinomials vs collapse_to_species, autonym collapse, and f./cv. rejection."""
        cases_preserve = (
            ("Pinus ponderosa var. scopulorum", "Pinus", "Pinus ponderosa var. scopulorum"),
            ("Pinus ponderosa Douglas ex C.Lawson var. scopulorum Engelm.", "Pinus", "Pinus ponderosa var. scopulorum"),
            ("Eriogonum umbellatum Torr. var. majus Hook.", "Eriogonum", "Eriogonum umbellatum var. majus"),
            ("Packera streptanthifolia ssp. borealis", "Packera", "Packera streptanthifolia subsp. borealis"),
            ("Quercus alba var. alba", "Quercus", "Quercus alba"),
            ("Quercus alba f. latiloba", "Quercus", "Quercus alba"),
            ("Quercus alba cv. Fastigiata", "Quercus", "Quercus alba"),
        )
        for raw_n, gen_n, exp_n in cases_preserve:
            self.assertEqual(standardize_taxon_name(raw_n, target_genus=gen_n, infraspecific_policy="preserve_accepted_trinomials"), exp_n)
        with patch.dict(os.environ, {"PLANTMATICS_INFRASPECIFIC_POLICY": "preserve_accepted_trinomials"}):
            self.assertEqual(resolve_infraspecific_rank_policy(), "preserve_accepted_trinomials")
            self.assertEqual(standardize_taxon_name("Pinus ponderosa var. scopulorum", target_genus="Pinus"), "Pinus ponderosa var. scopulorum")
        self.assertEqual(standardize_taxon_name("Pinus ponderosa var. scopulorum", target_genus="Pinus", infraspecific_policy="collapse_to_species"), "Pinus ponderosa")
        self.assertEqual(standardize_taxon_name("Eriogonum umbellatum var. majus", target_genus="Eriogonum", infraspecific_policy="collapse_to_species"), "Eriogonum umbellatum")

        mock_recs = [
            {
                "key": 701, "catalogNumber": "NCU_Q1", "institutionCode": "NCU",
                "scientificName": "Quercus alba var. alba", "stateProvince": "North Carolina",
                "decimalLatitude": 35.9, "decimalLongitude": -79.05, "year": 2020, "month": 5, "day": 1,
                "recordedBy": "Collector A", "recordNumber": "101",
                "media": [{"identifier": "https://example.org/q1.jpg", "type": "StillImage", "format": "image/jpeg"}],
            },
            {
                "key": 702, "catalogNumber": "NCU_Q2", "institutionCode": "NCU",
                "scientificName": "Quercus sinuata var. sinuata", "stateProvince": "North Carolina",
                "decimalLatitude": 35.91, "decimalLongitude": -79.06, "year": 2020, "month": 5, "day": 2,
                "recordedBy": "Collector B", "recordNumber": "102",
                "media": [{"identifier": "https://example.org/q2.jpg", "type": "StillImage", "format": "image/jpeg"}],
            },
            {
                "key": 703, "catalogNumber": "NCU_Q3", "institutionCode": "NCU",
                "scientificName": "Quercus alba subsp. montana", "stateProvince": "North Carolina",
                "decimalLatitude": 35.92, "decimalLongitude": -79.07, "year": 2020, "month": 5, "day": 3,
                "recordedBy": "Collector C", "recordNumber": "103",
                "media": [{"identifier": "https://example.org/q3.jpg", "type": "StillImage", "format": "image/jpeg"}],
            },
        ]
        with patch.dict(os.environ, {"PLANTMATICS_TARGET_GENUS": "Quercus"}), patch(
            "scripts.core.harvester.occ.search",
            side_effect=lambda **kw: {"results": [] if ("typeStatus" in kw or kw.get("country") == "CA") else mock_recs, "count": 3},
        ):
            harv = VoucherHarvester(taxa=["Quercus alba", "Quercus sinuata"], expand_synonyms=False, infraspecific_policy="preserve_accepted_trinomials")
            df_q = harv.harvest().set_index("catalogNumber")
            self.assertEqual(df_q.loc["NCU_Q1", "species_standardized"], "Quercus alba")
            self.assertTrue(pd.isna(df_q.loc["NCU_Q1", "infraspecific_rank"]))
            self.assertTrue(pd.isna(df_q.loc["NCU_Q1", "infraspecific_epithet"]))
            self.assertEqual(df_q.loc["NCU_Q3", "species_standardized"], "Quercus alba subsp. montana")
            self.assertEqual(df_q.loc["NCU_Q3", "infraspecific_rank"], "subsp.")
            self.assertEqual(df_q.loc["NCU_Q3", "infraspecific_epithet"], "montana")

            with tempfile.TemporaryDirectory() as td:
                csv_p, db_p = Path(td) / "curated.csv", Path(td) / "quercus.sqlite"
                export_curated_table(df_q.reset_index(), csv_p, db_path=db_p, infraspecific_policy="preserve_accepted_trinomials")
                csv_df = pd.read_csv(csv_p)
                for col in ("portal_source", "specimen_sex", "seasonal_shoot_phase", "infraspecific_rank", "infraspecific_epithet"):
                    self.assertIn(col, csv_df.columns)
                with sqlite3.connect(str(db_p)) as conn:
                    r3 = conn.execute("SELECT portal_source, infraspecific_rank, infraspecific_epithet FROM vouchers WHERE catalogNumber='NCU_Q3'").fetchone()
                    r1 = conn.execute("SELECT infraspecific_rank, infraspecific_epithet FROM vouchers WHERE catalogNumber='NCU_Q1'").fetchone()
                self.assertEqual(r3, ("GBIF", "subsp.", "montana"))
                self.assertEqual(r1, (None, None))

    def test_3_portal_dedup_backfill_sex_shoot_phase_and_cli(self) -> None:
        """Verify iDigBio/Symbiota deduplication, media backfilling, sex/shoot-phase regexes, and CLI wiring."""
        # 3A. Sex & seasonal shoot phase parser
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"sex": "male", "occurrenceRemarks": "Salix staminate catkins"}), ("staminate", None))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"reproductiveCondition": "pistillate flowers; female plant"}), ("pistillate", None))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"occurrenceRemarks": "perfect flowers, hermaphrodite"}), ("hermaphrodite", None))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"reproductiveCondition": "cleistogamous flowers"}), (None, "cleistogamous"))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"fieldNotes": "Viola chasmogamous corollas"}), (None, "chasmogamous"))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"dynamicProperties": '{"phase": "vernal culms"}'}), (None, "vernal"))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"occurrenceRemarks": "Dichanthelium autumnal phase"}), (None, "autumnal"))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"fieldNotes": "Osmunda sterile frond"}), (None, "sterile_frond"))
        self.assertEqual(parse_specimen_sex_and_shoot_phase({"reproductiveCondition": "fertile frond with sporophylls"}), (None, "fertile_frond"))

        # 3B. Multi-portal deduplication & backfill
        gbif_in = [
            {
                "key": 1, "catalogNumber": "NCU00012345", "institutionCode": "NCU", "occurrenceID": "urn:ncu:12345",
                "scientificName": "Packera dubia", "media": [{"identifier": "https://gbif.org/ncu1.jpg", "type": "StillImage"}],
            },
            {
                "key": 2, "catalogNumber": "DUKE000222", "institutionCode": "DUKE", "occurrenceID": "urn:duke:222",
                "scientificName": "Packera dubia", "media": [],
            },
        ]

        def _fake_idigbio(url: str) -> dict:
            if "view/mediarecords/m-uuid-9" in url:
                return {"data": {"ac:accessURI": "https://sernec.org/images/sernec9.jpg"}}
            return {
                "itemCount": 3,
                "items": [
                    {
                        "uuid": "idig-dup",
                        "data": {"dwc:catalogNumber": "NCU00012345", "dwc:institutionCode": "NCU", "dwc:scientificName": "Packera dubia", "ac:accessURI": "https://idigbio.org/dup.jpg"},
                    },
                    {
                        "uuid": "idig-backfill",
                        "data": {"dwc:catalogNumber": "DUKE000222", "dwc:institutionCode": "DUKE", "dwc:scientificName": "Packera dubia", "ac:accessURI": "https://sernec.org/duke_backfill.jpg"},
                    },
                    {
                        "uuid": "idig-new",
                        "indexTerms": {"geopoint": {"lat": 35.5, "lon": -79.2}, "mediarecords": ["m-uuid-9"]},
                        "data": {
                            "dwc:catalogNumber": "SERNEC_9001", "dwc:institutionCode": "CLEMS", "dwc:scientificName": "Packera dubia",
                            "dwc:stateProvince": "South Carolina", "dwc:countryCode": "US", "dwc:eventDate": "2018-05-12",
                            "dwc:occurrenceRemarks": "cleistogamous",
                        },
                    },
                ],
            }

        def _fake_symbiota(url: str) -> dict:
            return {
                "count": 1,
                "results": [
                    {
                        "occid": 8801, "catalogNumber": "SEINET_55", "institutionCode": "ARIZ", "occurrenceID": "urn:seinet:55",
                        "scientificName": "Packera dubia", "decimalLatitude": 35.8, "decimalLongitude": -82.5,
                        "stateProvince": "North Carolina", "countryCode": "US", "eventDate": "2019-04-22",
                        "images": [{"originalUrl": "https://seinet.org/orig/55.jpg"}],
                        "sex": "staminate", "reproductiveCondition": "vernal phase",
                    }
                ],
            }

        merged, stats = fetch_and_deduplicate_portal_records(
            taxa_queries=["Packera dubia"], gbif_records=gbif_in,
            portal_sources=["gbif", "idigbio", "symbiota"],
            symbiota_api_urls=["https://seinet.example.org/api/occurrences"],
            idigbio_http_fn=_fake_idigbio, symbiota_http_fn=_fake_symbiota,
        )
        self.assertEqual(stats, {
            "gbif_input": 2, "idigbio_fetched": 3, "symbiota_fetched": 1,
            "deduplicated_skipped": 2, "media_backfilled": 1, "net_new_added": 2, "total_merged": 4,
        })
        by_cat = {r["catalogNumber"]: r for r in merged}
        self.assertEqual(by_cat["NCU00012345"]["portal_source"], "GBIF")
        self.assertEqual(by_cat["DUKE000222"]["image_url"], "https://sernec.org/duke_backfill.jpg")
        self.assertEqual(by_cat["SERNEC_9001"]["portal_source"], "iDigBio")
        self.assertEqual(by_cat["SERNEC_9001"]["image_url"], "https://sernec.org/images/sernec9.jpg")
        self.assertEqual(by_cat["SEINET_55"]["portal_source"], "Symbiota")

        # 3C. End-to-end VoucherHarvester + pipeline CLI wiring
        pipe01 = importlib.import_module("scripts.pipeline.01_voucher_harvester")
        cli_args = pipe01.build_cli_parser().parse_args([
            "--taxon", "Packera dubia", "--portal-sources", "gbif", "idigbio", "symbiota",
            "--symbiota-api-url", "https://seinet.example.org/api/occurrences",
            "--infraspecific-policy", "preserve_accepted_trinomials",
        ])
        self.assertEqual(cli_args.portal_sources, ["gbif", "idigbio", "symbiota"])
        self.assertEqual(cli_args.symbiota_api_urls, ["https://seinet.example.org/api/occurrences"])
        self.assertEqual(cli_args.infraspecific_policy, "preserve_accepted_trinomials")

        with patch("scripts.core.harvester.occ.search", side_effect=lambda **kw: {"results": [] if ("typeStatus" in kw or kw.get("country") == "CA") else [
            {
                "key": 910, "catalogNumber": "NCU00012345", "institutionCode": "NCU", "occurrenceID": "urn:ncu:12345",
                "scientificName": "Packera dubia", "stateProvince": "North Carolina", "decimalLatitude": 35.9, "decimalLongitude": -79.05,
                "year": 2021, "month": 5, "media": [{"identifier": "https://gbif.org/ncu1.jpg", "type": "StillImage"}],
                "sex": "pistillate",
            }
        ], "count": 1}):
            harv_mp = VoucherHarvester(
                taxa=["Packera dubia"], expand_synonyms=False,
                portal_sources=["gbif", "idigbio", "symbiota"],
                symbiota_api_urls=["https://seinet.example.org/api/occurrences"],
                idigbio_http_fn=_fake_idigbio, symbiota_http_fn=_fake_symbiota,
            )
            df_mp = harv_mp.harvest().set_index("catalogNumber")
            self.assertEqual(set(df_mp.index), {"NCU00012345", "SERNEC_9001", "SEINET_55"})
            self.assertEqual(df_mp.loc["NCU00012345", "portal_source"], "GBIF")
            self.assertEqual(df_mp.loc["NCU00012345", "specimen_sex"], "pistillate")
            self.assertEqual(df_mp.loc["SERNEC_9001", "portal_source"], "iDigBio")
            self.assertEqual(df_mp.loc["SERNEC_9001", "seasonal_shoot_phase"], "cleistogamous")
            self.assertEqual(df_mp.loc["SEINET_55", "portal_source"], "Symbiota")
            self.assertEqual(df_mp.loc["SEINET_55", "specimen_sex"], "staminate")
            self.assertEqual(df_mp.loc["SEINET_55", "seasonal_shoot_phase"], "vernal")


if __name__ == "__main__":
    unittest.main()
